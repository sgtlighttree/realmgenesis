#!/usr/bin/env python3
"""Local map editor server (stdlib only). Run: tools/map/py tools/map/serve.py"""
import argparse, ipaddress, json, os, re, shutil, socket, subprocess, sys, tempfile, threading, time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

HOST, PORT = "127.0.0.1", 8765
ROOT = Path(__file__).resolve().parents[2]
EDITOR = Path(__file__).resolve().parent / "editor.html"
MAPS = ROOT / "library" / "maps"
STATE = MAPS / "map-state.json"
HIST = MAPS / "history"
QUAR = HIST / "quarantine"
KEEP = 50
MAX_BODY = 32 * 1024 * 1024
LOCK = threading.Lock()
COASTGEN = Path(os.environ.get("COASTGEN_PATH") or Path(__file__).resolve().parent / "terrain" / "coastgen.py")  # env override: test stubs
REGEN_TIMEOUT = 180
GROUP_RE = re.compile(r"^[a-z0-9][a-z0-9-]{0,63}$")
REGEN = {"running": False, "group": None, "started": None, "last_summary": None, "stderr_tail": ""}
REGEN_LOCK = threading.Lock()


GEN_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,63}\.json$")


def gen_dir():
    return MAPS / "generated"


def gen_rev(p):
    st = p.stat()
    return f"{st.st_mtime_ns}-{st.st_size}"


def gen_list():
    """Read-only listing of generated coast files: [{name, size, rev}]."""
    d, out = gen_dir(), []
    if d.is_dir():
        for p in sorted(d.iterdir()):
            if GEN_NAME.match(p.name) and p.is_file() and not p.is_symlink():
                out.append({"name": p.name[:-5], "size": p.stat().st_size, "rev": gen_rev(p)})
    return out


def gen_file(name):
    """Resolve /generated/<name>.json safely (strict name, no symlinks, must stay inside generated/), else None."""
    if not GEN_NAME.match(name):
        return None
    d = gen_dir()
    p = d / name
    try:
        if p.is_symlink() or not p.is_file() or p.resolve().parent != d.resolve():
            return None
    except OSError:
        return None
    return p


def rev_of():
    """Disk revision of the state file (mtime_ns-size), or None if it does not exist."""
    try:
        st = STATE.stat()
    except OSError:
        return None
    return f"{st.st_mtime_ns}-{st.st_size}"


def rotate_history():
    files = sorted(p for p in HIST.glob("map-state-*.json") if p.is_file())
    QUAR.mkdir(parents=True, exist_ok=True)
    for old in files[:-KEEP] if len(files) > KEEP else []:
        dest = QUAR / old.name
        n = 1
        while dest.exists():
            dest = QUAR / f"{old.stem}.{n}{old.suffix}"
            n += 1
        shutil.move(str(old), str(dest))


def write_state(text):
    MAPS.mkdir(parents=True, exist_ok=True)
    HIST.mkdir(parents=True, exist_ok=True)
    if STATE.exists():
        if STATE.read_text(encoding="utf-8") == text:
            return
        t = time.time()
        stamp = time.strftime("%Y%m%d-%H%M%S", time.localtime(t)) + f"-{int(t * 1000) % 1000:03d}"
        shutil.copy2(STATE, HIST / f"map-state-{stamp}.json")
        rotate_history()
    fd, tmp = tempfile.mkstemp(dir=MAPS, prefix=".map-state-", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, STATE)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def regen_groups():
    """Whitelist for /regen: landmass tags (tags on `land` features) present in the state file."""
    try:
        d = json.loads(STATE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return set()
    out = set()
    for f in d.get("features", []):
        if isinstance(f, dict) and f.get("type") == "land" and isinstance(f.get("groups"), list):
            out.update(t for t in f["groups"] if isinstance(t, str) and GROUP_RE.match(t))
    return out


def regen_run(group):
    """Worker thread: fixed argv, no shell, no request-supplied args beyond the whitelisted group."""
    summary, err = None, ""
    try:
        r = subprocess.run([sys.executable, str(COASTGEN), "--group", group, "--state", str(STATE)],
                           capture_output=True, text=True, timeout=REGEN_TIMEOUT, stdin=subprocess.DEVNULL, cwd=str(ROOT))
        err = (r.stderr or "")[-2000:]
        for line in reversed((r.stdout or "").splitlines()):
            line = line.strip()
            if line.startswith("{"):
                try:
                    summary = json.loads(line)
                except ValueError:
                    continue
                break
        if summary is None:
            summary = {"group": group, "error": f"no JSON summary (exit {r.returncode})"}
    except subprocess.TimeoutExpired as e:
        raw = e.stderr or b""
        err = (raw.decode("utf-8", "replace") if isinstance(raw, bytes) else raw)[-2000:]
        summary = {"group": group, "error": f"timeout after {REGEN_TIMEOUT}s"}
    except OSError as e:
        summary = {"group": group, "error": f"spawn failed: {e}"}
    with REGEN_LOCK:
        REGEN.update(running=False, last_summary=summary, stderr_tail=err)


class Handler(BaseHTTPRequestHandler):
    server_version = "MapEditor/1"

    def log_message(self, fmt, *args):
        if "/state" in (args[0] if args else "") and "POST" not in (args[0] if args else ""):
            return
        sys.stderr.write("%s %s\n" % (self.address_string(), fmt % args))

    def _host_ok(self):
        # DNS-rebinding guard: only IP literals on loopback/private/Tailscale ranges,
        # localhost, and *.local / *.ts.net names, on our port.
        name, _, port = self.headers.get("Host", "").rpartition(":")
        if port != str(PORT):
            return False
        name = name.strip("[]").lower()
        if name == "localhost" or name.endswith((".local", ".ts.net")):
            return True
        try:
            ip = ipaddress.ip_address(name)
        except ValueError:
            return False
        return ip.is_loopback or ip.is_private or ip in ipaddress.ip_network("100.64.0.0/10")

    def _send(self, code, body=b"", ctype="text/plain; charset=utf-8", rev=None):
        self.send_response(code)
        if rev:
            self.send_header("X-Map-Rev", rev)
            self.send_header("Access-Control-Expose-Headers", "X-Map-Rev")
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def do_GET(self):
        if not self._host_ok():
            return self._send(403, b"bad host")
        path = self.path.split("?", 1)[0]
        if path in ("/", "/index.html", "/editor.html"):
            return self._send(200, EDITOR.read_bytes(), "text/html; charset=utf-8")
        if path == "/state":
            if not STATE.exists():
                return self._send(404, b"no state file")
            rev = rev_of()
            return self._send(200, STATE.read_bytes(), "application/json; charset=utf-8", rev)
        if path in ("/generated", "/generated/"):
            body = json.dumps({"files": gen_list()}).encode()
            return self._send(200, body, "application/json; charset=utf-8")
        if path.startswith("/generated/"):
            p = gen_file(path[len("/generated/"):])
            if p is None:
                return self._send(404, b"not found")
            try:
                return self._send(200, p.read_bytes(), "application/json; charset=utf-8", gen_rev(p))
            except OSError:
                return self._send(404, b"not found")
        if path == "/regen/status":
            with REGEN_LOCK:
                body = json.dumps(REGEN).encode()
            return self._send(200, body, "application/json; charset=utf-8")
        if path == "/favicon.ico":
            return self._send(204)
        self._send(404, b"not found")

    do_HEAD = do_GET

    def do_POST(self):
        if not self._host_ok():
            return self._send(403, b"bad host")
        pth = self.path.split("?", 1)[0]
        if pth == "/regen":
            return self._regen()
        if pth != "/state":
            return self._send(404, b"not found")
        try:
            n = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            return self._send(400, b"bad length")
        if n <= 0 or n > MAX_BODY:
            return self._send(413, b"bad size")
        try:
            data = json.loads(self.rfile.read(n).decode("utf-8"))
            if not isinstance(data, dict) or not isinstance(data.get("features"), list):
                raise ValueError("expected object with features[]")
        except (ValueError, UnicodeDecodeError) as e:
            return self._send(400, f"invalid state: {e}".encode())
        text = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
        base = self.headers.get("X-Map-Rev")
        try:
            with LOCK:
                cur = rev_of()
                if cur is not None:
                    if base is None:
                        return self._send(428, b"X-Map-Rev header required (reload the page)", rev=cur)
                    if base != cur:
                        return self._send(409, b"state file changed on disk", rev=cur)
                write_state(text)
                new = rev_of()
        except OSError as e:
            return self._send(500, f"write failed: {e}".encode())
        self._send(200, b'{"ok":true}', "application/json", new)


def _regen(self):
    if not self.headers.get("Content-Type", "").lower().startswith("application/json"):
        return self._send(415, b"application/json required")
    try:
        n = int(self.headers.get("Content-Length", "0"))
        if n <= 0 or n > 4096:
            raise ValueError("bad size")
        data = json.loads(self.rfile.read(n).decode("utf-8"))
        group, rev = data.get("group"), data.get("rev")
        if not isinstance(group, str) or not isinstance(rev, str):
            raise ValueError("expected {group: str, rev: str}")
    except (ValueError, UnicodeDecodeError, AttributeError) as e:
        return self._send(400, f"bad request: {e}".encode())
    if group not in regen_groups():
        return self._send(400, b"unknown group (not a landmass tag in the state file)")
    if not COASTGEN.is_file():
        return self._send(501, b"coastgen not found")
    with REGEN_LOCK:
        if REGEN["running"]:
            return self._send(409, b"a regeneration is already running")
        cur = rev_of()
        if cur is None or rev != cur:
            return self._send(409, b"stale rev: reload the map, save, then retry", rev=cur)
        REGEN.update(running=True, group=group, started=time.time(), last_summary=None, stderr_tail="")
    threading.Thread(target=regen_run, args=(group,), daemon=True).start()
    self._send(202, json.dumps({"ok": True, "group": group}).encode(), "application/json")


Handler._regen = _regen


def lan_ip():
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
            s.connect(("192.0.2.1", 9))  # no packet sent; picks the outbound interface
            return s.getsockname()[0]
    except OSError:
        return None


def main():
    global HOST, PORT, STATE, MAPS, HIST, QUAR
    ap = argparse.ArgumentParser(description="Local map editor server")
    ap.add_argument("--host", default=HOST,
                    help="bind address; 0.0.0.0 to reach it from phone/LAN/Tailscale")
    ap.add_argument("--port", type=int, default=PORT)
    ap.add_argument("--state", type=Path, default=STATE, help="state file (history/ goes next to it)")
    a = ap.parse_args()
    try:
        import numpy, PIL  # noqa: F401  (not used here; /regen runs coastgen with this interpreter)
    except ImportError:
        print("warning: this interpreter lacks numpy/Pillow, so /regen will fail. Start with: tools/map/py tools/map/serve.py --host 0.0.0.0", file=sys.stderr)
    HOST, PORT, STATE = a.host, a.port, a.state.resolve()
    MAPS = STATE.parent
    HIST = MAPS / "history"
    QUAR = HIST / "quarantine"
    try:
        srv = ThreadingHTTPServer((HOST, PORT), Handler)
    except OSError as e:
        print(f"Port {PORT} on {HOST} is already in use ({e}). "
              f"The map server is probably already running: open http://{HOST}:{PORT}. "
              "Not killing anything.")
        sys.exit(1)
    print(f"Map editor: http://127.0.0.1:{PORT}  (state: {STATE})  Ctrl+C to stop")
    if HOST not in ("127.0.0.1", "localhost"):
        ip = lan_ip()
        print(f"  from phone/LAN: http://{ip or '<this-mac-ip>'}:{PORT}  or  http://{socket.gethostname().split('.')[0]}.local:{PORT}")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()


if __name__ == "__main__":
    main()
