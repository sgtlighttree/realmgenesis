#!/bin/sh
# (Re)build tools/map/.venv from requirements.txt. Uses Homebrew's stable python@3.14 path (Matt, 2026-10-11: stay on 3.14; 3.15 lacks some wheels) (not the `python3` link,
# which moves when brew changes the default); override with PY=/path/to/python3.x tools/map/setup-venv.sh
set -e
d=$(cd "$(dirname "$0")" && pwd)
PY=${PY:-/opt/homebrew/opt/python@3.14/bin/python3.14}
[ -x "$PY" ] || { echo "no $PY; set PY=/path/to/python3.14" >&2; exit 1; }   # never bare python3 (can be an old 3.10)
if [ -d "$d/.venv" ]; then mkdir -p "$d/.quarantine" && mv "$d/.venv" "$d/.quarantine/venv-$(date +%Y%m%d-%H%M%S)"; fi
"$PY" -m venv "$d/.venv"
"$d/.venv/bin/pip" install -q --only-binary=:all: -r "$d/requirements.txt"
"$d/.venv/bin/python" -c "import sys, numpy, PIL; print('ok', sys.version.split()[0], 'numpy', numpy.__version__, 'Pillow', PIL.__version__)"
