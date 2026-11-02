//! The world file container: a versioned header plus named binary chunks.
//!
//! This crate knows nothing about tectonics or climate. It only stores bytes. The pipeline
//! maps a `World` to and from chunks, so the file layout doesn't couple to the solvers
//! (ARCHITECTURE §3).
//!
//! Layout (all integers little-endian):
//! ```text
//! magic        8 bytes  "RGWORLD\0"
//! format       u32      container format version (FORMAT_VERSION)
//! header_len   u32
//! header       JSON     free-form metadata (params, generator versions)
//! chunk*:
//!   name_len   u16
//!   name       utf-8
//!   kind       u8       ChunkKind
//!   byte_len   u64
//!   bytes
//! ```

use std::io::{self, Read, Write};

pub const MAGIC: &[u8; 8] = b"RGWORLD\0";
/// Bump when the container layout changes (not when chunk contents change; those are
/// versioned in the header).
pub const FORMAT_VERSION: u32 = 1;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum ChunkKind {
    Bytes = 0,
    F32 = 1,
    U16 = 2,
    U8 = 3,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Chunk {
    pub name: String,
    pub kind: ChunkKind,
    pub bytes: Vec<u8>,
}

impl Chunk {
    pub fn f32s(name: &str, data: &[f32]) -> Self {
        Self { name: name.into(), kind: ChunkKind::F32, bytes: data.iter().flat_map(|v| v.to_le_bytes()).collect() }
    }
    pub fn u16s(name: &str, data: &[u16]) -> Self {
        Self { name: name.into(), kind: ChunkKind::U16, bytes: data.iter().flat_map(|v| v.to_le_bytes()).collect() }
    }
    pub fn u8s(name: &str, data: &[u8]) -> Self {
        Self { name: name.into(), kind: ChunkKind::U8, bytes: data.to_vec() }
    }
    pub fn as_f32s(&self) -> Option<Vec<f32>> {
        (self.kind == ChunkKind::F32 && self.bytes.len().is_multiple_of(4))
            .then(|| self.bytes.as_chunks::<4>().0.iter().map(|b| f32::from_le_bytes(*b)).collect())
    }
    pub fn as_u16s(&self) -> Option<Vec<u16>> {
        (self.kind == ChunkKind::U16 && self.bytes.len().is_multiple_of(2))
            .then(|| self.bytes.as_chunks::<2>().0.iter().map(|b| u16::from_le_bytes(*b)).collect())
    }
}

/// An in-memory world file.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct WorldFile {
    pub header: serde_json::Value,
    pub chunks: Vec<Chunk>,
}

impl WorldFile {
    pub fn chunk(&self, name: &str) -> Option<&Chunk> {
        self.chunks.iter().find(|c| c.name == name)
    }

    pub fn write_to(&self, mut w: impl Write) -> io::Result<()> {
        w.write_all(MAGIC)?;
        w.write_all(&FORMAT_VERSION.to_le_bytes())?;
        let header = serde_json::to_vec(&self.header).map_err(io::Error::other)?;
        w.write_all(&(header.len() as u32).to_le_bytes())?;
        w.write_all(&header)?;
        for c in &self.chunks {
            w.write_all(&(c.name.len() as u16).to_le_bytes())?;
            w.write_all(c.name.as_bytes())?;
            w.write_all(&[c.kind as u8])?;
            w.write_all(&(c.bytes.len() as u64).to_le_bytes())?;
            w.write_all(&c.bytes)?;
        }
        Ok(())
    }

    pub fn read_from(mut r: impl Read) -> io::Result<Self> {
        let bad = |m: &str| io::Error::new(io::ErrorKind::InvalidData, m.to_string());
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(bad("not a RealmGenesis world file"));
        }
        let format = read_u32(&mut r)?;
        if format != FORMAT_VERSION {
            return Err(bad(&format!("unsupported world file format {format} (this build reads {FORMAT_VERSION})")));
        }
        let hlen = read_u32(&mut r)? as usize;
        let mut header = vec![0u8; hlen];
        r.read_exact(&mut header)?;
        let header = serde_json::from_slice(&header).map_err(|e| bad(&format!("bad header: {e}")))?;
        let mut chunks = Vec::new();
        loop {
            let mut nl = [0u8; 2];
            match r.read_exact(&mut nl) {
                Ok(()) => {}
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e),
            }
            let mut name = vec![0u8; u16::from_le_bytes(nl) as usize];
            r.read_exact(&mut name)?;
            let name = String::from_utf8(name).map_err(|_| bad("chunk name is not utf-8"))?;
            let mut kind = [0u8; 1];
            r.read_exact(&mut kind)?;
            let kind = match kind[0] {
                0 => ChunkKind::Bytes,
                1 => ChunkKind::F32,
                2 => ChunkKind::U16,
                3 => ChunkKind::U8,
                k => return Err(bad(&format!("unknown chunk kind {k}"))),
            };
            let mut len = [0u8; 8];
            r.read_exact(&mut len)?;
            let mut bytes = vec![0u8; u64::from_le_bytes(len) as usize];
            r.read_exact(&mut bytes)?;
            chunks.push(Chunk { name, kind, bytes });
        }
        Ok(Self { header, chunks })
    }
}

fn read_u32(r: &mut impl Read) -> io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_is_byte_exact() {
        let f = WorldFile {
            header: serde_json::json!({"seed": 42, "name": "Test"}),
            chunks: vec![
                Chunk::f32s("height", &[0.0, -1.5, 8848.0, f32::MIN_POSITIVE]),
                Chunk::u16s("plate", &[0, 7, 65535]),
                Chunk::u8s("crust", &[0, 1, 1]),
            ],
        };
        let mut buf = Vec::new();
        f.write_to(&mut buf).unwrap();
        let g = WorldFile::read_from(&buf[..]).unwrap();
        assert_eq!(f, g);
        let mut buf2 = Vec::new();
        g.write_to(&mut buf2).unwrap();
        assert_eq!(buf, buf2);
        assert_eq!(g.chunk("height").unwrap().as_f32s().unwrap()[2], 8848.0);
    }

    #[test]
    fn rejects_garbage() {
        assert!(WorldFile::read_from(&b"NOTAWORLD_______"[..]).is_err());
    }
}
