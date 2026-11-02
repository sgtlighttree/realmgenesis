//! Save, load and export. Native builds use a fixed folder, `~/RealmGenesis/`, so no file
//! dialog dependency is needed for the MVP. The web preview can't write files from inside
//! its sandbox, so it reports that instead.

use crate::app::RealmApp;
use crate::commands::Command;

pub fn handle(app: &mut RealmApp, c: Command, rs: &egui_wgpu::RenderState) {
    let result = run(app, &c, rs);
    app.status = match result {
        Ok(msg) => msg,
        Err(e) => format!("⚠ {e}"),
    };
}

#[cfg(target_arch = "wasm32")]
fn run(_app: &mut RealmApp, _c: &Command, _rs: &egui_wgpu::RenderState) -> Result<String, String> {
    Err("Saving and exporting work in the desktop app; the web preview can't write files.".into())
}

#[cfg(not(target_arch = "wasm32"))]
fn run(app: &mut RealmApp, c: &Command, rs: &egui_wgpu::RenderState) -> Result<String, String> {
    use std::fs;
    let dir = folder()?;
    fs::create_dir_all(&dir).map_err(|e| format!("can't create {}: {e}", dir.display()))?;
    match c {
        Command::SaveWorld => {
            let w = app.world.as_ref().ok_or("nothing to save yet")?;
            let path = dir.join(format!("world-{}.rgworld", w.params.tectonics.seed));
            let file = fs::File::create(&path).map_err(|e| e.to_string())?;
            rg_pipeline::file::to_file(w)
                .write_to(std::io::BufWriter::new(file))
                .map_err(|e| e.to_string())?;
            Ok(format!("Saved {}", path.display()))
        }
        Command::LoadWorld => {
            let latest = fs::read_dir(&dir)
                .map_err(|e| e.to_string())?
                .filter_map(|e| e.ok())
                .filter(|e| e.path().extension().is_some_and(|x| x == "rgworld"))
                .max_by_key(|e| e.metadata().and_then(|m| m.modified()).ok())
                .ok_or(format!("no .rgworld files in {}", dir.display()))?
                .path();
            let file = fs::File::open(&latest).map_err(|e| e.to_string())?;
            let wf = rg_world::WorldFile::read_from(std::io::BufReader::new(file)).map_err(|e| e.to_string())?;
            let world = rg_pipeline::file::from_file(&wf)?;
            app.params = world.params.clone();
            app.world = Some(world);
            app.after_world_changed(rs, true);
            Ok(format!("Loaded {}", latest.display()))
        }
        Command::ExportHeightmap => {
            let w = app.world.as_ref().ok_or("nothing to export yet")?;
            let ex = rg_pipeline::export::heightmap_png(w, 4096)?;
            let base = dir.join(format!("heightmap-{}", w.params.tectonics.seed));
            fs::write(base.with_extension("png"), &ex.png).map_err(|e| e.to_string())?;
            fs::write(base.with_extension("json"), &ex.sidecar).map_err(|e| e.to_string())?;
            Ok(format!("Exported {}.png (+ .json)", base.display()))
        }
        _ => Ok(String::new()),
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn folder() -> Result<std::path::PathBuf, String> {
    let home = std::env::var_os("HOME")
        .or_else(|| std::env::var_os("USERPROFILE"))
        .ok_or("can't find the home folder")?;
    Ok(std::path::PathBuf::from(home).join("RealmGenesis"))
}
