//! Build the asset integrity manifest (`ancestry_manifest_<build>.json`) — run after the other
//! `panelbuild` steps. Hashes every `*_<build>.bin` in the ancestry directory so the app can verify
//! a loaded asset and refuse a corrupt / truncated CDN download.

use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::manifest::AssetManifest;

#[derive(Parser)]
pub struct ManifestArgs {
    /// The ancestry asset directory (where the `.bin` files live).
    #[arg(long)]
    pub dir: PathBuf,
    /// Build label; only files ending `_<build>.bin` are hashed.
    #[arg(long, default_value = "chm13v2.0")]
    pub build: String,
    /// Output manifest path (defaults to `<dir>/ancestry_manifest_<build>.json`).
    #[arg(long)]
    pub out: Option<PathBuf>,
}

pub fn build_manifest(args: ManifestArgs) -> Result<()> {
    let suffix = format!("_{}.bin", args.build);
    let mut manifest = AssetManifest { build: args.build.clone(), generated_at: String::new(), assets: Default::default() };
    for entry in fs::read_dir(&args.dir).with_context(|| format!("read dir {}", args.dir.display()))? {
        let path = entry?.path();
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else { continue };
        if !name.ends_with(&suffix) {
            continue;
        }
        let bytes = fs::read(&path).with_context(|| format!("read {}", path.display()))?;
        manifest.insert(name.to_string(), &bytes);
    }
    anyhow::ensure!(!manifest.assets.is_empty(), "no `*{suffix}` assets found in {}", args.dir.display());

    let out = args.out.unwrap_or_else(|| args.dir.join(format!("ancestry_manifest_{}.json", args.build)));
    let json = manifest.to_json().map_err(|e| anyhow::anyhow!(e.to_string()))?;
    File::create(&out)?.write_all(json.as_bytes())?;
    eprintln!("wrote {} ({} assets)", out.display(), manifest.assets.len());
    Ok(())
}
