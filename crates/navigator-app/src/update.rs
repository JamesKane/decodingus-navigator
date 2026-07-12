//! Checking GitHub Releases for a newer installer, and *notifying* the user — never auto-updating.
//!
//! Installers are published to the GitHub Releases of `JamesKane/decodingus-navigator` (`v*` tags;
//! Alpha/Beta/RC builds are marked *prerelease*). This module fetches that list, finds the highest
//! version, compares it against the running build ([`env!("CARGO_PKG_VERSION")`]), and — if it's
//! newer and the user hasn't chosen to skip it — returns an [`UpdateInfo`] pointing at the release
//! page and the platform-appropriate installer asset. The UI turns that into a dismissible prompt;
//! downloading/installing is entirely the user's choice.

use serde::Deserialize;

use crate::error::AppError;
use crate::settings::AppSettings;

/// The GitHub Releases API for the app repo. We list releases (rather than `/releases/latest`, which
/// excludes prereleases) so Alpha/Beta builds are considered too.
const RELEASES_URL: &str = "https://api.github.com/repos/JamesKane/decodingus-navigator/releases";

/// A newer installer is available. Serialized so it can cross the worker `Command`/`Event` channel.
#[derive(Debug, Clone, serde::Serialize)]
pub struct UpdateInfo {
    /// The running build's version (`CARGO_PKG_VERSION`).
    pub current_version: String,
    /// The newest published version (tag, without a leading `v`).
    pub latest_version: String,
    /// The release's display name (falls back to the tag).
    pub name: String,
    /// The GitHub release page — always present, used as the fallback download link.
    pub release_url: String,
    /// The direct download URL for this platform's installer asset, if one matched.
    pub download_url: Option<String>,
    /// ISO-8601 publish timestamp, if the API reported one.
    pub published_at: Option<String>,
    /// Whether the newest release is a prerelease (Alpha/Beta/RC).
    pub prerelease: bool,
    /// The release notes (Markdown); the UI truncates for display.
    pub notes: String,
}

#[derive(Debug, Deserialize)]
struct GhRelease {
    tag_name: String,
    #[serde(default)]
    name: Option<String>,
    #[serde(default)]
    prerelease: bool,
    #[serde(default)]
    draft: bool,
    html_url: String,
    #[serde(default)]
    body: Option<String>,
    #[serde(default)]
    published_at: Option<String>,
    #[serde(default)]
    assets: Vec<GhAsset>,
}

#[derive(Debug, Deserialize)]
struct GhAsset {
    name: String,
    browser_download_url: String,
}

impl crate::App {
    /// Check GitHub Releases for a newer installer than the running build. Returns `Ok(None)` when
    /// already current (or the newest release is one the user chose to skip); `Ok(Some(info))` when
    /// a newer version is available. Network/parse failures are surfaced as [`AppError::Update`] —
    /// callers treat a failed check as non-fatal.
    pub async fn check_for_update(&self) -> Result<Option<UpdateInfo>, AppError> {
        let current_str = env!("CARGO_PKG_VERSION");
        let current = Version::parse(current_str)
            .ok_or_else(|| AppError::Update(format!("unparseable build version {current_str}")))?;

        let releases = fetch_releases().await?;
        // Highest version among non-draft releases (prereleases included, so Alpha→Alpha upgrades
        // are offered). `max_by` needs the parsed version; releases with a non-version tag are skipped.
        let best = releases
            .into_iter()
            .filter(|r| !r.draft)
            .filter_map(|r| Version::parse(&r.tag_name).map(|v| (v, r)))
            .max_by(|(a, _), (b, _)| a.cmp(b));

        let Some((latest, rel)) = best else {
            return Ok(None);
        };
        if latest <= current {
            return Ok(None);
        }

        let latest_version = rel.tag_name.trim_start_matches(['v', 'V']).to_string();
        // Honor a "skip this version" choice — but a version *newer* than the skipped one still
        // notifies (the skip is keyed to the exact version string).
        if AppSettings::load().skip_update_version.as_deref() == Some(latest_version.as_str()) {
            return Ok(None);
        }

        Ok(Some(UpdateInfo {
            current_version: current_str.to_string(),
            latest_version,
            name: rel.name.clone().unwrap_or_else(|| rel.tag_name.clone()),
            release_url: rel.html_url,
            download_url: pick_installer_asset(&rel.assets),
            published_at: rel.published_at,
            prerelease: rel.prerelease,
            notes: rel.body.unwrap_or_default(),
        }))
    }
}

async fn fetch_releases() -> Result<Vec<GhRelease>, AppError> {
    // GitHub requires a User-Agent. reqwest is already a dependency (json + rustls-tls).
    let client = reqwest::Client::builder()
        .user_agent(concat!("DUNavigator/", env!("CARGO_PKG_VERSION")))
        .build()
        .map_err(|e| AppError::Update(e.to_string()))?;
    let resp = client
        .get(RELEASES_URL)
        .query(&[("per_page", "30")])
        .header("Accept", "application/vnd.github+json")
        .send()
        .await
        .map_err(|e| AppError::Update(e.to_string()))?;
    if !resp.status().is_success() {
        return Err(AppError::Update(format!("GitHub returned {}", resp.status())));
    }
    resp.json::<Vec<GhRelease>>()
        .await
        .map_err(|e| AppError::Update(e.to_string()))
}

/// Pick the installer asset for the current platform from a release's assets. macOS ships a single
/// universal2 `.dmg`; Windows an NSIS `*-setup.exe`; Linux an `.AppImage` / `.deb` per-arch. Returns
/// the first name-matching asset (preferring one whose name carries this arch or "universal").
fn pick_installer_asset(assets: &[GhAsset]) -> Option<String> {
    let exts: &[&str] = if cfg!(target_os = "macos") {
        &[".dmg"]
    } else if cfg!(target_os = "windows") {
        &["-setup.exe", ".msi", ".exe"]
    } else {
        &[".appimage", ".deb"]
    };
    let arch = std::env::consts::ARCH; // "x86_64" | "aarch64" | ...
    let matches: Vec<&GhAsset> = assets
        .iter()
        .filter(|a| {
            let n = a.name.to_ascii_lowercase();
            exts.iter().any(|e| n.ends_with(e))
        })
        .collect();
    // Prefer an arch- or universal-tagged asset; else the first match.
    matches
        .iter()
        .find(|a| {
            let n = a.name.to_ascii_lowercase();
            n.contains(arch) || n.contains("universal")
        })
        .or_else(|| matches.first())
        .map(|a| a.browser_download_url.clone())
}

/// A minimal `MAJOR.MINOR.PATCH[-prerelease]` version. Ordered so a release outranks its own
/// prereleases (`0.2.0` > `0.2.0-alpha.1`) and higher numbers win — sufficient for our `vX.Y.Z`
/// release tags; we deliberately don't implement full SemVer build-metadata precedence.
#[derive(Debug, PartialEq, Eq)]
struct Version {
    nums: (u64, u64, u64),
    pre: Option<String>,
}

impl Version {
    fn parse(s: &str) -> Option<Version> {
        let s = s.trim().trim_start_matches(['v', 'V']);
        let (core, pre) = match s.split_once('-') {
            Some((c, p)) => (c, Some(p.to_string())),
            None => (s, None),
        };
        let mut it = core.split('.');
        let major = it.next()?.parse::<u64>().ok()?;
        let minor = it.next().unwrap_or("0").parse::<u64>().ok()?;
        let patch = it.next().unwrap_or("0").parse::<u64>().ok()?;
        Some(Version {
            nums: (major, minor, patch),
            pre,
        })
    }
}

impl Ord for Version {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        match self.nums.cmp(&other.nums) {
            Ordering::Equal => match (&self.pre, &other.pre) {
                // No prerelease outranks a prerelease at the same version.
                (None, None) => Ordering::Equal,
                (None, Some(_)) => Ordering::Greater,
                (Some(_), None) => Ordering::Less,
                (Some(a), Some(b)) => a.cmp(b),
            },
            ord => ord,
        }
    }
}

impl PartialOrd for Version {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(s: &str) -> Version {
        Version::parse(s).unwrap()
    }

    #[test]
    fn version_ordering() {
        assert!(v("0.2.0") > v("0.1.0"));
        assert!(v("v0.2.0") > v("0.1.9"));
        assert!(v("1.0.0") > v("0.9.9"));
        assert!(v("0.2.1") > v("0.2.0"));
        // A release outranks its prereleases; prereleases compare lexically.
        assert!(v("0.2.0") > v("0.2.0-alpha"));
        assert!(v("0.2.0-beta") > v("0.2.0-alpha"));
        assert!(v("0.2.0-alpha.2") > v("0.2.0-alpha.1"));
        assert_eq!(v("0.1.0"), v("v0.1.0"));
    }

    #[test]
    fn parse_shapes() {
        assert_eq!(v("0.1.0").nums, (0, 1, 0));
        assert_eq!(v("v2").nums, (2, 0, 0));
        assert_eq!(v("1.5").nums, (1, 5, 0));
        assert_eq!(v("0.2.0-alpha.1").pre.as_deref(), Some("alpha.1"));
        assert!(Version::parse("not-a-version").is_none());
    }

    #[test]
    fn picks_platform_asset() {
        let assets = vec![
            GhAsset {
                name: "DUNavigator_0.2.0_universal.dmg".into(),
                browser_download_url: "https://example/dmg".into(),
            },
            GhAsset {
                name: "DUNavigator_0.2.0_x64-setup.exe".into(),
                browser_download_url: "https://example/exe".into(),
            },
            GhAsset {
                name: "SHA256SUMS".into(),
                browser_download_url: "https://example/sums".into(),
            },
        ];
        let picked = pick_installer_asset(&assets);
        assert!(picked.is_some());
        if cfg!(target_os = "macos") {
            assert_eq!(picked.as_deref(), Some("https://example/dmg"));
        } else if cfg!(target_os = "windows") {
            assert_eq!(picked.as_deref(), Some("https://example/exe"));
        }
    }
}
