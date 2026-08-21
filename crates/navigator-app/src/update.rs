//! A check of GitHub Releases for a newer installer. This module only *tells* the user. It never
//! updates the app without a command from the user.
//!
//! The team publishes each installer to the GitHub Releases of `JamesKane/decodingus-navigator`,
//! under a `v*` tag. An Alpha build, a Beta build, and an RC build each have the *prerelease*
//! mark.
//!
//! This module reads that list and finds the highest version. It then compares that version with
//! the version of this build, which is [`BUILD_VERSION`]. The module returns an [`UpdateInfo`]
//! value when the published version is newer and the user did not choose to skip it. The value
//! points to the release page and to the correct installer for the platform.
//!
//! The UI shows this information as a prompt that the user can close. Only the user can start the
//! download and the installation.

use serde::Deserialize;

use crate::error::AppError;
use crate::settings::AppSettings;

/// The GitHub Releases API for the repository of the app. The code lists all releases. It does not
/// use `/releases/latest`, because that endpoint hides a prerelease. So the check also sees an
/// Alpha build and a Beta build.
const RELEASES_URL: &str = "https://api.github.com/repos/JamesKane/decodingus-navigator/releases";

/// The version name of this build. The code compares this name with the published releases.
///
/// `CARGO_PKG_VERSION` alone does not work for that comparison. The reason is easy to miss. The
/// workspace version is a plain `0.1.0`, but each tag that the team ships is `v0.1.0-alpha.N`.
///
/// In SemVer, a release has a higher rank than its own prereleases. So this build always looked
/// *newer* than each alpha on the server, and the check always answered "up to date". The team
/// published sixteen alphas, and no user received a notification.
///
/// The workspace version stays numeric by design, because a Windows installer format needs `x.y.z`.
/// So the package step puts the full version here, from the tag of that build. A build on a
/// developer machine has no tag and uses the default value. That result is correct, because such a
/// build must not claim to be a release.
///
/// A note for a local test. Rust reads `option_env!` at compile time, and Cargo does not treat the
/// variable as an input. So a change to the variable does not start a rebuild. Touch this file, or
/// build clean. If not, you continue to see the value from the earlier build. CI always builds from
/// a new checkout, so this problem never reaches a packaged artifact.
pub(crate) const BUILD_VERSION: &str = match option_env!("NAVIGATOR_RELEASE_VERSION") {
    Some(v) => v,
    None => env!("CARGO_PKG_VERSION"),
};

/// A newer installer is available. Serialized so it can cross the worker `Command`/`Event` channel.
#[derive(Debug, Clone, serde::Serialize)]
pub struct UpdateInfo {
    /// The version of this build ([`BUILD_VERSION`]).
    pub current_version: String,
    /// The newest published version. The value is the tag with no `v` at the start.
    pub latest_version: String,
    /// The release's display name (falls back to the tag).
    pub name: String,
    /// The GitHub release page. This value is always present, and the UI uses it as the second
    /// download link.
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
    /// Check GitHub Releases for an installer that is newer than this build.
    ///
    /// The method returns `Ok(None)` when this build is current. It also returns `Ok(None)` when
    /// the user chose to skip the newest release. It returns `Ok(Some(info))` when a newer version
    /// exists.
    ///
    /// A network fault or a parse fault becomes an [`AppError::Update`] value. A caller must
    /// continue after a failed check.
    pub async fn check_for_update(&self) -> Result<Option<UpdateInfo>, AppError> {
        let current_str = BUILD_VERSION;
        let current = Version::parse(current_str)
            .ok_or_else(|| AppError::Update(format!("unparseable build version {current_str}")))?;

        let releases = fetch_releases().await?;
        // Find the highest version among the releases that are not a draft. The set holds each
        // prerelease, so the check offers an upgrade from one Alpha to the next Alpha. `max_by`
        // needs the parsed version. The code skips a release when its tag is not a version.
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
        // Obey a "skip this version" choice from the user. But a version that is *newer* than the
        // skipped version still gives a notification, because the skip holds one exact version
        // string.
        if AppSettings::load().skip_update_version.as_deref() == Some(latest_version.as_str()) {
            return Ok(None);
        }

        Ok(Some(UpdateInfo {
            // Remove the same characters that the code removes from `latest_version`. The
            // "current > latest" text in the UI is then consistent. The build step supplies a tag,
            // and a tag starts with `v`.
            current_version: current_str.trim_start_matches(['v', 'V']).to_string(),
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

/// Select the installer asset for this platform from the assets of a release.
///
/// macOS has one universal2 `.dmg` file. Windows has an NSIS `*-setup.exe` file. Linux has an
/// `.AppImage` file and a `.deb` file for each architecture.
///
/// The function returns the first asset with a name that matches. It selects an asset whose name
/// holds this architecture or the word "universal" before it selects another asset.
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

/// A small `MAJOR.MINOR.PATCH[-prerelease]` version.
///
/// The order puts a release above its own prereleases, so `0.2.0` is above `0.2.0-alpha.1`. A
/// higher number also wins. This order is enough for a `vX.Y.Z` release tag. The code does not
/// implement the full SemVer rule for build metadata, by design.
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
                (Some(a), Some(b)) => compare_prerelease(a, b),
            },
            ord => ord,
        }
    }
}

/// Compare two prerelease strings dot-part by dot-part, numerically where both parts are numbers.
///
/// A plain string compare is wrong when a counter reaches two digits. `"alpha.9"` then sorts
/// *above* `"alpha.16"`, because the character `9` is above the character `1`.
///
/// This fault occurred in this project. The version reached `alpha.16` with the old comparator, and
/// a search for the highest published tag returned `alpha.9`. The SemVer rule is the correction. A
/// numeric part compares as a number, and a numeric part ranks below a part with letters.
fn compare_prerelease(a: &str, b: &str) -> std::cmp::Ordering {
    use std::cmp::Ordering;

    let mut left = a.split('.');
    let mut right = b.split('.');
    loop {
        match (left.next(), right.next()) {
            (None, None) => return Ordering::Equal,
            // A prerelease with more parts outranks its own prefix: alpha.1 > alpha.
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (Some(x), Some(y)) => {
                let ord = match (x.parse::<u64>(), y.parse::<u64>()) {
                    (Ok(xn), Ok(yn)) => xn.cmp(&yn),
                    (Ok(_), Err(_)) => Ordering::Less,
                    (Err(_), Ok(_)) => Ordering::Greater,
                    (Err(_), Err(_)) => x.cmp(y),
                };
                if ord != Ordering::Equal {
                    return ord;
                }
            }
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
        // Two digits: the case that was wrong for sixteen releases. A string compare puts
        // "alpha.9" above "alpha.16" because '9' > '1'.
        assert!(v("0.1.0-alpha.16") > v("0.1.0-alpha.9"));
        assert!(v("0.1.0-alpha.17") > v("0.1.0-alpha.16"));
        assert!(v("0.1.0-alpha.100") > v("0.1.0-alpha.99"));
        // A longer prerelease outranks its own prefix, and numeric parts rank below alphanumeric.
        assert!(v("0.1.0-alpha.1") > v("0.1.0-alpha"));
        assert!(v("0.1.0-alpha.beta") > v("0.1.0-alpha.1"));
        assert_eq!(v("0.1.0"), v("v0.1.0"));
    }

    /// The fault that made the full feature inert. A build with the name `0.1.0` ranks above each
    /// `0.1.0-alpha.N` tag. So the newest alpha never looked newer, and no user received a
    /// notification. A release build must use its own prerelease name, or the comparison fails.
    #[test]
    fn a_release_build_compares_against_the_alphas_it_shipped_beside() {
        // What the bare workspace version did.
        assert!(
            v("v0.1.0-alpha.17") < v("0.1.0"),
            "a bare release version outranks its prereleases — this is why the check went silent"
        );
        // What an injected release version does instead.
        assert!(v("v0.1.0-alpha.17") > v("0.1.0-alpha.16"));
    }

    /// The parser must accept the name of this build. If not, the check stops with the message
    /// "unparseable build version" and does no useful work.
    #[test]
    fn the_build_version_is_a_version() {
        assert!(
            Version::parse(BUILD_VERSION).is_some(),
            "unparseable BUILD_VERSION: {BUILD_VERSION}"
        );
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
        // Add an asset for each platform, so the test is useful on any CI runner. The runners are
        // macOS, Windows, and Linux. On Linux the code looks for a .AppImage file or a .deb
        // file.
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
                name: "DUNavigator_0.2.0_x86_64.AppImage".into(),
                browser_download_url: "https://example/appimage".into(),
            },
            GhAsset {
                name: "DUNavigator_0.2.0_amd64.deb".into(),
                browser_download_url: "https://example/deb".into(),
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
        } else {
            // On Linux, the code selects a .AppImage file or a .deb file. It first looks for an
            // asset with this architecture or the word "universal" in its name. One example is the
            // x86_64 AppImage on an x86_64 runner. If it finds none, it takes the first file with
            // a correct extension.
            assert!(matches!(
                picked.as_deref(),
                Some("https://example/appimage") | Some("https://example/deb")
            ));
        }
    }
}
