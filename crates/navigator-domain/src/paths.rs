//! Where the data of each user lives: the one place that answers "what is home?".
//!
//! Everything the app persists hangs off `~/.decodingus` (workspace DB, references, liftover
//! chains, trees, panels, config). Six places used to derive that root on their own, each with
//! `std::env::var("HOME").unwrap_or(".")`. That is wrong on Windows, where `HOME` is normally
//! unset. Every one of those paths then resolved **relative to the current directory**, with no
//! message. The app made a fresh `.decodingus` tree wherever somebody started it, and two starts
//! did not always share one.
//!
//! [`home_dir`] resolves the real home of the platform instead, and the callers join their own
//! subpaths onto it. No `dirs` or `directories` dependency, on purpose: the rules are short enough
//! to state, and to test, directly.

use std::ffi::OsString;
use std::path::PathBuf;

/// The current user's home directory, or `None` if the platform's variables do not say.
///
/// Unix: `$HOME`. Windows: `%USERPROFILE%`, else `%HOMEDRIVE%%HOMEPATH%`.
pub fn home_dir() -> Option<PathBuf> {
    #[cfg(windows)]
    {
        windows_home(
            std::env::var_os("USERPROFILE"),
            std::env::var_os("HOMEDRIVE"),
            std::env::var_os("HOMEPATH"),
        )
    }
    #[cfg(not(windows))]
    {
        std::env::var_os("HOME").filter(|h| !h.is_empty()).map(PathBuf::from)
    }
}

/// The `~/.decodingus` root: every persisted artifact hangs off this. It falls back to a relative
/// `.decodingus` only when the platform reports no home at all. That is the previous behaviour,
/// and it stays. An environment with no home (some CI containers) then still runs, and does not
/// fail at startup.
pub fn decodingus_dir() -> PathBuf {
    home_dir().unwrap_or_else(|| PathBuf::from(".")).join(".decodingus")
}

/// Windows home resolution, split out from [`home_dir`] so its precedence is testable on any
/// platform.
///
/// `%USERPROFILE%` first, and `%HOMEDRIVE%` with `%HOMEPATH%` second, for a domain profile where
/// `%USERPROFILE%` is unset. This does **not** read `%HOME%`. MSYS2 and Git-Bash set it to a POSIX
/// path like `/c/Users/name`, which native Windows APIs read as a rooted path on the *current
/// drive*. To obey it would put user data in a wrong location that looks correct, and that is
/// worse than the fallback.
// Compiled on every platform so its precedence stays under test anywhere; only *called* on Windows.
#[cfg_attr(not(windows), allow(dead_code))]
fn windows_home(
    userprofile: Option<OsString>,
    homedrive: Option<OsString>,
    homepath: Option<OsString>,
) -> Option<PathBuf> {
    if let Some(p) = userprofile.filter(|p| !p.is_empty()) {
        return Some(PathBuf::from(p));
    }
    match (homedrive.filter(|d| !d.is_empty()), homepath.filter(|p| !p.is_empty())) {
        (Some(drive), Some(path)) => {
            let mut joined = drive;
            joined.push(path);
            Some(PathBuf::from(joined))
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn os(s: &str) -> Option<OsString> {
        Some(OsString::from(s))
    }

    #[test]
    fn windows_home_prefers_userprofile_then_homedrive_homepath() {
        assert_eq!(
            windows_home(os(r"C:\Users\ada"), os("D:"), os(r"\profiles\ada")),
            Some(PathBuf::from(r"C:\Users\ada"))
        );
        // Domain profile: no USERPROFILE, so this joins the drive and the path pair verbatim.
        assert_eq!(
            windows_home(None, os("D:"), os(r"\profiles\ada")),
            Some(PathBuf::from(r"D:\profiles\ada"))
        );
        // Empty is as good as unset (a set-but-blank variable must not win).
        assert_eq!(
            windows_home(os(""), os("D:"), os(r"\profiles\ada")),
            Some(PathBuf::from(r"D:\profiles\ada"))
        );
        // Half a pair is no answer. It is better to fall back than to build a half-formed path.
        assert_eq!(windows_home(None, os("D:"), None), None);
        assert_eq!(windows_home(None, None, os(r"\profiles\ada")), None);
        assert_eq!(windows_home(None, None, None), None);
    }

    /// The fallback stays relative, and does not panic: a container with no home still runs.
    #[test]
    fn decodingus_dir_ends_in_the_conventional_directory() {
        assert!(decodingus_dir().ends_with(".decodingus"));
    }
}
