//! Lightweight i18n for the egui Navigator: Play-style `key=value` catalogs embedded at
//! compile time, mirroring the AppView's (`decodingus/rust` du-web) approach so both Rust
//! front-ends share one catalog format. Dependency-free (no fluent).
//!
//! Lookup falls back from the active language → English → the key itself, so a partial
//! translation degrades to English rather than showing raw keys. Catalog values are
//! `&'static str`, so `tr()` returns `&'static str` and never borrows app state — convenient
//! inside egui closures.

use std::collections::HashMap;
use std::sync::OnceLock;

/// A supported UI language.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Lang {
    En,
    Es,
}

impl Lang {
    /// Language code (used to persist the chosen locale; see [`save_lang`]).
    pub fn code(self) -> &'static str {
        match self {
            Lang::En => "en",
            Lang::Es => "es",
        }
    }

    /// Native display name (for the language switcher).
    pub fn label(self) -> &'static str {
        match self {
            Lang::En => "English",
            Lang::Es => "Español",
        }
    }

    /// Resolve from a code (e.g. an env var or saved setting); `None` if unsupported.
    pub fn parse(s: &str) -> Option<Lang> {
        match s.get(0..2).map(str::to_ascii_lowercase).as_deref() {
            Some("en") => Some(Lang::En),
            Some("es") => Some(Lang::Es),
            _ => None,
        }
    }

    /// All languages, for rendering the switcher.
    pub fn all() -> &'static [Lang] {
        &[Lang::En, Lang::Es]
    }
}

/// Path of the persisted-language file (`~/.decodingus/navigator-lang`), matching the
/// `~/.decodingus` convention used for the workspace DB.
fn lang_file() -> std::path::PathBuf {
    let home = std::env::var("HOME").unwrap_or_else(|_| ".".to_string());
    std::path::PathBuf::from(home)
        .join(".decodingus")
        .join("navigator-lang")
}

/// The previously chosen UI language, if one was saved.
pub fn load_lang() -> Option<Lang> {
    std::fs::read_to_string(lang_file())
        .ok()
        .and_then(|s| Lang::parse(s.trim()))
}

/// Persist the chosen UI language so it survives a restart (best-effort; I/O errors ignored).
pub fn save_lang(lang: Lang) {
    let path = lang_file();
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let _ = std::fs::write(path, lang.code());
}

const EN_SRC: &str = include_str!("../locales/en.txt");
const ES_SRC: &str = include_str!("../locales/es.txt");

fn parse_catalog(src: &'static str) -> HashMap<&'static str, &'static str> {
    src.lines()
        .filter_map(|line| {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                return None;
            }
            line.split_once('=').map(|(k, v)| (k.trim(), v.trim()))
        })
        .collect()
}

fn catalog(lang: Lang) -> &'static HashMap<&'static str, &'static str> {
    static EN: OnceLock<HashMap<&str, &str>> = OnceLock::new();
    static ES: OnceLock<HashMap<&str, &str>> = OnceLock::new();
    match lang {
        Lang::En => EN.get_or_init(|| parse_catalog(EN_SRC)),
        Lang::Es => ES.get_or_init(|| parse_catalog(ES_SRC)),
    }
}

/// Translate `key` for `lang`, falling back to English then the key itself.
pub fn tr(lang: Lang, key: &'static str) -> &'static str {
    if let Some(v) = catalog(lang).get(key).copied() {
        return v;
    }
    if lang != Lang::En {
        if let Some(v) = catalog(Lang::En).get(key).copied() {
            return v;
        }
    }
    key
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translates_and_falls_back() {
        assert_eq!(tr(Lang::En, "nav.subjects"), "Subjects");
        assert_eq!(tr(Lang::Es, "nav.subjects"), "Sujetos");
        // Missing in Es → English fallback (assuming this key isn't translated).
        assert_eq!(tr(Lang::Es, "status.label"), tr(Lang::Es, "status.label"));
        // Unknown key → the key itself.
        assert_eq!(tr(Lang::En, "totally.unknown.key"), "totally.unknown.key");
    }

    #[test]
    fn every_es_key_exists_in_en() {
        // Catch typos: a translated key with no English source would never be reachable.
        let en = catalog(Lang::En);
        for k in catalog(Lang::Es).keys() {
            assert!(en.contains_key(k), "es.txt key not in en.txt: {k}");
        }
    }
}
