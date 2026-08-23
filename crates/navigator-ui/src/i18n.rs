//! The view the UI has of the shared i18n catalog. That catalog lives in `navigator-domain`, so
//! that the layers below the UI can localize too. Those layers are the prose of the Subject Brief,
//! and the HTML report export. See [`navigator_domain::i18n`]. This re-exports the catalog, and does not wrap
//! it, so that `crate::i18n::tr` and `NavigatorApp::tr` do not change.

// This does not re-export `tr_fmt` (positional interpolation). No UI string needs arguments yet,
// and a re-export that nothing uses is code with no purpose in a binary crate. Add it here when a
// string does need arguments.
pub use navigator_domain::i18n::{load_lang, save_lang, tr, Lang};
