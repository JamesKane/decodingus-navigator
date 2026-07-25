//! The UI's view of the shared i18n catalog, which lives in `navigator-domain` so that the layers
//! below the UI — the Subject Brief's prose, the HTML report export — can localize too. See
//! [`navigator_domain::i18n`]. Re-exported rather than wrapped so `crate::i18n::tr` and
//! `NavigatorApp::tr` keep working unchanged.

// `tr_fmt` (positional interpolation) is not re-exported: no UI string needs arguments yet,
// and an unused re-export is dead code in a binary crate. Add it here when one does.
pub use navigator_domain::i18n::{load_lang, save_lang, tr, Lang};
