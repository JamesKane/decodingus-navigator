//! The grounding context for the M4 "ask my results" chat (see
//! `docs/design/local-llm-expansion.md`). Narration (M1) stays grounded in the lean
//! [`SubjectBrief`](crate::brief::SubjectBrief) fact sheet; the chat needs to answer about *more*
//! signals (Y-STR panels, private-Y variants, mtDNA mutations, IBD matches, genetic sex), so it
//! grounds in a [`ResultsContext`] = the brief plus curated, summary-level facts for those signals.
//!
//! As with [`llm_prompt`](crate::llm_prompt), this layer is pure: the per-signal facts are plain,
//! already-vetted values the app fills in (the domain crate must not depend on analysis/store), and
//! [`results_fact_sheet`] is the unit-tested builder so the exact text we send stays reviewable. The
//! model remains a *rewriter, not a source of facts*: every section is summary-level, and the inline
//! notes keep STR/mtDNA/IBD framed as lineage facts, never health or trait claims.

use crate::brief::SubjectBrief;
use crate::llm_prompt::narrate_fact_sheet;

/// Genetic-sex call reduced to a label + confidence phrase.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SexFact {
    /// e.g. `"Male (XY)"`.
    pub label: String,
    /// `"high"` | `"medium"` | `"low"`.
    pub confidence: String,
}

/// One Y-STR panel: its name and how many markers it carries (summary only — never the raw values).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct YStrPanelFact {
    pub panel: String,
    pub markers: usize,
}

/// Private Y-variant counts, separated by confidence class (the `navigator-app` `PrivateBucket`
/// distinction: novel-in-unique-sequence vs off-path vs structural-region).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrivateYFact {
    /// Novel calls in unique sequence — the high-confidence new-branch candidates.
    pub novel_unique: usize,
    /// Variants off known branches — suggest previously-unknown branch depth.
    pub off_path: usize,
    /// Calls in structural/paralog-prone regions — suspect, to be treated with caution.
    pub structural: usize,
}

/// Above this many novel-in-unique private-Y calls, a single sample is almost certainly reporting
/// artifacts rather than real new-branch candidates. The de-novo tree pipeline's per-WGS-sample novel
/// count runs ~3–39 (median), so a count in the dozens is normal and the low hundreds is a red flag —
/// typically contamination, shallow/uneven coverage, or a reference-build mismatch.
pub const PRIVATE_Y_QC_WARN: usize = 50;

/// A one-line QC banner when the novel-in-unique private-Y count is implausibly high for one sample
/// (see [`PRIVATE_Y_QC_WARN`]), else `None`. Surfaced in reports and the `private-y` CLI so an
/// elevated count reads as "check this sample" rather than "you have this many new branches".
pub fn private_y_qc_banner(novel_unique: usize) -> Option<String> {
    (novel_unique >= PRIVATE_Y_QC_WARN).then(|| {
        format!(
            "⚠ elevated private-Y count ({novel_unique} novel in unique sequence) — unusually high for \
             one sample; check for contamination, low/uneven coverage, or a reference-build mismatch \
             before treating these as real new branches"
        )
    })
}

/// mtDNA differences from rCRS, summarized by region with a few example notations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MtMutationsFact {
    pub total: usize,
    pub hvr1: usize,
    pub hvr2: usize,
    pub coding: usize,
    /// A handful of example notations (e.g. `"263A>G"`), capped by the caller.
    pub examples: Vec<String>,
}

/// One IBD match, reduced to relationship band + sharing — **no identifying details** (the partner is
/// deliberately not named; this is the user's own workspace data described as a relationship, not a
/// person).
#[derive(Debug, Clone, PartialEq)]
pub struct IbdMatchFact {
    pub relationship: String,
    pub total_shared_cm: f64,
    pub segment_count: i64,
}

/// IBD/network matches summary: how many, and the closest by shared cM.
#[derive(Debug, Clone, PartialEq)]
pub struct IbdFact {
    pub match_count: usize,
    pub closest: Option<IbdMatchFact>,
}

/// The brief plus curated summaries of the other signals — the grounding context for the M4 chat.
/// Absent signals are `None` / empty and are simply omitted from the fact sheet (so the model can't
/// restate what isn't there), exactly like the brief's own optional sections.
#[derive(Debug, Clone, PartialEq)]
pub struct ResultsContext {
    pub brief: SubjectBrief,
    pub sex: Option<SexFact>,
    pub ystr: Vec<YStrPanelFact>,
    pub private_y: Option<PrivateYFact>,
    pub mt_mutations: Option<MtMutationsFact>,
    pub ibd: Option<IbdFact>,
}

/// Which result signal a per-tab "Explain this" narration (M5) targets — also the key used to pull a
/// single signal's section out of a [`ResultsContext`] via [`signal_section`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum SignalKind {
    Sex,
    YStr,
    PrivateY,
    MtMutations,
    Ibd,
}

impl SignalKind {
    /// Human label for the signal (used in the focused narration prompt and as a display heading).
    pub fn label(self) -> &'static str {
        match self {
            SignalKind::Sex => "genetic sex",
            SignalKind::YStr => "Y-STR markers",
            SignalKind::PrivateY => "private Y variants",
            SignalKind::MtMutations => "mtDNA mutations",
            SignalKind::Ibd => "DNA relatives (IBD matches)",
        }
    }
}

// --- Per-signal section builders -----------------------------------------------------------------
// Each returns the labelled, summary-level block for one signal (leading `\n`, so they concatenate
// cleanly after the brief's fact sheet), or `None` when the subject has nothing for that signal. The
// blocks are shared verbatim between the M4 chat grounding ([`results_fact_sheet`]) and the M5
// per-tab narration ([`signal_section`]).

fn sex_section(sex: &Option<SexFact>) -> Option<String> {
    let sex = sex.as_ref()?;
    Some(format!("\nGenetic sex:\n- {} ({} confidence)\n", sex.label, sex.confidence))
}

fn ystr_section(ystr: &[YStrPanelFact]) -> Option<String> {
    if ystr.is_empty() {
        return None;
    }
    let mut s = String::from("\nY-STR panels:\n");
    for p in ystr {
        s.push_str(&format!("- {} ({} markers)\n", p.panel, p.markers));
    }
    s.push_str(
        "- note: STR markers describe the paternal-lineage pattern only; they are not trait or \
         health information.\n",
    );
    Some(s)
}

fn private_y_section(private_y: &Option<PrivateYFact>) -> Option<String> {
    let p = private_y.as_ref()?;
    let mut s = String::from("\nPrivate Y variants (relative to the known tree):\n");
    s.push_str(&format!(
        "- {} novel variant(s) in unique sequence — candidates for a new branch\n",
        p.novel_unique
    ));
    s.push_str(&format!(
        "- {} variant(s) off known branches — suggest previously-unknown branch depth\n",
        p.off_path
    ));
    if p.structural > 0 {
        s.push_str(&format!(
            "- {} call(s) in structural/paralog-prone regions — uncertain, not confident new variants\n",
            p.structural
        ));
    }
    if let Some(warn) = private_y_qc_banner(p.novel_unique) {
        s.push_str(&warn);
        s.push('\n');
    }
    Some(s)
}

fn mt_section(mt: &Option<MtMutationsFact>) -> Option<String> {
    let m = mt.as_ref()?;
    let mut s = format!(
        "\nmtDNA mutations (differences from the rCRS reference): {} total — HVR1 {}, HVR2 {}, coding {}\n",
        m.total, m.hvr1, m.hvr2, m.coding
    );
    if !m.examples.is_empty() {
        s.push_str(&format!("- examples: {}\n", m.examples.join(", ")));
    }
    s.push_str("- note: these are maternal-lineage markers, not trait or health information.\n");
    Some(s)
}

fn ibd_section(ibd: &Option<IbdFact>) -> Option<String> {
    let ibd = ibd.as_ref()?;
    let mut s = format!(
        "\nGenetic relatives (IBD matches in this workspace): {}\n",
        ibd.match_count
    );
    if let Some(c) = &ibd.closest {
        s.push_str(&format!(
            "- closest: about {:.0} cM shared across {} segment(s), consistent with {}\n",
            c.total_shared_cm, c.segment_count, c.relationship
        ));
    }
    s.push_str(
        "- note: these are relationships inferred from shared DNA, not genealogically verified, \
         and are described without identifying the other person.\n",
    );
    Some(s)
}

/// The labelled section for a single signal, or `None` when the subject has nothing for it — the
/// grounding for an M5 per-tab "Explain this" narration of just that signal.
pub fn signal_section(ctx: &ResultsContext, kind: SignalKind) -> Option<String> {
    match kind {
        SignalKind::Sex => sex_section(&ctx.sex),
        SignalKind::YStr => ystr_section(&ctx.ystr),
        SignalKind::PrivateY => private_y_section(&ctx.private_y),
        SignalKind::MtMutations => mt_section(&ctx.mt_mutations),
        SignalKind::Ibd => ibd_section(&ctx.ibd),
    }
}

/// Build the chat's grounding fact sheet: the brief's fact sheet (unchanged from narration) plus a
/// labelled, summary-level block per present signal. Pure and unit-tested — this is the exact text
/// that goes in the chat system message as "your only source of facts".
pub fn results_fact_sheet(ctx: &ResultsContext) -> String {
    let mut s = narrate_fact_sheet(&ctx.brief);
    for section in [
        sex_section(&ctx.sex),
        ystr_section(&ctx.ystr),
        private_y_section(&ctx.private_y),
        mt_section(&ctx.mt_mutations),
        ibd_section(&ctx.ibd),
    ]
    .into_iter()
    .flatten()
    {
        s.push_str(&section);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ancestry::SuperPopulationSummary;
    use crate::brief::{AncestryBrief, Headline, LineageBrief, LineageKind, PackStatus, SubjectBrief, TestBrief};
    use crate::llm_prompt::mentions_health;

    #[test]
    fn qc_banner_only_fires_above_threshold() {
        // A normal WGS sample (single-/low-double-digit novels) is silent.
        assert!(private_y_qc_banner(7).is_none());
        assert!(private_y_qc_banner(PRIVATE_Y_QC_WARN - 1).is_none());
        // At/above the threshold, warn and name the count + the likely causes.
        let w = private_y_qc_banner(PRIVATE_Y_QC_WARN).expect("should warn at threshold");
        assert!(w.contains(&PRIVATE_Y_QC_WARN.to_string()) && w.contains("contamination"));
        // And it appears in the rendered section.
        let fact = PrivateYFact {
            novel_unique: 400,
            off_path: 3,
            structural: 10,
        };
        let section = private_y_section(&Some(fact)).unwrap();
        assert!(section.contains("elevated private-Y count (400"));
    }

    fn base_brief() -> SubjectBrief {
        SubjectBrief {
            headline: Headline {
                name: "James".into(),
                test_chip: "Whole Genome Sequencing".into(),
                summary: "summary".into(),
            },
            paternal: Some(LineageBrief {
                kind: LineageKind::Paternal,
                haplogroup: "R-FGC29071".into(),
                lineage_path: vec!["R".into()],
                matched_ancestor: None,
                age_phrase: None,
                origin_phrase: None,
                story: None,
                confidence_phrase: "tentative placement".into(),
                sources: vec![],
            }),
            maternal: None,
            ancestry: Some(AncestryBrief {
                summary_phrase: "Predominantly European".into(),
                super_populations: vec![SuperPopulationSummary {
                    super_population: "European".into(),
                    percentage: 98.0,
                    populations: vec![],
                }],
                fine_pops: vec![],
                ancient_pops: vec![],
                interpretation: None,
                method_note: "estimated from 400,000 markers".into(),
            }),
            test: TestBrief {
                test_name: "Whole Genome Sequencing".into(),
                what_it_tells: "Reads your whole genome.".into(),
                limitations: None,
                quality_phrase: "high-quality (30× average depth)".into(),
                quality_ok: true,
            },
            needs_analysis: false,
            caveats: vec![],
            pack_version: None,
            pack_status: PackStatus::Bundled,
            enriched: false,
        }
    }

    fn full_context() -> ResultsContext {
        ResultsContext {
            brief: base_brief(),
            sex: Some(SexFact {
                label: "Male (XY)".into(),
                confidence: "high".into(),
            }),
            ystr: vec![
                YStrPanelFact { panel: "Y-111".into(), markers: 111 },
                YStrPanelFact { panel: "Y-37".into(), markers: 37 },
            ],
            private_y: Some(PrivateYFact { novel_unique: 12, off_path: 3, structural: 2 }),
            mt_mutations: Some(MtMutationsFact {
                total: 41,
                hvr1: 5,
                hvr2: 4,
                coding: 32,
                examples: vec!["263A>G".into(), "315.1C".into()],
            }),
            ibd: Some(IbdFact {
                match_count: 3,
                closest: Some(IbdMatchFact {
                    relationship: "2nd–3rd cousin".into(),
                    total_shared_cm: 210.0,
                    segment_count: 9,
                }),
            }),
        }
    }

    #[test]
    fn sheet_includes_brief_facts_and_every_present_signal() {
        let s = results_fact_sheet(&full_context());
        // Brief grounding still present (we extend, not replace).
        assert!(s.contains("R-FGC29071"));
        assert!(s.contains("Predominantly European"));
        // Each new signal section appears.
        assert!(s.contains("Genetic sex:"));
        assert!(s.contains("Male (XY) (high confidence)"));
        assert!(s.contains("Y-111 (111 markers)"));
        assert!(s.contains("Y-37 (37 markers)"));
        assert!(s.contains("12 novel variant(s) in unique sequence"));
        assert!(s.contains("3 variant(s) off known branches"));
        assert!(s.contains("41 total — HVR1 5, HVR2 4, coding 32"));
        assert!(s.contains("263A>G"));
        assert!(s.contains("Genetic relatives (IBD matches in this workspace): 3"));
        assert!(s.contains("about 210 cM shared across 9 segment(s), consistent with 2nd–3rd cousin"));
    }

    #[test]
    fn absent_signals_are_omitted() {
        let ctx = ResultsContext {
            brief: base_brief(),
            sex: None,
            ystr: vec![],
            private_y: None,
            mt_mutations: None,
            ibd: None,
        };
        let s = results_fact_sheet(&ctx);
        assert!(!s.contains("Genetic sex:"));
        assert!(!s.contains("Y-STR panels:"));
        assert!(!s.contains("Private Y variants"));
        assert!(!s.contains("mtDNA mutations"));
        assert!(!s.contains("Genetic relatives"));
        // It still degrades to exactly the brief fact sheet.
        assert_eq!(s, narrate_fact_sheet(&base_brief()));
    }

    #[test]
    fn structural_line_only_when_nonzero() {
        let mut ctx = full_context();
        ctx.private_y = Some(PrivateYFact { novel_unique: 4, off_path: 1, structural: 0 });
        let s = results_fact_sheet(&ctx);
        assert!(!s.contains("structural/paralog-prone"));
    }

    #[test]
    fn full_sheet_stays_clear_of_health_language() {
        // The curated grounding itself must not trip the post-generation health guard.
        assert!(!mentions_health(&results_fact_sheet(&full_context())));
    }

    #[test]
    fn signal_section_returns_one_signal_or_none() {
        let ctx = full_context();
        let ystr = signal_section(&ctx, SignalKind::YStr).unwrap();
        assert!(ystr.contains("Y-111 (111 markers)"));
        // It is just that signal — not the private-Y or sex blocks.
        assert!(!ystr.contains("Private Y variants"));
        assert!(!ystr.contains("Genetic sex"));

        let pvt = signal_section(&ctx, SignalKind::PrivateY).unwrap();
        assert!(pvt.contains("12 novel variant(s) in unique sequence"));
        assert!(!pvt.contains("Y-STR panels"));

        // Absent signal → None (this context has no IBD-less variant; build one).
        let empty = ResultsContext {
            brief: base_brief(),
            sex: None,
            ystr: vec![],
            private_y: None,
            mt_mutations: None,
            ibd: None,
        };
        assert!(signal_section(&empty, SignalKind::YStr).is_none());
        assert!(signal_section(&empty, SignalKind::Ibd).is_none());
    }
}
