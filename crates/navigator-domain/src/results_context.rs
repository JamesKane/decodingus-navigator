//! The grounding context for the M4 "ask my results" chat (see
//! `documents/design/local-llm-expansion.md`). Narration (M1) keeps its grounding in the lean
//! [`SubjectBrief`](crate::brief::SubjectBrief) fact sheet. The chat must answer about *more*
//! signals: Y-STR panels, private-Y variants, mtDNA mutations, IBD matches, and genetic sex. So it
//! grounds in a [`ResultsContext`], which is the brief plus curated, summary-level facts for those
//! signals.
//!
//! Like [`llm_prompt`](crate::llm_prompt), this layer is pure. The facts for each signal are plain
//! values that the app already vetted and filled in. The domain crate must not depend on analysis
//! or the store. [`results_fact_sheet`] is the builder, and it has unit tests, so a person can
//! review the exact text we send.
//!
//! The model stays a *rewriter, not a source of facts*. Every section is summary-level, and the
//! inline notes keep STR, mtDNA and IBD as lineage facts, never as health or trait claims.

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

/// One Y-STR panel: its name, and how many markers it carries. Summary only, and never the raw
/// values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct YStrPanelFact {
    pub panel: String,
    pub markers: usize,
}

/// Private Y-variant counts, separated by confidence class (the `navigator-app` `PrivateBucket`
/// distinction: novel-in-unique-sequence vs off-path vs structural-region).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrivateYFact {
    /// Novel calls in unique sequence: the high-confidence new-branch candidates.
    pub novel_unique: usize,
    /// Variants off known branches. They suggest branch depth that nobody knew before.
    pub off_path: usize,
    /// Calls in structural regions, and regions where paralogs are common. Read them with
    /// caution.
    pub structural: usize,
}

/// Above this count of novel-in-unique private-Y calls, one sample almost certainly gives
/// artifacts, and not real new-branch candidates. In the de-novo tree pipeline, the novel count of
/// one WGS sample runs ~3–39 (median). So a count in the dozens is normal, and a count in the low
/// hundreds is a warning. The usual causes are contamination, shallow or uneven coverage, and a
/// mismatch of the reference build.
pub const PRIVATE_Y_QC_WARN: usize = 50;

/// A one-line QC banner when the novel-in-unique private-Y count is too high to be plausible for
/// one sample (see [`PRIVATE_Y_QC_WARN`]). `None` if not. Reports and the `private-y` CLI show it.
/// A high count then reads as "check this sample", and not as "you have this many new branches".
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

/// One IBD match, reduced to the relationship band and how much DNA the two share. It carries **no
/// detail that identifies anybody**. It does not name the partner, and that is deliberate. This is
/// the workspace data of the user, described as a relationship and not as a person.
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

/// The brief, plus curated summaries of the other signals: the grounding context for the M4 chat.
/// A signal that is absent is `None` or empty, and the fact sheet leaves it out. The model can
/// then not restate what is not there. This is exactly what the brief does with its own optional
/// sections.
#[derive(Debug, Clone, PartialEq)]
pub struct ResultsContext {
    pub brief: SubjectBrief,
    pub sex: Option<SexFact>,
    pub ystr: Vec<YStrPanelFact>,
    pub private_y: Option<PrivateYFact>,
    pub mt_mutations: Option<MtMutationsFact>,
    pub ibd: Option<IbdFact>,
}

/// Which result signal an "Explain this" narration on one tab (M5) points at. It is also the key
/// that takes the section of one signal out of a [`ResultsContext`], through [`signal_section`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum SignalKind {
    Sex,
    YStr,
    PrivateY,
    MtMutations,
    Ibd,
    Roh,
    Archaic,
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
            SignalKind::Roh => "shared ancestry (runs of homozygosity)",
            SignalKind::Archaic => "Neanderthal (archaic) markers",
        }
    }
}

// --- Section builders, one for each signal -------------------------------------------------------
// Each returns the labelled, summary-level block for one signal, or `None` when the subject has
// nothing for that signal. Each block starts with `\n`, so the blocks join cleanly after the fact
// sheet of the brief. The M4 chat grounding ([`results_fact_sheet`]) and the M5 narration on one
// tab ([`signal_section`]) share these blocks verbatim.

fn sex_section(sex: &Option<SexFact>) -> Option<String> {
    let sex = sex.as_ref()?;
    Some(format!(
        "\nGenetic sex:\n- {} ({} confidence)\n",
        sex.label, sex.confidence
    ))
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

/// Runs-of-homozygosity section, from the [`RohBrief`](crate::brief::RohBrief) of the brief. ROH is
/// a brief signal, and the others are not, so its source is `ctx.brief.roh`. The M5 "Explain this"
/// on the ROH card uses it. The M4 chat already carries ROH through the fact sheet of the brief
/// (see [`results_fact_sheet`]). So this does *not* add it there again, and that is deliberate.
fn roh_section(brief: &SubjectBrief) -> Option<String> {
    let r = brief.roh.as_ref()?;
    let mut s = String::from("\nShared ancestry (runs of homozygosity):\n");
    s.push_str(&format!("- pattern: {}\n", r.pattern));
    s.push_str(&format!(
        "- F_ROH: {:.4} (share of the genome in long identical runs)\n",
        r.f_roh
    ));
    s.push_str(&format!(
        "- {} run(s), about {:.0} Mb in total, longest {:.0} Mb\n",
        r.n_segments, r.total_mb, r.longest_mb
    ));
    s.push_str(&format!("- {}\n", r.summary_phrase));
    s.push_str(
        "- note: this describes shared ancestry between the parents' lines (endogamy / consanguinity), \
         not a health or trait result.\n",
    );
    Some(s)
}

fn archaic_section(brief: &SubjectBrief) -> Option<String> {
    let a = brief.archaic.as_ref()?;
    let mut s = String::from("\nNeanderthal (archaic) markers:\n");
    s.push_str(&format!(
        "- {} archaic-allele copies out of {} assayed\n",
        a.total_copies, a.possible_copies
    ));
    s.push_str(&format!(
        "- measured at {} of the {} marker sites in the panel ({:.1}% of it)\n",
        a.called_sites,
        a.panel_sites,
        if a.panel_sites == 0 {
            0.0
        } else {
            a.called_sites as f64 * 100.0 / a.panel_sites as f64
        }
    ));
    match (a.percentile, &a.cohort) {
        (Some(p), Some(c)) => s.push_str(&format!("- more than {p:.0}% of {c} reference samples\n")),
        _ => s.push_str("- no population comparison: the test covers too few marker sites to rank\n"),
    }
    s.push_str(&format!("- pattern: {}\n", a.pattern));
    // The model must not turn a count into a percentage, and must not make a Denisovan claim. The
    // design puts both out of scope for this signal (S1 and S7), and this text is the only
    // grounding the narration has.
    s.push_str(
        "- note: this is a COUNT of marker copies, not a percentage of the genome, and it is specific \
         to this panel — do not compare it to another company's number, and do not restate it as a \
         percent Neanderthal.\n",
    );
    s.push_str(
        "- note: no Denisovan result is reported; outside Oceania that signal is at the noise floor. \
         This is a genealogical curiosity, not a health or trait result.\n",
    );
    Some(s)
}

/// The labelled section for one signal, or `None` when the subject has nothing for it. This is the
/// grounding for an M5 "Explain this" narration of that one signal.
pub fn signal_section(ctx: &ResultsContext, kind: SignalKind) -> Option<String> {
    match kind {
        SignalKind::Sex => sex_section(&ctx.sex),
        SignalKind::YStr => ystr_section(&ctx.ystr),
        SignalKind::PrivateY => private_y_section(&ctx.private_y),
        SignalKind::MtMutations => mt_section(&ctx.mt_mutations),
        SignalKind::Ibd => ibd_section(&ctx.ibd),
        SignalKind::Roh => roh_section(&ctx.brief),
        SignalKind::Archaic => archaic_section(&ctx.brief),
    }
}

/// Build the grounding fact sheet of the chat. It is the fact sheet of the brief, which narration
/// does not change, plus a labelled, summary-level block for each signal that is there. Pure, and
/// unit tested. This is the exact text that goes into the chat system message as "your only source
/// of facts".
pub fn results_fact_sheet(ctx: &ResultsContext) -> String {
    let mut s = narrate_fact_sheet(&ctx.brief);
    for section in [
        sex_section(&ctx.sex),
        ystr_section(&ctx.ystr),
        private_y_section(&ctx.private_y),
        mt_section(&ctx.mt_mutations),
        ibd_section(&ctx.ibd),
        // ROH already reaches the sheet through narrate_fact_sheet, because it lives on the
        // brief. The archaic block does not, so add it here, or the chat can not answer about
        // it.
        archaic_section(&ctx.brief),
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
    use crate::brief::{
        AncestryBrief, Headline, LineageBrief, LineageKind, PackStatus, RohBrief, SubjectBrief, TestBrief,
    };
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
            roh: None,
            archaic: None,
            test: TestBrief {
                test_name: "Whole Genome Sequencing".into(),
                what_it_tells: "Reads your whole genome.".into(),
                limitations: None,
                quality_phrase: "high-quality (30× average depth)".into(),
                quality_ok: true,
            },
            needs_analysis: false,
            realign_offer: None,
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
                YStrPanelFact {
                    panel: "Y-111".into(),
                    markers: 111,
                },
                YStrPanelFact {
                    panel: "Y-37".into(),
                    markers: 37,
                },
            ],
            private_y: Some(PrivateYFact {
                novel_unique: 12,
                off_path: 3,
                structural: 2,
            }),
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
    fn roh_is_answerable_in_chat_and_per_tab() {
        let mut ctx = full_context();
        ctx.brief.roh = Some(RohBrief {
            f_roh: 0.0080,
            pattern: "Outbred".into(),
            summary_phrase: "Your parents' lines don't share a recent common ancestor.".into(),
            n_segments: 6,
            total_mb: 22.7,
            longest_mb: 7.1,
        });

        // M4 chat: the fact sheet carries the ROH facts, through the brief, with the run counts.
        let sheet = results_fact_sheet(&ctx);
        assert!(sheet.contains("Shared ancestry (runs of homozygosity)"));
        assert!(sheet.contains("pattern: Outbred"));
        assert!(sheet.contains("F_ROH: 0.0080"));
        assert!(sheet.contains("6 run(s)"));

        // M5 on one tab: the focused section draws, and stays on ancestry (no health language).
        let section = signal_section(&ctx, SignalKind::Roh).expect("roh section");
        assert!(section.contains("F_ROH: 0.0080"));
        assert!(section.contains("longest 7 Mb"));
        assert!(!mentions_health(&section), "ROH must not read as a health result");

        // Absent until something computes ROH.
        ctx.brief.roh = None;
        assert!(signal_section(&ctx, SignalKind::Roh).is_none());
    }

    #[test]
    fn archaic_section_is_a_count_never_a_percentage_or_a_denisovan_claim() {
        let mut ctx = full_context();
        ctx.brief.archaic = Some(crate::brief::ArchaicBrief {
            total_copies: 12126,
            possible_copies: 599864,
            called_sites: 299932,
            panel_sites: 299958,
            percentile: Some(12.0),
            cohort: Some("EUR".into()),
            pattern: "Fewer than most".into(),
            summary_phrase: "You carry fewer Neanderthal markers than most people with similar ancestry.".into(),
        });

        let section = signal_section(&ctx, SignalKind::Archaic).expect("archaic section");
        assert!(section.contains("12126 archaic-allele copies out of 599864"));
        assert!(section.contains("more than 12% of EUR"));

        // The grounding must push the model off the two forms the design forbids. It must not
        // restate a count as a percent-Neanderthal, and it must not report a Denisovan result.
        assert!(
            section.contains("not a percentage"),
            "must warn against percent framing"
        );
        assert!(section.contains("no Denisovan result is reported"));
        assert!(!mentions_health(&section), "archaic must not read as a health result");

        // And it must not itself state a percent-of-genome figure anywhere.
        assert!(!section.contains('%') || section.contains("% of EUR") || section.contains("% of it"));

        // The chat fact sheet carries it too.
        assert!(results_fact_sheet(&ctx).contains("Neanderthal (archaic) markers"));

        // Absent until something computes the count.
        ctx.brief.archaic = None;
        assert!(signal_section(&ctx, SignalKind::Archaic).is_none());
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
        ctx.private_y = Some(PrivateYFact {
            novel_unique: 4,
            off_path: 1,
            structural: 0,
        });
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
        // It is only that signal, and not the private-Y block or the sex block.
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
