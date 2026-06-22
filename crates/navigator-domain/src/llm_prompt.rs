//! Pure prompt construction + grounding for the local-LLM narration (see
//! `docs/design/local-llm-integration.md`). No I/O — given a [`SubjectBrief`], produce the exact
//! `system`/`user` message text we send, so the guardrails are reviewable and unit-tested (the same
//! discipline as the deterministic brief templating). The model is a **rewriter, not a source of
//! facts**: it only restates the already-curated, already-rounded strings in the fact sheet.

use crate::brief::{LineageBrief, SubjectBrief};

/// The fixed grounding/guardrail instructions for brief narration. Returned as an owned `String` so
/// callers (and tests) see the literal text we send.
pub fn narrate_system_prompt() -> String {
    "You are a genetic-genealogy guide writing a warm, insightful summary for a curious non-expert. \
     Your goal is to help them UNDERSTAND their results — not to repeat them. Weave the paternal \
     line, maternal line, and ancestry into one connected story; explain what the findings mean and \
     why they are interesting; note anything distinctive; and briefly define any term a beginner \
     wouldn't know (e.g. what a haplogroup is). Do NOT simply restate or list each fact one by one.\n\
     Stay grounded in the facts given in the user message. You may interpret, connect, and add \
     general context that follows directly from those facts, but do NOT introduce specific new \
     claims — no haplogroup ages, place names, peoples, dates, or percentages that are not in the \
     facts. NEVER make health, medical, disease, trait, or clinical statements; this covers ancestry \
     and lineage only. Preserve the stated confidence: if a placement is described as \"tentative\", \
     keep it tentative, and never overstate certainty. If something is not in the facts, leave it out \
     rather than guessing.\n\
     Write two to four short paragraphs in second person (\"your paternal line\"). Open with a single \
     sentence capturing who they are genetically, then move from deep ancestry toward their specific \
     lineages, and close with what this test can and cannot tell them. Do not address the reader as a \
     patient. No preamble, headings, bullet lists, or sign-off — return only the prose."
        .to_string()
}

/// The same facts-only / no-health / preserve-uncertainty rules, adapted for the M2 Q&A chat (added
/// here so both prompts share one reviewed guardrail source).
pub fn answer_system_prompt() -> String {
    format!(
        "{}\n\nYou are answering a question about these results. If the answer is not in the provided \
         context, say you don't know rather than guessing. If asked for medical, health, disease, or \
         clinical interpretation, reply that Navigator covers ancestry and lineage, not health.",
        narrate_system_prompt()
    )
}

fn lineage_lines(out: &mut String, label: &str, lb: &LineageBrief) {
    out.push_str(&format!("\n{label}:\n"));
    out.push_str(&format!("- haplogroup: {}\n", lb.haplogroup));
    if let Some(anc) = &lb.matched_ancestor {
        out.push_str(&format!("- description is for the ancestor: {anc}\n"));
    }
    if let Some(age) = &lb.age_phrase {
        out.push_str(&format!("- age: {age}\n"));
    }
    if let Some(origin) = &lb.origin_phrase {
        out.push_str(&format!("- origin: {origin}\n"));
    }
    if let Some(story) = &lb.story {
        out.push_str(&format!("- background: {story}\n"));
    }
    out.push_str(&format!("- confidence: {}\n", lb.confidence_phrase));
}

/// The user-message fact sheet built from the brief — only already-curated strings from the
/// deterministic pipeline. A missing section is simply absent (so the model can't restate it).
pub fn narrate_fact_sheet(b: &SubjectBrief) -> String {
    let mut s = String::from("FACTS:\n");
    s.push_str(&format!("Name: {}\n", b.headline.name));
    s.push_str(&format!("Test: {}\n", b.test.test_name));

    if let Some(p) = &b.paternal {
        lineage_lines(&mut s, "Paternal line (Y-DNA)", p);
    }
    if let Some(m) = &b.maternal {
        lineage_lines(&mut s, "Maternal line (mtDNA)", m);
    }

    if let Some(a) = &b.ancestry {
        s.push_str("\nAncestry:\n");
        s.push_str(&format!("- summary: {}\n", a.summary_phrase));
        for sp in a.super_populations.iter().filter(|p| p.percentage >= 0.5) {
            s.push_str(&format!("- {}: {:.1}%\n", sp.super_population, sp.percentage));
        }
        if let Some(interp) = &a.interpretation {
            s.push_str(&format!("- note: {interp}\n"));
        }
        for c in a.ancient_pops.iter() {
            s.push_str(&format!("- ancient source {} ({:.1}%)", c.name, c.percentage));
            match &c.blurb {
                Some(blurb) => s.push_str(&format!(": {blurb}\n")),
                None => s.push('\n'),
            }
        }
    }

    s.push_str("\nTest quality:\n");
    s.push_str(&format!("- {}\n", b.test.what_it_tells));
    if let Some(lim) = &b.test.limitations {
        s.push_str(&format!("- limitation: {lim}\n"));
    }
    s.push_str(&format!("- quality: {}\n", b.test.quality_phrase));

    s
}

/// The fixed reply for a health/medical question (the M2 scope guard) — keeps the assistant in the
/// ancestry/lineage lane instead of attempting a clinical answer.
pub fn health_deflection() -> &'static str {
    "I can only help with ancestry and lineage — Navigator doesn't provide health, medical, or \
     clinical interpretation. Ask me about your paternal or maternal line, your ancestry, or your test."
}

/// Conservative post-generation guard: does `text` contain clearly health/clinical language? Used to
/// reject (and fall back from) any model output that strays out of the genealogy/ancestry lane.
pub fn mentions_health(text: &str) -> bool {
    let t = text.to_lowercase();
    const TERMS: &[&str] = &[
        "disease",
        "diagnos", // diagnosis / diagnose / diagnostic
        "cancer",
        "tumor",
        "tumour",
        "clinical",
        "symptom",
        "prognos",
        "medication",
        "treatment",
        "medical condition",
        "health risk",
        "carrier status",
        "pathogenic",
    ];
    TERMS.iter().any(|w| t.contains(w))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ancestry::SuperPopulationSummary;
    use crate::brief::{
        AncestryBrief, AncientComponent, Headline, LineageBrief, LineageKind, PackStatus, SubjectBrief, TestBrief,
    };

    fn sample_brief() -> SubjectBrief {
        SubjectBrief {
            headline: Headline {
                name: "James".into(),
                test_chip: "Whole Genome Sequencing".into(),
                summary: "summary".into(),
            },
            paternal: Some(LineageBrief {
                kind: LineageKind::Paternal,
                haplogroup: "R-FGC29071".into(),
                lineage_path: vec!["R".into(), "R-M269".into()],
                matched_ancestor: Some("R-M269".into()),
                age_phrase: Some("formed roughly 6,400 years ago".into()),
                origin_phrase: Some("associated with the steppe".into()),
                story: Some("A Western European lineage.".into()),
                confidence_phrase: "tentative placement, from a single test".into(),
                sources: vec!["YFull".into()],
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
                ancient_pops: vec![AncientComponent {
                    code: "WHG".into(),
                    name: "Western Hunter-Gatherer".into(),
                    percentage: 50.0,
                    color: "#e15759".into(),
                    blurb: Some("Europe's oldest layer.".into()),
                }],
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
            caveats: vec![],
            pack_version: Some("2026.06-seed.2".into()),
            pack_status: PackStatus::Bundled,
            enriched: false,
        }
    }

    #[test]
    fn system_prompt_carries_the_guardrails() {
        let p = narrate_system_prompt();
        assert!(p.contains("grounded in the facts"));
        assert!(p.contains("do NOT introduce specific new"));
        assert!(p.contains("NEVER make health"));
        assert!(p.contains("tentative"));
        assert!(answer_system_prompt().contains("say you don't know"));
    }

    #[test]
    fn fact_sheet_contains_brief_fields_and_preserves_confidence() {
        let s = narrate_fact_sheet(&sample_brief());
        assert!(s.contains("Name: James"));
        assert!(s.contains("R-FGC29071"));
        assert!(s.contains("ancestor: R-M269"));
        assert!(s.contains("formed roughly 6,400 years ago"));
        assert!(s.contains("tentative placement"), "confidence must survive");
        assert!(s.contains("Predominantly European"));
        assert!(s.contains("Western Hunter-Gatherer"));
        assert!(s.contains("high-quality (30× average depth)"));
    }

    #[test]
    fn fact_sheet_omits_absent_sections() {
        let s = narrate_fact_sheet(&sample_brief());
        // No maternal line in the sample brief → it must not appear.
        assert!(!s.contains("Maternal line"));
    }

    #[test]
    fn health_guard_flags_clinical_language() {
        assert!(mentions_health("This suggests a higher disease risk."));
        assert!(mentions_health("a pathogenic variant"));
        assert!(!mentions_health("Your paternal line is common in Western Europe."));
    }
}
