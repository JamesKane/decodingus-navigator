//! Composition of a casual-reader [`SubjectBrief`]: pull the existing analysis signals for one
//! subject, load the narrative reference pack, and assemble the render-ready model via the pure
//! templating in `navigator_domain::brief`.
//!
//! The reference pack is loaded with **graceful fallback** (decided 2026-06-22): a bundled seed is
//! the always-available floor; a CDN-hosted pack refreshes/augments it when reachable; a stale cache
//! covers a failed refresh. A brief is never blocked by a missing pack — sections degrade to the
//! structured facts the analysis already provides, and [`SubjectBrief::pack_status`] records how
//! fresh the narrative is.

use crate::{App, AppError};
use navigator_domain::brief::{
    self, BriefPack, Headline, LineageBrief, LineageKind, PackStatus, SubjectBrief, TestBrief,
};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::reconciliation::{CompatibilityLevel, Consensus, DnaType};
use navigator_domain::testtype::{self, TargetType};
use navigator_refgenome::cache as refgenome_cache;

/// The bundled seed pack — the offline floor. Authored in `assets/brief-pack.seed.json`.
const SEED_PACK: &str = include_str!("../assets/brief-pack.seed.json");

/// Default CDN location of the refreshable reference pack. Override with `NAVIGATOR_BRIEF_PACK_URL`.
/// A 404 / unreachable host falls back gracefully to the cache, then the bundled seed.
const DEFAULT_BRIEF_PACK_URL: &str = "https://assets.decodingus.org/briefs/brief-pack.json";

/// How long a downloaded pack is trusted before a refresh is attempted (days).
const BRIEF_PACK_TTL_DAYS: u64 = 7;

fn brief_pack_url() -> String {
    std::env::var("NAVIGATOR_BRIEF_PACK_URL")
        .ok()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| DEFAULT_BRIEF_PACK_URL.to_string())
}

fn brief_pack_cache_path() -> std::path::PathBuf {
    refgenome_cache::base_dir().join("briefs").join("brief-pack.json")
}

/// Is the cached pack within its TTL? Unknown/unreadable mtime → not fresh (forces a refresh try).
fn brief_pack_is_fresh(path: &std::path::Path) -> bool {
    let ttl = std::time::Duration::from_secs(BRIEF_PACK_TTL_DAYS * 24 * 3600);
    std::fs::metadata(path)
        .and_then(|m| m.modified())
        .ok()
        .and_then(|mtime| std::time::SystemTime::now().duration_since(mtime).ok())
        .map(|age| age < ttl)
        .unwrap_or(false)
}

impl App {
    /// Build the plain-language brief for one subject. Pulls the consensus haplogroups, the best
    /// alignment's coverage, and the run's test type, and joins them to the reference pack. Always
    /// returns a brief (degrading per-section); only a store error propagates.
    pub async fn subject_brief(&self, biosample_guid: SampleGuid) -> Result<SubjectBrief, AppError> {
        let bio = navigator_store::biosample::get(self.store.pool(), biosample_guid)
            .await?
            .ok_or_else(|| AppError::Conflict(format!("unknown biosample {biosample_guid}")))?;

        // Best alignment + its run drive the test/quality section.
        let default_aln = self.default_alignment_for_subject(biosample_guid).await?;
        let (run, coverage) = match default_aln {
            Some((run_id, aln_id)) => {
                let run = navigator_store::sequence_run::get(self.store.pool(), run_id).await?;
                let coverage = self.cached_coverage(aln_id).await?;
                (run, coverage)
            }
            None => (None, None),
        };
        let test_code = run.as_ref().map(|r| r.test_type.clone());

        let (pack, pack_status) = self.load_brief_pack().await;

        // Consensus lineages (None when not placed yet, or N/A for the test).
        let cons_y = self.haplogroup_consensus(biosample_guid, DnaType::Y).await?;
        let cons_mt = self.haplogroup_consensus(biosample_guid, DnaType::Mt).await?;
        let paternal = cons_y
            .as_ref()
            .map(|c| build_lineage(LineageKind::Paternal, c, &pack, true));
        let maternal = cons_mt
            .as_ref()
            .map(|c| build_lineage(LineageKind::Maternal, c, &pack, false));

        let test = build_test(test_code.as_deref(), coverage.as_ref(), &pack);

        // Global caveats.
        let mut caveats = Vec::new();
        if matches!(pack_status, PackStatus::Bundled | PackStatus::Unavailable) {
            caveats.push(
                "Lineage descriptions are from the offline reference pack; connect to the internet for the latest."
                    .to_string(),
            );
        }
        if !test.quality_ok {
            caveats.push("This test's depth is limited, so some results are preliminary.".to_string());
        }

        let headline = Headline {
            name: bio.donor_identifier.clone(),
            test_chip: test.test_name.clone(),
            summary: headline_summary(&bio.donor_identifier, paternal.as_ref(), maternal.as_ref()),
        };

        Ok(SubjectBrief {
            headline,
            paternal,
            maternal,
            test,
            caveats,
            pack_version: (!pack.version.trim().is_empty()).then(|| pack.version.clone()),
            pack_status,
        })
    }

    /// Load the reference pack with graceful fallback: bundled seed (floor) → cached file (if fresh)
    /// → CDN refresh → stale cache. Never errors; the worst case is the seed (or an empty pack if
    /// even the seed fails to parse, flagged [`PackStatus::Unavailable`]).
    async fn load_brief_pack(&self) -> (BriefPack, PackStatus) {
        let (mut pack, mut status) = match serde_json::from_str::<BriefPack>(SEED_PACK) {
            Ok(p) => (p, PackStatus::Bundled),
            Err(e) => {
                eprintln!("brief: bundled seed pack failed to parse ({e}); descriptions unavailable");
                (BriefPack::default(), PackStatus::Unavailable)
            }
        };

        let cache_path = brief_pack_cache_path();
        let cached: Option<BriefPack> = std::fs::read_to_string(&cache_path)
            .ok()
            .and_then(|s| serde_json::from_str(&s).ok());

        // Fresh cache → use it without touching the network.
        if let Some(cp) = &cached {
            if brief_pack_is_fresh(&cache_path) {
                pack.merge(cp.clone());
                return (pack, PackStatus::Cached);
            }
        }

        // Stale / absent → try a refresh, falling back to the stale cache (then the seed).
        let url = brief_pack_url();
        let fetched: Result<BriefPack, String> = async {
            let resp = self
                .auth
                .http
                .get(&url)
                .send()
                .await
                .and_then(|r| r.error_for_status())
                .map_err(|e| format!("downloading {url}: {e}"))?;
            let body = resp.text().await.map_err(|e| format!("reading {url}: {e}"))?;
            serde_json::from_str::<BriefPack>(&body)
                .map(|p| (p, body))
                .map_err(|e| format!("parsing {url}: {e}"))
                .map(|(p, body)| {
                    if let Some(parent) = cache_path.parent() {
                        let _ = std::fs::create_dir_all(parent);
                    }
                    let _ = std::fs::write(&cache_path, &body);
                    p
                })
        }
        .await;

        match fetched {
            Ok(dp) => {
                pack.merge(dp);
                status = PackStatus::Downloaded;
            }
            Err(e) => {
                if let Some(cp) = cached {
                    eprintln!("brief: pack refresh failed ({e}); using the cached copy");
                    pack.merge(cp);
                    status = PackStatus::Cached;
                } else {
                    eprintln!("brief: pack refresh failed ({e}); using the bundled seed");
                }
            }
        }
        (pack, status)
    }
}

/// Assemble a lineage section from the consensus + pack content. `is_paternal` only chooses the
/// Y vs mt lookup.
fn build_lineage(kind: LineageKind, c: &Consensus, pack: &BriefPack, is_paternal: bool) -> LineageBrief {
    let matched = if is_paternal {
        pack.y_lookup(&c.haplogroup, &c.lineage)
    } else {
        pack.mt_lookup(&c.haplogroup, &c.lineage)
    };
    let (matched_name, entry) = match matched {
        Some((n, e)) => (Some(n), Some(e)),
        None => (None, None),
    };
    // Surface the ancestor only when the story is for an ancestor (not the terminal itself).
    let matched_ancestor = matched_name.filter(|n| n != &c.haplogroup);

    let conflict = matches!(
        c.compatibility,
        CompatibilityLevel::MajorDivergence | CompatibilityLevel::Incompatible
    );

    LineageBrief {
        kind,
        haplogroup: c.haplogroup.clone(),
        lineage_path: c.lineage.clone(),
        matched_ancestor,
        age_phrase: brief::age_phrase(entry.and_then(|e| e.formed_ybp)),
        origin_phrase: brief::origin_phrase(entry.and_then(|e| e.origin.as_deref())),
        story: entry.and_then(|e| e.story.clone()),
        confidence_phrase: brief::confidence_phrase(c.confidence, c.run_count, conflict),
        sources: entry.map(|e| e.sources.clone()).unwrap_or_default(),
    }
}

/// Assemble the test & quality section.
fn build_test(test_code: Option<&str>, coverage: Option<&crate::CoverageResult>, pack: &BriefPack) -> TestBrief {
    let code = test_code.unwrap_or("");
    let test_name = if code.is_empty() {
        "Unknown test".to_string()
    } else {
        testtype::display_name(code).to_string()
    };
    let target = testtype::by_code(code).map(|t| t.target).unwrap_or(TargetType::Mixed);

    // What it tells you + limits: pack description, else a target-derived fallback.
    let (what_it_tells, limitations) = match pack.test(code) {
        Some(e) => (e.what.clone(), e.limits.clone()),
        None => fallback_test_text(target),
    };

    let (quality_phrase, quality_ok) = match coverage {
        Some(c) => brief::quality_phrase(c.mean_coverage, target),
        None if code.starts_with("ARRAY") => brief::chip_quality_phrase(0),
        None => (
            "sequencing depth not yet measured — run analysis to see quality".to_string(),
            false,
        ),
    };

    TestBrief {
        test_name,
        what_it_tells,
        limitations,
        quality_phrase,
        quality_ok,
    }
}

/// Plain-language test description when the pack doesn't cover the code, derived from what the test
/// targets.
fn fallback_test_text(target: TargetType) -> (String, Option<String>) {
    match target {
        TargetType::WholeGenome => (
            "Reads your whole genome, covering your paternal line, maternal line and ancestry.".to_string(),
            None,
        ),
        TargetType::YChromosome => (
            "Tests your Y chromosome to determine your paternal line.".to_string(),
            Some("Covers only the Y chromosome — no maternal line or ancestry composition.".to_string()),
        ),
        TargetType::MtDna => (
            "Tests your mitochondrial DNA to determine your maternal line.".to_string(),
            Some("Covers only mitochondrial DNA — no paternal line or ancestry composition.".to_string()),
        ),
        TargetType::Autosomal | TargetType::Mixed => (
            "Genotypes markers across your genome for ancestry and a broad read of your lineages.".to_string(),
            None,
        ),
        TargetType::XChromosome => (
            "Tests your X chromosome.".to_string(),
            Some("Covers only the X chromosome.".to_string()),
        ),
    }
}

/// The one-line "who you are" headline summary.
fn headline_summary(name: &str, paternal: Option<&LineageBrief>, maternal: Option<&LineageBrief>) -> String {
    match (paternal, maternal) {
        (Some(p), Some(m)) => format!(
            "Your data places {name}'s paternal line at {} and maternal line at {}.",
            p.haplogroup, m.haplogroup
        ),
        (Some(p), None) => format!("Your data places {name}'s paternal line at {}.", p.haplogroup),
        (None, Some(m)) => format!("Your data places {name}'s maternal line at {}.", m.haplogroup),
        (None, None) => "Import or analyze this person's DNA to reveal their paternal and maternal lines.".to_string(),
    }
}
