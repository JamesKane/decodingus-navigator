//! Composition of a casual-reader [`SubjectBrief`]: pull the existing analysis signals for one
//! subject, load the narrative reference pack, and assemble the render-ready model via the pure
//! templating in `navigator_domain::brief`.
//!
//! The reference pack is loaded with **graceful fallback** (decided 2026-06-22): a bundled seed is
//! the always-available floor; a CDN-hosted pack refreshes/augments it when reachable; a stale cache
//! covers a failed refresh. A brief is never blocked by a missing pack — sections degrade to the
//! structured facts the analysis already provides, and [`SubjectBrief::pack_status`] records how
//! fresh the narrative is.

use crate::{decodingus_appview_url, App, AppError};
use navigator_domain::ancestry::AncestryResult;
use navigator_domain::brief::{
    self, AncestryBrief, BriefPack, Headline, LineageBrief, LineageKind, PackStatus, SubjectBrief, TestBrief,
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

/// Is the cached file within `ttl_days`? Unknown/unreadable mtime → not fresh (forces a refresh try).
fn cache_is_fresh(path: &std::path::Path, ttl_days: u64) -> bool {
    let ttl = std::time::Duration::from_secs(ttl_days * 24 * 3600);
    std::fs::metadata(path)
        .and_then(|m| m.modified())
        .ok()
        .and_then(|mtime| std::time::SystemTime::now().duration_since(mtime).ok())
        .map(|age| age < ttl)
        .unwrap_or(false)
}

/// How long a per-haplogroup enrichment record is trusted before a refresh is attempted (days).
const HAPLO_ENRICH_TTL_DAYS: u64 = 30;

/// Live haplogroup content fetched from the AppView, cached per (dna-type, name). `found = false` is
/// a negative-cache marker (the endpoint answered but had nothing) so a definitively-absent
/// haplogroup isn't re-requested every rebuild.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
struct HaploEnrichment {
    found: bool,
    #[serde(default)]
    formed_ybp: Option<i32>,
    #[serde(default)]
    tmrca_ybp: Option<i32>,
    #[serde(default)]
    origin: Option<String>,
    #[serde(default)]
    story: Option<String>,
}

impl HaploEnrichment {
    /// Does this carry any narrative/age content worth folding in?
    fn has_content(&self) -> bool {
        self.found && (self.formed_ybp.is_some() || self.origin.is_some() || self.story.is_some())
    }
}

fn haplo_enrich_cache_path(dna_type: DnaType, name: &str) -> std::path::PathBuf {
    // Sanitize the name for a filename (haplogroup names are SNP-ish but be defensive).
    let safe: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect();
    refgenome_cache::base_dir()
        .join("briefs")
        .join("haplo")
        .join(dna_type.as_str())
        .join(format!("{safe}.json"))
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

        // Consensus lineages (None when not placed yet, or N/A for the test). Each terminal is
        // enriched best-effort from the live haplogroup endpoint (cached); pack values stand offline.
        let cons_y = self.haplogroup_consensus(biosample_guid, DnaType::Y).await?;
        let cons_mt = self.haplogroup_consensus(biosample_guid, DnaType::Mt).await?;
        let mut enriched = false;
        let y_enrich = match &cons_y {
            Some(c) => self.enrich_haplogroup(&c.haplogroup, DnaType::Y).await,
            None => None,
        };
        let mt_enrich = match &cons_mt {
            Some(c) => self.enrich_haplogroup(&c.haplogroup, DnaType::Mt).await,
            None => None,
        };
        enriched |= y_enrich.is_some() || mt_enrich.is_some();
        let paternal = cons_y
            .as_ref()
            .map(|c| build_lineage(LineageKind::Paternal, c, &pack, true, y_enrich.as_ref()));
        let maternal = cons_mt
            .as_ref()
            .map(|c| build_lineage(LineageKind::Maternal, c, &pack, false, mt_enrich.as_ref()));

        let test = build_test(test_code.as_deref(), coverage.as_ref(), &pack);

        // Ancestry composition (from the persisted consensus estimate; None for Y/mt-only tests).
        let ancestry = match self.donor_ancestry(biosample_guid).await? {
            Some((_, result)) if !result.super_population_summary.is_empty() => {
                let fine = self
                    .consensus_ancestry(biosample_guid, "FINE_ADMIXTURE")
                    .await
                    .ok()
                    .flatten();
                // Ancient components are gated off (degenerate reference asset — see
                // [`crate::ANCIENT_ANCESTRY_ENABLED`]). Read nothing, so a *stale* row persisted by an
                // earlier build can't resurface in the brief, the DNA-story HTML export, or the LLM facts.
                let ancient = if crate::ANCIENT_ANCESTRY_ENABLED {
                    self.consensus_ancestry(biosample_guid, "PCA_PROJECTION_GMM")
                        .await
                        .ok()
                        .flatten()
                } else {
                    None
                };
                Some(build_ancestry(&result, fine.as_ref(), ancient.as_ref(), &pack))
            }
            _ => None,
        };

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
            ancestry,
            test,
            // Has a sequencing alignment but no coverage computed → offer the one-click Analyze.
            needs_analysis: default_aln.is_some() && coverage.is_none(),
            caveats,
            pack_version: (!pack.version.trim().is_empty()).then(|| pack.version.clone()),
            pack_status,
            enriched,
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
            if cache_is_fresh(&cache_path, BRIEF_PACK_TTL_DAYS) {
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

    /// Best-effort live enrichment for one haplogroup: cache-first (30-day TTL), else a short-timeout
    /// `GET {appview}/api/v1/haplogroup/{name}`. A definitive answer (200 / 404) is cached — including
    /// "not found" — so it isn't re-requested each rebuild; a transient network error is *not* cached,
    /// so enrichment self-heals once connectivity returns. Returns content only when there's something
    /// worth folding in (an age or narrative).
    async fn enrich_haplogroup(&self, name: &str, dna_type: DnaType) -> Option<HaploEnrichment> {
        if name.trim().is_empty() {
            return None;
        }
        let path = haplo_enrich_cache_path(dna_type, name);
        if cache_is_fresh(&path, HAPLO_ENRICH_TTL_DAYS) {
            if let Some(e) = std::fs::read_to_string(&path)
                .ok()
                .and_then(|s| serde_json::from_str::<HaploEnrichment>(&s).ok())
            {
                return e.has_content().then_some(e);
            }
        }

        let base = decodingus_appview_url();
        let url = format!("{base}/api/v1/haplogroup/{name}");
        let resp = self
            .auth
            .http
            .get(&url)
            .timeout(std::time::Duration::from_secs(4))
            .send()
            .await;

        let entry = match resp {
            Ok(r) if r.status().is_success() => {
                let body = r.text().await.unwrap_or_default();
                parse_haplo_enrichment(&body)
            }
            // The endpoint answered but had nothing (404 etc.) → cache a negative result.
            Ok(_) => HaploEnrichment::default(),
            // Network/timeout error → don't cache (retry next time).
            Err(_) => return None,
        };

        if let Some(parent) = path.parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        if let Ok(json) = serde_json::to_string(&entry) {
            let _ = std::fs::write(&path, json);
        }
        entry.has_content().then_some(entry)
    }
}

/// Parse the AppView haplogroup response into the enrichment subset. Tolerant of camelCase /
/// snake_case keys and a nested `provenance` blob; absent fields stay `None`.
fn parse_haplo_enrichment(body: &str) -> HaploEnrichment {
    let Ok(v) = serde_json::from_str::<serde_json::Value>(body) else {
        return HaploEnrichment::default();
    };
    let int = |keys: &[&str]| -> Option<i32> {
        keys.iter()
            .find_map(|k| v.get(*k).and_then(|x| x.as_i64()))
            .map(|n| n as i32)
    };
    let text = |keys: &[&str]| -> Option<String> {
        keys.iter()
            .find_map(|k| v.get(*k).and_then(|x| x.as_str()))
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
    };
    let formed_ybp = int(&["formed_ybp", "formedYbp"]);
    let tmrca_ybp = int(&["tmrca_ybp", "tmrcaYbp"]);
    let origin = text(&["origin"]).or_else(|| {
        v.get("provenance")
            .and_then(|p| p.get("origin"))
            .and_then(|x| x.as_str())
            .map(|s| s.to_string())
    });
    let story = text(&["story", "description", "summary"]);
    HaploEnrichment {
        found: true,
        formed_ybp,
        tmrca_ybp,
        origin,
        story,
    }
}

/// Assemble a lineage section from the consensus + pack content, overlaying live `enrich`ment when
/// present (it wins over pack values for age/origin/story). `is_paternal` only chooses the lookup.
fn build_lineage(
    kind: LineageKind,
    c: &Consensus,
    pack: &BriefPack,
    is_paternal: bool,
    enrich: Option<&HaploEnrichment>,
) -> LineageBrief {
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

    // Live enrichment wins over pack content for age/origin/story; pack fills the rest.
    let formed_ybp = enrich
        .and_then(|e| e.formed_ybp)
        .or_else(|| entry.and_then(|e| e.formed_ybp));
    let origin = enrich
        .and_then(|e| e.origin.clone())
        .or_else(|| entry.and_then(|e| e.origin.clone()));
    let story = enrich
        .and_then(|e| e.story.clone())
        .or_else(|| entry.and_then(|e| e.story.clone()));
    let mut sources = entry.map(|e| e.sources.clone()).unwrap_or_default();
    if enrich.is_some_and(|e| e.has_content()) {
        sources.push("DecodingUs (live)".to_string());
    }

    LineageBrief {
        kind,
        haplogroup: c.haplogroup.clone(),
        lineage_path: c.lineage.clone(),
        matched_ancestor,
        age_phrase: brief::age_phrase(formed_ybp),
        origin_phrase: brief::origin_phrase(origin.as_deref()),
        story,
        confidence_phrase: brief::confidence_phrase(c.confidence, c.run_count, conflict),
        sources,
    }
}

/// Assemble the ancestry section from the consensus estimate (+ optional fine-grained and ancient
/// estimates).
fn build_ancestry(
    result: &AncestryResult,
    fine: Option<&AncestryResult>,
    ancient: Option<&AncestryResult>,
    pack: &BriefPack,
) -> AncestryBrief {
    use navigator_domain::ancestry::{population_color, population_name, population_super};
    use navigator_domain::brief::AncientComponent;

    let super_populations = result.super_population_summary.clone();
    let summary_phrase = brief::ancestry_summary(&super_populations);
    let method_note = brief::ancestry_method_note(result.snps_with_genotype, &result.panel_type);

    let fine_pops = fine
        .map(|f| {
            f.components
                .iter()
                .map(|c| (c.population_name.clone(), c.percentage))
                .collect()
        })
        .unwrap_or_default();

    // Ancient components, biggest first, each with its palette color + a pack explanation (by code,
    // then display name).
    let ancient_pops: Vec<AncientComponent> = ancient
        .map(|a| {
            let mut comps: Vec<AncientComponent> = a
                .components
                .iter()
                .filter(|c| c.percentage >= 0.5)
                .map(|c| {
                    // Pack content (by code, then display name) supplies an optional friendlier name
                    // and the explanation — so a bare code like "ANF" reads as "Anatolian Farmer".
                    let direct = pack
                        .population(&c.population_code)
                        .or_else(|| pack.population(&c.population_name));
                    let name = direct
                        .and_then(|p| p.name.clone())
                        .unwrap_or_else(|| c.population_name.clone());
                    // The model's reference set mixes ancient and *modern* populations; the modern
                    // ones (e.g. Colombian/Puerto Rican standing in for Native American ancestry)
                    // rarely have their own blurb, so fall back to the continental description rather
                    // than leaving real non-European signal unexplained.
                    let blurb = direct.and_then(|p| p.blurb.clone()).or_else(|| {
                        population_super(&c.population_code)
                            .map(population_name)
                            .and_then(|sp| pack.population(&sp).and_then(|p| p.blurb.clone()))
                    });
                    AncientComponent {
                        code: c.population_code.clone(),
                        name,
                        percentage: c.percentage,
                        color: population_color(&c.population_code),
                        blurb,
                    }
                })
                .collect();
            comps.sort_by(|x, y| y.percentage.total_cmp(&x.percentage));
            comps
        })
        .unwrap_or_default();

    // Optional plain-language note for the dominant population (pack-supplied; tries code then name).
    let interpretation = super_populations
        .iter()
        .max_by(|a, b| a.percentage.total_cmp(&b.percentage))
        .and_then(|top| pack.population(&top.super_population).and_then(|p| p.blurb.clone()));

    AncestryBrief {
        summary_phrase,
        super_populations,
        fine_pops,
        ancient_pops,
        interpretation,
        method_note,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bundled_seed_pack_parses_and_has_content() {
        let pack: BriefPack = serde_json::from_str(SEED_PACK).expect("seed pack must be valid JSON");
        assert!(!pack.version.trim().is_empty());
        assert!(
            pack.y_haplogroups.contains_key("R-M269"),
            "expected a common Y haplogroup"
        );
        assert!(pack.mt_haplogroups.contains_key("H"), "expected a common mt haplogroup");
        assert!(pack.test_types.contains_key("WGS"), "expected the WGS test type");
        assert!(pack.populations.contains_key("European"), "expected a population blurb");
    }
}
