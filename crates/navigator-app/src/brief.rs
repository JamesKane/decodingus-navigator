//! This module makes a [`SubjectBrief`] for a reader who is not a specialist.
//!
//! It reads the analysis signals of one subject, reads the narrative reference pack, and builds the
//! model that the UI draws. The template code in `navigator_domain::brief` does the last step, and
//! that code is pure.
//!
//! The module reads the reference pack in three steps, and each step has a fallback. The team
//! decided this on 2026-06-22. The app always holds a seed pack, which is the lowest step. A pack
//! on the CDN refreshes and extends the seed when the network permits it. An old cache covers a
//! failed refresh.
//!
//! An absent pack never stops a brief. Each section then falls back to the structured facts of the
//! analysis. [`SubjectBrief::pack_status`] records the age of the narrative.

use crate::{decodingus_appview_url, App, AppError};
use navigator_domain::ancestry::AncestryResult;
use navigator_domain::brief::{
    self, AncestryBrief, BriefPack, Headline, LineageBrief, LineageKind, PackStatus, RealignOffer, SubjectBrief,
    TestBrief,
};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::i18n::{self, tr, tr_fmt, Lang};
use navigator_domain::reconciliation::{CompatibilityLevel, Consensus, DnaType};
use navigator_domain::testtype::{self, TargetType};
use navigator_refgenome::cache as refgenome_cache;

/// The seed pack in the application bundle. It is the lowest step, and it works offline. The file
/// is `assets/brief-pack.seed.json`.
const SEED_PACK: &str = include_str!("../assets/brief-pack.seed.json");

/// Default CDN location of the refreshable reference pack. Override with `NAVIGATOR_BRIEF_PACK_URL`.
/// A 404 / unreachable host falls back gracefully to the cache, then the bundled seed.
const DEFAULT_BRIEF_PACK_URL: &str = "https://assets.decodingus.org/briefs/brief-pack.json";

/// The count of days that the app trusts a downloaded pack. After this time, the app tries a
/// refresh.
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
pub(crate) fn cache_is_fresh(path: &std::path::Path, ttl_days: u64) -> bool {
    let ttl = std::time::Duration::from_secs(ttl_days * 24 * 3600);
    std::fs::metadata(path)
        .and_then(|m| m.modified())
        .ok()
        .and_then(|mtime| std::time::SystemTime::now().duration_since(mtime).ok())
        .map(|age| age < ttl)
        .unwrap_or(false)
}

/// The count of days that the app trusts the extra record of one haplogroup. After this time, the
/// app tries a refresh.
const HAPLO_ENRICH_TTL_DAYS: u64 = 30;

/// Haplogroup content from the AppView. The cache key is the DNA type together with the name.
///
/// A `found = false` value marks an absent record. The endpoint answered, but it held nothing. So
/// the app does not request that haplogroup again at each rebuild.
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
    /// Shows that this record holds narrative content or age content for the brief.
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
    /// Build the plain-language brief for one subject.
    ///
    /// The method reads the consensus haplogroups, the coverage of the best alignment, and the test
    /// type of the run. It then joins these values to the reference pack.
    ///
    /// The method always returns a brief. A section with no data falls back to a simpler form. Only
    /// a store error stops the method.
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
        // The brief is prose for the reader, so the code builds it in the language of that
        // reader. The code finds the language here, and the caller does not supply it, because
        // each caller needs the same language. The UI draws the brief. The HTML export writes it
        // to a file that the user keeps. The local-LLM prompt gives it to a model, and that model
        // must answer in the language that the user reads.
        let lang = i18n::load_lang().unwrap_or(Lang::En);

        // The consensus lineages. A value is None when the app did not place the subject, or when
        // the test does not cover that lineage. The code tries to add content for each terminal
        // node from the haplogroup endpoint, and it caches the result. Offline, the pack values
        // apply.
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
            .map(|c| build_lineage(lang, LineageKind::Paternal, c, &pack, true, y_enrich.as_ref()));
        let maternal = cons_mt
            .as_ref()
            .map(|c| build_lineage(lang, LineageKind::Maternal, c, &pack, false, mt_enrich.as_ref()));

        let test = build_test(lang, test_code.as_deref(), coverage.as_ref(), &pack);

        // Ancestry composition (from the persisted consensus estimate; None for Y/mt-only tests).
        let ancestry = match self.donor_ancestry(biosample_guid).await? {
            Some((_, result)) if !result.super_population_summary.is_empty() => {
                let fine = self
                    .consensus_ancestry(biosample_guid, "FINE_ADMIXTURE")
                    .await
                    .ok()
                    .flatten();
                // The deep, or ancient, components. The code reads *only* `ANCIENT_ADMIXTURE`.
                // An older build wrote incorrect numbers to a `PCA_PROJECTION_GMM` row and a
                // `G25_NMONTE` row, and that fault caused this rebuild. A read of only one source
                // keeps those old rows out of three places. They are the brief, the DNA-story HTML
                // export, and the facts for the LLM.
                //
                // The value is absent when the three ancient sources can not express the ancestry
                // of the sample. No card is better than a wrong card.
                let ancient = if crate::ANCIENT_ANCESTRY_ENABLED {
                    self.consensus_ancestry(biosample_guid, navigator_analysis::ancestry::ANCIENT_ADMIXTURE)
                        .await
                        .ok()
                        .flatten()
                } else {
                    None
                };
                Some(build_ancestry(lang, &result, fine.as_ref(), ancient.as_ref(), &pack))
            }
            _ => None,
        };

        // The runs of homozygosity, which show relatedness and endogamy. The code only reads
        // them. It shows a value only when an earlier run calculated it and wrote it to the cache.
        // The brief must stay fast, and the ROH calculation runs only at the request of the
        // user.
        let roh = self.cached_roh(biosample_guid).await?.map(|r| {
            brief::roh_brief(
                lang,
                r.summary.pattern,
                r.summary.f_roh,
                r.summary.n_segments,
                r.summary.total_roh_mb,
                r.summary.longest_mb,
            )
        });

        // The archaic markers for Neanderthal. The rule is the same as the rule for ROH. The code
        // only reads them, and it shows a value only after an earlier run calculated the count and
        // wrote it to the cache. The brief must stay fast.
        let archaic = self.cached_archaic(biosample_guid).await?.map(|a| {
            brief::archaic_brief(
                lang,
                a.total_copies,
                a.possible_copies,
                a.called_sites,
                a.panel_sites,
                a.percentile,
                a.cohort.clone(),
            )
        });

        // Global caveats.
        let mut caveats = Vec::new();
        if matches!(pack_status, PackStatus::Bundled | PackStatus::Unavailable) {
            caveats.push(tr(lang, "brief.caveatOfflinePack").to_string());
        }
        if !test.quality_ok {
            caveats.push(tr(lang, "brief.caveatShallowDepth").to_string());
        }

        let headline = Headline {
            name: bio.donor_identifier.clone(),
            test_chip: test.test_name.clone(),
            summary: headline_summary(lang, &bio.donor_identifier, paternal.as_ref(), maternal.as_ref()),
        };

        // The code calculates this value first, because the next step moves `paternal` into the
        // brief.
        let realign_offer = self
            .realign_offer(biosample_guid, paternal.is_some(), default_aln.map(|(_, aln)| aln))
            .await?;

        Ok(SubjectBrief {
            headline,
            paternal,
            maternal,
            ancestry,
            roh,
            archaic,
            test,
            // Has a sequencing alignment but no coverage computed → offer the one-click Analyze.
            needs_analysis: default_aln.is_some() && coverage.is_none(),
            realign_offer,
            caveats,
            pack_version: (!pack.version.trim().is_empty()).then(|| pack.version.clone()),
            pack_status,
            enriched,
        })
    }

    /// Shows whether to offer a realignment to this subject, and which alignment to offer.
    ///
    /// Each condition below stops an offer of many hours of work that would change nothing.
    ///
    /// - **The subject has a paternal line to improve.** A realignment gives discovery on the Y
    ///   chromosome and almost nothing more. The ancestry, the IBD, and the autosomes already work
    ///   on GRCh37 and GRCh38, and they give the same answer on CHM13. With no Y placement there is
    ///   no gain, so the app makes no offer.
    /// - **The subject has reads to map again.** A subject with only a chip or a VCF has no
    ///   alignment. An alignment row with no file has no reads.
    /// - **The subject has no CHM13 alignment.** The offer states that the app can not read part of
    ///   the paternal line of the user. That statement is false for a user who already has data on
    ///   the complete assembly, by any route. An older file of that user can still be without a
    ///   realignment, and the statement stays false.
    /// - **The job would act on the reads.** Three tests apply. The alignment must not be on
    ///   CHM13. It must not be a realignment. It must not have a realignment already.
    ///
    ///   This code does not repeat those tests. [`crate::realign::realignable_for_subject`] holds
    ///   them, and the batch count and the job call the same function. An offer of work that the
    ///   job then refuses is worse than no offer.
    ///
    /// Among the alignments that pass, the code selects the default alignment of the subject. That
    /// is the test with the largest breadth, and then the largest depth, from
    /// [`Self::default_alignment_for_subject`]. If there is none, the code takes the first
    /// alignment. This choice matters for a user with a whole genome and a Y-only test. The
    /// realignment of the whole genome answers more questions.
    async fn realign_offer(
        &self,
        biosample_guid: SampleGuid,
        has_paternal: bool,
        preferred: Option<i64>,
    ) -> Result<Option<RealignOffer>, AppError> {
        if !has_paternal {
            return Ok(None);
        }

        let alignments = navigator_store::alignment::list_for_biosample(self.store.pool(), biosample_guid).await?;

        // The app already reads the Y chromosome of this person against the complete assembly.
        // The route to that state does not matter.
        //
        // The rule below, which looks at one alignment, can still find an old GRCh37 file with no
        // realignment. The *job* would act on that file. But the offer makes a promise about the
        // paternal line of the subject, not about one file. For this reader the app already keeps
        // that promise.
        //
        // A donor with four CHM13 alignments saw this fault. The app told that donor that it could
        // read nothing of the line of their father.
        if alignments
            .iter()
            .any(|a| crate::realign::is_target_build(&a.reference_build))
        {
            return Ok(None);
        }

        let eligible = crate::realign::realignable_for_subject(&alignments, crate::realign::DEFAULT_TARGET_BUILD);

        let chosen = preferred
            .filter(|id| eligible.contains(id))
            .or_else(|| eligible.first().copied());

        Ok(chosen.and_then(|id| {
            alignments.iter().find(|a| a.id == id).map(|a| RealignOffer {
                alignment_id: a.id,
                current_build: a.reference_build.clone(),
            })
        }))
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

        // The cache is new enough. Use it, and do not use the network.
        if let Some(cp) = &cached {
            if cache_is_fresh(&cache_path, BRIEF_PACK_TTL_DAYS) {
                pack.merge(cp.clone());
                return (pack, PackStatus::Cached);
            }
        }

        // The cache is old or absent. Try a refresh. If the refresh fails, use the old cache, and
        // then the seed pack.
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

    /// Read the extra content for one haplogroup, if the network permits it.
    ///
    /// The method reads the cache first, and the cache entry is valid for 30 days. If the cache has
    /// no entry, the method sends `GET {appview}/api/v1/haplogroup/{name}` with a short timeout.
    ///
    /// The method caches a definite answer, which is a 200 response or a 404 response. It caches
    /// the "not found" result also, so it does not send the request again at each rebuild. It does
    /// not cache a temporary network error, so the content appears after the network returns.
    ///
    /// The method returns content only when that content has an age or a narrative.
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
            // Network/timeout error → do not cache (retry next time).
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

/// Build a lineage section from the consensus and the pack content. The `enrich` value, when it is
/// present, replaces the age, the origin, and the story of the pack. `is_paternal` only selects the
/// lookup.
fn build_lineage(
    lang: Lang,
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
        age_phrase: brief::age_phrase(lang, formed_ybp),
        origin_phrase: brief::origin_phrase(lang, origin.as_deref()),
        story,
        confidence_phrase: brief::confidence_phrase(lang, c.confidence, c.run_count, conflict),
        sources,
    }
}

/// Assemble the ancestry section from the consensus estimate (+ optional fine-grained and ancient
/// estimates).
fn build_ancestry(
    lang: Lang,
    result: &AncestryResult,
    fine: Option<&AncestryResult>,
    ancient: Option<&AncestryResult>,
    pack: &BriefPack,
) -> AncestryBrief {
    use navigator_domain::ancestry::{population_color, population_name, population_super};
    use navigator_domain::brief::AncientComponent;

    let super_populations = result.super_population_summary.clone();
    let summary_phrase = brief::ancestry_summary(lang, &super_populations);
    let method_note = brief::ancestry_method_note(lang, result.snps_with_genotype, &result.panel_type);

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
                    // The pack content gives a clearer name and the explanation. The code looks
                    // up the pack by the code first, and then by the display name. So a plain code
                    // such as "ANF" becomes "Anatolian Farmer".
                    let direct = pack
                        .population(&c.population_code)
                        .or_else(|| pack.population(&c.population_name));
                    let name = direct
                        .and_then(|p| p.name.clone())
                        .unwrap_or_else(|| c.population_name.clone());
                    // The reference set of the model holds ancient populations and *modern*
                    // populations. A modern population usually has no text of its own. One example
                    // is a Colombian or Puerto Rican population, which represents Native American
                    // ancestry. So the code uses the continental description. Without it, a real
                    // signal from outside Europe has no explanation.
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
fn build_test(
    lang: Lang,
    test_code: Option<&str>,
    coverage: Option<&crate::CoverageResult>,
    pack: &BriefPack,
) -> TestBrief {
    let code = test_code.unwrap_or("");
    let test_name = if code.is_empty() {
        tr(lang, "brief.testUnknown").to_string()
    } else {
        testtype::display_name(code).to_string()
    };
    let target = testtype::by_code(code).map(|t| t.target).unwrap_or(TargetType::Mixed);

    // What it tells you + limits: pack description, else a target-derived fallback.
    let (what_it_tells, limitations) = match pack.test(code) {
        Some(e) => (e.what.clone(), e.limits.clone()),
        None => fallback_test_text(lang, target),
    };

    let (quality_phrase, quality_ok) = match coverage {
        Some(c) => brief::quality_phrase(lang, c.mean_coverage, target),
        None if code.starts_with("ARRAY") => brief::chip_quality_phrase(lang, 0),
        None => (tr(lang, "brief.depthNotMeasured").to_string(), false),
    };

    TestBrief {
        test_name,
        what_it_tells,
        limitations,
        quality_phrase,
        quality_ok,
    }
}

/// Plain-language test description when the pack does not cover the code, derived from what the test
/// targets.
fn fallback_test_text(lang: Lang, target: TargetType) -> (String, Option<String>) {
    let (what, limits) = match target {
        TargetType::WholeGenome => ("brief.testWholeGenome", None),
        TargetType::YChromosome => ("brief.testY", Some("brief.testYLimits")),
        TargetType::MtDna => ("brief.testMt", Some("brief.testMtLimits")),
        TargetType::Autosomal | TargetType::Mixed => ("brief.testAutosomal", None),
        TargetType::XChromosome => ("brief.testX", Some("brief.testXLimits")),
    };
    (tr(lang, what).to_string(), limits.map(|k| tr(lang, k).to_string()))
}

/// The one-line "who you are" headline summary.
fn headline_summary(
    lang: Lang,
    name: &str,
    paternal: Option<&LineageBrief>,
    maternal: Option<&LineageBrief>,
) -> String {
    match (paternal, maternal) {
        (Some(p), Some(m)) => tr_fmt(lang, "brief.headlineBoth", &[name, &p.haplogroup, &m.haplogroup]),
        (Some(p), None) => tr_fmt(lang, "brief.headlinePaternal", &[name, &p.haplogroup]),
        (None, Some(m)) => tr_fmt(lang, "brief.headlineMaternal", &[name, &m.haplogroup]),
        (None, None) => tr(lang, "brief.headlineNone").to_string(),
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
