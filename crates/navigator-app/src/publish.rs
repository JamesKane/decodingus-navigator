//! `impl App` methods extracted from `lib.rs` (the `publish` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- publish -----------------------------------------------------------

    /// Build the alignment (coverage) record JSON for an alignment — the shared
    /// `com.decodingus.atmosphere.alignment` contract the AppView ingests (floats as strings).
    /// Links back to the subject's biosample + sequence-run records via their deterministic at://
    /// URIs in `did`'s repo, so the AppView can tie this coverage summary to its subject.
    pub(crate) async fn coverage_record(&self, did: &str, alignment_id: i64) -> Result<serde_json::Value, AppError> {
        let cov = self
            .cached_coverage(alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("coverage for alignment {alignment_id}"))))?;
        let aln = self.alignment_or_err(alignment_id).await?;
        // A whole-genome-labeled alignment whose reads are actually Y-scoped — a chrY-only extract,
        // or a Y test (Big Y / Y Elite) that came in mislabeled WGS — must not publish a coverage
        // summary. The AppView files an `alignment` record under whole-genome statistics, so its
        // near-zero autosomal depth and callable footprint would skew the aggregate WGS coverage
        // distributions. Genuine Y-targeted tests are exempt: their test type is published, so the
        // AppView cohorts their Y coverage separately from WGS.
        let is_wgs = matches!(
            sequence_run::get(self.store.pool(), aln.sequence_run_id)
                .await?
                .as_ref()
                .and_then(|r| navigator_domain::testtype::target_of(&r.test_type)),
            Some(navigator_domain::testtype::TargetType::WholeGenome)
        );
        if is_wgs
            && navigator_analysis::sex::is_y_scoped(
                cov.contig_coverage_stats.iter().map(|s| (s.contig.as_str(), s.num_reads)),
            )
        {
            return Err(AppError::Conflict(format!(
                "alignment {alignment_id} is a Y-scoped file labeled whole-genome — its coverage \
                 summary is withheld from the PDS so it can't skew AppView whole-genome statistics"
            )));
        }
        let guid = self.biosample_of_alignment(alignment_id).await?;
        let record = AlignmentRecord::new(
            aln.reference_build,
            Some(aln.aligner),
            cov.mean_coverage,
            cov.median_coverage,
            cov.sd_coverage,
            cov.pct_10x,
            cov.pct_20x,
            cov.pct_30x,
            cov.genome_territory,
            cov.callable_bases,
            Utc::now().to_rfc3339(),
        )
        .with_refs(
            Some(biosample_at_uri(did, guid)),
            Some(seqrun_at_uri(did, aln.sequence_run_id)),
        )
        .with_contigs(contig_metrics(&cov));
        Ok(serde_json::to_value(&record)?)
    }

    /// The subject's persisted **consensus** ancestry estimates ([`CONSENSUS_SOURCE_ID`]) — one per
    /// method (ADMIXTURE / FINE_ADMIXTURE), newest-first. Ancestry is estimated from the pooled
    /// autosomal consensus (not per alignment), so this is the subject's authoritative breakdown.
    /// Empty until the consensus ancestry has been estimated.
    ///
    /// The ancient methods (`PCA_PROJECTION_GMM`, `G25_NMONTE`) are **excluded while ancient ancestry
    /// is gated off** ([`crate::ANCIENT_ANCESTRY_ENABLED`]) — their reference asset is degenerate, and
    /// federating fabricated breakdowns to the PDS is far worse than showing them locally. The filter
    /// (not just the compute gate) is what keeps a *stale* row persisted by an earlier build from being
    /// published.
    pub(crate) async fn consensus_ancestry_results(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Vec<AncestryResult>, AppError> {
        const ANCIENT_METHODS: [&str; 2] = ["PCA_PROJECTION_GMM", "G25_NMONTE"];
        let all = ancestry_result::for_biosample(self.store.pool(), biosample_guid).await?;
        Ok(all
            .into_iter()
            .filter(|(id, _)| *id == CONSENSUS_SOURCE_ID)
            .map(|(_, r)| r)
            .filter(|r| crate::ANCIENT_ANCESTRY_ENABLED || !ANCIENT_METHODS.contains(&r.method.as_str()))
            .collect())
    }

    /// The populationBreakdown record JSON for each consensus ancestry estimate of a subject (one
    /// per method), linked to the biosample — the shared `com.decodingus.atmosphere.populationBreakdown`
    /// contract the AppView ingests (floats as strings). Empty if none computed.
    async fn consensus_ancestry_records(
        &self,
        did: &str,
        biosample_guid: SampleGuid,
    ) -> Result<Vec<serde_json::Value>, AppError> {
        let biosample_ref = biosample_at_uri(did, biosample_guid);
        self.consensus_ancestry_results(biosample_guid)
            .await?
            .iter()
            .map(|r| {
                let rec = population_breakdown_record(r).with_biosample_ref(Some(biosample_ref.clone()));
                serde_json::to_value(rec).map_err(AppError::from)
            })
            .collect()
    }

    /// Build the anonymized biosample record JSON — sex, center, and best-effort Y/mt
    /// haplogroup calls. Donor identifiers / accession / description are never carried.
    pub(crate) async fn biosample_record(
        &self,
        did: &str,
        biosample_guid: SampleGuid,
    ) -> Result<serde_json::Value, AppError> {
        let bio = biosample::get(self.store.pool(), biosample_guid)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("biosample {biosample_guid:?}"))))?;
        let y = self.consensus_haplogroup(biosample_guid, DnaType::Y).await?;
        let mt = self.consensus_haplogroup(biosample_guid, DnaType::Mt).await?;
        let runs = self.list_sequence_runs(biosample_guid).await?;
        // External identifiers (vendor kits + public catalog ids), a pure field rename onto the wire
        // shape. Published plaintext — the AppView keeps vendor ids off every public surface via its
        // `is_public` namespace policy; catalog ids (PGP/IGSR/ENA…) are already public. This is the
        // deterministic dedup anchor the AppView keys a re-published donor on.
        let external_ids = self
            .external_ids(biosample_guid)
            .await?
            .into_iter()
            .map(|e| du_domain::fed::ExternalId {
                namespace: e.source,
                value: e.external_id,
            })
            .collect();
        // Sequence-run refs are the runs' deterministic at:// URIs (not local ids), so the AppView
        // can follow them to the published sequence-run records.
        let record = BiosampleRecord::new(bio.sex, y, mt, bio.center_name, Utc::now().to_rfc3339())
            .with_refs(runs.iter().map(|r| seqrun_at_uri(did, r.id)).collect(), None, None)
            .with_external_ids(external_ids);
        Ok(serde_json::to_value(&record)?)
    }

    /// Build a sequence-run characterization record JSON (platform/instrument/test — no files).
    /// `instrument_id` (the sequencer serial inferred from read names) is published so the AppView
    /// can grow its crowd-sourced instrument→lab map (`fed.sequencerun.instrument_id` → the
    /// `instrument_observation`→proposal→accept consensus). It identifies the physical sequencer,
    /// not the donor — no PII, consistent with the anonymized fed-record posture.
    pub(crate) async fn sequence_run_record(
        &self,
        did: &str,
        run: &SequenceRun,
    ) -> Result<serde_json::Value, AppError> {
        let record = SequenceRunRecord::new(
            Some(biosample_at_uri(did, run.biosample_guid)),
            Some(run.platform_name.clone()),
            run.instrument_model.clone(),
            run.instrument_id.clone(),
            Some(run.test_type.clone()),
            run.library_layout.clone(),
            run.total_reads,
            run.mean_read_length.map(|l| l.round() as i32),
            run.mean_insert_size,
            Utc::now().to_rfc3339(),
        )
        // Publish the known lab so the AppView can display it (and learn the instrument→lab map —
        // many serials, e.g. PacBio, aren't in its dataset). See [`SequenceRun::sequencing_facility`].
        .with_facility(run.sequencing_facility.clone())
        // Exact sequenced yield + read chemistry back the standardized DTC test label the AppView
        // renders/groups by (`du_domain::testprofile`). Both `Option`al — older records omit them.
        .with_read_profile(run.total_bases, run.read_type.clone());
        Ok(serde_json::to_value(&record)?)
    }

    /// Best-effort consensus haplogroup for a subject arm, for the federated biosample record:
    /// manual override > genome-level placed terminal > per-run label reconciliation (all via
    /// [`haplogroup_consensus`](Self::haplogroup_consensus)). `None` when nothing has been called.
    async fn consensus_haplogroup(
        &self,
        biosample_guid: SampleGuid,
        dna_type: DnaType,
    ) -> Result<Option<String>, AppError> {
        Ok(self
            .haplogroup_consensus(biosample_guid, dna_type)
            .await?
            .map(|c| c.haplogroup))
    }

    /// Build the private-variants record JSON to publish for an alignment/contig.
    ///
    /// **chrY** publishes only the *filtered, publishable* private-Y set — the whole-chrY de-novo
    /// calls after backbone subtraction, callable masking, structural-region filtering, and the
    /// strict novel-marker [`PublishGate`] — tagged as unverified singleton candidates. It never
    /// publishes the raw de-novo flood (CHM13's Y is haplogroup J; an R sample's J-vs-R divergence
    /// plus paralog mismaps would otherwise drown AppView curators in non-viable SNPs).
    ///
    /// Other contigs (chrM) publish their raw de-novo calls — a small, well-behaved rCRS-relative
    /// set that needs no tree-relative filtering.
    pub(crate) async fn variants_record(&self, alignment_id: i64, contig: &str) -> Result<serde_json::Value, AppError> {
        let variants = if navigator_analysis::contig::is_chr_y(contig) {
            let bucket = self.private_y_variants_self_masked(alignment_id).await?;
            // QC gate: if the filtered novel count is implausibly high (contamination / low coverage /
            // reference-build mismatch — e.g. a GRCh38 alignment, whose chrY reference is far noisier
            // and whose shared-lineage variants the hs1-native tree can't fully resolve), the whole
            // set is suspect. Publish nothing rather than flood curators with candidates from a sample
            // we've already flagged; the variants still show in the in-app DISPLAY under the banner.
            if let Some(warn) = bucket.qc_banner() {
                eprintln!("private-variants publish skipped for alignment {alignment_id}: {warn}");
                Vec::new()
            } else {
                let gate = self.publish_gate_for_alignment(alignment_id).await?;
                bucket
                    .publishable(gate)
                    .into_iter()
                    .map(|v| {
                        VariantCallEntry::new(
                            v.position,
                            v.reference,
                            v.alternate,
                            v.depth,
                            v.alt_depth.min(v.depth),
                            v.allele_fraction,
                        )
                    })
                    .collect()
            }
        } else {
            let calls = self.cached_denovo(alignment_id, contig).await?.ok_or_else(|| {
                AppError::Store(StoreError::NotFound(format!(
                    "de-novo calls for alignment {alignment_id} {contig}"
                )))
            })?;
            calls
                .iter()
                .map(|c| {
                    VariantCallEntry::new(
                        c.position,
                        c.reference_allele,
                        c.alternate_allele,
                        c.depth,
                        c.alt_depth,
                        c.allele_fraction,
                    )
                })
                .collect()
        };
        let record = PrivateVariantsRecord::new(contig, caller::DENOVO_VERSION, Utc::now().to_rfc3339(), variants);
        Ok(serde_json::to_value(&record)?)
    }

    /// Publish an alignment's cached coverage summary using an explicit `client` (the
    /// testable core; production callers use [`publish_coverage`](Self::publish_coverage)).
    pub async fn publish_coverage_summary(&self, client: &PdsClient, alignment_id: i64) -> Result<RecordRef, AppError> {
        let value = self.coverage_record(client.did(), alignment_id).await?;
        Ok(client.create_record(NS_ALIGNMENT, value, None).await?)
    }

    /// Publish a subject's **consensus** ancestry estimates (one populationBreakdown per method)
    /// using an explicit `client` (the testable core; production callers use
    /// [`publish_ancestry`](Self::publish_ancestry)). Returns a ref per record.
    pub async fn publish_ancestry_with(
        &self,
        client: &PdsClient,
        biosample_guid: SampleGuid,
    ) -> Result<Vec<RecordRef>, AppError> {
        let mut refs = Vec::new();
        for value in self.consensus_ancestry_records(client.did(), biosample_guid).await? {
            refs.push(client.create_record(NS_POPULATION_BREAKDOWN, value, None).await?);
        }
        Ok(refs)
    }

    /// Publish the anonymized biosample summary using an explicit `client`.
    pub async fn publish_biosample_with(
        &self,
        client: &PdsClient,
        biosample_guid: SampleGuid,
    ) -> Result<RecordRef, AppError> {
        let value = self.biosample_record(client.did(), biosample_guid).await?;
        Ok(client.create_record(NS_BIOSAMPLE, value, Some(&biosample_rkey(biosample_guid))).await?)
    }

    /// Publish a sequence-run characterization using an explicit `client`.
    pub async fn publish_sequence_run_with(
        &self,
        client: &PdsClient,
        run: &SequenceRun,
    ) -> Result<RecordRef, AppError> {
        let value = self.sequence_run_record(client.did(), run).await?;
        Ok(client.create_record(NS_SEQUENCERUN, value, Some(&seqrun_rkey(run.id))).await?)
    }

    /// Publish an alignment's cached de-novo calls for `contig` using an explicit `client`
    /// (the testable core; production callers use [`publish_variants`](Self::publish_variants)).
    pub async fn publish_private_variants(
        &self,
        client: &PdsClient,
        alignment_id: i64,
        contig: &str,
    ) -> Result<RecordRef, AppError> {
        let value = self.variants_record(alignment_id, contig).await?;
        Ok(client.create_record(PRIVATE_VARIANTS_COLLECTION, value, None).await?)
    }
}

/// Fold a [`CoverageResult`]'s two per-contig views (samtools-style stats +
/// callable-state counts) into the shared lexicon's `contigs[]`, paired by contig
/// name — the same join `export::coverage_tsv` uses. Contigs present in the stats
/// but missing callable counts (shouldn't happen) fall back to zeros.
fn contig_metrics(cov: &CoverageResult) -> Vec<ContigMetrics> {
    cov.contig_coverage_stats
        .iter()
        .map(|s| {
            let c = cov.contig_callable.iter().find(|m| m.contig == s.contig);
            ContigMetrics {
                contig: s.contig.clone(),
                length: s.end_pos as i64,
                num_reads: s.num_reads as i64,
                mean_depth: s.mean_depth.into(),
                coverage_pct: s.coverage.into(),
                callable: c.map_or(0, |m| m.callable as i64),
                no_coverage: c.map_or(0, |m| m.no_coverage as i64),
                low_coverage: c.map_or(0, |m| m.low_coverage as i64),
                excessive_coverage: c.map_or(0, |m| m.excessive_coverage as i64),
                poor_mapping_quality: c.map_or(0, |m| m.poor_mapping_quality as i64),
                ref_n: c.map_or(0, |m| m.ref_n as i64),
                mean_base_q: s.mean_base_q.into(),
                mean_map_q: s.mean_map_q.into(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use navigator_analysis::coverage::{ContigCoverageStats, CoverageResult, COVERAGE_VERSION};
    use navigator_store::Store;

    fn cstat(contig: &str, num_reads: u64) -> ContigCoverageStats {
        ContigCoverageStats {
            contig: contig.into(),
            start_pos: 1,
            end_pos: 1,
            num_reads,
            cov_bases: 0,
            coverage: 0.0,
            mean_depth: 0.0,
            mean_base_q: 0.0,
            mean_map_q: 0.0,
            histogram: Vec::new(),
        }
    }

    /// chrY carrying millions of reads, autosomes/chrX only a trace of mismapped ones — the
    /// Y-only-extract shape.
    fn y_scoped_coverage() -> CoverageResult {
        CoverageResult {
            contig_coverage_stats: vec![cstat("chrY", 3_000_000), cstat("chr1", 30), cstat("chrX", 12)],
            ..Default::default()
        }
    }

    async fn alignment_with_test_type(app: &App, test_type: &str) -> i64 {
        let b = app.add_biosample(None, "yscoped", None, None).await.unwrap();
        let run = app
            .record_sequence_run(NewSequenceRun {
                biosample_guid: b.guid,
                platform_name: "ILLUMINA".into(),
                instrument_model: None,
                test_type: test_type.into(),
                library_layout: None,
                total_reads: None,
                pf_reads_aligned: None,
                mean_read_length: None,
                mean_insert_size: None,
            })
            .await
            .unwrap();
        app.record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "synthetic".into(),
            variant_caller: None,
            bam_path: Some("/nonexistent.cram".into()),
            reference_path: None,
            content_sha256: None,
        })
        .await
        .unwrap()
        .id
    }

    /// A WGS-labeled but Y-scoped alignment must not publish a coverage summary — it would poison
    /// the AppView's whole-genome statistics.
    #[tokio::test]
    async fn wgs_y_scoped_coverage_is_withheld() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let aln = alignment_with_test_type(&app, "WGS").await;
        app.save_analysis(aln, "coverage", COVERAGE_VERSION, &y_scoped_coverage())
            .await
            .unwrap();
        let err = app.coverage_record("did:plc:test", aln).await.unwrap_err();
        assert!(matches!(err, AppError::Conflict(_)), "expected Conflict, got {err:?}");
    }

    /// A Y-targeted test (Big Y) with the *same* Y-scoped shape publishes normally — its Y coverage
    /// is expected and the AppView cohorts it apart from WGS.
    #[tokio::test]
    async fn y_targeted_coverage_still_publishes() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let aln = alignment_with_test_type(&app, "BIG_Y_700").await;
        app.save_analysis(aln, "coverage", COVERAGE_VERSION, &y_scoped_coverage())
            .await
            .unwrap();
        app.coverage_record("did:plc:test", aln)
            .await
            .expect("a Y test's Y coverage must still publish");
    }

    /// A genuine whole-genome distribution publishes fine (the guard is specific to the Y-only shape).
    #[tokio::test]
    async fn normal_wgs_coverage_publishes() {
        let app = App::new(Store::open_in_memory().await.unwrap());
        let aln = alignment_with_test_type(&app, "WGS").await;
        let wgs = CoverageResult {
            contig_coverage_stats: vec![
                cstat("chr1", 200_000_000),
                cstat("chrX", 5_000_000),
                cstat("chrY", 3_000_000),
            ],
            ..Default::default()
        };
        app.save_analysis(aln, "coverage", COVERAGE_VERSION, &wgs).await.unwrap();
        app.coverage_record("did:plc:test", aln)
            .await
            .expect("normal WGS coverage should publish");
    }
}
