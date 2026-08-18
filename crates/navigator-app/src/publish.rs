//! `impl App` methods extracted from `lib.rs` (the `publish` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- publish -----------------------------------------------------------

    /// Build the JSON of the alignment record, which holds the coverage. The record follows the
    /// shared `com.decodingus.atmosphere.alignment` contract that the AppView reads, and each float
    /// is a string.
    ///
    /// The record links to the biosample record and the sequence-run records of the subject. It
    /// uses their fixed at:// URIs in the repository of `did`. So the AppView can join this
    /// coverage summary to its subject.
    pub(crate) async fn coverage_record(&self, did: &str, alignment_id: i64) -> Result<serde_json::Value, AppError> {
        let cov = self
            .cached_coverage(alignment_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("coverage for alignment {alignment_id}"))))?;
        let aln = self.alignment_or_err(alignment_id).await?;
        // An alignment can carry the label of a whole genome while its reads cover only the Y
        // chromosome. The file can be a chrY extract, or a Y test such as Big Y or Y Elite with a
        // wrong WGS label. The app must not publish a coverage summary for such an alignment.
        //
        // The AppView files an `alignment` record under the statistics of a whole genome. The
        // autosomal depth and the callable area of these files are almost zero. Those values move
        // the WGS coverage distribution of the full cohort.
        //
        // A true Y test does not have this problem. The app publishes its test type, so the AppView
        // puts its Y coverage in a group that is separate from WGS.
        let is_wgs = matches!(
            sequence_run::get(self.store.pool(), aln.sequence_run_id)
                .await?
                .as_ref()
                .and_then(|r| navigator_domain::testtype::target_of(&r.test_type)),
            Some(navigator_domain::testtype::TargetType::WholeGenome)
        );
        if is_wgs
            && navigator_analysis::sex::is_y_scoped(
                cov.contig_coverage_stats
                    .iter()
                    .map(|s| (s.contig.as_str(), s.num_reads)),
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

    /// The stored **consensus** ancestry estimates of the subject ([`CONSENSUS_SOURCE_ID`]). There
    /// is one estimate for each method, which is ADMIXTURE or FINE_ADMIXTURE, and the newest comes
    /// first.
    ///
    /// The code estimates the ancestry from the pooled autosomal consensus, and not from one
    /// alignment. So this list is the breakdown of the subject with authority. The list is empty
    /// until an estimate runs.
    ///
    /// Two filters control what leaves the machine. A wrong breakdown on the network is much worse
    /// than a wrong breakdown on the screen. A PDS record stays after the app corrects the fault.
    ///
    /// * The app **never** publishes a method in `RETIRED_METHODS`, with a flag or without one. The
    ///   ancient estimators that used a PCA centroid gave incorrect breakdowns. Those estimators are
    ///   gone, but a database from before the rebuild still holds their rows. This filter makes sure
    ///   that a build which can no longer *make* those numbers can not *publish* them.
    /// * The app publishes the current ancient method, `ANCIENT_ADMIXTURE`, only while ancient
    ///   ancestry is on, in [`crate::ANCIENT_ANCESTRY_ENABLED`]. That flag is then a true switch to
    ///   stop the feature.
    pub(crate) async fn consensus_ancestry_results(
        &self,
        biosample_guid: SampleGuid,
    ) -> Result<Vec<AncestryResult>, AppError> {
        const RETIRED_METHODS: [&str; 2] = ["PCA_PROJECTION_GMM", "G25_NMONTE"];
        let ancient = navigator_analysis::ancestry::ANCIENT_ADMIXTURE;
        let all = ancestry_result::for_biosample(self.store.pool(), biosample_guid).await?;
        Ok(all
            .into_iter()
            .filter(|(id, _)| *id == CONSENSUS_SOURCE_ID)
            .map(|(_, r)| r)
            .filter(|r| !RETIRED_METHODS.contains(&r.method.as_str()))
            .filter(|r| crate::ANCIENT_ANCESTRY_ENABLED || r.method != ancient)
            .collect())
    }

    /// The JSON of a populationBreakdown record for each consensus ancestry estimate of a subject.
    /// There is one record for each method, and each record links to the biosample.
    ///
    /// The record follows the shared `com.decodingus.atmosphere.populationBreakdown` contract that
    /// the AppView reads, and each float is a string. The list is empty when no estimate exists.
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

    /// Build the JSON of the anonymous biosample record. It holds the sex, the center, and the Y
    /// and mt haplogroup calls when they exist. The record never holds a donor identifier, an
    /// accession, or a description.
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
        // The external identifiers, which are the vendor kits and the public catalog ids. This
        // step only renames the fields for the wire format.
        //
        // The app publishes these values as plaintext. The `is_public` namespace policy of the
        // AppView keeps a vendor id off each public screen. A catalog id from PGP, IGSR, or ENA is
        // already public.
        //
        // These identifiers are the fixed anchor that the AppView uses to find a duplicate when a
        // user publishes the same donor again.
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

    /// Build the JSON of a sequence-run record. It holds the platform, the instrument, and the
    /// test. It holds no file.
    ///
    /// The app publishes `instrument_id`, which is the serial number of the sequencer that the code
    /// deduces from the read names. The AppView uses that value to build its map from an instrument
    /// to a laboratory. The path is `fed.sequencerun.instrument_id`, then an
    /// `instrument_observation`, then a proposal, and then the accepted consensus.
    ///
    /// The value names the physical sequencer. It does not name the donor. It holds no personal
    /// data, and it follows the rule for each anonymous federated record.
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
        // Publish the laboratory when the app knows it. The AppView then shows it, and the
        // AppView also learns the map from an instrument to a laboratory. Its dataset holds no
        // entry for many serial numbers, such as the PacBio numbers. See
        // [`SequenceRun::sequencing_facility`].
        .with_facility(run.sequencing_facility.clone())
        // The exact yield and the read chemistry support the standard DTC test label. The AppView
        // draws that label and groups by it, in `du_domain::testprofile`. Both fields are
        // `Option`, because an older record holds neither.
        .with_read_profile(run.total_bases, run.read_type.clone());
        Ok(serde_json::to_value(&record)?)
    }

    /// The consensus haplogroup of one lineage of a subject, for the federated biosample record.
    ///
    /// The method takes the first value that exists, in this order. First, a manual value from the
    /// user. Second, the terminal node of the genome-level placement. Third, the reconciled label
    /// of each run.
    ///
    /// [`haplogroup_consensus`](Self::haplogroup_consensus) gives all three values. The method
    /// returns `None` when no call exists.
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
    /// For **chrY**, the method publishes only the private-Y set that passed each filter.
    ///
    /// That set holds the de-novo calls across chrY after four steps. The code removes the backbone
    /// variants, applies the callable mask, removes the structural regions, and applies the strict
    /// novel-marker [`PublishGate`]. Each published variant carries the mark of a single unverified
    /// candidate.
    ///
    /// The method never publishes the full de-novo set for chrY. The Y chromosome of CHM13 belongs
    /// to haplogroup J. Take a sample in haplogroup R. The difference between J and R, and the
    /// reads that map to the wrong paralog, give many SNPs that no curator can use.
    ///
    /// For another contig, such as chrM, the method publishes the raw de-novo calls. That set is
    /// small, it behaves well, and it is relative to rCRS. It needs no filter against a tree.
    pub(crate) async fn variants_record(&self, alignment_id: i64, contig: &str) -> Result<serde_json::Value, AppError> {
        let variants = if navigator_analysis::contig::is_chr_y(contig) {
            let bucket = self.private_y_variants_self_masked(alignment_id).await?;
            // A quality gate. A count of new variants that is too high shows a problem with the
            // sample. The causes are contamination, low coverage, and a wrong reference build. A
            // GRCh38 alignment is one example: its chrY reference holds more noise, and the
            // hs1-native tree can not resolve each of its shared-lineage variants.
            //
            // In that case the full set is doubtful, and the app publishes nothing. A curator must
            // not receive many candidates from a sample that the app already marked. The app still
            // shows those variants on the screen, below the warning banner.
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

    /// Publish the **consensus** ancestry estimates of a subject with the `client` that the caller
    /// gives. The method writes one populationBreakdown record for each method.
    ///
    /// This function is the core, and a test can call it directly. In the app, callers use
    /// [`publish_ancestry`](Self::publish_ancestry). The method returns one ref for each record.
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
        Ok(client
            .create_record(NS_BIOSAMPLE, value, Some(&biosample_rkey(biosample_guid)))
            .await?)
    }

    /// Publish a sequence-run characterization using an explicit `client`.
    pub async fn publish_sequence_run_with(
        &self,
        client: &PdsClient,
        run: &SequenceRun,
    ) -> Result<RecordRef, AppError> {
        let value = self.sequence_run_record(client.did(), run).await?;
        Ok(client
            .create_record(NS_SEQUENCERUN, value, Some(&seqrun_rkey(run.id)))
            .await?)
    }

    /// Build an ancestral-origin record for one MDKA row. The method returns `None` when the app
    /// must not publish that row.
    ///
    /// [`AncestralOriginRecord::build`] holds the gate for each field. This method only gives the
    /// join keys.
    ///
    /// The method changes the lineage name to the form that the AppView uses, which is `Y_DNA` or
    /// `MT_DNA`. It never publishes an `Auto` lineage, because that lineage has no tree.
    pub(crate) async fn ancestral_origin_record(
        &self,
        did: &str,
        m: &Mdka,
    ) -> Result<Option<serde_json::Value>, AppError> {
        let Some(lineage) = (match m.lineage.as_str() {
            "Y" => Some("Y_DNA"),
            "Mt" => Some("MT_DNA"),
            _ => None,
        }) else {
            return Ok(None);
        };
        let external_ids = self
            .external_ids(m.biosample_guid)
            .await?
            .into_iter()
            .map(|e| OriginExternalId {
                namespace: e.source,
                value: e.external_id,
            })
            .collect();
        let coord = m.latitude.zip(m.longitude);
        let record = AncestralOriginRecord::build(
            Some(biosample_at_uri(did, m.biosample_guid)),
            external_ids,
            lineage,
            m.ancestor_name.as_deref(),
            m.origin_place.as_deref(),
            m.origin_country.as_deref(),
            m.birth_year,
            m.death_year,
            coord,
            Utc::now().to_rfc3339(),
        );
        record
            .map(|r| serde_json::to_value(&r))
            .transpose()
            .map_err(AppError::from)
    }

    /// Publish the ancestral origin of one MDKA with the `client` that the caller gives. A result
    /// of `Ok(None)` shows that a gate refused the row. That result is normal and is not an error.
    ///
    /// The method calls `putRecord` at a fixed rkey. So a user can correct an MDKA and run the
    /// method again, and the new record replaces the record of that ancestor. The repository does
    /// not collect duplicates of one man.
    pub async fn publish_ancestral_origin_with(
        &self,
        client: &PdsClient,
        m: &Mdka,
    ) -> Result<Option<RecordRef>, AppError> {
        let Some(value) = self.ancestral_origin_record(client.did(), m).await? else {
            return Ok(None);
        };
        let rkey = origin_rkey(m.biosample_guid, &m.lineage);
        Ok(Some(
            client.put_record(ANCESTRAL_ORIGIN_COLLECTION, &rkey, value).await?,
        ))
    }

    /// Publish the ancestral origins of every subject this workspace may publish for.
    ///
    /// The method puts each record in the outbox. It does not send a record directly. The outbox
    /// tries again after a failure, it continues after an offline period, and it maps a fixed rkey
    /// onto `putRecord`.
    ///
    /// So a user can correct an MDKA and run the batch again, and the new record replaces the
    /// record of that ancestor. The repository does not collect duplicates. For this reason the
    /// batch is always safe to run again, and a second run changes nothing on the AppView.
    ///
    /// With `dry_run`, the method builds each record and applies each gate, but it adds nothing to
    /// the outbox. A user can then read the counts before any genealogy leaves the machine.
    ///
    /// [`navigator_store::mdka::publishable`] holds the consent test.
    /// [`AncestralOriginRecord::build`] holds the gate for each field.
    pub async fn publish_ancestral_origins(
        &self,
        lineage: Lineage,
        dry_run: bool,
    ) -> Result<OriginPublishReport, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let rows = mdka::publishable(self.store.pool(), lineage.as_str()).await?;
        let mut report = OriginPublishReport {
            considered: rows.len(),
            ..Default::default()
        };
        for m in &rows {
            let Some(value) = self.ancestral_origin_record(&did, m).await? else {
                report.refused += 1;
                continue;
            };
            // The count of rows that pass each gate. A dry run then reports the coverage, and not
            // only a total.
            if value.get("originPlace").is_some() {
                report.with_place += 1;
            } else if value.get("originCountry").is_some() {
                report.country_only += 1;
            }
            if !dry_run {
                let rkey = origin_rkey(m.biosample_guid, &m.lineage);
                self.enqueue_publish(
                    "ancestralOrigin",
                    &m.biosample_guid.0.to_string(),
                    ANCESTRAL_ORIGIN_COLLECTION,
                    Some(&rkey),
                    value,
                )
                .await?;
            }
            report.publishable += 1;
        }
        Ok(report)
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

/// Join the two views that a [`CoverageResult`] holds for each contig, and write the result to the
/// `contigs[]` field of the shared lexicon.
///
/// The two views are the statistics in the samtools form and the counts of each callable state. The
/// key of the join is the contig name. `export::coverage_tsv` uses the same join.
///
/// A contig can appear in the statistics with no callable count. That state must not occur, and the
/// function then writes zeros.
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

    /// The shape of a chrY extract. The chrY contig holds millions of reads. Each autosome and the
    /// X chromosome hold only a few reads that map to the wrong place.
    fn y_scoped_coverage() -> CoverageResult {
        CoverageResult {
            contig_coverage_stats: vec![cstat("chrY", 3_000_000), cstat("chr1", 30), cstat("chrX", 12)],
            ..Default::default()
        }
    }

    async fn alignment_with_test_type(app: &App, test_type: &str) -> i64 {
        let b = app.add_biosample(None, "yscoped", None, None).await.unwrap();
        let run = app
            .record_sequence_run(NewSequenceRun::new(b.guid, "ILLUMINA", test_type))
            .await
            .unwrap();
        app.record_alignment(NewAlignment {
            bam_path: Some("/nonexistent.cram".into()),
            ..NewAlignment::new(run.id, "chm13v2.0", "synthetic")
        })
        .await
        .unwrap()
        .id
    }

    /// An alignment with a WGS label that holds only Y reads must not publish a coverage summary.
    /// Such a summary makes the whole-genome statistics of the AppView incorrect.
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

    /// A Y test, such as Big Y, has the *same* shape and publishes as usual. Its Y coverage is
    /// correct, and the AppView puts it in a group that is separate from WGS.
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
        app.save_analysis(aln, "coverage", COVERAGE_VERSION, &wgs)
            .await
            .unwrap();
        app.coverage_record("did:plc:test", aln)
            .await
            .expect("normal WGS coverage should publish");
    }
}
