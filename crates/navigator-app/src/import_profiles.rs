//! `impl App` methods extracted from `lib.rs` (the `import_profiles` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

impl App {
    // ---- panels + IBD ------------------------------------------------------

    /// Create a genotyping panel from explicit sites.
    pub async fn import_panel(&self, name: &str, sites: &[PanelSite]) -> Result<Panel, AppError> {
        Ok(panel::create(self.store.pool(), name, sites).await?)
    }

    /// Create a panel from a (plain-text) sites VCF — biallelic SNP rows only.
    pub async fn import_panel_from_vcf(&self, name: &str, vcf_path: &Path) -> Result<Panel, AppError> {
        let variants = navigator_analysis::parity::parse_truth_vcf(vcf_path)?;
        let sites: Vec<PanelSite> = variants
            .iter()
            .filter_map(|v| {
                let alt = v.alternate.first()?;
                (v.reference.len() == 1 && alt.len() == 1).then(|| PanelSite {
                    chrom: v.chrom.clone(),
                    position: v.pos,
                    reference_allele: v.reference.clone(),
                    alternate_allele: alt.clone(),
                    name: v
                        .ids
                        .first()
                        .cloned()
                        .unwrap_or_else(|| format!("{}:{}", v.chrom, v.pos)),
                })
            })
            .collect();
        self.import_panel(name, &sites).await
    }

    pub async fn list_panels(&self) -> Result<Vec<Panel>, AppError> {
        Ok(panel::list(self.store.pool()).await?)
    }

    // ---- STR profiles ------------------------------------------------------

    /// Import a Y-STR profile for a subject from an exported marker table (CSV/TSV).
    pub async fn import_str_profile_from_csv(
        &self,
        biosample_guid: SampleGuid,
        panel_name: &str,
        provider: Option<String>,
        source: Option<String>,
        csv_path: &Path,
    ) -> Result<StrProfile, AppError> {
        let text = std::fs::read_to_string(csv_path)?;
        let markers = strprofile::parse_csv(&text).map_err(AppError::Import)?;
        // Merge into an existing same-panel profile rather than creating a duplicate — e.g. a Big Y
        // CUSTOM (700/500) panel re-imported after the FTDNA project import already made one. Union
        // the markers, the freshly-imported value winning on a conflict.
        if let Some(existing) = str_profile::find_by_panel(self.store.pool(), biosample_guid, panel_name).await? {
            let mut merged = existing.markers.clone();
            for m in markers {
                match merged.iter_mut().find(|e| e.marker == m.marker) {
                    Some(e) => e.value = m.value,
                    None => merged.push(m),
                }
            }
            str_profile::replace_markers(self.store.pool(), existing.id, &merged).await?;
            self.assign_male_for_y_evidence(biosample_guid).await?;
            return Ok(StrProfile {
                markers: merged,
                ..existing
            });
        }
        let new = NewStrProfile {
            biosample_guid,
            panel_name: panel_name.to_string(),
            provider,
            source,
            markers,
        };
        let created = str_profile::create(self.store.pool(), &new).await?;
        self.assign_male_for_y_evidence(biosample_guid).await?;
        Ok(created)
    }

    /// All STR profiles for a subject.
    pub async fn list_str_profiles(&self, biosample_guid: SampleGuid) -> Result<Vec<StrProfile>, AppError> {
        Ok(str_profile::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    // ---- SNP variants ------------------------------------------------------

    /// Import a subject's SNP variant calls from a file. `.vcf` is parsed as a VCF (reusing
    /// the shared column parser); `.csv`/`.tsv` as a `contig,position,ref,alt[,rsid][,gt]`
    /// table (a YSEQ/Sanger panel export fits this). Indels/symbolic alleles are dropped
    /// (SNP-only). `source_type` sets the concordance weight (Sanger = gold standard).
    pub async fn import_variants_from_file(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
        source_type: SourceType,
    ) -> Result<VariantSet, AppError> {
        let label = path
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "variants".into());
        let is_vcf = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.eq_ignore_ascii_case("vcf"))
            .unwrap_or(false);

        let calls = if is_vcf {
            navigator_analysis::parity::parse_truth_vcf(path)?
                .into_iter()
                .filter_map(|v| {
                    let alt = v.alternate.first()?;
                    variants::snp_call(&v.chrom, v.pos, &v.reference, alt, v.ids.first().cloned(), None)
                })
                .collect()
        } else {
            let text = std::fs::read_to_string(path)?;
            variants::parse_csv(&text).map_err(AppError::Import)?
        };
        if calls.is_empty() {
            return Err(AppError::Import("no SNP variants found in file".into()));
        }

        // Vendor-aware tagging for VCFs: recognize FTDNA Big Y / Y Elite / YSEQ / mtFull from the
        // header + filename + sibling readme, and record the vendor label, a meaningful SourceType,
        // and the reference build (feeds Y/mt placement liftover). A generic VCF keeps the caller's
        // label/source_type. CSV imports are unchanged.
        let (source_label, source_type, reference_build) = if is_vcf {
            let (meta, contigs) = peek_vcf_header(path);
            let vendor =
                navigator_domain::vendorvcf::classify(&meta, &contigs, &label, sibling_readme(path).as_deref());
            let build = detect_vcf_build(&meta);
            if vendor.is_recognized() {
                (
                    format!("{} ({})", vendor.display(), vcf_label_context(path, &label)),
                    vendor.source_type(),
                    build,
                )
            } else {
                (label, source_type, build)
            }
        } else {
            (label, source_type, None)
        };

        let new = NewVariantSet {
            biosample_guid,
            source_label,
            source_type,
            reference_build,
            calls,
        };
        Ok(variant_set::create(self.store.pool(), &new).await?)
    }

    /// Add a manually-entered variant set — paste `contig,position,ref,alt` rows (e.g.
    /// Sanger/YSEQ confirmations). `source_type` sets the weight (Sanger = 1.0).
    pub async fn add_variants(
        &self,
        biosample_guid: SampleGuid,
        source_label: &str,
        source_type: SourceType,
        text: &str,
    ) -> Result<VariantSet, AppError> {
        let calls = variants::parse_csv(text).map_err(AppError::Import)?;
        let new = NewVariantSet {
            biosample_guid,
            source_label: source_label.to_string(),
            source_type,
            reference_build: None,
            calls,
        };
        Ok(variant_set::create(self.store.pool(), &new).await?)
    }

    /// The build to emit a subject's BISDNA calls on: the first of its alignments whose
    /// reference build maps to a dictionary key, else `"hs1"` (the project default).
    pub(crate) async fn bisdna_target_build(&self, biosample_guid: SampleGuid) -> String {
        if let Ok(aligns) = alignment::list_for_biosample(self.store.pool(), biosample_guid).await {
            for a in &aligns {
                if let Some(key) = decodingus_build_key(&a.reference_build) {
                    return key.to_string();
                }
            }
        }
        "hs1".to_string()
    }

    /// Annotate position-only Y variants with the catalogued Y-SNP **name** at that site, for the two
    /// Y-SNP tables (multi-source variant profile + private-Y union). Resolves the subject's Y build
    /// key (CHM13→`hs1`, else GRCh38/GRCh37 — same rule as the BISDNA importer), loads the Y-SNP
    /// dictionary (the small chromo2 panel manifest preferred — fast, ~14k SNPs), and returns
    /// `position → canonical name` for the requested positions only. Best-effort: a missing dictionary
    /// yields an empty map (not an error), so the tables simply show no extra names. Looking a position
    /// up against the wrong build just misses — there are no false labels, only possibly-absent ones.
    pub async fn y_snp_names_at(
        &self,
        biosample_guid: SampleGuid,
        positions: &[i64],
    ) -> Result<HashMap<i64, String>, AppError> {
        if positions.is_empty() {
            return Ok(HashMap::new());
        }
        let build = self.bisdna_target_build(biosample_guid).await;
        let Ok(dict) = YsnpDictionary::load(&ysnp_dict::asset_dir()) else {
            return Ok(HashMap::new()); // no dictionary installed — degrade gracefully
        };
        let idx = dict.position_index(&build);
        let names = positions
            .iter()
            .filter_map(|p| idx.get(p).map(|n| (*p, n.to_string())))
            .collect();
        Ok(names)
    }

    /// Import a BISDNA chromo2 Y-SNP export. Each named marker is resolved to a locus via the
    /// Y-SNP dictionary on `build` (when `None`, the subject's alignment build, else `"hs1"`).
    /// Only **positive** (derived) calls become variant calls: a negative is not a variant.
    /// `no_call`, back-mutated, and dictionary-unresolved markers are tallied but not emitted.
    /// The genotype is a QC cross-check only — the file's verdict (independent of the Illumina
    /// TOP strand) decides derived/ancestral. Stored as a `Chip`-weighted [`VariantSet`].
    pub async fn import_bisdna_from_file(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
        build: Option<&str>,
    ) -> Result<BisdnaImportSummary, AppError> {
        let text = std::fs::read_to_string(path)?;
        let calls = bisdna::parse(&text).map_err(AppError::Import)?;
        let build = match build {
            Some(b) => b.to_string(),
            None => self.bisdna_target_build(biosample_guid).await,
        };

        let dict_dir = ysnp_dict::asset_dir();
        let dict = YsnpDictionary::load(&dict_dir).map_err(|e| {
            AppError::Import(format!(
                "{e}. Build the Y-SNP dictionary with scripts/ysnp-dictionary (expected under {})",
                dict_dir.display()
            ))
        })?;

        const UNRESOLVED_SAMPLE_CAP: usize = 25;
        let outcome = bisdna::resolve_calls(&calls, &dict, &build, UNRESOLVED_SAMPLE_CAP);

        let derived_calls = outcome.calls.len();
        let label = path
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "BISDNA".into());

        // Also record an array QC summary so the chromo2 chip appears under Data Sources →
        // Chip / Array Profiles (the placeable per-SNP calls live in the variant set below; a
        // genotyping array legitimately has both a QC/provenance summary and its calls). BISDNA
        // is a Y-only haploid panel: every called marker is a Y marker, heterozygosity is n/a.
        let total = calls.len() as i64;
        let called = total - outcome.no_call as i64;
        let chip = NewChipProfile {
            biosample_guid,
            provider: "BISDNA".into(),
            chip_version: Some("chromo2".into()),
            summary: chipprofile::ChipSummary {
                total_markers_possible: total,
                total_markers_called: called,
                no_call_rate: if total > 0 {
                    outcome.no_call as f64 / total as f64
                } else {
                    0.0
                },
                het_rate: None,
                y_markers_called: called,
                mt_markers_called: 0,
                autosomal_markers_called: 0,
            },
            source_file_name: Some(label.clone()),
            source_path: None, // BISDNA is a Y-only panel — no autosomal genotypes for ancestry
        };
        chip_profile::create(self.store.pool(), &chip).await?;

        let new = NewVariantSet {
            biosample_guid,
            source_label: label,
            source_type: SourceType::Chip,
            reference_build: Some(build.clone()),
            calls: outcome.calls,
        };
        let variant_set = variant_set::create(self.store.pool(), &new).await?;

        Ok(BisdnaImportSummary {
            variant_set,
            build,
            total_markers: calls.len(),
            derived_calls,
            ancestral: outcome.ancestral,
            no_call: outcome.no_call,
            back_mutated: outcome.back_mutated,
            unresolved: outcome.unresolved,
            unresolved_names: outcome.unresolved_names,
            strand_mismatches: outcome.strand_mismatches,
        })
    }

    /// All variant sets for a subject.
    pub async fn list_variant_sets(&self, biosample_guid: SampleGuid) -> Result<Vec<VariantSet>, AppError> {
        Ok(variant_set::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    // ---- chip / array profiles ---------------------------------------------

    /// Import a genotyping-array raw-data export (CSV/TSV) and store its QC summary.
    /// `provider` overrides vendor detection when given; `chip_version` is optional.
    /// Import a genotyping-array raw-data export and (1) store its QC summary as a [`ChipProfile`],
    /// (2) store the haploid Y/MT genotype rows as a `Chip`-source [`VariantSet`], and (3)
    /// best-effort place the Y (and, where present, mtDNA) haplogroup on import — the consumer-array
    /// counterpart to BISDNA's chromo2 path. 23andMe carries both Y and MT rows; AncestryDNA carries
    /// Y but no usable mtDNA. The stored observed bases flow through the same
    /// [`assign_y_bisdna`](Self::assign_y_bisdna) / [`assign_mt_chip`](Self::assign_mt_chip) +
    /// `assemble_assignment_robust` placement as BISDNA, with plus-strand reconciliation to the tree.
    /// Placement is best-effort: an unreachable tree (offline) leaves the calls stored for a later
    /// manual "Assign … (panel)" — it does not fail the import.
    pub async fn import_chip_profile_from_csv(
        &self,
        biosample_guid: SampleGuid,
        provider: Option<String>,
        chip_version: Option<String>,
        path: &Path,
    ) -> Result<ChipProfile, AppError> {
        let text = std::fs::read_to_string(path)?;
        let (summary, detected) = chipprofile::summarize(&text).map_err(AppError::Import)?;
        let provider = provider.or(detected).unwrap_or_else(|| "OTHER".into());
        let source_file_name = path.file_name().map(|s| s.to_string_lossy().into_owned());
        let label = source_file_name.clone().unwrap_or_else(|| provider.clone());
        // Record the absolute path so ancestry-from-chip can re-read the autosomal genotypes later
        // (like alignments re-read bam_path). Canonicalize best-effort; fall back to the given path.
        let source_path = Some(
            std::fs::canonicalize(path)
                .unwrap_or_else(|_| path.to_path_buf())
                .to_string_lossy()
                .into_owned(),
        );
        let new = NewChipProfile {
            biosample_guid,
            provider: provider.clone(),
            chip_version,
            summary,
            source_file_name,
            source_path,
        };
        let profile = chip_profile::create(self.store.pool(), &new).await?;

        // Pull the haploid Y/MT genotype rows and store them as Chip-source variant calls so the
        // haplogroup placement (and later re-placement) has them without re-reading the file. The
        // observed allele goes in both `reference` and `alternate` (we don't know the ancestral);
        // the placement reads `alternate`.
        let haplo = chipprofile::haplo_calls(&text);
        if !haplo.is_empty() {
            let build = chipprofile::detect_build(&text);
            let (mut y_count, mut mt_count) = (0usize, 0usize);
            let mut variant_calls = Vec::with_capacity(haplo.len());
            for c in &haplo {
                let (contig, is_y) = match c.dna {
                    chipprofile::ChipDna::Y => ("chrY", true),
                    chipprofile::ChipDna::Mt => ("chrM", false),
                };
                let b = c.base.to_string();
                if let Some(call) =
                    variants::snp_call(contig, c.position, &b, &b, Some(c.rsid.clone()), Some("1".into()))
                {
                    if is_y {
                        y_count += 1;
                    } else {
                        mt_count += 1;
                    }
                    variant_calls.push(call);
                }
            }
            let set = NewVariantSet {
                biosample_guid,
                source_label: format!("{label} Y/MT calls"),
                source_type: SourceType::Chip,
                reference_build: Some(build.clone()),
                calls: variant_calls,
            };
            variant_set::create(self.store.pool(), &set).await?;

            // Compute the haplogroups on import (best-effort; an offline tree just leaves the calls).
            if y_count > 0 {
                if let Err(e) = self.assign_y_bisdna(biosample_guid, Some(&build)).await {
                    eprintln!("chip Y placement deferred ({e})");
                }
            }
            // AncestryDNA's stray MT rows aren't a usable mtDNA panel — only place mtDNA when the
            // array carries a real MT marker set (23andMe has thousands; the threshold filters noise).
            const MIN_MT_CALLS: usize = 20;
            if mt_count >= MIN_MT_CALLS {
                if let Err(e) = self.assign_mt_chip(biosample_guid).await {
                    eprintln!("chip mtDNA placement deferred ({e})");
                }
            }
        }

        Ok(profile)
    }

    /// All chip profiles for a subject.
    pub async fn list_chip_profiles(&self, biosample_guid: SampleGuid) -> Result<Vec<ChipProfile>, AppError> {
        Ok(chip_profile::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    // ---- mtDNA sequences ---------------------------------------------------

    /// Import a vendor mtDNA FASTA (~16,569 bp) for a subject. Validates the header,
    /// length, and bases; stores the sequence + N count.
    pub async fn import_mtdna_from_fasta(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
    ) -> Result<MtdnaSequence, AppError> {
        let text = std::fs::read_to_string(path)?;
        let parsed = mtdna::parse_fasta(&text).map_err(AppError::Import)?;
        let source_file_name = path.file_name().map(|s| s.to_string_lossy().into_owned());
        let new = NewMtdnaSequence {
            biosample_guid,
            defline: parsed.defline,
            sequence: parsed.sequence,
            n_count: parsed.n_count,
            source_file_name,
        };
        let seq = mtdna_store::create(self.store.pool(), &new).await?;

        // Derive rCRS-relative variants and persist them, so an mtDNA FASTA yields a variant set on
        // import (not only on the on-demand "show mutations" view) — like a chip/VCF import does.
        let derived = navigator_analysis::mtvariants::derive(navigator_analysis::mtvariants::rcrs(), &seq.sequence);
        if !derived.is_empty() {
            let label = mt_vendor_label(seq.source_file_name.as_deref(), seq.defline.as_deref());
            let calls = derived
                .iter()
                .map(|v| variants::VariantCall {
                    contig: "rCRS".to_string(),
                    position: v.position,
                    reference: v.reference.to_string(),
                    alternate: v.alternate.to_string(),
                    rs_id: None,
                    genotype: None,
                })
                .collect();
            let set = NewVariantSet {
                biosample_guid,
                source_label: format!("{label} ({} variants vs rCRS)", derived.len()),
                // A full-mtDNA consensus is authoritative for its calls (gold-standard weight).
                source_type: variants::SourceType::Sanger,
                reference_build: None, // calls are rCRS-relative (contig "rCRS"), not a nuclear build
                calls,
            };
            // Best-effort: a variant-set hiccup must not lose the stored sequence.
            let _ = variant_set::create(self.store.pool(), &set).await;
        }

        // Haplogroup placement is intentionally NOT run here: it needs the mt haplotree (network),
        // and coupling a deterministic import to a network fetch is what the alignment import
        // deliberately avoids too. The mtDNA tab's "Assign mtDNA haplogroup" places it on demand.
        Ok(seq)
    }

    /// All mtDNA sequences for a subject.
    pub async fn list_mtdna_sequences(&self, biosample_guid: SampleGuid) -> Result<Vec<MtdnaSequence>, AppError> {
        Ok(mtdna_store::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    /// Derive mtDNA variants for a stored sequence by comparing it to an rCRS reference
    /// FASTA, and save them as a variant set (contig `rCRS`) so they appear alongside the
    /// subject's other variants. The reference is validated as an mtDNA FASTA.
    /// The mtDNA mutation list for a stored sequence: variants relative to the **bundled** rCRS
    /// (NC_012920.1), via banded alignment — substitutions, insertions, and deletions in standard
    /// mtDNA notation. On-demand (one ~16.5 kb alignment), not stored. The classic mtDNA result.
    pub async fn mtdna_variants(&self, mtdna_id: i64) -> Result<Vec<MtVariant>, AppError> {
        let seq = mtdna_store::get(self.store.pool(), mtdna_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("mtDNA sequence {mtdna_id}"))))?;
        Ok(navigator_analysis::mtvariants::derive(
            navigator_analysis::mtvariants::rcrs(),
            &seq.sequence,
        ))
    }
}
