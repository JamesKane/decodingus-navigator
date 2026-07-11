//! `impl App` methods extracted from `lib.rs` (the `import_profiles` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

/// Process-wide memo of the parsed Y-SNP dictionary. Now that [`YsnpDictionary`] prefers the full
/// ~2M-row catalog, parsing it per resolve/annotate call (`y_snp_names_at` runs on every Y-SNP-table
/// view) would re-read ~200 MB each time; this parses once and reuses it. Keyed by the resolved
/// dictionary file's path + signature (mtime:size), so a refreshed dictionary is picked up.
type YsnpMemo = Mutex<Option<(String, Arc<YsnpDictionary>)>>;
static YSNP_MEMO: std::sync::OnceLock<YsnpMemo> = std::sync::OnceLock::new();

/// Load the Y-SNP dictionary from its asset dir, memoized process-wide (see [`YSNP_MEMO`]).
fn load_ysnp_dictionary_cached() -> Result<Arc<YsnpDictionary>, String> {
    let dir = ysnp_dict::asset_dir();
    let dict_path = YsnpDictionary::ASSET_FILENAMES
        .iter()
        .map(|f| dir.join(f))
        .find(|p| p.is_file())
        .ok_or_else(|| format!("no Y-SNP dictionary in {}", dir.display()))?;
    let key = format!("{}|{}", dict_path.display(), file_signature(&dict_path).unwrap_or_default());
    let memo = YSNP_MEMO.get_or_init(|| Mutex::new(None));
    if let Some((k, d)) = memo.lock().unwrap().as_ref() {
        if *k == key {
            return Ok(d.clone());
        }
    }
    let dict = Arc::new(YsnpDictionary::load(&dir)?);
    *memo.lock().unwrap() = Some((key, dict.clone()));
    Ok(dict)
}

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
        // Match `.vcf`, plus bgzipped/gzipped `.vcf.gz` / `.vcf.bgz` (extension() alone sees only
        // the trailing `.gz`, which would mis-route a compressed VCF to the CSV branch).
        let is_vcf = path
            .file_name()
            .and_then(|n| n.to_str())
            .map(|n| n.to_ascii_lowercase())
            .is_some_and(|n| n.ends_with(".vcf") || n.ends_with(".vcf.gz") || n.ends_with(".vcf.bgz"));

        let calls = if is_vcf {
            // Genotype-aware: a vendor VCF (FTDNA Big Y / YSEQ) reports reference sites too, so only
            // the genotype-selected ALT is kept (see parse_vcf_subject_snps). Sites-only VCFs keep
            // every listed variant. Handles a bgzipped `.vcf.gz` transparently.
            parse_vcf_subject_snps(path)?
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
        let set = variant_set::create(self.store.pool(), &new).await?;

        // Place a vendor Y-NGS VCF (FTDNA Big Y / YSEQ / Full Genomes / …) on import so it lands a
        // Y haplogroup without a manual Refresh — the VCF *is* the called Y-SNP set. Best-effort: an
        // offline tree or an autosomal/mt-only VCF just leaves the calls (no chrY → no-op).
        let has_chr_y = set
            .calls
            .iter()
            .any(|c| c.contig.eq_ignore_ascii_case("chrY") || c.contig.eq_ignore_ascii_case("y"));
        if has_chr_y {
            if let Err(e) = self.assign_y_vendor_vcfs(biosample_guid).await {
                eprintln!("vendor Y-VCF placement deferred ({e})");
            }
        }
        Ok(set)
    }

    /// Import a CompleteGenomics **masterVar** whole-genome variant table (`var-*-ASM.tsv[.bz2]`,
    /// the old CG sequencing service's `cgatools` output). The file is streamed and decompressed
    /// off-thread ([`navigator_analysis::mastervar`]) into SNP calls — each diploid het becomes a
    /// `0/1`, a homozygous/haploid call a `1/1` / `1`, indels and `ref`/`no-call` spans dropped
    /// (SNP-only, matching the VCF/CSV importer). Stored as a `WgsShortRead` set on GRCh37 (CG's
    /// only build; chrM = rCRS), then Y-placed on import like a vendor Y-NGS VCF. mtDNA falls out
    /// via the multi-source mt consensus (a non-chip set's chrM feeds `mt_source_calls`).
    pub async fn import_mastervar_from_file(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
    ) -> Result<VariantSet, AppError> {
        let label = path
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "masterVar".into());
        let file = path.to_path_buf();
        let parsed = tokio::task::spawn_blocking(move || navigator_analysis::mastervar::parse_file(&file))
            .await?
            .map_err(|e| AppError::Import(format!("reading masterVar {label}: {e}")))?;
        if parsed.calls.is_empty() {
            return Err(AppError::Import(format!(
                "no SNP variants found in masterVar {label} ({} loci scanned)",
                parsed.loci_seen
            )));
        }

        let new = NewVariantSet {
            biosample_guid,
            source_label: format!("CompleteGenomics masterVar ({label})"),
            source_type: SourceType::WgsShortRead,
            reference_build: Some(parsed.reference_build),
            calls: parsed.calls,
        };
        let set = variant_set::create(self.store.pool(), &new).await?;

        // Place Y from the derived chrY calls on import (the vendor-VCF path), so the haplogroup
        // lands without a manual Refresh. Best-effort: an offline tree just leaves the calls.
        let has_chr_y = set
            .calls
            .iter()
            .any(|c| c.contig.eq_ignore_ascii_case("chrY") || c.contig.eq_ignore_ascii_case("y"));
        if has_chr_y {
            if let Err(e) = self.assign_y_vendor_vcfs(biosample_guid).await {
                eprintln!("masterVar Y placement deferred ({e})");
            }
        }
        Ok(set)
    }

    /// Import an FTDNA Big Y CSV variant report (Named or Private Variants) — the data a project
    /// admin gets when their access tier exposes the browser CSVs but not the BAM/CRAM/VCF. The
    /// rows are GRCh38 chrY derived-allele calls, so they're stored as a `TargetedNgs` variant set
    /// on GRCh38 (FTDNA's native Y-tree build) and placed via the vendor path on import — the Named
    /// report lands a Y haplogroup directly (positions match the tree, no liftover). Private
    /// Variants are stored too (novel loci, off-tree) for the record.
    pub async fn import_ftdna_csv_variants(
        &self,
        biosample_guid: SampleGuid,
        path: &Path,
    ) -> Result<VariantSet, AppError> {
        let text = std::fs::read_to_string(path)?;
        let (report, calls) = navigator_domain::ftdna_csv::parse(&text).map_err(AppError::Import)?;
        let new = NewVariantSet {
            biosample_guid,
            source_label: report.label().to_string(),
            source_type: SourceType::TargetedNgs,
            reference_build: Some("GRCh38".to_string()),
            calls,
        };
        let set = variant_set::create(self.store.pool(), &new).await?;
        // Place Y from the vendor (non-Chip) sets — the Named report carries the tree-defining SNPs.
        if let Err(e) = self.assign_y_vendor_vcfs(biosample_guid).await {
            eprintln!("FTDNA CSV Y placement deferred ({e})");
        }
        Ok(set)
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
    /// dictionary (the full catalog, memoized), and returns `position → canonical name` for the
    /// requested positions only. Best-effort: a missing dictionary yields an empty map (not an error),
    /// so the tables simply show no extra names. Looking a position up against the wrong build just
    /// misses — there are no false labels, only possibly-absent ones.
    pub async fn y_snp_names_at(
        &self,
        biosample_guid: SampleGuid,
        positions: &[i64],
    ) -> Result<HashMap<i64, String>, AppError> {
        if positions.is_empty() {
            return Ok(HashMap::new());
        }
        let build = self.bisdna_target_build(biosample_guid).await;
        let Ok(dict) = load_ysnp_dictionary_cached() else {
            return Ok(HashMap::new()); // no dictionary installed — degrade gracefully
        };
        let idx = dict.position_index(&build);
        let names = positions
            .iter()
            .filter_map(|p| idx.get(p).map(|n| (*p, n.to_string())))
            .collect();
        Ok(names)
    }

    /// Ensure a Y-SNP dictionary is present, downloading the full catalog (`dictionary.tsv`,
    /// ~208 MB) from the asset release on first use — it's too big and too volatile (~weekly YBrowse
    /// refresh) to bundle in the installer. No-op when a dictionary (the chromo2 panel or the full
    /// catalog) is already installed, or the user pointed `NAVIGATOR_YSNP_DIR` at one. The download
    /// is verified against a small published manifest (`ysnp_manifest.json`, the ancestry
    /// [`AssetManifest`](navigator_analysis::manifest::AssetManifest) shape) so a rebuild is a
    /// re-publish, not a client change. Best-effort — the caller then loads, degrading clearly if the
    /// dictionary is still absent. Publish with `packaging/publish-assets.sh ysnp`.
    pub async fn ensure_ysnp_dictionary(&self) -> Result<(), AppError> {
        const YSNP_ASSET_BASE: &str =
            "https://github.com/JamesKane/decodingus-navigator/releases/download/assets-ysnp";

        let dir = ysnp_dict::asset_dir();
        if YsnpDictionary::ASSET_FILENAMES.iter().any(|f| dir.join(f).is_file()) {
            return Ok(());
        }
        std::fs::create_dir_all(&dir)?;

        let manifest_json = self
            .auth
            .http
            .get(format!("{YSNP_ASSET_BASE}/ysnp_manifest.json"))
            .send()
            .await
            .and_then(|r| r.error_for_status())
            .map_err(|e| AppError::Import(format!("fetching Y-SNP dictionary manifest: {e}")))?
            .text()
            .await
            .map_err(|e| AppError::Import(format!("reading Y-SNP dictionary manifest: {e}")))?;
        let manifest = navigator_analysis::manifest::AssetManifest::from_json(&manifest_json)
            .map_err(|e| AppError::Import(format!("parsing Y-SNP dictionary manifest: {e}")))?;

        let dest = dir.join("dictionary.tsv");
        let mut noop = |_: u64, _: Option<u64>| {};
        let got = navigator_refgenome::download::download(
            &self.auth.http,
            &format!("{YSNP_ASSET_BASE}/dictionary.tsv"),
            &dest,
            &mut noop,
        )
        .await?;
        // Verify the streamed digest against the manifest (no 208 MB re-read). A manifest without an
        // entry passes through advisory, matching `AssetManifest::verify`.
        if let Some(entry) = manifest.assets.get("dictionary.tsv") {
            if !got.eq_ignore_ascii_case(&entry.sha256) {
                let _ = std::fs::remove_file(&dest);
                return Err(AppError::Import(format!(
                    "Y-SNP dictionary failed its integrity check (manifest {}, download {got}) — re-try",
                    entry.sha256
                )));
            }
        }
        Ok(())
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

        // Fetch the full Y-SNP dictionary on first use (best-effort); the load below then finds it.
        if let Err(e) = self.ensure_ysnp_dictionary().await {
            eprintln!("Y-SNP dictionary download failed ({e}); trying any local copy");
        }
        let dict_dir = ysnp_dict::asset_dir();
        let dict = load_ysnp_dictionary_cached().map_err(|e| {
            AppError::Import(format!(
                "{e}. The Y-SNP dictionary downloads automatically on first import, or build it with \
                 scripts/ysnp-dictionary (expected under {})",
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
