//! `impl App` methods extracted from `lib.rs` (the `import_profiles` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

/// One copy of the parsed Y-SNP dictionary for the full process.
///
/// [`YsnpDictionary`] now selects the full catalog, which holds about 2 million rows. A parse of
/// that file at each call would read about 200 MB each time. The `y_snp_names_at` function runs at
/// each view of the Y-SNP table, so those calls are frequent. This value holds the result of one
/// parse.
///
/// The key is the path of the dictionary file together with its signature, which is the mtime and
/// the size. So the code reads a new dictionary after the user replaces the file.
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
    let key = format!(
        "{}|{}",
        dict_path.display(),
        file_signature(&dict_path).unwrap_or_default()
    );
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
        // Add the markers to a profile of the same panel, when one exists. Do not make a second
        // profile. One example is a Big Y CUSTOM panel, of 700 or 500 markers, that the user
        // imports after the FTDNA project import made the profile. The code joins the two marker
        // sets. On a conflict, the value from the new import wins.
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

    /// Import the SNP variant calls of a subject from a file.
    ///
    /// The code parses a `.vcf` file as a VCF, with the shared column parser. It parses a `.csv`
    /// file or a `.tsv` file as a `contig,position,ref,alt[,rsid][,gt]` table. A YSEQ panel export
    /// and a Sanger panel export have that shape.
    ///
    /// The code keeps only SNPs. It removes each indel and each symbolic allele. `source_type` sets
    /// the weight of the source in the concordance calculation, and a Sanger source has the highest
    /// weight.
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
        // Match `.vcf`, and also `.vcf.gz` and `.vcf.bgz` from bgzip or gzip. The `extension()`
        // function reads only the last `.gz` part, and that value would send a compressed VCF to
        // the CSV branch.
        let is_vcf = path
            .file_name()
            .and_then(|n| n.to_str())
            .map(|n| n.to_ascii_lowercase())
            .is_some_and(|n| n.ends_with(".vcf") || n.ends_with(".vcf.gz") || n.ends_with(".vcf.bgz"));

        let calls = if is_vcf {
            // The parser reads the genotype. A vendor VCF from FTDNA Big Y or YSEQ also reports
            // a reference site. So the code keeps only the ALT value that the genotype selects.
            // See parse_vcf_subject_snps. For a VCF with sites only, the code keeps each listed
            // variant. The parser also reads a `.vcf.gz` file from bgzip.
            parse_vcf_subject_snps(path)?
        } else {
            let text = std::fs::read_to_string(path)?;
            variants::parse_csv(&text).map_err(AppError::Import)?
        };
        if calls.is_empty() {
            return Err(AppError::Import("no SNP variants found in file".into()));
        }

        // Find the vendor of a VCF and mark the record. The code recognizes FTDNA Big Y, Y Elite,
        // YSEQ, and mtFull. It reads the header, the file name, and a readme file in the same
        // directory.
        //
        // The code then records the vendor label, a correct SourceType, and the reference build.
        // The Y placement and the mt placement use that build for the liftover.
        //
        // A VCF with no vendor keeps the label and the source_type of the caller. A CSV import does
        // not change.
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
            // The code records this path, so it can read the VCF again and genotype at the
            // positions of the tree. The `alignment.bam_path` field has the same role for a CRAM.
            // See `App::vset_base_calls`.
            source_path: Some(path.to_string_lossy().into_owned()),
        };
        let set = variant_set::create(self.store.pool(), &new).await?;

        // Place a vendor Y-NGS VCF at the import, so the subject gets a Y haplogroup and the user
        // does not press Refresh. Such a VCF comes from FTDNA Big Y, YSEQ, Full Genomes, or a
        // similar test, and the file *is* the set of Y-SNP calls.
        //
        // The step is optional. The cache can hold no tree, and a VCF can hold only autosomal
        // data or mt data. In each case the code writes the calls and does nothing more. A file
        // with no chrY data gives no placement.
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

    /// Import a CompleteGenomics **masterVar** whole-genome variant table. The file name is
    /// `var-*-ASM.tsv` or `var-*-ASM.tsv.bz2`, and the `cgatools` program of the old CG sequencing
    /// service wrote it.
    ///
    /// [`navigator_analysis::mastervar`] reads and decompresses the file on another thread, and it
    /// makes SNP calls. A diploid heterozygous call becomes `0/1`. A homozygous call becomes `1/1`,
    /// and a haploid call becomes `1`. The code removes each indel, each `ref` span, and each
    /// `no-call` span. It keeps only SNPs, as the VCF importer and the CSV importer do.
    ///
    /// The code stores the result as a `WgsShortRead` set on GRCh37, which is the only build of CG.
    /// The chrM contig uses rCRS. The code then places the Y haplogroup at the import, as it does
    /// for a vendor Y-NGS VCF.
    ///
    /// The mtDNA result comes from the mt consensus of many sources. The chrM data of a set that is
    /// not a chip feeds `mt_source_calls`.
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
            source_path: Some(path.to_string_lossy().into_owned()),
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

    /// Import a Big Y CSV variant report from FTDNA. The report is the Named report or the Private
    /// Variants report. A project administrator receives these files when the access level gives
    /// the browser CSV files but no BAM file, CRAM file, or VCF file.
    ///
    /// Each row is a derived-allele call on chrY in GRCh38. So the code stores the rows as a
    /// `TargetedNgs` variant set on GRCh38, which is the native build of the Y tree of FTDNA. The
    /// code then places the subject with the vendor path at the import.
    ///
    /// The Named report gives a Y haplogroup directly, because its positions match the tree and
    /// need no liftover. The code also stores the Private Variants. Those loci are new and are not
    /// on the tree, and the store keeps them as a record.
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
            source_path: Some(path.to_string_lossy().into_owned()),
        };
        let set = variant_set::create(self.store.pool(), &new).await?;
        // Place the Y haplogroup from the vendor sets, which are the sets that are not a chip.
        // The Named report holds the SNPs that define a node of the tree.
        if let Err(e) = self.assign_y_vendor_vcfs(biosample_guid).await {
            eprintln!("FTDNA CSV Y placement deferred ({e})");
        }
        Ok(set)
    }

    /// Add a variant set that the user typed. The user pastes `contig,position,ref,alt` rows. One
    /// example is a set of confirmations from Sanger or YSEQ. `source_type` sets the weight, and a
    /// Sanger source has the weight 1.0.
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
            source_path: None, // parsed from text already in hand, not a file we can re-read
        };
        Ok(variant_set::create(self.store.pool(), &new).await?)
    }

    /// The build for the BISDNA calls of a subject. The method takes the first alignment whose
    /// reference build has a dictionary key. If there is none, it returns `"hs1"`, which is the
    /// default of the project.
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

    /// Add the catalogued Y-SNP **name** to each Y variant that has a position and no name. Two
    /// tables use this map: the variant profile with many sources, and the union of the private-Y
    /// sets.
    ///
    /// The method finds the Y build key of the subject. A CHM13 build gives `hs1`, and the other
    /// builds give GRCh38 or GRCh37. The BISDNA importer uses the same rule.
    ///
    /// The method then reads the Y-SNP dictionary, which is the full catalog and stays in memory.
    /// It returns a map from a position to a canonical name, for the requested positions only.
    ///
    /// The step is optional. An absent dictionary gives an empty map and no error, and the tables
    /// then show no extra name. A lookup against the wrong build finds nothing. So the table can
    /// hold an absent name, but it never holds a wrong name.
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

    /// Make sure that a Y-SNP dictionary is on the machine. At the first use, the method downloads
    /// the full catalog, `dictionary.tsv`, which is about 208 MB.
    ///
    /// The installer does not hold that file. The file is too large, and YBrowse refreshes it about
    /// once each week.
    ///
    /// The method does nothing when the machine already holds a dictionary. That dictionary can be
    /// the chromo2 panel or the full catalog. It also does nothing when `NAVIGATOR_YSNP_DIR` points
    /// to one.
    ///
    /// The method checks the download against a small published manifest,
    /// `ysnp_manifest.json`. That file has the shape of the ancestry
    /// [`AssetManifest`](navigator_analysis::manifest::AssetManifest). So a rebuild of the catalog
    /// is a new publish and not a change to the client.
    ///
    /// The step is optional. The caller then reads the dictionary, and it reports the state clearly
    /// when the file is still absent. Publish the file with `packaging/publish-assets.sh ysnp`.
    pub async fn ensure_ysnp_dictionary(&self) -> Result<(), AppError> {
        const YSNP_ASSET_BASE: &str = "https://github.com/JamesKane/decodingus-navigator/releases/download/assets-ysnp";

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
        // Check the digest from the stream against the manifest. The code does not read the 208 MB
        // file again. A manifest with no entry for the file gives a warning only, as
        // `AssetManifest::verify` does.
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

    /// Import a chromo2 Y-SNP export from BISDNA.
    ///
    /// The Y-SNP dictionary changes each marker name into a locus on `build`. When `build` is
    /// `None`, the method uses the alignment build of the subject, and then `"hs1"`.
    ///
    /// Only a **positive**, or derived, call becomes a variant call. A negative call is not a
    /// variant. The method counts a `no_call` marker, a back-mutated marker, and a marker that the
    /// dictionary does not hold. It writes none of those three.
    ///
    /// The genotype is a quality cross-check only. The verdict in the file decides between derived
    /// and ancestral, and that verdict does not depend on the Illumina TOP strand.
    ///
    /// The method stores the result as a [`VariantSet`] with the `Chip` weight.
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

        // Also write a quality summary for the array. The chromo2 chip then appears under
        // Data Sources, in the Chip and Array Profiles list. The variant set below holds the SNP
        // calls that the code can place. An array correctly has both a quality summary with its
        // provenance and a set of calls.
        //
        // BISDNA is a haploid Y panel. Each called marker is a Y marker, and heterozygosity does
        // not apply.
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
            source_path: None, // resolved through the Y-SNP dictionary; the panel is not a VCF
        };
        let variant_set = variant_set::create(self.store.pool(), &new).await?;

        // Calculate the Y haplogroup at the import. The step is optional, and with no tree the
        // code writes the calls only. The array path in `import_chip_profile_from_csv` does the
        // same.
        //
        // Without this step, a chromo2 panel from BISDNA imports its calls and never gets a
        // placement. Such a panel has no alignment genotypes in the cache, so `rebuild-signatures`
        // can not place it later either. The Y value of the subject then stays at <none>.
        if derived_calls > 0 {
            if let Err(e) = self.assign_y_bisdna(biosample_guid, Some(&build)).await {
                eprintln!("BISDNA Y placement deferred ({e})");
            }
        }

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

    /// Import the raw-data export of a genotyping array, as a CSV file or a TSV file.
    ///
    /// The method does three things. It writes the quality summary as a [`ChipProfile`]. It writes
    /// the haploid Y rows and MT rows as a [`VariantSet`] with the `Chip` source. It then tries to
    /// place the Y haplogroup, and the mtDNA haplogroup when the file holds one.
    ///
    /// This method is the consumer-array form of the chromo2 path of BISDNA. A 23andMe file holds
    /// both Y rows and MT rows. An AncestryDNA file holds Y rows and no mtDNA rows that the app can
    /// use.
    ///
    /// The stored bases go through the same placement as BISDNA. That path is
    /// [`assign_y_bisdna`](Self::assign_y_bisdna) or [`assign_mt_chip`](Self::assign_mt_chip),
    /// followed by `assemble_assignment_robust`, and it reconciles each call to the plus strand of
    /// the tree.
    ///
    /// The placement is optional. With no network, the code stores the calls, and the user can
    /// press "Assign … (panel)" later. A failed placement does not fail the import.
    ///
    /// `provider` replaces the vendor that the code finds, when the caller gives it. `chip_version`
    /// is optional.
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

        // Read the haploid Y rows and MT rows, and store them as variant calls with the Chip
        // source. The haplogroup placement then has them, and a later placement also has them,
        // with no second read of the file.
        //
        // The observed allele goes in `reference` and in `alternate`, because the app does not know
        // the ancestral allele. The placement reads `alternate`.
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
                source_path: None, // extracted from a chip export, not a re-readable call file
            };
            variant_set::create(self.store.pool(), &set).await?;

            // Compute the haplogroups on import (best-effort; an offline tree just leaves the calls).
            if y_count > 0 {
                if let Err(e) = self.assign_y_bisdna(biosample_guid, Some(&build)).await {
                    eprintln!("chip Y placement deferred ({e})");
                }
            }
            // The few MT rows of an AncestryDNA file are not an mtDNA panel that the app can use.
            // Place mtDNA only when the array holds a true MT marker set. A 23andMe file holds
            // some thousands of such markers, and the limit below removes the noise.
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

        // Derive the variants against rCRS and write them to the store. An mtDNA FASTA then gives
        // a variant set at the import. Before this step, the set appeared only in the "show
        // mutations" view. A chip import and a VCF import behave in the same way.
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
                    // These calls come from a comparison with rCRS and not from a source VCF. So
                    // there is no evidence to store.
                    evidence: Default::default(),
                })
                .collect();
            let set = NewVariantSet {
                biosample_guid,
                source_label: format!("{label} ({} variants vs rCRS)", derived.len()),
                // A full-mtDNA consensus is authoritative for its calls (gold-standard weight).
                source_type: variants::SourceType::Sanger,
                reference_build: None, // calls are rCRS-relative (contig "rCRS"), not a nuclear build
                calls,
                source_path: None, // derived by diffing the consensus against rCRS
            };
            // Best-effort: a variant-set hiccup must not lose the stored sequence.
            let _ = variant_set::create(self.store.pool(), &set).await;
        }

        // This method does NOT place the haplogroup, by design. A placement needs the mt
        // haplotree, and a read of that tree needs the network. An import must stay deterministic,
        // so it must not depend on the network. The alignment import follows the same rule. The
        // user presses "Assign mtDNA haplogroup" on the mtDNA tab to place the subject.
        Ok(seq)
    }

    /// All mtDNA sequences for a subject.
    pub async fn list_mtdna_sequences(&self, biosample_guid: SampleGuid) -> Result<Vec<MtdnaSequence>, AppError> {
        Ok(mtdna_store::list_for_biosample(self.store.pool(), biosample_guid).await?)
    }

    /// Derive the mtDNA variants of a stored sequence. The method compares that sequence with an
    /// rCRS reference FASTA and writes the result as a variant set on the contig `rCRS`. The
    /// variants then appear with the other variants of the subject. The method checks that the
    /// reference file is an mtDNA FASTA.
    ///
    /// The mutation list holds the variants against the **bundled** rCRS sequence, NC_012920.1. A
    /// banded alignment gives them. The list holds substitutions, insertions, and deletions, in the
    /// standard mtDNA notation.
    ///
    /// The method runs at the request of the user, and it does one alignment of about 16.5 kb. It
    /// stores nothing. This list is the classic mtDNA result.
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
