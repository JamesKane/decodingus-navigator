//! `impl App` methods extracted from `lib.rs` (the `gvcf` cluster). Split out in the
//! 2026-06 simplification round; `use super::*` reaches the crate-root types + free helpers.
use super::*;

/// Read a chrY position mask, or blocklist, from a BED file in the application bundle, for one
/// build. The step is optional.
///
/// The variable `env_var` gives the path when the user sets it. If not, the code reads
/// `<cache base>/masks/<stem>.<build_token>.bed`. It tries the `.bed.gz` file first, because the
/// bundle holds that form, and then a plain `.bed` file.
///
/// The function returns `None` when the file is absent, when the parser refuses it, and when it
/// holds no row. An absent cohort asset then removes that filter, and it stops no work.
fn load_y_position_bed(env_var: &str, stem: &str, build_token: &str) -> Option<navigator_analysis::mask::RegionMask> {
    let candidates: Vec<PathBuf> = if let Ok(p) = std::env::var(env_var) {
        vec![PathBuf::from(p)]
    } else {
        let dir = refgenome_cache::base_dir().join("masks");
        let file = format!("{stem}.{build_token}.bed");
        vec![dir.join(format!("{file}.gz")), dir.join(file)]
    };
    candidates.into_iter().find_map(|path| {
        navigator_analysis::mask::RegionMask::from_bed(&path, "chrY")
            .ok()
            .filter(|m| !m.is_empty())
    })
}

/// The output directories of a caller, beside the alignment, that usually hold the GVCF file of one
/// sample.
///
/// A pipeline with more than one caller puts the output of each caller in its own directory, and not
/// beside the CRAM file. A search of the directory of the alignment alone finds no such file.
const CALLER_SUBDIRS: [&str; 3] = ["gatk4", "gatk3", "gvcf"];

/// Find the GVCF file of one sample for `contig`, beside an alignment.
///
/// The function reads the directory of the alignment first, where the ytree layout puts a
/// `*.chrY.g.vcf.gz` file. It then reads the known directories of each caller, where the usual name
/// is a plain `chrY.g.vcf.gz`.
///
/// Both names matter. The flat ytree layout writes `<sample>.chrY.g.vcf.gz` beside the CRAM file. A
/// pipeline that works on one run writes `gatk4/chrY.g.vcf.gz`, and that name holds no sample
/// prefix.
///
/// An earlier version matched the dotted name only, and it found no file of the second kind. The
/// GVCF file is what lets the placement skip the CRAM decode. So a search that fails changes a read
/// of some seconds into a walk of the full chromosome, which needs some minutes.
fn gvcf_beside_alignment(aln: &Alignment, contig_token: &str) -> Option<PathBuf> {
    let dotted = format!(".{contig_token}.g.vcf.gz");
    let bare = format!("{contig_token}.g.vcf.gz");
    let matches = |name: &str| {
        let n = name.to_ascii_lowercase();
        n.ends_with(&dotted) || n == bare
    };
    let dir = Path::new(aln.bam_path.as_ref()?).parent()?;
    let scan = |d: &Path| -> Option<PathBuf> {
        std::fs::read_dir(d)
            .ok()?
            .flatten()
            .find_map(|e| matches(&e.file_name().to_string_lossy()).then(|| e.path()))
    };
    scan(dir).or_else(|| CALLER_SUBDIRS.iter().find_map(|sub| scan(&dir.join(sub))))
}

/// Find the chrY GVCF file of one sample, for an alignment. The `NAVIGATOR_Y_GVCF` variable gives
/// the path when the user sets it. If not, the code calls [`gvcf_beside_alignment`].
///
/// The function returns `None` when it finds no file. The private-Y path then uses the pileup
/// caller, and the placement walks the full CRAM file.
pub(crate) fn chr_y_gvcf_for_alignment(aln: &Alignment) -> Option<PathBuf> {
    if let Ok(p) = std::env::var("NAVIGATOR_Y_GVCF") {
        let p = PathBuf::from(p);
        if p.is_file() {
            return Some(p);
        }
    }
    gvcf_beside_alignment(aln, "chry")
}

/// Find the chrM GVCF file of one sample, for an alignment. The `NAVIGATOR_M_GVCF` variable gives
/// the path when the user sets it. If not, the code calls [`gvcf_beside_alignment`]. This function
/// is the mtDNA form of [`chr_y_gvcf_for_alignment`].
pub(crate) fn chr_m_gvcf_for_alignment(aln: &Alignment) -> Option<PathBuf> {
    if let Ok(p) = std::env::var("NAVIGATOR_M_GVCF") {
        let p = PathBuf::from(p);
        if p.is_file() {
            return Some(p);
        }
    }
    gvcf_beside_alignment(aln, "chrm")
}

/// The token in the mask file name for the reference build of an alignment. The function returns
/// `None` when the bundle holds no chrY mask for that build.
///
/// The CHM13 masks are native, on hs1. CrossMap moves them from hs1 to hg38, and those files are the
/// GRCh38 masks. There is no mask for GRCh37 yet. That build names its contig `Y`, and no code moved
/// the masks to it.
fn y_mask_build_token(build: &str) -> Option<&'static str> {
    match canonical_build(build) {
        Some(ReferenceBuild::Chm13v2 | ReferenceBuild::Chm13v2MaskedRcrs) => Some("chm13v2"),
        Some(ReferenceBuild::Grch38) => Some("grch38"),
        _ => None,
    }
}

/// A shared, process-wide copy of one build's chrY structural regions (see the memo in
/// [`App::y_structural_regions_for`]).
type YRegionsHandle = std::sync::Arc<navigator_analysis::mask::YStructuralRegions>;

impl App {
    // ---- fast path: place haplogroups from precomputed pipeline GVCFs ---------

    /// Build the base call of each tree position for an alignment, from a **GVCF file that another
    /// tool made**. This path is the fast path, and it does no pileup on the CRAM file.
    ///
    /// The method moves each tree position onto the build of the GVCF file when the two coordinate
    /// spaces differ. One example is an mt tree in rCRS against a CHM13 `chrM` contig. The CRAM path
    /// does the same step. The method then reads the GVCF file, and it reads no alignment record.
    pub(crate) async fn gvcf_base_calls(
        &self,
        alignment_id: i64,
        contig: &str,
        gvcf: &Path,
        tree: &navigator_analysis::haplo::HaploTree,
        tree_source_build: Option<&str>,
    ) -> Result<HashMap<i64, char>, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        // The method needs the reference. A hom-ref site in a GVCF file states that the base of
        // the sample equals the base of the reference.
        //
        // The reference is itself deep in the tree. The CHM13 reference holds the Y chromosome of
        // HG002, which is in haplogroup J1. So its base at a tree position is often the *derived*
        // allele and not the ancestral one.
        //
        // So the code reads the reference base at each callable tree position, which is the value
        // that `call_bases_at` also reads.
        let reference = match aln.reference_path {
            Some(p) => PathBuf::from(p),
            None => {
                self.gateway
                    .resolve_reference(&aln.reference_build, &mut |_, _| {})
                    .await?
            }
        };
        let targets: HashSet<i64> = tree
            .nodes
            .values()
            .flat_map(|n| n.loci.iter().map(|l| l.position))
            .collect();
        if targets.is_empty() {
            return Ok(HashMap::new());
        }
        let params = gvcf::GvcfReadParams::default();

        let lifted = self
            .lifted_targets(
                &aln.reference_build,
                Some(&reference),
                contig,
                &targets,
                tree_source_build,
            )
            .await?;

        match lifted {
            // Native: tree positions are already in the GVCF's coordinates → direct read, then
            // resolve hom-ref bases from the reference at the same positions.
            None => {
                let gvcf = gvcf.to_path_buf();
                let contig_s = contig.to_string();
                let targets2 = targets.clone();
                let called =
                    tokio::task::spawn_blocking(move || gvcf::read_called_bases(&gvcf, &contig_s, &targets2, &params))
                        .await??;
                let ref_base = self.reference_bases(&reference, contig, &called.callable).await?;
                Ok(gvcf::assemble_calls(&called, &ref_base))
            }
            // The code moved the positions. It reads the GVCF file at each new contig, and it
            // reads the reference bases there. It then maps each observation back to a tree
            // position. For a position on the minus strand, it takes the reverse complement.
            Some(lifted) => {
                let mut by_contig: HashMap<String, HashSet<i64>> = HashMap::new();
                for lp in &lifted {
                    by_contig.entry(lp.contig.clone()).or_default().insert(lp.pos);
                }
                let mut all = gvcf::CalledBases::default();
                let mut ref_base: HashMap<i64, char> = HashMap::new();
                for (qcontig, set) in by_contig {
                    let gvcf = gvcf.to_path_buf();
                    let qc = qcontig.clone();
                    let set2 = set.clone();
                    let called =
                        tokio::task::spawn_blocking(move || gvcf::read_called_bases(&gvcf, &qc, &set2, &params))
                            .await??;
                    ref_base.extend(self.reference_bases(&reference, &qcontig, &called.callable).await?);
                    all.variant_bases.extend(called.variant_bases);
                    all.callable.extend(called.callable);
                }
                Ok(assemble_calls_lifted(&all, &lifted, &ref_base))
            }
        }
    }

    /// The reference genome bases at each of the `positions` on `contig`. Each base is an upper-case
    /// A, C, G, or T.
    ///
    /// The method reads the contig sequence one time, on another thread. Each position is 1-based.
    /// The result holds no position with another base, and no position outside the contig.
    ///
    /// The GVCF fast path calls this method. It needs the real base at a hom-ref tree site.
    async fn reference_bases(
        &self,
        reference: &Path,
        contig: &str,
        positions: &HashSet<i64>,
    ) -> Result<HashMap<i64, char>, AppError> {
        if positions.is_empty() {
            return Ok(HashMap::new());
        }
        let reference = reference.to_path_buf();
        let contig = contig.to_string();
        let positions: Vec<i64> = positions.iter().copied().collect();
        let map = tokio::task::spawn_blocking(
            move || -> Result<HashMap<i64, char>, navigator_analysis::AnalysisError> {
                let seq = navigator_analysis::reader::read_contig_sequence(&reference, &contig)?;
                let mut m = HashMap::with_capacity(positions.len());
                for p in positions {
                    if p >= 1 && (p as usize) <= seq.len() {
                        let b = seq[p as usize - 1].to_ascii_uppercase();
                        if matches!(b, b'A' | b'C' | b'G' | b'T') {
                            m.insert(p, b as char);
                        }
                    }
                }
                Ok(m)
            },
        )
        .await??;
        Ok(map)
    }

    /// The fingerprint of a placement that came from a GVCF file. The value joins the content hash
    /// of that file with the hash of the tree.
    ///
    /// The value is different from the fingerprint of [`Self::y_score_fingerprint`], which the CRAM
    /// path writes. This one starts with `gv:`, and that one starts with `f:`. So a later deep
    /// analysis can see that the call came from a sidecar file, and it can then skip a step.
    async fn gvcf_fingerprint(&self, gvcf: &Path, tree_json: &str, tag: &str) -> Result<String, AppError> {
        let h = sha256_file_async(gvcf.to_path_buf()).await?;
        Ok(format!("gv:{}|{}:{}", &h[..16], tag, &sha256_str(tree_json)[..16]))
    }

    /// Assign a Y haplogroup from a chrY GVCF file that another tool made. The method walks no CRAM
    /// file.
    ///
    /// It places the sample against the DecodingUs tree, at the native build of the alignment, and
    /// it needs no liftover. It writes the call under the same source key as the CRAM path, which is
    /// `aln:{id}`, with a fingerprint that starts with `gv:`.
    ///
    /// The method fails when the DecodingUs tree has no coordinates for that build, and when it can
    /// not read the tree. The caller is `ingest_sidecars`, and it then leaves the Y haplogroup for
    /// the deep pass.
    pub async fn assign_y_from_gvcf(&self, alignment_id: i64, gvcf: &Path) -> Result<HaploAssignment, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let build_key = decodingus_build_key(&aln.reference_build).ok_or_else(|| {
            AppError::Import(format!(
                "no DecodingUs tree coordinates for build {}",
                aln.reference_build
            ))
        })?;
        let tree_json = self.fetch_decodingus_y_tree().await?;
        let tree = navigator_analysis::haplo::parse_decodingus_json(&tree_json, build_key).map_err(AppError::Import)?;
        let calls = self.gvcf_base_calls(alignment_id, "chrY", gvcf, &tree, None).await?;
        // Use the proportional-top selection, and not the strict guard that the alignment path
        // needs.
        //
        // A GVCF file from a joint genotype step gives confident calls. A few of those calls
        // contradict the deep backbone with an ancestral state. The causes are a recurrent site,
        // the CHM13 reference, which is in haplogroup J1, and the hard filters of the joint step.
        //
        // The strict `path_admissible` rule then refuses the true deep lineage and takes a node
        // near the root. Sample HG00096 gave A1b, and its true terminal is R1b1a1b1a1a. The `score`
        // function ranks that terminal first, at 344 of 364.
        //
        // The data has the same shape as BISDNA chip data: confident, with a few contradictions.
        // See [`assemble_assignment_robust`].
        let assignment = assemble_assignment_robust(&tree, &calls);
        if let Ok(bio) = self.biosample_of_alignment(alignment_id).await {
            let fp = self.gvcf_fingerprint(gvcf, &tree_json, "yt").await.ok();
            self.record_call_fp(
                bio,
                DnaType::Y,
                &external_y_source_key(alignment_id),
                format!("aln #{alignment_id} Y (pipeline GVCF)"),
                &assignment,
                CallProvenance::External,
                fp.as_deref(),
            )
            .await?;
        }
        Ok(assignment)
    }

    /// Assign an mtDNA haplogroup from a chrM GVCF file that another tool made. The method walks no
    /// CRAM file.
    ///
    /// It places the sample against the mt tree of FTDNA. On GRCh38, it reads the rCRS positions of
    /// that tree directly. On CHM13, it moves those positions onto the `chrM` contig. The code makes
    /// that map itself, at a low cost.
    ///
    /// The method writes the call under the mt source key of the CRAM path, which is `aln:{id}:mt`,
    /// with a fingerprint that starts with `gv:`.
    pub async fn assign_mt_from_gvcf(&self, alignment_id: i64, gvcf: &Path) -> Result<HaploAssignment, AppError> {
        let tree_json = self.fetch_ftdna_mt_tree().await?;
        let tree = navigator_analysis::haplo::parse_ftdna_json(&tree_json).map_err(AppError::Import)?;
        let source_build = tree_build_for_contig("chrM"); // None → rCRS-direct / chrM lift
        let calls = self
            .gvcf_base_calls(alignment_id, "chrM", gvcf, &tree, source_build)
            .await?;
        // Use the proportional-top selection, as the Y path does. See assign_y_from_gvcf. The
        // confident calls of a GVCF file fit that rule better than the strict alignment guard.
        let assignment = assemble_assignment_robust(&tree, &calls);
        if let Ok(bio) = self.biosample_of_alignment(alignment_id).await {
            let fp = self.gvcf_fingerprint(gvcf, &tree_json, "mt").await.ok();
            self.record_call_fp(
                bio,
                DnaType::Mt,
                &external_mt_source_key(alignment_id),
                format!("aln #{alignment_id} mtDNA (pipeline GVCF)"),
                &assignment,
                CallProvenance::External,
                fp.as_deref(),
            )
            .await?;
        }
        Ok(assignment)
    }

    /// The sidecar paths that this alignment came from. [`Self::ingest_sidecars`] writes them.
    ///
    /// The method reads the value and does not compare the mtime of the source. This value records
    /// *what the app used*. It is not a derived result, so a change to the CRAM file does not make
    /// it wrong.
    ///
    /// The method returns `None` for an alignment that never used the fast path. Such an alignment
    /// arrived before this record existed, or it had no sidecar file.
    pub async fn recorded_sidecars(&self, alignment_id: i64) -> Result<Option<SampleSidecars>, AppError> {
        match artifact::get(self.store.pool(), alignment_id, SIDECARS_KIND, SIDECARS_VERSION).await? {
            Some(a) => Ok(serde_json::from_str(&a.payload).ok()),
            None => Ok(None),
        }
    }

    /// Read the pipeline sidecar files of a sample onto one alignment, on the fast path.
    ///
    /// The method places the Y haplogroup and the mt haplogroup from the GVCF files. It fills the
    /// sex, the read metrics, and a small coverage result from the text sidecar files. It reads no
    /// CRAM file.
    ///
    /// Each step is independent, and each one is optional. A failure goes into the report that the
    /// method returns, and the other steps continue. An absent sidecar file, or one that does not
    /// match, leaves that result for the deep pass.
    ///
    /// The method returns the values that it filled.
    pub async fn ingest_sidecars(
        &self,
        alignment_id: i64,
        sidecars: &SampleSidecars,
    ) -> Result<SidecarIngest, AppError> {
        let mut out = SidecarIngest::default();

        // Record the files that this alignment came from, before the code reads them.
        //
        // The code finds those files with one directory scan, at the import. Without this record,
        // the fast path runs one time only. A Y placement from a GVCF file, against the tree of
        // that day, could never run again. The `haplogroup_call` row then stayed after each later
        // tree.
        //
        // `App::replace_against_current_tree` reads this record and runs the placement again.
        //
        // The step is optional. A workspace that can not write the paths must still receive the
        // data.
        let _ = self
            .save_analysis_with_provenance(
                alignment_id,
                SIDECARS_KIND,
                SIDECARS_VERSION,
                sidecars,
                "sidecar",
                "full",
            )
            .await;

        if let Some(gvcf) = &sidecars.chr_y_gvcf {
            match self.assign_y_from_gvcf(alignment_id, gvcf).await {
                Ok(a) => out.y_haplogroup = a.ranked.first().map(|r| r.name.clone()),
                Err(e) => out.errors.push(format!("Y from GVCF: {e}")),
            }
        }
        if let Some(gvcf) = &sidecars.chr_m_gvcf {
            match self.assign_mt_from_gvcf(alignment_id, gvcf).await {
                Ok(a) => out.mt_haplogroup = a.ranked.first().map(|r| r.name.clone()),
                Err(e) => out.errors.push(format!("mt from GVCF: {e}")),
            }
        }
        if let Some(path) = &sidecars.sex {
            match self.ingest_sex_sidecar(alignment_id, path).await {
                Ok(Some(s)) => out.sex = Some(s),
                Ok(None) => {} // kept an existing full result (reimport) — not overwritten
                Err(e) => out.errors.push(format!("sex: {e}")),
            }
        }
        // The read metrics. The source with the most data wins. The order is samtools `stats`,
        // which holds the full data with each histogram, then Picard AlignmentSummaryMetrics, then
        // samtools `flagstat`, which holds counts only.
        match self.ingest_read_metrics(alignment_id, sidecars).await {
            Ok(true) => out.read_metrics = true,
            Ok(false) => {}
            Err(e) => out.errors.push(format!("read metrics: {e}")),
        }
        // The coverage. The samtools `coverage` output holds the statistics of each contig. The
        // Picard CollectWgsMetrics output holds the depth distribution of the full genome. That
        // distribution is the median, the sd, the MAD, the fractions that the tool excluded, and
        // the pct_Nx values.
        //
        // Use each file that exists, and write the distribution onto the table of contigs.
        if sidecars.coverage.is_some() || sidecars.wgs_metrics.is_some() {
            match self.ingest_coverage_sidecar(alignment_id, sidecars).await {
                Ok(wrote) => out.lite_coverage = wrote,
                Err(e) => out.errors.push(format!("coverage: {e}")),
            }
        }
        Ok(out)
    }

    /// Ingest inferred sex from the sidecar. `Ok(None)` when an equal-or-fuller result is already
    /// stored (reimport) so we neither overwrite the artifact nor re-stamp the sequence run.
    async fn ingest_sex_sidecar(&self, alignment_id: i64, path: &Path) -> Result<Option<String>, AppError> {
        let text = tokio::fs::read_to_string(path)
            .await
            .map_err(|e| AppError::Import(format!("{}: {e}", path.display())))?;
        let result = sidecar::parse_sex(&text);
        let wrote = self
            .save_analysis_no_downgrade(alignment_id, "sex", "1", &result, "pipeline-sidecar", "full")
            .await?;
        if !wrote {
            return Ok(None);
        }
        self.write_back_inferred_sex(alignment_id, &result).await?;
        Ok(Some(
            match result.inferred_sex {
                InferredSex::Male => "M",
                InferredSex::Female => "F",
                InferredSex::Unknown => "U",
            }
            .to_string(),
        ))
    }

    /// Ingest read metrics from the best available sidecar (priority: samtools `stats` →
    /// Picard AlignmentSummaryMetrics → samtools `flagstat`). Returns whether one was found.
    async fn ingest_read_metrics(&self, alignment_id: i64, sidecars: &SampleSidecars) -> Result<bool, AppError> {
        let read = |p: &Path| {
            let p = p.to_path_buf();
            async move {
                tokio::fs::read_to_string(&p)
                    .await
                    .map_err(|e| AppError::Import(format!("{}: {e}", p.display())))
            }
        };
        // (metrics, completeness): samtools stats is full (carries histograms); the others are
        // counts/scalars only, so `partial` lets a deep read-metrics walk upgrade them later.
        let (metrics, completeness) = if let Some(p) = &sidecars.stats {
            (sidecar::parse_samtools_stats(&read(p).await?), "full")
        } else if let Some(p) = &sidecars.alignment_summary {
            match sidecar::parse_alignment_summary(&read(p).await?) {
                Some(m) => (m, "partial"),
                None => return Ok(false),
            }
        } else if let Some(p) = &sidecars.flagstat {
            (sidecar::parse_flagstat(&read(p).await?), "partial")
        } else {
            return Ok(false);
        };
        // A second import must not replace a full deep walk with a smaller result. Keep the stored
        // result when it is the same or better.
        let wrote = self
            .save_analysis_no_downgrade(
                alignment_id,
                "read_metrics",
                "1",
                &metrics,
                "pipeline-sidecar",
                completeness,
            )
            .await?;
        Ok(wrote)
    }

    /// Read the small coverage result from the sidecar files. The method returns `true` when it
    /// wrote the result.
    ///
    /// It returns `false` when the store already holds a coverage artifact that is the same or
    /// better. A deep walk from an earlier run gives such an artifact.
    async fn ingest_coverage_sidecar(&self, alignment_id: i64, sidecars: &SampleSidecars) -> Result<bool, AppError> {
        let read = |p: &Path| {
            let p = p.to_path_buf();
            async move {
                tokio::fs::read_to_string(&p)
                    .await
                    .map_err(|e| AppError::Import(format!("{}: {e}", p.display())))
            }
        };
        // The statistics of each contig, and the callable counts, from the samtools coverage
        // output. The code starts from an empty value when that file is absent.
        let lite = match &sidecars.coverage {
            Some(cp) => {
                let cov = read(cp).await?;
                let summary = match &sidecars.callable_summary {
                    Some(p) => Some(read(p).await?),
                    None => None,
                };
                sidecar::lite_coverage(&cov, summary.as_deref())
            }
            None => CoverageResult::default(),
        };
        // Write the genome-wide depth distribution of Picard onto the table of contigs. Start from
        // the Picard result, which holds the median, the sd, the MAD, the fractions that the tool
        // excluded, and the pct_Nx values. Then add the statistics of each contig.
        let result = match &sidecars.wgs_metrics {
            Some(wp) => match sidecar::parse_wgs_metrics(&read(wp).await?) {
                Some(mut w) => {
                    w.contig_coverage_stats = lite.contig_coverage_stats;
                    w.contig_callable = lite.contig_callable;
                    w.callable_bases = lite.callable_bases;
                    if w.genome_territory == 0 {
                        w.genome_territory = lite.genome_territory;
                    }
                    if w.mean_coverage == 0.0 {
                        w.mean_coverage = lite.mean_coverage;
                    }
                    w
                }
                None => lite,
            },
            None => lite,
        };
        // The result keeps the `partial` mark. It holds no depth histogram for each base, because
        // only the deep walk makes one. So the deep pass still replaces this result.
        //
        // The store holds it under the standard coverage key. A second import must never replace a
        // full deep-walk result with this one. Keep the stored result when it exists.
        let wrote = self
            .save_analysis_no_downgrade(
                alignment_id,
                "coverage",
                coverage::COVERAGE_VERSION,
                &result,
                "pipeline-sidecar",
                "partial",
            )
            .await?;
        Ok(wrote)
    }

    /// The callable intervals of `contig`, from the reads of the alignment itself. The intervals
    /// are in the BED form, which is 0-based and half open.
    ///
    /// The parameters change with the sample. A long read from a HiFi test becomes callable at a
    /// lower depth. The limit on a CALLABLE run also grows with the length of the molecule, at
    /// `f` times the fragment length. So a long molecule passes that limit across much more of
    /// chrY.
    ///
    /// The method needs the BAM file.
    pub async fn callable_chr_intervals(&self, alignment_id: i64, contig: &str) -> Result<Vec<(i64, i64)>, AppError> {
        // Find the reference through the gateway when the alignment holds no path. No reader can
        // decode a CRAM file without a reference, and most imported alignments hold a NULL
        // `reference_path` value. The row then records the build only. The de-novo caller finds the
        // reference in the same way.
        let (bam, reference) = self.alignment_bam_reference(alignment_id).await?;
        let reference = Some(reference);
        let contig = contig.to_string();
        tokio::task::spawn_blocking(move || {
            let (read_len, frag_len) = coverage::estimate_molecule_lengths(&bam, reference.as_deref())?;
            let molecule = frag_len.max(read_len);
            let mut params = CallableLociParams::default();
            // Long, accurate reads (HiFi) are callable from a single read (see adaptive_min_depth).
            params.min_depth = adaptive_min_depth(params.min_depth, read_len);
            let min_run_len = molecule.round().max(1.0) as u32; // f = 1.0
            coverage::callable_intervals(&bam, &contig, &params, min_run_len, reference.as_deref())
        })
        .await?
        .map_err(Into::into)
    }

    /// The **private bucket**. It holds the de-novo SNP calls on chrY that the Y placement does not
    /// explain. Those calls are not on the backbone that the code assigned.
    ///
    /// The method puts each call in one of two groups. A known call off the path marks a finer FTDNA
    /// branch, or a branch beside the assigned one. A new call is a candidate for a new branch.
    ///
    /// With `callable_bed`, such as the Poznik file `b38_sites.bed` from 1KG, the method removes
    /// each call outside a reliable region.
    pub async fn private_y_variants(
        &self,
        alignment_id: i64,
        callable_bed: Option<&Path>,
    ) -> Result<PrivateBucket, AppError> {
        let mask = match callable_bed {
            Some(p) => Some(navigator_analysis::mask::RegionMask::from_bed(p, "chrY")?),
            None => None,
        };
        self.private_y_core(alignment_id, mask).await
    }

    /// The work of [`private_y_variants`], with the callable-Y BED of the sample itself as the
    /// mask. That mask changes with the depth and the read technology of the sample, and it needs no
    /// other file.
    ///
    /// With a GVCF sidecar for the sample, the method **does not** apply that mask. The confidence
    /// values in the GVCF file are the evidence that a site is callable.
    ///
    /// A second depth limit from Navigator would remove GATK calls, and the purpose of this path is
    /// to trust those calls. The method also avoids a CRAM walk, so the GVCF fast path stays fast.
    ///
    /// The reliability then comes from the callable mask of the cohort and the GQ value of the GVCF
    /// file.
    pub async fn private_y_variants_self_masked(&self, alignment_id: i64) -> Result<PrivateBucket, AppError> {
        let aln = self.alignment_or_err(alignment_id).await?;
        let mask = if chr_y_gvcf_for_alignment(&aln).is_some() {
            None
        } else {
            let intervals = self.callable_chr_intervals(alignment_id, "chrY").await?;
            Some(navigator_analysis::mask::RegionMask::from_intervals(intervals))
        };
        let bucket = self.private_y_core(alignment_id, mask).await?;
        // Write the masked bucket to the store, so the next session reads it and calculates
        // nothing.
        //
        // Version 3 takes the GVCF sidecar of the sample as its source of derived calls. Version 2
        // used the pileup only. So the code must calculate a version 2 value again, and it must not
        // read it.
        //
        // Version 4 classifies each private variant against the structural masks in the build of
        // the alignment. A version 3 bucket on a GRCh38 alignment used no mask.
        self.save_analysis(alignment_id, "private_y", "4", &bucket).await?;
        Ok(bucket)
    }

    /// Cached self-masked private-Y bucket for an alignment, if previously computed.
    pub async fn cached_private_y(&self, alignment_id: i64) -> Result<Option<PrivateBucket>, AppError> {
        self.load_analysis(alignment_id, "private_y", "4").await
    }

    /// The shared core. It runs five steps. It assigns the Y haplogroup and calls de-novo variants
    /// on chrY. It then removes the backbone calls, applies a mask when the caller asks for one, and
    /// classifies each remaining call.
    ///
    /// It also reads the curated structural regions of chrY on CHM13, which are the palindromes, the
    /// amplicons, and the AZF-DYZ regions. It finds and caches the three BED files at the first use.
    /// Each step is optional, and a failed download or a failed parse gives `None`. So this
    /// annotation never stops the analysis.
    ///
    /// It also reads the genome-region metadata of a build, which is the centromere, the telomeres,
    /// the cytobands, and the PAR regions. Those values come from the two-layer cache of the
    /// gateway, and the gateway reads the UCSC cytoBand table when the cache holds nothing. The UI
    /// uses them for quality checks and for context.
    pub async fn genome_regions(&self, build: &str) -> Result<std::sync::Arc<GenomeRegions>, AppError> {
        Ok(self.gateway.genome_regions(build, &mut |_, _| {}).await?)
    }

    /// The region of a 1-based `position` on `contig` in `build`. The value states whether the
    /// position is in a centromere, a telomere, or a PAR region, and it gives the cytoband name. The
    /// method reads the cache only, and it returns `None` when the cache holds nothing.
    pub fn region_annotation(&self, build: &str, contig: &str, position: i64) -> Option<RegionAnnotation> {
        self.gateway
            .cached_genome_regions(build)
            .map(|r| r.annotate(contig, position))
    }

    /// A memory for [`y_structural_regions_for`]. A liftover parses the full chain file. A project
    /// pass across thousands of subjects would repeat that parse for each subject, and the tree
    /// fetch had the same fault. The masks do not change inside one process, so the code resolves
    /// each build one time.
    fn y_regions_memo() -> &'static std::sync::Mutex<HashMap<String, Option<YRegionsHandle>>> {
        static MEMO: std::sync::OnceLock<std::sync::Mutex<HashMap<String, Option<YRegionsHandle>>>> =
            std::sync::OnceLock::new();
        MEMO.get_or_init(|| std::sync::Mutex::new(HashMap::new()))
    }

    /// The curated chrY structural regions **in `build`'s coordinates**.
    ///
    /// The three BED files are native to CHM13, so the code moves them to any other build.
    ///
    /// That step is important. Without it, a GRCh38 source or a GRCh37 source has no structural
    /// mask. Each call in a palindrome and in an amplicon then counts as unique sequence, and the
    /// count of private variants grows into the hundreds. One donor gave a mean of 4 with the mask
    /// and a mean of 661 without it.
    ///
    /// Each step is optional. A failed download, a failed liftover, and a failed parse each give
    /// `None`, so this annotation never stops the analysis.
    async fn y_structural_regions_for(&self, build: &str) -> Option<YRegionsHandle> {
        // The key is the *canonical* build. The builds `hs1`, `CHM13v2.0`, and the masked build
        // use the same coordinates. So they must share one entry, and the code must not do three
        // liftovers.
        let key = canonical_build(build)?.as_str().to_string();
        if let Some(hit) = Self::y_regions_memo().lock().unwrap().get(&key) {
            return hit.clone();
        }
        let built = self.build_y_structural_regions(build).await.map(std::sync::Arc::new);
        if built.is_none() {
            // The code caches this state, so a batch does not try a failed download again for
            // each subject. It also writes a message.
            //
            // The state "no structural mask" grew the count of private variants into the hundreds.
            // The app must never reach that state again with no message.
            eprintln!(
                "no chrY structural mask available for {key} — private-variant counts will include \
                 palindromic and amplicon calls"
            );
        }
        Self::y_regions_memo().lock().unwrap().insert(key, built.clone());
        built
    }

    async fn build_y_structural_regions(&self, build: &str) -> Option<navigator_analysis::mask::YStructuralRegions> {
        use navigator_analysis::mask::{RegionMask, YStructuralRegions};

        let amplicon = self
            .gateway
            .resolve_mask("chm13v2.0Y_amplicons_v1", &mut |_, _| {})
            .await
            .ok()?;
        let palindrome = self
            .gateway
            .resolve_mask("chm13v2.0Y_inverted_repeats_v1", &mut |_, _| {})
            .await
            .ok()?;
        let azf_dyz = self
            .gateway
            .resolve_mask("chm13v2.0Y_AZF_DYZ_v1", &mut |_, _| {})
            .await
            .ok()?;
        let native = YStructuralRegions::from_beds(&amplicon, &palindrome, &azf_dyz).ok()?;

        let target = canonical_build(build)?;
        if matches!(target, ReferenceBuild::Chm13v2 | ReferenceBuild::Chm13v2MaskedRcrs) {
            return Some(native);
        }
        // The code reads the PAR regions and the heterochromatin natively for each build. It does
        // not move them from CHM13. A chain file is least reliable in those places. Both chrX and
        // chrY hold the PAR regions, and Yq12 is satellite sequence.
        let landmarks = navigator_analysis::mask::y_landmarks(build)?;

        self.gateway
            .resolve_chain(ReferenceBuild::Chm13v2.as_str(), target.as_str(), &mut |_, _| {})
            .await
            .ok()?;
        let lift = |m: &RegionMask, what: &str| -> Option<RegionMask> {
            let (iv, dropped) = self
                .gateway
                .lift_intervals(ReferenceBuild::Chm13v2.as_str(), target.as_str(), "chrY", m.intervals())
                .ok()?;
            if iv.is_empty() {
                // A mask that gives no interval after the liftover is not a mask. The code then
                // writes no annotation. It must not report "no structural regions here" as a
                // result of a real check.
                eprintln!("chrY {what} mask lifted CHM13→{} to nothing; skipping", target.as_str());
                return None;
            }
            if dropped > 0 {
                eprintln!(
                    "chrY {what} mask CHM13→{}: {} of {} intervals dropped as unliftable",
                    target.as_str(),
                    dropped,
                    m.intervals().len()
                );
            }
            Some(RegionMask::from_intervals(iv))
        };
        let (palindrome_m, amplicon_m) = native.structural_masks();
        Some(YStructuralRegions::from_masks(
            RegionMask::from_intervals(landmarks.par.to_vec()),
            lift(palindrome_m, "palindrome")?,
            lift(amplicon_m, "amplicon")?,
            // A satellite array rarely survives a chain file. So the heterochromatin bound of the
            // build carries this test, and the AZF and DYZ intervals from the liftover only make it
            // more exact.
            lift(native.heterochromatin_mask(), "AZF/DYZ")
                .unwrap_or_else(|| RegionMask::from_intervals(vec![]))
                .union(&[landmarks.heterochromatin]),
        ))
    }

    /// Find the private-variant candidates on chrY from the GVCF file of one sample.
    ///
    /// The method returns the same [`VariantCall`] shape that the pileup de-novo path returns. So
    /// the classification in `private_y_core` is the same for both paths.
    ///
    /// The reassembly step of GATK finds SNVs that the pileup caller does not find.
    async fn run_denovo_from_gvcf(&self, gvcf: &Path) -> Result<Vec<VariantCall>, AppError> {
        let gvcf = gvcf.to_path_buf();
        let snvs = tokio::task::spawn_blocking(move || {
            // Use min_dp 4, and not the default of 2 that the reader allows. A real private SNV
            // has 4 reads or more. A group of misaligned reads gives 2 or 3 reads across many false
            // SNVs that are near each other. So this depth limit removes each artefact group, and
            // it keeps each true call with a DP of 4 or more.
            let params = navigator_analysis::gvcf::GvcfReadParams { min_dp: 4, min_gq: 20 };
            navigator_analysis::gvcf::read_derived_snvs(&gvcf, "chrY", &params)
        })
        .await??;
        Ok(snvs
            .into_iter()
            .map(|s| VariantCall {
                contig: "chrY".to_string(),
                position: s.position,
                reference_allele: s.reference,
                alternate_allele: s.alternate,
                depth: s.depth,
                alt_depth: s.alt_depth,
                allele_fraction: s.allele_fraction,
                quality: None,
            })
            .collect())
    }

    async fn private_y_core(
        &self,
        alignment_id: i64,
        mask: Option<navigator_analysis::mask::RegionMask>,
    ) -> Result<PrivateBucket, AppError> {
        // Classify each new call against the **DecodingUs** tree. That tree has the authority for
        // a placement in this app, and it holds the branches that the de-novo tree pipeline found
        // in the cohort.
        //
        // A variant of a shared lineage has a name in that tree, so the code marks it OffPathKnown
        // and not "novel". A variant that the tree does not hold, and that the cohort shares, is
        // doubtful.
        //
        // The FTDNA tree is the second choice. It keeps the report correct when the code can not
        // read the AppView tree, and when the build has no DecodingUs coordinates.
        let (tree, tree_calls) = match self.y_decodingus_tree_calls(alignment_id).await {
            Ok(tc) => tc,
            // An absent alignment file is not a fault of the tree. The second path reads the same
            // absent file, so it also fails. It then writes a log entry that names a tree provider
            // with no fault.
            Err(e) if e.is_missing_alignment_file() => return Err(e),
            Err(e) => {
                eprintln!("DecodingUs Y tree unavailable ({e}); private-Y classifying against FTDNA");
                let tree_json = self.fetch_ftdna_y_tree().await?;
                self.tree_base_calls(alignment_id, "chrY", &tree_json).await?
            }
        };
        let assignment = assemble_assignment(&tree, &tree_calls);
        let terminal = assignment
            .ranked
            .first()
            .ok_or_else(|| AppError::Import("no Y haplogroup match".into()))?;
        let path = navigator_analysis::haplo::path_positions(&tree, terminal.id);
        let known = navigator_analysis::haplo::tree_positions(&tree);

        // The structural BED files hold CHM13 chrY coordinates, so they annotate a CHM13 alignment
        // only. The cohort masks apply to each build: the CHM13 files are native, and CrossMap
        // moved the GRCh38 files from hs1 to hg38.
        let aln = self.alignment_or_err(alignment_id).await?;
        let regions = self.y_structural_regions_for(&aln.reference_build).await;
        // Layer 2 is the **callable mask** of the cohort, in the Poznik form. A position is in that
        // mask when it is CALLABLE in 90% or more of a cohort of about 3,000 men. Only about 25% of
        // chrY outside the PAR regions is reliably callable across a cohort.
        //
        // Layer 3 is a blocklist of the **sites that the cohort shares**. It holds each position
        // that varies with two carriers or more across the cohort. It also holds the homoplasy
        // hotspots.
        //
        // A real variant of a shared lineage is in the DecodingUs tree, and the step above marks it
        // off-path-known. A variant that the cohort shares and the tree does not hold is a
        // recurrent artefact. It is not a private SNP.
        //
        // A true private variant has one carrier in the cohort, so it passes this filter. This
        // layer takes the place of the cohort carrier filter of the de-novo pipeline, for one
        // sample.
        //
        // The bundle holds a file for each build. The CHM13 file is native, and CrossMap moved the
        // GRCh38 file. The code skips this layer when the file is absent.
        let mask_token = y_mask_build_token(&aln.reference_build);
        let cohort_mask =
            mask_token.and_then(|t| load_y_position_bed("NAVIGATOR_Y_CALLABLE_MASK", "chrY_callable_mask", t));
        let cohort_shared =
            mask_token.and_then(|t| load_y_position_bed("NAVIGATOR_Y_COHORT_SHARED", "chrY_cohort_shared_sites", t));

        // The source of the derived calls. Take the chrY GVCF file of the sample first.
        //
        // GATK HaplotypeCaller builds the local haplotypes again. Take a site with about 50%
        // reference reads, where the mapper placed the reference reads wrongly. That step gives an
        // answer at such a site, and the pileup caller removes it. See the recall gap of WGS229.
        //
        // When the sample has no sidecar file, use the de-novo pileup caller of Navigator.
        let denovo = match chr_y_gvcf_for_alignment(&aln) {
            Some(gvcf) => {
                eprintln!("private-Y: sourcing chrY calls from GVCF sidecar {}", gvcf.display());
                self.run_denovo_from_gvcf(&gvcf).await?
            }
            None => self.run_denovo_for_alignment(alignment_id, "chrY".to_string()).await?,
        };
        // Then keep only off-backbone, callable, non-shared calls.
        let mut variants: Vec<PrivateVariant> = denovo
            .iter()
            .filter(|c| !path.contains(&c.position))
            .filter(|c| mask.as_ref().map_or(true, |m| m.contains(c.position))) // self-callable
            .filter(|c| cohort_mask.as_ref().map_or(true, |m| m.contains(c.position))) // L3 cohort callable mask
            .filter(|c| cohort_shared.as_ref().map_or(true, |m| !m.contains(c.position))) // L4 cohort-shared exclude
            .map(|c| PrivateVariant {
                position: c.position,
                reference: c.reference_allele,
                alternate: c.alternate_allele,
                depth: c.depth,
                alt_depth: c.alt_depth,
                allele_fraction: c.allele_fraction,
                class: match known.get(&c.position) {
                    Some(name) => PrivateClass::OffPathKnown(name.clone()),
                    None => PrivateClass::Novel,
                },
                region: regions.as_ref().and_then(|r| r.classify(c.position)),
            })
            .collect();
        variants.sort_by_key(|v| v.position);
        Ok(PrivateBucket {
            terminal: terminal.name.clone(),
            variants,
        })
    }
}

#[cfg(test)]
mod gvcf_discovery_tests {
    use super::*;

    fn alignment_at(bam: &Path) -> Alignment {
        Alignment {
            id: 1,
            sequence_run_id: 1,
            bam_path: Some(bam.to_string_lossy().into_owned()),
            reference_path: None,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            content_sha256: None,
            derived_from_alignment_id: None,
            derivation: None,
        }
    }

    /// Each case has its own directory, because these tests run at the same time.
    fn scratch(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("nav-gvcf-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn finds_the_ytree_sibling_sidecar() {
        let d = scratch("sibling");
        std::fs::write(d.join("chrYM.cram"), "").unwrap();
        std::fs::write(d.join("HG002.chrY.g.vcf.gz"), "").unwrap();
        let got = chr_y_gvcf_for_alignment(&alignment_at(&d.join("chrYM.cram"))).unwrap();
        assert_eq!(got.file_name().unwrap(), "HG002.chrY.g.vcf.gz");
    }

    #[test]
    fn finds_a_bare_named_gvcf_in_a_caller_subdirectory() {
        // The D2C layout for one run: `<run>/CP086569.2/gatk4/chrY.g.vcf.gz`. That name holds no
        // sample prefix, and the file is one directory below the alignment.
        //
        // An earlier version matched `*.chry.g.vcf.gz` beside the CRAM file only, and it found no
        // file of this kind. So each subject decoded the full chromosome.
        let d = scratch("subdir");
        std::fs::write(d.join("chrYM.cram"), "").unwrap();
        std::fs::create_dir_all(d.join("gatk4")).unwrap();
        std::fs::write(d.join("gatk4").join("chrY.g.vcf.gz"), "").unwrap();
        let got = chr_y_gvcf_for_alignment(&alignment_at(&d.join("chrYM.cram"))).unwrap();
        assert_eq!(got.file_name().unwrap(), "chrY.g.vcf.gz");
        assert_eq!(got.parent().unwrap().file_name().unwrap(), "gatk4");
    }

    #[test]
    fn the_alignments_own_directory_wins_over_a_subdirectory() {
        let d = scratch("precedence");
        std::fs::write(d.join("chrYM.cram"), "").unwrap();
        std::fs::write(d.join("HG002.chrY.g.vcf.gz"), "").unwrap();
        std::fs::create_dir_all(d.join("gatk4")).unwrap();
        std::fs::write(d.join("gatk4").join("chrY.g.vcf.gz"), "").unwrap();
        let got = chr_y_gvcf_for_alignment(&alignment_at(&d.join("chrYM.cram"))).unwrap();
        assert_eq!(got.file_name().unwrap(), "HG002.chrY.g.vcf.gz");
    }

    #[test]
    fn a_chr_m_gvcf_is_not_mistaken_for_chr_y() {
        let d = scratch("contig");
        std::fs::write(d.join("chrYM.cram"), "").unwrap();
        std::fs::create_dir_all(d.join("gatk4")).unwrap();
        std::fs::write(d.join("gatk4").join("chrM.g.vcf.gz"), "").unwrap();
        let aln = alignment_at(&d.join("chrYM.cram"));
        assert!(
            chr_y_gvcf_for_alignment(&aln).is_none(),
            "chrM must not satisfy a chrY lookup"
        );
        assert!(chr_m_gvcf_for_alignment(&aln).is_some());
    }

    #[test]
    fn absent_sidecars_yield_none() {
        let d = scratch("absent");
        std::fs::write(d.join("chrYM.cram"), "").unwrap();
        assert!(chr_y_gvcf_for_alignment(&alignment_at(&d.join("chrYM.cram"))).is_none());
    }
}

/// Read depth a VCF call must show before it can be a private-variant candidate. Matches the GVCF
/// path's floor: below this, a call is far more likely a misaligned-read cluster than a real SNV.
const VCF_PRIVATE_MIN_DP: u32 = 4;

/// The lowest genotype quality that a call can have. The value equals the `min_gq` value of the GVCF
/// path.
///
/// This value makes the gates of the two sources as similar as the evidence permits. It does not
/// make their results comparable. [`App::private_y_from_variant_set`] gives the reason: the call set
/// of a vendor comes from a different instrument.
const VCF_PRIVATE_MIN_GQ: u32 = 20;

/// The fraction of derived reads that a call needs to be **deterministic** on a haploid chromosome.
///
/// A cell holds one copy of chrY. So a true call has almost no reference read. A fraction between
/// the two states marks a locus with two answers, and such a call can not support a claim about a
/// private variant.
const VCF_PRIVATE_MIN_AF: f64 = 0.95;

/// The maximum depth of a call, as a factor of the usual depth of that donor at a good call.
///
/// A cell holds one copy of chrY. A locus with many more reads than the rest of the chromosome takes
/// those reads from another place. That place is a collapsed repeat.
///
/// A review of one candidate branch found this fault. Its two carriers had a DP of 413 and 504,
/// against a median of 57. Each one also held about 5% reference reads. That shape is a pile of
/// reads from a paralog, and the call passed each other gate.
///
/// A limit of three times the median keeps about 91% of the good calls and removes that group.
const VCF_PRIVATE_MAX_DEPTH_RATIO: u32 = 3;

/// The count of good calls that the code needs before it can trust a depth ratio. Below this count,
/// the median does not describe the coverage of the donor. The rule then does nothing, and it makes
/// no comparison against noise.
const VCF_PRIVATE_MIN_CALLS_FOR_RATIO: usize = 20;

impl App {
    /// The **private bucket of a variant set**. It is the VCF form of
    /// [`Self::private_y_variants`].
    ///
    /// The key of the private-Y data was always an alignment. The code walks a BAM file or a CRAM
    /// file, or its GVCF sidecar, and it caches the result under `alignment_id`.
    ///
    /// A subject whose Y data arrived as a VCF file from another tool has no alignment. So the app
    /// never offered this option to that subject. On R1b-CTS4466Plus, about 1,600 of the 1,881
    /// members are in that group. For this reason, each cohort feature that needs private variants
    /// had almost no data.
    ///
    /// The classification is the same as the classification of the alignment path, by design. The
    /// code removes the backbone that it placed. It then removes each call outside the callable mask
    /// of the cohort, and each call on the blocklist of the cohort. It then separates the known
    /// off-path calls from the new ones.
    ///
    /// Two things are different: the source of the evidence, and the test of the reliability of the
    /// donor.
    ///
    /// - **The placement** uses [`Self::vset_base_calls`]. So the terminal comes from the genotypes
    ///   at each tree position, and those genotypes include the hom-ref calls. It does not come from
    ///   the few derived calls alone.
    /// - **There is no self-callable mask**, because a VCF file holds no coverage track. The
    ///   evidence of each call takes its place. The code removes a call with a `FILTER` flag, a call
    ///   below [`VCF_PRIVATE_MIN_DP`], and each heterozygous call on chrY. A cell holds one copy of
    ///   chrY, so such a call comes from a paralog or from a read that the mapper placed wrongly.
    ///   Those calls are about two thirds of the chrY rows of a Big Y test.
    /// - **The depth and the allele fraction come from the source**, so [`PublishGate`] judges these
    ///   calls on real read evidence. A set that the app imported before it stored the evidence has
    ///   a `call_schema` of 1. Such a set gives nothing that the app can publish, and that result is
    ///   correct.
    ///
    /// **These counts do not compare with the counts of the alignment path, and they are not yet
    /// good enough for branch inference.**
    ///
    /// This path gives a median of about 175 new calls for each donor. The GVCF path gives 3 to 13.
    /// The difference is the instrument and not a fault here. The alignment path reads GATK
    /// HaplotypeCaller at ploidy 1. A vendor export comes from a diploid caller, and that caller
    /// writes many more chrY calls. A match of the DP gates and the GQ gates removes only about 10%
    /// of the difference.
    ///
    /// Note also that `Novel` means "this variant defines no branch in *this* tree". The tree holds
    /// the branches that FTDNA supports, and the splits that the cohort resolved. It is not a
    /// catalogue of each known Y variant. So a real, well-known variant that defines no branch is
    /// `Novel` here.
    ///
    /// The block tree read these buckets for its candidate detection. On CTS4466 the count of
    /// candidates went from 3 to 20, with 39 conflicts and 105 recurrent positions removed. Only 111
    /// of about 1,600 sets took part. For this reason,
    /// [`Self::private_y_for_biosamples`] does not yet join them.
    pub async fn private_y_from_variant_set(&self, set: &VariantSet) -> Result<PrivateBucket, AppError> {
        use navigator_analysis::haplo;

        // With no evidence for each call, every quality gate below does nothing. The result is
        // then the list that the caller of the vendor wrote. On a real set that list holds 400 to
        // 550 "new" calls. The set of the same donor with evidence holds about 70.
        //
        // A call that the app can not judge is the least deterministic call of all. So the method
        // refuses, and it does not publish a number that looks like a result.
        //
        // A second import of the source writes the `CallEvidence` rows. See migration 0042.
        if !set.has_evidence() {
            return Err(AppError::Import(format!(
                "variant set {} carries no per-call evidence (call_schema {}); re-import it to enable private-Y",
                set.id, set.call_schema
            )));
        }
        let build = set.reference_build.clone().unwrap_or_else(|| "GRCh38".to_string());
        let tree = self.chip_y_tree(&build).await?;

        let cache_key = {
            let targets: HashSet<i64> = tree
                .nodes
                .values()
                .flat_map(|n| n.loci.iter().map(|l| l.position))
                .collect();
            // The `pv2` key. The code now moves the chrY structural masks to the build of the
            // set.
            //
            // A `pv1` bucket used **no** structural mask on any build except CHM13, so its counts
            // are too high. The version number removes the old value from the cache, and the
            // `--force` option can not reach this cache.
            format!("pv2:{}", crate::haplogroup::genotype_cache_key("chrY", None, &targets))
        };
        if let Ok(Some(json)) = variant_set_private_y::get(self.store.pool(), set.id, &cache_key).await {
            if let Ok(bucket) = serde_json::from_str::<PrivateBucket>(&json) {
                return Ok(bucket);
            }
        }

        // Place the set to get its backbone. Tree-position genotypes (with ancestral evidence) place
        // far deeper than the derived calls alone would.
        let genotypes = self.vset_base_calls(set, "chrY", &tree).await;
        let assignment = assemble_assignment(&tree, &genotypes);
        let Some(terminal) = assignment.ranked.first() else {
            return Err(AppError::Import("no Y haplogroup match for this call set".into()));
        };
        let path = haplo::path_positions(&tree, terminal.id);
        let known = haplo::tree_positions(&tree);

        let mask_token = y_mask_build_token(&build);
        let cohort_mask =
            mask_token.and_then(|t| load_y_position_bed("NAVIGATOR_Y_CALLABLE_MASK", "chrY_callable_mask", t));
        let cohort_shared =
            mask_token.and_then(|t| load_y_position_bed("NAVIGATOR_Y_COHORT_SHARED", "chrY_cohort_shared_sites", t));
        // The structural BED files are native to CHM13, and the code moves them to the build of
        // this set. Without that step, a GRCh38 set has no structural mask, and its count of private
        // variants grows into the hundreds.
        let regions = self.y_structural_regions_for(&build).await;

        // Take the calls that pass the quality gates first. The depth limit below then compares
        // against the good coverage of the donor. Without this order, it compares against a median
        // that the calls below the gates make smaller.
        let passing: Vec<&navigator_domain::variants::VariantCall> = set
            .calls
            .iter()
            .filter(|c| c.contig.eq_ignore_ascii_case("chrY") || c.contig.eq_ignore_ascii_case("y"))
            .filter(|c| !path.contains(&c.position))
            .filter(|c| cohort_mask.as_ref().map_or(true, |m| m.contains(c.position)))
            .filter(|c| cohort_shared.as_ref().map_or(true, |m| !m.contains(c.position)))
            .filter(|c| is_hemizygous(c.genotype.as_deref()))
            .filter(|c| !c.evidence.is_filtered())
            .filter(|c| !c.evidence.dp.is_some_and(|dp| dp < VCF_PRIVATE_MIN_DP))
            .filter(|c| !c.evidence.gq.is_some_and(|gq| gq < VCF_PRIVATE_MIN_GQ))
            .filter(|c| !c.evidence.allele_fraction().is_some_and(|af| af < VCF_PRIVATE_MIN_AF))
            .collect();
        let depth_ceiling = median_depth(&passing).map(|m| m * VCF_PRIVATE_MAX_DEPTH_RATIO);

        let mut variants: Vec<PrivateVariant> = passing
            .into_iter()
            .filter(|c| match (depth_ceiling, c.evidence.dp) {
                (Some(ceiling), Some(dp)) => dp <= ceiling,
                _ => true,
            })
            .filter_map(|c| {
                let (reference, alternate) = (c.reference.chars().next()?, c.alternate.chars().next()?);
                Some(PrivateVariant {
                    position: c.position,
                    reference,
                    alternate,
                    // The numbers of the source. The value is absent when the source gave none.
                    // The publish gate then refuses that call, and that decision is correct. It
                    // must not read an absent value as evidence.
                    depth: c.evidence.dp.unwrap_or(0),
                    alt_depth: c.evidence.ad_alt.unwrap_or(0),
                    allele_fraction: c.evidence.allele_fraction().unwrap_or(0.0),
                    class: match known.get(&c.position) {
                        Some(name) => PrivateClass::OffPathKnown(name.clone()),
                        None => PrivateClass::Novel,
                    },
                    region: regions.as_ref().and_then(|r| r.classify(c.position)),
                })
            })
            .collect();
        variants.sort_by_key(|v| v.position);

        let bucket = PrivateBucket {
            terminal: terminal.name.clone(),
            variants,
        };
        if let Ok(json) = serde_json::to_string(&bucket) {
            let _ = variant_set_private_y::upsert(self.store.pool(), set.id, &cache_key, &json).await;
        }
        Ok(bucket)
    }
}

/// Shows whether a genotype holds one allele. Such a call is hemizygous or homozygous for the
/// alternate allele.
///
/// A cell holds one copy of chrY. So a heterozygous call there is not possible in biology. It marks
/// a paralog, or a locus where the mapper placed the reads wrongly.
///
/// In a real Big Y export, those calls are about two thirds of the chrY rows. With them, most of the
/// private set is an artefact.
///
/// The function accepts a call with no genotype. A source that writes no GT field makes no statement
/// about heterozygosity.
fn is_hemizygous(gt: Option<&str>) -> bool {
    let Some(gt) = gt else { return true };
    let alleles: Vec<&str> = gt.split(['/', '|']).filter(|a| *a != ".").collect();
    !alleles.is_empty() && alleles.windows(2).all(|w| w[0] == w[1])
}

#[cfg(test)]
mod vcf_private_y_tests {
    use super::is_hemizygous;

    #[test]
    fn a_chr_y_heterozygote_is_rejected() {
        // A cell holds one copy of chrY. So a heterozygous call marks a paralog, or a read that
        // the mapper placed wrongly. Those calls are about two thirds of the chrY rows of a Big Y
        // test. With them, most of the private set is an artefact.
        assert!(!is_hemizygous(Some("0/1")));
        assert!(!is_hemizygous(Some("1|2")));
        assert!(!is_hemizygous(Some("1/2")));
    }

    #[test]
    fn hemizygous_and_homozygous_alt_calls_are_kept() {
        assert!(is_hemizygous(Some("1")));
        assert!(is_hemizygous(Some("1/1")));
        assert!(is_hemizygous(Some("1|1")));
        assert!(is_hemizygous(Some("2/2")));
    }

    #[test]
    fn a_source_that_reports_no_genotype_is_not_treated_as_heterozygous() {
        // An absent GT field makes no statement about the ploidy. A rule that refused it would
        // remove each set with sites only, and each set from a CSV file, with no message.
        assert!(is_hemizygous(None));
        assert!(is_hemizygous(Some("1/.")), "a partial call still carries one allele");
    }
}

/// The median read depth across `calls`. The function returns `None` when too few calls hold a depth
/// to describe the coverage of the donor.
///
/// The function takes the median and not the mean. The piles of reads that this code must find make
/// a mean larger, and they then hide behind that larger value.
fn median_depth(calls: &[&navigator_domain::variants::VariantCall]) -> Option<u32> {
    let mut depths: Vec<u32> = calls.iter().filter_map(|c| c.evidence.dp).collect();
    if depths.len() < VCF_PRIVATE_MIN_CALLS_FOR_RATIO {
        return None;
    }
    depths.sort_unstable();
    Some(depths[depths.len() / 2]).filter(|&m| m > 0)
}

#[cfg(test)]
mod depth_ratio_tests {
    use super::{median_depth, VCF_PRIVATE_MIN_CALLS_FOR_RATIO};
    use navigator_domain::variants::{CallEvidence, VariantCall};

    fn call(dp: Option<u32>) -> VariantCall {
        VariantCall {
            contig: "chrY".into(),
            position: 1,
            reference: "A".into(),
            alternate: "G".into(),
            rs_id: None,
            genotype: Some("1/1".into()),
            evidence: CallEvidence {
                dp,
                ..Default::default()
            },
        }
    }

    #[test]
    fn the_median_ignores_the_pile_ups_it_exists_to_find() {
        // The extreme values make a mean larger, and they then hide behind that larger value.
        let mut calls: Vec<VariantCall> = (0..VCF_PRIVATE_MIN_CALLS_FOR_RATIO).map(|_| call(Some(50))).collect();
        calls.push(call(Some(2584)));
        calls.push(call(Some(1191)));
        let refs: Vec<&VariantCall> = calls.iter().collect();
        assert_eq!(median_depth(&refs), Some(50));
    }

    #[test]
    fn it_abstains_on_too_few_calls_to_describe_coverage() {
        let calls: Vec<VariantCall> = (0..VCF_PRIVATE_MIN_CALLS_FOR_RATIO - 1)
            .map(|_| call(Some(50)))
            .collect();
        let refs: Vec<&VariantCall> = calls.iter().collect();
        assert_eq!(median_depth(&refs), None, "no ceiling rather than one built on noise");
    }

    #[test]
    fn a_source_reporting_no_depth_yields_no_ceiling() {
        // Absent depth must not read as zero and reject everything.
        let calls: Vec<VariantCall> = (0..40).map(|_| call(None)).collect();
        let refs: Vec<&VariantCall> = calls.iter().collect();
        assert_eq!(median_depth(&refs), None);
    }
}
