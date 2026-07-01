//! Headless command-line interface for the Navigator workbench. The same binary launches
//! the egui GUI when run with no subcommand; with a subcommand it opens the *same* workspace
//! database and runs scripted ingestion or read-only probes, then exits.
//!
//! This makes the workbench scriptable: bulk-load an assortment of files into a subject via
//! the unified auto-detect importer (`app.add_data`), then query the resulting rows for
//! verification (`subjects` / `show` / `projects`, with `--json` for machine consumption).
//!
//!   navigator ingest --subject "James Kane" --project mine /Volumes/nas/Genomics/mine/*
//!   navigator show --subject "James Kane" --json

use std::path::{Path, PathBuf};

use clap::{Args, Parser, Subcommand};
use navigator_app::{App, DnaType};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::workspace::NewProject;

#[derive(Parser)]
#[command(
    name = "navigator",
    version,
    about = "DUNavigator workbench (launches the GUI when run with no subcommand)"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,
}

#[derive(Subcommand)]
pub enum Command {
    /// Ingest files/directories into a subject via auto-detection (BAM/CRAM, VCF, chip, STR, mtDNA FASTA).
    Ingest(IngestArgs),
    /// List every subject with its data-source counts.
    Subjects(ProbeArgs),
    /// Show one subject's runs, alignments, profiles, and haplogroup consensus.
    Show(ShowArgs),
    /// List projects with their subject counts.
    Projects(ProbeArgs),
    /// De-novo diploid variant calling → VCF (whole-genome, or a single `--contig`).
    Call(CallArgs),
    /// Lift a VCF from one reference build to another (chain-based; the GATK LiftoverVcf replacement).
    LiftVcf(LiftVcfArgs),
    /// Run the per-alignment analysis steps with per-step wall-clock timing (the GUI Full Analysis
    /// path), to profile where time goes. Mutates the workspace (caches results).
    Analyze(AnalyzeArgs),
    /// Rebuild the genome-consensus Y signature (variant profile + descent) for subjects, e.g. to
    /// re-place existing profiles on the current tree provider. Reuses cached genotypes, so it is
    /// cheap for already-analyzed subjects. By default only subjects that already have a Y profile
    /// are rebuilt (`--all` rebuilds every subject).
    RebuildSignatures(RebuildArgs),
}

#[derive(Args)]
pub struct IngestArgs {
    /// Subject donor identifier (found by exact match, or created if absent). Mutually exclusive
    /// with --external-id; one of the two is required.
    #[arg(long, short, required_unless_present = "external_id")]
    subject: Option<String>,
    /// Resolve the subject by a vendor id `(--id-source, ID)` instead of a donor identifier — e.g.
    /// an FTDNA kit number. The subject must already exist (never created); pair with
    /// --skip-unmatched to skip unknown ids quietly.
    #[arg(long, conflicts_with = "subject")]
    external_id: Option<String>,
    /// Vendor source for --external-id (default FTDNA).
    #[arg(long, default_value = navigator_domain::identity::IdSource::FTDNA)]
    id_source: String,
    /// With --external-id: if no subject matches the id, skip quietly (exit 0) instead of erroring.
    #[arg(long)]
    skip_unmatched: bool,
    /// Force the sequencing-run test type for alignment files (e.g. "Big Y") instead of inferring
    /// it. Useful for bulk imports where the directory layout names the test; CRAMs have no `.bai`
    /// for the coverage-shape detector, so they otherwise fall back to WGS. Ignored for non-BAM/CRAM.
    #[arg(long)]
    test_type: Option<String>,
    /// Optional project name to assign the subject to (found or created).
    #[arg(long, short)]
    project: Option<String>,
    /// Sex recorded only when the subject is created (e.g. male / female).
    #[arg(long)]
    sex: Option<String>,
    /// Recurse into directories (default: only their immediate files).
    #[arg(long, short)]
    recursive: bool,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
    /// Files and/or directories to ingest.
    #[arg(required = true)]
    paths: Vec<PathBuf>,
}

#[derive(Args)]
pub struct ProbeArgs {
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
    /// Emit JSON instead of a human-readable table.
    #[arg(long)]
    json: bool,
}

#[derive(Args)]
pub struct ShowArgs {
    /// Subject donor identifier to show.
    #[arg(long, short)]
    subject: String,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
    /// Emit JSON instead of a human-readable summary.
    #[arg(long)]
    json: bool,
}

#[derive(Args)]
pub struct CallArgs {
    /// Subject donor identifier (used to resolve the alignment when `--alignment` is omitted).
    #[arg(long, short)]
    subject: Option<String>,
    /// Alignment id to call (from `show --json`). If omitted, the subject's sole alignment is used.
    #[arg(long, short)]
    alignment: Option<i64>,
    /// Restrict to a single contig (e.g. chrM, chr21). Default: every primary chromosome.
    #[arg(long, short)]
    contig: Option<String>,
    /// Write the VCF here instead of stdout.
    #[arg(long, short)]
    out: Option<PathBuf>,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
}

#[derive(Args)]
pub struct AnalyzeArgs {
    /// Alignment id to analyze (from `show --json`).
    #[arg(long, short)]
    alignment: i64,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
}

#[derive(Args)]
pub struct RebuildArgs {
    /// Rebuild every subject's Y signature, including those without one yet (default: only subjects
    /// that already have a Y profile — the ones a placement change leaves stale).
    #[arg(long)]
    all: bool,
    /// Restrict to subjects in this project (by exact name).
    #[arg(long, short)]
    project: Option<String>,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
}

#[derive(Args)]
pub struct LiftVcfArgs {
    /// Input VCF (`.vcf` or `.vcf.gz`).
    #[arg(long, short)]
    r#in: PathBuf,
    /// Target reference build to lift to (e.g. chm13v2.0, GRCh38, GRCh37).
    #[arg(long, short)]
    to: String,
    /// Source build of the input VCF (e.g. GRCh38). Inferred from the header when omitted.
    #[arg(long, short)]
    from: Option<String>,
    /// Output VCF path (`.vcf` or `.vcf.gz`).
    #[arg(long, short)]
    out: PathBuf,
    /// Drop variants landing in the target chrY PAR.
    #[arg(long)]
    filter_par: bool,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
}

/// Run a CLI subcommand to completion, returning a process exit code. Spins its own tokio
/// runtime so `main` (which must keep the GUI on the main thread) stays sync.
pub fn run(command: Command) -> i32 {
    // 16 MiB stacks: Y/mt tree parse + placement recurse to the haplotree depth and overflow
    // tokio's default 2 MiB worker stack on deep lineages (matches the GUI worker runtime).
    let rt = match tokio::runtime::Builder::new_multi_thread()
        .enable_all()
        .thread_stack_size(16 * 1024 * 1024)
        .build()
    {
        Ok(rt) => rt,
        Err(e) => {
            eprintln!("error: could not start runtime: {e}");
            return 1;
        }
    };
    rt.block_on(async move {
        match command {
            Command::Ingest(a) => ingest(a).await,
            Command::Subjects(a) => subjects(a).await,
            Command::Show(a) => show(a).await,
            Command::Projects(a) => projects(a).await,
            Command::Call(a) => call(a).await,
            Command::LiftVcf(a) => lift_vcf(a).await,
            Command::Analyze(a) => analyze(a).await,
            Command::RebuildSignatures(a) => rebuild_signatures(a).await,
        }
    })
}

/// Rebuild the genome-consensus Y signature for a set of subjects. Used to re-place existing
/// profiles after a placement/tree-provider change: the batch analyzer only *creates* a signature
/// when one is missing, so profiles built on an older tree stay stale until rebuilt here.
async fn rebuild_signatures(args: RebuildArgs) -> i32 {
    use std::time::Instant;
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };

    // Resolve the optional project filter to an id.
    let project_id = match &args.project {
        Some(name) => {
            let overview = app.project_overview().await.unwrap_or_default();
            match overview.iter().find(|o| o.project.name == *name) {
                Some(o) => Some(o.project.id),
                None => {
                    eprintln!("error: no project named \"{name}\"");
                    return 1;
                }
            }
        }
        None => None,
    };

    let bios = match app.list_all_biosamples().await {
        Ok(v) => v,
        Err(e) => return report(e),
    };

    let (mut rebuilt, mut skipped, mut failed) = (0usize, 0usize, 0usize);
    for b in &bios {
        if let Some(pid) = project_id {
            if b.project_id != Some(pid) {
                continue;
            }
        }
        // Default: only refresh subjects that already carry a Y profile (the stale ones). With
        // --all, build for every subject (skips those with no Y evidence, which just yield empty).
        if !args.all {
            match app.cached_y_profile(b.guid).await {
                Ok(Some(_)) => {}
                Ok(None) => {
                    skipped += 1;
                    continue;
                }
                Err(e) => {
                    eprintln!("FAIL {:<24} reading profile: {e}", truncate(&b.donor_identifier, 24));
                    failed += 1;
                    continue;
                }
            }
        }
        let t = Instant::now();
        match app.build_y_profile(b.guid).await {
            Ok(p) => {
                rebuilt += 1;
                println!(
                    "OK   {:<24} {:<12} {:>3} variants  [{:.1?}]",
                    truncate(&b.donor_identifier, 24),
                    p.terminal.as_deref().unwrap_or("(none)"),
                    p.variants.len(),
                    t.elapsed()
                );
            }
            Err(e) => {
                failed += 1;
                eprintln!("FAIL {:<24} {e}", truncate(&b.donor_identifier, 24));
            }
        }
    }
    println!("\nrebuilt {rebuilt} signature(s), {skipped} skipped (no profile), {failed} failed");
    if failed > 0 {
        1
    } else {
        0
    }
}

/// Time the per-alignment analysis steps (the GUI Full Analysis path) to profile where time goes.
async fn analyze(args: AnalyzeArgs) -> i32 {
    use std::time::Instant;
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    let id = args.alignment;
    eprintln!("analyzing alignment #{id} (cold — run_* recompute, exercising localize + scoping)…");

    let t = Instant::now();
    match app.run_unified_metrics(id).await {
        Ok(r) => {
            eprintln!(
                "  [{:>8.1?}]  Step 1: unified metrics (coverage + sex + read-metrics, incl. localize copy) — coverage mean {:.2}x",
                t.elapsed(),
                r.coverage.mean_coverage
            );
            // Mirror the batch path: clear any stale failure marker now that the walk succeeded.
            app.clear_analysis_error(id).await;
        }
        Err(e) => {
            eprintln!("  Step 1 (unified metrics) FAILED: {e}");
            // Persist the failure (corrupt/undecodable file) so the project report surfaces it.
            app.record_analysis_error(id, "metrics", &e.to_string()).await;
            return 1;
        }
    }

    let t = Instant::now();
    match app.assign_y_haplogroup(id).await {
        Ok(a) => eprintln!(
            "  [{:>8.1?}]  Step 3: Y haplogroup placement — {}",
            t.elapsed(),
            a.ranked.first().map(|h| h.name.as_str()).unwrap_or("(none)")
        ),
        Err(e) => eprintln!("  Step 3 (Y placement) FAILED: {e}"),
    }

    // Step 4: build the genome-consensus Y signature (deep placement + variant profile → descent
    // report), mirroring batch analysis so the descent report is ready without an explicit click.
    if let Ok(bio) = app.biosample_of_alignment(id).await {
        let t = Instant::now();
        match app.build_y_profile(bio).await {
            Ok(p) => eprintln!(
                "  [{:>8.1?}]  Step 4: Y signature — terminal {} · {} variants",
                t.elapsed(),
                p.terminal.as_deref().unwrap_or("(none)"),
                p.variants.len()
            ),
            Err(e) => eprintln!("  Step 4 (Y signature) FAILED: {e}"),
        }
    }
    eprintln!("done.");
    0
}

fn db_path(over: Option<PathBuf>) -> PathBuf {
    over.unwrap_or_else(crate::default_db_path)
}

async fn open(db: Option<PathBuf>) -> Result<App, i32> {
    let path = db_path(db);
    App::open(&path).await.map_err(|e| {
        eprintln!("error: could not open workspace {}: {e}", path.display());
        1
    })
}

/// Find a subject by exact donor identifier, returning its guid if present.
async fn find_subject(app: &App, donor: &str) -> Result<Option<SampleGuid>, i32> {
    let all = app.list_all_biosamples().await.map_err(report)?;
    Ok(all.into_iter().find(|b| b.donor_identifier == donor).map(|b| b.guid))
}

/// Find a project id by exact name, or create it.
async fn find_or_create_project(app: &App, name: &str) -> Result<i64, i32> {
    let overview = app.project_overview().await.map_err(report)?;
    if let Some(o) = overview.iter().find(|o| o.project.name == name) {
        return Ok(o.project.id);
    }
    let p = app
        .create_project(NewProject {
            name: name.to_string(),
            description: None,
            administrator: "cli".into(),
        })
        .await
        .map_err(report)?;
    Ok(p.id)
}

fn report(e: navigator_app::AppError) -> i32 {
    eprintln!("error: {e}");
    1
}

/// Expand the path list into a sorted file list, descending into directories (one level, or
/// fully with `recursive`). Hidden dotfiles are skipped.
fn collect_files(paths: &[PathBuf], recursive: bool) -> Vec<PathBuf> {
    let mut out = Vec::new();
    for p in paths {
        push_path(p, recursive, &mut out);
    }
    out.sort();
    out.dedup();
    out
}

fn is_hidden(p: &Path) -> bool {
    p.file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|n| n.starts_with('.'))
}

fn push_path(p: &Path, recursive: bool, out: &mut Vec<PathBuf>) {
    if is_hidden(p) {
        return;
    }
    if p.is_dir() {
        let Ok(entries) = std::fs::read_dir(p) else {
            eprintln!("warning: cannot read directory {}", p.display());
            return;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if is_hidden(&path) {
                continue;
            }
            if path.is_dir() {
                if recursive {
                    push_path(&path, recursive, out);
                }
            } else {
                out.push(path);
            }
        }
    } else if p.is_file() {
        out.push(p.to_path_buf());
    } else {
        eprintln!("warning: skipping {} (not a file or directory)", p.display());
    }
}

async fn ingest(args: IngestArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };

    // Resolve the project first so subject creation can attach to it.
    let project_id = match &args.project {
        Some(name) => match find_or_create_project(&app, name).await {
            Ok(id) => Some(id),
            Err(c) => return c,
        },
        None => None,
    };

    // Resolve the target subject. Two modes:
    //   --external-id : match an existing subject by vendor id; NEVER create (skip or error).
    //   --subject     : find-or-create by donor identifier (the original behaviour).
    // `label` is the human identifier used in the summary line.
    let (guid, label) = if let Some(ext) = &args.external_id {
        match app.find_biosample_by_external_id(&args.id_source, ext).await {
            Ok(Some(g)) => (g, format!("{}:{ext}", args.id_source)),
            Ok(None) => {
                if args.skip_unmatched {
                    println!("SKIP no subject for {}:{ext}", args.id_source);
                    return 0;
                }
                eprintln!("error: no subject with {} id \"{ext}\"", args.id_source);
                return 1;
            }
            Err(e) => return report(e),
        }
    } else {
        // --subject is required-unless-present(external_id), so it is Some here.
        let subject = args.subject.clone().unwrap_or_default();
        match find_subject(&app, &subject).await {
            Ok(Some(g)) => (g, subject),
            Ok(None) => match app
                .add_biosample(project_id, subject.clone(), None, args.sex.clone())
                .await
            {
                Ok(b) => {
                    println!("created subject {subject} ({})", b.guid.0);
                    (b.guid, subject)
                }
                Err(e) => return report(e),
            },
            Err(c) => return c,
        }
    };

    // Assign to the named project (applies to both resolution modes).
    if let Some(pid) = project_id {
        if let Err(e) = app.add_biosample_to_project(guid, Some(pid)).await {
            return report(e);
        }
    }

    let files = collect_files(&args.paths, args.recursive);
    if files.is_empty() {
        eprintln!("error: no files found in the given paths");
        return 1;
    }

    let (mut ok, mut failed, mut ysnp_panels) = (0usize, 0usize, 0usize);
    for path in &files {
        match app.add_data_with_test_type(guid, path, args.test_type.as_deref()).await {
            Ok(detected) => {
                ok += 1;
                if detected == navigator_app::DetectedData::YSnpPanel {
                    ysnp_panels += 1;
                }
                println!("OK   {:<18} {}", detected.description(), path.display());
            }
            Err(e) => {
                failed += 1;
                eprintln!("FAIL {:<18} {}: {e}", "—", path.display());
            }
        }
    }
    println!("\ningested {ok} file(s), {failed} failed, into subject \"{label}\"");

    // A Y-SNP panel (BISDNA) was imported — place a Y haplogroup from its derived calls and
    // report the terminal (the call is recorded for the donor consensus).
    if ysnp_panels > 0 {
        match app.assign_y_bisdna(guid, None).await {
            Ok(a) => match a.ranked.first() {
                Some(top) => println!(
                    "Y-DNA (panel): {} (score {:.3}, {}/{} defining SNPs)",
                    top.name, top.score, top.matched, top.expected
                ),
                None => println!("Y-DNA (panel): no haplogroup match"),
            },
            Err(e) => eprintln!("warning: Y-SNP panel placement failed: {e}"),
        }
    }
    if failed > 0 {
        1
    } else {
        0
    }
}

async fn subjects(args: ProbeArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    let bios = match app.list_all_biosamples().await {
        Ok(v) => v,
        Err(e) => return report(e),
    };
    let overview = app.project_overview().await.unwrap_or_default();
    let project_name = |id: Option<i64>| -> Option<String> {
        id.and_then(|pid| {
            overview
                .iter()
                .find(|o| o.project.id == pid)
                .map(|o| o.project.name.clone())
        })
    };

    let mut rows = Vec::new();
    for b in &bios {
        let runs = app.list_sequence_runs(b.guid).await.map(|v| v.len()).unwrap_or(0);
        let mut alns = 0usize;
        if let Ok(rs) = app.list_sequence_runs(b.guid).await {
            for r in rs {
                alns += app.list_alignments(r.id).await.map(|v| v.len()).unwrap_or(0);
            }
        }
        let strs = app.list_str_profiles(b.guid).await.map(|v| v.len()).unwrap_or(0);
        let vars = app.list_variant_sets(b.guid).await.map(|v| v.len()).unwrap_or(0);
        let chips = app.list_chip_profiles(b.guid).await.map(|v| v.len()).unwrap_or(0);
        let mt = app.list_mtdna_sequences(b.guid).await.map(|v| v.len()).unwrap_or(0);
        rows.push((b, runs, alns, strs, vars, chips, mt));
    }

    if args.json {
        let arr: Vec<_> = rows
            .iter()
            .map(|(b, runs, alns, strs, vars, chips, mt)| {
                serde_json::json!({
                    "guid": b.guid.0.to_string(),
                    "donor_identifier": b.donor_identifier,
                    "sample_accession": b.sample_accession,
                    "sex": b.sex,
                    "project": project_name(b.project_id),
                    "runs": runs, "alignments": alns,
                    "str_profiles": strs, "variant_sets": vars,
                    "chip_profiles": chips, "mtdna_sequences": mt,
                })
            })
            .collect();
        println!("{}", serde_json::to_string_pretty(&arr).unwrap());
    } else if rows.is_empty() {
        println!("(no subjects)");
    } else {
        println!("{:<24} {:<16} runs aln str var chip mt", "SUBJECT", "PROJECT");
        for (b, runs, alns, strs, vars, chips, mt) in &rows {
            println!(
                "{:<24} {:<16} {:>4} {:>3} {:>3} {:>3} {:>4} {:>2}",
                truncate(&b.donor_identifier, 24),
                truncate(project_name(b.project_id).as_deref().unwrap_or("—"), 16),
                runs,
                alns,
                strs,
                vars,
                chips,
                mt
            );
        }
    }
    0
}

async fn show(args: ShowArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    let guid = match find_subject(&app, &args.subject).await {
        Ok(Some(g)) => g,
        Ok(None) => {
            eprintln!("error: no subject with identifier \"{}\"", args.subject);
            return 1;
        }
        Err(c) => return c,
    };

    let runs = app.list_sequence_runs(guid).await.unwrap_or_default();
    let strs = app.list_str_profiles(guid).await.unwrap_or_default();
    let vars = app.list_variant_sets(guid).await.unwrap_or_default();
    let chips = app.list_chip_profiles(guid).await.unwrap_or_default();
    let mt = app.list_mtdna_sequences(guid).await.unwrap_or_default();
    let y_cons = app.haplogroup_consensus(guid, DnaType::Y).await.ok().flatten();
    let mt_cons = app.haplogroup_consensus(guid, DnaType::Mt).await.ok().flatten();

    if args.json {
        let mut runs_json = Vec::new();
        for r in &runs {
            let alns = app.list_alignments(r.id).await.unwrap_or_default();
            runs_json.push(serde_json::json!({
                "id": r.id, "test_type": r.test_type, "platform": r.platform_name,
                "instrument": r.instrument_model, "library_layout": r.library_layout,
                "alignments": alns.iter().map(|a| serde_json::json!({
                    "id": a.id, "reference_build": a.reference_build, "aligner": a.aligner,
                    "variant_caller": a.variant_caller, "bam_path": a.bam_path,
                })).collect::<Vec<_>>(),
            }));
        }
        let out = serde_json::json!({
            "guid": guid.0.to_string(),
            "donor_identifier": args.subject,
            "runs": runs_json,
            "str_profiles": strs.iter().map(|p| serde_json::json!({"id": p.id, "panel": p.panel_name, "markers": p.markers.len()})).collect::<Vec<_>>(),
            "variant_sets": vars.iter().map(|s| serde_json::json!({"id": s.id, "source": s.source_label, "calls": s.calls.len()})).collect::<Vec<_>>(),
            "chip_profiles": chips.iter().map(|c| serde_json::json!({"id": c.id, "provider": c.provider})).collect::<Vec<_>>(),
            "mtdna_sequences": mt.iter().map(|m| serde_json::json!({"id": m.id, "length": m.length()})).collect::<Vec<_>>(),
            "y_haplogroup": y_cons.as_ref().map(|c| c.haplogroup.clone()),
            "mt_haplogroup": mt_cons.as_ref().map(|c| c.haplogroup.clone()),
        });
        println!("{}", serde_json::to_string_pretty(&out).unwrap());
        return 0;
    }

    println!("Subject: {} ({})", args.subject, guid.0);
    if let Some(c) = &y_cons {
        println!("  Y-DNA : {}", c.haplogroup);
    }
    if let Some(c) = &mt_cons {
        println!("  mtDNA : {}", c.haplogroup);
    }
    println!("\nSequencing runs ({}):", runs.len());
    for r in &runs {
        println!(
            "  #{} {} · {} · {}",
            r.id,
            r.test_type,
            r.platform_name,
            r.instrument_model.as_deref().unwrap_or("—")
        );
        let alns = app.list_alignments(r.id).await.unwrap_or_default();
        for a in &alns {
            println!(
                "      aln #{} {} / {}{}",
                a.id,
                a.reference_build,
                a.aligner,
                a.bam_path.as_deref().map(|p| format!("  [{p}]")).unwrap_or_default()
            );
        }
    }
    println!(
        "\nProfiles: {} STR, {} variant-set, {} chip, {} mtDNA",
        strs.len(),
        vars.len(),
        chips.len(),
        mt.len()
    );
    for p in &strs {
        println!("  STR  #{} {} ({} markers)", p.id, p.panel_name, p.markers.len());
    }
    for s in &vars {
        println!("  VAR  #{} {} ({} calls)", s.id, s.source_label, s.calls.len());
    }
    for c in &chips {
        println!("  CHIP #{} {}", c.id, c.provider);
    }
    for m in &mt {
        println!("  MT   #{} ({} bp)", m.id, m.length());
    }
    0
}

/// Resolve the alignment id to call: the explicit `--alignment`, else the subject's sole alignment.
async fn resolve_alignment(app: &App, subject: Option<&str>, explicit: Option<i64>) -> Result<i64, i32> {
    if let Some(id) = explicit {
        return Ok(id);
    }
    let Some(donor) = subject else {
        eprintln!("error: pass --alignment <id> or --subject <donor>");
        return Err(1);
    };
    let guid = match find_subject(app, donor).await? {
        Some(g) => g,
        None => {
            eprintln!("error: no subject with identifier \"{donor}\"");
            return Err(1);
        }
    };
    let runs = app.list_sequence_runs(guid).await.map_err(report)?;
    let mut alns = Vec::new();
    for r in &runs {
        alns.extend(app.list_alignments(r.id).await.map_err(report)?);
    }
    match alns.as_slice() {
        [] => {
            eprintln!("error: subject \"{donor}\" has no alignments");
            Err(1)
        }
        [a] => Ok(a.id),
        many => {
            eprintln!(
                "error: subject \"{donor}\" has {} alignments — pass --alignment <id> (one of: {})",
                many.len(),
                many.iter().map(|a| a.id.to_string()).collect::<Vec<_>>().join(", ")
            );
            Err(1)
        }
    }
}

async fn call(args: CallArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    let alignment_id = match resolve_alignment(&app, args.subject.as_deref(), args.alignment).await {
        Ok(id) => id,
        Err(c) => return c,
    };

    let scope = args.contig.clone().unwrap_or_else(|| "whole genome".into());
    eprintln!("calling de-novo diploid variants on alignment #{alignment_id} ({scope})…");
    let vcf = match args.contig {
        Some(contig) => app.diploid_vcf(alignment_id, contig).await,
        None => app.diploid_vcf_genome(alignment_id).await,
    };
    let vcf = match vcf {
        Ok(v) => v,
        Err(e) => return report(e),
    };

    // Summary to stderr (records, of which multiallelic) so a redirected stdout stays pure VCF.
    let records: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
    let multiallelic = records
        .iter()
        .filter(|l| l.split('\t').nth(4).is_some_and(|alt| alt.contains(',')))
        .count();
    eprintln!("{} variant record(s), {multiallelic} multiallelic", records.len());

    match &args.out {
        Some(path) => {
            if let Err(e) = std::fs::write(path, &vcf) {
                eprintln!("error: writing {}: {e}", path.display());
                return 1;
            }
            eprintln!("wrote {}", path.display());
        }
        None => print!("{vcf}"),
    }
    0
}

async fn lift_vcf(args: LiftVcfArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    let source = match args
        .from
        .clone()
        .or_else(|| navigator_app::infer_vcf_source_build(&args.r#in))
    {
        Some(s) => s,
        None => {
            eprintln!(
                "error: could not infer the source build from {} — pass --from <build>",
                args.r#in.display()
            );
            return 1;
        }
    };
    eprintln!("lifting {} : {} → {} …", args.r#in.display(), source, args.to);
    let opts = navigator_app::VcfLiftOpts {
        filter_par: args.filter_par,
    };
    let mut last = 0u64;
    let mut progress = |received: u64, _total: Option<u64>| {
        // Coarse byte ticks during any chain/reference download.
        if received >= last + 50_000_000 {
            last = received;
            eprint!(".");
        }
    };
    match app
        .lift_vcf(
            &source,
            &args.to,
            args.r#in.clone(),
            args.out.clone(),
            opts,
            &mut progress,
        )
        .await
    {
        Ok(s) => {
            eprintln!(
                "\nlifted {}/{} ({} unmapped, {} split, {} ref-mismatch, {} swap-ambiguous, {} complex-rev, {} PAR) → {}",
                s.lifted, s.total, s.unmapped, s.split, s.ref_mismatch, s.swap_ambiguous, s.complex_reverse, s.par,
                args.out.display()
            );
            0
        }
        Err(e) => report(e),
    }
}

async fn projects(args: ProbeArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    let overview = match app.project_overview().await {
        Ok(v) => v,
        Err(e) => return report(e),
    };
    if args.json {
        let arr: Vec<_> = overview
            .iter()
            .map(|o| {
                serde_json::json!({
                    "id": o.project.id, "name": o.project.name,
                    "administrator": o.project.administrator,
                    "description": o.project.description,
                    "subjects": o.sample_count,
                })
            })
            .collect();
        println!("{}", serde_json::to_string_pretty(&arr).unwrap());
    } else if overview.is_empty() {
        println!("(no projects)");
    } else {
        println!("{:<6} {:<24} {:<16} subjects", "ID", "NAME", "ADMIN");
        for o in &overview {
            println!(
                "{:<6} {:<24} {:<16} {}",
                o.project.id,
                truncate(&o.project.name, 24),
                truncate(&o.project.administrator, 16),
                o.sample_count
            );
        }
    }
    0
}

fn truncate(s: &str, max: usize) -> String {
    if s.chars().count() <= max {
        s.to_string()
    } else {
        format!("{}…", s.chars().take(max - 1).collect::<String>())
    }
}
