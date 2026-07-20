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

use std::path::PathBuf;

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
    /// Diagnostic: trace the genome-consensus Y placement for one subject (candidates + lineage tally).
    DebugPlace(ShowArgs),
    /// Diagnostic: dump the Y descent SNP-by-SNP — state + observed base vs the tree's per-build polarity.
    DebugDescent(ShowArgs),
    /// Diagnostic: dump the raw read pileup (ref + A/C/G/T) behind each lineage call for one alignment.
    DebugCalls(DebugCallsArgs),
    /// Diagnostic: the filtered private-Y bucket for an alignment — DISPLAY vs PUBLISH counts.
    PrivateY(DebugCallsArgs),
    /// Per-marker branch report: the sample's genotype at every defining marker of a Y/mtDNA tree
    /// node's descendant subtree (observed base + derived/ancestral status + evidence). For
    /// spot-checking placement and exchanging observations. Table by default; `--tsv` / `--json`.
    BranchReport(BranchReportArgs),
    /// Diagnostic: explain why an alignment cannot be read. Probes the BAM/CRAM, its coordinate
    /// index, the reference FASTA and that FASTA's `.fai` **separately**, so a failure names the
    /// file actually at fault instead of whichever path the failing call happened to be handed —
    /// and reports the raw errno, which on macOS is the only thing distinguishing a privacy (TCC)
    /// denial from a Unix permission denial. Prints a report meant for pasting into a bug report.
    /// Use `--file` to check a file that was never imported. Exits non-zero if a check failed.
    Doctor(DoctorArgs),
    /// List projects with their subject counts.
    Projects(ProbeArgs),
    /// De-novo diploid variant calling → VCF (whole-genome, or a single `--contig`).
    Call(CallArgs),
    /// Lift a VCF from one reference build to another (chain-based; the GATK LiftoverVcf replacement).
    LiftVcf(LiftVcfArgs),
    /// Run the per-alignment analysis steps with per-step wall-clock timing (the GUI Full Analysis
    /// path), to profile where time goes. Mutates the workspace (caches results).
    Analyze(AnalyzeArgs),
    /// Rebuild the genome-consensus Y and mtDNA signatures (variant profile + descent) for subjects,
    /// e.g. to re-place existing profiles on the current tree provider. Reuses cached genotypes, so it
    /// is cheap for already-analyzed subjects. By default only subjects that already have a Y or mt
    /// profile are rebuilt (`--all` rebuilds every subject).
    RebuildSignatures(RebuildArgs),
    /// Backfill the standardized-test-label read-profile fields (`total_bases`, `read_type`) on runs
    /// imported before those fields existed. `total_bases` is recovered for free from cached
    /// read-metrics; `read_type` is inferred from platform/test-type, with `--rescan` reading a
    /// bounded prefix of the alignment to tell HiFi from CLR on generic-WGS PacBio runs.
    BackfillProfiles(BackfillArgs),
    /// Delete orphaned alignment (coverage-summary) records from the signed-in account's PDS — the
    /// duplicates left by the old create-race (two records for one alignment). Dry-run by default;
    /// pass `--apply` to actually delete. Requires being signed in.
    PruneOrphans(PruneArgs),
    /// Sign in to a PDS account via OAuth (opens a browser, waits for the loopback callback) and
    /// persist the session so the other subcommands (publish, prune-orphans) can authenticate.
    Login(LoginArgs),
    /// Attach public-catalog external ids (IGSR/HGDP/INSDC) derivable from each subject's local
    /// provenance, so bulk-imported public datasets publish ids that match their AppView catalog
    /// rows. Dry-run by default; pass `--apply` to write. Local-only (no PDS writes).
    BackfillCatalogIds(CatalogArgs),
    /// Resolve subjects against the AppView samples API and attach, in one pass, the catalog name id
    /// (IGSR/HGDP) plus the authoritative INSDC accession it returns (`SAMN…`→BIOSAMPLE, `ERS…`→ENA,
    /// `SRS…`→SRA), correcting the local placeholder. Dry-run by default; `--apply` to write. Only
    /// queries recognizable catalog aliases unless `--all`. Local writes only (no PDS).
    BackfillAccessions(AccessionArgs),
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
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
    /// Files and/or directories to ingest. A directory is one staged sample (sidecar fast path);
    /// a file is imported on its own.
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
pub struct DebugCallsArgs {
    /// Subject donor identifier (used to pick an alignment when `--alignment` is omitted — prefers a
    /// CHM13/HiFi alignment, else the first).
    #[arg(long, short)]
    subject: Option<String>,
    /// Alignment id to genotype (from `show --json`). Takes precedence over `--subject`.
    #[arg(long, short)]
    alignment: Option<i64>,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
}

#[derive(Args)]
pub struct BranchReportArgs {
    /// Subject donor identifier (used to pick a Y/mt alignment when `--alignment` is omitted —
    /// prefers a CHM13/HiFi alignment, else the first).
    #[arg(long, short)]
    subject: Option<String>,
    /// Alignment id to genotype (from `show --json`). Takes precedence over `--subject`.
    #[arg(long, short)]
    alignment: Option<i64>,
    /// Node to report: a haplogroup name (`R-FGC29071`) or a defining marker (`FGC29071`). The
    /// report covers this node's descendant subtree.
    #[arg(long, short)]
    node: String,
    /// Which tree to read: `y` or `mt`.
    #[arg(long, short, default_value = "y")]
    tree: String,
    /// Limit descent to this many levels below the node (default: the whole subtree).
    #[arg(long)]
    depth: Option<usize>,
    /// Write the TSV here instead of printing a table.
    #[arg(long)]
    tsv: Option<PathBuf>,
    /// Emit JSON instead of a table. Mutually exclusive with `--tsv`.
    #[arg(long, conflicts_with = "tsv")]
    json: bool,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
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
    /// Rebuild every subject's signatures, including those without one yet (default: only subjects
    /// that already have a Y or mt profile — the ones a placement change leaves stale).
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
pub struct BackfillArgs {
    /// Also read a bounded prefix of the alignment file to resolve `read_type` on runs the cheap
    /// platform/test-type inference can't (generic-`WGS` PacBio: HiFi vs CLR). Touches the files.
    #[arg(long)]
    rescan: bool,
    /// Restrict to subjects in this project (by exact name).
    #[arg(long, short)]
    project: Option<String>,
    /// Emit the per-field counts as JSON.
    #[arg(long)]
    json: bool,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
}

#[derive(Args)]
pub struct AccessionArgs {
    /// Actually attach the accessions and correct the local `sample_accession`. Without this it's a
    /// dry run (queries the API read-only, writes nothing).
    #[arg(long)]
    apply: bool,
    /// Query every subject, not just those whose name is a recognizable catalog alias (IGSR/HGDP).
    #[arg(long)]
    all: bool,
    /// Restrict to subjects in this project (by exact name).
    #[arg(long, short)]
    project: Option<String>,
    /// Cap how many subjects are queried (for a bounded test run).
    #[arg(long)]
    limit: Option<usize>,
    /// Emit the outcome as JSON.
    #[arg(long)]
    json: bool,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
}

#[derive(Args)]
pub struct CatalogArgs {
    /// Actually write the derived ids. Without this flag the command is a dry run (reports counts,
    /// writes nothing).
    #[arg(long)]
    apply: bool,
    /// Restrict to subjects in this project (by exact name).
    #[arg(long, short)]
    project: Option<String>,
    /// Emit the outcome as JSON.
    #[arg(long)]
    json: bool,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
}

#[derive(Args)]
pub struct LoginArgs {
    /// Account handle or DID (e.g. `jameskane.blog`).
    #[arg(required = true)]
    handle: String,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
}

#[derive(Args)]
pub struct PruneArgs {
    /// Actually delete the orphans. Without this flag the command is a dry run (lists what it would
    /// remove and touches nothing) — a PDS delete is irreversible, so it's opt-in.
    #[arg(long)]
    apply: bool,
    /// Emit the outcome as JSON.
    #[arg(long)]
    json: bool,
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
    // 64 MiB stacks: Y/mt tree parse + placement recurse to the haplotree depth, and noodles' CRAM
    // decoder recurses on `spawn_blocking` decode paths (deepest on CRAM 3.1) — either overflows
    // tokio's default 2 MiB stack and aborts the process (matches the GUI worker runtime). See
    // `NAVIGATOR_DECODE_STACK_MB`.
    let rt = match tokio::runtime::Builder::new_multi_thread()
        .enable_all()
        .thread_stack_size(64 * 1024 * 1024)
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
            Command::DebugPlace(a) => debug_place(a).await,
            Command::DebugDescent(a) => debug_descent(a).await,
            Command::DebugCalls(a) => debug_calls(a).await,
            Command::PrivateY(a) => private_y(a).await,
            Command::BranchReport(a) => branch_report(a).await,
            Command::Doctor(a) => doctor(a).await,
            Command::Projects(a) => projects(a).await,
            Command::Call(a) => call(a).await,
            Command::LiftVcf(a) => lift_vcf(a).await,
            Command::Analyze(a) => analyze(a).await,
            Command::RebuildSignatures(a) => rebuild_signatures(a).await,
            Command::BackfillProfiles(a) => backfill_profiles(a).await,
            Command::PruneOrphans(a) => prune_orphans(a).await,
            Command::Login(a) => login(a).await,
            Command::BackfillCatalogIds(a) => backfill_catalog_ids(a).await,
            Command::BackfillAccessions(a) => backfill_accessions(a).await,
        }
    })
}

/// Resolve subjects against the AppView samples API and attach the authoritative INSDC accessions.
async fn backfill_accessions(args: AccessionArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
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
    let r = match app.backfill_accessions(project_id, args.apply, args.all, args.limit).await {
        Ok(r) => r,
        Err(e) => return report(e),
    };
    if args.json {
        println!("{}", serde_json::to_string_pretty(&r).unwrap_or_default());
    } else {
        println!(
            "Queried {} subject(s): {} resolved, {} not found, {} error(s).",
            r.examined, r.resolved, r.not_found, r.errors
        );
        if r.applied {
            println!("  ids attached (name + accession): {}", r.ids_added);
            println!("  local accession fixed:           {}", r.accession_updated);
            if r.conflicts > 0 {
                println!("  conflicts:                       {} (id already owned by another subject)", r.conflicts);
            }
        } else {
            println!("  ids to attach (name + accession): {} (dry run — pass --apply)", r.ids_to_add);
        }
        for ex in &r.examples {
            println!("  e.g. {ex}");
        }
    }
    0
}

/// Attach public-catalog external ids derivable from provenance across the workspace (or a project).
async fn backfill_catalog_ids(args: CatalogArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
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
    let r = match app.backfill_catalog_ids(project_id, args.apply).await {
        Ok(r) => r,
        Err(e) => return report(e),
    };
    if args.json {
        println!("{}", serde_json::to_string_pretty(&r).unwrap_or_default());
    } else {
        println!(
            "Examined {} subject(s); {} have derivable catalog ids.",
            r.subjects_examined, r.subjects_matched
        );
        if r.applied {
            println!("  ids added:     {}", r.ids_added);
            if r.conflicts > 0 {
                println!("  conflicts:     {} (id already owned by another subject — skipped)", r.conflicts);
            }
        } else {
            println!("  ids to add:    {} (dry run — pass --apply)", r.ids_to_add);
        }
    }
    0
}

/// Sign in via OAuth (browser + loopback callback) and persist the session for later subcommands.
async fn login(args: LoginArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    eprintln!("Opening browser to sign in as {}…", args.handle);
    match app.login(&args.handle).await {
        Ok(did) => {
            println!("Signed in: {did}");
            0
        }
        Err(e) => report(e),
    }
}

/// Delete orphaned alignment records from the signed-in account's PDS (dry-run unless `--apply`).
async fn prune_orphans(args: PruneArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    let report = match app.prune_orphan_alignments(args.apply).await {
        Ok(r) => r,
        Err(navigator_app::AppError::NotAuthenticated) => {
            eprintln!("error: not signed in — run the GUI and sign in first (the session is reused here)");
            return 1;
        }
        Err(e) => return report(e),
    };

    if args.json {
        println!("{}", serde_json::to_string_pretty(&report).unwrap_or_default());
    } else if report.applied {
        println!("Examined {} alignment record(s); deleted {} orphan(s).", report.examined, report.deleted);
        for rk in &report.orphans {
            println!("  deleted {rk}");
        }
    } else {
        println!(
            "Examined {} alignment record(s); {} orphan(s) would be deleted (dry run — pass --apply):",
            report.examined,
            report.orphans.len()
        );
        for rk in &report.orphans {
            println!("  {rk}");
        }
    }
    0
}

/// Backfill the standardized-test-label read-profile fields across the workspace (or one project).
async fn backfill_profiles(args: BackfillArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
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

    let r = match app.backfill_read_profiles(project_id, args.rescan).await {
        Ok(r) => r,
        Err(e) => return report(e),
    };

    if args.json {
        println!("{}", serde_json::to_string_pretty(&r).unwrap_or_default());
    } else {
        println!("Examined {} run(s):", r.runs_examined);
        println!("  total_bases filled:   {}", r.total_bases_filled);
        println!("  read_type inferred:   {}", r.read_type_filled);
        if args.rescan {
            println!("  read_type rescanned:  {}", r.read_type_rescanned);
        }
        println!(
            "  read_type unresolved: {}{}",
            r.read_type_unresolved,
            if !args.rescan && r.read_type_unresolved > 0 {
                "  (retry with --rescan to resolve PacBio HiFi/CLR)"
            } else {
                ""
            }
        );
    }
    0
}

/// Rebuild the genome-consensus Y **and** mtDNA signatures for a set of subjects. Used to re-place
/// existing profiles after a placement/tree-provider change (e.g. the FTDNA→DecodingUs switch): the
/// batch analyzer only *creates* a signature when one is missing, so profiles built on an older tree
/// stay stale until rebuilt here.
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
        let label = truncate(&b.donor_identifier, 24);
        // Default: only refresh subjects that already carry a Y or mt profile (the stale ones). With
        // --all, build for every subject (those with no evidence just yield an empty profile).
        if !args.all {
            let has_y = matches!(app.cached_y_profile(b.guid).await, Ok(Some(_)));
            let has_mt = matches!(app.cached_mt_profile(b.guid).await, Ok(Some(_)));
            if !has_y && !has_mt {
                skipped += 1;
                continue;
            }
        }
        // Rebuild both signatures; a source-less DNA type just produces an empty profile.
        let t = Instant::now();
        let y = app.build_y_profile(b.guid).await;
        let m = app.build_mt_profile(b.guid).await;
        match (&y, &m) {
            (Err(e), _) | (_, Err(e)) => {
                failed += 1;
                eprintln!("FAIL {label:<24} {e}");
            }
            (Ok(yp), Ok(mp)) => {
                rebuilt += 1;
                println!(
                    "OK   {label:<24} Y {:<12} mt {:<10} ({}Y/{}mt variants) [{:.1?}]",
                    yp.terminal.as_deref().unwrap_or("(none)"),
                    mp.terminal.as_deref().unwrap_or("(none)"),
                    yp.variants.len(),
                    mp.variants.len(),
                    t.elapsed()
                );
            }
        }
    }
    println!("\nrebuilt {rebuilt} subject(s), {skipped} skipped (no profile), {failed} failed");
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

    // Partition the top-level inputs. A **directory** is one staged sample: it takes the sidecar
    // fast path (Y/mt haplogroup from the GVCF + sex/read-metrics/coverage from text sidecars, no
    // CRAM decode) via `add_sample_dir`, which groups the sidecars to their alignment — something
    // per-file `add_data` can't do (and it would mis-route a `*.g.vcf.gz` through the plain-VCF
    // reader). Individual **files** keep the per-file detect-and-import path.
    let mut sample_dirs: Vec<PathBuf> = Vec::new();
    let mut files: Vec<PathBuf> = Vec::new();
    for p in &args.paths {
        if p.is_dir() {
            sample_dirs.push(p.clone());
        } else if p.is_file() {
            files.push(p.clone());
        } else {
            eprintln!("warning: skipping {} (not a file or directory)", p.display());
        }
    }
    if sample_dirs.is_empty() && files.is_empty() {
        eprintln!("error: no files or directories found in the given paths");
        return 1;
    }

    let (mut ok, mut failed, mut ysnp_panels, mut ftdna_csv) = (0usize, 0usize, 0usize, 0usize);

    // Directories: the fast-path sidecar ingest onto the resolved subject.
    for dir in &sample_dirs {
        match app.add_sample_dir(guid, dir, true).await {
            Ok(s) => {
                ok += 1;
                println!(
                    "OK   {:<18} {} ({} new, {} existing alignment(s))",
                    "sample dir",
                    dir.display(),
                    s.alignments_created,
                    s.alignments_skipped
                );
                if let Some(y) = &s.y_haplogroup {
                    println!("     Y-DNA (GVCF): {y}");
                }
                if let Some(mt) = &s.mt_haplogroup {
                    println!("     mtDNA (GVCF): {mt}");
                }
                let mut filled = Vec::new();
                if s.sex.is_some() {
                    filled.push("sex");
                }
                if s.read_metrics {
                    filled.push("read-metrics");
                }
                if s.lite_coverage {
                    filled.push("coverage");
                }
                if !filled.is_empty() {
                    println!("     sidecars: {}", filled.join(", "));
                }
                for (name, desc) in &s.imported {
                    println!("     imported {desc}: {name}");
                }
                for (name, why) in &s.skipped {
                    eprintln!("     skip {name}: {why}");
                }
                for e in &s.errors {
                    eprintln!("     warning: {e}");
                }
            }
            Err(e) => {
                failed += 1;
                eprintln!("FAIL {:<18} {}: {e}", "sample dir", dir.display());
            }
        }
    }

    // Files: per-file detect + import.
    for path in &files {
        match app.add_data_with_test_type(guid, path, args.test_type.as_deref()).await {
            Ok(detected) => {
                ok += 1;
                if detected == navigator_app::DetectedData::YSnpPanel {
                    ysnp_panels += 1;
                }
                if detected == navigator_app::DetectedData::FtdnaCsvVariants {
                    ftdna_csv += 1;
                }
                println!("OK   {:<18} {}", detected.description(), path.display());
            }
            Err(e) => {
                failed += 1;
                eprintln!("FAIL {:<18} {}: {e}", "—", path.display());
            }
        }
    }
    println!("\ningested {ok} item(s), {failed} failed, into subject \"{label}\"");

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

    // FTDNA Big Y CSV variant report(s) imported — the importer placed Y from the named (on-tree)
    // calls; report the donor's reconciled Y terminal so the admin sees it land.
    if ftdna_csv > 0 {
        match app.haplogroup_consensus(guid, DnaType::Y).await {
            Ok(Some(c)) => println!("Y-DNA (FTDNA CSV): {}", c.haplogroup),
            Ok(None) => println!("Y-DNA (FTDNA CSV): no haplogroup placed"),
            Err(e) => eprintln!("warning: reading Y consensus failed: {e}"),
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

async fn debug_place(args: ShowArgs) -> i32 {
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
    match app.debug_y_placement(guid).await {
        Ok(trace) => {
            println!("{trace}");
            0
        }
        Err(e) => {
            eprintln!("error: {e}");
            1
        }
    }
}

async fn debug_descent(args: ShowArgs) -> i32 {
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
    match app.debug_y_descent(guid).await {
        Ok(trace) => {
            println!("{trace}");
            0
        }
        Err(e) => {
            eprintln!("error: {e}");
            1
        }
    }
}

async fn debug_calls(args: DebugCallsArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    // Resolve the target alignment: explicit --alignment wins; otherwise pick from the subject
    // (prefer a CHM13/HiFi alignment, else the first).
    let alignment_id = match args.alignment {
        Some(id) => id,
        None => {
            let Some(subject) = args.subject.as_deref() else {
                eprintln!("error: provide --alignment <id> or --subject <id>");
                return 2;
            };
            let guid = match find_subject(&app, subject).await {
                Ok(Some(g)) => g,
                Ok(None) => {
                    eprintln!("error: no subject with identifier \"{subject}\"");
                    return 1;
                }
                Err(c) => return c,
            };
            match app.pick_y_debug_alignment(guid).await {
                Ok(Some(id)) => id,
                Ok(None) => {
                    eprintln!("error: subject \"{subject}\" has no alignments");
                    return 1;
                }
                Err(e) => {
                    eprintln!("error: {e}");
                    return 1;
                }
            }
        }
    };
    match app.debug_y_calls(alignment_id).await {
        Ok(trace) => {
            println!("{trace}");
            0
        }
        Err(e) => {
            eprintln!("error: {e}");
            1
        }
    }
}

async fn private_y(args: DebugCallsArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    let alignment_id = match args.alignment {
        Some(id) => id,
        None => {
            let Some(subject) = args.subject.as_deref() else {
                eprintln!("error: provide --alignment <id> or --subject <id>");
                return 2;
            };
            let guid = match find_subject(&app, subject).await {
                Ok(Some(g)) => g,
                Ok(None) => {
                    eprintln!("error: no subject with identifier \"{subject}\"");
                    return 1;
                }
                Err(c) => return c,
            };
            match app.pick_y_debug_alignment(guid).await {
                Ok(Some(id)) => id,
                Ok(None) => {
                    eprintln!("error: subject \"{subject}\" has no alignments");
                    return 1;
                }
                Err(e) => {
                    eprintln!("error: {e}");
                    return 1;
                }
            }
        }
    };
    if std::env::var("NAVIGATOR_DUMP_DENOVO").is_ok() {
        match app.run_denovo_for_alignment(alignment_id, "chrY".to_string()).await {
            Ok(calls) => {
                eprintln!("raw de-novo chrY calls: {}", calls.len());
                for c in &calls {
                    println!("DENOVO\t{}\t{}\t{}\t{}\t{:.2}", c.position, c.depth, c.alt_depth, c.alternate_allele, c.allele_fraction);
                }
                return 0;
            }
            Err(e) => {
                eprintln!("error: {e}");
                return 1;
            }
        }
    }
    let bucket = match app.private_y_variants_self_masked(alignment_id).await {
        Ok(b) => b,
        Err(e) => {
            eprintln!("error: {e}");
            return 1;
        }
    };
    let gate = app
        .publish_gate_for_alignment(alignment_id)
        .await
        .unwrap_or_default();
    println!("alignment {alignment_id} — terminal {}", bucket.terminal);
    println!("  DISPLAY (filtered) total: {}", bucket.variants.len());
    println!("    off-path known:          {}", bucket.off_path());
    println!("    novel:                   {}", bucket.novel());
    println!("    novel in structural zone:{}", bucket.in_structural_region());
    println!("    novel in unique seq:     {}", bucket.novel_in_unique_sequence());
    println!(
        "  PUBLISH (gate af>={}, alt>={}): {}",
        gate.min_allele_fraction,
        gate.min_alt_depth,
        bucket.publishable_count(gate)
    );
    if let Some(warn) = bucket.qc_banner() {
        println!("  {warn}");
    }
    // Per-variant detail (pos class region depth altDepth af publishable) for diagnosis.
    println!("  pos\tclass\tregion\tdepth\talt\taf\tpublish");
    for v in &bucket.variants {
        let class = match &v.class {
            navigator_app::PrivateClass::Novel => "novel".to_string(),
            navigator_app::PrivateClass::OffPathKnown(n) => format!("known:{n}"),
        };
        let region = v.region.map(|r| r.as_str()).unwrap_or("-");
        println!(
            "  {}\t{}\t{}\t{}\t{}\t{:.2}\t{}",
            v.position,
            class,
            region,
            v.depth,
            v.alt_depth,
            v.allele_fraction,
            gate.admits(v)
        );
    }
    0
}

/// Per-marker branch report over a Y/mtDNA node's descendant subtree — table / `--tsv` / `--json`.
async fn branch_report(args: BranchReportArgs) -> i32 {
    let app = match open(args.db).await {
        Ok(a) => a,
        Err(c) => return c,
    };
    let dna = match args.tree.to_ascii_lowercase().as_str() {
        "y" | "ydna" | "y-dna" => DnaType::Y,
        "mt" | "mtdna" | "mt-dna" => DnaType::Mt,
        other => {
            eprintln!("error: unknown --tree \"{other}\" (use \"y\" or \"mt\")");
            return 2;
        }
    };
    let alignment_id = match args.alignment {
        Some(id) => id,
        None => {
            let Some(subject) = args.subject.as_deref() else {
                eprintln!("error: provide --alignment <id> or --subject <id>");
                return 2;
            };
            let guid = match find_subject(&app, subject).await {
                Ok(Some(g)) => g,
                Ok(None) => {
                    eprintln!("error: no subject with identifier \"{subject}\"");
                    return 1;
                }
                Err(c) => return c,
            };
            // Y and mt want different alignments — a Big-Y run carries no chrM reads.
            match app.pick_alignment_for(guid, dna).await {
                Ok(Some(id)) => id,
                Ok(None) => {
                    eprintln!("error: subject \"{subject}\" has no alignments");
                    return 1;
                }
                Err(e) => {
                    eprintln!("error: {e}");
                    return 1;
                }
            }
        }
    };

    let report = match app.branch_report(alignment_id, dna, &args.node, args.depth).await {
        Ok(r) => r,
        Err(e) => {
            eprintln!("error: {e}");
            return 1;
        }
    };

    if let Some(path) = args.tsv {
        let body = navigator_app::export::branch_report_tsv(&report);
        if let Err(e) = std::fs::write(&path, body) {
            eprintln!("error: writing {}: {e}", path.display());
            return 1;
        }
        eprintln!("wrote {} ({} markers)", path.display(), report.rows.len());
        return 0;
    }

    let status = |s: navigator_app::CallState| match s {
        navigator_app::CallState::Derived => "derived",
        navigator_app::CallState::Ancestral => "ancestral",
        navigator_app::CallState::NoCall => "nocall",
    };
    let gt = |s: navigator_app::CallState| match s {
        navigator_app::CallState::Derived => "1",
        navigator_app::CallState::Ancestral => "0",
        navigator_app::CallState::NoCall => ".",
    };
    let (d, a, n) = report.counts();

    if args.json {
        let markers: Vec<_> = report
            .rows
            .iter()
            .map(|r| {
                serde_json::json!({
                    "node": r.node, "parent": r.parent, "marker": r.marker,
                    "chrom": report.contig, "pos": r.position,
                    "ancestral": r.ancestral, "derived": r.derived,
                    "observed_base": r.observed_base.map(|c| c.to_string()),
                    "status": status(r.state),
                    "ad": r.ad.map(|(rf, al)| format!("{rf},{al}")),
                    "dp": r.dp, "gq": r.gq, "source": r.source, "note": r.note,
                })
            })
            .collect();
        let out = serde_json::json!({
            "node": report.root, "contig": report.contig,
            "dna": match dna { DnaType::Y => "Y", DnaType::Mt => "mt" },
            "gvcf_backed": report.gvcf_backed,
            "derived": d, "ancestral": a, "nocall": n,
            "markers": markers,
        });
        println!("{}", serde_json::to_string_pretty(&out).unwrap_or_default());
        return 0;
    }

    println!(
        "branch report: node {} ({}, {}) — {} markers: {d} derived / {a} ancestral / {n} no-call",
        report.root,
        report.contig,
        if report.gvcf_backed { "gVCF" } else { "pileup" },
        report.rows.len(),
    );
    println!("  node\tparent\tmarker\tpos\tanc>der\tobs\tstatus\tGT\tAD\tDP\tGQ\tsource\tnote");
    for r in &report.rows {
        let opt = |o: Option<u32>| o.map(|v| v.to_string()).unwrap_or_else(|| ".".to_string());
        println!(
            "  {}\t{}\t{}\t{}\t{}>{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.node,
            r.parent,
            r.marker,
            r.position,
            r.ancestral,
            r.derived,
            r.observed_base.map(|c| c.to_string()).unwrap_or_else(|| ".".to_string()),
            status(r.state),
            gt(r.state),
            r.ad.map(|(rf, al)| format!("{rf},{al}")).unwrap_or_else(|| ".".to_string()),
            opt(r.dp),
            opt(r.gq),
            r.source,
            r.note,
        );
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
#[derive(Args)]
pub struct DoctorArgs {
    /// Alignment id to diagnose.
    #[arg(long)]
    alignment: Option<i64>,
    /// Subject donor identifier — used when `--alignment` is omitted and the subject has exactly one.
    #[arg(long, short)]
    subject: Option<String>,
    /// Diagnose a BAM/CRAM path directly, bypassing the workspace (for a file that was never imported).
    #[arg(long)]
    file: Option<PathBuf>,
    /// Reference FASTA to pair with `--file`. Required to decode a CRAM; ignored for a BAM.
    #[arg(long)]
    reference: Option<PathBuf>,
    /// Workspace database path (defaults to the GUI's ~/.decodingus/navigator-rs.db).
    #[arg(long)]
    db: Option<PathBuf>,
    /// Emit JSON instead of the human-readable report.
    #[arg(long)]
    json: bool,
}

/// Run the alignment preflight and print it. Exits 1 when a check failed, so this is usable as a
/// gate in a script and not just by eye.
///
/// `--file` deliberately skips opening the workspace: the file being undiagnosable is often *why*
/// the user cannot import it, so requiring a workspace record first would make the diagnostic
/// unavailable in the case it exists for.
async fn doctor(args: DoctorArgs) -> i32 {
    let diagnosis = if let Some(file) = args.file {
        let reference = args.reference;
        match tokio::task::spawn_blocking(move || {
            navigator_app::diagnose_alignment_file(&file, reference.as_deref())
        })
        .await
        {
            Ok(r) => r,
            Err(e) => {
                eprintln!("error: diagnosis task failed: {e}");
                return 1;
            }
        }
    } else {
        let app = match open(args.db).await {
            Ok(a) => a,
            Err(c) => return c,
        };
        let id = match resolve_alignment(&app, args.subject.as_deref(), args.alignment).await {
            Ok(id) => id,
            Err(c) => return c,
        };
        match app.diagnose_alignment(id).await {
            Ok(r) => r,
            Err(e) => return report(e),
        }
    };

    if args.json {
        println!("{}", serde_json::to_string_pretty(&diagnosis).unwrap());
    } else {
        print!("{diagnosis}");
    }
    i32::from(diagnosis.failed())
}

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
        Some(contig) => app.diploid_vcf(alignment_id, contig, navigator_app::CancelToken::none())
            .await,
        None => app.diploid_vcf_genome(alignment_id, navigator_app::CancelToken::none())
            .await,
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
