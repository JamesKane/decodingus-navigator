//! Sync↔async bridge. egui runs the immediate-mode loop on the main thread; the
//! `App` (tokio + sqlx) runs on a dedicated worker thread with its own runtime. The UI
//! sends [`Command`]s and drains [`Event`]s each frame — no DB calls or domain
//! decisions on the UI thread (plan §5).
//!
//! Each command is handled on its own task so a long analysis run never blocks quick
//! queries. The command→event mapping ([`handle`]) is pure and unit-tested; [`spawn`]
//! is the thread/runtime/channel glue.

use std::path::PathBuf;
use std::sync::mpsc::{Receiver, Sender};
use std::sync::Arc;

use navigator_app::{App, Coverage, DenovoCall, IbdComparison, IbdDetectorConfig, PanelGenotype, ProjectOverview};
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::chipprofile::ChipProfile;
use navigator_domain::mtdna::MtdnaSequence;
use navigator_domain::strprofile::StrProfile;
use navigator_domain::variants::VariantSet;
use navigator_domain::workspace::{
    Alignment, Biosample, NewAlignment, NewProject, NewSequenceRun, Panel, Project, SequenceRun,
};
use tokio::sync::mpsc::{unbounded_channel, UnboundedSender};

/// Fields for adding a biosample (the app assigns its `SampleGuid`). `project_id` is
/// optional — biosamples are first-class and need not belong to a project.
#[derive(Debug, Clone)]
pub struct NewBiosample {
    pub project_id: Option<i64>,
    pub donor_identifier: String,
    pub sample_accession: Option<String>,
    pub sex: Option<String>,
}

/// A request from the UI to the worker.
#[derive(Debug, Clone)]
pub enum Command {
    LoadOverview,
    CreateProject(NewProject),
    LoadSamples(i64),
    /// Load every biosample (subjects list), regardless of project.
    LoadAllBiosamples,
    AddBiosample(NewBiosample),
    LoadRuns(SampleGuid),
    AddRun(NewSequenceRun),
    LoadStrProfiles(SampleGuid),
    ImportStrProfile {
        biosample_guid: SampleGuid,
        panel_name: String,
        provider: Option<String>,
        source: Option<String>,
        path: PathBuf,
    },
    LoadVariantSets(SampleGuid),
    ImportVariants { biosample_guid: SampleGuid, path: PathBuf },
    LoadChipProfiles(SampleGuid),
    ImportChipProfile { biosample_guid: SampleGuid, provider: Option<String>, path: PathBuf },
    LoadMtdna(SampleGuid),
    ImportMtdna { biosample_guid: SampleGuid, path: PathBuf },
    /// Unified import: detect the file's type and route it to the right importer.
    AddData { biosample_guid: SampleGuid, path: PathBuf },
    LoadAlignments(i64),
    AddAlignment(NewAlignment),
    LoadCoverage(i64),
    RunCoverage(i64),
    LoadDenovo { alignment_id: i64, contig: String },
    RunDenovo { alignment_id: i64, contig: String },
    LoadPanels,
    ImportPanel { name: String, path: PathBuf },
    LoadAllAlignments,
    GenotypePanel { alignment_id: i64, panel_id: i64, ploidy: u8 },
    LoadPanelGenotypes { alignment_id: i64, panel_id: i64, ploidy: u8 },
    CompareIbd { a: i64, b: i64, panel_id: i64, ploidy: u8 },
    /// Report who's signed in (no side effects) — sent on startup.
    AuthStatus,
    /// Report the current online/offline state (no side effects).
    SyncStatus,
    /// Sign in to a PDS via OAuth (opens a browser); `handle` is a handle or DID.
    Login { handle: String },
    Logout,
    PublishCoverage(i64),
    PublishVariants { alignment_id: i64, contig: String },
}

/// A panel with its site count, for the panel list.
#[derive(Debug, Clone)]
pub struct PanelInfo {
    pub panel: Panel,
    pub site_count: i64,
}

/// A result/notification from the worker to the UI.
#[derive(Debug, Clone)]
pub enum Event {
    Overview(Vec<ProjectOverview>),
    ProjectCreated(Project),
    Samples { project_id: i64, samples: Vec<Biosample> },
    /// All biosamples (the project-independent subjects list).
    AllBiosamples(Vec<Biosample>),
    /// A biosample was added/changed; reload the subjects list (and any open project view).
    BiosamplesChanged,
    Runs { biosample_guid: SampleGuid, runs: Vec<SequenceRun> },
    RunsChanged(SampleGuid),
    StrProfiles { biosample_guid: SampleGuid, profiles: Vec<StrProfile> },
    StrProfilesChanged(SampleGuid),
    VariantSets { biosample_guid: SampleGuid, sets: Vec<VariantSet> },
    VariantSetsChanged(SampleGuid),
    ChipProfiles { biosample_guid: SampleGuid, profiles: Vec<ChipProfile> },
    ChipProfilesChanged(SampleGuid),
    MtdnaSequences { biosample_guid: SampleGuid, sequences: Vec<MtdnaSequence> },
    MtdnaChanged(SampleGuid),
    /// A unified import succeeded; `label` describes the detected type. The UI should
    /// reload the subject's data sections.
    DataImported { biosample_guid: SampleGuid, label: String },
    Alignments { sequence_run_id: i64, alignments: Vec<Alignment> },
    AlignmentsChanged(i64),
    Coverage { alignment_id: i64, result: Option<Coverage> },
    Denovo { alignment_id: i64, contig: String, result: Option<Vec<DenovoCall>> },
    Panels(Vec<PanelInfo>),
    PanelImported,
    AllAlignments(Vec<Alignment>),
    PanelGenotypes { alignment_id: i64, panel_id: i64, ploidy: u8, genotypes: Vec<PanelGenotype> },
    Ibd(IbdComparison),
    /// Current signed-in account (DID), or `None` when signed out.
    Authenticated(Option<String>),
    /// A record was published; `kind` is a human label, `uri` the `at://` URI.
    Published { kind: String, uri: String },
    /// Whether the last PDS write reached the server (offline indicator).
    SyncOnline(bool),
    Error(String),
}

/// Execute one command against the app, mapping success/failure to an [`Event`].
pub async fn handle(app: &App, cmd: Command) -> Event {
    match cmd {
        Command::LoadOverview => match app.project_overview().await {
            Ok(v) => Event::Overview(v),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::CreateProject(new) => match app.create_project(new).await {
            Ok(p) => Event::ProjectCreated(p),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadSamples(project_id) => match app.list_biosamples(project_id).await {
            Ok(samples) => Event::Samples { project_id, samples },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAllBiosamples => match app.list_all_biosamples().await {
            Ok(samples) => Event::AllBiosamples(samples),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AddBiosample(b) => {
            match app
                .add_biosample(b.project_id, b.donor_identifier, b.sample_accession, b.sex)
                .await
            {
                Ok(_) => Event::BiosamplesChanged,
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadRuns(biosample_guid) => match app.list_sequence_runs(biosample_guid).await {
            Ok(runs) => Event::Runs { biosample_guid, runs },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AddRun(new) => match app.record_sequence_run(new).await {
            Ok(run) => Event::RunsChanged(run.biosample_guid),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadStrProfiles(guid) => match app.list_str_profiles(guid).await {
            Ok(profiles) => Event::StrProfiles { biosample_guid: guid, profiles },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportStrProfile { biosample_guid, panel_name, provider, source, path } => {
            match app
                .import_str_profile_from_csv(biosample_guid, &panel_name, provider, source, &path)
                .await
            {
                Ok(_) => Event::StrProfilesChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadVariantSets(guid) => match app.list_variant_sets(guid).await {
            Ok(sets) => Event::VariantSets { biosample_guid: guid, sets },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportVariants { biosample_guid, path } => {
            match app.import_variants_from_file(biosample_guid, &path).await {
                Ok(_) => Event::VariantSetsChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadChipProfiles(guid) => match app.list_chip_profiles(guid).await {
            Ok(profiles) => Event::ChipProfiles { biosample_guid: guid, profiles },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportChipProfile { biosample_guid, provider, path } => {
            match app.import_chip_profile_from_csv(biosample_guid, provider, None, &path).await {
                Ok(_) => Event::ChipProfilesChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadMtdna(guid) => match app.list_mtdna_sequences(guid).await {
            Ok(sequences) => Event::MtdnaSequences { biosample_guid: guid, sequences },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportMtdna { biosample_guid, path } => {
            match app.import_mtdna_from_fasta(biosample_guid, &path).await {
                Ok(_) => Event::MtdnaChanged(biosample_guid),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::AddData { biosample_guid, path } => match app.add_data(biosample_guid, &path).await {
            Ok(detected) => Event::DataImported { biosample_guid, label: detected.description().to_string() },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAlignments(sequence_run_id) => match app.list_alignments(sequence_run_id).await {
            Ok(alignments) => Event::Alignments { sequence_run_id, alignments },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::AddAlignment(new) => match app.record_alignment(new).await {
            Ok(a) => Event::AlignmentsChanged(a.sequence_run_id),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadCoverage(alignment_id) => match app.cached_coverage(alignment_id).await {
            Ok(result) => Event::Coverage { alignment_id, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunCoverage(alignment_id) => match app.run_coverage_for_alignment(alignment_id).await {
            Ok(result) => Event::Coverage { alignment_id, result: Some(result) },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadDenovo { alignment_id, contig } => match app.cached_denovo(alignment_id, &contig).await {
            Ok(result) => Event::Denovo { alignment_id, contig, result },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::RunDenovo { alignment_id, contig } => {
            match app.run_denovo_for_alignment(alignment_id, contig.clone()).await {
                Ok(result) => Event::Denovo { alignment_id, contig, result: Some(result) },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadPanels => match app.list_panels().await {
            Ok(panels) => {
                let mut infos = Vec::with_capacity(panels.len());
                for panel in panels {
                    let site_count = app.panel_site_count(panel.id).await.unwrap_or(0);
                    infos.push(PanelInfo { panel, site_count });
                }
                Event::Panels(infos)
            }
            Err(e) => Event::Error(e.to_string()),
        },
        Command::ImportPanel { name, path } => match app.import_panel_from_vcf(&name, &path).await {
            Ok(_) => Event::PanelImported,
            Err(e) => Event::Error(e.to_string()),
        },
        Command::LoadAllAlignments => match app.list_all_alignments().await {
            Ok(alns) => Event::AllAlignments(alns),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::GenotypePanel { alignment_id, panel_id, ploidy } => {
            match app.genotype_panel(alignment_id, panel_id, ploidy).await {
                Ok(genotypes) => Event::PanelGenotypes { alignment_id, panel_id, ploidy, genotypes },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::LoadPanelGenotypes { alignment_id, panel_id, ploidy } => {
            match app.cached_panel_genotypes(alignment_id, panel_id, ploidy).await {
                Ok(genotypes) => Event::PanelGenotypes {
                    alignment_id,
                    panel_id,
                    ploidy,
                    genotypes: genotypes.unwrap_or_default(),
                },
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::CompareIbd { a, b, panel_id, ploidy } => {
            match app.compare_ibd(a, b, panel_id, ploidy, IbdDetectorConfig::default()).await {
                Ok(cmp) => Event::Ibd(cmp),
                Err(e) => Event::Error(e.to_string()),
            }
        }
        Command::AuthStatus => Event::Authenticated(app.current_account()),
        Command::SyncStatus => Event::SyncOnline(app.is_online()),
        Command::Login { handle } => match app.login(&handle).await {
            Ok(did) => Event::Authenticated(Some(did)),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::Logout => match app.logout().await {
            Ok(()) => Event::Authenticated(None),
            Err(e) => Event::Error(e.to_string()),
        },
        Command::PublishCoverage(alignment_id) => match app.publish_coverage(alignment_id).await {
            Ok(r) => Event::Published { kind: "coverage summary".into(), uri: r.uri },
            Err(e) => Event::Error(e.to_string()),
        },
        Command::PublishVariants { alignment_id, contig } => {
            match app.publish_variants(alignment_id, &contig).await {
                Ok(r) => Event::Published { kind: format!("{contig} variants"), uri: r.uri },
                Err(e) => Event::Error(e.to_string()),
            }
        }
    }
}

/// Spawn the worker thread: open the workspace at `db_path` inside the worker's runtime
/// (so the connection pool lives there), then serve commands. `wake` is called after
/// each event so the UI can `request_repaint`. Returns the command sender and event
/// receiver the UI holds.
pub fn spawn(
    db_path: PathBuf,
    wake: impl Fn() + Send + Sync + 'static,
) -> (UnboundedSender<Command>, Receiver<Event>) {
    let (cmd_tx, mut cmd_rx) = unbounded_channel::<Command>();
    let (evt_tx, evt_rx) = std::sync::mpsc::channel::<Event>();
    let wake = Arc::new(wake);

    std::thread::Builder::new()
        .name("navigator-worker".into())
        .spawn(move || {
            let rt = match tokio::runtime::Runtime::new() {
                Ok(rt) => rt,
                Err(e) => {
                    let _ = evt_tx.send(Event::Error(format!("runtime: {e}")));
                    wake();
                    return;
                }
            };
            rt.block_on(async move {
                let app = match App::open(&db_path).await {
                    Ok(app) => app,
                    Err(e) => {
                        let _ = evt_tx.send(Event::Error(format!("open workspace: {e}")));
                        wake();
                        return;
                    }
                };
                while let Some(cmd) = cmd_rx.recv().await {
                    let app = app.clone();
                    let evt_tx: Sender<Event> = evt_tx.clone();
                    let wake = wake.clone();
                    tokio::spawn(async move {
                        let event = handle(&app, cmd).await;
                        let _ = evt_tx.send(event);
                        wake();
                    });
                }
            });
        })
        .expect("spawn worker thread");

    (cmd_tx, evt_rx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use navigator_store::Store;

    async fn app() -> App {
        App::new(Store::open_in_memory().await.unwrap())
    }

    #[tokio::test]
    async fn handle_maps_commands_to_events() {
        let app = app().await;

        // create a project
        let created = handle(&app, Command::CreateProject(NewProject {
            name: "Trio".into(),
            description: None,
            administrator: "jk".into(),
        }))
        .await;
        let pid = match created {
            Event::ProjectCreated(p) => p.id,
            other => panic!("expected ProjectCreated, got {other:?}"),
        };

        // overview reflects it
        match handle(&app, Command::LoadOverview).await {
            Event::Overview(v) => {
                assert_eq!(v.len(), 1);
                assert_eq!(v[0].sample_count, 0);
            }
            other => panic!("expected Overview, got {other:?}"),
        }

        // samples for the project (empty)
        match handle(&app, Command::LoadSamples(pid)).await {
            Event::Samples { project_id, samples } => {
                assert_eq!(project_id, pid);
                assert!(samples.is_empty());
            }
            other => panic!("expected Samples, got {other:?}"),
        }
    }

    #[tokio::test]
    async fn handle_runs_coverage_for_a_stored_alignment() {
        use navigator_domain::workspace::{NewAlignment, NewSequenceRun};
        use std::path::PathBuf;

        let app = app().await;
        let b = app.add_biosample(None, "HG002", None, None).await.unwrap();
        let run = app
            .record_sequence_run(NewSequenceRun {
                biosample_guid: b.guid,
                platform_name: "ILLUMINA".into(),
                instrument_model: None,
                test_type: "WGS".into(),
                library_layout: None,
                total_reads: None,
                pf_reads_aligned: None,
                mean_read_length: None,
                mean_insert_size: None,
            })
            .await
            .unwrap();
        let fixtures = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../navigator-analysis/tests/fixtures");
        let aln = app
            .record_alignment(NewAlignment {
                sequence_run_id: run.id,
                reference_build: "chrM".into(),
                aligner: "synthetic".into(),
                variant_caller: None,
                bam_path: Some(fixtures.join("coverage.bam").to_string_lossy().into_owned()),
                reference_path: Some(fixtures.join("ref.fa").to_string_lossy().into_owned()),
            })
            .await
            .unwrap();

        // alignments query (keyed by run)
        match handle(&app, Command::LoadAlignments(run.id)).await {
            Event::Alignments { alignments, .. } => assert_eq!(alignments, vec![aln.clone()]),
            other => panic!("expected Alignments, got {other:?}"),
        }

        // cold cache
        match handle(&app, Command::LoadCoverage(aln.id)).await {
            Event::Coverage { result, .. } => assert!(result.is_none()),
            other => panic!("expected Coverage(None), got {other:?}"),
        }

        // run + persist (uses the alignment's stored paths, via spawn_blocking)
        match handle(&app, Command::RunCoverage(aln.id)).await {
            Event::Coverage { alignment_id, result } => {
                assert_eq!(alignment_id, aln.id);
                assert_eq!(result.unwrap().genome_territory, 50);
            }
            other => panic!("expected Coverage(Some), got {other:?}"),
        }

        // now cached
        match handle(&app, Command::LoadCoverage(aln.id)).await {
            Event::Coverage { result, .. } => assert_eq!(result.unwrap().callable_bases, 10),
            other => panic!("expected cached Coverage, got {other:?}"),
        }

        // de-novo on the fixture contig (cold -> run -> cached), per-contig keyed
        match handle(&app, Command::LoadDenovo { alignment_id: aln.id, contig: "chrM".into() }).await {
            Event::Denovo { result, .. } => assert!(result.is_none()),
            other => panic!("expected Denovo(None), got {other:?}"),
        }
        match handle(&app, Command::RunDenovo { alignment_id: aln.id, contig: "chrM".into() }).await {
            Event::Denovo { contig, result, .. } => {
                assert_eq!(contig, "chrM");
                let calls = result.unwrap();
                assert_eq!(calls.iter().map(|c| c.position).collect::<Vec<_>>(), vec![2, 3, 4, 6, 7, 8, 10]);
            }
            other => panic!("expected Denovo(Some), got {other:?}"),
        }
        match handle(&app, Command::LoadDenovo { alignment_id: aln.id, contig: "chrM".into() }).await {
            Event::Denovo { result, .. } => assert_eq!(result.unwrap().len(), 7),
            other => panic!("expected cached Denovo, got {other:?}"),
        }
    }

    #[tokio::test]
    async fn add_commands_create_and_signal_reload() {
        use navigator_domain::workspace::{NewAlignment, NewSequenceRun};

        let app = app().await;
        let pid = match handle(&app, Command::CreateProject(NewProject {
            name: "P".into(),
            description: None,
            administrator: "jk".into(),
        }))
        .await
        {
            Event::ProjectCreated(p) => p.id,
            other => panic!("got {other:?}"),
        };

        // add a sample (tagged to the project) -> BiosamplesChanged
        match handle(&app, Command::AddBiosample(NewBiosample {
            project_id: Some(pid),
            donor_identifier: "HG002".into(),
            sample_accession: None,
            sex: Some("male".into()),
        }))
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        // it shows up both under the project and in the all-subjects list
        let guid = match handle(&app, Command::LoadSamples(pid)).await {
            Event::Samples { samples, .. } => samples[0].guid,
            other => panic!("got {other:?}"),
        };
        match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => assert_eq!(all.len(), 1),
            other => panic!("got {other:?}"),
        }

        // a project-less subject is also allowed
        match handle(&app, Command::AddBiosample(NewBiosample {
            project_id: None,
            donor_identifier: "NA12878".into(),
            sample_accession: None,
            sex: None,
        }))
        .await
        {
            Event::BiosamplesChanged => {}
            other => panic!("got {other:?}"),
        }
        match handle(&app, Command::LoadAllBiosamples).await {
            Event::AllBiosamples(all) => assert_eq!(all.len(), 2),
            other => panic!("got {other:?}"),
        }

        // add a run -> RunsChanged(sample)
        match handle(&app, Command::AddRun(NewSequenceRun {
            biosample_guid: guid,
            platform_name: "ILLUMINA".into(),
            instrument_model: None,
            test_type: "WGS".into(),
            library_layout: None,
            total_reads: None,
            pf_reads_aligned: None,
            mean_read_length: None,
            mean_insert_size: None,
        }))
        .await
        {
            Event::RunsChanged(g) => assert_eq!(g, guid),
            other => panic!("got {other:?}"),
        }
        let run_id = match handle(&app, Command::LoadRuns(guid)).await {
            Event::Runs { runs, .. } => runs[0].id,
            other => panic!("got {other:?}"),
        };

        // add an alignment -> AlignmentsChanged(run)
        match handle(&app, Command::AddAlignment(NewAlignment {
            sequence_run_id: run_id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
        }))
        .await
        {
            Event::AlignmentsChanged(r) => assert_eq!(r, run_id),
            other => panic!("got {other:?}"),
        }
    }
}
