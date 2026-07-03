//! App command/query layer tests against an in-memory store.

use std::path::PathBuf;

use navigator_analysis::caller::HaploidCallerParams;
use navigator_analysis::coverage::CallableLociParams;
use navigator_app::{App, AppError};
use navigator_domain::workspace::{NewAlignment, NewProject, NewSequenceRun};
use navigator_store::Store;
use serde::{Deserialize, Serialize};

async fn app() -> App {
    App::new(Store::open_in_memory().await.unwrap())
}

/// Serializes tests that mutate the process-global `NAVIGATOR_TREE_DIR`: one test's `remove_var`
/// would otherwise yank the seeded tree dir out from under another running concurrently. Held for
/// the whole test body; ignores poisoning so a panicking test doesn't wedge the rest.
static TREE_DIR_ENV_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// Reuse the analysis crate's committed fixtures (workspace-relative).
fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../navigator-analysis/tests/fixtures")
}

/// Serializes the `NAVIGATOR_REFGENOME_DIR` env write (read once in `App::new`) so
/// parallel tests pointing the gateway cache at different temp dirs don't race.
static REF_ENV_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// An `App` whose reference-gateway cache is `cache`. The store is opened first (async), then
/// the env write + `App::new` happen synchronously under the lock so the gateway captures the
/// right base dir without racing other tests.
async fn app_with_ref_cache(cache: &std::path::Path) -> App {
    let store = Store::open_in_memory().await.unwrap();
    let _g = REF_ENV_LOCK.lock().unwrap();
    std::env::set_var("NAVIGATOR_REFGENOME_DIR", cache);
    App::new(store)
}

#[tokio::test]
async fn import_str_profile_from_csv_round_trips() {
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap(); // no project needed

    // A small FTDNA-style marker export.
    let path = std::env::temp_dir().join(format!("str-{}.csv", subject.guid.0));
    std::fs::write(&path, "Marker,Value\nDYS393,13\nDYS390,24\nDYS385,11-14\n").unwrap();

    let profile = app
        .import_str_profile_from_csv(
            subject.guid,
            "Y-37",
            Some("FTDNA".into()),
            Some("DIRECT_TEST".into()),
            &path,
        )
        .await
        .unwrap();
    assert_eq!(profile.panel_name, "Y-37");
    assert_eq!(profile.markers.len(), 3);

    let listed = app.list_str_profiles(subject.guid).await.unwrap();
    assert_eq!(listed.len(), 1);
    assert_eq!(listed[0].provider.as_deref(), Some("FTDNA"));
    assert_eq!(listed[0].markers[2].value, "11-14"); // multi-copy preserved through the store

    let _ = std::fs::remove_file(&path);
}

#[tokio::test]
async fn import_variants_from_csv_keeps_only_snps() {
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();

    let path = std::env::temp_dir().join(format!("variants-{}.csv", subject.guid.0));
    // header layout; one indel row that must be dropped (SNP-only)
    std::fs::write(
        &path,
        "contig,position,ref,alt,rsid,genotype\nchr1,1000,A,G,rs1,0/1\nchr1,2000,A,AT,rs2,0/1\nchrM,73,G,A,.,1/1\n",
    )
    .unwrap();

    let set = app
        .import_variants_from_file(subject.guid, &path, navigator_app::SourceType::Imported)
        .await
        .unwrap();
    assert_eq!(set.source_label, path.file_name().unwrap().to_string_lossy());
    assert_eq!(set.source_type, navigator_app::SourceType::Imported);
    assert_eq!(set.calls.len(), 2); // the A>AT indel dropped

    let listed = app.list_variant_sets(subject.guid).await.unwrap();
    assert_eq!(listed.len(), 1);
    assert_eq!(listed[0].calls[0].rs_id.as_deref(), Some("rs1"));
    assert_eq!(listed[0].calls[0].genotype.as_deref(), Some("0/1"));
    assert_eq!(listed[0].calls[1].rs_id, None); // "." normalized away

    let _ = std::fs::remove_file(&path);
}

#[tokio::test]
async fn y_profile_build_persists_and_reloads() {
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();

    // No Y sources yet → cached is empty until built.
    assert!(app.cached_y_profile(subject.guid).await.unwrap().is_none());

    // Build (no alignments/chip/private → an empty but valid, persisted snapshot).
    let built = app.build_y_profile(subject.guid).await.unwrap();
    assert!(built.variants.is_empty());

    // The snapshot reloads cheaply and equals the built profile (round-trips through the table).
    let cached = app
        .cached_y_profile(subject.guid)
        .await
        .unwrap()
        .expect("snapshot persisted");
    assert_eq!(cached, built);
}

#[tokio::test]
async fn import_vendor_big_y_vcf_is_tagged() {
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();

    // A Big Y bundle: a generically-named variants.vcf with the FTDNA aengine signature + a
    // sibling readme, in a per-sample directory (the parent dir disambiguates the label).
    let dir = std::env::temp_dir().join(format!("bigy-{}", subject.guid.0));
    std::fs::create_dir_all(&dir).unwrap();
    let vcf = dir.join("variants.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.1\n##reference=ucsc.hg38.fasta\n##source=aengine\n\
         ##contig=<ID=chrY,length=57227415,assembly=ucsc.hg38>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n\
         chrY\t2781339\t.\tC\tT\t.\tPASS\t.\tGT\t1\n",
    )
    .unwrap();
    std::fs::write(dir.join("readme.txt"), "…the BigY raw data files…").unwrap();

    let set = app
        .import_variants_from_file(subject.guid, &vcf, navigator_app::SourceType::Imported)
        .await
        .unwrap();
    assert!(
        set.source_label.starts_with("FTDNA Big Y"),
        "label was {}",
        set.source_label
    );
    assert!(
        set.source_label.contains(&format!("bigy-{}", subject.guid.0)),
        "label disambiguated by dir"
    );
    assert_eq!(set.source_type, navigator_app::SourceType::TargetedNgs);
    assert_eq!(set.reference_build.as_deref(), Some("GRCh38"));
    assert_eq!(set.calls.len(), 1);

    let _ = std::fs::remove_dir_all(&dir);
}

#[tokio::test]
async fn import_mtdna_fasta_derives_variants() {
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();

    // Bundled rCRS with two substitutions → an mtDNA FASTA with exactly those diffs.
    let mut seq: Vec<u8> = navigator_analysis::mtvariants::rcrs().bytes().collect();
    seq[262] = if seq[262] == b'G' { b'A' } else { b'G' }; // 1-based 263 (classic 263A>G site)
    seq[749] = if seq[749] == b'G' { b'A' } else { b'G' }; // 1-based 750
    let body = String::from_utf8(seq).unwrap();
    let path = std::env::temp_dir().join(format!("mt-{}.fasta", subject.guid.0));
    std::fs::write(&path, format!(">sample mtDNA\n{body}\n")).unwrap();

    app.import_mtdna_from_fasta(subject.guid, &path).await.unwrap();

    // The import derived + persisted an rCRS-relative variant set (haplogroup placement needs the
    // network, so it's best-effort and not asserted here).
    let sets = app.list_variant_sets(subject.guid).await.unwrap();
    let mt = sets
        .iter()
        .find(|s| s.source_label.contains("vs rCRS"))
        .expect("mtDNA variant set");
    assert_eq!(mt.source_type, navigator_app::SourceType::Sanger);
    assert!(
        mt.calls.iter().any(|c| c.position == 263),
        "expected the 263 substitution"
    );
    assert!(mt.calls.iter().all(|c| c.contig == "rCRS"));

    let _ = std::fs::remove_file(&path);
}

#[tokio::test]
async fn import_chip_profile_detects_vendor_and_summarizes() {
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();

    let path = std::env::temp_dir().join(format!("chip-{}.txt", subject.guid.0));
    std::fs::write(
        &path,
        "# This data file generated by 23andMe\nrsid\tchromosome\tposition\tgenotype\n\
         rs1\t1\t100\tAA\nrs2\t1\t200\tAG\nrs3\t1\t300\t--\nrs4\tY\t400\tG\n",
    )
    .unwrap();

    // provider=None -> auto-detect from the header
    let profile = app
        .import_chip_profile_from_csv(subject.guid, None, None, &path)
        .await
        .unwrap();
    assert_eq!(profile.provider, "23andMe");
    assert_eq!(profile.summary.total_markers_possible, 4);
    assert_eq!(profile.summary.total_markers_called, 3); // "--" no-call
    assert_eq!(profile.summary.y_markers_called, 1);
    assert_eq!(profile.summary.het_rate, Some(0.5)); // 1 het of 2 autosomal

    let listed = app.list_chip_profiles(subject.guid).await.unwrap();
    assert_eq!(listed.len(), 1);
    assert_eq!(
        listed[0].source_file_name.as_deref(),
        Some(path.file_name().unwrap().to_str().unwrap())
    );

    let _ = std::fs::remove_file(&path);
}

#[tokio::test]
async fn import_mtdna_fasta_round_trips() {
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();

    // 16,569 bp with two Ns.
    let mut body = "A".repeat(16_567);
    body.insert_str(100, "NN");
    let path = std::env::temp_dir().join(format!("mtdna-{}.fasta", subject.guid.0));
    std::fs::write(&path, format!(">sample mtDNA\n{body}\n")).unwrap();

    let seq = app.import_mtdna_from_fasta(subject.guid, &path).await.unwrap();
    assert_eq!(seq.length(), 16_569);
    assert_eq!(seq.n_count, 2);
    assert_eq!(seq.defline.as_deref(), Some("sample mtDNA"));

    let listed = app.list_mtdna_sequences(subject.guid).await.unwrap();
    assert_eq!(listed.len(), 1);
    assert_eq!(listed[0].sequence.len(), 16_569);

    // a too-short sequence is rejected
    let bad = std::env::temp_dir().join(format!("mtdna-bad-{}.fasta", subject.guid.0));
    std::fs::write(&bad, ">x\nACGT\n").unwrap();
    assert!(matches!(
        app.import_mtdna_from_fasta(subject.guid, &bad).await,
        Err(AppError::Import(_))
    ));

    let _ = std::fs::remove_file(&path);
    let _ = std::fs::remove_file(&bad);
}

#[tokio::test]
async fn derive_mtdna_variants_vs_rcrs() {
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();
    let dir = std::env::temp_dir();

    // A synthetic 16,569 bp "rCRS"; the sample differs at two positions.
    let reference = "A".repeat(16_569);
    let mut sample = reference.clone().into_bytes();
    sample[262] = b'G'; // position 263 A>G
    sample[749] = b'C'; // position 750 A>C
    let sample = String::from_utf8(sample).unwrap();

    let ref_path = dir.join(format!("rcrs-{}.fasta", subject.guid.0));
    let samp_path = dir.join(format!("mt-{}.fasta", subject.guid.0));
    std::fs::write(&ref_path, format!(">rCRS\n{reference}\n")).unwrap();
    std::fs::write(&samp_path, format!(">sample\n{sample}\n")).unwrap();

    let mt = app.import_mtdna_from_fasta(subject.guid, &samp_path).await.unwrap();
    let set = app.derive_mtdna_variants(mt.id, &ref_path).await.unwrap();
    assert_eq!(set.calls.len(), 2);
    assert_eq!(set.calls[0].contig, "rCRS");
    assert_eq!(
        (
            set.calls[0].position,
            set.calls[0].reference.as_str(),
            set.calls[0].alternate.as_str()
        ),
        (263, "A", "G")
    );
    assert_eq!((set.calls[1].position, set.calls[1].alternate.as_str()), (750, "C"));

    // Two sets now: import auto-derives one vs the bundled rCRS, plus this explicit
    // derive_mtdna_variants against the provided reference.
    assert_eq!(app.list_variant_sets(subject.guid).await.unwrap().len(), 2);

    for p in [ref_path, samp_path] {
        let _ = std::fs::remove_file(p);
    }
}

#[tokio::test]
async fn derive_mtdna_variants_detects_a_deletion() {
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();
    let dir = std::env::temp_dir();

    // rCRS: A-runs around a 9-base landmark at positions 301-309; the sample lacks it. The
    // landmark is A-free so the flanking A-runs can't absorb any of it — the 9-base
    // deletion is unambiguous.
    let landmark = "CCCCCCCCC";
    let reference = format!("{}{}{}", "A".repeat(300), landmark, "A".repeat(16_569 - 309));
    let sample = "A".repeat(16_560);
    assert_eq!(reference.len(), 16_569);
    assert_eq!(sample.len(), 16_560);

    let ref_path = dir.join(format!("rcrs-del-{}.fasta", subject.guid.0));
    let samp_path = dir.join(format!("mt-del-{}.fasta", subject.guid.0));
    std::fs::write(&ref_path, format!(">rCRS\n{reference}\n")).unwrap();
    std::fs::write(&samp_path, format!(">sample\n{sample}\n")).unwrap();

    let mt = app.import_mtdna_from_fasta(subject.guid, &samp_path).await.unwrap();
    let set = app.derive_mtdna_variants(mt.id, &ref_path).await.unwrap();
    assert_eq!(set.calls.len(), 1);
    let del = &set.calls[0];
    assert_eq!(del.position, 301);
    assert_eq!(del.reference, landmark);
    assert_eq!(del.alternate, ""); // deletion: empty alt

    for p in [ref_path, samp_path] {
        let _ = std::fs::remove_file(p);
    }
}

#[tokio::test]
async fn assign_mtdna_haplogroup_ranks_best() {
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();
    let dir = std::env::temp_dir();

    // Sample carries the tree's derived G at 263 and 750 (RSRS-anchored, no reference).
    let mut sample = "A".repeat(16_569).into_bytes();
    sample[262] = b'G'; // position 263
    sample[749] = b'G'; // position 750
    let sample = String::from_utf8(sample).unwrap();
    let samp_path = dir.join(format!("mt-hg-{}.fasta", subject.guid.0));
    std::fs::write(&samp_path, format!(">sample\n{sample}\n")).unwrap();
    let mt = app.import_mtdna_from_fasta(subject.guid, &samp_path).await.unwrap();

    let tree = r#"{"allNodes":{
        "1":{"haplogroupId":1,"name":"root","isRoot":true,"variants":[],"children":[2]},
        "2":{"haplogroupId":2,"name":"H","isRoot":false,"variants":[{"variant":"A263G","position":263,"ancestral":"A","derived":"G"}],"children":[3]},
        "3":{"haplogroupId":3,"name":"H2","isRoot":false,"variants":[{"variant":"A750G","position":750,"ancestral":"A","derived":"G"}],"children":[]}
    }}"#;

    let ranked = app.assign_mtdna_haplogroup_with_tree(mt.id, tree).await.unwrap().ranked;
    assert_eq!(ranked[0].name, "H2"); // carries both derived alleles -> deepest node
    assert_eq!(ranked[0].matched, 2);
    assert!((ranked[0].score - 1.0).abs() < 1e-9);
    assert_eq!(ranked[0].lineage, vec!["root", "H", "H2"]);

    let _ = std::fs::remove_file(&samp_path);
}

/// Real-data validation: assign mt + Y haplogroups from a GRCh38-aligned HG002 BAM
/// (chrM = rCRS, chrY = GRCh38 — matching the FTDNA trees). Needs network (live FTDNA
/// fetch). Run: `HG002_B38_BAM=/path/HG002.b38.bam cargo test -p navigator-app --test app \
/// validate_hg002 -- --ignored --nocapture`.
#[tokio::test]
#[ignore = "requires HG002_B38_BAM (GRCh38) + network"]
async fn validate_hg002_haplogroups() {
    let Ok(bam) = std::env::var("HG002_B38_BAM") else {
        eprintln!("HG002_B38_BAM unset — skipping HG002 validation");
        return;
    };
    let app = app().await;
    let b = app
        .add_biosample(None, "HG002", None, Some("male".into()))
        .await
        .unwrap();
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "GRCh38".into(),
            aligner: "bwa-mem2".into(),
            variant_caller: None,
            bam_path: Some(bam),
            reference_path: std::env::var("B38_REF").ok(), // needed for the private (de-novo) bucket
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;

    let mt = app
        .assign_mtdna_haplogroup_from_alignment(aln)
        .await
        .expect("mt assign")
        .ranked;
    let top = &mt[0];
    eprintln!(
        "HG002 mtDNA: {}  ({}/{} mutations, score {:.3})",
        top.name, top.matched, top.expected, top.score
    );
    eprintln!("  lineage: {}", top.lineage.join(" › "));
    for r in mt.iter().skip(1).take(3) {
        eprintln!("  alt: {} ({:.3})", r.name, r.score);
    }
    assert!(top.depth > 0 && top.matched > 0, "mt should resolve below root");

    let y = app.assign_y_haplogroup(aln).await.expect("Y assign");
    let top = &y.ranked[0];
    eprintln!(
        "HG002 Y: {}  ({}/{} mutations, score {:.3})",
        top.name, top.matched, top.expected, top.score
    );
    eprintln!("  lineage: {}", top.lineage.join(" › "));
    for r in y.ranked.iter().skip(1).take(3) {
        eprintln!("  alt: {} ({:.3})", r.name, r.score);
    }
    // Why descent stopped: child branches and their defining-SNP states.
    use navigator_app::CallState;
    for b in &y.branches {
        eprintln!("  child {} — {}/{} derived:", b.name, b.derived, b.snps.len());
        for s in b.snps.iter().take(8) {
            let st = match s.state {
                CallState::Derived => "DERIVED",
                CallState::Ancestral => "ancestral",
                CallState::NoCall => "no-call",
            };
            eprintln!("      {} {}{}>{}  {}", s.name, s.position, s.ancestral, s.derived, st);
        }
    }
    assert!(top.depth > 0 && top.matched > 0, "Y should resolve below root");

    // Private bucket (de-novo chrY off the backbone) — gated separately (slow + needs ref).
    if std::env::var("PRIVATE_Y").is_ok() {
        use navigator_app::PrivateClass;
        // Y_MASK_BED=path -> external mask; SELF_MASK set -> self-referential; else none.
        let bucket = if std::env::var("SELF_MASK").is_ok() {
            let ivs = app
                .callable_chr_intervals(aln, "chrY")
                .await
                .expect("callable intervals");
            let cov: i64 = ivs.iter().map(|(s, e)| e - s).sum();
            eprintln!("Self-referential callable chrY: {} intervals, {} bp", ivs.len(), cov);
            app.private_y_variants_self_masked(aln).await.expect("private Y (self)")
        } else {
            let mask = std::env::var("Y_MASK_BED").ok().map(std::path::PathBuf::from);
            app.private_y_variants(aln, mask.as_deref()).await.expect("private Y")
        };
        eprintln!(
            "Private Y below {}: {} novel, {} off-path",
            bucket.terminal,
            bucket.novel(),
            bucket.off_path()
        );
        for v in bucket
            .variants
            .iter()
            .filter(|v| matches!(v.class, PrivateClass::OffPathKnown(_)))
            .take(12)
        {
            if let PrivateClass::OffPathKnown(n) = &v.class {
                eprintln!(
                    "  off-path {} {}>{} d{} -> {}",
                    v.position, v.reference, v.alternate, v.depth, n
                );
            }
        }
        for v in bucket
            .variants
            .iter()
            .filter(|v| v.class == PrivateClass::Novel)
            .take(12)
        {
            eprintln!("  novel    {} {}>{} d{}", v.position, v.reference, v.alternate, v.depth);
        }
        assert!(!bucket.variants.is_empty(), "expected some off-backbone calls");
    }
}

/// Real-data validation that liftover gives the right answer: a CHM13-aligned HiFi BAM of a
/// donor whose GRCh38 terminals are known (Y: R-FGC29071, mt: U5a1b1g). Y is assigned by
/// lifting the GRCh38 tree positions onto CHM13 via the cached chain (auto-downloaded);
/// mtDNA is a direct chrM query (this BAM's chrM is 16,569 bp = rCRS). The calls should match
/// the GRCh38 result. Needs network (FTDNA tree + the GRCh38→CHM13 chain). Run:
///   GFX_CHM13_BAM=/Users/jkane/Genomics/GFX0457637/GFX0457637.pbmm2.chm13v2.bam \
///   cargo test -p navigator-app --test app validate_gfx_chm13 -- --ignored --nocapture
#[tokio::test]
#[ignore = "requires GFX_CHM13_BAM (CHM13) + network (FTDNA tree + liftover chain)"]
async fn validate_gfx_chm13_haplogroups() {
    let Ok(bam) = std::env::var("GFX_CHM13_BAM") else {
        eprintln!("GFX_CHM13_BAM unset — skipping CHM13 liftover validation");
        return;
    };
    let app = app().await; // default ~/.decodingus cache → the chain persists across runs
    let b = app
        .add_biosample(None, "GFX0457637", None, Some("male".into()))
        .await
        .unwrap();
    let run = app
        .record_sequence_run(NewSequenceRun {
            biosample_guid: b.guid,
            platform_name: "PACBIO".into(),
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
    // mt now needs the CHM13 reference (to self-generate the rCRS↔chrM map): resolve it
    // (downloads ~1 GB on first run, then cached), or take GFX_CHM13_REF if provided.
    let reference = match std::env::var("GFX_CHM13_REF") {
        Ok(p) => p,
        Err(_) => app
            .resolve_reference("chm13v2.0", &mut |_, _| {})
            .await
            .expect("resolve chm13v2.0 reference")
            .to_string_lossy()
            .into_owned(),
    };
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(), // triggers GRCh38→CHM13 liftover for chrY; rCRS↔chrM map for mt
            aligner: "pbmm2".into(),
            variant_caller: None,
            bam_path: Some(bam),
            reference_path: Some(reference),
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;

    // mtDNA: now lifted via the self-generated rCRS↔CHM13-chrM map — expect the U5a1b1g lineage,
    // matching the GRCh38 result.
    let mt = app
        .assign_mtdna_haplogroup_from_alignment(aln)
        .await
        .expect("mt assign")
        .ranked;
    let top = &mt[0];
    eprintln!(
        "GFX0457637 mtDNA: {}  ({}/{} mutations, score {:.3})",
        top.name, top.matched, top.expected, top.score
    );
    eprintln!("  lineage: {}", top.lineage.join(" › "));
    assert!(top.matched > 0, "mt should resolve below root");
    assert!(
        top.lineage.iter().any(|h| h.starts_with("U5")),
        "expected the U5 clade (known U5a1b1g) via the rCRS↔chrM map, got {}",
        top.lineage.join(" › ")
    );

    // Y: lifted GRCh38 tree → CHM13 chrY — expect the R-FGC29071 clade.
    let y = app.assign_y_haplogroup(aln).await.expect("Y assign");
    let top = &y.ranked[0];
    eprintln!(
        "GFX0457637 Y: {}  ({}/{} mutations, score {:.3})",
        top.name, top.matched, top.expected, top.score
    );
    eprintln!("  lineage: {}", top.lineage.join(" › "));
    for r in y.ranked.iter().skip(1).take(3) {
        eprintln!("  alt: {} ({:.3})", r.name, r.score);
    }
    assert!(top.matched > 0, "Y should resolve below root");
    assert!(
        top.lineage.iter().any(|h| h.starts_with("R-")),
        "expected the R clade (known terminal R-FGC29071), got {}",
        top.lineage.join(" › ")
    );

    // Private-Y bucket with curated CHM13 structural annotation: novel calls in palindrome /
    // amplicon / AZF-DYZ regions are paralog-prone (down-weight), the rest are unique-sequence
    // new-branch candidates. (Heavy: runs the de-novo chrY sweep; gate already requires the BAM.)
    if std::env::var("NAVIGATOR_VALIDATE_PRIVATE_Y").is_ok() {
        let bucket = app.private_y_variants_self_masked(aln).await.expect("private Y");
        use navigator_app::YRegionClass;
        let count = |c: YRegionClass| bucket.variants.iter().filter(|v| v.region == Some(c)).count();
        eprintln!(
            "private Y: {} calls — {} novel, {} off-path; structural {} (amp {}, palindrome {}, azf/dyz {}); novel-in-unique {}",
            bucket.variants.len(), bucket.novel(), bucket.off_path(), bucket.in_structural_region(),
            count(YRegionClass::Amplicon), count(YRegionClass::Palindrome), count(YRegionClass::Heterochromatin),
            bucket.novel_in_unique_sequence(),
        );
        assert!(
            bucket.in_structural_region() > 0,
            "expected some calls flagged in CHM13 Y structural regions"
        );
    }
}

/// End-to-end DecodingUs Y-tree provider against a locally-running AppView, using the CHM13
/// alignment's **native `hs1` coordinates** (no liftover). Verifies the integration places the
/// GFX sample deep onto the decoding-us backbone (the K2b clade, en route to its known
/// R-FGC29071 terminal). Reaching the R tips requires the AppView to enrich `hs1` coords for the
/// FTDNA-grafted variants (today `hs1` covers the backbone only); until then deep CHM13 placement
/// stops at the backbone. Gated on a reachable AppView. Run (AppView up on :9000, default URL):
///   GFX_CHM13_BAM=/Users/jkane/Genomics/GFX0457637/GFX0457637.pbmm2.chm13v2.bam \
///   GFX_CHM13_REF=/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa \
///   DECODINGUS_APPVIEW_URL=http://localhost:9000 \
///   cargo test -p navigator-app --test app validate_gfx_decodingus_y -- --ignored --nocapture
#[tokio::test]
#[ignore = "requires GFX_CHM13_BAM + a running DecodingUs AppView (DECODINGUS_APPVIEW_URL)"]
async fn validate_gfx_decodingus_y() {
    let (Ok(bam), Ok(reference)) = (std::env::var("GFX_CHM13_BAM"), std::env::var("GFX_CHM13_REF")) else {
        eprintln!("set GFX_CHM13_BAM + GFX_CHM13_REF (+ DECODINGUS_APPVIEW_URL) to run this");
        return;
    };
    // Force the DecodingUs provider (it's the default, but be explicit); host from env or :9000.
    std::env::set_var("NAVIGATOR_Y_TREE_PROVIDER", "decodingus");

    let app = app().await;
    let b = app
        .add_biosample(None, "GFX0457637", None, Some("male".into()))
        .await
        .unwrap();
    let run = app
        .record_sequence_run(NewSequenceRun {
            biosample_guid: b.guid,
            platform_name: "PACBIO".into(),
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(), // DecodingUs native hs1 coords → direct, no liftover
            aligner: "pbmm2".into(),
            variant_caller: None,
            bam_path: Some(bam),
            reference_path: Some(reference),
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;

    let y = app.assign_y_haplogroup(aln).await.expect("Y assign via DecodingUs");
    let top = &y.ranked[0];
    eprintln!(
        "GFX0457637 Y (DecodingUs): {}  ({}/{} mutations, score {:.3})",
        top.name, top.matched, top.expected, top.score
    );
    eprintln!("  lineage: {}", top.lineage.join(" › "));
    // Native hs1 coords place GFX deep on the decoding-us backbone (K2b, toward R-FGC29071).
    // Substantial match count + reaching the K backbone confirms the end-to-end provider works;
    // the R tips need AppView hs1 enrichment (see fn docs).
    assert!(
        top.matched >= 50,
        "expected a substantial native-hs1 match count, got {}",
        top.matched
    );
    assert!(
        top.lineage
            .iter()
            .any(|h| h == "K" || h.starts_with("K2") || h.starts_with("K-")),
        "expected to reach the K backbone via DecodingUs native hs1, got {}",
        top.lineage.join(" › ")
    );
}

/// Fast-path-only smoke test: place HG00096's Y from the precomputed GVCF (no CRAM walk),
/// against the DecodingUs tree (warm cache or AppView). Prints the terminal + lineage so we
/// can eyeball that the GVCF path produces a sensible deep placement.
///   GVCF_SMOKE_Y_GVCF (+ optional GVCF_SMOKE_M_GVCF)
#[tokio::test]
#[ignore = "requires a chrY GVCF + DecodingUs tree (AppView/cache)"]
async fn gvcf_y_placement_smoke() {
    let Ok(y_gvcf) = std::env::var("GVCF_SMOKE_Y_GVCF") else {
        eprintln!("set GVCF_SMOKE_Y_GVCF to run this");
        return;
    };
    std::env::set_var("NAVIGATOR_Y_TREE_PROVIDER", "decodingus");
    let app = app().await;
    let b = app
        .add_biosample(None, "SMOKE", None, Some("male".into()))
        .await
        .unwrap();
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa-mem".into(),
            variant_caller: None,
            bam_path: Some("/nonexistent.cram".into()), // native Y path never reads the CRAM
            reference_path: None,
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;
    let y = app
        .assign_y_from_gvcf(aln, std::path::Path::new(&y_gvcf))
        .await
        .expect("Y from GVCF");
    let top = &y.ranked[0];
    eprintln!(
        "HG00096 Y (GVCF): {} ({}/{} mutations, score {:.3})",
        top.name, top.matched, top.expected, top.score
    );
    eprintln!("  lineage: {}", top.lineage.join(" › "));
    // HG00096 is a 1000G GBR sample → deep R1b. Confirm the GVCF places it deep on the R
    // backbone (not a shallow veto), with a substantial match count.
    assert!(
        top.matched >= 100,
        "expected a deep match count from the GVCF, got {}",
        top.matched
    );
    assert!(
        top.lineage.iter().any(|h| h == "R" || h.starts_with("R1")),
        "expected HG00096 to place on the R lineage, got {}",
        top.lineage.join(" › ")
    );
}

/// **The fast-path correctness gate.** Placing a sample's Y (and mt) from the precomputed
/// pipeline GVCF must reach the same terminal as walking the CRAM — otherwise the fast path
/// is silently wrong. Set the env to a sample dir that has BOTH the CRAM and the GVCFs
/// (the ytree layout), with a running DecodingUs AppView (or a warm tree cache).
///   GVCF_PARITY_CRAM, GVCF_PARITY_REF, GVCF_PARITY_Y_GVCF[, GVCF_PARITY_M_GVCF]
#[tokio::test]
#[ignore = "requires a ytree sample dir (CRAM + GVCFs) + DecodingUs tree (AppView/cache)"]
async fn gvcf_fast_path_matches_cram_walk() {
    let (Ok(cram), Ok(reference), Ok(y_gvcf)) = (
        std::env::var("GVCF_PARITY_CRAM"),
        std::env::var("GVCF_PARITY_REF"),
        std::env::var("GVCF_PARITY_Y_GVCF"),
    ) else {
        eprintln!("set GVCF_PARITY_CRAM + GVCF_PARITY_REF + GVCF_PARITY_Y_GVCF to run this");
        return;
    };
    std::env::set_var("NAVIGATOR_Y_TREE_PROVIDER", "decodingus");

    let app = app().await;
    let b = app
        .add_biosample(None, "PARITY", None, Some("male".into()))
        .await
        .unwrap();
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa-mem".into(),
            variant_caller: None,
            bam_path: Some(cram),
            reference_path: Some(reference),
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;

    // Fast path (reads the ~MB GVCF) vs. the CRAM walk (genotypes the multi-GB CRAM).
    let fast = app
        .assign_y_from_gvcf(aln, std::path::Path::new(&y_gvcf))
        .await
        .expect("Y from GVCF");
    let slow = app.assign_y_haplogroup(aln).await.expect("Y from CRAM");
    let (ft, st) = (&fast.ranked[0], &slow.ranked[0]);
    eprintln!(
        "Y  GVCF: {} ({}/{})   CRAM: {} ({}/{})",
        ft.name, ft.matched, ft.expected, st.name, st.matched, st.expected
    );
    // Same-lineage consistency, not exact-terminal equality: the GVCF path uses robust
    // (proportional-top) selection and the CRAM path the strict guard, so they can stop at
    // different depths on the *same* path. The gate is that neither places on a different
    // branch — one lineage must contain the other's terminal.
    assert!(
        ft.lineage.contains(&st.name) || st.lineage.contains(&ft.name),
        "GVCF and CRAM placed on different Y branches: {} vs {}",
        ft.lineage.join(">"),
        st.lineage.join(">")
    );

    if let Ok(m_gvcf) = std::env::var("GVCF_PARITY_M_GVCF") {
        let fast_mt = app
            .assign_mt_from_gvcf(aln, std::path::Path::new(&m_gvcf))
            .await
            .expect("mt from GVCF");
        let slow_mt = app
            .assign_mtdna_haplogroup_from_alignment(aln)
            .await
            .expect("mt from CRAM");
        let (fm, sm) = (&fast_mt.ranked[0], &slow_mt.ranked[0]);
        eprintln!("mt GVCF: {}   CRAM: {}", fm.name, sm.name);
        assert!(
            fm.lineage.contains(&sm.name) || sm.lineage.contains(&fm.name),
            "GVCF and CRAM placed on different mt branches: {} vs {}",
            fm.lineage.join(">"),
            sm.lineage.join(">")
        );
    }
}

#[tokio::test]
async fn analysis_provenance_roundtrips_and_defaults_full_walk() {
    let app = app().await;
    let b = app.add_biosample(None, "PROV", None, None).await.unwrap();
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "x".into(),
            variant_caller: None,
            bam_path: Some("/x.cram".into()),
            reference_path: None,
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;

    // No artifact yet → no provenance.
    assert_eq!(app.analysis_provenance(aln, "coverage", "v1").await.unwrap(), None);

    // A fast-path sidecar result is partial.
    app.save_analysis_with_provenance(
        aln,
        "coverage",
        "v1",
        &serde_json::json!({"m": 1}),
        "pipeline-sidecar",
        "partial",
    )
    .await
    .unwrap();
    assert_eq!(
        app.analysis_provenance(aln, "coverage", "v1").await.unwrap(),
        Some(("pipeline-sidecar".into(), "partial".into()))
    );

    // The deep walk overwrites it → full / navigator-walk (the default save_analysis).
    app.save_analysis(aln, "coverage", "v1", &serde_json::json!({"m": 2}))
        .await
        .unwrap();
    assert_eq!(
        app.analysis_provenance(aln, "coverage", "v1").await.unwrap(),
        Some(("navigator-walk".into(), "full".into()))
    );
}

#[tokio::test]
async fn save_analysis_no_downgrade_keeps_the_fuller_result() {
    let app = app().await;
    let b = app.add_biosample(None, "NODG", None, None).await.unwrap();
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "x".into(),
            variant_caller: None,
            bam_path: Some("/x.cram".into()),
            reference_path: None,
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;

    // No artifact yet → the sidecar write goes through.
    let wrote = app
        .save_analysis_no_downgrade(aln, "coverage", "v1", &serde_json::json!({"m": 1}), "pipeline-sidecar", "partial")
        .await
        .unwrap();
    assert!(wrote, "first sidecar write with nothing present");

    // A deep walk upgrades partial → full.
    app.save_analysis(aln, "coverage", "v1", &serde_json::json!({"m": 2}))
        .await
        .unwrap();

    // Reimport: a partial sidecar must NOT clobber the full deep walk.
    let wrote = app
        .save_analysis_no_downgrade(aln, "coverage", "v1", &serde_json::json!({"m": 3}), "pipeline-sidecar", "partial")
        .await
        .unwrap();
    assert!(!wrote, "partial must not downgrade a full result");
    assert_eq!(
        app.analysis_provenance(aln, "coverage", "v1").await.unwrap(),
        Some(("navigator-walk".into(), "full".into())),
        "the full deep-walk result is preserved"
    );
    let kept: serde_json::Value = app.load_analysis(aln, "coverage", "v1").await.unwrap().unwrap();
    assert_eq!(kept["m"], 2, "payload is still the deep walk's, not the sidecar's");
}

#[tokio::test]
async fn haplogroup_consensus_combines_recorded_calls() {
    use navigator_app::{CompatibilityLevel, DnaType};
    use navigator_domain::reconciliation::RunHaplogroupCall;
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();

    // Two sources on one path: confident short-read at R-FGC29067, tentative deeper HiFi.
    for (key, label, hg, lineage, score) in [
        ("aln:1", "wgs", "R-FGC29067", vec!["root", "R", "R-FGC29067"], 0.75),
        (
            "aln:2",
            "hifi",
            "R-FGC29071",
            vec!["root", "R", "R-FGC29067", "R-FGC29071"],
            0.54,
        ),
    ] {
        let call = RunHaplogroupCall {
            source_label: label.into(),
            haplogroup: hg.into(),
            lineage: lineage.iter().map(|s| s.to_string()).collect(),
            score,
            matched: 0,
            expected: 0,
        };
        app.record_haplogroup_call(subject.guid, DnaType::Y, key, &call)
            .await
            .unwrap();
    }

    let c = app
        .haplogroup_consensus(subject.guid, DnaType::Y)
        .await
        .unwrap()
        .unwrap();
    assert_eq!(c.haplogroup, "R-FGC29067"); // confident node, not the tentative deeper one
    assert_eq!(c.compatibility, CompatibilityLevel::Compatible);
    assert_eq!(c.run_count, 2);
    assert_eq!(c.warnings.len(), 1); // flags the deeper HiFi placement

    // re-recording the same source key replaces (no duplicate)
    let calls = app.haplogroup_calls(subject.guid, DnaType::Y).await.unwrap();
    assert_eq!(calls.len(), 2);
    // mt has nothing recorded
    assert!(app
        .haplogroup_consensus(subject.guid, DnaType::Mt)
        .await
        .unwrap()
        .is_none());

    // manual override replaces the computed consensus and is flagged + audited.
    app.set_manual_override(subject.guid, DnaType::Y, "R-FGC29071", Some("Sanger-confirmed"))
        .await
        .unwrap();
    let o = app
        .haplogroup_consensus(subject.guid, DnaType::Y)
        .await
        .unwrap()
        .unwrap();
    assert_eq!(o.haplogroup, "R-FGC29071");
    assert!(o.overridden);
    assert!(o.warnings.iter().any(|w| w.contains("Sanger-confirmed")));

    app.clear_manual_override(subject.guid, DnaType::Y).await.unwrap();
    let back = app
        .haplogroup_consensus(subject.guid, DnaType::Y)
        .await
        .unwrap()
        .unwrap();
    assert_eq!(back.haplogroup, "R-FGC29067"); // back to computed
    assert!(!back.overridden);

    // audit recorded the recordings + override + clear
    let audit = app.reconciliation_audit(subject.guid, DnaType::Y).await.unwrap();
    assert!(audit.iter().any(|e| e.action == "MANUAL_OVERRIDE"));
    assert!(audit.iter().any(|e| e.action == "OVERRIDE_CLEARED"));
    assert!(audit.iter().filter(|e| e.action == "RUN_RECORDED").count() >= 2);
}

#[tokio::test]
async fn assign_haplogroup_from_alignment_calls_and_ranks() {
    // Uses the chrM coverage fixture as a stand-in alignment; the reference is ACGTACGT…,
    // so the consensus base at positions 1 and 5 is 'A' (non-variant sites).
    let app = app().await;
    let dir = fixtures();
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chrM".into(),
            aligner: "synthetic".into(),
            variant_caller: None,
            bam_path: Some(dir.join("coverage.bam").to_string_lossy().into_owned()),
            reference_path: Some(dir.join("ref.fa").to_string_lossy().into_owned()),
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;

    // root --A@1--> N1 --A@5--> N2 ; the sample carries 'A' at both -> N2 is best.
    let tree = r#"{"allNodes":{
        "1":{"haplogroupId":1,"name":"root","isRoot":true,"variants":[],"children":[2]},
        "2":{"haplogroupId":2,"name":"N1","isRoot":false,"variants":[{"variant":"x1A","position":1,"ancestral":"C","derived":"A"}],"children":[3]},
        "3":{"haplogroupId":3,"name":"N2","isRoot":false,"variants":[{"variant":"x5A","position":5,"ancestral":"C","derived":"A"}],"children":[]}
    }}"#;

    let ranked = app
        .assign_haplogroup_from_alignment(aln, "chrM", tree)
        .await
        .unwrap()
        .ranked;
    assert_eq!(ranked[0].name, "N2");
    assert_eq!(ranked[0].matched, 2);
    assert!((ranked[0].score - 1.0).abs() < 1e-9);
    assert_eq!(ranked[0].lineage, vec!["root", "N1", "N2"]);
}

#[tokio::test]
async fn add_data_detects_and_routes() {
    use navigator_app::DetectedData;
    let app = app().await;
    let subject = app.add_biosample(None, "HG002", None, None).await.unwrap();
    let dir = std::env::temp_dir();

    // An STR table -> StrProfile, stored under the subject's STR profiles.
    let str_path = dir.join(format!("data-str-{}.csv", subject.guid.0));
    std::fs::write(
        &str_path,
        "Marker,Value\nDYS393,13\nDYS390,24\nDYS19,14\nDYS391,11\nDYS385,11-14\n",
    )
    .unwrap();
    assert_eq!(
        app.add_data(subject.guid, &str_path).await.unwrap(),
        DetectedData::StrProfile
    );
    assert_eq!(app.list_str_profiles(subject.guid).await.unwrap().len(), 1);

    // A 23andMe-style export -> ChipData.
    let chip_path = dir.join(format!("data-genome-{}.txt", subject.guid.0));
    std::fs::write(
        &chip_path,
        "# 23andMe\nrsid\tchromosome\tposition\tgenotype\nrs4477212\t1\t82154\tAA\nrs3094315\t1\t752566\tAG\n",
    )
    .unwrap();
    assert_eq!(
        app.add_data(subject.guid, &chip_path).await.unwrap(),
        DetectedData::ChipData
    );
    assert_eq!(app.list_chip_profiles(subject.guid).await.unwrap().len(), 1);

    // A BAM/CRAM auto-imports: it creates a sequencing run + alignment (header probed
    // best-effort; here the bytes aren't a real BAM so detection falls back to defaults).
    let bam = dir.join(format!("data-{}.bam", subject.guid.0));
    std::fs::write(&bam, b"\x1f\x8b").unwrap();
    assert_eq!(app.add_data(subject.guid, &bam).await.unwrap(), DetectedData::Alignment);
    let runs = app.list_sequence_runs(subject.guid).await.unwrap();
    assert_eq!(runs.len(), 1);
    let alns = app.list_alignments(runs[0].id).await.unwrap();
    assert_eq!(alns.len(), 1);
    // The content hash is deferred (not computed at import) so a multi-GB alignment imports
    // instantly; it's filled in lazily on the first analysis that needs it.
    assert_eq!(alns[0].content_sha256, None, "content hash is deferred at import");
    // Idempotent: re-adding the same path doesn't duplicate the run/alignment.
    assert_eq!(app.add_data(subject.guid, &bam).await.unwrap(), DetectedData::Alignment);
    assert_eq!(app.list_sequence_runs(subject.guid).await.unwrap().len(), 1);

    for p in [str_path, chip_path, bam] {
        let _ = std::fs::remove_file(p);
    }
}

/// Create a sample → run → alignment chain and return the alignment id.
async fn alignment_id(app: &App) -> i64 {
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
    app.record_alignment(NewAlignment {
        sequence_run_id: run.id,
        reference_build: "chrM-fixture".into(),
        aligner: "synthetic".into(),
        variant_caller: None,
        bam_path: None,
        reference_path: None,
        content_sha256: None,
    })
    .await
    .unwrap()
    .id
}

#[tokio::test]
async fn command_flow_and_overview() {
    let app = app().await;
    let p = app
        .create_project(NewProject {
            name: "Trio".into(),
            description: None,
            administrator: "jk".into(),
        })
        .await
        .unwrap();

    let b1 = app
        .add_biosample(Some(p.id), "HG002", Some("SAMEA1".into()), Some("male".into()))
        .await
        .unwrap();
    app.add_biosample(Some(p.id), "HG003", None, Some("female".into()))
        .await
        .unwrap();

    let overview = app.project_overview().await.unwrap();
    assert_eq!(overview.len(), 1);
    assert_eq!(overview[0].project, p);
    assert_eq!(overview[0].sample_count, 2);

    // chain a run + alignment off the first sample
    let run = app
        .record_sequence_run(NewSequenceRun {
            biosample_guid: b1.guid,
            platform_name: "ILLUMINA".into(),
            instrument_model: None,
            test_type: "WGS".into(),
            library_layout: Some("PAIRED".into()),
            total_reads: Some(8_000_000),
            pf_reads_aligned: Some(7_956_881),
            mean_read_length: Some(148.0),
            mean_insert_size: Some(580.7),
        })
        .await
        .unwrap();
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        })
        .await
        .unwrap();
    assert_eq!(aln.sequence_run_id, run.id);
}

#[tokio::test]
async fn add_biosample_to_missing_project_is_not_found() {
    let app = app().await;
    let err = app.add_biosample(Some(123), "HG002", None, None).await;
    assert!(
        matches!(err, Err(AppError::Store(navigator_store::StoreError::NotFound(_)))),
        "got {err:?}"
    );
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
struct CoverageSummary {
    mean_coverage: f64,
    callable_bases: u64,
}

#[tokio::test]
async fn typed_analysis_artifact_round_trips_and_versions() {
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "bwa".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        })
        .await
        .unwrap();

    let summary = CoverageSummary {
        mean_coverage: 178.81,
        callable_bases: 16_292,
    };
    app.save_analysis(aln.id, "coverage", "walker-v1", &summary)
        .await
        .unwrap();

    let loaded: CoverageSummary = app
        .load_analysis(aln.id, "coverage", "walker-v1")
        .await
        .unwrap()
        .unwrap();
    assert_eq!(loaded, summary);

    // a different version is absent until written
    let v2: Option<CoverageSummary> = app.load_analysis(aln.id, "coverage", "walker-v2").await.unwrap();
    assert!(v2.is_none());
}

#[tokio::test]
async fn run_coverage_persists_and_reads_back_from_cache() {
    let app = app().await;
    let aln = alignment_id(&app).await;
    let dir = fixtures();

    assert!(app.cached_coverage(aln).await.unwrap().is_none()); // cold

    let result = app
        .run_coverage(
            aln,
            dir.join("coverage.bam"),
            dir.join("ref.fa"),
            None,
            CallableLociParams::default(),
        )
        .await
        .unwrap();
    assert_eq!(result.genome_territory, 50); // chrM fixture
    assert_eq!(result.callable_bases, 10);

    // now cached for this version (integer fields exact; floats survive round-trip to
    // ~1 ULP, so compare those approximately rather than with fragile float ==)
    let cached = app.cached_coverage(aln).await.unwrap().unwrap();
    assert_eq!(cached.genome_territory, result.genome_territory);
    assert_eq!(cached.callable_bases, result.callable_bases);
    assert_eq!(cached.contig_callable, result.contig_callable);
    assert_eq!(cached.coverage_histogram, result.coverage_histogram);
    assert!((cached.mean_coverage - result.mean_coverage).abs() < 1e-9);

    // re-running is idempotent (upsert in place; store-layer test covers row count)
    let rerun = app
        .run_coverage(
            aln,
            dir.join("coverage.bam"),
            dir.join("ref.fa"),
            None,
            CallableLociParams::default(),
        )
        .await
        .unwrap();
    assert_eq!(rerun, result);
}

#[tokio::test]
async fn run_denovo_caller_persists_snp_calls() {
    let app = app().await;
    let aln = alignment_id(&app).await;
    let dir = fixtures();

    let calls = app
        .run_denovo_caller(
            aln,
            dir.join("coverage.bam"),
            dir.join("ref.fa"),
            "chrM".into(),
            HaploidCallerParams::default(),
        )
        .await
        .unwrap();
    // fixture: ref ACGT.. with all-A reads -> SNPs where ref != A at depth>=4
    assert_eq!(
        calls.iter().map(|c| c.position).collect::<Vec<_>>(),
        vec![2, 3, 4, 6, 7, 8, 10]
    );

    let cached = app.cached_denovo(aln, "chrM").await.unwrap().unwrap();
    assert_eq!(cached, calls);
    // a different contig has no cached result
    assert!(app.cached_denovo(aln, "chrY").await.unwrap().is_none());
}

/// Build a sample → run → alignment chain whose BAM is the diploid fixture, return id.
async fn diploid_alignment(app: &App) -> i64 {
    let b = app.add_biosample(None, "diploid", None, None).await.unwrap();
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
    let bam = fixtures().join("diploid.bam").to_string_lossy().into_owned();
    app.record_alignment(NewAlignment {
        sequence_run_id: run.id,
        reference_build: "chr1".into(),
        aligner: "synthetic".into(),
        variant_caller: None,
        bam_path: Some(bam),
        reference_path: None,
        content_sha256: None,
    })
    .await
    .unwrap()
    .id
}

#[tokio::test]
async fn publish_coverage_summary_requires_cached_coverage() {
    let app = app().await;
    let aln = diploid_alignment(&app).await; // has a BAM but no coverage run
                                             // Bearer client is never reached — the missing-coverage check fails first.
    let client = navigator_app::PdsClient::bearer(reqwest::Client::new(), "http://127.0.0.1:1", "did:plc:x", "tok");
    let err = app.publish_coverage_summary(&client, aln).await;
    assert!(
        matches!(
            err,
            Err(navigator_app::AppError::Store(navigator_store::StoreError::NotFound(_)))
        ),
        "got {err:?}"
    );
}

/// Full path: run coverage on the fixture → publish the summary (a real CoverageResult,
/// floats encoded as strings) to a live PDS via a throwaway Bearer account.
#[tokio::test]
#[ignore = "requires PDS_TEST_URL (local atproto PDS container)"]
async fn publish_coverage_summary_to_live_pds() {
    let Ok(pds) = std::env::var("PDS_TEST_URL") else {
        eprintln!("PDS_TEST_URL unset — skipping live publish test");
        return;
    };
    let pds = pds.trim_end_matches('/').to_string();
    let app = app().await;
    let dir = fixtures();
    let aln = alignment_id(&app).await; // bam = coverage.bam, ref = ref.fa, build = chrM-fixture
    app.run_coverage(
        aln,
        dir.join("coverage.bam"),
        dir.join("ref.fa"),
        None,
        CallableLociParams::default(),
    )
    .await
    .expect("run_coverage");

    let client = live_bearer_client(&pds).await;
    let r = app
        .publish_coverage_summary(&client, aln)
        .await
        .expect("publish coverage summary");
    assert!(r.uri.starts_with("at://"), "uri: {}", r.uri);
    eprintln!("✓ published coverage summary {}", r.uri);
}

/// Create a throwaway PDS account and return a Bearer client for it (live tests).
async fn live_bearer_client(pds: &str) -> navigator_app::PdsClient {
    let http = reqwest::Client::new();
    let n = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos()
        % 1_000_000_000;
    let acct: serde_json::Value = http
        .post(format!("{pds}/xrpc/com.atproto.server.createAccount"))
        .json(&serde_json::json!({ "handle": format!("navp{n}.pds.test"), "email": format!("navp{n}@example.test"), "password": "navp-pw-123456" }))
        .send()
        .await
        .expect("createAccount")
        .json()
        .await
        .expect("createAccount json");
    navigator_app::PdsClient::bearer(
        http,
        pds,
        acct["did"].as_str().unwrap(),
        acct["accessJwt"].as_str().unwrap(),
    )
}

#[tokio::test]
async fn publish_without_login_is_not_authenticated() {
    let app = app().await;
    assert_eq!(app.current_account(), None);
    // Both convenience publishers refuse before sign-in (no keychain access needed).
    assert!(matches!(
        app.publish_coverage(1).await,
        Err(navigator_app::AppError::NotAuthenticated)
    ));
    assert!(matches!(
        app.publish_variants(1, "chrM").await,
        Err(navigator_app::AppError::NotAuthenticated)
    ));
}

#[tokio::test]
async fn publish_private_variants_requires_cached_calls() {
    let app = app().await;
    let aln = diploid_alignment(&app).await; // no de-novo run
    let client = navigator_app::PdsClient::bearer(reqwest::Client::new(), "http://127.0.0.1:1", "did:plc:x", "tok");
    let err = app.publish_private_variants(&client, aln, "chrM").await;
    assert!(
        matches!(
            err,
            Err(navigator_app::AppError::Store(navigator_store::StoreError::NotFound(_)))
        ),
        "got {err:?}"
    );
}

/// Full path: run de-novo on the fixture → publish the private-variants record (with
/// allele_fraction as a string) to a live PDS, read it back.
#[tokio::test]
#[ignore = "requires PDS_TEST_URL (local atproto PDS container)"]
async fn publish_private_variants_to_live_pds() {
    let Ok(pds) = std::env::var("PDS_TEST_URL") else {
        eprintln!("PDS_TEST_URL unset — skipping live private-variants publish test");
        return;
    };
    let pds = pds.trim_end_matches('/').to_string();
    let app = app().await;
    let dir = fixtures();
    let aln = alignment_id(&app).await;
    app.run_denovo_caller(
        aln,
        dir.join("coverage.bam"),
        dir.join("ref.fa"),
        "chrM".into(),
        HaploidCallerParams::default(),
    )
    .await
    .expect("run de-novo");

    let client = live_bearer_client(&pds).await;
    let r = app
        .publish_private_variants(&client, aln, "chrM")
        .await
        .expect("publish private variants");
    assert!(r.uri.starts_with("at://"), "uri: {}", r.uri);

    let got = client
        .get_record(navigator_app::PRIVATE_VARIANTS_COLLECTION, r.rkey())
        .await
        .expect("getRecord");
    assert_eq!(got["value"]["contig"], "chrM");
    assert_eq!(got["value"]["variants"].as_array().unwrap().len(), 7); // fixture de-novo
    assert!(got["value"]["variants"][0]["alleleFraction"].is_string());
    eprintln!("✓ published private variants {}", r.uri);
}

#[tokio::test]
async fn panel_genotyping_then_ibd_compare() {
    use navigator_analysis::ibd::IbdDetectorConfig;
    use navigator_domain::workspace::PanelSite;

    let app = app().await;
    let site = |pos, r: &str, a: &str| PanelSite {
        chrom: "chr1".into(),
        position: pos,
        reference_allele: r.into(),
        alternate_allele: a.into(),
        name: format!("s{pos}"),
    };
    // The four informative sites in the diploid fixture: hom-ref, het, het, hom-alt.
    let sites = vec![
        site(1, "A", "G"),
        site(2, "C", "G"),
        site(5, "A", "T"),
        site(8, "T", "A"),
    ];

    let panel = app.import_panel("test-panel", &sites).await.unwrap();
    assert_eq!(app.panel_site_count(panel.id).await.unwrap(), 4);

    let aln = diploid_alignment(&app).await;
    let genos = app.genotype_panel(aln, panel.id, 2).await.unwrap();
    let dosages: Vec<i32> = genos.iter().map(|g| g.dosage).collect();
    assert_eq!(dosages, vec![0, 1, 1, 2]);

    // cached read-back
    assert_eq!(
        app.cached_panel_genotypes(aln, panel.id, 2).await.unwrap().unwrap(),
        genos
    );

    // IBD self-compare with relaxed thresholds (only 4 sites): one fully-shared segment.
    let cfg = IbdDetectorConfig {
        min_snp_count: 3,
        window_size: 3,
        min_segment_cm: 0.0,
        ..IbdDetectorConfig::default()
    };
    let cmp = app.compare_ibd(aln, aln, panel.id, 2, cfg).await.unwrap();
    assert_eq!(cmp.segments.len(), 1);
    assert!(cmp.summary.total_shared_cm >= 0.0);

    // comparing against an un-genotyped alignment errors clearly.
    let other = diploid_alignment(&app).await;
    assert!(matches!(
        app.compare_ibd(aln, other, panel.id, 2, cfg).await,
        Err(AppError::NotGenotyped(_))
    ));
}

#[tokio::test]
async fn import_project_dir_creates_rows_is_idempotent_and_coverage_runs_on_cram() {
    let app = app().await;
    let fx = fixtures();

    // Build a temp project tree: <root>/HG00096/HG00096.chm13.cram(+.crai), reusing the
    // committed CRAM fixture (the .crai is index-by-offset, so the rename is fine).
    let root = std::env::temp_dir().join(format!("dun-import-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    let sample = root.join("HG00096");
    std::fs::create_dir_all(&sample).unwrap();
    std::fs::copy(fx.join("coverage.cram"), sample.join("HG00096.chm13.cram")).unwrap();
    std::fs::copy(fx.join("coverage.cram.crai"), sample.join("HG00096.chm13.cram.crai")).unwrap();

    let reference = fx.join("ref.fa");
    let summary = app
        .import_project_dir(&root, Some(reference.clone()), "tester".into(), false)
        .await
        .unwrap();
    assert_eq!(summary.samples_total, 1);
    assert_eq!(summary.samples_created, 1);
    assert_eq!(summary.alignments_created, 1);
    assert_eq!(summary.alignments_skipped, 0);
    assert!(summary.missing_index.is_empty());

    let bios = app.list_biosamples(summary.project.id).await.unwrap();
    assert_eq!(bios.len(), 1);
    assert_eq!(bios[0].donor_identifier, "HG00096");

    // Re-import: project/sample/alignment reused, nothing new created.
    let again = app
        .import_project_dir(&root, Some(reference), "tester".into(), false)
        .await
        .unwrap();
    assert_eq!(again.project.id, summary.project.id);
    assert_eq!(again.samples_created, 0);
    assert_eq!(again.alignments_created, 0);
    assert_eq!(again.alignments_skipped, 1);
    assert_eq!(app.project_overview().await.unwrap().len(), 1);

    // Coverage recompute works on the imported CRAM (reference_path was stamped).
    let aln = app.list_all_alignments().await.unwrap();
    assert_eq!(aln.len(), 1);
    assert_eq!(aln[0].reference_build, "chm13v2.0");
    let cov = app.run_coverage_for_alignment(aln[0].id).await.unwrap();
    assert_eq!(cov.genome_territory, 50); // the fixture chrM is 50 bp

    let _ = std::fs::remove_dir_all(&root);
}

#[tokio::test]
async fn reimport_under_different_project_name_reuses_subject() {
    // A person is one subject across projects: re-importing the same sample folder under a
    // different project name must reuse the subject (join it to the new project), not duplicate it.
    let app = app().await;
    let fx = fixtures();
    let reference = fx.join("ref.fa");

    let stage = |tag: &str| {
        let root = std::env::temp_dir().join(format!("dun-reimport-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let sample = root.join("HG00096");
        std::fs::create_dir_all(&sample).unwrap();
        std::fs::copy(fx.join("coverage.cram"), sample.join("HG00096.chm13.cram")).unwrap();
        std::fs::copy(fx.join("coverage.cram.crai"), sample.join("HG00096.chm13.cram.crai")).unwrap();
        root
    };

    let a = stage("a");
    let s1 = app.import_project_dir(&a, Some(reference.clone()), "t".into(), false).await.unwrap();
    assert_eq!(s1.samples_created, 1);

    // Same sample, different folder name → a distinct project, but the SAME person.
    let b = stage("b");
    let s2 = app.import_project_dir(&b, Some(reference), "t".into(), false).await.unwrap();
    assert_ne!(s2.project.id, s1.project.id, "a different folder name is a different project");
    assert_eq!(s2.samples_created, 0, "the subject is reused, not duplicated");

    // Exactly one subject in the workspace, and it's a roster member of BOTH projects.
    assert_eq!(app.list_all_biosamples().await.unwrap().len(), 1);
    assert_eq!(app.list_biosamples(s1.project.id).await.unwrap().len(), 1);
    assert_eq!(app.list_biosamples(s2.project.id).await.unwrap().len(), 1);

    let _ = std::fs::remove_dir_all(&a);
    let _ = std::fs::remove_dir_all(&b);
}

#[tokio::test]
async fn delete_project_detaches_members_and_keeps_subjects() {
    // A project is a grouping — deleting a non-empty one must succeed by detaching its members,
    // not refuse ("N subjects still belong to it"). The subjects themselves survive.
    let app = app().await;
    let p = app
        .create_project(NewProject { name: "P".into(), description: None, administrator: "t".into() })
        .await
        .unwrap();
    let b = app.add_biosample(Some(p.id), "S1", None, None).await.unwrap();
    app.add_biosample_to_project(b.guid, Some(p.id)).await.unwrap();

    app.delete_project(p.id).await.unwrap();

    assert!(app.project_overview().await.unwrap().iter().all(|o| o.project.id != p.id), "project removed");
    assert!(
        app.list_all_biosamples().await.unwrap().iter().any(|x| x.guid == b.guid),
        "subject survives the project deletion"
    );
}

#[tokio::test]
async fn project_report_rolls_up_coverage_and_csv_round_trips() {
    let app = app().await;
    let fx = fixtures();

    let root = std::env::temp_dir().join(format!("dun-report-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    let sample = root.join("HG00096");
    std::fs::create_dir_all(&sample).unwrap();
    std::fs::copy(fx.join("coverage.cram"), sample.join("HG00096.chm13.cram")).unwrap();
    std::fs::copy(fx.join("coverage.cram.crai"), sample.join("HG00096.chm13.cram.crai")).unwrap();

    let summary = app
        .import_project_dir(&root, Some(fx.join("ref.fa")), "tester".into(), false)
        .await
        .unwrap();
    let pid = summary.project.id;

    // Before coverage runs: report has the sample, coverage cells empty, haplogroups none.
    let before = app.project_report(pid).await.unwrap();
    assert_eq!(before.len(), 1);
    assert_eq!(before[0].alignment_count, 1);
    assert!(before[0].mean_coverage.is_none());
    assert!(before[0].y_haplogroup.is_none() && before[0].mt_haplogroup.is_none());

    let aln = app.list_all_alignments().await.unwrap();

    // A lite (sidecar) coverage is flagged `partial` in the report so the UI can badge it.
    let lite = app.run_coverage_for_alignment(aln[0].id).await.unwrap();
    app.save_analysis_with_provenance(
        aln[0].id,
        "coverage",
        navigator_analysis::coverage::COVERAGE_VERSION,
        &lite,
        "pipeline-sidecar",
        "partial",
    )
    .await
    .unwrap();
    let partial = app.project_report(pid).await.unwrap();
    assert!(partial[0].coverage_partial, "sidecar coverage shows as lite/partial");

    // Run the full coverage walk, then the report fills in and the partial flag clears.
    app.run_coverage_for_alignment(aln[0].id).await.unwrap();
    let after = app.project_report(pid).await.unwrap();
    assert!(after[0].mean_coverage.is_some());
    assert!(!after[0].coverage_partial, "a full walk upgrades the partial flag");
    assert_eq!(after[0].callable_bases, Some(10)); // fixture: 10 callable bases

    // CSV: header + one data row, sample id present.
    let csv = navigator_app::report_csv(&after);
    let lines: Vec<&str> = csv.lines().collect();
    assert_eq!(lines.len(), 2);
    assert!(lines[0].starts_with("sample_id,alignment_count,mean_coverage"));
    assert!(lines[1].starts_with("HG00096,1,"));

    let _ = std::fs::remove_dir_all(&root);
}

#[tokio::test]
async fn import_without_reference_resolves_from_cache_else_reports_needed() {
    use navigator_app::AppError;

    // Point the gateway cache at a temp dir; no network involved.
    let cache = std::env::temp_dir().join(format!("dun-refcache-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&cache);
    let app = app_with_ref_cache(&cache).await;
    let fx = fixtures();
    let root = std::env::temp_dir().join(format!("dun-import-noref-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    let sample = root.join("HG00096");
    std::fs::create_dir_all(&sample).unwrap();
    std::fs::copy(fx.join("coverage.cram"), sample.join("HG00096.chm13.cram")).unwrap();
    std::fs::copy(fx.join("coverage.cram.crai"), sample.join("HG00096.chm13.cram.crai")).unwrap();

    // Empty cache → import (no explicit reference) reports the chm13v2.0 build is needed, no writes.
    match app.import_project_dir(&root, None, "tester".into(), false).await {
        Err(AppError::ReferenceNeeded(needs)) => {
            assert_eq!(needs.len(), 1);
            assert_eq!(needs[0].build, "chm13v2.0");
            assert!(needs[0].url.ends_with("chm13v2.0.fa.gz"));
        }
        other => panic!("expected ReferenceNeeded, got {other:?}"),
    }
    assert!(
        app.project_overview().await.unwrap().is_empty(),
        "no writes on ReferenceNeeded"
    );

    // Seed the cache as if a download had completed (reuse the fixture ref as chm13v2.0).
    let refs = cache.join("references");
    std::fs::create_dir_all(&refs).unwrap();
    std::fs::copy(fx.join("ref.fa"), refs.join("chm13v2.0.fa")).unwrap();
    std::fs::copy(fx.join("ref.fa.fai"), refs.join("chm13v2.0.fa.fai")).unwrap();

    // Now import (no explicit reference) resolves from the cache and creates rows.
    let summary = app
        .import_project_dir(&root, None, "tester".into(), false)
        .await
        .unwrap();
    assert_eq!(summary.alignments_created, 1);
    let aln = app.list_all_alignments().await.unwrap();
    assert_eq!(
        aln[0].reference_path.as_deref(),
        Some(refs.join("chm13v2.0.fa").to_string_lossy().as_ref())
    );
    // Coverage runs on the CRAM using the cache-resolved reference.
    app.run_coverage_for_alignment(aln[0].id).await.unwrap();

    std::env::remove_var("NAVIGATOR_REFGENOME_DIR");
    let _ = std::fs::remove_dir_all(&root);
    let _ = std::fs::remove_dir_all(&cache);
}

#[tokio::test]
async fn assign_y_haplogroup_lifts_grch38_tree_onto_chm13_alignment() {
    // Seed a GRCh38→chm13v2.0 chain: chrY t[100,200) → chrY q[0,100). So 1-based tree
    // positions 101 & 105 lift (lift p-1, +1 back) to ychr.bam positions 1 & 5, where the
    // fixture has 'A' callable. The chain is pre-seeded so resolve_chain is a cache hit.
    let cache = std::env::temp_dir().join(format!("dun-lift-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&cache);
    let chains = cache.join("liftover");
    std::fs::create_dir_all(&chains).unwrap();
    std::fs::write(
        chains.join("GRCh38-to-chm13v2.0.chain"),
        "chain 1 chrY 1000 + 100 200 chrY 1000 + 0 100 1\n100\n",
    )
    .unwrap();

    let app = app_with_ref_cache(&cache).await;
    let dir = fixtures();
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(), // triggers GRCh38→CHM13 liftover for chrY
            aligner: "synthetic".into(),
            variant_caller: None,
            bam_path: Some(dir.join("ychr.bam").to_string_lossy().into_owned()),
            reference_path: None,
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;

    // GRCh38-coordinate tree: root --A@101--> N1 --A@105--> N2.
    let tree = r#"{"allNodes":{
        "1":{"haplogroupId":1,"name":"root","isRoot":true,"variants":[],"children":[2]},
        "2":{"haplogroupId":2,"name":"N1","isRoot":false,"variants":[{"variant":"y101","position":101,"ancestral":"C","derived":"A"}],"children":[3]},
        "3":{"haplogroupId":3,"name":"N2","isRoot":false,"variants":[{"variant":"y105","position":105,"ancestral":"C","derived":"A"}],"children":[]}
    }}"#;

    let ranked = app
        .assign_haplogroup_from_alignment(aln, "chrY", tree)
        .await
        .unwrap()
        .ranked;
    assert_eq!(
        ranked[0].name, "N2",
        "lifted GRCh38 positions found the derived alleles on CHM13 chrY"
    );
    assert_eq!(ranked[0].matched, 2);
    assert!((ranked[0].score - 1.0).abs() < 1e-9);

    let _ = std::fs::remove_dir_all(&cache);
}

#[tokio::test]
// Holds TREE_DIR_ENV_LOCK across awaits on purpose — it serializes tests that mutate the
// process-global NAVIGATOR_TREE_DIR env var; the guard must outlive the async body.
#[allow(clippy::await_holding_lock)]
async fn analyze_project_runs_coverage_and_attempts_y_per_sample() {
    let _env = TREE_DIR_ENV_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    // Seed the Y-tree cache so assign_y is offline; a root-only tree (no loci) means no query
    // targets and no chain — exercises the orchestration without network. Force the FTDNA provider
    // so the seeded tree is used (the default DecodingUs provider would reach out to the AppView).
    let trees = std::env::temp_dir().join(format!("dun-trees-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&trees);
    std::fs::create_dir_all(&trees).unwrap();
    std::fs::write(
        trees.join("ftdna-ytree.json"),
        r#"{"allNodes":{"1":{"haplogroupId":1,"name":"root","isRoot":true,"variants":[],"children":[]}}}"#,
    )
    .unwrap();
    std::env::set_var("NAVIGATOR_TREE_DIR", &trees);
    std::env::set_var("NAVIGATOR_Y_TREE_PROVIDER", "ftdna");

    let app = app().await;
    let dir = fixtures();
    let p = app
        .create_project(NewProject {
            name: "Proj".into(),
            description: None,
            administrator: "jk".into(),
        })
        .await
        .unwrap();
    let b = app.add_biosample(Some(p.id), "S1", None, None).await.unwrap();
    let run = app
        .record_sequence_run(NewSequenceRun {
            biosample_guid: b.guid,
            platform_name: "X".into(),
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
    app.record_alignment(NewAlignment {
        sequence_run_id: run.id,
        reference_build: "chm13v2.0".into(),
        aligner: "x".into(),
        variant_caller: None,
        bam_path: Some(dir.join("coverage.cram").to_string_lossy().into_owned()),
        reference_path: Some(dir.join("ref.fa").to_string_lossy().into_owned()),
        content_sha256: None,
    })
    .await
    .unwrap();

    let s = app.analyze_project(p.id).await.unwrap();
    assert_eq!(s.samples, 1);
    assert_eq!(s.coverage_done, 1, "coverage computed on the CRAM");
    // Y was attempted: recorded, or (here) errored on the chrM-only fixture lacking chrY.
    assert_eq!(s.y_done + s.errors.iter().filter(|e| e.contains("Y:")).count(), 1);

    // The report now shows coverage filled for the sample.
    let report = app.project_report(p.id).await.unwrap();
    assert!(report[0].mean_coverage.is_some());

    std::env::remove_var("NAVIGATOR_TREE_DIR");
    std::env::remove_var("NAVIGATOR_Y_TREE_PROVIDER");
    let _ = std::fs::remove_dir_all(&trees);
}

/// The AppView instrument→lab lookup (D8): a seeded `sequencer-lab-instruments.json` cache stands
/// in for the live endpoint (a fresh cache short-circuits the network). The returned lab name is
/// normalized to the local labs catalog's canonical display name when it matches; unknown labs
/// pass through; an unassociated instrument resolves to `None`.
#[tokio::test]
#[allow(clippy::await_holding_lock)] // see analyze_project_* — env-var serialization guard held across awaits
async fn lookup_lab_by_instrument_resolves_and_normalizes_from_cache() {
    let _env = TREE_DIR_ENV_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let trees = std::env::temp_dir().join(format!("dun-labs-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&trees);
    std::fs::create_dir_all(&trees).unwrap();
    std::fs::write(
        trees.join("sequencer-lab-instruments.json"),
        r#"[
            {"instrument_id":"A00182","lab_name":"FTDNA","is_d2c":true,"manufacturer":"Illumina","model_name":"NovaSeq 6000","website_url":null},
            {"instrument_id":"m84005","lab_name":"Acme Genomics","is_d2c":false}
        ]"#,
    )
    .unwrap();
    std::env::set_var("NAVIGATOR_TREE_DIR", &trees);

    let app = app().await;
    // "FTDNA" is a catalog alias → normalized to the canonical display name.
    assert_eq!(
        app.lookup_lab_by_instrument("A00182").await.as_deref(),
        Some("FamilyTreeDNA")
    );
    // An unlisted lab passes through unchanged.
    assert_eq!(
        app.lookup_lab_by_instrument("m84005").await.as_deref(),
        Some("Acme Genomics")
    );
    // No association → None (best-effort; the caller leaves the facility unset).
    assert_eq!(app.lookup_lab_by_instrument("UNASSOCIATED").await, None);

    std::env::remove_var("NAVIGATOR_TREE_DIR");
    let _ = std::fs::remove_dir_all(&trees);
}

/// A 23andMe import stores the haploid Y/MT genotype rows as a `Chip` variant set and places
/// BOTH a Y and an mtDNA haplogroup on import (best-effort), offline against seeded FTDNA trees.
/// The file declares build 38 so Y placement uses the FTDNA fallback (no DecodingUs AppView).
#[tokio::test]
#[allow(clippy::await_holding_lock)] // see analyze_project_* — env-var serialization guard held across awaits
async fn import_23andme_stores_calls_and_places_y_and_mt() {
    use navigator_app::DnaType;
    let _env = TREE_DIR_ENV_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    // Seed real (subset) FTDNA Y + mt trees so both placements run with NO network/integration —
    // small connected subtrees of the FTDNA trees, committed under tests/fixtures: Y is
    // R-L761 → R-L389 → R-P297 → R-M269 (GRCh38 coords); mt is H2a → H2a2 → H2a2a → H2a2a1 (rCRS).
    // Force the FTDNA provider so Y placement uses the seeded tree, not the (mutable) DecodingUs
    // instance — a hard assert on a terminal label can't depend on live curated data.
    let trees = std::env::temp_dir().join(format!("dun-chip-trees-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&trees);
    std::fs::create_dir_all(&trees).unwrap();
    std::fs::write(
        trees.join("ftdna-ytree.json"),
        include_str!("fixtures/ftdna-ytree.json"),
    )
    .unwrap();
    std::fs::write(
        trees.join("ftdna-mttree.json"),
        include_str!("fixtures/ftdna-mttree.json"),
    )
    .unwrap();
    std::env::set_var("NAVIGATOR_TREE_DIR", &trees);
    std::env::set_var("NAVIGATOR_Y_TREE_PROVIDER", "ftdna");

    // A synthetic 23andMe export (GRCh38): one autosomal row (ignored), the four informative Y rows
    // derived along the R-M269 lineage (L761 G, L389 G, P297 C, M269 C), the four informative MT
    // rows derived along the H2a2a1 lineage (4769/750/8860/263 → A), plus filler MT rows so the MT
    // marker count clears the "real mt panel" threshold (≥20).
    let mut file = String::from(
        "# This data file generated by 23andMe at human assembly build 38\n\
         rsid\tchromosome\tposition\tgenotype\n\
         rsauto\t1\t100\tAG\n\
         rsY1\tY\t14661990\tG\n\
         rsY2\tY\t26586954\tG\n\
         rsY3\tY\t16544628\tC\n\
         rsY4\tY\t20577481\tC\n\
         rsM1\tMT\t4769\tA\n\
         rsM2\tMT\t750\tA\n\
         rsM3\tMT\t8860\tA\n\
         rsM4\tMT\t263\tA\n",
    );
    for i in 0..16 {
        file.push_str(&format!("rsMf{i}\tMT\t{}\tA\n", 3000 + i));
    }
    let dir = std::env::temp_dir().join(format!("dun-chip-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("genome_23andMe.txt");
    std::fs::write(&path, &file).unwrap();

    let app = app().await;
    let b = app.add_biosample(None, "ARRAY1", None, None).await.unwrap();
    let profile = app
        .import_chip_profile_from_csv(b.guid, None, None, &path)
        .await
        .unwrap();
    assert_eq!(profile.provider, "23andMe");

    // The haploid Y/MT rows are stored as a Chip variant set on the vendor build (GRCh38 here).
    let sets = app.list_variant_sets(b.guid).await.unwrap();
    assert_eq!(sets.len(), 1);
    assert_eq!(sets[0].source_type, navigator_app::SourceType::Chip);
    assert_eq!(sets[0].reference_build.as_deref(), Some("GRCh38"));
    assert_eq!(
        sets[0].calls.len(),
        4 + 20,
        "4 Y + 20 MT haploid calls (autosomal row dropped)"
    );

    // Both haplogroups are placed on import, against the seeded FTDNA subset trees.
    let y = app
        .haplogroup_consensus(b.guid, DnaType::Y)
        .await
        .unwrap()
        .expect("Y placed");
    assert_eq!(y.haplogroup, "R-M269", "derived down the R-M269 lineage");
    let mt = app
        .haplogroup_consensus(b.guid, DnaType::Mt)
        .await
        .unwrap()
        .expect("mt placed");
    assert_eq!(mt.haplogroup, "H2a2a1", "derived down the H2a2a1 lineage");

    std::env::remove_var("NAVIGATOR_TREE_DIR");
    std::env::remove_var("NAVIGATOR_Y_TREE_PROVIDER");
    let _ = std::fs::remove_dir_all(&trees);
    let _ = std::fs::remove_dir_all(&dir);
}

/// Exact GRCh38-vs-CHM13 mtDNA comparison on the SAME donor (GFX0457637): the GRCh38 BAM
/// queries chrM directly (rCRS), the CHM13 BAM lifts via the self-generated rCRS↔chrM map.
/// Prints both terminals and the per-SNP lineage states, and diffs them position-by-position
/// so any difference is attributable (NoCall = coverage/unmapped vs Ancestral = different base).
///   GFX_B38_BAM=/Users/jkane/Genomics/GFX0457637/GFX0457637.b38.bam \
///   GFX_CHM13_BAM=/Users/jkane/Genomics/GFX0457637/GFX0457637.pbmm2.chm13v2.bam \
///   cargo test -p navigator-app --test app compare_mt_grch38_vs_chm13 -- --ignored --nocapture
#[tokio::test]
#[ignore = "requires GFX_B38_BAM + GFX_CHM13_BAM (+ network for the mt tree / CHM13 reference)"]
async fn compare_mt_grch38_vs_chm13() {
    let (Ok(b38), Ok(chm13)) = (std::env::var("GFX_B38_BAM"), std::env::var("GFX_CHM13_BAM")) else {
        eprintln!("GFX_B38_BAM / GFX_CHM13_BAM unset — skipping");
        return;
    };
    let app = app().await;

    async fn aln_for(app: &App, bam: String, build: &str, reference: Option<String>) -> i64 {
        let b = app
            .add_biosample(None, "GFX0457637", None, Some("male".into()))
            .await
            .unwrap();
        let run = app
            .record_sequence_run(NewSequenceRun {
                biosample_guid: b.guid,
                platform_name: "X".into(),
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
        app.record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: build.into(),
            aligner: "x".into(),
            variant_caller: None,
            bam_path: Some(bam),
            reference_path: reference,
            content_sha256: None,
        })
        .await
        .unwrap()
        .id
    }

    // GRCh38: chrM == rCRS → direct query (no reference needed for the mt path).
    let a_b38 = aln_for(&app, b38, "GRCh38", std::env::var("GFX_B38_REF").ok()).await;
    // CHM13: needs the reference to self-generate the rCRS↔chrM map.
    let chm_ref = match std::env::var("GFX_CHM13_REF") {
        Ok(p) => p,
        Err(_) => app
            .resolve_reference("chm13v2.0", &mut |_, _| {})
            .await
            .unwrap()
            .to_string_lossy()
            .into_owned(),
    };
    let a_chm = aln_for(&app, chm13, "chm13v2.0", Some(chm_ref)).await;

    let (g, g_lin, g_calls) = app.assign_mtdna_haplogroup_detail(a_b38).await.expect("b38 mt");
    let (c, c_lin, c_calls) = app.assign_mtdna_haplogroup_detail(a_chm).await.expect("chm13 mt");
    eprintln!(
        "GRCh38 mtDNA: {}  ({}/{} matched)",
        g.ranked[0].name, g.ranked[0].matched, g.ranked[0].expected
    );
    eprintln!(
        "CHM13  mtDNA: {}  ({}/{} matched)",
        c.ranked[0].name, c.ranked[0].matched, c.ranked[0].expected
    );

    // Compare the lineage element-wise (same tree + terminal → same ordered path). A position
    // can RECUR with opposite polarity (e.g. C182T then a T182C reversal), so a position-keyed
    // map would falsely diff the two occurrences against each other — match 1:1 by order.
    let _ = (&g_calls, &c_calls);
    assert_eq!(
        g_lin.len(),
        c_lin.len(),
        "lineage length differs → different terminal/path"
    );
    let mut diffs = 0;
    eprintln!("--- per-SNP lineage differences (GRCh38 vs CHM13), matched 1:1 ---");
    for (gs, cs) in g_lin.iter().zip(c_lin.iter()) {
        assert_eq!(
            (gs.position, &gs.derived),
            (cs.position, &cs.derived),
            "lineage diverged"
        );
        if gs.state != cs.state {
            diffs += 1;
            let gb = g_calls.get(&gs.position).copied().unwrap_or('-');
            let cb = c_calls.get(&cs.position).copied().unwrap_or('-');
            eprintln!(
                "  {} {}{}>{}  GRCh38={:?}(read {})  CHM13={:?}(read {})",
                gs.name, gs.position, gs.ancestral, gs.derived, gs.state, gb, cs.state, cb
            );
        }
    }
    eprintln!(
        "total TRUE differing lineage SNPs: {diffs}  (lineage length {})",
        g_lin.len()
    );
    assert_eq!(
        diffs, 0,
        "GRCh38 and CHM13 mtDNA calls should be identical on the same reads"
    );
}

// ---- sex + read metrics ----------------------------------------------------

/// Offline: run_read_metrics + run_sex persist and round-trip. Uses sex.bam (autosome + chrX;
/// sex inference requires chrX, which the chr1-only diploid fixture lacks).
#[tokio::test]
async fn sex_and_read_metrics_persist_and_reload() {
    let app = app().await;
    let b = app.add_biosample(None, "sx", None, None).await.unwrap();
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "synthetic".into(),
            variant_caller: None,
            bam_path: Some(fixtures().join("sex.bam").to_string_lossy().into_owned()),
            reference_path: None,
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;

    let m = app.run_read_metrics(aln).await.expect("run_read_metrics");
    assert!(m.total_reads > 0, "expected reads");
    assert_eq!(app.cached_read_metrics(aln).await.unwrap(), Some(m));

    let s = app.run_sex(aln).await.expect("run_sex");
    assert_eq!(app.cached_sex(aln).await.unwrap(), Some(s));
}

/// Live: GFX0457637 carries a Y haplogroup (R-FGC29071), so sex inference should call Male.
/// Uses the BAI fast-path, so it's quick. Requires GFX_CHM13_BAM.
#[tokio::test]
#[ignore = "requires GFX_CHM13_BAM"]
async fn gfx_sex_is_male() {
    let Ok(bam) = std::env::var("GFX_CHM13_BAM") else {
        eprintln!("GFX_CHM13_BAM unset — skipping");
        return;
    };
    let app = app().await;
    let b = app.add_biosample(None, "GFX0457637", None, None).await.unwrap();
    let run = app
        .record_sequence_run(NewSequenceRun {
            biosample_guid: b.guid,
            platform_name: "PACBIO_SMRT".into(),
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
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "chm13v2.0".into(),
            aligner: "pbmm2".into(),
            variant_caller: None,
            bam_path: Some(bam),
            reference_path: std::env::var("GFX_CHM13_REF").ok(),
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;
    let s = app.run_sex(aln).await.expect("run_sex");
    eprintln!(
        "sex={:?} ratio={:.3} conf={:?}",
        s.inferred_sex, s.x_autosome_ratio, s.confidence
    );
    assert_eq!(s.inferred_sex, navigator_app::InferredSex::Male);
}

/// Analysis-cache staleness: a cached artifact is reused while the source file is unchanged, and
/// invalidated (recomputed) once the file's signature changes (BAM-mtime invalidation, §6).
#[tokio::test]
async fn cached_artifact_invalidated_when_source_file_changes() {
    let app = app().await;
    let b = app
        .add_biosample(None, "HG002", None, Some("male".into()))
        .await
        .unwrap();
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

    // A stand-in "BAM" file whose mtime/size we control.
    let bam = std::env::temp_dir().join(format!("dun-cache-{}.bam", std::process::id()));
    std::fs::write(&bam, b"original").unwrap();
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "GRCh38".into(),
            aligner: "bwa-mem2".into(),
            variant_caller: None,
            bam_path: Some(bam.to_string_lossy().into_owned()),
            reference_path: None,
            content_sha256: None,
        })
        .await
        .unwrap()
        .id;

    // Save → load round-trips while the source is unchanged.
    app.save_analysis(aln, "testkind", "v1", &vec![1u32, 2, 3])
        .await
        .unwrap();
    let got: Option<Vec<u32>> = app.load_analysis(aln, "testkind", "v1").await.unwrap();
    assert_eq!(got, Some(vec![1, 2, 3]), "fresh cache is served");

    // Change the source file's content (size differs → signature differs) → cache goes stale.
    std::fs::write(&bam, b"re-aligned, different content").unwrap();
    let stale: Option<Vec<u32>> = app.load_analysis(aln, "testkind", "v1").await.unwrap();
    assert_eq!(stale, None, "changed source invalidates the cached artifact");

    // Recomputing re-stamps the new signature → served again.
    app.save_analysis(aln, "testkind", "v1", &vec![9u32]).await.unwrap();
    let fresh: Option<Vec<u32>> = app.load_analysis(aln, "testkind", "v1").await.unwrap();
    assert_eq!(fresh, Some(vec![9]), "recomputed cache is fresh again");

    let _ = std::fs::remove_file(&bam);
}

/// FTDNA project import — the B5163↔GFX merge scenario plus a new subject and an orphan
/// (ancestry without a roster row). Exercises plan → resolve fuzzy → commit end to end.
#[tokio::test]
async fn ftdna_project_import_plans_and_commits_merge_new_and_orphan() {
    use navigator_app::{DnaType, FtdnaImportOptions, FtdnaResolution, MatchKind};
    use navigator_domain::reconciliation::RunHaplogroupCall;
    use std::collections::BTreeMap;

    let ftdna = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ftdna");
    let app = app().await;
    let project = app
        .create_project(NewProject {
            name: "R1b-CTS4466Plus".into(),
            description: None,
            administrator: "admin".into(),
        })
        .await
        .unwrap();

    // The existing WGS Subject: GFX, already placed at R-FGC29071 — the same person as kit B5163.
    let gfx = app
        .add_biosample(Some(project.id), "GFX0457637", None, None)
        .await
        .unwrap();
    let call = RunHaplogroupCall {
        source_label: "wgs".into(),
        haplogroup: "R-FGC29071".into(),
        lineage: vec!["R".into(), "R-FGC29071".into()],
        score: 0.9,
        matched: 0,
        expected: 0,
    };
    app.record_haplogroup_call(gfx.guid, DnaType::Y, "aln:1", &call)
        .await
        .unwrap();

    // Dry-run plan over the roster + paternal ancestry fixtures (no maternal file).
    let plan = app
        .plan_ftdna_import(
            Some(project.id),
            None,
            Some(ftdna.join("Member_Information.csv")),
            Some(ftdna.join("Paternal_Ancestry.csv")),
            None,
            Some(ftdna.join("YDNA_Results_Overview.csv")),
            FtdnaImportOptions::default(),
        )
        .await
        .unwrap();

    // Recognized-input stats: 2 roster rows, 3 paternal, 1 Y-STR, 1 workspace subject scanned (GFX).
    assert_eq!(plan.stats.roster, 2);
    assert_eq!(plan.stats.paternal, 3);
    assert_eq!(plan.stats.ystr, 1);
    assert_eq!(plan.stats.scanned_subjects, 1);

    // Three kits: B5163 (roster+ancestry), K000002 (roster+ancestry), K000003 (ancestry only = orphan).
    assert_eq!(plan.rows.len(), 3);
    let row = |kit: &str| plan.rows.iter().find(|r| r.kit_number == kit).unwrap();

    // B5163 fuzzy-matches GFX on the Y terminal (no auto-merge: GFX had no kit#).
    let b = row("B5163");
    assert!(b.in_roster);
    assert_eq!(b.y_terminal.as_deref(), Some("FGC29071"));
    assert!(
        b.ystr_count >= 4,
        "Y-STR markers joined from the overview (got {})",
        b.ystr_count
    );
    match &b.kind {
        MatchKind::NeedsConfirm { candidates } => {
            assert!(
                candidates.iter().any(|c| c.guid == gfx.guid),
                "GFX offered as a candidate"
            );
        }
        other => panic!("expected B5163 NeedsConfirm, got {other:?}"),
    }
    // K000002: free-text Sub Group → no Y terminal → New.
    assert!(matches!(row("K000002").kind, MatchKind::New));
    assert!(row("K000002").in_roster);
    // K000003: ancestry only → New + orphan.
    assert!(matches!(row("K000003").kind, MatchKind::New));
    assert!(!row("K000003").in_roster);

    // Resolve B5163 → merge into GFX, then commit.
    let mut res = BTreeMap::new();
    res.insert("B5163".to_string(), FtdnaResolution::Merge(gfx.guid));
    let summary = app.commit_ftdna_import(&plan, &res).await.unwrap();
    assert_eq!(summary.merged, 1);
    assert_eq!(summary.created, 2);
    assert_eq!(summary.orphans, 1);
    assert_eq!(
        summary.str_profiles, 0,
        "merge adds metadata only — no duplicate Y-STR profile"
    );
    assert!(summary.errors.is_empty(), "{:?}", summary.errors);

    // Merge attaches identity/membership/MDKA but NOT a Y-STR profile (GFX keeps its own sources).
    assert!(
        app.list_str_profiles(gfx.guid).await.unwrap().is_empty(),
        "no Y-STR profile attached on merge"
    );

    // GFX gained the FTDNA kit identity, member labels, and MDKA — without a duplicate Subject.
    let ids = app.external_ids(gfx.guid).await.unwrap();
    assert_eq!(ids.len(), 1);
    assert_eq!(ids[0].external_id, "B5163");
    let member = app.ftdna_member(gfx.guid).await.unwrap().unwrap();
    assert_eq!(member.access_granted.as_deref(), Some("Limited"));
    assert_eq!(member.publicly_shares, Some(true));
    let mdkas = app.mdka_for(gfx.guid).await.unwrap();
    assert_eq!(mdkas.len(), 1);
    assert_eq!(mdkas[0].lineage, "Y");
    assert_eq!(mdkas[0].ancestor_name.as_deref(), Some("Thomas Michael Kane"));
    assert_eq!(mdkas[0].birth_year, Some(1830));
    assert_eq!(mdkas[0].origin_country.as_deref(), Some("Ireland"));
    assert_eq!(mdkas[0].latitude, Some(52.75));

    // The detail-card bundle composes all three (and a never-imported subject is empty).
    let bundle = app.subject_genealogy(gfx.guid).await.unwrap();
    assert!(!bundle.is_empty());
    assert_eq!(bundle.external_ids.len(), 1);
    assert!(bundle.member.is_some());
    assert_eq!(bundle.mdka.len(), 1);
    assert!(app
        .project_membership_ids(gfx.guid)
        .await
        .unwrap()
        .contains(&project.id));

    // The exact-kit path now auto-merges on a re-plan (the kit# is attached).
    let replan = app
        .plan_ftdna_import(
            Some(project.id),
            None,
            Some(ftdna.join("Member_Information.csv")),
            Some(ftdna.join("Paternal_Ancestry.csv")),
            None,
            None,
            FtdnaImportOptions::default(),
        )
        .await
        .unwrap();
    match &replan.rows.iter().find(|r| r.kit_number == "B5163").unwrap().kind {
        MatchKind::AutoMerge { guid, .. } => assert_eq!(*guid, gfx.guid),
        other => panic!("expected B5163 AutoMerge on re-plan, got {other:?}"),
    }
}

/// FTDNA import with no pre-selected project creates one (named from the caller) at commit — the
/// fix for the "dead Import button". A cancelled dry-run (no commit) creates nothing.
#[tokio::test]
async fn ftdna_import_into_new_project_creates_it_at_commit() {
    use navigator_app::{FtdnaImportOptions, MatchKind};
    use std::collections::BTreeMap;

    let ftdna = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ftdna");
    let app = app().await;

    // Plan into a NEW project (no project_id) — read-only, so nothing is created yet.
    let plan = app
        .plan_ftdna_import(
            None,
            Some("R1b-CTS4466Plus".into()),
            Some(ftdna.join("Member_Information.csv")),
            Some(ftdna.join("Paternal_Ancestry.csv")),
            None,
            None,
            FtdnaImportOptions::default(),
        )
        .await
        .unwrap();
    assert!(plan.project_id.is_none());
    assert_eq!(plan.project_name, "R1b-CTS4466Plus");
    // No GFX in the workspace → every kit is New.
    assert!(plan.rows.iter().all(|r| matches!(r.kind, MatchKind::New)));
    assert!(
        app.project_overview().await.unwrap().is_empty(),
        "dry-run created no project"
    );

    // Commit creates the project and imports into it.
    let summary = app.commit_ftdna_import(&plan, &BTreeMap::new()).await.unwrap();
    assert!(summary.project_id > 0);
    assert_eq!(summary.created, 3);
    let overview = app.project_overview().await.unwrap();
    assert_eq!(overview.len(), 1);
    assert_eq!(overview[0].project.name, "R1b-CTS4466Plus");
    assert_eq!(overview[0].project.id, summary.project_id);
}

/// FTDNA matching via Y-STR genetic distance: an existing subject whose Y haplogroup is an ISOGG
/// long-form label (no SNP terminal to compare) but which carries the same Y-STR profile still
/// surfaces as a fuzzy candidate (the real KANE-0001 = GFX case).
#[tokio::test]
async fn ftdna_matches_existing_subject_by_ystr_distance() {
    use navigator_app::{DnaType, FtdnaImportOptions, MatchKind};
    use navigator_domain::reconciliation::RunHaplogroupCall;

    let ftdna = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ftdna");
    let app = app().await;

    // Existing subject: ISOGG-form Y label that does NOT reduce to the FGC29071 SNP.
    let kane = app.add_biosample(None, "KANE-0001", None, None).await.unwrap();
    let call = RunHaplogroupCall {
        source_label: "wgs".into(),
        haplogroup: "R1b1a1b1a1a2c1a3a".into(),
        lineage: vec!["R".into(), "R1b1a1b1a1a2c1a3a".into()],
        score: 0.9,
        matched: 0,
        expected: 0,
    };
    app.record_haplogroup_call(kane.guid, DnaType::Y, "aln:1", &call)
        .await
        .unwrap();

    // Give KANE-0001 B5163's Y-STR markers (from the overview fixture) via a tall CSV.
    let ydna = std::fs::read_to_string(ftdna.join("YDNA_Results_Overview.csv")).unwrap();
    let per_kit = navigator_domain::ftdna::parse_ydna_overview(&ydna).unwrap();
    let (_, markers) = per_kit.iter().find(|(k, _)| k == "B5163").unwrap();
    assert!(
        markers.len() >= 67,
        "fixture must carry enough markers for an exact-haplotype match"
    );
    let mut tall = String::from("marker,value\n");
    for m in markers {
        tall.push_str(&format!("{},{}\n", m.marker, m.value));
    }
    let tmp = std::env::temp_dir().join("ftdna_ystr_match_test.csv");
    std::fs::write(&tmp, tall).unwrap();
    app.import_str_profile_from_csv(
        kane.guid,
        "FTDNA Y-700",
        Some("FTDNA".into()),
        Some("IMPORTED".into()),
        &tmp,
    )
    .await
    .unwrap();

    // Plan with roster + Y-STR. B5163's SNP terminal won't match the ISOGG label, but the Y-STR
    // genetic distance (GD 0) must surface KANE-0001 as a candidate.
    let plan = app
        .plan_ftdna_import(
            None,
            Some("R1b-CTS4466Plus".into()),
            Some(ftdna.join("Member_Information.csv")),
            None,
            None,
            Some(ftdna.join("YDNA_Results_Overview.csv")),
            FtdnaImportOptions::default(),
        )
        .await
        .unwrap();
    let b = plan.rows.iter().find(|r| r.kit_number == "B5163").unwrap();
    match &b.kind {
        MatchKind::NeedsConfirm { candidates } => {
            let c = candidates
                .iter()
                .find(|c| c.guid == kane.guid)
                .expect("KANE-0001 candidate");
            assert!(
                c.reasons.iter().any(|r| r.contains("Y-STR")),
                "matched on Y-STR: {:?}",
                c.reasons
            );
        }
        other => panic!("expected B5163 NeedsConfirm via Y-STR, got {other:?}"),
    }

    // Commit the merge into the (new) project. KANE-0001 has NO home project (`project_id` is NULL) —
    // the merge adds an M:N membership row only. The project report must still surface it (regression
    // for "matched samples don't appear in the Project report" — it reads membership ∪ home column).
    let mut res = std::collections::BTreeMap::new();
    res.insert("B5163".to_string(), navigator_app::FtdnaResolution::Merge(kane.guid));
    let summary = app.commit_ftdna_import(&plan, &res).await.unwrap();
    assert_eq!(summary.merged, 1, "{:?}", summary.errors);
    let pid = summary.project_id;

    let report = app.project_report(pid).await.unwrap();
    assert!(
        report.iter().any(|r| r.biosample.guid == kane.guid),
        "merged subject (membership-only, no home project) must appear in the project report"
    );
    // The projects-list badge counts membership-only members too.
    let overview = app.project_overview().await.unwrap();
    let ov = overview.iter().find(|o| o.project.id == pid).unwrap();
    assert!(
        ov.sample_count >= 1,
        "membership-only member counted in the project badge"
    );

    let _ = std::fs::remove_file(&tmp);
}

/// Deleting a sequencing run purges the haplogroup calls + consensus placement derived from its
/// alignments, so a wrong haplogroup doesn't linger after the run is removed.
#[tokio::test]
async fn deleting_run_purges_derived_haplogroup_and_consensus() {
    use navigator_app::DnaType;
    use navigator_domain::reconciliation::RunHaplogroupCall;

    let app = app().await;
    let b = app.add_biosample(None, "103589", None, None).await.unwrap();
    let run = app
        .record_sequence_run(NewSequenceRun {
            biosample_guid: b.guid,
            platform_name: "ILLUMINA".into(),
            instrument_model: None,
            test_type: "Targeted Y".into(),
            library_layout: None,
            total_reads: None,
            pf_reads_aligned: None,
            mean_read_length: None,
            mean_insert_size: None,
        })
        .await
        .unwrap();
    let aln = app
        .record_alignment(NewAlignment {
            sequence_run_id: run.id,
            reference_build: "GRCh38".into(),
            aligner: "unknown".into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
        })
        .await
        .unwrap();

    // The alignment-derived Y call (source_key `aln:<id>`) is what placement records.
    let call = RunHaplogroupCall {
        source_label: format!("aln #{} unknown", aln.id),
        haplogroup: "R-BY30544".into(),
        lineage: vec!["R".into(), "R-BY30544".into()],
        score: 0.72,
        matched: 0,
        expected: 0,
    };
    app.record_haplogroup_call(b.guid, DnaType::Y, &format!("aln:{}", aln.id), &call)
        .await
        .unwrap();
    assert!(
        app.haplogroup_consensus(b.guid, DnaType::Y).await.unwrap().is_some(),
        "precondition: a Y haplogroup is shown before deletion"
    );

    // Delete the run → the derived Y call + consensus must be gone.
    app.delete_sequence_run(run.id).await.unwrap();
    assert!(
        app.haplogroup_consensus(b.guid, DnaType::Y).await.unwrap().is_none(),
        "Y haplogroup must not linger after its only source run is deleted"
    );
    assert!(app.haplogroup_calls(b.guid, DnaType::Y).await.unwrap().is_empty());
}
