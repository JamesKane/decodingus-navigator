# Sequencing lab + platform inference (read-name crowd-source)

Status: implemented (Navigator side, phases 1–4); live lab lookup blocked on AppView D8.
Scope: `navigator-domain` (SequenceRun fields, `labs` catalog), `navigator-store` (mig 0018),
`navigator-analysis` (`library_stats`), `navigator-app` (import wiring), `navigator-ui` (lab
chip + edit-run dropdown). Branch `rust-rewrite`.

## Why

The Rust rewrite dropped the Scala `SequenceRun` lab/instrument block and the
`LibraryStatsProcessor` read-name inference. Concretely it lost: the **lab**
(`sequencing_facility` — FGC/FTDNA/YSEQ/Dante/Nebula…), the **instrument serial**
(`instrument_id`, the crowd-source key), `@RG SM/LB/PU`, the flowcell, and the
**specific platform/instrument inference** (the probe only read `@RG PL/PM`, which vendor
BAMs often leave sparse). Restoring these also re-enables crowd-sourcing the lab from read
names.

## What the Scala did (ported)

`LibraryStatsProcessor` scanned ≤10k reads and per read classified the platform from the
qname, extracted the instrument + flowcell, then mapped the instrument-id prefix to a model
(Illumina `A`→NovaSeq, MGI `V300`→DNBSEQ, PacBio `m84`→Revio, Nanopore UUID, …). The lab
was resolved by `GET /api/v1/sequencer/lab?instrument_id=…` against a community-curated
instrument→lab map (`LabsConfig`/`labs.conf` was the static display catalog).

## Design

1. **Model** (`workspace::SequenceRun`) — added `sequencing_facility`, `instrument_id`,
   `sample_name`, `library_id`, `platform_unit`, `flowcell_id` (mig 0018). `NewSequenceRun`
   is unchanged; the block is written post-create by `sequence_run::set_library_stats` (33
   `NewSequenceRun` literals vs 2 read-model literals → far less churn). `update` also sets
   `sequencing_facility` (manual lab edit).
2. **Inference** (`navigator_analysis::library_stats`) — `scan_library_stats(path, ref,
   max_reads)` reads `@RG SM/LB/PU` + scans ≤10k records, `detect_platform_from_qname` /
   `parse_instrument_and_flowcell` / `infer_model` (structural matchers, no regex dep).
   Returns the most-frequent instrument/flowcell/platform + inferred model.
3. **Catalog** (`navigator_domain::labs`) — the 25 labs from `labs.conf` as static Rust data,
   with id/display/alias lookup, ≤6-char abbreviations, and `sequence_run_lab_names()` for the
   dropdown.
4. **Import + UI** — `import_alignment_file` resolves the reference first (CRAM decode), scans
   (best-effort), prefers header `@RG` for platform/model and falls back to inference, then
   persists the identity block. The run card shows a lab chip + the instrument/flowcell; the
   edit-run dialog has a lab dropdown.

## Validated

Real 43 GB BAM (`WGS229.b38.bam`) via the CLI, ~1.6 s: `platform=ILLUMINA`,
`instrument_model=NovaSeq` (inferred from `A00182`), `instrument_id=A00182`,
`sample_name=WGS229`, `library_id=WGS229_Lib1`, `platform_unit=CeGaT_NovaSeq`,
`flowcell_id=H5WLTDMXX`. Unit tests cover the qname matchers/model map + the labs catalog;
store + worker tests cover the column round-trip + the manual lab edit.

## Deferred — AppView D8 (backlog)

`sequencing_facility` is **manual** until the AppView ships the read endpoints from
`sequencer-lab-inference-system.md`: `GET /api/v1/sequencer/lab?instrument_id=…` (+ the bulk
`/lab-instruments` cache seed). Navigator already collects + stores `instrument_id`; it must
also **publish it on the `sequencerun` fed record** to feed the
`instrument_observation`→proposal→accept consensus (`fed.sequencerun.instrument_id`). Pinned
as a cross-repo contract in the AppView roadmap (`design-roadmap-rust-rewrite.md` §6, D8).
Once live, resolve the facility from `instrument_id` at import / on a backfill pass.
