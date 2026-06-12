# Sequencing lab + platform inference (read-name crowd-source)

Status: implemented (Navigator side, phases 1–4 + live AppView lab lookup). AppView D8 endpoints
shipped (decodingus 9c28b6d); resolution wired Navigator-side (commit c0dfae7).
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

## Live lab lookup — DONE (AppView D8 shipped)

The AppView shipped the endpoints (`decodingus` 9c28b6d): `GET /api/v1/sequencer/lab?instrument_id=…`
→ `SequencerLabDto {instrument_id, lab_name, is_d2c, manufacturer, model_name, website_url}`, and the
bulk `GET /api/v1/sequencer/lab-instruments` (join `genomics.sequencer_instrument.lab_id →
sequencing_lab` — the redesign re-added a preseeded `lab_id` FK). Navigator wiring (commit c0dfae7):

- `App::fetch_lab_instruments` GETs the bulk list, on-disk cached via the same `fetch_tree` path
  (7-day TTL + offline fallback) → one network call per batch.
- `App::lookup_lab_by_instrument` resolves + normalizes the name to the local catalog's canonical
  display name (unlisted labs pass through).
- `import_alignment_file` auto-resolves the lab inline after inferring `instrument_id` (best-effort:
  an unreachable AppView leaves it unset — no error, no blocking).
- `App::backfill_run_labs` (worker `Command::BackfillLabs`, run on startup) fills any run with an
  `instrument_id` but no facility, so runs imported before D8 landed pick up associations.

**Still open — the contribution side:** Navigator does not yet **publish `instrument_id` on the
`sequencerun` fed record**, which is the AppView's `instrument_observation`→proposal→accept consensus
source (`fed.sequencerun.instrument_id`). Today the AppView map is preseeded/curator-driven; closing
the crowd-source loop means emitting that field when the sequencerun record is published. Pinned in
the AppView roadmap (`design-roadmap-rust-rewrite.md` §6, D8).
