# Example: link FTDNA Big Y CRAMs to workspace subjects

`link-bigy-crams.sh` is a worked example that bulk-links FTDNA Big Y CRAM exports to the biosamples
they belong to, matched by **FTDNA kit number**. Nothing in it is tied to a particular machine —
point `--root` at your own Big Y export tree (and `--db` at your workspace if it isn't in the
default location). It also serves as a reference for the `navigator ingest --external-id` workflow.

## Layout it expects

The Big Y export tree is one directory per FTDNA sample id, named by the kit number:

```
<ROOT>/<KIT>/Big_Y-700/CP086569.2/chrYM.cram     # preferred
<ROOT>/<KIT>/Big_Y-500/CP086569.2/chrYM.cram     # fallback when there is no 700
```

The directory name `<KIT>` is the FTDNA kit number. Each CRAM is aligned to `CP086569.2`
(T2T/CHM13v2 chrY + rCRS chrM). When both `Big_Y-700` and `Big_Y-500` exist only the **700** is
taken — the 500 data is merged into it.

## Matching

A directory is linked only when its name is an `external_id` of `source = FTDNA` in the workspace
DB. Everything else is **skipped** — e.g. anonymous Big Y UUID folders (which carry no kit number)
and any kit not loaded into this workspace. In a typical export only a subset of the directories
will match the subjects you actually have loaded; that is expected, and the run prints the counts.

## Idempotent

A CRAM already recorded as an alignment is skipped — checked here against the DB before invoking,
and enforced again inside `navigator ingest` (it no-ops on a duplicate `bam_path`). Safe to re-run.

## Usage

`--root` is required; the default run is a dry run.

```bash
# Dry run: report what would be linked, no writes.
./link-bigy-crams.sh --root /path/to/FTDNA

# Apply: link the matched CRAMs.
./link-bigy-crams.sh --root /path/to/FTDNA --apply

# Other options
./link-bigy-crams.sh \
  --root /path/to/FTDNA \                 # Big Y export root (required)
  --apply \
  --db   ~/.decodingus/navigator-rs.db \  # workspace DB (default shown)
  --navigator ../../target/release/navigator \
  --project "My Project"                  # also add matched subjects to this project
```

The run is recorded with `test_type = "Big Y"` (forced via `navigator ingest --test-type`, since
CRAMs ship no `.bai` for the coverage-shape detector) and `reference_build = chm13v2.0`.

## Rollback

The links are additive. To remove every alignment created from a given export root:

```sql
DELETE FROM alignment WHERE bam_path LIKE '/path/to/FTDNA/%';
-- then prune the now-orphaned sequence_run rows it created.
```

## CLI enhancements this relies on

`navigator ingest` gained (in `crates/navigator-ui/src/cli.rs` + `navigator-app`):

- `--external-id <ID>` / `--id-source <SOURCE>` (default `FTDNA`) — resolve the target subject by a
  vendor id instead of donor identifier. The subject must already exist (never created).
- `--skip-unmatched` — when the id is unknown, skip quietly (exit 0) instead of erroring.
- `--test-type <TYPE>` — force the run test type for alignment files instead of inferring it.
