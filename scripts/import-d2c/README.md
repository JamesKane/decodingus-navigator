# Bulk-load the D2C sequencing collection

`load_d2c.py` loads the D2C collection described by `biosample_map.tsv` into a Navigator workspace,
driving `navigator ingest` once per sample. It is the multi-lab counterpart to
[`../import-ftdna-bigy/`](../import-ftdna-bigy/), which links a single vendor's CRAMs.

```bash
# Dry run (the default): resolve everything, write the plan, touch nothing.
./load_d2c.py

# Apply: back the workspace up, then ingest.
./load_d2c.py --apply
```

## What it does

Each manifest row is resolved along two axes — **which files** back it, and **which subject** owns
them — and then handed to `navigator ingest`.

### Subject resolution (the matching rule)

1. An existing subject carrying this lab's kit as an `external_id` (`source = FTDNA`, …) is
   **matched** and ingested into via `ingest --external-id`. The subject is never re-created.
2. Failing that — **for FTDNA rows only** — an existing subject whose `donor_identifier` *is* the
   raw kit id is matched. The members loaded from GAP are named this way. The fallback is scoped to
   FTDNA deliberately: `donor_identifier` carries no vendor namespace, so an unscoped match reads a
   bare `21418` as the same person whether it is an FGC kit or an FTDNA kit.
3. Otherwise a **new** subject is created, named `[lab]-[kit]`.

The `[lab]-[kit]` name is taken from the manifest's own `name` column rather than rebuilt from the
`lab` and `kit` fields. The column already renders that pattern, and it sanitizes the ~80 rows whose
`kit` is a UUID or a stray file path — yielding `Dante-0e74a433`, not
`Dante-184b9e17-19f0-4f4c-8ee3-9a962e36dc1b`.

A matched subject that **already has sequencing data** is left alone: the files are present, so
there is nothing to add.

### Lab precedence

When several rows resolve to the same subject, one wins and the rest are skipped as
`skip-superseded`. The order is **FTDNA > FGC > YSEQ > Dante > Nebula > other**, then richer
artifact (alignment+GVCF > VCF > STR/SNP), then subject uuid as a stable tie-break.

### Artifact resolution

The manifest records **server** paths (`/mnt/md0/Repo/…`), which are rewritten onto `--d2c-root`.
The `y_artifact` column needs more care:

| tier | artifact | ingested as |
|---|---|---|
| `gvcf_haploid`, `derived_rescue` | `<subject>/<run>/CP086569.2/` | **directory** (sidecar fast path) |
| `vcf_aengine_native` | `b38/aengine/variants.vcf.gz` | file |
| `str_only` | `b38/{DYS_Results,strs,STR}.csv` | file |
| `snp_panel` | `b38/SNPs.csv` | **skipped** — see below |
| `none` | — | skipped (no data on disk) |

The `ytree/d2c/flat/…` and `ytree/d2c/work/recover/…` Y artifacts are **server-side derivations
with no local copy**. Those rows fall back to the per-run GATK output they were derived from
(`<subject>/<run>/CP086569.2/gatk4/chrY.g.vcf.gz`) and ingest that directory instead.

Passing a **directory** is what makes this cheap: `navigator ingest` routes a directory through
`add_sample_dir`, the sidecar fast path that places Y/mt from the GVCF and fills sex, read-metrics
and coverage from the text sidecars **without decoding the CRAM**. FTDNA CRAMs are forced to
`--test-type "Big Y"` — a CRAM ships no `.bai`, so the coverage-shape detector would otherwise fall
back to `WGS`.

Because the collection lives on an external volume where ~60k individual `stat` calls cost minutes,
the tree is walked **once** and cached (`--index-cache`, default `d2c_files.txt`). Pass
`--refresh-index` after the collection changes.

## Known gaps

**`snp_panel` (11 rows) is skipped, not loaded.** FTDNA's `SNPs.csv` is a *comma*-separated Y-SNP
panel with the verdict in column 2:

```
SNP Name,Test Results,Test Type
A2578,Positive,Big Y-700
```

`looks_like_ysnp_panel` ([`navigator-domain/src/filetype.rs`](../../crates/navigator-domain/src/filetype.rs))
splits only on **tab** and only reads a verdict from **column 3**, so this never matches and `ingest`
fails with *"could not recognize the data"*. The loader skips these rows rather than create an empty
subject on a guaranteed failure. Teaching that detector the comma form would unlock them; the
conversion is not done here because the verdict column does not line up, and guessing it would risk
writing wrong genotypes.

**A handful of `variants.vcf.gz` are truncated at source** (`gzip -t` fails on them independently of
Navigator). They are reported as failures and skipped; re-copy them from the server and re-run — the
ledger will pick up only what is missing.

## Safety

- **Dry run by default.** `--apply` is required to write.
- **Backs up first.** `--apply` takes a consistent copy via SQLite's backup API (not a `cp`, which
  would tear the `.db` away from its `-wal`) and `quick_check`s it before ingesting.
- **Resumable.** Every ingested subject is appended to a ledger (`--state`). Re-running skips them,
  so an interrupted load picks up where it stopped. This also guards the one non-idempotent step:
  variant-set import is not content-idempotent, so a re-ingested VCF would *duplicate* its set.
- **Idempotent on alignments.** A CRAM already recorded as an alignment is skipped.
- **Staged rollout.** `--lab`, `--mode`, `--match` and `--limit` narrow what runs.

Every row's decision — including why it was skipped — is written to `--report`
(`d2c_load_report.tsv`), whether or not you `--apply`.

## Rollback

Restore the backup the run printed:

```bash
cp ~/.decodingus/navigator-rs.db.bak-<stamp> ~/.decodingus/navigator-rs.db
```

To unpick just this load without reverting the whole workspace, the subjects it created are the
members of the `D2C` project (`--project`, default `D2C`); the alignments it linked all live under
the D2C root:

```sql
DELETE FROM alignment WHERE bam_path LIKE '/Volumes/Genomics/D2C/%';
-- then prune the now-orphaned sequence_run rows.
```

## Options

```
--manifest PATH      biosample_map.tsv            (default: <d2c-root>/_manifests/biosample_map.tsv)
--d2c-root PATH      collection root              (default: /Volumes/Genomics/D2C)
--db PATH            workspace                    (default: ~/.decodingus/navigator-rs.db)
--navigator PATH     binary                       (default: the repo's release build)
--project NAME       file subjects under it       (default: D2C; "" to disable)
--apply              write (default: dry run)
--no-backup          skip the pre-flight backup   (NOT recommended)
--limit N            ingest at most N samples
--lab LAB            only this lab (repeatable)
--mode dir|file      only directory / single-file ingests
--match ext|donor|new  only matched / new subjects
--state PATH         resumable ledger             (default: alongside the report)
--report PATH        per-row plan                 (default: d2c_load_report.tsv)
--index-cache PATH   cached file walk             (default: d2c_files.txt)
--refresh-index      re-walk the tree
```
