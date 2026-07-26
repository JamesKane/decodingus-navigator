# Academic / Public-Dataset Import (ENA) — Design

**Status:** Design only (no implementation). Drafted 2026-06-06.
**Companion to:** `ftdna-project-import.md` — this doc is the **IRB- and
publishing-facing** view of the *Academic / ENA import profile* (that doc's §9).
The shared engine (vendor-neutral Subject, `external_id`, M:N projects, STR/SNP/
haplogroup stores, matching, assertions) is described there and **not** repeated
here; this doc covers only what makes the academic path ethically and editorially
defensible.

## 1. Purpose & audience

Enable researchers to assemble, analyze, and collaborate over **public genomic
datasets** (ENA/SRA open-access studies, 1000 Genomes, HGDP, SGDP, etc.) in
DUNavigator, producing results that **withstand IRB review and meet journal data/
ethics requirements**. Written so an **IRB reviewer or co-author** can read it
standalone to understand the data handling.

## 2. Defining constraint — PII is out of scope, entirely

The academic profile **never ingests, derives, or stores personally identifiable
information**. Subjects are identified **solely by public accession** (e.g.
biosample `SAMEA…`/`ERS…`, run `ERR…`, study `PRJEB…`). There is no name, email,
kit number, MDKA, ancestor, or contact field — the FTDNA-side `ftdna_member` and
`mdka` tables are simply **not written** in this profile.

**Why this clears IRB:** open-access repositories like ENA/1000G/HGDP hold data
that is already de-identified and consented for unrestricted research use. A tool
that (a) only consumes such accessions, (b) stores no identifiers, and (c) makes no
attempt to re-identify is, for these datasets, typically **non-human-subjects /
exempt** work. This doc gives the reviewer the concrete basis for that
determination.

### 2.1 How the no-PII guarantee is enforced (policy + tested convention)

Chosen model (vs. structural enforcement): **same schema, disciplined importer,
proven by tests** —

- The ENA importer code path has **no inputs** that carry PII (accession metadata
  only) and **never** calls the `ftdna_member`/`mdka`/name/email writers.
- A **workspace-level `academic` mode flag** is set on academic projects; the UI
  badges it and hides PII-only surfaces.
- A **regression test** asserts that after an academic import, the PII columns
  (`biosample.donor_identifier` set to accession only; no `ftdna_member`/`mdka`
  rows; no email/name anywhere) are **null/absent**. This test is the artifact
  cited to the IRB as evidence the guarantee holds.
- *Trade-off vs. structural:* lighter to build, but the guarantee is procedural —
  if a future code path violated it, only the test would catch it. Revisit
  structural enforcement (separate PII-free schema / write-block) if a target IRB
  demands a physical guarantee. (Main doc §11.10.)

## 3. Identity & ingest (reuses existing infra)

- **Identity:** `external_id(source='ENA', external_id=<accession>)`, with
  `Biosample.sample_accession` carrying the canonical accession. The
  `(source, external_id)` UNIQUE is the dedup anchor — same as FTDNA, different
  source.
- **Ingest paths already present** (memory `project-import`):
  - `EnaClient` — ENA **Portal API** + **BioSamples API** metadata resolution.
  - `scan.rs` + `import_project_dir` — NAS directory scan for downloaded
    BAM/CRAM/VCF, keyed by accession-named subdirectories.
- **Acquisition is legitimately batch** here (unlike FTDNA §3.4): public APIs and
  bulk downloads are designed for programmatic, whole-study retrieval.

## 4. Provenance & reproducibility metadata (for publication)

Where the FTDNA profile records member metadata, the academic profile records
**citation- and reproducibility-grade provenance** (no PII), so a methods section
can be regenerated from the workspace:

```
-- per academic Subject / study (sketch; finalize with first real study)
study_accession        -- PRJEB…/PRJNA…
sample_accession       -- ERS…/SAMEA…
run_accession(s)       -- ERR…/SRR…
population / cohort     -- e.g. 'GBR', 'YRI' (1000G), 'HGDP00521'
data_use_conditions    -- declared access class (open / consented-research)
reference_build        -- already on Alignment (GRCh38/CHM13)
pipeline + params       -- aligner, caller, versions (for methods repro)
retrieved_at / source_url
```

- **Citation:** every result traces to study + sample + run accessions and the
  reference build — the journal "data availability" statement writes itself.
- **Reproducibility:** record pipeline + parameters + versions so a reviewer can
  reproduce from the same public inputs.
- **Data-use conditions:** capture each study's declared access/consent class so
  the workspace can flag anything that *isn't* open-access (see §6).

## 5. Collaboration & federation — the private tier collapses

Because there is no PII, the three-tier visibility model (main doc §8.1) reduces to
**local + federated (public)** — the project-shared *private* tier is unnecessary.
Consequences:

- Assertions (§8.4) over academic Subjects are **public** by default — research
  over open data, attributed and citable.
- Academic results are **freely federatable/publishable** (ancestry, coverage,
  haplogroup, branch assignments) — the *inverse* of the FTDNA PII constraint.
- Cross-researcher Subject resolution is by **accession + genetic** similarity (the
  reused IBD backbone); there is **no name-fuzzing** because no names exist.

## 6. Scope boundaries / risks

1. **Open-access only for v1.** **Controlled-access** sources (**EGA, dbGaP**)
   carry Data Use Agreements, consent limitations, and often *are* identifiable —
   they re-introduce exactly the IRB/PII complexity this profile avoids. **Out of
   scope for v1;** if pursued, they need DUA tracking, access-class enforcement, and
   likely the structural PII guarantee — a separate effort.
2. **Aggregate re-identification risk.** Even de-identified genomic data can be
   re-identified in principle (beacon/aggregate attacks). The profile must not
   publish per-individual genotypes beyond what the source already exposes; keep
   federated outputs at the same granularity/aggregation as the public source.
3. **Mixed workspaces.** If an academic project and an FTDNA project share one
   workspace, the `academic` flag + the null-PII test guard the academic project,
   but **publishing from a mixed workspace needs care** — only academic-flagged
   Subjects are publication-clean. Consider warning on cross-project publish.
4. **License/attribution.** 1000G/HGDP/ENA each have citation + acknowledgement
   expectations; surface them with the dataset so publications credit correctly.

## 7. Relationship to the FTDNA doc

| Concern | FTDNA profile | Academic profile (this doc) |
| --- | --- | --- |
| Shared engine | — defined in `ftdna-project-import.md` §§4–5, 8 — | reused as-is |
| Identity source | kit# | ENA accession |
| PII | present, private tier | **none, out of scope** |
| Acquisition | batch report CSV + per-member manual | batch via public APIs/NAS |
| Federation | gated by member consent | unrestricted (no PII) |
| New work | parsers, MDKA, collaboration ACL | **constrain** importer + provenance/citation metadata + no-PII test |

## 8. Next step

1. Pick a **reference study** (e.g. a 1000G subset on CHM13 already in the ancestry
   pipeline) to pin the provenance schema and write the **null-PII regression
   test** against a real academic import.
2. Specify the `academic` workspace flag + UI badge.
3. Defer EGA/dbGaP (controlled-access) to a separate design.
