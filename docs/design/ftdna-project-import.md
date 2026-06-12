# FTDNA Project Import — Backlog Design

**Status:** Design only (no implementation). Drafted 2026-06-06.
**Branch context:** `rust-rewrite`. Targets the Rust crates (`navigator-domain`, `navigator-store`, `navigator-app`, `navigator-ui`).

## 1. Goal

Let an FTDNA **group-project administrator** bulk-import a project they manage —
its member roster plus the associated genetic results (Y-STR, Y-SNP/Big Y, mtDNA,
and Family Finder where available) — into a DUNavigator project as Subjects with
their data sources attached.

The import must:

1. **Batch-create Subjects** from a member roster.
2. **Fill in metadata** FTDNA provides that we don't model yet (kit number,
   paternal/maternal ancestor, country of origin, FTDNA-reported haplogroups).
3. **Attach sequencing/marker data** (STR panels, SNP calls, mtDNA mutations).
4. **Match against existing Subjects** before creating duplicates — admins run
   **many projects** and the **same kit appears in several of them**.

### The larger frame

This is not just an importer. It **bootstraps a vendor-neutral genealogical
research platform** on top of FTDNA data — FTDNA is the *onboarding ramp*, not the
foundation. The durable object is a vendor-neutral research **Subject** (FTDNA
kit# is just one of its identifiers); co-admins on a shared project **collaborate
on the research** through the App Layer (AppView), reusing the existing
IBD-matching / aggregation backbone. The import (§§2–7) and the collaboration
layer (§8) are two halves of the same feature; the import seeds the data the
collaboration layer coordinates over.

### Product decisions (locked 2026-06-06)

| Decision | Choice |
| --- | --- |
| Subject ↔ project relationship | **One Subject, many projects (M:N)** — the FTDNA kit is the canonical person. |
| Match behavior | **Auto-merge on exact kit-number match; queue fuzzy matches for admin confirmation.** |
| Ingest scope | **Everything FTDNA exports** (members, Y-STR, Y-SNP/Big Y, mtDNA, Family Finder). Phased delivery below. |

## 2. What exists today (grounding)

| Concept | Type / table | Key | Notes |
| --- | --- | --- | --- |
| Subject | `Biosample` (`navigator-domain/src/workspace.rs:28`) / `biosample` | `guid` (UUID) | Single `project_id` FK today. Fields: `donor_identifier`, `sample_accession`, `description`, `center_name`, `sex`, `project_id`. |
| Y-STR | `StrProfile`/`StrMarker` (`strprofile.rs`) / `str_profile`+`str_marker` | `biosample_guid` | `panel_name` (Y-12…Y-700), `provider` (FTDNA…), `source`. Parser is **long-format** marker/value. |
| SNP | `VariantSet`/`VariantCall` / `variant_set`+`variant_call` | `biosample_guid` | contig/pos/ref/alt/rs_id/genotype. SNP-only. |
| Haplogroup | `RunHaplogroupCall` / `haplogroup_call` | `(biosample_guid, dna_type, source_key)` **UPSERT** | `dna_type` ∈ {Y, Mt}. Already idempotent per source — ideal for "FTDNA-reported" calls. |
| Chip/array | `ChipProfile` / `chip_profile` | `biosample_guid` | Provider auto-detect already includes FTDNA. Family Finder lands here. |
| mtDNA seq | `MtdnaSequence` / `mtdna_sequence` | `biosample_guid` | FASTA sequence. |
| Sequence run | `SequenceRun` / `sequence_run` | `biosample_guid` | platform/instrument/test_type/read metrics. |
| Batch import | `import_project_dir` (`navigator-app/src/lib.rs:2170`) | — | Dedupes biosample by `donor_identifier == sample_id` **within one project**. Append-only, idempotent by path. |

**Gaps this feature must close:**

- No kit number, no cross-project identity, no FTDNA member metadata.
- STR parser is long-format; FTDNA exports are **wide** (one row per person, one
  column per marker).
- `Biosample.project_id` is 1:1; the locked decision needs **M:N**.
- No "review/confirm before write" path — current importers write directly.

## 3. FTDNA export formats (VERIFIED against real exports)

Confirmed 2026-06-06 against three GAP exports from the **R1b-CTS4466Plus**
project (~1,879 members): `Member_Information`, `Paternal_Ancestry`, and
`YDNA_Results_Overview`. Exact headers below. (Redaction utility:
`/Users/jkane/Downloads/redact_ftdna.py` — blanks PII columns, CSV-aware.)

**Format gotchas observed (parser must handle):**
- Fields are quoted and **contain commas** (e.g. ancestor
  `"Joseph Abbett, b. 19 Mar 1819 and d. 2 Nov 1852"`) — must use a real CSV
  parser, never a naïve split.
- Quoting is **inconsistent across files**: Member/Ancestry exports quote every
  field; `YDNA_Results_Overview` is unquoted in the header but quoted per-cell in
  data rows, and **every STR value has a leading space** (`" 13"`).
- HTML entities appear in headers/values: `Name&darr;` (a "Name↓" sort arrow),
  `&gt;` for `>` in haplogroup paths, `&amp;`. **Normalize/unescape on import.**
- The `YDNA_Results_Overview` has **two leading non-member rows** before real
  data: a panel/min row (`" 00000. R-…"`) and a `" MIN"` row of minimum marker
  values. **Skip rows whose Kit Number isn't a real kit.**
- Null marker value = blank, `0`, or `-`.

### 3.1 Member_Information — the roster (spine of the import)

12 columns: `Kit Number(1)`, `Family Tree(2)` (YES/NO), `Name(3)` **[PII]**,
`Email(4)` **[PII]**, `Note(5)` **[PII, free-form]**, `Release(6)` (YES/NO),
`Kit Back(7)`, `Last Sign In(8)`, `Access Granted(9)` (Limited/Full/…),
`Allows MyHeritage Connection(10)`, **`Publicly Share DNA Results(11)`** (YES/NO),
`Remove From Group(12)`.

- **`Publicly Share DNA Results(11)` is a real per-member consent flag** — maps
  directly onto the visibility tiers (§8.1): it gates whether a Subject may move
  past project-shared toward anything federated. Plus `Access Granted(9)`.
- **`Note(5)` is free-form admin text** — in this sample admins use it for
  research coordination ("Administrator is …", "Old email address …", test-status
  notes). It frequently embeds **third-party** names/emails, so it **cannot be
  heuristically scrubbed** and is fully redacted for sharing. *Design relevance:*
  this is exactly the informal coordination the assertion/notes model (§8.4)
  should replace with structured, attributed, scoped annotations.

### 3.2 Paternal_Ancestry — MDKA source (maternal export is identical layout)

13 columns: `Kit Number(1)`, `Name(2)` **[PII]**, `Sub Group(3)`, `Email(4)`
(empty in this export), `Country(5)`, `Comment(6)` **[free-form]**,
**`Paternal Ancestor Name(7)`**, `Map Location(8)`, `Latitude(9)`,
`Longitude(10)`, `Family Tree(11)`, `Family Tree(12)` (dup header; tree links like
"MyHeritage WikiTree …"), `Remove From Group(13)`.

- **`Sub Group(3)` carries the clade/branch path**, e.g.
  `"20270. CTS4466>S1115>Z3023>FGC84010>…>BY34724>FT19839>FT22709"` (with `&gt;`).
  This is the project's branch assignment → feeds `belongs_to_branch` assertions
  (§8.4) and the shared clade tree. Leading number is the group's sort/label.
- **MDKA is fully present and richer than assumed:** `Paternal Ancestor Name(7)`
  embeds **birth/death inline** ("…, b. 19 Mar 1819 and d. 2 Nov 1852"), plus
  `Country(5)`, `Map Location(8)`, `Latitude(9)`, `Longitude(10)`. → populates the
  `mdka` table (§4.3), `lineage='Y'`. *Caveat:* in this sample most rows have
  `Map Location="No Location Saved"` and `Lat/Long="0"/"0"`, so geocoded fields
  are sparsely populated — parse name+dates reliably, coordinates opportunistically.
- The **maternal-ancestry export has the same column layout** (with
  `Maternal Ancestor Name`) → same parser, `lineage='Mt'`.

### 3.3 YDNA_Results_Overview — the wide Y-STR chart

109 columns: identity block `Kit Number(1)`, `Name(2)` **[PII]**,
`Paternal Ancestor Name(3)`, `Country(4)`, `Haplogroup(5)`, `Test(6)`,
`Subgroup(7)`, then **102 STR marker columns (8–109)**.

- Marker headers are **plain DYS names with multi-copy as ONE column**, not
  suffixed: `DYS385`, `DYS459`, `DYS464`, `CDY`, `YCAII`, `DYF395S1`, `DYS413`
  hold dash-joined palindromic values in a single cell (`" 10-14"`, `" 8-9-10"`,
  `" 12-15-15-17-17-17"`). The wide parser (§4.4) must **split multi-allele cells
  into a/b/c/d `StrMarker`s** rather than expect separate columns. Full observed
  marker order is in the header row of the file.
- `Haplogroup(5)` = FTDNA terminal label; `Subgroup(7)` repeats the clade path.
- Panel inference from count of populated markers (this export carries the Y-700
  superset, 102 columns).
- Skip the two leading non-member rows (panel/`MIN`).

### 3.4 What is NOT batch-extractable — the hard constraint

**FTDNA gives admins NO batch export of per-member testing data.** The three files
above (roster, ancestry, Y-STR overview) are **project-level report CSVs** — that
is the *entire* batchable surface. Everything deeper is **manual, one member at a
time**: the admin must **"pose as" the member** in GAP (subject to the per-member
**`Access Granted(9)`** level — Limited/Full) and download that member's results
individually. There is no API, no bulk endpoint.

Consequences for the design:

- **Deep data is single-file, per-member, manual:**
  - **Y-SNP / Big Y** (derived/ancestral/no-call + **novel variants**, hg38 coords),
  - **mtDNA** (HVR1/HVR2/coding mutations vs rCRS + haplogroup),
  - **Family Finder** raw autosomal,
  - any **BAM/VCF**,

  are obtained by posing-as-member and imported **one file at a time** — which the
  existing single-file path (`import_file`, `navigator-app/src/lib.rs:2082`)
  already handles. The "FTDNA import" feature does **not** bulk-pull these; its job
  is to **organize** the manually-retrieved files and attach them to the right
  Subject. Phases 3–5 are *per-member import + normalization*, not batch parsers.
- **`Access Granted(9)` is an acquisition gate, not just a consent label** — it
  determines whether the admin can pose-as and at what depth. **Store it per
  Subject** so the UI can show "deep data retrievable? Full / Limited / No," drive
  a per-member **work-list** of what's still manually pullable, and reflect that we
  only ever hold what the member granted. (Reinforces the consent basis, §11.7.)
- Terminal Y/mt haplogroup labels we already get from the batch files (`Haplogroup`,
  `Sub Group`); the *deep* SNP/variant detail is the manual part.

**Big Y formats are now confirmed** (§3.5). mtDNA and Family Finder per-member
examples are still useful to lock those single-file parsers (Phases 4–5).

### 3.5 Big Y artifacts (VERIFIED) — and access level *is* the data tier

Confirmed 2026-06-06 against a real Big Y. **Which artifacts an admin can pull is
exactly the `access_granted` level** — so that field doubles as "what Big Y data is
reachable":

| `access_granted` | Big Y artifacts | Format |
| --- | --- | --- |
| **Limited** | `<kit>_BigY_Named_Variants.csv`, `<kit>_BigY_Private_Variants.csv` | derived/ancestral SNP **lists** (CSV) |
| **Advanced** | `bigy2-<sampleUUID>/` archive: `variants.vcf` + `regions.bed` + `readme.txt`; **BAM** (separate, pending) | full per-site **VCF** + callable **BED** |

**Build = hg38 — and the readme lies.** `readme.txt` is stale boilerplate claiming
"build 37.3"; the VCF header is authoritative: `##reference=ucsc.hg38.fasta`,
`##contig=<ID=chrY,length=57227415,assembly=ucsc.hg38>` (+ `chrY_KI270740v1_random`
decoy). **Parser must read build from the VCF header, never the readme**, then
lift hg38→CHM13 (partially resolves §11.2). CSV positions are bare (no contig) but
**match the VCF exactly → also hg38, chrY**.

**Named_Variants.csv** — `SNP_Name, Position, On_Haplotree, Ancestral, Derived`.
996 known Y-SNPs (420 `On_Haplotree=Yes` = FTDNA tree-recognized, 576 No). These
are the branch-defining calls → drive haplogroup placement.

**Private_Variants.csv** — `Position, Ancestral, Derived`. Novel/unnamed derived
variants (7 here) → the "private SNPs" a project hunts for new branches.

**variants.vcf** — VCFv4.1, single sample (column name = the sample UUID), chrY +
decoy. 263,632 records / **206,999 PASS**. **`FILTER` is `PASS` or a fail-reason
string** (`QUAL=…`, `GTL=…;QUAL=…`) — keep **PASS only** (readme confirms). Per-site
genotype: `GT 1/1` = **derived** (~384 PASS, the real positives), `0/0` =
**ancestral-confirmed** (~205.5k), `ALT="."` = ref-only ancestral; `0/1`/`1/2`
(~1k) are het/multiallelic **artifacts on haploid Y — flag/drop**. Rich FORMAT
(AD/DP/GQ/AB/…) available for quality gating.

**regions.bed** — standard 0-based half-open; chrY(+decoy); the Big Y **callable
footprint**. Use to distinguish *no-call* from *ancestral* when running our own
caller. **Footprint scales with test version** (below).

**BY500 (`bigy2-…`) vs BY700 (`bigy3-…`) — same format, different size.** Both
verified; the VCF/BED schema is **identical** (VCFv4.1, `ucsc.hg38.fasta`, same
FORMAT/INFO, sample column = UUID, readme still mislabels build) — **one parser
handles both**. They differ only in footprint: BY500 ≈ 9.7 Mb / 12,171 BED
intervals / 263k VCF records; BY700 ≈ 14.8 Mb / 22,136 intervals / 233k records.
So **don't hardcode region counts**; **the BED footprint is the reliable version
signal** (the `@RG LB` tag is not — see BAM below).

**Advanced-tier BAM** (`<kit>-<sampleUUID>.bam` + `.bai`, ~624 MB) — verified:
- **Aligned to full hg38**: 455 `@SQ` contigs, `AS:ucsc.hg38`, **with `M5`
  checksums** (use them to verify/resolve the exact reference). `chrM`=16569,
  `chrY`=57227415.
- **Reads are *usually* chrY-only — but detect, don't assume**: in all three
  verified BAMs `idxstats` shows reads only on `chrY` (~8.9 M here), **0 elsewhere
  (incl. `chrM`)** — the whole-genome `@SQ` header is kept but **off-target reads
  are deliberately stripped** so Big Y can't double as a "free" (mostly-accurate)
  mtDNA test. **Exception:** some **first-generation Big Y** BAMs predate that
  policy and **retain `chrM` (and other off-target) reads** — rare, but real. So the
  importer must **probe `idxstats` for actually-covered contigs** rather than
  hardcode Y-only: treat as Y-restricted by default, but **if `chrM` carries reads,
  opportunistically offer mtDNA extraction** (§4.5). Coverage/sex/SV analysis should
  likewise scope to whatever contigs actually have reads.
- **Coordinate-sorted** (`@HD SO:coordinate`) + `.bai` present → directly usable.
- **`@RG LB` is an UNRELIABLE version signal** (confirmed across two kits): the
  BY700 BAM has `LB:unknown-library-Big Y-700`, but a **BY500 BAM has a bare
  `LB:unknown-library`** with no suffix. So treat `LB` as *best-effort* test_type
  (use it if it carries a "Big Y-…" suffix, else `test_type="Big Y (version
  unknown)"`). **Definitive version comes from the paired VCF/BED footprint**, or
  the admin/order context — **not the BAM alone**. (`SM` is `GRC…<kit>` but the kit
  isn't a clean fixed-prefix substring — don't parse kit# from it; the dir name is
  the kit.)
- **Original aligner is GONE**: FTDNA **reheadered/stripped** the BAM, so `@PG`
  lists only `samtools reheader`, not the original mapper. **`probe.rs` aligner
  auto-detect will fail → record `aligner = "unknown"`** for FTDNA BAMs (don't
  guess). Platform/test still inferable from `LB`.
- **Join key:** the **sample UUID in the filename** (`5a291bcb-…`) matches the
  `bigy3-<UUID>/` package and the VCF sample column → that UUID, **not** the `@RG
  SM` (an internal `GRC…` id) nor the dir name, links BAM ↔ variants for one kit.

**Delivery = ZIP.** Both the variant package (`bigy*-<UUID>/`) and the BAM package
(`<kit>/`) arrive as **ZIP archives**. The importer must **accept a ZIP, sniff its
contents** (`variants.vcf`+`regions.bed` ⇒ Big Y package; `*.bam`+`*.bai` ⇒ BAM
package), **extract, and route** — see §6.

**Coordinate/build note:** FTDNA reports on **hg19/hg38**. SNP positions and Big Y
novel variants (§3.4) must be lifted to our working build (CHM13/T2T) via the
existing `navigator-refgenome` liftover before storage, or stored with an explicit
`build` tag and lifted lazily. Decide per data type (open question §11.2).

## 4. Data model changes

### 4.1 M:N Subject ↔ Project (new join table)

Today `biosample.project_id` is a single FK. Introduce membership:

```
-- migration 0014_biosample_project_membership
CREATE TABLE biosample_project (
    biosample_guid TEXT NOT NULL REFERENCES biosample(guid),
    project_id     INTEGER NOT NULL REFERENCES project(id),
    role           TEXT,          -- optional: subgroup/branch label within the project
    added_at       TEXT NOT NULL, -- ISO-8601
    PRIMARY KEY (biosample_guid, project_id)
);
```

**Migration strategy:** backfill one membership row per existing
`biosample.project_id`. Keep `project_id` column for one release as a
"primary/home project" (nullable) to avoid a big-bang UI refactor, OR drop it and
make all reads go through the join. **Recommendation:** keep `project_id` as a
nullable "home project" pointer, treat `biosample_project` as the source of truth
for membership. Decide during implementation.

Store API additions (`navigator-store`):
- `biosample_project::add(pool, guid, project_id, role)` (idempotent upsert)
- `biosample_project::list_projects_for(pool, guid)`
- `biosample::list_for_project` switches to a join through `biosample_project`.

### 4.2 Subject identity — vendor-neutral (decision: generalize now)

The platform is **not** locked to FTDNA — FTDNA is the onboarding ramp. So a
Subject's vendor identifiers are modeled generically, and the FTDNA kit number is
**one identifier among many** (YSEQ, Nebula, direct-WGS sample id, manual). The
`(source, external_id)` pair is the **global cross-project dedup anchor**:

```
-- migration 0015_subject_identity
CREATE TABLE external_id (
    id             INTEGER PRIMARY KEY,
    biosample_guid TEXT NOT NULL REFERENCES biosample(guid),
    source         TEXT NOT NULL,        -- 'FTDNA' | 'YSEQ' | 'NEBULA' | 'WGS' | 'MANUAL'
    external_id    TEXT NOT NULL,        -- kit number / vendor id
    UNIQUE (source, external_id)         -- one Subject per (vendor, id)
);
CREATE INDEX idx_external_id_lookup ON external_id(source, external_id);

-- FTDNA-specific extras (only the labels FTDNA reports; ancestry → MDKA §4.3)
CREATE TABLE ftdna_member (
    biosample_guid      TEXT PRIMARY KEY REFERENCES biosample(guid),
    member_name         TEXT,            -- member/contact name as in GAP
    y_haplogroup_ftdna  TEXT,            -- as reported by FTDNA (label only)
    mt_haplogroup_ftdna TEXT,
    haplo_status        TEXT,            -- predicted | confirmed
    access_granted      TEXT,            -- 'Advanced'|'Limited'|'None' — pose-as gate; also the Big Y data tier (§3.5)
    publicly_shares     INTEGER          -- consent flag (col 11), 0/1; gates federation
);
```

Rationale: `UNIQUE(source, external_id)` enforces "one Subject per vendor id" and
is the key the matching engine (§5) and the cross-admin resolver (§8) both use.
FTDNA-reported haplogroups are stored as **labels** in `ftdna_member`; any
*computed* haplogroup from our own analysis stays in `haplogroup_call` (different
provenance — keep them distinct). New vendors only add `source` values + their own
optional details table; no schema churn.

### 4.3 MDKA — Most Distant Known Ancestor (PRIVATE, vendor-agnostic)

MDKA (a.k.a. FTDNA's "Earliest Known Ancestor") is genealogy, not genetics: the
earliest documented ancestor on a given lineage, with name, dates, place of
origin, and often **geocoded coordinates** so it can be pinned on a map. It is
**not** FTDNA-specific — a Subject may acquire MDKA data from FTDNA import, manual
entry, or another vendor — so it gets its own Subject-level table rather than
living on `ftdna_kit`.

A Subject can have one MDKA **per lineage** (paternal/Y, maternal/mt, and
optionally an autosomal/"general" entry):

```
-- migration 0016_mdka
CREATE TABLE mdka (
    id              INTEGER PRIMARY KEY,
    biosample_guid  TEXT NOT NULL REFERENCES biosample(guid),
    lineage         TEXT NOT NULL,        -- 'Y' | 'Mt' | 'Auto'
    ancestor_name   TEXT,
    birth_year      INTEGER,              -- nullable; FTDNA often gives only a year
    death_year      INTEGER,
    origin_place    TEXT,                 -- free-text place as entered (e.g. "Cork, Ireland")
    origin_country  TEXT,                 -- normalized country (for grouping/maps)
    latitude        REAL,                 -- geocoded; nullable
    longitude       REAL,
    source          TEXT,                 -- 'FTDNA' | 'MANUAL' | ...
    notes           TEXT,
    updated_at      TEXT NOT NULL,        -- ISO-8601
    UNIQUE (biosample_guid, lineage)      -- one MDKA per line per Subject
);
CREATE INDEX idx_mdka_biosample ON mdka(biosample_guid);
```

**Privacy posture (important).** MDKA is the most sensitive data in this feature —
it names living-adjacent real people and places. Its visibility tier is
**project-shared (private)**, never **federated (public)** — see the three-tier
model in §8.1:

- It is **never** published as a world-readable PDS / `du-domain::fed` record
  (contrast with the ancestry/coverage/biosample/seqrun aggregate records that
  *are* federated — see memory `fed-report-records`). AT Proto is public-by-default,
  so PII must not ride that rail.
- It **may** be shared with the **co-admin team of a project**, but **only over the
  encrypted Edge-to-Edge (P2P) channel** (§8.2) — **never stored in AppView**, which
  is anonymized-only. AppView merely brokers the consented connection; the MDKA
  ciphertext flows admin-to-admin and is decrypted locally. Continuous with how
  co-admins see member data in FTDNA GAP, but without a central custodian.
- Mark `mdka`/`ftdna_member`/`external_id` in code as **PII / never-leaves-encrypted**
  so neither the public `fed` record-derivation **nor any AppView sync** reads them;
  they may only enter the P2P encrypted payload (§8.2).

Store API additions: `mdka::upsert(pool, guid, lineage, fields)`,
`mdka::list_for(pool, guid)`. Roster import (§3.1) populates MDKA; the Subject
detail UI gets a small editable MDKA panel (paternal/maternal), and the coordinates
feed a future origins-map view.

> **FTDNA mapping (confirmed §3.2):** the `Paternal_Ancestry` export →
> `mdka(lineage='Y')`, the maternal export (identical layout) → `mdka(lineage='Mt')`.
> Column map: `Paternal/Maternal Ancestor Name(7)` → `ancestor_name` **+ parse
> inline `b.`/`d.` dates** into `birth_year`/`death_year`; `Country(5)` →
> `origin_country`; `Map Location(8)` → `origin_place` (drop the
> `"No Location Saved"` sentinel); `Latitude(9)`/`Longitude(10)` →
> `latitude`/`longitude` (drop `0/0` sentinel). Coordinates are sparsely
> populated, so name+dates are the reliable fields.

### 4.4 STR — wide-format parser (new), existing storage

`StrProfile`/`StrMarker` storage is fine as-is. Add a **wide-format CSV parser**
alongside the existing long-format `strprofile::parse_csv`. Confirmed against
`YDNA_Results_Overview` (§3.3):

- Input: header row of identity columns (1–7) + 102 marker columns; one data row
  per kit. **Skip the two leading non-member rows** (panel / `MIN`) — gate on a
  real `Kit Number`.
- Per data row → `(identity, Vec<StrMarker>)`. Identity from cols 1–7
  (`Kit Number`, `Name`, `Paternal Ancestor Name`, `Country`, `Haplogroup`,
  `Test`, `Subgroup`).
- **Trim leading spaces** on every value (`" 13"` → `13`); null = blank / `0` / `-`.
- **Multi-copy markers arrive as ONE dash-joined cell**, not suffixed columns:
  split `DYS385 " 10-14"` → `DYS385a=10, DYS385b=14`; same for `DYS459`, `DYS464`
  (a–d), `CDY`, `YCAII`, `DYF395S1`, `DYS413`, etc. Normalization map keyed on the
  plain DYS header.
- **Unescape HTML entities** (`&gt;` etc.) in `Haplogroup`/`Subgroup`.
- Panel inference from count of populated markers (12/25/37/67/111/700) →
  `panel_name`; `provider = "FTDNA"`, `source = "IMPORTED"`.

### 4.5 SNP / Big Y, mtDNA

**Big Y (formats confirmed §3.5).** All variants land in `variant_set`/
`variant_call` (`contig='chrY'`, `position`, `reference=Ancestral`,
`alternate=Derived`, `rs_id=SNP_Name`, `genotype` from GT), tagged with the source
tier, **after hg38→CHM13 liftover** (`navigator-refgenome`). The two access tiers
import differently but normalize to the same store:

- **Limited tier (CSVs):** `Named_Variants.csv` → `variant_set` (source_label
  `"FTDNA BigY Named <kit#>"`); `rs_id=SNP_Name`, keep `On_Haplotree`. `Private_
  Variants.csv` → `variant_set` (`"FTDNA BigY Private <kit#>"`), `rs_id` null. Both
  are **derived-only** lists.
- **Advanced tier (VCF/BED + BAM) — store *positives + BED* (decided):** parse
  `variants.vcf` keeping **`FILTER=="PASS"` only**; persist **only derived
  positives** (`GT 1/1`) to `variant_set`, and **drop `0/1`/`1/2` het artifacts**.
  **Do NOT persist the ~190k ancestral `0/0` calls.** Instead store `regions.bed`
  as a callable-regions artifact: a site **inside the BED but not among our stored
  positives ⇒ ancestral-confirmed**; a site **outside the BED ⇒ no-call/uncovered**.
  BED + positives reconstructs ancestral-vs-no-call **losslessly without ~50× the
  rows per subject**. Reuse the existing VCF reader (`parity::parse_truth_vcf`); add
  a thin Big Y adapter for the FTDNA FILTER/FORMAT quirks. **Read build from the VCF
  header, not the readme.** The **BAM** (when present) imports via the existing BAM
  path (`import_file` → `probe.rs` → run + alignment): `reference_build=hg38`,
  `test_type` from `@RG LB` **if it carries a "Big Y-…" suffix, else "Big Y (version
  unknown)"** (§3.5), **`aligner="unknown"`** (header stripped); it is **chrY-only**,
  so downstream analysis must treat it as Y-restricted.
- **Terminal haplogroup, two provenances (keep distinct):** FTDNA's reported label
  (from the batch Y-STR `Haplogroup`, §3.3) → `haplogroup_call` upsert
  `source_key="ftdna:kit:<kit#>"`. Our **own** call computed from the derived named
  SNPs against the tree (`haplo.rs`) → its own `source_key` (e.g. `"bigy:<kit#>"`).
  Both visible via the existing reconciliation.

**mtDNA** — **usually not obtainable from the Big Y BAM** (off-target reads stripped,
§3.5), so the primary source is a separate FTDNA **mtDNA / Full Mitochondrial
Sequence (FMS)** test. **Exception:** rare **first-gen Big Y** BAMs retain `chrM`
reads — when `idxstats` shows `chrM` coverage, **opportunistically call mtDNA from
the BAM** via the existing mtDNA path (rCRS-relative variants + sequence into
`mtdna_sequence`), flagged lower-confidence than a true FMS. Either source →
mutation list as a variant set vs rCRS and/or sequence into `mtdna_sequence`; mtDNA
haplogroup → `haplogroup_call` (`dna_type = Mt`). *(Per-member FMS format still to
confirm — §3.4.)*

## 5. Matching / dedup engine

The heart of the feature. Runs **before** any write, producing a **plan** the
admin can review.

**Inputs:** parsed roster rows (each with kit#, name, ancestors, country,
haplogroups) + the target project.

**Match precedence per incoming row:**

1. **Exact vendor id** → look up `external_id WHERE source='FTDNA' AND
   external_id=<kit#>`. If found → **auto-merge**: reuse that biosample, add
   membership to the target project, update/fill metadata, attach data.
   *(Locked behavior.)* This is the **local** half of the same resolution the
   cross-admin resolver runs at AppView (§8.3).
2. **No id match → fuzzy candidates** within workspace, scored on:
   - normalized name similarity (paternal-ancestor name is often the de-facto
     name in FTDNA),
   - same Y or mt haplogroup label,
   - same country.
   Above a threshold → **queue as "possible match"** for admin confirm/reject.
   Below → **treat as new Subject**.
3. **New Subject** → create `Biosample` (`donor_identifier = kit#` or member
   name), `external_id` row, `ftdna_member` row, MDKA rows, project membership,
   then attach data.

**Conflict handling on merge:** when filling metadata that already exists and
differs (e.g., different country), don't silently overwrite — record both and
flag in the import report (prefer existing non-null unless admin chose "prefer
incoming"). Append-only data (STR/SNP/mtDNA) attaches as additional sources;
existing reconciliation handles multi-source consensus.

**Cross-file join key:** within one import session, the **kit number** joins the
roster row to its STR/SNP/mtDNA rows across the multiple files. Rows whose kit#
isn't in the roster → report as "data without a member" (still importable as a
bare Subject, flagged).

**Big Y exception — UUID-only, admin-attached.** The Big Y *variant package* and
*BAM* carry only an internal **sample UUID**, no kit# (only the Limited-tier CSVs
are named `<kit>_…`). There is no in-file UUID→kit# map. So raw Big Y artifacts are
**associated to a Subject by the admin at import** (they downloaded them while posed
as that member, §3.4) — a "attach this Big Y package/BAM to <Subject>" step, not an
automatic kit-match. The **UUID then links BAM ↔ variant package** to each other
within that member, and we may persist it as `external_id(source='FTDNA_BIGY_UUID',
…)` for re-import idempotency.

## 6. Import flow (proposed)

```
app::import_ftdna_project(
    project_id: i64,                 // target DUNavigator project
    files: Vec<FtdnaFile>,           // {path, declared_kind?} — kind auto-detected if absent
    options: FtdnaImportOptions,     // fuzzy threshold, prefer-incoming/existing, dry_run
) -> FtdnaImportPlan | FtdnaImportSummary
```

0. **Unpack archives.** Per-member Big Y / BAM deliverables arrive as **ZIPs**
   (§3.5). Accept a `.zip`, **sniff contents to classify**: `variants.vcf`+
   `regions.bed` ⇒ Big Y variant package; `*.bam`+`*.bai` ⇒ BAM package; extract to
   a temp dir and feed the inner files into step 1. The **sample UUID** in the dir/
   filename is the join key linking a BAM to its variant package.
1. **Classify** each file (roster / Y-STR / Y-SNP CSV / Big Y package / BAM /
   mtDNA / FF) by header/content sniffing (mirrors `chipprofile` provider detection).
2. **Parse** each into typed rows keyed by kit# (or sample UUID → kit# for Big Y
   artifacts that carry only the UUID).
3. **Join** rows across files by kit#.
4. **Match** each kit against the workspace (§5) → build a **plan**:
   `{ new: [...], auto_merge: [...], needs_confirm: [...], orphan_data: [...] }`.
5. **Present plan** in UI (dry-run by default). Admin resolves the
   `needs_confirm` queue.
6. **Commit** in a transaction-ish batch: create/merge Subjects, add memberships,
   upsert metadata + haplogroup calls, create STR/SNP/mtDNA/chip records.
7. **Return summary**: counts of created/merged/skipped, conflicts, orphans,
   liftover failures, unrecognized markers.

`FtdnaImportSummary` mirrors the existing `ProjectImportSummary` shape
(`navigator-app/src/lib.rs:452`) for UI consistency.

## 7. UI (egui Workbench)

Extends the in-progress Workbench redesign (see memory `ui-workbench-redesign`).

- Entry point: project action bar → **"Import FTDNA Project…"** (multi-file
  picker).
- **Review screen** (the key new surface): three sections — *New Subjects*,
  *Auto-merged (kit match)*, *Needs confirmation (fuzzy)* — the last with
  side-by-side incoming-vs-existing diff and Approve/Reject/"It's new" per row.
- Post-import **summary panel** reusing the Full Analysis modal styling.

## 8. Collaboration / App Layer architecture

This is the part that turns the importer into a **genealogical research platform**.
FTDNA is the onboarding ramp; the durable object is a **vendor-neutral research
Subject** that co-admins collaborate on.

**Decision (2026-06-06): AppView holds NO PII — it is a pure broker.** AppView keeps
its deliberate *anonymized/aggregate-only* posture (it never stores member names,
MDKA, or kit linkage). The collaboration layer **reuses the IBD-matching
architecture wholesale**: AppView coordinates **discovery, consent, and key
exchange** and persists **match/assertion *state*** (PII-free, pseudonymous), while
the **PII itself flows admin-to-admin over an encrypted Edge-to-Edge (P2P) channel**
— exactly the mechanism the IBD system uses for genetic comparison
(`decodingus/documents/planning/ibd-matching-system.md`: ECDH X25519 session key +
AES-256-GCM, AT-Proto-brokered handshake, P2P/relay transport). One encrypted
exchange substrate serves both IBD comparison and genealogy-PII exchange.

### 8.1 Three visibility tiers

Every datum sits in exactly one tier. The middle tier is private *to a consented
circle* and **never touches AppView storage** — it moves peer-to-peer.

| Tier | Examples | Where it lives | Who sees it |
| --- | --- | --- | --- |
| **Local** | raw alignments, scratch notes | local SQLite only | the one admin |
| **Project-shared (private)** | Subject roster, `external_id` links, MDKA, branch/clade assignments, annotations, tasks | **exchanged Edge-to-Edge, E2E-encrypted; held only in each participating admin's local store.** AppView brokers (discovery/consent/keys) but **stores no PII** | the project's **admin team** (membership brokered + audited by AppView; payload AppView can't read) |
| **Federated (public)** | ancestry estimate, coverage summary (aggregate, non-PII) | PDS `du-domain::fed` records → AppView ingest | world-readable |

**Why neither AppView nor the PDS rail holds PII:** AT Proto records are
public-by-default (so PII can't be a PDS record), **and** AppView is anonymized-only
by policy (so PII can't be an AppView row either). The only place PII lives off the
originating machine is **inside the encrypted P2P payload** between consented admins.

### 8.2 Substrate — broker + encrypted P2P

- **Non-PII attributed assertions** (branch assignments by opaque subject id,
  pseudonymous `same_person` links, clade hypotheses, aggregate stats) → published
  as **PDS records** in the asserting admin's repo and ingested by AppView du-jobs —
  same path as today's fed records, new record types. AppView may fold these into a
  PII-free `current_view`.
- **PII-bearing data** (member names, MDKA, contact info, kit↔identity linkage, raw
  STR/SNP for comparison) → **never sent to AppView**. AppView brokers an
  introduction between two project admins (verifies both are members, relays the
  **ECDH key-exchange** messages and consent, records that an exchange occurred);
  the **ciphertext flows Edge-to-Edge** (P2P, or via a blind relay that sees only
  ciphertext), decrypted only in each admin's Navigator. **Same crypto + handshake
  as IBD** — build once, use for both.

So each admin's local store is the working copy; **there is no shared server-side
copy of the private tier** — co-admins reach a shared view by exchanging encrypted
assertions P2P. AppView holds only the PII-free coordination state. Sync is one
project at a time, opt-in, per consented peer.

### 8.3 Cross-admin Subject resolution (reuses IBD backbone)

> **Mechanism specified in AppView design D2**
> (`decodingus/documents/planning/d2-research-subject-registry.md`), which **corrects
> an earlier sketch here** (AppView-stored salted `id_hashes` — rejected: kit numbers
> are enumerable, so any hash the broker can see is brute-forceable).

Two admins in a **shared project** → the same person. Resolution runs **without
AppView ever seeing an identifier or a hash of one**:

1. **Deterministic (exact):** co-admins exchange their `(source, external_id)` lists
   **over the encrypted D1 channel** (within the consented project scope —
   GAP-equivalent) and compute the intersection **locally**; matches ⇒ same person →
   agree a shared `research_subject_id`. AppView learns only that a pseudonymous
   subject gained a membership.
2. **Genetic:** ids differ but the **IBD / haplotype comparison** (Edge-to-Edge, §8.2,
   the IBD pipeline) says same-person/close-kin → *suggested* merge, never auto.
3. **Assertion-mediated:** a pseudonymous `same_person(subject_a, subject_b)` assertion
   (no kit#); the group accepts/rejects; provenance retained.

`ResearchSubject` is the AppView-side **PII-free** node: `{research_subject_id,
custody_did}` + project memberships — **no names, no kit numbers, no MDKA, and no
hashes** (D2). The clear-text identity and the `external_id↔research_subject_id` map
live only in each admin's local store and cross between admins via the D1 channel.
Cross-*project* linking is **member-claim only**, never a silent AppView merge (D2 §4.4).

### 8.4 The assertion model (the collaboration primitive)

Co-admin research is modeled as **attributed, scoped assertions** rather than
direct mutation of shared rows. One shape, reused everywhere:

```
Assertion {
  id,
  subject:    research_subject_id,
  predicate:  'same_person_as' | 'belongs_to_branch' | 'mdka_is' | 'haplogroup_is' | 'note' | ...
  value:      <predicate-specific payload>,
  author_did: <admin's DID>,        -- attribution
  scope:      'project:<id>' ,       -- visibility/consent boundary
  evidence:   optional (STR distance, SNP, doc citation),
  created_at,
  retracted_by / supersedes:        -- assertions are append-only + retractable
}
```

This one model buys four things at once:

- **Multi-admin conflict handling with provenance** — two admins disagree → two
  assertions, both visible with author + time; the group resolves, nothing is
  silently overwritten. **PII-free** assertions fold into an AppView `current_view`;
  **PII-bearing** assertions (`mdka_is`, `note` with a name, …) are exchanged
  **P2P-encrypted** and folded **locally** in each admin's Navigator — AppView never
  sees them.
- **Scope = consent + transport** — the `scope` field is the privacy boundary, and
  it **decides the rail**: `public` → PDS record; `project:<id>` non-PII → AppView
  coordination state; `project:<id>` **PII → encrypted P2P only**, never persisted
  server-side.
- **Vendor-independence** — assertions reference `research_subject_id`, never a kit
  number.
- **Future member self-claim** (below) — custody transfer is just another assertion.

### 8.5 Governance & the member-claim path

- **Scope = co-admin team now** (locked): on import, PII may be shared to the
  project's admin team — continuous with FTDNA GAP — but **via encrypted P2P
  exchange, not a shared server copy** (§8.2). AppView brokers membership/consent and
  audits that an exchange happened; it never holds the data.
- **Custody / stewardship:** an imported Subject is **admin-stewarded** (the
  importing admin is custodian). Future path: the actual member onboards, proves
  control of the kit, and **claims** the Subject — custody transfers to their DID,
  and *they* decide federation/sharing of their own data. This is the consent
  endgame and a platform selling point; schema leaves room (`custody` on
  `ResearchSubject`) but the claim flow is out of scope for the first cuts.

### 8.6 New AppView surfaces (sketch) — all PII-free

- **PII-free** `ResearchSubject` registry: opaque ids + `id_hashes[]` + DIDs +
  custody (no names/kit#/MDKA).
- **Non-PII** `Assertion` store + `current_view`; PII assertions are **not** stored
  here.
- Project **admin team** membership + roles (owner/editor/viewer) + audit log of
  *coordination events* (who exchanged with whom, when) — not the payloads.
- **Broker endpoints**: peer discovery within a project, ECDH key-exchange relay,
  consent handshake, exchange-occurred attestation — modeled on the IBD
  request/consent/attestation flow.
- Reuse existing **messaging** for non-PII discussion; PII-bearing messages go P2P.
- Branch/clade tree as a shared structure (assignments are pseudonymous assertions).

> These are AppView (separate repo) changes — see the companion roadmap
> `decodingus/documents/planning/design-roadmap-rust-rewrite.md`. This doc scopes the
> **Navigator-side** contract: Navigator runs the **encrypted P2P endpoint** (key
> exchange + ciphertext exchange + local fold), publishes **PII-free** assertions/
> hashes to AppView, and consumes the brokered `current_view` + suggested merges.
> **No PII payload is ever sent to AppView.**

## 9. Import profiles (sources)

The same engine serves more than FTDNA. An **import profile** is a thin adapter
over the shared spine (vendor-neutral Subject + `external_id`, M:N projects, STR/
SNP/haplogroup stores, matching engine §5, assertions §8). Two profiles are in
scope; more (YSEQ, Nebula, direct-WGS) slot in later as just a new `source` value
plus a parser.

| | **FTDNA profile** (§§1–8) | **Academic / ENA profile** |
| --- | --- | --- |
| Identity | `external_id(source='FTDNA', kit#)` | `external_id(source='ENA', accession)` — ERS/SAMEA/ERR; reuse `Biosample.sample_accession` |
| Metadata | `ftdna_member` + **MDKA** (PII) | study/cohort/population/data-use accession metadata — **no PII** |
| Acquisition | roster batch + per-member manual (§3.4) | public repos: **ENA Portal API + BioSamples** (existing `EnaClient`), NAS scan (`scan.rs`, `import_project_dir`) |
| Matching | kit# exact + fuzzy name/ancestor (§5) | accession exact + genetic; **no name fuzzing** (no names exist) |
| Privacy tiers (§8.1) | local / **project-shared-private** / federated | local / **federated** only — the private PII tier **collapses** |
| Collaboration | assertions, scoped to admin team | assertions, **public** — research over open data |

**The academic profile's defining rule: PII is out of scope entirely** — to survive
IRB review and journal publishing requirements. The ENA importer **never** writes
`ftdna_member`, `mdka`, member name, or email; a Subject is identified solely by its
public accession. Per the chosen guarantee model (**policy + tested convention**):
same schema, but the academic importer never populates PII columns, enforced by a
test asserting those columns stay `NULL` for academic-profile imports, plus a
workspace-level "academic" flag surfaced in the UI. Because nothing is PII, the
academic profile's data is **freely federatable/publishable** (the inverse of the
FTDNA constraint).

> **Much of this profile already exists** — `EnaClient` (Portal + BioSamples
> resolution), `scan.rs`, and `import_project_dir` are the ENA path's foundation
> (memory `project-import`). The academic profile is mostly *constraining* the
> existing importer (no-PII guarantee, provenance/citation metadata) rather than
> new ingest. Its IRB-/publishing-facing specifics live in a companion doc:
> **`docs/design/academic-ena-import.md`**.

## 10. Phasing (delivery slices)

Even though scope is "everything," ship in safe increments. **Phases 0–5 are
local/single-admin** (the importer); **Phases 6–8 light up collaboration** and can
proceed in parallel once the identity model (Phase 0) lands.

- **Phase 0 — schema & identity:** migrations `0014` (M:N) + `0015`
  (`external_id` + `ftdna_member`) + `0016` (`mdka`); store APIs; backfill; mark
  `external_id`/`ftdna_member`/`mdka` as **PII / non-federated** so the public
  `fed` layer skips them. Vendor-neutral from day one. No UI yet.
- **Phase 1 — roster + matching engine:** parse members CSV (including MDKA
  fields → paternal/maternal `mdka` rows, `access_granted`/`publicly_shares`),
  build/commit the match plan keyed on `external_id`, M:N membership, editable MDKA
  panel in Subject detail, dry-run + review UI. *Spine and riskiest logic — land it
  first.*
- **Phase 2 — Y-STR wide-format:** wide parser + marker normalization + panel
  inference; attach to matched Subjects.

  > **Phases 3–5 are per-member, NOT batch** (§3.4): FTDNA exposes no bulk export
  > of deep test data — the admin manually poses-as-member (gated by
  > `access_granted`) and pulls one file at a time. These phases harden the
  > **single-file** path (`import_file`) per data type + a per-Subject "deep data
  > retrievable?" work-list, rather than building batch parsers.

- **Phase 3 — Y-SNP / Big Y (per-member, formats confirmed §3.5):** start with the
  **Limited tier** (Named/Private Variant CSVs → `variant_set`) — smallest, no VCF
  machinery — then **Advanced tier**: ZIP unpack + sniff (§6 step 0), PASS-filtered
  `variants.vcf` + callable `regions.bed`, and the **chrY-only BAM** via the existing
  BAM path (`aligner=unknown`, `test_type` from `@RG LB`, UUID-attached §5). One
  parser covers BY500/BY700. hg38→CHM13 liftover; own-vs-FTDNA haplogroup both
  recorded (§4.5).
- **Phase 4 — mtDNA (per-member):** single-file mtDNA results → mutation lists +
  haplogroup.
- **Phase 5 — Family Finder (per-member):** single-file FF raw → `ChipProfile`
  (FTDNA provider already auto-detected).
- **Phase 6 — private-tier sync:** Navigator ↔ AppView project-shared sync
  (Subjects, `external_id`, MDKA, memberships) over the authenticated
  access-controlled API; audit log. *Depends on a parallel AppView design doc.*
- **Phase 7 — assertions & cross-admin resolution:** publish/consume assertion
  records; AppView folds `current_view`; suggested merges via the reused
  IBD-resolution backbone; conflict UI with provenance.
- **Phase 8 — coordination surfaces:** shared branch/clade tree, discussion
  threads (reuse messaging), research tasks. Member self-claim deferred beyond
  this set.

## 11. Open questions / risks

1. **File formats** — Member / Paternal-Ancestry / Y-STR-Overview **confirmed**
   against real exports (§3). Still need samples for **Y-SNP/Big Y, mtDNA, and
   Family Finder** (§3.4) before Phases 3–5.
2. **Build/liftover for SNP & Big Y** — Big Y build **confirmed hg38** (§3.5, read
   from VCF header not the stale readme). Remaining decision: lift hg38→CHM13 on
   import vs. store-as-reported + lazy lift. Affects immediate comparability to our
   CHM13 analyses.
3. **`project_id` retirement** — keep as nullable "home project" or fully migrate
   reads to the join table now? (§4.1)
4. **Fuzzy threshold & signals** — needs a small labeled sample to tune; start
   conservative (high precision, more rows land in the confirm queue).
5. **AppView design doc** — the collaboration tier (§8) is a separate, substantial
   AppView/Postgres design (assertion store, `ResearchSubject` registry, ACL +
   audit, sync API). Phases 6–8 are blocked on it. This doc only fixes the
   Navigator-side contract.
6. **Privacy boundary** — confirmed: **AppView holds NO PII** (anonymized-only
   stance retained); PII never a public PDS record either; the private tier moves
   **encrypted Edge-to-Edge (P2P)**, AppView brokers only (§8.1/8.2), reusing the IBD
   crypto. Watch for leakage in the `fed` record-derivation **and** ensure no PII
   field is ever put in an AppView-bound payload — tables marked
   `PII / never-leaves-encrypted`.
7. **Consent basis** — co-admin-team sharing rests on the existing FTDNA admin
   relationship as the legal/consent basis. Validate against FTDNA TOS + applicable
   data law before Phase 6; the member-claim path (§8.5) is the longer-term fix.
8. **IBD-pipeline coupling** — cross-admin resolution (§8.3) reuses the IBD
   match-suggestion backbone; confirm that pipeline's API/maturity can carry the
   extra same-person/merge-suggestion load before Phase 7.
9. **No batch deep-data export (resolved fact, §3.4)** — FTDNA admins cannot bulk
   export Big Y/mtDNA/FF; it's manual pose-as-member, gated by `access_granted`.
   Phases 3–5 are per-member single-file imports, not batch. Mainly affects UX
   (work-list) and throughput expectations, not feasibility.
10. **Academic profile guarantee** — no-PII is policy+tested-convention, not
    structural (§9). Confirm the null-PII test + "academic" workspace flag are
    sufficient for the target IRB; controlled-access sources (EGA/dbGaP) are out of
    scope for v1 (see companion doc).

## 12. Next step

1. **Importer:** roster, ancestry, Y-STR, **and Big Y** formats are locked (§3,
   §3.5); Phase 0/1/2/3 parser specs are ready to groom into tickets. Still need
   **mtDNA + Family Finder** per-member examples to unblock Phases 4–5.
2. **Unblock collaboration:** spin up a companion **AppView research-layer design
   doc** (assertion store, `ResearchSubject` registry, ACL/audit, sync API) so
   Phases 6–8 have a target. The Navigator-side contract here is the input to it.
3. **Academic profile:** see companion **`docs/design/academic-ena-import.md`**
   (IRB-/publishing-facing); the shared engine (§9) is the same.
```

