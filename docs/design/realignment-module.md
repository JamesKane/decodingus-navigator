# Realignment module — design & options

Status: **design / specification only** (no code). Branch context: `rust-rewrite`.
Scope if built: `navigator-analysis` (revert + align + post-process), `navigator-refgenome`
(aligner-index cache), `navigator-app` (job orchestration + provenance), `navigator-store`
(alignment provenance migration), `navigator-ui` (opt-in background job + warnings).

## Problem

Most vendor whole-genome data arrives aligned to **GRCh37 or GRCh38** (Dante, Nebula,
Sequencing.com, YSEQ, clinical labs). Navigator's modern analyses increasingly want
**CHM13v2 / hs1** (`Build::Chm13v2`):

- **Ancestry (admixture/PCA) is build-locked.** A panel is built for one reference and a
  build mismatch is a hard error (`AppError::AncestryPanelBuildMismatch` — the panel's
  `build` must equal the alignment's `reference_build`). A GRCh38 sample therefore **cannot
  run ancestry against the hs1 panel at all** today.
- **IBD panel genotyping** is likewise defined on a fixed build.
- **Coverage / callable / SV** computed on GRCh38 describe the *old* reference — they miss
  the ~200 Mbp of sequence CHM13 adds and the ~thousands of collapsed-duplication and gap
  regions T2T resolves.

We already have **liftover** (UCSC chains via `ReferenceGateway::load_liftover` /
`lift_positions`) and it is the right tool for *coordinates of a fixed site list* — it is how
Y-haplogroup placement reaches hs1 positions. But liftover **cannot** do what realignment
does:

- It only moves positions that have a 1:1 chain block; novel/!restructured CHM13 sequence has
  **no source coordinate**, so nothing lifts into it.
- It cannot recover **reads that were unmapped or mismapped** on GRCh38 but belong in
  CHM13-resolved regions — those reads simply aren't where a lifted query looks.
- Lifting *genotypes* (rather than re-mapping reads) bakes in GRCh38 reference bias and
  indel-representation artifacts.

So for genome-wide, read-level analyses (ancestry, IBD, native coverage/SV on hs1) of an
off-build vendor sample, the correct path is to **re-map the reads to CHM13v2**, not to lift.

This doc lays out the options for doing that and recommends a path. It is a **heavyweight,
opt-in** capability — hours of compute and tens of GB of scratch per WGS — not part of the
import fast path.

## Goals

1. Take a vendor `Alignment` on GRCh37/GRCh38 and produce a new `Alignment` on `Build::Chm13v2`
   (or `Chm13v2MaskedRcrs`) whose reads are *actually mapped* to that reference.
2. Make the realigned alignment a first-class workspace row so every existing analysis runs on
   it natively (no liftover): ancestry, IBD, coverage, SV, haplogroups.
3. Preserve provenance — the realigned alignment must point back to its source.
4. Be correct for the read technologies Navigator already classifies (Illumina short read;
   PacBio HiFi / ONT long read).
5. Degrade honestly where realignment is unavailable (notably Windows — see platform options).

## Non-goals

- Re-deriving joint panels, trees, or masks (off-box pipeline work).
- Replacing liftover — liftover stays the path for fixed-site coordinate mapping (Y/mt).
- Pangenome/graph realignment — that's the separate, later `pangenome-gam-data-sources.md`
  track; this module is the linear-reference bridge that coexists with and predates it.
- Realigning *within* a build (e.g. GRCh38→GRCh38 to change aligner) — out of scope, though
  the machinery would support it.

## Where this sits relative to existing coordinate strategies

| Strategy | What it does | Cost | When |
|----------|--------------|------|------|
| **Native build** | Analyze on the alignment's own reference | free | sample already on the wanted build |
| **Liftover** (today) | Map a fixed site list's coordinates across builds via chains | ms | Y/mt haplogroups; any fixed-panel coordinate hop with 1:1 chain |
| **Realignment** (this doc) | Re-map all reads to a new reference | hours / tens of GB | unlock build-locked genome-wide analyses (ancestry, IBD) for off-build samples |
| **Pangenome** (future) | Map to a graph; project to linear | external, future | supersedes liftover for hard regions; separate track |

---

## Pipeline overview

```
source Alignment (GRCh37/38 BAM/CRAM)
        │
   (A) revert ──────────►  collate by read name → reset/clean → paired FASTQ (+ singletons)
        │                  [includes unmapped reads + their mates]
        │
   (B) align ───────────►  map FASTQ to CHM13v2 with the build's aligner index
        │
   (C) post-process ────►  sort by coord → mark duplicates (short read) → CRAM + .crai
        │
   (D) register ────────►  new Alignment row (reference_build=chm13v2.0, aligner=…,
                            derived_from=<source id>, content_sha256 computed)
        │
   downstream analyses run natively on hs1 (ancestry, IBD, coverage, SV, haplogroups)
```

### Stage A — read extraction & revert (the hard part)

To re-map correctly we must reconstruct the *original unaligned reads* from the source BAM/CRAM.
This is the GATK `RevertSam` + `SamToFastq` / `samtools collate | fastq` job, in Rust on
noodles. The non-obvious requirements:

- **Collate by read name.** Mates in a coordinate-sorted BAM can be arbitrarily far apart, so
  pairing requires grouping by name. A read-name→record hash is infeasible at WGS scale
  (~10⁹ records). Options: (a) an **external merge sort by name** to scratch, then stream
  mates out in order; (b) two passes with an on-disk index. This sort is the dominant cost of
  Stage A and must be disk-backed and cancellable.
- **Keep primaries only.** Drop `secondary` and `supplementary` records. Primary records from
  mainstream aligners are *soft-clipped* (full SEQ/QUAL retained); supplementaries are
  *hard-clipped* (sequence lost) — so keeping primaries preserves the full read. Flag the rare
  hard-clipped-primary case and either skip or best-effort.
- **Restore orientation.** Reverse-strand alignments store SEQ/QUAL reverse-complemented;
  revert them so FASTQ carries the original read.
- **Restore base qualities.** If an `OQ` tag is present (original qualities before BQSR),
  prefer it.
- **Strip alignment state & tags.** Clear position, CIGAR, MAPQ, mate fields, the
  duplicate/QC-fail flags, and aligner tags (`NM`, `MD`, `AS`, `XS`, …).
- **Include unmapped reads and their mates.** This is a *feature*, not an afterthought — reads
  unmapped on GRCh38 are exactly the ones that may map into CHM13-resolved sequence. They are
  the realignment payoff and must flow into the FASTQ.
- **Output shape.** Paired FASTQ (R1/R2 in sync) + a singletons file; or an unaligned BAM
  (uBAM) preserving read-group metadata. uBAM better preserves `@RG`; FASTQ is simpler for
  every aligner. Recommend uBAM-or-FASTQ behind one writer; default FASTQ.

noodles already provides the readers (`reader.rs` `open_seq`/`records`/`records_lazy`) and a
FASTQ writer; the collate-sort and the revert transform are new.

### Stage B — alignment

Map the reverted reads to the CHM13v2 FASTA resolved by `ReferenceGateway::resolve_reference`.
Aligner choice is dictated by read technology (Stage-by-tech below). This is the stage with the
external-tool decision (Decision 1).

### Stage C — post-processing

- **Sort** by coordinate (external merge sort; noodles bam/cram writers).
- **Mark duplicates** for short-read data (coordinate+orientation+UMI-less grouping,
  samblaster-style). **Skip dup-marking for long reads** (HiFi/ONT) — standard practice.
- **Compress to CRAM** against the CHM13v2 reference (smaller; Navigator already reads CRAM
  with a reference repository) and index (`.crai`).

### Stage D — registration & provenance

Insert a new `Alignment` under the **same `SequenceRun`** as the source (same physical library;
only the mapping changed), with `reference_build = "chm13v2.0"`, `aligner = "<backend>"`,
`bam_path` = the new CRAM, `reference_path` set, and `content_sha256` computed. See the data
model section for the provenance column.

---

## Decision 1 — aligner backend integration

The project's defining constraint: **"no external tools."** Today that is *strictly* true — the
only `std::process::Command` in the workspace is a best-effort browser launch for OAuth. A
production short-read/long-read mapper is ~100k lines of hand-tuned SIMD C/C++; re-implementing
one in pure Rust to production accuracy is not realistic. So this is the one place we must choose
how to break, bend, or preserve that constraint.

**The decision is to settle on a single backend: minimap2, linked in-process via FFI.** The
alternatives below are recorded for why, not as parallel paths we'll ship.

| Option | How | Verdict |
|--------|-----|---------|
| **1a. Shell out to a minimap2 binary** | `Command::new("minimap2")`, discover on PATH or a managed tools dir | Rejected — adds a PATH/version dependency and per-OS binary shipping for no benefit over linking the same code in-process |
| **1b. minimap2 via FFI (static)** | Link `libminimap2` in-process via the `minimap2-rs` crate (`minimap2` + `minimap2-sys`), `static` + `simde` features | **Chosen** — compiled into Navigator; no separate binary, no PATH dependency; `static` embeds the lib, `simde` covers non-x86 (Apple Silicon). Cost: a C toolchain at build time and an `unsafe` FFI surface |
| **1c. Pure-Rust aligner** | Use/extend a Rust mapper (rust-bio primitives; no production WGS mapper exists) | Rejected — not accuracy-competitive for WGS; multi-year effort |

Ship **minimap2 via the `minimap2-rs` FFI crate** with `static` + `simde`, so the mapper is
compiled *into* Navigator — the closest we can get to the single-artifact / no-external-tools
spirit while using a real mapper. One library maps short reads (`sr`) and long reads
(`map-hifi` / `map-ont`); the builder API is
`Aligner::builder().sr()/.map_hifi()/.map_ont().with_threads(n).with_index(fasta, None)` then
`aligner.map(seq, …)`. Record the backend in `Alignment.aligner` (the `probe.rs` aligner list
already knows `minimap2`) so downstream and provenance stay honest.

**Why minimap2 is the only backend we need — accuracy is "good enough," and that's not a hedge:**

- **It's already the production aligner for consumer WGS.** **Nebula Genomics delivers customer
  CRAMs aligned with minimap2** ("similar accuracy for variant detection while providing a
  significant runtime speedup compared to bwa-mem"). The exact off-build vendor data this module
  ingests is, in a large share of cases, *already minimap2 output* upstream.
- **Small-variant accuracy is comparable to the BWA-MEM gold standard.** On real human data the
  `sr` preset shows SNP FN 2.6% vs 2.3% with *fewer* false positives, and near-identical indel
  rates; an independent somatic-WGS comparison concluded "it looks pretty safe to migrate to
  Minimap2" (recall even higher; the one precision wrinkle was a caller/EVS-score interaction,
  not a mapping error).
- **It's generally faster** — typically ~3–4× on >100 bp Illumina reads. (Known edge case:
  `-ax sr` can slow on pathologically repetitive WGS — minimap2 issue #1180 — a perf note, not a
  blocker.)
- **The headline business case forgives any margin anyway.** Ancestry/IBD panels are common,
  well-behaved SNPs; minimap2-`sr` is comfortably adequate there. This is the 80% case, and it's
  the *only* case we're targeting.

**Why not also offer bwa-mem2 (descoped):** even as an opt-in it isn't worth carrying. Its index
is **not buildable by home users** — constructing the human bwa-mem2 index has been clocked at
**~85 GB RAM**, far beyond a desktop — and the prebuilt index is far too large for the project to
distribute. That makes it intractable for the audience this module serves. minimap2's `.mmi`, by
contrast, builds in minutes within a few GB of RAM from the FASTA we already cache. Advanced users
who specifically want bwa-mem2 already have their own pipelines and don't need Navigator to
provide it.

## Decision 2 — aligner by read technology

Navigator already classifies technology (`SequenceRun.platform_name` from `@RG PL`,
`SequenceRun.test_type` ∈ {`WGS`, `WGS_HIFI`, `WGS_NANOPORE`, `WES`, `BIG_Y_700`, …}, with
`testtype.rs` inferring HiFi from PacBio/long mean read length and Nanopore from ONT). Map that
to mapper presets:

One backend, one preset switch:

| Source technology | minimap2 preset |
|-------------------|-----------------|
| Illumina / short-read WGS/WES | `sr` |
| PacBio HiFi (`WGS_HIFI`) | `map-hifi` |
| Oxford Nanopore (`WGS_NANOPORE`) | `map-ont` |
| Targeted Y/mt panels | match the panel's underlying chemistry (usually `sr`) |

Pick the preset from `test_type`/`platform_name`; let the user override. Refuse (or warn loudly)
on mixed/unknown technology rather than guessing.

## Decision 3 — cross-platform strategy (the Windows gap)

minimap2 is **not officially supported on Windows** (upstream recommends WSL), and the
`minimap2-rs` FFI crate documents testing only on **x86_64 and aarch64** with no Windows/macOS
coverage ("other platforms may work"). So minimap2 is *not* a guaranteed Windows win; it must be
treated as "unproven on Windows until we validate the build."

What minimap2 FFI **does** cleanly cover, via `static` (embed the lib) + `simde` (portable SIMD):

- **Linux** x86_64 and arm64
- **macOS** Intel *and* Apple Silicon (the `simde` feature is what makes arm64 work)

That is the realignment audience's 80%+. Windows is the residual gap. Options for it:

| Option | Behavior on Windows | Trade-off |
|--------|---------------------|-----------|
| **3a. Validate minimap2 FFI on Windows** | If the C lib + `minimap2-rs` build under MSVC/mingw, realignment works natively | Upstream is unsupported on Windows and the crate is untested there; needs real validation and possibly source patches (cf. the nanoporetech `msvc14` fork) — don't assume it |
| **3b. POSIX/Apple-only feature (degrade gracefully)** | Realignment disabled on Windows; UI explains the sample stays on its native build and can't run hs1-locked analyses (ancestry/IBD) | Simple and honest; Windows users lose off-build realignment until 3a lands |
| **3c. WSL2 delegation** | Run the minimap2 backend inside WSL2 | Heavy install burden; brittle path translation; poor desktop UX — only as an explicit advanced path |
| **3d. Cloud/off-box realignment** | Hand the job to a remote service | Defeats the local-privacy premise; only as an explicit, consented option |

**Recommendation: ship 3b first (macOS + Linux, the 80% case), pursue 3a as a follow-up.** Treat
realignment as a POSIX/Apple capability for v1 and have Windows degrade with a clear "realignment
isn't available on Windows yet — this sample stays on GRCh38" message rather than a mid-job
failure. In parallel, spike the minimap2 FFI build on Windows (3a); if it validates, Windows
joins the native path with no design change. Offer WSL guidance (3c) only as documentation for
advanced users; avoid 3d for v1 — it breaks the single-artifact / local-privacy posture.

## Decision 4 — aligner index management

Realignment needs a minimap2 index (`.mmi`) of CHM13v2, which **does not exist in the cache today**
(`refgenome` caches FASTA + `.fai` + chains + masks only).

- **Cache layout:** extend the refgenome cache with `<base>/minimap2_index/<build>/…` next to
  `references/` and `liftover/`.
- **Build on demand, no download needed.** A minimap2 `.mmi` builds in minutes within a few GB of
  RAM from the FASTA we already resolve — well within a desktop's means — so we generate it
  lazily on first realignment against a build and cache it, surfaced through the same
  `ReferenceGateway` progress-callback pattern as `resolve_reference`. (This cheap, home-buildable
  index is a core reason the project can ship realignment at all; see Decision 1 on why a
  ~85 GB-RAM bwa-mem2 index was descoped.) The `.mmi` is preset-specific, so key the cache by
  preset (`sr` / `map-hifi` / `map-ont`) as well as build.
- **Preflight:** check free disk for index + scratch + output CRAM before starting; refuse with
  a clear message rather than filling the disk mid-run.

## Decision 5 — realignment scope

| Scope | What | Verdict |
|-------|------|---------|
| **Whole-genome** | Re-map every read | **Recommended for v1.** The only scope that unlocks genome-wide ancestry/IBD/coverage on hs1; conceptually simple. |
| **Y/mt-only** | Re-map only chrY+chrM reads (+ unmapped) | **Useful add-on.** Cheap; good for Big-Y / mt-only products and structurally complex Y. Small enough to run anywhere. |
| **Targeted (panel ± flanks)** | Extract reads near lifted AIM/IBD sites and realign just those | **Not recommended.** AIM/IBD panels are genome-wide (tens of thousands of sites across all chromosomes), so "targeted" still touches most reads while adding edge-effect risk and missing reads that *moved* between builds — the very signal realignment exists to recover. |

Recommend: ship **whole-genome** as the core, with **Y/mt-only** as a lightweight mode for
uniparental products and for users who can't afford a full WGS realignment.

---

## Data model & provenance

The current `Alignment` has **no parent/derived-from field** — alignments are independent rows
keyed only to a `SequenceRun`. A realigned alignment must record its lineage:

```rust
// proposed addition (navigator-domain::workspace::Alignment) + store migration
pub derived_from_alignment_id: Option<i64>,  // source alignment this was realigned from
pub derivation: Option<String>,              // e.g. "realign:minimap2-sr" | "realign:minimap2-map-hifi"
```

- New nullable columns via a `navigator-store` migration (follows the `0017`/`0018` numbering);
  existing rows default to `NULL` (= not derived).
- The realigned row sits under the **same `SequenceRun`** (same library), with
  `reference_build = "chm13v2.0"`, `aligner` = the backend, `derived_from_alignment_id` = source,
  fresh `content_sha256`.
- UI: badge realigned alignments ("realigned to hs1 from GRCh38 alignment #N"), and let the
  user pick which alignment a given analysis runs against. Never silently delete the source.

This mirrors how the sidecar fast path stayed *additive* — realignment **adds** an alignment;
it never mutates or replaces the vendor's original.

## Resource profile & UX

- **Heavy and long:** WGS realignment is hours of CPU and tens of GB of scratch (revert sort +
  index + sorted output). It must be an **opt-in, cancellable background job**, modeled on the
  existing streaming deep-analyze / `RunFullAnalysis` spawn-loop (per-sample progress events,
  honor the `CancelAnalysis` flag, `await` between stages so the UI stays responsive).
- **Preflight checks:** disk free (index + 2–3× the source size for scratch + output), and a
  RAM estimate for the chosen backend; refuse early with guidance.
- **Threads:** reuse/extend `NAVIGATOR_ANALYSIS_THREADS` (or a dedicated
  `NAVIGATOR_REALIGN_THREADS`) for both the sort and the mapper.
- **Scratch:** a managed temp dir under the cache, cleaned on success/cancel; resumable stages
  are a nice-to-have, not v1.
- **Triggering:** an explicit "Realign to hs1" action on an off-build alignment (and a batch
  "realign project" that queues sequentially). Never auto-queued at import — same discipline as
  the deep pass.

## Read-technology edge cases

- **Long reads (HiFi/ONT):** no duplicate marking; minimap2 long-read presets; expect larger
  per-read compute but far fewer reads. CRAM compression of long reads is fine.
- **Hard-clipped primaries:** rare but real (some pipelines emit them); detect and skip-or-warn
  rather than emit truncated reads.
- **Read groups:** preserve `@RG` across the revert (uBAM path preserves it best); the realigned
  header should carry the original RGs plus a new `@PG` for the realignment step.
- **Already-CHM13 input:** no-op / refuse — there's nothing to realign (offer a same-build
  re-map only behind an explicit flag).

## Validation plan

- **Concordance vs liftover where both apply:** Y/mt terminal haplogroup from a realigned hs1
  alignment must match the liftover-based call on the validated donor (GFX → R-FGC29071 +
  U5a1b1g) and the HG00096 fast-path result (R1b1a1b1a1a). Realignment must not regress the
  uniparental calls.
- **Ancestry enablement:** a GRCh38 vendor sample that *cannot* run ancestry today should, after
  realignment, produce a stable estimate consistent with a truth/expectation (e.g. the
  EUR-100% GFX donor).
- **Read-recovery sanity:** measure reads that were unmapped on GRCh38 and now map into
  CHM13-resolved regions — the realignment payoff should be visible and non-trivial.
- **Coverage parity:** genome-wide coverage on the realigned hs1 alignment vs the
  pipeline/native expectation for the same sample.
- **Determinism:** fixed thread count + fixed backend version → reproducible terminal calls.

## Phasing / rollout

1. **Revert in pure Rust** (collate-by-name external sort → cleaned paired FASTQ/uBAM, unmapped
   included) + unit tests on small fixtures. This is the hard, backend-agnostic core.
2. **minimap2 FFI backend** (`sr` + long presets) + the `minimap2_index` cache; end-to-end on a
   small genome/region fixture.
3. **Stage C** (sort, short-read markdup, CRAM emit/index) + Stage D registration with the new
   provenance columns (store migration).
4. **App orchestration + UI**: opt-in cancellable background job, preflight, progress, badges;
   wire realigned alignments into the analysis selectors.
5. **Windows FFI spike** — confirm whether minimap2/`minimap2-rs` builds there (3a) or Windows
   stays a graceful POSIX/Apple-only degradation (3b).
6. **Y/mt-only mode** as a lightweight scope.

Target macOS + Linux (incl. Apple Silicon, via the `simde` feature) for phases 1–4 — that's the
80% case and ships realignment to most users; Windows is handled by graceful degradation until
the spike in phase 5 resolves it. Phases 1–2 prove correctness in isolation before any UI;
phase 4 is the first user-visible delivery (ancestry unlocked for off-build samples).

## Open questions

- Does minimap2 / `minimap2-rs` build on Windows (MSVC or mingw)? Upstream is unsupported there
  and the crate is untested on Windows, so this is a **validation spike**, not an assumption — it
  decides whether Windows gets native realignment (3a) or graceful degradation (3b) for v1.
- uBAM vs FASTQ as the revert interchange default (RG fidelity vs simplicity)?
- Mark-duplicates: implement in Rust, or fold into the chosen backend's ecosystem? (Affects
  whether the no-external-tools posture holds for the *whole* pipeline or just the mapper.)
- Should realignment target `Chm13v2` or the analysis-tuned `Chm13v2MaskedRcrs` (PAR-masked +
  rCRS) by default? The masked variant is the recommended short-read calling reference and
  shares CHM13 chains, but ancestry/IBD panel builds must match whichever we pick.

## Sources

Backend evidence gathered for Decision 1 / Decision 3:

- minimap2 paper (accuracy, `sr` preset, ~3–4× speed vs BWA-MEM): <https://academic.oup.com/bioinformatics/article/34/18/3094/4994778>
- "Minimap2 and the future of BWA" (lh3): <https://lh3.github.io/2018/04/02/minimap2-and-the-future-of-bwa>
- UMCCR — BWA-MEM vs minimap2 for WGS variant calling ("pretty safe to migrate to Minimap2"): <https://umccr.org/blog/bwa-mem-vs-minimap2/>
- Nebula Genomics ships minimap2-aligned consumer CRAMs (ecseq, inspecting consumer WGS): <https://www.ecseq.com/blog/2023/Inspecting-Consumer-Whole-Genome-Sequencing-Data>
- `minimap2-rs` crate (FFI; `static`/`simde` features; tested x86_64 + aarch64): <https://github.com/jguhlin/minimap2-rs> · <https://crates.io/crates/minimap2>
- minimap2 Windows status (unsupported upstream; WSL recommended; mingw needs patches): <https://github.com/lh3/minimap2/issues/19>
- `-ax sr` perf edge case on repetitive WGS: <https://github.com/lh3/minimap2/issues/1180>
