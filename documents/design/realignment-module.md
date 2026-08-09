# Realignment module — design & options

Status: **design / specification only** (no code), with a **phase 0 validation spike run
2026-08-08** whose measurements are folded in below. Branch context: `rust-rewrite`.
Scope if built: `navigator-analysis` (revert + post-process), a new `navigator-align` crate
(the mapper — see Decision 1), `navigator-refgenome` (aligner-index cache), `navigator-app`
(job orchestration + provenance), `navigator-store` (alignment provenance migration),
`navigator-ui` (opt-in background job + warnings).

> **The spike invalidated three of this document's load-bearing claims** — the motivating
> ancestry premise, the aligner-backend decision, and the Windows strategy — and put a
> condition on a fourth, the index cost, which holds only if the index batch size (`-I`) is set
> deliberately. Each is corrected in place below and the measurements are in
> [Phase 0 spike results](#phase-0-spike-results-2026-08-08). Read that section before treating
> any decision here as settled.

## Problem

Most vendor whole-genome data arrives aligned to **GRCh37 or GRCh38** (Dante, Nebula,
Sequencing.com, YSEQ, clinical labs). Navigator's modern analyses increasingly want
**CHM13v2 / hs1** (`Build::Chm13v2`):

- **Coverage / callable / SV** computed on GRCh38 describe the *old* reference — they miss
  the ~200 Mbp of sequence CHM13 adds and the ~thousands of collapsed-duplication and gap
  regions T2T resolves.
- **Reads unmapped or mismapped on GRCh38** are invisible to every analysis. The ones that
  belong in CHM13-resolved sequence can only be recovered by re-mapping.
- **Genotypes carry GRCh38 reference bias** and that build's indel representation, whatever
  coordinate system they are later expressed in.

> ### ⚠️ Correction (2026-08-08): ancestry is *not* build-locked
>
> Earlier revisions of this document claimed that "a panel is built for one reference and a
> build mismatch is a hard error (`AppError::AncestryPanelBuildMismatch`)", so that a GRCh38
> sample **cannot run ancestry against the hs1 panel at all**. That was the headline
> justification for this module. **It is false against the current tree**, and the design must
> not be re-argued from it:
>
> - `AppError::AncestryPanelBuildMismatch` (`navigator-app/src/error.rs:28`) is declared and
>   **never constructed anywhere in the workspace**. It is vestigial.
> - `navigator-analysis/src/ibd_panel.rs` is a **multi-build** panel: each site carries
>   `chm13`, `grch38`, and `grch37` loci, mapped offline by allele-aware GATK liftover, with
>   the same biological alleles on each build so the dosage is build-agnostic.
> - `App::ibd_panel_dosages` (`navigator-app/src/import_unified.rs:1305`) branches on build.
>   A CHM13 alignment is genotyped at native loci; a **GRCh37/GRCh38 alignment is genotyped at
>   that build's own coordinates** and the dosages are re-keyed to canonical CHM13.
> - `build_autosomal_profile_inner` (`navigator-app/src/haplogroup.rs:2113`) pools those
>   per-alignment dosages into the autosomal consensus, and `estimate_ancestry_from_consensus`
>   runs the estimators off that consensus with the build pinned to `Chm13v2`.
>
> Off-build vendor samples therefore **already produce ancestry estimates today**, via the
> multi-build panel and the consensus. The architecture routed around the problem this module
> was conceived to solve.
>
> What survives is the read-level case listed above — recovered reads, native hs1
> coverage/callable/SV, and bias-free genotypes. Those are real and liftover cannot deliver
> them. But they are a **narrower and less urgent** payoff than "ancestry is impossible without
> this," and the cost/benefit in [Resource profile](#resource-profile--ux) should be judged
> against them, not against the retracted claim.

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

So for genome-wide, read-level analyses of an off-build vendor sample — native coverage,
callable, and SV on hs1, and any analysis that needs reads GRCh38 lost — the correct path is to
**re-map the reads to CHM13v2**, not to lift. (Ancestry and IBD are *not* in that list; see the
correction above.)

This doc lays out the options for doing that and recommends a path. It is a **heavyweight,
opt-in** capability — hours of compute and tens of GB of scratch per WGS — not part of the
import fast path.

## Goals

1. Take a vendor `Alignment` on GRCh37/GRCh38 and produce a new `Alignment` on `Build::Chm13v2`
   (or `Chm13v2MaskedRcrs`) whose reads are *actually mapped* to that reference.
2. Make the realigned alignment a first-class workspace row so every existing analysis runs on
   it natively (no liftover): coverage, callable, SV, haplogroups — and ancestry/IBD at their
   native loci rather than through the multi-build panel's GRCh38 coordinates.
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
| **Realignment** (this doc) | Re-map all reads to a new reference | hours / tens of GB disk / ~6–12 GB RAM | recover reads GRCh38 lost; native hs1 coverage/callable/SV; bias-free genotypes |
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

## Phase 0 spike results (2026-08-08)

Run on an Apple Silicon macOS box (128 GB RAM, 16 cores — so nothing below was
memory-constrained), against the already-cached `~/.decodingus/references/chm13v2.0.fa`.
Spike crates were standalone, outside the workspace; nothing here touched `Cargo.lock`.

**S1 — the minimap2 FFI builds, but `static` smuggles in htslib.** `minimap2 =
"0.1.31+minimap2.2.30"` with `static` + `simde` compiles clean on arm64 and maps correctly.
However the crate declares `static = [minimap2-sys/static, rust-htslib/static]` with no `?`
sigil, so **enabling `static` activates the optional `rust-htslib` dependency and compiles
`hts-sys`** — silently linking htslib, the C library noodles exists to replace, even with
`default-features = false`. Depending on `minimap2` with `["simde"]` and on `minimap2-sys`
directly with `["static", "simde"]` gets identical behavior with htslib gone (10 packages
instead of ~50; 5.5s instead of 18.3s). Any FFI dependency spelled the obvious way is wrong.

**S2 — index RAM is bounded by `-I`, but only outside the library bindings.** CHM13v2, `sr`.
The first measurements used the *library* builder APIs, which construct one monolithic in-memory
index; the later ones used the CLI's `-d` dump path, which builds a part, writes it, and frees
it — the same thing `minimap2 -d ref.mmi ref.fa` has always done:

| Path | `-I` (parts) | Wall | Peak RSS |
|---|---|---:|---:|
| Library builder (C FFI) | default (1) | 35.2 s | 18.8 GiB |
| Library builder (C FFI), `batch_size` = 1 GB | 1 GB | 23.1 s | 17.2 GiB |
| Library builder (pure Rust) | default (1) | 20.8 s | 18.4 GiB |
| **CLI `-d` dump** | default (1) | 24.8 s | 19.2 GiB |
| **CLI `-d` dump** | 1 GB (4) | 22.2 s | **11.7 GiB** |
| **CLI `-d` dump** | 400 MB | 24.1 s | **8.7 GiB** |
| **CLI `-d` dump**, `-t 4` | 200 MB | 26.4 s | **7.5 GiB** |

Mapping against the resulting index, same reads, `-t 4`:

| Index | Peak RSS |
|---|---:|
| Single-part `.mmi` | 10.25 GiB |
| 4-part `.mmi` (`-I 1G`) + `--split-prefix` | **5.4 GiB** |

**Wall time is flat across every `-I`**, so bounding the memory is close to free. A 16 GB desktop
runs this comfortably at `-I 1G`; 8 GB is plausible at `-I 400M`. The `.mmi` is 8.93 GB on disk
either way.

Why the library numbers mislead: `minimap2-rs` retains every part
(`idx_parts: Vec<Arc<MmIdx>>`), so lowering `batch_size` *through the binding* saved only ~8% and
made the cost look algorithmic. It is not — it is an artifact of asking a binding for one
resident index. **Any implementation here must build and map part-by-part rather than
materialising a whole-genome index in memory.** Note that per-part mapping plus result merging is
CLI-level logic in minimap2 (`--split-prefix`), not something either library API does for us, so
`navigator-align` has to implement it.

**The trade is MAPQ, not placement.** Mapping 5,000 reads against a ~5-part split of chr21 versus
a single-part index: **all 5,045 alignment records had identical target, start, end, and strand**.
Seven records (0.14%) differed, in MAPQ only, always revised *upward* (20→38, 19→29, 0→30, …),
because a read's second-best hit can fall in a different part and go uncounted in `s2`. The
worst case is the dangerous direction — an ambiguous MAPQ-0 read reported as MAPQ 30 — so the
effect is small but not benign, and it may be larger genome-wide than in this single-chromosome
proxy, since homologous repeats are spread across chromosomes and split apart more readily. Pick
`-I` as large as the RAM budget allows rather than as small as possible, and validate MAPQ
distributions at the chosen setting.

**S3 — the ancestry premise is false.** See the correction in [Problem](#problem).

**S4 — a credible pure-Rust minimap2 now exists.** `minimap2-pure-rs 0.5.3` is a faithful
translation of minimap2 v2.31 pinned to an upstream commit, with the full preset set,
paired-end support, CIGAR/SAM/PAF output, `.mmi` I/O, and a `tracehash` feature built for
cross-language stage-by-stage hash comparison against C minimap2. Independent parity test —
5,000 simulated 150 bp reads from CHM13 chr21, 1% substitutions, both strands, `sr` preset,
both aligners configured identically with CIGAR:

| Outcome | Count |
|---|---:|
| Byte-identical output lines | **4,987 / 5,000 (99.74%)** |
| Differ, both report MAPQ 0 (multi-mapping; locus arbitrary by design) | 11 |
| Differ in MAPQ only, at identical coordinates (30↔35, 53↔60) | 2 |
| **Disagreements at MAPQ > 0** | **0** |

Scored against simulated truth both got the **same 4,503 reads right (90.06%)** — the same
set, not merely the same count. (The 10% miss is chr21's acrocentric satellite/rDNA arrays
where reads genuinely multi-map; identical on both sides, which is the point.) Dependency tree
is **zero `-sys` crates, no `cc`, no bindgen, no C-toolchain build script**.

Caveats that must travel with that result: 5,000 reads on one 45 Mb chromosome is not a
whole-genome validation; the two MAPQ-only divergences are small but real and Navigator's
callers filter on MAPQ; and the project self-describes as an **"LLM-mediated faithful
translation"** that is **"experimental"**, is not endorsed by minimap2's author, and
*deliberately replicates upstream bugs* for reproducibility.

## Decision 1 — aligner backend integration

The project's defining constraint: **"no external tools."** Today that is *strictly* true — the
only `std::process::Command` in the workspace is a best-effort browser launch for OAuth. It is
also stricter than "no subprocesses": the workspace has repeatedly paid to stay **MSVC-clean**,
picking `bzip2 = "0.6"` for its pure-Rust backend and `bio = "4"` over C-binding POA/WFA crates
precisely because those "don't build under MSVC" (see `navigator-analysis/Cargo.toml`). Any C
FFI dropped into `navigator-analysis` would break `cargo build` for that whole crate on Windows,
not merely disable a feature.

**Revised decision (2026-08-08): the default backend is `minimap2-pure-rs`, with the C FFI
retained behind an off-by-default feature as a cross-check.** The previous revision chose the
FFI on the grounds that "re-implementing [a mapper] in pure Rust to production accuracy is not
realistic" and that no pure-Rust production mapper existed. Spike S4 falsified that.

| Option | How | Verdict |
|--------|-----|---------|
| **1a. Shell out to a minimap2 binary** | `Command::new("minimap2")`, discover on PATH or a managed tools dir | Rejected — adds a PATH/version dependency and per-OS binary shipping for no benefit over linking the same code in-process |
| **1b. minimap2 via FFI (static)** | Link `libminimap2` in-process via the `minimap2-rs` crate (`minimap2` + `minimap2-sys`), `simde` + `minimap2-sys/static` | **Retained as an off-by-default validation backend.** Works (S1), but needs a C toolchain, is unproven on Windows, carries an `unsafe` surface, and the obvious feature spelling silently links htslib |
| **1c. Pure-Rust aligner** (`minimap2-pure-rs`) | A faithful translation of minimap2 v2.31, pinned to an upstream commit; full preset set, paired-end, CIGAR/SAM/PAF, `.mmi` I/O | **Chosen as default** — 99.74% byte-identical to the C implementation with **zero MAPQ>0 disagreements** and identical accuracy against truth (S4); no `-sys` crates, no `cc`, no build script; indexes CHM13 ~40% faster at the same RAM |

Choosing 1c buys three things at once: **Windows works** with no spike and no graceful
degradation (Decision 3 mostly dissolves), the **single-artifact / no-external-tools posture is
preserved intact** rather than bent, and the htslib trap in S1 never arises. Its cost is the
maturity risk below.

**Put the mapper in a new `navigator-align` crate**, not in `navigator-analysis`. The surface is
small — build/load an index, map reads — so a thin internal trait with two implementations is
cheap. That trait is also the *mitigation* for the maturity risk: it lets a suspicious
realignment be re-run through the C backend and diffed, which is what makes a pure-Rust default
defensible rather than a leap of faith. Isolating it in its own crate also means that if the FFI
feature is ever switched on, only that crate needs a C toolchain.

**Maturity risk, stated plainly.** `minimap2-pure-rs` self-describes as an **"LLM-mediated
faithful translation"** that is **"experimental"** and asks users to "stay vigilant to possible
bugs." It is not endorsed by minimap2's author, and it *deliberately reproduces upstream bugs*
for study-to-study comparability. It is v0.5.x from a single lab. Our own parity evidence (S4)
is 5,000 reads on one 45 Mb chromosome, not a genome. Before phase 4 ships to users, parity must
be re-established at WGS scale on a real donor, and the two MAPQ-only divergences chased down —
MAPQ feeds Navigator's caller filters, so that is the one path by which a silent translation
divergence could reach a persisted result.

Whichever backend is active, one library maps short reads (`sr`) and long reads (`map-hifi` /
`map-ont`). Record the backend *and its version* in `Alignment.aligner` (the `probe.rs` aligner
list already knows `minimap2`) so provenance distinguishes a C-mapped from a Rust-mapped
alignment. Note the FFI builder's `with_threads` is deprecated in favour of `with_index_threads`.

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
- **The remaining business case forgives any margin anyway.** ~~Ancestry/IBD panels are common,
  well-behaved SNPs.~~ *(Retracted 2026-08-08 — ancestry/IBD are not blocked; see the correction
  in [Problem](#problem).)* What is left — recovering reads GRCh38 could not place, and
  describing coverage/callable/SV against the reference that actually resolves those regions —
  is served fine by minimap2's accuracy. Neither payoff turns on a fraction of a percent of
  small-variant recall.

**Why not also offer bwa-mem2 (descoped):** even as an opt-in it isn't worth carrying. Its index
is **not buildable by home users** — constructing the human bwa-mem2 index has been clocked at
**~85 GB RAM**, far beyond a desktop — and the prebuilt index is far too large for the project to
distribute. Advanced users who specifically want bwa-mem2 already have their own pipelines and
don't need Navigator to provide it.

The contrast holds up, with one correction: minimap2's `.mmi` was described here as building "in
minutes within a few GB of RAM." It builds in **~25 seconds**, and RAM depends entirely on the
index batch size — **~19 GB if built monolithically, ~7.5–11.7 GB with `-I` set** (S2). So
"a few GB" was optimistic but the conclusion is intact: unlike bwa-mem2's ~85 GB, this is
tunable to fit a 16 GB desktop at no meaningful time cost. The artifact is 8.93 GB on disk.
See Decision 4.

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

**This decision is largely dissolved by Decision 1.** It existed because minimap2 is **not
officially supported on Windows** (upstream recommends WSL) and the `minimap2-rs` FFI crate
documents testing only on x86_64 and aarch64. Both facts still hold — they are now simply
irrelevant to the default path, because the default backend is pure Rust.

`minimap2-pure-rs` depends only on `clap`, `rayon`, `hashbrown`, `bitflags`, `log`,
`env_logger`, `thiserror`, `flate2`, and `memchr` — **no `-sys` crates, no `cc`, no bindgen, no
build script needing a C toolchain** (S4). Its SIMD is runtime-dispatched with a scalar
fallback, which is why it ran unmodified on Apple Silicon. That is the same "compiles everywhere
Rust compiles" property the rest of the workspace has been protecting.

So realignment targets **Windows, macOS (Intel and Apple Silicon), and Linux (x86_64 and arm64)
alike**, with no per-platform capability gate and no degradation message. The historical options
are kept only to explain why they are no longer needed:

| Option | Status |
|--------|--------|
| **3a. Validate minimap2 FFI on Windows** | Only needed if someone enables the off-by-default FFI validation backend there. Not on the shipping path; no longer a release blocker |
| **3b. POSIX/Apple-only feature (degrade gracefully)** | **Withdrawn** — a pure-Rust backend makes the Windows gap it was papering over disappear |
| **3c. WSL2 delegation** | Withdrawn; nothing left to delegate |
| **3d. Cloud/off-box realignment** | Still rejected for v1 on privacy/single-artifact grounds, unchanged |

The residual platform-shaped constraint is **RAM rather than OS** — but per S2 it is a tuning
parameter, not a wall: sizing `-I` to the machine keeps a 16 GB desktop of any OS inside budget.
Decision 4 covers the sizing and the preflight.

## Decision 4 — aligner index management

Realignment needs a minimap2 index (`.mmi`) of CHM13v2, which **does not exist in the cache today**
(`refgenome` caches FASTA + `.fai` + chains + masks only).

- **Cache layout:** extend the refgenome cache with `<base>/minimap2_index/<build>/…` next to
  `references/` and `liftover/`. The `.mmi` is preset-specific, so key the cache by preset
  (`sr` / `map-hifi` / `map-ont`) as well as build.
- **Build on demand.** Generate the `.mmi` lazily on first realignment against a build and cache
  it, surfaced through the same `ReferenceGateway` progress-callback pattern as
  `resolve_reference`.

- **Size `-I` to the machine — this is the memory control, and it is the whole ballgame.** The
  original "a few GB of RAM" was optimistic but the spirit was right, *provided the index is
  built and used part-by-part*. Build monolithically and CHM13 costs ~19 GB; set `-I` and the
  same index costs 11.7 GB at 1 GB batches or 7.5 GB at 200 MB batches, with **no meaningful
  change in wall time** (S2). Mapping follows the same curve: 10.25 GB against a single-part
  index, 5.4 GB against a 4-part one.

  Pick the batch size from detected physical RAM — largest that fits the budget, not smallest
  that fits the machine, because a coarser split costs MAPQ fidelity (below). A first cut:
  `-I 1G` at ≥16 GB RAM, `-I 400M` at 8 GB. Record the chosen `-I` in the cache key and in the
  alignment's provenance, since it is not a purely cosmetic knob.

- **Implement per-part mapping and merging.** `navigator-align` owns the orchestration, but not
  the hard part: `minimap2-pure-rs` exposes `index::split::{create_split_tmp,
  write_split_query_record, read_split_query_record, merge_split_query_records}`, so the
  cross-part MAPQ recomputation — the piece that is easy to get subtly wrong — is reused rather
  than reimplemented. What we must own is the loop and the output, because the crate's own
  file-level split entry points (`map_file_pe_sam_split` and friends) write to **stdout** and take
  `parts: &[MmIdx]`, i.e. they hold every part resident and so give up the memory bound this whole
  decision exists to buy. `index::reader::IdxReader::read_next` is the part-at-a-time source.
  Building one whole-genome resident index instead is the failure mode that produced this
  document's original wrong estimate; don't reintroduce it.

- **Part boundaries overshoot, and `mini_batch_size` is not the control.** Measured: a part
  accumulates whole sequences until the running total *exceeds* `batch_size`, so parts land at or
  above the batch and a reference only slightly larger than the batch still yields one part.
  `mini_batch_size` does not affect the part count at all. Size progress bars from an upper bound,
  not an equality.

- **Known cost: MAPQ inflation, not misplacement.** Against a ~5-part split, 5,045 alignment
  records kept **identical target/start/end/strand**; 7 (0.14%) differed in MAPQ only, always
  revised upward, because a read's second-best hit can sit in another part and go uncounted in
  `s2`. Worst observed case was a genuinely ambiguous MAPQ-0 read reported at MAPQ 30. Navigator's
  callers filter on MAPQ, so validate the MAPQ distribution at the chosen `-I` before shipping,
  and prefer larger batches when RAM allows.

- **Preflight:** detect physical RAM and pick `-I` from it (or refuse if even the smallest
  sensible batch won't fit) *and* check free disk for the 8.93 GB `.mmi` + scratch + output CRAM.
  Refuse early with a clear message rather than filling the disk or thrashing mid-run.

## Decision 5 — realignment scope

| Scope | What | Verdict |
|-------|------|---------|
| **Whole-genome** | Re-map every read | **Recommended for v1.** The only scope that delivers genome-wide coverage/callable/SV on hs1 and recovers reads across the whole genome; conceptually simple. |
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
- **Preflight checks:** disk free (the 8.93 GB `.mmi` + 2–3× the source size for scratch and
  output), and physical RAM — which selects the index batch size rather than gating the feature
  (Decision 4). Refuse early with guidance.
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
- **Ancestry non-regression** (reframed — ancestry is not blocked today, see
  [Problem](#problem)): a GRCh38 vendor sample's ancestry estimate via the multi-build panel and
  the estimate from its realigned hs1 alignment must **agree**. Realignment must not move a
  result users have already seen; if it does, that is a finding about one of the two paths.
- **Read-recovery sanity:** measure reads that were unmapped on GRCh38 and now map into
  CHM13-resolved regions — the realignment payoff should be visible and non-trivial. With the
  ancestry claim retracted this is now the module's **primary** justification, so quantify it on
  a real donor before phase 4 rather than assuming it.
- **Coverage parity:** genome-wide coverage on the realigned hs1 alignment vs the
  pipeline/native expectation for the same sample.
- **Backend parity at WGS scale:** re-run the S4 comparison (pure Rust vs C FFI) on a whole
  genome, not a single chromosome, and diff placements and MAPQ distributions. The phase 0
  evidence is 5,000 reads on chr21; that is enough to justify the backend choice, not enough to
  ship it. Chase the MAPQ-only divergences specifically.
- **Split-index MAPQ validation:** at the `-I` the preflight actually picks, confirm the MAPQ
  inflation stays at the ~0.1% level seen in S2 and that no MAPQ-0 repeat read is promoted into
  a caller's confident band.
- **Determinism:** fixed thread count + fixed `-I` + fixed backend version → reproducible
  terminal calls.

## Phasing / rollout

0. **Phase 0 spikes — done 2026-08-08.** Backend viability, index cost, and the ancestry premise;
   results above. The one open pre-work item they created is deciding whether the retracted
   ancestry justification leaves enough value to proceed (see [Problem](#problem)).
1. **Revert in pure Rust** (collate-by-name external sort → cleaned paired FASTQ/uBAM, unmapped
   included) + unit tests on small fixtures. This is the hard, backend-agnostic core, and it is
   worth building regardless of how the backend question settles.
2. **`navigator-align` crate**: the mapper trait, `minimap2-pure-rs` as default backend, the C
   FFI behind an off-by-default feature, part-by-part index build/map with `-I` sized from RAM,
   and the `minimap2_index` cache; end-to-end on a small genome/region fixture.
3. **Stage C** (sort, short-read markdup, CRAM emit/index) + Stage D registration with the new
   provenance columns (store migration).
4. **App orchestration + UI**: opt-in cancellable background job, preflight, progress, badges;
   wire realigned alignments into the analysis selectors.
5. **WGS-scale backend parity + MAPQ validation** (replaces the former Windows FFI spike, which
   Decision 1 made unnecessary) — the gate before this is exposed to users.
6. **Y/mt-only mode** as a lightweight scope.

All five desktop platforms are in scope from phase 2 — Windows, macOS (Intel and Apple Silicon),
and Linux (x86_64 and arm64) — because the default backend needs no C toolchain. Phases 1–2 prove
correctness in isolation before any UI; phase 4 is the first user-visible delivery.

## Open questions

- **Does the remaining justification carry the cost?** With ancestry/IBD retracted, the payoff is
  read recovery, native hs1 coverage/callable/SV, and bias-free genotypes, against hours of
  compute and a substantial new subsystem. This is the *first* question now, and it is a product
  call rather than a technical one.
- ~~Does minimap2 / `minimap2-rs` build on Windows?~~ **Moot** — the default backend is pure Rust
  (Decision 1). Only relevant if someone enables the FFI validation backend on Windows.
- **Is `minimap2-pure-rs` trustworthy enough to be the default?** Phase 0 says yes at small scale
  (S4); phase 5 has to say yes at WGS scale. The fallback if it doesn't: flip the default to the
  FFI and accept the Windows gap and htslib linkage after all.
- uBAM vs FASTQ as the revert interchange default (RG fidelity vs simplicity)?
- Mark-duplicates: implement in Rust, or fold into the chosen backend's ecosystem? (With a
  pure-Rust mapper the no-external-tools posture now holds for the *whole* pipeline, so this is
  the last place it could leak.)
- Should realignment target `Chm13v2` or the analysis-tuned `Chm13v2MaskedRcrs` (PAR-masked +
  rCRS) by default? Two facts narrow this: `registry.rs:49` normalizes the masked variant to
  `Chm13v2` for chains/coordinates, so the choice only affects chrM and the asset name; and every
  ancestry/IBD asset currently published under `~/.decodingus/ancestry/` is keyed `chm13v2.0`.
  That points at `Chm13v2` unless someone intends to publish masked-variant panels.

## Sources

Backend evidence gathered for Decision 1 / Decision 3:

- minimap2 paper (accuracy, `sr` preset, ~3–4× speed vs BWA-MEM): <https://academic.oup.com/bioinformatics/article/34/18/3094/4994778>
- "Minimap2 and the future of BWA" (lh3): <https://lh3.github.io/2018/04/02/minimap2-and-the-future-of-bwa>
- UMCCR — BWA-MEM vs minimap2 for WGS variant calling ("pretty safe to migrate to Minimap2"): <https://umccr.org/blog/bwa-mem-vs-minimap2/>
- Nebula Genomics ships minimap2-aligned consumer CRAMs (ecseq, inspecting consumer WGS): <https://www.ecseq.com/blog/2023/Inspecting-Consumer-Whole-Genome-Sequencing-Data>
- `minimap2-rs` crate (FFI; `static`/`simde` features; tested x86_64 + aarch64): <https://github.com/jguhlin/minimap2-rs> · <https://crates.io/crates/minimap2>
- minimap2 Windows status (unsupported upstream; WSL recommended; mingw needs patches): <https://github.com/lh3/minimap2/issues/19>
- `-ax sr` perf edge case on repetitive WGS: <https://github.com/lh3/minimap2/issues/1180>
- `minimap2-pure-rs` — pure-Rust translation of minimap2 v2.31; the chosen default backend
  (Decision 1): <https://github.com/henriksson-lab/minimap2-pure-rs> · <https://crates.io/crates/minimap2-pure-rs>
- `rammap` — a more ambitious pure-Rust aligner the above README points at, if the translation
  approach ever proves untenable: <https://github.com/jwanglab/rammap>

Phase 0 measurements were taken with throwaway spike crates (FFI, pure-Rust, and a simulated-read
parity harness) on Apple Silicon macOS, 128 GB RAM, 16 cores, against the cached
`chm13v2.0.fa`. They are reproducible from the tables in
[Phase 0 spike results](#phase-0-spike-results-2026-08-08); nothing from the spikes was committed
to the workspace.
