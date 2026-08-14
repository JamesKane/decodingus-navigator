# Realignment module — design & options

Status: **built and validated — phases 1–5 complete (2026-08-14).** Branch: `worktree-realignment`,
not yet merged. The measurements from the phase 0 spike (2026-08-08), the backend parity work
(2026-08-10), the first whole-genome run (2026-08-12), the second (2026-08-13) and the third — which
passed every acceptance criterion — are folded in below. See
[Phase 5 result](#phase-5-result--wgs229-end-to-end-2026-08-14).
Scope if built: `navigator-analysis` (revert + post-process), a new `navigator-align` crate
(the mapper — see Decision 1), `navigator-refgenome` (aligner-index cache), `navigator-app`
(job orchestration + provenance), `navigator-store` (alignment provenance migration),
`navigator-ui` (opt-in background job + warnings).

> **The purpose of this module is Y-chromosome variant discovery** (2026-08-10). Earlier revisions
> argued it from autosomal ancestry, and an intermediate one from generic read recovery; both were
> wrong, and Decision 5's scope ranking is being re-argued as a result. See [Problem](#problem).
>
> Separately, the phase 0 spike invalidated three technical claims — the aligner-backend decision
> and the Windows strategy — and put a condition on the index cost, which holds only if the index
> batch size (`-I`) is set deliberately. Each is corrected in place below and the measurements are
> in [Phase 0 spike results](#phase-0-spike-results-2026-08-08).

## Problem

**This module exists for Y-chromosome variant discovery.** Everything else it enables is a
side benefit; if Y discovery did not need it, it would not be worth building.

Most vendor whole-genome data arrives aligned to **GRCh37 or GRCh38** (Dante, Nebula,
Sequencing.com, YSEQ, clinical labs). Discovering a sample's private and novel Y-SNPs — placing
it on the tree, and producing claims an AppView curator can accept — requires it to be on
**CHM13v2 / hs1** instead, for three reasons that compound:

- **GRCh38's chrY is not a usable discovery substrate.** It is gap-ridden and partly collapsed,
  while CHM13v2's Y is a complete T2T assembly (HG002's). Discovery means finding variants at
  positions nobody enumerated in advance, so the parts of the Y a reference *fails to represent*
  are exactly the parts that never yield a call.
- **Every curated Y asset is CHM13-defined.** The cohort callable mask
  (`chrY.callable_mask.chm13v2.bed`), the non-PAR restriction (`chrY_nonPAR.chm13v2.bed`), the
  recurrent-site blocklist, and the de-novo tree itself all come from a pipeline that
  joint-genotyped ~3,352 **CHM13-aligned** males (see
  [`private-y-variant-filtering.md`](private-y-variant-filtering.md)). A GRCh38-aligned sample
  cannot be filtered by any of them, and the private-Y doc's own evidence is that *without* that
  filter stack a WGS sample leaks hundreds of fake privates per tip — so an unfiltered call set
  is not a weaker result, it is an unusable one.
- **Liftover cannot substitute, specifically here.** It moves a *fixed site list* across builds,
  which is why Y-haplogroup placement already uses it. Discovery is the opposite problem: the
  variant is not in any list yet, and sequence T2T resolves that GRCh38 lacks has no source
  coordinate to lift *from*. Reads that mismapped in GRCh38's collapsed ampliconic and
  palindromic regions are likewise not where a lifted query looks.

Two consequences follow that are worth stating plainly, because they are easy to over-claim:

- Realignment puts a sample **into the coordinate system where the filters and the tree live**.
  It does not by itself solve reference polarity — CHM13's Y is haplogroup J, so the reference
  base is the *derived* allele at many Y-SNP sites (`registry.rs::reference_polarity`), and an
  R sample still needs the carrier filter to avoid reporting the whole J-vs-R divergence.
- The same argument applies to **mtDNA**, which travels with the Y in every uniparental product,
  and only weakly to the autosomes.

Genome-wide benefits do come along — native coverage/callable/SV on hs1, and reads GRCh38 could
not place — but they are consequences, not the reason.

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
> multi-build panel and the consensus.
>
> Correcting the correction (2026-08-10): the retraction above is accurate, but an intermediate
> revision then re-argued the module from generic read recovery, which was also wrong. **Fixed-site
> autosomal matching was never the point.** The module is for Y-chromosome variant discovery — see
> [Problem](#problem). Ancestry is mentioned here only so nobody reinstates it as a justification.

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

1. **Make an off-build sample eligible for Y variant discovery** — on the reference the callable
   mask, the non-PAR restriction, the recurrent blocklist, and the de-novo tree are all defined
   on. This is the goal the module is measured against.
2. Take a vendor `Alignment` on GRCh37/GRCh38 and produce a new `Alignment` on `Build::Chm13v2`
   (or `Chm13v2MaskedRcrs`) whose reads are *actually mapped* to that reference.
3. Make the realigned alignment a first-class workspace row so every existing analysis runs on
   it natively (no liftover): private-Y and haplogroup placement first, then coverage, callable,
   and SV.
4. Preserve provenance — the realigned alignment must point back to its source.
5. Be correct for the read technologies Navigator already classifies (Illumina short read;
   PacBio HiFi / ONT long read).
6. Degrade honestly where realignment is unavailable.

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
   (C) post-process ────►  sort by coord → mark duplicates (short read) → BAM + .bai
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

- **Sort** by coordinate (external merge sort; noodles bam/cram writers). Runs spill as ordinary
  BAM files rather than a bespoke encoding — an alignment record carries a CIGAR, a tag
  dictionary, and per-base sequence and quality, so inventing a serialization means re-deriving
  BAM's encoding badly. Unplaced reads sort last, which matters more here than usual: recovering
  reads the old reference could not place is much of why this module exists.
- **Mark duplicates** for short-read data (coordinate+orientation+UMI-less grouping,
  samblaster-style). **Skip dup-marking for long reads** (HiFi/ONT) — standard practice.
  Group on the **unclipped** 5' position: two copies of one molecule can be soft-clipped
  differently, which moves the alignment start without moving the fragment, so grouping on the
  alignment start looks right and silently misses them. Both ends of a template must receive the
  same verdict — one end marked and the other not shows consumers half a pair — which holds
  because the signature is symmetric and the sort is deterministic.
- **Finalise and index.** Duplicate marking has already written exactly the bytes that belong in
  the output, so this is a rename plus a `.bai` rather than another compression pass.

  This was **CRAM emit** — reference-compressed, ~17 GB against 60–80 GB for a 30x WGS, which is
  the better container on paper. Two defects in `noodles-cram` 0.94 make it unreachable, and both
  appear only at real scale: writing panics on a secondary alignment with `SEQ: *` (legal SAM, and
  what minimap2 emits), and `cram::fs::index` then panics on any *multi-reference* slice because it
  decodes records against an empty `fasta::Repository` — `// TODO` still in the upstream source.
  With 25 contigs every slice straddling a contig boundary is multi-reference, so a whole-genome
  CRAM cannot be indexed at all. A desktop user holds a handful of whole genomes; disk is the
  cheaper side of that trade. The CRAM path is kept and still tested — the defects are upstream and
  fixable, and this should be revisitable without rebuilding the stage.

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

## Phase 5 results — backend parity on real reads (2026-08-10)

**Verdict: the pure-Rust backend does not pass the gate on this evidence.** It is not wrong about
*where* reads go; it is systematically more confident about it, concentrated on chrY, in a way
that would change the callable set the private-Y filter stack depends on.

### What was run

168,184 real reads from **WGS229** (`WGS229.bwa-mem2.b38.cram`, GRCh38, bwa-mem2) — the donor
[`private-y-variant-filtering.md`](private-y-variant-filtering.md) measures — extracted from
`chrY:2,800,000-6,000,000` through the shipped revert stage, then mapped against CHM13v2 with the
`sr` preset by both `minimap2-pure-rs` and C minimap2 via FFI. Single-end, R1 only.

Stage A on real vendor data for the first time: 338,511 records in, 168,184 pairs, 1,304
singletons, 142 unmapped carried through, 839 secondary and 0 supplementary dropped, 0
hard-clipped primaries.

### Agreement

| | reads | share |
|---|---:|---:|
| Identical | 166,837 | 99.20% |
| Differ in MAPQ only (same locus) | 1,204 | 0.72% |
| Differ in locus, both MAPQ 0 | 92 | 0.055% |
| **Differ in locus, MAPQ > 0** | **51** | **0.030%** |

Neither implementation mapped a read the other could not. Phase 0's simulated comparison found
*zero* MAPQ>0 disagreements; real reads find 51.

### The finding: MAPQ is systematically inflated, on chrY

Of the 1,204 same-locus MAPQ differences, **1,157 have the Rust value higher and only 47 lower** —
median `+10`, mean `+8.55`. That is a one-directional bias, not noise.

It is concentrated where the module lives. 89.0% of these reads map to chrY, but the
MAPQ-difference *rate* is **0.80% on chrY against 0.053% elsewhere — a 15x enrichment**, so the
concentration is real rather than an artifact of having sampled a chrY region.

**111 reads cross the MQ>=20 callable threshold upward** (3 downward). The cohort callable mask is
defined as `depth >= 4, MQ >= 20`, so those are reads that would newly enter the callable set on
the strength of a confidence the reference implementation does not share. Admitting artifact reads
to the callable set is the documented mechanism of the fake-private flood this module's whole
filter stack exists to prevent — the private-Y doc measures it at hundreds of spurious privates
per tip when filtering is inadequate.

A plausible mechanism, consistent with both halves of the data: the two differ in how the
second-best alignment is scored. Where there is no real competitor the Rust version rates the hit
higher; where there is a strong one it finds it and the C version does not. chrY's ampliconic and
palindromic structure produces near-equal competitors constantly, which is why the divergence
concentrates there.

### The second finding: chrX/chrY swaps in the X-transposed region

Of the 51 confident placement disagreements, 32 involve chrY and **17 are chrY<->chrX swaps**, all
in the X-transposed region (chrY ~2.8-5.7 Mb against chrX ~88-91 Mb, where the two are ~99%
identical). The pattern runs both ways, and in each case one implementation is confident (MAPQ up
to 60) while the other reports MAPQ 1. For Y variant discovery these are the reads whose assignment
decides whether a Y variant has support at all.

### Paired-end: the gap widens

The single-end result left open whether pairing would close it — mate rescue and the proper-pair
bonus both feed MAPQ. It does not. Re-run as `-ax sr` against **upstream minimap2 2.31-r1302**, the
exact release `minimap2-pure-rs` translates, on the same reverted FASTQ pair:

| | single-end | paired-end |
|---|---:|---:|
| Records compared | 168,184 | 336,368 |
| Identical (locus + MAPQ) | 99.20% | **98.59%** |
| MAPQ differs at same locus | 0.72% | **1.09%** |
| Rust higher / lower | 1,157 / 47 | **3,545 / 124** |
| Median delta | +10 | +7 |
| Crosses MQ>=20 upward | 111 (0.066%) | 219 (0.065%) |
| chrY<->chrX swaps | 17 | 53 |
| Different locus | 143 | 110 |

Pairing made per-record agreement *worse* — 1.09% of records differ in MAPQ against 0.72%
single-end — while leaving the callable-threshold crossing rate essentially unchanged at ~0.065%.
The direction is unchanged and remains overwhelming: 96.6% of MAPQ differences are the Rust value
being higher. Neither implementation produced a primary record the other lacked, and the
proper-pair flag disagreed on only 12 records, so the pairing logic itself is not the problem —
the MAPQ arithmetic underneath it is.

This also **does not reproduce the upstream crate's own claim**. Its README reports exact PAF
parity against C minimap2 for `sr` on "HG002 WGS 1M paired reads". On real chrY reads at the same
version, with the same invocation, we do not see it. The likeliest explanation is that a
genome-wide benchmark is dominated by unique sequence, where the two agree, and chrY's ampliconic
and palindromic structure is where they do not — which is precisely the sequence this module
exists to map.

### What this does and does not establish

- It is **one donor and one 3.2 Mb window**. The direction and concentration are consistent enough
  to act on, but the magnitude is not yet a genome-wide figure.
- It says nothing against *placement* accuracy: 98.6% fully identical, no primary record present
  in one and missing from the other, and confident locus disagreements are 0.03%.
- The failure is specific and bounded: **MAPQ, on repeat-structured sequence**. That is unlucky,
  because it is the one number the Y filter stack gates on and the one chromosome it gates.

### The measurement that decides it: do the calls change?

The MAPQ divergence only matters if it reaches the output. It does not.

Both SAMs — the same 168k read pairs, mapped by upstream 2.31 and by
`minimap2-pure-rs` — were put through the shipped stage C (sort, mark duplicates) and the shipped
de-novo caller over `chrY:2,800,000-6,000,000`. That window is deliberately the unfavourable one:
it contains the X-transposed region and carries the highest density of mapping disagreement found
anywhere in the parity runs.

| | calls |
|---|---:|
| Upstream minimap2 | 232 |
| `minimap2-pure-rs` | 233 |
| **Shared** | **232** |
| Only upstream | **0** |
| Only pure-Rust | 1 |
| Shared, differing depth | 6 |
| Shared, differing allele fraction | 1 |

**Every call upstream makes, the Rust backend also makes.** The one extra call is
`chrY:5,120,338 T>G` at depth 3 — below the cohort callable mask's `depth >= 4` gate, so the
filter stack the private-Y path already applies removes it. Duplicate marking was identical on
both (19,101 records), so stage C is stable across the backends too.

So the divergence is real, measurable, and confined to a quantity the downstream thresholds
absorb. A MAPQ that is systematically ~7 points high does not move a call across a `MQ >= 20`
gate often enough to matter once `depth >= 4` is also required — the reads it promotes are the
thin, ambiguous ones that the depth gate was there to catch regardless.

### Verdict

**Ship the pure-Rust backend, with the divergence documented.** This is Decision 1 option 4, chosen
because options 2 and 3 both fail the project's constraints — a desktop application for
non-specialists cannot require minimap2 on `PATH`, and linking the C library reintroduces the
toolchain, `unsafe` surface, htslib trap, and unproven Windows build that Decision 1 rejected.

What this rests on, stated so a later reader can judge it:

- One donor (WGS229), one 3.2 Mb window, ~168k read pairs. The window was chosen to be the worst
  case, not a representative one, which strengthens a null result but does not make it genome-wide.
- SNV calls under default `HaploidCallerParams`. Indels were not separately compared.
- The six shared calls with differing depth are the visible edge of the effect: the evidence behind
  a call can shift even when the call does not. At a different depth or coverage that margin could
  tip.

Re-run this comparison if the backend is upgraded, the caller's thresholds change, or the callable
mask is rebuilt — those are the three things that could turn an absorbed divergence into a visible
one.

### Consequence

**Decision 1 was reopened and re-settled on the same backend**, for a better reason than the one
that originally chose it. It chose the pure-Rust backend on phase 0 evidence of "zero
MAPQ>0 disagreements", measured on simulated reads over a pseudo-random reference — sequence with
no repeat structure, and therefore no second-best alignments to disagree about. That evidence does
not survive contact with chrY. The options:

1. **Find and fix the divergence.** `minimap2-pure-rs` ships a `tracehash` feature built for
   exactly this: stage-by-stage hash comparison against C minimap2, enabled on both sides. Its
   existence means upstream anticipated translation divergences and built the tool to localise
   them. Best outcome if it works — the pure-Rust backend is what makes this module buildable on
   every platform without a C toolchain — but it is open-ended work on someone else's code.
2. **Shell out to an installed minimap2** (Decision 1a, previously rejected). Correct by
   construction, and the reference implementation is already on this machine. Costs a PATH
   dependency, per-OS binary shipping, and the single-artifact property the project values.
3. **Link the C library** (Decision 1b, previously rejected). Correct by construction and no PATH
   dependency, but reintroduces a C toolchain, an `unsafe` surface, the htslib trap, and an
   unproven Windows build.
4. **Ship anyway with the divergence documented.** ✅ **Chosen** — the experiment was run and the
   divergence does not propagate to calls. See "the measurement that decides it" above.

## The first WGS-scale run (2026-08-12)

WGS229's GRCh38 CRAM (17.3 GB, 615.6M reads) through stages A–D on 16 cores. It **did not
complete** — it died in stage 7 of 8 after 10.7 hours — but it is the first measurement of the
pipeline at the scale it exists for, and the failure is more informative than the timings.

| stage | wall clock | share |
|---|---:|---:|
| Revert | 1 h 05 m | 10% |
| Index | 0 s (cached) | — |
| Map | 3 h 40 m | 34% |
| Sort | **4 h 44 m** | **44%** |
| Mark duplicates | 1 h 09 m | 11% |
| CRAM | died | — |
| **Total** | **10 h 41 m** | |

### What it found

**A secondary alignment has no SEQ, and CRAM cannot encode it.** The panic — `range end index 65
out of range for slice of length 0`, from inside noodles — is the *read's* sequence being empty,
not the reference's. CRAM stores a read as its difference from the reference, so encoding walks the
CIGAR comparing bases; a record with an aligned CIGAR and `SEQ: *` has nothing to compare, and
noodles indexes the empty slice rather than checking. This is legal SAM — a secondary alignment may
omit SEQ, and minimap2 uses that permission, since only the primary carries the bases — so the
pipeline was always going to produce it. Nothing smaller than a real mapping run would have shown
it: every fixture in the suite built records with sequences. Non-primary records of that shape are
now dropped and counted; the same shape on a primary is a read going missing and fails loudly.

**The sort is the pipeline's most expensive stage**, at 44% of wall clock — more than mapping,
which is the stage everyone expects to dominate. Worth profiling before this is exposed to users.

**Stage A is single-threaded, and on a CRAM source that shows.** `reader::open_seq` gives the BAM
path a multithreaded bgzf reader and the CRAM path a plain one, so the revert decodes 17 GB of
reference-compressed data on one core while fifteen sit idle. The gzip FASTQ write is serial too
(already at `Compression::fast()`). The fix is not to thread the decoder but to **revert
per-contig in parallel** — the codebase already does this elsewhere via `reader::decode_pool`, a
coordinate-sorted CRAM with a `.crai` supports it, and collation is a sort, so per-worker output
order does not matter: each worker produces sorted runs the existing k-way merge already consumes.
Unmapped reads need one extra pass over the `*` bin. That should take stage A from ~65 minutes to
roughly ten.

### Where the time actually went (2026-08-13)

Three hypotheses were tried against the mapping stage before one was measured, which is the part
worth recording: the stage alternates read → map → write per batch, the CPU-load graph shows a
valley between peaks, and *both* plausible readings of that valley were wrong.

Profiling was blocked by the build, not by the code. `[profile.release]` sets `strip = true`, so
every sample resolves to `???`. A `profiling` profile now inherits release and keeps symbols, and
`navigator-align/examples/map_profile` runs stage B alone against an already-reverted FASTQ pair so
the profile attributes the stall to a function rather than to a stage. With that, the main thread:

| | share of the serial phase |
|---|---:|
| zlib deflate (`longest_match` alone a third) | **~60%** |
| malloc | 10.7% |
| BAM record encoding | ~5.7% |
| inflate + FASTQ parse | ~3.8% |

**The valley is the output BAM being compressed on one thread.** Not read delivery, which was the
standing assumption and is under 4%. Giving the mapper's writer the same worker pool stage C got
took the stage from **21,500 to 26,500 reads/s (+23%)** on identical input, and removed deflate
from the main-thread profile entirely.

What remains serial there is ~15% BAM record encoding and ~9% inflate, plus malloc churn from the
SAM-text→`RecordBuf` round trip that [`output`](../../crates/navigator-align/src/output.rs)
documents as a deliberate trade. Moving encoding into the pool and prefetching the next batch on
its own thread are the next two steps, in that order.

One reading correction for anyone profiling this again: main-thread samples inside minimap2 are not
overhead. `rayon::ThreadPool::install` makes the calling thread a pool worker, so the main thread
does its share of the mapping and shows minimap2 frames while doing it.

### Resumability is not a nice-to-have

[Resource profile & UX](#resource-profile--ux) records "resumable stages are a nice-to-have, not
v1". This run is the counter-evidence: a bug in stage 7 discarded the seven stages that worked,
because [`JobScratch`] deletes every intermediate however the job ends. Ten hours of correct
compute was thrown away to reclaim disk from a job that had already failed.

The scratch guard's default is still right for a shipped desktop application — hundreds of GB must
not outlive a failure the user did not cause — so the immediate change is narrow:
`NAVIGATOR_REALIGN_KEEP_SCRATCH=1` inverts it for anyone iterating on the pipeline, and a stage-7
failure then costs one hour to retry instead of eleven. Real stage-level resume, keyed on the
source's `source_sig` and the stage inputs already on disk, should be reconsidered before this
ships: the cost of not having it scales with the length of the job, and this job is eleven hours.

> **Built, 2026-08-14** — the second WGS run made the case unarguable; see
> [The second WGS-scale run](#the-second-wgs-scale-run-2026-08-13). `RealignParams::resume` keys on
> the scratch path, which is already derived from the source alignment and the target build, rather
> than on `source_sig`: anything in that directory belongs to the job about to run. The UI opts in.

## The second WGS-scale run (2026-08-13)

The run that followed the BAM switch and the mapper's overlap work was killed at 20:24, five hours
and forty-four minutes in, part-way through the sort's merge. Its cause of death is worth recording
carefully, because the obvious reading of it was wrong in a way that would have sent the next
change to the wrong place.

**It was not a reboot.** `kern.boottime` was unbroken across the whole run — the machine's last
boot predated the run's start by a day. **It was not an out-of-memory kill either**, which was the
working hypothesis: at the moment of death the compressor held 105 pages, no jetsam report was
filed against the process, and anonymous memory sat at 32 GB of 128, most of the rest being
reclaimable file cache. The run log ends after "stage 5/8: Sorting by position" with no panic
appended, which by itself rules out the failure mode of the first run.

What happened is in `WindowServer-2026-08-13-202436.ips`: a `WATCHDOG` termination, "40 seconds
since last successful checkin", against WindowServer's main thread. Killing WindowServer tears down
the graphical login session and every process in it, this job included. To anyone watching the
screen it is indistinguishable from a reboot.

**The pressure that starved it was I/O, not memory.** macOS filed a disk-writes resource notice
against the run for dirtying **549.76 GB of file-backed memory** at 8758 KB/s sustained, against a
limit of 6362 KB/s — the sort writing as fast as it could into the page cache and leaving write-back
to the operating system. That is the right division of labour until the volume gets this far out of
scale, at which point the machine has a debt it cannot settle inside a 40-second watchdog window.

Three things came out of it:

- **`bamio::PacedFile`** flushes on a byte cadence (`NAVIGATOR_IO_SYNC_MB`, 256 MB by default), so
  the write path pays for its own I/O in instalments. It sits under `bamio::create`, the one choke
  point every stage-C write already goes through.
- **`navigator_analysis::resource::ResourceWatch`** samples memory *and* the write rate every 30
  seconds — deliberately inside the 40-second watchdog window it is trying to catch the shadow of.
  It reports and never intervenes. Nothing was recording the number that turned out to matter.
- **Resume**, above. The killed run left 59 GB of complete `mapped.bam` on disk: the revert and the
  mapping, 3 h 58 m, intact and unusable. Resuming from it started the next attempt at the sort.

The sort buffer is worth revisiting separately: at the default 512 MB it spilled **688 runs**, which
the merge then opens at once. That is bounded memory by design and it works, but on a 128 GB machine
it is a lot of fan-in bought for no reason. `NAVIGATOR_SORT_MB` already exists; sizing its default
from installed RAM is not yet done.

## Phase 5 result — WGS229 end to end (2026-08-14)

**Phase 5 passes.** The third WGS run completed in **245.6 min** and every acceptance criterion in
the [Validation plan](#validation-plan) is met. Alignment **#5807**, 43.0 GB,
`chm13v2.0` / minimap2, registered against the same sequence run as its GRCh38 source (#5), which
was read and never written.

The donor (`huF98AFD`) already held native CHM13 alignments, so this is measured against the same
person aligned independently — not against expectation alone.

### Acceptance

| Criterion | Result | |
|---|---|---|
| **Private-Y count** | **11** (0 off-path known, 11 novel, 4 above the publish gate) | dozens, not hundreds — the doc expects ~12 against a 3–39 median |
| Y concordance | R-FGC29071 | ground truth; 1798 markers vs 1811 for native CHM13 (#9) |
| mtDNA concordance | U5a1b1g | ground truth; score and matched count identical to all three native CHM13 alignments |
| Coverage parity | 26.67x mean, 96.1% ≥10x | native CHM13 is 26.95x / 27.64x |
| Ancestry non-regression | European ~100% | matches the validated EUR expectation |
| Y read recovery | chrY breadth **41.28% → 98.14%** | +96k chrY reads, callable 15.28 M → 15.61 M |

That last row is the module's purpose stated as a measurement. On GRCh38 only 41% of chrY receives
coverage at all, most of the remainder being reference gap; the same reads on CHM13 cover 98%,
slightly ahead of both native CHM13 alignments. The gain in *callable* bases is far more modest
(+2.2%) because much of what CHM13 adds is heterochromatin that maps ambiguously — visible as
`poor_mapping_quality` rising from 8.3 M to 45.4 M, which the native alignments show too (44.3 M,
43.9 M) and which is therefore a property of the reference, not of this pipeline.

The de-novo chrY SNP artifact is 655,286 bytes against native CHM13's 655,331 — the two disagree by
45 bytes on a 640 KB call set.

Two caveats worth keeping with the result. Ancestry non-regression rests on the estimate matching
the independently validated EUR expectation, not on a strict before/after diff — `ancestry_result`
was not snapshotted before the run. And the coverage SD is the one place the realignment *beats*
every alternative (22.2 against GRCh38's 190.5), which is collapsed-repeat pileup on GRCh38 rather
than anything the realignment does well.

### What the run cost, and why it is now affordable

| Stage | First run (2026-08-12) | This run | |
|---|---|---|---|
| Revert | 65 min | 64 min | |
| Mapping | 3 h 40 m | 2 h 16 m | mapper overlap work |
| **Sort** | **4 h 44 m** | **29.3 min** | merge readers + `NAVIGATOR_SORT_MB=8192` |
| Mark duplicates | 68.6 min | 11.3 min | not touched — the machine had stopped thrashing |
| Index | crashed | 3.1 min | |
| **Total** | **10 h 41 m** | **4 h 5.6 m** | |

The sort's ten-fold improvement has two causes changed together — the per-run reader fix and a 16x
larger sort budget — and this run cannot separate them. Duplicate marking got six times faster with
no change to its code at all, which is the clearest evidence that the earlier runs were being
starved rather than being slow.

`NAVIGATOR_SORT_MB` was set to 8192 for this run only; the shipped default is still 512 MB, and
sizing it from installed RAM remains open.

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
| **1b. minimap2 via FFI (static)** | Link `libminimap2` in-process via the `minimap2-rs` crate (`minimap2` + `minimap2-sys`), `simde` + `minimap2-sys/static` | **Rejected, including as a bundled cross-check.** Works (S1), but needs a C toolchain, is unproven on Windows, carries an `unsafe` surface, and the obvious feature spelling silently links htslib. Validating against the original does not require shipping it — see below |
| **1c. Pure-Rust aligner** (`minimap2-pure-rs`) | A faithful translation of minimap2 v2.31, pinned to an upstream commit; full preset set, paired-end, CIGAR/SAM/PAF, `.mmi` I/O | **Chosen as default** — 99.74% byte-identical to the C implementation with **zero MAPQ>0 disagreements** and identical accuracy against truth (S4); no `-sys` crates, no `cc`, no build script; indexes CHM13 ~40% faster at the same RAM |

Choosing 1c buys three things at once: **Windows works** with no spike and no graceful
degradation (Decision 3 mostly dissolves), the **single-artifact / no-external-tools posture is
preserved intact** rather than bent, and the htslib trap in S1 never arises. Its cost is the
maturity risk below.

**Put the mapper in a new `navigator-align` crate**, not in `navigator-analysis`: it is a distinct
concern with its own dependency set, and a leaf crate keeps the layering rule intact.

**The cross-check is out-of-band, not a bundled backend.** An earlier revision proposed carrying
the C FFI behind an off-by-default feature so a suspicious realignment could be re-run through it.
That was the wrong shape. Checking a translation against its original is a *development* activity;
building it into the shipped crate would put a C toolchain, an `unsafe` surface, and a
Windows-unproven dependency into an artifact no user ever exercises. Upstream minimap2 is a
perfectly good reference implementation run the ordinary way — over the same FASTQ, diffing the
output — and that is how phase 5 validates. (It was also, in practice, broken: `minimap2-pure-rs`
publishes its lib target as `minimap2`, the same extern name the FFI crate uses, so enabling both
did not resolve. That is a symptom of the shape being wrong, not a reason to work around it.)

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
  that fits the machine, because a coarser split costs MAPQ fidelity (below). The shipped table:
  200 Mbase under 8 GiB, 400 Mbase at 8–15, 1 Gbase at 16–31, unsplit at 32+. Record the chosen
  batch in the cache key and in the alignment's provenance, since it is not a cosmetic knob.

  **Detected, not asked.** The module's audience clicks a button; "bases per index part" is not a
  question they can answer, and a wrong answer is an out-of-memory failure rather than a
  preference. `navigator-align` takes `sysinfo` (`default-features = false`, `system` only) to
  read physical memory — Windows-clean, since it binds through the pure-Rust `windows` crate on
  Windows, `objc2-*` on macOS, and `libc` on Linux, with no C toolchain. Sizing uses **total**
  memory rather than currently-available: the `.mmi` is cached and reused for every later job, so
  a machine that happens to be busy at the first click would otherwise bake a more-split index —
  and its permanent MAPQ cost — into the cache. Available memory is reported separately, for the
  preflight's "can this start now" question. `NAVIGATOR_ALIGN_BATCH_MBASE` overrides both.

- **Implement per-part mapping and merging.** ✅ Built (`navigator-align::map`).
  `navigator-align` owns the orchestration, but not the hard part: `minimap2-pure-rs` exposes
  `index::split::{create_split_tmp, write_split_query_record, read_split_query_record,
  merge_split_query_records}`, so the cross-part MAPQ recomputation — the piece that is easy to
  get subtly wrong — is reused rather than reimplemented. What we own is the loop and the output,
  because the crate's own file-level split entry points (`map_file_pe_sam_split` and friends)
  write to **stdout** and take `parts: &[MmIdx]`, i.e. they hold every part resident and so give
  up the memory bound this whole decision exists to buy.
  `index::reader::IdxReader::read_next` is the part-at-a-time source. Building one whole-genome
  resident index instead is the failure mode that produced this document's original wrong
  estimate; don't reintroduce it.

  **The equivalence is now tested, not assumed.** `a_split_index_places_reads_exactly_where_a
  _whole_index_does` maps the same reads against a 1-part and a multi-part index of the same
  reference and requires byte-identical placements and MAPQ. That test is what makes the memory
  table above safe to rely on; if it ever fails, splitting has become an accuracy decision and
  Decision 4 has to be re-argued.

  Two traps found while building it, both worth knowing before touching this code. Mapping must
  set `MapFlags::CIGAR`: without it `map_query` stops at chaining, so records carry coordinates
  but no CIGAR *and* the block that assigns primary/secondary never runs, leaving every record
  flagged supplementary. And per-part scratch is joined to reads **positionally**, so any
  batching must preserve input order.

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
| **Whole-genome** | Re-map every read | **Built, and the only correct scope.** See below: a read subset cannot be chosen without already knowing the answer. Costs hours per sample. |
| **Y/mt-only** | Re-map only chrY+chrM reads (+ unmapped) | **Withdrawn.** Would return the reads GRCh38 already agreed were Y — the ones needing realignment least — and miss the mismapped ones that are the point. |
| **Targeted (panel ± flanks)** | Extract reads near lifted AIM/IBD sites and realign just those | **Not recommended.** AIM/IBD panels are genome-wide (tens of thousands of sites across all chromosomes), so "targeted" still touches most reads while adding edge-effect risk and missing reads that *moved* between builds — the very signal realignment exists to recover. |

> ### Withdrawn (2026-08-10): there is no useful Y-only scope
>
> A previous revision of this section argued that, since the module is for Y variant discovery,
> Y/mt-only should be the primary scope. That was wrong, and the reason is worth recording so it
> is not proposed again.
>
> **You cannot choose the read subset without already knowing the answer.** The reads that belong
> on the T2T Y are exactly the ones GRCh38 placed wrongly — on chrX, on Y-homologous autosomal
> regions, or nowhere. Identifying them requires letting every read compete against the whole
> reference, which is the computation being asked for. Selecting by source coordinate can only
> return reads GRCh38 already agreed were Y, which is the subset that needed realigning least.
>
> The costs a subset would avoid do not exist either. Revert collates by read name across the
> entire file, so the source BAM is read in full whatever the scope. The index is
> per-`(build, preset)` and cached, so it is paid once for the machine regardless. What remains is
> the mapping pass, and that is the part that cannot be safely narrowed.
>
> **Whole-genome is the only correct scope.** The genome-wide results — coverage, callable, SV on
> hs1 — come along with it rather than being a reason for it.

---

## Data model & provenance

The current `Alignment` has **no parent/derived-from field** — alignments are independent rows
keyed only to a `SequenceRun`. A realigned alignment must record its lineage:

```rust
// navigator-domain::workspace::Alignment — added by migration 0045_alignment_derivation
pub derived_from_alignment_id: Option<i64>,  // source alignment this was realigned from
pub derivation: Option<String>,              // e.g. "realign:minimap2-sr" | "realign:minimap2-map-hifi"
```

Built. `Alignment::is_derived()` answers the UI's question directly, and
`App::{register_realigned_alignment, derived_alignments, derivation_source}` are the write and
read paths. The derivation string is assembled inside `register_realigned_alignment` from a
backend and a preset rather than passed in ready-made, so callers cannot invent a spelling that
later queries fail to recognise.

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
- **Private-Y discovery is the acceptance test.** A GRCh38 vendor sample that cannot be filtered
  by the CHM13 callable mask today should, after realignment, produce a private-Y set of the right
  order of magnitude — **dozens, not hundreds** — once the standard filter stack is applied. The
  private-Y doc measures WGS229 at ~12 privates against a de-novo tree that expects a 3–39 median;
  a realigned sample landing far outside that range means the realignment is manufacturing
  artifacts, not recovering signal. This is the number the module lives or dies by.
- **Ancestry non-regression** (ancestry is not blocked today, and was never the point — see
  [Problem](#problem)): a GRCh38 vendor sample's ancestry estimate via the multi-build panel and
  the estimate from its realigned hs1 alignment must **agree**. Realignment must not move a
  result users have already seen; if it does, that is a finding about one of the two paths.
- **Read-recovery sanity, scoped to the Y:** measure reads that were unmapped or mismapped on
  GRCh38 and now place into chrY — particularly the ampliconic and palindromic regions GRCh38
  collapses. Genome-wide recovery is worth reporting too, but the Y number is the one tied to the
  purpose.
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
2. **`navigator-align` crate** ✅ built: `minimap2-pure-rs` as default backend with the C FFI
   behind an off-by-default feature, preset selection from `test_type`/`platform_name`,
   part-by-part index build and map with the batch size detected from the machine's RAM, the
   `minimap2_index` cache resolved against the shared cache root, single-end and paired-end
   mapping, and SAM/BAM/CRAM output through noodles (BAM by default). Split-vs-whole index
   equivalence is tested for both single and paired reads.
3. **Stage C + Stage D** ✅ built. Stage C is `navigator-analysis::postprocess`: coordinate sort,
   short-read-only duplicate marking on unclipped 5' positions, and a finalise step that moves the
   marked BAM into place and indexes it (`.bai`). Both the sort and the marking verify they are
   lossless. CRAM emit exists and is tested but is not what the job runs — see Stage C above for
   the two `noodles-cram` defects that rule it out at whole-genome scale. Stage D is migration `0045_alignment_derivation` plus
   `navigator-app::realign`: the realigned row is inserted under the source's `sequence_run_id`
   with `derived_from_alignment_id` and `derivation` set, the source untouched, and realigning to
   a build the sample is already on is refused.
4. **App orchestration + UI** ✅ built. `navigator-app::realign_job` drives stages A–D as one
   cancellable job with a disk preflight; the UI is **cards, not a modal** — a job running for
   hours (or days, for a project) must not own the screen. A per-alignment card offers, reports
   progress, and states the outcome, and on a derived row shows its provenance; a project card
   does the same for a sequential batch. Realigned rows are badged in the alignment list, and
   `default_alignment_for_subject` prefers a realignment over its source as a tie-break after
   breadth and depth.
5. **WGS-scale backend parity + MAPQ validation** (replaces the former Windows FFI spike, which
   Decision 1 made unnecessary) — the gate before this is exposed to users.
6. ~~**Y/mt-only mode**~~ — withdrawn; see Decision 5. Every read has to be realigned anyway.

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
  (S4); phase 5 has to say yes at WGS scale, on real reads, against upstream minimap2 run
  normally. The fallback if it doesn't is to shell out to a real minimap2 — which reopens
  Decision 1a's PATH-and-binary-shipping problems, and is the reason phase 5 is the gate.
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
