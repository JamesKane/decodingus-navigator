# Haploid reassembly caller — design & specification

Status: **design / specification** (prove-out done, no production code). Branch context:
`feat/reassembly-caller` (POC commits `64f8721`, `6e59b12`).
Scope if built: `navigator-analysis` (`caller.rs` de-novo path, new `reassembly.rs`), with the
existing private-Y pipeline (`navigator-app::fastpath`) and downstream filters unchanged.

This is **Option B** from [`private-y-variant-filtering.md`](private-y-variant-filtering.md) §5a/§6 —
the general fix for the single-sample recall gap. Option A (source calls from a GATK gVCF sidecar)
is already shipped, but it only helps samples that arrive with a per-sample chrY gVCF — the
batch-imported ytree minority. **Every ordinary user imports a plain BAM/CRAM with no gVCF, so the
caller itself must resolve these sites.** That is what this module does.

## Problem

Navigator's haploid caller (`caller.rs`) is a **pileup-consensus caller**: at each position it
tallies the passing bases, takes the consensus, and calls a variant when the consensus differs from
the reference and clears an allele-fraction / paralog gate. This is correct and fast for the ~99% of
haploid positions that are cleanly monoallelic, but it has no way to resolve a position whose pileup
is *ambiguous because the reads are misaligned* — precisely the sites GATK's HaplotypeCaller recovers
by local reassembly.

Measured on WGS229 (= `huF98AFD`, the ground-truth donor; ytree truth = 12 private SNPs at
R-FGC29071), §5a of the private-Y doc catalogued exactly where the 12 truth privates are lost:

| ytree private | pileup caller | why |
|---|---|---|
| 3318203 · 4665675 · 7062156 · 16652092 | called | clean af 1.0, depth 11–19 |
| 11008394 · 11913711 | **not called** | q20 depth 3 (< min-depth 4) — genuinely low evidence |
| **4284195 · 11191589 · 20973395 · 21149865** | **not called** | **~50% alt fraction → paralog/allele-balance gate rejects them** |

The last group is the failure this module targets. At these loci, reads from a **paralogous Y region**
(segmental duplication / ampliconic / palindrome) mismap onto the position carrying the *reference*
base, so the true derived allele sits at ~50% and the `is_paralogous` gate (caller.rs) throws the
site out as a suspected bi-allelic artifact. GATK sees the same raw pileup but **reassembles the local
haplotypes and scores each read against them with a base-quality-aware PairHMM**, which down-weights
the low-quality / low-mapping-quality paralog reads and resolves the site.

The root cause is explicit in the code: `realign.rs` does *indel-only* local realignment (homopolymer
smear correction); there is **no SNV haplotype reassembly** (`caller.rs:7` — "no diploid local
reassembly"). So the pileup caller cannot separate a true derived SNV from a paralog artifact when
both present as ~50/50.

**Proven with GATK as oracle** (private-Y §5a): GATK's *single-sample* ploidy-1 gVCF calls **all 12**
truth privates — including the four ~50/50 sites — so the evidence is present in one sample; the gap
is purely the caller. This module closes it.

## Goals

1. Recover the misaligned-ref haploid SNVs the pileup caller drops, on **BAM/CRAM-only** samples
   (no gVCF), so private-Y recall no longer depends on the Option-A sidecar.
2. Match GATK's single-sample ploidy-1 gVCF on the WGS229 truth set as the accuracy bar (recover the
   four ~50/50 sites; the two q20-depth-3 sites are genuinely sub-threshold and out of scope).
3. Emit the **same `VariantCall`** shape the pileup path emits, so every downstream filter (self-mask,
   DecodingUs-tree classification, cohort-shared exclude, region-class, publish gate, QC banner) runs
   unchanged.
4. Add **negligible overhead** to de-novo calling: reassembly runs only over the small set of
   *ambiguous* windows; the cheap pileup path still handles the clean 99%.
5. Stay **pure-Rust and Windows/MSVC-clean** (90% of the user base) — no C-binding bioinformatics
   dependency.
6. Be honest about confidence: carry a per-call quality so borderline recoveries (GATK's own GQ 6–12
   sites) are flagged, not silently published.

## Non-goals

- **Diploid / autosomal calling.** This is a *haploid* reassembler for Y (and later mtDNA). The
  separate `call_denovo_diploid` path is untouched. Genome-wide diploid reassembly is a different,
  much larger problem.
- **Replacing the pileup caller.** Reassembly is a *targeted* pass over active regions; clean
  monoallelic positions stay on the existing fast path.
- **Replacing Option A.** Where a per-sample gVCF exists, the gVCF fast path stays the cheapest route
  (GATK already did the reassembly); this module is the fallback when it doesn't.
- **A general read re-mapper.** Reads are realigned *locally* to candidate haplotypes over a small
  window — not re-mapped genome-wide. Whole-genome re-mapping is the separate
  [`realignment-module.md`](realignment-module.md) track.
- **Structural variants / large indels.** v1 targets SNVs (and short indels as a v2 extension); SV
  detection stays in `sv.rs`.

## Where this sits relative to the existing calling strategies

| Strategy | What it does | Cost | When |
|----------|--------------|------|------|
| **Pileup consensus** (today) | Tally bases, call consensus vs reference | cheap, whole-contig | the clean monoallelic majority of positions |
| **Indel local realign** (today, `realign.rs`) | Re-fit reads over indel windows to un-smear homopolymers | cheap, per-indel-window | ambiguous homopolymer indels |
| **gVCF fast-path** (Option A, shipped) | Read GATK's already-reassembled derived SNVs from a sidecar | ~free | sample arrived with a `*.chry.g.vcf.gz` |
| **Haploid reassembly** (this doc, Option B) | Assemble local haplotypes, score reads with base-quality PairHMM, genotype | modest, per-active-region | ambiguous misaligned-ref pileups the pileup gate rejects, **BAM-only samples** |

The design principle: **escalate only where the cheap path is provably ambiguous.** Active-region
detection (below) is the escalation trigger; everywhere else the pileup consensus already agrees with
GATK, so we don't pay for reassembly.

## Decision 1 — library stack (the constraint that inverts `realignment-module.md`)

Unlike whole-genome re-mapping — where the only production-accuracy option is a ~100k-line C mapper
(minimap2) and we accept an FFI dependency — **local haploid reassembly is small enough to do in pure
Rust to full accuracy.** The whole algorithm is POA assembly + pairwise realignment + a PairHMM over
~200 bp windows, all of which exist as mature, pure-Rust primitives.

**The gate is Windows/MSVC compilability** (per the project's audience: ~90% Windows). C-binding
bioinformatics libraries (abPOA, WFA2, htslib) lean on autotools/Make/POSIX toolchains that do not
build under MSVC. So the decision is the opposite of the realignment module's: **use pure Rust.**

| Option | Provides | Verdict |
|--------|----------|---------|
| **1a. rust-bio (`bio` 4.0.1)** | `bio::alignment::poa` (POA assembly), `bio::alignment::pairwise` (SW/affine realign), `bio::stats::pairhmm` (base-quality PairHMM) | **Chosen** — one MIT dependency, 30 non-optional deps all standard Rust (no `*-sys`/`cc`/`cmake`), compiles clean in-tree, MSVC-clean. Covers the entire algorithm. |
| **1b. abpoa-rs / rust-spoa / libwfa2** | C-binding POA / WFA | Rejected — autotools/Make; won't build under MSVC. |
| **1c. lorikeet-genome** | A GATK-HaplotypeCaller **port** in Rust | **Reference only** — depends on rust-htslib (C) + gkl and unmaintained since 2023. Study its active-region + genotyping logic; don't depend on it. |
| **1d. `rust_wfa`** | Pure-Rust WFA | Fallback only if `bio::alignment::pairwise` semiglobal proves too slow at scale (not observed in the POC). |

CI already builds `x86_64-pc-windows-msvc` on `windows-latest` (`release.yml`) — the real Windows
gate. **Add that job to this branch's PR checks** so the pure-Rust claim is enforced pre-merge, not
just at release.

## Decision 2 — assembly granularity (per-SNV vs POA multi-haplotype)

There are two ways to turn an active region into candidate haplotypes:

- **2a. Per-candidate-SNV** — for each position in the window carrying ≥2 non-reference reads, form
  `H_alt = H_ref` with that single substitution and score `P(read | H_ref)` vs `P(read | H_alt)` per
  read. This is exactly the POC probe, and it already recovers 4/5 of the diagnostic sites. Simple,
  fast, no assembly ambiguity.
- **2b. POA multi-haplotype** — assemble the reads into candidate haplotype(s) with
  `bio::alignment::poa`, align each back to the reference window to enumerate its variants, and score
  every read against every haplotype (a GATK-style read-likelihood matrix). Handles **linked/phased
  variants and short indels** in one pass, which 2a cannot.

**Decision: ship 2a for v1, add 2b as v2.** 2a is proven on real data and covers the private-Y driver
(isolated misaligned-ref SNVs). 2b is the correct generalization for indels and multi-variant windows
but adds assembly-consensus edge cases (POA of ragged low-coverage reads was the POC's first dead end —
see §POC). POA still earns its place in v1 as a **cross-check** (assemble the survivors; confirm the
consensus base at the site), which is cheap and catches per-SNV mistakes.

## Algorithm / pipeline

```
per haploid contig (chrY; chrM in a later phase)
        │
  (A) active-region scan ──►  one pileup pass (we already do this) flags windows where a position is
        │                      ambiguous: minor-allele fraction above the error floor, OR clustered
        │                      soft-clips / indels. Merge+pad adjacent flags into ~200 bp windows.
        │                      Everything NOT flagged → unchanged fast pileup consensus.
        │
  (B) read selection ───────►  pull reads overlapping the window; drop MAPQ < 20 (paralog exclusion);
        │                      collapse overlapping mate pairs to one fragment (dedup); keep per-base
        │                      Phred qualities. This is the layer the POC lacked (the 20973395 miss).
        │
  (C) candidate haplotypes ─►  v1: H_ref + one H_alt per candidate SNV position.
        │                      v2: POA-assembled haplotype(s) (bio::alignment::poa), each aligned back
        │                      to H_ref (bio::alignment::pairwise) to enumerate its variants.
        │
  (D) read-likelihood ──────►  P(read | hap) via base-quality PairHMM (bio::stats::pairhmm) for every
        │                      (read, haplotype): err = 10^(-q/10), Match(1-err)/Mismatch(err/3),
        │                      affine gap-open 1e-4 / extend 0.1, semiglobal. (POC-validated params.)
        │
  (E) genotype ─────────────►  haploid: called haplotype = argmax_h Σ_r ln P(read_r | h). Per site,
        │                      log-odds = Σ ln P(read|alt) − Σ ln P(read|ref); DERIVED if > threshold.
        │                      Phred-scale the log-odds → GQ. Derive depth / alt_depth / AF from the
        │                      per-read best-haplotype assignment.
        │
  (F) emit ─────────────────►  VariantCall { contig, position, ref, alt, depth, alt_depth,
                               allele_fraction, (new) quality } — identical shape to the pileup path,
                               so subtract_known + all app-side filters run unchanged.
```

### Stage A — active-region detection

The escalation trigger, and the key to keeping this cheap. During the existing de-novo tally
(`denovo_chunk` already builds a per-position `[u32; 4]` count array), flag a position as *active*
when either:

- the second-most-common base clears a small read floor **and** its fraction exceeds an error-model
  threshold (i.e. the `is_paralogous` condition — the very sites we currently drop), or
- there is a local cluster of soft-clips or candidate indels (evidence the local alignment is wrong).

Merge flags within a gap tolerance and pad to a fixed window (~150–250 bp, sized so short reads span
it — the POC used ±40 bp around a single site; a multi-SNV window wants the wider span). Positions
that never flag stay on the current pileup consensus and cost nothing extra. On chrY this makes
reassembly run over a **small fraction** of the contig — the ambiguous loci only.

### Stage B — read selection (the ingredient the POC still lacked)

For each active window, query the indexed BAM/CRAM (`reader::open_indexed`, decode-safe
`reader::decode_pool` for CRAM 3.1) and keep reads that:

- are primary, non-dup, non-QC-fail (existing `passes` filter), **and MAPQ ≥ 20** (GATK default) —
  this excludes ambiguously-placed paralog reads;
- **collapse overlapping mate pairs to a single fragment** — where a read and its mate both cover the
  site, count the fragment once (take the higher-quality base on disagreement). This is what reduces
  GATK's own DP at 20973395 from 7→5, and it is exactly what the POC omitted (its extra reference
  reads outvoted the truth). Fragment dedup + MAPQ gate together are expected to close the 5th site.

Carry each read's window bases **and per-base Phred qualities** (the PairHMM's whole point).

### Stage C — candidate haplotypes

v1: `H_ref` = reference window (uppercased — the FASTA is soft-masked lowercase, a POC gotcha), plus
one `H_alt` per position with ≥2 non-reference fragments (the alt = the majority non-ref base there).

v2: POA-assemble the selected reads (`bio::alignment::poa`, heaviest-bundle consensus), then
semiglobally align each assembled haplotype back to `H_ref` (`bio::alignment::pairwise`) to read off
its SNVs/indels. Bound the haplotype count (haploid ⇒ expect 1–2 real haplotypes; extra POA bundles
are noise). POA also serves as the v1 cross-check.

### Stage D — read↔haplotype likelihood (PairHMM)

For each (read, haplotype) compute `P(read | haplotype)` with `bio::stats::pairhmm::PairHMM`, using
the POC-validated emission model:

- emission: `err = 10^(−q/10)` (Phred clamped to Q2–Q60); `read[i] == hap[j]` → `Match(1−err)`, else
  `Mismatch(err/3)`;
- gap: affine, open `1e-4`, extend `0.1`; semiglobal in the read (free window-edge offset).

This is the base-quality-aware tie-breaker the crude match/mismatch score lacks.

### Stage E — haploid genotyping

Haploid ⇒ one true haplotype. Called haplotype = the one maximizing `Σ_r ln P(read_r | h)`. Per
candidate site, `log-odds = Σ ln P(read|H_alt) − Σ ln P(read|H_ref)`; call **DERIVED** when
`log-odds > τ` (POC τ ≈ 2 nats ≈ Phred ~9; calibrate on the truth set). Assign each read to its
max-likelihood haplotype to derive `depth` / `alt_depth` / `allele_fraction`, and **Phred-scale the
log-odds into a `quality`/GQ** so the publish gate and QC banner can treat low-confidence recoveries
(GATK's GQ 6–12 sites) honestly.

### Stage F — emit & integrate

Emit `VariantCall`s at the called-haplotype's differing positions. `subtract_known` (tree-position
removal) and the entire app-side private-Y filter stack (`private_y_core`: self-mask, DecodingUs-tree
classification, cohort-shared exclude, region-class, publish gate, QC banner) consume these
unchanged. **Bump `DENOVO_VERSION`** to e.g. `haploid-denovo-3` so cached de-novo artifacts
invalidate (the source-signature cache keys on it).

## Integration into `caller.rs`

The reassembly pass slots into `denovo_chunk`, after the tally, as a replacement for *dropping*
active positions:

- Today: `tally_region` → (optional indel `realign_region`) → per-position consensus with
  `is_paralogous` **rejecting** the ambiguous sites.
- Proposed: `tally_region` → active-region detection → for each active window, run the reassembly
  resolver (new `reassembly.rs`); for non-active positions, the current consensus path is unchanged.

New module `crates/navigator-analysis/src/reassembly.rs` (pure, unit-tested, mirroring how
`realign.rs` is a pure module `caller.rs` drives). It owns Stages B–E over a single window; `caller.rs`
owns Stage A (it already has the counts) and Stage F (it already builds `VariantCall`s). A
`HaploidCallerParams` field (e.g. `reassembly: bool`, defaulting on for chrY) gates the pass so it can
be toggled for parity testing against the pileup-only baseline.

`VariantCall` gains an optional `quality: Option<f64>` (Phred GQ) — additive, `#[serde(default)]`, so
existing cached artifacts and the pileup path (which can leave it `None`) stay compatible.

## POC evidence (branch `feat/reassembly-caller`)

`crates/navigator-analysis/examples/reassembly_probe.rs`, run on the real WGS229 CHM13 CRAM, scored
against GATK's own chrY gVCF (which calls all five sites `GT=1`/DERIVED):

| site | pileup | GATK gVCF | crude realign | **base-qual PairHMM** | result |
|---|---|---|---|---|---|
| 3318203 (control) | C@1.00 | AD 0,17 GQ99 | 0/17 → alt | +142.3 | ✅ DERIVED |
| 16652092 (control) | T@1.00 | AD 0,10 GQ99 | 0/11 → alt | +103.5 | ✅ DERIVED |
| **4284195** | T@0.50 | AD 9,10 GQ44 MQRankSum −3.55 | **10/10 tie** | **+3.8** | ✅ **DERIVED** |
| **11191589** | T@0.43 | AD 2,3 GQ25 | 3/3 tie | +9.7 | ✅ **DERIVED** |
| 20973395 | A@0.43 | AD 2,3 GQ41 | 4/3 | −4.4 | ❌ ancestral |

Findings that shaped this design:

- **Base-quality PairHMM is the load-bearing tie-breaker.** Crude match/mismatch scoring *ties* the
  misaligned-ref sites (4284195 at 10/10); the PairHMM breaks them exactly as GATK does, recovering
  4/5.
- **The 5th site (20973395) is a read-selection problem, not a scoring one.** Its paralog reference
  reads clear the MQ≥20 gate; GATK excludes them via active-region **fragment/read selection** (its DP
  falls 7→5). Hence Stage B's fragment dedup — the POC omitted it and was outvoted.
- **POA of ragged low-coverage reads is fragile** (the POC's first dead end: full-span POA of partial
  reads produced garbage consensus). Hence v1 uses per-SNV PairHMM with POA as a *cross-check*, and
  defers full POA multi-haplotype assembly to v2.
- **FASTA soft-masking gotcha:** `reader::read_contig_sequence` returns lowercase soft-masked bases;
  uppercase the reference window before scoring.

## Performance & resource profile

- Reassembly runs only over **active windows** — a small fraction of chrY — so the added cost over the
  current de-novo pass is dominated by the active-region scan, which reuses the tally we already
  compute. Net overhead: modest.
- PairHMM is `O(read_len × hap_len)` per (read, haplotype); over ~200 bp windows with a handful of
  haplotypes and tens of fragments, each window is sub-millisecond. The POC scored five windows on a
  full-CRAM run in ~50 s wall — and essentially all of that was the CRAM open/seek per site, not the
  HMM; the production path amortizes the open across the whole contig.
- Reuses the existing per-chunk rayon parallelism (`decode_pool`) and CRAM-decode stack safety
  (`NAVIGATOR_DECODE_STACK_MB`).

## Validation plan

- **Truth set:** WGS229's 12 ytree privates at R-FGC29071. Target = recover the 4 misaligned-ref
  sites (4284195, 11191589, 20973395, 21149865) that the pileup gate drops, matching GATK's
  single-sample gVCF; the 2 q20-depth-3 sites (11008394, 11913711) stay out of scope (genuinely
  sub-threshold). Success = **10/12** from a BAM-only run, parity with the Option-A gVCF path.
- **No-regression:** the clean sites (3318203 · 4665675 · 7062156 · 16652092) and every currently
  correct call must be unchanged; run the pileup-only baseline vs reassembly and diff.
- **Specificity:** reassembly must not *add* false privates — the whole point is to reject paralog
  artifacts while keeping true derived SNVs. Confirm the post-reassembly DISPLAY/PUBLISH counts on
  WGS229 stay at the validated single-digit level (private-Y §5a: DISPLAY 15, PUBLISH 4/4 truth via
  gVCF; BAM reassembly should land in the same neighborhood).
- **HiFi cross-check:** the GFX HiFi donor (private-Y: PUBLISH 2 @ R-FGC29071) must not regress under
  the longer-read regime (the publish gate already special-cases HiFi alt-depth).
- **Determinism:** fixed thread count + fixed params → reproducible calls (same discipline as the rest
  of the analysis engine).
- **Unit tests** in `reassembly.rs` on small synthetic fixtures: a clean derived site, a balanced
  paralog artifact (must stay rejected), an overlapping-mate double-count (dedup must collapse it).

## Phasing / rollout

1. **`reassembly.rs` core (Stages B–E, per-SNV / v1)** + unit tests on synthetic windows — pure,
   backend-agnostic, no `caller.rs` wiring yet. This is the algorithm, proven in isolation.
   **DONE** (commit `eb4270c`).
2. **Active-region detection + `caller.rs` wiring (Stage A/F)** behind a `HaploidCallerParams`
   toggle; `DENOVO_VERSION` bump; end-to-end on the WGS229 CRAM. **DONE** — see the measured result
   below.
3. **Windows CI on the branch** — confirm the pure-Rust claim under MSVC before merge.
4. **POA multi-haplotype (v2)** — linked variants + short indels via `bio::alignment::poa`; extends
   coverage beyond isolated SNVs. **This is also what closes `20973395`** (active-region read
   selection beyond the MQ gate — see the measured result).
5. **mtDNA (v3)** — apply to chrM, but heteroplasmy needs *fractional* genotyping (not the haploid
   argmax), so it's a genuine extension, not a config flip.

Phases 1–3 deliver the private-Y driver (misaligned-ref SNV recovery for BAM-only samples) and are
the merge target; 4–5 are follow-ons.

### Measured phase-2 result (WGS229 CHM13 CRAM, BAM-only, reassembly OFF vs ON)

`examples/reassembly_validate.rs` runs the bounded de-novo caller with reassembly off (pileup) vs on:

| truth private | pileup | reassembly | GATK gVCF | outcome |
|---|---|---|---|---|
| 3318203 · 16652092 (+ 4665675 · 7062156) | called | called | GQ99 | no regression |
| **4284195** | dropped | **C>T 10/20 GQ16** | GQ44 | ✅ recovered |
| **11191589** | dropped | **C>T 2/4 GQ40** | GQ25 | ✅ recovered |
| 20973395 | dropped | dropped | GQ41 | deferred to v2 (paralog reads clear MQ≥20; needs active-region read selection — logodds −4.4) |
| 21149865 | dropped | dropped | **GQ8 / QUAL 0.18** | marginal even for GATK (3-allele 4T/4A/1G) |
| 11008394 · 11913711 | dropped | dropped | (low) | out of scope (q20 depth 3 < min-depth 4) |

v1 recovers the two misaligned-ref sites resolvable by the MQ gate + base-quality PairHMM (matching
the POC exactly), with **no regression** on the clean calls — which is structural: reassembly only
*escalates positions the paralog gate already dropped* and *appends* recoveries, so it can never
remove an existing pileup call. `20973395` is the concrete v2 driver; `21149865` is GATK-marginal.

## Open questions

- **Active-region window size & merge tolerance** — POC used ±40 bp around one site; multi-SNV windows
  want ~150–250 bp. Calibrate against read length and the density of real linked Y variants.
- **Genotype threshold τ & GQ calibration** — the POC's `log-odds > 2` recovered the clean cases;
  where exactly to set τ (and how to Phred-scale it into a GQ the publish gate trusts) wants tuning on
  the full truth set, ideally cross-checked against GATK's GQ at the same sites.
- **Fragment dedup on disagreement** — take the higher-quality base, or drop the fragment as
  discordant? GATK collapses to a consensus; confirm which rule flips 20973395 without new FPs.
- **Interaction with the existing indel `realign.rs`** — should active-region reassembly *subsume* the
  homopolymer indel realignment (one unified local reassembly), or stay a separate SNV-only pass over
  `realign.rs`'s output? v1 keeps them separate; v2/POA is the natural point to unify.
- **mtDNA heteroplasmy** — fractional allele genotyping vs the haploid argmax; defer to phase 5 but
  keep the `reassembly.rs` API from hard-coding "one haplotype wins."

## Sources

- rust-bio (`bio` 4.0.1; POA + pairwise + `stats::pairhmm`; MIT, pure Rust):
  <https://docs.rs/bio/latest/bio/> · <https://crates.io/crates/bio>
- `bio::stats::pairhmm` (the base-quality PairHMM used for read↔haplotype likelihood):
  <https://docs.rs/bio/latest/bio/stats/pairhmm/index.html>
- GATK HaplotypeCaller (local reassembly + PairHMM, the algorithm this ports to the haploid case):
  <https://gatk.broadinstitute.org/hc/en-us/articles/360035531412>
- lorikeet-genome (Rust GATK-HC port; reference only — rust-htslib/C, unmaintained):
  <https://github.com/rhysnewell/Lorikeet>
- POC + measured evidence: this repo, branch `feat/reassembly-caller`
  (`crates/navigator-analysis/examples/reassembly_probe.rs`, commits `64f8721`, `6e59b12`) and
  [`private-y-variant-filtering.md`](private-y-variant-filtering.md) §5a.
