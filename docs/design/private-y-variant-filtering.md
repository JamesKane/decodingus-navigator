# Private-Y Variant Filtering Plan

**Status:** design / not yet implemented
**Trigger:** `huF98AFD` (= `WGS229`, my own genome) reports far too many "Private Y" variants. We are
surfacing SNPs that should be pre-filtered — sending users on a wild goose chase and threatening to
flood AppView curators with non-viable novel-SNP claims.

## 1. The evidence: what "sensible" looks like

The `/Users/jkane/Genomics/ytree` de-novo tree pipeline joint-genotyped **~3,352 CHM13-aligned males**
(1000G + HGDP + PGP + ancient + `mine/WGS229`) and, after its full filter stack, places **WGS229 at
`R-FGC29071` with exactly 12 private/novel defining SNPs** (`chrY.ftdna.refined.ingest.json`, node
`WGS229`). Across all WGS samples the pipeline's per-tip novel-SNP median is **3–39** (BigY, which
*discovers* novels, is higher at ~85).

So the correct order of magnitude for a modern short-read WGS sample's private Y set is **dozens, not
hundreds or thousands.** If Navigator shows more than that, it is emitting artifacts.

Why the flood happens: CHM13v2's Y is HG002's Y — **haplogroup J**. Our sample is **R**. A single-sample
caller therefore sees the *entire J-vs-R divergence* plus every mismapped ampliconic/palindromic/
heterochromatic read as a "derived, not-in-tree" call. The ytree pipeline learned this the hard way:
without its carrier filter, WGS samples leaked **641 fake privates per tip** (4,969 spurious sub-clades);
with it, **53/tip** — and after the callable mask + recurrent exclusion, ~12 for WGS229.

## 2. What the ytree pipeline does to get viable SNPs

The pipeline's filter stack, in order of impact:

| # | Filter | Artifact / parameter | Effect |
|---|--------|----------------------|--------|
| F1 | **Non-PAR restriction** | `chrY:2,458,320–62,122,809` (`chrY_nonPAR.chm13v2.bed`) | PAR recombines with X → never a Y marker |
| F2 | **Cohort callable mask (Poznik ≥90%)** | `chrY.callable_mask.chm13v2.bed` — CALLABLE (depth ≥4, MQ ≥20) in **≥90% of ~3,000 males** | **14.96 Mbp = only 25.1% of non-PAR.** 75% of the Y is *structurally* unreliable. **The single biggest lever.** |
| F3 | **Per-call gates** | `min-depth 4, MQ 20, min-base-qual`, GQ ≥50 for `jointConfirmed` | drops low-evidence calls |
| F4 | **Cohort carrier filter** *(the key fix)* | keep an off-tree SNP only if **all** cohort carriers sit at one terminal | kills backbone/divergence leakage (641→53/tip). **Inherently cohort-wide.** |
| F5 | **Recurrent / homoplasy exclusion** | `branch_recurrent_exclude.chm13v2.bed` (56 positions on ≥3 branches) + `low_confidence_defining.chm13v2.tsv` (**80,720 positive-only + 12,295 homoplasic** positions) | drops recurrent-mutation hotspots and never-joint-confirmed positions |
| F6 | **Indels/spanning-deletions → N**, multiallelic decoded per-sample | build_phylip | avoids STR/homopolymer indel noise |
| F7 | **≥3-carrier rule for a real marker** | a per-sample singleton "can't be a ≥3-carrier marker" | a lone private is a *candidate*, not a branch |
| F8 | **aDNA damage trim** (ancient only) | 7 bp terminal qual→0 + transversion-restriction | not relevant to modern WGS |

## 3. Gap analysis — Navigator today

Traced through `caller.rs`, `fastpath.rs`, `queries.rs`, `publish.rs`, `mask.rs`:

**What we already have (good):**
- Caller intrinsic gates (`caller.rs` `HaploidCallerParams`): MAPQ ≥20, base-qual ≥20, depth ≥4
  (≥1 HiFi), allele-fraction ≥0.5, and a paralog/allele-balance filter (`is_paralogous`). ≈ F3.
- In-app `private_y_core` (`fastpath.rs`): backbone subtraction + **optional** callable-mask +
  region **annotation**. Splits `OffPathKnown` vs `Novel`.
- `mask.rs` `YRegionClass` (PAR/palindrome/amplicon/XTR/STR/centromere/heterochromatin + modifiers)
  and `YStructuralRegions::classify`, gated to CHM13v2.

**The gaps that cause the flood:**
1. **The published record bypasses every private-path filter.** `publish.rs::variants_record` builds
   `PrivateVariantsRecord` **directly from `cached_denovo`** (raw whole-chrY de-novo set) — no backbone
   subtraction, no mask, no region logic. `ibd_exchange.rs::publish_variants` does the same. **This is
   what floods curators.** (Highest-priority fix.)
2. **No cohort callable mask by default.** `private_y_variants(None)` applies no mask; the
   `_self_masked` path uses the *sample's own* callable intervals — which still include paralog-prone Y
   that maps in one sample. We never apply the cohort **≥90%** mask, so ~75% of ineligible territory
   stays in scope.
3. **`YRegionClass` is annotation-only — it never drops a variant.** A novel in a palindrome/amplicon/
   heterochromatin block is reported as if it were unique sequence. (Reporting *counts*
   `novel_in_unique_sequence()` but nothing gates on it.)
4. **No recurrent-position / low-confidence blocklist** exists in Navigator.
5. **No "publishable" tier.** Placement gates (depth ≥4, AF ≥0.5) are correct for *confirming a
   backbone SNP* but far too loose for *proposing a novel branch marker to a curator*, where a haploid
   call should be near-homozygous (AF ≈1.0) at solid depth.
6. **Single-sample cannot reproduce F4/F7** (cohort carrier / ≥3-carrier). Those are AppView-side by
   nature — Navigator must publish novels as *candidates*, not as proposed branches.

## 4. The plan

Two principles:
- **Ship the cohort's verdict as bundled assets.** Navigator is single-sample and can't recompute F2/F4/F5
  on the fly — but the ytree pipeline *already produced* those cohort-derived artifacts. Bundle them.
- **Two tiers.** *Display* everything (with honest flags) so a user can explore; *publish/escalate* only
  the high-confidence, unique-sequence, in-mask subset so curators aren't flooded.

### 4.1 Layered filter (applied on BOTH the in-app and the publish path)

```
raw de-novo chrY calls
  └─ L0  caller intrinsic gates            (exists: MAPQ≥20, BQ≥20, depth, AF≥0.5, paralog)
  └─ L1  backbone subtraction vs the       ★ classify against the DecodingUs tree, not FTDNA's
         DecodingUs tree (!path & known)      → shared lineage variants read as off-path-known, not novel
  └─ L2  cohort callable mask (≥90%)        (F2 — BUNDLE; CHM13 only)
  └─ L3  cohort-shared-sites exclude        ★ F4 distilled: drop positions with ≥2 cohort carriers
         (AC≥2 ∪ homoplasy)                    (a real shared variant belongs in the tree; shared-but-
                                                untree'd = suspect artifact; a true private is AC=1)
  └─ L4  structural-region class            (YRegionClass — annotated in DISPLAY, excluded from PUBLISH)
  └─ L5  publishable-tier gate              (novel-marker gates: AF≥0.9, alt-depth≥10 SR / ≥3 HiFi)
        ├─ DISPLAY  = L0–L4 (region kept as a labeled bucket, not published)
        └─ PUBLISH  = L0–L5 AND class==Novel AND region.is_none()  → "singleton candidate"
```

- **★ L1 (classify against the DecodingUs tree) and L3 (cohort-shared exclude) are the real levers —
  not the callable mask.** Empirically (WGS229), the callable mask alone removed only ~11% (764→701
  publishable) because the excess "novels" are *not* in non-callable regions — they are shared-R
  lineage variants. The two changes that matter:
  1. **L1** — `private_y_core` was classifying against the *FTDNA* tree while the rest of the app
     places against the *DecodingUs* tree (which folds in the de-novo pipeline's cohort branches).
     Switching it made shared-R variants resolve as `OffPathKnown` (off-path-known jumped 49→553),
     collapsing publishable 701→151.
  2. **L3** — the residual 151 are "shared across the cohort yet absent from the tree" = suspect. A
     bundled exclusion of every ≥2-carrier cohort position (from the joint-VCF `AC`) drops them,
     because a *true* private has exactly one cohort carrier (`AC=1`) and survives. This is the
     single-sample stand-in for the de-novo pipeline's cohort **carrier filter** (F4), which no
     single-sample analysis can compute live. Publishable 151→**4**; DISPLAY total 1481→**9**.
- **L4:** a `Novel` call inside PAR/palindrome/amplicon/XTR/STR/centromere/heterochromatin is not
  viable. Kept in the DISPLAY view under a `structural` bucket (`in_structural_region()`) but never
  published (the L5 gate requires `region.is_none()`).
- **L5 publishable gate:** haploid chrY should be effectively homozygous. Require **AF ≥ 0.9** (0.5–0.9
  is a mixed/paralog/contamination signal) and **alt-depth ≥ 10** for short-read (≥3 for HiFi). This is
  separate from the placement params — placement stays permissive; *publishing a claim* is strict.

### 4.2 Bundled assets (env-overridable; live under `<cache base>/masks/`, CHM13-only)

Committed (gzipped) under repo `assets/masks/`, seeded to `<cache base>/masks/` on first run:

| Asset (env override) | Build | Source | Purpose |
|-------|-------|--------|---------|
| `chrY_callable_mask.{chm13v2,grch38}.bed.gz` — `NAVIGATOR_Y_CALLABLE_MASK` | CHM13 native / GRCh38 lifted | ytree `results/chrY.callable_mask.chm13v2.bed` | L2 cohort callable mask |
| `chrY_cohort_shared_sites.{chm13v2,grch38}.bed.gz` — `NAVIGATOR_Y_COHORT_SHARED` | CHM13 native / GRCh38 lifted | joint-VCF `AC≥2` ∪ ytree `branch_recurrent_exclude` homoplasy | L3 cohort-shared exclude (the distilled carrier filter) |

Loaded best-effort per build via `y_mask_build_token` (CHM13 → `chm13v2`, GRCh38 → `grch38`, else none);
absent ⇒ that filter is skipped. **Regeneration:** the `AC≥2` list comes from the ytree joint VCF
(`bcftools query -f '%CHROM\t%POS\t%INFO/AC\n' chrY.joint.vcf.gz`, keep max-AC ≥ 2 ∪ homoplasy). The
**GRCh38** files are the two CHM13 BEDs lifted with CrossMap and the UCSC `hs1ToHg38.over.chain`
(`CrossMap bed <chain> in.bed out.bed`, then chrY-only + sort/merge). See `assets/masks/README.md`.
**Follow-up:** sha256 `AssetManifest` entries for the masks (they currently seed unverified).

### 4.3 Code changes — IMPLEMENTED

1. **Publish path fixed.** `publish.rs::variants_record` now routes **chrY** through the filtered,
   publish-gated bucket (`private_y_variants_self_masked` → `publishable(gate)`), not raw
   `cached_denovo`. chrM keeps the raw path. `ibd_exchange.rs::publish_variants` inherits it. This
   stops the curator flood on its own.
2. **Classify against the DecodingUs tree.** `private_y_core` switched from `fetch_ftdna_y_tree` to
   `y_decodingus_tree_calls` (FTDNA fallback), so shared-lineage variants resolve as `OffPathKnown`.
3. **Cohort assets bundled + applied** in `private_y_core` (`load_y_position_bed`): L2 callable mask +
   L3 cohort-shared exclude, CHM13-gated. Private-Y cache version bumped `1→2`.
4. **Publishable tier** = `PublishGate { min_allele_fraction: 0.9, min_alt_depth: 10 (3 HiFi) }` +
   `PrivateBucket::publishable()/publishable_count()`; `PrivateVariant` gained `alt_depth`. Structural
   region stays annotation in DISPLAY, excluded from PUBLISH (`region.is_none()`).
5. **`callable_chr_intervals` hardened** to resolve the reference via the gateway (CRAMs with no stored
   `reference_path` previously errored on the self-masked path).
6. **`private-y` CLI** diagnostic (`navigator private-y --alignment N`) reports DISPLAY vs PUBLISH counts.
7. **QC banner** (`private_y_qc_banner`, threshold `PRIVATE_Y_QC_WARN = 50`) in the report section +
   `PrivateBucket::qc_banner()` + CLI. **Publishing is skipped entirely when it fires** (an
   artifact-laden sample federates zero candidates; DISPLAY still shows them under the banner).
8. **Masks bundled as committed repo assets** (`assets/masks/`, gzipped; `RegionMask::from_bed` gunzips
   transparently), seeded on first run (`seed_bundled_masks`) + staged for packaging.
9. **GRCh38 lifted masks** (CrossMap `hs1ToHg38`), selected per build by `y_mask_build_token`.
10. **GVCF-sourced calls** (§5a): `private_y_core` prefers a per-sample chrY gVCF sidecar
    (`gvcf::read_derived_snvs`) over the pileup caller — GATK's reassembly recovers SNVs the pileup
    caller misses. Self-mask skipped on this path; cache v2→v3.

**Still open:** Option B — SNV local reassembly in Navigator's pileup caller (for samples with only a
BAM, no gVCF); sha256 `AssetManifest` entries for the masks; GRCh37 masks.

### 4.4 AppView-side (F4/F7 — out of Navigator's scope, stated for completeness)

Navigator publishes filtered **candidates**. The cohort verdict — a novel is a real marker only if it
recurs (**≥3 carriers**) and its carriers are monophyletic — is AppView-side, exactly like the
cross-analyst dedup already agreed for external identifiers. Curators should see candidates promoted to
markers only after AppView corroboration, never a raw per-sample singleton dump.

## 5. Validation (measured 2026-07-07)

`huF98AFD` = `WGS229`, both alignments place at **R-FGC29071** (ytree ground truth: 12 private SNPs).

| Stage | WGS229 Illumina/CHM13 (aln 9) | GFX0457637 HiFi/CHM13 (aln 6) |
|-------|------------------------------|------------------------------|
| DISPLAY total (was raw ~thousands) | **9** | 144 (124 in structural zones) |
| novel in unique sequence | 7 | 19 |
| **PUBLISH** (old path published raw denovo) | **4** | **2** |

Progression on WGS229 as each lever was added: self-mask only → PUBLISH **764**; + callable mask →
**701** (mask alone is *not* the lever); + DecodingUs-tree classification → **151**; + cohort-shared
exclude → **4**. This tracks the ytree "dozens, not thousands" target and single-digit publishable.

- Both same-person alignments stay at R-FGC29071; the HiFi DISPLAY is higher only because HiFi reads
  reach further into ampliconic/structural zones — those are flagged `structural` and excluded from
  PUBLISH, so publishable stays at 2.

**GRCh38 (aln 5, same person, WGS229.b38.bam) — the lifted masks help but GRCh38 stays noisier:**
without masks PUBLISH **190** / DISPLAY 441; with the lifted GRCh38 masks PUBLISH **148** / DISPLAY 337.
It does *not* reach CHM13's cleanliness — GRCh38's chrY reference is intrinsically noisier, and the
DecodingUs tree's de-novo cohort nodes carry hs1 coords, so shared-R leaks as novel (off-path-known is
only ~25 vs CHM13's 553). Two things make this safe rather than a regression: (1) the **QC banner fires**
(318 novel ≥ 50), and (2) **publishing is skipped entirely when the banner fires** — a GRCh38 sample
flagged as artifact-laden federates *zero* candidates, though the calls still appear in the in-app
DISPLAY under the banner. Practical guidance: **use a CHM13 alignment for private-Y**; GRCh38 is
filtered best-effort and honestly flagged.

## 5a. Recall ceiling — why WGS229 shows fewer than its 12 ytree privates (measured 2026-07-07)

Checked all 12 of WGS229's ytree-truth privates against every stage (`navigator private-y` +
`NAVIGATOR_DUMP_DENOVO`). Result: the mask/tree/cohort **post-filters keep 100 % of what the caller
emits** — the shortfall is entirely upstream, in single-sample calling:

| ytree private | raw de-novo? | reason |
|---|---|---|
| 3318203 · 4665675 · 7062156 · 16652092 | called | clean af 1.0, depth 11–19 → **published** |
| 14650090 | called | af 1.0, alt-depth 8 → dropped by the publish gate (≥10) |
| 21626529 | called | af 1.0, but AZF/amplicon region → region-excluded |
| 11008394 · 11913711 | **not called** | q20 depth **3** (< caller min-depth 4) |
| 4284195 · 11191589 · 20973395 · 21149865 | **not called** | **~50 % alt fraction** → paralog/allele-balance gate |

**This was first (wrongly) attributed to a single-sample-vs-joint ceiling. The test below refuted that.**
Tested against GATK's *per-sample* `WGS229.chm13.chrY.g.vcf.gz` (HaplotypeCaller ploidy-1, **pre-joint**):
GATK single-sample calls **all 12 derived** — including the 6 Navigator missed (`4284195` DP19/GQ44; the
low-depth DP4–5 ones). So the 6 are single-sample-callable; the gap is **Navigator's pileup caller**, which
has no SNV haplotype reassembly (`local_realign` is indel-only, caller.rs:7) and so rejects the misaligned-ref
~50/50 pileup that GATK's reassembly resolves. (The served tree is *not* the issue — it is complete; per-sample
private leaves are correctly not served, so off-path-known = 0 for a sample's own privates is by design.)

**Fix — Option A (implemented): source chrY calls from a per-sample GVCF when present** (§4.3 #10). GATK's
gVCF already did the reassembly, so `private_y_core` reads its derived SNVs instead of the pileup caller;
downstream filters are identical. `gvcf::read_derived_snvs` (min_dp 4, min_gq 20) + `chr_y_gvcf_for_alignment`
(`NAVIGATOR_Y_GVCF` override or a `*.chry.g.vcf.gz` sibling). The self-mask is skipped on this path (the gVCF
is the callable evidence; also avoids a CRAM walk). **Result: 9/12 recovered** (was 6/12), DISPLAY 15, PUBLISH
4/4 truth; the min_dp-4 floor drops the misaligned DP2–3 artifact clusters without touching the DP≥4 truth.
The 3 still missing (GQ 6/8/12) are the ones GATK's own single-sample GQ flags as needing cohort confirmation.
Option B (SNV local reassembly in Navigator's caller, for samples with only a BAM) remains the general fix.

## 6. Limitations

- Single-sample Navigator **cannot** do the cohort carrier filter (F4) or the ≥3-carrier rule (F7) — it
  ships the cohort's precomputed mask/blocklist instead and defers corroboration to AppView.
- **Recall is bounded by the caller, not by single-sample data** (§5a): GATK's single-sample gVCF calls all
  12; Navigator's pileup caller (no SNV reassembly) gets 6. A gVCF sidecar closes most of the gap (9/12);
  the general fix for BAM-only samples is SNV local reassembly (Option B). Post-filters are lossless over
  whatever the source emits.
- **GRCh38 private-Y is intrinsically noisier than CHM13** and the lifted masks only partly close the gap
  (its chrY reference is worse, and the tree's de-novo cohort nodes are hs1-native so shared-R leaks as
  novel). The QC banner + publish-skip make this safe, but CHM13 is the recommended reference.
- Masks ship for **CHM13v2 + GRCh38**; GRCh37 has none (bare-`Y` naming and the de-novo path is chrY-only).
- The mask reflects the ytree cohort as of its build; refreshing it means re-running stage 3 and
  bumping the asset manifest (same cadence as the tree/ancestry assets).
