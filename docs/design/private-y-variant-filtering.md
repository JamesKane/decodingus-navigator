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

| Asset (env override) | Source | Purpose |
|-------|--------|---------|
| `chrY_callable_mask.chm13v2.bed` (14.96 Mbp) — `NAVIGATOR_Y_CALLABLE_MASK` | ytree `results/chrY.callable_mask.chm13v2.bed` | L2 cohort callable mask |
| `chrY_cohort_shared_sites.chm13v2.bed` (~323k positions) — `NAVIGATOR_Y_COHORT_SHARED` | joint-VCF positions with `AC≥2` ∪ ytree `branch_recurrent_exclude` homoplasy hotspots | L3 cohort-shared exclude (the distilled carrier filter) |

Both are loaded best-effort (absent ⇒ that filter is skipped) and gated to **CHM13v2 (hs1)**
alignments — they are in hs1 coordinates, matching the existing `YStructuralRegions` CHM13 gate.
GRCh38/GRCh37 samples fall back to DecodingUs-tree classification + self-callable + region-class only.
The `AC≥2` list is regenerated from the ytree joint VCF with
`bcftools query -f '%CHROM\t%POS\t%INFO/AC\n' chrY.joint.vcf.gz` (keep max-AC ≥ 2). **Follow-up:** ship
via the offline packaging bundle + a sha256 `AssetManifest` entry (same as ancestry/IBD assets).

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

**Still open:** QC banner (warn when novel-in-unique ≫ ~50); AssetManifest entries + packaging-bundle
staging for the two masks; GRCh38 lifted callable mask.

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
- Follow-up: cross-check the 4/2 surviving positions against ytree node `WGS229`'s 12 defining
  variants; confirm a GRCh38 alignment degrades gracefully (DecodingUs GRCh38 coords, no CHM13 masks).

## 6. Limitations

- Single-sample Navigator **cannot** do the cohort carrier filter (F4) or the ≥3-carrier rule (F7) — it
  ships the cohort's precomputed mask/blocklist instead and defers corroboration to AppView.
- Bundled masks are **CHM13v2-only** today; GRCh38/GRCh37 need a lifted callable mask (future).
- The mask reflects the ytree cohort as of its build; refreshing it means re-running stage 3 and
  bumping the asset manifest (same cadence as the tree/ancestry assets).
