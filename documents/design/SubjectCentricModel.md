# Subject-centric analysis model (design)

Status: **P1–P3 implemented** (commits ba6ffc2 / 9259f7a / d928af8). The remaining open item is
true genotype-level pooling for ancestry (see Phasing/Open decisions). Captures the shift from
run-centric tabs to donor-level aggregate reports.

## Problem

The detail tabs (Y-DNA / mtDNA / Ancestry / IBD) are driven by a **selected alignment** — one
sequencing run. But a subject (donor) commonly has **many** sources: a WGS BAM, a Big-Y, an
mtDNA full sequence, several STR panels, a chip array. The genealogically-meaningful answer for a
person is the **consensus across all of their data**, not "the result for this one run."

Symptoms today:
- You must drill into Data Sources and pick an alignment before any analysis tab works.
- Per-run results are shown as if they were the subject's answer.
- Re-running analysis re-does work per alignment rather than reusing a donor-level result.

## What's already donor-level vs run-level

| Data | Today | Donor-level machinery exists? |
|------|-------|------------------------------|
| Y haplogroup | per-alignment assign; **consensus_y** reconciles all recorded Y calls | **Yes** — `reconciliation::reconcile` over `haplogroup_call` rows (any source) |
| mtDNA haplogroup | per-alignment assign; **consensus_mt** likewise | **Yes** — same path |
| Ancestry | per-alignment estimate (`ancestry_result` keyed by alignment + method) | No — no donor rollup |
| STR profile | per-panel rows (`str_profile`) | Partial — multiple profiles stored, no consensus marker set |
| Coverage / sequencing QC | per-alignment (correct — it *is* a per-run property) | n/a (stays per-run) |
| Private-Y | per-alignment bucket | No — no cross-run union |

So Y/mt are **already** reconciled donor-level; the UI just doesn't lean on it. Ancestry, STR,
and private-Y need a donor rollup.

## Target model

Each tab presents the **donor aggregate**, with per-source contributions as drill-down. Analysis
*actions* (assign/estimate) operate on the subject and pick the appropriate source automatically.

- **Y-DNA tab** — donor Y consensus (haplogroup + lineage + confidence) as the headline; the
  per-source calls that fed it listed below; private-Y as the union across Y-bearing sources,
  re-classified against the consensus backbone. No alignment picking.
- **mtDNA tab** — donor mt consensus; contributing sources; heteroplasmy aggregated.
- **Ancestry tab** — the donor estimate. Default policy: use the **best autosomal source** (highest
  mean-coverage WGS); later, pool genotypes across sources. Show which source it came from.
- **IBD tab** — uses the donor's best genotyped source.
- **Data Sources tab** — unchanged: the per-run/per-panel detail view (the one place run-level
  belongs), where you add data and inspect each source's coverage/calls.

### The "consensus call pool"

A donor's reconciled view is built from the pool of all calls across their sources:
- **Haplogroup calls** (`haplogroup_call`, keyed by `(biosample, dna_type, source)`) → `reconcile`
  → consensus. Already implemented; just surface it as the tab headline.
- **Ancestry**: pick-best now; pool-genotypes later (union the panel genotypes across sources,
  then estimate once). A `biosample`-keyed `ancestry_result` (best/aggregate) alongside the
  existing alignment-keyed rows.
- **STR**: a consensus marker profile = per-marker best call across panels (highest-quality /
  modal repeat), with disagreements flagged.
- **Private-Y**: union of off-backbone calls across Y sources, masked + classified against the
  consensus terminal.

## UI / data-flow changes

1. Each subject gets a **default source per modality** (e.g. best-coverage WGS for autosomal/Y/mt;
   the most complete STR panel). Tabs use it for actions without the user navigating Data Sources.
2. Tabs **read the donor aggregate** (consensus / rollup), not `selected_alignment`. The
   `selected_alignment` becomes a Data-Sources-only concept.
3. Analysis events update the **donor** view (consensus reload), which they already trigger.

## Phasing

- **Phase 1 — default-source + consensus headline (small, no schema change).** Compute a default
  alignment per subject (highest mean coverage; the project report already has a notion of a
  "drive" alignment). Tabs use it for actions and show the donor consensus as the headline. Removes
  the "pick an alignment first" friction; Y/mt already reconcile.
- **Phase 2 — ancestry + STR donor rollups.** `biosample`-keyed ancestry (best source now; show
  provenance) and a consensus STR profile. Small store additions.
- **Phase 3 — full pooled calls.** Genotype-pool ancestry across sources; private-Y union; a single
  "consensus call pool" read path the tabs consume.

## Open decisions

- **Ancestry across sources**: pick-best vs genotype-pooling. Pooling is more correct for a donor
  with several partial sources but needs careful site-union + dedup; pick-best is trivial and fine
  when one WGS dominates. Recommend pick-best in P2, pooling in P3.
- **STR consensus rule**: modal repeat vs highest-quality-source per marker; how to surface
  disagreements (a curator-judgable conflict, like haplogroup reconciliation).
- **Default-source policy**: pure max-coverage, or test-type aware (e.g. prefer a Big-Y for the Y
  tab, an mtDNA full sequence for the mt tab, WGS for ancestry).

## Non-goals

Coverage/sequencing QC stays per-run (it's intrinsically a run property). Data Sources stays the
run-level workspace. This is about the *analysis tabs* presenting the donor, not removing per-run
detail.
