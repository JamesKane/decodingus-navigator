# Multi-source reconciliation — Y & mtDNA across runs + Sanger confirmation

Status: **built** (all 6 phases landed on `rust-rewrite`). Scope: combine a donor's uniparental (Y, mtDNA)
results from multiple sources — sequencing runs across platforms (short-read, HiFi),
chip/array, STR panels, and **supplemental Sanger confirmations** — into a donor-level
consensus with conflict detection, source weighting, identity verification, and a manual
override / audit trail.

This is the natural consequence of supporting >1 run per donor (e.g. WGS229 short-read +
GFX0457637 HiFi are the same individual). Today the Rust app assigns haplogroups
*per alignment / per mtDNA sequence*; nothing combines them at the donor level.

## The Scala reference (what exists today)

The Scala app implements this fully — this design ports it, **simplified** per Rewrite
Plan §6 ("7-table Y-profile schema + manual concordance → one `YDnaProfile {snps, strs,
sources, reconciliation}` aggregate"). Key pieces to preserve:

- **Source types with quality weights** (`YProfileSourceType`): `SANGER` (SNP 1.0 — gold
  standard; STR 0.9), `CAPILLARY_ELECTROPHORESIS` (STR 1.0, SNP 0.5), `WGS_LONG_READ`
  (0.95/0.7), `WGS_SHORT_READ` (0.85/0.5), `TARGETED_NGS`, `CHIP` (0.5/0.3), `MANUAL`
  (0.3/0.2). Final per-call weight = method weight × callable-state weight (`CALLABLE` 1.0,
  `LOW_COVERAGE` 0.5, `EXCESSIVE`/`POOR_MAPQ` 0.3, `NO_COVERAGE`/`REF_N` 0.0).
- **Variant concordance** (`YVariantConcordance`): quality-weighted voting per position →
  consensus allele + state (`DERIVED`/`ANCESTRAL`/`HETEROPLASMY`/`NO_CALL`); a position is
  a **CONFLICT** when >30% of weighted votes disagree. Per-variant status: `CONFIRMED`
  (concordant, in tree), `NOVEL` (private — branch candidate), `CONFLICT`, `NO_COVERAGE`.
- **Haplogroup reconciliation** (`HaplogroupReconciliation`, per biosample per `DnaType`):
  `runCalls: [RunHaplogroupCall]` (per source: haplogroup, confidence, method, score,
  supporting/conflicting SNPs, technology, mean coverage, tree provider+version, lineage),
  `snpConflicts`, `heteroplasmyObservations` (mtDNA), `identityVerification`,
  `manualOverride`, `auditLog`, and a `ReconciliationStatus` (consensus haplogroup,
  `CompatibilityLevel` = COMPATIBLE / MINOR_DIVERGENCE / MAJOR_DIVERGENCE / INCOMPATIBLE,
  confidence, divergence point, SNP concordance, run count, warnings).
- **Identity verification** (`IdentityVerification`): are these runs the *same individual*?
  via autosomal kinship coefficient, Y-STR distance, and/or fingerprint-SNP concordance →
  VERIFIED_SAME … VERIFIED_DIFFERENT. Gates whether runs may be merged.
- These already map to the **Atmosphere Lexicon** (`com.decodingus.atmosphere.haplogroup
  Reconciliation`, `…#runHaplogroupCall`, `…#reconciliationStatus`, `…#heteroplasmy
  Observation`, `…#identityVerification`), so the consensus is PDS-syncable.

## What the rewrite already has to build on

- **Per-source haplogroup calls** — `assign_y_haplogroup` / `assign_mtdna_haplogroup*`
  return a `HaploAssignment` (terminal, score, lineage, matched/expected, branch evidence).
  Each run/source produces exactly the `RunHaplogroupCall` fields.
- **Variant evidence** — `private_y_variants` already classifies de-novo calls as
  off-path-known / novel; `call_bases_at` gives per-position bases. These are the
  per-source variant inputs to concordance.
- **Callable weighting** — self-referential `callable_intervals` (§4e) gives the per-source,
  per-position callable state that scales each source's concordance weight.
- **Identity verification, mostly free** — the IBD engine computes kinship/relationship
  (autosomal), and STR profiles give Y-STR distance. We can answer "same donor?" already.
- **STR profiles, SNP variant sets, chip profiles, mtDNA sequences** — the per-subject
  imports are all in place as sources.

## Proposed Rust design (simplified)

One aggregate per `(biosample, dna_type)`, not the 7-table Y-profile:

```
Reconciliation {
  biosample_guid, dna_type: Y | Mt,
  run_calls: [RunHaplogroupCall],          // one per contributing source
  consensus: { haplogroup, lineage, confidence,
               compatibility: Compatible | MinorDivergence | MajorDivergence | Incompatible,
               divergence_point, snp_concordance, run_count, warnings },
  conflicts: [VariantConflict],            // positions where sources disagree
  heteroplasmy: [HeteroplasmyObs],         // mtDNA only
  identity: Option<IdentityVerification>,  // are the runs the same donor?
  manual_override: Option<ManualOverride>,
  audit: [AuditEntry],
}
```

Sources (`SourceType` with `(snp_weight, str_weight)`, ported incl. **Sanger**), a
`Source` row referencing the originating run/profile + its type, and per-source variant
calls feeding concordance.

**Algorithms (all pure, testable):**
1. **Haplogroup consensus** — combine `run_calls` by tree topology: two calls are
   *compatible* iff one lineage is a prefix of the other (same branch, different depth);
   the consensus is the **deepest call supported by sufficient weighted agreement**, the
   divergence point is the LCA on disagreement, and `CompatibilityLevel` follows from how
   deep the LCA sits (tip vs backbone vs root). Confidence = weighted agreement fraction.
   (We already have lineage paths + the tree, so this is graph logic over existing data.)
2. **Variant concordance** — port `YVariantConcordance`: per position, weighted vote
   (source weight × callable-state weight) → consensus state; CONFLICT if >30% weighted
   disagreement; status from tree membership (CONFIRMED/NOVEL/CONFLICT/NO_COVERAGE). Sanger
   calls, at weight 1.0, decisively confirm or overturn an NGS call.
3. **Identity verification** — reuse IBD kinship + Y-STR distance to set
   VERIFIED_SAME…DIFFERENT; refuse to merge (or warn loudly) below a threshold.
4. **mtDNA heteroplasmy** — per position, major/minor allele + frequency from the
   per-source allele fractions we already compute.

**Sanger confirmation flow:** a manual `Source` of type `SANGER` carrying a small set of
confirmed variant calls (position, allele, derived/ancestral) entered via the UI. Because
its SNP weight is 1.0, it dominates concordance at those positions — confirming a tentative
NGS branch SNP or flagging an NGS error. No new calling; it's a high-trust manual source.

**Storage:** one `reconciliation` row (consensus + identity + override as columns/JSON);
`reconciliation_source` and `reconciliation_call` child tables (or a single JSON blob for
the edge app, matching the Atmosphere record shape for sync). Recompute on any source
add/remove; append an `AuditEntry` each time (INITIAL / RUN_ADDED / RUN_REMOVED /
MANUAL_OVERRIDE / CONFLICT_RESOLVED / RECOMPUTED).

**Atmosphere alignment:** the aggregate serializes to the existing
`com.decodingus.atmosphere.haplogroupReconciliation` record (floats-as-strings per the
no-float rule) so it publishes through the same `navigator-sync` path as coverage/private
variants.

## Phasing (build order) — ✅ all phases complete

1. ✅ **Source model + per-source run calls.** `SourceType` (+ weights, incl. Sanger), a
   `Source` table, and capture each `assign_*` result as a `RunHaplogroupCall`.
2. ✅ **Haplogroup consensus + compatibility.** Pure topology combine over `run_calls`;
   surface consensus + divergence in the subject view. (Highest value, smallest lift —
   directly answers "combine the two runs.")
3. ✅ **Identity verification.** Wire IBD kinship + Y-STR distance; gate/flag merges.
4. ✅ **Variant concordance + status.** Weighted voting; CONFIRMED/NOVEL/CONFLICT; feeds the
   private/branch-proposal bucket with multi-source confidence.
5. ✅ **Sanger / manual source entry.** UI to add confirmed calls at weight 1.0.
6. ✅ **mtDNA heteroplasmy + manual override + audit + PDS record.** Heteroplasmy via a chrM
   pileup scan (`navigator_analysis::heteroplasmy`); override + audit persisted in
   `reconciliation_override`/`reconciliation_audit` (migration 0010); the aggregate
   publishes as `com.decodingus.atmosphere.haplogroupReconciliation` (floats-as-strings).
   UI: per-DNA-type override/clear, audit-log expander, heteroplasmy scan, publish button.

## Open questions

- **Autosomal kinship across runs** needs autosomal genotyping per run (the panel
  genotyper exists, but identity verification wants a shared fingerprint-SNP set) — define
  the fingerprint panel, or rely on Y-STR + uniparental concordance when autosomal is thin.
- **Consensus when sources truly diverge** (INCOMPATIBLE) — surface as "likely different
  individuals," do not merge; this is the identity-verification guardrail.
- Keep the simplified aggregate from regrowing into the 7-table Scala shape — resist
  per-region/per-variant tables; prefer a compact record + child calls.
