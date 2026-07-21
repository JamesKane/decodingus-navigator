# External-caller precedence + autosomal sidecar fast path

Status: **Phases 1–3 landed** (provenance + reconcile precedence + guards + backfill migration +
`reingest-external` CLI; observation-pooling now GVCF-sourced; Preferences toggle + "Compare
callers" diagnostic). Phases 4–5 (autosomal ingest) not started.
Scope: `navigator-store` (provenance/migration), `navigator-analysis` (call-set reader,
diploid gVCF), `navigator-app` (reconcile precedence, ingest, guards, observation pooling),
`navigator-ui` (Preferences toggle, "Compare callers" action).
Supersedes the "additive haplogroups" assumption in
[`pipeline-artifact-import.md`](pipeline-artifact-import.md) §4, which does not actually hold.

## 1. Motivation

Advanced users run an established external calling workflow (GATK4 for chrY/chrM GVCFs today;
GATK4 / `pileupCaller` for autosomal). Navigator already ingests the Y/mt GVCFs through the
sidecar fast path. Two problems block treating that external caller as authoritative:

1. **Navigator's internal caller silently overwrites the external call.** Observed on the
   ancient-DNA collection **PRJEB37976**: re-running analysis changed placed haplogroups,
   because the damaged-CRAM walk replaced the clean GATK4 GVCF placement. This is a
   correctness bug, not a preference gap.
2. **There is no autosomal fast path.** Autosomal ancestry/IBD still require a full CRAM
   decode even when the user already has an external autosomal call set — and for ancient
   DNA the correct representation (pseudo-haploid **1240K**, the AADR/qpAdm space) is exactly
   what their pipeline emits and what our validated deep-ancestry estimator consumes.

This design (a) makes external calls first-class and never-overwritten, with a user preference
to make them authoritative, and (b) extends the sidecar fast path to autosomal via a targeted
1240K call set, reusing the existing 1240K frontend/backend.

## 2. Root-cause analysis — the overwrite

Y/mt haplogroup calls persist in `haplogroup_call`, keyed
`UNIQUE(biosample_guid, dna_type, source_key)`
(`navigator-store/migrations/0008_haplogroup_calls.up.sql`; `source_fingerprint` added by 0016).
The writer `haplogroup_call::upsert` (`navigator-store/src/haplogroup_call.rs:40`) is
`INSERT … ON CONFLICT(biosample_guid, dna_type, source_key) DO UPDATE SET haplogroup = excluded.haplogroup, …`.

Four compounding defects:

1. **Shared `source_key`.** The sidecar path (`fastpath.rs:188 assign_y_from_gvcf` →
   `aln:{id}`; `fastpath.rs:226 assign_mt_from_gvcf` → `aln:{id}:mt`) and the internal walk
   (`haplogroup.rs:3080 assign_y_haplogroup` → `aln:{id}`; `haplogroup.rs:2471
   assign_mtdna_haplogroup_from_alignment` → `aln:{id}:mt`) write the **same** `source_key`.
   The upsert clobbers.
2. **The skip-guard cannot fire across the boundary.** The fast path stamps a `gv:`-prefixed
   fingerprint (`fastpath.rs:178 gvcf_fingerprint`); the walk compares against an `f:`-prefixed
   `y_score_fingerprint` (`haplogroup.rs:3099-3109`). `stored == fp` is never true → the walk
   always re-genotypes and overwrites (`record_call_fp`, `haplogroup.rs:3113`).
3. **`RunFullAnalysis` is unguarded.** `run_full_analysis_streaming`
   (`navigator-ui/src/worker.rs:2414`) unconditionally enqueues `AssignYHaplogroup`
   (worker.rs:2557) + `AssignMtdnaHaplogroupFromAlignment` (worker.rs:2562). This is the
   primary trigger. Batch `analyze_biosample` (`queries.rs:666`) is *partially* guarded for Y
   (`haplogroup_consensus(...).is_some()`, queries.rs:761) but has **no mt step/guard at all**.
4. **Reconcile is source-blind.** `reconciliation::reconcile`
   (`navigator-domain/src/reconciliation.rs:179`) returns the highest-`score` call. Ancient
   DNA deamination (C→T / G→A) manufactures false-derived alleles, so a damaged walk can
   out-score the clean external call even where both rows survive.

`haplogroup_call` has **no** `source`/`provenance`/`locked`/`priority` column. The only
provenance-aware protection in the codebase — `save_analysis_no_downgrade` /
`completeness_rank` (`commands.rs:462`) over `analysis_artifact.source/completeness`
(migration 0017) — covers coverage/read-metrics/sex/SV and is **bypassed** by haplogroups.

## 3. Decisions (locked)

- **Precedence control = one global setting.** `AppSettings.prefer_external_calls: bool`
  (Preferences, alongside the haplogroup-tree provider). When on, external calls win the
  consensus wherever present. Default: **on** (the safe direction — never dilute a call the
  user deliberately produced).
- **Row-level no-clobber is always on**, independent of the setting: external and internal
  calls are distinct rows; neither ever silently overwrites the other.
- **When external is preferred, the internal walk is skipped entirely** for that modality — no
  CRAM decode. An explicit **"Compare callers"** action runs the walk on demand to audit
  divergence.
- **Autosomal external ingest, first increment = a 1240K call set** (EIGENSTRAT
  `.geno/.snp/.ind` and/or PLINK; `pileupCaller` pseudo-haploid). Autosomal gVCF (GATK4
  diploid) is the second increment, same dosage sink.

## 4. Part 1 — call provenance & precedence (Y / mt / autosomal, unified)

### 4.1 Provenance dimension

A single provenance enum reused everywhere calls are pooled:

```rust
// navigator-domain
pub enum CallProvenance { Manual, External, NavigatorWalk }
```

Precedence (highest first): **Manual > External > NavigatorWalk** when
`prefer_external_calls`; when off, External and Walk are peers ranked by `score`. Manual (the
existing `override_consensus`, `haplogroup.rs:390`) always wins. External is **never**
demoted below Walk — the toggle only chooses between "external auto-wins" and "score breaks
the tie"; it never lets a walk silently supersede an external row.

### 4.2 Storage (`navigator-store`)

Migration `00NN_haplogroup_call_provenance`:
- Add `provenance TEXT NOT NULL DEFAULT 'navigator-walk'` to `haplogroup_call`.
- **Distinct `source_key` per provenance** so the upsert can never cross the boundary:
  external → `aln:{id}:ext` / `aln:{id}:ext:mt`; walk → `aln:{id}:walk` / `aln:{id}:walk:mt`;
  manual unchanged. (`source_key` distinctness is what actually stops the clobber; the
  `provenance` column is the typed key reconcile ranks on, cheaper than parsing the string.)
- **Backfill** existing rows (§7): rows with a `gv:` fingerprint → `provenance='external'`,
  `source_key` re-suffixed `:ext`; `f:` → `'navigator-walk'`, `:walk`. This heals the
  PRJEB37976 workspace: any external call that a walk had already overwritten is *gone* and
  must be re-ingested (the fast path re-runs cheaply from the GVCF), but future re-analysis
  can no longer clobber.

### 4.3 Precedence-aware reconcile

`reconciliation::reconcile` gains provenance + a `prefer_external: bool` policy. It first
partitions by tier; if `prefer_external`, it reconciles **within the top non-empty tier**
(Manual, else External, else Walk) and reports lower tiers only as a divergence warning
(reusing the existing "per-run calls vary" warning surface). Within a tier, the current
score/path logic (reconciliation.rs:187-212) is unchanged. `haplogroup_consensus`
(`haplogroup.rs:285`) passes the setting through.

### 4.4 Guards — no walk when external is preferred

The single behavior "skip the internal walk for a modality that has a preferred external
call" must hold at **every** entry point (the agents found three, only one partially guarded):

- Add `has_preferred_external_call(biosample, dna_type) -> bool` (checks for an `External` row
  when `prefer_external_calls`).
- `assign_y_haplogroup` / `assign_mtdna_haplogroup_from_alignment`: early-return the existing
  external call instead of walking, when preferred. (Also add the missing mt guard to
  `analyze_biosample`.)
- `run_full_analysis_streaming` (worker.rs:2557-2566): gate the `AssignY*`/`AssignMtdna*`
  enqueues on `!has_preferred_external_call`. This closes the primary PRJEB37976 trigger.
- The walk writer additionally refuses to downgrade: a walk `record_call_fp` never writes onto
  a higher-tier `source_key` (belt-and-suspenders mirror of `completeness_rank`).

### 4.5 Observation pooling must respect provenance

`place_y_consensus` / `place_mt_consensus` (`haplogroup.rs:821` / mt sibling) pool genotypes
**by position** across every WGS alignment. For an alignment with a preferred external call,
the pooled placement must take that alignment's observations from the **GVCF** (the genotype
the fast path already cached — `place_*_consensus` already "reuses the exact genotype the
assignment cached"), never re-decoding the CRAM. Without this, aDNA damage in the CRAM still
dilutes the pooled `ObservedProfile` even though the row-level call is protected. Concretely:
`place_*_consensus`'s per-alignment genotype source branches on
`has_preferred_external_call` → `gvcf_base_calls` (fastpath.rs:62) vs `base_calls`.

## 5. Part 2 — autosomal external ingest (1240K first)

The 1240K frontend/backend already exists (deep-ancestry §7.16-7.18): `IbdPanel` is the
full-1240K multi-build genotyping frontend (`resolve_alignment`/`resolve_chip` re-key any
build/chip to canonical CHM13); `reconcile_diploid` pools **per-source dosages** into the
`DiploidProfile` consensus and is an incremental reducer; the consensus feeds modern/fine/deep
(qpAdm) ancestry and IBD. The extension is one new "resolve external call set → per-source
1240K panel dosages" path that feeds `reconcile_diploid` as an `External` source — parallel to
the Y/mt GVCF sidecar, **no CRAM decode**.

### 5.1 Call-set reader (`navigator-analysis`, new `callset.rs`)

Input: a targeted 1240K call set. First formats: **EIGENSTRAT** (`.geno`/`.snp`/`.ind`) and
**PLINK** (`.bed`/`.bim`/`.fam`); both carry (build, pos, ref/alt or rsID) + genotype.
Output: `Vec<(canonical_site, PanelDosage)>` in the exact shape `ibd_panel_dosages` produces.

- Map each call-set SNP to a 1240K panel locus by **rsID** and by **(build, pos)** (the panel
  carries per-build loci + rsIDs — same join deep-ancestry uses). Orientation-check against the
  panel's canonical CHM13 allele (the §7.16 orientation fix): flip dosage where ref/alt reversed;
  drop sites matching neither allele.
- **Pseudo-haploid** (`pileupCaller`; `.geno` values 0/2 only) → dosage {0,2}; carry a
  `pseudo_haploid` flag on the source so `reconcile_diploid` / qpAdm treat it correctly (do not
  synthesize hets). **Diploid** call sets → {0,1,2} directly.

### 5.2 Sidecar discovery + ingest wiring

- Extend `SampleSidecars` (`scan.rs:69`) + `detect_sidecars` (`scan.rs:105`): add
  `autosomal_callset: Option<CallSetRef>` matched on an EIGENSTRAT/PLINK triplet or a
  `*.1240k.*` name. (Autosomal gVCF `*.autosome*.g.vcf.gz` reserved for increment 2.)
  Extend `has_haplogroup_gvcf`'s sibling gate with `has_autosomal_callset`.
- `ingest_sidecars` (`fastpath.rs:256`) gains an autosomal step: call-set → panel dosages →
  persist as a per-source dosage row (`External` provenance) → `refresh_autosomal_consensus`
  (the §7.18 reducer, which never decodes an uncached alignment). Cheap; folds into the
  existing `FastPathSummary`.

### 5.3 Provenance in the diploid consensus

`reconcile_diploid` already pools per-source dosages; tag each source with `CallProvenance`.
Under `prefer_external_calls`, an alignment's **CRAM-genotyped** panel dosages are not computed
(and if cached, not pooled) for a subject that has an `External` autosomal source — same
skip-when-preferred rule as Y/mt, so a damaged aDNA CRAM never dilutes the external 1240K
consensus. Modern/fine/deep ancestry and IBD read the consensus unchanged.

## 6. "Compare callers" action

An explicit per-subject action (report button + CLI `compare-callers`) that, ignoring the
preference, runs the internal walk into its `:walk` rows / CRAM panel dosages and renders the
divergence (external terminal vs Navigator terminal; per-site autosomal dosage disagreement).
This is the opt-in audit that replaces the automatic secondary walk we chose *not* to run.

## 7. Migration / backfill

- The provenance migration backfills existing `haplogroup_call` rows by fingerprint prefix
  (§4.2). Rows already clobbered by a walk are lost and re-ingested from the GVCF.
- Add `navigator reingest-external <project>` (or fold into `rebuild-signatures`) to re-run the
  sidecar fast path over a collection, restoring external calls the old code had overwritten —
  the operational fix for the current PRJEB37976 state.

## 8. Phasing

1. **Provenance storage + reconcile precedence + guards** (Part 1, no autosomal). Fixes the
   PRJEB37976 overwrite; independently shippable. Backfill migration + reingest command. **DONE** —
   see §10.
2. **Observation-pooling provenance** (§4.5) — pooled placement stops re-walking preferred-
   external alignments. **DONE** — see §11.
3. **Preferences toggle + "Compare callers"** UI (§3, §6). **DONE** — see §12.
4. **Autosomal call-set reader + ingest** (Part 2, 1240K first) → consensus → ancestry/IBD.
5. **Autosomal gVCF** (GATK4 diploid) as a second call-set source into the same sink.

Phase 1 alone resolves the reported bug; 4 delivers the autosomal fast path.

## 9. Validation gates

- **PRJEB37976 idempotence (the reported bug).** Ingest external Y/mt → run Full Analysis and
  the deep pass → assert the placed haplogroups are **unchanged** (external row survives, no
  walk row wins). With `prefer_external_calls=false`, assert both rows exist and the divergence
  warning renders.
- **No-clobber unit test.** Sidecar then walk on the same alignment → two distinct rows; walk
  never mutates the `:ext` row.
- **Observation pooling.** A preferred-external subject's pooled placement is byte-identical
  whether or not the CRAM is present (no re-walk dilution).
- **Autosomal parity.** A subject genotyped from an external 1240K call set vs the internal
  CRAM panel genotyping agree at shared sites (reuse the 99.84% same-person concordance harness
  from deep-ancestry §7.13); deep qpAdm from the external-sourced consensus reproduces the WGS
  fit within the stability band.
- **aDNA end-to-end.** A PRJEB37976 sample with an external 1240K call set produces a deep-
  ancestry/consensus result with **no CRAM decode**, and the internal walk (via "Compare
  callers") shows the expected deamination-driven divergence.

## 10. Phase 1 — as built

Landed on `main` (working tree); `cargo clippy --all-targets -- -D warnings` clean.

- **Provenance model** — `navigator_domain::reconciliation::CallProvenance`
  {`NavigatorWalk`,`External`,`Manual`} (`as_str`/`from_token`/`rank`) +
  `reconcile_with_provenance(calls, prefer_external)`: under the policy, only the top tier present
  reconciles; lower tiers that differ become a warning. Source-blind (== `reconcile`) when off.
  3 unit tests.
- **Storage** — migration `0036_haplogroup_call_provenance`: adds `provenance` (default
  `navigator-walk`), stamps existing `gv:`-fingerprint rows `external`, and re-keys them
  (`aln:{id}` → `aln:{id}:ext`, `…:mt` → `…:ext:mt`) so a walk can never share their key.
  `haplogroup_call::upsert` takes `CallProvenance`; new `list_for_with_provenance`; `list_all`
  carries it.
- **Writers** — the fast path (`assign_{y,mt}_from_gvcf`) records `External` on the `:ext` keys via
  `external_{y,mt}_source_key`; the CRAM walk records `NavigatorWalk` on the unchanged `aln:{id}` /
  `aln:{id}:mt` keys. External and internal now coexist as distinct rows — no upsert can clobber.
- **Consensus** — `haplogroup_consensus` + `haplogroup_terminals` use `reconcile_with_provenance`
  and **skip the CRAM-pooled placed label** (`consensus_profile.consensus_label`) on a
  preferred-external subject, so the external terminal is authoritative there. (Provenance-aware
  *placement* — sourcing `place_{y,mt}_consensus` from the GVCF instead of walking — is Phase 2.)
- **Guards** — `assign_y_haplogroup` / `assign_mtdna_haplogroup_from_alignment` early-return the
  external call (via `preferred_external_call`) instead of walking; the UI worker's Full-Analysis
  step-list skips the Y/mt enqueues when `has_preferred_external_call` (belt-and-suspenders over the
  assign-level guard). `analyze_biosample`'s Y step was already consensus-guarded; it has no mt step.
- **Setting** — `AppSettings.prefer_external_calls: Option<bool>` (+ env
  `NAVIGATOR_PREFER_EXTERNAL_CALLS`), resolver `prefer_external_calls()`, **default on**. The
  Settings modal preserves it (no toggle UI yet — Phase 3).
- **Operational fix** — `App::reingest_external_for_biosample` re-runs the sidecar fast path from the
  still-present GVCFs (no CRAM); CLI `navigator reingest-external [--project] [--db]` restores
  external calls a pre-provenance build had overwritten in a live workspace.
- **Tests** — domain `reconcile_with_provenance` (×3); store migration+upsert round-trip; app
  `external_call_is_not_clobbered_and_wins_consensus` (the PRJEB37976 idempotence gate).

## 11. Phase 2 — as built

Landed on `feat/external-caller-precedence`; `cargo clippy --all-targets -- -D warnings` clean; the
app suite (lib + integration) still green (the fallback path is behavior-preserving).

The genome-consensus placement no longer walks the CRAM for a preferred-external alignment — it
sources that alignment's tree-locus genotype from the sidecar GVCF, so the pooled placement, the Y/mt
variant profile, and the descent report are all GVCF-derived (no ancient-DNA dilution, no decode).

- **`App::consensus_base_calls(aln, contig, tree, tree_source_build)`** — drop-in for the per-
  alignment genotype in the consensus placements. When `prefer_external_calls()` and the alignment
  has a sidecar GVCF (`chr_{y,m}_gvcf_for_alignment`), it calls the existing (tested) `gvcf_base_calls`
  (native read, or lift for a non-native tree build); otherwise the cached CRAM walk `base_calls`.
  Same signature/shape as `base_calls`, so behavior is identical whenever no GVCF is present.
- **Wired at all three consensus call sites:** `place_y_consensus_decodingus` (native, `None` build),
  `mt_source_calls` (rCRS, `None` build), and `place_y_consensus_ftdna` (lifts native→GRCh38 via
  `tree_build_for_contig("chrY")`, falling back to `assign_haplogroup_detail` when no GVCF — its
  prior path is untouched).
- **Added** `chr_m_gvcf_for_alignment` (the chrM sibling finder / `NAVIGATOR_M_GVCF`), mirroring the
  existing chrY one.

**Interaction with the Phase 1 skip.** `haplogroup_consensus` still skips the placed
`consensus_label` on a preferred-external subject — now not for correctness (a freshly built label is
GVCF-sourced and agrees with the external call) but so a **stale** label from a pre-Phase-2 (CRAM-
pooled) build cannot resurface before the profile is rebuilt. After `rebuild-signatures` the placed
label and the external reconcile agree.

**Validation.** Fallback path unchanged → existing consensus/placement tests pass. GVCF-vs-CRAM
parity for the substituted genotype is the existing `gvcf_fast_path_matches_cram_walk` /
`gvcf_y_placement_smoke` gates (same `gvcf_base_calls` the placement now calls). A placement-level
parity gate (place-with-GVCF == place-with-CRAM on the running fixture) is the natural addition when a
committed GVCF+ref+tree fixture is added.

**Not in Phases 1–2:** autosomal external ingest (§5, Phases 4–5), the Preferences toggle +
"Compare callers" UI (§6, Phase 3).

## 12. Phase 3 — as built

Landed on `feat/external-caller-precedence`; clippy clean; i18n parity + app suite green.

- **Preferences toggle** — a checkbox in the Settings modal ("Prefer external caller (GATK4 /
  1240K)") backed by `SettingsForm.prefer_external_calls` → `AppSettings.prefer_external_calls`
  (default on). New i18n keys `settings.preferExternalCalls[/Hint]` in `en`+`es` (parity test green).
  The env override `NAVIGATOR_PREFER_EXTERNAL_CALLS` still wins over the setting.
- **"Compare callers" diagnostic** — `App::compare_callers(alignment_id) -> Vec<CallerComparison>`
  (`{dna_type, external, navigator}`, `.agree()`), which **forces** the internal walk regardless of
  the policy. The walk was extracted into `assign_y_haplogroup_walk` / `assign_mtdna_haplogroup_walk`
  (the guarded public `assign_*` now = guard → walk), so compare can run it without the guard; the
  walk records its own `aln:{id}` / `:mt` (`NavigatorWalk`) rows and never touches the external `:ext`
  row (non-destructive). CLI `navigator compare-callers <subject>` prints, per alignment, the external
  vs Navigator Y/mt terminal and flags real divergences.
- **Not done (optional):** a GUI "Compare callers" button (worker `Command`/`Event` + detail-view
  rendering). The CLI + app method deliver the capability; the button is a thin follow-up.

**Not in Phases 1–3:** autosomal external ingest — the 1240K EIGENSTRAT/pileupCaller call-set reader
→ `reconcile_diploid` → consensus → ancestry/IBD (§5, Phases 4–5).
```
