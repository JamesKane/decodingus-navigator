# Atmosphere Lexicon Alignment Plan

Navigator Desktop model updates to align with the Web Application's Atmosphere Lexicon
(`documents/atmosphere/`). This plan covers every model, field, file, repository, codec,
UI component, and migration required.

**Date:** 2026-03-07
**Baseline:** Web App Lexicon v1.9 (2025-12-09)

---

## Summary

| Change Area | Scope | Files Affected | Priority |
|:---|:---|---:|:---|
| 1. ChipProfile field additions + vendor/provider rename | Small | ~12 | High |
| 2. AncestryResult -> PopulationBreakdown promotion | Medium | ~14 | High |
| 3. HaplogroupReconciliation enrichment | Small | ~6 | Medium |
| 4. Sync-time validation & test type mapping | Small | ~3 | Medium |

---

## 1. ChipProfile Field Additions & Vendor/Provider Rename

**Goal:** Align `ChipProfile` with `com.decodingus.atmosphere.genotype` schema.

### 1.1 Field Changes to `ChipProfile`

| Field | Current | Target | Action |
|:---|:---|:---|:---|
| `vendor` | `String` | rename to `provider` | Rename field |
| `yMarkersTotal` | missing | `Option[Int]` | Add |
| `mtMarkersCalled` | missing | `Option[Int]` | Add |
| `mtMarkersTotal` | missing | `Option[Int]` | Add |
| `testDate` | missing | `Option[LocalDateTime]` | Add |
| `processedAt` | missing | `Option[LocalDateTime]` | Add |
| `buildVersion` | missing | `Option[String]` | Add (GRCh37/GRCh38) |
| `derivedHaplogroups` | missing | `Option[HaplogroupAssignments]` | Add |
| `populationBreakdownRef` | missing | `Option[String]` | Add |
| `imputationRef` | missing | `Option[String]` | Add |
| `sourceFileName` | `Option[String]` | keep | Navigator extension (not in web schema) |

### 1.2 Files to Modify

**Model:**
- `src/main/scala/com/decodingus/workspace/model/ChipProfile.scala`
  - Rename `vendor: String` -> `provider: String`
  - Add new optional fields listed above
  - Update `KnownVendors` -> `KnownProviders` (same values)
  - Update companion object methods that reference `vendor`

**Repository:**
- `src/main/scala/com/decodingus/repository/ChipProfileRepository.scala`
  - Rename `vendor` column -> `provider` in `ChipProfileEntity`
  - Add columns: `y_markers_total`, `mt_markers_called`, `mt_markers_total`,
    `test_date`, `processed_at`, `build_version`, `derived_haplogroups` (JSONB),
    `population_breakdown_ref`, `imputation_ref`
  - Update `ChipProfileCodecs`
  - Update `findByVendor()` -> `findByProvider()`

**Service:**
- `src/main/scala/com/decodingus/service/EntityConversions.scala`
  - Update `toChipProfileEntity()` and `fromChipProfileEntity()` for renamed/new fields

- `src/main/scala/com/decodingus/service/H2WorkspaceService.scala`
  - No structural changes; flows through EntityConversions

- `src/main/scala/com/decodingus/workspace/services/WorkspaceOperations.scala`
  - No structural changes

**Parsers:**
- `src/main/scala/com/decodingus/genotype/parser/ChipDataParser.scala`
  - Update all parser implementations: `vendor = "23andMe"` -> `provider = "23andMe"`, etc.
  - All 6 parsers (Parser23andMe, ParserAncestryDna, ParserFtdna, ParserMyHeritage,
    ParserLivingDna, ParserBisdna)

**UI:**
- `src/main/scala/com/decodingus/ui/v2/SubjectDetailView.scala`
  - Update references from `chip.vendor` -> `chip.provider`
  - Update dialog display labels

**ViewModel:**
- `src/main/scala/com/decodingus/workspace/WorkbenchViewModel.scala`
  - Update `parseChipFile` method: field name change and populate new fields
    (`processedAt = Some(LocalDateTime.now())`, `buildVersion` from parser detection)

**Codecs:**
- `src/main/scala/com/decodingus/pds/PdsClient.scala`
  - Codecs use `deriveEncoder`/`deriveDecoder` so they auto-update with field rename
  - Verify serialized JSON field name is `provider` (not `vendor`)

**Sync:**
- `src/main/scala/com/decodingus/repository/SyncQueueRepository.scala`
  - No changes (uses EntityType enum, not field names)

**Tests:**
- `src/test/scala/com/decodingus/service/Phase2EntityIntegrationSpec.scala`
  - Update test fixtures: `vendor` -> `provider`
- `src/test/scala/com/decodingus/repository/ChipProfileRepositorySpec.scala` (if exists)
  - Update test fixtures

**Database Migration:**
- New migration: `V0XX__align_chip_profile_with_lexicon.sql`
  - `ALTER TABLE chip_profiles RENAME COLUMN vendor TO provider;`
  - `ALTER TABLE chip_profiles ADD COLUMN y_markers_total INTEGER;`
  - `ALTER TABLE chip_profiles ADD COLUMN mt_markers_called INTEGER;`
  - `ALTER TABLE chip_profiles ADD COLUMN mt_markers_total INTEGER;`
  - `ALTER TABLE chip_profiles ADD COLUMN test_date TIMESTAMP;`
  - `ALTER TABLE chip_profiles ADD COLUMN processed_at TIMESTAMP;`
  - `ALTER TABLE chip_profiles ADD COLUMN build_version VARCHAR(20);`
  - `ALTER TABLE chip_profiles ADD COLUMN derived_haplogroups CLOB;` (JSON)
  - `ALTER TABLE chip_profiles ADD COLUMN population_breakdown_ref VARCHAR(500);`
  - `ALTER TABLE chip_profiles ADD COLUMN imputation_ref VARCHAR(500);`

### 1.3 Behavioral Changes

- `ChipProfile.isAcceptableForAncestry()` - no change (uses marker counts, not vendor)
- `ChipProfile.hasSufficientYCoverage()` - can now use `yMarkersTotal` if available
  for more accurate threshold checking
- `ChipProfile.hasSufficientMtCoverage()` - can now use `mtMarkersTotal` similarly
- Chip parsers should populate `buildVersion` based on vendor detection (23andMe v5 = GRCh37,
  AncestryDNA v2 = GRCh37, etc.)

---

## 2. AncestryResult -> PopulationBreakdown Promotion

**Goal:** Promote `AncestryResult` from a standalone analysis result to a first-class
Atmosphere Lexicon record (`com.decodingus.atmosphere.populationBreakdown`) with AT URI,
RecordMeta, and biosample reference. Align field names with web app schema.

### 2.1 Model Changes

**New wrapper record: `PopulationBreakdown`**

```scala
// New file or wraps existing AncestryResult
case class PopulationBreakdown(
  atUri: Option[String],
  meta: RecordMeta,
  biosampleRef: String,
  analysisMethod: String,               // "PCA_PROJECTION_GMM", "ADMIXTURE", etc.
  panelType: String,                     // "aims" or "genome-wide"
  referencePopulations: Option[String],  // "1000G_HGDP_v1"
  snpsAnalyzed: Int,
  snpsWithGenotype: Int,
  snpsMissing: Int,
  confidenceLevel: Double,
  components: List[PopulationComponent], // renamed from percentages
  superPopulationSummary: List[SuperPopulationSummary], // renamed
  pcaCoordinates: Option[List[Double]],
  analysisDate: Option[LocalDateTime],
  pipelineVersion: Option[String],       // renamed from analysisVersion
  referenceVersion: Option[String]
)
```

**Rename `PopulationPercentage` -> `PopulationComponent`**

| Current Field | Target Field | Action |
|:---|:---|:---|
| `populationCode` | `populationCode` | Keep |
| `populationName` | `populationName` | Keep |
| (implicit from Population.scala) | `superPopulation` | Add explicitly |
| `percentage` | `percentage` | Keep |
| `confidenceLow` | `confidenceInterval.lower` | Restructure |
| `confidenceHigh` | `confidenceInterval.upper` | Restructure |
| `rank` | `rank` | Keep |

New nested type:
```scala
case class ConfidenceInterval(lower: Double, upper: Double)
```

**Rename `SuperPopulationPercentage` -> `SuperPopulationSummary`**

| Current Field | Target Field | Action |
|:---|:---|:---|
| `superPopulation` | `superPopulation` | Keep |
| `percentage` | `percentage` | Keep |
| `populations` | `populations` | Keep |

### 2.2 Files to Modify

**Model (core changes):**
- `src/main/scala/com/decodingus/ancestry/model/AncestryResult.scala`
  - Rename `PopulationPercentage` -> `PopulationComponent`
  - Add `superPopulation: String` field to `PopulationComponent`
  - Replace `confidenceLow`/`confidenceHigh` with `confidenceInterval: ConfidenceInterval`
  - Add `ConfidenceInterval` case class
  - Rename `SuperPopulationPercentage` -> `SuperPopulationSummary`
  - Rename `analysisVersion` -> `pipelineVersion`
  - Keep `AncestryResult` as the internal analysis return type (processor output)
  - Update `fromProbabilities()` factory to build `PopulationComponent` with new structure
  - Update `calculateCiWidth()` to produce `ConfidenceInterval`

- `src/main/scala/com/decodingus/workspace/model/PopulationBreakdown.scala` (NEW FILE)
  - `PopulationBreakdown` case class wrapping the analysis result with AT Protocol fields
  - Factory method: `fromAncestryResult(result: AncestryResult, biosampleRef: String)`
  - Re-export `PopulationComponent`, `SuperPopulationSummary`, `ConfidenceInterval`

- `src/main/scala/com/decodingus/workspace/model/Workspace.scala`
  - Add `populationBreakdowns: List[PopulationBreakdown]` to `WorkspaceContent`
  - Add `getPopulationBreakdownForBiosample()` method

**Processor pipeline (return type changes):**
- `src/main/scala/com/decodingus/ancestry/processor/AncestryEstimator.scala`
  - Update `estimate()` return: use renamed types (`PopulationComponent`, etc.)

- `src/main/scala/com/decodingus/ancestry/processor/AncestryProcessor.scala`
  - Return `Either[String, AncestryResult]` stays the same (processor output)
  - The wrapping into `PopulationBreakdown` happens at the workspace level

- `src/main/scala/com/decodingus/genotype/processor/ChipAncestryAdapter.scala`
  - Same: returns `AncestryResult`, wrapping happens at workspace level

**Reporting (field name updates):**
- `src/main/scala/com/decodingus/ancestry/report/AncestryReportWriter.scala`
  - Update field references: `percentages` -> `components`
  - Update `PopulationPercentage` -> `PopulationComponent`
  - Update `confidenceLow`/`confidenceHigh` -> `confidenceInterval.lower`/`.upper`
  - Update JSON output structure

**UI:**
- `src/main/scala/com/decodingus/ui/components/AncestryResultDialog.scala`
  - Update type references: `PopulationPercentage` -> `PopulationComponent`
  - Update confidence interval display code

- `src/main/scala/com/decodingus/ui/v2/SubjectDetailView.scala`
  - Update references to `PopulationPercentage` type

**ViewModel:**
- `src/main/scala/com/decodingus/workspace/WorkbenchViewModel.scala`
  - `analyzeChipDataForAncestry()`: Wrap `AncestryResult` in `PopulationBreakdown`
  - Create proper AT URI: `local:populationBreakdown:{sampleAccession}:{uuid}`
  - Save `PopulationBreakdown` record to workspace
  - Update `biosample.populationBreakdownRef` with the AT URI

**Persistence (new repository):**
- `src/main/scala/com/decodingus/repository/PopulationBreakdownRepository.scala` (NEW FILE)
  - `PopulationBreakdownEntity` with all fields
  - `components` stored as JSONB (List[PopulationComponent])
  - `superPopulationSummary` stored as JSONB
  - `pcaCoordinates` stored as JSONB
  - Standard CRUD + `findByBiosample()`

- `src/main/scala/com/decodingus/service/EntityConversions.scala`
  - Add `toPopulationBreakdownEntity()` / `fromPopulationBreakdownEntity()`

- `src/main/scala/com/decodingus/service/H2WorkspaceService.scala`
  - Add PopulationBreakdown CRUD methods
  - Load PopulationBreakdowns in workspace load

- `src/main/scala/com/decodingus/workspace/services/WorkspaceOperations.scala`
  - Add `addPopulationBreakdown()`, `getPopulationBreakdownForBiosample()`

**Codecs:**
- `src/main/scala/com/decodingus/pds/PdsClient.scala`
  - Add `PopulationBreakdown`, `PopulationComponent`, `SuperPopulationSummary`,
    `ConfidenceInterval` encoder/decoder

**Sync:**
- `src/main/scala/com/decodingus/repository/SyncQueueRepository.scala`
  - Add `PopulationBreakdown` to `EntityType` enum

**Database Migration:**
- New migration: `V0XX__add_population_breakdown_table.sql`

```sql
CREATE TABLE population_breakdowns (
  id UUID PRIMARY KEY,
  biosample_id UUID NOT NULL REFERENCES biosamples(id),
  at_uri VARCHAR(500),
  analysis_method VARCHAR(50) NOT NULL,
  panel_type VARCHAR(20) NOT NULL,
  reference_populations VARCHAR(50),
  snps_analyzed INTEGER NOT NULL,
  snps_with_genotype INTEGER NOT NULL,
  snps_missing INTEGER NOT NULL,
  confidence_level DOUBLE NOT NULL,
  components CLOB NOT NULL,          -- JSON
  super_population_summary CLOB,     -- JSON
  pca_coordinates CLOB,              -- JSON
  analysis_date TIMESTAMP,
  pipeline_version VARCHAR(50),
  reference_version VARCHAR(50),
  -- EntityMeta columns
  version INTEGER NOT NULL DEFAULT 1,
  created_at TIMESTAMP NOT NULL,
  updated_at TIMESTAMP,
  last_modified_field VARCHAR(100),
  sync_status VARCHAR(20) NOT NULL DEFAULT 'PENDING'
);

CREATE INDEX idx_pop_breakdown_biosample ON population_breakdowns(biosample_id);
CREATE INDEX idx_pop_breakdown_at_uri ON population_breakdowns(at_uri);
```

**Tests:**
- Update all test fixtures using `PopulationPercentage` -> `PopulationComponent`
- Update confidence interval structure in test assertions
- Add repository integration test for `PopulationBreakdownRepository`

---

## 3. HaplogroupReconciliation Enrichment

**Goal:** Add missing Atmosphere Lexicon types to the reconciliation model:
`HeteroplasmyObservation`, `IdentityVerification`, `ManualOverride`, `AuditEntry`.

### 3.1 New Types

```scala
case class HeteroplasmyObservation(
  position: Int,
  majorAllele: String,
  minorAllele: String,
  majorAlleleFrequency: Double,      // 0.5-1.0
  depth: Option[Int],
  isDefiningSnp: Option[Boolean],
  affectedHaplogroup: Option[String]
)

case class IdentityVerification(
  kinshipCoefficient: Option[Double],
  fingerprintSnpConcordance: Option[Double],
  yStrDistance: Option[Int],
  verificationStatus: Option[String], // VERIFIED_SAME, LIKELY_SAME, UNCERTAIN, etc.
  verificationMethod: Option[String]  // AUTOSOMAL_KINSHIP, Y_STR, FINGERPRINT_SNPS, COMBINED
)

case class ManualOverride(
  overriddenHaplogroup: String,
  reason: Option[String],
  overriddenAt: Option[LocalDateTime],
  overriddenBy: Option[String]        // DID
)

case class AuditEntry(
  timestamp: LocalDateTime,
  action: String,                     // INITIAL_RECONCILIATION, RUN_ADDED, etc.
  previousConsensus: Option[String],
  newConsensus: Option[String],
  runRef: Option[String],
  notes: Option[String]
)
```

### 3.2 Field Additions to `HaplogroupReconciliation`

| Field | Type | Default |
|:---|:---|:---|
| `heteroplasmyObservations` | `List[HeteroplasmyObservation]` | `List.empty` |
| `identityVerification` | `Option[IdentityVerification]` | `None` |
| `manualOverride` | `Option[ManualOverride]` | `None` |
| `auditLog` | `List[AuditEntry]` | `List.empty` |

### 3.3 Files to Modify

**Model:**
- `src/main/scala/com/decodingus/workspace/model/HaplogroupReconciliation.scala`
  - Add the 4 new case classes above
  - Add 4 new optional fields to `HaplogroupReconciliation`
  - Default all new fields (backward compatible - existing records deserialize fine)

**Repository:**
- `src/main/scala/com/decodingus/repository/HaplogroupReconciliationRepository.scala`
  - Add Circe codecs for new types
  - Existing JSON column storage handles new fields automatically (additive JSON change)
  - No entity structure change needed if reconciliation data is stored as single JSONB blob
  - If columns are separate: add `heteroplasmy_observations`, `identity_verification`,
    `manual_override`, `audit_log` columns

**Codecs:**
- `src/main/scala/com/decodingus/pds/PdsClient.scala`
  - Add encoder/decoder for `HeteroplasmyObservation`, `IdentityVerification`,
    `ManualOverride`, `AuditEntry`

**Service:**
- `src/main/scala/com/decodingus/workspace/services/WorkspaceOperations.scala`
  - Update `addHaplogroupCall()` to append an `AuditEntry` with action `RUN_ADDED`
  - Update `removeHaplogroupCall()` to append `AuditEntry` with action `RUN_REMOVED`

**UI (optional, low priority):**
- `src/main/scala/com/decodingus/ui/components/ReconciliationDetailDialog.scala`
  - Display heteroplasmy observations if present
  - Display identity verification metrics if present
  - Display audit log entries

**Database Migration (if columns are separate):**
- New migration: `V0XX__enrich_haplogroup_reconciliation.sql`
  - Add JSONB columns for new fields (or no migration needed if stored as single JSON blob)

**Tests:**
- Update `HaplogroupReconciliationRepositorySpec` with new fields
- Update `Phase2EntityIntegrationSpec` fixtures

---

## 4. Sync-Time Validation & Test Type Mapping

**Goal:** Validate required fields before PDS submission and map Navigator's detailed
test type taxonomy to the web app's generic categories.

### 4.1 PDS Submission Validation

Before syncing a record to PDS, validate that fields marked as **required** in the
Atmosphere Lexicon are populated:

**Biosample required fields:**
- `citizenDid` - must be non-empty (currently `Option[String]`)
- `centerName` - must be non-empty (currently `Option[String]`)
- `atUri` - must be assigned (currently `Option[String]`)

**SequenceRun required fields:**
- `atUri`, `biosampleRef`, `platformName`, `testType`, `files` (non-empty)

**Alignment required fields:**
- `atUri`, `sequenceRunRef`, `referenceBuild`, `aligner`

**Genotype required fields:**
- `atUri`, `biosampleRef`, `testTypeCode`, `provider`

**PopulationBreakdown required fields:**
- `atUri`, `biosampleRef`, `analysisMethod`, `panelType`, `components` (non-empty)

### 4.2 Test Type Mapping

Navigator uses a detailed taxonomy; the web app uses broader categories for `sequencerun.testType`:

| Navigator Code | Web App testType | Notes |
|:---|:---|:---|
| `WGS` | `WGS` | Direct match |
| `WGS_LOW_PASS` | `WGS` | Subtype of WGS |
| `WGS_HIFI` | `WGS` | PacBio HiFi is WGS |
| `WGS_NANOPORE` | `WGS` | Nanopore is WGS |
| `WGS_CLR` | `WGS` | PacBio CLR is WGS |
| `EXOME` / `WES` | `EXOME` | Direct match |
| `BIG_Y_500` | `TARGETED` | Y-chromosome targeted |
| `BIG_Y_700` | `TARGETED` | Y-chromosome targeted |
| `Y_ELITE` | `TARGETED` | Y-chromosome targeted |
| `Y_PRIME` | `TARGETED` | Y-chromosome targeted |
| `MT_FULL_SEQUENCE` | `TARGETED` | mtDNA targeted |
| `MT_PLUS` | `TARGETED` | mtDNA targeted |
| `MT_CR_ONLY` | `TARGETED` | mtDNA control region |
| `AMPLICON` | `AMPLICON` | Direct match |
| `RNA_SEQ` | `RNA_SEQ` | Direct match |
| `TARGETED` | `TARGETED` | Direct match |

### 4.3 Files to Modify

- `src/main/scala/com/decodingus/pds/PdsClient.scala` (or new validation module)
  - Add `validateForSync(record: Any): Either[List[String], Unit]` per record type
  - Add `mapTestTypeForSync(navigatorCode: String): String` mapping function

- `src/main/scala/com/decodingus/workspace/model/SequenceRun.scala`
  - Add `toSyncTestType: String` method to companion object
  - Maps detailed Navigator codes to web app categories

### 4.4 Approach

**Do NOT change the internal test type codes.** Navigator's detailed taxonomy is
more useful for local analysis (e.g., knowing `BIG_Y_700` vs `WGS` affects haplogroup
analysis strategy). The mapping is applied only at the PDS sync boundary.

---

## 5. Records Already Aligned (No Changes)

These Navigator models already match or exceed the web app schema. Confirmed no action needed:

| Navigator Model | Lexicon Record | Status |
|:---|:---|:---|
| `RecordMeta` | `defs#recordMeta` | Field-for-field match |
| `FileInfo` | `defs#fileInfo` | Field-for-field match |
| `VariantCall` | `defs#variantCall` | Field-for-field match |
| `PrivateVariantData` | `defs#privateVariantData` | Field-for-field match |
| `HaplogroupAssignments` | `defs#haplogroupAssignments` | Match |
| `HaplogroupResult` | `defs#haplogroupResult` | Navigator superset (adds source, sourceRef, treeProvider, treeVersion, analyzedAt) |
| `StrProfile` | `strProfile` | Match |
| `StrMarkerValue` | `defs#strMarkerValue` | Navigator superset (adds startPosition, endPosition, orderedDate) |
| `StrPanel` | `defs#strPanel` | Match |
| `StrValue` hierarchy | `defs#strValue` union | Match |
| `SnpConflict` | `defs#snpConflict` | Match |
| `SnpCallFromRun` | `defs#snpCallFromRun` | Match |
| `Workspace` | `workspace` | Navigator superset (denormalized records) |
| `Project` | `project` | Match |
| `Alignment` | `alignment` | Match (Navigator already has biosampleRef, variantCaller) |
| `AlignmentMetrics` | `defs#alignmentMetrics` | Navigator superset (adds callableBases, VCF status, sex inference, SV calling) |
| `ContigMetrics` | `defs#contigMetrics` | Navigator superset (adds excessiveCoverage, refN) |
| `SequenceRun` | `sequencerun` | Navigator superset (adds GATK-derived metrics, BAM header fields, fingerprint) |
| `ReconciliationStatus` | `defs#reconciliationStatus` | Match (Navigator already has divergencePoint, branchCompatibilityScore, snpConcordance) |
| `RunHaplogroupCall` | `defs#runHaplogroupCall` | Navigator superset (adds lineagePath, treeProvider) |

## 6. Navigator Extensions (Not in Web Schema)

These are Navigator-specific features that extend beyond the Atmosphere Lexicon.
They should be preserved locally and excluded from PDS sync (or proposed to web schema later):

| Extension | Location | Notes |
|:---|:---|:---|
| `Biosample.ySnpPanelRefs` | Biosample.scala | Y-DNA SNP pack result references |
| `YDnaSnpPanelResult` | YDnaSnpResult.scala | Full SNP pack data model |
| `Biosample.strProfileRefs` (plural) | Biosample.scala | Multiple STR profiles; web has singular `strProfileRef` |
| `SequenceRun` GATK fields | SequenceRun.scala | pfReads, pfReadsAligned, readsPaired, etc. |
| `SequenceRun` BAM header fields | SequenceRun.scala | sampleName, libraryId, platformUnit, runFingerprint |
| `SequenceRun.sequencingFacility` | SequenceRun.scala | Inferred from instrumentId |
| `AlignmentMetrics` VCF/SV/sex fields | AlignmentMetrics.scala | vcfPath, svVcfPath, inferredSex, etc. |
| `ContigMetrics.excessiveCoverage` | ContigMetrics.scala | Not in web schema |
| `HaplogroupResult.source/sourceRef` | HaplogroupResult.scala | Provenance tracking |
| Detailed test type taxonomy | SequenceRun, TestType | BIG_Y_500, WGS_HIFI, MT_PLUS, etc. |

## 7. Future Records (No Action Now)

These exist in the web app schema but are not needed in Navigator until their
features are implemented:

| Record | NSID | When Needed |
|:---|:---|:---|
| `imputation` | `com.decodingus.atmosphere.imputation` | When imputation support added |
| `matchConsent` | `com.decodingus.atmosphere.matchConsent` | When IBD matching added |
| `matchList` | `com.decodingus.atmosphere.matchList` | When IBD matching added |
| `matchRequest` | `com.decodingus.atmosphere.matchRequest` | When IBD matching added |
| `instrumentObservation` | `com.decodingus.atmosphere.instrumentObservation` | When lab inference added |
| `haplogroupAncestralStr` | `com.decodingus.atmosphere.haplogroupAncestralStr` | AppView-computed only |

---

## Implementation Order

### Phase A: ChipProfile alignment (Section 1)
1. DB migration for column rename + new columns
2. Update `ChipProfile` model
3. Update `ChipProfileRepository` entity + codecs
4. Update `EntityConversions`
5. Update all 6 chip parsers (vendor -> provider)
6. Update `WorkbenchViewModel.parseChipFile()`
7. Update `SubjectDetailView` UI references
8. Update PdsClient codecs (verify)
9. Update tests

### Phase B: PopulationBreakdown promotion (Section 2)
1. Create `PopulationBreakdown` model file
2. Rename `PopulationPercentage` -> `PopulationComponent` + restructure CI
3. Rename `SuperPopulationPercentage` -> `SuperPopulationSummary`
4. Update `AncestryResult.fromProbabilities()` factory
5. DB migration for population_breakdowns table
6. Create `PopulationBreakdownRepository`
7. Update `EntityConversions`
8. Update `H2WorkspaceService` + `WorkspaceOperations`
9. Add to `WorkspaceContent`
10. Update `AncestryReportWriter`
11. Update `AncestryResultDialog`
12. Update `WorkbenchViewModel` ancestry handling
13. Add PdsClient codecs
14. Update tests

### Phase C: Reconciliation enrichment (Section 3)
1. Add new case classes to `HaplogroupReconciliation.scala`
2. Add optional fields to `HaplogroupReconciliation`
3. Update repository codecs
4. Update PdsClient codecs
5. Update `WorkspaceOperations` to append audit entries
6. Update tests

### Phase D: Sync validation (Section 4)
1. Add test type mapping to `SequenceRun` companion
2. Add validation module/methods
3. Wire into PDS sync flow

---

## Risk Notes

- **Backward compatibility**: All new fields default to `None` / `List.empty` so existing
  H2 databases deserialize without issues. DB migrations add nullable columns.
- **vendor -> provider rename**: This is a breaking change for the database column.
  Migration must run before app starts. PDS JSON format changes from `"vendor"` to
  `"provider"`. If existing PDS records use `"vendor"`, add a Circe decoder fallback.
- **AncestryResult restructure**: Internal processors can continue returning `AncestryResult`.
  The `PopulationBreakdown` wrapper is applied at the workspace persistence boundary.
  This limits blast radius to persistence/UI layers.
- **Test type mapping**: Only applied at sync boundary. Internal analysis logic continues
  using Navigator's detailed codes.
