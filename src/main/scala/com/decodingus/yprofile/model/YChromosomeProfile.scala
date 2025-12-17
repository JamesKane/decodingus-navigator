package com.decodingus.yprofile.model

import com.decodingus.repository.{Entity, EntityMeta}

import java.time.LocalDateTime
import java.util.UUID

// ============================================
// Enums
// ============================================

/**
 * Source type for Y chromosome profile data.
 * Different testing technologies have different quality characteristics.
 *
 * Each source type carries its quality weights for SNP and STR concordance:
 * - snpWeight: Reliability for SNP/INDEL/MNP detection (0.0-1.0)
 * - strWeight: Reliability for STR repeat counting (0.0-1.0)
 *
 * @param snpWeight Weight for SNP/INDEL/MNP concordance (higher = more reliable)
 * @param strWeight Weight for STR concordance (higher = more reliable)
 */
enum YProfileSourceType(val snpWeight: Double, val strWeight: Double):
  case SANGER extends YProfileSourceType(1.0, 0.9) // Gold standard for SNPs, good for STRs
  case CAPILLARY_ELECTROPHORESIS extends YProfileSourceType(0.5, 1.0) // Not for SNPs, gold standard for STRs
  case WGS_LONG_READ extends YProfileSourceType(0.95, 0.7) // Excellent for SNPs, good for repeats
  case WGS_SHORT_READ extends YProfileSourceType(0.85, 0.5) // Good for SNPs, repeat estimation error-prone
  case TARGETED_NGS extends YProfileSourceType(0.75, 0.4) // Good but limited regions
  case CHIP extends YProfileSourceType(0.5, 0.3) // Probe-based, limited
  case MANUAL extends YProfileSourceType(0.3, 0.2) // User-provided, lowest confidence

  /** Get method tier (0-5 integer) for SNP concordance. */
  def snpTier: Int = math.round(snpWeight * 5).toInt

  /** Get method tier (0-5 integer) for STR concordance. */
  def strTier: Int = math.round(strWeight * 5).toInt

object YProfileSourceType:
  def fromString(s: String): YProfileSourceType = s match
    case "SANGER" => SANGER
    case "CAPILLARY_ELECTROPHORESIS" => CAPILLARY_ELECTROPHORESIS
    case "WGS_LONG_READ" => WGS_LONG_READ
    case "WGS_SHORT_READ" => WGS_SHORT_READ
    case "TARGETED_NGS" => TARGETED_NGS
    case "CHIP" => CHIP
    case "MANUAL" => MANUAL
    case other => throw new IllegalArgumentException(s"Unknown YProfileSourceType: $other")

/**
 * Variant type classification.
 */
enum YVariantType:
  case SNP // Single nucleotide polymorphism
  case INDEL // Insertion/deletion
  case MNP // Multi-nucleotide polymorphism
  case STR // Short tandem repeat

object YVariantType:
  def fromString(s: String): YVariantType = s match
    case "SNP" => SNP
    case "INDEL" => INDEL
    case "MNP" => MNP
    case "STR" => STR
    case other => throw new IllegalArgumentException(s"Unknown YVariantType: $other")

/**
 * Consensus state for a variant call.
 */
enum YConsensusState:
  case DERIVED // Mutant/derived state (positive for haplogroup)
  case ANCESTRAL // Ancestral/reference state
  case HETEROPLASMY // Mixed signal (rare on Y chromosome)
  case NO_CALL // Unable to make a call

object YConsensusState:
  def fromString(s: String): YConsensusState = s match
    case "DERIVED" => DERIVED
    case "ANCESTRAL" => ANCESTRAL
    case "HETEROPLASMY" => HETEROPLASMY
    case "NO_CALL" => NO_CALL
    case other => throw new IllegalArgumentException(s"Unknown YConsensusState: $other")

/**
 * Status of a variant in the unified profile.
 */
enum YVariantStatus:
  case CONFIRMED // Concordant across sources, in reference tree
  case NOVEL // High confidence but not in reference tree (private)
  case CONFLICT // Discordant calls across sources
  case NO_COVERAGE // No source has data at this position
  case PENDING // Awaiting reconciliation

object YVariantStatus:
  def fromString(s: String): YVariantStatus = s match
    case "CONFIRMED" => CONFIRMED
    case "NOVEL" => NOVEL
    case "CONFLICT" => CONFLICT
    case "NO_COVERAGE" => NO_COVERAGE
    case "PENDING" => PENDING
    case other => throw new IllegalArgumentException(s"Unknown YVariantStatus: $other")

/**
 * Callable state for a genomic region.
 *
 * Each state carries a weight factor for concordance calculations:
 * - CALLABLE regions get full weight (1.0)
 * - Problem regions get reduced weight based on reliability
 *
 * @param weight Factor applied to concordance weight (0.0-1.0)
 */
enum YCallableState(val weight: Double):
  case CALLABLE extends YCallableState(1.0) // Full confidence
  case NO_COVERAGE extends YCallableState(0.0) // Zero confidence
  case LOW_COVERAGE extends YCallableState(0.5) // Reduced confidence
  case EXCESSIVE_COVERAGE extends YCallableState(0.3) // Likely artifact
  case POOR_MAPPING_QUALITY extends YCallableState(0.3) // Unreliable mapping
  case REF_N extends YCallableState(0.0) // Can't call
  case SUMMARY extends YCallableState(0.5) // Aggregated

object YCallableState:
  def fromString(s: String): YCallableState = s match
    case "CALLABLE" => CALLABLE
    case "NO_COVERAGE" => NO_COVERAGE
    case "LOW_COVERAGE" => LOW_COVERAGE
    case "EXCESSIVE_COVERAGE" => EXCESSIVE_COVERAGE
    case "POOR_MAPPING_QUALITY" => POOR_MAPPING_QUALITY
    case "REF_N" => REF_N
    case "SUMMARY" => SUMMARY
    case other => throw new IllegalArgumentException(s"Unknown YCallableState: $other")

/**
 * Audit action for manual overrides.
 */
enum YAuditAction:
  case OVERRIDE // Manual override of consensus
  case CONFIRM // Manual confirmation of call
  case REJECT // Manual rejection of call
  case ANNOTATE // Add annotation without changing call
  case REVERT // Revert a previous override

object YAuditAction:
  def fromString(s: String): YAuditAction = s match
    case "OVERRIDE" => OVERRIDE
    case "CONFIRM" => CONFIRM
    case "REJECT" => REJECT
    case "ANNOTATE" => ANNOTATE
    case "REVERT" => REVERT
    case other => throw new IllegalArgumentException(s"Unknown YAuditAction: $other")

/**
 * Naming status for variants.
 */
enum YNamingStatus:
  case UNNAMED // Novel variant, not yet named
  case PENDING_REVIEW // Submitted for naming review
  case NAMED // Has canonical name in the tree

object YNamingStatus:
  def fromString(s: String): YNamingStatus = s match
    case "UNNAMED" => UNNAMED
    case "PENDING_REVIEW" => PENDING_REVIEW
    case "NAMED" => NAMED
    case other => throw new IllegalArgumentException(s"Unknown YNamingStatus: $other")

// ============================================
// STR Metadata
// ============================================

/**
 * Metadata for STR variants, stored as JSON in the database.
 *
 * @param repeatMotif The repeat unit sequence (e.g., "GATA")
 * @param repeatUnit  Length of repeat unit in bp
 * @param copies      For multi-copy markers like DYS385a/b
 * @param rawNotation Original notation if complex (e.g., "22t-25c-26.1t")
 */
case class StrMetadata(
                        repeatMotif: Option[String] = None,
                        repeatUnit: Option[Int] = None,
                        copies: Option[List[Int]] = None,
                        rawNotation: Option[String] = None
                      )

// ============================================
// Entity Classes
// ============================================

/**
 * A contributing test source for the Y chromosome profile.
 *
 * @param id                 Unique identifier
 * @param yProfileId         Parent profile ID
 * @param sourceType         Testing technology used
 * @param sourceRef          Reference identifier (e.g., file path, test ID)
 * @param vendor             Testing provider (FTDNA, YSEQ, etc.)
 * @param testName           Name of the test (Big Y-700, Y-Elite, etc.)
 * @param testDate           When the test was performed
 * @param methodTier         Quality tier (computed from sourceType + variantType)
 * @param meanReadDepth      Average read depth (for sequencing sources)
 * @param meanMappingQuality Average mapping quality (for sequencing sources)
 * @param coveragePct        Percentage of Y chromosome covered
 * @param variantCount       Number of SNP/INDEL/MNP variants contributed
 * @param strMarkerCount     Number of STR markers contributed
 * @param novelVariantCount  Number of novel/private variants
 * @param alignmentId        Optional link to alignment entity
 * @param referenceBuild     Reference genome build (GRCh38, GRCh37)
 * @param importedAt         When this source was imported
 */
case class YProfileSourceEntity(
                                 id: UUID,
                                 yProfileId: UUID,
                                 sourceType: YProfileSourceType,
                                 sourceRef: Option[String],
                                 vendor: Option[String],
                                 testName: Option[String],
                                 testDate: Option[LocalDateTime],
                                 methodTier: Int,
                                 meanReadDepth: Option[Double],
                                 meanMappingQuality: Option[Double],
                                 coveragePct: Option[Double],
                                 variantCount: Int,
                                 strMarkerCount: Int,
                                 novelVariantCount: Int,
                                 alignmentId: Option[UUID],
                                 referenceBuild: Option[String],
                                 importedAt: LocalDateTime
                               ) extends Entity[UUID]

object YProfileSourceEntity:
  def create(
              yProfileId: UUID,
              sourceType: YProfileSourceType,
              sourceRef: Option[String] = None,
              vendor: Option[String] = None,
              testName: Option[String] = None,
              testDate: Option[LocalDateTime] = None,
              methodTier: Int = 1,
              meanReadDepth: Option[Double] = None,
              meanMappingQuality: Option[Double] = None,
              coveragePct: Option[Double] = None,
              variantCount: Int = 0,
              strMarkerCount: Int = 0,
              novelVariantCount: Int = 0,
              alignmentId: Option[UUID] = None,
              referenceBuild: Option[String] = None
            ): YProfileSourceEntity = YProfileSourceEntity(
    id = UUID.randomUUID(),
    yProfileId = yProfileId,
    sourceType = sourceType,
    sourceRef = sourceRef,
    vendor = vendor,
    testName = testName,
    testDate = testDate,
    methodTier = methodTier,
    meanReadDepth = meanReadDepth,
    meanMappingQuality = meanMappingQuality,
    coveragePct = coveragePct,
    variantCount = variantCount,
    strMarkerCount = strMarkerCount,
    novelVariantCount = novelVariantCount,
    alignmentId = alignmentId,
    referenceBuild = referenceBuild,
    importedAt = LocalDateTime.now()
  )

/**
 * A callable region from a source in the Y chromosome profile.
 *
 * @param id                  Unique identifier
 * @param yProfileId          Parent profile ID
 * @param sourceId            Contributing source ID
 * @param contig              Contig name (chrY)
 * @param startPosition       Region start (GRCh38)
 * @param endPosition         Region end (GRCh38)
 * @param callableState       Coverage/quality state
 * @param meanCoverage        Average coverage in region
 * @param meanMappingQuality  Average mapping quality in region
 * @param callableLociCacheId Link to cached callable_loci analysis
 */
case class YProfileRegionEntity(
                                 id: UUID,
                                 yProfileId: UUID,
                                 sourceId: UUID,
                                 contig: String,
                                 startPosition: Long,
                                 endPosition: Long,
                                 callableState: YCallableState,
                                 meanCoverage: Option[Double],
                                 meanMappingQuality: Option[Double],
                                 callableLociCacheId: Option[UUID]
                               ) extends Entity[UUID]

object YProfileRegionEntity:
  def create(
              yProfileId: UUID,
              sourceId: UUID,
              startPosition: Long,
              endPosition: Long,
              callableState: YCallableState,
              contig: String = "chrY",
              meanCoverage: Option[Double] = None,
              meanMappingQuality: Option[Double] = None,
              callableLociCacheId: Option[UUID] = None
            ): YProfileRegionEntity = YProfileRegionEntity(
    id = UUID.randomUUID(),
    yProfileId = yProfileId,
    sourceId = sourceId,
    contig = contig,
    startPosition = startPosition,
    endPosition = endPosition,
    callableState = callableState,
    meanCoverage = meanCoverage,
    meanMappingQuality = meanMappingQuality,
    callableLociCacheId = callableLociCacheId
  )

/**
 * A per-source call for a variant.
 *
 * @param id                     Unique identifier
 * @param variantId              Parent variant ID
 * @param sourceId               Contributing source ID
 * @param calledAllele           The allele called by this source
 * @param callState              DERIVED, ANCESTRAL, NO_CALL, or HETEROPLASMY
 * @param calledRepeatCount      For STRs: the repeat count called
 * @param readDepth              Read depth at this position
 * @param qualityScore           Quality/confidence score
 * @param mappingQuality         Mapping quality
 * @param variantAlleleFrequency VAF for heteroplasmy detection
 * @param callableState          Callable state at this position
 * @param concordanceWeight      Calculated weight for voting
 * @param calledAt               When this call was made
 */
case class YVariantSourceCallEntity(
                                     id: UUID,
                                     variantId: UUID,
                                     sourceId: UUID,
                                     calledAllele: String,
                                     callState: YConsensusState,
                                     calledRepeatCount: Option[Int],
                                     readDepth: Option[Int],
                                     qualityScore: Option[Double],
                                     mappingQuality: Option[Double],
                                     variantAlleleFrequency: Option[Double],
                                     callableState: Option[YCallableState],
                                     concordanceWeight: Double,
                                     calledAt: LocalDateTime
                                   ) extends Entity[UUID]

object YVariantSourceCallEntity:
  def create(
              variantId: UUID,
              sourceId: UUID,
              calledAllele: String,
              callState: YConsensusState,
              calledRepeatCount: Option[Int] = None,
              readDepth: Option[Int] = None,
              qualityScore: Option[Double] = None,
              mappingQuality: Option[Double] = None,
              variantAlleleFrequency: Option[Double] = None,
              callableState: Option[YCallableState] = None,
              concordanceWeight: Double = 1.0
            ): YVariantSourceCallEntity = YVariantSourceCallEntity(
    id = UUID.randomUUID(),
    variantId = variantId,
    sourceId = sourceId,
    calledAllele = calledAllele,
    callState = callState,
    calledRepeatCount = calledRepeatCount,
    readDepth = readDepth,
    qualityScore = qualityScore,
    mappingQuality = mappingQuality,
    variantAlleleFrequency = variantAlleleFrequency,
    callableState = callableState,
    concordanceWeight = concordanceWeight,
    calledAt = LocalDateTime.now()
  )

/**
 * Coordinate representation of a source call in a specific reference build.
 *
 * Key insight: The same source call (one piece of evidence) can be represented
 * in multiple reference builds. These are NOT separate pieces of evidence for
 * concordance - they're just different coordinate representations of the same call.
 *
 * @param id                     Unique identifier
 * @param sourceCallId           Parent source call
 * @param referenceBuild         Reference build (GRCh38, GRCh37, hs1, T2T-CHM13)
 * @param contig                 Contig in this reference
 * @param position               Position in this reference
 * @param refAllele              Reference allele in this reference (may differ due to strand)
 * @param altAllele              Alternate allele in this reference
 * @param calledAllele           The allele called at this position
 * @param readDepth              Read depth for this alignment
 * @param mappingQuality         Mapping quality for this alignment
 * @param baseQuality            Base quality at this position
 * @param variantAlleleFrequency VAF for this alignment
 * @param graphNode              For pangenome: node ID
 * @param graphOffset            For pangenome: offset within node
 * @param alignmentId            Optional link to alignment entity
 * @param createdAt              When this record was created
 */
case class YSourceCallAlignmentEntity(
                                       id: UUID,
                                       sourceCallId: UUID,
                                       referenceBuild: String,
                                       contig: String,
                                       position: Long,
                                       refAllele: String,
                                       altAllele: String,
                                       calledAllele: String,
                                       readDepth: Option[Int],
                                       mappingQuality: Option[Double],
                                       baseQuality: Option[Double],
                                       variantAlleleFrequency: Option[Double],
                                       graphNode: Option[String],
                                       graphOffset: Option[Int],
                                       alignmentId: Option[UUID],
                                       createdAt: LocalDateTime
                                     ) extends Entity[UUID]

object YSourceCallAlignmentEntity:
  def create(
              sourceCallId: UUID,
              referenceBuild: String,
              position: Long,
              refAllele: String,
              altAllele: String,
              calledAllele: String,
              contig: String = "chrY",
              readDepth: Option[Int] = None,
              mappingQuality: Option[Double] = None,
              baseQuality: Option[Double] = None,
              variantAlleleFrequency: Option[Double] = None,
              graphNode: Option[String] = None,
              graphOffset: Option[Int] = None,
              alignmentId: Option[UUID] = None
            ): YSourceCallAlignmentEntity = YSourceCallAlignmentEntity(
    id = UUID.randomUUID(),
    sourceCallId = sourceCallId,
    referenceBuild = referenceBuild,
    contig = contig,
    position = position,
    refAllele = refAllele,
    altAllele = altAllele,
    calledAllele = calledAllele,
    readDepth = readDepth,
    mappingQuality = mappingQuality,
    baseQuality = baseQuality,
    variantAlleleFrequency = variantAlleleFrequency,
    graphNode = graphNode,
    graphOffset = graphOffset,
    alignmentId = alignmentId,
    createdAt = LocalDateTime.now()
  )

/**
 * Coordinates for a novel variant in a specific reference build.
 * Stored as JSON in novel_coordinates field.
 */
case class NovelCoordinates(
                             position: Long,
                             ref: String,
                             alt: String,
                             contig: String = "chrY"
                           )

/**
 * A variant in the Y chromosome profile with concordance information.
 *
 * Key concept: Variant identity is (canonical_name, defining_haplogroup), NOT position.
 * Position varies by reference build; the same variant has different coordinates in
 * GRCh37, GRCh38, hs1, etc. Coordinates are stored at the source_call_alignment level.
 *
 * @param id                    Unique identifier
 * @param yProfileId            Parent profile ID
 * @param canonicalName         Primary variant name (M269, L21) - NULL for unnamed
 * @param namingStatus          UNNAMED, PENDING_REVIEW, or NAMED
 * @param novelCoordinates      For unnamed variants: GRCh38 coordinates (JSON)
 * @param contig                Contig name (chrY) - deprecated, use alignments
 * @param position              Position (GRCh38) - deprecated, use alignments
 * @param endPosition           End position for INDELs
 * @param refAllele             Reference allele - deprecated, use alignments
 * @param altAllele             Alternate allele - deprecated, use alignments
 * @param variantType           SNP, INDEL, MNP, or STR
 * @param variantName           Marker name (deprecated, use canonicalName)
 * @param rsId                  dbSNP rsID if available
 * @param markerName            STR marker name (DYS393, etc.)
 * @param repeatCount           Consensus repeat count for STRs
 * @param strMetadata           Additional STR metadata (JSON)
 * @param consensusAllele       Consensus allele from voting
 * @param consensusState        DERIVED, ANCESTRAL, etc.
 * @param status                CONFIRMED, NOVEL, CONFLICT, etc.
 * @param sourceCount           Number of sources with data (NOT alignments!)
 * @param concordantCount       Sources agreeing with consensus
 * @param discordantCount       Sources disagreeing with consensus
 * @param confidenceScore       Weighted confidence score
 * @param maxReadDepth          Best read depth across sources
 * @param maxQualityScore       Best quality score across sources
 * @param definingHaplogroup    Haplogroup this variant defines
 * @param haplogroupBranchDepth Tree depth of defining haplogroup
 * @param lastUpdatedAt         Last update timestamp
 */
case class YProfileVariantEntity(
                                  id: UUID,
                                  yProfileId: UUID,
                                  // Variant identity (new)
                                  canonicalName: Option[String],
                                  namingStatus: YNamingStatus,
                                  novelCoordinates: Option[Map[String, NovelCoordinates]], // {"GRCh38": {...}}
                                  // Legacy coordinates (kept for backward compatibility, deprecated)
                                  contig: String,
                                  position: Long,
                                  endPosition: Option[Long],
                                  refAllele: String,
                                  altAllele: String,
                                  variantType: YVariantType,
                                  variantName: Option[String],
                                  rsId: Option[String],
                                  markerName: Option[String],
                                  repeatCount: Option[Int],
                                  strMetadata: Option[StrMetadata],
                                  consensusAllele: Option[String],
                                  consensusState: YConsensusState,
                                  status: YVariantStatus,
                                  sourceCount: Int,
                                  concordantCount: Int,
                                  discordantCount: Int,
                                  confidenceScore: Double,
                                  maxReadDepth: Option[Int],
                                  maxQualityScore: Option[Double],
                                  definingHaplogroup: Option[String],
                                  haplogroupBranchDepth: Option[Int],
                                  lastUpdatedAt: LocalDateTime
                                ) extends Entity[UUID]

object YProfileVariantEntity:
  def create(
              yProfileId: UUID,
              position: Long,
              refAllele: String,
              altAllele: String,
              variantType: YVariantType = YVariantType.SNP,
              contig: String = "chrY",
              endPosition: Option[Long] = None,
              canonicalName: Option[String] = None,
              namingStatus: YNamingStatus = YNamingStatus.UNNAMED,
              novelCoordinates: Option[Map[String, NovelCoordinates]] = None,
              variantName: Option[String] = None,
              rsId: Option[String] = None,
              markerName: Option[String] = None,
              repeatCount: Option[Int] = None,
              strMetadata: Option[StrMetadata] = None,
              consensusAllele: Option[String] = None,
              consensusState: YConsensusState = YConsensusState.NO_CALL,
              status: YVariantStatus = YVariantStatus.PENDING,
              sourceCount: Int = 0,
              concordantCount: Int = 0,
              discordantCount: Int = 0,
              confidenceScore: Double = 0.0,
              maxReadDepth: Option[Int] = None,
              maxQualityScore: Option[Double] = None,
              definingHaplogroup: Option[String] = None,
              haplogroupBranchDepth: Option[Int] = None
            ): YProfileVariantEntity =
    // Auto-populate canonical_name from variant_name if provided
    val effectiveCanonicalName = canonicalName.orElse(variantName)
    val effectiveNamingStatus = if effectiveCanonicalName.isDefined then YNamingStatus.NAMED else namingStatus
    // Auto-create novel_coordinates for unnamed variants
    val effectiveNovelCoordinates = if effectiveCanonicalName.isEmpty && novelCoordinates.isEmpty then
      Some(Map("GRCh38" -> NovelCoordinates(position, refAllele, altAllele, contig)))
    else
      novelCoordinates

    YProfileVariantEntity(
      id = UUID.randomUUID(),
      yProfileId = yProfileId,
      canonicalName = effectiveCanonicalName,
      namingStatus = effectiveNamingStatus,
      novelCoordinates = effectiveNovelCoordinates,
      contig = contig,
      position = position,
      endPosition = endPosition,
      refAllele = refAllele,
      altAllele = altAllele,
      variantType = variantType,
      variantName = variantName,
      rsId = rsId,
      markerName = markerName,
      repeatCount = repeatCount,
      strMetadata = strMetadata,
      consensusAllele = consensusAllele,
      consensusState = consensusState,
      status = status,
      sourceCount = sourceCount,
      concordantCount = concordantCount,
      discordantCount = discordantCount,
      confidenceScore = confidenceScore,
      maxReadDepth = maxReadDepth,
      maxQualityScore = maxQualityScore,
      definingHaplogroup = definingHaplogroup,
      haplogroupBranchDepth = haplogroupBranchDepth,
      lastUpdatedAt = LocalDateTime.now()
    )

/**
 * Audit trail entry for manual variant overrides.
 *
 * @param id                      Unique identifier
 * @param variantId               Variant being audited
 * @param action                  Type of audit action
 * @param previousConsensusAllele Previous consensus allele
 * @param previousConsensusState  Previous consensus state
 * @param previousStatus          Previous status
 * @param previousConfidence      Previous confidence score
 * @param newConsensusAllele      New consensus allele
 * @param newConsensusState       New consensus state
 * @param newStatus               New status
 * @param newConfidence           New confidence score
 * @param userId                  User who made the change
 * @param reason                  Required reason for change
 * @param supportingEvidence      Optional supporting evidence
 * @param createdAt               When the audit entry was created
 */
case class YVariantAuditEntity(
                                id: UUID,
                                variantId: UUID,
                                action: YAuditAction,
                                previousConsensusAllele: Option[String],
                                previousConsensusState: Option[YConsensusState],
                                previousStatus: Option[YVariantStatus],
                                previousConfidence: Option[Double],
                                newConsensusAllele: Option[String],
                                newConsensusState: Option[YConsensusState],
                                newStatus: Option[YVariantStatus],
                                newConfidence: Option[Double],
                                userId: Option[String],
                                reason: String,
                                supportingEvidence: Option[String],
                                createdAt: LocalDateTime
                              ) extends Entity[UUID]

object YVariantAuditEntity:
  def create(
              variantId: UUID,
              action: YAuditAction,
              reason: String,
              previousConsensusAllele: Option[String] = None,
              previousConsensusState: Option[YConsensusState] = None,
              previousStatus: Option[YVariantStatus] = None,
              previousConfidence: Option[Double] = None,
              newConsensusAllele: Option[String] = None,
              newConsensusState: Option[YConsensusState] = None,
              newStatus: Option[YVariantStatus] = None,
              newConfidence: Option[Double] = None,
              userId: Option[String] = None,
              supportingEvidence: Option[String] = None
            ): YVariantAuditEntity = YVariantAuditEntity(
    id = UUID.randomUUID(),
    variantId = variantId,
    action = action,
    previousConsensusAllele = previousConsensusAllele,
    previousConsensusState = previousConsensusState,
    previousStatus = previousStatus,
    previousConfidence = previousConfidence,
    newConsensusAllele = newConsensusAllele,
    newConsensusState = newConsensusState,
    newStatus = newStatus,
    newConfidence = newConfidence,
    userId = userId,
    reason = reason,
    supportingEvidence = supportingEvidence,
    createdAt = LocalDateTime.now()
  )

/**
 * Main Y chromosome profile entity.
 * One per biosample, aggregating all Y-DNA test results.
 *
 * @param id                     Unique identifier
 * @param biosampleId            Parent biosample
 * @param consensusHaplogroup    Determined haplogroup from unified variants
 * @param haplogroupConfidence   Confidence in haplogroup determination
 * @param haplogroupTreeProvider Tree provider used (ftdna, decodingus)
 * @param haplogroupTreeVersion  Tree version used
 * @param totalVariants          Total SNP/INDEL/MNP variants
 * @param confirmedCount         Variants with CONFIRMED status
 * @param novelCount             Variants with NOVEL status (private)
 * @param conflictCount          Variants with CONFLICT status
 * @param noCoverageCount        Variants with NO_COVERAGE status
 * @param strMarkerCount         Total STR markers
 * @param strConfirmedCount      STR markers with CONFIRMED status
 * @param overallConfidence      Weighted overall profile confidence
 * @param callableRegionPct      Percentage of Y with callable coverage
 * @param meanCoverage           Mean coverage across sources
 * @param sourceCount            Number of contributing test sources
 * @param primarySourceType      Primary/best source type
 * @param lastReconciledAt       Last reconciliation timestamp
 * @param meta                   Entity metadata (sync, version)
 */
case class YChromosomeProfileEntity(
                                     id: UUID,
                                     biosampleId: UUID,
                                     consensusHaplogroup: Option[String],
                                     haplogroupConfidence: Option[Double],
                                     haplogroupTreeProvider: Option[String],
                                     haplogroupTreeVersion: Option[String],
                                     totalVariants: Int,
                                     confirmedCount: Int,
                                     novelCount: Int,
                                     conflictCount: Int,
                                     noCoverageCount: Int,
                                     strMarkerCount: Int,
                                     strConfirmedCount: Int,
                                     overallConfidence: Option[Double],
                                     callableRegionPct: Option[Double],
                                     meanCoverage: Option[Double],
                                     sourceCount: Int,
                                     primarySourceType: Option[YProfileSourceType],
                                     lastReconciledAt: Option[LocalDateTime],
                                     meta: EntityMeta
                                   ) extends Entity[UUID]

object YChromosomeProfileEntity:
  def create(
              biosampleId: UUID,
              consensusHaplogroup: Option[String] = None,
              haplogroupConfidence: Option[Double] = None,
              haplogroupTreeProvider: Option[String] = None,
              haplogroupTreeVersion: Option[String] = None,
              totalVariants: Int = 0,
              confirmedCount: Int = 0,
              novelCount: Int = 0,
              conflictCount: Int = 0,
              noCoverageCount: Int = 0,
              strMarkerCount: Int = 0,
              strConfirmedCount: Int = 0,
              overallConfidence: Option[Double] = None,
              callableRegionPct: Option[Double] = None,
              meanCoverage: Option[Double] = None,
              sourceCount: Int = 0,
              primarySourceType: Option[YProfileSourceType] = None
            ): YChromosomeProfileEntity = YChromosomeProfileEntity(
    id = UUID.randomUUID(),
    biosampleId = biosampleId,
    consensusHaplogroup = consensusHaplogroup,
    haplogroupConfidence = haplogroupConfidence,
    haplogroupTreeProvider = haplogroupTreeProvider,
    haplogroupTreeVersion = haplogroupTreeVersion,
    totalVariants = totalVariants,
    confirmedCount = confirmedCount,
    novelCount = novelCount,
    conflictCount = conflictCount,
    noCoverageCount = noCoverageCount,
    strMarkerCount = strMarkerCount,
    strConfirmedCount = strConfirmedCount,
    overallConfidence = overallConfidence,
    callableRegionPct = callableRegionPct,
    meanCoverage = meanCoverage,
    sourceCount = sourceCount,
    primarySourceType = primarySourceType,
    lastReconciledAt = None,
    meta = EntityMeta.create()
  )
