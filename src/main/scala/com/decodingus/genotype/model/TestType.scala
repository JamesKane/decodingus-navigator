package com.decodingus.genotype.model

import com.typesafe.config.{Config, ConfigFactory}
import io.circe.{Codec, Decoder, Encoder}
import org.slf4j.LoggerFactory

import scala.jdk.CollectionConverters.*
import scala.util.Try

/**
 * Data generation method - sequencing vs genotyping arrays.
 */
enum DataGenerationMethod derives Codec.AsObject:
  case Sequencing
  case Genotyping

/**
 * Target type for the test.
 */
enum TargetType derives Codec.AsObject:
  case WholeGenome
  case YChromosome
  case MtDna
  case Autosomal
  case XChromosome
  case Mixed

/**
 * Test type definition matching the server-side taxonomy.
 *
 * @see multi-test-type-roadmap.md for full taxonomy
 */
case class TestTypeDefinition(
                               code: String,
                               displayName: String,
                               category: DataGenerationMethod,
                               vendor: Option[String],
                               targetType: TargetType,
                               expectedMinDepth: Option[Double],
                               expectedTargetDepth: Option[Double],
                               expectedMarkerCount: Option[Int],
                               supportsHaplogroupY: Boolean,
                               supportsHaplogroupMt: Boolean,
                               supportsAutosomalIbd: Boolean,
                               supportsAncestry: Boolean,
                               typicalFileFormats: List[String]
                             ) derives Codec.AsObject

/**
 * Known test types - sequencing and genotyping.
 *
 * Test types are loaded from HOCON configuration (test_types.conf) at startup.
 * This allows adding or modifying test types without recompiling.
 */
object TestTypes {

  private val log = LoggerFactory.getLogger(getClass)

  // ============================================================================
  // Configuration Loading
  // ============================================================================

  private def loadTestTypesFromConfig(): List[TestTypeDefinition] = {
    try {
      val config = ConfigFactory.load("test_types.conf")
      if (!config.hasPath("test-types")) {
        log.warn("No test-types configuration found, using defaults")
        return defaultTestTypes
      }

      val testTypesConfig = config.getConfig("test-types")
      val codes = testTypesConfig.root().keySet().asScala.toList

      codes.flatMap { code =>
        Try {
          val typeConfig = testTypesConfig.getConfig(code)
          parseTestType(code, typeConfig)
        }.toOption.orElse {
          log.warn(s"Failed to parse test type: $code")
          None
        }
      }
    } catch {
      case e: Exception =>
        log.error(s"Failed to load test types from config: ${e.getMessage}", e)
        defaultTestTypes
    }
  }

  private def parseTestType(code: String, config: Config): TestTypeDefinition = {
    val category = config.getString("category").toLowerCase match {
      case "sequencing" => DataGenerationMethod.Sequencing
      case "genotyping" => DataGenerationMethod.Genotyping
      case other => throw new IllegalArgumentException(s"Unknown category: $other")
    }

    val targetType = config.getString("target-type").toLowerCase match {
      case "whole-genome" => TargetType.WholeGenome
      case "y-chromosome" => TargetType.YChromosome
      case "mt-dna" => TargetType.MtDna
      case "autosomal" => TargetType.Autosomal
      case "x-chromosome" => TargetType.XChromosome
      case "mixed" => TargetType.Mixed
      case other => throw new IllegalArgumentException(s"Unknown target type: $other")
    }

    val supportsConfig = config.getConfig("supports")

    TestTypeDefinition(
      code = code,
      displayName = config.getString("display-name"),
      category = category,
      vendor = if (config.hasPath("vendor")) Some(config.getString("vendor")) else None,
      targetType = targetType,
      expectedMinDepth = if (config.hasPath("expected-min-depth")) Some(config.getDouble("expected-min-depth")) else None,
      expectedTargetDepth = if (config.hasPath("expected-target-depth")) Some(config.getDouble("expected-target-depth")) else None,
      expectedMarkerCount = if (config.hasPath("expected-marker-count")) Some(config.getInt("expected-marker-count")) else None,
      supportsHaplogroupY = supportsConfig.getBoolean("haplogroup-y"),
      supportsHaplogroupMt = supportsConfig.getBoolean("haplogroup-mt"),
      supportsAutosomalIbd = supportsConfig.getBoolean("autosomal-ibd"),
      supportsAncestry = supportsConfig.getBoolean("ancestry"),
      typicalFileFormats = config.getStringList("typical-file-formats").asScala.toList
    )
  }

  // ============================================================================
  // Default Test Types (Fallback)
  // ============================================================================

  private def defaultTestTypes: List[TestTypeDefinition] = List(
    TestTypeDefinition("WGS", "Whole Genome Sequencing", DataGenerationMethod.Sequencing, None,
      TargetType.WholeGenome, Some(10.0), Some(30.0), None, true, true, true, true, List("BAM", "CRAM", "VCF")),
    TestTypeDefinition("WGS_LOW_PASS", "Low-Pass WGS", DataGenerationMethod.Sequencing, None,
      TargetType.WholeGenome, Some(0.5), Some(4.0), None, true, true, true, true, List("BAM", "CRAM", "VCF")),
    TestTypeDefinition("BIG_Y_700", "FTDNA Big Y-700", DataGenerationMethod.Sequencing, Some("FamilyTreeDNA"),
      TargetType.YChromosome, Some(30.0), Some(50.0), None, true, false, false, false, List("BAM", "VCF", "BED")),
    TestTypeDefinition("Y_ELITE", "Full Genomes Y Elite", DataGenerationMethod.Sequencing, Some("Full Genomes"),
      TargetType.YChromosome, Some(20.0), Some(30.0), None, true, true, false, false, List("BAM", "CRAM", "VCF")),
    TestTypeDefinition("MT_FULL_SEQUENCE", "mtDNA Full Sequence", DataGenerationMethod.Sequencing, None,
      TargetType.MtDna, Some(500.0), Some(1000.0), None, false, true, false, false, List("BAM", "FASTA", "VCF"))
  )

  // ============================================================================
  // Loaded Test Types
  // ============================================================================

  /**
   * All known test types loaded from configuration.
   */
  lazy val all: List[TestTypeDefinition] = {
    val loaded = loadTestTypesFromConfig()
    if (loaded.nonEmpty) {
      log.info(s"Loaded ${loaded.size} test types from configuration")
      loaded
    } else {
      log.warn("No test types loaded, using defaults")
      defaultTestTypes
    }
  }

  // Named accessors for commonly used test types (for backward compatibility)
  lazy val WGS: TestTypeDefinition = byCode("WGS").getOrElse(defaultTestTypes.head)
  lazy val WGS_LOW_PASS: TestTypeDefinition = byCode("WGS_LOW_PASS").getOrElse(WGS)
  lazy val WGS_HIFI: TestTypeDefinition = byCode("WGS_HIFI").getOrElse(WGS)
  lazy val WGS_NANOPORE: TestTypeDefinition = byCode("WGS_NANOPORE").getOrElse(WGS)
  lazy val WGS_CLR: TestTypeDefinition = byCode("WGS_CLR").getOrElse(WGS)
  lazy val WES: TestTypeDefinition = byCode("WES").getOrElse(WGS)
  lazy val BIG_Y_500: TestTypeDefinition = byCode("BIG_Y_500").getOrElse(byCode("BIG_Y_700").getOrElse(WGS))
  lazy val BIG_Y_700: TestTypeDefinition = byCode("BIG_Y_700").getOrElse(WGS)
  lazy val Y_ELITE: TestTypeDefinition = byCode("Y_ELITE").getOrElse(BIG_Y_700)
  lazy val Y_PRIME: TestTypeDefinition = byCode("Y_PRIME").getOrElse(BIG_Y_700)
  lazy val MT_FULL_SEQUENCE: TestTypeDefinition = byCode("MT_FULL_SEQUENCE").getOrElse(WGS)
  lazy val MT_PLUS: TestTypeDefinition = byCode("MT_PLUS").getOrElse(MT_FULL_SEQUENCE)
  lazy val MT_CR_ONLY: TestTypeDefinition = byCode("MT_CR_ONLY").getOrElse(MT_FULL_SEQUENCE)
  lazy val YDNA_SNP_PACK_FTDNA: TestTypeDefinition = byCode("YDNA_SNP_PACK_FTDNA").getOrElse(BIG_Y_700)
  lazy val YDNA_PANEL_YSEQ: TestTypeDefinition = byCode("YDNA_PANEL_YSEQ").getOrElse(BIG_Y_700)
  lazy val YDNA_SNP_PANEL: TestTypeDefinition = byCode("YDNA_SNP_PANEL").getOrElse(BIG_Y_700)
  lazy val ARRAY_BISDNA: TestTypeDefinition = byCode("ARRAY_BISDNA").getOrElse(WGS)
  lazy val ARRAY_23ANDME_V5: TestTypeDefinition = byCode("ARRAY_23ANDME_V5").getOrElse(WGS)
  lazy val ARRAY_23ANDME_V4: TestTypeDefinition = byCode("ARRAY_23ANDME_V4").getOrElse(WGS)
  lazy val ARRAY_ANCESTRY_V2: TestTypeDefinition = byCode("ARRAY_ANCESTRY_V2").getOrElse(WGS)
  lazy val ARRAY_FTDNA_FF: TestTypeDefinition = byCode("ARRAY_FTDNA_FF").getOrElse(WGS)
  lazy val ARRAY_MYHERITAGE: TestTypeDefinition = byCode("ARRAY_MYHERITAGE").getOrElse(WGS)
  lazy val ARRAY_LIVINGDNA: TestTypeDefinition = byCode("ARRAY_LIVINGDNA").getOrElse(WGS)

  // ============================================================================
  // Collections
  // ============================================================================

  /**
   * All sequencing test types.
   */
  def sequencing: List[TestTypeDefinition] =
    all.filter(_.category == DataGenerationMethod.Sequencing)

  /**
   * Short-read WGS types (Illumina, BGI, Element, Ultima).
   */
  def shortReadWgs: List[TestTypeDefinition] =
    all.filter(t => t.code == "WGS" || t.code == "WGS_LOW_PASS")

  /**
   * Long-read WGS types (PacBio HiFi, Nanopore, PacBio CLR).
   * These excel at resolving complex/repetitive regions including Y palindromes.
   */
  def longReadWgs: List[TestTypeDefinition] =
    all.filter(t => t.code == "WGS_HIFI" || t.code == "WGS_NANOPORE" || t.code == "WGS_CLR")

  /**
   * All genotyping array types.
   */
  def genotypingArrays: List[TestTypeDefinition] =
    all.filter(_.category == DataGenerationMethod.Genotyping)

  /**
   * Targeted Y-DNA sequencing types (Big Y, Y Elite, etc.).
   */
  def targetedYDnaSequencing: List[TestTypeDefinition] =
    sequencing.filter(_.targetType == TargetType.YChromosome)

  /**
   * Y-DNA SNP packs and panels (probe-based genotyping).
   * These are NOT sequencing - they test specific named SNPs.
   */
  def yDnaSnpPanels: List[TestTypeDefinition] =
    genotypingArrays.filter(_.targetType == TargetType.YChromosome)

  /**
   * All Y-DNA capable test types (sequencing + SNP panels).
   */
  def allYDnaTests: List[TestTypeDefinition] =
    targetedYDnaSequencing ++ yDnaSnpPanels

  /**
   * Targeted mtDNA sequencing types.
   */
  def targetedMtDna: List[TestTypeDefinition] =
    sequencing.filter(_.targetType == TargetType.MtDna)

  /**
   * Get test type by code.
   */
  def byCode(code: String): Option[TestTypeDefinition] =
    all.find(_.code.equalsIgnoreCase(code))

  /**
   * Get test types that support Y-DNA haplogroup analysis.
   */
  def yDnaCapable: List[TestTypeDefinition] =
    all.filter(_.supportsHaplogroupY)

  /**
   * Get test types that support mtDNA haplogroup analysis.
   */
  def mtDnaCapable: List[TestTypeDefinition] =
    all.filter(_.supportsHaplogroupMt)

  /**
   * Get test types that support ancestry analysis.
   */
  def ancestryCapable: List[TestTypeDefinition] =
    all.filter(_.supportsAncestry)

  /**
   * Infer test type from coverage statistics and platform info.
   *
   * Uses a heuristic based on:
   * - Platform and read length (for long-read distinction) - checked FIRST
   * - Whether only Y or MT has coverage (targeted tests)
   * - Coverage depth patterns
   *
   * Note: Platform takes precedence over coverage for WGS classification.
   * A 4x HiFi is still HiFi (not low-pass) because the long-read technology
   * still provides value in resolving complex regions.
   *
   * @param yCoverage         Mean Y chromosome coverage (None if no reads)
   * @param mtCoverage        Mean mtDNA coverage (None if no reads)
   * @param autosomalCoverage Mean autosomal coverage (None if no reads)
   * @param totalReads        Total read count
   * @param vendor            Optional vendor hint from file/header
   * @param platform          Optional platform (PacBio, Nanopore, Illumina, etc.)
   * @param meanReadLength    Optional mean read length for long-read detection
   * @return Best-guess test type, defaults to WGS
   */
  def inferFromCoverage(
                         yCoverage: Option[Double],
                         mtCoverage: Option[Double],
                         autosomalCoverage: Option[Double],
                         totalReads: Long,
                         vendor: Option[String] = None,
                         platform: Option[String] = None,
                         meanReadLength: Option[Int] = None
                       ): TestTypeDefinition = {

    val hasYCoverage = yCoverage.exists(_ > 1.0)
    val hasMtCoverage = mtCoverage.exists(_ > 10.0)
    val hasAutosomalCoverage = autosomalCoverage.exists(_ > 1.0)

    // Check for long-read platform FIRST - platform determines capabilities
    // regardless of coverage depth (a 4x HiFi is still HiFi, not low-pass)
    val isLongReadPlatform = platform.exists { p =>
      val pl = p.toLowerCase
      pl.contains("pacbio") || pl.contains("nanopore") || pl.contains("ont")
    }
    val isLongReadByLength = meanReadLength.exists(_ > 1000)

    if (isLongReadPlatform || isLongReadByLength) {
      // Long-read WGS (unless it's targeted Y or MT only)
      if (hasYCoverage && !hasAutosomalCoverage) {
        // Targeted Y on long-read - still use Y-specific type
        vendor.map(_.toLowerCase) match {
          case Some(v) if v.contains("ftdna") || v.contains("familytreedna") => BIG_Y_700
          case Some(v) if v.contains("full genomes") => Y_ELITE
          case Some(v) if v.contains("yseq") => Y_PRIME
          case _ => BIG_Y_700
        }
      } else if (hasMtCoverage && !hasAutosomalCoverage && !hasYCoverage) {
        MT_FULL_SEQUENCE
      } else {
        // Full WGS on long-read platform
        inferWgsType(platform, meanReadLength)
      }
    }
    // Targeted Y test: Y coverage but minimal/no autosomal
    else if (hasYCoverage && !hasAutosomalCoverage) {
      vendor.map(_.toLowerCase) match {
        case Some(v) if v.contains("ftdna") || v.contains("familytreedna") => BIG_Y_700
        case Some(v) if v.contains("full genomes") => Y_ELITE
        case Some(v) if v.contains("yseq") => Y_PRIME
        case _ => BIG_Y_700 // Default targeted Y
      }
    }
    // Targeted MT test: MT coverage but minimal/no autosomal or Y
    else if (hasMtCoverage && !hasAutosomalCoverage && !hasYCoverage) {
      vendor.map(_.toLowerCase) match {
        case Some(v) if v.contains("ftdna") || v.contains("familytreedna") => MT_PLUS
        case _ => MT_FULL_SEQUENCE
      }
    }
    // Exome: high autosomal but sparse overall (no Y/MT)
    else if (hasAutosomalCoverage && !hasYCoverage && !hasMtCoverage) {
      autosomalCoverage match {
        case Some(cov) if cov > 50.0 => WES
        case _ => WGS // Could be low-pass short-read
      }
    }
    // Low-pass WGS: very low autosomal coverage (short-read only)
    else if (autosomalCoverage.exists(_ < 5.0) && hasAutosomalCoverage) {
      WGS_LOW_PASS
    }
    // Default WGS (short-read)
    else {
      WGS
    }
  }

  /**
   * Infer WGS subtype based on platform and read length.
   *
   * @param platform       Sequencing platform (PacBio, Nanopore, Illumina, etc.)
   * @param meanReadLength Mean read length in bp
   * @return Appropriate WGS test type
   */
  def inferWgsType(
                    platform: Option[String],
                    meanReadLength: Option[Int]
                  ): TestTypeDefinition = {
    val platformLower = platform.map(_.toLowerCase)
    val isLongRead = meanReadLength.exists(_ > 1000)

    platformLower match {
      case Some(p) if p.contains("pacbio") =>
        // HiFi reads are typically 10-25kb with high accuracy
        // CLR reads can be longer but lower accuracy
        if (meanReadLength.exists(_ > 5000)) WGS_HIFI
        else WGS_CLR

      case Some(p) if p.contains("nanopore") || p.contains("ont") =>
        WGS_NANOPORE

      case Some(p) if isLongRead =>
        // Unknown long-read platform - default to HiFi as most common now
        WGS_HIFI

      case _ =>
        // Short-read or unknown - default to generic WGS (Illumina-style)
        WGS
    }
  }
}
