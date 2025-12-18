package com.decodingus.genotype.model

import io.circe.{Codec, Decoder, Encoder}

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
 */
object TestTypes {

  // ============ SEQUENCING: Whole Genome - Short Read ============

  val WGS: TestTypeDefinition = TestTypeDefinition(
    code = "WGS",
    displayName = "Whole Genome Sequencing",
    category = DataGenerationMethod.Sequencing,
    vendor = None,
    targetType = TargetType.WholeGenome,
    expectedMinDepth = Some(10.0),
    expectedTargetDepth = Some(30.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = true,
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = true,
    supportsAncestry = true,
    typicalFileFormats = List("BAM", "CRAM", "VCF")
  )

  val WGS_LOW_PASS: TestTypeDefinition = TestTypeDefinition(
    code = "WGS_LOW_PASS",
    displayName = "Low-Pass WGS",
    category = DataGenerationMethod.Sequencing,
    vendor = None,
    targetType = TargetType.WholeGenome,
    expectedMinDepth = Some(0.5),
    expectedTargetDepth = Some(4.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = true, // Limited but possible
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = true, // Via imputation
    supportsAncestry = true,
    typicalFileFormats = List("BAM", "CRAM", "VCF")
  )

  // ============ SEQUENCING: Whole Genome - Long Read ============
  // Long-read technologies excel at resolving complex/repetitive regions
  // including palindromic regions on the Y chromosome.

  val WGS_HIFI: TestTypeDefinition = TestTypeDefinition(
    code = "WGS_HIFI",
    displayName = "PacBio HiFi WGS",
    category = DataGenerationMethod.Sequencing,
    vendor = Some("PacBio"),
    targetType = TargetType.WholeGenome,
    expectedMinDepth = Some(15.0),
    expectedTargetDepth = Some(30.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = true, // Excellent - resolves complex Y regions
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = true,
    supportsAncestry = true,
    typicalFileFormats = List("BAM", "CRAM", "VCF")
  )

  val WGS_NANOPORE: TestTypeDefinition = TestTypeDefinition(
    code = "WGS_NANOPORE",
    displayName = "Nanopore WGS",
    category = DataGenerationMethod.Sequencing,
    vendor = Some("Oxford Nanopore"),
    targetType = TargetType.WholeGenome,
    expectedMinDepth = Some(15.0),
    expectedTargetDepth = Some(30.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = true, // Good - ultra-long reads span complex regions
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = true,
    supportsAncestry = true,
    typicalFileFormats = List("BAM", "CRAM", "VCF", "FASTQ")
  )

  // PacBio CLR (Continuous Long Read) - older technology, lower accuracy than HiFi
  val WGS_CLR: TestTypeDefinition = TestTypeDefinition(
    code = "WGS_CLR",
    displayName = "PacBio CLR WGS",
    category = DataGenerationMethod.Sequencing,
    vendor = Some("PacBio"),
    targetType = TargetType.WholeGenome,
    expectedMinDepth = Some(20.0),
    expectedTargetDepth = Some(40.0), // Higher depth needed due to lower accuracy
    expectedMarkerCount = None,
    supportsHaplogroupY = true,
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = true,
    supportsAncestry = true,
    typicalFileFormats = List("BAM", "VCF")
  )

  // ============ SEQUENCING: Exome ============

  val WES: TestTypeDefinition = TestTypeDefinition(
    code = "WES",
    displayName = "Whole Exome Sequencing",
    category = DataGenerationMethod.Sequencing,
    vendor = None,
    targetType = TargetType.Autosomal, // Primarily exonic regions
    expectedMinDepth = Some(50.0),
    expectedTargetDepth = Some(100.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = false, // No Y coverage
    supportsHaplogroupMt = false, // No mtDNA coverage
    supportsAutosomalIbd = false, // Sparse coverage
    supportsAncestry = false,
    typicalFileFormats = List("BAM", "CRAM", "VCF")
  )

  // ============ GENOTYPING: Y-DNA SNP Packs/Panels ============
  // These are probe-based or Sanger sequencing tests for specific Y-SNPs.
  // Results typically delivered as positive/negative calls for named SNPs.

  val YDNA_SNP_PACK_FTDNA: TestTypeDefinition = TestTypeDefinition(
    code = "YDNA_SNP_PACK_FTDNA",
    displayName = "FTDNA SNP Pack",
    category = DataGenerationMethod.Genotyping, // Probe-based, not sequencing
    vendor = Some("FamilyTreeDNA"),
    targetType = TargetType.YChromosome,
    expectedMinDepth = None,
    expectedTargetDepth = None,
    expectedMarkerCount = Some(100), // Varies by pack, typically 50-200 SNPs
    supportsHaplogroupY = true,
    supportsHaplogroupMt = false,
    supportsAutosomalIbd = false,
    supportsAncestry = false,
    typicalFileFormats = List("CSV", "TXT")
  )

  val YDNA_PANEL_YSEQ: TestTypeDefinition = TestTypeDefinition(
    code = "YDNA_PANEL_YSEQ",
    displayName = "YSEQ Panel",
    category = DataGenerationMethod.Genotyping,
    vendor = Some("YSEQ"),
    targetType = TargetType.YChromosome,
    expectedMinDepth = None,
    expectedTargetDepth = None,
    expectedMarkerCount = None, // Continuous delivery - grows over time
    supportsHaplogroupY = true,
    supportsHaplogroupMt = false,
    supportsAutosomalIbd = false,
    supportsAncestry = false,
    typicalFileFormats = List("CSV", "YSEQ") // Custom YSEQ format for incremental updates
  )

  // BISDNA was a genotyping array with Y-DNA emphasis.
  // The underlying chip included autosomal markers, but only Y-DNA results
  // were delivered to customers as raw data.
  val ARRAY_BISDNA: TestTypeDefinition = TestTypeDefinition(
    code = "ARRAY_BISDNA",
    displayName = "BISDNA Array",
    category = DataGenerationMethod.Genotyping,
    vendor = Some("BISDNA"),
    targetType = TargetType.YChromosome, // Only Y-DNA delivered to customers
    expectedMinDepth = None,
    expectedTargetDepth = None,
    expectedMarkerCount = Some(15000), // Y-DNA markers delivered
    supportsHaplogroupY = true,
    supportsHaplogroupMt = false,
    supportsAutosomalIbd = false, // Autosomal not delivered
    supportsAncestry = false, // Autosomal not delivered
    typicalFileFormats = List("CSV", "TXT")
  )

  // Generic Y-SNP panel for other vendors or custom panels
  val YDNA_SNP_PANEL: TestTypeDefinition = TestTypeDefinition(
    code = "YDNA_SNP_PANEL",
    displayName = "Y-DNA SNP Panel",
    category = DataGenerationMethod.Genotyping,
    vendor = None,
    targetType = TargetType.YChromosome,
    expectedMinDepth = None,
    expectedTargetDepth = None,
    expectedMarkerCount = None,
    supportsHaplogroupY = true,
    supportsHaplogroupMt = false,
    supportsAutosomalIbd = false,
    supportsAncestry = false,
    typicalFileFormats = List("CSV", "TXT", "VCF")
  )

  // ============ SEQUENCING: Targeted Y-DNA ============

  val BIG_Y_500: TestTypeDefinition = TestTypeDefinition(
    code = "BIG_Y_500",
    displayName = "FTDNA Big Y-500",
    category = DataGenerationMethod.Sequencing,
    vendor = Some("FamilyTreeDNA"),
    targetType = TargetType.YChromosome,
    expectedMinDepth = Some(30.0),
    expectedTargetDepth = Some(50.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = true,
    supportsHaplogroupMt = false,
    supportsAutosomalIbd = false,
    supportsAncestry = false,
    typicalFileFormats = List("BAM", "VCF", "BED")
  )

  val BIG_Y_700: TestTypeDefinition = TestTypeDefinition(
    code = "BIG_Y_700",
    displayName = "FTDNA Big Y-700",
    category = DataGenerationMethod.Sequencing,
    vendor = Some("FamilyTreeDNA"),
    targetType = TargetType.YChromosome,
    expectedMinDepth = Some(30.0),
    expectedTargetDepth = Some(50.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = true,
    supportsHaplogroupMt = false,
    supportsAutosomalIbd = false,
    supportsAncestry = false,
    typicalFileFormats = List("BAM", "VCF", "BED")
  )

  val Y_ELITE: TestTypeDefinition = TestTypeDefinition(
    code = "Y_ELITE",
    displayName = "Full Genomes Y Elite",
    category = DataGenerationMethod.Sequencing,
    vendor = Some("Full Genomes"),
    targetType = TargetType.YChromosome,
    expectedMinDepth = Some(20.0),
    expectedTargetDepth = Some(30.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = true,
    supportsHaplogroupMt = true, // Off-target mtDNA reads provide sufficient coverage for precise haplogroup
    supportsAutosomalIbd = false,
    supportsAncestry = false,
    typicalFileFormats = List("BAM", "CRAM", "VCF")
  )

  val Y_PRIME: TestTypeDefinition = TestTypeDefinition(
    code = "Y_PRIME",
    displayName = "YSEQ Y-Prime",
    category = DataGenerationMethod.Sequencing,
    vendor = Some("YSEQ"),
    targetType = TargetType.YChromosome,
    expectedMinDepth = Some(20.0),
    expectedTargetDepth = Some(30.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = true,
    supportsHaplogroupMt = false,
    supportsAutosomalIbd = false,
    supportsAncestry = false,
    typicalFileFormats = List("BAM", "VCF")
  )

  // ============ SEQUENCING: Targeted mtDNA ============

  val MT_FULL_SEQUENCE: TestTypeDefinition = TestTypeDefinition(
    code = "MT_FULL_SEQUENCE",
    displayName = "mtDNA Full Sequence",
    category = DataGenerationMethod.Sequencing,
    vendor = None,
    targetType = TargetType.MtDna,
    expectedMinDepth = Some(500.0),
    expectedTargetDepth = Some(1000.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = false,
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = false,
    supportsAncestry = false,
    typicalFileFormats = List("BAM", "FASTA", "VCF")
  )

  val MT_PLUS: TestTypeDefinition = TestTypeDefinition(
    code = "MT_PLUS",
    displayName = "FTDNA mtDNA Plus",
    category = DataGenerationMethod.Sequencing,
    vendor = Some("FamilyTreeDNA"),
    targetType = TargetType.MtDna,
    expectedMinDepth = Some(500.0),
    expectedTargetDepth = Some(1000.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = false,
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = false,
    supportsAncestry = false,
    typicalFileFormats = List("FASTA", "VCF")
  )

  val MT_CR_ONLY: TestTypeDefinition = TestTypeDefinition(
    code = "MT_CR_ONLY",
    displayName = "mtDNA Control Region (HVR1/HVR2)",
    category = DataGenerationMethod.Sequencing,
    vendor = None,
    targetType = TargetType.MtDna,
    expectedMinDepth = Some(100.0),
    expectedTargetDepth = Some(500.0),
    expectedMarkerCount = None,
    supportsHaplogroupY = false,
    supportsHaplogroupMt = true, // Limited resolution
    supportsAutosomalIbd = false,
    supportsAncestry = false,
    typicalFileFormats = List("FASTA", "VCF")
  )

  // Genotyping array types
  val ARRAY_23ANDME_V5: TestTypeDefinition = TestTypeDefinition(
    code = "ARRAY_23ANDME_V5",
    displayName = "23andMe v5 Chip",
    category = DataGenerationMethod.Genotyping,
    vendor = Some("23andMe"),
    targetType = TargetType.Mixed,
    expectedMinDepth = None,
    expectedTargetDepth = None,
    expectedMarkerCount = Some(640000),
    supportsHaplogroupY = true,
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = true,
    supportsAncestry = true,
    typicalFileFormats = List("TXT")
  )

  val ARRAY_23ANDME_V4: TestTypeDefinition = TestTypeDefinition(
    code = "ARRAY_23ANDME_V4",
    displayName = "23andMe v4 Chip",
    category = DataGenerationMethod.Genotyping,
    vendor = Some("23andMe"),
    targetType = TargetType.Mixed,
    expectedMinDepth = None,
    expectedTargetDepth = None,
    expectedMarkerCount = Some(570000),
    supportsHaplogroupY = true,
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = true,
    supportsAncestry = true,
    typicalFileFormats = List("TXT")
  )

  val ARRAY_ANCESTRY_V2: TestTypeDefinition = TestTypeDefinition(
    code = "ARRAY_ANCESTRY_V2",
    displayName = "AncestryDNA v2",
    category = DataGenerationMethod.Genotyping,
    vendor = Some("AncestryDNA"),
    targetType = TargetType.Mixed,
    expectedMinDepth = None,
    expectedTargetDepth = None,
    expectedMarkerCount = Some(700000),
    supportsHaplogroupY = true,
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = true,
    supportsAncestry = true,
    typicalFileFormats = List("TXT")
  )

  val ARRAY_FTDNA_FF: TestTypeDefinition = TestTypeDefinition(
    code = "ARRAY_FTDNA_FF",
    displayName = "FTDNA Family Finder",
    category = DataGenerationMethod.Genotyping,
    vendor = Some("FamilyTreeDNA"),
    targetType = TargetType.Autosomal,
    expectedMinDepth = None,
    expectedTargetDepth = None,
    expectedMarkerCount = Some(700000),
    supportsHaplogroupY = false,
    supportsHaplogroupMt = false,
    supportsAutosomalIbd = true,
    supportsAncestry = true,
    typicalFileFormats = List("CSV")
  )

  val ARRAY_MYHERITAGE: TestTypeDefinition = TestTypeDefinition(
    code = "ARRAY_MYHERITAGE",
    displayName = "MyHeritage DNA",
    category = DataGenerationMethod.Genotyping,
    vendor = Some("MyHeritage"),
    targetType = TargetType.Mixed,
    expectedMinDepth = None,
    expectedTargetDepth = None,
    expectedMarkerCount = Some(700000),
    supportsHaplogroupY = true,
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = true,
    supportsAncestry = true,
    typicalFileFormats = List("CSV")
  )

  val ARRAY_LIVINGDNA: TestTypeDefinition = TestTypeDefinition(
    code = "ARRAY_LIVINGDNA",
    displayName = "LivingDNA",
    category = DataGenerationMethod.Genotyping,
    vendor = Some("LivingDNA"),
    targetType = TargetType.Mixed,
    expectedMinDepth = None,
    expectedTargetDepth = None,
    expectedMarkerCount = Some(630000),
    supportsHaplogroupY = true,
    supportsHaplogroupMt = true,
    supportsAutosomalIbd = true,
    supportsAncestry = true,
    typicalFileFormats = List("CSV", "TXT")
  )

  // ============ Collections ============

  /**
   * All known test types.
   */
  val all: List[TestTypeDefinition] = List(
    // Whole genome sequencing - short read (Illumina, BGI, Element, Ultima)
    WGS, WGS_LOW_PASS,
    // Whole genome sequencing - long read (PacBio, Nanopore)
    WGS_HIFI, WGS_NANOPORE, WGS_CLR,
    // Exome sequencing
    WES,
    // Targeted Y-DNA sequencing
    BIG_Y_500, BIG_Y_700, Y_ELITE, Y_PRIME,
    // Targeted mtDNA sequencing
    MT_FULL_SEQUENCE, MT_PLUS, MT_CR_ONLY,
    // Y-DNA SNP packs/panels (genotyping - probe-based or Sanger)
    YDNA_SNP_PACK_FTDNA, YDNA_PANEL_YSEQ, YDNA_SNP_PANEL,
    // Y-DNA focused arrays (genotyping arrays with Y emphasis)
    ARRAY_BISDNA,
    // Genotyping arrays (autosomal + mixed)
    ARRAY_23ANDME_V5, ARRAY_23ANDME_V4, ARRAY_ANCESTRY_V2,
    ARRAY_FTDNA_FF, ARRAY_MYHERITAGE, ARRAY_LIVINGDNA
  )

  /**
   * All sequencing test types.
   */
  val sequencing: List[TestTypeDefinition] =
    all.filter(_.category == DataGenerationMethod.Sequencing)

  /**
   * Short-read WGS types (Illumina, BGI, Element, Ultima).
   */
  val shortReadWgs: List[TestTypeDefinition] =
    List(WGS, WGS_LOW_PASS)

  /**
   * Long-read WGS types (PacBio HiFi, Nanopore, PacBio CLR).
   * These excel at resolving complex/repetitive regions including Y palindromes.
   */
  val longReadWgs: List[TestTypeDefinition] =
    List(WGS_HIFI, WGS_NANOPORE, WGS_CLR)

  /**
   * All genotyping array types.
   */
  val genotypingArrays: List[TestTypeDefinition] =
    all.filter(_.category == DataGenerationMethod.Genotyping)

  /**
   * Targeted Y-DNA sequencing types (Big Y, Y Elite, etc.).
   */
  val targetedYDnaSequencing: List[TestTypeDefinition] =
    sequencing.filter(_.targetType == TargetType.YChromosome)

  /**
   * Y-DNA SNP packs and panels (probe-based genotyping).
   * These are NOT sequencing - they test specific named SNPs.
   */
  val yDnaSnpPanels: List[TestTypeDefinition] =
    genotypingArrays.filter(_.targetType == TargetType.YChromosome)

  /**
   * All Y-DNA capable test types (sequencing + SNP panels).
   */
  val allYDnaTests: List[TestTypeDefinition] =
    targetedYDnaSequencing ++ yDnaSnpPanels

  /**
   * Targeted mtDNA sequencing types.
   */
  val targetedMtDna: List[TestTypeDefinition] =
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
