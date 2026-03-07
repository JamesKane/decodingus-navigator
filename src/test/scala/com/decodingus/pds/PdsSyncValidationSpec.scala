package com.decodingus.pds

import com.decodingus.workspace.model.*
import munit.FunSuite

class PdsSyncValidationSpec extends FunSuite:

  // ============================================
  // Test Type Mapping
  // ============================================

  test("toSyncTestType maps WGS variants to WGS") {
    assertEquals(SequenceRun.toSyncTestType("WGS"), "WGS")
    assertEquals(SequenceRun.toSyncTestType("WGS_LOW_PASS"), "WGS")
    assertEquals(SequenceRun.toSyncTestType("WGS_HIFI"), "WGS")
    assertEquals(SequenceRun.toSyncTestType("WGS_NANOPORE"), "WGS")
    assertEquals(SequenceRun.toSyncTestType("WGS_CLR"), "WGS")
  }

  test("toSyncTestType maps exome variants to EXOME") {
    assertEquals(SequenceRun.toSyncTestType("WES"), "EXOME")
    assertEquals(SequenceRun.toSyncTestType("EXOME"), "EXOME")
  }

  test("toSyncTestType maps targeted Y-DNA to TARGETED") {
    assertEquals(SequenceRun.toSyncTestType("BIG_Y_500"), "TARGETED")
    assertEquals(SequenceRun.toSyncTestType("BIG_Y_700"), "TARGETED")
    assertEquals(SequenceRun.toSyncTestType("Y_ELITE"), "TARGETED")
    assertEquals(SequenceRun.toSyncTestType("Y_PRIME"), "TARGETED")
  }

  test("toSyncTestType maps targeted mtDNA to TARGETED") {
    assertEquals(SequenceRun.toSyncTestType("MT_FULL_SEQUENCE"), "TARGETED")
    assertEquals(SequenceRun.toSyncTestType("MT_PLUS"), "TARGETED")
    assertEquals(SequenceRun.toSyncTestType("MT_CR_ONLY"), "TARGETED")
  }

  test("toSyncTestType passes through direct matches") {
    assertEquals(SequenceRun.toSyncTestType("AMPLICON"), "AMPLICON")
    assertEquals(SequenceRun.toSyncTestType("RNA_SEQ"), "RNA_SEQ")
    assertEquals(SequenceRun.toSyncTestType("TARGETED"), "TARGETED")
  }

  test("toSyncTestType passes through unknown codes") {
    assertEquals(SequenceRun.toSyncTestType("CUSTOM_TYPE"), "CUSTOM_TYPE")
  }

  // ============================================
  // Biosample Validation
  // ============================================

  test("validateBiosample passes with all required fields") {
    val b = makeBiosample(atUri = Some("at://did:plc:test/bio/1"), citizenDid = Some("did:plc:test"), centerName = Some("Lab"))
    assert(PdsSyncValidation.validateBiosample(b).isRight)
  }

  test("validateBiosample fails when missing citizenDid") {
    val b = makeBiosample(atUri = Some("at://test"), citizenDid = None, centerName = Some("Lab"))
    val Left(errors) = PdsSyncValidation.validateBiosample(b): @unchecked
    assert(errors.exists(_.contains("citizenDid")))
  }

  test("validateBiosample fails when missing atUri") {
    val b = makeBiosample(atUri = None, citizenDid = Some("did"), centerName = Some("Lab"))
    val Left(errors) = PdsSyncValidation.validateBiosample(b): @unchecked
    assert(errors.exists(_.contains("atUri")))
  }

  test("validateBiosample collects multiple errors") {
    val b = makeBiosample(atUri = None, citizenDid = None, centerName = None)
    val Left(errors) = PdsSyncValidation.validateBiosample(b): @unchecked
    assertEquals(errors.size, 3)
  }

  // ============================================
  // SequenceRun Validation
  // ============================================

  test("validateSequenceRun passes with all required fields") {
    val sr = makeSequenceRun(atUri = Some("at://test"), files = List(FileInfo("test.bam", None, "BAM", None)))
    assert(PdsSyncValidation.validateSequenceRun(sr).isRight)
  }

  test("validateSequenceRun fails when files empty") {
    val sr = makeSequenceRun(atUri = Some("at://test"), files = List.empty)
    val Left(errors) = PdsSyncValidation.validateSequenceRun(sr): @unchecked
    assert(errors.exists(_.contains("files")))
  }

  // ============================================
  // PopulationBreakdown Validation
  // ============================================

  test("validatePopulationBreakdown passes with required fields") {
    val pb = makePopulationBreakdown(components = List(makeComponent()))
    assert(PdsSyncValidation.validatePopulationBreakdown(pb).isRight)
  }

  test("validatePopulationBreakdown fails when components empty") {
    val pb = makePopulationBreakdown(components = List.empty)
    val Left(errors) = PdsSyncValidation.validatePopulationBreakdown(pb): @unchecked
    assert(errors.exists(_.contains("components")))
  }

  // ============================================
  // Helpers
  // ============================================

  private def makeBiosample(
                             atUri: Option[String] = Some("at://test"),
                             citizenDid: Option[String] = Some("did:plc:test"),
                             centerName: Option[String] = Some("TestLab")
                           ): Biosample =
    Biosample(
      atUri = atUri,
      meta = RecordMeta.initial,
      sampleAccession = "TEST001",
      donorIdentifier = "DONOR001",
      citizenDid = citizenDid,
      centerName = centerName
    )

  private def makeSequenceRun(
                                atUri: Option[String] = Some("at://test"),
                                files: List[FileInfo] = List(FileInfo("test.bam", None, "BAM", None))
                              ): SequenceRun =
    SequenceRun(
      atUri = atUri,
      meta = RecordMeta.initial,
      biosampleRef = "at://test/bio/1",
      platformName = "ILLUMINA",
      testType = "WGS",
      files = files
    )

  private def makePopulationBreakdown(
                                       components: List[com.decodingus.ancestry.model.PopulationComponent] = List.empty
                                     ): PopulationBreakdown =
    PopulationBreakdown(
      meta = RecordMeta.initial,
      biosampleRef = "at://test/bio/1",
      analysisMethod = "PCA_PROJECTION_GMM",
      panelType = "aims",
      referencePopulations = "1000G_HGDP_v1",
      snpsAnalyzed = 5000,
      snpsWithGenotype = 4500,
      snpsMissing = 500,
      confidenceLevel = 0.9,
      components = components,
      superPopulationSummary = List.empty,
      pcaCoordinates = None,
      analysisDate = None,
      pipelineVersion = "1.0.0",
      referenceVersion = "v1"
    )

  private def makeComponent(): com.decodingus.ancestry.model.PopulationComponent =
    com.decodingus.ancestry.model.PopulationComponent(
      populationCode = "CEU",
      populationName = "Northwestern European",
      superPopulation = "European",
      percentage = 85.0,
      confidenceInterval = com.decodingus.ancestry.model.ConfidenceInterval(80.0, 90.0),
      rank = 1
    )
