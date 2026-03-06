package com.decodingus.workspace.services

import com.decodingus.model.LibraryStats
import com.decodingus.workspace.model.{RecordMeta, SequenceRun}
import munit.FunSuite

/**
 * Tests for FingerprintMatchService - verifies matching logic for
 * detecting when alignment files belong to the same sequencing run
 * (e.g., same library aligned to different reference builds).
 */
class FingerprintMatchServiceSpec extends FunSuite:

  val service = new FingerprintMatchService()

  // ============================================
  // Test Data Helpers
  // ============================================

  private def createSequenceRun(
    runFingerprint: Option[String] = None,
    platformUnit: Option[String] = None,
    libraryId: Option[String] = None,
    sampleName: Option[String] = None
  ): SequenceRun = SequenceRun(
    atUri = Some("local://sequencerun/test"),
    meta = RecordMeta.initial,
    biosampleRef = "local://biosample/test",
    platformName = "ILLUMINA",
    testType = "WGS",
    runFingerprint = runFingerprint,
    platformUnit = platformUnit,
    libraryId = libraryId,
    sampleName = sampleName
  )

  private def createLibraryStats(
    platformUnit: Option[String] = None,
    libraryId: String = "Unknown",
    sampleName: String = "Unknown"
  ): LibraryStats = LibraryStats(
    readCount = 1000000,
    pairedReads = 800000,
    lengthDistribution = Map(150 -> 1000000),
    insertSizeDistribution = Map.empty,
    aligner = "BWA-MEM2",
    referenceBuild = "GRCh38",
    sampleName = sampleName,
    libraryId = libraryId,
    platformUnit = platformUnit,
    flowCells = Map.empty,
    instruments = Map.empty,
    mostFrequentInstrumentId = "NovaSeq",
    mostFrequentInstrument = "NovaSeq 6000",
    inferredPlatform = "ILLUMINA",
    platformCounts = Map.empty
  )

  // ============================================
  // Tier 1: Exact Fingerprint Match (HIGH)
  // ============================================

  test("findMatch returns HIGH confidence for exact fingerprint match") {
    val fingerprint = "LB:lib123_PU:flowcell.1.A01_SM:sample001"
    val existingRun = createSequenceRun(runFingerprint = Some(fingerprint))
    val candidateRuns = List((existingRun, 0))
    val libraryStats = createLibraryStats()

    val result = service.findMatch(candidateRuns, fingerprint, libraryStats)

    result match
      case FingerprintMatchResult.MatchFound(run, index, confidence) =>
        assertEquals(confidence, "HIGH")
        assertEquals(index, 0)
        assertEquals(run.runFingerprint, Some(fingerprint))
      case FingerprintMatchResult.NoMatch =>
        fail("Expected MatchFound but got NoMatch")
  }

  test("findMatch returns HIGH confidence when fingerprint matches - different reference builds") {
    // Scenario: User has GRCh38 BAM, now adding GRCh37 BAM from same sequencing run
    val fingerprint = "LB:WGS_LIB_001_PU:HVMK5DSX2.1.AAACGGCG_SM:Sample_A"
    val grch38Run = createSequenceRun(
      runFingerprint = Some(fingerprint),
      platformUnit = Some("HVMK5DSX2.1.AAACGGCG"),
      libraryId = Some("WGS_LIB_001"),
      sampleName = Some("Sample_A")
    )
    val candidateRuns = List((grch38Run, 0))

    // LibraryStats from the new GRCh37 file
    val grch37Stats = createLibraryStats(
      platformUnit = Some("HVMK5DSX2.1.AAACGGCG"),
      libraryId = "WGS_LIB_001",
      sampleName = "Sample_A"
    )

    val result = service.findMatch(candidateRuns, fingerprint, grch37Stats)

    result match
      case FingerprintMatchResult.MatchFound(_, _, confidence) =>
        assertEquals(confidence, "HIGH")
      case FingerprintMatchResult.NoMatch =>
        fail("Should match same sequencing run with different reference")
  }

  // ============================================
  // Tier 2: Platform Unit Match (HIGH)
  // ============================================

  test("findMatch returns HIGH confidence for PU match when fingerprint doesn't match") {
    val existingRun = createSequenceRun(
      runFingerprint = Some("old_fingerprint"),
      platformUnit = Some("HVMK5DSX2.1.AAACGGCG")
    )
    val candidateRuns = List((existingRun, 0))
    val libraryStats = createLibraryStats(
      platformUnit = Some("HVMK5DSX2.1.AAACGGCG")
    )

    val result = service.findMatch(candidateRuns, "new_fingerprint", libraryStats)

    result match
      case FingerprintMatchResult.MatchFound(_, _, confidence) =>
        assertEquals(confidence, "HIGH")
      case FingerprintMatchResult.NoMatch =>
        fail("Expected HIGH confidence match on PU")
  }

  test("findMatch by PU works for same run aligned to T2T-CHM13 vs GRCh38") {
    // Scenario: First file was GRCh38, second is T2T-CHM13 alignment
    val existingRun = createSequenceRun(
      runFingerprint = Some("fingerprint_grch38"),
      platformUnit = Some("FLOWCELL123.2.ACGTACGT")
    )
    val candidateRuns = List((existingRun, 0))

    // T2T alignment has different fingerprint hash but same PU
    val t2tStats = createLibraryStats(
      platformUnit = Some("FLOWCELL123.2.ACGTACGT"),
      libraryId = "LIB_001",
      sampleName = "SAMPLE_001"
    )

    val result = service.findMatch(candidateRuns, "fingerprint_t2t", t2tStats)

    result match
      case FingerprintMatchResult.MatchFound(_, _, confidence) =>
        assertEquals(confidence, "HIGH")
      case FingerprintMatchResult.NoMatch =>
        fail("Should match by PU for different reference builds")
  }

  // ============================================
  // Tier 3: Library ID + Sample Name Match (MEDIUM)
  // ============================================

  test("findMatch returns MEDIUM confidence for LB+SM match when no PU match") {
    val existingRun = createSequenceRun(
      libraryId = Some("LIB_001"),
      sampleName = Some("SAMPLE_001")
      // No PU or fingerprint set
    )
    val candidateRuns = List((existingRun, 0))
    val libraryStats = createLibraryStats(
      libraryId = "LIB_001",
      sampleName = "SAMPLE_001"
      // No PU
    )

    val result = service.findMatch(candidateRuns, "some_fingerprint", libraryStats)

    result match
      case FingerprintMatchResult.MatchFound(_, _, confidence) =>
        assertEquals(confidence, "MEDIUM")
      case FingerprintMatchResult.NoMatch =>
        fail("Expected MEDIUM confidence match on LB+SM")
  }

  test("findMatch returns NoMatch when LB matches but SM doesn't") {
    val existingRun = createSequenceRun(
      libraryId = Some("LIB_001"),
      sampleName = Some("SAMPLE_001")
    )
    val candidateRuns = List((existingRun, 0))
    val libraryStats = createLibraryStats(
      libraryId = "LIB_001",
      sampleName = "DIFFERENT_SAMPLE"  // Different sample
    )

    val result = service.findMatch(candidateRuns, "fingerprint", libraryStats)

    assertEquals(result, FingerprintMatchResult.NoMatch)
  }

  // ============================================
  // No Match Scenarios
  // ============================================

  test("findMatch returns NoMatch for completely different run") {
    val existingRun = createSequenceRun(
      runFingerprint = Some("fingerprint_A"),
      platformUnit = Some("FLOWCELL_A.1.AAAA"),
      libraryId = Some("LIB_A"),
      sampleName = Some("SAMPLE_A")
    )
    val candidateRuns = List((existingRun, 0))
    val libraryStats = createLibraryStats(
      platformUnit = Some("FLOWCELL_B.1.TTTT"),
      libraryId = "LIB_B",
      sampleName = "SAMPLE_B"
    )

    val result = service.findMatch(candidateRuns, "fingerprint_B", libraryStats)

    assertEquals(result, FingerprintMatchResult.NoMatch)
  }

  test("findMatch returns NoMatch when libraryStats has Unknown values") {
    val existingRun = createSequenceRun(
      libraryId = Some("LIB_001"),
      sampleName = Some("SAMPLE_001")
    )
    val candidateRuns = List((existingRun, 0))
    val libraryStats = createLibraryStats(
      libraryId = "Unknown",  // Unknown - won't match
      sampleName = "Unknown"
    )

    val result = service.findMatch(candidateRuns, "fingerprint", libraryStats)

    assertEquals(result, FingerprintMatchResult.NoMatch)
  }

  test("findMatch returns NoMatch when no candidate runs exist") {
    val libraryStats = createLibraryStats(
      platformUnit = Some("FLOWCELL.1.AAAA"),
      libraryId = "LIB_001",
      sampleName = "SAMPLE_001"
    )

    val result = service.findMatch(List.empty, "fingerprint", libraryStats)

    assertEquals(result, FingerprintMatchResult.NoMatch)
  }

  // ============================================
  // Multiple Candidates
  // ============================================

  test("findMatch returns correct index when matching second candidate") {
    val run1 = createSequenceRun(
      runFingerprint = Some("fingerprint_1"),
      platformUnit = Some("FLOWCELL_A.1.AAAA")
    )
    val run2 = createSequenceRun(
      runFingerprint = Some("fingerprint_2"),
      platformUnit = Some("FLOWCELL_B.1.TTTT")
    )
    val candidateRuns = List((run1, 0), (run2, 1))
    val libraryStats = createLibraryStats(
      platformUnit = Some("FLOWCELL_B.1.TTTT")
    )

    val result = service.findMatch(candidateRuns, "new_fingerprint", libraryStats)

    result match
      case FingerprintMatchResult.MatchFound(_, index, _) =>
        assertEquals(index, 1)  // Should match run2 at index 1
      case FingerprintMatchResult.NoMatch =>
        fail("Expected match on run2")
  }

  // ============================================
  // Confidence Level Helpers
  // ============================================

  test("canAutoGroup returns true for HIGH confidence") {
    assert(service.canAutoGroup("HIGH"))
  }

  test("canAutoGroup returns true for MEDIUM confidence") {
    assert(service.canAutoGroup("MEDIUM"))
  }

  test("canAutoGroup returns false for LOW confidence") {
    assert(!service.canAutoGroup("LOW"))
  }

  test("requiresUserConfirmation returns true only for LOW confidence") {
    assert(service.requiresUserConfirmation("LOW"))
    assert(!service.requiresUserConfirmation("MEDIUM"))
    assert(!service.requiresUserConfirmation("HIGH"))
  }
