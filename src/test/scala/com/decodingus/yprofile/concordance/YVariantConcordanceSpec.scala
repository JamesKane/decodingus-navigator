package com.decodingus.yprofile.concordance

import com.decodingus.yprofile.model.*
import munit.FunSuite
import java.util.UUID

class YVariantConcordanceSpec extends FunSuite:

  import YVariantConcordance.*

  // ============================================
  // Weight Calculation Tests
  // ============================================

  test("SNP weight: Sanger has highest base weight") {
    val sangerWeight = calculateWeight(YProfileSourceType.SANGER, YVariantType.SNP)
    val wgsWeight = calculateWeight(YProfileSourceType.WGS_SHORT_READ, YVariantType.SNP)
    val chipWeight = calculateWeight(YProfileSourceType.CHIP, YVariantType.SNP)

    assert(sangerWeight > wgsWeight)
    assert(wgsWeight > chipWeight)
    assertEquals(sangerWeight, 1.0)
    assertEquals(wgsWeight, 0.85)
    assertEquals(chipWeight, 0.5)
  }

  test("STR weight: Capillary electrophoresis has highest base weight") {
    val ceWeight = calculateWeight(YProfileSourceType.CAPILLARY_ELECTROPHORESIS, YVariantType.STR)
    val wgsWeight = calculateWeight(YProfileSourceType.WGS_SHORT_READ, YVariantType.STR)
    val chipWeight = calculateWeight(YProfileSourceType.CHIP, YVariantType.STR)

    assert(ceWeight > wgsWeight)
    assert(wgsWeight > chipWeight)
    assertEquals(ceWeight, 1.0)
    assertEquals(wgsWeight, 0.5)
    assertEquals(chipWeight, 0.3)
  }

  test("SNP vs STR: WGS is better for SNPs, CE is better for STRs") {
    val wgsSnpWeight = calculateWeight(YProfileSourceType.WGS_SHORT_READ, YVariantType.SNP)
    val wgsStrWeight = calculateWeight(YProfileSourceType.WGS_SHORT_READ, YVariantType.STR)
    val ceSnpWeight = calculateWeight(YProfileSourceType.CAPILLARY_ELECTROPHORESIS, YVariantType.SNP)
    val ceStrWeight = calculateWeight(YProfileSourceType.CAPILLARY_ELECTROPHORESIS, YVariantType.STR)

    // WGS should be better for SNPs
    assert(wgsSnpWeight > wgsStrWeight)
    // CE should be better for STRs
    assert(ceStrWeight > ceSnpWeight)
  }

  test("depth bonus increases weight for sequencing sources") {
    val noDepth = calculateWeight(YProfileSourceType.WGS_SHORT_READ, YVariantType.SNP)
    val lowDepth = calculateWeight(YProfileSourceType.WGS_SHORT_READ, YVariantType.SNP, readDepth = Some(10))
    val highDepth = calculateWeight(YProfileSourceType.WGS_SHORT_READ, YVariantType.SNP, readDepth = Some(100))

    assert(highDepth > lowDepth)
    assert(lowDepth > noDepth)
  }

  test("mapping quality affects weight") {
    val highMQ = calculateWeight(YProfileSourceType.WGS_SHORT_READ, YVariantType.SNP, mappingQuality = Some(60))
    val lowMQ = calculateWeight(YProfileSourceType.WGS_SHORT_READ, YVariantType.SNP, mappingQuality = Some(30))

    assert(highMQ > lowMQ)
    assertEqualsDouble(highMQ / lowMQ, 2.0, 0.01)
  }

  test("callable state affects weight") {
    val callable = calculateWeight(
      YProfileSourceType.WGS_SHORT_READ, YVariantType.SNP,
      callableState = Some(YCallableState.CALLABLE)
    )
    val lowCoverage = calculateWeight(
      YProfileSourceType.WGS_SHORT_READ, YVariantType.SNP,
      callableState = Some(YCallableState.LOW_COVERAGE)
    )
    val noCoverage = calculateWeight(
      YProfileSourceType.WGS_SHORT_READ, YVariantType.SNP,
      callableState = Some(YCallableState.NO_COVERAGE)
    )

    assert(callable > lowCoverage)
    assertEquals(noCoverage, 0.0)
  }

  // ============================================
  // Consensus Calculation Tests
  // ============================================

  test("unanimous DERIVED calls yield CONFIRMED status") {
    val calls = List(
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.WGS_SHORT_READ, "A", YConsensusState.DERIVED),
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.TARGETED_NGS, "A", YConsensusState.DERIVED),
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.CHIP, "A", YConsensusState.DERIVED)
    )

    val result = calculateConsensus(calls, YVariantType.SNP, isInTree = true)

    assertEquals(result.consensusAllele, Some("A"))
    assertEquals(result.consensusState, YConsensusState.DERIVED)
    assertEquals(result.status, YVariantStatus.CONFIRMED)
    assertEquals(result.concordantCount, 3)
    assertEquals(result.discordantCount, 0)
    assertEquals(result.confidenceScore, 1.0)
  }

  test("unanimous ANCESTRAL calls yield CONFIRMED status") {
    val calls = List(
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.WGS_SHORT_READ, "G", YConsensusState.ANCESTRAL),
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.TARGETED_NGS, "G", YConsensusState.ANCESTRAL)
    )

    val result = calculateConsensus(calls, YVariantType.SNP, isInTree = true)

    assertEquals(result.consensusAllele, Some("G"))
    assertEquals(result.consensusState, YConsensusState.ANCESTRAL)
    assertEquals(result.status, YVariantStatus.CONFIRMED)
  }

  test("novel variant (not in tree) yields NOVEL status") {
    val calls = List(
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.WGS_SHORT_READ, "A", YConsensusState.DERIVED),
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.TARGETED_NGS, "A", YConsensusState.DERIVED)
    )

    val result = calculateConsensus(calls, YVariantType.SNP, isInTree = false)

    assertEquals(result.status, YVariantStatus.NOVEL)
  }

  test("empty calls yield NO_COVERAGE status") {
    val result = calculateConsensus(List.empty, YVariantType.SNP)

    assertEquals(result.status, YVariantStatus.NO_COVERAGE)
    assertEquals(result.sourceCount, 0)
    assertEquals(result.confidenceScore, 0.0)
  }

  test("all NO_CALL states yield NO_COVERAGE status") {
    val calls = List(
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.WGS_SHORT_READ, ".", YConsensusState.NO_CALL),
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.CHIP, ".", YConsensusState.NO_CALL)
    )

    val result = calculateConsensus(calls, YVariantType.SNP)

    assertEquals(result.status, YVariantStatus.NO_COVERAGE)
    assertEquals(result.sourceCount, 2)  // Tracks total sources
  }

  test("high discordance yields CONFLICT status") {
    val calls = List(
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.WGS_SHORT_READ, "A", YConsensusState.DERIVED),
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.TARGETED_NGS, "G", YConsensusState.ANCESTRAL)
    )

    val result = calculateConsensus(calls, YVariantType.SNP)

    assertEquals(result.status, YVariantStatus.CONFLICT)
    assertEquals(result.discordantCount, 1)
  }

  // ============================================
  // User Scenario: 2 WGS vs 1 CE for SNP
  // ============================================

  test("SNP: 2 WGS DERIVED outweigh 1 CE ANCESTRAL") {
    // This is the exact scenario from the user's requirements:
    // "There is a SNP in my own profile that exists in short-read NGS and HiFi reads.
    //  The capillary results are ancestral. Will the two WGS approaches outrank here?"

    val calls = List(
      // Short-read WGS: DERIVED (weight ~0.85)
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.WGS_SHORT_READ, "A", YConsensusState.DERIVED,
        readDepth = Some(30), mappingQuality = Some(60)
      ),
      // Long-read WGS (HiFi): DERIVED (weight ~0.95)
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.WGS_LONG_READ, "A", YConsensusState.DERIVED,
        readDepth = Some(15), mappingQuality = Some(60)
      ),
      // CE: ANCESTRAL (weight ~0.5 for SNPs)
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.CAPILLARY_ELECTROPHORESIS, "G", YConsensusState.ANCESTRAL
      )
    )

    val result = calculateConsensus(calls, YVariantType.SNP, isInTree = true)

    // The two WGS sources should outweigh the CE source for SNPs
    assertEquals(result.consensusAllele, Some("A"))
    assertEquals(result.consensusState, YConsensusState.DERIVED)
    assertEquals(result.concordantCount, 2)
    assertEquals(result.discordantCount, 1)

    // Verify the weights make sense
    val derivedWeight = result.weightedCalls
      .filter(_._1.calledAllele == "A")
      .map(_._2).sum
    val ancestralWeight = result.weightedCalls
      .filter(_._1.calledAllele == "G")
      .map(_._2).sum

    // WGS weights combined should exceed CE weight for SNPs
    assert(derivedWeight > ancestralWeight,
      s"DERIVED weight ($derivedWeight) should exceed ANCESTRAL weight ($ancestralWeight)")
  }

  test("STR: CE outweighs WGS for repeat counting") {
    // For STRs, capillary electrophoresis should be trusted more
    val calls = List(
      // CE: 13 repeats (weight 1.0 for STRs)
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.CAPILLARY_ELECTROPHORESIS, "(GATA)13", YConsensusState.DERIVED,
        calledRepeatCount = Some(13)
      ),
      // WGS: 14 repeats (weight 0.5 for STRs)
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.WGS_SHORT_READ, "(GATA)14", YConsensusState.DERIVED,
        readDepth = Some(30), calledRepeatCount = Some(14)
      )
    )

    val result = calculateStrConsensus(calls, isInTree = true)

    // CE should win for STRs
    assertEquals(result.consensusAllele, Some("(GATA)13"))
  }

  // ============================================
  // STR Concordance Tests
  // ============================================

  test("STR: exact repeat count matches are concordant") {
    val calls = List(
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.CAPILLARY_ELECTROPHORESIS, "(GATA)13", YConsensusState.DERIVED,
        calledRepeatCount = Some(13)
      ),
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.CAPILLARY_ELECTROPHORESIS, "(GATA)13", YConsensusState.DERIVED,
        calledRepeatCount = Some(13)
      )
    )

    val result = calculateStrConsensus(calls)

    assertEquals(result.status, YVariantStatus.CONFIRMED)
    assertEquals(result.concordantCount, 2)
    assertEquals(result.discordantCount, 0)
  }

  test("STR: off-by-one repeat counts are concordant") {
    // Small differences in repeat count can occur due to stutter
    val calls = List(
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.CAPILLARY_ELECTROPHORESIS, "(GATA)13", YConsensusState.DERIVED,
        calledRepeatCount = Some(13)
      ),
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.CAPILLARY_ELECTROPHORESIS, "(GATA)14", YConsensusState.DERIVED,
        calledRepeatCount = Some(14)
      )
    )

    val result = calculateStrConsensus(calls)

    // Off-by-one should be concordant
    assertEquals(result.concordantCount, 2)
    assertEquals(result.discordantCount, 0)
  }

  test("STR: large repeat count differences are discordant") {
    val calls = List(
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.CAPILLARY_ELECTROPHORESIS, "(GATA)13", YConsensusState.DERIVED,
        calledRepeatCount = Some(13)
      ),
      SourceCallInput(
        UUID.randomUUID(), YProfileSourceType.WGS_SHORT_READ, "(GATA)16", YConsensusState.DERIVED,
        calledRepeatCount = Some(16)
      )
    )

    val result = calculateStrConsensus(calls)

    // Difference of 3 should be discordant
    assertEquals(result.discordantCount, 1)
  }

  // ============================================
  // Heteroplasmy Detection Tests
  // ============================================

  test("heteroplasmy detection for intermediate VAF values") {
    assert(!isLikelyHeteroplasmy(Some(0.0)))   // Homozygous ref
    assert(!isLikelyHeteroplasmy(Some(0.10)))  // Below threshold
    assert(isLikelyHeteroplasmy(Some(0.20)))   // Mixed signal
    assert(isLikelyHeteroplasmy(Some(0.50)))   // Mixed signal
    assert(isLikelyHeteroplasmy(Some(0.80)))   // Mixed signal
    assert(!isLikelyHeteroplasmy(Some(0.90)))  // Above threshold
    assert(!isLikelyHeteroplasmy(Some(1.0)))   // Homozygous alt
    assert(!isLikelyHeteroplasmy(None))        // No data
  }

  // ============================================
  // Profile Confidence Tests
  // ============================================

  test("profile confidence calculation") {
    // All confirmed
    assertEquals(calculateProfileConfidence(100, 0, 0, 100), 1.0)

    // All novel
    assertEqualsDouble(calculateProfileConfidence(0, 100, 0, 100), 0.7, 0.01)

    // Mixed confirmed and novel
    assertEqualsDouble(calculateProfileConfidence(50, 50, 0, 100), 0.85, 0.01)

    // With conflicts
    assert(calculateProfileConfidence(90, 5, 5, 100) < 1.0)

    // Empty profile
    assertEquals(calculateProfileConfidence(0, 0, 0, 0), 0.0)
  }

  // ============================================
  // Edge Cases
  // ============================================

  test("single source call yields consensus") {
    val calls = List(
      SourceCallInput(UUID.randomUUID(), YProfileSourceType.WGS_SHORT_READ, "A", YConsensusState.DERIVED)
    )

    val result = calculateConsensus(calls, YVariantType.SNP)

    assertEquals(result.consensusAllele, Some("A"))
    assertEquals(result.sourceCount, 1)
    assertEquals(result.concordantCount, 1)
    // Single source can still be CONFIRMED if high confidence
    assert(result.confidenceScore == 1.0)
  }

  test("weight calculation with all quality metrics") {
    val weight = calculateWeight(
      YProfileSourceType.WGS_SHORT_READ,
      YVariantType.SNP,
      readDepth = Some(100),      // sqrt(100)/10 = 1.0 depth bonus
      mappingQuality = Some(60),  // 60/60 = 1.0 factor
      callableState = Some(YCallableState.CALLABLE)  // 1.0 factor
    )

    // 0.85 * (1 + 1.0) * 1.0 * 1.0 = 1.7
    assertEqualsDouble(weight, 1.7, 0.01)
  }
