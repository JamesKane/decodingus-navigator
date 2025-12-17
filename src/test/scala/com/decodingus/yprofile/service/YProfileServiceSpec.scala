package com.decodingus.yprofile.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.{BiosampleRepository, BiosampleEntity}
import com.decodingus.yprofile.repository.{
  YChromosomeProfileRepository, YProfileSourceRepository, YProfileRegionRepository,
  YProfileVariantRepository, YVariantSourceCallRepository, YVariantAuditRepository,
  YSourceCallAlignmentRepository
}
import com.decodingus.yprofile.model.*
import com.decodingus.yprofile.service.VariantCallInput
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class YProfileServiceSpec extends FunSuite with DatabaseTestSupport:

  def createService(tx: Transactor): YProfileService =
    YProfileService(
      tx,
      YChromosomeProfileRepository(),
      YProfileSourceRepository(),
      YProfileRegionRepository(),
      YProfileVariantRepository(),
      YVariantSourceCallRepository(),
      YVariantAuditRepository(),
      YSourceCallAlignmentRepository()
    )

  val biosampleRepo = BiosampleRepository()

  def createTestBiosample()(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = s"TEST${UUID.randomUUID().toString.take(8)}",
      donorIdentifier = "DONOR001"
    ))

  // ============================================
  // Profile Management Tests
  // ============================================

  testTransactor.test("getOrCreateProfile creates new profile") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get

    val result = service.getOrCreateProfile(biosampleId)
    assert(result.isRight)
    assertEquals(result.toOption.get.biosampleId, biosampleId)
  }

  testTransactor.test("getOrCreateProfile returns existing profile") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get

    val first = service.getOrCreateProfile(biosampleId)
    val second = service.getOrCreateProfile(biosampleId)

    assertEquals(first.toOption.get.id, second.toOption.get.id)
  }

  // ============================================
  // Source Management Tests
  // ============================================

  testTransactor.test("addSource creates source and updates profile count") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get

    val source = service.addSource(
      profile.id,
      YProfileSourceType.WGS_SHORT_READ,
      vendor = Some("Nebula"),
      testName = Some("30x WGS")
    )

    assert(source.isRight)
    assertEquals(source.toOption.get.vendor, Some("Nebula"))

    // Check profile was updated
    val updatedProfile = service.getProfile(profile.id).toOption.flatten.get
    assertEquals(updatedProfile.sourceCount, 1)
  }

  testTransactor.test("addSource sets correct method tier for WGS") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get

    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ)

    assertEquals(source.toOption.get.methodTier, 4) // WGS_SHORT_READ = tier 4 (0.85 * 5 = 4.25 -> 4)
  }

  testTransactor.test("removeSource deletes source and updates profile") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.CHIP).toOption.get

    val removed = service.removeSource(source.id)
    assert(removed.toOption.get)

    val updatedProfile = service.getProfile(profile.id).toOption.flatten.get
    assertEquals(updatedProfile.sourceCount, 0)
  }

  // ============================================
  // Variant Management Tests
  // ============================================

  testTransactor.test("addVariantCall creates variant and source call") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    val result = service.addVariantCall(
      profile.id, source.id,
      position = 2787994L,
      refAllele = "G",
      altAllele = "A",
      calledAllele = "A",
      callState = YConsensusState.DERIVED,
      variantName = Some("M343")
    )

    assert(result.isRight)
    val (variant, call) = result.toOption.get
    assertEquals(variant.position, 2787994L)
    assertEquals(variant.variantName, Some("M343"))
    assertEquals(call.calledAllele, "A")
    assertEquals(call.callState, YConsensusState.DERIVED)
  }

  testTransactor.test("addVariantCall adds call to existing variant") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source1 = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get
    val source2 = service.addSource(profile.id, YProfileSourceType.TARGETED_NGS).toOption.get

    // Add first call
    val result1 = service.addVariantCall(
      profile.id, source1.id,
      position = 2787994L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED
    )

    // Add second call for same position
    val result2 = service.addVariantCall(
      profile.id, source2.id,
      position = 2787994L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED
    )

    // Should reuse the same variant
    assertEquals(result1.toOption.get._1.id, result2.toOption.get._1.id)

    // Should have two source calls
    val calls = service.getVariantCalls(result1.toOption.get._1.id).toOption.get
    assertEquals(calls.size, 2)
  }

  // ============================================
  // Reconciliation Tests
  // ============================================

  testTransactor.test("reconcileVariant calculates consensus from source calls") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source1 = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get
    val source2 = service.addSource(profile.id, YProfileSourceType.WGS_LONG_READ).toOption.get

    // Add concordant calls
    val (variant, _) = service.addVariantCall(
      profile.id, source1.id,
      position = 2787994L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED
    ).toOption.get

    service.addVariantCall(
      profile.id, source2.id,
      position = 2787994L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED
    )

    val reconciled = service.reconcileVariant(variant.id, isInTree = true)

    assert(reconciled.isRight)
    val result = reconciled.toOption.get
    assertEquals(result.consensusAllele, Some("A"))
    assertEquals(result.consensusState, YConsensusState.DERIVED)
    assertEquals(result.status, YVariantStatus.CONFIRMED)
    assertEquals(result.concordantCount, 2)
    assertEquals(result.discordantCount, 0)
  }

  testTransactor.test("reconcileVariant detects conflicts") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source1 = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get
    val source2 = service.addSource(profile.id, YProfileSourceType.CHIP).toOption.get

    // Add conflicting calls
    val (variant, _) = service.addVariantCall(
      profile.id, source1.id,
      position = 2787994L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED
    ).toOption.get

    service.addVariantCall(
      profile.id, source2.id,
      position = 2787994L, refAllele = "G", altAllele = "A",
      calledAllele = "G", callState = YConsensusState.ANCESTRAL
    )

    val reconciled = service.reconcileVariant(variant.id)

    assert(reconciled.isRight)
    val result = reconciled.toOption.get
    // WGS should win over CHIP for SNPs
    assertEquals(result.consensusAllele, Some("A"))
    // But with discordance, status depends on threshold
    assertEquals(result.discordantCount, 1)
  }

  testTransactor.test("reconcileProfile updates all variants and statistics") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    // Add multiple variants
    service.addVariantCall(profile.id, source.id, 1000L, "G", "A", "A", YConsensusState.DERIVED)
    service.addVariantCall(profile.id, source.id, 2000L, "C", "T", "T", YConsensusState.DERIVED)
    service.addVariantCall(profile.id, source.id, 3000L, "A", "G", "G", YConsensusState.DERIVED)

    val count = service.reconcileProfile(profile.id)

    assertEquals(count.toOption.get, 3)

    val updatedProfile = service.getProfile(profile.id).toOption.flatten.get
    assertEquals(updatedProfile.totalVariants, 3)
    assert(updatedProfile.lastReconciledAt.isDefined)
  }

  // ============================================
  // Manual Curation Tests
  // ============================================

  testTransactor.test("overrideVariant creates audit entry and updates variant") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    val (variant, _) = service.addVariantCall(
      profile.id, source.id,
      position = 1000L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED
    ).toOption.get

    // Override the variant
    val overridden = service.overrideVariant(
      variant.id,
      newConsensusAllele = "G",
      newConsensusState = YConsensusState.ANCESTRAL,
      newStatus = YVariantStatus.CONFIRMED,
      reason = "IGV inspection shows ancestral state",
      userId = Some("user@example.com")
    )

    assert(overridden.isRight)
    val result = overridden.toOption.get
    assertEquals(result.consensusAllele, Some("G"))
    assertEquals(result.consensusState, YConsensusState.ANCESTRAL)
    assertEquals(result.status, YVariantStatus.CONFIRMED)
    assertEquals(result.confidenceScore, 1.0)

    // Check audit trail
    val audits = service.getAuditHistory(variant.id).toOption.get
    assertEquals(audits.size, 1)
    assertEquals(audits.head.action, YAuditAction.OVERRIDE)
    assertEquals(audits.head.reason, "IGV inspection shows ancestral state")
  }

  testTransactor.test("revertOverride recalculates consensus") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    val (variant, _) = service.addVariantCall(
      profile.id, source.id,
      position = 1000L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED
    ).toOption.get

    // Override
    service.overrideVariant(
      variant.id, "G", YConsensusState.ANCESTRAL, YVariantStatus.CONFIRMED,
      "Wrong override"
    )

    // Revert
    val reverted = service.revertOverride(
      variant.id,
      reason = "Override was incorrect",
      isInTree = true
    )

    assert(reverted.isRight)
    val result = reverted.toOption.get
    // Should be back to original consensus from source call
    assertEquals(result.consensusAllele, Some("A"))
    assertEquals(result.consensusState, YConsensusState.DERIVED)

    // Check audit trail has both entries
    val audits = service.getAuditHistory(variant.id).toOption.get
    assertEquals(audits.size, 2)
    assertEquals(audits.head.action, YAuditAction.REVERT)
  }

  // ============================================
  // User Scenario: 2 WGS vs 1 CE
  // ============================================

  testTransactor.test("SNP: 2 WGS sources outweigh 1 CE source") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get

    // Add sources
    val wgsShort = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ,
      vendor = Some("Nebula")).toOption.get
    val wgsLong = service.addSource(profile.id, YProfileSourceType.WGS_LONG_READ,
      vendor = Some("PacBio")).toOption.get
    val ce = service.addSource(profile.id, YProfileSourceType.CAPILLARY_ELECTROPHORESIS,
      vendor = Some("FTDNA")).toOption.get

    // WGS sources say DERIVED
    val (variant, _) = service.addVariantCall(
      profile.id, wgsShort.id,
      position = 2787994L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED,
      readDepth = Some(30), mappingQuality = Some(60)
    ).toOption.get

    service.addVariantCall(
      profile.id, wgsLong.id,
      position = 2787994L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED,
      readDepth = Some(15), mappingQuality = Some(60)
    )

    // CE says ANCESTRAL
    service.addVariantCall(
      profile.id, ce.id,
      position = 2787994L, refAllele = "G", altAllele = "A",
      calledAllele = "G", callState = YConsensusState.ANCESTRAL
    )

    // Reconcile
    val reconciled = service.reconcileVariant(variant.id, isInTree = true)

    assert(reconciled.isRight)
    val result = reconciled.toOption.get

    // WGS should win for SNPs
    assertEquals(result.consensusAllele, Some("A"),
      "Two WGS sources should outweigh one CE source for SNP")
    assertEquals(result.consensusState, YConsensusState.DERIVED)
    assertEquals(result.concordantCount, 2)
    assertEquals(result.discordantCount, 1)
  }

  testTransactor.test("STR: CE source outweighs WGS source") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get

    // Add sources
    val ce = service.addSource(profile.id, YProfileSourceType.CAPILLARY_ELECTROPHORESIS,
      vendor = Some("FTDNA")).toOption.get
    val wgs = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ,
      vendor = Some("Nebula")).toOption.get

    // CE says 13 repeats
    val (variant, _) = service.addVariantCall(
      profile.id, ce.id,
      position = 12500000L, refAllele = "(GATA)13", altAllele = "(GATA)14",
      calledAllele = "(GATA)13", callState = YConsensusState.DERIVED,
      variantType = YVariantType.STR,
      markerName = Some("DYS393"),
      repeatCount = Some(13)
    ).toOption.get

    // WGS says 14 repeats
    service.addVariantCall(
      profile.id, wgs.id,
      position = 12500000L, refAllele = "(GATA)13", altAllele = "(GATA)14",
      calledAllele = "(GATA)14", callState = YConsensusState.DERIVED,
      variantType = YVariantType.STR,
      markerName = Some("DYS393"),
      repeatCount = Some(14),
      readDepth = Some(30)
    )

    // Reconcile
    val reconciled = service.reconcileVariant(variant.id, isInTree = true)

    assert(reconciled.isRight)
    val result = reconciled.toOption.get

    // CE should win for STRs
    assertEquals(result.consensusAllele, Some("(GATA)13"),
      "CE source should outweigh WGS source for STR")
  }

  // ============================================
  // Callable Loci Integration Tests
  // ============================================

  testTransactor.test("importCallableIntervals stores regions and updates profile pct") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    // Import some callable intervals
    val intervals = List(
      (1000000L, 5000000L, YCallableState.CALLABLE),      // 4M bases callable
      (5000001L, 6000000L, YCallableState.LOW_COVERAGE),  // 1M low coverage
      (6000001L, 10000000L, YCallableState.CALLABLE)      // 4M bases callable
    )

    val result = service.importCallableIntervals(profile.id, source.id, intervals)

    assert(result.isRight)
    assertEquals(result.toOption.get, 3)

    // Check profile has callable region pct updated
    val updatedProfile = service.getProfile(profile.id).toOption.flatten.get
    assert(updatedProfile.callableRegionPct.isDefined)
    assert(updatedProfile.callableRegionPct.get > 0)
  }

  testTransactor.test("queryCallableState returns stored state for position") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    // Import intervals
    val intervals = List(
      (1000000L, 2000000L, YCallableState.CALLABLE),
      (2000001L, 3000000L, YCallableState.LOW_COVERAGE),
      (3000001L, 4000000L, YCallableState.NO_COVERAGE)
    )
    service.importCallableIntervals(profile.id, source.id, intervals)

    // Query positions in different regions
    val callable = service.queryCallableState(profile.id, 1500000L)
    val lowCov = service.queryCallableState(profile.id, 2500000L)
    val noCov = service.queryCallableState(profile.id, 3500000L)
    val unknown = service.queryCallableState(profile.id, 100L) // Outside any region

    assertEquals(callable.toOption.get, YCallableState.CALLABLE)
    assertEquals(lowCov.toOption.get, YCallableState.LOW_COVERAGE)
    assertEquals(noCov.toOption.get, YCallableState.NO_COVERAGE)
    assertEquals(unknown.toOption.get, YCallableState.NO_COVERAGE) // Default for unknown
  }

  testTransactor.test("queryCallableStates returns map for multiple positions") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    val intervals = List(
      (1000000L, 2000000L, YCallableState.CALLABLE),
      (2000001L, 3000000L, YCallableState.LOW_COVERAGE)
    )
    service.importCallableIntervals(profile.id, source.id, intervals)

    val positions = List(1500000L, 2500000L, 5000000L)
    val states = service.queryCallableStates(profile.id, positions)

    assert(states.isRight)
    val statesMap = states.toOption.get
    assertEquals(statesMap(1500000L), YCallableState.CALLABLE)
    assertEquals(statesMap(2500000L), YCallableState.LOW_COVERAGE)
    assertEquals(statesMap(5000000L), YCallableState.NO_COVERAGE)
  }

  testTransactor.test("getCallableSummary returns correct statistics") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    val intervals = List(
      (1000000L, 5000000L, YCallableState.CALLABLE),      // 4M bases
      (5000001L, 6000000L, YCallableState.LOW_COVERAGE),  // 1M bases
      (6000001L, 7000000L, YCallableState.NO_COVERAGE),   // 1M bases
      (7000001L, 10000000L, YCallableState.CALLABLE)      // 3M bases
    )
    service.importCallableIntervals(profile.id, source.id, intervals)

    val summary = service.getCallableSummary(source.id)

    assert(summary.isRight)
    val s = summary.toOption.get
    assertEquals(s.regionCount, 4)
    assertEquals(s.callableRegionCount, 2)
    assertEquals(s.lowCoverageRegionCount, 1)
    assertEquals(s.noCoverageRegionCount, 1)
    assertEquals(s.callableBases, 7000001L) // 4000001 + 3000000 from inclusive range
    assertEquals(s.totalBases, 9000001L)    // 4000001 + 1000000 + 1000000 + 3000000
  }

  testTransactor.test("callable state affects variant concordance weight") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    // Add variant with callable state = LOW_COVERAGE
    val (variant, sourceCall) = service.addVariantCall(
      profile.id, source.id,
      position = 2000000L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED,
      readDepth = Some(30), mappingQuality = Some(60),
      callableState = Some(YCallableState.LOW_COVERAGE)
    ).toOption.get

    // The weight should be reduced due to LOW_COVERAGE
    // (callableFactor = 0.5 instead of 1.0)
    assert(sourceCall.concordanceWeight < 0.85) // Base WGS weight without penalty
  }

  // ============================================
  // Alignment Record Tests
  // ============================================

  testTransactor.test("addVariantCall creates alignment record") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(
      profile.id, YProfileSourceType.WGS_SHORT_READ,
      referenceBuild = Some("GRCh38")
    ).toOption.get

    val result = service.addVariantCall(
      profile.id, source.id,
      position = 2887824L,
      refAllele = "G",
      altAllele = "A",
      calledAllele = "A",
      callState = YConsensusState.DERIVED,
      variantName = Some("M269"),
      readDepth = Some(30),
      mappingQuality = Some(60.0)
    )

    assert(result.isRight)
    val (variant, sourceCall) = result.toOption.get

    // Verify alignment was created
    val alignments = service.getAlignments(sourceCall.id)
    assert(alignments.isRight)
    assertEquals(alignments.toOption.get.size, 1)

    val alignment = alignments.toOption.get.head
    assertEquals(alignment.referenceBuild, "GRCh38")
    assertEquals(alignment.position, 2887824L)
    assertEquals(alignment.refAllele, "G")
    assertEquals(alignment.altAllele, "A")
    assertEquals(alignment.calledAllele, "A")
    assertEquals(alignment.readDepth, Some(30))
  }

  testTransactor.test("addAlignmentToSourceCall adds multi-reference alignment") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    // Add variant (creates GRCh38 alignment by default)
    val (variant, sourceCall) = service.addVariantCall(
      profile.id, source.id,
      position = 2887824L,
      refAllele = "G", altAllele = "A",
      calledAllele = "A",
      callState = YConsensusState.DERIVED,
      referenceBuild = Some("GRCh38")
    ).toOption.get

    // Add GRCh37 alignment (same data, different coordinates)
    val grch37Result = service.addAlignmentToSourceCall(
      sourceCallId = sourceCall.id,
      referenceBuild = "GRCh37",
      position = 2793009L,  // Different position in GRCh37
      refAllele = "G",
      altAllele = "A",
      calledAllele = "A"
    )

    assert(grch37Result.isRight)

    // Verify both alignments exist
    val alignments = service.getAlignments(sourceCall.id)
    assertEquals(alignments.toOption.get.size, 2)

    // Verify specific build queries work
    val grch38 = service.getAlignmentForBuild(sourceCall.id, "GRCh38")
    val grch37 = service.getAlignmentForBuild(sourceCall.id, "GRCh37")

    assert(grch38.toOption.get.isDefined)
    assert(grch37.toOption.get.isDefined)
    assertEquals(grch38.toOption.get.get.position, 2887824L)
    assertEquals(grch37.toOption.get.get.position, 2793009L)
  }

  testTransactor.test("multiple alignments still count as one evidence") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    // Add variant with GRCh38 alignment
    val (variant, sourceCall) = service.addVariantCall(
      profile.id, source.id,
      position = 2887824L, refAllele = "G", altAllele = "A",
      calledAllele = "A", callState = YConsensusState.DERIVED,
      referenceBuild = Some("GRCh38")
    ).toOption.get

    // Add GRCh37 and hs1 alignments
    service.addAlignmentToSourceCall(
      sourceCall.id, "GRCh37", 2793009L, "G", "A", "A"
    )
    service.addAlignmentToSourceCall(
      sourceCall.id, "hs1", 2912345L, "C", "T", "T"  // Reverse complement
    )

    // Should have 3 alignments but only 1 source call
    val alignments = service.getAlignments(sourceCall.id)
    assertEquals(alignments.toOption.get.size, 3)

    val calls = service.getVariantCalls(variant.id)
    assertEquals(calls.toOption.get.size, 1) // Only ONE piece of evidence

    // Reconcile - should have source_count = 1
    val reconciled = service.reconcileVariant(variant.id, isInTree = true)
    assertEquals(reconciled.toOption.get.sourceCount, 1)
  }

  // ============================================
  // Source Import Tests
  // ============================================

  testTransactor.test("importVariantCalls imports list of variants") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    val calls = List(
      VariantCallInput(
        position = 2887824L,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A",
        derived = true,
        variantName = Some("M269")
      ),
      VariantCallInput(
        position = 2790184L,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A",
        derived = true,
        variantName = Some("L21")
      ),
      VariantCallInput(
        position = 15571716L,
        refAllele = "C",
        altAllele = "T",
        calledAllele = "C",
        derived = false  // Ancestral
      )
    )

    val result = service.importVariantCalls(profile.id, source.id, calls)

    assert(result.isRight)
    val importResult = result.toOption.get
    assertEquals(importResult.snpCallsImported, 3)
    assertEquals(importResult.errorsEncountered, 0)

    // Verify variants were created
    val variants = service.getVariants(profile.id)
    assertEquals(variants.toOption.get.size, 3)

    // Verify alignments were created
    val variantCalls = variants.toOption.get.flatMap { v =>
      service.getVariantCalls(v.id).toOption.get
    }
    assertEquals(variantCalls.size, 3)

    variantCalls.foreach { call =>
      val alignments = service.getAlignments(call.id)
      assertEquals(alignments.toOption.get.size, 1)
    }
  }

  testTransactor.test("importVariantCalls handles readDepth and quality") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = tx.readWrite { createTestBiosample().id }.toOption.get
    val profile = service.getOrCreateProfile(biosampleId).toOption.get
    val source = service.addSource(profile.id, YProfileSourceType.WGS_SHORT_READ).toOption.get

    val calls = List(
      VariantCallInput(
        position = 2887824L,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A",
        derived = true,
        variantName = Some("M269"),
        readDepth = Some(45),
        qualityScore = Some(99.5),
        mappingQuality = Some(60.0)
      )
    )

    val result = service.importVariantCalls(profile.id, source.id, calls)
    assert(result.isRight)

    // Verify quality metrics were captured
    val variants = service.getVariants(profile.id).toOption.get
    val variantCalls = service.getVariantCalls(variants.head.id).toOption.get

    assertEquals(variantCalls.head.readDepth, Some(45))
    assertEquals(variantCalls.head.qualityScore, Some(99.5))
    assertEquals(variantCalls.head.mappingQuality, Some(60.0))
  }
