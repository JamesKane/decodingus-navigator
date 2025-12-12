package com.decodingus.yprofile.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.{
  BiosampleRepository, BiosampleEntity,
  YChromosomeProfileRepository, YProfileSourceRepository, YProfileRegionRepository,
  YProfileVariantRepository, YVariantSourceCallRepository, YVariantAuditRepository
}
import com.decodingus.yprofile.model.*
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
      YVariantAuditRepository()
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

    assertEquals(source.toOption.get.methodTier, 3) // WGS_SHORT_READ = tier 3
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
