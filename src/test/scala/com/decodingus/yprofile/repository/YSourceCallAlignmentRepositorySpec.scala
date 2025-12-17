package com.decodingus.yprofile.repository

import com.decodingus.db.DatabaseTestSupport
import com.decodingus.repository.{BiosampleRepository, BiosampleEntity}
import com.decodingus.yprofile.model.*
import munit.FunSuite
import java.util.UUID

class YSourceCallAlignmentRepositorySpec extends FunSuite with DatabaseTestSupport:

  val alignmentRepo = new YSourceCallAlignmentRepository()
  val sourceCallRepo = new YVariantSourceCallRepository()
  val variantRepo = new YProfileVariantRepository()
  val sourceRepo = new YProfileSourceRepository()
  val profileRepo = new YChromosomeProfileRepository()
  val biosampleRepo = new BiosampleRepository()

  def createTestBiosample()(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = s"TEST-${UUID.randomUUID().toString.take(8)}",
      donorIdentifier = "test-donor",
      sex = Some("male")
    ))

  def createTestProfile(biosampleId: UUID)(using java.sql.Connection): YChromosomeProfileEntity =
    profileRepo.insert(YChromosomeProfileEntity.create(biosampleId))

  def createTestSource(profileId: UUID)(using java.sql.Connection): YProfileSourceEntity =
    sourceRepo.insert(YProfileSourceEntity.create(
      yProfileId = profileId,
      sourceType = YProfileSourceType.WGS_SHORT_READ,
      referenceBuild = Some("GRCh38")
    ))

  def createTestVariant(profileId: UUID)(using java.sql.Connection): YProfileVariantEntity =
    variantRepo.insert(YProfileVariantEntity.create(
      yProfileId = profileId,
      position = 2887824,
      refAllele = "G",
      altAllele = "A",
      variantName = Some("M269")
    ))

  def createTestSourceCall(variantId: UUID, sourceId: UUID)(using java.sql.Connection): YVariantSourceCallEntity =
    sourceCallRepo.insert(YVariantSourceCallEntity.create(
      variantId = variantId,
      sourceId = sourceId,
      calledAllele = "A",
      callState = YConsensusState.DERIVED
    ))

  // ============================================
  // Basic CRUD Tests
  // ============================================

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val sourceCall = createTestSourceCall(variant.id, source.id)

      val alignment = alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A",
        readDepth = Some(30),
        mappingQuality = Some(60.0)
      ))

      val found = alignmentRepo.findById(alignment.id)
      assert(found.isDefined)
      assertEquals(found.get.sourceCallId, sourceCall.id)
      assertEquals(found.get.referenceBuild, "GRCh38")
      assertEquals(found.get.position, 2887824L)
      assertEquals(found.get.readDepth, Some(30))
    }
  }

  testTransactor.test("update modifies entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val sourceCall = createTestSourceCall(variant.id, source.id)

      val alignment = alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))

      val updated = alignmentRepo.update(alignment.copy(readDepth = Some(50)))
      val found = alignmentRepo.findById(alignment.id)

      assertEquals(found.get.readDepth, Some(50))
    }
  }

  testTransactor.test("delete removes entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val sourceCall = createTestSourceCall(variant.id, source.id)

      val alignment = alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))

      assert(alignmentRepo.delete(alignment.id))
      assert(alignmentRepo.findById(alignment.id).isEmpty)
    }
  }

  // ============================================
  // Multi-Reference Tests
  // ============================================

  testTransactor.test("same source call can have multiple reference alignments") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val sourceCall = createTestSourceCall(variant.id, source.id)

      // Same source call aligned to different references
      val grch38 = alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))

      val grch37 = alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh37",
        position = 2793009,  // Different position in GRCh37!
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))

      val alignments = alignmentRepo.findBySourceCall(sourceCall.id)
      assertEquals(alignments.size, 2)

      // Verify different positions for same variant
      val grch38Alignment = alignments.find(_.referenceBuild == "GRCh38").get
      val grch37Alignment = alignments.find(_.referenceBuild == "GRCh37").get
      assertEquals(grch38Alignment.position, 2887824L)
      assertEquals(grch37Alignment.position, 2793009L)
    }
  }

  testTransactor.test("findBySourceCallAndBuild returns specific alignment") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val sourceCall = createTestSourceCall(variant.id, source.id)

      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))

      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh37",
        position = 2793009,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))

      val found = alignmentRepo.findBySourceCallAndBuild(sourceCall.id, "GRCh37")
      assert(found.isDefined)
      assertEquals(found.get.position, 2793009L)
    }
  }

  testTransactor.test("findByPosition returns alignments at position") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val sourceCall = createTestSourceCall(variant.id, source.id)

      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))

      val found = alignmentRepo.findByPosition("GRCh38", "chrY", 2887824)
      assertEquals(found.size, 1)
      assertEquals(found.head.calledAllele, "A")
    }
  }

  testTransactor.test("findByPositionRange returns alignments in range") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant1 = variantRepo.insert(YProfileVariantEntity.create(
        yProfileId = profile.id,
        position = 1000000,
        refAllele = "G",
        altAllele = "A"
      ))
      val variant2 = variantRepo.insert(YProfileVariantEntity.create(
        yProfileId = profile.id,
        position = 2000000,
        refAllele = "C",
        altAllele = "T"
      ))
      val variant3 = variantRepo.insert(YProfileVariantEntity.create(
        yProfileId = profile.id,
        position = 3000000,
        refAllele = "A",
        altAllele = "G"
      ))

      val sc1 = createTestSourceCall(variant1.id, source.id)
      val sc2 = createTestSourceCall(variant2.id, source.id)
      val sc3 = createTestSourceCall(variant3.id, source.id)

      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sc1.id, referenceBuild = "GRCh38", position = 1000000,
        refAllele = "G", altAllele = "A", calledAllele = "A"
      ))
      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sc2.id, referenceBuild = "GRCh38", position = 2000000,
        refAllele = "C", altAllele = "T", calledAllele = "T"
      ))
      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sc3.id, referenceBuild = "GRCh38", position = 3000000,
        refAllele = "A", altAllele = "G", calledAllele = "G"
      ))

      // Query range that includes only 2 of 3
      val found = alignmentRepo.findByPositionRange("GRCh38", "chrY", 1500000, 2500000)
      assertEquals(found.size, 1)
      assertEquals(found.head.position, 2000000L)
    }
  }

  testTransactor.test("upsert creates or updates based on source_call + build") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val sourceCall = createTestSourceCall(variant.id, source.id)

      // First upsert creates
      val created = alignmentRepo.upsert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A",
        readDepth = Some(30)
      ))

      // Second upsert with same source_call + build updates
      val updated = alignmentRepo.upsert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A",
        readDepth = Some(50)
      ))

      // Should only have one alignment
      val alignments = alignmentRepo.findBySourceCall(sourceCall.id)
      assertEquals(alignments.size, 1)
      assertEquals(alignments.head.readDepth, Some(50))

      // ID should be preserved from original
      assertEquals(updated.id, created.id)
    }
  }

  testTransactor.test("deleteBySourceCall removes all alignments for source call") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val sourceCall = createTestSourceCall(variant.id, source.id)

      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))
      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh37",
        position = 2793009,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))

      assertEquals(alignmentRepo.findBySourceCall(sourceCall.id).size, 2)

      val deleted = alignmentRepo.deleteBySourceCall(sourceCall.id)
      assertEquals(deleted, 2)
      assertEquals(alignmentRepo.findBySourceCall(sourceCall.id).size, 0)
    }
  }

  testTransactor.test("getDistinctReferenceBuilds returns all unique builds") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val sourceCall = createTestSourceCall(variant.id, source.id)

      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))
      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh37",
        position = 2793009,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))
      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "hs1",
        position = 2912345,
        refAllele = "C",  // Reverse complemented
        altAllele = "T",
        calledAllele = "T"
      ))

      val builds = alignmentRepo.getDistinctReferenceBuilds()
      assertEquals(builds.size, 3)
      assert(builds.contains("GRCh38"))
      assert(builds.contains("GRCh37"))
      assert(builds.contains("hs1"))
    }
  }

  testTransactor.test("cascades delete on source call deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val sourceCall = createTestSourceCall(variant.id, source.id)

      alignmentRepo.insert(YSourceCallAlignmentEntity.create(
        sourceCallId = sourceCall.id,
        referenceBuild = "GRCh38",
        position = 2887824,
        refAllele = "G",
        altAllele = "A",
        calledAllele = "A"
      ))

      // Delete the source call
      sourceCallRepo.delete(sourceCall.id)

      // Alignment should be gone due to cascade
      assertEquals(alignmentRepo.findBySourceCall(sourceCall.id).size, 0)
    }
  }
