package com.decodingus.yprofile.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.{BiosampleRepository, BiosampleEntity}
import com.decodingus.yprofile.model.*
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class YVariantSourceCallRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = YVariantSourceCallRepository()
  val variantRepo = YProfileVariantRepository()
  val sourceRepo = YProfileSourceRepository()
  val profileRepo = YChromosomeProfileRepository()
  val biosampleRepo = BiosampleRepository()

  def createTestBiosample()(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = s"TEST${UUID.randomUUID().toString.take(8)}",
      donorIdentifier = "DONOR001"
    ))

  def createTestProfile(biosampleId: UUID)(using java.sql.Connection): YChromosomeProfileEntity =
    profileRepo.insert(YChromosomeProfileEntity.create(biosampleId = biosampleId))

  def createTestSource(profileId: UUID, sourceType: YProfileSourceType = YProfileSourceType.WGS_SHORT_READ)(using java.sql.Connection): YProfileSourceEntity =
    sourceRepo.insert(YProfileSourceEntity.create(profileId, sourceType))

  def createTestVariant(profileId: UUID, position: Long = 1000L)(using java.sql.Connection): YProfileVariantEntity =
    variantRepo.insert(YProfileVariantEntity.create(profileId, position, "G", "A"))

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)

      val entity = YVariantSourceCallEntity.create(
        variantId = variant.id,
        sourceId = source.id,
        calledAllele = "A",
        callState = YConsensusState.DERIVED,
        readDepth = Some(45),
        qualityScore = Some(99.0),
        mappingQuality = Some(58.0),
        variantAlleleFrequency = Some(1.0),
        callableState = Some(YCallableState.CALLABLE),
        concordanceWeight = 0.85
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.calledAllele, "A")
      assertEquals(found.get.callState, YConsensusState.DERIVED)
      assertEquals(found.get.readDepth, Some(45))
      assertEquals(found.get.concordanceWeight, 0.85)
      assertEquals(found.get.callableState, Some(YCallableState.CALLABLE))
    }
  }

  testTransactor.test("insert STR call with repeat count") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id, YProfileSourceType.CAPILLARY_ELECTROPHORESIS)
      val variant = variantRepo.insert(YProfileVariantEntity.create(
        profile.id, 1000L, "(GATA)13", "(GATA)14",
        variantType = YVariantType.STR,
        markerName = Some("DYS393")
      ))

      val entity = YVariantSourceCallEntity.create(
        variantId = variant.id,
        sourceId = source.id,
        calledAllele = "(GATA)13",
        callState = YConsensusState.DERIVED,
        calledRepeatCount = Some(13),
        concordanceWeight = 1.0
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.calledRepeatCount, Some(13))
    }
  }

  testTransactor.test("findByVariant returns all calls for variant") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source1 = createTestSource(profile.id, YProfileSourceType.WGS_SHORT_READ)
      val source2 = createTestSource(profile.id, YProfileSourceType.TARGETED_NGS)
      val variant = createTestVariant(profile.id)

      repo.insert(YVariantSourceCallEntity.create(variant.id, source1.id, "A", YConsensusState.DERIVED, concordanceWeight = 0.85))
      repo.insert(YVariantSourceCallEntity.create(variant.id, source2.id, "A", YConsensusState.DERIVED, concordanceWeight = 0.75))

      val calls = repo.findByVariant(variant.id)
      assertEquals(calls.size, 2)
      // Should be ordered by concordance_weight DESC
      assertEquals(calls.head.concordanceWeight, 0.85)
    }
  }

  testTransactor.test("findBySource returns all calls from source") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant1 = createTestVariant(profile.id, 1000L)
      val variant2 = createTestVariant(profile.id, 2000L)

      repo.insert(YVariantSourceCallEntity.create(variant1.id, source.id, "A", YConsensusState.DERIVED))
      repo.insert(YVariantSourceCallEntity.create(variant2.id, source.id, "T", YConsensusState.ANCESTRAL))

      val calls = repo.findBySource(source.id)
      assertEquals(calls.size, 2)
    }
  }

  testTransactor.test("findByVariantAndSource returns specific call") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source1 = createTestSource(profile.id)
      val source2 = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)

      val call1 = repo.insert(YVariantSourceCallEntity.create(variant.id, source1.id, "A", YConsensusState.DERIVED))
      repo.insert(YVariantSourceCallEntity.create(variant.id, source2.id, "G", YConsensusState.ANCESTRAL))

      val found = repo.findByVariantAndSource(variant.id, source1.id)
      assert(found.isDefined)
      assertEquals(found.get.id, call1.id)
      assertEquals(found.get.calledAllele, "A")
    }
  }

  testTransactor.test("countByVariant returns call count") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source1 = createTestSource(profile.id)
      val source2 = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)

      repo.insert(YVariantSourceCallEntity.create(variant.id, source1.id, "A", YConsensusState.DERIVED))
      repo.insert(YVariantSourceCallEntity.create(variant.id, source2.id, "A", YConsensusState.DERIVED))

      assertEquals(repo.countByVariant(variant.id), 2L)
    }
  }

  testTransactor.test("sumWeightsForAllele calculates weight sum") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source1 = createTestSource(profile.id)
      val source2 = createTestSource(profile.id)
      val source3 = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)

      // Two sources call A with weights 0.85 and 0.75
      repo.insert(YVariantSourceCallEntity.create(variant.id, source1.id, "A", YConsensusState.DERIVED, concordanceWeight = 0.85))
      repo.insert(YVariantSourceCallEntity.create(variant.id, source2.id, "A", YConsensusState.DERIVED, concordanceWeight = 0.75))
      // One source calls G with weight 0.5
      repo.insert(YVariantSourceCallEntity.create(variant.id, source3.id, "G", YConsensusState.ANCESTRAL, concordanceWeight = 0.5))

      val aWeight = repo.sumWeightsForAllele(variant.id, "A")
      val gWeight = repo.sumWeightsForAllele(variant.id, "G")

      assertEquals(aWeight, 1.6)
      assertEquals(gWeight, 0.5)
    }
  }

  testTransactor.test("deleteByVariant removes all calls for variant") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)

      repo.insert(YVariantSourceCallEntity.create(variant.id, source.id, "A", YConsensusState.DERIVED))

      val deleted = repo.deleteByVariant(variant.id)
      assertEquals(deleted, 1)
      assertEquals(repo.countByVariant(variant.id), 0L)
    }
  }

  testTransactor.test("deleteBySource removes all calls from source") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant1 = createTestVariant(profile.id, 1000L)
      val variant2 = createTestVariant(profile.id, 2000L)

      repo.insert(YVariantSourceCallEntity.create(variant1.id, source.id, "A", YConsensusState.DERIVED))
      repo.insert(YVariantSourceCallEntity.create(variant2.id, source.id, "T", YConsensusState.ANCESTRAL))

      val deleted = repo.deleteBySource(source.id)
      assertEquals(deleted, 2)
      assertEquals(repo.findBySource(source.id).size, 0)
    }
  }

  testTransactor.test("update modifies entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)

      val entity = repo.insert(YVariantSourceCallEntity.create(
        variant.id, source.id, "A", YConsensusState.DERIVED,
        concordanceWeight = 0.5
      ))

      val updated = entity.copy(concordanceWeight = 0.95, readDepth = Some(50))
      repo.update(updated)

      val found = repo.findById(entity.id)
      assertEquals(found.get.concordanceWeight, 0.95)
      assertEquals(found.get.readDepth, Some(50))
    }
  }

  testTransactor.test("cascades delete on variant deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val call = repo.insert(YVariantSourceCallEntity.create(variant.id, source.id, "A", YConsensusState.DERIVED))

      variantRepo.delete(variant.id)

      assertEquals(repo.findById(call.id), None)
    }
  }

  testTransactor.test("cascades delete on source deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val variant = createTestVariant(profile.id)
      val call = repo.insert(YVariantSourceCallEntity.create(variant.id, source.id, "A", YConsensusState.DERIVED))

      sourceRepo.delete(source.id)

      assertEquals(repo.findById(call.id), None)
    }
  }
