package com.decodingus.yprofile.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.{BiosampleRepository, BiosampleEntity}
import com.decodingus.yprofile.model.*
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class YProfileRegionRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = YProfileRegionRepository()
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

  def createTestSource(profileId: UUID)(using java.sql.Connection): YProfileSourceEntity =
    sourceRepo.insert(YProfileSourceEntity.create(profileId, YProfileSourceType.WGS_SHORT_READ))

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)

      val entity = YProfileRegionEntity.create(
        yProfileId = profile.id,
        sourceId = source.id,
        startPosition = 1000000L,
        endPosition = 2000000L,
        callableState = YCallableState.CALLABLE,
        contig = "chrY",
        meanCoverage = Some(45.5),
        meanMappingQuality = Some(58.0)
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.startPosition, 1000000L)
      assertEquals(found.get.endPosition, 2000000L)
      assertEquals(found.get.callableState, YCallableState.CALLABLE)
      assertEquals(found.get.meanCoverage, Some(45.5))
    }
  }

  testTransactor.test("findByProfile returns all regions for profile") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)

      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 1000L, 2000L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 3000L, 4000L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 5000L, 6000L, YCallableState.LOW_COVERAGE))

      val regions = repo.findByProfile(profile.id)
      assertEquals(regions.size, 3)
      // Should be ordered by start_position
      assertEquals(regions.map(_.startPosition), List(1000L, 3000L, 5000L))
    }
  }

  testTransactor.test("findBySource returns regions for specific source") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source1 = createTestSource(profile.id)
      val source2 = createTestSource(profile.id)

      repo.insert(YProfileRegionEntity.create(profile.id, source1.id, 1000L, 2000L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source1.id, 3000L, 4000L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source2.id, 1000L, 2000L, YCallableState.LOW_COVERAGE))

      val source1Regions = repo.findBySource(source1.id)
      assertEquals(source1Regions.size, 2)
    }
  }

  testTransactor.test("findByState filters by callable state") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)

      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 1000L, 2000L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 3000L, 4000L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 5000L, 6000L, YCallableState.NO_COVERAGE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 7000L, 8000L, YCallableState.LOW_COVERAGE))

      val callable = repo.findByState(profile.id, YCallableState.CALLABLE)
      assertEquals(callable.size, 2)

      val noCoverage = repo.findByState(profile.id, YCallableState.NO_COVERAGE)
      assertEquals(noCoverage.size, 1)
    }
  }

  testTransactor.test("findOverlapping returns regions containing position") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)

      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 1000L, 2000L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 1500L, 2500L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 3000L, 4000L, YCallableState.CALLABLE))

      // Position 1800 should overlap with first two regions
      val overlapping = repo.findOverlapping(profile.id, 1800L)
      assertEquals(overlapping.size, 2)

      // Position 3500 should overlap with third region only
      val overlapping2 = repo.findOverlapping(profile.id, 3500L)
      assertEquals(overlapping2.size, 1)
    }
  }

  testTransactor.test("findOverlappingRange returns regions overlapping range") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)

      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 1000L, 2000L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 2500L, 3500L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 5000L, 6000L, YCallableState.CALLABLE))

      // Range 1500-3000 should overlap with first two regions
      val overlapping = repo.findOverlappingRange(profile.id, 1500L, 3000L)
      assertEquals(overlapping.size, 2)

      // Range 4000-5500 should overlap with third region only
      val overlapping2 = repo.findOverlappingRange(profile.id, 4000L, 5500L)
      assertEquals(overlapping2.size, 1)
    }
  }

  testTransactor.test("deleteByProfile removes all regions for profile") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)

      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 1000L, 2000L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 3000L, 4000L, YCallableState.CALLABLE))

      val deleted = repo.deleteByProfile(profile.id)
      assertEquals(deleted, 2)
      assertEquals(repo.findByProfile(profile.id).size, 0)
    }
  }

  testTransactor.test("deleteBySource removes all regions for source") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)

      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 1000L, 2000L, YCallableState.CALLABLE))
      repo.insert(YProfileRegionEntity.create(profile.id, source.id, 3000L, 4000L, YCallableState.CALLABLE))

      val deleted = repo.deleteBySource(source.id)
      assertEquals(deleted, 2)
      assertEquals(repo.findBySource(source.id).size, 0)
    }
  }

  testTransactor.test("insertBatch inserts multiple regions") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)

      val regions = List(
        YProfileRegionEntity.create(profile.id, source.id, 1000L, 2000L, YCallableState.CALLABLE),
        YProfileRegionEntity.create(profile.id, source.id, 3000L, 4000L, YCallableState.CALLABLE),
        YProfileRegionEntity.create(profile.id, source.id, 5000L, 6000L, YCallableState.LOW_COVERAGE)
      )

      val inserted = repo.insertBatch(regions)
      assertEquals(inserted.size, 3)
      assertEquals(repo.findByProfile(profile.id).size, 3)
    }
  }

  testTransactor.test("update modifies entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)

      val entity = repo.insert(YProfileRegionEntity.create(
        profile.id, source.id, 1000L, 2000L,
        YCallableState.LOW_COVERAGE,
        meanCoverage = Some(5.0)
      ))

      val updated = entity.copy(
        callableState = YCallableState.CALLABLE,
        meanCoverage = Some(45.0)
      )
      repo.update(updated)

      val found = repo.findById(entity.id)
      assertEquals(found.get.callableState, YCallableState.CALLABLE)
      assertEquals(found.get.meanCoverage, Some(45.0))
    }
  }

  testTransactor.test("delete removes entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val entity = repo.insert(YProfileRegionEntity.create(
        profile.id, source.id, 1000L, 2000L, YCallableState.CALLABLE
      ))

      assert(repo.exists(entity.id))
      repo.delete(entity.id)
      assert(!repo.exists(entity.id))
    }
  }

  testTransactor.test("cascades delete on profile deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val region = repo.insert(YProfileRegionEntity.create(
        profile.id, source.id, 1000L, 2000L, YCallableState.CALLABLE
      ))

      profileRepo.delete(profile.id)

      assertEquals(repo.findById(region.id), None)
    }
  }

  testTransactor.test("cascades delete on source deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = createTestSource(profile.id)
      val region = repo.insert(YProfileRegionEntity.create(
        profile.id, source.id, 1000L, 2000L, YCallableState.CALLABLE
      ))

      sourceRepo.delete(source.id)

      assertEquals(repo.findById(region.id), None)
    }
  }
