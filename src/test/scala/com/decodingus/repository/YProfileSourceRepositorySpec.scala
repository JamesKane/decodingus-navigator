package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.yprofile.model.*
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class YProfileSourceRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = YProfileSourceRepository()
  val profileRepo = YChromosomeProfileRepository()
  val biosampleRepo = BiosampleRepository()
  val alignmentRepo = AlignmentRepository()
  val seqRunRepo = SequenceRunRepository()

  def createTestBiosample()(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = s"TEST${UUID.randomUUID().toString.take(8)}",
      donorIdentifier = "DONOR001"
    ))

  def createTestProfile(biosampleId: UUID)(using java.sql.Connection): YChromosomeProfileEntity =
    profileRepo.insert(YChromosomeProfileEntity.create(biosampleId = biosampleId))

  def createTestAlignment(biosampleId: UUID)(using java.sql.Connection): AlignmentEntity =
    val seqRun = seqRunRepo.insert(SequenceRunEntity.create(
      biosampleId = biosampleId,
      platform = "ILLUMINA",
      testType = "WGS"
    ))
    alignmentRepo.insert(AlignmentEntity.create(
      sequenceRunId = seqRun.id,
      referenceBuild = "GRCh38",
      aligner = "BWA-MEM2"
    ))

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      val entity = YProfileSourceEntity.create(
        yProfileId = profile.id,
        sourceType = YProfileSourceType.WGS_SHORT_READ,
        sourceRef = Some("/path/to/bam"),
        vendor = Some("Nebula"),
        testName = Some("30x WGS"),
        methodTier = 3,
        meanReadDepth = Some(35.5),
        meanMappingQuality = Some(58.0),
        coveragePct = Some(0.92),
        variantCount = 450,
        strMarkerCount = 0,
        novelVariantCount = 15
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.sourceType, YProfileSourceType.WGS_SHORT_READ)
      assertEquals(found.get.vendor, Some("Nebula"))
      assertEquals(found.get.methodTier, 3)
      assertEquals(found.get.meanReadDepth, Some(35.5))
      assertEquals(found.get.variantCount, 450)
    }
  }

  testTransactor.test("findByProfile returns all sources for profile") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.WGS_SHORT_READ, methodTier = 3))
      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.TARGETED_NGS, methodTier = 2))
      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.CHIP, methodTier = 1))

      val sources = repo.findByProfile(profile.id)
      assertEquals(sources.size, 3)
      // Should be ordered by method_tier DESC
      assertEquals(sources.head.methodTier, 3)
    }
  }

  testTransactor.test("findByType returns sources by type") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.WGS_SHORT_READ))
      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.WGS_SHORT_READ))
      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.CHIP))

      val wgs = repo.findByType(YProfileSourceType.WGS_SHORT_READ)
      assertEquals(wgs.size, 2)
    }
  }

  testTransactor.test("findByAlignment returns sources linked to alignment") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val alignment = createTestAlignment(biosample.id)

      val withAlignment = repo.insert(YProfileSourceEntity.create(
        profile.id,
        YProfileSourceType.WGS_SHORT_READ,
        alignmentId = Some(alignment.id)
      ))
      repo.insert(YProfileSourceEntity.create(
        profile.id,
        YProfileSourceType.CHIP,
        alignmentId = None
      ))

      val sources = repo.findByAlignment(alignment.id)
      assertEquals(sources.size, 1)
      assertEquals(sources.head.id, withAlignment.id)
    }
  }

  testTransactor.test("findByVendor filters by vendor") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.TARGETED_NGS, vendor = Some("FTDNA")))
      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.TARGETED_NGS, vendor = Some("FTDNA")))
      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.TARGETED_NGS, vendor = Some("YSEQ")))

      val ftdna = repo.findByVendor("FTDNA")
      assertEquals(ftdna.size, 2)
    }
  }

  testTransactor.test("deleteByProfile removes all sources for profile") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.WGS_SHORT_READ))
      repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.CHIP))

      assertEquals(repo.countByProfile(profile.id), 2L)

      val deleted = repo.deleteByProfile(profile.id)
      assertEquals(deleted, 2)
      assertEquals(repo.countByProfile(profile.id), 0L)
    }
  }

  testTransactor.test("update modifies entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val entity = repo.insert(YProfileSourceEntity.create(
        profile.id,
        YProfileSourceType.WGS_SHORT_READ,
        variantCount = 100
      ))

      val updated = entity.copy(variantCount = 500, strMarkerCount = 50)
      repo.update(updated)

      val found = repo.findById(entity.id)
      assertEquals(found.get.variantCount, 500)
      assertEquals(found.get.strMarkerCount, 50)
    }
  }

  testTransactor.test("delete removes entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val entity = repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.CHIP))

      assert(repo.exists(entity.id))
      repo.delete(entity.id)
      assert(!repo.exists(entity.id))
    }
  }

  testTransactor.test("cascades delete on profile deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val source = repo.insert(YProfileSourceEntity.create(profile.id, YProfileSourceType.CHIP))

      profileRepo.delete(profile.id)

      assertEquals(repo.findById(source.id), None)
    }
  }
