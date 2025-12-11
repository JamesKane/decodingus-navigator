package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.workspace.model.FileInfo
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class ChipProfileRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = ChipProfileRepository()
  val biosampleRepo = BiosampleRepository()

  def createTestBiosample(accession: String = s"TEST${UUID.randomUUID().toString.take(8)}")(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = accession,
      donorIdentifier = "DONOR001"
    ))

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = ChipProfileEntity.create(
        biosampleId = biosample.id,
        vendor = "23andMe",
        testTypeCode = "V5",
        chipVersion = Some("v5.2"),
        totalMarkersCalled = 650000,
        totalMarkersPossible = 700000,
        noCallRate = 0.071,
        autosomalMarkersCalled = 640000,
        importDate = LocalDateTime.now(),
        yMarkersCalled = Some(4000),
        mtMarkersCalled = Some(3000)
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.vendor, "23andMe")
      assertEquals(found.get.testTypeCode, "V5")
      assertEquals(found.get.chipVersion, Some("v5.2"))
      assertEquals(found.get.totalMarkersCalled, 650000)
      assertEquals(found.get.noCallRate, 0.071)
      assertEquals(found.get.yMarkersCalled, Some(4000))
      assertEquals(found.get.meta.syncStatus, SyncStatus.Local)
    }
  }

  testTransactor.test("findByBiosample returns all profiles for biosample") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val other = createTestBiosample()

      repo.insert(createChipProfile(biosample.id, "23andMe", "V5"))
      repo.insert(createChipProfile(biosample.id, "AncestryDNA", "V2"))
      repo.insert(createChipProfile(other.id, "FTDNA", "V3"))

      val profiles = repo.findByBiosample(biosample.id)
      assertEquals(profiles.size, 2)
      assert(profiles.forall(_.biosampleId == biosample.id))
    }
  }

  testTransactor.test("findByVendor filters by vendor") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()

      repo.insert(createChipProfile(biosample.id, "23andMe", "V5"))
      repo.insert(createChipProfile(biosample.id, "23andMe", "V4"))
      repo.insert(createChipProfile(biosample.id, "AncestryDNA", "V2"))

      val results = repo.findByVendor("23andMe")
      assertEquals(results.size, 2)
      assert(results.forall(_.vendor == "23andMe"))
    }
  }

  testTransactor.test("findByTestType filters by test type code") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()

      repo.insert(createChipProfile(biosample.id, "23andMe", "V5"))
      repo.insert(createChipProfile(biosample.id, "AncestryDNA", "V5"))
      repo.insert(createChipProfile(biosample.id, "23andMe", "V4"))

      val v5Results = repo.findByTestType("V5")
      assertEquals(v5Results.size, 2)
    }
  }

  testTransactor.test("findBySourceFileHash finds by hash for deduplication") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = createChipProfile(biosample.id, "23andMe", "V5")
        .copy(sourceFileHash = Some("abc123hash"))
      repo.insert(entity)

      val found = repo.findBySourceFileHash("abc123hash")
      assert(found.isDefined)

      val notFound = repo.findBySourceFileHash("nonexistent")
      assertEquals(notFound, None)
    }
  }

  testTransactor.test("update modifies entity and increments version") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = createChipProfile(biosample.id, "23andMe", "V5")
      val saved = repo.insert(entity)

      val updated = saved.copy(
        vendor = "Updated Vendor",
        noCallRate = 0.05
      )
      repo.update(updated)

      val found = repo.findById(saved.id)
      assert(found.isDefined)
      assertEquals(found.get.vendor, "Updated Vendor")
      assertEquals(found.get.noCallRate, 0.05)
      assertEquals(found.get.meta.version, 2)
    }
  }

  testTransactor.test("delete removes entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(createChipProfile(biosample.id, "23andMe", "V5"))

      assert(repo.exists(entity.id))
      val deleted = repo.delete(entity.id)
      assert(deleted)
      assert(!repo.exists(entity.id))
    }
  }

  testTransactor.test("markSynced updates sync status and AT fields") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(createChipProfile(biosample.id, "23andMe", "V5"))

      repo.markSynced(entity.id, "at://did:plc:test/chipprofile/1", "bafycid456")

      val found = repo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.meta.syncStatus, SyncStatus.Synced)
      assertEquals(found.get.meta.atUri, Some("at://did:plc:test/chipprofile/1"))
    }
  }

  testTransactor.test("cascades delete on biosample deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = repo.insert(createChipProfile(biosample.id, "23andMe", "V5"))

      biosampleRepo.delete(biosample.id)

      assertEquals(repo.findById(profile.id), None)
    }
  }

  private def createChipProfile(biosampleId: UUID, vendor: String, testType: String): ChipProfileEntity =
    ChipProfileEntity.create(
      biosampleId = biosampleId,
      vendor = vendor,
      testTypeCode = testType,
      totalMarkersCalled = 600000,
      totalMarkersPossible = 700000,
      noCallRate = 0.14,
      autosomalMarkersCalled = 590000,
      importDate = LocalDateTime.now()
    )
