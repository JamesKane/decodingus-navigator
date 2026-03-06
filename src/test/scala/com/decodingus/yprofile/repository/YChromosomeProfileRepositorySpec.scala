package com.decodingus.yprofile.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.{BiosampleRepository, BiosampleEntity, SyncStatus}
import com.decodingus.yprofile.model.*
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class YChromosomeProfileRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = YChromosomeProfileRepository()
  val biosampleRepo = BiosampleRepository()

  def createTestBiosample(accession: String = s"TEST${UUID.randomUUID().toString.take(8)}")(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = accession,
      donorIdentifier = "DONOR001"
    ))

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = YChromosomeProfileEntity.create(
        biosampleId = biosample.id,
        consensusHaplogroup = Some("R-BY140757"),
        haplogroupConfidence = Some(0.99),
        haplogroupTreeProvider = Some("ftdna"),
        haplogroupTreeVersion = Some("2024.1"),
        totalVariants = 500,
        confirmedCount = 450,
        novelCount = 30,
        conflictCount = 5,
        noCoverageCount = 15,
        strMarkerCount = 111,
        strConfirmedCount = 108,
        overallConfidence = Some(0.95),
        callableRegionPct = Some(0.87),
        meanCoverage = Some(45.5),
        sourceCount = 3,
        primarySourceType = Some(YProfileSourceType.WGS_SHORT_READ)
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.consensusHaplogroup, Some("R-BY140757"))
      assertEquals(found.get.haplogroupConfidence, Some(0.99))
      assertEquals(found.get.totalVariants, 500)
      assertEquals(found.get.confirmedCount, 450)
      assertEquals(found.get.novelCount, 30)
      assertEquals(found.get.strMarkerCount, 111)
      assertEquals(found.get.primarySourceType, Some(YProfileSourceType.WGS_SHORT_READ))
      assertEquals(found.get.meta.syncStatus, SyncStatus.Local)
    }
  }

  testTransactor.test("findByBiosample returns profile for biosample") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val other = createTestBiosample()

      repo.insert(YChromosomeProfileEntity.create(
        biosampleId = biosample.id,
        consensusHaplogroup = Some("R-M269")
      ))

      val found = repo.findByBiosample(biosample.id)
      assert(found.isDefined)
      assertEquals(found.get.consensusHaplogroup, Some("R-M269"))

      val notFound = repo.findByBiosample(other.id)
      assert(notFound.isEmpty)
    }
  }

  testTransactor.test("findByHaplogroup returns exact matches") { case (db, tx) =>
    tx.readWrite {
      val b1 = createTestBiosample()
      val b2 = createTestBiosample()
      val b3 = createTestBiosample()

      repo.insert(YChromosomeProfileEntity.create(b1.id, consensusHaplogroup = Some("R-M269")))
      repo.insert(YChromosomeProfileEntity.create(b2.id, consensusHaplogroup = Some("R-L21")))
      repo.insert(YChromosomeProfileEntity.create(b3.id, consensusHaplogroup = Some("R-M269")))

      val m269 = repo.findByHaplogroup("R-M269")
      assertEquals(m269.size, 2)
    }
  }

  testTransactor.test("findByHaplogroupBranch finds prefix matches") { case (db, tx) =>
    tx.readWrite {
      val b1 = createTestBiosample()
      val b2 = createTestBiosample()
      val b3 = createTestBiosample()
      val b4 = createTestBiosample()

      repo.insert(YChromosomeProfileEntity.create(b1.id, consensusHaplogroup = Some("R-M269")))
      repo.insert(YChromosomeProfileEntity.create(b2.id, consensusHaplogroup = Some("R-L21")))
      repo.insert(YChromosomeProfileEntity.create(b3.id, consensusHaplogroup = Some("R-BY12345")))
      repo.insert(YChromosomeProfileEntity.create(b4.id, consensusHaplogroup = Some("I1-M253")))

      val rBranch = repo.findByHaplogroupBranch("R-")
      assertEquals(rBranch.size, 3)

      val iBranch = repo.findByHaplogroupBranch("I1")
      assertEquals(iBranch.size, 1)
    }
  }

  testTransactor.test("findWithConflicts returns profiles with conflict_count > 0") { case (db, tx) =>
    tx.readWrite {
      val b1 = createTestBiosample()
      val b2 = createTestBiosample()

      repo.insert(YChromosomeProfileEntity.create(b1.id, conflictCount = 5))
      repo.insert(YChromosomeProfileEntity.create(b2.id, conflictCount = 0))

      val conflicts = repo.findWithConflicts()
      assertEquals(conflicts.size, 1)
      assertEquals(conflicts.head.conflictCount, 5)
    }
  }

  testTransactor.test("update modifies entity and increments version") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = YChromosomeProfileEntity.create(
        biosampleId = biosample.id,
        consensusHaplogroup = Some("R-M269")
      )
      val saved = repo.insert(entity)

      val updated = saved.copy(
        consensusHaplogroup = Some("R-L21"),
        confirmedCount = 100
      )
      repo.update(updated)

      val found = repo.findById(saved.id)
      assert(found.isDefined)
      assertEquals(found.get.consensusHaplogroup, Some("R-L21"))
      assertEquals(found.get.confirmedCount, 100)
      assertEquals(found.get.meta.version, 2)
    }
  }

  testTransactor.test("markReconciled updates timestamp") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(YChromosomeProfileEntity.create(biosample.id))

      assert(entity.lastReconciledAt.isEmpty)

      repo.markReconciled(entity.id)

      val found = repo.findById(entity.id)
      assert(found.get.lastReconciledAt.isDefined)
    }
  }

  testTransactor.test("delete removes entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(YChromosomeProfileEntity.create(biosample.id))

      assert(repo.exists(entity.id))
      val deleted = repo.delete(entity.id)
      assert(deleted)
      assert(!repo.exists(entity.id))
    }
  }

  testTransactor.test("markSynced updates sync status and AT fields") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(YChromosomeProfileEntity.create(biosample.id))

      repo.markSynced(entity.id, "at://did:plc:test/yprofile/1", "bafycid123")

      val found = repo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.meta.syncStatus, SyncStatus.Synced)
      assertEquals(found.get.meta.atUri, Some("at://did:plc:test/yprofile/1"))
    }
  }

  testTransactor.test("cascades delete on biosample deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = repo.insert(YChromosomeProfileEntity.create(biosample.id))

      biosampleRepo.delete(biosample.id)

      assertEquals(repo.findById(profile.id), None)
    }
  }

  testTransactor.test("unique constraint prevents multiple profiles per biosample") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      repo.insert(YChromosomeProfileEntity.create(biosample.id))

      val caught = intercept[Exception] {
        repo.insert(YChromosomeProfileEntity.create(biosample.id))
      }
      assert(caught.getMessage.contains("constraint") || caught.getMessage.contains("Unique"))
    }
  }
