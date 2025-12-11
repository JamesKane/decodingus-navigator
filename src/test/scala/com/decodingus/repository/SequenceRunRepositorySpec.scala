package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import munit.FunSuite
import java.util.UUID

class SequenceRunRepositorySpec extends FunSuite with DatabaseTestSupport:

  val biosampleRepo = BiosampleRepository()
  val seqRunRepo = SequenceRunRepository()

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      val entity = SequenceRunEntity.create(
        biosampleId = biosample.id,
        platform = "ILLUMINA",
        testType = "WGS",
        instrumentModel = Some("NovaSeq 6000"),
        libraryId = Some("LIB001"),
        libraryLayout = Some("PAIRED")
      )

      val saved = seqRunRepo.insert(entity)
      val found = seqRunRepo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.biosampleId, biosample.id)
      assertEquals(found.get.platform, "ILLUMINA")
      assertEquals(found.get.testType, "WGS")
      assertEquals(found.get.instrumentModel, Some("NovaSeq 6000"))
      assertEquals(found.get.libraryId, Some("LIB001"))
      assertEquals(found.get.libraryLayout, Some("PAIRED"))
      assertEquals(found.get.meta.syncStatus, SyncStatus.Local)
    }
  }

  testTransactor.test("findAll returns all sequence runs") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS"))
      seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "PACBIO", "WGS"))

      val all = seqRunRepo.findAll()
      assertEquals(all.size, 2)
    }
  }

  testTransactor.test("update modifies entity correctly") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))
      val entity = seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS"))

      val updated = entity.copy(
        instrumentModel = Some("NextSeq 2000"),
        totalReads = Some(1000000L),
        readLength = Some(150)
      )
      seqRunRepo.update(updated)

      val found = seqRunRepo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.instrumentModel, Some("NextSeq 2000"))
      assertEquals(found.get.totalReads, Some(1000000L))
      assertEquals(found.get.readLength, Some(150))
      assertEquals(found.get.meta.version, 2)
    }
  }

  testTransactor.test("delete removes sequence run") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))
      val entity = seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS"))

      assert(seqRunRepo.exists(entity.id))

      val deleted = seqRunRepo.delete(entity.id)
      assert(deleted)
      assert(!seqRunRepo.exists(entity.id))
    }
  }

  testTransactor.test("findByBiosample returns runs for biosample") { case (db, tx) =>
    tx.readWrite {
      val bs1 = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))
      val bs2 = biosampleRepo.insert(BiosampleEntity.create("BS002", "D2"))

      seqRunRepo.insert(SequenceRunEntity.create(bs1.id, "ILLUMINA", "WGS"))
      seqRunRepo.insert(SequenceRunEntity.create(bs1.id, "ILLUMINA", "WES"))
      seqRunRepo.insert(SequenceRunEntity.create(bs2.id, "PACBIO", "WGS"))

      val bs1Runs = seqRunRepo.findByBiosample(bs1.id)
      assertEquals(bs1Runs.size, 2)
      assert(bs1Runs.forall(_.biosampleId == bs1.id))

      val bs2Runs = seqRunRepo.findByBiosample(bs2.id)
      assertEquals(bs2Runs.size, 1)
    }
  }

  testTransactor.test("findByPlatform filters by platform") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS"))
      seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WES"))
      seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "PACBIO", "WGS"))

      val illuminaRuns = seqRunRepo.findByPlatform("ILLUMINA")
      assertEquals(illuminaRuns.size, 2)

      val pacbioRuns = seqRunRepo.findByPlatform("PACBIO")
      assertEquals(pacbioRuns.size, 1)
    }
  }

  testTransactor.test("updateMetrics updates sequencing metrics") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))
      val entity = seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS"))

      seqRunRepo.updateMetrics(
        entity.id,
        totalReads = Some(2000000L),
        pfReads = Some(1800000L),
        readLength = Some(150),
        meanInsertSize = Some(350.5)
      )

      val found = seqRunRepo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.totalReads, Some(2000000L))
      assertEquals(found.get.pfReads, Some(1800000L))
      assertEquals(found.get.readLength, Some(150))
      assertEquals(found.get.meanInsertSize, Some(350.5))
    }
  }

  testTransactor.test("findPendingSync returns Local and Modified runs") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))

      val local = seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS"))
      val synced = seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "PACBIO", "WGS"))

      seqRunRepo.markSynced(synced.id, "at://test/1", "cid1")

      val pending = seqRunRepo.findPendingSync()
      assertEquals(pending.size, 1)
      assertEquals(pending.head.id, local.id)
    }
  }

  testTransactor.test("foreign key to biosample is enforced") { case (db, tx) =>
    val result = tx.readWrite {
      val fakeId = UUID.randomUUID()
      seqRunRepo.insert(SequenceRunEntity.create(fakeId, "ILLUMINA", "WGS"))
    }

    assert(result.isLeft, "Should fail with foreign key violation")
  }

  testTransactor.test("cascade delete removes sequence runs when biosample deleted") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))
      val seqRun = seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS"))

      assert(seqRunRepo.exists(seqRun.id))

      biosampleRepo.delete(biosample.id)

      assert(!seqRunRepo.exists(seqRun.id), "Sequence run should be cascade deleted")
    }
  }
