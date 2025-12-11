package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.workspace.model.{
  DnaType, CompatibilityLevel, HaplogroupTechnology, CallMethod, ConflictResolution,
  ReconciliationStatus, RunHaplogroupCall, SnpCallFromRun, SnpConflict
}
import munit.FunSuite
import java.time.Instant
import java.util.UUID

class HaplogroupReconciliationRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = HaplogroupReconciliationRepository()
  val biosampleRepo = BiosampleRepository()

  def createTestBiosample(accession: String = s"TEST${UUID.randomUUID().toString.take(8)}")(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = accession,
      donorIdentifier = "DONOR001"
    ))

  def createTestStatus(haplogroup: String = "R-M269"): ReconciliationStatus =
    ReconciliationStatus(
      compatibilityLevel = CompatibilityLevel.COMPATIBLE,
      consensusHaplogroup = haplogroup,
      confidence = 0.95,
      runCount = 1
    )

  def createTestRunCall(haplogroup: String = "R-M269"): RunHaplogroupCall =
    RunHaplogroupCall(
      sourceRef = "at://test/source/1",
      haplogroup = haplogroup,
      confidence = 0.95,
      callMethod = CallMethod.SNP_PHYLOGENETIC,
      technology = Some(HaplogroupTechnology.WGS)
    )

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = createTestStatus(),
        runCalls = List(createTestRunCall()),
        lastReconciliationAt = Some(Instant.now())
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.dnaType, DnaType.Y_DNA)
      assertEquals(found.get.status.consensusHaplogroup, "R-M269")
      assertEquals(found.get.status.confidence, 0.95)
      assertEquals(found.get.runCalls.size, 1)
      assertEquals(found.get.meta.syncStatus, SyncStatus.Local)
    }
  }

  testTransactor.test("findByBiosampleAndDnaType returns unique reconciliation") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()

      val yDna = HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = createTestStatus("R-M269")
      )
      val mtDna = HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.MT_DNA,
        status = createTestStatus("H1a1")
      )

      repo.insert(yDna)
      repo.insert(mtDna)

      val foundY = repo.findByBiosampleAndDnaType(biosample.id, DnaType.Y_DNA)
      assert(foundY.isDefined)
      assertEquals(foundY.get.status.consensusHaplogroup, "R-M269")

      val foundMt = repo.findByBiosampleAndDnaType(biosample.id, DnaType.MT_DNA)
      assert(foundMt.isDefined)
      assertEquals(foundMt.get.status.consensusHaplogroup, "H1a1")
    }
  }

  testTransactor.test("findByBiosample returns all reconciliations for biosample") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()

      repo.insert(HaplogroupReconciliationEntity.create(biosample.id, DnaType.Y_DNA, createTestStatus()))
      repo.insert(HaplogroupReconciliationEntity.create(biosample.id, DnaType.MT_DNA, createTestStatus("H1")))

      val all = repo.findByBiosample(biosample.id)
      assertEquals(all.size, 2)
    }
  }

  testTransactor.test("findByDnaType filters by DNA type") { case (db, tx) =>
    tx.readWrite {
      val b1 = createTestBiosample()
      val b2 = createTestBiosample()

      repo.insert(HaplogroupReconciliationEntity.create(b1.id, DnaType.Y_DNA, createTestStatus()))
      repo.insert(HaplogroupReconciliationEntity.create(b1.id, DnaType.MT_DNA, createTestStatus("H1")))
      repo.insert(HaplogroupReconciliationEntity.create(b2.id, DnaType.Y_DNA, createTestStatus("I1")))

      val yDnas = repo.findByDnaType(DnaType.Y_DNA)
      assertEquals(yDnas.size, 2)

      val mtDnas = repo.findByDnaType(DnaType.MT_DNA)
      assertEquals(mtDnas.size, 1)
    }
  }

  testTransactor.test("unique constraint on biosample_id + dna_type") { case (db, tx) =>
    val result = tx.readWrite {
      val biosample = createTestBiosample()
      repo.insert(HaplogroupReconciliationEntity.create(biosample.id, DnaType.Y_DNA, createTestStatus()))
      // Second insert with same biosample + dna_type should fail
      repo.insert(HaplogroupReconciliationEntity.create(biosample.id, DnaType.Y_DNA, createTestStatus("R-L21")))
    }
    assert(result.isLeft, "Should fail on duplicate biosample_id + dna_type")
  }

  testTransactor.test("upsert inserts when not exists") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = createTestStatus()
      )

      val result = repo.upsert(entity)
      assertEquals(result.status.consensusHaplogroup, "R-M269")

      val found = repo.findByBiosampleAndDnaType(biosample.id, DnaType.Y_DNA)
      assert(found.isDefined)
    }
  }

  testTransactor.test("upsert updates when exists") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val initial = HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = createTestStatus("R-M269")
      )
      repo.insert(initial)

      val updated = HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = createTestStatus("R-L21")
      )
      val result = repo.upsert(updated)

      // Should use the same ID as the existing record
      val found = repo.findByBiosampleAndDnaType(biosample.id, DnaType.Y_DNA)
      assert(found.isDefined)
      assertEquals(found.get.status.consensusHaplogroup, "R-L21")
      // Count should still be 1
      assertEquals(repo.findByBiosample(biosample.id).filter(_.dnaType == DnaType.Y_DNA).size, 1)
    }
  }

  testTransactor.test("stores and retrieves run calls correctly") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val calls = List(
        createTestRunCall("R-M269"),
        RunHaplogroupCall(
          sourceRef = "at://test/source/2",
          haplogroup = "R-L21",
          confidence = 0.92,
          callMethod = CallMethod.SNP_PHYLOGENETIC,
          technology = Some(HaplogroupTechnology.BIG_Y)
        )
      )

      val entity = HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = createTestStatus(),
        runCalls = calls
      )
      repo.insert(entity)

      val found = repo.findByBiosampleAndDnaType(biosample.id, DnaType.Y_DNA)
      assert(found.isDefined)
      assertEquals(found.get.runCalls.size, 2)
      val haplogroups = found.get.runCalls.map(_.haplogroup).toSet
      assert(haplogroups.contains("R-M269"))
      assert(haplogroups.contains("R-L21"))
    }
  }

  testTransactor.test("stores and retrieves SNP conflicts correctly") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val conflicts = List(
        SnpConflict(
          position = 14722131,
          snpName = Some("L21"),
          calls = List(
            SnpCallFromRun("at://source/1", "A", Some(30.0), Some(100)),
            SnpCallFromRun("at://source/2", "G", Some(25.0), Some(80))
          ),
          resolution = Some(ConflictResolution.ACCEPT_HIGHER_QUALITY),
          resolvedValue = Some("A")
        )
      )

      val entity = HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = createTestStatus().copy(compatibilityLevel = CompatibilityLevel.MINOR_DIVERGENCE),
        snpConflicts = conflicts
      )
      repo.insert(entity)

      val found = repo.findByBiosampleAndDnaType(biosample.id, DnaType.Y_DNA)
      assert(found.isDefined)
      assertEquals(found.get.snpConflicts.size, 1)
      assertEquals(found.get.snpConflicts.head.snpName, Some("L21"))
      assertEquals(found.get.snpConflicts.head.resolution, Some(ConflictResolution.ACCEPT_HIGHER_QUALITY))
    }
  }

  testTransactor.test("delete removes entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = createTestStatus()
      ))

      assert(repo.exists(entity.id))
      val deleted = repo.delete(entity.id)
      assert(deleted)
      assert(!repo.exists(entity.id))
    }
  }

  testTransactor.test("cascades delete on biosample deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = createTestStatus()
      ))

      biosampleRepo.delete(biosample.id)

      assertEquals(repo.findById(entity.id), None)
    }
  }
