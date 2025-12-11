package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.workspace.model.{FileInfo, SimpleStrValue, StrMarkerValue, StrPanel}
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class StrProfileRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = StrProfileRepository()
  val biosampleRepo = BiosampleRepository()

  def createTestBiosample(accession: String = s"TEST${UUID.randomUUID().toString.take(8)}")(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = accession,
      donorIdentifier = "DONOR001"
    ))

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = StrProfileEntity.create(
        biosampleId = biosample.id,
        source = Some("DIRECT_TEST"),
        importedFrom = Some("FTDNA"),
        totalMarkers = Some(37),
        panels = List(StrPanel("Y-37", 37, Some("FTDNA"), None)),
        markers = List(
          StrMarkerValue("DYS393", SimpleStrValue(13), panel = Some("Y12")),
          StrMarkerValue("DYS390", SimpleStrValue(24), panel = Some("Y12"))
        )
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.source, Some("DIRECT_TEST"))
      assertEquals(found.get.importedFrom, Some("FTDNA"))
      assertEquals(found.get.totalMarkers, Some(37))
      assertEquals(found.get.panels.size, 1)
      assertEquals(found.get.markers.size, 2)
      assertEquals(found.get.meta.syncStatus, SyncStatus.Local)
    }
  }

  testTransactor.test("stores and retrieves STR markers with genomic positions") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val orderDate = LocalDateTime.of(2017, 2, 22, 0, 0)
      val entity = StrProfileEntity.create(
        biosampleId = biosample.id,
        source = Some("BIG_Y_DERIVED"),
        markers = List(
          // DYR33 from the sample file: 17388563-17388824, value=14
          StrMarkerValue(
            marker = "DYR33",
            value = SimpleStrValue(14),
            startPosition = Some(17388563L),
            endPosition = Some(17388824L),
            orderedDate = Some(orderDate),
            panel = Some("Big Y-700")
          ),
          // DYR76: 12080003-12080116, value=13
          StrMarkerValue(
            marker = "DYR76",
            value = SimpleStrValue(13),
            startPosition = Some(12080003L),
            endPosition = Some(12080116L),
            orderedDate = Some(orderDate),
            panel = Some("Big Y-700")
          )
        )
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.markers.size, 2)

      val dyr33 = found.get.markers.find(_.marker == "DYR33")
      assert(dyr33.isDefined)
      assertEquals(dyr33.get.startPosition, Some(17388563L))
      assertEquals(dyr33.get.endPosition, Some(17388824L))
      assertEquals(dyr33.get.orderedDate, Some(orderDate))
      assertEquals(dyr33.get.regionSpan, Some(262L))  // 17388824 - 17388563 + 1

      val dyr76 = found.get.markers.find(_.marker == "DYR76")
      assert(dyr76.isDefined)
      assertEquals(dyr76.get.startPosition, Some(12080003L))
      assertEquals(dyr76.get.endPosition, Some(12080116L))
      assertEquals(dyr76.get.regionSpan, Some(114L))  // 12080116 - 12080003 + 1
    }
  }

  testTransactor.test("findByBiosample returns all profiles for biosample") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val other = createTestBiosample("OTHER001")

      repo.insert(StrProfileEntity.create(biosample.id, source = Some("DIRECT_TEST")))
      repo.insert(StrProfileEntity.create(biosample.id, source = Some("WGS_DERIVED")))
      repo.insert(StrProfileEntity.create(other.id, source = Some("IMPORTED")))

      val profiles = repo.findByBiosample(biosample.id)
      assertEquals(profiles.size, 2)
      assert(profiles.forall(_.biosampleId == biosample.id))
    }
  }

  testTransactor.test("findBySource filters by source type") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()

      repo.insert(StrProfileEntity.create(biosample.id, source = Some("DIRECT_TEST")))
      repo.insert(StrProfileEntity.create(biosample.id, source = Some("WGS_DERIVED")))
      repo.insert(StrProfileEntity.create(biosample.id, source = Some("DIRECT_TEST")))

      val direct = repo.findBySource("DIRECT_TEST")
      assertEquals(direct.size, 2)

      val wgs = repo.findBySource("WGS_DERIVED")
      assertEquals(wgs.size, 1)
    }
  }

  testTransactor.test("findByImportedFrom filters by provider") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()

      repo.insert(StrProfileEntity.create(biosample.id, importedFrom = Some("FTDNA")))
      repo.insert(StrProfileEntity.create(biosample.id, importedFrom = Some("YSEQ")))
      repo.insert(StrProfileEntity.create(biosample.id, importedFrom = Some("FTDNA")))

      val ftdna = repo.findByImportedFrom("FTDNA")
      assertEquals(ftdna.size, 2)

      val yseq = repo.findByImportedFrom("YSEQ")
      assertEquals(yseq.size, 1)
    }
  }

  testTransactor.test("update modifies entity and increments version") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = StrProfileEntity.create(biosample.id, source = Some("DIRECT_TEST"))
      val saved = repo.insert(entity)

      val updated = saved.copy(
        source = Some("WGS_DERIVED"),
        totalMarkers = Some(500)
      )
      repo.update(updated)

      val found = repo.findById(saved.id)
      assert(found.isDefined)
      assertEquals(found.get.source, Some("WGS_DERIVED"))
      assertEquals(found.get.totalMarkers, Some(500))
      assertEquals(found.get.meta.version, 2)
    }
  }

  testTransactor.test("delete removes entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(StrProfileEntity.create(biosample.id))

      assert(repo.exists(entity.id))
      val deleted = repo.delete(entity.id)
      assert(deleted)
      assert(!repo.exists(entity.id))
    }
  }

  testTransactor.test("markSynced updates sync status and AT fields") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(StrProfileEntity.create(biosample.id))

      repo.markSynced(entity.id, "at://did:plc:test/strprofile/1", "bafycid123")

      val found = repo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.meta.syncStatus, SyncStatus.Synced)
      assertEquals(found.get.meta.atUri, Some("at://did:plc:test/strprofile/1"))
      assertEquals(found.get.meta.atCid, Some("bafycid123"))
    }
  }

  testTransactor.test("findPendingSync returns Local and Modified entities") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val local = repo.insert(StrProfileEntity.create(biosample.id))
      val synced = repo.insert(StrProfileEntity.create(biosample.id))

      repo.markSynced(synced.id, "at://test/1", "cid1")

      val pending = repo.findPendingSync()
      assertEquals(pending.size, 1)
      assertEquals(pending.head.id, local.id)
    }
  }

  testTransactor.test("cascades delete on biosample deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = repo.insert(StrProfileEntity.create(biosample.id))

      biosampleRepo.delete(biosample.id)

      assertEquals(repo.findById(profile.id), None)
    }
  }
