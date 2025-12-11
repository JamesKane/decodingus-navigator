package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import munit.FunSuite
import java.util.UUID

class BiosampleRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = BiosampleRepository()

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val entity = BiosampleEntity.create(
        sampleAccession = "TEST001",
        donorIdentifier = "DONOR001",
        description = Some("Test biosample"),
        centerName = Some("Test Center"),
        sex = Some("Male")
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.sampleAccession, "TEST001")
      assertEquals(found.get.donorIdentifier, "DONOR001")
      assertEquals(found.get.description, Some("Test biosample"))
      assertEquals(found.get.centerName, Some("Test Center"))
      assertEquals(found.get.sex, Some("Male"))
      assertEquals(found.get.meta.syncStatus, SyncStatus.Local)
      assertEquals(found.get.meta.version, 1)
    }
  }

  testTransactor.test("findAll returns all entities") { case (db, tx) =>
    tx.readWrite {
      repo.insert(BiosampleEntity.create("ACC001", "DONOR001"))
      repo.insert(BiosampleEntity.create("ACC002", "DONOR001"))
      repo.insert(BiosampleEntity.create("ACC003", "DONOR002"))

      val all = repo.findAll()
      assertEquals(all.size, 3)
    }
  }

  testTransactor.test("update modifies entity correctly") { case (db, tx) =>
    tx.readWrite {
      val entity = BiosampleEntity.create(
        sampleAccession = "UPDATE001",
        donorIdentifier = "DONOR001"
      )
      val saved = repo.insert(entity)

      // Update
      val updated = saved.copy(
        description = Some("Updated description"),
        sex = Some("Female")
      )
      repo.update(updated)

      // Verify
      val found = repo.findById(saved.id)
      assert(found.isDefined)
      assertEquals(found.get.description, Some("Updated description"))
      assertEquals(found.get.sex, Some("Female"))
      // Version should be incremented
      assertEquals(found.get.meta.version, 2)
    }
  }

  testTransactor.test("delete removes entity") { case (db, tx) =>
    tx.readWrite {
      val entity = repo.insert(BiosampleEntity.create("DELETE001", "DONOR001"))

      assert(repo.exists(entity.id), "Entity should exist before delete")

      val deleted = repo.delete(entity.id)
      assert(deleted, "Delete should return true")

      assert(!repo.exists(entity.id), "Entity should not exist after delete")
      assertEquals(repo.findById(entity.id), None)
    }
  }

  testTransactor.test("count returns correct number") { case (db, tx) =>
    tx.readWrite {
      assertEquals(repo.count(), 0L)

      repo.insert(BiosampleEntity.create("COUNT001", "DONOR001"))
      assertEquals(repo.count(), 1L)

      repo.insert(BiosampleEntity.create("COUNT002", "DONOR001"))
      assertEquals(repo.count(), 2L)
    }
  }

  testTransactor.test("findByAccession finds by unique accession") { case (db, tx) =>
    tx.readWrite {
      repo.insert(BiosampleEntity.create("UNIQUE001", "DONOR001"))
      repo.insert(BiosampleEntity.create("UNIQUE002", "DONOR002"))

      val found = repo.findByAccession("UNIQUE001")
      assert(found.isDefined)
      assertEquals(found.get.donorIdentifier, "DONOR001")

      val notFound = repo.findByAccession("NONEXISTENT")
      assertEquals(notFound, None)
    }
  }

  testTransactor.test("findByDonor returns all samples for donor") { case (db, tx) =>
    tx.readWrite {
      repo.insert(BiosampleEntity.create("DONOR_A_1", "DONOR_A"))
      repo.insert(BiosampleEntity.create("DONOR_A_2", "DONOR_A"))
      repo.insert(BiosampleEntity.create("DONOR_B_1", "DONOR_B"))

      val donorASamples = repo.findByDonor("DONOR_A")
      assertEquals(donorASamples.size, 2)
      assert(donorASamples.forall(_.donorIdentifier == "DONOR_A"))

      val donorBSamples = repo.findByDonor("DONOR_B")
      assertEquals(donorBSamples.size, 1)
    }
  }

  testTransactor.test("findByStatus filters by sync status") { case (db, tx) =>
    tx.readWrite {
      val local1 = repo.insert(BiosampleEntity.create("LOCAL1", "D1"))
      val local2 = repo.insert(BiosampleEntity.create("LOCAL2", "D2"))

      // Mark one as synced
      repo.markSynced(local1.id, "at://did:plc:test/biosample/1", "cid123")

      val localOnly = repo.findByStatus(SyncStatus.Local)
      assertEquals(localOnly.size, 1)
      assertEquals(localOnly.head.sampleAccession, "LOCAL2")

      val synced = repo.findByStatus(SyncStatus.Synced)
      assertEquals(synced.size, 1)
      assertEquals(synced.head.sampleAccession, "LOCAL1")
    }
  }

  testTransactor.test("markSynced updates status and AT Protocol fields") { case (db, tx) =>
    tx.readWrite {
      val entity = repo.insert(BiosampleEntity.create("SYNC001", "DONOR001"))
      assertEquals(entity.meta.syncStatus, SyncStatus.Local)

      repo.markSynced(entity.id, "at://did:plc:test/biosample/1", "bafycid123")

      val updated = repo.findById(entity.id)
      assert(updated.isDefined)
      assertEquals(updated.get.meta.syncStatus, SyncStatus.Synced)
      assertEquals(updated.get.meta.atUri, Some("at://did:plc:test/biosample/1"))
      assertEquals(updated.get.meta.atCid, Some("bafycid123"))
    }
  }

  testTransactor.test("updateStatus changes only sync status") { case (db, tx) =>
    tx.readWrite {
      val entity = repo.insert(BiosampleEntity.create("STATUS001", "DONOR001"))

      repo.updateStatus(entity.id, SyncStatus.Modified)

      val updated = repo.findById(entity.id)
      assert(updated.isDefined)
      assertEquals(updated.get.meta.syncStatus, SyncStatus.Modified)
    }
  }

  testTransactor.test("findPendingSync returns Local and Modified entities") { case (db, tx) =>
    tx.readWrite {
      val local = repo.insert(BiosampleEntity.create("PENDING1", "D1"))
      val synced = repo.insert(BiosampleEntity.create("PENDING2", "D2"))
      val modified = repo.insert(BiosampleEntity.create("PENDING3", "D3"))

      // Mark one as synced, then modify it
      repo.markSynced(synced.id, "at://test/1", "cid1")
      repo.markSynced(modified.id, "at://test/2", "cid2")
      repo.updateStatus(modified.id, SyncStatus.Modified)

      val pending = repo.findPendingSync()
      assertEquals(pending.size, 2)
      val accessions = pending.map(_.sampleAccession).toSet
      assert(accessions.contains("PENDING1"))
      assert(accessions.contains("PENDING3"))
      assert(!accessions.contains("PENDING2"))
    }
  }

  testTransactor.test("searchByAccession finds by prefix") { case (db, tx) =>
    tx.readWrite {
      repo.insert(BiosampleEntity.create("SEARCH_ABC_001", "D1"))
      repo.insert(BiosampleEntity.create("SEARCH_ABC_002", "D2"))
      repo.insert(BiosampleEntity.create("SEARCH_XYZ_001", "D3"))

      val abcResults = repo.searchByAccession("SEARCH_ABC")
      assertEquals(abcResults.size, 2)

      val xyzResults = repo.searchByAccession("SEARCH_XYZ")
      assertEquals(xyzResults.size, 1)

      val noResults = repo.searchByAccession("SEARCH_NONE")
      assertEquals(noResults.size, 0)
    }
  }

  testTransactor.test("unique constraint on sample_accession") { case (db, tx) =>
    val result = tx.readWrite {
      repo.insert(BiosampleEntity.create("DUPE001", "D1"))
      repo.insert(BiosampleEntity.create("DUPE001", "D2")) // Duplicate accession
    }

    assert(result.isLeft, "Should fail on duplicate accession")
  }
