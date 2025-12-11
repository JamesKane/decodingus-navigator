package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import io.circe.Json
import io.circe.syntax.*
import munit.FunSuite
import java.util.UUID

class SyncConflictRepositorySpec extends FunSuite with DatabaseTestSupport:

  val conflictRepo = SyncConflictRepository()

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()

      val entity = SyncConflictEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = entityId,
        localVersion = 3,
        remoteVersion = 5,
        atUri = Some("at://test/biosample/1"),
        suggestedResolution = Some(ConflictResolution.AcceptRemote),
        resolutionReason = Some("Remote version is newer")
      )

      val saved = conflictRepo.insert(entity)
      val found = conflictRepo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.entityType, SyncEntityType.Biosample)
      assertEquals(found.get.entityId, entityId)
      assertEquals(found.get.localVersion, 3)
      assertEquals(found.get.remoteVersion, 5)
      assertEquals(found.get.status, ConflictStatus.Unresolved)
      assertEquals(found.get.suggestedResolution, Some(ConflictResolution.AcceptRemote))
    }
  }

  testTransactor.test("findUnresolved returns only unresolved conflicts") { case (db, tx) =>
    tx.readWrite {
      val c1 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))
      val c2 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Project, UUID.randomUUID(), 1, 2))
      val c3 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Alignment, UUID.randomUUID(), 1, 2))

      conflictRepo.resolveKeepLocal(c1.id, "user")

      val unresolved = conflictRepo.findUnresolved()
      assertEquals(unresolved.size, 2)
      assert(unresolved.forall(_.status == ConflictStatus.Unresolved))
    }
  }

  testTransactor.test("findByEntity returns conflicts for entity") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()

      conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, entityId, 1, 2))
      conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, entityId, 2, 3))
      conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))

      val conflicts = conflictRepo.findByEntity(SyncEntityType.Biosample, entityId)
      assertEquals(conflicts.size, 2)
    }
  }

  testTransactor.test("findActiveConflict returns unresolved conflict for entity") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()

      val c1 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, entityId, 1, 2))
      conflictRepo.resolveKeepLocal(c1.id, "user")

      val c2 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, entityId, 2, 4))

      val active = conflictRepo.findActiveConflict(SyncEntityType.Biosample, entityId)
      assert(active.isDefined)
      assertEquals(active.get.id, c2.id)
    }
  }

  testTransactor.test("resolveKeepLocal resolves with KEPT_LOCAL action") { case (db, tx) =>
    tx.readWrite {
      val conflict = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))

      val resolved = conflictRepo.resolveKeepLocal(conflict.id, "test_user")
      assert(resolved)

      val found = conflictRepo.findById(conflict.id)
      assert(found.isDefined)
      assertEquals(found.get.status, ConflictStatus.Resolved)
      assertEquals(found.get.resolutionAction, Some(ResolutionAction.KeptLocal))
      assertEquals(found.get.resolvedBy, Some("test_user"))
      assert(found.get.resolvedAt.isDefined)
    }
  }

  testTransactor.test("resolveAcceptRemote resolves with ACCEPTED_REMOTE action") { case (db, tx) =>
    tx.readWrite {
      val conflict = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))

      val resolved = conflictRepo.resolveAcceptRemote(conflict.id, "test_user")
      assert(resolved)

      val found = conflictRepo.findById(conflict.id)
      assertEquals(found.get.status, ConflictStatus.Resolved)
      assertEquals(found.get.resolutionAction, Some(ResolutionAction.AcceptedRemote))
    }
  }

  testTransactor.test("resolveMerge resolves with MERGED action") { case (db, tx) =>
    tx.readWrite {
      val conflict = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))

      val resolved = conflictRepo.resolveMerge(conflict.id, "test_user")
      assert(resolved)

      val found = conflictRepo.findById(conflict.id)
      assertEquals(found.get.status, ConflictStatus.Resolved)
      assertEquals(found.get.resolutionAction, Some(ResolutionAction.Merged))
    }
  }

  testTransactor.test("resolveManual resolves with MANUAL_EDIT action") { case (db, tx) =>
    tx.readWrite {
      val conflict = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))

      val resolved = conflictRepo.resolveManual(conflict.id, "test_user")
      assert(resolved)

      val found = conflictRepo.findById(conflict.id)
      assertEquals(found.get.status, ConflictStatus.Resolved)
      assertEquals(found.get.resolutionAction, Some(ResolutionAction.ManualEdit))
    }
  }

  testTransactor.test("dismiss sets status to DISMISSED") { case (db, tx) =>
    tx.readWrite {
      val conflict = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))

      val dismissed = conflictRepo.dismiss(conflict.id, "test_user")
      assert(dismissed)

      val found = conflictRepo.findById(conflict.id)
      assertEquals(found.get.status, ConflictStatus.Dismissed)
      assertEquals(found.get.resolvedBy, Some("test_user"))
    }
  }

  testTransactor.test("resolve operations only affect unresolved conflicts") { case (db, tx) =>
    tx.readWrite {
      val conflict = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))

      // Resolve once
      conflictRepo.resolveKeepLocal(conflict.id, "user1")

      // Try to resolve again - should not change
      val resolved = conflictRepo.resolveAcceptRemote(conflict.id, "user2")
      assert(!resolved)

      val found = conflictRepo.findById(conflict.id)
      // Should still be KEPT_LOCAL from first resolution
      assertEquals(found.get.resolutionAction, Some(ResolutionAction.KeptLocal))
      assertEquals(found.get.resolvedBy, Some("user1"))
    }
  }

  testTransactor.test("countUnresolved returns correct count") { case (db, tx) =>
    tx.readWrite {
      val c1 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))
      conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Project, UUID.randomUUID(), 1, 2))
      conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Alignment, UUID.randomUUID(), 1, 2))

      conflictRepo.resolveKeepLocal(c1.id, "user")

      assertEquals(conflictRepo.countUnresolved(), 2L)
    }
  }

  testTransactor.test("countByStatus returns correct counts") { case (db, tx) =>
    tx.readWrite {
      val c1 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))
      val c2 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Project, UUID.randomUUID(), 1, 2))
      val c3 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Alignment, UUID.randomUUID(), 1, 2))

      conflictRepo.resolveKeepLocal(c1.id, "user")
      conflictRepo.dismiss(c2.id, "user")

      val counts = conflictRepo.countByStatus()
      assertEquals(counts.getOrElse(ConflictStatus.Unresolved, 0L), 1L)
      assertEquals(counts.getOrElse(ConflictStatus.Resolved, 0L), 1L)
      assertEquals(counts.getOrElse(ConflictStatus.Dismissed, 0L), 1L)
    }
  }

  testTransactor.test("findByStatus returns conflicts with matching status") { case (db, tx) =>
    tx.readWrite {
      val c1 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))
      val c2 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Project, UUID.randomUUID(), 1, 2))

      conflictRepo.resolveKeepLocal(c1.id, "user")

      val resolved = conflictRepo.findByStatus(ConflictStatus.Resolved)
      assertEquals(resolved.size, 1)
      assertEquals(resolved.head.id, c1.id)

      val unresolved = conflictRepo.findByStatus(ConflictStatus.Unresolved)
      assertEquals(unresolved.size, 1)
      assertEquals(unresolved.head.id, c2.id)
    }
  }

  testTransactor.test("cleanupResolved removes old resolved conflicts") { case (db, tx) =>
    tx.readWrite {
      val c1 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))
      val c2 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Project, UUID.randomUUID(), 1, 2))
      val c3 = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Alignment, UUID.randomUUID(), 1, 2))

      conflictRepo.resolveKeepLocal(c1.id, "user")
      conflictRepo.dismiss(c2.id, "user")
      // c3 stays unresolved

      // Cleanup with 0 days should remove resolved and dismissed
      val removed = conflictRepo.cleanupResolved(0)
      assertEquals(removed, 2)

      // Only unresolved should remain
      assertEquals(conflictRepo.findAll().size, 1)
      assertEquals(conflictRepo.findUnresolved().size, 1)
    }
  }

  testTransactor.test("conflict with JSON snapshots stores and retrieves correctly") { case (db, tx) =>
    tx.readWrite {
      val localSnapshot = Json.obj("name" -> Json.fromString("Local Name"), "version" -> Json.fromInt(3))
      val remoteSnapshot = Json.obj("name" -> Json.fromString("Remote Name"), "version" -> Json.fromInt(5))
      val localChanges = Json.arr(Json.fromString("name"))
      val remoteChanges = Json.arr(Json.fromString("name"), Json.fromString("description"))
      val overlapping = Json.arr(Json.fromString("name"))

      val entity = SyncConflictEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = UUID.randomUUID(),
        localVersion = 3,
        remoteVersion = 5,
        localChanges = Some(localChanges),
        remoteChanges = Some(remoteChanges),
        overlappingFields = Some(overlapping),
        localSnapshot = Some(localSnapshot),
        remoteSnapshot = Some(remoteSnapshot)
      )

      val saved = conflictRepo.insert(entity)
      val found = conflictRepo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.localSnapshot, Some(localSnapshot))
      assertEquals(found.get.remoteSnapshot, Some(remoteSnapshot))
      assertEquals(found.get.localChanges, Some(localChanges))
      assertEquals(found.get.remoteChanges, Some(remoteChanges))
      assertEquals(found.get.overlappingFields, Some(overlapping))
    }
  }

  testTransactor.test("delete removes conflict") { case (db, tx) =>
    tx.readWrite {
      val conflict = conflictRepo.insert(SyncConflictEntity.create(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2))

      assert(conflictRepo.findById(conflict.id).isDefined)

      val deleted = conflictRepo.delete(conflict.id)
      assert(deleted)
      assertEquals(conflictRepo.findById(conflict.id), None)
    }
  }
