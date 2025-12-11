package com.decodingus.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.*
import munit.FunSuite

import java.time.LocalDateTime
import java.util.UUID

/**
 * Integration tests for full conflict resolution workflow.
 * Covers: create → update → sync → conflict → resolve
 * Based on design doc integration test requirements.
 */
class ConflictResolutionWorkflowSpec extends FunSuite with DatabaseTestSupport:

  val biosampleRepo = BiosampleRepository()
  val syncQueueRepo = SyncQueueRepository()
  val syncHistoryRepo = SyncHistoryRepository()
  val syncConflictRepo = SyncConflictRepository()

  // Full workflow: create → update → sync → conflict detection
  testTransactor.test("full workflow: create biosample, update, sync, detect conflict") { case (db, tx) =>
    tx.readWrite {
      // Step 1: CREATE - User creates a new biosample
      val biosample = biosampleRepo.insert(BiosampleEntity.create(
        sampleAccession = "CONFLICT-TEST-001",
        donorIdentifier = "DONOR-001"
      ))
      assertEquals(biosample.meta.syncStatus, SyncStatus.Local)

      // Queue for sync
      syncQueueRepo.enqueue(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        operation = SyncOperation.Create
      )

      // Step 2: Simulate sync completion - mark as synced
      biosampleRepo.markSynced(biosample.id, "at://did:plc:user/biosample/1", "bafycid123")
      val syncedBiosample = biosampleRepo.findById(biosample.id).get
      assertEquals(syncedBiosample.meta.syncStatus, SyncStatus.Synced)

      // Record sync history
      val startTime = LocalDateTime.now()
      syncHistoryRepo.insert(SyncHistoryEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        operation = SyncOperation.Create,
        direction = SyncDirection.Push,
        status = SyncResultStatus.Success,
        startedAt = startTime,
        completedAt = LocalDateTime.now(),
        atUri = Some("at://did:plc:user/biosample/1")
      ))

      // Step 3: UPDATE - User makes a local change
      val updated = syncedBiosample.copy(
        donorIdentifier = "DONOR-001-UPDATED",
        meta = syncedBiosample.meta.copy(
          syncStatus = SyncStatus.Modified,
          version = syncedBiosample.meta.version + 1
        )
      )
      biosampleRepo.update(updated)

      // Queue update for sync
      syncQueueRepo.enqueue(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        operation = SyncOperation.Update
      )

      // Verify modified status
      val modifiedBiosample = biosampleRepo.findById(biosample.id).get
      assertEquals(modifiedBiosample.meta.syncStatus, SyncStatus.Modified)

      // Step 4: CONFLICT - Simulate remote change detected during sync
      val conflict = syncConflictRepo.insert(SyncConflictEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        localVersion = modifiedBiosample.meta.version,
        remoteVersion = modifiedBiosample.meta.version + 1
      ))

      // Update entity to conflict status
      biosampleRepo.update(modifiedBiosample.copy(
        meta = modifiedBiosample.meta.copy(syncStatus = SyncStatus.Conflict)
      ))

      // Verify conflict state
      val conflictedBiosample = biosampleRepo.findById(biosample.id).get
      assertEquals(conflictedBiosample.meta.syncStatus, SyncStatus.Conflict)

      val pendingConflicts = syncConflictRepo.findUnresolved()
      assertEquals(pendingConflicts.size, 1)
      assertEquals(pendingConflicts.head.entityId, biosample.id)
    }
  }

  // Conflict resolution: keep local
  testTransactor.test("conflict resolution: keep local changes") { case (db, tx) =>
    tx.readWrite {
      // Setup: create a conflicted biosample
      val biosample = biosampleRepo.insert(BiosampleEntity.create(
        sampleAccession = "KEEP-LOCAL-001",
        donorIdentifier = "LOCAL-VALUE"
      ))
      biosampleRepo.markSynced(biosample.id, "at://test/1", "cid1")

      // Update locally
      val updated = biosample.copy(
        donorIdentifier = "LOCAL-UPDATED",
        meta = biosample.meta.copy(
          syncStatus = SyncStatus.Conflict,
          version = 2
        )
      )
      biosampleRepo.update(updated)

      // Create conflict record
      val conflict = syncConflictRepo.insert(SyncConflictEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        localVersion = 2,
        remoteVersion = 3
      ))

      // RESOLVE: Keep local
      val resolved = syncConflictRepo.resolveKeepLocal(conflict.id, "test-user")
      assert(resolved)

      // Re-queue for sync with local version
      syncQueueRepo.enqueue(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        operation = SyncOperation.Update,
        priority = 1 // High priority for conflict resolution
      )

      // Update entity to modified status (ready for re-sync)
      val current = biosampleRepo.findById(biosample.id).get
      biosampleRepo.update(current.copy(
        meta = current.meta.copy(syncStatus = SyncStatus.Modified)
      ))

      // Verify resolution
      val resolvedConflict = syncConflictRepo.findById(conflict.id).get
      assertEquals(resolvedConflict.status, ConflictStatus.Resolved)
      assertEquals(resolvedConflict.resolutionAction, Some(ResolutionAction.KeptLocal))

      val finalBiosample = biosampleRepo.findById(biosample.id).get
      assertEquals(finalBiosample.donorIdentifier, "LOCAL-UPDATED")
      assertEquals(finalBiosample.meta.syncStatus, SyncStatus.Modified)
    }
  }

  // Conflict resolution: accept remote
  testTransactor.test("conflict resolution: accept remote changes") { case (db, tx) =>
    tx.readWrite {
      // Setup: create a conflicted biosample
      val biosample = biosampleRepo.insert(BiosampleEntity.create(
        sampleAccession = "ACCEPT-REMOTE-001",
        donorIdentifier = "LOCAL-VALUE"
      ))
      biosampleRepo.markSynced(biosample.id, "at://test/2", "cid2")

      // Update locally
      val updated = biosample.copy(
        donorIdentifier = "LOCAL-UPDATED",
        meta = biosample.meta.copy(
          syncStatus = SyncStatus.Conflict,
          version = 2
        )
      )
      biosampleRepo.update(updated)

      // Create conflict record
      val conflict = syncConflictRepo.insert(SyncConflictEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        localVersion = 2,
        remoteVersion = 3
      ))

      // RESOLVE: Accept remote
      val resolved = syncConflictRepo.resolveAcceptRemote(conflict.id, "test-user")
      assert(resolved)

      // Apply remote value (simulating what sync service would do)
      val applying = biosampleRepo.findById(biosample.id).get
      biosampleRepo.update(applying.copy(
        donorIdentifier = "REMOTE-VALUE",
        meta = applying.meta.copy(version = 3)
      ))
      // Mark as synced to properly set sync status
      biosampleRepo.markSynced(biosample.id, "at://test/2", "cid2-new")

      // Verify resolution
      val resolvedConflict = syncConflictRepo.findById(conflict.id).get
      assertEquals(resolvedConflict.status, ConflictStatus.Resolved)
      assertEquals(resolvedConflict.resolutionAction, Some(ResolutionAction.AcceptedRemote))

      val finalBiosample = biosampleRepo.findById(biosample.id).get
      assertEquals(finalBiosample.donorIdentifier, "REMOTE-VALUE")
      assertEquals(finalBiosample.meta.syncStatus, SyncStatus.Synced)
    }
  }

  // Multiple conflicts resolution batch
  testTransactor.test("batch conflict resolution across multiple entities") { case (db, tx) =>
    tx.readWrite {
      // Create multiple conflicted biosamples
      val biosamples = (1 to 5).map { i =>
        val b = biosampleRepo.insert(BiosampleEntity.create(
          sampleAccession = s"BATCH-CONFLICT-$i",
          donorIdentifier = s"LOCAL-$i"
        ))
        biosampleRepo.markSynced(b.id, s"at://test/$i", s"cid$i")
        biosampleRepo.update(b.copy(
          meta = b.meta.copy(syncStatus = SyncStatus.Conflict, version = 2)
        ))
        b
      }

      // Create conflict records for each
      val conflicts = biosamples.map { b =>
        syncConflictRepo.insert(SyncConflictEntity.create(
          entityType = SyncEntityType.Biosample,
          entityId = b.id,
          localVersion = 2,
          remoteVersion = 3
        ))
      }

      // Verify pending conflicts
      assertEquals(syncConflictRepo.findUnresolved().size, 5)

      // Batch resolve all with KeepLocal
      conflicts.foreach { c =>
        syncConflictRepo.resolveKeepLocal(c.id, "batch-user")
      }

      // Verify all resolved
      assertEquals(syncConflictRepo.findUnresolved().size, 0)
      assertEquals(syncConflictRepo.findAll().count(_.status == ConflictStatus.Resolved), 5)
    }
  }

  // Conflict resolution with merge
  testTransactor.test("conflict resolution: merge non-overlapping changes") { case (db, tx) =>
    tx.readWrite {
      // Setup: biosample with description
      val biosample = biosampleRepo.insert(BiosampleEntity.create(
        sampleAccession = "AUTO-MERGE-001",
        donorIdentifier = "DONOR-001",
        description = Some("Original description")
      ))
      biosampleRepo.markSynced(biosample.id, "at://test/merge", "cid-merge")

      // Local change: update description
      val localUpdate = biosample.copy(
        description = Some("Updated local description"),
        meta = biosample.meta.copy(
          syncStatus = SyncStatus.Modified,
          version = 2
        )
      )
      biosampleRepo.update(localUpdate)

      // Create conflict record - non-overlapping changes
      val conflict = syncConflictRepo.insert(SyncConflictEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        localVersion = 2,
        remoteVersion = 2,
        suggestedResolution = Some(ConflictResolution.Merge)
      ))

      // Resolve with merge
      val resolved = syncConflictRepo.resolveMerge(conflict.id, "merge-user")
      assert(resolved)

      // Verify resolution
      val resolvedConflict = syncConflictRepo.findById(conflict.id).get
      assertEquals(resolvedConflict.status, ConflictStatus.Resolved)
      assertEquals(resolvedConflict.resolutionAction, Some(ResolutionAction.Merged))
    }
  }

  // Sync history tracking through full workflow
  testTransactor.test("sync history tracks full workflow lifecycle") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create(
        sampleAccession = "HISTORY-TEST-001",
        donorIdentifier = "DONOR-001"
      ))

      val now = LocalDateTime.now()

      // Record create push
      syncHistoryRepo.insert(SyncHistoryEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        operation = SyncOperation.Create,
        direction = SyncDirection.Push,
        status = SyncResultStatus.Success,
        startedAt = now,
        completedAt = now,
        atUri = Some("at://test/history")
      ))

      // Record update push that failed
      syncHistoryRepo.insert(SyncHistoryEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        operation = SyncOperation.Update,
        direction = SyncDirection.Push,
        status = SyncResultStatus.Failed,
        startedAt = now,
        completedAt = now,
        errorMessage = Some("Conflict detected"),
        atUri = Some("at://test/history")
      ))

      // Record conflict resolution
      syncHistoryRepo.insert(SyncHistoryEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        operation = SyncOperation.Update,
        direction = SyncDirection.Push,
        status = SyncResultStatus.Success,
        startedAt = now,
        completedAt = now,
        atUri = Some("at://test/history")
      ))

      // Verify history
      val history = syncHistoryRepo.findByEntity(SyncEntityType.Biosample, biosample.id)
      assertEquals(history.size, 3)

      val successCount = history.count(_.status == SyncResultStatus.Success)
      val failedCount = history.count(_.status == SyncResultStatus.Failed)
      assertEquals(successCount, 2)
      assertEquals(failedCount, 1)
    }
  }

  // Conflict dismissal
  testTransactor.test("conflicts can be dismissed") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create(
        sampleAccession = "DISMISS-001",
        donorIdentifier = "DONOR-001"
      ))

      // Create conflict
      val conflict = syncConflictRepo.insert(SyncConflictEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        localVersion = 1,
        remoteVersion = 2
      ))

      assertEquals(syncConflictRepo.findUnresolved().size, 1)

      // Dismiss the conflict
      val dismissed = syncConflictRepo.dismiss(conflict.id, "dismiss-user")
      assert(dismissed)

      // Verify dismissed
      val dismissedConflict = syncConflictRepo.findById(conflict.id).get
      assertEquals(dismissedConflict.status, ConflictStatus.Dismissed)

      // Should not appear in unresolved
      assertEquals(syncConflictRepo.findUnresolved().size, 0)
    }
  }

  // Active conflict check
  testTransactor.test("findActiveConflict finds unresolved conflict for entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create(
        sampleAccession = "ACTIVE-001",
        donorIdentifier = "DONOR-001"
      ))

      // Initially no active conflict
      val noConflict = syncConflictRepo.findActiveConflict(SyncEntityType.Biosample, biosample.id)
      assert(noConflict.isEmpty)

      // Create conflict
      syncConflictRepo.insert(SyncConflictEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = biosample.id,
        localVersion = 1,
        remoteVersion = 2
      ))

      // Now should find active conflict
      val activeConflict = syncConflictRepo.findActiveConflict(SyncEntityType.Biosample, biosample.id)
      assert(activeConflict.isDefined)
      assertEquals(activeConflict.get.entityId, biosample.id)
    }
  }
