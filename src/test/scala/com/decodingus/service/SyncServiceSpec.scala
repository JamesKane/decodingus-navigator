package com.decodingus.service

import com.decodingus.db.{Database, DatabaseTestSupport, Migrator, Transactor}
import com.decodingus.repository.*
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class SyncServiceSpec extends FunSuite with DatabaseTestSupport:

  private def createSyncService(tx: Transactor): H2SyncService =
    H2SyncService(
      transactor = tx,
      queueRepo = SyncQueueRepository(),
      historyRepo = SyncHistoryRepository(),
      conflictRepo = SyncConflictRepository()
    )

  testTransactor.test("enqueuePush creates queue entry") { case (db, tx) =>
    val service = createSyncService(tx)
    val entityId = UUID.randomUUID()

    val result = service.enqueuePush(
      entityType = SyncEntityType.Biosample,
      entityId = entityId,
      operation = SyncOperation.Create,
      priority = 3
    )

    assert(result.isRight)
    result.foreach { entry =>
      assertEquals(entry.entityType, SyncEntityType.Biosample)
      assertEquals(entry.entityId, entityId)
      assertEquals(entry.status, QueueStatus.Pending)
      assertEquals(entry.priority, 3)
    }
  }

  testTransactor.test("getPendingCount returns correct count") { case (db, tx) =>
    val service = createSyncService(tx)

    service.enqueuePush(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)
    service.enqueuePush(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create)

    val count = service.getPendingCount()
    assertEquals(count, Right(2L))
  }

  testTransactor.test("getNextBatch returns ordered batch") { case (db, tx) =>
    val service = createSyncService(tx)

    service.enqueuePush(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create, priority = 5)
    service.enqueuePush(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create, priority = 1)

    val batch = service.getNextBatch(10)
    assert(batch.isRight)
    batch.foreach { entries =>
      assertEquals(entries.size, 2)
      assertEquals(entries.head.priority, 1) // Higher priority first
    }
  }

  testTransactor.test("markSynced records success in history") { case (db, tx) =>
    val service = createSyncService(tx)
    val entityId = UUID.randomUUID()

    // Enqueue and start processing
    val entry = service.enqueuePush(SyncEntityType.Biosample, entityId, SyncOperation.Create).toOption.get
    service.startProcessing(entry.id)

    // Mark as synced
    val result = service.markSynced(entry.id, Some("at://test/1"), Some("cid123"))
    assert(result.isRight)

    // Check history
    val history = service.getEntityHistory(SyncEntityType.Biosample, entityId)
    assert(history.isRight)
    history.foreach { entries =>
      assertEquals(entries.size, 1)
      assertEquals(entries.head.status, SyncResultStatus.Success)
    }
  }

  testTransactor.test("markFailed records failure and schedules retry") { case (db, tx) =>
    val service = createSyncService(tx)
    val entityId = UUID.randomUUID()

    val entry = service.enqueuePush(SyncEntityType.Biosample, entityId, SyncOperation.Create).toOption.get
    service.startProcessing(entry.id)

    val result = service.markFailed(entry.id, "Connection refused")
    assert(result.isRight)

    // Check history
    val history = service.getEntityHistory(SyncEntityType.Biosample, entityId)
    history.foreach { entries =>
      assertEquals(entries.size, 1)
      assertEquals(entries.head.status, SyncResultStatus.Failed)
    }

    // Should still be pending (for retry)
    val pending = service.getPendingCount()
    assertEquals(pending, Right(1L))
  }

  testTransactor.test("recordConflict and resolution flow") { case (db, tx) =>
    val service = createSyncService(tx)
    val entityId = UUID.randomUUID()

    // Record a conflict
    val conflictResult = service.recordConflict(
      entityType = SyncEntityType.Biosample,
      entityId = entityId,
      localVersion = 3,
      remoteVersion = 5,
      atUri = Some("at://test/1")
    )
    assert(conflictResult.isRight)
    val conflict = conflictResult.toOption.get

    // Should have 1 unresolved conflict
    assertEquals(service.getUnresolvedConflictCount(), Right(1L))

    // Resolve it
    val resolveResult = service.resolveKeepLocal(conflict.id)
    assert(resolveResult.isRight)

    // Should have 0 unresolved conflicts
    assertEquals(service.getUnresolvedConflictCount(), Right(0L))
  }

  testTransactor.test("getSyncStatus returns comprehensive status") { case (db, tx) =>
    val service = createSyncService(tx)

    // Create some queue entries
    val e1 = service.enqueuePush(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create).toOption.get
    service.enqueuePush(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create)
    service.startProcessing(e1.id)

    // Create a conflict
    service.recordConflict(SyncEntityType.Alignment, UUID.randomUUID(), 1, 2)

    val status = service.getSyncStatus()
    assert(status.isRight)
    status.foreach { info =>
      assertEquals(info.pendingCount, 1L)
      assertEquals(info.inProgressCount, 1L)
      assertEquals(info.unresolvedConflicts, 1L)
      // Not healthy because of conflict
      assert(!info.isHealthy)
    }
  }

  testTransactor.test("isSyncHealthy returns false when conflicts exist") { case (db, tx) =>
    val service = createSyncService(tx)

    // Initially healthy
    assertEquals(service.isSyncHealthy(), Right(true))

    // Add a conflict
    service.recordConflict(SyncEntityType.Biosample, UUID.randomUUID(), 1, 2)

    // Now unhealthy
    assertEquals(service.isSyncHealthy(), Right(false))
  }

  testTransactor.test("cancelSync cancels pending entry") { case (db, tx) =>
    val service = createSyncService(tx)

    val entry = service.enqueuePush(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create).toOption.get

    val result = service.cancelSync(entry.id)
    assert(result.isRight)

    assertEquals(service.getPendingCount(), Right(0L))
  }

  testTransactor.test("cancelAllForEntity cancels all entries for entity") { case (db, tx) =>
    val service = createSyncService(tx)
    val entityId = UUID.randomUUID()

    service.enqueuePush(SyncEntityType.Biosample, entityId, SyncOperation.Create)
    service.enqueuePush(SyncEntityType.Biosample, entityId, SyncOperation.Update)
    service.enqueuePush(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create)

    val result = service.cancelAllForEntity(SyncEntityType.Biosample, entityId)
    assertEquals(result, Right(2))
    assertEquals(service.getPendingCount(), Right(1L))
  }

  testTransactor.test("getRecentHistory returns ordered history") { case (db, tx) =>
    val service = createSyncService(tx)
    val entityId = UUID.randomUUID()

    // Create entries and sync them
    val e1 = service.enqueuePush(SyncEntityType.Biosample, entityId, SyncOperation.Create).toOption.get
    service.startProcessing(e1.id)
    service.markSynced(e1.id)

    val e2 = service.enqueuePush(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create).toOption.get
    service.startProcessing(e2.id)
    service.markFailed(e2.id, "Error")

    val history = service.getRecentHistory(10)
    assert(history.isRight)
    history.foreach { entries =>
      assertEquals(entries.size, 2)
    }
  }

  testTransactor.test("getRecentFailures returns only failures") { case (db, tx) =>
    val service = createSyncService(tx)

    val e1 = service.enqueuePush(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create).toOption.get
    service.startProcessing(e1.id)
    service.markSynced(e1.id)

    val e2 = service.enqueuePush(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create).toOption.get
    service.startProcessing(e2.id)
    service.markPermanentlyFailed(e2.id, "Fatal error")

    val failures = service.getRecentFailures()
    assert(failures.isRight)
    failures.foreach { entries =>
      assertEquals(entries.size, 1)
      assertEquals(entries.head.status, SyncResultStatus.Failed)
    }
  }

  testTransactor.test("cleanup operations work correctly") { case (db, tx) =>
    val service = createSyncService(tx)

    // Create completed queue entry
    val e1 = service.enqueuePush(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create).toOption.get
    service.startProcessing(e1.id)
    service.markSynced(e1.id)

    // Create history
    val e2 = service.enqueuePush(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create).toOption.get
    service.startProcessing(e2.id)
    service.markSynced(e2.id)

    // Create resolved conflict
    val conflict = service.recordConflict(SyncEntityType.Alignment, UUID.randomUUID(), 1, 2).toOption.get
    service.resolveKeepLocal(conflict.id)

    // Cleanup (0 days removes all)
    val queueCleanup = service.cleanupQueue(0)
    assert(queueCleanup.isRight)

    val historyCleanup = service.cleanupHistory(0)
    assert(historyCleanup.isRight)

    val conflictCleanup = service.cleanupConflicts(0)
    assert(conflictCleanup.isRight)
  }
