package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import munit.FunSuite
import java.util.UUID

class SyncQueueRepositorySpec extends FunSuite with DatabaseTestSupport:

  val queueRepo = SyncQueueRepository()

  testTransactor.test("enqueue creates new queue entry") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()
      val entry = queueRepo.enqueue(
        entityType = SyncEntityType.Biosample,
        entityId = entityId,
        operation = SyncOperation.Create,
        priority = 3
      )

      assertEquals(entry.entityType, SyncEntityType.Biosample)
      assertEquals(entry.entityId, entityId)
      assertEquals(entry.operation, SyncOperation.Create)
      assertEquals(entry.status, QueueStatus.Pending)
      assertEquals(entry.priority, 3)
      assertEquals(entry.attemptCount, 0)
    }
  }

  testTransactor.test("enqueue updates existing pending entry") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()

      // First enqueue
      val first = queueRepo.enqueue(SyncEntityType.Biosample, entityId, SyncOperation.Create, priority = 5)

      // Second enqueue with higher priority should update
      val second = queueRepo.enqueue(SyncEntityType.Biosample, entityId, SyncOperation.Create, priority = 1)

      // Should still have only one entry
      val all = queueRepo.findByEntity(SyncEntityType.Biosample, entityId)
      assertEquals(all.size, 1)

      // Priority should be the minimum (1)
      assertEquals(all.head.priority, 1)
    }
  }

  testTransactor.test("findPendingBatch returns ordered batch") { case (db, tx) =>
    tx.readWrite {
      // Create entries with different priorities
      queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create, priority = 5)
      queueRepo.enqueue(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create, priority = 1)
      queueRepo.enqueue(SyncEntityType.Alignment, UUID.randomUUID(), SyncOperation.Create, priority = 3)

      val batch = queueRepo.findPendingBatch(10)
      assertEquals(batch.size, 3)

      // Should be ordered by priority (ascending)
      assertEquals(batch.head.priority, 1)
      assertEquals(batch(1).priority, 3)
      assertEquals(batch(2).priority, 5)
    }
  }

  testTransactor.test("markInProgress updates status and increments attempt count") { case (db, tx) =>
    tx.readWrite {
      val entry = queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)

      val marked = queueRepo.markInProgress(entry.id)
      assert(marked)

      val found = queueRepo.findById(entry.id)
      assert(found.isDefined)
      assertEquals(found.get.status, QueueStatus.InProgress)
      assertEquals(found.get.attemptCount, 1)
      assert(found.get.startedAt.isDefined)
    }
  }

  testTransactor.test("markCompleted updates status") { case (db, tx) =>
    tx.readWrite {
      val entry = queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)
      queueRepo.markInProgress(entry.id)

      val marked = queueRepo.markCompleted(entry.id)
      assert(marked)

      val found = queueRepo.findById(entry.id)
      assert(found.isDefined)
      assertEquals(found.get.status, QueueStatus.Completed)
      assert(found.get.completedAt.isDefined)
    }
  }

  testTransactor.test("markFailedWithRetry sets next retry time with backoff") { case (db, tx) =>
    tx.readWrite {
      val entry = queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)
      queueRepo.markInProgress(entry.id)

      queueRepo.markFailedWithRetry(entry.id, "Connection timeout")

      val found = queueRepo.findById(entry.id)
      assert(found.isDefined)
      assertEquals(found.get.status, QueueStatus.Pending) // Back to pending
      assertEquals(found.get.lastError, Some("Connection timeout"))
      assert(found.get.nextRetryAt.isDefined)
    }
  }

  testTransactor.test("markFailed sets permanent failure") { case (db, tx) =>
    tx.readWrite {
      val entry = queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)
      queueRepo.markInProgress(entry.id)

      queueRepo.markFailed(entry.id, "Entity not found")

      val found = queueRepo.findById(entry.id)
      assert(found.isDefined)
      assertEquals(found.get.status, QueueStatus.Failed)
      assertEquals(found.get.lastError, Some("Entity not found"))
      assert(found.get.completedAt.isDefined)
    }
  }

  testTransactor.test("cancel changes status to Cancelled") { case (db, tx) =>
    tx.readWrite {
      val entry = queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)

      val cancelled = queueRepo.cancel(entry.id)
      assert(cancelled)

      val found = queueRepo.findById(entry.id)
      assertEquals(found.get.status, QueueStatus.Cancelled)
    }
  }

  testTransactor.test("cancelByEntity cancels all pending for entity") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()

      queueRepo.enqueue(SyncEntityType.Biosample, entityId, SyncOperation.Create)
      queueRepo.enqueue(SyncEntityType.Biosample, entityId, SyncOperation.Update)

      // Different entity should not be affected
      queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)

      val cancelled = queueRepo.cancelByEntity(SyncEntityType.Biosample, entityId)
      assertEquals(cancelled, 2)

      val remaining = queueRepo.findByStatus(QueueStatus.Pending)
      assertEquals(remaining.size, 1)
    }
  }

  testTransactor.test("findByEntity returns all entries for entity") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()

      queueRepo.enqueue(SyncEntityType.Biosample, entityId, SyncOperation.Create)
      queueRepo.enqueue(SyncEntityType.Biosample, entityId, SyncOperation.Update)
      queueRepo.enqueue(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create)

      val entries = queueRepo.findByEntity(SyncEntityType.Biosample, entityId)
      assertEquals(entries.size, 2)
    }
  }

  testTransactor.test("findByStatus returns entries with matching status") { case (db, tx) =>
    tx.readWrite {
      val e1 = queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)
      val e2 = queueRepo.enqueue(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create)

      queueRepo.markInProgress(e1.id)

      val pending = queueRepo.findByStatus(QueueStatus.Pending)
      assertEquals(pending.size, 1)

      val inProgress = queueRepo.findByStatus(QueueStatus.InProgress)
      assertEquals(inProgress.size, 1)
    }
  }

  testTransactor.test("countByStatus returns correct counts") { case (db, tx) =>
    tx.readWrite {
      val e1 = queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)
      val e2 = queueRepo.enqueue(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create)
      val e3 = queueRepo.enqueue(SyncEntityType.Alignment, UUID.randomUUID(), SyncOperation.Create)

      queueRepo.markInProgress(e1.id)
      queueRepo.markCompleted(e1.id)
      queueRepo.markInProgress(e2.id)

      val counts = queueRepo.countByStatus()
      assertEquals(counts.getOrElse(QueueStatus.Pending, 0L), 1L)
      assertEquals(counts.getOrElse(QueueStatus.InProgress, 0L), 1L)
      assertEquals(counts.getOrElse(QueueStatus.Completed, 0L), 1L)
    }
  }

  testTransactor.test("countPending returns pending count") { case (db, tx) =>
    tx.readWrite {
      queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)
      queueRepo.enqueue(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create)
      val e3 = queueRepo.enqueue(SyncEntityType.Alignment, UUID.randomUUID(), SyncOperation.Create)

      queueRepo.markInProgress(e3.id)

      assertEquals(queueRepo.countPending(), 2L)
    }
  }

  testTransactor.test("cleanupCompleted removes old completed entries") { case (db, tx) =>
    tx.readWrite {
      val e1 = queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)
      val e2 = queueRepo.enqueue(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create)

      queueRepo.markInProgress(e1.id)
      queueRepo.markCompleted(e1.id)
      queueRepo.markInProgress(e2.id)
      queueRepo.markCompleted(e2.id)

      // Cleanup with 0 days should remove completed entries
      val removed = queueRepo.cleanupCompleted(0)
      assertEquals(removed, 2)
    }
  }

  testTransactor.test("delete removes queue entry") { case (db, tx) =>
    tx.readWrite {
      val entry = queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)

      assert(queueRepo.findById(entry.id).isDefined)

      val deleted = queueRepo.delete(entry.id)
      assert(deleted)
      assertEquals(queueRepo.findById(entry.id), None)
    }
  }
