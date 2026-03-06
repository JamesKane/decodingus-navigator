package com.decodingus.sync

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.*
import munit.FunSuite
import java.util.UUID
import scala.concurrent.{Await, ExecutionContext}
import scala.concurrent.duration.*

class AsyncSyncServiceSpec extends FunSuite with DatabaseTestSupport:

  given ExecutionContext = ExecutionContext.global

  testTransactor.test("queueForSync creates queue entry") { case (db, tx) =>
    val service = createService(tx)

    val entityId = UUID.randomUUID()
    val future = service.queueForSync(
      SyncEntityType.Biosample,
      entityId,
      SyncOperation.Create,
      priority = 3
    )

    val entry = Await.result(future, 5.seconds)
    assertEquals(entry.entityType, SyncEntityType.Biosample)
    assertEquals(entry.entityId, entityId)
    assertEquals(entry.operation, SyncOperation.Create)
    assertEquals(entry.priority, 3)
    assertEquals(entry.status, QueueStatus.Pending)

    // Clean up
    service.shutdown()
  }

  testTransactor.test("queueBiosampleSync creates biosample queue entry") { case (db, tx) =>
    val service = createService(tx)

    val biosampleId = UUID.randomUUID()
    val future = service.queueBiosampleSync(biosampleId, SyncOperation.Update)

    val entry = Await.result(future, 5.seconds)
    assertEquals(entry.entityType, SyncEntityType.Biosample)
    assertEquals(entry.entityId, biosampleId)
    assertEquals(entry.operation, SyncOperation.Update)

    service.shutdown()
  }

  testTransactor.test("queueProjectSync creates project queue entry") { case (db, tx) =>
    val service = createService(tx)

    val projectId = UUID.randomUUID()
    val future = service.queueProjectSync(projectId, SyncOperation.Delete)

    val entry = Await.result(future, 5.seconds)
    assertEquals(entry.entityType, SyncEntityType.Project)
    assertEquals(entry.entityId, projectId)
    assertEquals(entry.operation, SyncOperation.Delete)

    service.shutdown()
  }

  testTransactor.test("queueSequenceRunSync creates sequence run queue entry") { case (db, tx) =>
    val service = createService(tx)

    val seqRunId = UUID.randomUUID()
    val future = service.queueSequenceRunSync(seqRunId, SyncOperation.Create)

    val entry = Await.result(future, 5.seconds)
    assertEquals(entry.entityType, SyncEntityType.SequenceRun)
    assertEquals(entry.entityId, seqRunId)

    service.shutdown()
  }

  testTransactor.test("queueAlignmentSync creates alignment queue entry") { case (db, tx) =>
    val service = createService(tx)

    val alignmentId = UUID.randomUUID()
    val future = service.queueAlignmentSync(alignmentId, SyncOperation.Update)

    val entry = Await.result(future, 5.seconds)
    assertEquals(entry.entityType, SyncEntityType.Alignment)
    assertEquals(entry.entityId, alignmentId)

    service.shutdown()
  }

  testTransactor.test("getQueueStats returns correct counts") { case (db, tx) =>
    val service = createService(tx)
    val queueRepo = SyncQueueRepository()

    // Add some queue entries directly
    tx.readWrite {
      queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)
      queueRepo.enqueue(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create)
      val e3 = queueRepo.enqueue(SyncEntityType.Alignment, UUID.randomUUID(), SyncOperation.Create)
      queueRepo.markInProgress(e3.id)
      queueRepo.markCompleted(e3.id)
    }

    val stats = Await.result(service.getQueueStats, 5.seconds)
    assertEquals(stats.pending, 2L)
    assertEquals(stats.completed, 1L)
    assertEquals(stats.total, 3L)

    service.shutdown()
  }

  testTransactor.test("processQueueNow returns 0 when no user logged in") { case (db, tx) =>
    val service = createService(tx)
    val queueRepo = SyncQueueRepository()

    // Add queue entry
    tx.readWrite {
      queueRepo.enqueue(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create)
    }

    // No user set, so processing should skip
    val processed = Await.result(service.processQueueNow(), 5.seconds)
    assertEquals(processed, 0)

    service.shutdown()
  }

  testTransactor.test("setIncomingSyncEnabled updates flag") { case (db, tx) =>
    val service = createService(tx)

    // Just verify it doesn't throw
    service.setIncomingSyncEnabled(false)
    service.setIncomingSyncEnabled(true)

    service.shutdown()
  }

  testTransactor.test("setUser updates current user") { case (db, tx) =>
    val service = createService(tx)

    // Just verify it doesn't throw
    service.setUser(None)

    service.shutdown()
  }

  testTransactor.test("start and shutdown complete without error") { case (db, tx) =>
    val service = createService(tx)

    // Start should schedule tasks
    service.start(None)

    // Give it a moment to initialize
    Thread.sleep(100)

    // Shutdown should clean up
    service.shutdown()
  }

  private def createService(tx: Transactor): AsyncSyncService =
    AsyncSyncService(
      transactor = tx,
      syncQueueRepo = SyncQueueRepository(),
      syncConflictRepo = SyncConflictRepository(),
      syncHistoryRepo = SyncHistoryRepository(),
      biosampleRepo = BiosampleRepository(),
      projectRepo = ProjectRepository(),
      sequenceRunRepo = SequenceRunRepository(),
      alignmentRepo = AlignmentRepository(),
      conflictNotifier = ConflictNotifier()
    )
