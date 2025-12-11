package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class SyncHistoryRepositorySpec extends FunSuite with DatabaseTestSupport:

  val historyRepo = SyncHistoryRepository()

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()
      val startedAt = LocalDateTime.now().minusSeconds(5)

      val entity = SyncHistoryEntity.create(
        entityType = SyncEntityType.Biosample,
        entityId = entityId,
        operation = SyncOperation.Create,
        direction = SyncDirection.Push,
        status = SyncResultStatus.Success,
        startedAt = startedAt,
        completedAt = LocalDateTime.now(),
        atUri = Some("at://did:plc:test/biosample/1"),
        remoteCid = Some("bafycid123")
      )

      val saved = historyRepo.insert(entity)
      val found = historyRepo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.entityType, SyncEntityType.Biosample)
      assertEquals(found.get.entityId, entityId)
      assertEquals(found.get.operation, SyncOperation.Create)
      assertEquals(found.get.direction, SyncDirection.Push)
      assertEquals(found.get.status, SyncResultStatus.Success)
      assertEquals(found.get.atUri, Some("at://did:plc:test/biosample/1"))
      assert(found.get.durationMs.isDefined)
    }
  }

  testTransactor.test("recordSuccess creates success entry") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()
      val startedAt = LocalDateTime.now().minusSeconds(2)

      val entry = historyRepo.recordSuccess(
        entityType = SyncEntityType.Project,
        entityId = entityId,
        operation = SyncOperation.Update,
        direction = SyncDirection.Push,
        startedAt = startedAt,
        atUri = Some("at://test/1")
      )

      assertEquals(entry.status, SyncResultStatus.Success)
      assertEquals(entry.direction, SyncDirection.Push)
      assert(entry.errorMessage.isEmpty)
    }
  }

  testTransactor.test("recordFailure creates failure entry with error message") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()
      val startedAt = LocalDateTime.now().minusSeconds(1)

      val entry = historyRepo.recordFailure(
        entityType = SyncEntityType.Alignment,
        entityId = entityId,
        operation = SyncOperation.Create,
        direction = SyncDirection.Push,
        startedAt = startedAt,
        errorMessage = "Connection refused"
      )

      assertEquals(entry.status, SyncResultStatus.Failed)
      assertEquals(entry.errorMessage, Some("Connection refused"))
    }
  }

  testTransactor.test("recordConflict creates conflict entry") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()
      val startedAt = LocalDateTime.now().minusSeconds(1)

      val entry = historyRepo.recordConflict(
        entityType = SyncEntityType.Biosample,
        entityId = entityId,
        operation = SyncOperation.Update,
        direction = SyncDirection.Pull,
        startedAt = startedAt,
        localVersionBefore = 3,
        remoteVersionBefore = 5
      )

      assertEquals(entry.status, SyncResultStatus.Conflict)
      assertEquals(entry.localVersionBefore, Some(3))
      assertEquals(entry.remoteVersionBefore, Some(5))
    }
  }

  testTransactor.test("findByEntity returns history for entity") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()
      val startedAt = LocalDateTime.now()

      historyRepo.recordSuccess(SyncEntityType.Biosample, entityId, SyncOperation.Create, SyncDirection.Push, startedAt)
      historyRepo.recordSuccess(SyncEntityType.Biosample, entityId, SyncOperation.Update, SyncDirection.Push, startedAt)
      historyRepo.recordSuccess(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)

      val history = historyRepo.findByEntity(SyncEntityType.Biosample, entityId)
      assertEquals(history.size, 2)
    }
  }

  testTransactor.test("findByStatus returns entries with matching status") { case (db, tx) =>
    tx.readWrite {
      val startedAt = LocalDateTime.now()

      historyRepo.recordSuccess(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)
      historyRepo.recordSuccess(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)
      historyRepo.recordFailure(SyncEntityType.Alignment, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt, "Error")

      val success = historyRepo.findByStatus(SyncResultStatus.Success)
      assertEquals(success.size, 2)

      val failed = historyRepo.findByStatus(SyncResultStatus.Failed)
      assertEquals(failed.size, 1)
    }
  }

  testTransactor.test("findByDirection returns entries with matching direction") { case (db, tx) =>
    tx.readWrite {
      val startedAt = LocalDateTime.now()

      historyRepo.recordSuccess(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)
      historyRepo.recordSuccess(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Update, SyncDirection.Pull, startedAt)
      historyRepo.recordSuccess(SyncEntityType.Alignment, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)

      val pushes = historyRepo.findByDirection(SyncDirection.Push)
      assertEquals(pushes.size, 2)

      val pulls = historyRepo.findByDirection(SyncDirection.Pull)
      assertEquals(pulls.size, 1)
    }
  }

  testTransactor.test("findRecentFailures returns only failures") { case (db, tx) =>
    tx.readWrite {
      val startedAt = LocalDateTime.now()

      historyRepo.recordSuccess(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)
      historyRepo.recordFailure(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt, "Error 1")
      historyRepo.recordFailure(SyncEntityType.Alignment, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt, "Error 2")

      val failures = historyRepo.findRecentFailures()
      assertEquals(failures.size, 2)
      assert(failures.forall(_.status == SyncResultStatus.Failed))
    }
  }

  testTransactor.test("getLastSyncForEntity returns most recent") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()
      val startedAt = LocalDateTime.now()

      historyRepo.recordSuccess(SyncEntityType.Biosample, entityId, SyncOperation.Create, SyncDirection.Push, startedAt)
      Thread.sleep(10) // Ensure different timestamps
      historyRepo.recordSuccess(SyncEntityType.Biosample, entityId, SyncOperation.Update, SyncDirection.Push, startedAt)

      val last = historyRepo.getLastSyncForEntity(SyncEntityType.Biosample, entityId)
      assert(last.isDefined)
      assertEquals(last.get.operation, SyncOperation.Update)
    }
  }

  testTransactor.test("getLastSuccessfulSync returns most recent success") { case (db, tx) =>
    tx.readWrite {
      val entityId = UUID.randomUUID()
      val startedAt = LocalDateTime.now()

      historyRepo.recordSuccess(SyncEntityType.Biosample, entityId, SyncOperation.Create, SyncDirection.Push, startedAt)
      Thread.sleep(10)
      historyRepo.recordFailure(SyncEntityType.Biosample, entityId, SyncOperation.Update, SyncDirection.Push, startedAt, "Error")

      val lastSuccess = historyRepo.getLastSuccessfulSync(SyncEntityType.Biosample, entityId)
      assert(lastSuccess.isDefined)
      assertEquals(lastSuccess.get.operation, SyncOperation.Create)
      assertEquals(lastSuccess.get.status, SyncResultStatus.Success)
    }
  }

  testTransactor.test("countByStatus returns correct counts") { case (db, tx) =>
    tx.readWrite {
      val startedAt = LocalDateTime.now()

      historyRepo.recordSuccess(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)
      historyRepo.recordSuccess(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)
      historyRepo.recordFailure(SyncEntityType.Alignment, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt, "Error")

      val counts = historyRepo.countByStatus()
      assertEquals(counts.getOrElse(SyncResultStatus.Success, 0L), 2L)
      assertEquals(counts.getOrElse(SyncResultStatus.Failed, 0L), 1L)
    }
  }

  testTransactor.test("getStatsForPeriod calculates statistics") { case (db, tx) =>
    tx.readWrite {
      val now = LocalDateTime.now()
      val startedAt = now.minusMinutes(1)

      historyRepo.recordSuccess(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)
      historyRepo.recordSuccess(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Update, SyncDirection.Pull, startedAt)
      historyRepo.recordFailure(SyncEntityType.Alignment, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt, "Error")

      val stats = historyRepo.getStatsForPeriod(now.minusHours(1), now.plusHours(1))
      assertEquals(stats.total, 3)
      assertEquals(stats.successful, 2)
      assertEquals(stats.failed, 1)
      assertEquals(stats.pushes, 2)
      assertEquals(stats.pulls, 1)
    }
  }

  testTransactor.test("cleanupOlderThan removes old entries") { case (db, tx) =>
    tx.readWrite {
      val startedAt = LocalDateTime.now()

      historyRepo.recordSuccess(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)
      historyRepo.recordSuccess(SyncEntityType.Project, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)

      // Cleanup with 0 days should remove all
      val removed = historyRepo.cleanupOlderThan(0)
      assertEquals(removed, 2)

      assertEquals(historyRepo.findAll().size, 0)
    }
  }

  testTransactor.test("delete removes history entry") { case (db, tx) =>
    tx.readWrite {
      val startedAt = LocalDateTime.now()
      val entry = historyRepo.recordSuccess(SyncEntityType.Biosample, UUID.randomUUID(), SyncOperation.Create, SyncDirection.Push, startedAt)

      assert(historyRepo.findById(entry.id).isDefined)

      val deleted = historyRepo.delete(entry.id)
      assert(deleted)
      assertEquals(historyRepo.findById(entry.id), None)
    }
  }
