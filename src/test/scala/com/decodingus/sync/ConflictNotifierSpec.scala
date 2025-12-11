package com.decodingus.sync

import com.decodingus.repository.{SyncConflictEntity, SyncEntityType, ConflictStatus}
import io.circe.Json
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class ConflictNotifierSpec extends FunSuite:

  // Note: ScalaFX Platform operations are bypassed in tests since we're not on FX thread
  // The observable properties still work, just without thread marshalling

  test("initial state has no conflicts") {
    val notifier = ConflictNotifier()

    assertEquals(notifier.hasConflicts.value, false)
    assertEquals(notifier.conflictCount.value, 0)
    assertEquals(notifier.pendingCount.value, 0)
    assertEquals(notifier.isSyncing.value, false)
    assertEquals(notifier.isOnline.value, true)
    assertEquals(notifier.lastError.value, None)
    assert(notifier.conflicts.isEmpty)
  }

  test("updateCounts updates pending and conflict counts") {
    val notifier = ConflictNotifier()

    // Note: In test context without FX thread, runOnFxThread executes directly
    notifier.updateCounts(5, 2)

    // Values should be updated (may need Platform.runLater in real FX context)
    assertEquals(notifier.pendingCount.value, 5)
    assertEquals(notifier.conflictCount.value, 2)
    assertEquals(notifier.hasConflicts.value, true)
  }

  test("notifyConflicts adds conflicts to list") {
    val notifier = ConflictNotifier()

    val conflicts = List(
      createTestConflict(UUID.randomUUID(), SyncEntityType.Biosample),
      createTestConflict(UUID.randomUUID(), SyncEntityType.Project)
    )

    notifier.notifyConflicts(conflicts)

    assertEquals(notifier.conflicts.size, 2)
    assertEquals(notifier.conflictCount.value, 2)
    assertEquals(notifier.hasConflicts.value, true)
  }

  test("addConflict adds single conflict") {
    val notifier = ConflictNotifier()
    val conflict = createTestConflict(UUID.randomUUID(), SyncEntityType.Biosample)

    notifier.addConflict(conflict)

    assertEquals(notifier.conflicts.size, 1)
    assertEquals(notifier.conflictCount.value, 1)
    assertEquals(notifier.hasConflicts.value, true)
  }

  test("addConflict replaces existing conflict for same entity") {
    val notifier = ConflictNotifier()
    val entityId = UUID.randomUUID()

    val conflict1 = createTestConflict(entityId, SyncEntityType.Biosample)
    val conflict2 = createTestConflict(entityId, SyncEntityType.Biosample)

    notifier.addConflict(conflict1)
    notifier.addConflict(conflict2)

    // Should still be 1 (replaced, not added)
    assertEquals(notifier.conflicts.size, 1)
    assertEquals(notifier.conflictCount.value, 1)
  }

  test("clearConflict removes conflict by entity ID") {
    val notifier = ConflictNotifier()
    val entityId1 = UUID.randomUUID()
    val entityId2 = UUID.randomUUID()

    notifier.addConflict(createTestConflict(entityId1, SyncEntityType.Biosample))
    notifier.addConflict(createTestConflict(entityId2, SyncEntityType.Project))

    assertEquals(notifier.conflicts.size, 2)

    notifier.clearConflict(entityId1)

    assertEquals(notifier.conflicts.size, 1)
    assertEquals(notifier.conflictCount.value, 1)
    assertEquals(notifier.conflicts.head.entityId, entityId2)
  }

  test("clearAllConflicts removes all conflicts") {
    val notifier = ConflictNotifier()

    notifier.addConflict(createTestConflict(UUID.randomUUID(), SyncEntityType.Biosample))
    notifier.addConflict(createTestConflict(UUID.randomUUID(), SyncEntityType.Project))
    notifier.addConflict(createTestConflict(UUID.randomUUID(), SyncEntityType.Alignment))

    assertEquals(notifier.conflicts.size, 3)

    notifier.clearAllConflicts()

    assertEquals(notifier.conflicts.size, 0)
    assertEquals(notifier.conflictCount.value, 0)
    assertEquals(notifier.hasConflicts.value, false)
  }

  test("setSyncStatus updates syncing and online flags") {
    val notifier = ConflictNotifier()

    notifier.setSyncStatus(syncing = true, online = false)

    assertEquals(notifier.isSyncing.value, true)
    assertEquals(notifier.isOnline.value, false)

    notifier.setSyncStatus(syncing = false, online = true)

    assertEquals(notifier.isSyncing.value, false)
    assertEquals(notifier.isOnline.value, true)
  }

  test("setError and clearError manage error message") {
    val notifier = ConflictNotifier()

    assertEquals(notifier.lastError.value, None)

    notifier.setError("Connection failed")

    assertEquals(notifier.lastError.value, Some("Connection failed"))

    notifier.clearError()

    assertEquals(notifier.lastError.value, None)
  }

  test("incrementPending increases count") {
    val notifier = ConflictNotifier()

    assertEquals(notifier.pendingCount.value, 0)

    notifier.incrementPending()
    assertEquals(notifier.pendingCount.value, 1)

    notifier.incrementPending()
    assertEquals(notifier.pendingCount.value, 2)
  }

  test("decrementPending decreases count but not below zero") {
    val notifier = ConflictNotifier()

    notifier.incrementPending()
    notifier.incrementPending()
    assertEquals(notifier.pendingCount.value, 2)

    notifier.decrementPending()
    assertEquals(notifier.pendingCount.value, 1)

    notifier.decrementPending()
    assertEquals(notifier.pendingCount.value, 0)

    // Should not go negative
    notifier.decrementPending()
    assertEquals(notifier.pendingCount.value, 0)
  }

  private def createTestConflict(entityId: UUID, entityType: SyncEntityType): SyncConflictEntity =
    val now = LocalDateTime.now()
    SyncConflictEntity(
      id = UUID.randomUUID(),
      entityType = entityType,
      entityId = entityId,
      atUri = None,
      detectedAt = now,
      localVersion = 1,
      remoteVersion = 2,
      localChanges = Some(Json.obj()),
      remoteChanges = Some(Json.obj()),
      overlappingFields = Some(Json.arr()),
      suggestedResolution = None,
      resolutionReason = None,
      status = ConflictStatus.Unresolved,
      resolvedAt = None,
      resolvedBy = None,
      resolutionAction = None,
      localSnapshot = None,
      remoteSnapshot = None,
      createdAt = now,
      updatedAt = now
    )
