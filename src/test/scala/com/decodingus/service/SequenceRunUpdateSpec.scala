package com.decodingus.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.*
import com.decodingus.workspace.model.*
import munit.FunSuite
import java.util.UUID

class SequenceRunUpdateSpec extends FunSuite with DatabaseTestSupport:

  private def createWorkspaceService(tx: Transactor): H2WorkspaceService =
    H2WorkspaceService(
      transactor = tx,
      biosampleRepo = BiosampleRepository(),
      projectRepo = ProjectRepository(),
      sequenceRunRepo = SequenceRunRepository(),
      alignmentRepo = AlignmentRepository()
    )

  testTransactor.test("updateSequenceRun persists libraryId and platformUnit") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    // 1. Create Biosample
    val biosample = service.createBiosample(
      Biosample(None, RecordMeta.initial, "UPDATE-TEST-001", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    // 2. Create initial SequenceRun (empty libraryId/platformUnit)
    val initialRun = SequenceRun(
      atUri = None,
      meta = RecordMeta.initial,
      biosampleRef = biosample.atUri.get,
      platformName = "Unknown",
      testType = "Unknown",
      libraryId = None,
      platformUnit = None
    )

    val createdRun = service.createSequenceRun(initialRun, biosampleId).toOption.get
    assert(createdRun.libraryId.isEmpty)
    assert(createdRun.platformUnit.isEmpty)

    // 3. Update SequenceRun with new info (mimicking WorkbenchViewModel update)
    val updatedRun = createdRun.copy(
      platformName = "ILLUMINA",
      testType = "WGS",
      libraryId = Some("LIB-123"),
      platformUnit = Some("FLOWCELL.1.BARCODE")
    )

    val result = service.updateSequenceRun(updatedRun)
    assert(result.isRight)

    // 4. Verify returned object
    val saved = result.toOption.get
    assertEquals(saved.libraryId, Some("LIB-123"))
    assertEquals(saved.platformUnit, Some("FLOWCELL.1.BARCODE"))
    assertEquals(saved.platformName, "ILLUMINA")

    // 5. Verify persistence by reloading from DB
    val reloaded = service.getSequenceRun(EntityConversions.parseIdFromRef(saved.atUri.get).get).toOption.flatten.get
    assertEquals(reloaded.libraryId, Some("LIB-123"))
    assertEquals(reloaded.platformUnit, Some("FLOWCELL.1.BARCODE"))
    assertEquals(reloaded.platformName, "ILLUMINA")
  }
