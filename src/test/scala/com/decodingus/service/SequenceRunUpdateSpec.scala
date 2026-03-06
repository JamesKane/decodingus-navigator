package com.decodingus.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.*
import com.decodingus.service.EntityConversions
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
      alignmentRepo = AlignmentRepository(),
      strProfileRepo = StrProfileRepository(),
      chipProfileRepo = ChipProfileRepository()
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

  testTransactor.test("updateSequenceRun preserves alignmentRefs in workspace") { case (db, tx) =>
    // Test that updating a sequence run doesn't clear its alignmentRefs
    // This is a regression test for the bug where editing metadata would cause
    // alignment files to no longer render under their run.

    val biosampleRepo = BiosampleRepository()
    val sequenceRunRepo = SequenceRunRepository()
    val alignmentRepo = AlignmentRepository()

    val manager = new SequenceDataManager(tx, biosampleRepo, sequenceRunRepo, alignmentRepo)

    // Set up workspace tracking
    var currentWorkspace = Workspace.empty
    manager.setWorkspaceCallbacks(
      getter = () => currentWorkspace,
      updater = ws => currentWorkspace = ws
    )

    // 1. Create Biosample directly in DB
    val biosampleEntity = tx.readWrite {
      biosampleRepo.insert(BiosampleEntity(
        id = UUID.randomUUID(),
        sampleAccession = "ALIGNREF-TEST-001",
        donorIdentifier = "DONOR-001",
        description = None,
        centerName = None,
        sex = None,
        citizenDid = None,
        haplogroups = None,
        meta = EntityMeta.create()
      ))
    }.toOption.get

    // 2. Create SequenceRun via manager
    val initialRun = SequenceRun(
      atUri = None,
      meta = RecordMeta.initial,
      biosampleRef = EntityConversions.localUri("biosample", biosampleEntity.id),
      platformName = "ILLUMINA",
      testType = "WGS"
    )
    val fileInfo = FileInfo("test.bam", Some(1000L), "BAM", None, None, None)
    val createResult = manager.createSequenceRun(biosampleEntity.sampleAccession, initialRun, fileInfo)
    assert(createResult.isRight, s"Failed to create sequence run: ${createResult.left.getOrElse("")}")
    val createdRun = createResult.toOption.get.sequenceRun
    val sequenceRunUri = createdRun.atUri.get

    // 3. Create Alignment for the run
    val alignment = Alignment(
      atUri = None,
      meta = RecordMeta.initial,
      sequenceRunRef = sequenceRunUri,
      biosampleRef = None,
      referenceBuild = "GRCh38",
      aligner = "BWA-MEM"
    )
    val alignResult = manager.createAlignment(sequenceRunUri, alignment)
    assert(alignResult.isRight, s"Failed to create alignment: ${alignResult.left.getOrElse("")}")
    val createdAlignment = alignResult.toOption.get.alignment
    val alignmentUri = createdAlignment.atUri.get

    // 4. Verify workspace has the alignmentRef
    val runBeforeUpdate = currentWorkspace.main.sequenceRuns.find(_.atUri.contains(sequenceRunUri))
    assert(runBeforeUpdate.isDefined, "Sequence run not in workspace")
    assert(runBeforeUpdate.get.alignmentRefs.contains(alignmentUri),
      s"AlignmentRef missing before update: ${runBeforeUpdate.get.alignmentRefs}")

    // 5. Update SequenceRun metadata (like editing testType in the UI)
    val updatedRun = runBeforeUpdate.get.copy(
      testType = "Y_ELITE",
      sequencingFacility = Some("Test Lab")
    )
    val updateResult = manager.updateSequenceRun(updatedRun)
    assert(updateResult.isRight, s"Failed to update sequence run: ${updateResult.left.getOrElse("")}")

    // 6. Verify alignmentRefs are preserved in workspace
    val runAfterUpdate = currentWorkspace.main.sequenceRuns.find(_.atUri.contains(sequenceRunUri))
    assert(runAfterUpdate.isDefined, "Sequence run missing after update")
    assert(runAfterUpdate.get.alignmentRefs.contains(alignmentUri),
      s"AlignmentRef lost after update! Was: ${runBeforeUpdate.get.alignmentRefs}, Now: ${runAfterUpdate.get.alignmentRefs}")

    // 7. Verify other fields were actually updated
    assertEquals(runAfterUpdate.get.testType, "Y_ELITE")
    assertEquals(runAfterUpdate.get.sequencingFacility, Some("Test Lab"))
  }
