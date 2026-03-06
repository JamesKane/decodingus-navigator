package com.decodingus.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.*
import com.decodingus.workspace.model.*
import munit.FunSuite
import java.util.UUID

class H2WorkspaceServiceSpec extends FunSuite with DatabaseTestSupport:

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

  // ============================================
  // Biosample Operations
  // ============================================

  testTransactor.test("createBiosample creates and returns biosample with ID") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = Biosample(
      atUri = None,
      meta = RecordMeta.initial,
      sampleAccession = "TEST-001",
      donorIdentifier = "DONOR-001",
      description = Some("Test biosample"),
      centerName = Some("Test Center"),
      sex = Some("Male")
    )

    val result = service.createBiosample(biosample)

    assert(result.isRight)
    result.foreach { saved =>
      assert(saved.atUri.isDefined, "Should have atUri assigned")
      assertEquals(saved.sampleAccession, "TEST-001")
      assertEquals(saved.donorIdentifier, "DONOR-001")
      assertEquals(saved.description, Some("Test biosample"))
    }
  }

  testTransactor.test("createBiosample rejects duplicate accession") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = Biosample(None, RecordMeta.initial, "DUPLICATE-001", "DONOR-001")
    service.createBiosample(biosample)

    val result = service.createBiosample(biosample)
    assert(result.isLeft, "Should reject duplicate accession")
  }

  testTransactor.test("getBiosample returns existing biosample") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = Biosample(None, RecordMeta.initial, "GET-TEST-001", "DONOR-001")
    val created = service.createBiosample(biosample).toOption.get
    val id = EntityConversions.parseIdFromRef(created.atUri.get).get

    val result = service.getBiosample(id)

    assert(result.isRight)
    assert(result.toOption.flatten.isDefined)
    assertEquals(result.toOption.flatten.get.sampleAccession, "GET-TEST-001")
  }

  testTransactor.test("getBiosampleByAccession returns matching biosample") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    service.createBiosample(Biosample(None, RecordMeta.initial, "ACCESSION-FIND", "DONOR-001"))

    val result = service.getBiosampleByAccession("ACCESSION-FIND")

    assert(result.isRight)
    assert(result.toOption.flatten.isDefined)
    assertEquals(result.toOption.flatten.get.sampleAccession, "ACCESSION-FIND")
  }

  testTransactor.test("updateBiosample modifies existing biosample") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = Biosample(None, RecordMeta.initial, "UPDATE-TEST", "DONOR-001", description = Some("Original"))
    val created = service.createBiosample(biosample).toOption.get

    val updated = created.copy(description = Some("Updated"))
    val result = service.updateBiosample(updated)

    assert(result.isRight)
    result.foreach { saved =>
      assertEquals(saved.description, Some("Updated"))
    }
  }

  testTransactor.test("deleteBiosample removes biosample") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val created = service.createBiosample(
      Biosample(None, RecordMeta.initial, "DELETE-TEST", "DONOR-001")
    ).toOption.get
    val id = EntityConversions.parseIdFromRef(created.atUri.get).get

    val result = service.deleteBiosample(id)

    assertEquals(result, Right(true))
    assertEquals(service.getBiosample(id), Right(None))
  }

  testTransactor.test("getAllBiosamples returns all biosamples") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    service.createBiosample(Biosample(None, RecordMeta.initial, "LIST-001", "DONOR-001"))
    service.createBiosample(Biosample(None, RecordMeta.initial, "LIST-002", "DONOR-002"))
    service.createBiosample(Biosample(None, RecordMeta.initial, "LIST-003", "DONOR-003"))

    val result = service.getAllBiosamples()

    assert(result.isRight)
    assertEquals(result.toOption.get.size, 3)
  }

  testTransactor.test("updateBiosampleHaplogroups updates haplogroup assignments") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val created = service.createBiosample(
      Biosample(None, RecordMeta.initial, "HAPLO-TEST", "DONOR-001")
    ).toOption.get
    val id = EntityConversions.parseIdFromRef(created.atUri.get).get

    val haplogroups = HaplogroupAssignments(
      yDna = Some(HaplogroupResult(
        haplogroupName = "R-M269",
        score = 0.99,
        treeProvider = Some("FTDNA")
      ))
    )

    val updateResult = service.updateBiosampleHaplogroups(id, haplogroups)
    assertEquals(updateResult, Right(true))

    val found = service.getBiosample(id).toOption.flatten.get
    assertEquals(found.haplogroups.flatMap(_.yDna).map(_.haplogroupName), Some("R-M269"))
  }

  // ============================================
  // Project Operations
  // ============================================

  testTransactor.test("createProject creates and returns project with ID") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val project = Project(
      atUri = None,
      meta = RecordMeta.initial,
      projectName = "Test Project",
      description = Some("A test project"),
      administrator = "did:plc:admin123"
    )

    val result = service.createProject(project)

    assert(result.isRight)
    result.foreach { saved =>
      assert(saved.atUri.isDefined)
      assertEquals(saved.projectName, "Test Project")
    }
  }

  testTransactor.test("createProject rejects duplicate name") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val project = Project(None, RecordMeta.initial, "Duplicate Project", None, "did:plc:admin")
    service.createProject(project)

    val result = service.createProject(project)
    assert(result.isLeft)
  }

  testTransactor.test("getProjectByName returns matching project") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    service.createProject(Project(None, RecordMeta.initial, "Named Project", None, "did:plc:admin"))

    val result = service.getProjectByName("Named Project")

    assert(result.isRight)
    assert(result.toOption.flatten.isDefined)
    assertEquals(result.toOption.flatten.get.projectName, "Named Project")
  }

  testTransactor.test("addProjectMember and getProjectMembers work together") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample1 = service.createBiosample(
      Biosample(None, RecordMeta.initial, "MEMBER-001", "DONOR-001")
    ).toOption.get
    val biosample2 = service.createBiosample(
      Biosample(None, RecordMeta.initial, "MEMBER-002", "DONOR-002")
    ).toOption.get

    val project = service.createProject(
      Project(None, RecordMeta.initial, "Members Project", None, "did:plc:admin")
    ).toOption.get

    val projectId = EntityConversions.parseIdFromRef(project.atUri.get).get
    val biosampleId1 = EntityConversions.parseIdFromRef(biosample1.atUri.get).get
    val biosampleId2 = EntityConversions.parseIdFromRef(biosample2.atUri.get).get

    service.addProjectMember(projectId, biosampleId1)
    service.addProjectMember(projectId, biosampleId2)

    val members = service.getProjectMembers(projectId)

    assert(members.isRight)
    assertEquals(members.toOption.get.size, 2)
  }

  testTransactor.test("removeProjectMember removes member from project") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = service.createBiosample(
      Biosample(None, RecordMeta.initial, "REMOVE-MEMBER", "DONOR-001")
    ).toOption.get
    val project = service.createProject(
      Project(None, RecordMeta.initial, "Remove Member Project", None, "did:plc:admin")
    ).toOption.get

    val projectId = EntityConversions.parseIdFromRef(project.atUri.get).get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    service.addProjectMember(projectId, biosampleId)
    assertEquals(service.getProjectMembers(projectId).toOption.get.size, 1)

    service.removeProjectMember(projectId, biosampleId)
    assertEquals(service.getProjectMembers(projectId).toOption.get.size, 0)
  }

  // ============================================
  // SequenceRun Operations
  // ============================================

  testTransactor.test("createSequenceRun creates sequence run linked to biosample") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = service.createBiosample(
      Biosample(None, RecordMeta.initial, "SEQRUN-BIOSAMPLE", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    val sequenceRun = SequenceRun(
      atUri = None,
      meta = RecordMeta.initial,
      biosampleRef = biosample.atUri.get,
      platformName = "ILLUMINA",
      instrumentModel = Some("NovaSeq 6000"),
      testType = "WGS"
    )

    val result = service.createSequenceRun(sequenceRun, biosampleId)

    assert(result.isRight)
    result.foreach { saved =>
      assert(saved.atUri.isDefined)
      assertEquals(saved.platformName, "ILLUMINA")
      assertEquals(saved.testType, "WGS")
    }
  }

  testTransactor.test("createSequenceRun fails for non-existent biosample") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val sequenceRun = SequenceRun(None, RecordMeta.initial, "fake:ref", "ILLUMINA", testType = "WGS")
    val result = service.createSequenceRun(sequenceRun, UUID.randomUUID())

    assert(result.isLeft)
  }

  testTransactor.test("getSequenceRunsForBiosample returns linked sequence runs") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = service.createBiosample(
      Biosample(None, RecordMeta.initial, "MULTI-SEQRUN", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    service.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "ILLUMINA", testType = "WGS"),
      biosampleId
    )
    service.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "PACBIO", testType = "WGS_HIFI"),
      biosampleId
    )

    val result = service.getSequenceRunsForBiosample(biosampleId)

    assert(result.isRight)
    assertEquals(result.toOption.get.size, 2)
  }

  // ============================================
  // Alignment Operations
  // ============================================

  testTransactor.test("createAlignment creates alignment linked to sequence run") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = service.createBiosample(
      Biosample(None, RecordMeta.initial, "ALIGN-BIOSAMPLE", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    val sequenceRun = service.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "ILLUMINA", testType = "WGS"),
      biosampleId
    ).toOption.get
    val sequenceRunId = EntityConversions.parseIdFromRef(sequenceRun.atUri.get).get

    val alignment = Alignment(
      atUri = None,
      meta = RecordMeta.initial,
      sequenceRunRef = sequenceRun.atUri.get,
      referenceBuild = "GRCh38",
      aligner = "BWA-MEM2",
      variantCaller = Some("DeepVariant")
    )

    val result = service.createAlignment(alignment, sequenceRunId)

    assert(result.isRight)
    result.foreach { saved =>
      assert(saved.atUri.isDefined)
      assertEquals(saved.referenceBuild, "GRCh38")
      assertEquals(saved.aligner, "BWA-MEM2")
    }
  }

  testTransactor.test("createAlignment fails for non-existent sequence run") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val alignment = Alignment(None, RecordMeta.initial, "fake:ref", referenceBuild = "GRCh38", aligner = "BWA")
    val result = service.createAlignment(alignment, UUID.randomUUID())

    assert(result.isLeft)
  }

  testTransactor.test("getAlignmentsForBiosample returns all alignments for biosample") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = service.createBiosample(
      Biosample(None, RecordMeta.initial, "MULTI-ALIGN", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    val seqRun1 = service.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "ILLUMINA", testType = "WGS"),
      biosampleId
    ).toOption.get
    val seqRunId1 = EntityConversions.parseIdFromRef(seqRun1.atUri.get).get

    val seqRun2 = service.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "PACBIO", testType = "WGS_HIFI"),
      biosampleId
    ).toOption.get
    val seqRunId2 = EntityConversions.parseIdFromRef(seqRun2.atUri.get).get

    service.createAlignment(
      Alignment(None, RecordMeta.initial, seqRun1.atUri.get, referenceBuild = "GRCh38", aligner = "BWA-MEM2"),
      seqRunId1
    )
    service.createAlignment(
      Alignment(None, RecordMeta.initial, seqRun2.atUri.get, referenceBuild = "GRCh38", aligner = "minimap2"),
      seqRunId2
    )

    val result = service.getAlignmentsForBiosample(biosampleId)

    assert(result.isRight)
    assertEquals(result.toOption.get.size, 2)
  }

  testTransactor.test("updateAlignmentMetrics updates metrics for alignment") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = service.createBiosample(
      Biosample(None, RecordMeta.initial, "METRICS-TEST", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    val seqRun = service.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "ILLUMINA", testType = "WGS"),
      biosampleId
    ).toOption.get
    val seqRunId = EntityConversions.parseIdFromRef(seqRun.atUri.get).get

    val alignment = service.createAlignment(
      Alignment(None, RecordMeta.initial, seqRun.atUri.get, referenceBuild = "GRCh38", aligner = "BWA-MEM2"),
      seqRunId
    ).toOption.get
    val alignmentId = EntityConversions.parseIdFromRef(alignment.atUri.get).get

    val metrics = AlignmentMetrics(
      genomeTerritory = Some(3000000000L),
      meanCoverage = Some(30.5),
      medianCoverage = Some(30.0),
      pctExcDupe = Some(5.2)
    )

    val updateResult = service.updateAlignmentMetrics(alignmentId, metrics)
    assertEquals(updateResult, Right(true))

    val found = service.getAlignment(alignmentId).toOption.flatten.get
    assertEquals(found.metrics.flatMap(_.genomeTerritory), Some(3000000000L))
    assertEquals(found.metrics.flatMap(_.meanCoverage), Some(30.5))
  }

  // ============================================
  // Cascade Delete Operations
  // ============================================

  testTransactor.test("deleteBiosample cascades to sequence runs and alignments") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = service.createBiosample(
      Biosample(None, RecordMeta.initial, "CASCADE-TEST", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    val seqRun = service.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "ILLUMINA", testType = "WGS"),
      biosampleId
    ).toOption.get
    val seqRunId = EntityConversions.parseIdFromRef(seqRun.atUri.get).get

    val alignment = service.createAlignment(
      Alignment(None, RecordMeta.initial, seqRun.atUri.get, referenceBuild = "GRCh38", aligner = "BWA"),
      seqRunId
    ).toOption.get
    val alignmentId = EntityConversions.parseIdFromRef(alignment.atUri.get).get

    // Verify they exist
    assert(service.getSequenceRun(seqRunId).toOption.flatten.isDefined)
    assert(service.getAlignment(alignmentId).toOption.flatten.isDefined)

    // Delete biosample
    service.deleteBiosample(biosampleId)

    // Verify cascade
    assertEquals(service.getBiosample(biosampleId), Right(None))
    assertEquals(service.getSequenceRun(seqRunId), Right(None))
    assertEquals(service.getAlignment(alignmentId), Right(None))
  }

  testTransactor.test("deleteSequenceRun cascades to alignments") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = service.createBiosample(
      Biosample(None, RecordMeta.initial, "SEQRUN-CASCADE", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    val seqRun = service.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "ILLUMINA", testType = "WGS"),
      biosampleId
    ).toOption.get
    val seqRunId = EntityConversions.parseIdFromRef(seqRun.atUri.get).get

    val alignment = service.createAlignment(
      Alignment(None, RecordMeta.initial, seqRun.atUri.get, referenceBuild = "GRCh38", aligner = "BWA"),
      seqRunId
    ).toOption.get
    val alignmentId = EntityConversions.parseIdFromRef(alignment.atUri.get).get

    // Delete sequence run
    service.deleteSequenceRun(seqRunId)

    // Biosample should remain
    assert(service.getBiosample(biosampleId).toOption.flatten.isDefined)
    // Sequence run and alignment should be gone
    assertEquals(service.getSequenceRun(seqRunId), Right(None))
    assertEquals(service.getAlignment(alignmentId), Right(None))
  }

  // ============================================
  // Sync Status Operations
  // ============================================

  testTransactor.test("getSyncStatusSummary returns correct counts") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    // Create some entities (all start as Local)
    service.createBiosample(Biosample(None, RecordMeta.initial, "SYNC-001", "DONOR-001"))
    service.createBiosample(Biosample(None, RecordMeta.initial, "SYNC-002", "DONOR-002"))
    service.createProject(Project(None, RecordMeta.initial, "Sync Project", None, "did:plc:admin"))

    val result = service.getSyncStatusSummary()

    assert(result.isRight)
    result.foreach { summary =>
      assertEquals(summary.localCount, 3) // 2 biosamples + 1 project
      assertEquals(summary.syncedCount, 0)
      assertEquals(summary.modifiedCount, 0)
      assertEquals(summary.conflictCount, 0)
    }
  }

  testTransactor.test("getPendingSyncEntities returns local entities") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    service.createBiosample(Biosample(None, RecordMeta.initial, "PENDING-001", "DONOR-001"))
    service.createProject(Project(None, RecordMeta.initial, "Pending Project", None, "did:plc:admin"))

    val result = service.getPendingSyncEntities()

    assert(result.isRight)
    result.foreach { pending =>
      assertEquals(pending.biosamples.size, 1)
      assertEquals(pending.projects.size, 1)
    }
  }

  // ============================================
  // Bulk Operations
  // ============================================

  testTransactor.test("loadWorkspaceContent returns full workspace") { case (db, tx) =>
    val service = createWorkspaceService(tx)

    val biosample = service.createBiosample(
      Biosample(None, RecordMeta.initial, "WORKSPACE-001", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    service.createProject(Project(None, RecordMeta.initial, "Workspace Project", None, "did:plc:admin"))

    val seqRun = service.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "ILLUMINA", testType = "WGS"),
      biosampleId
    ).toOption.get
    val seqRunId = EntityConversions.parseIdFromRef(seqRun.atUri.get).get

    service.createAlignment(
      Alignment(None, RecordMeta.initial, seqRun.atUri.get, referenceBuild = "GRCh38", aligner = "BWA"),
      seqRunId
    )

    val result = service.loadWorkspaceContent()

    assert(result.isRight)
    result.foreach { content =>
      assertEquals(content.samples.size, 1)
      assertEquals(content.projects.size, 1)
      assertEquals(content.sequenceRuns.size, 1)
      assertEquals(content.alignments.size, 1)
    }
  }
