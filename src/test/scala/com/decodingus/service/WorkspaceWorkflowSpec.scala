package com.decodingus.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.*
import com.decodingus.workspace.model.*
import munit.FunSuite
import java.util.UUID

/**
 * End-to-end workflow tests that verify integration between services.
 *
 * These tests simulate real user workflows:
 * - Creating samples and running analysis
 * - Sync queue integration
 * - Cache management during analysis
 */
class WorkspaceWorkflowSpec extends FunSuite with DatabaseTestSupport:

  // ============================================
  // Service Factory
  // ============================================

  private def createServices(tx: Transactor): WorkspaceServices =
    WorkspaceServices(
      workspace = H2WorkspaceService(
        tx, BiosampleRepository(), ProjectRepository(),
        SequenceRunRepository(), AlignmentRepository(),
        StrProfileRepository(), ChipProfileRepository()
      ),
      sync = H2SyncService(
        tx, SyncQueueRepository(), SyncHistoryRepository(), SyncConflictRepository()
      ),
      cache = H2CacheService(
        tx, AnalysisArtifactRepository(), SourceFileRepository()
      )
    )

  case class WorkspaceServices(
    workspace: H2WorkspaceService,
    sync: H2SyncService,
    cache: H2CacheService
  )

  // ============================================
  // Workflow: Sample Import and Analysis
  // ============================================

  testTransactor.test("workflow: import sample, create sequence run, add alignment, run analysis") { case (db, tx) =>
    val services = createServices(tx)

    // Step 1: Create biosample from imported data
    val biosample = services.workspace.createBiosample(
      Biosample(
        atUri = None,
        meta = RecordMeta.initial,
        sampleAccession = "WGS-SAMPLE-001",
        donorIdentifier = "DONOR-001",
        description = Some("Imported WGS sample"),
        sex = Some("Male")
      )
    ).toOption.get

    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    // Step 2: Create sequence run for the BAM file
    val sequenceRun = services.workspace.createSequenceRun(
      SequenceRun(
        atUri = None,
        meta = RecordMeta.initial,
        biosampleRef = biosample.atUri.get,
        platformName = "ILLUMINA",
        instrumentModel = Some("NovaSeq 6000"),
        testType = "WGS",
        totalReads = Some(500000000L)
      ),
      biosampleId
    ).toOption.get

    val sequenceRunId = EntityConversions.parseIdFromRef(sequenceRun.atUri.get).get

    // Step 3: Create alignment
    val alignment = services.workspace.createAlignment(
      Alignment(
        atUri = None,
        meta = RecordMeta.initial,
        sequenceRunRef = sequenceRun.atUri.get,
        referenceBuild = "GRCh38",
        aligner = "BWA-MEM2",
        variantCaller = Some("DeepVariant")
      ),
      sequenceRunId
    ).toOption.get

    val alignmentId = EntityConversions.parseIdFromRef(alignment.atUri.get).get

    // Step 4: Register source file
    val sourceFile = services.cache.registerSourceFile(
      filePath = "/data/samples/WGS-SAMPLE-001.bam",
      fileChecksum = "sha256:abc123def456",
      fileSize = Some(50000000000L),
      fileFormat = Some(SourceFileFormat.Bam)
    ).toOption.get

    // Step 5: Start analysis artifact (simulating WGS metrics analysis)
    val artifact = services.cache.startArtifact(
      alignmentId = alignmentId,
      artifactType = ArtifactType.WgsMetrics,
      cachePath = s"wgs/${alignmentId}/metrics.txt",
      generatorVersion = Some("GATK-4.6.2"),
      dependsOnSourceChecksum = Some(sourceFile.fileChecksum),
      dependsOnReferenceBuild = Some("GRCh38")
    ).toOption.get

    // Step 6: Complete analysis and update metrics
    services.cache.completeArtifact(artifact.id, 5000L, "metrics-checksum", Some("TXT"))

    val metrics = AlignmentMetrics(
      genomeTerritory = Some(3000000000L),
      meanCoverage = Some(32.5),
      medianCoverage = Some(32.0),
      pctExcDupe = Some(4.2)
    )
    services.workspace.updateAlignmentMetrics(alignmentId, metrics)

    // Step 7: Update haplogroups after analysis
    val haplogroups = HaplogroupAssignments(
      yDna = Some(HaplogroupResult(
        haplogroupName = "R-M269",
        score = 0.99,
        treeProvider = Some("FTDNA")
      )),
      mtDna = Some(HaplogroupResult(
        haplogroupName = "H1a",
        score = 0.95
      ))
    )
    services.workspace.updateBiosampleHaplogroups(biosampleId, haplogroups)

    // Verify final state
    val finalBiosample = services.workspace.getBiosample(biosampleId).toOption.flatten.get
    assert(finalBiosample.haplogroups.flatMap(_.yDna).isDefined)
    assertEquals(finalBiosample.haplogroups.flatMap(_.yDna).map(_.haplogroupName), Some("R-M269"))

    val finalAlignment = services.workspace.getAlignment(alignmentId).toOption.flatten.get
    assertEquals(finalAlignment.metrics.flatMap(_.meanCoverage), Some(32.5))

    val artifactResult = services.cache.getArtifact(alignmentId, ArtifactType.WgsMetrics)
    assert(artifactResult.toOption.flatten.exists(_.status == ArtifactStatus.Available))
  }

  // ============================================
  // Workflow: Sync Queue Integration
  // ============================================

  testTransactor.test("workflow: create entities and queue for sync") { case (db, tx) =>
    val services = createServices(tx)

    // Create biosample
    val biosample = services.workspace.createBiosample(
      Biosample(None, RecordMeta.initial, "SYNC-WORKFLOW-001", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    // Queue for sync
    services.sync.enqueuePush(SyncEntityType.Biosample, biosampleId, SyncOperation.Create)

    // Create project and add member
    val project = services.workspace.createProject(
      Project(None, RecordMeta.initial, "Sync Workflow Project", None, "did:plc:admin")
    ).toOption.get
    val projectId = EntityConversions.parseIdFromRef(project.atUri.get).get

    services.workspace.addProjectMember(projectId, biosampleId)
    services.sync.enqueuePush(SyncEntityType.Project, projectId, SyncOperation.Create)

    // Check sync status
    val pendingCount = services.sync.getPendingCount()
    assertEquals(pendingCount, Right(2L))

    // Get sync batch
    val batch = services.sync.getNextBatch(10)
    assert(batch.isRight)
    assertEquals(batch.toOption.get.size, 2)
  }

  // ============================================
  // Workflow: Cache Invalidation on Source Change
  // ============================================

  testTransactor.test("workflow: source file change triggers cache invalidation") { case (db, tx) =>
    val services = createServices(tx)

    // Setup: Create biosample → sequence run → alignment
    val biosample = services.workspace.createBiosample(
      Biosample(None, RecordMeta.initial, "CACHE-INVALIDATE-001", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    val sequenceRun = services.workspace.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "ILLUMINA", testType = "WGS"),
      biosampleId
    ).toOption.get
    val sequenceRunId = EntityConversions.parseIdFromRef(sequenceRun.atUri.get).get

    val alignment = services.workspace.createAlignment(
      Alignment(None, RecordMeta.initial, sequenceRun.atUri.get, referenceBuild = "GRCh38", aligner = "BWA-MEM2"),
      sequenceRunId
    ).toOption.get
    val alignmentId = EntityConversions.parseIdFromRef(alignment.atUri.get).get

    val originalChecksum = "sha256:original123"

    // Create artifacts dependent on source checksum
    val wgsArtifact = services.cache.startArtifact(
      alignmentId, ArtifactType.WgsMetrics, "wgs/metrics.txt",
      dependsOnSourceChecksum = Some(originalChecksum)
    ).toOption.get
    services.cache.completeArtifact(wgsArtifact.id, 1000L, "wgs-checksum", None)

    val callableLociArtifact = services.cache.startArtifact(
      alignmentId, ArtifactType.CallableLoci, "callable/loci.bed",
      dependsOnSourceChecksum = Some(originalChecksum)
    ).toOption.get
    services.cache.completeArtifact(callableLociArtifact.id, 2000L, "callable-checksum", None)

    // Verify artifacts are available
    assert(services.cache.getArtifact(alignmentId, ArtifactType.WgsMetrics)
      .toOption.flatten.exists(_.status == ArtifactStatus.Available))

    // Simulate source file change by invalidating checksum
    val invalidatedCount = services.cache.invalidateBySourceChecksum(
      originalChecksum,
      "Source BAM file was modified"
    )

    assertEquals(invalidatedCount, Right(2))

    // Verify artifacts are now stale
    val wgsResult = services.cache.getArtifact(alignmentId, ArtifactType.WgsMetrics)
    assert(wgsResult.toOption.flatten.exists(_.status == ArtifactStatus.Stale))

    val callableResult = services.cache.getArtifact(alignmentId, ArtifactType.CallableLoci)
    assert(callableResult.toOption.flatten.exists(_.status == ArtifactStatus.Stale))
  }

  // ============================================
  // Workflow: Project Management
  // ============================================

  testTransactor.test("workflow: create project, add members, delete project") { case (db, tx) =>
    val services = createServices(tx)

    // Create biosamples
    val biosample1 = services.workspace.createBiosample(
      Biosample(None, RecordMeta.initial, "PROJECT-MEMBER-001", "DONOR-001")
    ).toOption.get
    val biosample2 = services.workspace.createBiosample(
      Biosample(None, RecordMeta.initial, "PROJECT-MEMBER-002", "DONOR-002")
    ).toOption.get

    val id1 = EntityConversions.parseIdFromRef(biosample1.atUri.get).get
    val id2 = EntityConversions.parseIdFromRef(biosample2.atUri.get).get

    // Create project
    val project = services.workspace.createProject(
      Project(
        atUri = None,
        meta = RecordMeta.initial,
        projectName = "Research Project Alpha",
        description = Some("Testing haplogroup distribution"),
        administrator = "did:plc:researcher"
      )
    ).toOption.get
    val projectId = EntityConversions.parseIdFromRef(project.atUri.get).get

    // Add members
    services.workspace.addProjectMember(projectId, id1)
    services.workspace.addProjectMember(projectId, id2)

    // Verify membership
    val members = services.workspace.getProjectMembers(projectId).toOption.get
    assertEquals(members.size, 2)

    // Verify project shows up with members in getAllProjects
    val allProjects = services.workspace.getAllProjects().toOption.get
    val foundProject = allProjects.find(_.projectName == "Research Project Alpha")
    assert(foundProject.isDefined)
    assertEquals(foundProject.get.memberRefs.size, 2)

    // Delete project
    services.workspace.deleteProject(projectId)

    // Verify project deleted but biosamples remain
    assertEquals(services.workspace.getProject(projectId), Right(None))
    assert(services.workspace.getBiosample(id1).toOption.flatten.isDefined)
    assert(services.workspace.getBiosample(id2).toOption.flatten.isDefined)
  }

  // ============================================
  // Workflow: Multi-Alignment Analysis
  // ============================================

  testTransactor.test("workflow: single sample with multiple alignments") { case (db, tx) =>
    val services = createServices(tx)

    // One biosample
    val biosample = services.workspace.createBiosample(
      Biosample(None, RecordMeta.initial, "MULTI-ALIGN-001", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    // One sequence run
    val seqRun = services.workspace.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "ILLUMINA", testType = "WGS"),
      biosampleId
    ).toOption.get
    val seqRunId = EntityConversions.parseIdFromRef(seqRun.atUri.get).get

    // Multiple alignments to different references
    val grch38Alignment = services.workspace.createAlignment(
      Alignment(None, RecordMeta.initial, seqRun.atUri.get, referenceBuild = "GRCh38", aligner = "BWA-MEM2"),
      seqRunId
    ).toOption.get

    val t2tAlignment = services.workspace.createAlignment(
      Alignment(None, RecordMeta.initial, seqRun.atUri.get, referenceBuild = "T2T-CHM13", aligner = "BWA-MEM2"),
      seqRunId
    ).toOption.get

    val grch38Id = EntityConversions.parseIdFromRef(grch38Alignment.atUri.get).get
    val t2tId = EntityConversions.parseIdFromRef(t2tAlignment.atUri.get).get

    // Create artifacts for each alignment
    val grch38Artifact = services.cache.startArtifact(
      grch38Id, ArtifactType.HaplogroupVcf, "grch38/haplogroup.vcf",
      dependsOnReferenceBuild = Some("GRCh38")
    ).toOption.get
    services.cache.completeArtifact(grch38Artifact.id, 1000L, "grch38-vcf", Some("VCF"))

    val t2tArtifact = services.cache.startArtifact(
      t2tId, ArtifactType.HaplogroupVcf, "t2t/haplogroup.vcf",
      dependsOnReferenceBuild = Some("T2T-CHM13")
    ).toOption.get
    services.cache.completeArtifact(t2tArtifact.id, 1100L, "t2t-vcf", Some("VCF"))

    // Verify alignments for biosample
    val alignments = services.workspace.getAlignmentsForBiosample(biosampleId).toOption.get
    assertEquals(alignments.size, 2)
    assert(alignments.exists(_.referenceBuild == "GRCh38"))
    assert(alignments.exists(_.referenceBuild == "CHM13v2"))

    // Invalidate by reference build (e.g., reference genome update)
    val invalidated = services.cache.invalidateByReferenceBuild("GRCh38", "GRCh38 reference updated")
    assertEquals(invalidated, Right(1))

    // Verify only GRCh38 artifact is stale
    val grch38Result = services.cache.getArtifact(grch38Id, ArtifactType.HaplogroupVcf)
    assert(grch38Result.toOption.flatten.exists(_.status == ArtifactStatus.Stale))

    val t2tResult = services.cache.getArtifact(t2tId, ArtifactType.HaplogroupVcf)
    assert(t2tResult.toOption.flatten.exists(_.status == ArtifactStatus.Available))
  }

  // ============================================
  // Workflow: Sync Conflict Resolution
  // ============================================

  testTransactor.test("workflow: detect and resolve sync conflict") { case (db, tx) =>
    val services = createServices(tx)

    // Create and queue entity
    val biosample = services.workspace.createBiosample(
      Biosample(None, RecordMeta.initial, "CONFLICT-001", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    val queueEntry = services.sync.enqueuePush(
      SyncEntityType.Biosample, biosampleId, SyncOperation.Create
    ).toOption.get

    services.sync.startProcessing(queueEntry.id)

    // Simulate conflict detection during sync
    val conflict = services.sync.recordConflict(
      entityType = SyncEntityType.Biosample,
      entityId = biosampleId,
      localVersion = 1,
      remoteVersion = 2,
      atUri = Some("at://did:plc:remote/biosample/123")
    ).toOption.get

    // Check sync health
    val isHealthy = services.sync.isSyncHealthy()
    assertEquals(isHealthy, Right(false))

    val unresolvedCount = services.sync.getUnresolvedConflictCount()
    assertEquals(unresolvedCount, Right(1L))

    // Resolve conflict by keeping local
    services.sync.resolveKeepLocal(conflict.id)

    // Sync should be healthy again
    assertEquals(services.sync.isSyncHealthy(), Right(true))
    assertEquals(services.sync.getUnresolvedConflictCount(), Right(0L))
  }

  // ============================================
  // Workflow: Full Sync Status
  // ============================================

  testTransactor.test("workflow: comprehensive sync status tracking") { case (db, tx) =>
    val services = createServices(tx)

    // Create multiple entities
    val biosample1 = services.workspace.createBiosample(
      Biosample(None, RecordMeta.initial, "STATUS-001", "DONOR-001")
    ).toOption.get
    val biosample2 = services.workspace.createBiosample(
      Biosample(None, RecordMeta.initial, "STATUS-002", "DONOR-002")
    ).toOption.get

    val project = services.workspace.createProject(
      Project(None, RecordMeta.initial, "Status Project", None, "did:plc:admin")
    ).toOption.get

    // Check workspace sync summary (all should be Local)
    val summary = services.workspace.getSyncStatusSummary().toOption.get
    assertEquals(summary.localCount, 3)
    assertEquals(summary.syncedCount, 0)

    // Check pending entities
    val pending = services.workspace.getPendingSyncEntities().toOption.get
    assertEquals(pending.biosamples.size, 2)
    assertEquals(pending.projects.size, 1)

    // Queue items for sync
    val id1 = EntityConversions.parseIdFromRef(biosample1.atUri.get).get
    val id2 = EntityConversions.parseIdFromRef(biosample2.atUri.get).get
    val projectId = EntityConversions.parseIdFromRef(project.atUri.get).get

    services.sync.enqueuePush(SyncEntityType.Biosample, id1, SyncOperation.Create)
    services.sync.enqueuePush(SyncEntityType.Biosample, id2, SyncOperation.Create)
    services.sync.enqueuePush(SyncEntityType.Project, projectId, SyncOperation.Create)

    // Check sync service status
    val syncStatus = services.sync.getSyncStatus().toOption.get
    assertEquals(syncStatus.pendingCount, 3L)
    assertEquals(syncStatus.inProgressCount, 0L)
    assert(syncStatus.isHealthy)
  }

  // ============================================
  // Workflow: Multi-Reference Build Alignment Handling
  // ============================================

  testTransactor.test("workflow: same sequencing run with alignments to multiple reference builds") { case (db, tx) =>
    val services = createServices(tx)

    // Step 1: Create biosample
    val biosample = services.workspace.createBiosample(
      Biosample(
        atUri = None,
        meta = RecordMeta.initial,
        sampleAccession = "MULTI-REF-001",
        donorIdentifier = "DONOR-001",
        description = Some("Sample with multiple reference alignments"),
        sex = Some("Male")
      )
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    // Step 2: Create ONE sequence run (represents the physical sequencing)
    // This is the key concept: one sequencing event, multiple alignments
    val sequenceRun = services.workspace.createSequenceRun(
      SequenceRun(
        atUri = None,
        meta = RecordMeta.initial,
        biosampleRef = biosample.atUri.get,
        platformName = "ILLUMINA",
        instrumentModel = Some("NovaSeq 6000"),
        testType = "WGS",
        runFingerprint = Some("LB:WGS_LIB_001_PU:HVMK5DSX2.1.AAACGGCG_SM:MULTI-REF-001"),
        platformUnit = Some("HVMK5DSX2.1.AAACGGCG"),
        libraryId = Some("WGS_LIB_001"),
        sampleName = Some("MULTI-REF-001"),
        totalReads = Some(500000000L)
      ),
      biosampleId
    ).toOption.get
    val sequenceRunId = EntityConversions.parseIdFromRef(sequenceRun.atUri.get).get

    // Step 3: Add GRCh38 alignment (first file)
    val grch38Alignment = services.workspace.createAlignment(
      Alignment(
        atUri = None,
        meta = RecordMeta.initial,
        sequenceRunRef = sequenceRun.atUri.get,
        biosampleRef = Some(biosample.atUri.get),
        referenceBuild = "GRCh38",
        aligner = "BWA-MEM2",
        variantCaller = Some("DeepVariant"),
        files = List(FileInfo(
          fileName = "sample_grch38.bam",
          fileSizeBytes = Some(50000000000L),
          fileFormat = "BAM",
          checksum = Some("sha256:grch38checksum123"),
          location = Some("/data/samples/sample_grch38.bam")
        ))
      ),
      sequenceRunId
    ).toOption.get

    // Step 4: Add GRCh37 alignment (second file, SAME sequencing run)
    val grch37Alignment = services.workspace.createAlignment(
      Alignment(
        atUri = None,
        meta = RecordMeta.initial,
        sequenceRunRef = sequenceRun.atUri.get,  // Same sequence run!
        biosampleRef = Some(biosample.atUri.get),
        referenceBuild = "GRCh37",
        aligner = "BWA-MEM2",
        variantCaller = Some("GATK HaplotypeCaller"),
        files = List(FileInfo(
          fileName = "sample_grch37.bam",
          fileSizeBytes = Some(50000000000L),
          fileFormat = "BAM",
          checksum = Some("sha256:grch37checksum456"),
          location = Some("/data/samples/sample_grch37.bam")
        ))
      ),
      sequenceRunId
    ).toOption.get

    // Step 5: Add T2T-CHM13 alignment (third file, SAME sequencing run)
    val t2tAlignment = services.workspace.createAlignment(
      Alignment(
        atUri = None,
        meta = RecordMeta.initial,
        sequenceRunRef = sequenceRun.atUri.get,  // Same sequence run!
        biosampleRef = Some(biosample.atUri.get),
        referenceBuild = "T2T-CHM13",
        aligner = "BWA-MEM2",
        files = List(FileInfo(
          fileName = "sample_t2t.bam",
          fileSizeBytes = Some(50000000000L),
          fileFormat = "BAM",
          checksum = Some("sha256:t2tchecksum789"),
          location = Some("/data/samples/sample_t2t.bam")
        ))
      ),
      sequenceRunId
    ).toOption.get

    // Verify: All alignments belong to the SAME sequence run
    val alignmentsForRun = services.workspace.getAlignmentsForSequenceRun(sequenceRunId).toOption.get
    assertEquals(alignmentsForRun.size, 3, "Should have 3 alignments for 1 sequence run")

    // Verify: All alignments are for the same biosample
    val alignmentsForBiosample = services.workspace.getAlignmentsForBiosample(biosampleId).toOption.get
    assertEquals(alignmentsForBiosample.size, 3)

    // Verify: Reference builds are different
    val refBuilds = alignmentsForBiosample.map(_.referenceBuild).toSet
    assertEquals(refBuilds, Set("GRCh38", "GRCh37", "CHM13v2"))

    // Verify: Still only ONE sequence run for this biosample
    val sequenceRunsForBiosample = services.workspace.getSequenceRunsForBiosample(biosampleId).toOption.get
    assertEquals(sequenceRunsForBiosample.size, 1, "Should have only 1 sequence run, not 3")

    // Verify: Sequence run has correct metadata
    val run = sequenceRunsForBiosample.head
    assertEquals(run.platformUnit, Some("HVMK5DSX2.1.AAACGGCG"))
    assertEquals(run.libraryId, Some("WGS_LIB_001"))
  }

  testTransactor.test("workflow: two different sequencing runs should NOT be merged") { case (db, tx) =>
    val services = createServices(tx)

    // Create biosample
    val biosample = services.workspace.createBiosample(
      Biosample(None, RecordMeta.initial, "TWO-RUNS-001", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    // First sequencing run (e.g., initial WGS)
    val run1 = services.workspace.createSequenceRun(
      SequenceRun(
        atUri = None,
        meta = RecordMeta.initial,
        biosampleRef = biosample.atUri.get,
        platformName = "ILLUMINA",
        testType = "WGS",
        runFingerprint = Some("LB:LIB_001_PU:FLOWCELL_A.1.AAAA_SM:TWO-RUNS-001"),
        platformUnit = Some("FLOWCELL_A.1.AAAA"),
        libraryId = Some("LIB_001"),
        totalReads = Some(300000000L)
      ),
      biosampleId
    ).toOption.get
    val run1Id = EntityConversions.parseIdFromRef(run1.atUri.get).get

    // Second sequencing run (e.g., re-sequencing at higher depth)
    val run2 = services.workspace.createSequenceRun(
      SequenceRun(
        atUri = None,
        meta = RecordMeta.initial,
        biosampleRef = biosample.atUri.get,
        platformName = "ILLUMINA",
        testType = "WGS",
        runFingerprint = Some("LB:LIB_002_PU:FLOWCELL_B.1.TTTT_SM:TWO-RUNS-001"),
        platformUnit = Some("FLOWCELL_B.1.TTTT"),  // Different flowcell!
        libraryId = Some("LIB_002"),  // Different library!
        totalReads = Some(600000000L)  // Higher depth
      ),
      biosampleId
    ).toOption.get
    val run2Id = EntityConversions.parseIdFromRef(run2.atUri.get).get

    // Add alignment to first run
    services.workspace.createAlignment(
      Alignment(None, RecordMeta.initial, run1.atUri.get, None, "GRCh38", "BWA-MEM2"),
      run1Id
    )

    // Add alignment to second run
    services.workspace.createAlignment(
      Alignment(None, RecordMeta.initial, run2.atUri.get, None, "GRCh38", "BWA-MEM2"),
      run2Id
    )

    // Verify: Should have TWO sequence runs (not merged)
    val runs = services.workspace.getSequenceRunsForBiosample(biosampleId).toOption.get
    assertEquals(runs.size, 2, "Different sequencing runs should NOT be merged")

    // Verify: Each run has different fingerprint data
    val fingerprints = runs.flatMap(_.runFingerprint).toSet
    assertEquals(fingerprints.size, 2)

    val platformUnits = runs.flatMap(_.platformUnit).toSet
    assertEquals(platformUnits, Set("FLOWCELL_A.1.AAAA", "FLOWCELL_B.1.TTTT"))
  }

  testTransactor.test("workflow: alignments track their reference build correctly for cache invalidation") { case (db, tx) =>
    val services = createServices(tx)

    // Setup biosample → sequence run → multiple alignments
    val biosample = services.workspace.createBiosample(
      Biosample(None, RecordMeta.initial, "CACHE-MULTI-REF-001", "DONOR-001")
    ).toOption.get
    val biosampleId = EntityConversions.parseIdFromRef(biosample.atUri.get).get

    val seqRun = services.workspace.createSequenceRun(
      SequenceRun(None, RecordMeta.initial, biosample.atUri.get, "ILLUMINA", testType = "WGS"),
      biosampleId
    ).toOption.get
    val seqRunId = EntityConversions.parseIdFromRef(seqRun.atUri.get).get

    // Two alignments with different reference builds
    val grch38Align = services.workspace.createAlignment(
      Alignment(None, RecordMeta.initial, seqRun.atUri.get, None, "GRCh38", "BWA-MEM2"),
      seqRunId
    ).toOption.get
    val grch38Id = EntityConversions.parseIdFromRef(grch38Align.atUri.get).get

    val grch37Align = services.workspace.createAlignment(
      Alignment(None, RecordMeta.initial, seqRun.atUri.get, None, "GRCh37", "BWA-MEM2"),
      seqRunId
    ).toOption.get
    val grch37Id = EntityConversions.parseIdFromRef(grch37Align.atUri.get).get

    // Create artifacts for each alignment (dependent on their specific reference)
    val grch38Artifact = services.cache.startArtifact(
      grch38Id, ArtifactType.HaplogroupVcf, "grch38/haplogroup.vcf",
      dependsOnReferenceBuild = Some("GRCh38")
    ).toOption.get
    services.cache.completeArtifact(grch38Artifact.id, 1000L, "grch38-vcf", Some("VCF"))

    val grch37Artifact = services.cache.startArtifact(
      grch37Id, ArtifactType.HaplogroupVcf, "grch37/haplogroup.vcf",
      dependsOnReferenceBuild = Some("GRCh37")
    ).toOption.get
    services.cache.completeArtifact(grch37Artifact.id, 1000L, "grch37-vcf", Some("VCF"))

    // Invalidate only GRCh38 artifacts (e.g., reference genome patch)
    val invalidated = services.cache.invalidateByReferenceBuild("GRCh38", "GRCh38 patch applied")
    assertEquals(invalidated, Right(1))

    // Verify: GRCh38 artifact is stale
    val grch38Result = services.cache.getArtifact(grch38Id, ArtifactType.HaplogroupVcf)
    assert(grch38Result.toOption.flatten.exists(_.status == ArtifactStatus.Stale))

    // Verify: GRCh37 artifact is still available (different reference)
    val grch37Result = services.cache.getArtifact(grch37Id, ArtifactType.HaplogroupVcf)
    assert(grch37Result.toOption.flatten.exists(_.status == ArtifactStatus.Available))
  }
