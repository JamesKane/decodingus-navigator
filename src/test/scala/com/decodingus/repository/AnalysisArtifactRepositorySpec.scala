package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import munit.FunSuite
import java.util.UUID

class AnalysisArtifactRepositorySpec extends FunSuite with DatabaseTestSupport:

  val biosampleRepo = BiosampleRepository()
  val seqRunRepo = SequenceRunRepository()
  val alignmentRepo = AlignmentRepository()
  val artifactRepo = AnalysisArtifactRepository()

  private def createTestAlignment(tx: Transactor): AlignmentEntity =
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))
      val seqRun = seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS"))
      alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh38", "BWA-MEM2"))
    }.getOrElse(throw new RuntimeException("Failed to create test data"))

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      val entity = AnalysisArtifactEntity.create(
        alignmentId = alignment.id,
        artifactType = ArtifactType.WgsMetrics,
        cachePath = "cache/wgs_metrics_123.txt",
        generatorVersion = Some("1.0.0"),
        dependsOnSourceChecksum = Some("abc123")
      )

      val saved = artifactRepo.insert(entity)
      val found = artifactRepo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.alignmentId, alignment.id)
      assertEquals(found.get.artifactType, ArtifactType.WgsMetrics)
      assertEquals(found.get.cachePath, "cache/wgs_metrics_123.txt")
      assertEquals(found.get.status, ArtifactStatus.InProgress)
      assertEquals(found.get.generatorVersion, Some("1.0.0"))
      assertEquals(found.get.dependsOnSourceChecksum, Some("abc123"))
    }
  }

  testTransactor.test("findByAlignment returns all artifacts for alignment") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/wgs.txt"))
      artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.CallableLoci, "cache/callable.bed"))
      artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.CoverageSummary, "cache/coverage.txt"))

      val artifacts = artifactRepo.findByAlignment(alignment.id)
      assertEquals(artifacts.size, 3)
    }
  }

  testTransactor.test("findByAlignmentAndType returns specific artifact") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/wgs.txt"))
      artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.CallableLoci, "cache/callable.bed"))

      val wgs = artifactRepo.findByAlignmentAndType(alignment.id, ArtifactType.WgsMetrics)
      assert(wgs.isDefined)
      assertEquals(wgs.get.artifactType, ArtifactType.WgsMetrics)

      val callable = artifactRepo.findByAlignmentAndType(alignment.id, ArtifactType.CallableLoci)
      assert(callable.isDefined)
      assertEquals(callable.get.artifactType, ArtifactType.CallableLoci)

      val missing = artifactRepo.findByAlignmentAndType(alignment.id, ArtifactType.HaplogroupVcf)
      assertEquals(missing, None)
    }
  }

  testTransactor.test("markAvailable updates status and file info") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      val entity = artifactRepo.insert(
        AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/wgs.txt")
      )
      assertEquals(entity.status, ArtifactStatus.InProgress)

      artifactRepo.markAvailable(entity.id, fileSize = 12345L, fileChecksum = "sha256:abc", fileFormat = Some("TSV"))

      val found = artifactRepo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.status, ArtifactStatus.Available)
      assertEquals(found.get.fileSize, Some(12345L))
      assertEquals(found.get.fileChecksum, Some("sha256:abc"))
      assertEquals(found.get.fileFormat, Some("TSV"))
    }
  }

  testTransactor.test("markStale updates status with reason") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      val entity = artifactRepo.insert(
        AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/wgs.txt")
      )
      artifactRepo.markAvailable(entity.id, 1000L, "checksum")

      artifactRepo.markStale(entity.id, "Source file changed")

      val found = artifactRepo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.status, ArtifactStatus.Stale)
      assertEquals(found.get.staleReason, Some("Source file changed"))
    }
  }

  testTransactor.test("markError updates status with error reason") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      val entity = artifactRepo.insert(
        AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/wgs.txt")
      )

      artifactRepo.markError(entity.id, "Out of memory")

      val found = artifactRepo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.status, ArtifactStatus.Error)
      assertEquals(found.get.staleReason, Some("Out of memory"))
    }
  }

  testTransactor.test("findByStatus returns artifacts with matching status") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      val a1 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/1.txt"))
      val a2 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.CallableLoci, "cache/2.txt"))
      val a3 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.CoverageSummary, "cache/3.txt"))

      // Mark one as available, one as error
      artifactRepo.markAvailable(a1.id, 1000L, "checksum1")
      artifactRepo.markError(a2.id, "Failed")

      val inProgress = artifactRepo.findByStatus(ArtifactStatus.InProgress)
      assertEquals(inProgress.size, 1)

      val available = artifactRepo.findByStatus(ArtifactStatus.Available)
      assertEquals(available.size, 1)

      val error = artifactRepo.findByStatus(ArtifactStatus.Error)
      assertEquals(error.size, 1)
    }
  }

  testTransactor.test("findStale returns stale artifacts") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      val a1 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/1.txt"))
      val a2 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.CallableLoci, "cache/2.txt"))

      artifactRepo.markAvailable(a1.id, 1000L, "checksum1")
      artifactRepo.markAvailable(a2.id, 2000L, "checksum2")

      artifactRepo.markStale(a1.id, "Source changed")

      val stale = artifactRepo.findStale()
      assertEquals(stale.size, 1)
      assertEquals(stale.head.id, a1.id)
    }
  }

  testTransactor.test("markStaleBySourceChecksum marks all matching artifacts") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      val a1 = artifactRepo.insert(
        AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/1.txt", dependsOnSourceChecksum = Some("checksum_A"))
      )
      val a2 = artifactRepo.insert(
        AnalysisArtifactEntity.create(alignment.id, ArtifactType.CallableLoci, "cache/2.txt", dependsOnSourceChecksum = Some("checksum_A"))
      )
      val a3 = artifactRepo.insert(
        AnalysisArtifactEntity.create(alignment.id, ArtifactType.CoverageSummary, "cache/3.txt", dependsOnSourceChecksum = Some("checksum_B"))
      )

      // Mark all as available
      artifactRepo.markAvailable(a1.id, 1000L, "c1")
      artifactRepo.markAvailable(a2.id, 2000L, "c2")
      artifactRepo.markAvailable(a3.id, 3000L, "c3")

      // Mark stale by checksum
      val count = artifactRepo.markStaleBySourceChecksum("checksum_A", "Source file updated")
      assertEquals(count, 2)

      val stale = artifactRepo.findStale()
      assertEquals(stale.size, 2)
      assert(stale.forall(_.dependsOnSourceChecksum.contains("checksum_A")))
    }
  }

  testTransactor.test("countByStatus returns correct counts") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      val a1 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/1.txt"))
      val a2 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.CallableLoci, "cache/2.txt"))
      val a3 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.CoverageSummary, "cache/3.txt"))

      artifactRepo.markAvailable(a1.id, 1000L, "c1")
      artifactRepo.markAvailable(a2.id, 2000L, "c2")
      artifactRepo.markError(a3.id, "Failed")

      val counts = artifactRepo.countByStatus()
      assertEquals(counts.getOrElse(ArtifactStatus.Available, 0L), 2L)
      assertEquals(counts.getOrElse(ArtifactStatus.Error, 0L), 1L)
    }
  }

  testTransactor.test("totalCacheSize sums available artifact sizes") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      val a1 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/1.txt"))
      val a2 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.CallableLoci, "cache/2.txt"))
      val a3 = artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.CoverageSummary, "cache/3.txt"))

      artifactRepo.markAvailable(a1.id, 1000L, "c1")
      artifactRepo.markAvailable(a2.id, 2000L, "c2")
      // a3 stays in progress, shouldn't be counted

      val totalSize = artifactRepo.totalCacheSize()
      assertEquals(totalSize, 3000L)
    }
  }

  testTransactor.test("delete removes artifact") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    tx.readWrite {
      val entity = artifactRepo.insert(
        AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/wgs.txt")
      )

      assert(artifactRepo.findById(entity.id).isDefined)

      val deleted = artifactRepo.delete(entity.id)
      assert(deleted)
      assertEquals(artifactRepo.findById(entity.id), None)
    }
  }

  testTransactor.test("unique constraint on alignment_id and artifact_type") { case (db, tx) =>
    val alignment = createTestAlignment(tx)

    val result = tx.readWrite {
      artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/1.txt"))
      artifactRepo.insert(AnalysisArtifactEntity.create(alignment.id, ArtifactType.WgsMetrics, "cache/2.txt")) // Duplicate
    }

    assert(result.isLeft, "Should fail on duplicate alignment_id + artifact_type")
  }
