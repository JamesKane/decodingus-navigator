package com.decodingus.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.*
import munit.FunSuite
import java.nio.file.{Files, Path, Paths}
import java.util.UUID

class CacheServiceSpec extends FunSuite with DatabaseTestSupport:

  // Use the actual cache directory for tests that need file existence
  private val CacheDir: Path = Paths.get(System.getProperty("user.home"), ".decodingus", "cache")

  override def beforeAll(): Unit =
    super.beforeAll()
    Files.createDirectories(CacheDir)

  private def createCacheService(tx: Transactor): H2CacheService =
    H2CacheService(
      transactor = tx,
      artifactRepo = AnalysisArtifactRepository(),
      sourceFileRepo = SourceFileRepository()
    )

  // Helper to create an alignment for testing
  private def createAlignment(tx: Transactor): UUID =
    val biosampleRepo = BiosampleRepository()
    val seqRunRepo = SequenceRunRepository()
    val alignmentRepo = AlignmentRepository()

    tx.readWrite {
      val biosample = BiosampleEntity.create("TEST-001", "DONOR-001")
      biosampleRepo.insert(biosample)

      val seqRun = SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS")
      seqRunRepo.insert(seqRun)

      val alignment = AlignmentEntity.create(seqRun.id, "GRCh38", "BWA-MEM2")
      alignmentRepo.insert(alignment)
      alignment.id
    }.getOrElse(throw new RuntimeException("Failed to create test alignment"))

  // ============================================
  // Source File Tests
  // ============================================

  testTransactor.test("registerSourceFile creates new source file") { case (db, tx) =>
    val service = createCacheService(tx)

    val result = service.registerSourceFile(
      filePath = "/test/sample.bam",
      fileChecksum = "abc123",
      fileSize = Some(1000L),
      fileFormat = Some(SourceFileFormat.Bam)
    )

    assert(result.isRight)
    result.foreach { sf =>
      assertEquals(sf.fileChecksum, "abc123")
      assertEquals(sf.filePath, Some("/test/sample.bam"))
      assertEquals(sf.fileFormat, Some(SourceFileFormat.Bam))
    }
  }

  testTransactor.test("registerSourceFile updates existing by checksum") { case (db, tx) =>
    val service = createCacheService(tx)

    // First registration
    service.registerSourceFile("/old/path.bam", "checksum123", Some(1000L), None)

    // Second registration with same checksum but different path
    val result = service.registerSourceFile("/new/path.bam", "checksum123", Some(1000L), None)

    assert(result.isRight)
    result.foreach { sf =>
      assertEquals(sf.filePath, Some("/new/path.bam"))
    }
  }

  testTransactor.test("getSourceFileByChecksum finds existing file") { case (db, tx) =>
    val service = createCacheService(tx)

    service.registerSourceFile("/test/file.bam", "unique-checksum", None, None)

    val result = service.getSourceFileByChecksum("unique-checksum")
    assert(result.isRight)
    assert(result.toOption.flatten.isDefined)
  }

  // ============================================
  // Artifact Management Tests
  // ============================================

  testTransactor.test("startArtifact creates artifact in progress") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    val result = service.startArtifact(
      alignmentId = alignmentId,
      artifactType = ArtifactType.WgsMetrics,
      cachePath = "wgs/metrics.txt",
      generatorVersion = Some("GATK-4.6"),
      dependsOnSourceChecksum = Some("source-checksum"),
      dependsOnReferenceBuild = Some("GRCh38")
    )

    assert(result.isRight)
    result.foreach { artifact =>
      assertEquals(artifact.artifactType, ArtifactType.WgsMetrics)
      assertEquals(artifact.status, ArtifactStatus.InProgress)
      assertEquals(artifact.dependsOnSourceChecksum, Some("source-checksum"))
      assertEquals(artifact.dependsOnReferenceBuild, Some("GRCh38"))
    }
  }

  testTransactor.test("completeArtifact marks artifact as available") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    val artifact = service.startArtifact(
      alignmentId, ArtifactType.WgsMetrics, "wgs/metrics.txt"
    ).toOption.get

    val result = service.completeArtifact(artifact.id, 500L, "file-checksum", Some("TXT"))

    assert(result.isRight)
    assertEquals(result, Right(true))

    // Verify status changed
    val retrieved = service.getArtifact(alignmentId, ArtifactType.WgsMetrics)
    assert(retrieved.toOption.flatten.exists(_.status == ArtifactStatus.Available))
  }

  testTransactor.test("isArtifactAvailable returns false for missing file") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    val artifact = service.startArtifact(alignmentId, ArtifactType.WgsMetrics, "nonexistent.txt").toOption.get
    service.completeArtifact(artifact.id, 100L, "checksum", None)

    // File doesn't exist, so should return false
    val result = service.isArtifactAvailable(alignmentId, ArtifactType.WgsMetrics)
    assertEquals(result, Right(false))
  }

  // ============================================
  // Cache Invalidation Tests
  // ============================================

  testTransactor.test("invalidateArtifactsForAlignment marks artifacts stale") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    // Create and complete an artifact
    val artifact = service.startArtifact(alignmentId, ArtifactType.WgsMetrics, "metrics.txt").toOption.get
    service.completeArtifact(artifact.id, 100L, "checksum", None)

    // Invalidate
    val result = service.invalidateArtifactsForAlignment(alignmentId, "Test invalidation")

    assertEquals(result, Right(1))

    // Verify it's stale
    val retrieved = service.getArtifact(alignmentId, ArtifactType.WgsMetrics)
    assert(retrieved.toOption.flatten.exists(_.status == ArtifactStatus.Stale))
  }

  testTransactor.test("invalidateBySourceChecksum marks matching artifacts stale") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)
    val sourceChecksum = "source-abc123"

    // Create artifact with source dependency
    val artifact = service.startArtifact(
      alignmentId, ArtifactType.CallableLoci, "callable.bed",
      dependsOnSourceChecksum = Some(sourceChecksum)
    ).toOption.get
    service.completeArtifact(artifact.id, 200L, "file-checksum", None)

    // Invalidate by source checksum
    val result = service.invalidateBySourceChecksum(sourceChecksum, "Source file changed")

    assertEquals(result, Right(1))

    // Verify it's stale
    val retrieved = service.getArtifact(alignmentId, ArtifactType.CallableLoci)
    assert(retrieved.toOption.flatten.exists(_.status == ArtifactStatus.Stale))
  }

  testTransactor.test("invalidateByReferenceBuild marks matching artifacts stale") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    // Create artifact with reference dependency
    val artifact = service.startArtifact(
      alignmentId, ArtifactType.HaplogroupVcf, "haplogroup.vcf",
      dependsOnReferenceBuild = Some("GRCh38")
    ).toOption.get
    service.completeArtifact(artifact.id, 300L, "vcf-checksum", Some("VCF"))

    // Invalidate by reference build
    val result = service.invalidateByReferenceBuild("GRCh38", "Reference updated")

    assertEquals(result, Right(1))

    // Verify it's stale
    val retrieved = service.getArtifact(alignmentId, ArtifactType.HaplogroupVcf)
    assert(retrieved.toOption.flatten.exists(_.status == ArtifactStatus.Stale))
  }

  testTransactor.test("validateArtifact returns NotFound for missing artifact") { case (db, tx) =>
    val service = createCacheService(tx)

    val result = service.validateArtifact(UUID.randomUUID())

    assertEquals(result, Right(ArtifactValidationResult.NotFound))
  }

  testTransactor.test("validateArtifact returns AlreadyStale for stale artifact") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    // Create and make stale
    val artifact = service.startArtifact(alignmentId, ArtifactType.WgsMetrics, "metrics.txt").toOption.get
    service.completeArtifact(artifact.id, 100L, "checksum", None)
    service.invalidateArtifactsForAlignment(alignmentId, "Manual invalidation")

    val result = service.validateArtifact(artifact.id)

    assertEquals(result, Right(ArtifactValidationResult.AlreadyStale))
  }

  testTransactor.test("validateArtifact returns FileNotFound for missing cache file") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    // Create artifact pointing to non-existent file
    val artifact = service.startArtifact(alignmentId, ArtifactType.WgsMetrics, "nonexistent.txt").toOption.get
    service.completeArtifact(artifact.id, 100L, "checksum", None)

    val result = service.validateArtifact(artifact.id)

    assertEquals(result, Right(ArtifactValidationResult.FileNotFound))
  }

  testTransactor.test("validateArtifact returns Valid for valid artifact with file") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    // Create actual file in cache dir
    val cachePath = s"test-valid-metrics-${UUID.randomUUID()}.txt"
    Files.writeString(CacheDir.resolve(cachePath), "test content")

    try
      val artifact = service.startArtifact(alignmentId, ArtifactType.WgsMetrics, cachePath).toOption.get
      service.completeArtifact(artifact.id, 100L, "checksum", None)

      val result = service.validateArtifact(artifact.id)

      assertEquals(result, Right(ArtifactValidationResult.Valid))
    finally
      // Cleanup test file
      Files.deleteIfExists(CacheDir.resolve(cachePath))
  }

  testTransactor.test("getStaleArtifacts returns stale artifacts") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    // Create and make stale
    val artifact = service.startArtifact(alignmentId, ArtifactType.WgsMetrics, "metrics.txt").toOption.get
    service.completeArtifact(artifact.id, 100L, "checksum", None)
    service.invalidateArtifactsForAlignment(alignmentId, "Test")

    val result = service.getStaleArtifacts()

    assert(result.isRight)
    assertEquals(result.toOption.get.size, 1)
  }

  testTransactor.test("cleanupMissingArtifacts marks missing files as deleted") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    // Create artifact with non-existent file
    val artifact = service.startArtifact(alignmentId, ArtifactType.WgsMetrics, "missing.txt").toOption.get
    service.completeArtifact(artifact.id, 100L, "checksum", None)

    val result = service.cleanupMissingArtifacts()

    assertEquals(result, Right(1))

    // Verify status changed to deleted
    val retrieved = service.getArtifact(alignmentId, ArtifactType.WgsMetrics)
    assert(retrieved.toOption.flatten.exists(_.status == ArtifactStatus.Deleted))
  }

  // ============================================
  // Statistics Tests
  // ============================================

  testTransactor.test("getCacheStats returns accurate counts") { case (db, tx) =>
    val service = createCacheService(tx)
    val alignmentId = createAlignment(tx)

    // Create some artifacts in different states
    val a1 = service.startArtifact(alignmentId, ArtifactType.WgsMetrics, "m1.txt").toOption.get
    service.completeArtifact(a1.id, 100L, "c1", None)

    val a2 = service.startArtifact(alignmentId, ArtifactType.CallableLoci, "c1.bed").toOption.get
    // Leave in progress

    val a3 = service.startArtifact(alignmentId, ArtifactType.CoverageSummary, "cov.txt").toOption.get
    service.completeArtifact(a3.id, 200L, "c3", None)
    service.invalidateArtifactsForAlignment(alignmentId, "Test") // Makes a1 and a3 stale

    val result = service.getCacheStats()

    assert(result.isRight)
    result.foreach { stats =>
      assertEquals(stats.totalArtifacts, 3)
      assertEquals(stats.inProgressArtifacts, 1)
      assertEquals(stats.staleArtifacts, 2)
    }
  }
