package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import munit.FunSuite
import java.util.UUID

class SourceFileRepositorySpec extends FunSuite with DatabaseTestSupport:

  val sourceFileRepo = SourceFileRepository()

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val entity = SourceFileEntity.create(
        filePath = "/data/sample.bam",
        fileChecksum = "sha256:abc123",
        fileSize = Some(1000000L),
        fileFormat = Some(SourceFileFormat.Bam)
      )

      val saved = sourceFileRepo.insert(entity)
      val found = sourceFileRepo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.filePath, Some("/data/sample.bam"))
      assertEquals(found.get.fileChecksum, "sha256:abc123")
      assertEquals(found.get.fileSize, Some(1000000L))
      assertEquals(found.get.fileFormat, Some(SourceFileFormat.Bam))
      assert(found.get.isAccessible)
      assert(!found.get.hasBeenAnalyzed)
    }
  }

  testTransactor.test("findByChecksum finds by stable identifier") { case (db, tx) =>
    tx.readWrite {
      sourceFileRepo.insert(SourceFileEntity.create("/path1/sample.bam", "checksum_A"))
      sourceFileRepo.insert(SourceFileEntity.create("/path2/sample.bam", "checksum_B"))

      val found = sourceFileRepo.findByChecksum("checksum_A")
      assert(found.isDefined)
      assertEquals(found.get.filePath, Some("/path1/sample.bam"))

      val notFound = sourceFileRepo.findByChecksum("nonexistent")
      assertEquals(notFound, None)
    }
  }

  testTransactor.test("findByPath finds by file path") { case (db, tx) =>
    tx.readWrite {
      sourceFileRepo.insert(SourceFileEntity.create("/data/sample1.bam", "checksum_A"))
      sourceFileRepo.insert(SourceFileEntity.create("/data/sample2.bam", "checksum_B"))

      val found = sourceFileRepo.findByPath("/data/sample1.bam")
      assert(found.isDefined)
      assertEquals(found.get.fileChecksum, "checksum_A")
    }
  }

  testTransactor.test("existsByChecksum checks existence") { case (db, tx) =>
    tx.readWrite {
      sourceFileRepo.insert(SourceFileEntity.create("/data/sample.bam", "existing_checksum"))

      assert(sourceFileRepo.existsByChecksum("existing_checksum"))
      assert(!sourceFileRepo.existsByChecksum("nonexistent_checksum"))
    }
  }

  testTransactor.test("markAccessible updates accessibility") { case (db, tx) =>
    tx.readWrite {
      val entity = sourceFileRepo.insert(SourceFileEntity.create("/data/sample.bam", "checksum"))
      // Make inaccessible first
      sourceFileRepo.markInaccessible(entity.id)

      val inaccessible = sourceFileRepo.findById(entity.id)
      assert(!inaccessible.get.isAccessible)

      sourceFileRepo.markAccessible(entity.id)

      val accessible = sourceFileRepo.findById(entity.id)
      assert(accessible.get.isAccessible)
      assert(accessible.get.lastVerifiedAt.isDefined)
    }
  }

  testTransactor.test("markInaccessible updates accessibility") { case (db, tx) =>
    tx.readWrite {
      val entity = sourceFileRepo.insert(SourceFileEntity.create("/data/sample.bam", "checksum"))
      assert(entity.isAccessible)

      sourceFileRepo.markInaccessible(entity.id)

      val found = sourceFileRepo.findById(entity.id)
      assert(!found.get.isAccessible)
    }
  }

  testTransactor.test("markAnalyzed updates analysis status") { case (db, tx) =>
    tx.readWrite {
      val entity = sourceFileRepo.insert(SourceFileEntity.create("/data/sample.bam", "checksum"))
      assert(!entity.hasBeenAnalyzed)
      assertEquals(entity.analysisCompletedAt, None)

      sourceFileRepo.markAnalyzed(entity.id)

      val found = sourceFileRepo.findById(entity.id)
      assert(found.get.hasBeenAnalyzed)
      assert(found.get.analysisCompletedAt.isDefined)
    }
  }

  testTransactor.test("updatePath changes file path") { case (db, tx) =>
    tx.readWrite {
      val entity = sourceFileRepo.insert(SourceFileEntity.create("/old/path/sample.bam", "checksum"))

      sourceFileRepo.updatePath(entity.id, "/new/path/sample.bam")

      val found = sourceFileRepo.findById(entity.id)
      assertEquals(found.get.filePath, Some("/new/path/sample.bam"))
      assert(found.get.isAccessible)
    }
  }

  testTransactor.test("linkToAlignment links source file to alignment") { case (db, tx) =>
    tx.readWrite {
      val entity = sourceFileRepo.insert(SourceFileEntity.create("/data/sample.bam", "checksum"))
      assertEquals(entity.alignmentId, None)

      val alignmentId = UUID.randomUUID()
      sourceFileRepo.linkToAlignment(entity.id, alignmentId)

      val found = sourceFileRepo.findById(entity.id)
      assertEquals(found.get.alignmentId, Some(alignmentId))
    }
  }

  testTransactor.test("findAccessible returns only accessible files") { case (db, tx) =>
    tx.readWrite {
      val f1 = sourceFileRepo.insert(SourceFileEntity.create("/data/1.bam", "checksum1"))
      val f2 = sourceFileRepo.insert(SourceFileEntity.create("/data/2.bam", "checksum2"))
      val f3 = sourceFileRepo.insert(SourceFileEntity.create("/data/3.bam", "checksum3"))

      sourceFileRepo.markInaccessible(f2.id)

      val accessible = sourceFileRepo.findAccessible()
      assertEquals(accessible.size, 2)
      assert(accessible.forall(_.isAccessible))
    }
  }

  testTransactor.test("findInaccessible returns only inaccessible files") { case (db, tx) =>
    tx.readWrite {
      val f1 = sourceFileRepo.insert(SourceFileEntity.create("/data/1.bam", "checksum1"))
      val f2 = sourceFileRepo.insert(SourceFileEntity.create("/data/2.bam", "checksum2"))

      sourceFileRepo.markInaccessible(f1.id)

      val inaccessible = sourceFileRepo.findInaccessible()
      assertEquals(inaccessible.size, 1)
      assertEquals(inaccessible.head.id, f1.id)
    }
  }

  testTransactor.test("findNotAnalyzed returns accessible files not yet analyzed") { case (db, tx) =>
    tx.readWrite {
      val f1 = sourceFileRepo.insert(SourceFileEntity.create("/data/1.bam", "checksum1"))
      val f2 = sourceFileRepo.insert(SourceFileEntity.create("/data/2.bam", "checksum2"))
      val f3 = sourceFileRepo.insert(SourceFileEntity.create("/data/3.bam", "checksum3"))

      sourceFileRepo.markAnalyzed(f1.id)
      sourceFileRepo.markInaccessible(f3.id)

      val notAnalyzed = sourceFileRepo.findNotAnalyzed()
      assertEquals(notAnalyzed.size, 1)
      assertEquals(notAnalyzed.head.id, f2.id)
    }
  }

  testTransactor.test("upsertByChecksum creates new if not exists") { case (db, tx) =>
    tx.readWrite {
      val entity = sourceFileRepo.upsertByChecksum(
        filePath = "/data/new.bam",
        fileChecksum = "new_checksum",
        fileSize = Some(5000L),
        fileFormat = Some(SourceFileFormat.Bam)
      )

      assertEquals(entity.filePath, Some("/data/new.bam"))
      assertEquals(entity.fileChecksum, "new_checksum")

      val found = sourceFileRepo.findByChecksum("new_checksum")
      assert(found.isDefined)
    }
  }

  testTransactor.test("upsertByChecksum updates path if exists") { case (db, tx) =>
    tx.readWrite {
      // First insert
      sourceFileRepo.insert(SourceFileEntity.create("/old/path.bam", "existing_checksum"))

      // Upsert with different path
      val entity = sourceFileRepo.upsertByChecksum(
        filePath = "/new/path.bam",
        fileChecksum = "existing_checksum",
        fileSize = None,
        fileFormat = None
      )

      assertEquals(entity.filePath, Some("/new/path.bam"))

      // Should not create duplicate
      val all = sourceFileRepo.findAll()
      assertEquals(all.size, 1)
    }
  }

  testTransactor.test("countByAccessibility returns correct counts") { case (db, tx) =>
    tx.readWrite {
      val f1 = sourceFileRepo.insert(SourceFileEntity.create("/data/1.bam", "checksum1"))
      val f2 = sourceFileRepo.insert(SourceFileEntity.create("/data/2.bam", "checksum2"))
      val f3 = sourceFileRepo.insert(SourceFileEntity.create("/data/3.bam", "checksum3"))

      sourceFileRepo.markInaccessible(f3.id)

      val (accessible, inaccessible) = sourceFileRepo.countByAccessibility()
      assertEquals(accessible, 2L)
      assertEquals(inaccessible, 1L)
    }
  }

  testTransactor.test("totalFileSize sums accessible file sizes") { case (db, tx) =>
    tx.readWrite {
      val f1 = sourceFileRepo.insert(SourceFileEntity.create("/data/1.bam", "checksum1", fileSize = Some(1000L)))
      val f2 = sourceFileRepo.insert(SourceFileEntity.create("/data/2.bam", "checksum2", fileSize = Some(2000L)))
      val f3 = sourceFileRepo.insert(SourceFileEntity.create("/data/3.bam", "checksum3", fileSize = Some(3000L)))

      sourceFileRepo.markInaccessible(f3.id) // Shouldn't be counted

      val totalSize = sourceFileRepo.totalFileSize()
      assertEquals(totalSize, 3000L)
    }
  }

  testTransactor.test("delete removes source file") { case (db, tx) =>
    tx.readWrite {
      val entity = sourceFileRepo.insert(SourceFileEntity.create("/data/sample.bam", "checksum"))

      assert(sourceFileRepo.findById(entity.id).isDefined)

      val deleted = sourceFileRepo.delete(entity.id)
      assert(deleted)
      assertEquals(sourceFileRepo.findById(entity.id), None)
    }
  }

  testTransactor.test("unique constraint on file_checksum") { case (db, tx) =>
    val result = tx.readWrite {
      sourceFileRepo.insert(SourceFileEntity.create("/path1/file.bam", "same_checksum"))
      sourceFileRepo.insert(SourceFileEntity.create("/path2/file.bam", "same_checksum")) // Duplicate
    }

    assert(result.isLeft, "Should fail on duplicate checksum")
  }
