package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.workspace.model.{AlignmentMetrics, FileInfo}
import munit.FunSuite
import java.util.UUID

class AlignmentRepositorySpec extends FunSuite with DatabaseTestSupport:

  val biosampleRepo = BiosampleRepository()
  val seqRunRepo = SequenceRunRepository()
  val alignmentRepo = AlignmentRepository()

  private def createTestSequenceRun(tx: Transactor): (BiosampleEntity, SequenceRunEntity) =
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))
      val seqRun = seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS"))
      (biosample, seqRun)
    }.getOrElse(throw new RuntimeException("Failed to create test data"))

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    val (_, seqRun) = createTestSequenceRun(tx)

    tx.readWrite {
      val entity = AlignmentEntity.create(
        sequenceRunId = seqRun.id,
        referenceBuild = "GRCh38",
        aligner = "BWA-MEM2",
        variantCaller = Some("GATK HaplotypeCaller")
      )

      val saved = alignmentRepo.insert(entity)
      val found = alignmentRepo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.sequenceRunId, seqRun.id)
      assertEquals(found.get.referenceBuild, "GRCh38")
      assertEquals(found.get.aligner, "BWA-MEM2")
      assertEquals(found.get.variantCaller, Some("GATK HaplotypeCaller"))
      assertEquals(found.get.meta.syncStatus, SyncStatus.Local)
    }
  }

  testTransactor.test("findAll returns all alignments") { case (db, tx) =>
    val (_, seqRun) = createTestSequenceRun(tx)

    tx.readWrite {
      alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh38", "BWA-MEM2"))
      alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh37", "Bowtie2"))

      val all = alignmentRepo.findAll()
      assertEquals(all.size, 2)
    }
  }

  testTransactor.test("update modifies entity correctly") { case (db, tx) =>
    val (_, seqRun) = createTestSequenceRun(tx)

    tx.readWrite {
      val entity = alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh38", "BWA-MEM2"))

      val updated = entity.copy(
        variantCaller = Some("DeepVariant"),
        aligner = "Minimap2"
      )
      alignmentRepo.update(updated)

      val found = alignmentRepo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.aligner, "Minimap2")
      assertEquals(found.get.variantCaller, Some("DeepVariant"))
      assertEquals(found.get.meta.version, 2)
    }
  }

  testTransactor.test("delete removes alignment") { case (db, tx) =>
    val (_, seqRun) = createTestSequenceRun(tx)

    tx.readWrite {
      val entity = alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh38", "BWA-MEM2"))

      assert(alignmentRepo.exists(entity.id))

      val deleted = alignmentRepo.delete(entity.id)
      assert(deleted)
      assert(!alignmentRepo.exists(entity.id))
    }
  }

  testTransactor.test("findBySequenceRun returns alignments for sequence run") { case (db, tx) =>
    tx.readWrite {
      val biosample = biosampleRepo.insert(BiosampleEntity.create("BS001", "D1"))
      val seqRun1 = seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "ILLUMINA", "WGS"))
      val seqRun2 = seqRunRepo.insert(SequenceRunEntity.create(biosample.id, "PACBIO", "WGS"))

      alignmentRepo.insert(AlignmentEntity.create(seqRun1.id, "GRCh38", "BWA-MEM2"))
      alignmentRepo.insert(AlignmentEntity.create(seqRun1.id, "GRCh37", "BWA-MEM2"))
      alignmentRepo.insert(AlignmentEntity.create(seqRun2.id, "GRCh38", "Minimap2"))

      val seqRun1Alignments = alignmentRepo.findBySequenceRun(seqRun1.id)
      assertEquals(seqRun1Alignments.size, 2)

      val seqRun2Alignments = alignmentRepo.findBySequenceRun(seqRun2.id)
      assertEquals(seqRun2Alignments.size, 1)
    }
  }

  testTransactor.test("findByReferenceBuild filters by reference") { case (db, tx) =>
    val (_, seqRun) = createTestSequenceRun(tx)

    tx.readWrite {
      alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh38", "BWA-MEM2"))
      alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh38", "Bowtie2"))
      alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh37", "BWA-MEM2"))

      val grch38 = alignmentRepo.findByReferenceBuild("GRCh38")
      assertEquals(grch38.size, 2)

      val grch37 = alignmentRepo.findByReferenceBuild("GRCh37")
      assertEquals(grch37.size, 1)
    }
  }

  testTransactor.test("updateMetrics updates alignment metrics") { case (db, tx) =>
    val (_, seqRun) = createTestSequenceRun(tx)

    tx.readWrite {
      val entity = alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh38", "BWA-MEM2"))

      val metrics = AlignmentMetrics(
        meanCoverage = Some(30.5),
        sdCoverage = Some(5.2),
        medianCoverage = Some(30.0),
        pct10x = Some(95.0),
        pct20x = Some(85.0),
        pctExcDupe = Some(10.5)
      )

      alignmentRepo.updateMetrics(entity.id, metrics)

      val found = alignmentRepo.findById(entity.id)
      assert(found.isDefined)
      assert(found.get.metrics.isDefined)
      assertEquals(found.get.metrics.get.meanCoverage, Some(30.5))
      assertEquals(found.get.metrics.get.pct10x, Some(95.0))
    }
  }

  testTransactor.test("update preserves file list in entity") { case (db, tx) =>
    val (_, seqRun) = createTestSequenceRun(tx)

    tx.readWrite {
      val files = List(
        FileInfo(
          fileName = "sample.bam",
          fileFormat = "BAM",
          checksum = Some("abc123"),
          fileSizeBytes = Some(1000000L),
          location = Some("/data/sample.bam")
        ),
        FileInfo(
          fileName = "sample.bam.bai",
          fileFormat = "BAI",
          checksum = Some("def456"),
          fileSizeBytes = Some(10000L),
          location = Some("/data/sample.bam.bai")
        )
      )

      val entity = AlignmentEntity.create(seqRun.id, "GRCh38", "BWA-MEM2").copy(files = files)
      alignmentRepo.insert(entity)

      val found = alignmentRepo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.files.size, 2)
      assertEquals(found.get.files.head.fileName, "sample.bam")
    }
  }

  testTransactor.test("findPendingSync returns Local and Modified alignments") { case (db, tx) =>
    val (_, seqRun) = createTestSequenceRun(tx)

    tx.readWrite {
      val local = alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh38", "BWA-MEM2"))
      val synced = alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh37", "BWA-MEM2"))

      alignmentRepo.markSynced(synced.id, "at://test/1", "cid1")

      val pending = alignmentRepo.findPendingSync()
      assertEquals(pending.size, 1)
      assertEquals(pending.head.id, local.id)
    }
  }

  testTransactor.test("foreign key to sequence_run is enforced") { case (db, tx) =>
    val result = tx.readWrite {
      val fakeId = UUID.randomUUID()
      alignmentRepo.insert(AlignmentEntity.create(fakeId, "GRCh38", "BWA-MEM2"))
    }

    assert(result.isLeft, "Should fail with foreign key violation")
  }

  testTransactor.test("cascade delete removes alignments when sequence run deleted") { case (db, tx) =>
    val (biosample, seqRun) = createTestSequenceRun(tx)

    tx.readWrite {
      val alignment = alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "GRCh38", "BWA-MEM2"))

      assert(alignmentRepo.exists(alignment.id))

      seqRunRepo.delete(seqRun.id)

      assert(!alignmentRepo.exists(alignment.id), "Alignment should be cascade deleted")
    }
  }

  testTransactor.test("reference build constraint is enforced") { case (db, tx) =>
    val (_, seqRun) = createTestSequenceRun(tx)

    val result = tx.readWrite {
      alignmentRepo.insert(AlignmentEntity.create(seqRun.id, "InvalidBuild", "BWA-MEM2"))
    }

    assert(result.isLeft, "Should fail with invalid reference build")
  }
