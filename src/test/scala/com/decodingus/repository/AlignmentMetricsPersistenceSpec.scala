package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.workspace.model.{AlignmentMetrics, ContigMetrics}
import munit.FunSuite
import java.util.UUID

class AlignmentMetricsPersistenceSpec extends FunSuite with DatabaseTestSupport:

  testTransactor.test("AlignmentRepository persists and retrieves metrics correctly") { case (db, tx) =>
    val repo = new AlignmentRepository()
    val seqRunRepo = new SequenceRunRepository()
    
    tx.readWrite {
      // Setup dependency: SequenceRun
      val seqRun = SequenceRunEntity.create(
        biosampleId = UUID.randomUUID(),
        platform = "ILLUMINA",
        testType = "WGS"
      )
      seqRunRepo.insert(seqRun)

      // 1. Create Alignment
      val alignment = AlignmentEntity.create(
        sequenceRunId = seqRun.id,
        referenceBuild = "GRCh38",
        aligner = "BWA-MEM2"
      )
      repo.insert(alignment)

      // 2. Create Metrics
      val metrics = AlignmentMetrics(
        genomeTerritory = Some(3000000000L),
        meanCoverage = Some(30.5),
        pctExcDupe = Some(0.05),
        // List field
        contigs = List(
          ContigMetrics("chr1", 1000, 10, 5, 2, 1),
          ContigMetrics("chr2", 2000, 20, 10, 4, 2)
        )
      )

      // 3. Update Metrics using updateMetrics
      val updateSuccess = repo.updateMetrics(alignment.id, metrics)
      assert(updateSuccess, "updateMetrics should return true")

      // 4. Retrieve and Verify
      val retrieved = repo.findById(alignment.id).get
      assert(retrieved.metrics.isDefined, "Metrics should be defined")
      
      val loadedMetrics = retrieved.metrics.get
      assertEquals(loadedMetrics.genomeTerritory, Some(3000000000L))
      assertEquals(loadedMetrics.meanCoverage, Some(30.5))
      assertEquals(loadedMetrics.contigs.size, 2)
      assertEquals(loadedMetrics.contigs.head.contigName, "chr1")

      // 5. Update using generic update method
      val updatedMetrics2 = metrics.copy(meanCoverage = Some(45.0))
      val entityToUpdate = retrieved.copy(metrics = Some(updatedMetrics2))
      repo.update(entityToUpdate)

      val retrieved2 = repo.findById(alignment.id).get
      assertEquals(retrieved2.metrics.get.meanCoverage, Some(45.0))
    }
  }
