package com.decodingus.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.*
import munit.FunSuite

import java.util.UUID
import java.util.concurrent.{CountDownLatch, Executors, TimeUnit}
import scala.concurrent.{Await, ExecutionContext, Future}
import scala.concurrent.duration.*

/**
 * Integration tests for performance and concurrency requirements.
 * Based on design doc requirements:
 * - Large dataset performance (1000+ biosamples)
 * - Concurrent access simulation
 * - Query latency targets: <100ms for common operations
 * - Bulk import: 1000 samples in <10 seconds
 */
class PerformanceIntegrationSpec extends FunSuite with DatabaseTestSupport:

  val biosampleRepo = BiosampleRepository()
  val projectRepo = ProjectRepository()
  val sequenceRunRepo = SequenceRunRepository()
  val alignmentRepo = AlignmentRepository()
  val syncQueueRepo = SyncQueueRepository()

  // Performance test: bulk insert 1000 biosamples
  testTransactor.test("bulk insert 1000 biosamples completes in under 10 seconds") { case (db, tx) =>
    val startTime = System.currentTimeMillis()

    tx.readWrite {
      (1 to 1000).foreach { i =>
        val entity = BiosampleEntity.create(
          sampleAccession = f"PERF-$i%04d",
          donorIdentifier = f"DONOR-$i%04d"
        )
        biosampleRepo.insert(entity)
      }
    }

    val elapsed = System.currentTimeMillis() - startTime
    assert(elapsed < 10000, s"Bulk insert took ${elapsed}ms, expected <10000ms")

    // Verify all inserted
    tx.readOnly {
      val count = biosampleRepo.count()
      assertEquals(count, 1000L)
    }
  }

  // Performance test: query latency for large dataset
  testTransactor.test("findAll query on 1000 biosamples completes in under 100ms") { case (db, tx) =>
    // Setup: insert 1000 biosamples
    tx.readWrite {
      (1 to 1000).foreach { i =>
        biosampleRepo.insert(BiosampleEntity.create(
          sampleAccession = f"QUERY-$i%04d",
          donorIdentifier = f"DONOR-$i%04d"
        ))
      }
    }

    // Measure query time
    val startTime = System.currentTimeMillis()
    tx.readOnly {
      val results = biosampleRepo.findAll()
      val elapsed = System.currentTimeMillis() - startTime
      assertEquals(results.size, 1000)
      assert(elapsed < 100, s"Query took ${elapsed}ms, expected <100ms")
    }
  }

  // Performance test: findById with large dataset
  testTransactor.test("findById query with 1000 biosamples completes in under 10ms") { case (db, tx) =>
    // Setup: insert 1000 biosamples, capture one ID
    var targetId: UUID = null
    tx.readWrite {
      (1 to 1000).foreach { i =>
        val entity = biosampleRepo.insert(BiosampleEntity.create(
          sampleAccession = f"FIND-$i%04d",
          donorIdentifier = f"DONOR-$i%04d"
        ))
        if (i == 500) targetId = entity.id
      }
    }

    // Measure query time
    val startTime = System.currentTimeMillis()
    tx.readOnly {
      val result = biosampleRepo.findById(targetId)
      val elapsed = System.currentTimeMillis() - startTime
      assert(result.isDefined)
      assert(elapsed < 10, s"FindById took ${elapsed}ms, expected <10ms")
    }
  }

  // Performance test: findByAccession with large dataset
  testTransactor.test("findByAccession query with 1000 biosamples completes in under 10ms") { case (db, tx) =>
    // Setup: insert 1000 biosamples
    tx.readWrite {
      (1 to 1000).foreach { i =>
        biosampleRepo.insert(BiosampleEntity.create(
          sampleAccession = f"ACC-$i%04d",
          donorIdentifier = f"DONOR-$i%04d"
        ))
      }
    }

    // Measure query time
    val startTime = System.currentTimeMillis()
    tx.readOnly {
      val result = biosampleRepo.findByAccession("ACC-0500")
      val elapsed = System.currentTimeMillis() - startTime
      assert(result.isDefined)
      assert(elapsed < 10, s"FindByAccession took ${elapsed}ms, expected <10ms")
    }
  }

  // Performance test: project with many members
  testTransactor.test("project with 500 members queries in under 100ms") { case (db, tx) =>
    tx.readWrite {
      // Create project
      val project = projectRepo.insert(ProjectEntity.create(
        projectName = "Large Project",
        administratorDid = "did:plc:admin"
      ))

      // Add 500 biosamples as members
      (1 to 500).foreach { i =>
        val biosample = biosampleRepo.insert(BiosampleEntity.create(
          sampleAccession = f"MEMBER-$i%04d",
          donorIdentifier = f"DONOR-$i%04d"
        ))
        projectRepo.addMember(project.id, biosample.id)
      }

      // Measure query time
      val startTime = System.currentTimeMillis()
      val memberIds = projectRepo.getMemberIds(project.id)
      val elapsed = System.currentTimeMillis() - startTime

      assertEquals(memberIds.size, 500)
      assert(elapsed < 100, s"GetMemberIds took ${elapsed}ms, expected <100ms")
    }
  }

  // Performance test: cascade delete performance
  testTransactor.test("cascade delete of biosample with sequence runs and alignments completes in under 100ms") { case (db, tx) =>
    tx.readWrite {
      // Create biosample with multiple sequence runs and alignments
      val biosample = biosampleRepo.insert(BiosampleEntity.create(
        sampleAccession = "CASCADE-TEST",
        donorIdentifier = "DONOR-CASCADE"
      ))

      // Add 10 sequence runs, each with 5 alignments
      (1 to 10).foreach { i =>
        val seqRun = sequenceRunRepo.insert(SequenceRunEntity.create(
          biosampleId = biosample.id,
          platform = "ILLUMINA",
          testType = "WGS"
        ))

        (1 to 5).foreach { j =>
          alignmentRepo.insert(AlignmentEntity.create(
            sequenceRunId = seqRun.id,
            referenceBuild = "GRCh38",
            aligner = "BWA-MEM2"
          ))
        }
      }

      // Verify setup
      assertEquals(sequenceRunRepo.findByBiosample(biosample.id).size, 10)

      // Measure cascade delete time
      val startTime = System.currentTimeMillis()
      biosampleRepo.delete(biosample.id)
      val elapsed = System.currentTimeMillis() - startTime

      // Verify cascade
      assertEquals(sequenceRunRepo.findByBiosample(biosample.id).size, 0)
      assert(elapsed < 100, s"Cascade delete took ${elapsed}ms, expected <100ms")
    }
  }

  // Concurrency test: concurrent reads
  testTransactor.test("concurrent reads do not block each other") { case (db, tx) =>
    // Setup: insert test data
    tx.readWrite {
      (1 to 100).foreach { i =>
        biosampleRepo.insert(BiosampleEntity.create(
          sampleAccession = f"CONCURRENT-$i%04d",
          donorIdentifier = f"DONOR-$i%04d"
        ))
      }
    }

    // Run concurrent reads with fresh transactors
    val executor = Executors.newFixedThreadPool(10)
    implicit val ec: ExecutionContext = ExecutionContext.fromExecutor(executor)

    val startTime = System.currentTimeMillis()
    val futures = (1 to 10).map { _ =>
      Future {
        tx.readOnly { biosampleRepo.findAll().size }
      }
    }

    val results = Await.result(Future.sequence(futures), 5.seconds)
    val elapsed = System.currentTimeMillis() - startTime

    executor.shutdown()

    // All should return same results (handle Either)
    results.foreach {
      case Right(count) => assertEquals(count, 100)
      case Left(err) => fail(s"Read failed: $err")
    }
    // Concurrent reads should complete relatively quickly
    assert(elapsed < 2000, s"Concurrent reads took ${elapsed}ms")
  }

  // Sync queue performance test
  testTransactor.test("sync queue handles 100 items efficiently") { case (db, tx) =>
    tx.readWrite {
      // Create 100 biosamples and queue them for sync
      (1 to 100).foreach { i =>
        val biosample = biosampleRepo.insert(BiosampleEntity.create(
          sampleAccession = f"SYNC-$i%04d",
          donorIdentifier = f"DONOR-$i%04d"
        ))

        syncQueueRepo.enqueue(
          entityType = SyncEntityType.Biosample,
          entityId = biosample.id,
          operation = SyncOperation.Create
        )
      }
    }

    // Measure batch fetch time
    val startTime = System.currentTimeMillis()
    tx.readOnly {
      val batch = syncQueueRepo.findPendingBatch(100)
      val elapsed = System.currentTimeMillis() - startTime
      assertEquals(batch.size, 100)
      assert(elapsed < 100, s"Batch fetch took ${elapsed}ms, expected <100ms")
    }
  }

  // Memory efficiency test: large result set handling
  testTransactor.test("handles large result sets without excessive memory") { case (db, tx) =>
    // Insert 2000 biosamples
    tx.readWrite {
      (1 to 2000).foreach { i =>
        biosampleRepo.insert(BiosampleEntity.create(
          sampleAccession = f"MEMORY-$i%05d",
          donorIdentifier = f"DONOR-$i%05d"
        ))
      }
    }

    // Force GC before measurement
    System.gc()
    val memoryBefore = Runtime.getRuntime.totalMemory() - Runtime.getRuntime.freeMemory()

    // Load all biosamples
    tx.readOnly {
      val results = biosampleRepo.findAll()
      assertEquals(results.size, 2000)
    }

    val memoryAfter = Runtime.getRuntime.totalMemory() - Runtime.getRuntime.freeMemory()
    val memoryUsed = memoryAfter - memoryBefore

    // Should use less than 50MB for 2000 biosamples
    assert(memoryUsed < 50 * 1024 * 1024, s"Memory usage ${memoryUsed / 1024 / 1024}MB exceeds 50MB limit")
  }
