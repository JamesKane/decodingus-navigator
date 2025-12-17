package com.decodingus.yprofile.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.{
  BiosampleRepository, BiosampleEntity, AlignmentRepository, AlignmentEntity,
  SequenceRunRepository, SequenceRunEntity, SyncStatus
}
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class YSnpPanelRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = YSnpPanelRepository()
  val biosampleRepo = BiosampleRepository()
  val alignmentRepo = AlignmentRepository()
  val seqRunRepo = SequenceRunRepository()

  def createTestBiosample(accession: String = s"TEST${UUID.randomUUID().toString.take(8)}")(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = accession,
      donorIdentifier = "DONOR001"
    ))

  def createTestAlignment(biosampleId: UUID)(using java.sql.Connection): AlignmentEntity =
    val seqRun = seqRunRepo.insert(SequenceRunEntity.create(
      biosampleId = biosampleId,
      platform = "ILLUMINA",
      testType = "WGS"
    ))
    alignmentRepo.insert(AlignmentEntity.create(
      sequenceRunId = seqRun.id,
      referenceBuild = "GRCh38",
      aligner = "BWA-MEM2"
    ))

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = YSnpPanelEntity.create(
        biosampleId = biosample.id,
        panelName = Some("Big Y-700"),
        provider = Some("FTDNA"),
        testDate = Some(LocalDateTime.now()),
        totalSnpsTested = Some(35000),
        derivedCount = Some(480),
        ancestralCount = Some(34000),
        noCallCount = Some(520),
        terminalHaplogroup = Some("R-BY140757"),
        confidence = Some(0.99),
        snpCalls = List(
          YSnpCall("M343", 2787994, None, "A", true, Some(YVariantType.SNP), None, Some(99.0)),
          YSnpCall("M269", 22739368, None, "C", true, Some(YVariantType.SNP), None, Some(98.5))
        ),
        privateVariants = List(
          YPrivateVariant(12345678, "G", "A", Some("FGC99999"), Some(95.0), Some(45))
        )
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined, "Should find inserted entity")
      assertEquals(found.get.panelName, Some("Big Y-700"))
      assertEquals(found.get.provider, Some("FTDNA"))
      assertEquals(found.get.totalSnpsTested, Some(35000))
      assertEquals(found.get.terminalHaplogroup, Some("R-BY140757"))
      assertEquals(found.get.snpCalls.size, 2)
      assertEquals(found.get.privateVariants.size, 1)
      assertEquals(found.get.meta.syncStatus, SyncStatus.Local)
    }
  }

  testTransactor.test("findByBiosample returns all panels for biosample") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val other = createTestBiosample()

      repo.insert(YSnpPanelEntity.create(biosample.id, panelName = Some("Big Y-700")))
      repo.insert(YSnpPanelEntity.create(biosample.id, panelName = Some("WGS Analysis")))
      repo.insert(YSnpPanelEntity.create(other.id, panelName = Some("Big Y-500")))

      val panels = repo.findByBiosample(biosample.id)
      assertEquals(panels.size, 2)
      assert(panels.forall(_.biosampleId == biosample.id))
    }
  }

  testTransactor.test("findByAlignment returns panels linked to alignment") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val alignment = createTestAlignment(biosample.id)

      val withAlignment = repo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        alignmentId = Some(alignment.id),
        panelName = Some("WGS-derived")
      ))
      repo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        alignmentId = None,
        panelName = Some("Big Y-700")
      ))

      val aligned = repo.findByAlignment(alignment.id)
      assertEquals(aligned.size, 1)
      assertEquals(aligned.head.id, withAlignment.id)
    }
  }

  testTransactor.test("findByProvider filters by provider") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()

      repo.insert(YSnpPanelEntity.create(biosample.id, provider = Some("FTDNA")))
      repo.insert(YSnpPanelEntity.create(biosample.id, provider = Some("FTDNA")))
      repo.insert(YSnpPanelEntity.create(biosample.id, provider = Some("YSEQ")))

      val ftdna = repo.findByProvider("FTDNA")
      assertEquals(ftdna.size, 2)

      val yseq = repo.findByProvider("YSEQ")
      assertEquals(yseq.size, 1)
    }
  }

  testTransactor.test("findByHaplogroup finds exact match") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()

      repo.insert(YSnpPanelEntity.create(biosample.id, terminalHaplogroup = Some("R-M269")))
      repo.insert(YSnpPanelEntity.create(biosample.id, terminalHaplogroup = Some("R-L21")))
      repo.insert(YSnpPanelEntity.create(biosample.id, terminalHaplogroup = Some("I1-M253")))

      val m269 = repo.findByHaplogroup("R-M269")
      assertEquals(m269.size, 1)
    }
  }

  testTransactor.test("findByHaplogroupBranch finds prefix matches") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()

      repo.insert(YSnpPanelEntity.create(biosample.id, terminalHaplogroup = Some("R-M269")))
      repo.insert(YSnpPanelEntity.create(biosample.id, terminalHaplogroup = Some("R-L21")))
      repo.insert(YSnpPanelEntity.create(biosample.id, terminalHaplogroup = Some("R-BY12345")))
      repo.insert(YSnpPanelEntity.create(biosample.id, terminalHaplogroup = Some("I1-M253")))

      val rBranch = repo.findByHaplogroupBranch("R-")
      assertEquals(rBranch.size, 3)

      val iBranch = repo.findByHaplogroupBranch("I1")
      assertEquals(iBranch.size, 1)
    }
  }

  testTransactor.test("stores and retrieves SNP calls correctly") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val snpCalls = List(
        YSnpCall("M343", 2787994, None, "A", true, Some(YVariantType.SNP), None, Some(99.0)),
        YSnpCall("M269", 22739368, None, "C", true, Some(YVariantType.SNP), None, Some(98.5)),
        YSnpCall("P312", 15579988, None, "C", true, None, None, Some(97.0)),
        YSnpCall("L21", 14722131, None, "A", true, None, None, Some(96.5))
      )

      val entity = YSnpPanelEntity.create(
        biosampleId = biosample.id,
        snpCalls = snpCalls
      )
      repo.insert(entity)

      val found = repo.findByBiosample(biosample.id).head
      assertEquals(found.snpCalls.size, 4)

      val m343 = found.snpCalls.find(_.name == "M343")
      assert(m343.isDefined)
      assertEquals(m343.get.startPosition, 2787994L)
      assertEquals(m343.get.endPosition, None)
      assertEquals(m343.get.allele, "A")
      assertEquals(m343.get.derived, true)
      assertEquals(m343.get.variantType, Some(YVariantType.SNP))
      assertEquals(m343.get.quality, Some(99.0))
      assert(!m343.get.isIndel)
    }
  }

  testTransactor.test("stores and retrieves INDEL calls correctly") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val orderDate = LocalDateTime.of(2014, 11, 20, 0, 0)
      val snpCalls = List(
        // SNP
        YSnpCall("A81", 8718573, None, "T", true, Some(YVariantType.SNP), Some(orderDate), None),
        // INDEL with range
        YSnpCall("A1133", 17199439, Some(17199443L), "ins", false, Some(YVariantType.INDEL), Some(orderDate), None)
      )

      val entity = YSnpPanelEntity.create(
        biosampleId = biosample.id,
        snpCalls = snpCalls
      )
      repo.insert(entity)

      val found = repo.findByBiosample(biosample.id).head
      assertEquals(found.snpCalls.size, 2)

      val indel = found.snpCalls.find(_.name == "A1133")
      assert(indel.isDefined)
      assertEquals(indel.get.startPosition, 17199439L)
      assertEquals(indel.get.endPosition, Some(17199443L))
      assertEquals(indel.get.effectiveEndPosition, 17199443L)
      assertEquals(indel.get.allele, "ins")
      assertEquals(indel.get.derived, false)
      assertEquals(indel.get.variantType, Some(YVariantType.INDEL))
      assertEquals(indel.get.orderedDate, Some(orderDate))
      assert(indel.get.isIndel)

      val snp = found.snpCalls.find(_.name == "A81")
      assert(snp.isDefined)
      assertEquals(snp.get.effectiveEndPosition, snp.get.startPosition)
      assert(!snp.get.isIndel)
    }
  }

  testTransactor.test("stores and retrieves private variants correctly") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val privateVariants = List(
        YPrivateVariant(12345678, "G", "A", Some("FGC99999"), Some(95.0), Some(45)),
        YPrivateVariant(23456789, "C", "T", None, Some(88.0), Some(32))
      )

      val entity = YSnpPanelEntity.create(
        biosampleId = biosample.id,
        privateVariants = privateVariants
      )
      repo.insert(entity)

      val found = repo.findByBiosample(biosample.id).head
      assertEquals(found.privateVariants.size, 2)

      val named = found.privateVariants.find(_.snpName.isDefined)
      assert(named.isDefined)
      assertEquals(named.get.position, 12345678L)
      assertEquals(named.get.snpName, Some("FGC99999"))
    }
  }

  testTransactor.test("update modifies entity and increments version") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = YSnpPanelEntity.create(
        biosampleId = biosample.id,
        terminalHaplogroup = Some("R-M269")
      )
      val saved = repo.insert(entity)

      val updated = saved.copy(
        terminalHaplogroup = Some("R-L21"),
        confidence = Some(0.98)
      )
      repo.update(updated)

      val found = repo.findById(saved.id)
      assert(found.isDefined)
      assertEquals(found.get.terminalHaplogroup, Some("R-L21"))
      assertEquals(found.get.confidence, Some(0.98))
      assertEquals(found.get.meta.version, 2)
    }
  }

  testTransactor.test("delete removes entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(YSnpPanelEntity.create(biosample.id))

      assert(repo.exists(entity.id))
      val deleted = repo.delete(entity.id)
      assert(deleted)
      assert(!repo.exists(entity.id))
    }
  }

  testTransactor.test("markSynced updates sync status and AT fields") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val entity = repo.insert(YSnpPanelEntity.create(biosample.id))

      repo.markSynced(entity.id, "at://did:plc:test/ysnppanel/1", "bafycid789")

      val found = repo.findById(entity.id)
      assert(found.isDefined)
      assertEquals(found.get.meta.syncStatus, SyncStatus.Synced)
      assertEquals(found.get.meta.atUri, Some("at://did:plc:test/ysnppanel/1"))
    }
  }

  testTransactor.test("cascades delete on biosample deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val panel = repo.insert(YSnpPanelEntity.create(biosample.id))

      biosampleRepo.delete(biosample.id)

      assertEquals(repo.findById(panel.id), None)
    }
  }

  testTransactor.test("sets alignment_id to null on alignment deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val alignment = createTestAlignment(biosample.id)
      val panel = repo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        alignmentId = Some(alignment.id)
      ))

      // Delete the alignment
      alignmentRepo.delete(alignment.id)

      // Panel should still exist but with null alignment_id
      val found = repo.findById(panel.id)
      assert(found.isDefined)
      assertEquals(found.get.alignmentId, None)
    }
  }
