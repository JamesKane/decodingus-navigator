package com.decodingus.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.*
import com.decodingus.workspace.model.{
  FileInfo, SimpleStrValue, StrMarkerValue, StrPanel,
  DnaType, CompatibilityLevel, HaplogroupTechnology, CallMethod, ConflictResolution,
  ReconciliationStatus, RunHaplogroupCall, SnpCallFromRun, SnpConflict
}
import munit.FunSuite
import java.time.{Instant, LocalDateTime}
import java.util.UUID

/**
 * Integration tests for Phase 2 entities:
 * - STR Profiles
 * - Chip Profiles
 * - Y-SNP Panels
 * - Haplogroup Reconciliation
 *
 * Tests cross-repository operations and relationships.
 */
class Phase2EntityIntegrationSpec extends FunSuite with DatabaseTestSupport:

  // ============================================
  // Repository Setup
  // ============================================

  val biosampleRepo = BiosampleRepository()
  val sequenceRunRepo = SequenceRunRepository()
  val alignmentRepo = AlignmentRepository()
  val strProfileRepo = StrProfileRepository()
  val chipProfileRepo = ChipProfileRepository()
  val ySnpPanelRepo = YSnpPanelRepository()
  val haplogroupReconciliationRepo = HaplogroupReconciliationRepository()

  def createTestBiosample(accession: String = s"TEST${UUID.randomUUID().toString.take(8)}")(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = accession,
      donorIdentifier = "DONOR001"
    ))

  def createTestAlignment(biosampleId: UUID)(using java.sql.Connection): AlignmentEntity =
    val seqRun = sequenceRunRepo.insert(SequenceRunEntity.create(
      biosampleId = biosampleId,
      platform = "ILLUMINA",
      testType = "WGS"
    ))
    alignmentRepo.insert(AlignmentEntity.create(
      sequenceRunId = seqRun.id,
      referenceBuild = "GRCh38",
      aligner = "BWA-MEM2"
    ))

  // ============================================
  // STR Profile Integration Tests
  // ============================================

  testTransactor.test("STR profile workflow: create profile, query by biosample") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample("STR-WORKFLOW-001")

      // Create STR profile with markers
      val strProfile = strProfileRepo.insert(StrProfileEntity.create(
        biosampleId = biosample.id,
        source = Some("DIRECT_TEST"),
        importedFrom = Some("FTDNA"),
        totalMarkers = Some(111),
        panels = List(
          StrPanel("Y-37", 37, Some("FTDNA"), None),
          StrPanel("Y-67", 67, Some("FTDNA"), None),
          StrPanel("Y-111", 111, Some("FTDNA"), None)
        ),
        markers = List(
          StrMarkerValue("DYS393", SimpleStrValue(13), panel = Some("Y12")),
          StrMarkerValue("DYS390", SimpleStrValue(24), panel = Some("Y12")),
          StrMarkerValue("DYS19", SimpleStrValue(14), panel = Some("Y25")),
          StrMarkerValue("DYS385a", SimpleStrValue(11), panel = Some("Y37")),
          StrMarkerValue("DYS385b", SimpleStrValue(14), panel = Some("Y37"))
        )
      ))

      // Query by biosample
      val profiles = strProfileRepo.findByBiosample(biosample.id)
      assertEquals(profiles.size, 1)
      assertEquals(profiles.head.markers.size, 5)
      assertEquals(profiles.head.panels.size, 3)
    }
  }

  testTransactor.test("STR profile: multiple profiles per biosample from different sources") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample("MULTI-STR-001")

      // Direct test profile
      strProfileRepo.insert(StrProfileEntity.create(
        biosampleId = biosample.id,
        source = Some("DIRECT_TEST"),
        importedFrom = Some("FTDNA"),
        totalMarkers = Some(37)
      ))

      // WGS-derived profile
      strProfileRepo.insert(StrProfileEntity.create(
        biosampleId = biosample.id,
        source = Some("WGS_DERIVED"),
        derivationMethod = Some("HIPSTR"),
        totalMarkers = Some(500)
      ))

      val profiles = strProfileRepo.findByBiosample(biosample.id)
      assertEquals(profiles.size, 2)

      val directProfiles = strProfileRepo.findBySource("DIRECT_TEST")
      assertEquals(directProfiles.size, 1)

      val wgsProfiles = strProfileRepo.findBySource("WGS_DERIVED")
      assertEquals(wgsProfiles.size, 1)
    }
  }

  testTransactor.test("STR profile: cascade delete on biosample removal") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample("STR-CASCADE-001")

      val profile = strProfileRepo.insert(StrProfileEntity.create(
        biosampleId = biosample.id,
        source = Some("DIRECT_TEST")
      ))

      // Verify profile exists
      assert(strProfileRepo.exists(profile.id))

      // Delete biosample
      biosampleRepo.delete(biosample.id)

      // Profile should be deleted via cascade
      assertEquals(strProfileRepo.findById(profile.id), None)
    }
  }

  // ============================================
  // Chip Profile Integration Tests
  // ============================================

  testTransactor.test("Chip profile workflow: create profile with marker counts") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample("CHIP-WORKFLOW-001")

      val chipProfile = chipProfileRepo.insert(ChipProfileEntity.create(
        biosampleId = biosample.id,
        vendor = "ILLUMINA",
        testTypeCode = "ANCESTRY_V2",
        totalMarkersCalled = 650000,
        totalMarkersPossible = 700000,
        noCallRate = 0.07,
        autosomalMarkersCalled = 600000,
        importDate = LocalDateTime.now(),
        chipVersion = Some("Illumina Global Screening Array v2"),
        files = List(FileInfo(
          fileName = "chip_data.csv",
          fileSizeBytes = Some(50000000L),
          fileFormat = "CSV",
          checksum = None,
          location = Some("/data/chip/chip_data.csv")
        ))
      ))

      val profiles = chipProfileRepo.findByBiosample(biosample.id)
      assertEquals(profiles.size, 1)
      assertEquals(profiles.head.vendor, "ILLUMINA")
      assertEquals(profiles.head.totalMarkersCalled, 650000)
    }
  }

  testTransactor.test("Chip profile: find by vendor") { case (db, tx) =>
    tx.readWrite {
      val biosample1 = createTestBiosample("CHIP-VENDOR-001")
      val biosample2 = createTestBiosample("CHIP-VENDOR-002")

      chipProfileRepo.insert(ChipProfileEntity.create(
        biosampleId = biosample1.id,
        vendor = "ILLUMINA",
        testTypeCode = "ANCESTRY_V2",
        totalMarkersCalled = 650000,
        totalMarkersPossible = 700000,
        noCallRate = 0.07,
        autosomalMarkersCalled = 600000,
        importDate = LocalDateTime.now()
      ))
      chipProfileRepo.insert(ChipProfileEntity.create(
        biosampleId = biosample2.id,
        vendor = "ILLUMINA",
        testTypeCode = "GSA_V3",
        totalMarkersCalled = 700000,
        totalMarkersPossible = 750000,
        noCallRate = 0.06,
        autosomalMarkersCalled = 650000,
        importDate = LocalDateTime.now()
      ))
      chipProfileRepo.insert(ChipProfileEntity.create(
        biosampleId = biosample1.id,
        vendor = "23ANDME",
        testTypeCode = "CUSTOM_V5",
        totalMarkersCalled = 650000,
        totalMarkersPossible = 680000,
        noCallRate = 0.04,
        autosomalMarkersCalled = 600000,
        importDate = LocalDateTime.now()
      ))

      val illuminaProfiles = chipProfileRepo.findByVendor("ILLUMINA")
      assertEquals(illuminaProfiles.size, 2)

      val andMeProfiles = chipProfileRepo.findByVendor("23ANDME")
      assertEquals(andMeProfiles.size, 1)
    }
  }

  // ============================================
  // Y-SNP Panel Integration Tests
  // ============================================

  testTransactor.test("Y-SNP panel workflow: Big Y import with SNP calls and private variants") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample("BIG-Y-001")
      val orderDate = LocalDateTime.of(2023, 6, 15, 0, 0)

      val panel = ySnpPanelRepo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        panelName = Some("Big Y-700"),
        provider = Some("FTDNA"),
        testDate = Some(orderDate),
        totalSnpsTested = Some(35000),
        derivedCount = Some(485),
        ancestralCount = Some(34000),
        noCallCount = Some(515),
        terminalHaplogroup = Some("R-BY140757"),
        confidence = Some(0.99),
        snpCalls = List(
          YSnpCall("M343", 2787994, None, "A", true, Some(YVariantType.SNP), Some(orderDate), Some(99.0)),
          YSnpCall("M269", 22739368, None, "C", true, Some(YVariantType.SNP), Some(orderDate), Some(98.5)),
          YSnpCall("P312", 15579988, None, "C", true, None, None, Some(97.0)),
          YSnpCall("L21", 14722131, None, "A", true, None, None, Some(96.5)),
          // INDEL example
          YSnpCall("A1133", 17199439, Some(17199443L), "ins", false, Some(YVariantType.INDEL), Some(orderDate), None)
        ),
        privateVariants = List(
          YPrivateVariant(12345678, "G", "A", Some("FGC99999"), Some(95.0), Some(45)),
          YPrivateVariant(23456789, "C", "T", None, Some(88.0), Some(32))
        )
      ))

      val found = ySnpPanelRepo.findById(panel.id)
      assert(found.isDefined)
      assertEquals(found.get.terminalHaplogroup, Some("R-BY140757"))
      assertEquals(found.get.snpCalls.size, 5)
      assertEquals(found.get.privateVariants.size, 2)

      // Test INDEL detection
      val indel = found.get.snpCalls.find(_.name == "A1133")
      assert(indel.isDefined)
      assert(indel.get.isIndel)
      assertEquals(indel.get.endPosition, Some(17199443L))
    }
  }

  testTransactor.test("Y-SNP panel: link to WGS alignment") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample("WGS-SNP-001")
      val alignment = createTestAlignment(biosample.id)

      // Create panel derived from WGS alignment
      val panel = ySnpPanelRepo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        alignmentId = Some(alignment.id),
        panelName = Some("WGS-derived Y-SNP"),
        provider = Some("INTERNAL")
      ))

      // Find by alignment
      val alignmentPanels = ySnpPanelRepo.findByAlignment(alignment.id)
      assertEquals(alignmentPanels.size, 1)
      assertEquals(alignmentPanels.head.id, panel.id)

      // Verify alignment deletion sets alignment_id to null but keeps panel
      alignmentRepo.delete(alignment.id)

      val updatedPanel = ySnpPanelRepo.findById(panel.id)
      assert(updatedPanel.isDefined)
      assertEquals(updatedPanel.get.alignmentId, None)
    }
  }

  testTransactor.test("Y-SNP panel: haplogroup branch queries") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample("HAPLO-BRANCH-001")

      ySnpPanelRepo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        terminalHaplogroup = Some("R-M269")
      ))
      ySnpPanelRepo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        terminalHaplogroup = Some("R-L21")
      ))
      ySnpPanelRepo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        terminalHaplogroup = Some("R-BY140757")
      ))
      ySnpPanelRepo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        terminalHaplogroup = Some("I1-M253")
      ))

      // Find all R branch
      val rBranch = ySnpPanelRepo.findByHaplogroupBranch("R-")
      assertEquals(rBranch.size, 3)

      // Find exact haplogroup
      val exactMatch = ySnpPanelRepo.findByHaplogroup("R-L21")
      assertEquals(exactMatch.size, 1)
    }
  }

  // ============================================
  // Haplogroup Reconciliation Integration Tests
  // ============================================

  testTransactor.test("Haplogroup reconciliation: basic create and query") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample("RECON-001")
      val sourceRef = s"local:biosample:${biosample.id}"

      val reconciliation = haplogroupReconciliationRepo.insert(HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = ReconciliationStatus(
          compatibilityLevel = CompatibilityLevel.COMPATIBLE,
          consensusHaplogroup = "R-BY140757",
          confidence = 0.95,
          runCount = 1
        ),
        runCalls = List(
          RunHaplogroupCall(
            sourceRef = sourceRef,
            haplogroup = "R-BY140757",
            confidence = 0.99,
            callMethod = CallMethod.SNP_PHYLOGENETIC,
            technology = Some(HaplogroupTechnology.BIG_Y),
            treeProvider = Some("FTDNA"),
            treeVersion = Some("2024.1"),
            supportingSnps = Some(485)
          )
        )
      ))

      val found = haplogroupReconciliationRepo.findById(reconciliation.id)
      assert(found.isDefined)
      assertEquals(found.get.status.consensusHaplogroup, "R-BY140757")
      assertEquals(found.get.runCalls.size, 1)
    }
  }

  testTransactor.test("Haplogroup reconciliation: query by DNA type and biosample") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample("DNA-TYPE-001")

      // Create Y-DNA reconciliation
      haplogroupReconciliationRepo.insert(HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = ReconciliationStatus(
          compatibilityLevel = CompatibilityLevel.COMPATIBLE,
          consensusHaplogroup = "R-M269",
          confidence = 0.98,
          runCount = 1
        )
      ))

      // Create mtDNA reconciliation
      haplogroupReconciliationRepo.insert(HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.MT_DNA,
        status = ReconciliationStatus(
          compatibilityLevel = CompatibilityLevel.COMPATIBLE,
          consensusHaplogroup = "H1a",
          confidence = 0.96,
          runCount = 1
        )
      ))

      // Find by biosample
      val biosampleReconciliations = haplogroupReconciliationRepo.findByBiosample(biosample.id)
      assertEquals(biosampleReconciliations.size, 2)

      // Find by DNA type
      val yDnaReconciliations = haplogroupReconciliationRepo.findByDnaType(DnaType.Y_DNA)
      assertEquals(yDnaReconciliations.size, 1)

      val mtDnaReconciliations = haplogroupReconciliationRepo.findByDnaType(DnaType.MT_DNA)
      assertEquals(mtDnaReconciliations.size, 1)
    }
  }

  // ============================================
  // Cross-Entity Workflow Tests
  // ============================================

  testTransactor.test("workflow: complete sample import with all Phase 2 entities") { case (db, tx) =>
    tx.readWrite {
      // 1. Create biosample
      val biosample = createTestBiosample("FULL-IMPORT-001")

      // 2. Create alignment from WGS
      val alignment = createTestAlignment(biosample.id)

      // 3. Create STR profile from panel test
      val strProfile = strProfileRepo.insert(StrProfileEntity.create(
        biosampleId = biosample.id,
        source = Some("DIRECT_TEST"),
        importedFrom = Some("FTDNA"),
        totalMarkers = Some(111),
        markers = List(
          StrMarkerValue("DYS393", SimpleStrValue(13), panel = Some("Y12")),
          StrMarkerValue("DYS390", SimpleStrValue(24), panel = Some("Y12"))
        )
      ))

      // 4. Create Y-SNP panel from Big Y
      val ySnpPanel = ySnpPanelRepo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        panelName = Some("Big Y-700"),
        provider = Some("FTDNA"),
        terminalHaplogroup = Some("R-BY140757"),
        snpCalls = List(
          YSnpCall("M343", 2787994, None, "A", true, Some(YVariantType.SNP), None, None),
          YSnpCall("M269", 22739368, None, "C", true, Some(YVariantType.SNP), None, None)
        )
      ))

      // 5. Create chip profile from autosomal test
      val chipProfile = chipProfileRepo.insert(ChipProfileEntity.create(
        biosampleId = biosample.id,
        vendor = "ILLUMINA",
        testTypeCode = "GSA_V3",
        totalMarkersCalled = 650000,
        totalMarkersPossible = 700000,
        noCallRate = 0.07,
        autosomalMarkersCalled = 600000,
        importDate = LocalDateTime.now()
      ))

      // 6. Create haplogroup reconciliation
      val reconciliation = haplogroupReconciliationRepo.insert(HaplogroupReconciliationEntity.create(
        biosampleId = biosample.id,
        dnaType = DnaType.Y_DNA,
        status = ReconciliationStatus(
          compatibilityLevel = CompatibilityLevel.COMPATIBLE,
          consensusHaplogroup = "R-BY140757",
          confidence = 0.99,
          runCount = 1
        ),
        runCalls = List(
          RunHaplogroupCall(
            sourceRef = s"local:biosample:${biosample.id}",
            haplogroup = "R-BY140757",
            confidence = 0.99,
            callMethod = CallMethod.SNP_PHYLOGENETIC,
            technology = Some(HaplogroupTechnology.BIG_Y)
          )
        )
      ))

      // Verify all entities exist and are linked to biosample
      val strProfiles = strProfileRepo.findByBiosample(biosample.id)
      assertEquals(strProfiles.size, 1)

      val ySnpPanels = ySnpPanelRepo.findByBiosample(biosample.id)
      assertEquals(ySnpPanels.size, 1)

      val chipProfiles = chipProfileRepo.findByBiosample(biosample.id)
      assertEquals(chipProfiles.size, 1)

      val reconciliations = haplogroupReconciliationRepo.findByBiosample(biosample.id)
      assertEquals(reconciliations.size, 1)

      // Test cascade delete
      biosampleRepo.delete(biosample.id)

      // All Phase 2 entities should be deleted
      assertEquals(strProfileRepo.findById(strProfile.id), None)
      assertEquals(ySnpPanelRepo.findById(ySnpPanel.id), None)
      assertEquals(chipProfileRepo.findById(chipProfile.id), None)
      assertEquals(haplogroupReconciliationRepo.findById(reconciliation.id), None)
    }
  }

  testTransactor.test("workflow: sync status tracking across Phase 2 entities") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample("SYNC-PHASE2-001")

      // Create entities (all start with Local status)
      val strProfile = strProfileRepo.insert(StrProfileEntity.create(
        biosampleId = biosample.id,
        source = Some("DIRECT_TEST")
      ))

      val ySnpPanel = ySnpPanelRepo.insert(YSnpPanelEntity.create(
        biosampleId = biosample.id,
        panelName = Some("Big Y-700")
      ))

      // Mark STR profile as synced
      strProfileRepo.markSynced(strProfile.id, "at://did:plc:test/strprofile/1", "cid123")

      // Verify sync statuses
      val foundStr = strProfileRepo.findById(strProfile.id).get
      assertEquals(foundStr.meta.syncStatus, SyncStatus.Synced)

      val foundYSnp = ySnpPanelRepo.findById(ySnpPanel.id).get
      assertEquals(foundYSnp.meta.syncStatus, SyncStatus.Local)

      // Find pending sync entities
      val pendingStr = strProfileRepo.findPendingSync()
      assertEquals(pendingStr.size, 0) // STR is synced

      val pendingYSnp = ySnpPanelRepo.findPendingSync()
      assertEquals(pendingYSnp.size, 1) // Y-SNP is still local
    }
  }
