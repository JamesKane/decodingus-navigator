package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.yprofile.model.*
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class YVariantAuditRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = YVariantAuditRepository()
  val variantRepo = YProfileVariantRepository()
  val profileRepo = YChromosomeProfileRepository()
  val biosampleRepo = BiosampleRepository()

  def createTestBiosample()(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = s"TEST${UUID.randomUUID().toString.take(8)}",
      donorIdentifier = "DONOR001"
    ))

  def createTestProfile(biosampleId: UUID)(using java.sql.Connection): YChromosomeProfileEntity =
    profileRepo.insert(YChromosomeProfileEntity.create(biosampleId = biosampleId))

  def createTestVariant(profileId: UUID, position: Long = 1000L)(using java.sql.Connection): YProfileVariantEntity =
    variantRepo.insert(YProfileVariantEntity.create(profileId, position, "G", "A"))

  testTransactor.test("insert and findById returns correct entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)

      val entity = YVariantAuditEntity.create(
        variantId = variant.id,
        action = YAuditAction.OVERRIDE,
        reason = "IGV inspection shows clear derived signal",
        previousConsensusAllele = Some("G"),
        previousConsensusState = Some(YConsensusState.ANCESTRAL),
        previousStatus = Some(YVariantStatus.CONFLICT),
        previousConfidence = Some(0.5),
        newConsensusAllele = Some("A"),
        newConsensusState = Some(YConsensusState.DERIVED),
        newStatus = Some(YVariantStatus.CONFIRMED),
        newConfidence = Some(0.99),
        userId = Some("user@example.com"),
        supportingEvidence = Some("Screenshot attached: igv_pos1000.png")
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.action, YAuditAction.OVERRIDE)
      assertEquals(found.get.reason, "IGV inspection shows clear derived signal")
      assertEquals(found.get.previousConsensusState, Some(YConsensusState.ANCESTRAL))
      assertEquals(found.get.newConsensusState, Some(YConsensusState.DERIVED))
      assertEquals(found.get.userId, Some("user@example.com"))
    }
  }

  testTransactor.test("findByVariant returns all audit entries for variant") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)

      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "First override"))
      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.ANNOTATE, "Added note"))
      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.REVERT, "Reverted override"))

      val audits = repo.findByVariant(variant.id)
      assertEquals(audits.size, 3)
    }
  }

  testTransactor.test("findByAction filters by action type") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)

      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "Override 1"))
      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "Override 2"))
      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.ANNOTATE, "Note"))

      val overrides = repo.findByAction(YAuditAction.OVERRIDE)
      assertEquals(overrides.size, 2)
    }
  }

  testTransactor.test("findByUser filters by user") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)

      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "By user1", userId = Some("user1@example.com")))
      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "By user1 again", userId = Some("user1@example.com")))
      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "By user2", userId = Some("user2@example.com")))

      val user1 = repo.findByUser("user1@example.com")
      assertEquals(user1.size, 2)
    }
  }

  testTransactor.test("findRecent returns limited results") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)

      for i <- 1 to 10 do
        repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.ANNOTATE, s"Note $i"))

      val recent = repo.findRecent(5)
      assertEquals(recent.size, 5)
    }
  }

  testTransactor.test("findOverrides returns only OVERRIDE actions") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)

      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "Override"))
      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.ANNOTATE, "Note"))

      val overrides = repo.findOverrides(variant.id)
      assertEquals(overrides.size, 1)
      assertEquals(overrides.head.action, YAuditAction.OVERRIDE)
    }
  }

  testTransactor.test("hasOverride returns true when override exists") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)

      assert(!repo.hasOverride(variant.id))

      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "Override"))

      assert(repo.hasOverride(variant.id))
    }
  }

  testTransactor.test("findMostRecent returns latest audit") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)

      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "First"))
      Thread.sleep(10) // Ensure different timestamps
      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.ANNOTATE, "Second"))
      Thread.sleep(10)
      val last = repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.REVERT, "Third"))

      val mostRecent = repo.findMostRecent(variant.id)
      assert(mostRecent.isDefined)
      assertEquals(mostRecent.get.reason, "Third")
    }
  }

  testTransactor.test("deleteByVariant removes all audit entries for variant") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)

      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "Override"))
      repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.ANNOTATE, "Note"))

      val deleted = repo.deleteByVariant(variant.id)
      assertEquals(deleted, 2)
      assertEquals(repo.countByVariant(variant.id), 0L)
    }
  }

  testTransactor.test("delete removes entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)
      val entity = repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.ANNOTATE, "Note"))

      assert(repo.exists(entity.id))
      repo.delete(entity.id)
      assert(!repo.exists(entity.id))
    }
  }

  testTransactor.test("cascades delete on variant deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)
      val audit = repo.insert(YVariantAuditEntity.create(variant.id, YAuditAction.OVERRIDE, "Override"))

      variantRepo.delete(variant.id)

      assertEquals(repo.findById(audit.id), None)
    }
  }

  testTransactor.test("audit entries are immutable except for supporting evidence") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = createTestVariant(profile.id)

      val entity = repo.insert(YVariantAuditEntity.create(
        variant.id,
        YAuditAction.OVERRIDE,
        "Original reason",
        supportingEvidence = Some("Initial evidence")
      ))

      // Update with new supporting evidence
      val updated = entity.copy(
        reason = "Attempted change", // This should NOT change
        supportingEvidence = Some("Updated evidence with screenshot")
      )
      repo.update(updated)

      val found = repo.findById(entity.id)
      // Reason should remain unchanged (audit immutability)
      assertEquals(found.get.reason, "Original reason")
      // Supporting evidence can be updated
      assertEquals(found.get.supportingEvidence, Some("Updated evidence with screenshot"))
    }
  }
