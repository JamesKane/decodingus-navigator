package com.decodingus.repository

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.yprofile.model.*
import munit.FunSuite
import java.time.LocalDateTime
import java.util.UUID

class YProfileVariantRepositorySpec extends FunSuite with DatabaseTestSupport:

  val repo = YProfileVariantRepository()
  val profileRepo = YChromosomeProfileRepository()
  val biosampleRepo = BiosampleRepository()

  def createTestBiosample()(using java.sql.Connection): BiosampleEntity =
    biosampleRepo.insert(BiosampleEntity.create(
      sampleAccession = s"TEST${UUID.randomUUID().toString.take(8)}",
      donorIdentifier = "DONOR001"
    ))

  def createTestProfile(biosampleId: UUID)(using java.sql.Connection): YChromosomeProfileEntity =
    profileRepo.insert(YChromosomeProfileEntity.create(biosampleId = biosampleId))

  testTransactor.test("insert and findById returns correct SNP entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      val entity = YProfileVariantEntity.create(
        yProfileId = profile.id,
        position = 2787994L,
        refAllele = "G",
        altAllele = "A",
        variantType = YVariantType.SNP,
        variantName = Some("M343"),
        rsId = Some("rs9786184"),
        consensusAllele = Some("A"),
        consensusState = YConsensusState.DERIVED,
        status = YVariantStatus.CONFIRMED,
        sourceCount = 3,
        concordantCount = 3,
        discordantCount = 0,
        confidenceScore = 0.99,
        maxReadDepth = Some(45),
        maxQualityScore = Some(99.0),
        definingHaplogroup = Some("R1b"),
        haplogroupBranchDepth = Some(5)
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.position, 2787994L)
      assertEquals(found.get.variantName, Some("M343"))
      assertEquals(found.get.variantType, YVariantType.SNP)
      assertEquals(found.get.consensusState, YConsensusState.DERIVED)
      assertEquals(found.get.status, YVariantStatus.CONFIRMED)
      assertEquals(found.get.confidenceScore, 0.99)
    }
  }

  testTransactor.test("insert and findById returns correct STR entity with metadata") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      val strMeta = StrMetadata(
        repeatMotif = Some("GATA"),
        repeatUnit = Some(4),
        copies = None,
        rawNotation = None
      )

      val entity = YProfileVariantEntity.create(
        yProfileId = profile.id,
        position = 12500000L,
        refAllele = "(GATA)13",
        altAllele = "(GATA)14",
        variantType = YVariantType.STR,
        markerName = Some("DYS393"),
        repeatCount = Some(13),
        strMetadata = Some(strMeta),
        consensusState = YConsensusState.DERIVED,
        status = YVariantStatus.CONFIRMED
      )

      val saved = repo.insert(entity)
      val found = repo.findById(saved.id)

      assert(found.isDefined)
      assertEquals(found.get.variantType, YVariantType.STR)
      assertEquals(found.get.markerName, Some("DYS393"))
      assertEquals(found.get.repeatCount, Some(13))
      assert(found.get.strMetadata.isDefined)
      assertEquals(found.get.strMetadata.get.repeatMotif, Some("GATA"))
      assertEquals(found.get.strMetadata.get.repeatUnit, Some(4))
    }
  }

  testTransactor.test("findByProfile returns all variants for profile") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileVariantEntity.create(profile.id, 1000L, "G", "A"))
      repo.insert(YProfileVariantEntity.create(profile.id, 2000L, "C", "T"))
      repo.insert(YProfileVariantEntity.create(profile.id, 3000L, "A", "G"))

      val variants = repo.findByProfile(profile.id)
      assertEquals(variants.size, 3)
      // Should be ordered by position
      assertEquals(variants.map(_.position), List(1000L, 2000L, 3000L))
    }
  }

  testTransactor.test("findByType filters by variant type") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileVariantEntity.create(profile.id, 1000L, "G", "A", variantType = YVariantType.SNP))
      repo.insert(YProfileVariantEntity.create(profile.id, 2000L, "C", "T", variantType = YVariantType.SNP))
      repo.insert(YProfileVariantEntity.create(profile.id, 3000L, "(GATA)13", "(GATA)14", variantType = YVariantType.STR))

      val snps = repo.findByType(profile.id, YVariantType.SNP)
      assertEquals(snps.size, 2)

      val strs = repo.findByType(profile.id, YVariantType.STR)
      assertEquals(strs.size, 1)
    }
  }

  testTransactor.test("findByStatus filters by status") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileVariantEntity.create(profile.id, 1000L, "G", "A", status = YVariantStatus.CONFIRMED))
      repo.insert(YProfileVariantEntity.create(profile.id, 2000L, "C", "T", status = YVariantStatus.NOVEL))
      repo.insert(YProfileVariantEntity.create(profile.id, 3000L, "A", "G", status = YVariantStatus.CONFLICT))

      val confirmed = repo.findByStatus(profile.id, YVariantStatus.CONFIRMED)
      assertEquals(confirmed.size, 1)

      val conflicts = repo.findByStatus(profile.id, YVariantStatus.CONFLICT)
      assertEquals(conflicts.size, 1)
    }
  }

  testTransactor.test("findDerivedVariants returns only DERIVED state") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileVariantEntity.create(profile.id, 1000L, "G", "A", consensusState = YConsensusState.DERIVED))
      repo.insert(YProfileVariantEntity.create(profile.id, 2000L, "C", "T", consensusState = YConsensusState.ANCESTRAL))
      repo.insert(YProfileVariantEntity.create(profile.id, 3000L, "A", "G", consensusState = YConsensusState.DERIVED))

      val derived = repo.findDerivedVariants(profile.id)
      assertEquals(derived.size, 2)
    }
  }

  testTransactor.test("findNovelVariants returns only NOVEL status") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileVariantEntity.create(profile.id, 1000L, "G", "A", status = YVariantStatus.CONFIRMED))
      repo.insert(YProfileVariantEntity.create(profile.id, 2000L, "C", "T", status = YVariantStatus.NOVEL))

      val novel = repo.findNovelVariants(profile.id)
      assertEquals(novel.size, 1)
      assertEquals(novel.head.status, YVariantStatus.NOVEL)
    }
  }

  testTransactor.test("findConflictingVariants returns only CONFLICT status") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileVariantEntity.create(
        profile.id, 1000L, "G", "A",
        status = YVariantStatus.CONFLICT,
        discordantCount = 2
      ))
      repo.insert(YProfileVariantEntity.create(profile.id, 2000L, "C", "T", status = YVariantStatus.CONFIRMED))

      val conflicts = repo.findConflictingVariants(profile.id)
      assertEquals(conflicts.size, 1)
      assertEquals(conflicts.head.discordantCount, 2)
    }
  }

  testTransactor.test("findStrMarkers returns only STR type") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileVariantEntity.create(
        profile.id, 1000L, "(GATA)13", "(GATA)14",
        variantType = YVariantType.STR,
        markerName = Some("DYS393")
      ))
      repo.insert(YProfileVariantEntity.create(
        profile.id, 2000L, "(TAGA)9", "(TAGA)10",
        variantType = YVariantType.STR,
        markerName = Some("DYS19")
      ))
      repo.insert(YProfileVariantEntity.create(profile.id, 3000L, "G", "A", variantType = YVariantType.SNP))

      val strs = repo.findStrMarkers(profile.id)
      assertEquals(strs.size, 2)
    }
  }

  testTransactor.test("findByVariantName searches by marker name") { case (db, tx) =>
    tx.readWrite {
      val b1 = createTestBiosample()
      val b2 = createTestBiosample()
      val p1 = createTestProfile(b1.id)
      val p2 = createTestProfile(b2.id)

      repo.insert(YProfileVariantEntity.create(p1.id, 2787994L, "G", "A", variantName = Some("M343")))
      repo.insert(YProfileVariantEntity.create(p2.id, 2787994L, "G", "A", variantName = Some("M343")))
      repo.insert(YProfileVariantEntity.create(p1.id, 22739368L, "T", "C", variantName = Some("M269")))

      val m343 = repo.findByVariantName("M343")
      assertEquals(m343.size, 2)
    }
  }

  testTransactor.test("countByStatus returns status counts") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileVariantEntity.create(profile.id, 1000L, "G", "A", status = YVariantStatus.CONFIRMED))
      repo.insert(YProfileVariantEntity.create(profile.id, 2000L, "C", "T", status = YVariantStatus.CONFIRMED))
      repo.insert(YProfileVariantEntity.create(profile.id, 3000L, "A", "G", status = YVariantStatus.NOVEL))

      val counts = repo.countByStatus(profile.id)
      assertEquals(counts.get(YVariantStatus.CONFIRMED), Some(2L))
      assertEquals(counts.get(YVariantStatus.NOVEL), Some(1L))
    }
  }

  testTransactor.test("update modifies entity") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val entity = repo.insert(YProfileVariantEntity.create(
        profile.id, 1000L, "G", "A",
        status = YVariantStatus.PENDING,
        confidenceScore = 0.5
      ))

      val updated = entity.copy(
        status = YVariantStatus.CONFIRMED,
        confidenceScore = 0.99,
        concordantCount = 3
      )
      repo.update(updated)

      val found = repo.findById(entity.id)
      assertEquals(found.get.status, YVariantStatus.CONFIRMED)
      assertEquals(found.get.confidenceScore, 0.99)
      assertEquals(found.get.concordantCount, 3)
    }
  }

  testTransactor.test("deleteByProfile removes all variants for profile") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)

      repo.insert(YProfileVariantEntity.create(profile.id, 1000L, "G", "A"))
      repo.insert(YProfileVariantEntity.create(profile.id, 2000L, "C", "T"))

      val deleted = repo.deleteByProfile(profile.id)
      assertEquals(deleted, 2)
      assertEquals(repo.findByProfile(profile.id).size, 0)
    }
  }

  testTransactor.test("cascades delete on profile deletion") { case (db, tx) =>
    tx.readWrite {
      val biosample = createTestBiosample()
      val profile = createTestProfile(biosample.id)
      val variant = repo.insert(YProfileVariantEntity.create(profile.id, 1000L, "G", "A"))

      profileRepo.delete(profile.id)

      assertEquals(repo.findById(variant.id), None)
    }
  }
