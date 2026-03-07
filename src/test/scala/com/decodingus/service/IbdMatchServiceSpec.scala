package com.decodingus.service

import com.decodingus.db.{DatabaseTestSupport, Transactor}
import com.decodingus.repository.*
import com.decodingus.workspace.model.*
import munit.FunSuite

import java.util.UUID

class IbdMatchServiceSpec extends FunSuite with DatabaseTestSupport:

  private def createService(tx: Transactor): IbdMatchService =
    IbdMatchService(
      transactor = tx,
      consentRepo = MatchConsentRepository(),
      requestRepo = MatchRequestRepository(),
      resultRepo = MatchResultRepository()
    )

  private def createBiosample(tx: Transactor): UUID =
    val biosampleRepo = BiosampleRepository()
    val entity = BiosampleEntity.create(
      sampleAccession = s"TEST-${UUID.randomUUID().toString.take(8)}",
      donorIdentifier = "DONOR-001",
      citizenDid = Some("did:plc:test")
    )
    tx.readWrite { biosampleRepo.insert(entity)(using summon) }
    entity.id

  // ============================================
  // Consent Operations
  // ============================================

  testTransactor.test("grantConsent creates new consent") { case (_, tx) =>
    val service = createService(tx)
    val biosampleId = createBiosample(tx)

    val result = service.grantConsent(biosampleId, s"local://$biosampleId", ConsentLevel.Full)
    assert(result.isRight)
    result.foreach { consent =>
      assertEquals(consent.consentLevel, ConsentLevel.Full)
      assert(consent.allowedMatchTypes.contains(MatchType.Ibd))
    }
  }

  testTransactor.test("grantConsent updates existing consent") { case (_, tx) =>
    val service = createService(tx)
    val biosampleId = createBiosample(tx)

    service.grantConsent(biosampleId, s"local://$biosampleId", ConsentLevel.Full)
    val result = service.grantConsent(biosampleId, s"local://$biosampleId", ConsentLevel.Anonymous)
    assert(result.isRight)
    result.foreach { consent =>
      assertEquals(consent.consentLevel, ConsentLevel.Anonymous)
    }
  }

  testTransactor.test("hasConsent returns true when consent exists") { case (_, tx) =>
    val service = createService(tx)
    val biosampleId = createBiosample(tx)

    service.grantConsent(biosampleId, s"local://$biosampleId", ConsentLevel.Full)
    val result = service.hasConsent(biosampleId)
    assertEquals(result, Right(true))
  }

  testTransactor.test("hasConsent returns false when no consent") { case (_, tx) =>
    val service = createService(tx)
    val biosampleId = createBiosample(tx)

    val result = service.hasConsent(biosampleId)
    assertEquals(result, Right(false))
  }

  testTransactor.test("revokeConsent removes consent") { case (_, tx) =>
    val service = createService(tx)
    val biosampleId = createBiosample(tx)

    service.grantConsent(biosampleId, s"local://$biosampleId", ConsentLevel.Full)
    val revokeResult = service.revokeConsent(biosampleId)
    assertEquals(revokeResult, Right(true))

    val checkResult = service.hasConsent(biosampleId)
    assertEquals(checkResult, Right(false))
  }

  // ============================================
  // Match Request Operations
  // ============================================

  testTransactor.test("sendMatchRequest creates pending request") { case (_, tx) =>
    val service = createService(tx)

    val result = service.sendMatchRequest(
      fromBiosampleRef = "at://did:plc:a/bio/1",
      toBiosampleRef = "at://did:plc:b/bio/1",
      message = Some("I think we share a common ancestor")
    )
    assert(result.isRight)
    result.foreach { request =>
      assertEquals(request.status, RequestStatus.Pending)
      assertEquals(request.fromBiosampleRef, "at://did:plc:a/bio/1")
      assertEquals(request.toBiosampleRef, "at://did:plc:b/bio/1")
      assertEquals(request.message, Some("I think we share a common ancestor"))
    }
  }

  testTransactor.test("getOutgoingRequests returns requests from biosample") { case (_, tx) =>
    val service = createService(tx)

    service.sendMatchRequest("at://did:plc:a/bio/1", "at://did:plc:b/bio/1")
    service.sendMatchRequest("at://did:plc:a/bio/1", "at://did:plc:c/bio/1")
    service.sendMatchRequest("at://did:plc:other/bio/1", "at://did:plc:d/bio/1")

    val result = service.getOutgoingRequests("at://did:plc:a/bio/1")
    assert(result.isRight)
    assertEquals(result.map(_.size), Right(2))
  }

  testTransactor.test("getPendingRequests returns only pending requests") { case (_, tx) =>
    val service = createService(tx)

    val Right(req1) = service.sendMatchRequest("at://did:plc:x/bio/1", "at://did:plc:target/bio/1"): @unchecked
    service.sendMatchRequest("at://did:plc:y/bio/1", "at://did:plc:target/bio/1")

    // Accept first request
    req1.atUri.flatMap(uri => scala.util.Try(UUID.fromString(uri.split("/").last)).toOption).foreach { id =>
      service.respondToRequest(id, accept = true)
    }

    val result = service.getPendingRequests("at://did:plc:target/bio/1")
    assert(result.isRight)
    // Should have at most 1 pending (the second one)
    result.foreach { pending =>
      assert(pending.forall(_.status == RequestStatus.Pending))
    }
  }

  testTransactor.test("respondToRequest changes status") { case (_, tx) =>
    val service = createService(tx)
    val requestRepo = MatchRequestRepository()

    val Right(_) = service.sendMatchRequest("at://did:plc:a/bio/1", "at://did:plc:b/bio/1"): @unchecked

    // Extract the actual entity ID from the repository
    val Right(entityId) = tx.readOnly {
      requestRepo.findAll()(using summon).head.id
    }: @unchecked

    val result = service.respondToRequest(entityId, accept = false)
    assertEquals(result, Right(true))

    // Verify status changed
    val Right(updated) = tx.readOnly {
      requestRepo.findById(entityId)(using summon)
    }: @unchecked
    updated.foreach { entity =>
      assertEquals(entity.status, "DECLINED")
      assert(entity.respondedAt.isDefined)
    }
  }

  // ============================================
  // Match Result Operations
  // ============================================

  testTransactor.test("getMatchResults returns results for biosample") { case (_, tx) =>
    val service = createService(tx)
    val biosampleId = createBiosample(tx)
    val resultRepo = MatchResultRepository()

    // Insert a match result directly
    tx.readWrite {
      resultRepo.insert(MatchResultEntity.create(
        biosampleId = biosampleId,
        matchedBiosampleRef = "at://did:plc:b/bio/1",
        totalSharedCm = 150.5,
        segmentCount = 5,
        longestSegmentCm = Some(45.2),
        relationshipEstimate = Some("2ND_COUSIN")
      ))(using summon)
    }

    val result = service.getMatchResults(biosampleId)
    assert(result.isRight)
    result.foreach { results =>
      assertEquals(results.size, 1)
      assertEquals(results.head.totalSharedCm, 150.5)
      assertEquals(results.head.segmentCount, 5)
      assertEquals(results.head.relationshipEstimate, Some(RelationshipEstimate.SecondCousin))
    }
  }

  testTransactor.test("getMatchResultsAboveThreshold filters by cM") { case (_, tx) =>
    val service = createService(tx)
    val biosampleId = createBiosample(tx)
    val resultRepo = MatchResultRepository()

    tx.readWrite {
      resultRepo.insert(MatchResultEntity.create(
        biosampleId = biosampleId,
        matchedBiosampleRef = "at://did:plc:b/bio/1",
        totalSharedCm = 150.5,
        segmentCount = 5
      ))(using summon)
      resultRepo.insert(MatchResultEntity.create(
        biosampleId = biosampleId,
        matchedBiosampleRef = "at://did:plc:c/bio/1",
        totalSharedCm = 5.2,
        segmentCount = 1
      ))(using summon)
    }

    val allResults = service.getMatchResults(biosampleId)
    assertEquals(allResults.map(_.size), Right(2))

    val filteredResults = service.getMatchResultsAboveThreshold(biosampleId, 10.0)
    assertEquals(filteredResults.map(_.size), Right(1))
    filteredResults.foreach { results =>
      assertEquals(results.head.totalSharedCm, 150.5)
    }
  }
