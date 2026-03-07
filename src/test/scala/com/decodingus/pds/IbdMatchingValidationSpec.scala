package com.decodingus.pds

import com.decodingus.workspace.model.*
import munit.FunSuite

import java.time.LocalDateTime

class IbdMatchingValidationSpec extends FunSuite:

  // ============================================
  // MatchConsent Validation
  // ============================================

  test("validateMatchConsent passes with required fields") {
    val mc = makeConsent()
    assert(PdsSyncValidation.validateMatchConsent(mc).isRight)
  }

  test("validateMatchConsent fails when missing biosampleRef") {
    val mc = makeConsent(biosampleRef = "")
    val Left(errors) = PdsSyncValidation.validateMatchConsent(mc): @unchecked
    assert(errors.exists(_.contains("biosampleRef")))
  }

  test("validateMatchConsent fails when missing consentLevel") {
    val mc = makeConsent(consentLevel = "")
    val Left(errors) = PdsSyncValidation.validateMatchConsent(mc): @unchecked
    assert(errors.exists(_.contains("consentLevel")))
  }

  test("validateMatchConsent fails with invalid consentLevel") {
    val mc = makeConsent(consentLevel = "INVALID")
    val Left(errors) = PdsSyncValidation.validateMatchConsent(mc): @unchecked
    assert(errors.exists(_.contains("invalid consentLevel")))
  }

  test("validateMatchConsent accepts all valid consent levels") {
    for level <- MatchConsent.ConsentLevels do
      val mc = makeConsent(consentLevel = level)
      assert(PdsSyncValidation.validateMatchConsent(mc).isRight, s"Failed for $level")
  }

  // ============================================
  // MatchRequest Validation
  // ============================================

  test("validateMatchRequest passes with required fields") {
    val mr = makeRequest()
    assert(PdsSyncValidation.validateMatchRequest(mr).isRight)
  }

  test("validateMatchRequest fails when missing atUri") {
    val mr = makeRequest(atUri = None)
    val Left(errors) = PdsSyncValidation.validateMatchRequest(mr): @unchecked
    assert(errors.exists(_.contains("atUri")))
  }

  test("validateMatchRequest fails when missing fromBiosampleRef") {
    val mr = makeRequest(fromBiosampleRef = "")
    val Left(errors) = PdsSyncValidation.validateMatchRequest(mr): @unchecked
    assert(errors.exists(_.contains("fromBiosampleRef")))
  }

  test("validateMatchRequest fails when missing toBiosampleRef") {
    val mr = makeRequest(toBiosampleRef = "")
    val Left(errors) = PdsSyncValidation.validateMatchRequest(mr): @unchecked
    assert(errors.exists(_.contains("toBiosampleRef")))
  }

  test("validateMatchRequest collects multiple errors") {
    val mr = makeRequest(atUri = None, fromBiosampleRef = "", toBiosampleRef = "", status = "")
    val Left(errors) = PdsSyncValidation.validateMatchRequest(mr): @unchecked
    assertEquals(errors.size, 4)
  }

  // ============================================
  // IbdSegment Model
  // ============================================

  test("IbdSegment stores chromosome segment data") {
    val segment = IbdSegment(
      chromosome = "1",
      startPosition = 1000000,
      endPosition = 5000000,
      lengthCm = 12.5,
      snpCount = Some(500),
      isHalfIdentical = Some(true)
    )
    assertEquals(segment.chromosome, "1")
    assertEquals(segment.lengthCm, 12.5)
    assert(segment.isHalfIdentical.contains(true))
  }

  // ============================================
  // MatchConsent Known Values
  // ============================================

  test("MatchConsent.ConsentLevels contains expected values") {
    assert(MatchConsent.ConsentLevels.contains("FULL"))
    assert(MatchConsent.ConsentLevels.contains("ANONYMOUS"))
    assert(MatchConsent.ConsentLevels.contains("PROJECT_ONLY"))
  }

  test("MatchConsent.MatchTypes contains expected values") {
    assert(MatchConsent.MatchTypes.contains("IBD"))
    assert(MatchConsent.MatchTypes.contains("Y_STR"))
    assert(MatchConsent.MatchTypes.contains("MT_SEQUENCE"))
    assert(MatchConsent.MatchTypes.contains("AUTOSOMAL"))
  }

  test("MatchRequest.Statuses contains expected values") {
    assert(MatchRequest.Statuses.contains("PENDING"))
    assert(MatchRequest.Statuses.contains("ACCEPTED"))
    assert(MatchRequest.Statuses.contains("DECLINED"))
    assert(MatchRequest.Statuses.contains("EXPIRED"))
    assert(MatchRequest.Statuses.contains("WITHDRAWN"))
  }

  test("MatchResult.RelationshipEstimates contains expected values") {
    assert(MatchResult.RelationshipEstimates.contains("PARENT_CHILD"))
    assert(MatchResult.RelationshipEstimates.contains("2ND_COUSIN"))
    assert(MatchResult.RelationshipEstimates.contains("DISTANT"))
  }

  // ============================================
  // Helpers
  // ============================================

  private def makeConsent(
                           biosampleRef: String = "at://did:plc:test/bio/1",
                           consentLevel: String = "FULL"
                         ): MatchConsent =
    MatchConsent(
      atUri = Some("at://did:plc:test/matchconsent/1"),
      meta = RecordMeta.initial,
      biosampleRef = biosampleRef,
      consentLevel = consentLevel,
      consentedAt = LocalDateTime.now()
    )

  private def makeRequest(
                            atUri: Option[String] = Some("at://did:plc:test/matchrequest/1"),
                            fromBiosampleRef: String = "at://did:plc:a/bio/1",
                            toBiosampleRef: String = "at://did:plc:b/bio/1",
                            status: String = "PENDING"
                          ): MatchRequest =
    MatchRequest(
      atUri = atUri,
      meta = RecordMeta.initial,
      fromBiosampleRef = fromBiosampleRef,
      toBiosampleRef = toBiosampleRef,
      status = status,
      createdAt = LocalDateTime.now()
    )
