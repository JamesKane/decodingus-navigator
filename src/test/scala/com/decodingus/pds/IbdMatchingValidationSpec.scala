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
    val mr = makeRequest(atUri = None, fromBiosampleRef = "", toBiosampleRef = "")
    val Left(errors) = PdsSyncValidation.validateMatchRequest(mr): @unchecked
    assertEquals(errors.size, 3)
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
  // Enum Completeness
  // ============================================

  test("ConsentLevel has all expected values") {
    assertEquals(ConsentLevel.values.length, 3)
    assertEquals(ConsentLevel.fromString("FULL"), ConsentLevel.Full)
    assertEquals(ConsentLevel.fromString("ANONYMOUS"), ConsentLevel.Anonymous)
    assertEquals(ConsentLevel.fromString("PROJECT_ONLY"), ConsentLevel.ProjectOnly)
  }

  test("MatchType has all expected values") {
    assertEquals(MatchType.values.length, 4)
    assertEquals(MatchType.fromString("IBD"), MatchType.Ibd)
    assertEquals(MatchType.fromString("Y_STR"), MatchType.YStr)
    assertEquals(MatchType.fromString("MT_SEQUENCE"), MatchType.MtSequence)
    assertEquals(MatchType.fromString("AUTOSOMAL"), MatchType.Autosomal)
  }

  test("RequestStatus has all expected values") {
    assertEquals(RequestStatus.values.length, 5)
    assertEquals(RequestStatus.fromString("PENDING"), RequestStatus.Pending)
    assertEquals(RequestStatus.fromString("ACCEPTED"), RequestStatus.Accepted)
    assertEquals(RequestStatus.fromString("DECLINED"), RequestStatus.Declined)
    assertEquals(RequestStatus.fromString("EXPIRED"), RequestStatus.Expired)
    assertEquals(RequestStatus.fromString("WITHDRAWN"), RequestStatus.Withdrawn)
  }

  test("RelationshipEstimate has all expected values") {
    assertEquals(RelationshipEstimate.values.length, 14)
    assertEquals(RelationshipEstimate.fromString("PARENT_CHILD"), RelationshipEstimate.ParentChild)
    assertEquals(RelationshipEstimate.fromString("2ND_COUSIN"), RelationshipEstimate.SecondCousin)
    assertEquals(RelationshipEstimate.fromString("DISTANT"), RelationshipEstimate.Distant)
  }

  test("ConsentLevel round-trips through toDbString/fromString") {
    for level <- ConsentLevel.values do
      assertEquals(ConsentLevel.fromString(level.toDbString), level)
  }

  test("RelationshipEstimate round-trips through toDbString/fromString") {
    for est <- RelationshipEstimate.values do
      assertEquals(RelationshipEstimate.fromString(est.toDbString), est)
  }

  // ============================================
  // Helpers
  // ============================================

  private def makeConsent(
                           biosampleRef: String = "at://did:plc:test/bio/1",
                           consentLevel: ConsentLevel = ConsentLevel.Full
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
                            toBiosampleRef: String = "at://did:plc:b/bio/1"
                          ): MatchRequest =
    MatchRequest(
      atUri = atUri,
      meta = RecordMeta.initial,
      fromBiosampleRef = fromBiosampleRef,
      toBiosampleRef = toBiosampleRef,
      status = RequestStatus.Pending,
      createdAt = LocalDateTime.now()
    )
