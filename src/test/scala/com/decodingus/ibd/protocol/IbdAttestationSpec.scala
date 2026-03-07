package com.decodingus.ibd.protocol

import com.decodingus.ibd.crypto.IbdCryptoService
import com.decodingus.ibd.engine.MatchSummary
import com.decodingus.workspace.model.{IbdSegment, RelationshipEstimate}
import io.circe.parser.decode
import io.circe.syntax.*
import munit.FunSuite

class IbdAttestationSpec extends FunSuite:

  private def makeSummary(): MatchSummary =
    val segments = List(
      IbdSegment("1", 5000000, 10000000, 15.23),
      IbdSegment("3", 20000000, 30000000, 28.76)
    )
    MatchSummary.fromSegments(segments)

  test("create and verify attestation") {
    val signingKeyPair = IbdCryptoService.generateEd25519KeyPair()
    val summary = makeSummary()

    val attestation = IbdAttestation.create(
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "session-001",
      attestingDid = "did:plc:a",
      attestingSampleRef = "at://did:plc:a/bio/1",
      partnerSampleRef = "at://did:plc:b/bio/1",
      summary = summary,
      partnerSummaryHash = summary.summaryHash,
      signingKeyPair = signingKeyPair
    )

    assertEquals(attestation.attestingDid, "did:plc:a")
    assertEquals(attestation.totalSharedCm, summary.totalSharedCm)
    assertEquals(attestation.segmentCount, summary.segmentCount)
    assertEquals(attestation.summaryHash, summary.summaryHash)
    assertEquals(attestation.partnerSummaryHash, summary.summaryHash)
    assert(attestation.signature.nonEmpty)
    assert(attestation.signingPublicKey.nonEmpty)

    assert(IbdAttestation.verify(attestation), "Attestation should verify with its own key")
  }

  test("tampered attestation fails verification") {
    val signingKeyPair = IbdCryptoService.generateEd25519KeyPair()
    val summary = makeSummary()

    val attestation = IbdAttestation.create(
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "session-001",
      attestingDid = "did:plc:a",
      attestingSampleRef = "at://did:plc:a/bio/1",
      partnerSampleRef = "at://did:plc:b/bio/1",
      summary = summary,
      partnerSummaryHash = summary.summaryHash,
      signingKeyPair = signingKeyPair
    )

    // Tamper with the total shared cM
    val tampered = attestation.copy(totalSharedCm = 999.0)
    assert(!IbdAttestation.verify(tampered), "Tampered attestation should fail verification")
  }

  test("attestation signed by different key fails verification") {
    val signingKeyPair1 = IbdCryptoService.generateEd25519KeyPair()
    val signingKeyPair2 = IbdCryptoService.generateEd25519KeyPair()
    val summary = makeSummary()

    val attestation = IbdAttestation.create(
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "session-001",
      attestingDid = "did:plc:a",
      attestingSampleRef = "at://did:plc:a/bio/1",
      partnerSampleRef = "at://did:plc:b/bio/1",
      summary = summary,
      partnerSummaryHash = summary.summaryHash,
      signingKeyPair = signingKeyPair1
    )

    // Replace the public key with a different one
    val wrongKey = attestation.copy(
      signingPublicKey = IbdCryptoService.encodePublicKey(signingKeyPair2.getPublic)
    )
    assert(!IbdAttestation.verify(wrongKey), "Wrong public key should fail verification")
  }

  test("attestation JSON round trip") {
    val signingKeyPair = IbdCryptoService.generateEd25519KeyPair()
    val summary = makeSummary()

    val attestation = IbdAttestation.create(
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "session-001",
      attestingDid = "did:plc:a",
      attestingSampleRef = "at://did:plc:a/bio/1",
      partnerSampleRef = "at://did:plc:b/bio/1",
      summary = summary,
      partnerSummaryHash = summary.summaryHash,
      signingKeyPair = signingKeyPair
    )

    val json = attestation.asJson.noSpaces
    val decoded = decode[IbdAttestation](json)

    assert(decoded.isRight, s"Failed to decode: ${decoded.left.getOrElse("")}")
    val roundTripped = decoded.toOption.get

    assertEquals(roundTripped.matchRequestUri, attestation.matchRequestUri)
    assertEquals(roundTripped.sessionId, attestation.sessionId)
    assertEquals(roundTripped.attestingDid, attestation.attestingDid)
    assertEquals(roundTripped.totalSharedCm, attestation.totalSharedCm)
    assertEquals(roundTripped.summaryHash, attestation.summaryHash)
    assertEquals(roundTripped.signature, attestation.signature)
    assertEquals(roundTripped.relationshipEstimate, attestation.relationshipEstimate)

    // Verify still works after round trip
    assert(IbdAttestation.verify(roundTripped), "Round-tripped attestation should verify")
  }

  test("two parties produce matching hashes for same segments") {
    val aliceKey = IbdCryptoService.generateEd25519KeyPair()
    val bobKey = IbdCryptoService.generateEd25519KeyPair()
    val summary = makeSummary()

    val aliceAttestation = IbdAttestation.create(
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "session-001",
      attestingDid = "did:plc:a",
      attestingSampleRef = "at://did:plc:a/bio/1",
      partnerSampleRef = "at://did:plc:b/bio/1",
      summary = summary,
      partnerSummaryHash = summary.summaryHash,
      signingKeyPair = aliceKey
    )

    val bobAttestation = IbdAttestation.create(
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "session-001",
      attestingDid = "did:plc:b",
      attestingSampleRef = "at://did:plc:b/bio/1",
      partnerSampleRef = "at://did:plc:a/bio/1",
      summary = summary,
      partnerSummaryHash = summary.summaryHash,
      signingKeyPair = bobKey
    )

    // Both should have the same summary hash
    assertEquals(aliceAttestation.summaryHash, bobAttestation.summaryHash)
    // But different signatures (different keys)
    assert(aliceAttestation.signature != bobAttestation.signature)
    // Both should verify independently
    assert(IbdAttestation.verify(aliceAttestation))
    assert(IbdAttestation.verify(bobAttestation))
  }

  test("attestation includes relationship estimate") {
    val signingKeyPair = IbdCryptoService.generateEd25519KeyPair()
    val summary = makeSummary()

    val attestation = IbdAttestation.create(
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "session-001",
      attestingDid = "did:plc:a",
      attestingSampleRef = "at://did:plc:a/bio/1",
      partnerSampleRef = "at://did:plc:b/bio/1",
      summary = summary,
      partnerSummaryHash = summary.summaryHash,
      signingKeyPair = signingKeyPair
    )

    assertEquals(attestation.relationshipEstimate, summary.relationshipEstimate)
  }
