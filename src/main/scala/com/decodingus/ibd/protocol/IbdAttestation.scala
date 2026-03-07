package com.decodingus.ibd.protocol

import com.decodingus.ibd.crypto.IbdCryptoService
import com.decodingus.ibd.engine.MatchSummary
import com.decodingus.workspace.model.RelationshipEstimate
import io.circe.*

import java.nio.charset.StandardCharsets
import java.security.{KeyPair, PrivateKey, PublicKey}
import java.time.Instant
import java.util.Base64

/**
 * An attestation that a Navigator instance computed an IBD match.
 *
 * Both parties produce an attestation; when both hashes agree,
 * the AppView "stamps" the match as confirmed.
 *
 * @param matchRequestUri    AT URI of the match request that triggered this comparison
 * @param sessionId          Unique session ID for the comparison exchange
 * @param attestingDid       DID of the citizen whose Navigator produced this attestation
 * @param attestingSampleRef AT URI of the biosample that was compared
 * @param partnerSampleRef   AT URI of the partner biosample
 * @param totalSharedCm      Total centiMorgans shared
 * @param segmentCount       Number of shared segments
 * @param longestSegmentCm   Longest segment in cM
 * @param relationshipEstimate Estimated relationship category
 * @param summaryHash        SHA-256 of the canonical match summary
 * @param partnerSummaryHash SHA-256 hash received from the partner
 * @param signature          Ed25519 signature over the canonical attestation content
 * @param signingPublicKey   Base64-encoded Ed25519 public key for verification
 * @param attestedAt         When this attestation was created
 */
case class IbdAttestation(
                           matchRequestUri: String,
                           sessionId: String,
                           attestingDid: String,
                           attestingSampleRef: String,
                           partnerSampleRef: String,
                           totalSharedCm: Double,
                           segmentCount: Int,
                           longestSegmentCm: Double,
                           relationshipEstimate: RelationshipEstimate,
                           summaryHash: String,
                           partnerSummaryHash: String,
                           signature: String,
                           signingPublicKey: String,
                           attestedAt: Instant
                         )

object IbdAttestation:

  /**
   * Create and sign an attestation from a match summary.
   *
   * The signature covers the canonical attestation string to prevent tampering.
   */
  def create(
              matchRequestUri: String,
              sessionId: String,
              attestingDid: String,
              attestingSampleRef: String,
              partnerSampleRef: String,
              summary: MatchSummary,
              partnerSummaryHash: String,
              signingKeyPair: KeyPair
            ): IbdAttestation =
    val attestedAt = Instant.now()
    val canonical = canonicalString(
      matchRequestUri, sessionId, attestingDid, attestingSampleRef, partnerSampleRef,
      summary.totalSharedCm, summary.segmentCount, summary.longestSegmentCm,
      summary.relationshipEstimate, summary.summaryHash, partnerSummaryHash, attestedAt
    )

    val signatureBytes = IbdCryptoService.signAttestation(
      canonical.getBytes(StandardCharsets.UTF_8),
      signingKeyPair.getPrivate
    )

    IbdAttestation(
      matchRequestUri = matchRequestUri,
      sessionId = sessionId,
      attestingDid = attestingDid,
      attestingSampleRef = attestingSampleRef,
      partnerSampleRef = partnerSampleRef,
      totalSharedCm = summary.totalSharedCm,
      segmentCount = summary.segmentCount,
      longestSegmentCm = summary.longestSegmentCm,
      relationshipEstimate = summary.relationshipEstimate,
      summaryHash = summary.summaryHash,
      partnerSummaryHash = partnerSummaryHash,
      signature = Base64.getEncoder.encodeToString(signatureBytes),
      signingPublicKey = IbdCryptoService.encodePublicKey(signingKeyPair.getPublic),
      attestedAt = attestedAt
    )

  /**
   * Verify an attestation signature using the embedded public key.
   */
  def verify(attestation: IbdAttestation): Boolean =
    val canonical = canonicalString(
      attestation.matchRequestUri, attestation.sessionId, attestation.attestingDid,
      attestation.attestingSampleRef, attestation.partnerSampleRef,
      attestation.totalSharedCm, attestation.segmentCount, attestation.longestSegmentCm,
      attestation.relationshipEstimate, attestation.summaryHash,
      attestation.partnerSummaryHash, attestation.attestedAt
    )

    val publicKey = IbdCryptoService.decodeEd25519PublicKey(attestation.signingPublicKey)
    val signatureBytes = Base64.getDecoder.decode(attestation.signature)
    IbdCryptoService.verifyAttestation(
      canonical.getBytes(StandardCharsets.UTF_8),
      signatureBytes,
      publicKey
    )

  private def canonicalString(
                               matchRequestUri: String, sessionId: String, attestingDid: String,
                               attestingSampleRef: String, partnerSampleRef: String,
                               totalSharedCm: Double, segmentCount: Int, longestSegmentCm: Double,
                               relationshipEstimate: RelationshipEstimate, summaryHash: String,
                               partnerSummaryHash: String, attestedAt: Instant
                             ): String =
    s"$matchRequestUri|$sessionId|$attestingDid|$attestingSampleRef|$partnerSampleRef|" +
      s"$totalSharedCm|$segmentCount|$longestSegmentCm|${relationshipEstimate.toDbString}|" +
      s"$summaryHash|$partnerSummaryHash|${attestedAt.getEpochSecond}"

  given Encoder[IbdAttestation] = Encoder.instance { a =>
    Json.obj(
      "matchRequestUri" -> Json.fromString(a.matchRequestUri),
      "sessionId" -> Json.fromString(a.sessionId),
      "attestingDid" -> Json.fromString(a.attestingDid),
      "attestingSampleRef" -> Json.fromString(a.attestingSampleRef),
      "partnerSampleRef" -> Json.fromString(a.partnerSampleRef),
      "totalSharedCm" -> Json.fromDoubleOrNull(a.totalSharedCm),
      "segmentCount" -> Json.fromInt(a.segmentCount),
      "longestSegmentCm" -> Json.fromDoubleOrNull(a.longestSegmentCm),
      "relationshipEstimate" -> Json.fromString(a.relationshipEstimate.toDbString),
      "summaryHash" -> Json.fromString(a.summaryHash),
      "partnerSummaryHash" -> Json.fromString(a.partnerSummaryHash),
      "signature" -> Json.fromString(a.signature),
      "signingPublicKey" -> Json.fromString(a.signingPublicKey),
      "attestedAt" -> Json.fromString(a.attestedAt.toString)
    )
  }

  given Decoder[IbdAttestation] = Decoder.instance { c =>
    for
      matchRequestUri <- c.get[String]("matchRequestUri")
      sessionId <- c.get[String]("sessionId")
      attestingDid <- c.get[String]("attestingDid")
      attestingSampleRef <- c.get[String]("attestingSampleRef")
      partnerSampleRef <- c.get[String]("partnerSampleRef")
      totalSharedCm <- c.get[Double]("totalSharedCm")
      segmentCount <- c.get[Int]("segmentCount")
      longestSegmentCm <- c.get[Double]("longestSegmentCm")
      relationshipEstimate <- c.get[String]("relationshipEstimate").map(RelationshipEstimate.fromString)
      summaryHash <- c.get[String]("summaryHash")
      partnerSummaryHash <- c.get[String]("partnerSummaryHash")
      signature <- c.get[String]("signature")
      signingPublicKey <- c.get[String]("signingPublicKey")
      attestedAt <- c.get[String]("attestedAt").map(Instant.parse)
    yield IbdAttestation(
      matchRequestUri, sessionId, attestingDid, attestingSampleRef, partnerSampleRef,
      totalSharedCm, segmentCount, longestSegmentCm, relationshipEstimate,
      summaryHash, partnerSummaryHash, signature, signingPublicKey, attestedAt
    )
  }
