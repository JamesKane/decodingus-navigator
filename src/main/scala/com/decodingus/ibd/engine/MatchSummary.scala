package com.decodingus.ibd.engine

import com.decodingus.ibd.crypto.IbdCryptoService
import com.decodingus.workspace.model.{IbdSegment, RelationshipEstimate}

import java.nio.charset.StandardCharsets

/**
 * Summary of an IBD comparison between two individuals.
 *
 * @param totalSharedCm        Total centiMorgans shared across all segments
 * @param segmentCount         Number of shared segments
 * @param longestSegmentCm     Length of the longest shared segment
 * @param segments             All detected IBD segments
 * @param relationshipEstimate Estimated relationship category
 * @param summaryHash          SHA-256 hash of canonical summary for attestation verification
 */
case class MatchSummary(
                         totalSharedCm: Double,
                         segmentCount: Int,
                         longestSegmentCm: Double,
                         segments: List[IbdSegment],
                         relationshipEstimate: RelationshipEstimate,
                         summaryHash: String
                       )

object MatchSummary:

  /**
   * Compute a match summary from detected IBD segments.
   * The summary hash is deterministic — two parties computing the same
   * segments will produce the same hash.
   */
  def fromSegments(segments: List[IbdSegment]): MatchSummary =
    if segments.isEmpty then
      MatchSummary(
        totalSharedCm = 0.0,
        segmentCount = 0,
        longestSegmentCm = 0.0,
        segments = Nil,
        relationshipEstimate = RelationshipEstimate.Unknown,
        summaryHash = IbdCryptoService.sha256Hex("no_segments".getBytes(StandardCharsets.UTF_8))
      )
    else
      val totalCm = roundCm(segments.map(_.lengthCm).sum)
      val longestCm = roundCm(segments.map(_.lengthCm).max)
      val relationship = RelationshipEstimator.estimate(totalCm)
      val hash = computeHash(segments)

      MatchSummary(
        totalSharedCm = totalCm,
        segmentCount = segments.size,
        longestSegmentCm = longestCm,
        segments = segments,
        relationshipEstimate = relationship,
        summaryHash = hash
      )

  /**
   * Compute a deterministic SHA-256 hash of the match summary.
   *
   * The canonical form sorts segments by (chromosome numerically, startPosition)
   * and rounds cM to 2 decimal places to ensure both parties produce identical hashes
   * despite potential floating-point differences.
   */
  def computeHash(segments: List[IbdSegment]): String =
    val canonical = segments
      .sortBy(s => (chrSortKey(s.chromosome), s.startPosition))
      .map { s =>
        s"${s.chromosome}:${s.startPosition}-${s.endPosition}:${roundCm(s.lengthCm)}"
      }
      .mkString("|")

    IbdCryptoService.sha256Hex(canonical.getBytes(StandardCharsets.UTF_8))

  private def roundCm(cm: Double): Double =
    math.round(cm * 100.0) / 100.0

  private def chrSortKey(chr: String): Int =
    chr.toIntOption.getOrElse(chr match
      case "X" => 23
      case "Y" => 24
      case _ => 99
    )
