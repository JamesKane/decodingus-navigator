package com.decodingus.workspace.model

import java.time.LocalDateTime

/**
 * A confirmed IBD match with another biosample.
 *
 * @param atUri                The AT URI of this match result record
 * @param meta                 Record metadata for tracking changes and sync
 * @param biosampleRef         AT URI of the local biosample
 * @param matchedBiosampleRef  AT URI of the matched biosample
 * @param matchedCitizenDid    DID of the matched citizen (if they consent to sharing)
 * @param relationshipEstimate Estimated relationship category
 * @param totalSharedCm        Total centiMorgans shared across all segments
 * @param longestSegmentCm     Length of the longest shared segment in cM
 * @param segmentCount         Number of shared segments
 * @param sharedSegments       Detailed segment information
 * @param xMatchSharedCm       cM shared on X chromosome
 * @param matchedAt            When this match was confirmed
 * @param attestationHash      SHA-256 hash of the match attestation
 */
case class MatchResult(
                        atUri: Option[String],
                        meta: RecordMeta,
                        biosampleRef: String,
                        matchedBiosampleRef: String,
                        matchedCitizenDid: Option[String] = None,
                        relationshipEstimate: Option[String] = None,
                        totalSharedCm: Double,
                        longestSegmentCm: Option[Double] = None,
                        segmentCount: Int,
                        sharedSegments: List[IbdSegment] = List.empty,
                        xMatchSharedCm: Option[Double] = None,
                        matchedAt: LocalDateTime,
                        attestationHash: Option[String] = None
                      )

object MatchResult:
  val RelationshipEstimates: Set[String] = Set(
    "PARENT_CHILD", "FULL_SIBLING", "HALF_SIBLING", "GRANDPARENT", "AUNT_UNCLE",
    "1ST_COUSIN", "1ST_COUSIN_1R", "2ND_COUSIN", "2ND_COUSIN_1R", "3RD_COUSIN",
    "4TH_COUSIN", "5TH_COUSIN", "DISTANT", "UNKNOWN"
  )
