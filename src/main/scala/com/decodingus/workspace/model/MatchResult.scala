package com.decodingus.workspace.model

import java.time.LocalDateTime

enum RelationshipEstimate:
  case ParentChild, FullSibling, HalfSibling, Grandparent, AuntUncle
  case FirstCousin, FirstCousinOnceRemoved
  case SecondCousin, SecondCousinOnceRemoved
  case ThirdCousin, FourthCousin, FifthCousin
  case Distant, Unknown

object RelationshipEstimate:
  def fromString(s: String): RelationshipEstimate = s match
    case "PARENT_CHILD" => ParentChild
    case "FULL_SIBLING" => FullSibling
    case "HALF_SIBLING" => HalfSibling
    case "GRANDPARENT" => Grandparent
    case "AUNT_UNCLE" => AuntUncle
    case "1ST_COUSIN" => FirstCousin
    case "1ST_COUSIN_1R" => FirstCousinOnceRemoved
    case "2ND_COUSIN" => SecondCousin
    case "2ND_COUSIN_1R" => SecondCousinOnceRemoved
    case "3RD_COUSIN" => ThirdCousin
    case "4TH_COUSIN" => FourthCousin
    case "5TH_COUSIN" => FifthCousin
    case "DISTANT" => Distant
    case "UNKNOWN" => Unknown
    case other => throw new IllegalArgumentException(s"Unknown relationship estimate: $other")

  extension (re: RelationshipEstimate)
    def toDbString: String = re match
      case ParentChild => "PARENT_CHILD"
      case FullSibling => "FULL_SIBLING"
      case HalfSibling => "HALF_SIBLING"
      case Grandparent => "GRANDPARENT"
      case AuntUncle => "AUNT_UNCLE"
      case FirstCousin => "1ST_COUSIN"
      case FirstCousinOnceRemoved => "1ST_COUSIN_1R"
      case SecondCousin => "2ND_COUSIN"
      case SecondCousinOnceRemoved => "2ND_COUSIN_1R"
      case ThirdCousin => "3RD_COUSIN"
      case FourthCousin => "4TH_COUSIN"
      case FifthCousin => "5TH_COUSIN"
      case Distant => "DISTANT"
      case Unknown => "UNKNOWN"

    def label: String = re match
      case ParentChild => "Parent/Child"
      case FullSibling => "Full Sibling"
      case HalfSibling => "Half Sibling"
      case Grandparent => "Grandparent/Grandchild"
      case AuntUncle => "Aunt/Uncle/Niece/Nephew"
      case FirstCousin => "1st Cousin"
      case FirstCousinOnceRemoved => "1st Cousin Once Removed"
      case SecondCousin => "2nd Cousin"
      case SecondCousinOnceRemoved => "2nd Cousin Once Removed"
      case ThirdCousin => "3rd Cousin"
      case FourthCousin => "4th Cousin"
      case FifthCousin => "5th Cousin"
      case Distant => "Distant Relative"
      case Unknown => "Unknown"

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
                        relationshipEstimate: Option[RelationshipEstimate] = None,
                        totalSharedCm: Double,
                        longestSegmentCm: Option[Double] = None,
                        segmentCount: Int,
                        sharedSegments: List[IbdSegment] = List.empty,
                        xMatchSharedCm: Option[Double] = None,
                        matchedAt: LocalDateTime,
                        attestationHash: Option[String] = None
                      )
