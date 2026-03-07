package com.decodingus.workspace.model

import java.time.LocalDateTime

enum ConsentLevel:
  case Full, Anonymous, ProjectOnly

object ConsentLevel:
  def fromString(s: String): ConsentLevel = s match
    case "FULL" => Full
    case "ANONYMOUS" => Anonymous
    case "PROJECT_ONLY" => ProjectOnly
    case other => throw new IllegalArgumentException(s"Unknown consent level: $other")

  extension (cl: ConsentLevel)
    def toDbString: String = cl match
      case Full => "FULL"
      case Anonymous => "ANONYMOUS"
      case ProjectOnly => "PROJECT_ONLY"

enum MatchType:
  case Ibd, YStr, MtSequence, Autosomal

object MatchType:
  def fromString(s: String): MatchType = s match
    case "IBD" => Ibd
    case "Y_STR" => YStr
    case "MT_SEQUENCE" => MtSequence
    case "AUTOSOMAL" => Autosomal
    case other => throw new IllegalArgumentException(s"Unknown match type: $other")

  extension (mt: MatchType)
    def toDbString: String = mt match
      case Ibd => "IBD"
      case YStr => "Y_STR"
      case MtSequence => "MT_SEQUENCE"
      case Autosomal => "AUTOSOMAL"

/**
 * Consent record for IBD matching participation.
 * NSID: com.decodingus.atmosphere.matchConsent
 *
 * Presence enables matching for this biosample; deletion revokes consent.
 *
 * @param atUri             The AT URI of this consent record
 * @param meta              Record metadata for tracking changes and sync
 * @param biosampleRef      AT URI of the biosample for which consent is granted
 * @param consentLevel      Level of matching participation
 * @param allowedMatchTypes Types of matching allowed
 * @param minimumSegmentCm  Minimum segment size (cM) for matches to be shared
 * @param shareContactInfo  Whether to share contact information with matches
 * @param consentedAt       When consent was granted
 * @param expiresAt         Optional expiration date for consent
 */
case class MatchConsent(
                         atUri: Option[String],
                         meta: RecordMeta,
                         biosampleRef: String,
                         consentLevel: ConsentLevel,
                         allowedMatchTypes: List[MatchType] = List(MatchType.Ibd),
                         minimumSegmentCm: Double = 7.0,
                         shareContactInfo: Boolean = false,
                         consentedAt: LocalDateTime,
                         expiresAt: Option[LocalDateTime] = None
                       )
