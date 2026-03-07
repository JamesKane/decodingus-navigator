package com.decodingus.workspace.model

import java.time.LocalDateTime

enum RequestStatus:
  case Pending, Accepted, Declined, Expired, Withdrawn

object RequestStatus:
  def fromString(s: String): RequestStatus = s match
    case "PENDING" => Pending
    case "ACCEPTED" => Accepted
    case "DECLINED" => Declined
    case "EXPIRED" => Expired
    case "WITHDRAWN" => Withdrawn
    case other => throw new IllegalArgumentException(s"Unknown request status: $other")

  extension (rs: RequestStatus)
    def toDbString: String = rs match
      case Pending => "PENDING"
      case Accepted => "ACCEPTED"
      case Declined => "DECLINED"
      case Expired => "EXPIRED"
      case Withdrawn => "WITHDRAWN"

enum RequestType:
  case Autosomal, YChromosome, MtDna, Full

object RequestType:
  def fromString(s: String): RequestType = s match
    case "AUTOSOMAL" => Autosomal
    case "Y_CHROMOSOME" => YChromosome
    case "MT_DNA" => MtDna
    case "FULL" => Full
    case other => throw new IllegalArgumentException(s"Unknown request type: $other")

  extension (rt: RequestType)
    def toDbString: String = rt match
      case Autosomal => "AUTOSOMAL"
      case YChromosome => "Y_CHROMOSOME"
      case MtDna => "MT_DNA"
      case Full => "FULL"

/**
 * A request to initiate contact and IBD comparison with a genetic match.
 * NSID: com.decodingus.atmosphere.matchRequest
 *
 * @param atUri              The AT URI of this match request record
 * @param meta               Record metadata for tracking changes and sync
 * @param fromBiosampleRef   AT URI of the requesting biosample
 * @param toBiosampleRef     AT URI of the target biosample
 * @param status             Current status
 * @param requestType        Type of comparison
 * @param message            Optional message to the match
 * @param sharedAncestorHint Suspected common ancestor or family line
 * @param discoveryReason    Why this match was suggested (JSON-serializable)
 * @param createdAt          When the request was created
 * @param expiresAt          When this request expires if not responded to
 * @param respondedAt        When the request was responded to
 */
case class MatchRequest(
                         atUri: Option[String],
                         meta: RecordMeta,
                         fromBiosampleRef: String,
                         toBiosampleRef: String,
                         status: RequestStatus,
                         requestType: RequestType = RequestType.Autosomal,
                         message: Option[String] = None,
                         sharedAncestorHint: Option[String] = None,
                         discoveryReason: Option[String] = None,
                         createdAt: LocalDateTime,
                         expiresAt: Option[LocalDateTime] = None,
                         respondedAt: Option[LocalDateTime] = None
                       )
