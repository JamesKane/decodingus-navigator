package com.decodingus.workspace.model

import java.time.LocalDateTime

/**
 * A request to initiate contact and IBD comparison with a genetic match.
 * NSID: com.decodingus.atmosphere.matchRequest
 *
 * @param atUri              The AT URI of this match request record
 * @param meta               Record metadata for tracking changes and sync
 * @param fromBiosampleRef   AT URI of the requesting biosample
 * @param toBiosampleRef     AT URI of the target biosample
 * @param status             Current status (PENDING, ACCEPTED, DECLINED, EXPIRED, WITHDRAWN)
 * @param requestType        Type of comparison (AUTOSOMAL, Y_CHROMOSOME, MT_DNA, FULL)
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
                         status: String,
                         requestType: String = "AUTOSOMAL",
                         message: Option[String] = None,
                         sharedAncestorHint: Option[String] = None,
                         discoveryReason: Option[String] = None,
                         createdAt: LocalDateTime,
                         expiresAt: Option[LocalDateTime] = None,
                         respondedAt: Option[LocalDateTime] = None
                       )

object MatchRequest:
  val Statuses: Set[String] = Set("PENDING", "ACCEPTED", "DECLINED", "EXPIRED", "WITHDRAWN")
  val RequestTypes: Set[String] = Set("AUTOSOMAL", "Y_CHROMOSOME", "MT_DNA", "FULL")
