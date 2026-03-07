package com.decodingus.workspace.model

import java.time.LocalDateTime

/**
 * Consent record for IBD matching participation.
 * NSID: com.decodingus.atmosphere.matchConsent
 *
 * Presence enables matching for this biosample; deletion revokes consent.
 *
 * @param atUri             The AT URI of this consent record
 * @param meta              Record metadata for tracking changes and sync
 * @param biosampleRef      AT URI of the biosample for which consent is granted
 * @param consentLevel      Level of matching participation (FULL, ANONYMOUS, PROJECT_ONLY)
 * @param allowedMatchTypes Types of matching allowed (IBD, Y_STR, MT_SEQUENCE, AUTOSOMAL)
 * @param minimumSegmentCm  Minimum segment size (cM) for matches to be shared
 * @param shareContactInfo  Whether to share contact information with matches
 * @param consentedAt       When consent was granted
 * @param expiresAt         Optional expiration date for consent
 */
case class MatchConsent(
                         atUri: Option[String],
                         meta: RecordMeta,
                         biosampleRef: String,
                         consentLevel: String,
                         allowedMatchTypes: List[String] = List("IBD"),
                         minimumSegmentCm: Double = 7.0,
                         shareContactInfo: Boolean = false,
                         consentedAt: LocalDateTime,
                         expiresAt: Option[LocalDateTime] = None
                       )

object MatchConsent:
  val ConsentLevels: Set[String] = Set("FULL", "ANONYMOUS", "PROJECT_ONLY")
  val MatchTypes: Set[String] = Set("IBD", "Y_STR", "MT_SEQUENCE", "AUTOSOMAL")
