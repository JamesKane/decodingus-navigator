package com.decodingus.service

import com.decodingus.client.DecodingUsClient
import com.decodingus.config.FeatureToggles
import com.decodingus.db.Transactor
import com.decodingus.repository.{MatchConsentRepository, MatchRequestRepository, MatchResultRepository}
import com.decodingus.service.EntityConversions.*
import com.decodingus.util.Logger
import com.decodingus.workspace.model.*

import java.time.LocalDateTime
import java.util.UUID
import scala.concurrent.{ExecutionContext, Future}

/**
 * Service layer for IBD match operations.
 *
 * Coordinates consent management, match suggestions from AppView,
 * match requests, and confirmed match results.
 */
class IbdMatchService(
                       transactor: Transactor,
                       consentRepo: MatchConsentRepository,
                       requestRepo: MatchRequestRepository,
                       resultRepo: MatchResultRepository
                     ):
  private val log = Logger[IbdMatchService]

  // ============================================
  // Feature Gate
  // ============================================

  def isEnabled: Boolean = FeatureToggles.ibdMatchingEnabled

  // ============================================
  // Consent Operations
  // ============================================

  def getConsent(biosampleId: UUID): Either[String, Option[MatchConsent]] =
    transactor.readOnly {
      consentRepo.findByBiosample(biosampleId).map { entity =>
        fromMatchConsentEntity(entity, localUri("biosample", entity.biosampleId))
      }
    }

  def hasConsent(biosampleId: UUID): Either[String, Boolean] =
    transactor.readOnly {
      consentRepo.findByBiosample(biosampleId).isDefined
    }

  def grantConsent(biosampleId: UUID, biosampleRef: String, consentLevel: ConsentLevel,
                   allowedMatchTypes: List[MatchType] = List(MatchType.Ibd),
                   minimumSegmentCm: Double = 7.0,
                   shareContactInfo: Boolean = false): Either[String, MatchConsent] =
    transactor.readWrite {
      // Upsert: update existing or create new
      consentRepo.findByBiosample(biosampleId) match
        case Some(existing) =>
          val updated = existing.copy(
            consentLevel = consentLevel.toDbString,
            allowedMatchTypes = allowedMatchTypes.map(_.toDbString),
            minimumSegmentCm = minimumSegmentCm,
            shareContactInfo = shareContactInfo
          )
          fromMatchConsentEntity(consentRepo.update(updated), biosampleRef)
        case None =>
          val entity = com.decodingus.repository.MatchConsentEntity.create(
            biosampleId = biosampleId,
            consentLevel = consentLevel.toDbString,
            allowedMatchTypes = allowedMatchTypes.map(_.toDbString),
            minimumSegmentCm = minimumSegmentCm,
            shareContactInfo = shareContactInfo
          )
          fromMatchConsentEntity(consentRepo.insert(entity), biosampleRef)
    }

  def revokeConsent(biosampleId: UUID): Either[String, Boolean] =
    transactor.readWrite {
      consentRepo.deleteByBiosample(biosampleId)
    }

  // ============================================
  // Match Suggestions (from AppView API)
  // ============================================

  def fetchSuggestions(biosampleRef: String, limit: Int = 20)(implicit ec: ExecutionContext): Future[Either[String, List[MatchSuggestion]]] =
    if !isEnabled then Future.successful(Left("IBD matching is not enabled"))
    else DecodingUsClient.getMatchSuggestions(biosampleRef, limit)

  def dismissSuggestion(suggestionId: String)(implicit ec: ExecutionContext): Future[Either[String, Boolean]] =
    if !isEnabled then Future.successful(Left("IBD matching is not enabled"))
    else DecodingUsClient.dismissMatchSuggestion(suggestionId)

  // ============================================
  // Match Request Operations
  // ============================================

  def getOutgoingRequests(biosampleRef: String): Either[String, List[MatchRequest]] =
    transactor.readOnly {
      requestRepo.findByFromBiosample(biosampleRef).map(fromMatchRequestEntity)
    }

  def getIncomingRequests(biosampleRef: String): Either[String, List[MatchRequest]] =
    transactor.readOnly {
      requestRepo.findByToBiosample(biosampleRef).map(fromMatchRequestEntity)
    }

  def getPendingRequests(biosampleRef: String): Either[String, List[MatchRequest]] =
    transactor.readOnly {
      requestRepo.findPending(biosampleRef).map(fromMatchRequestEntity)
    }

  def sendMatchRequest(fromBiosampleRef: String, toBiosampleRef: String,
                       requestType: RequestType = RequestType.Autosomal,
                       message: Option[String] = None,
                       sharedAncestorHint: Option[String] = None,
                       discoveryReason: Option[String] = None): Either[String, MatchRequest] =
    transactor.readWrite {
      val entity = com.decodingus.repository.MatchRequestEntity.create(
        fromBiosampleRef = fromBiosampleRef,
        toBiosampleRef = toBiosampleRef,
        requestType = requestType.toDbString,
        message = message,
        sharedAncestorHint = sharedAncestorHint,
        discoveryReason = discoveryReason
      )
      fromMatchRequestEntity(requestRepo.insert(entity))
    }

  def respondToRequest(requestId: UUID, accept: Boolean): Either[String, Boolean] =
    transactor.readWrite {
      val newStatus = if accept then RequestStatus.Accepted else RequestStatus.Declined
      requestRepo.updateRequestStatus(requestId, newStatus.toDbString)
    }

  def withdrawRequest(requestId: UUID): Either[String, Boolean] =
    transactor.readWrite {
      requestRepo.updateRequestStatus(requestId, RequestStatus.Withdrawn.toDbString)
    }

  // ============================================
  // Match Result Operations
  // ============================================

  def getMatchResults(biosampleId: UUID): Either[String, List[MatchResult]] =
    transactor.readOnly {
      val biosampleRef = localUri("biosample", biosampleId)
      resultRepo.findByBiosample(biosampleId).map(entity => fromMatchResultEntity(entity, biosampleRef))
    }

  def getMatchResultsAboveThreshold(biosampleId: UUID, minCm: Double): Either[String, List[MatchResult]] =
    transactor.readOnly {
      val biosampleRef = localUri("biosample", biosampleId)
      resultRepo.findByBiosampleAboveThreshold(biosampleId, minCm).map(entity => fromMatchResultEntity(entity, biosampleRef))
    }

  def getMatchResultCount(biosampleId: UUID): Either[String, Int] =
    transactor.readOnly {
      resultRepo.findByBiosample(biosampleId).size
    }
