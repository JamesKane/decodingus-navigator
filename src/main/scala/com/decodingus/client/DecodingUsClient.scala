package com.decodingus.client

import com.decodingus.model.LibraryStats
import com.decodingus.refgenome.model.{GenomeRegions, GenomeRegionsCodecs}
import com.decodingus.workspace.model.MatchSuggestion
import io.circe.*
import io.circe.generic.auto.*
import sttp.client3.*
import sttp.client3.circe.*

import scala.concurrent.{ExecutionContext, Future}

case class PdsRegistrationRequest(did: String, token: String, pdsUrl: String)

/**
 * Response from the match suggestions discovery endpoint.
 */
case class MatchSuggestionsResponse(
                                     data: List[MatchSuggestionDto],
                                     count: Int
                                   )

/**
 * DTO for a match suggestion from the AppView API.
 */
case class MatchSuggestionDto(
                               suggestionId: String,
                               biosampleRef: String,
                               matchedBiosampleRef: String,
                               matchedDid: Option[String],
                               matchedLabel: String,
                               score: Double,
                               reasonType: String,
                               reasonDetail: Option[String],
                               populationOverlap: Option[Double]
                             )

/**
 * Response from the population overlap endpoint.
 */
case class PopulationOverlapResponse(overlap: Double)

/**
 * Information about a sequencing instrument and its associated lab.
 *
 * @param instrumentId Unique instrument identifier (matches @RG PU/PM fields)
 * @param labName      Name of the sequencing facility
 * @param isD2c        Whether this is a direct-to-consumer lab
 * @param manufacturer Instrument manufacturer (e.g., "Illumina", "PacBio")
 * @param model        Instrument model (e.g., "NovaSeq 6000", "Sequel II")
 * @param websiteUrl   Lab website URL if available
 */
case class SequencerLabInfo(
                             instrumentId: String,
                             labName: String,
                             isD2c: Boolean,
                             manufacturer: Option[String] = None,
                             model: Option[String] = None,
                             websiteUrl: Option[String] = None
                           )

/**
 * Response from the lab-instruments endpoint.
 */
case class SequencerLabInstrumentsResponse(
                                            data: List[SequencerLabInfo],
                                            count: Int
                                          )

object DecodingUsClient {

  private val backend = HttpClientFutureBackend()
  private val BaseUrl = uri"https://decoding-us.com/api/v1"

  // In-memory cache for lab instruments (refreshed on demand)
  @volatile private var labInstrumentsCache: Option[Map[String, SequencerLabInfo]] = None

  /**
   * Registers the user's PDS with the DecodingUs platform.
   *
   * @param did    The user's DID.
   * @param token  The authentication token (R_Token).
   * @param pdsUrl The URL of the PDS.
   * @param ec     Execution context.
   * @return A Future completing on success.
   */
  def registerPds(did: String, token: String, pdsUrl: String)(implicit ec: ExecutionContext): Future[Unit] = {
    val request = basicRequest
      .post(BaseUrl.addPath("registerPDS"))
      .body(PdsRegistrationRequest(did, token, pdsUrl))
      .response(asString)

    request.send(backend).flatMap { response =>
      if (response.code.isSuccess) {
        println(s"Successfully registered PDS for $did")
        Future.successful(())
      } else {
        Future.failed(new RuntimeException(s"PDS Registration failed: ${response.code} ${response.body}"))
      }
    }
  }

  /**
   * Stubs the retrieval of a Biosample ID from the DecodingUs platform.
   * This ID uniquely identifies a specific sequencing event and alignment version for a donor.
   *
   * @param userId       The user's ID.
   * @param libraryStats Metadata about the library (sample name, reference, platform).
   * @param ec           Execution context.
   * @return A Future containing the Biosample ID string.
   */
  def resolveBiosampleId(userId: String, libraryStats: LibraryStats)(implicit ec: ExecutionContext): Future[String] = {
    // In a real implementation, this would POST metadata to the platform to find or register the biosample.
    // val request = basicRequest
    //   .post(BaseUrl.addPath("biosamples", "resolve"))
    //   .body(...)
    //   .response(asJson[String])

    Future {
      // simulating network delay
      Thread.sleep(500)

      // Mock logic to generate a consistent ID based on input for testing, 
      // or just a random one if it were a real 'new' sample.
      // For stubbing purposes, we'll make it look like a UUID but deterministic for the same sample name.
      val seed = s"$userId-${libraryStats.sampleName}-${libraryStats.referenceBuild}"
      java.util.UUID.nameUUIDFromBytes(seed.getBytes).toString
    }
  }

  /**
   * Fetches all lab-instrument associations from the API.
   *
   * @param ec Execution context
   * @return Future containing the list of lab instrument info
   */
  def getLabInstruments()(implicit ec: ExecutionContext): Future[List[SequencerLabInfo]] = {
    val request = basicRequest
      .get(BaseUrl.addPath("sequencer", "lab-instruments"))
      .response(asJson[SequencerLabInstrumentsResponse])

    request.send(backend).flatMap { response =>
      response.body match {
        case Right(labResponse) =>
          // Update cache
          labInstrumentsCache = Some(labResponse.data.map(info => info.instrumentId -> info).toMap)
          Future.successful(labResponse.data)
        case Left(error) =>
          Future.failed(new RuntimeException(s"Failed to fetch lab instruments: $error"))
      }
    }
  }

  /**
   * Looks up lab information for a given instrument ID.
   * Uses cached data if available, otherwise fetches from API.
   *
   * @param instrumentId The instrument identifier to look up
   * @param ec           Execution context
   * @return Future containing optional lab info if found
   */
  def lookupLabByInstrument(instrumentId: String)(implicit ec: ExecutionContext): Future[Option[SequencerLabInfo]] = {
    labInstrumentsCache match {
      case Some(cache) =>
        Future.successful(cache.get(instrumentId))
      case None =>
        getLabInstruments().map { instruments =>
          instruments.find(_.instrumentId == instrumentId)
        }
    }
  }

  /**
   * Clears the lab instruments cache, forcing a refresh on next lookup.
   */
  def clearLabInstrumentsCache(): Unit = {
    labInstrumentsCache = None
  }

  /**
   * Fetches genome region metadata for a reference build.
   * Returns centromeres, telomeres, cytobands, and Y-specific region annotations.
   *
   * @param build The reference genome build (GRCh38, GRCh37, CHM13v2)
   * @param ec    Execution context
   * @return Future containing Either error message or GenomeRegions data
   */
  def getGenomeRegions(build: String)(implicit ec: ExecutionContext): Future[Either[String, GenomeRegions]] = {
    import GenomeRegionsCodecs.given

    val request = basicRequest
      .get(BaseUrl.addPath("genome-regions", build))
      .response(asJson[GenomeRegions])

    request.send(backend).map { response =>
      response.body match {
        case Right(regions) =>
          Right(regions)
        case Left(error) =>
          Left(s"Failed to fetch genome regions for $build: ${error.getMessage}")
      }
    }.recover { case e: Exception =>
      Left(s"Network error fetching genome regions: ${e.getMessage}")
    }
  }

  // ============================================
  // IBD Match Discovery API
  // ============================================

  /**
   * Fetches match suggestions for a biosample from the AppView discovery engine.
   *
   * @param biosampleRef AT URI of the biosample to get suggestions for
   * @param limit        Maximum number of suggestions to return
   * @param ec           Execution context
   * @return Future containing either error message or list of suggestions
   */
  def getMatchSuggestions(biosampleRef: String, limit: Int = 20)(implicit ec: ExecutionContext): Future[Either[String, List[MatchSuggestion]]] = {
    val request = basicRequest
      .get(BaseUrl.addPath("ibd", "suggestions")
        .addParam("biosampleRef", biosampleRef)
        .addParam("limit", limit.toString))
      .response(asJson[MatchSuggestionsResponse])

    request.send(backend).map { response =>
      response.body match {
        case Right(suggestionsResp) =>
          Right(suggestionsResp.data.map(dto => MatchSuggestion(
            suggestionId = dto.suggestionId,
            biosampleRef = dto.biosampleRef,
            matchedBiosampleRef = dto.matchedBiosampleRef,
            matchedDid = dto.matchedDid,
            matchedLabel = dto.matchedLabel,
            score = dto.score,
            reasonType = dto.reasonType,
            reasonDetail = dto.reasonDetail,
            populationOverlap = dto.populationOverlap
          )))
        case Left(error) =>
          Left(s"Failed to fetch match suggestions: ${error.getMessage}")
      }
    }.recover { case e: Exception =>
      Left(s"Network error fetching match suggestions: ${e.getMessage}")
    }
  }

  /**
   * Dismisses a match suggestion so it won't appear again.
   *
   * @param suggestionId The suggestion ID to dismiss
   * @param ec           Execution context
   * @return Future containing either error message or success
   */
  def dismissMatchSuggestion(suggestionId: String)(implicit ec: ExecutionContext): Future[Either[String, Boolean]] = {
    val request = basicRequest
      .post(BaseUrl.addPath("ibd", "suggestions", suggestionId, "dismiss"))
      .response(asString)

    request.send(backend).map { response =>
      if (response.code.isSuccess) Right(true)
      else Left(s"Failed to dismiss suggestion: ${response.code}")
    }.recover { case e: Exception =>
      Left(s"Network error dismissing suggestion: ${e.getMessage}")
    }
  }

  /**
   * Gets the estimated population overlap between two biosamples.
   *
   * @param biosampleRef1 AT URI of the first biosample
   * @param biosampleRef2 AT URI of the second biosample
   * @param ec            Execution context
   * @return Future containing either error message or overlap percentage (0.0-100.0)
   */
  def getPopulationOverlap(biosampleRef1: String, biosampleRef2: String)(implicit ec: ExecutionContext): Future[Either[String, Double]] = {
    val request = basicRequest
      .get(BaseUrl.addPath("ibd", "population-overlap")
        .addParam("biosample1", biosampleRef1)
        .addParam("biosample2", biosampleRef2))
      .response(asJson[PopulationOverlapResponse])

    request.send(backend).map { response =>
      response.body match {
        case Right(overlapResp) => Right(overlapResp.overlap)
        case Left(error) => Left(s"Failed to fetch population overlap: ${error.getMessage}")
      }
    }.recover { case e: Exception =>
      Left(s"Network error fetching population overlap: ${e.getMessage}")
    }
  }
}
