package com.decodingus.client

import com.decodingus.model.LibraryStats
import com.decodingus.refgenome.model.{GenomeRegions, GenomeRegionsCodecs}
import sttp.client3._
import sttp.client3.circe._
import io.circe.generic.auto._

import scala.concurrent.{ExecutionContext, Future}

case class PdsRegistrationRequest(did: String, token: String, pdsUrl: String)

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
   * @return       A Future completing on success.
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
   * @return             A Future containing the Biosample ID string.
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
}
