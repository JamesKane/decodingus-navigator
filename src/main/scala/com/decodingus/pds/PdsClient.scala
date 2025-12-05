package com.decodingus.pds

import com.decodingus.auth.User
import com.decodingus.model.{ContigSummary, CoverageSummary, LibraryStats, WgsMetrics}
import sttp.client3._
import sttp.client3.circe._
import io.circe.{Encoder, Json}
import io.circe.generic.semiauto.deriveEncoder
import io.circe.syntax._

import scala.concurrent.{ExecutionContext, Future}

object PdsClient {

  implicit val libraryStatsEncoder: Encoder[LibraryStats] = deriveEncoder
  implicit val wgsMetricsEncoder: Encoder[WgsMetrics] = deriveEncoder
  implicit val contigSummaryEncoder: Encoder[ContigSummary] = deriveEncoder
  implicit val coverageSummaryEncoder: Encoder[CoverageSummary] = deriveEncoder

  // Using the Future backend
  private val backend = HttpClientFutureBackend()

  // Placeholder endpoint - would normally be in config
  // In a real scenario, this URL would handle the ingestion of the JSON payload.
  private val PdsEndpoint = uri"https://decoding.us.com/api/v1/pds/ingest"

  /**
   * Legacy Upload: Transmits the summary data to the user's PDS data vault via the DecodingUs ingestion API.
   *
   * @param summary The CoverageSummary to upload.
   * @param ec      The execution context for the future.
   * @return A Future that completes when the upload is finished.
   */
  def uploadSummary(summary: CoverageSummary)(implicit ec: ExecutionContext): Future[Unit] = {
    val request = basicRequest
      .post(PdsEndpoint)
      .body(summary)
      .response(asString)

    println(s"Initiating PDS upload for user ${summary.pdsUserId} to $PdsEndpoint...")

    request.send(backend).flatMap { response =>
      response.body match {
        case Right(body) =>
          println(s"Successfully uploaded summary. Server response: $body")
          Future.successful(())
        case Left(error) =>
          println(s"Failed to upload summary. Status: ${response.code}, Error: $error")
          Future.failed(new RuntimeException(s"PDS Upload failed with status ${response.code}: $error"))
      }
    }.recoverWith { case e: Exception =>
      println(s"Exception during PDS upload: ${e.getMessage}")
      Future.failed(e)
    }
  }

  /**
   * AT Protocol Upload: Transmits the summary data to the user's PDS data vault via the AT Protocol.
   * Uses com.atproto.repo.createRecord to store the data.
   *
   * @param user    The authenticated User.
   * @param summary The CoverageSummary to upload.
   * @param ec      The execution context for the future.
   * @return        A Future that completes when the upload is finished.
   */
  def uploadSummaryAtProto(user: User, summary: CoverageSummary)(implicit ec: ExecutionContext): Future[Unit] = {
    // Check if PDS URL is valid (it should be if user is logged in)
    val pdsUrl = if (user.pdsUrl.endsWith("/")) user.pdsUrl.dropRight(1) else user.pdsUrl
    val endpoint = uri"$pdsUrl/xrpc/com.atproto.repo.createRecord"
    
    val collection = "com.decodingus.genome.summary"
    // Generate a random Record Key (rkey)
    val rkey = java.util.UUID.randomUUID().toString

    val payload = Json.obj(
      "repo" -> Json.fromString(user.did),
      "collection" -> Json.fromString(collection),
      "rkey" -> Json.fromString(rkey),
      "record" -> summary.asJson
    )

    val request = basicRequest
      .post(endpoint)
      .header("Authorization", s"Bearer ${user.token}")
      .body(payload)
      .response(asString)

    println(s"Initiating PDS upload for user ${user.did} to $endpoint...")

    request.send(backend).flatMap { response =>
      response.body match {
        case Right(body) =>
          println(s"Successfully uploaded summary to PDS. Response: $body")
          Future.successful(())
        case Left(error) =>
          println(s"Failed to upload summary. Status: ${response.code}, Error: $error")
          Future.failed(new RuntimeException(s"PDS Upload failed with status ${response.code}: $error"))
      }
    }.recoverWith { case e: Exception =>
      println(s"Exception during PDS upload: ${e.getMessage}")
      Future.failed(e)
    }
  }
}