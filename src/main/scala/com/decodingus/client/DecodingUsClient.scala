package com.decodingus.client

import com.decodingus.model.LibraryStats
import sttp.client3._
import sttp.client3.circe._
import io.circe.generic.auto._

import scala.concurrent.{ExecutionContext, Future}

case class PdsRegistrationRequest(did: String, token: String, pdsUrl: String)

object DecodingUsClient {

  private val backend = HttpClientFutureBackend()
  private val BaseUrl = uri"https://decoding.us.com/api/v1"

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
}
