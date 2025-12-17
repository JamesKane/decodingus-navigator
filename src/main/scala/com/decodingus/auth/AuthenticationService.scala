package com.decodingus.auth

import com.decodingus.util.Logger
import io.circe.generic.auto.*
import sttp.client3.*
import sttp.client3.circe.*

import scala.concurrent.{ExecutionContext, Future}

case class User(
                 id: String,
                 username: String,
                 token: String,
                 did: String = "",
                 pdsUrl: String = ""
               )

case class AtpSessionResponse(did: String, handle: String, accessJwt: String, refreshJwt: String)

object AuthenticationService {
  private val log = Logger("AuthenticationService")
  private val backend = HttpClientFutureBackend()

  /**
   * Legacy Mock Login
   */
  def login(username: String, password: String): Option[User] = {
    if (username.nonEmpty && password.nonEmpty) {
      Some(User("60820188481374", username, "mock-session-token"))
    } else {
      None
    }
  }

  /**
   * AT Protocol Login
   */
  def loginAtProto(handle: String, password: String, pdsUrl: String)(implicit ec: ExecutionContext): Future[Option[User]] = {
    val cleanPdsUrl = if (pdsUrl.endsWith("/")) pdsUrl.dropRight(1) else pdsUrl
    val sessionUri = uri"$cleanPdsUrl/xrpc/com.atproto.server.createSession"

    val request = basicRequest
      .post(sessionUri)
      .body(Map("identifier" -> handle, "password" -> password))
      .response(asJson[AtpSessionResponse])

    log.info(s"Attempting login to $cleanPdsUrl for $handle...")

    request.send(backend).map { response =>
      response.body match {
        case Right(session) =>
          log.info(s"Login successful for ${session.handle} (DID: ${session.did})")
          Some(User(
            id = session.did,
            username = session.handle,
            token = session.accessJwt,
            did = session.did,
            pdsUrl = cleanPdsUrl
          ))
        case Left(error) =>
          log.warn(s"Login failed. Status: ${response.code}, Error: $error")
          None
      }
    }.recover {
      case e: Exception =>
        log.error(s"Login exception: ${e.getMessage}", e)
        None
    }
  }
}
