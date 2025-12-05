package com.decodingus.auth

import sttp.client3._
import sttp.client3.circe._
import io.circe.generic.auto._
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

    println(s"Attempting login to $cleanPdsUrl for $handle...")

    request.send(backend).map { response =>
      response.body match {
        case Right(session) =>
          println(s"Login successful for ${session.handle} (DID: ${session.did})")
          Some(User(
            id = session.did,
            username = session.handle,
            token = session.accessJwt,
            did = session.did,
            pdsUrl = cleanPdsUrl
          ))
        case Left(error) =>
          println(s"Login failed. Status: ${response.code}, Error: $error")
          None
      }
    }.recover {
      case e: Exception =>
        println(s"Login exception: ${e.getMessage}")
        None
    }
  }
}
