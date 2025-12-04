package com.decodingus.auth

case class User(id: String, username: String, token: String)

object AuthenticationService {
  /**
   * Authenticates a user.
   * For this MVP, we accept any non-empty username/password and return a fixed mock user ID.
   * In a real application, this would call an identity provider.
   */
  def login(username: String, password: String): Option[User] = {
    if (username.nonEmpty && password.nonEmpty) {
      // Simulating a successful login with a fixed PDS ID for consistency with previous placeholders
      Some(User("60820188481374", username, "mock-session-token"))
    } else {
      None
    }
  }
}
