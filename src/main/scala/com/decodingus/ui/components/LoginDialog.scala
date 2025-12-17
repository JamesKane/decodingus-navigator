package com.decodingus.ui.components

import com.decodingus.auth.{AuthenticationService, User}
import com.decodingus.config.FeatureToggles
import scalafx.Includes.*
import scalafx.application.Platform
import scalafx.geometry.Insets
import scalafx.scene.control.*
import scalafx.scene.layout.GridPane
import scalafx.stage.Window

import java.util.concurrent.Executors
import scala.concurrent.duration.*
import scala.concurrent.{Await, ExecutionContext}

object LoginDialog {
  def show(ownerWindow: Window): Option[User] = {
    val dialog = new Dialog[User]() {
      initOwner(ownerWindow)
      title = "Login to Decoding-Us"
      headerText = if (FeatureToggles.atProtocolEnabled) "Sign in to access PDS features" else "Login"
    }

    val loginButtonType = new ButtonType("Login", ButtonBar.ButtonData.OKDone)
    dialog.dialogPane().buttonTypes = Seq(loginButtonType, ButtonType.Cancel)

    val usernameField = new TextField() {
      promptText = if (FeatureToggles.atProtocolEnabled) "Handle" else "Username"
    }
    val passwordField = new PasswordField() {
      promptText = "Password"
    }

    // PDS URL Field (only if feature enabled)
    val pdsUrlField = new TextField() {
      promptText = "PDS URL (e.g. https://bsky.social)"
      text = "https://bsky.social"
    }

    val grid = new GridPane() {
      hgap = 10
      vgap = 10
      padding = Insets(20, 150, 10, 10)
      add(new Label(if (FeatureToggles.atProtocolEnabled) "Handle:" else "Username:"), 0, 0)
      add(usernameField, 1, 0)
      add(new Label("Password:"), 0, 1)
      add(passwordField, 1, 1)

      if (FeatureToggles.atProtocolEnabled) {
        add(new Label("PDS URL:"), 0, 2)
        add(pdsUrlField, 1, 2)
      }
    }

    dialog.dialogPane().content = grid

    Platform.runLater(usernameField.requestFocus())

    dialog.resultConverter = dialogButton => {
      if (dialogButton == loginButtonType) {
        if (FeatureToggles.atProtocolEnabled) {
          implicit val ec = ExecutionContext.fromExecutor(Executors.newSingleThreadExecutor())
          try {
            val future = AuthenticationService.loginAtProto(
              usernameField.text.value,
              passwordField.text.value,
              pdsUrlField.text.value
            )
            Await.result(future, 15.seconds).orNull
          } catch {
            case e: Exception =>
              e.printStackTrace()
              null
          }
        } else {
          // Legacy sync login
          AuthenticationService.login(usernameField.text.value, passwordField.text.value).orNull
        }
      } else {
        null
      }
    }

    val result = dialog.showAndWait()

    result match {
      case Some(u) =>
        val user = u.asInstanceOf[com.decodingus.auth.User]
        if (user != null) Some(user) else None
      case _ => None
    }
  }
}
