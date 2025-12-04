package com.decodingus.ui.components

import com.decodingus.auth.{AuthenticationService, User}
import scalafx.Includes._
import scalafx.application.Platform
import scalafx.geometry.Insets
import scalafx.scene.control._
import scalafx.scene.layout.GridPane
import scalafx.stage.Window

object LoginDialog {
  def show(ownerWindow: Window): Option[User] = {
    val dialog = new Dialog[User]() {
      initOwner(ownerWindow)
      title = "Login to Decoding-Us"
      headerText = "Sign in to access PDS features"
    }

    val loginButtonType = new ButtonType("Login", ButtonBar.ButtonData.OKDone)
    dialog.dialogPane().buttonTypes = Seq(loginButtonType, ButtonType.Cancel)

    val usernameField = new TextField() { promptText = "Username" }
    val passwordField = new PasswordField() { promptText = "Password" }

    val grid = new GridPane() {
      hgap = 10
      vgap = 10
      padding = Insets(20, 150, 10, 10)
      add(new Label("Username:"), 0, 0)
      add(usernameField, 1, 0)
      add(new Label("Password:"), 0, 1)
      add(passwordField, 1, 1)
    }

    dialog.dialogPane().content = grid

    Platform.runLater(usernameField.requestFocus())

    dialog.resultConverter = dialogButton => {
      if (dialogButton == loginButtonType) {
        AuthenticationService.login(usernameField.text.value, passwordField.text.value).orNull
      } else {
        null
      }
    }

    val result = dialog.showAndWait()

    result match {
      case Some(u) =>
        // Explicit cast to bypass ScalaFX DConvert type issue
        val user = u.asInstanceOf[com.decodingus.auth.User]
        if (user != null) Some(user) else None
      case _ => None
    }
  }
}
