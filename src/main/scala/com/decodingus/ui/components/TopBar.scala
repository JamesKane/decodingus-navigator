package com.decodingus.ui.components

import com.decodingus.auth.User
import com.decodingus.config.FeatureToggles
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.{Button, Label}
import scalafx.scene.layout.HBox

class TopBar(onLogin: () => Unit, onLogout: () => Unit) extends HBox {
  alignment = Pos.CenterRight
  padding = Insets(10)
  style = "-fx-background-color: #333333;"

  def update(currentUser: Option[User]): Unit = {
    children.clear()
    currentUser match {
      case Some(user) =>
        val userLabel = new Label(s"Logged in as: ${user.username}") {
          styleClass.add("info-label")
          style = "-fx-text-fill: #E0E0E0; -fx-padding: 0 10 0 0;"
        }
        val logoutBtn = new Button("Logout") {
          styleClass.add("button-select")
          onAction = _ => onLogout()
        }
        children.addAll(userLabel, logoutBtn)
      case None =>
        if (FeatureToggles.authEnabled) {
          val loginBtn = new Button("Login") {
            styleClass.add("button-select")
            onAction = _ => onLogin()
          }
          children.add(loginBtn)
        }
    }
  }
}
