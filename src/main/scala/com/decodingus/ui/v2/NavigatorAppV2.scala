package com.decodingus.ui.v2

import com.decodingus.i18n.I18n
import com.decodingus.service.DatabaseInitializer
import com.decodingus.ui.components.{SettingsDialog, StatusBar}
import com.decodingus.ui.theme.Theme
import com.decodingus.util.Logger
import com.decodingus.workspace.WorkbenchViewModel
import scalafx.application.JFXApp3
import scalafx.application.JFXApp3.PrimaryStage
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.Scene
import scalafx.scene.control.{Button, Label, Tooltip}
import scalafx.scene.layout.{BorderPane, HBox}

/**
 * Main entry point for DUNavigator application.
 *
 * Run with: sbt run
 */
object NavigatorAppV2 extends JFXApp3 {

  private val log = Logger("DUNavigator")

  override def start(): Unit = {
    // Initialize i18n from user preferences
    I18n.initializeFromPreferences()

    // Initialize theme from user preferences
    Theme.initializeFromPreferences()

    // Initialize database
    log.info("Initializing H2 database...")
    val databaseContext = DatabaseInitializer.initialize() match {
      case Right(context) =>
        log.info(s"Database initialized. Schema version: ${context.schemaVersion}")
        context
      case Left(error) =>
        throw new RuntimeException(s"Failed to initialize database: $error")
    }

    // Create ViewModel
    val viewModel = new WorkbenchViewModel(databaseContext)

    // Create the new MainShell
    val mainShell = new MainShell(viewModel)

    // Create status bar bound to sync notifier
    val statusBar = StatusBar(viewModel.syncNotifier)

    // Top bar with theme toggle
    val topBar = createTopBar()

    // Layout
    val mainLayout = new BorderPane {
      top = topBar
      center = mainShell
      bottom = statusBar
    }

    // Create stage with theme-aware stylesheet
    stage = new PrimaryStage {
      title = I18n.t("app.title")
      scene = new Scene(1400, 900) {
        root = mainLayout
        stylesheets.add(getClass.getResource("/style.css").toExternalForm)
        // Add theme-specific stylesheet
        val themeStylesheet = if (Theme.isDark) "/theme-dark.css" else "/theme-light.css"
        stylesheets.add(getClass.getResource(themeStylesheet).toExternalForm)
      }
    }

    log.info(s"Application started with ${if (Theme.isDark) "dark" else "light"} theme")
  }

  private def createTopBar(): BorderPane = {
    val colors = Theme.current

    val themeToggleButton = new Button {
      text = if (Theme.isDark) "☀" else "☾"
      style = s"-fx-background-color: transparent; -fx-text-fill: ${colors.textPrimary}; -fx-font-size: 16px; -fx-cursor: hand;"
      onAction = _ => {
        Theme.toggle()
        // Note: Full theme switch requires app restart for now
        // Future: implement live theme switching
        log.info(s"Theme changed to ${if (Theme.isDark) "dark" else "light"}. Restart app to apply.")
      }
    }

    new BorderPane {
      style = s"-fx-background-color: ${colors.surface}; -fx-padding: 10;"

      left = new HBox {
        alignment = Pos.CenterLeft
        children = Seq(
          new Label(I18n.t("app.title")) {
            style = s"-fx-font-size: 16px; -fx-font-weight: bold; -fx-text-fill: ${colors.textPrimary};"
          }
        )
      }

      right = new HBox(10) {
        alignment = Pos.CenterRight
        padding = Insets(0, 10, 0, 0)
        children = Seq(
          settingsButton,
          themeToggleButton
        )
      }
    }
  }

  private def settingsButton: Button = {
    val colors = Theme.current
    new Button("⚙") {
      style = s"-fx-background-color: transparent; -fx-text-fill: ${colors.textPrimary}; -fx-font-size: 16px; -fx-cursor: hand;"
      tooltip = new Tooltip(I18n.t("settings.title"))
      onAction = _ => openSettings()
    }
  }

  private def openSettings(): Unit = {
    val dialog = new SettingsDialog()
    Option(stage.scene.value).flatMap(s => Option(s.getWindow)).foreach { window =>
      dialog.initOwner(window)
    }
    dialog.showAndWait()
  }
}
