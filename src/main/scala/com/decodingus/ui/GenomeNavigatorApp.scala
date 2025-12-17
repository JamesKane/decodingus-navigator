package com.decodingus.ui

import com.decodingus.analysis.*
import com.decodingus.analysis.AnalysisCache.{contigSummaryEncoder, coverageSummaryEncoder, libraryStatsEncoder, wgsMetricsEncoder}
import com.decodingus.auth.*
import com.decodingus.client.DecodingUsClient
import com.decodingus.config.FeatureToggles
import com.decodingus.haplogroup.model.Haplogroup
import com.decodingus.haplogroup.tree.{TreeProvider, TreeProviderType, TreeType}
import com.decodingus.haplogroup.vendor.{DecodingUsTreeProvider, FtdnaTreeProvider}
import com.decodingus.model.*
import com.decodingus.pds.PdsClient
import com.decodingus.refgenome.ReferenceGateway
import com.decodingus.service.{DatabaseContext, DatabaseInitializer}
import com.decodingus.ui.components.*
import com.decodingus.workspace.model.*
import com.decodingus.workspace.{H2WorkspaceAdapter, WorkbenchViewModel, WorkspaceService}
import htsjdk.samtools.SamReaderFactory
import io.circe.syntax.*
import javafx.concurrent as jfxc
import javafx.stage.FileChooser
import scalafx.Includes.*
import scalafx.application.JFXApp3.PrimaryStage
import scalafx.application.{JFXApp3, Platform}
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.Scene
import scalafx.scene.control.*
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.control.TableColumn.sfxTableColumn2jfx
import scalafx.scene.input.{DragEvent, TransferMode}
import scalafx.scene.layout.*
import scalafx.scene.text.{Text, TextAlignment}
import scalafx.scene.web.WebView

import java.io.File
import java.nio.file.{Files, Paths}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import scala.collection.mutable
import scala.concurrent.Await
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration.*


case class ContigAnalysisRow(
                              contig: String,
                              callableBases: Long,
                              noCoverage: Long,
                              lowCoverage: Long,
                              poorMq: Long,
                              refN: Long,
                              svgFile: String
                            )

object GenomeNavigatorApp extends JFXApp3 {
  private val mainLayout = new StackPane()

  // Removed old state variables, now managing a single Workspace object
  // private var currentWorkspace: Workspace = Workspace(List.empty, List.empty) // ViewModel now owns this

  private var currentUser: Option[User] = None // Keep currentUser for login

  // Database context for H2 persistence (initialized lazily)
  private lazy val databaseContext: DatabaseContext = {
    println("[App] Initializing H2 database...")
    DatabaseInitializer.initialize() match {
      case Right(context) =>
        println(s"[App] H2 database initialized successfully. Schema version: ${context.schemaVersion}")
        context
      case Left(error) =>
        throw new RuntimeException(s"Failed to initialize H2 database: $error")
    }
  }

  // ViewModel is created early so topBar can reference it
  // Uses DatabaseContext directly for atomic CRUD operations
  private lazy val viewModel = new WorkbenchViewModel(databaseContext)

  private lazy val topBar: TopBar = new TopBar(
    onLogin = () => {
      LoginDialog.show(stage).foreach { user =>
        currentUser = Some(user)
        topBar.update(currentUser)

        // Update ViewModel with the logged-in user
        viewModel.currentUser.value = Some(user)

        if (FeatureToggles.atProtocolEnabled) {
          // Register PDS with the main server asynchronously if AT Protocol is enabled
          DecodingUsClient.registerPds(user.did, user.token, user.pdsUrl).failed.foreach { e =>
            Platform.runLater {
              new Alert(AlertType.Error) {
                initOwner(stage)
                title = "PDS Registration Failed"
                headerText = "Could not register PDS with DecodingUs"
                contentText = e.getMessage
              }.showAndWait()
            }
          }

          // Trigger PDS sync now that user is logged in
          viewModel.syncFromPdsIfAvailable()
        }
      }
    },
    onLogout = () => {
      currentUser = None
      topBar.update(currentUser)

      // Clear user from ViewModel
      viewModel.currentUser.value = None
    }
  )

  override def start(): Unit = {
    // Create status bar bound to ViewModel's sync notifier
    val statusBar = StatusBar(viewModel.syncNotifier)

    val borderPaneRoot = new BorderPane {
      top = topBar
      bottom = statusBar
    }

    stage = new PrimaryStage {
      title = "Decoding-Us Navigator - Workbench"
      scene = new Scene(1200, 850) {
        root = borderPaneRoot
        stylesheets.add(getClass.getResource("/style.css").toExternalForm)
      }
    }

    topBar.update(currentUser)
    // Instantiate the WorkbenchView with the ViewModel
    val workbenchView = new WorkbenchView(viewModel)
    borderPaneRoot.center = workbenchView
  }
}
