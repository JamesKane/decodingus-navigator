package com.decodingus.ui

import com.decodingus.analysis._
import com.decodingus.auth._
import com.decodingus.config.FeatureToggles
import com.decodingus.haplogroup.model.Haplogroup
import com.decodingus.haplogroup.tree.{TreeProvider, TreeProviderType, TreeType}
import com.decodingus.haplogroup.vendor.{DecodingUsTreeProvider, FtdnaTreeProvider}
import com.decodingus.model._
import com.decodingus.pds.PdsClient
import com.decodingus.client.DecodingUsClient
import com.decodingus.refgenome.ReferenceGateway
import com.decodingus.ui.components._
import htsjdk.samtools.SamReaderFactory
import io.circe.syntax._
import com.decodingus.analysis.AnalysisCache.{coverageSummaryEncoder, libraryStatsEncoder, wgsMetricsEncoder, contigSummaryEncoder} // Import implicits
import javafx.concurrent as jfxc
import scalafx.Includes.*
import scalafx.application.JFXApp3.PrimaryStage
import scalafx.application.{JFXApp3, Platform}
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.Scene
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.control.*
import scalafx.scene.control.TableColumn.sfxTableColumn2jfx
import scalafx.scene.input.{DragEvent, TransferMode}
import scalafx.scene.layout.{BorderPane, GridPane, HBox, Priority, Region, StackPane, VBox}
import scalafx.scene.text.{Text, TextAlignment}
import scalafx.scene.web.WebView

import com.decodingus.workspace.model.{Workspace, Biosample, Project, SequenceRun, Alignment, AlignmentMetrics, ContigMetrics, HaplogroupResult, HaplogroupAssignments, FileInfo, RecordMeta} // Explicitly import all workspace models
import com.decodingus.workspace.{WorkspaceService, LiveWorkspaceService, WorkbenchViewModel} // Import the workspace service trait and its live implementation, and the ViewModel
import com.decodingus.ui.components.WorkbenchView // Import the new WorkbenchView

import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import javafx.stage.FileChooser
import scala.collection.mutable
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.Await
import scala.concurrent.duration._


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

  // ViewModel is created early so topBar can reference it
  private lazy val viewModel = new WorkbenchViewModel(LiveWorkspaceService)

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
    val borderPaneRoot = new BorderPane {
      top = topBar
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
