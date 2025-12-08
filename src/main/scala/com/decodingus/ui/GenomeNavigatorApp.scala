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
  // Initialize Log4j workaround for GATK compatibility with Java 17+
  // This MUST be done before any GATK classes are loaded to avoid
  // "No class provided" errors from Log4j's stack-based caller lookup
  System.setProperty("log4j2.loggerContextFactory", "org.apache.logging.log4j.simple.SimpleLoggerContextFactory")
  System.setProperty("log4j2.disable.jmx", "true")

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

  // --- Commenting out old analysis methods for refactoring into WorkbenchView ---
  /*
  private def selectFile(): Unit = {
    println("File selection dialog would open here.")
    startInitialAnalysis("path/to/mock/file.bam") // Renamed
  }

  private def startInitialAnalysis(filePath: String): Unit = {
    currentFilePath = filePath
    val progressLabel = new Label("Initializing analysis...") {
      styleClass.add("progress-label")
    }
    val shaLabel = new Label("Calculating file signature...") {
      styleClass.add("info-label")
    }
    val progressBar = new ProgressBar {
      prefWidth = 400
    }
    val progressIndicator = new ProgressIndicator

    val progressScreen = new VBox(20) {
      alignment = Pos.Center
      styleClass.add("root-pane")
      children = Seq(
        progressLabel,
        shaLabel,
        new HBox(20) {
          alignment = Pos.Center
          children = Seq(progressBar, progressIndicator)
        }
      )
    }

    mainLayout.children = progressScreen

    // Background Task for SHA-256
    val shaTask = new jfxc.Task[String]() {
      override def call(): String = {
        AnalysisCache.calculateSha256(new File(filePath))
      }
    }

    shaTask.setOnSucceeded(_ => {
      currentSha256 = Some(shaTask.getValue)
      Platform.runLater {
        shaLabel.text = s"File Signature: ${shaTask.getValue.take(8)}..."
        // If we are still on the progress screen, this will just update the label
        // The check for cache existence happens in showInitialResultsAndChoices
        // However, if the main analysis finishes first, we need to update the choices screen dynamically.
        // We'll handle this by checking currentSha256 in showInitialResultsAndChoices
      }
    })

    new Thread(shaTask).start()


    val jfxTask = new jfxc.Task[(LibraryStats, String)]() { // Now returns LibraryStats and referencePath
      override def call(): (LibraryStats, String) = {
        try {
          // Step 1: Detect Reference Build
          Platform.runLater { progressLabel.text = "Detecting reference build from BAM/CRAM header..." }
          val header = SamReaderFactory.makeDefault().open(new File(filePath)).getFileHeader
          val libraryStatsProcessor = new LibraryStatsProcessor()
          val referenceBuild = libraryStatsProcessor.detectReferenceBuild(header)
          if (referenceBuild == "Unknown") {
            throw new IllegalStateException("Could not determine reference build from BAM/CRAM header.")
          }

          // Step 2: Resolve Reference Path
          Platform.runLater { progressLabel.text = s"Resolving reference: $referenceBuild" }
          val referenceGateway = new ReferenceGateway((done, total) => {
            Platform.runLater {
              progressLabel.text = s"Downloading reference: $done / $total bytes"
              updateProgress(done, total)
            }
          })
          val referencePath = referenceGateway.resolve(referenceBuild) match {
            case Right(path) => path.toString
            case Left(error) => throw new Exception(error)
          }

          // Phase 1: Library Stats (quick scan)
          val libraryStats = libraryStatsProcessor.process(filePath, referencePath, (message, current, total) => {
            Platform.runLater { progressLabel.text = s"Library Stats: $message" }
            updateProgress(current, total)
          })
          
          (libraryStats, referencePath) // Return only library stats and reference path
        } catch {
          case e: Exception =>
            e.printStackTrace()
            cancel()
            throw e
        }
      }
    }

    progressBar.progress <== jfxTask.progressProperty
    progressIndicator.progress <== jfxTask.progressProperty

    jfxTask.setOnSucceeded(_ => {
      val (libraryStats, referencePath) = jfxTask.getValue
      currentLibraryStats = Some(libraryStats)
      currentReferencePath = Some(referencePath)
      // Wait for SHA if it's not done? Or just proceed and update later?
      // Let's proceed, and in showInitialResultsAndChoices we can add a listener or check
      showInitialResultsAndChoices(libraryStats, referencePath)
    })

    jfxTask.setOnFailed(_ => {
      val errorMessage = Option(jfxTask.getException).map(_.getMessage).getOrElse("Unknown error during initial analysis.")
      val errorScreen = new VBox(20) {
        alignment = Pos.Center
        styleClass.add("root-pane")
        children = Seq(
          new Label("Initial Analysis Failed!") { // Changed title for clarity
            styleClass.add("error-label")
          },
          new Label(s"Reason: $errorMessage") { // Display the specific error message
            styleClass.add("info-label")
            wrapText = true
            textAlignment = TextAlignment.Center
          },
          new Button("Back to Welcome") {
            onAction = _ => mainLayout.children = new WelcomeScreen(
              onFileSelected = filePath => startInitialAnalysis(filePath),
              onSelectFileClicked = () => selectFile()
            )
          }
        )
      }
      mainLayout.children = errorScreen
    })

    new Thread(jfxTask).start()
  }

  // New method to show initial results and choices
  private def showInitialResultsAndChoices(libraryStats: LibraryStats, referencePath: String): Unit = {
    val resultsTitle = new Label("Initial Analysis Results") {
      styleClass.add("title-label")
    }

    val resultsVBox = new VBox(5) { // Use VBox to stack lines, with a small vertical gap
      padding = Insets(20)
      styleClass.add("stats-container") // Add a style class for potential styling
      children = Seq(
        new HBox(5) { // HBox for the first line
          children = Seq(
            new Text("Sample Name:") { styleClass.add("stat-name") },
            new Text(libraryStats.sampleName) { styleClass.add("stat-value") },
            new Text(s"(${f"${libraryStats.readCount}%,d"} reads sampled)") { styleClass.add("stat-value") }
          )
        },
        new HBox(20) { // HBox for the second line, with more horizontal gap
          children = Seq(
            new Text("Platform:") { styleClass.add("stat-name") },
            new Text(libraryStats.inferredPlatform) { styleClass.add("stat-value") },
            new Text("Instrument Type:") { styleClass.add("stat-name") },
            new Text(libraryStats.mostFrequentInstrument) { styleClass.add("stat-value") }
          )
        },
        new HBox(20) { // HBox for the third line, with more horizontal gap
          children = Seq(
            new Text("Aligner:") { styleClass.add("stat-name") },
            new Text(libraryStats.aligner) { styleClass.add("stat-value") },
            new Text("Reference:") { styleClass.add("stat-name") },
            new Text(libraryStats.referenceBuild) { styleClass.add("stat-value") }
          )
        }
      )
    }

    val deepAnalysisHelpText = new Text("Performs a comprehensive WGS (Whole Genome Sequencing) analysis, including detailed coverage metrics and callable loci, generating plots and tables for visual inspection. This can take significant time for large genomes.") {
      wrappingWidth = 400
      styleClass.add("info-label")
    }

    val deepAnalysisButton = new Button("Perform Deep Coverage Analysis") {
      onAction = _ => startDeepCoverageAnalysis(currentFilePath, libraryStats, referencePath) // Calls the renamed method
    }

    val loadCacheButton = new Button("Load Cached Analysis") {
      styleClass.add("button-success")
      visible = false
      managed = false // Don't take up space if not visible
      onAction = _ => {
        currentSha256.flatMap(AnalysisCache.load).foreach {
          coverageSummary = Some(cachedSummary)
          // We don't have the SVG strings stored in the summary object, which is a limitation of the current cache.
          // For now, we will pass an empty list for SVGs or handle it gracefully.
          // Ideally, we should cache SVGs too or regenerate them (fast) or store them in the summary.
          showResults(cachedSummary, List.empty) 
        }
      }
    }

    // Check cache availability
    def updateCacheButton(): Unit = {
      currentSha256.foreach {
        if (AnalysisCache.exists(sha)) {
          loadCacheButton.visible = true
          loadCacheButton.managed = true
          loadCacheButton.text = s"Load Cached Analysis (${sha.take(8)}...)"
        }
      }
    }
    
    // Initial check
    updateCacheButton()
    
    // Poll for SHA completion if not yet done (simple hack since we don't have a bound property easily here)
    if (currentSha256.isEmpty) {
      val timer = new java.util.Timer()
      timer.schedule(new java.util.TimerTask {
        override def run(): Unit = {
          Platform.runLater {
            if (currentSha256.isDefined) {
              updateCacheButton()
              timer.cancel()
            }
          }
        }
      }, 1000, 1000)
    }


    val haplogroupProviderChoice = new ChoiceBox[TreeProviderType] {
      items = ObservableBuffer(TreeProviderType.FTDNA, TreeProviderType.DECODINGUS)
      value = TreeProviderType.DECODINGUS // Default to DecodingUs
    }

    val haplogroupTypeChoice = new ChoiceBox[TreeType] {
      items = ObservableBuffer(TreeType.YDNA, TreeType.MTDNA)
      value = TreeType.YDNA // Default to YDNA
    }

    val analyzeHaplogroupBtn = new Button("Analyze Haplogroup") {
      onAction = _ => startHaplogroupAnalysis(haplogroupTypeChoice.value(), haplogroupProviderChoice.value())
    }

    val privateSnpButton = new Button("Find Private SNPs") {
      id = "privateSnpButton"
      onAction = _ => startPrivateSnpAnalysis(analyzedHaplogroupType.getOrElse(TreeType.YDNA)) // Use the type of the analyzed haplogroup
      disable = true // Disabled until haplogroup is determined and is YDNA
    }

    // Logic to disable DecodingUs mtDNA
    haplogroupProviderChoice.value.onChange { (_, _, newProvider) =>
      if (newProvider == TreeProviderType.DECODINGUS) {
        haplogroupTypeChoice.items.value.clear()
        haplogroupTypeChoice.items.value.add(TreeType.YDNA)
        haplogroupTypeChoice.value = TreeType.YDNA // Ensure YDNA is selected
      } else {
        haplogroupTypeChoice.items.value.clear()
        haplogroupTypeChoice.items.value.add(TreeType.YDNA)
        haplogroupTypeChoice.items.value.add(TreeType.MTDNA)
        haplogroupTypeChoice.value = TreeType.YDNA // Default to YDNA
      }
    }

    // Logic to enable/disable private SNP button - this will be called from startHaplogroupAnalysis
    // when a haplogroup analysis completes.
    // Initial state is disable, and is updated on a successful haplogroup analysis.

    val haplogroupControls = new HBox(10) {
      alignment = Pos.CenterLeft
      padding = Insets(10, 0, 10, 0)
      children = Seq(
        new Label("Provider:"),
        haplogroupProviderChoice,
        new Label("Type:"),
        haplogroupTypeChoice,
        analyzeHaplogroupBtn,
        privateSnpButton
      )
    }

    val resultsVBox = new VBox(20) {
      alignment = Pos.Center
      styleClass.add("root-pane")
      padding = Insets(20)
      children = Seq(resultsTitle, statsGrid, new Separator(), contigBreakdownTitle, contigTable, webView, new Separator(), haplogroupControls)
    }

    if (FeatureToggles.pdsSubmissionEnabled) {
      val pdsBox = new VBox(15) {
        alignment = Pos.Center
        padding = Insets(20)
        children = Seq(
          new Label("Help advance research by securely contributing your anonymized summary data to your Personal Data Store.") {
            wrapText = true
            textAlignment = TextAlignment.Center
            styleClass.add("info-label")
          },
          new CheckBox("I agree to upload my anonymized summary data.") {
            selected = true
            style = "-fx-text-fill: #E0E0E0; -fx-font-size: 14px;"
          },
          new Button("Upload to PDS") {
            styleClass.add("button-upload")
            onAction = _ => {
              // If AT Protocol enabled, ensure login flow handles registration
              val userToUse = if (FeatureToggles.atProtocolEnabled) {
                 currentUser.orElse {
                    val loggedIn = LoginDialog.show(stage)
                    loggedIn.foreach { u =>
                      currentUser = Some(u)
                      topBar.update(currentUser)
                      DecodingUsClient.registerPds(u.did, u.token, u.pdsUrl)
                    }
                    loggedIn
                 }
              } else {
                // Legacy: Login is optional for upload context locally but PDS ID needed? 
                // Actually legacy upload used summary.pdsUserId. 
                // We'll assume current logic is sufficient or prompts login.
                // But let's stick to the pattern: ensure user is logged in.
                currentUser.orElse {
                    val loggedIn = LoginDialog.show(stage)
                    loggedIn.foreach { u =>
                      currentUser = Some(u)
                      topBar.update(currentUser)
                    }
                    loggedIn
                 }
              }

              userToUse match {
                case Some(user) =>
                  text = "Uploading..."
                  disable = true
                  
                  val uploadFuture = if (FeatureToggles.atProtocolEnabled) {
                    PdsClient.uploadSummaryAtProto(user, summary)
                  } else {
                    PdsClient.uploadSummary(summary)
                  }

                  uploadFuture.onComplete {
                    case scala.util.Success(_) =>
                      Platform.runLater {
                        text = "Upload Complete!"
                        styleClass.remove("button-upload")
                        styleClass.add("button-success")
                      }
                    case scala.util.Failure(ex) =>
                      Platform.runLater {
                        text = "Upload to PDS" // Reset text
                        disable = false
                        new Alert(AlertType.Error) {
                          initOwner(stage)
                          title = "Upload Failed"
                          headerText = "Could not upload to PDS"
                          contentText = ex.getMessage
                        }.showAndWait()
                      }
                  }
                case None =>
                // User cancelled login
              }
            }
          }
        )
      }
      resultsVBox.children.addAll(new Separator(), pdsBox)
    }

    // Developer features section
    if (FeatureToggles.developerFeatures.saveJsonEnabled) {
      val saveJsonButton = new Button("Save Analysis JSON") {
        styleClass.add("button-select")
        onAction = _ => {
          import io.circe.syntax._
          import io.circe.generic.semiauto._ // Needed for implicit Encoder. 
          import java.time.LocalDateTime
          import java.time.format.DateTimeFormatter
          import javafx.stage.FileChooser
          import java.nio.file.Files
          import java.nio.file.Paths

          // Ensure implicit encoders are available for Summary components if not already globally present
          val jsonString = summary.asJson.spaces2

          val now = LocalDateTime.now()
          val formatter = DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss")
          val defaultFileName = s"coverage_summary_${now.format(formatter)}.json"

          val fileChooser = new FileChooser()
          fileChooser.title = "Save Analysis Summary JSON"
          fileChooser.initialFileName = defaultFileName
          fileChooser.extensionFilters.add(new FileChooser.ExtensionFilter("JSON files (*.json)", "*.json"))

          Option(fileChooser.showSaveDialog(stage)).foreach {
            try {
              Files.write(file.toPath, jsonString.getBytes)
              new Alert(AlertType.Information) {
                initOwner(stage)
                title = "JSON Saved"
                headerText = "Analysis summary saved successfully."
                contentText = s"File: ${file.getAbsolutePath}"
              }.showAndWait()
            } catch {
              case e: Exception =>
                new Alert(AlertType.Error) {
                  initOwner(stage)
                  title = "Save Failed"
                  headerText = "Could not save analysis summary JSON."
                  contentText = e.getMessage
                }.showAndWait()
            }
          }
        }
      }
      resultsVBox.children.addAll(new Separator(), saveJsonButton)
    }

    val resultsScreen = new ScrollPane {
      content = resultsVBox
      fitToWidth = true
    }

    mainLayout.children = resultsScreen
  }

  private def startHaplogroupAnalysis(treeType: TreeType, providerType: TreeProviderType): Unit = {
    // Ensure currentLibraryStats and currentReferencePath are available
    for {
      summary <- currentLibraryStats // Use currentLibraryStats
      tree <- haplogroupTree
      bestHg <- bestHaplogroup
      tpInstance <- treeProviderInstance // Get the stored TreeProvider instance
      // Ensure private SNP analysis is only run for YDNA
      if treeType == TreeType.YDNA
    } {
      val progressDialog = new Dialog[Unit]() {
        initOwner(stage)
        title = "Haplogroup Analysis"
        headerText = "Running haplogroup analysis..."
        dialogPane().content = new ProgressIndicator()
      }

      val haplogroupTask = new jfxc.Task[Either[String, List[com.decodingus.haplogroup.model.HaplogroupResult]]]() {
        override def call(): Either[String, List[com.decodingus.haplogroup.model.HaplogroupResult]] = {
          val treeProvider: TreeProvider = providerType match {
            case TreeProviderType.FTDNA => new FtdnaTreeProvider(treeType)
            case TreeProviderType.DECODINGUS => new DecodingUsTreeProvider(treeType)
          }
          treeProviderInstance = Some(treeProvider) // Store the instance
          
          haplogroupTree = treeProvider.loadTree(summary.referenceBuild).toOption
          
          val processor = new HaplogroupProcessor()
          processor.analyze(currentFilePath, summary, treeType, providerType, (message, current, total) => {
            Platform.runLater {
              progressDialog.headerText = message
            }
          })
        }
      }

      haplogroupTask.setOnSucceeded(_ => {
        progressDialog.close()
        haplogroupTask.getValue match {
          case Right(results) =>
            bestHaplogroup = results.headOption
            analyzedHaplogroupType = Some(treeType) // Store the analyzed tree type
            val topResult = bestHaplogroup.map(_.name).getOrElse("Not found")
            new Alert(AlertType.Information) {
              initOwner(stage)
              title = "Haplogroup Analysis Complete"
              headerText = "Top Haplogroup Result:"
              contentText = topResult
            }.showAndWait()
            // Enable private SNP button only if YDNA was analyzed
            mainLayout.scene().lookup("#privateSnpButton").setDisable(analyzedHaplogroupType != Some(TreeType.YDNA))
          case Left(error) =>
            new Alert(AlertType.Error) {
              initOwner(stage)
              title = "Haplogroup Analysis Failed"
              headerText = "An error occurred during haplogroup analysis."
              contentText = error
            }.showAndWait()
        }
      })

      haplogroupTask.setOnFailed(_ => {
        progressDialog.close()
        new Alert(AlertType.Error) {
          initOwner(stage)
          title = "Haplogroup Analysis Failed"
          headerText = "A critical error occurred during haplogroup analysis."
          contentText = haplogroupTask.getException.getMessage
        }.showAndWait()
      })

      progressDialog.show()
      new Thread(haplogroupTask).start()
    }
  }

  private def startPrivateSnpAnalysis(treeType: TreeType): Unit = {
    for {
      libraryStats <- currentLibraryStats // Use currentLibraryStats
      tree <- haplogroupTree
      bestHg <- bestHaplogroup
      tpInstance <- treeProviderInstance // Get the stored TreeProvider instance
      // Ensure private SNP analysis is only run for YDNA
      if treeType == TreeType.YDNA
    } {
      val progressDialog = new Dialog[Unit]() {
        initOwner(stage)
        title = "Private SNP Analysis"
        headerText = "Finding private SNPs..."
        dialogPane().content = new ProgressIndicator()
      }

      val privateSnpTask = new jfxc.Task[File]() {
        override def call(): File = {
          val referenceGateway = new ReferenceGateway((_, _) => {})
          val referencePath = referenceGateway.resolve(libraryStats.referenceBuild).toOption.get.toString

          val nameToHaplogroup = tree.flatMap(root => {
            def flatten(h: Haplogroup): List[(String, Haplogroup)] = (h.name -> h) :: h.children.flatMap(flatten)
            flatten(root)
          }).toMap

          val pathHgs = mutable.ListBuffer[Haplogroup]()
          var currentNameOpt: Option[String] = Some(bestHg.name)
          while (currentNameOpt.isDefined) {
            val currentName = currentNameOpt.get
            nameToHaplogroup.get(currentName) match {
              case Some(hg) =>
                pathHgs += hg
                currentNameOpt = hg.parent
              case None =>
                currentNameOpt = None
            }
          }

          val knownLoci = pathHgs.flatMap(_.loci).toSet

          val processor = new PrivateSnpProcessor()
          val contig = if (treeType == TreeType.YDNA) "chrY" else "chrM"
          val privateSnps = processor.findPrivateSnps(
            currentFilePath,
            referencePath,
            contig,
            knownLoci,
            tpInstance.sourceBuild, // Pass treeSourceBuild
            libraryStats.referenceBuild, // Pass referenceBuild
            (message, _, _) => {
              Platform.runLater {
                progressDialog.headerText = message
              }
            }
          )

          val reportFile = new File(s"private_snps_${libraryStats.sampleName}.txt")
          processor.writeReport(privateSnps, reportFile)
          reportFile
        }
      }

      privateSnpTask.setOnSucceeded(_ => {
        progressDialog.close()
        val reportFile = privateSnpTask.getValue
        new Alert(AlertType.Information) {
          initOwner(stage)
          title = "Private SNP Analysis Complete"
          headerText = "Private SNP report generated."
          contentText = s"Report saved to: ${reportFile.getAbsolutePath}"
        }.showAndWait()
      })

      privateSnpTask.setOnFailed(_ => {
        progressDialog.close()
        new Alert(AlertType.Error) {
          initOwner(stage)
          title = "Private SNP Analysis Failed"
          headerText = "A critical error occurred during private SNP analysis."
          contentText = privateSnpTask.getException.getMessage
        }.showAndWait()
      })

      progressDialog.show()
      new Thread(privateSnpTask).start()
    }
  }
  */ // Explicitly closing the comment block
}
