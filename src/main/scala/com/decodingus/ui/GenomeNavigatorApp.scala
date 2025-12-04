package com.decodingus.ui

import com.decodingus.analysis._
import com.decodingus.auth._
import com.decodingus.config.FeatureToggles
import com.decodingus.haplogroup.model.Haplogroup
import com.decodingus.haplogroup.tree.{TreeProvider, TreeProviderType, TreeType}
import com.decodingus.haplogroup.vendor.{DecodingUsTreeProvider, FtdnaTreeProvider}
import com.decodingus.model._
import com.decodingus.pds.PdsClient
import com.decodingus.refgenome.ReferenceGateway
import com.decodingus.ui.components._
import htsjdk.samtools.SamReaderFactory
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

import java.io.File
import scala.collection.mutable
import scala.concurrent.ExecutionContext.Implicits.global

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

  private var currentFilePath: String = ""
  private var currentLibraryStats: Option[LibraryStats] = None // Store initial library stats
  private var currentReferencePath: Option[String] = None     // Store resolved reference path
  private var coverageSummary: Option[CoverageSummary] = None // This will be set after deep analysis
  private var haplogroupTree: Option[List[Haplogroup]] = None
  private var bestHaplogroup: Option[com.decodingus.haplogroup.model.HaplogroupResult] = None
  private var treeProviderInstance: Option[TreeProvider] = None
  private var analyzedHaplogroupType: Option[TreeType] = None // To store the type of haplogroup that was analyzed

  private var currentSha256: Option[String] = None
  private var currentUser: Option[User] = None

  private lazy val topBar: TopBar = new TopBar(
    onLogin = () => {
      LoginDialog.show(stage).foreach { user =>
        currentUser = Some(user)
        topBar.update(currentUser)
      }
    },
    onLogout = () => {
      currentUser = None
      topBar.update(currentUser)
    }
  )

  override def start(): Unit = {
    stage = new PrimaryStage {
      title = "Decoding-Us Navigator"
      scene = new Scene(800, 850) {
        root = new BorderPane {
          top = topBar
          center = mainLayout
        }
        stylesheets.add(getClass.getResource("/style.css").toExternalForm)
      }
    }

    topBar.update(currentUser)
    val welcomeScreen = new WelcomeScreen(
      onFileSelected = filePath => {
        println(s"Dropped file: $filePath")
        startInitialAnalysis(filePath)
      },
      onSelectFileClicked = () => selectFile()
    )
    mainLayout.children = welcomeScreen
  }

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
      val errorScreen = new VBox(20) {
        alignment = Pos.Center
        styleClass.add("root-pane")
        children = Seq(
          new Label("Analysis Failed!") {
            styleClass.add("error-label")
          },
          new Label("Please check the console for more details and ensure the reference and input files are correct.") {
            styleClass.add("info-label")
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
        currentSha256.flatMap(AnalysisCache.load).foreach { cachedSummary =>
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
      currentSha256.foreach { sha =>
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

    val choicesBox = new VBox(10) {
      alignment = Pos.CenterLeft
      padding = Insets(20)
      children = Seq(
        new Label("Select an analysis to perform:") { styleClass.add("sub-title-label") },
        loadCacheButton, // Add the cache button here
        deepAnalysisHelpText,
        deepAnalysisButton,
        new Separator(),
        new Label("Haplogroup Analysis:") { styleClass.add("sub-title-label") },
        haplogroupControls
      )
    }

    val initialResultsScreen = new VBox(20) {
      alignment = Pos.TopCenter
      styleClass.add("root-pane")
      padding = Insets(20)
      children = Seq(resultsTitle, resultsVBox, new Separator(), choicesBox)
    }

    mainLayout.children = initialResultsScreen
  }

  private def startDeepCoverageAnalysis(filePath: String, libraryStats: LibraryStats, referencePath: String): Unit = {
    val progressLabel = new Label("Deep analysis in progress...") {
      styleClass.add("progress-label")
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
        new HBox(20) {
          alignment = Pos.Center
          children = Seq(progressBar, progressIndicator)
        }
      )
    }

    mainLayout.children = progressScreen

    val jfxTask = new jfxc.Task[(CoverageSummary, List[String])]() { // Returns CoverageSummary and svgStrings
      override def call(): (CoverageSummary, List[String]) = {
        try {
          val wgsMetricsProcessor = new WgsMetricsProcessor()
          val callableLociProcessor = new CallableLociProcessor()

          // Phase 2: WGS Metrics
          val wgsMetrics = wgsMetricsProcessor.process(filePath, referencePath, (message, current, total) => {
            Platform.runLater { progressLabel.text = message }
            updateProgress(current, total * 2) // 0-50%
          })

          // Phase 3: Callable Loci Analysis
          val (callableLociResult, svgStrings) = callableLociProcessor.process(filePath, referencePath, (message, current, total) => {
            Platform.runLater { progressLabel.text = s"Callable Loci: $message" }
            updateProgress(total + current, total * 2) // 50-100%
          })

          val summary = CoverageSummary(
            pdsUserId = currentUser.map(_.id).getOrElse("Anonymous"),
            libraryStats = libraryStats,
            wgsMetrics = wgsMetrics,
            callableBases = callableLociResult.callableBases,
            contigAnalysis = callableLociResult.contigAnalysis
          )
          
          // Save to cache if SHA is available
          // Note: We access currentSha256 from the main thread usually, but here we are in background.
          // However, currentSha256 is a simple var. If it's set (which it should be by now for a large file analysis), we use it.
          // If the analysis was faster than the SHA calculation (unlikely for deep analysis), we might miss it.
          // Ideally we should wait for it, but for MVP this is acceptable best-effort.
          if (currentSha256.isDefined) {
             AnalysisCache.save(currentSha256.get, summary)
          }
          
          coverageSummary = Some(summary)
          (summary, svgStrings)
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
      val (summary, svgStrings) = jfxTask.getValue
      showResults(summary, svgStrings) // Show detailed results after deep analysis
    })

    jfxTask.setOnFailed(_ => {
      val errorScreen = new VBox(20) {
        alignment = Pos.Center
        styleClass.add("root-pane")
        children = Seq(
          new Label("Deep Analysis Failed!") {
            styleClass.add("error-label")
          },
          new Label("Please check the console for more details.") {
            styleClass.add("info-label")
          },
          new Button("Back to Choices") {
            onAction = _ => showInitialResultsAndChoices(libraryStats, referencePath) // Go back to choices
          }
        )
      }
      mainLayout.children = errorScreen
    })

    new Thread(jfxTask).start()
  }

  private def showResults(summary: CoverageSummary, svgStrings: List[String]): Unit = {
    val resultsTitle = new Label("Deep Coverage Analysis Results") { // Updated title
      styleClass.add("title-label")
    }

    val statsGrid = new GridPane {
      hgap = 20
      vgap = 10
      padding = Insets(20)
      styleClass.add("stats-grid")
    }

    def addStat(name: String, value: String, row: Int): Unit = {
      statsGrid.add(new Text(name) { styleClass.add("stat-name") }, 0, row)
      statsGrid.add(new Text(value) { styleClass.add("stat-value") }, 1, row)
    }

    addStat("PDS User ID:", summary.pdsUserId, 0)
    addStat("Sample Name:", summary.libraryStats.sampleName, 1)
    addStat("Aligner:", summary.libraryStats.aligner, 2)
    addStat("Platform:", summary.libraryStats.inferredPlatform, 3)
    addStat("Instrument:", summary.libraryStats.mostFrequentInstrument, 4)
    addStat("Reference:", summary.libraryStats.referenceBuild, 5)
    addStat("Genome Size:", f"${summary.wgsMetrics.genomeTerritory}%,d", 6)
    addStat("Mean Coverage:", f"${summary.wgsMetrics.meanCoverage}%.2fx", 7)
    addStat("Callable Bases:", f"${summary.callableBases}%,d", 8)
    val callablePercent = if (summary.wgsMetrics.genomeTerritory > 0) (summary.callableBases.toDouble / summary.wgsMetrics.genomeTerritory * 100) else 0.0
    addStat("Callable Percentage:", f"$callablePercent%.2f%%", 9)

    val contigBreakdownTitle = new Label("🧬 Contig Breakdown") {
      styleClass.add("title-label")
      padding = Insets(20, 0, 10, 0)
    }

    val contigData = ObservableBuffer.from(summary.contigAnalysis.map { contig =>
      ContigAnalysisRow(
        contig.contigName,
        contig.callable,
        contig.noCoverage,
        contig.lowCoverage,
        contig.poorMappingQuality,
        contig.refN,
        s"${contig.contigName}.callable.svg"
      )
    })

    val contigTable = new TableView[ContigAnalysisRow](contigData) {
      columns ++= Seq(
        new TableColumn[ContigAnalysisRow, String]("Contig") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(c.value.contig)
        },
        new TableColumn[ContigAnalysisRow, String]("Callable Bases") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(f"${c.value.callableBases}%,d")
        },
        new TableColumn[ContigAnalysisRow, String]("No Coverage") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(f"${c.value.callableBases}%,d")
        },
        new TableColumn[ContigAnalysisRow, String]("Low Coverage") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(f"${c.value.lowCoverage}%,d")
        },
        new TableColumn[ContigAnalysisRow, String]("Poor MQ") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(f"${c.value.poorMq}%,d")
        },
        new TableColumn[ContigAnalysisRow, String]("REF N") {
          cellValueFactory = c => new scalafx.beans.property.StringProperty(f"${c.value.refN}%,d")
        }
      )
    }

    val webView = new WebView {
      prefHeight = 400
    }

    contigTable.selectionModel.value.selectedItem.onChange { (_, _, newValue) =>
      if (newValue != null) {
        val index = contigTable.items.value.indexOf(newValue)
        if (index >= 0 && index < svgStrings.length) {
          webView.engine.loadContent(svgStrings(index))
        }
      }
    }

    // Haplogroup Analysis Controls
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
              text = "Uploading..."
              disable = true
              PdsClient.uploadSummary(summary).onComplete {
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
            }
          }
        )
      }
      resultsVBox.children.addAll(new Separator(), pdsBox)
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
      summary <- currentLibraryStats
      referencePath <- currentReferencePath
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
}