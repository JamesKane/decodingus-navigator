package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{TableView, TableColumn, Button, Alert, ButtonType, ContextMenu, MenuItem, Tooltip, Dialog, Label, TextArea}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.layout.{VBox, HBox, Priority}
import scalafx.geometry.{Insets, Pos}
import scalafx.collections.ObservableBuffer
import scalafx.beans.property.StringProperty
import scalafx.application.Platform
import scalafx.stage.FileChooser
import com.decodingus.config.FeatureToggles
import com.decodingus.haplogroup.tree.TreeType
import com.decodingus.workspace.model.{Biosample, ChipProfile}
import com.decodingus.workspace.WorkbenchViewModel

/**
 * Table component displaying chip/array genotype data for a subject.
 * Supports importing raw data from 23andMe, AncestryDNA, FTDNA, MyHeritage, and LivingDNA.
 */
class ChipDataTable(
  viewModel: WorkbenchViewModel,
  subject: Biosample,
  chipProfiles: List[ChipProfile],
  onRemove: (String) => Unit  // Callback when remove is clicked, passes profile URI
) extends VBox(10) {

  padding = Insets(10, 0, 0, 0)

  // Convert the profiles to an observable buffer
  case class ChipProfileRow(profile: ChipProfile)

  private val tableData: ObservableBuffer[ChipProfileRow] = ObservableBuffer.from(
    chipProfiles.map(ChipProfileRow.apply)
  )

  private val table = new TableView[ChipProfileRow](tableData) {
    prefHeight = 120
    columnResizePolicy = TableView.ConstrainedResizePolicy

    // Vendor column
    columns += new TableColumn[ChipProfileRow, String] {
      text = "Vendor"
      cellValueFactory = { row =>
        StringProperty(row.value.profile.vendor)
      }
      prefWidth = 90
    }

    // Test Type column
    columns += new TableColumn[ChipProfileRow, String] {
      text = "Test"
      cellValueFactory = { row =>
        val version = row.value.profile.chipVersion.map(v => s" v$v").getOrElse("")
        StringProperty(row.value.profile.vendor + version)
      }
      prefWidth = 110
    }

    // Markers column
    columns += new TableColumn[ChipProfileRow, String] {
      text = "Markers"
      cellValueFactory = { row =>
        val count = row.value.profile.totalMarkersCalled
        val display = if (count >= 1000) f"${count / 1000.0}%.0fk" else count.toString
        StringProperty(display)
      }
      prefWidth = 60
    }

    // Call Rate column
    columns += new TableColumn[ChipProfileRow, String] {
      text = "Call Rate"
      cellValueFactory = { row =>
        val rate = row.value.profile.callRate * 100
        StringProperty(f"$rate%.1f%%")
      }
      prefWidth = 70
    }

    // Y Markers column
    columns += new TableColumn[ChipProfileRow, String] {
      text = "Y"
      cellValueFactory = { row =>
        val yMarkers = row.value.profile.yMarkersCalled.map(_.toString).getOrElse("—")
        StringProperty(yMarkers)
      }
      prefWidth = 40
    }

    // MT Markers column
    columns += new TableColumn[ChipProfileRow, String] {
      text = "MT"
      cellValueFactory = { row =>
        val mtMarkers = row.value.profile.mtMarkersCalled.map(_.toString).getOrElse("—")
        StringProperty(mtMarkers)
      }
      prefWidth = 40
    }

    // Status column
    columns += new TableColumn[ChipProfileRow, String] {
      text = "Status"
      cellValueFactory = { row =>
        StringProperty(row.value.profile.status)
      }
      prefWidth = 70
    }

    // Context menu for row actions
    rowFactory = { _ =>
      val row = new javafx.scene.control.TableRow[ChipProfileRow]()
      val contextMenu = new ContextMenu(
        new MenuItem("Details") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              showDetailsDialog(item.profile)
            }
          }
        },
        new MenuItem("Run Ancestry Analysis") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              if (item.profile.isAcceptableForAncestry) {
                handleAncestryAnalysis(item.profile)
              } else {
                new Alert(AlertType.Warning) {
                  title = "Quality Warning"
                  headerText = "Chip data may not be suitable for ancestry analysis"
                  contentText = s"Call rate: ${f"${item.profile.callRate * 100}%.1f"}%, " +
                    s"Autosomal markers: ${item.profile.autosomalMarkersCalled}. " +
                    "Ancestry analysis requires >95% call rate and >100k autosomal markers."
                }.showAndWait()
              }
            }
          }
        },
        new MenuItem("Y-DNA Haplogroup") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              val yMarkers = item.profile.yMarkersCalled.getOrElse(0)
              if (yMarkers >= FeatureToggles.chipData.minYMarkers) {
                handleHaplogroupAnalysis(item.profile, TreeType.YDNA)
              } else {
                new Alert(AlertType.Warning) {
                  title = "Insufficient Y-DNA Markers"
                  headerText = "Not enough Y-DNA markers for haplogroup analysis"
                  contentText = s"Found $yMarkers Y-DNA markers, minimum ${FeatureToggles.chipData.minYMarkers} required."
                }.showAndWait()
              }
            }
          }
        },
        new MenuItem("mtDNA Haplogroup") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              val mtMarkers = item.profile.mtMarkersCalled.getOrElse(0)
              if (mtMarkers >= FeatureToggles.chipData.minMtMarkers) {
                handleHaplogroupAnalysis(item.profile, TreeType.MTDNA)
              } else {
                new Alert(AlertType.Warning) {
                  title = "Insufficient mtDNA Markers"
                  headerText = "Not enough mtDNA markers for haplogroup analysis"
                  contentText = s"Found $mtMarkers mtDNA markers, minimum ${FeatureToggles.chipData.minMtMarkers} required."
                }.showAndWait()
              }
            }
          }
        },
        new MenuItem("Remove") {
          onAction = _ => {
            Option(row.getItem).foreach { item =>
              val confirm = new Alert(AlertType.Confirmation) {
                title = "Remove Chip Data"
                headerText = "Remove this chip/array data?"
                contentText = s"Vendor: ${item.profile.vendor}, ${item.profile.totalMarkersCalled} markers"
              }
              confirm.showAndWait() match {
                case Some(ButtonType.OK) =>
                  item.profile.atUri.foreach(onRemove)
                case _ =>
              }
            }
          }
        }
      )
      row.contextMenu = contextMenu
      row
    }
  }

  /** Shows a dialog with chip profile details */
  private def showDetailsDialog(profile: ChipProfile): Unit = {
    val dialog = new Dialog[Unit] {
      title = "Chip Data Details"
      headerText = s"${profile.vendor} - ${profile.testTypeCode}"
      dialogPane().buttonTypes = Seq(ButtonType.Close)
      dialogPane().setPrefSize(400, 350)

      val detailsText =
        s"""Vendor: ${profile.vendor}
           |Test Type: ${profile.testTypeCode}
           |Chip Version: ${profile.chipVersion.getOrElse("Unknown")}
           |
           |Total Markers: ${profile.totalMarkersCalled} / ${profile.totalMarkersPossible}
           |Call Rate: ${f"${profile.callRate * 100}%.2f"}%
           |No-Call Rate: ${f"${profile.noCallRate * 100}%.2f"}%
           |
           |Autosomal Markers: ${profile.autosomalMarkersCalled}
           |Y-DNA Markers: ${profile.yMarkersCalled.getOrElse("N/A")}
           |mtDNA Markers: ${profile.mtMarkersCalled.getOrElse("N/A")}
           |Heterozygosity Rate: ${profile.hetRate.map(r => f"${r * 100}%.2f%%").getOrElse("N/A")}
           |
           |Status: ${profile.status}
           |Suitable for Ancestry: ${if (profile.isAcceptableForAncestry) "Yes" else "No"}
           |Sufficient Y Coverage: ${if (profile.hasSufficientYCoverage) "Yes" else "No"}
           |Sufficient MT Coverage: ${if (profile.hasSufficientMtCoverage) "Yes" else "No"}
           |
           |Import Date: ${profile.importDate.toLocalDate}
           |Source File: ${profile.sourceFileName.getOrElse("Unknown")}
           |File Hash: ${profile.sourceFileHash.map(_.take(16) + "...").getOrElse("N/A")}
         """.stripMargin

      val textArea = new TextArea(detailsText) {
        editable = false
        wrapText = true
        style = "-fx-font-family: monospace; -fx-font-size: 12px;"
      }
      VBox.setVgrow(textArea, Priority.Always)

      dialogPane().content = new VBox(10) {
        padding = Insets(10)
        children = Seq(textArea)
      }
    }
    dialog.showAndWait()
  }

  /** Handles launching ancestry analysis */
  private def handleAncestryAnalysis(profile: ChipProfile): Unit = {
    import com.decodingus.ancestry.model.AncestryPanelType

    // Recommend panel type based on marker count
    val recommendedPanel = if (profile.autosomalMarkersCalled >= 500000) {
      AncestryPanelType.GenomeWide
    } else {
      AncestryPanelType.Aims
    }

    val panelLabel = recommendedPanel match {
      case AncestryPanelType.Aims => "AIMs (~5k markers, faster)"
      case AncestryPanelType.GenomeWide => "Genome-wide (~500k markers, detailed)"
    }

    // Confirm with user
    val confirm = new Alert(AlertType.Confirmation) {
      title = "Run Ancestry Analysis"
      headerText = s"Analyze ${profile.vendor} chip data for ancestry"
      contentText = s"""This will estimate population percentages using the $panelLabel panel.

Markers: ${profile.autosomalMarkersCalled}
Call Rate: ${f"${profile.callRate * 100}%.1f"}%

Note: Reference data download may be required on first run."""
    }

    confirm.showAndWait() match {
      case Some(ButtonType.OK) =>
        profile.atUri match {
          case Some(profileUri) =>
            // Show progress dialog
            val progressDialog = new AnalysisProgressDialog(
              "Ancestry Analysis",
              viewModel.analysisProgress,
              viewModel.analysisProgressPercent,
              viewModel.analysisInProgress
            )

            viewModel.runChipAncestryAnalysis(
              subject.sampleAccession,
              profileUri,
              recommendedPanel,
              onComplete = {
                case Right(ancestryResult) =>
                  Platform.runLater {
                    // Show results dialog
                    val resultDialog = new AncestryResultDialog(ancestryResult)
                    resultDialog.showAndWait()
                  }
                case Left(error) =>
                  Platform.runLater {
                    new Alert(AlertType.Error) {
                      title = "Ancestry Analysis Failed"
                      headerText = "Could not complete ancestry analysis"
                      contentText = error
                    }.showAndWait()
                  }
              }
            )

            progressDialog.show()

          case None =>
            new Alert(AlertType.Error) {
              title = "Error"
              headerText = "Invalid chip profile"
              contentText = "Profile has no AT URI."
            }.showAndWait()
        }
      case _ => // User cancelled
    }
  }

  /** Handles launching haplogroup analysis (Y-DNA or mtDNA) */
  private def handleHaplogroupAnalysis(profile: ChipProfile, treeType: TreeType): Unit = {
    import com.decodingus.genotype.processor.ChipHaplogroupAdapter

    val typeName = if (treeType == TreeType.YDNA) "Y-DNA" else "mtDNA"
    val markerCount = treeType match {
      case TreeType.YDNA => profile.yMarkersCalled.getOrElse(0)
      case TreeType.MTDNA => profile.mtMarkersCalled.getOrElse(0)
    }

    // Confirm with user
    val confirm = new Alert(AlertType.Confirmation) {
      title = s"Run $typeName Haplogroup Analysis"
      headerText = s"Analyze ${profile.vendor} chip data for $typeName haplogroup"
      contentText = s"""This will score chip genotypes against the $typeName haplogroup tree.

$typeName Markers: $markerCount

Note: Chip-based haplogroup estimation has limited resolution compared to WGS.
The terminal haplogroup may be upstream of the true assignment."""
    }

    confirm.showAndWait() match {
      case Some(ButtonType.OK) =>
        profile.atUri match {
          case Some(profileUri) =>
            // Show progress dialog
            val progressDialog = new AnalysisProgressDialog(
              s"$typeName Haplogroup Analysis",
              viewModel.analysisProgress,
              viewModel.analysisProgressPercent,
              viewModel.analysisInProgress
            )

            viewModel.runChipHaplogroupAnalysis(
              subject.sampleAccession,
              profileUri,
              treeType,
              onComplete = {
                case Right(haplogroupResult) =>
                  Platform.runLater {
                    // Show results dialog
                    val confidenceDesc = ChipHaplogroupAdapter.confidenceDescription(haplogroupResult.confidence)
                    new Alert(AlertType.Information) {
                      title = s"$typeName Haplogroup Result"
                      headerText = s"$typeName: ${haplogroupResult.topHaplogroup}"
                      contentText = s"""Confidence: $confidenceDesc (${f"${haplogroupResult.confidence * 100}%.0f"}%)
SNPs Matched: ${haplogroupResult.snpsMatched} / ${haplogroupResult.snpsTotal}
Tree Depth: ${haplogroupResult.results.headOption.map(_.depth).getOrElse(0)}

Note: Chip data covers ~${f"${haplogroupResult.snpsMatched.toDouble / haplogroupResult.snpsTotal * 100}%.0f"}% of tree positions.
For higher resolution, consider WGS analysis."""
                    }.showAndWait()
                  }
                case Left(error) =>
                  Platform.runLater {
                    new Alert(AlertType.Error) {
                      title = s"$typeName Haplogroup Analysis Failed"
                      headerText = "Could not complete haplogroup analysis"
                      contentText = error
                    }.showAndWait()
                  }
              }
            )

            progressDialog.show()

          case None =>
            new Alert(AlertType.Error) {
              title = "Error"
              headerText = "Invalid chip profile"
              contentText = "Profile has no AT URI."
            }.showAndWait()
        }
      case _ => // User cancelled
    }
  }

  // Action buttons
  private val importButton = new Button("Import Chip Data") {
    tooltip = Tooltip("Import raw chip data from 23andMe, AncestryDNA, FTDNA, etc.")
    onAction = _ => handleImportChipData()
  }

  /** Handles importing chip data file */
  private def handleImportChipData(): Unit = {
    val fileChooser = new FileChooser() {
      title = "Import Chip Data"
      extensionFilters.addAll(
        new FileChooser.ExtensionFilter("Chip Data Files", Seq("*.txt", "*.csv")),
        new FileChooser.ExtensionFilter("23andMe (*.txt)", Seq("*.txt")),
        new FileChooser.ExtensionFilter("AncestryDNA (*.txt)", Seq("*.txt")),
        new FileChooser.ExtensionFilter("FTDNA CSV (*.csv)", Seq("*.csv")),
        new FileChooser.ExtensionFilter("All Files", Seq("*.*"))
      )
    }

    Option(fileChooser.showOpenDialog(this.scene().getWindow)).foreach { file =>
      // Show progress dialog
      val progressDialog = new AnalysisProgressDialog(
        "Importing Chip Data",
        viewModel.analysisProgress,
        viewModel.analysisProgressPercent,
        viewModel.analysisInProgress
      )

      viewModel.importChipData(
        subject.sampleAccession,
        file,
        onComplete = {
          case Right(chipProfile) =>
            Platform.runLater {
              tableData += ChipProfileRow(chipProfile)
              new Alert(AlertType.Information) {
                title = "Import Complete"
                headerText = s"Successfully imported ${chipProfile.vendor} chip data"
                contentText = s"""Markers: ${chipProfile.totalMarkersCalled}
                                 |Call Rate: ${f"${chipProfile.callRate * 100}%.1f"}%
                                 |Y Markers: ${chipProfile.yMarkersCalled.getOrElse("N/A")}
                                 |MT Markers: ${chipProfile.mtMarkersCalled.getOrElse("N/A")}
                                 |Status: ${chipProfile.status}""".stripMargin
              }.showAndWait()
            }
          case Left(error) =>
            Platform.runLater {
              new Alert(AlertType.Error) {
                title = "Import Failed"
                headerText = "Could not import chip data"
                contentText = error
              }.showAndWait()
            }
        }
      )

      progressDialog.show()
    }
  }

  private val buttonBar = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(importButton)
  }

  children = Seq(
    new scalafx.scene.control.Label("Chip/Array Data:") { style = "-fx-font-weight: bold;" },
    table,
    buttonBar
  )

  VBox.setVgrow(table, Priority.Always)
}
