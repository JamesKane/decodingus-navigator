package com.decodingus.ui.components

import com.decodingus.yprofile.model.*
import com.decodingus.yprofile.repository.YSnpPanelEntity
import com.decodingus.yprofile.service.YProfileService
import scalafx.Includes.*
import scalafx.application.Platform
import scalafx.beans.property.{ObjectProperty, StringProperty}
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.{GridPane, HBox, Priority, Region, VBox}
import scalafx.scene.Node

import java.time.format.DateTimeFormatter
import java.util.UUID

/**
 * Result of profile management actions.
 */
sealed trait YProfileManagementResult
case class ProfileCreated(profileId: UUID) extends YProfileManagementResult
case class ProfileUpdated(profileId: UUID) extends YProfileManagementResult
case class SourceAdded(sourceId: UUID) extends YProfileManagementResult
case class SourceRemoved(sourceId: UUID) extends YProfileManagementResult
case class ReconciliationRun(profileId: UUID) extends YProfileManagementResult
case object DialogClosed extends YProfileManagementResult

/**
 * Dialog for managing Y Chromosome Profile - creation, sources, imports, and reconciliation.
 * Works both when a profile exists and when one needs to be created.
 *
 * @param biosampleId     The biosample UUID
 * @param biosampleName   Display name for the biosample
 * @param yProfileService Service for YProfile CRUD operations
 * @param existingProfile Optional existing profile (None = show create UI)
 * @param snpPanels       Available SNP panels for import
 * @param onRefresh       Callback to refresh profile data after mutations
 */
class YProfileManagementDialog(
                                biosampleId: UUID,
                                biosampleName: String,
                                yProfileService: YProfileService,
                                existingProfile: Option[YChromosomeProfileEntity],
                                sources: List[YProfileSourceEntity],
                                variants: List[YProfileVariantEntity],
                                snpPanels: List[YSnpPanelEntity],
                                onRefresh: () => Unit
                              ) extends Dialog[YProfileManagementResult] {

  private val dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm")

  // Track current profile state (may change during dialog lifetime)
  private val currentProfile = ObjectProperty(existingProfile)
  private val currentSources = ObservableBuffer.from(sources)
  private val statusMessage = StringProperty("")

  title = "Y Chromosome Profile Management"
  headerText = s"Manage Y Profile for $biosampleName"

  dialogPane().buttonTypes = Seq(ButtonType.Close)
  dialogPane().setPrefSize(900, 650)

  // Build dialog content based on whether profile exists
  private val dialogContent = new VBox(10) {
    padding = Insets(15)
    children = buildContent()
  }
  VBox.setVgrow(dialogContent, Priority.Always)

  dialogPane().content = dialogContent

  // Make dialog resizable
  dialogPane().getScene.getWindow match {
    case stage: javafx.stage.Stage => stage.setResizable(true)
    case _ =>
  }

  // Result converter
  resultConverter = _ => DialogClosed

  // --- Content Building ---

  private def buildContent(): Seq[Node] = {
    currentProfile.value match {
      case Some(profile) => buildProfileManagementUI(profile)
      case None => buildCreateProfileUI()
    }
  }

  private def buildCreateProfileUI(): Seq[Node] = {
    val infoBox = new VBox(15) {
      padding = Insets(30)
      alignment = Pos.Center
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 8;"

      children = Seq(
        new Label("No Y Chromosome Profile") {
          style = "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #888;"
        },
        new Label("A Y Profile allows you to combine Y-DNA data from multiple sources\nand reconcile conflicting calls for accurate haplogroup determination.") {
          style = "-fx-font-size: 14px; -fx-text-fill: #aaa; -fx-text-alignment: center;"
          wrapText = true
        },
        new Region { prefHeight = 20 },
        new Button("Create Y Profile") {
          style = "-fx-font-size: 16px; -fx-padding: 10 30;"
          onAction = _ => handleCreateProfile()
        },
        new Label {
          text <== statusMessage
          style = "-fx-text-fill: #ff6666;"
        }
      )
    }

    Seq(infoBox)
  }

  private def buildProfileManagementUI(profile: YChromosomeProfileEntity): Seq[Node] = {
    val summaryPanel = createSummaryPanel(profile)

    val overviewTab = createOverviewTab(profile)
    val sourcesTab = createSourcesTab(profile)
    val importTab = createImportTab(profile)
    val reconciliationTab = createReconciliationTab(profile)

    val tabPane = new TabPane {
      tabs = Seq(overviewTab, sourcesTab, importTab, reconciliationTab)
    }
    VBox.setVgrow(tabPane, Priority.Always)

    Seq(summaryPanel, tabPane)
  }

  // --- Summary Panel ---

  private def createSummaryPanel(profile: YChromosomeProfileEntity): VBox = {
    val haplogroupColor = if (profile.consensusHaplogroup.isDefined) "#2d5a2d" else "#4a4a4a"

    new VBox(8) {
      padding = Insets(15)
      style = s"-fx-background-color: linear-gradient(to bottom, $haplogroupColor, #1a3a1a); -fx-background-radius: 8;"

      val haplogroupDisplay = profile.consensusHaplogroup match {
        case Some(hg) =>
          val confidenceText = profile.haplogroupConfidence.map(c => f"${c * 100}%.0f%%").getOrElse("")
          new HBox(20) {
            alignment = Pos.CenterLeft
            children = Seq(
              new Label(hg) {
                style = "-fx-font-size: 28px; -fx-font-weight: bold; -fx-text-fill: white;"
              },
              if (confidenceText.nonEmpty)
                new Label(s"Confidence: $confidenceText") {
                  style = "-fx-font-size: 14px; -fx-text-fill: #88ff00; -fx-font-weight: bold;"
                }
              else new Region()
            )
          }
        case None =>
          new Label("Haplogroup Pending") {
            style = "-fx-font-size: 24px; -fx-text-fill: #aaa; -fx-font-style: italic;"
          }
      }

      val statsBox = new HBox(30) {
        children = Seq(
          createStatBox("Variants", profile.totalVariants.toString),
          createStatBox("Confirmed", profile.confirmedCount.toString),
          createStatBox("Conflicts", profile.conflictCount.toString),
          createStatBox("Sources", profile.sourceCount.toString)
        )
      }

      children = Seq(haplogroupDisplay, statsBox)
    }
  }

  private def createStatBox(label: String, value: String): VBox = {
    new VBox(2) {
      alignment = Pos.Center
      children = Seq(
        new Label(value) {
          style = "-fx-font-size: 18px; -fx-font-weight: bold; -fx-text-fill: white;"
        },
        new Label(label) {
          style = "-fx-font-size: 11px; -fx-text-fill: #aaa;"
        }
      )
    }
  }

  // --- Tab Creation ---

  private def createOverviewTab(profile: YChromosomeProfileEntity): Tab = {
    val overviewContent = new VBox(15) {
      padding = Insets(15)

      // Profile metadata grid
      val metadataGrid = new GridPane {
        hgap = 15
        vgap = 8
        padding = Insets(10)
      }

      def addRow(row: Int, labelText: String, valueText: String): Unit = {
        val label = new Label(s"$labelText:") {
          style = "-fx-font-weight: bold; -fx-text-fill: #888;"
        }
        val value = new Label(valueText) {
          style = "-fx-text-fill: #ccc;"
        }
        metadataGrid.add(label, 0, row)
        metadataGrid.add(value, 1, row)
      }

      addRow(0, "Profile ID", profile.id.toString.take(8) + "...")
      addRow(1, "Tree Provider", profile.haplogroupTreeProvider.getOrElse("Not set"))
      addRow(2, "Tree Version", profile.haplogroupTreeVersion.getOrElse("Not set"))
      addRow(3, "Primary Source", profile.primarySourceType.map(_.toString).getOrElse("None"))
      addRow(4, "Created", profile.meta.createdAt.format(dateFormatter))
      addRow(5, "Last Updated", profile.meta.updatedAt.format(dateFormatter))
      addRow(6, "Last Reconciled", profile.lastReconciledAt.map(_.format(dateFormatter)).getOrElse("Never"))

      // Status message
      val statusLabel = new Label {
        text <== statusMessage
        style = "-fx-text-fill: #88ff00;"
      }

      children = Seq(
        new Label("Profile Details") {
          style = "-fx-font-size: 16px; -fx-font-weight: bold;"
        },
        metadataGrid,
        new Region { prefHeight = 20 },
        statusLabel
      )
    }

    new Tab {
      text = "Overview"
      closable = false
      content = new ScrollPane {
        fitToWidth = true
        content = overviewContent
      }
    }
  }

  private def createSourcesTab(profile: YChromosomeProfileEntity): Tab = {
    case class SourceRow(
                          id: UUID,
                          vendor: String,
                          testName: String,
                          sourceType: String,
                          tier: String,
                          variants: Int,
                          importedAt: String
                        )

    val tableData = ObservableBuffer.from(currentSources.map { s =>
      SourceRow(
        s.id,
        s.vendor.getOrElse("-"),
        s.testName.getOrElse("-"),
        s.sourceType.toString,
        s"Tier ${s.methodTier}",
        s.variantCount,
        s.importedAt.format(dateFormatter)
      )
    }.toList)

    val table = new TableView[SourceRow](tableData) {
      columnResizePolicy = TableView.ConstrainedResizePolicy

      columns ++= Seq(
        new TableColumn[SourceRow, String] {
          text = "Vendor"
          cellValueFactory = { r => StringProperty(r.value.vendor) }
          prefWidth = 120
        },
        new TableColumn[SourceRow, String] {
          text = "Test"
          cellValueFactory = { r => StringProperty(r.value.testName) }
          prefWidth = 120
        },
        new TableColumn[SourceRow, String] {
          text = "Type"
          cellValueFactory = { r => StringProperty(r.value.sourceType) }
          prefWidth = 120
        },
        new TableColumn[SourceRow, String] {
          text = "Tier"
          cellValueFactory = { r => StringProperty(r.value.tier) }
          prefWidth = 60
        },
        new TableColumn[SourceRow, Int] {
          text = "Variants"
          cellValueFactory = { r => ObjectProperty(r.value.variants) }
          prefWidth = 70
        },
        new TableColumn[SourceRow, String] {
          text = "Imported"
          cellValueFactory = { r => StringProperty(r.value.importedAt) }
          prefWidth = 130
        }
      )
    }
    VBox.setVgrow(table, Priority.Always)

    // Action buttons
    val addButton = new Button("Add Source...") {
      onAction = _ => handleAddSource(profile)
    }

    val removeButton = new Button("Remove Selected") {
      disable = true
      onAction = _ => {
        Option(table.selectionModel.value.getSelectedItem).foreach { row =>
          handleRemoveSource(profile, row.id)
        }
      }
    }

    // Enable remove button when selection changes
    table.selectionModel.value.selectedItemProperty.onChange { (_, _, selected) =>
      removeButton.disable = selected == null
    }

    val buttonBar = new HBox(10) {
      padding = Insets(5)
      children = Seq(addButton, removeButton, new Region { hgrow = Priority.Always })
    }

    new Tab {
      text = s"Sources (${currentSources.size})"
      closable = false
      content = new VBox(5) {
        children = Seq(buttonBar, table)
        VBox.setVgrow(table, Priority.Always)
      }
    }
  }

  private def createImportTab(profile: YChromosomeProfileEntity): Tab = {
    val importContent = new VBox(20) {
      padding = Insets(15)

      // Import from SNP Panel section
      val snpPanelSection = new VBox(10) {
        val panelCombo = new ComboBox[String] {
          promptText = "Select an SNP Panel..."
          items = ObservableBuffer.from(
            snpPanels.map(p => s"${p.provider.getOrElse("Unknown")} - ${p.panelName.getOrElse("Panel")} (${p.snpCalls.size} SNPs)")
          )
          prefWidth = 400
        }

        val importPanelButton = new Button("Import Panel") {
          disable = true
          onAction = _ => {
            val selectedIdx = panelCombo.selectionModel.value.getSelectedIndex
            if (selectedIdx >= 0 && selectedIdx < snpPanels.size) {
              handleImportSnpPanel(profile, snpPanels(selectedIdx))
            }
          }
        }

        panelCombo.selectionModel.value.selectedIndexProperty.onChange { (_, _, idx) =>
          importPanelButton.disable = idx.intValue() < 0
        }

        children = Seq(
          new Label("Import from SNP Panel") {
            style = "-fx-font-size: 14px; -fx-font-weight: bold;"
          },
          new Label("Import Y-DNA SNP calls from an existing panel test.") {
            style = "-fx-text-fill: #888;"
          },
          new HBox(10) {
            children = Seq(panelCombo, importPanelButton)
          }
        )
      }

      // Import from haplogroup analysis section
      val analysisSection = new VBox(10) {
        children = Seq(
          new Label("Import from Haplogroup Analysis") {
            style = "-fx-font-size: 14px; -fx-font-weight: bold;"
          },
          new Label("Y-DNA variants are automatically imported when you run haplogroup analysis.") {
            style = "-fx-text-fill: #888;"
          },
          new Label("Use the Analysis menu to run Y-DNA haplogroup determination from a BAM/CRAM file.") {
            style = "-fx-text-fill: #aaa; -fx-font-style: italic;"
          }
        )
      }

      // Status message
      val statusLabel = new Label {
        text <== statusMessage
        style = "-fx-text-fill: #88ff00;"
      }

      children = Seq(
        snpPanelSection,
        new Separator(),
        analysisSection,
        new Region { vgrow = Priority.Always },
        statusLabel
      )
    }

    new Tab {
      text = "Import"
      closable = false
      content = importContent
    }
  }

  private def createReconciliationTab(profile: YChromosomeProfileEntity): Tab = {
    // Show conflicts
    val conflicts = variants.filter(_.status == YVariantStatus.CONFLICT)

    val reconcileContent = new VBox(15) {
      padding = Insets(15)

      // Reconciliation status
      val statusSection = new VBox(10) {
        val conflictCountLabel = if (conflicts.isEmpty) {
          new Label("No conflicts detected") {
            style = "-fx-font-size: 14px; -fx-text-fill: #88ff00;"
          }
        } else {
          new Label(s"${conflicts.size} conflicts requiring review") {
            style = "-fx-font-size: 14px; -fx-text-fill: #ff6666;"
          }
        }

        children = Seq(
          new Label("Reconciliation Status") {
            style = "-fx-font-size: 16px; -fx-font-weight: bold;"
          },
          conflictCountLabel
        )
      }

      // Run reconciliation button
      val runButton = new Button("Run Reconciliation") {
        style = "-fx-padding: 10 20;"
        onAction = _ => handleRunReconciliation(profile)
      }

      // Conflict list
      val conflictList = if (conflicts.nonEmpty) {
        val listView = new ListView[String] {
          items = ObservableBuffer.from(conflicts.map { v =>
            s"${v.variantName.getOrElse(s"pos:${v.position}")} at ${v.position}"
          })
          prefHeight = 200
        }
        Some(new VBox(10) {
          children = Seq(
            new Label("Conflicting Variants") {
              style = "-fx-font-size: 14px; -fx-font-weight: bold;"
            },
            listView
          )
        })
      } else None

      // Status message
      val statusLabel = new Label {
        text <== statusMessage
        style = "-fx-text-fill: #88ff00;"
      }

      children = Seq(statusSection, runButton) ++ conflictList.toSeq ++ Seq(
        new Region { vgrow = Priority.Always },
        statusLabel
      )
    }

    new Tab {
      text = s"Reconciliation${if (conflicts.nonEmpty) s" (${conflicts.size})" else ""}"
      closable = false
      content = reconcileContent
    }
  }

  // --- Action Handlers ---

  private def handleCreateProfile(): Unit = {
    statusMessage.value = "Creating profile..."
    yProfileService.getOrCreateProfile(biosampleId) match {
      case Right(profile) =>
        statusMessage.value = ""
        currentProfile.value = Some(profile)
        // Rebuild dialog content
        dialogContent.children.setAll(buildContent().map(_.delegate): _*)
        onRefresh()
      case Left(error) =>
        statusMessage.value = s"Error: $error"
    }
  }

  private def handleAddSource(profile: YChromosomeProfileEntity): Unit = {
    // Show add source dialog
    val dialog = new AddYProfileSourceDialog(profile.id, yProfileService)
    dialog.showAndWait() match {
      case Some(Some(source: YProfileSourceEntity)) =>
        currentSources.add(source)
        statusMessage.value = s"Added source: ${source.vendor.getOrElse("Unknown")}"
        onRefresh()
      case _ => // Cancelled
    }
  }

  private def handleRemoveSource(profile: YChromosomeProfileEntity, sourceId: UUID): Unit = {
    val confirm = new Alert(Alert.AlertType.Confirmation) {
      title = "Remove Source"
      headerText = "Remove this source from the profile?"
      contentText = "All variant calls from this source will be removed. This cannot be undone."
    }

    confirm.showAndWait() match {
      case Some(ButtonType.OK) =>
        yProfileService.removeSource(sourceId) match {
          case Right(_) =>
            currentSources.removeIf(_.id == sourceId)
            statusMessage.value = "Source removed"
            onRefresh()
          case Left(error) =>
            statusMessage.value = s"Error: $error"
        }
      case _ => // Cancelled
    }
  }

  private def handleImportSnpPanel(profile: YChromosomeProfileEntity, panel: YSnpPanelEntity): Unit = {
    statusMessage.value = "Importing SNP panel..."

    yProfileService.importFromSnpPanel(profile.id, panel) match {
      case Right(result) =>
        statusMessage.value = s"Imported ${result.totalImported} variants from panel"
        onRefresh()
      case Left(error) =>
        statusMessage.value = s"Error: $error"
    }
  }

  private def handleRunReconciliation(profile: YChromosomeProfileEntity): Unit = {
    statusMessage.value = "Running reconciliation..."

    yProfileService.reconcileProfile(profile.id) match {
      case Right(count) =>
        statusMessage.value = s"Reconciliation complete. Processed $count variants."
        onRefresh()
      case Left(error) =>
        statusMessage.value = s"Error: $error"
    }
  }
}

/**
 * Simple dialog for adding a Y Profile source manually.
 */
class AddYProfileSourceDialog(
                               profileId: UUID,
                               yProfileService: YProfileService
                             ) extends Dialog[Option[YProfileSourceEntity]] {

  title = "Add Y Profile Source"
  headerText = "Add a new data source to the Y Profile"

  dialogPane().buttonTypes = Seq(ButtonType.OK, ButtonType.Cancel)

  // Form fields
  private val sourceTypeCombo = new ComboBox[String] {
    items = ObservableBuffer(
      "WGS Short Read",
      "WGS Long Read",
      "Targeted NGS",
      "Chip/Array",
      "Sanger Sequencing",
      "Manual Entry"
    )
    selectionModel.value.selectFirst()
    prefWidth = 200
  }

  private val vendorField = new TextField {
    promptText = "e.g., FTDNA, 23andMe, Internal"
    prefWidth = 200
  }

  private val testNameField = new TextField {
    promptText = "e.g., Big Y-700, WGS 30x"
    prefWidth = 200
  }

  private val referenceBuildCombo = new ComboBox[String] {
    items = ObservableBuffer("GRCh38", "GRCh37", "CHM13v2")
    selectionModel.value.selectFirst()
    prefWidth = 200
  }

  // Layout
  private val formGrid = new GridPane {
    hgap = 10
    vgap = 10
    padding = Insets(20, 100, 10, 10)

    add(new Label("Source Type:"), 0, 0)
    add(sourceTypeCombo, 1, 0)
    add(new Label("Vendor:"), 0, 1)
    add(vendorField, 1, 1)
    add(new Label("Test Name:"), 0, 2)
    add(testNameField, 1, 2)
    add(new Label("Reference Build:"), 0, 3)
    add(referenceBuildCombo, 1, 3)
  }

  dialogPane().content = formGrid

  resultConverter = btn => {
    if (btn == ButtonType.OK) {
      val sourceType = sourceTypeCombo.value.value match {
        case "WGS Short Read" => YProfileSourceType.WGS_SHORT_READ
        case "WGS Long Read" => YProfileSourceType.WGS_LONG_READ
        case "Targeted NGS" => YProfileSourceType.TARGETED_NGS
        case "Chip/Array" => YProfileSourceType.CHIP
        case "Sanger Sequencing" => YProfileSourceType.SANGER
        case _ => YProfileSourceType.MANUAL
      }

      yProfileService.addSource(
        profileId = profileId,
        sourceType = sourceType,
        vendor = Option(vendorField.text.value).filter(_.nonEmpty),
        testName = Option(testNameField.text.value).filter(_.nonEmpty),
        referenceBuild = Some(referenceBuildCombo.value.value)
      ) match {
        case Right(source) => Some(source)
        case Left(_) => None
      }
    } else None
  }
}
