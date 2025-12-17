package com.decodingus.ui.v2

import com.decodingus.i18n.I18n.{t, bind}
import com.decodingus.i18n.Formatters
import com.decodingus.ui.components.{AddSequenceDataDialog, ConfirmDialog, EditSubjectDialog, SequenceDataInput}
import com.decodingus.ui.v2.BiosampleExtensions.*
import com.decodingus.util.Logger
import com.decodingus.workspace.WorkbenchViewModel
import com.decodingus.workspace.model.Biosample
import scalafx.Includes.*
import scalafx.beans.property.{ObjectProperty, StringProperty}
import scalafx.geometry.{Insets, Pos, Side}
import scalafx.scene.control.*
import scalafx.scene.layout.*

/**
 * Subject detail view with tabbed sections.
 *
 * Tabs:
 * - Overview: Summary of all genetic findings
 * - Y-DNA: Y-chromosome analysis
 * - mtDNA: Mitochondrial analysis
 * - Ancestry: Ancestry composition (future)
 * - IBD Matches: IBD matching results (future)
 * - Data Sources: Raw data management
 */
class SubjectDetailView(viewModel: WorkbenchViewModel) extends VBox {

  private val log = Logger[SubjectDetailView]

  spacing = 0
  styleClass += "subject-detail-view"
  style = "-fx-background-color: #1e1e1e;"

  // ============================================================================
  // State
  // ============================================================================

  private val currentSubject: ObjectProperty[Option[Biosample]] = ObjectProperty(None)

  // ============================================================================
  // Header Section
  // ============================================================================

  private val subjectNameLabel = new Label {
    styleClass += "subject-name"
    style = "-fx-font-size: 20px; -fx-font-weight: bold; -fx-text-fill: #ffffff;"
  }

  private val subjectIdLabel = new Label {
    styleClass += "subject-id"
    style = "-fx-font-size: 12px; -fx-text-fill: #b0b0b0;"
  }

  private val editButton = new Button {
    text = t("action.edit")
    onAction = _ => handleEdit()
  }

  private val deleteButton = new Button {
    text = t("action.delete")
    styleClass += "button-danger"
    onAction = _ => handleDelete()
  }

  private val headerSection = new HBox(15) {
    alignment = Pos.CenterLeft
    padding = Insets(15)
    style = "-fx-background-color: #2a2a2a;"
    children = Seq(
      new VBox(5) {
        children = Seq(subjectNameLabel, subjectIdLabel)
      },
      new Region { hgrow = Priority.Always },
      editButton,
      deleteButton
    )
  }

  // ============================================================================
  // Tab Content Views
  // ============================================================================

  private val overviewTab = createTab("subject.tab.overview", createOverviewContent())
  private val ydnaTab = createTab("subject.tab.ydna", createYdnaContent())
  private val mtdnaTab = createTab("subject.tab.mtdna", createMtdnaContent())
  private val ancestryTab = createTab("subject.tab.ancestry", createAncestryContent())
  private val ibdTab = createTab("subject.tab.ibd", createIbdContent())
  private val dataTab = createTab("subject.tab.data", createDataSourcesContent())

  private val tabPane = new TabPane {
    tabClosingPolicy = TabPane.TabClosingPolicy.Unavailable
    side = Side.Top
    styleClass += "subject-tab-pane"
    tabs = Seq(overviewTab, ydnaTab, mtdnaTab, ancestryTab, ibdTab, dataTab)
  }

  // ============================================================================
  // Layout
  // ============================================================================

  vgrow = Priority.Always
  children = Seq(headerSection, tabPane)

  VBox.setVgrow(tabPane, Priority.Always)

  // ============================================================================
  // Tab Creation
  // ============================================================================

  private def createTab(i18nKey: String, tabContent: javafx.scene.Node): Tab = {
    val tab = new Tab {
      text <== bind(i18nKey)
      closable = false
    }
    tab.content = tabContent
    tab
  }

  // ============================================================================
  // Overview Tab Content
  // ============================================================================

  private def createOverviewContent(): ScrollPane = {
    val ydnaCard = createHaplogroupCard("haplogroup.ydna.title", "#2d3a2d", "ydna")
    val mtdnaCard = createHaplogroupCard("haplogroup.mtdna.title", "#2d2d3a", "mtdna")

    val haplogroupSection = new HBox(20) {
      padding = Insets(0, 0, 20, 0)
      children = Seq(ydnaCard, mtdnaCard)
    }

    val ancestryCard = createPlaceholderCard("ancestry.title", "ancestry.not_analyzed", "action.analyze")
    val ibdCard = createPlaceholderCard("ibd.matches", "ibd.no_matches", "ibd.run_match")

    val secondarySection = new HBox(20) {
      padding = Insets(0, 0, 20, 0)
      children = Seq(ancestryCard, ibdCard)
    }

    val dataSummarySection = createDataSummarySection()

    new ScrollPane {
      fitToWidth = true
      content = new VBox(20) {
        padding = Insets(20)
        children = Seq(
          createSectionLabel("subject.genetic_summary"),
          haplogroupSection,
          secondarySection,
          createSectionLabel("subject.data_sources_summary"),
          dataSummarySection
        )
      }
    }
  }

  private def createHaplogroupCard(titleKey: String, bgColor: String, dataType: String): VBox = {
    val haplogroupLabel = new Label("-") {
      styleClass += "haplogroup-value"
      style = "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #ffffff;"
      id = s"${dataType}-haplogroup"
    }

    val confidenceLabel = new Label {
      text = ""
      style = "-fx-font-size: 11px; -fx-text-fill: #b0b0b0;"
      id = s"${dataType}-confidence"
    }

    val viewDetailsButton = new Button {
      text = t(if (dataType == "ydna") "haplogroup.view_profile" else "haplogroup.view_details")
      styleClass += "button-link"
      onAction = _ => {
        // Switch to the appropriate tab
        if (dataType == "ydna") tabPane.selectionModel.value.select(ydnaTab)
        else tabPane.selectionModel.value.select(mtdnaTab)
      }
    }

    new VBox(10) {
      padding = Insets(15)
      prefWidth = 220
      style = s"-fx-background-color: $bgColor; -fx-background-radius: 10;"
      children = Seq(
        new Label { text <== bind(titleKey); style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" },
        haplogroupLabel,
        confidenceLabel,
        viewDetailsButton
      )
    }
  }

  private def createPlaceholderCard(titleKey: String, placeholderKey: String, actionKey: String): VBox = {
    new VBox(10) {
      padding = Insets(15)
      prefWidth = 220
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"
      children = Seq(
        new Label { text <== bind(titleKey); style = "-fx-font-weight: bold;" },
        new Label { text <== bind(placeholderKey); style = "-fx-text-fill: #666666;" },
        new Button {
          text <== bind(actionKey)
          styleClass += "button-secondary"
        }
      )
    }
  }

  private def createDataSummarySection(): HBox = {
    new HBox(15) {
      id = "data-summary"
      children = Seq(
        createDataCountBadge("data.sequencing_runs", 0),
        createDataCountBadge("data.chip_profiles", 0),
        createDataCountBadge("data.str_profiles", 0)
      )
    }
  }

  private def createDataCountBadge(labelKey: String, count: Int): HBox = {
    new HBox(5) {
      alignment = Pos.CenterLeft
      padding = Insets(8, 12, 8, 12)
      style = "-fx-background-color: #333333; -fx-background-radius: 5;"
      children = Seq(
        new Label(count.toString) { style = "-fx-font-weight: bold;" },
        new Label { text <== bind(labelKey); style = "-fx-text-fill: #888888;" }
      )
    }
  }

  // ============================================================================
  // Y-DNA Tab Content
  // ============================================================================

  private def createYdnaContent(): ScrollPane = {
    val terminalHaplogroupSection = new VBox(10) {
      padding = Insets(15)
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"
      children = Seq(
        new Label { text <== bind("haplogroup.terminal"); style = "-fx-font-weight: bold;" },
        new Label("-") {
          id = "ydna-terminal"
          style = "-fx-font-size: 28px; -fx-font-weight: bold; -fx-text-fill: #4ade80;"
        },
        new Label {
          id = "ydna-path"
          text = ""
          style = "-fx-text-fill: #888888;"
        },
        new HBox(20) {
          children = Seq(
            createStatLabel("haplogroup.derived", "ydna-derived", "-"),
            createStatLabel("haplogroup.ancestral", "ydna-ancestral", "-"),
            createStatLabel("haplogroup.callable", "ydna-callable", "-")
          )
        }
      )
    }

    val analyzeButton = new Button {
      text <== bind("analysis.run")
      styleClass += "button-primary"
      onAction = _ => handleRunYdnaAnalysis()
    }

    new ScrollPane {
      fitToWidth = true
      content = new VBox(20) {
        padding = Insets(20)
        children = Seq(
          new HBox(10) {
            alignment = Pos.CenterLeft
            children = Seq(
              new Label { text <== bind("haplogroup.ydna.title"); style = "-fx-font-size: 18px; -fx-font-weight: bold;" },
              new Region { hgrow = Priority.Always },
              analyzeButton
            )
          },
          terminalHaplogroupSection
        )
      }
    }
  }

  // ============================================================================
  // mtDNA Tab Content
  // ============================================================================

  private def createMtdnaContent(): ScrollPane = {
    val haplogroupSection = new VBox(10) {
      padding = Insets(15)
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"
      children = Seq(
        new Label { text <== bind("haplogroup.title"); style = "-fx-font-weight: bold;" },
        new Label("-") {
          id = "mtdna-terminal"
          style = "-fx-font-size: 28px; -fx-font-weight: bold; -fx-text-fill: #60a5fa;"
        },
        new Label {
          id = "mtdna-path"
          text = ""
          style = "-fx-text-fill: #888888;"
        }
      )
    }

    val analyzeButton = new Button {
      text <== bind("analysis.run")
      styleClass += "button-primary"
      onAction = _ => handleRunMtdnaAnalysis()
    }

    new ScrollPane {
      fitToWidth = true
      content = new VBox(20) {
        padding = Insets(20)
        children = Seq(
          new HBox(10) {
            alignment = Pos.CenterLeft
            children = Seq(
              new Label { text <== bind("haplogroup.mtdna.title"); style = "-fx-font-size: 18px; -fx-font-weight: bold;" },
              new Region { hgrow = Priority.Always },
              analyzeButton
            )
          },
          haplogroupSection
        )
      }
    }
  }

  // ============================================================================
  // Ancestry Tab Content
  // ============================================================================

  private def createAncestryContent(): VBox = {
    new VBox(20) {
      padding = Insets(20)
      alignment = Pos.Center
      children = Seq(
        new Label { text <== bind("ancestry.title"); style = "-fx-font-size: 18px; -fx-font-weight: bold;" },
        new Label { text <== bind("ancestry.not_analyzed"); style = "-fx-text-fill: #666666;" },
        new Button {
          text <== bind("analysis.run")
          styleClass += "button-primary"
        }
      )
    }
  }

  // ============================================================================
  // IBD Tab Content
  // ============================================================================

  private def createIbdContent(): VBox = {
    new VBox(20) {
      padding = Insets(20)
      alignment = Pos.Center
      children = Seq(
        new Label { text <== bind("ibd.title"); style = "-fx-font-size: 18px; -fx-font-weight: bold;" },
        new Label { text <== bind("ibd.no_matches"); style = "-fx-text-fill: #666666;" },
        new Button {
          text <== bind("ibd.run_match")
          styleClass += "button-primary"
        }
      )
    }
  }

  // ============================================================================
  // Data Sources Tab Content
  // ============================================================================

  private def createDataSourcesContent(): ScrollPane = {
    val sequencingSection = createDataSection("data.sequencing_runs", "sequencing-list")
    val chipSection = createDataSection("data.chip_profiles", "chip-list")
    val strSection = createDataSection("data.str_profiles", "str-list")

    val addDataButton = new Button {
      text <== bind("data.add")
      styleClass += "button-primary"
      onAction = _ => handleAddData()
    }

    new ScrollPane {
      fitToWidth = true
      content = new VBox(20) {
        padding = Insets(20)
        children = Seq(
          new HBox(10) {
            alignment = Pos.CenterLeft
            children = Seq(
              new Label { text <== bind("data.title"); style = "-fx-font-size: 18px; -fx-font-weight: bold;" },
              new Region { hgrow = Priority.Always },
              addDataButton
            )
          },
          sequencingSection,
          chipSection,
          strSection
        )
      }
    }
  }

  private def createDataSection(titleKey: String, listId: String): VBox = {
    val placeholder = new Label {
      text <== bind(titleKey match {
        case "data.sequencing_runs" => "data.no_sequencing"
        case "data.chip_profiles" => "data.no_chip"
        case _ => "data.no_str"
      })
      style = "-fx-text-fill: #666666;"
    }

    new VBox(10) {
      padding = Insets(15)
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"
      children = Seq(
        new Label { text <== bind(titleKey); style = "-fx-font-weight: bold;" },
        new VBox {
          id = listId
          children = Seq(placeholder)
        }
      )
    }
  }

  // ============================================================================
  // Helper Methods
  // ============================================================================

  private def createSectionLabel(key: String): Label = {
    new Label {
      text <== bind(key)
      style = "-fx-font-size: 14px; -fx-font-weight: bold; -fx-text-fill: #aaaaaa;"
    }
  }

  private def createStatLabel(labelKey: String, valueId: String, defaultValue: String): VBox = {
    new VBox(2) {
      children = Seq(
        new Label { text <== bind(labelKey); style = "-fx-font-size: 11px; -fx-text-fill: #888888;" },
        new Label(defaultValue) { id = valueId; style = "-fx-font-weight: bold;" }
      )
    }
  }

  private def updateContent(subject: Biosample): Unit = {
    // Update header
    subjectNameLabel.text = subject.donorId.getOrElse(subject.accession)
    subjectIdLabel.text = s"ID: ${subject.accession} • ${formatSex(subject.sex)}"

    // Update Overview tab - Y-DNA card
    lookupLabel("ydna-haplogroup").foreach(_.text = subject.yHaplogroup.getOrElse("-"))

    // Update Overview tab - mtDNA card
    lookupLabel("mtdna-haplogroup").foreach(_.text = subject.mtHaplogroup.getOrElse("-"))

    // Update Y-DNA tab
    lookupLabel("ydna-terminal").foreach(_.text = subject.yHaplogroup.getOrElse("-"))

    // Update mtDNA tab
    lookupLabel("mtdna-terminal").foreach(_.text = subject.mtHaplogroup.getOrElse("-"))

    // TODO: Load and display sequence runs, chip profiles, STR profiles from workspace
  }

  private def lookupLabel(labelId: String): Option[Label] = {
    val node = this.lookup(s"#$labelId")
    if (node != null) {
      node.delegate match {
        case l: javafx.scene.control.Label => Some(new Label(l))
        case _ => None
      }
    } else None
  }

  private def formatSex(sex: Option[String]): String = {
    sex.map(_.toUpperCase) match {
      case Some("M") | Some("MALE") => t("sex.male")
      case Some("F") | Some("FEMALE") => t("sex.female")
      case _ => t("sex.unknown")
    }
  }

  // ============================================================================
  // Action Handlers
  // ============================================================================

  private def handleEdit(): Unit = {
    currentSubject.value.foreach { subject =>
      val dialog = new EditSubjectDialog(subject)
      dialog.showAndWait() match {
        case Some(Some(updatedSubject: Biosample)) =>
          viewModel.updateSubject(updatedSubject)
          currentSubject.value = Some(updatedSubject)
          updateContent(updatedSubject)
          log.info(s"Updated subject: ${updatedSubject.sampleAccession}")
        case _ =>
          log.debug("Edit subject cancelled")
      }
    }
  }

  private def handleDelete(): Unit = {
    currentSubject.value.foreach { subject =>
      val subjectName = subject.donorId.getOrElse(subject.accession)
      if (ConfirmDialog.confirmRemoval("subject", subjectName)) {
        viewModel.deleteSubject(subject.accession)
        currentSubject.value = None
        log.info(s"Deleted subject: ${subject.accession}")
      }
    }
  }

  private def handleRunYdnaAnalysis(): Unit = {
    currentSubject.value.foreach { subject =>
      // Analysis requires sequence data - check if available
      if (!subject.hasSequenceData) {
        showInfoDialog(
          t("analysis.title"),
          t("data.no_sequencing"),
          t("data.add_sequence_first")
        )
      } else {
        // TODO: Integrate with existing HaplogroupAnalysisDialog and analysis flow
        log.debug(s"Run Y-DNA analysis for: ${subject.accession} - not yet integrated")
        showInfoDialog(
          t("haplogroup.ydna.title"),
          t("analysis.not_integrated"),
          t("analysis.use_main_workflow")
        )
      }
    }
  }

  private def handleRunMtdnaAnalysis(): Unit = {
    currentSubject.value.foreach { subject =>
      if (!subject.hasSequenceData) {
        showInfoDialog(
          t("analysis.title"),
          t("data.no_sequencing"),
          t("data.add_sequence_first")
        )
      } else {
        // TODO: Integrate with existing HaplogroupAnalysisDialog and analysis flow
        log.debug(s"Run mtDNA analysis for: ${subject.accession} - not yet integrated")
        showInfoDialog(
          t("haplogroup.mtdna.title"),
          t("analysis.not_integrated"),
          t("analysis.use_main_workflow")
        )
      }
    }
  }

  private def handleAddData(): Unit = {
    currentSubject.value.foreach { subject =>
      // For now, pass empty set - duplicate detection can be enhanced later
      // when we have full sequence run resolution from refs
      val existingChecksums: Set[String] = Set.empty

      val dialog = new AddSequenceDataDialog(existingChecksums)
      dialog.showAndWait() match {
        case Some(Some(dataInput: SequenceDataInput)) =>
          // TODO: Process the data input and add to subject
          log.info(s"Adding data for ${subject.accession}: ${dataInput.fileInfo.fileName}")
          showInfoDialog(
            t("data.title"),
            t("data.processing"),
            t("data.processing.detail")
          )
        case _ =>
          log.debug("Add data cancelled")
      }
    }
  }

  private def showInfoDialog(dialogTitle: String, dialogHeader: String, dialogContent: String): Unit = {
    import scalafx.scene.control.Alert
    import scalafx.scene.control.Alert.AlertType
    new Alert(AlertType.Information) {
      title = dialogTitle
      headerText = dialogHeader
      contentText = dialogContent
    }.showAndWait()
  }

  // ============================================================================
  // Public API
  // ============================================================================

  /**
   * Set the subject to display.
   */
  def setSubject(subject: Biosample): Unit = {
    currentSubject.value = Some(subject)
    updateContent(subject)
  }

  /**
   * Clear the current subject.
   */
  def clearSubject(): Unit = {
    currentSubject.value = None
    subjectNameLabel.text = ""
    subjectIdLabel.text = ""
  }
}
