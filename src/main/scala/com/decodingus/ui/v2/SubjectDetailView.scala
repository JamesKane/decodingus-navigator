package com.decodingus.ui.v2

import com.decodingus.i18n.I18n.{t, bind}
import com.decodingus.i18n.Formatters
import com.decodingus.str.StrCsvParser
import com.decodingus.config.FeatureToggles
import com.decodingus.haplogroup.tree.TreeType
import com.decodingus.ui.components.{AddDataDialog, AddSequenceDataDialog, AnalysisProgressDialog, AncestryResultDialog, ConfirmDialog, DataInput, DataType, EditSubjectDialog, InfoDialog, SequenceDataInput, VcfMetadata, VcfMetadataDialog}
import com.decodingus.ui.v2.BiosampleExtensions.*
import com.decodingus.util.Logger
import com.decodingus.workspace.WorkbenchViewModel
import com.decodingus.workspace.model.{Alignment, Biosample, ChipProfile, HaplogroupResult, SequenceRun, StrProfile}
import scalafx.scene.control.Alert.AlertType
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
  // UI Labels (must be declared before tab creation)
  // ============================================================================

  // Overview tab labels
  private val overviewYdnaHaplogroupLabel = new Label("-") {
    styleClass += "haplogroup-value"
    style = "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #ffffff;"
  }
  private val overviewYdnaConfidenceLabel = new Label {
    text = ""
    style = "-fx-font-size: 11px; -fx-text-fill: #b0b0b0;"
  }
  private val overviewMtdnaHaplogroupLabel = new Label("-") {
    styleClass += "haplogroup-value"
    style = "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #ffffff;"
  }
  private val overviewMtdnaConfidenceLabel = new Label {
    text = ""
    style = "-fx-font-size: 11px; -fx-text-fill: #b0b0b0;"
  }
  private val sequencingCountLabel = new Label("0") { style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" }
  private val chipCountLabel = new Label("0") { style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" }
  private val strCountLabel = new Label("0") { style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" }

  // Y-DNA tab labels
  private val ydnaTerminalLabel = new Label("-") {
    id = "ydna-terminal"
    style = "-fx-font-size: 28px; -fx-font-weight: bold; -fx-text-fill: #4ade80;"
  }
  private val ydnaPathLabel = new Label {
    id = "ydna-path"
    text = ""
    style = "-fx-text-fill: #b0b0b0; -fx-font-size: 12px;"
    wrapText = true
  }
  private val ydnaDerivedLabel = new Label("-") { id = "ydna-derived"; style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" }
  private val ydnaAncestralLabel = new Label("-") { id = "ydna-ancestral"; style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" }
  private val ydnaConfidenceLabel = new Label("-") { id = "ydna-confidence"; style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" }
  private val ydnaSourceLabel = new Label("-") { id = "ydna-source"; style = "-fx-text-fill: #888888;" }
  private val ydnaQualityLabel = new Label("-") { id = "ydna-quality"; style = "-fx-font-weight: bold;" }
  private val ydnaNotAnalyzedPane = new VBox(10) {
    alignment = Pos.Center
    padding = Insets(40)
    children = Seq(
      new Label { text <== bind("haplogroup.not_determined"); style = "-fx-font-size: 16px; -fx-text-fill: #888888;" },
      new Label { text <== bind("data.add_sequence_first"); style = "-fx-text-fill: #666666;" }
    )
  }
  private val ydnaResultPane = new VBox(15) {
    padding = Insets(0)
  }

  // mtDNA tab labels
  private val mtdnaTerminalLabel = new Label("-") {
    id = "mtdna-terminal"
    style = "-fx-font-size: 28px; -fx-font-weight: bold; -fx-text-fill: #60a5fa;"
  }
  private val mtdnaPathLabel = new Label {
    id = "mtdna-path"
    text = ""
    style = "-fx-text-fill: #b0b0b0; -fx-font-size: 12px;"
    wrapText = true
  }
  private val mtdnaConfidenceLabel = new Label("-") { id = "mtdna-confidence"; style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" }
  private val mtdnaSourceLabel = new Label("-") { id = "mtdna-source"; style = "-fx-text-fill: #888888;" }
  private val mtdnaQualityLabel = new Label("-") { id = "mtdna-quality"; style = "-fx-font-weight: bold;" }
  private val mtdnaNotAnalyzedPane = new VBox(10) {
    alignment = Pos.Center
    padding = Insets(40)
    children = Seq(
      new Label { text <== bind("haplogroup.not_determined"); style = "-fx-font-size: 16px; -fx-text-fill: #888888;" },
      new Label { text <== bind("data.add_sequence_first"); style = "-fx-text-fill: #666666;" }
    )
  }
  private val mtdnaResultPane = new VBox(15) {
    padding = Insets(0)
  }

  // Data Sources tab containers
  private val sequencingListContainer = new VBox(8) {
    id = "sequencing-list"
  }
  private val chipListContainer = new VBox(8) {
    id = "chip-list"
  }
  private val strListContainer = new VBox(8) {
    id = "str-list"
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
    val ydnaCard = createHaplogroupCard("haplogroup.ydna.title", "#2d3a2d", "ydna",
      overviewYdnaHaplogroupLabel, overviewYdnaConfidenceLabel)
    val mtdnaCard = createHaplogroupCard("haplogroup.mtdna.title", "#2d2d3a", "mtdna",
      overviewMtdnaHaplogroupLabel, overviewMtdnaConfidenceLabel)

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
      style = "-fx-background: #1e1e1e; -fx-background-color: #1e1e1e;"
      content = new VBox(20) {
        padding = Insets(20)
        style = "-fx-background-color: #1e1e1e;"
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

  private def createHaplogroupCard(
    titleKey: String,
    bgColor: String,
    dataType: String,
    haplogroupLabel: Label,
    confidenceLabel: Label
  ): VBox = {
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
        new Label { text <== bind(titleKey); style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" },
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
        createDataCountBadge("data.sequencing_runs", sequencingCountLabel),
        createDataCountBadge("data.chip_profiles", chipCountLabel),
        createDataCountBadge("data.str_profiles", strCountLabel)
      )
    }
  }

  private def createDataCountBadge(labelKey: String, countLabel: Label): HBox = {
    new HBox(5) {
      alignment = Pos.CenterLeft
      padding = Insets(8, 12, 8, 12)
      style = "-fx-background-color: #333333; -fx-background-radius: 5;"
      children = Seq(
        countLabel,
        new Label { text <== bind(labelKey); style = "-fx-text-fill: #888888;" }
      )
    }
  }

  // ============================================================================
  // Y-DNA Tab Content
  // ============================================================================

  private def createYdnaContent(): ScrollPane = {
    val terminalHaplogroupSection = new VBox(15) {
      padding = Insets(20)
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"
      children = Seq(
        new Label { text <== bind("haplogroup.terminal"); style = "-fx-font-weight: bold; -fx-text-fill: #aaaaaa;" },
        ydnaTerminalLabel,
        new VBox(5) {
          children = Seq(
            new Label { text <== bind("haplogroup.phylogenetic_path"); style = "-fx-font-size: 11px; -fx-text-fill: #888888;" },
            ydnaPathLabel
          )
        },
        new Separator { style = "-fx-background-color: #444444;" },
        new HBox(30) {
          alignment = Pos.CenterLeft
          children = Seq(
            new VBox(2) {
              children = Seq(
                new Label { text <== bind("haplogroup.derived"); style = "-fx-font-size: 11px; -fx-text-fill: #888888;" },
                ydnaDerivedLabel
              )
            },
            new VBox(2) {
              children = Seq(
                new Label { text <== bind("haplogroup.ancestral"); style = "-fx-font-size: 11px; -fx-text-fill: #888888;" },
                ydnaAncestralLabel
              )
            },
            new VBox(2) {
              children = Seq(
                new Label { text <== bind("haplogroup.confidence"); style = "-fx-font-size: 11px; -fx-text-fill: #888888;" },
                ydnaConfidenceLabel
              )
            },
            new VBox(2) {
              children = Seq(
                new Label { text <== bind("analysis.quality"); style = "-fx-font-size: 11px; -fx-text-fill: #888888;" },
                ydnaQualityLabel
              )
            }
          )
        },
        ydnaSourceLabel
      )
    }

    ydnaResultPane.children = Seq(terminalHaplogroupSection)

    val analyzeButton = new Button {
      text <== bind("analysis.run")
      styleClass += "button-primary"
      onAction = _ => handleRunYdnaAnalysis()
    }

    val viewProfileButton = new Button {
      text <== bind("haplogroup.view_profile")
      onAction = _ => log.debug("View full Y profile - not yet implemented")
    }

    new ScrollPane {
      fitToWidth = true
      style = "-fx-background: #1e1e1e; -fx-background-color: #1e1e1e;"
      content = new VBox(20) {
        padding = Insets(20)
        style = "-fx-background-color: #1e1e1e;"
        children = Seq(
          new HBox(10) {
            alignment = Pos.CenterLeft
            children = Seq(
              new Label { text <== bind("haplogroup.ydna.title"); style = "-fx-font-size: 18px; -fx-font-weight: bold; -fx-text-fill: #ffffff;" },
              new Region { hgrow = Priority.Always },
              viewProfileButton,
              analyzeButton
            )
          },
          new StackPane {
            children = Seq(ydnaNotAnalyzedPane, ydnaResultPane)
          }
        )
      }
    }
  }

  // ============================================================================
  // mtDNA Tab Content
  // ============================================================================

  private def createMtdnaContent(): ScrollPane = {
    val haplogroupSection = new VBox(15) {
      padding = Insets(20)
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"
      children = Seq(
        new Label { text <== bind("haplogroup.title"); style = "-fx-font-weight: bold; -fx-text-fill: #aaaaaa;" },
        mtdnaTerminalLabel,
        new VBox(5) {
          children = Seq(
            new Label { text <== bind("haplogroup.phylogenetic_path"); style = "-fx-font-size: 11px; -fx-text-fill: #888888;" },
            mtdnaPathLabel
          )
        },
        new Separator { style = "-fx-background-color: #444444;" },
        new HBox(30) {
          alignment = Pos.CenterLeft
          children = Seq(
            new VBox(2) {
              children = Seq(
                new Label { text <== bind("haplogroup.confidence"); style = "-fx-font-size: 11px; -fx-text-fill: #888888;" },
                mtdnaConfidenceLabel
              )
            },
            new VBox(2) {
              children = Seq(
                new Label { text <== bind("analysis.quality"); style = "-fx-font-size: 11px; -fx-text-fill: #888888;" },
                mtdnaQualityLabel
              )
            }
          )
        },
        mtdnaSourceLabel
      )
    }

    mtdnaResultPane.children = Seq(haplogroupSection)

    val analyzeButton = new Button {
      text <== bind("analysis.run")
      styleClass += "button-primary"
      onAction = _ => handleRunMtdnaAnalysis()
    }

    val viewDetailsButton = new Button {
      text <== bind("haplogroup.view_details")
      onAction = _ => log.debug("View mtDNA details - not yet implemented")
    }

    new ScrollPane {
      fitToWidth = true
      style = "-fx-background: #1e1e1e; -fx-background-color: #1e1e1e;"
      content = new VBox(20) {
        padding = Insets(20)
        style = "-fx-background-color: #1e1e1e;"
        children = Seq(
          new HBox(10) {
            alignment = Pos.CenterLeft
            children = Seq(
              new Label { text <== bind("haplogroup.mtdna.title"); style = "-fx-font-size: 18px; -fx-font-weight: bold; -fx-text-fill: #ffffff;" },
              new Region { hgrow = Priority.Always },
              viewDetailsButton,
              analyzeButton
            )
          },
          new StackPane {
            children = Seq(mtdnaNotAnalyzedPane, mtdnaResultPane)
          }
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
    val sequencingSection = createDataSection("data.sequencing_runs", sequencingListContainer, "data.no_sequencing")
    val chipSection = createDataSection("data.chip_profiles", chipListContainer, "data.no_chip")
    val strSection = createDataSection("data.str_profiles", strListContainer, "data.no_str")

    val addDataButton = new Button {
      text <== bind("data.add")
      styleClass += "button-primary"
      onAction = _ => handleAddData()
    }

    new ScrollPane {
      fitToWidth = true
      style = "-fx-background: #1e1e1e; -fx-background-color: #1e1e1e;"
      content = new VBox(20) {
        padding = Insets(20)
        style = "-fx-background-color: #1e1e1e;"
        children = Seq(
          new HBox(10) {
            alignment = Pos.CenterLeft
            children = Seq(
              new Label { text <== bind("data.title"); style = "-fx-font-size: 18px; -fx-font-weight: bold; -fx-text-fill: #ffffff;" },
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

  private def createDataSection(titleKey: String, container: VBox, emptyKey: String): VBox = {
    new VBox(10) {
      padding = Insets(15)
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"
      children = Seq(
        new Label { text <== bind(titleKey); style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" },
        container
      )
    }
  }

  // ============================================================================
  // Data Source Item Creation
  // ============================================================================

  private def createSequenceRunItem(seqRun: SequenceRun, index: Int): VBox = {
    val testTypeDisplay = SequenceRun.testTypeDisplayName(seqRun.testType)
    val readsDisplay = seqRun.totalReads.map(r => formatReadCount(r)).getOrElse("-")
    val alignedPct = seqRun.pctPfReadsAligned.map(p => f"${p * 100}%.1f%%").getOrElse("-")

    // Get alignments for this sequence run
    val alignments = viewModel.workspace.value.main.getAlignmentsForSequenceRun(seqRun)

    val runInfoBox = new HBox(15) {
      alignment = Pos.CenterLeft
      padding = Insets(10)
      style = "-fx-background-color: #333333; -fx-background-radius: 5 5 0 0;"
      children = Seq(
        // Icon/type indicator
        new Label {
          text = seqRun.platformName.take(3)
          prefWidth = 40
          style = "-fx-font-weight: bold; -fx-text-fill: #4ade80; -fx-font-family: monospace;"
        },
        // Main info
        new VBox(3) {
          hgrow = Priority.Always
          children = Seq(
            new HBox(8) {
              alignment = Pos.CenterLeft
              children = Seq(
                new Label(testTypeDisplay) { style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" },
                new Label(s"• ${seqRun.platformName}") { style = "-fx-text-fill: #888888;" },
                seqRun.instrumentModel.map(m => new Label(s"• $m") { style = "-fx-text-fill: #888888;" }).getOrElse(new Region)
              )
            },
            new HBox(15) {
              children = Seq(
                new Label(s"${t("data.reads")}: $readsDisplay") { style = "-fx-text-fill: #b0b0b0; -fx-font-size: 11px;" },
                new Label(s"${t("data.aligned")}: $alignedPct") { style = "-fx-text-fill: #b0b0b0; -fx-font-size: 11px;" },
                seqRun.libraryLayout.map(l => new Label(l) { style = "-fx-text-fill: #666666; -fx-font-size: 11px;" }).getOrElse(new Region)
              )
            }
          )
        },
        // Capabilities badges
        new HBox(5) {
          alignment = Pos.CenterRight
          children = {
            val badges = scala.collection.mutable.ArrayBuffer[scalafx.scene.Node]()
            if (SequenceRun.supportsYDna(seqRun.testType)) {
              badges += new Label("Y") {
                style = "-fx-background-color: #2d3a2d; -fx-text-fill: #4ade80; -fx-padding: 2 6; -fx-background-radius: 3; -fx-font-size: 10px;"
              }
            }
            if (SequenceRun.supportsMtDna(seqRun.testType)) {
              badges += new Label("mt") {
                style = "-fx-background-color: #2d2d3a; -fx-text-fill: #60a5fa; -fx-padding: 2 6; -fx-background-radius: 3; -fx-font-size: 10px;"
              }
            }
            badges.toSeq
          }
        },
        // Action menu button
        new MenuButton("⋮") {
          style = "-fx-background-color: transparent; -fx-text-fill: #888888; -fx-font-size: 16px;"
          items = createSequenceRunMenuItems(seqRun, index)
        }
      )
    }

    // Alignments section (collapsed under the run info)
    val alignmentsContainer = new VBox(5) {
      padding = Insets(5, 10, 10, 50) // Indented from the left
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 0 0 5 5;"
      visible = alignments.nonEmpty
      managed = alignments.nonEmpty
      children = alignments.zipWithIndex.map { case (alignment, alignIdx) =>
        createAlignmentItem(alignment, index, alignIdx)
      }
    }

    new VBox(0) {
      style = "-fx-background-color: #333333; -fx-background-radius: 5;"
      children = if (alignments.nonEmpty) {
        Seq(runInfoBox, alignmentsContainer)
      } else {
        // Adjust runInfoBox style when no alignments
        runInfoBox.style = "-fx-background-color: #333333; -fx-background-radius: 5;"
        Seq(runInfoBox)
      }
    }
  }

  /** Creates a display item for an alignment with analysis actions */
  private def createAlignmentItem(alignmentData: Alignment, seqRunIndex: Int, alignIndex: Int): HBox = {
    val coverageDisplay = alignmentData.metrics.flatMap(_.meanCoverage).map(c => f"$c%.1fx").getOrElse("-")
    val callableDisplay = alignmentData.metrics.flatMap(_.callableBases).map(b => Formatters.formatNumber(b)).getOrElse("-")

    val refBuild = alignmentData.referenceBuild
    val alignerName = alignmentData.aligner
    val hasCoverage = alignmentData.metrics.flatMap(_.meanCoverage).isDefined
    val hasCallable = alignmentData.metrics.flatMap(_.callableBases).isDefined
    val hasVcf = alignmentData.metrics.exists(_.hasVcf)

    new HBox(10) {
      alignment = Pos.CenterLeft
      padding = Insets(8)
      style = "-fx-background-color: #3a3a3a; -fx-background-radius: 3;"
      children = Seq(
        // Reference badge
        new Label(refBuild) {
          prefWidth = 80
          style = "-fx-font-weight: bold; -fx-text-fill: #60a5fa; -fx-font-size: 11px;"
        },
        // Aligner
        new Label(alignerName) {
          style = "-fx-text-fill: #888888; -fx-font-size: 11px;"
        },
        // Metrics info
        new HBox(15) {
          hgrow = Priority.Always
          children = Seq(
            new Label(s"Coverage: $coverageDisplay") {
              style = s"-fx-text-fill: ${if (coverageDisplay == "-") "#666666" else "#b0b0b0"}; -fx-font-size: 11px;"
            },
            new Label(s"Callable: $callableDisplay") {
              style = s"-fx-text-fill: ${if (callableDisplay == "-") "#666666" else "#b0b0b0"}; -fx-font-size: 11px;"
            }
          )
        },
        // Status indicators
        new HBox(5) {
          children = {
            val indicators = scala.collection.mutable.ArrayBuffer[scalafx.scene.Node]()
            if (hasCoverage) {
              indicators += new Label("WGS") {
                style = "-fx-background-color: #2d3a2d; -fx-text-fill: #4ade80; -fx-padding: 1 4; -fx-background-radius: 2; -fx-font-size: 9px;"
              }
            }
            if (hasCallable) {
              indicators += new Label("CL") {
                style = "-fx-background-color: #2d3a2d; -fx-text-fill: #4ade80; -fx-padding: 1 4; -fx-background-radius: 2; -fx-font-size: 9px;"
              }
            }
            if (hasVcf) {
              indicators += new Label("VCF") {
                style = "-fx-background-color: #3a2d3a; -fx-text-fill: #c084fc; -fx-padding: 1 4; -fx-background-radius: 2; -fx-font-size: 9px;"
              }
            }
            indicators.toSeq
          }
        },
        // Action menu
        new MenuButton("⋮") {
          style = "-fx-background-color: transparent; -fx-text-fill: #666666; -fx-font-size: 14px;"
          items = createAlignmentMenuItems(alignmentData, seqRunIndex, alignIndex)
        }
      )
    }
  }

  /** Creates context menu items for alignment actions */
  private def createAlignmentMenuItems(alignment: Alignment, seqRunIndex: Int, alignIndex: Int): Seq[MenuItem] = {
    val items = scala.collection.mutable.ArrayBuffer[MenuItem]()

    // Run WGS Metrics
    items += new MenuItem("Run WGS Metrics") {
      disable = alignment.metrics.flatMap(_.meanCoverage).isDefined
      onAction = _ => handleRunWgsMetrics(seqRunIndex, alignIndex)
    }

    // Run Callable Loci
    items += new MenuItem("Run Callable Loci") {
      disable = alignment.metrics.flatMap(_.callableBases).isDefined
      onAction = _ => handleRunCallableLoci(seqRunIndex, alignIndex)
    }

    items += new SeparatorMenuItem()

    // Details
    items += new MenuItem("Details") {
      onAction = _ => showAlignmentDetailsDialog(alignment)
    }

    items.toSeq
  }

  /** Handle running WGS metrics analysis for a specific alignment */
  private def handleRunWgsMetrics(seqRunIndex: Int, alignIndex: Int): Unit = {
    currentSubject.value.foreach { subject =>
      val progressDialog = new AnalysisProgressDialog(
        "WGS Metrics Analysis",
        viewModel.analysisProgress,
        viewModel.analysisProgressPercent,
        viewModel.analysisInProgress
      )
      Option(SubjectDetailView.this.getScene).flatMap(s => Option(s.getWindow)).foreach { window =>
        progressDialog.initOwner(window)
      }

      viewModel.runWgsMetricsAnalysisForAlignment(
        subject.accession,
        seqRunIndex,
        alignIndex,
        {
          case Right(metrics) =>
            scalafx.application.Platform.runLater {
              updateDataSources(subject)
              showInfoDialog(
                "WGS Metrics Complete",
                "Analysis finished successfully",
                s"Mean coverage: ${f"${metrics.meanCoverage}%.1f"}x\nMedian coverage: ${f"${metrics.medianCoverage}%.1f"}x"
              )
            }
          case Left(error) =>
            scalafx.application.Platform.runLater {
              showInfoDialog(t("error.title"), "WGS Metrics Failed", error)
            }
        }
      )

      progressDialog.show()
    }
  }

  /** Handle running callable loci analysis for a specific alignment */
  private def handleRunCallableLoci(seqRunIndex: Int, alignIndex: Int): Unit = {
    currentSubject.value.foreach { subject =>
      val progressDialog = new AnalysisProgressDialog(
        "Callable Loci Analysis",
        viewModel.analysisProgress,
        viewModel.analysisProgressPercent,
        viewModel.analysisInProgress
      )
      Option(SubjectDetailView.this.getScene).flatMap(s => Option(s.getWindow)).foreach { window =>
        progressDialog.initOwner(window)
      }

      viewModel.runCallableLociAnalysisForAlignment(
        subject.accession,
        seqRunIndex,
        alignIndex,
        {
          case Right((result, _)) =>
            scalafx.application.Platform.runLater {
              updateDataSources(subject)
              showInfoDialog(
                "Callable Loci Complete",
                "Analysis finished successfully",
                s"Callable bases: ${Formatters.formatNumber(result.callableBases)}\nContigs analyzed: ${result.contigAnalysis.size}"
              )
            }
          case Left(error) =>
            scalafx.application.Platform.runLater {
              showInfoDialog(t("error.title"), "Callable Loci Failed", error)
            }
        }
      )

      progressDialog.show()
    }
  }

  /** Show alignment details dialog */
  private def showAlignmentDetailsDialog(alignment: Alignment): Unit = {
    val metrics = alignment.metrics.getOrElse(com.decodingus.workspace.model.AlignmentMetrics())
    val detailsText =
      s"""Reference: ${alignment.referenceBuild}
         |Aligner: ${alignment.aligner}
         |Variant Caller: ${alignment.variantCaller.getOrElse("N/A")}
         |Files: ${alignment.files.size}
         |
         |--- WGS Metrics ---
         |Mean Coverage: ${metrics.meanCoverage.map(c => f"$c%.2fx").getOrElse("Not analyzed")}
         |Median Coverage: ${metrics.medianCoverage.map(c => f"$c%.2fx").getOrElse("N/A")}
         |SD Coverage: ${metrics.sdCoverage.map(c => f"$c%.2f").getOrElse("N/A")}
         |% at 10x: ${metrics.pct10x.map(p => f"${p * 100}%.1f%%").getOrElse("N/A")}
         |% at 20x: ${metrics.pct20x.map(p => f"${p * 100}%.1f%%").getOrElse("N/A")}
         |% at 30x: ${metrics.pct30x.map(p => f"${p * 100}%.1f%%").getOrElse("N/A")}
         |
         |--- Callable Loci ---
         |Callable Bases: ${metrics.callableBases.map(b => Formatters.formatNumber(b)).getOrElse("Not analyzed")}
         |Analysis Complete: ${metrics.callableLociComplete.map(_.toString).getOrElse("N/A")}
         |
         |--- VCF Status ---
         |VCF Generated: ${if (metrics.hasVcf) "Yes" else "No"}
         |Variant Count: ${metrics.vcfVariantCount.map(v => Formatters.formatNumber(v)).getOrElse("N/A")}
         |
         |--- Sex Inference ---
         |Inferred Sex: ${metrics.inferredSex.getOrElse("Not determined")}
         |Confidence: ${metrics.sexInferenceConfidence.getOrElse("N/A")}
       """.stripMargin

    InfoDialog.showCode(
      "Alignment Details",
      s"${alignment.referenceBuild} - ${alignment.aligner}",
      detailsText,
      dialogWidth = 450,
      dialogHeight = 500
    )
  }

  /** Creates context menu items for sequence run actions */
  private def createSequenceRunMenuItems(seqRun: SequenceRun, index: Int): Seq[MenuItem] = {
    import com.decodingus.ui.components.{MergeSequenceRunsDialog, MergeDecision}

    val items = scala.collection.mutable.ArrayBuffer[MenuItem]()

    // Details
    items += new MenuItem("Details") {
      onAction = _ => showSequenceRunDetailsDialog(seqRun)
    }

    items += new SeparatorMenuItem()

    // Merge with another run
    items += new MenuItem("Merge with another run...") {
      onAction = _ => {
        currentSubject.value.foreach { subject =>
          val allRuns = viewModel.workspace.value.main.getSequenceRunsForBiosample(subject)
          if (allRuns.size < 2) {
            showInfoDialog("Cannot Merge", "Not enough runs", "You need at least two sequence runs to merge.")
          } else {
            // Find another run to merge with (pick the first one that isn't this one)
            val otherRuns = allRuns.zipWithIndex.filterNot(_._2 == index)
            if (otherRuns.isEmpty) {
              showInfoDialog("Cannot Merge", "Not enough runs", "No other sequence runs available to merge with.")
            } else {
              // For simplicity, show merge dialog with first other run
              // Could enhance to let user pick which run to merge with
              val (otherRun, otherIndex) = otherRuns.head
              val mergeDialog = new MergeSequenceRunsDialog(seqRun, index, otherRun, otherIndex)
              // Set owner window from this view if available
              Option(SubjectDetailView.this.getScene).flatMap(s => Option(s.getWindow)).foreach { window =>
                mergeDialog.initOwner(window)
              }

              mergeDialog.showAndWait() match {
                case Some(Some(decision: MergeDecision)) =>
                  viewModel.mergeSequenceRuns(
                    subject.sampleAccession,
                    decision.primaryIndex,
                    decision.secondaryIndex
                  ) match {
                    case Right(movedCount) =>
                      updateDataSources(subject)
                      showInfoDialog(
                        "Merge Complete",
                        "Sequence runs merged",
                        s"Moved $movedCount alignment${if (movedCount != 1) "s" else ""} to the primary run."
                      )
                    case Left(error) =>
                      showInfoDialog(t("error.title"), "Merge Failed", error)
                  }
                case _ =>
                  log.debug("Merge dialog cancelled")
              }
            }
          }
        }
      }
    }

    items += new SeparatorMenuItem()

    // Remove
    items += new MenuItem("Remove") {
      onAction = _ => {
        val details = s"${seqRun.platformName} - ${SequenceRun.testTypeDisplayName(seqRun.testType)}"
        if (ConfirmDialog.confirmRemoval("Sequence Run", details)) {
          currentSubject.value.foreach { subject =>
            viewModel.removeSequenceData(subject.accession, index)
            updateDataSources(subject)
          }
        }
      }
    }

    items.toSeq
  }

  /** Shows sequence run details dialog */
  private def showSequenceRunDetailsDialog(seqRun: SequenceRun): Unit = {
    val alignmentCount = seqRun.alignmentRefs.size
    val fileCount = seqRun.files.size

    val detailsText =
      s"""Platform: ${seqRun.platformName}
         |Instrument: ${seqRun.instrumentModel.getOrElse("Unknown")}
         |Test Type: ${SequenceRun.testTypeDisplayName(seqRun.testType)}
         |Library Layout: ${seqRun.libraryLayout.getOrElse("Unknown")}
         |
         |Total Reads: ${seqRun.totalReads.map(r => Formatters.formatNumber(r)).getOrElse("N/A")}
         |PF Reads Aligned: ${seqRun.pctPfReadsAligned.map(p => f"${p * 100}%.1f%%").getOrElse("N/A")}
         |Read Length: ${seqRun.readLength.map(_.toString).getOrElse("N/A")}
         |Mean Insert Size: ${seqRun.meanInsertSize.map(s => f"$s%.0f").getOrElse("N/A")}
         |
         |Sample Name: ${seqRun.sampleName.getOrElse("N/A")}
         |Library ID: ${seqRun.libraryId.getOrElse("N/A")}
         |Platform Unit: ${seqRun.platformUnit.getOrElse("N/A")}
         |Fingerprint: ${seqRun.runFingerprint.map(_.take(16) + "...").getOrElse("Not computed")}
         |
         |Alignments: $alignmentCount
         |Files: $fileCount
       """.stripMargin

    InfoDialog.showCode(
      "Sequence Run Details",
      s"${seqRun.platformName} - ${SequenceRun.testTypeDisplayName(seqRun.testType)}",
      detailsText,
      dialogWidth = 420,
      dialogHeight = 400
    )
  }

  private def createChipProfileItem(chip: ChipProfile, index: Int): HBox = {
    val callRatePct = f"${chip.callRate * 100}%.1f%%"
    val statusStyle = chip.status match {
      case "Good" => "-fx-text-fill: #4ade80;"
      case "Acceptable" => "-fx-text-fill: #fbbf24;"
      case _ => "-fx-text-fill: #f87171;"
    }

    val itemBox = new HBox(15) {
      alignment = Pos.CenterLeft
      padding = Insets(10)
      style = "-fx-background-color: #333333; -fx-background-radius: 5;"
      children = Seq(
        // Icon/vendor indicator
        new Label {
          text = chip.vendor.take(3).toUpperCase
          prefWidth = 40
          style = "-fx-font-weight: bold; -fx-text-fill: #fbbf24; -fx-font-family: monospace;"
        },
        // Main info
        new VBox(3) {
          hgrow = Priority.Always
          children = Seq(
            new HBox(8) {
              alignment = Pos.CenterLeft
              children = Seq(
                new Label(chip.vendor) { style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" },
                chip.chipVersion.map(v => new Label(s"• $v") { style = "-fx-text-fill: #888888;" }).getOrElse(new Region)
              )
            },
            new HBox(15) {
              children = Seq(
                new Label(s"${t("data.markers")}: ${Formatters.formatNumber(chip.totalMarkersCalled)}") { style = "-fx-text-fill: #b0b0b0; -fx-font-size: 11px;" },
                new Label(s"${t("data.call_rate")}: $callRatePct") { style = "-fx-text-fill: #b0b0b0; -fx-font-size: 11px;" },
                new Label(chip.status) { style = s"-fx-font-size: 11px; $statusStyle" }
              )
            }
          )
        },
        // Capabilities badges
        new HBox(5) {
          alignment = Pos.CenterRight
          children = {
            val badges = scala.collection.mutable.ArrayBuffer[scalafx.scene.Node]()
            if (chip.hasSufficientYCoverage) {
              badges += new Label("Y") {
                style = "-fx-background-color: #2d3a2d; -fx-text-fill: #4ade80; -fx-padding: 2 6; -fx-background-radius: 3; -fx-font-size: 10px;"
              }
            }
            if (chip.hasSufficientMtCoverage) {
              badges += new Label("mt") {
                style = "-fx-background-color: #2d2d3a; -fx-text-fill: #60a5fa; -fx-padding: 2 6; -fx-background-radius: 3; -fx-font-size: 10px;"
              }
            }
            if (chip.isAcceptableForAncestry) {
              badges += new Label("Anc") {
                style = "-fx-background-color: #3a2d3a; -fx-text-fill: #c084fc; -fx-padding: 2 6; -fx-background-radius: 3; -fx-font-size: 10px;"
              }
            }
            badges.toSeq
          }
        },
        // Action menu button
        new MenuButton("⋮") {
          style = "-fx-background-color: transparent; -fx-text-fill: #888888; -fx-font-size: 16px;"
          items = createChipProfileMenuItems(chip)
        }
      )
    }
    itemBox
  }

  /** Creates context menu items for chip profile actions */
  private def createChipProfileMenuItems(chip: ChipProfile): Seq[MenuItem] = {
    val items = scala.collection.mutable.ArrayBuffer[MenuItem]()

    // Details
    items += new MenuItem("Details") {
      onAction = _ => showChipDetailsDialog(chip)
    }

    items += new SeparatorMenuItem()

    // Y-DNA Haplogroup
    val yMenuItem = new MenuItem("Y-DNA Haplogroup") {
      onAction = _ => handleChipHaplogroupAnalysis(chip, TreeType.YDNA)
    }
    if (!chip.hasSufficientYCoverage) {
      yMenuItem.disable = true
    }
    items += yMenuItem

    // mtDNA Haplogroup
    val mtMenuItem = new MenuItem("mtDNA Haplogroup") {
      onAction = _ => handleChipHaplogroupAnalysis(chip, TreeType.MTDNA)
    }
    if (!chip.hasSufficientMtCoverage) {
      mtMenuItem.disable = true
    }
    items += mtMenuItem

    // Ancestry Analysis
    val ancMenuItem = new MenuItem("Ancestry Analysis") {
      onAction = _ => handleChipAncestryAnalysis(chip)
    }
    if (!chip.isAcceptableForAncestry) {
      ancMenuItem.disable = true
    }
    items += ancMenuItem

    items += new SeparatorMenuItem()

    // Remove
    items += new MenuItem("Remove") {
      onAction = _ => {
        val details = s"${chip.vendor} - ${Formatters.formatNumber(chip.totalMarkersCalled)} markers"
        if (ConfirmDialog.confirmRemoval("Chip Data", details)) {
          chip.atUri.foreach { uri =>
            currentSubject.value.foreach { subject =>
              viewModel.deleteChipProfile(subject.accession, uri)
              updateDataSources(subject)
            }
          }
        }
      }
    }

    items.toSeq
  }

  /** Shows chip profile details dialog */
  private def showChipDetailsDialog(chip: ChipProfile): Unit = {
    val detailsText =
      s"""Vendor: ${chip.vendor}
         |Test Type: ${chip.testTypeCode}
         |Chip Version: ${chip.chipVersion.getOrElse("Unknown")}
         |
         |Total Markers: ${Formatters.formatNumber(chip.totalMarkersCalled)} / ${Formatters.formatNumber(chip.totalMarkersPossible)}
         |Call Rate: ${f"${chip.callRate * 100}%.2f"}%
         |No-Call Rate: ${f"${chip.noCallRate * 100}%.2f"}%
         |
         |Autosomal Markers: ${Formatters.formatNumber(chip.autosomalMarkersCalled)}
         |Y-DNA Markers: ${chip.yMarkersCalled.map(n => Formatters.formatNumber(n)).getOrElse("N/A")}
         |mtDNA Markers: ${chip.mtMarkersCalled.map(n => Formatters.formatNumber(n)).getOrElse("N/A")}
         |Heterozygosity Rate: ${chip.hetRate.map(r => f"${r * 100}%.2f%%").getOrElse("N/A")}
         |
         |Status: ${chip.status}
         |Suitable for Ancestry: ${if (chip.isAcceptableForAncestry) "Yes" else "No"}
         |Sufficient Y Coverage: ${if (chip.hasSufficientYCoverage) "Yes" else "No"}
         |Sufficient MT Coverage: ${if (chip.hasSufficientMtCoverage) "Yes" else "No"}
         |
         |Import Date: ${chip.importDate.toLocalDate}
         |Source File: ${chip.sourceFileName.getOrElse("Unknown")}
       """.stripMargin

    InfoDialog.showCode(
      "Chip Data Details",
      s"${chip.vendor} - ${chip.testTypeCode}",
      detailsText,
      dialogWidth = 420,
      dialogHeight = 380
    )
  }

  /** Handles Y-DNA or mtDNA haplogroup analysis from chip data */
  private def handleChipHaplogroupAnalysis(chip: ChipProfile, treeType: TreeType): Unit = {
    val typeName = if (treeType == TreeType.YDNA) "Y-DNA" else "mtDNA"
    val markerCount = treeType match {
      case TreeType.YDNA => chip.yMarkersCalled.getOrElse(0)
      case TreeType.MTDNA => chip.mtMarkersCalled.getOrElse(0)
    }

    val confirm = new Alert(AlertType.Confirmation) {
      title = s"Run $typeName Haplogroup Analysis"
      headerText = s"Analyze ${chip.vendor} chip data for $typeName haplogroup"
      contentText =
        s"""This will score chip genotypes against the $typeName haplogroup tree.

$typeName Markers: ${Formatters.formatNumber(markerCount)}

Note: Chip-based haplogroup estimation has limited resolution compared to WGS.
The terminal haplogroup may be upstream of the true assignment."""
    }

    confirm.showAndWait() match {
      case Some(ButtonType.OK) =>
        chip.atUri match {
          case Some(profileUri) =>
            currentSubject.value.foreach { subject =>
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
                    scalafx.application.Platform.runLater {
                      import com.decodingus.genotype.processor.ChipHaplogroupAdapter
                      val confidenceDesc = ChipHaplogroupAdapter.confidenceDescription(haplogroupResult.confidence)
                      new Alert(AlertType.Information) {
                        title = s"$typeName Haplogroup Result"
                        headerText = s"$typeName: ${haplogroupResult.topHaplogroup}"
                        contentText =
                          s"""Confidence: $confidenceDesc (${f"${haplogroupResult.confidence * 100}%.0f"}%)
SNPs Matched: ${haplogroupResult.snpsMatched} / ${haplogroupResult.snpsTotal}
Tree Depth: ${haplogroupResult.results.headOption.map(_.depth).getOrElse(0)}

Note: Chip data covers ~${f"${haplogroupResult.snpsMatched.toDouble / haplogroupResult.snpsTotal * 100}%.0f"}% of tree positions.
For higher resolution, consider WGS analysis."""
                      }.showAndWait()
                    }
                  case Left(error) =>
                    scalafx.application.Platform.runLater {
                      new Alert(AlertType.Error) {
                        title = s"$typeName Haplogroup Analysis Failed"
                        headerText = "Could not complete haplogroup analysis"
                        contentText = error
                      }.showAndWait()
                    }
                }
              )

              progressDialog.show()
            }
          case None =>
            showInfoDialog("Error", "Invalid chip profile", "Profile has no AT URI.")
        }
      case _ => // User cancelled
    }
  }

  /** Handles ancestry analysis from chip data */
  private def handleChipAncestryAnalysis(chip: ChipProfile): Unit = {
    import com.decodingus.ancestry.model.AncestryPanelType

    val recommendedPanel = if (chip.autosomalMarkersCalled >= 500000) {
      AncestryPanelType.GenomeWide
    } else {
      AncestryPanelType.Aims
    }

    val panelLabel = recommendedPanel match {
      case AncestryPanelType.Aims => "AIMs (~5k markers, faster)"
      case AncestryPanelType.GenomeWide => "Genome-wide (~500k markers, detailed)"
    }

    val confirm = new Alert(AlertType.Confirmation) {
      title = "Run Ancestry Analysis"
      headerText = s"Analyze ${chip.vendor} chip data for ancestry"
      contentText =
        s"""This will estimate population percentages using the $panelLabel panel.

Markers: ${Formatters.formatNumber(chip.autosomalMarkersCalled)}
Call Rate: ${f"${chip.callRate * 100}%.1f"}%

Note: Reference data download may be required on first run."""
    }

    confirm.showAndWait() match {
      case Some(ButtonType.OK) =>
        chip.atUri match {
          case Some(profileUri) =>
            currentSubject.value.foreach { subject =>
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
                    scalafx.application.Platform.runLater {
                      val resultDialog = new AncestryResultDialog(ancestryResult)
                      resultDialog.showAndWait()
                    }
                  case Left(error) =>
                    scalafx.application.Platform.runLater {
                      new Alert(AlertType.Error) {
                        title = "Ancestry Analysis Failed"
                        headerText = "Could not complete ancestry analysis"
                        contentText = error
                      }.showAndWait()
                    }
                }
              )

              progressDialog.show()
            }
          case None =>
            showInfoDialog("Error", "Invalid chip profile", "Profile has no AT URI.")
        }
      case _ => // User cancelled
    }
  }

  private def createStrProfileItem(strProfile: StrProfile, index: Int): HBox = {
    val markerCount = strProfile.totalMarkers.getOrElse(strProfile.markers.size)
    val panelNames = strProfile.panels.map(_.panelName).mkString(", ")
    val sourceDisplay = strProfile.source.getOrElse(strProfile.importedFrom.getOrElse("-"))

    new HBox(15) {
      alignment = Pos.CenterLeft
      padding = Insets(10)
      style = "-fx-background-color: #333333; -fx-background-radius: 5;"
      children = Seq(
        // Icon/source indicator
        new Label {
          text = strProfile.importedFrom.map(_.take(3).toUpperCase).getOrElse("STR")
          prefWidth = 40
          style = "-fx-font-weight: bold; -fx-text-fill: #c084fc; -fx-font-family: monospace;"
        },
        // Main info
        new VBox(3) {
          hgrow = Priority.Always
          children = Seq(
            new HBox(8) {
              alignment = Pos.CenterLeft
              children = Seq(
                new Label(if (panelNames.nonEmpty) panelNames else t("data.str_profile")) { style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;" },
                strProfile.importedFrom.map(s => new Label(s"• $s") { style = "-fx-text-fill: #888888;" }).getOrElse(new Region)
              )
            },
            new HBox(15) {
              children = Seq(
                new Label(s"${t("data.markers")}: $markerCount") { style = "-fx-text-fill: #b0b0b0; -fx-font-size: 11px;" },
                new Label(s"${t("data.source")}: $sourceDisplay") { style = "-fx-text-fill: #b0b0b0; -fx-font-size: 11px;" }
              )
            }
          )
        }
      )
    }
  }

  private def formatReadCount(reads: Long): String = {
    if (reads >= 1_000_000_000) f"${reads / 1_000_000_000.0}%.1fB"
    else if (reads >= 1_000_000) f"${reads / 1_000_000.0}%.1fM"
    else if (reads >= 1_000) f"${reads / 1_000.0}%.1fK"
    else reads.toString
  }

  private def createEmptyPlaceholder(messageKey: String): Label = {
    new Label {
      text <== bind(messageKey)
      style = "-fx-text-fill: #666666; -fx-font-style: italic;"
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
    overviewYdnaHaplogroupLabel.text = subject.yHaplogroup.getOrElse("-")
    subject.yHaplogroupResult match {
      case Some(result) =>
        overviewYdnaConfidenceLabel.text = s"${t("haplogroup.confidence")}: ${t(s"haplogroup.confidence.${result.confidenceLevel.toLowerCase}")}"
        overviewYdnaConfidenceLabel.style = confidenceStyle(result.confidenceLevel)
      case None =>
        overviewYdnaConfidenceLabel.text = ""
    }

    // Update Overview tab - mtDNA card
    overviewMtdnaHaplogroupLabel.text = subject.mtHaplogroup.getOrElse("-")
    subject.mtHaplogroupResult match {
      case Some(result) =>
        overviewMtdnaConfidenceLabel.text = s"${t("haplogroup.confidence")}: ${t(s"haplogroup.confidence.${result.confidenceLevel.toLowerCase}")}"
        overviewMtdnaConfidenceLabel.style = confidenceStyle(result.confidenceLevel)
      case None =>
        overviewMtdnaConfidenceLabel.text = ""
    }

    // Update Y-DNA tab
    subject.yHaplogroupResult match {
      case Some(result) =>
        ydnaTerminalLabel.text = result.haplogroupName
        ydnaPathLabel.text = result.formattedPath
        ydnaDerivedLabel.text = result.derivedCount.toString
        ydnaAncestralLabel.text = result.ancestralCount.toString
        ydnaConfidenceLabel.text = s"${result.confidencePercent} (${t(s"haplogroup.confidence.${result.confidenceLevel.toLowerCase}")})"
        ydnaConfidenceLabel.style = s"-fx-font-weight: bold; ${confidenceTextColor(result.confidenceLevel)}"
        ydnaSourceLabel.text = s"${t("data.platform")}: ${result.sourceDisplay}"
        ydnaQualityLabel.text = t(s"analysis.quality.${result.qualityRating.toLowerCase}")
        ydnaQualityLabel.style = s"-fx-font-weight: bold; ${qualityTextColor(result.qualityRating)}"
        ydnaNotAnalyzedPane.visible = false
        ydnaResultPane.visible = true
      case None =>
        ydnaTerminalLabel.text = "-"
        ydnaPathLabel.text = ""
        ydnaDerivedLabel.text = "-"
        ydnaAncestralLabel.text = "-"
        ydnaConfidenceLabel.text = "-"
        ydnaSourceLabel.text = ""
        ydnaQualityLabel.text = "-"
        ydnaNotAnalyzedPane.visible = true
        ydnaResultPane.visible = false
    }

    // Update mtDNA tab
    subject.mtHaplogroupResult match {
      case Some(result) =>
        mtdnaTerminalLabel.text = result.haplogroupName
        mtdnaPathLabel.text = result.formattedPath
        mtdnaConfidenceLabel.text = s"${result.confidencePercent} (${t(s"haplogroup.confidence.${result.confidenceLevel.toLowerCase}")})"
        mtdnaConfidenceLabel.style = s"-fx-font-weight: bold; ${confidenceTextColor(result.confidenceLevel)}"
        mtdnaSourceLabel.text = s"${t("data.platform")}: ${result.sourceDisplay}"
        mtdnaQualityLabel.text = t(s"analysis.quality.${result.qualityRating.toLowerCase}")
        mtdnaQualityLabel.style = s"-fx-font-weight: bold; ${qualityTextColor(result.qualityRating)}"
        mtdnaNotAnalyzedPane.visible = false
        mtdnaResultPane.visible = true
      case None =>
        mtdnaTerminalLabel.text = "-"
        mtdnaPathLabel.text = ""
        mtdnaConfidenceLabel.text = "-"
        mtdnaSourceLabel.text = ""
        mtdnaQualityLabel.text = "-"
        mtdnaNotAnalyzedPane.visible = true
        mtdnaResultPane.visible = false
    }

    // Update data counts and Data Sources tab
    updateDataSources(subject)
  }

  private def updateDataSources(subject: Biosample): Unit = {
    // Get actual data from ViewModel
    val sequenceRuns = viewModel.workspace.value.main.getSequenceRunsForBiosample(subject)
    val chipProfiles = viewModel.getChipProfilesForBiosample(subject.accession)
    val strProfiles = viewModel.getStrProfilesForBiosample(subject.accession)

    // Update counts on Overview tab
    sequencingCountLabel.text = sequenceRuns.size.toString
    chipCountLabel.text = chipProfiles.size.toString
    strCountLabel.text = strProfiles.size.toString

    // Update sequencing runs container
    sequencingListContainer.children.clear()
    if (sequenceRuns.isEmpty) {
      sequencingListContainer.children += createEmptyPlaceholder("data.no_sequencing")
    } else {
      sequenceRuns.zipWithIndex.foreach { case (seqRun, idx) =>
        sequencingListContainer.children += createSequenceRunItem(seqRun, idx)
      }
    }

    // Update chip profiles container
    chipListContainer.children.clear()
    if (chipProfiles.isEmpty) {
      chipListContainer.children += createEmptyPlaceholder("data.no_chip")
    } else {
      chipProfiles.zipWithIndex.foreach { case (chip, idx) =>
        chipListContainer.children += createChipProfileItem(chip, idx)
      }
    }

    // Update STR profiles container
    strListContainer.children.clear()
    if (strProfiles.isEmpty) {
      strListContainer.children += createEmptyPlaceholder("data.no_str")
    } else {
      strProfiles.zipWithIndex.foreach { case (strProfile, idx) =>
        strListContainer.children += createStrProfileItem(strProfile, idx)
      }
    }
  }

  private def confidenceStyle(level: String): String = level match {
    case "HIGH" => "-fx-font-size: 11px; -fx-text-fill: #4ade80;"
    case "MEDIUM" => "-fx-font-size: 11px; -fx-text-fill: #fbbf24;"
    case _ => "-fx-font-size: 11px; -fx-text-fill: #f87171;"
  }

  private def confidenceTextColor(level: String): String = level match {
    case "HIGH" => "-fx-text-fill: #4ade80;"
    case "MEDIUM" => "-fx-text-fill: #fbbf24;"
    case _ => "-fx-text-fill: #f87171;"
  }

  private def qualityTextColor(quality: String): String = quality match {
    case "Excellent" => "-fx-text-fill: #4ade80;"
    case "Good" => "-fx-text-fill: #60a5fa;"
    case "Fair" => "-fx-text-fill: #fbbf24;"
    case _ => "-fx-text-fill: #f87171;"
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
    currentSubject.value match {
      case None =>
        log.warn("handleAddData called but no subject selected")
        showInfoDialog(
          t("error.title"),
          t("error.no_subject"),
          t("error.select_subject_first")
        )
      case Some(subject) =>
        try {
          // Get existing checksums for duplicate detection
          val existingChecksums = viewModel.getExistingChecksums(subject.accession)

          val dialog = new AddDataDialog(existingChecksums)
          // Set owner window for proper modal behavior
          Option(this.getScene).flatMap(s => Option(s.getWindow)).foreach { window =>
            dialog.initOwner(window)
          }

          dialog.showAndWait() match {
            case Some(Some(dataInput: DataInput)) =>
              log.info(s"Adding ${dataInput.dataType} data for ${subject.accession}: ${dataInput.fileInfo.fileName}")

              // Handle based on data type
              val file = new java.io.File(dataInput.fileInfo.location.getOrElse(""))

              dataInput.dataType match {
                case DataType.Alignment =>
                  // Use addFileAndAnalyze for BAM/CRAM to get proper fingerprint matching
                  // This detects if the file is a re-alignment of existing data to a different reference
                  val progressDialog = new AnalysisProgressDialog(
                    "Analyzing Alignment",
                    viewModel.analysisProgress,
                    viewModel.analysisProgressPercent,
                    viewModel.analysisInProgress
                  )

                  viewModel.addFileAndAnalyze(
                    subject.accession,
                    dataInput.fileInfo,
                    onProgress = (message, percent) => {
                      scalafx.application.Platform.runLater {
                        viewModel.analysisProgress.value = message
                        viewModel.analysisProgressPercent.value = percent
                      }
                    },
                    onComplete = {
                      case Right((index, libraryStats)) =>
                        scalafx.application.Platform.runLater {
                          log.info(s"Added alignment at index $index for ${subject.accession} (ref: ${libraryStats.referenceBuild})")
                          updateDataSources(subject)
                          showInfoDialog(
                            t("data.title"),
                            t("data.added.success"),
                            s"${dataInput.fileInfo.fileName}\n${libraryStats.inferredPlatform} - ${libraryStats.referenceBuild}"
                          )
                        }
                      case Left(error) =>
                        scalafx.application.Platform.runLater {
                          log.error(s"Failed to add alignment: $error")
                          showInfoDialog(t("error.title"), t("data.alignment"), error)
                        }
                    }
                  )
                  progressDialog.show()
                  // Alignment analysis is async, return early
                  return

                case DataType.Variants =>
                  // VCF files - show metadata dialog to get test type info
                  val metadataDialog = new VcfMetadataDialog(dataInput.fileInfo)
                  Option(this.getScene).flatMap(s => Option(s.getWindow)).foreach { window =>
                    metadataDialog.initOwner(window)
                  }

                  metadataDialog.showAndWait() match {
                    case Some(Some(metadata: VcfMetadata)) =>
                      // Create a SequenceRun with the VCF file using existing method
                      val vcfFileInfo = dataInput.fileInfo.copy(
                        fileFormat = if (dataInput.fileInfo.fileName.toLowerCase.endsWith(".vcf.gz")) "VCF.GZ" else "VCF"
                      )

                      // Add the file first (creates sequence run with placeholder values)
                      val newIndex = viewModel.addSequenceRunFromFile(subject.accession, vcfFileInfo)

                      // Get the new sequence run and update it with proper metadata
                      viewModel.getSequenceRun(subject.accession, newIndex).foreach { seqRun =>
                        val updatedRun = seqRun.copy(
                          platformName = metadata.platform,
                          instrumentModel = metadata.testType.vendor,
                          testType = metadata.testType.code
                        )
                        viewModel.updateSequenceRun(subject.accession, newIndex, updatedRun)
                        log.info(s"Created sequence run for VCF (${metadata.testType.displayName}): ${dataInput.fileInfo.fileName}")
                      }

                    case _ =>
                      log.debug("VCF metadata dialog cancelled")
                      return
                  }

                case DataType.StrProfile =>
                  // STR CSV import using StrCsvParser
                  val biosampleRef = subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}")
                  StrCsvParser.parse(file, biosampleRef) match {
                    case Right(parseResult) =>
                      viewModel.addStrProfile(subject.accession, parseResult.profile) match {
                        case Right(profileUri) =>
                          log.info(s"Added STR profile with ${parseResult.profile.markers.size} markers for ${subject.accession}")
                          if (parseResult.warnings.nonEmpty) {
                            log.warn(s"STR import warnings: ${parseResult.warnings.mkString("; ")}")
                          }
                        case Left(error) =>
                          log.error(s"Failed to add STR profile: $error")
                          showInfoDialog(t("error.title"), t("data.str_import"), error)
                          return
                      }
                    case Left(error) =>
                      log.error(s"Failed to parse STR file: $error")
                      showInfoDialog(t("error.title"), t("data.str_import"), error)
                      return
                  }

                case DataType.ChipData(_) =>
                  // Chip/array data import using existing parser
                  viewModel.importChipData(subject.accession, file, {
                    case Right(chipProfile) =>
                      log.info(s"Imported chip data: ${chipProfile.vendor} with ${chipProfile.totalMarkersCalled} markers")
                      scalafx.application.Platform.runLater {
                        updateDataSources(subject)
                        showInfoDialog(
                          t("data.title"),
                          t("data.added.success"),
                          s"${chipProfile.vendor} - ${Formatters.formatNumber(chipProfile.totalMarkersCalled)} ${t("data.markers")}"
                        )
                      }
                    case Left(error) =>
                      log.error(s"Failed to import chip data: $error")
                      scalafx.application.Platform.runLater {
                        showInfoDialog(t("error.title"), t("data.chip_import"), error)
                      }
                  })
                  // Chip import is async, so return here - callback will handle UI update
                  return

                case DataType.Unknown =>
                  // Should not happen as dialog disables Add button for unknown types
                  log.warn("Unknown file type - this should not happen")
                  return
              }

              // Refresh the data sources display
              updateDataSources(subject)

              // Show success message
              showInfoDialog(
                t("data.title"),
                t("data.added.success"),
                s"${dataInput.fileInfo.fileName} ${t("data.added.to_subject")}"
              )
            case _ =>
              log.debug("Add data cancelled")
          }
        } catch {
          case e: Exception =>
            log.error(s"Error in handleAddData: ${e.getMessage}", e)
            showInfoDialog(
              t("error.title"),
              t("error.unexpected"),
              e.getMessage
            )
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
