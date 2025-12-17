package com.decodingus.ui.v2

import com.decodingus.i18n.I18n.{t, bind}
import com.decodingus.i18n.Formatters
import com.decodingus.ui.v2.BiosampleExtensions.*
import com.decodingus.util.Logger
import com.decodingus.workspace.WorkbenchViewModel
import com.decodingus.workspace.model.Biosample
import scalafx.Includes.*
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.*
import scalafx.scene.text.{Font, FontWeight, Text}

/**
 * Dashboard view providing a workspace overview.
 *
 * Displays:
 * - Summary statistics (subjects, projects, pending analyses)
 * - Haplogroup distributions
 * - Pending work queue
 * - Recent activity
 */
class DashboardView(viewModel: WorkbenchViewModel) extends ScrollPane {

  private val log = Logger[DashboardView]

  fitToWidth = true
  hbarPolicy = ScrollPane.ScrollBarPolicy.Never
  styleClass += "dashboard-view"
  style = "-fx-background: #1e1e1e; -fx-background-color: #1e1e1e;"

  // ============================================================================
  // Summary Cards
  // ============================================================================

  private val subjectCountCard = createStatCard(
    "dashboard.subjects.count",
    () => viewModel.samples.size.toString,
    "stat-card-subjects"
  )

  private val projectCountCard = createStatCard(
    "dashboard.projects.count",
    () => viewModel.projects.size.toString,
    "stat-card-projects"
  )

  private val ydnaCountCard = createStatCard(
    "dashboard.subjects.with_ydna",
    () => countWithYdna.toString,
    "stat-card-ydna"
  )

  private val pendingCountCard = createStatCard(
    "dashboard.pending_analyses",
    () => countPendingAnalyses.toString,
    "stat-card-pending"
  )

  private val summaryCardsBox = new HBox(20) {
    alignment = Pos.CenterLeft
    padding = Insets(0, 0, 20, 0)
    children = Seq(subjectCountCard, projectCountCard, ydnaCountCard, pendingCountCard)
  }

  // ============================================================================
  // Haplogroup Distribution Charts
  // ============================================================================

  private val ydnaDistributionBox = createHaplogroupDistributionBox(
    "dashboard.ydna_distribution",
    () => getYdnaDistribution,
    "distribution-ydna"
  )

  private val mtdnaDistributionBox = createHaplogroupDistributionBox(
    "dashboard.mtdna_distribution",
    () => getMtdnaDistribution,
    "distribution-mtdna"
  )

  private val distributionsBox = new HBox(30) {
    alignment = Pos.TopLeft
    padding = Insets(0, 0, 20, 0)
    children = Seq(ydnaDistributionBox, mtdnaDistributionBox)
  }

  // ============================================================================
  // Pending Work Section
  // ============================================================================

  private val pendingWorkList = new ListView[String] {
    prefHeight = 150
    placeholder = new Label(t("dashboard.no_pending"))
    styleClass += "pending-work-list"
  }

  private val runAllButton = new Button {
    text = t("dashboard.run_all")
    styleClass += "button-primary"
    disable = true // Enable when there's pending work
    onAction = _ => runAllPendingAnalyses()
  }

  private val pendingWorkSection = new VBox(10) {
    styleClass += "dashboard-section"
    children = Seq(
      createSectionHeader("dashboard.pending_work", Some(runAllButton)),
      pendingWorkList
    )
  }

  // ============================================================================
  // Recent Activity Section
  // ============================================================================

  private val activityList = new ListView[String] {
    prefHeight = 150
    placeholder = new Label(t("info.no_data"))
    styleClass += "activity-list"
  }

  private val activitySection = new VBox(10) {
    styleClass += "dashboard-section"
    children = Seq(
      createSectionHeader("dashboard.recent_activity", None),
      activityList
    )
  }

  // ============================================================================
  // Main Layout
  // ============================================================================

  private val mainContent = new VBox(20) {
    padding = Insets(20)
    style = "-fx-background-color: #1e1e1e;"
    children = Seq(
      createPageHeader("dashboard.title"),
      summaryCardsBox,
      distributionsBox,
      pendingWorkSection,
      activitySection
    )
  }

  content = mainContent

  // ============================================================================
  // Data Binding
  // ============================================================================

  // Update cards when data changes
  viewModel.samples.onChange { (_, _) => updateStats() }
  viewModel.projects.onChange { (_, _) => updateStats() }

  // Initial data load
  updateStats()
  updatePendingWork()

  // ============================================================================
  // Helper Methods
  // ============================================================================

  private def createPageHeader(titleKey: String): HBox = {
    new HBox {
      alignment = Pos.CenterLeft
      children = Seq(
        new Label {
          text <== bind(titleKey)
          styleClass += "page-title"
          style = "-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: #ffffff;"
        }
      )
    }
  }

  private def createSectionHeader(titleKey: String, actionButton: Option[Button]): HBox = {
    new HBox {
      alignment = Pos.CenterLeft
      spacing = 10
      children = Seq(
        new Label {
          text <== bind(titleKey)
          styleClass += "section-title"
          style = "-fx-font-size: 16px; -fx-font-weight: bold; -fx-text-fill: #ffffff;"
        },
        new Region { hgrow = Priority.Always }
      ) ++ actionButton.toSeq
    }
  }

  private def createStatCard(labelKey: String, valueProvider: () => String, styleClassName: String): VBox = {
    val valueLabel = new Label {
      text = valueProvider()
      styleClass += "stat-value"
      style = "-fx-font-size: 32px; -fx-font-weight: bold; -fx-text-fill: #ffffff;"
    }

    val titleLabel = new Label {
      // Extract just the label part (remove the {0} placeholder)
      text = t(labelKey, "").trim
      styleClass += "stat-label"
      style = "-fx-font-size: 12px; -fx-text-fill: #b0b0b0;"
    }

    new VBox(5) {
      alignment = Pos.Center
      padding = Insets(20)
      prefWidth = 140
      prefHeight = 100
      styleClass ++= Seq("stat-card", styleClassName)
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"
      children = Seq(valueLabel, titleLabel)

      // Store reference for updates
      userData = (() => valueLabel.text = valueProvider(), valueProvider)
    }
  }

  private def createHaplogroupDistributionBox(
    titleKey: String,
    dataProvider: () => Seq[(String, Int)],
    styleClassName: String
  ): VBox = {
    val contentBox = new VBox(5) {
      styleClass += "distribution-content"
    }

    def updateContent(): Unit = {
      val data = dataProvider()
      val maxCount = data.map(_._2).maxOption.getOrElse(1)

      contentBox.children = data.take(5).map { case (haplogroup, count) =>
        createDistributionBar(haplogroup, count, maxCount)
      }

      if (data.isEmpty) {
        contentBox.children = Seq(new Label(t("info.no_data")) {
          style = "-fx-text-fill: #999999;"
        })
      }
    }

    updateContent()

    new VBox(10) {
      padding = Insets(15)
      prefWidth = 280
      styleClass ++= Seq("distribution-box", styleClassName)
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"
      children = Seq(
        new Label {
          text <== bind(titleKey)
          styleClass += "distribution-title"
          style = "-fx-font-size: 14px; -fx-font-weight: bold; -fx-text-fill: #ffffff;"
        },
        contentBox
      )

      // Store update function for refresh
      userData = updateContent _
    }
  }

  private def createDistributionBar(label: String, count: Int, maxCount: Int): HBox = {
    val barWidth = (count.toDouble / maxCount * 150).toInt.max(10)

    new HBox(10) {
      alignment = Pos.CenterLeft
      children = Seq(
        new Label(label) {
          prefWidth = 80
          style = "-fx-font-family: monospace; -fx-text-fill: #e0e0e0;"
        },
        new Region {
          prefWidth = barWidth
          prefHeight = 16
          style = "-fx-background-color: #4a9eff; -fx-background-radius: 3;"
        },
        new Label(count.toString) {
          style = "-fx-text-fill: #b0b0b0; -fx-font-size: 11px;"
        }
      )
    }
  }

  private def updateStats(): Unit = {
    // Update stat cards
    val statCards: Seq[VBox] = Seq(subjectCountCard, projectCountCard, ydnaCountCard, pendingCountCard)
    statCards.foreach { card =>
      card.userData match {
        case (updateFn: (() => Unit) @unchecked, _) => updateFn()
        case _ =>
      }
    }

    // Update distribution charts
    val distBoxes: Seq[VBox] = Seq(ydnaDistributionBox, mtdnaDistributionBox)
    distBoxes.foreach { box =>
      box.userData match {
        case updateFn: (() => Unit) @unchecked => updateFn()
        case _ =>
      }
    }
  }

  private def updatePendingWork(): Unit = {
    // For now, show subjects without haplogroup results as pending
    val pending = viewModel.samples.filter { s =>
      s.yHaplogroup.isEmpty && s.mtHaplogroup.isEmpty
    }.map { s =>
      s"${s.donorId.getOrElse(s.accession)} - ${t("analysis.pending")}"
    }.toSeq

    pendingWorkList.items = ObservableBuffer.from(pending)
    runAllButton.disable = pending.isEmpty
  }

  private def runAllPendingAnalyses(): Unit = {
    // TODO: Implement batch analysis
    log.debug("Run all pending analyses")
  }

  // ============================================================================
  // Data Calculations
  // ============================================================================

  private def countWithYdna: Int = {
    viewModel.samples.count(_.yHaplogroup.isDefined)
  }

  private def countPendingAnalyses: Int = {
    viewModel.samples.count { s =>
      s.yHaplogroup.isEmpty || s.mtHaplogroup.isEmpty
    }
  }

  private def getYdnaDistribution: Seq[(String, Int)] = {
    viewModel.samples
      .flatMap(_.yHaplogroup)
      .map(extractTopLevelHaplogroup)
      .groupBy(identity)
      .view.mapValues(_.size)
      .toSeq
      .sortBy(-_._2)
  }

  private def getMtdnaDistribution: Seq[(String, Int)] = {
    viewModel.samples
      .flatMap(_.mtHaplogroup)
      .map(extractTopLevelHaplogroup)
      .groupBy(identity)
      .view.mapValues(_.size)
      .toSeq
      .sortBy(-_._2)
  }

  /**
   * Extract top-level haplogroup (e.g., "R1b-P312" -> "R1b", "H1a1" -> "H")
   */
  private def extractTopLevelHaplogroup(haplogroup: String): String = {
    // For Y-DNA: take up to first hyphen or first digit
    // For mtDNA: take first letter(s) before digits
    val cleaned = haplogroup.trim
    if (cleaned.contains("-")) {
      cleaned.split("-").head
    } else {
      cleaned.takeWhile(c => c.isLetter || c == '*')
    }
  }

  // ============================================================================
  // Public API
  // ============================================================================

  /**
   * Refresh dashboard data.
   */
  def refresh(): Unit = {
    updateStats()
    updatePendingWork()
  }
}
