package com.decodingus.ui.components

import com.decodingus.i18n.I18n.t
import com.decodingus.ui.v2.BiosampleExtensions.*
import com.decodingus.workspace.model.*
import scalafx.Includes.*
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.*

/**
 * Dialog for comparing 2-3 subjects side-by-side.
 * Shows Y-STR, Y-DNA, and mtDNA comparisons with genetic distance calculations.
 */
class SubjectCompareDialog(
  subjects: List[Biosample],
  strProfiles: Map[String, StrProfile] // accession -> profile
) extends Dialog[Unit] {

  require(subjects.size >= 2 && subjects.size <= 3, "Must compare 2-3 subjects")

  title = t("compare.title")
  headerText = t("compare.subjects", subjects.size.toString)
  dialogPane().buttonTypes = Seq(ButtonType.Close)
  resizable = true
  dialogPane().setPrefSize(900, 700)

  // Subject names for header
  private val subjectNames = subjects.map(s => s.donorId.getOrElse(s.accession))

  // Tab pane for different comparison views
  private val tabPane = new TabPane {
    tabClosingPolicy = TabPane.TabClosingPolicy.Unavailable
    tabs = Seq(
      createYStrTab(),
      createYDnaTab(),
      createMtDnaTab()
    )
  }

  dialogPane().content = tabPane
  resultConverter = _ => ()

  // ============================================================================
  // Y-STR Comparison Tab
  // ============================================================================

  private def createYStrTab(): Tab = {
    val tab = new Tab {
      text = "Y-STR"
      closable = false
    }

    // Get STR profiles for each subject
    val profiles: List[Option[StrProfile]] = subjects.map(s => strProfiles.get(s.accession))
    val definedProfiles = profiles.collect { case Some(p) => p }

    if (definedProfiles.isEmpty) {
      // No STR data available
      tab.content = createNoDataPane(t("compare.no_str_data"))
      return tab
    }

    // Collect all unique markers across all profiles
    val allMarkers = definedProfiles.flatMap(_.markers.map(_.marker)).distinct.sorted
    val markerCount = allMarkers.size

    // Create comparison rows
    val comparisonRows = allMarkers.map { marker =>
      val values = profiles.map { profileOpt =>
        profileOpt.flatMap(_.markers.find(_.marker == marker)).map(formatStrValue)
      }
      StrComparisonRow(marker, values)
    }

    // Calculate genetic distance (mismatches)
    val (matches, mismatches) = calculateGeneticDistance(comparisonRows)
    val geneticDistance = mismatches

    // Header with subject names
    val headerBox = createSubjectHeader()

    // Summary panel
    val summaryPanel = new HBox(30) {
      alignment = Pos.Center
      padding = Insets(10)
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 8;"
      children = Seq(
        createStatBox(t("compare.markers_compared"), markerCount.toString),
        createStatBox(t("compare.matches"), matches.toString, "#4CAF50"),
        createStatBox(t("compare.mismatches"), mismatches.toString, if (mismatches > 0) "#F44336" else "#4CAF50"),
        createStatBox(t("compare.genetic_distance"), geneticDistance.toString, if (geneticDistance > 5) "#FF9800" else "#4CAF50")
      )
    }

    // Create table
    val tableData = ObservableBuffer.from(comparisonRows)
    val table = new TableView[StrComparisonRow](tableData) {
      columnResizePolicy = TableView.ConstrainedResizePolicy
      style = "-fx-background-color: #333333;"

      // Marker column
      columns += new TableColumn[StrComparisonRow, String] {
        text = t("compare.marker")
        prefWidth = 100
        cellValueFactory = p => StringProperty(p.value.marker)
        cellFactory = { (_: TableColumn[StrComparisonRow, String]) =>
          new TableCell[StrComparisonRow, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                style = "-fx-text-fill: #ffffff; -fx-font-weight: bold; -fx-font-size: 11px;"
              } else {
                text = ""
              }
            }
          }
        }
      }

      // Subject value columns
      subjects.zipWithIndex.foreach { case (subject, idx) =>
        columns += new TableColumn[StrComparisonRow, String] {
          text = subjectNames(idx)
          prefWidth = 120
          cellValueFactory = p => StringProperty(p.value.values(idx).getOrElse("-"))
          cellFactory = { (_: TableColumn[StrComparisonRow, String]) =>
            new TableCell[StrComparisonRow, String] {
              item.onChange { (_, _, newValue) =>
                if (newValue != null) {
                  text = newValue
                  style = "-fx-text-fill: #b0b0b0; -fx-font-family: monospace; -fx-font-size: 12px;"
                } else {
                  text = ""
                }
              }
            }
          }
        }
      }

      // Match status column
      columns += new TableColumn[StrComparisonRow, String] {
        text = t("compare.status")
        prefWidth = 120
        cellValueFactory = p => StringProperty(p.value.matchStatus)
        cellFactory = { (_: TableColumn[StrComparisonRow, String]) =>
          new TableCell[StrComparisonRow, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                val color = if (newValue.contains("Match")) "#4CAF50"
                else if (newValue.contains("1 step")) "#FF9800"
                else if (newValue.contains("step")) "#F44336"
                else if (newValue.contains("-")) "#888888"
                else "#F44336"
                style = s"-fx-text-fill: $color; -fx-font-size: 11px;"
              } else {
                text = ""
              }
            }
          }
        }
      }
    }

    tab.content = new VBox(15) {
      padding = Insets(15)
      children = Seq(headerBox, summaryPanel, table)
      VBox.setVgrow(table, Priority.Always)
    }

    tab
  }

  // ============================================================================
  // Y-DNA Comparison Tab
  // ============================================================================

  private def createYDnaTab(): Tab = {
    val tab = new Tab {
      text = "Y-DNA"
      closable = false
    }

    val yDnaResults: List[Option[HaplogroupResult]] = subjects.map(_.yHaplogroupResult)

    if (yDnaResults.flatten(identity).isEmpty) {
      tab.content = createNoDataPane(t("compare.no_ydna_data"))
      return tab
    }

    val headerBox = createSubjectHeader()

    val haplogroupCards = new HBox(20) {
      alignment = Pos.TopCenter
      padding = Insets(20)
      children = subjects.zipWithIndex.map { case (subject, idx) =>
        createHaplogroupCard(
          subjectNames(idx),
          subject.yHaplogroup,
          subject.yHaplogroupResult.map(_.formattedPath),
          "#4ade80" // Y-DNA green
        )
      }
    }

    // Comparison summary
    val haplogroups = subjects.flatMap(_.yHaplogroup).distinct
    val (summaryText, summaryColor) = if (haplogroups.size == 1) {
      (t("compare.haplogroups_match"), "#4CAF50")
    } else {
      (t("compare.haplogroups_differ"), "#F44336")
    }

    val summaryLabel = new Label(summaryText) {
      style = s"-fx-font-size: 16px; -fx-font-weight: bold; -fx-text-fill: $summaryColor;"
    }

    tab.content = new VBox(20) {
      padding = Insets(15)
      alignment = Pos.TopCenter
      children = Seq(headerBox, haplogroupCards, summaryLabel)
    }

    tab
  }

  // ============================================================================
  // mtDNA Comparison Tab
  // ============================================================================

  private def createMtDnaTab(): Tab = {
    val tab = new Tab {
      text = "mtDNA"
      closable = false
    }

    val mtDnaResults: List[Option[HaplogroupResult]] = subjects.map(_.mtHaplogroupResult)

    if (mtDnaResults.flatten(identity).isEmpty) {
      tab.content = createNoDataPane(t("compare.no_mtdna_data"))
      return tab
    }

    val headerBox = createSubjectHeader()

    val haplogroupCards = new HBox(20) {
      alignment = Pos.TopCenter
      padding = Insets(20)
      children = subjects.zipWithIndex.map { case (subject, idx) =>
        createHaplogroupCard(
          subjectNames(idx),
          subject.mtHaplogroup,
          subject.mtHaplogroupResult.map(_.formattedPath),
          "#60a5fa" // mtDNA blue
        )
      }
    }

    // Comparison summary
    val haplogroups = subjects.flatMap(_.mtHaplogroup).distinct
    val (summaryText, summaryColor) = if (haplogroups.size == 1) {
      (t("compare.haplogroups_match"), "#4CAF50")
    } else {
      (t("compare.haplogroups_differ"), "#F44336")
    }

    val summaryLabel = new Label(summaryText) {
      style = s"-fx-font-size: 16px; -fx-font-weight: bold; -fx-text-fill: $summaryColor;"
    }

    tab.content = new VBox(20) {
      padding = Insets(15)
      alignment = Pos.TopCenter
      children = Seq(headerBox, haplogroupCards, summaryLabel)
    }

    tab
  }

  // ============================================================================
  // Helper Methods
  // ============================================================================

  private def createSubjectHeader(): HBox = {
    new HBox(20) {
      alignment = Pos.Center
      padding = Insets(10)
      children = subjectNames.map { name =>
        new Label(name) {
          style = "-fx-font-size: 14px; -fx-font-weight: bold; -fx-text-fill: #ffffff; -fx-padding: 8 16; -fx-background-color: #3a3a3a; -fx-background-radius: 4;"
        }
      }
    }
  }

  private def createNoDataPane(message: String): VBox = {
    new VBox(15) {
      alignment = Pos.Center
      padding = Insets(40)
      children = Seq(
        new Label(message) {
          style = "-fx-font-size: 14px; -fx-text-fill: #888888;"
        }
      )
    }
  }

  private def createStatBox(label: String, value: String, color: String = "#ffffff"): VBox = {
    new VBox(4) {
      alignment = Pos.Center
      children = Seq(
        new Label(value) {
          style = s"-fx-font-size: 24px; -fx-font-weight: bold; -fx-text-fill: $color;"
        },
        new Label(label) {
          style = "-fx-font-size: 11px; -fx-text-fill: #888888;"
        }
      )
    }
  }

  private def createHaplogroupCard(name: String, haplogroup: Option[String], path: Option[String], color: String): VBox = {
    new VBox(10) {
      alignment = Pos.TopCenter
      padding = Insets(20)
      prefWidth = 250
      style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"
      children = Seq(
        new Label(name) {
          style = "-fx-font-size: 12px; -fx-text-fill: #888888;"
        },
        new Label(haplogroup.getOrElse("-")) {
          style = s"-fx-font-size: 28px; -fx-font-weight: bold; -fx-text-fill: $color;"
        },
        new Label(path.getOrElse("")) {
          style = "-fx-font-size: 10px; -fx-text-fill: #666666;"
          wrapText = true
          maxWidth = 230
        }
      )
    }
  }

  private def formatStrValue(mv: StrMarkerValue): String = {
    mv.value match {
      case SimpleStrValue(repeats) => repeats.toString
      case MultiCopyStrValue(copies) => copies.mkString("-")
      case ComplexStrValue(alleles, raw) => raw.getOrElse(
        alleles.map(a => s"${a.repeats}${a.designation.getOrElse("")}").mkString("-")
      )
    }
  }

  private def calculateGeneticDistance(rows: Seq[StrComparisonRow]): (Int, Int) = {
    var matches = 0
    var mismatches = 0

    rows.foreach { row =>
      val definedValues = row.values.flatten
      if (definedValues.size >= 2) {
        if (definedValues.distinct.size == 1) {
          matches += 1
        } else {
          mismatches += 1
        }
      }
    }

    (matches, mismatches)
  }

  /**
   * Row model for STR comparison table.
   */
  private case class StrComparisonRow(
    marker: String,
    values: Seq[Option[String]]
  ) {
    def matchStatus: String = {
      val definedValues = values.flatten
      if (definedValues.isEmpty) {
        "-"
      } else if (definedValues.size == 1) {
        t("compare.single_value")
      } else if (definedValues.distinct.size == 1) {
        s"[x] ${t("compare.match")}"
      } else {
        // Try to calculate step difference for numeric values
        val numericValues = definedValues.flatMap(s => s.toIntOption)
        if (numericValues.size == definedValues.size && numericValues.nonEmpty) {
          val diff = numericValues.max - numericValues.min
          if (diff == 1) s"[!] 1 ${t("compare.step")}"
          else if (diff == 2) s"[!] 2 ${t("compare.steps")}"
          else s"[X] $diff ${t("compare.steps")}"
        } else {
          s"[X] ${t("compare.diff")}"
        }
      }
    }
  }
}

object SubjectCompareDialog {
  /**
   * Creates a comparison dialog for the given subjects.
   * STR profiles should be pre-loaded and passed as a map.
   */
  def apply(
    subjects: List[Biosample],
    strProfiles: Map[String, StrProfile]
  ): SubjectCompareDialog = new SubjectCompareDialog(subjects, strProfiles)
}
