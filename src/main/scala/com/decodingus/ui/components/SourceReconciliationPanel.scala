package com.decodingus.ui.components

import com.decodingus.i18n.I18n.{t, bind}
import com.decodingus.workspace.model.*
import scalafx.Includes.*
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.*

/**
 * Panel showing multi-run source reconciliation for haplogroup calls.
 * Displays a table of all run calls with their haplogroups, confidence, and
 * a summary of the overall compatibility status.
 */
class SourceReconciliationPanel extends VBox {

  spacing = 10
  padding = Insets(15)
  style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"

  private val tableData = ObservableBuffer.empty[RunHaplogroupCall]
  private var currentReconciliation: Option[HaplogroupReconciliation] = None

  // Header with title and status
  private val titleLabel = new Label(t("reconciliation.source_title")) {
    style = "-fx-font-size: 14px; -fx-font-weight: bold; -fx-text-fill: #ffffff;"
  }

  private val statusIndicator = new Label {
    style = "-fx-font-size: 12px; -fx-font-weight: bold;"
  }

  private val headerBox = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(
      titleLabel,
      new Region { hgrow = Priority.Always },
      statusIndicator
    )
  }

  // Table of run calls
  private val runCallsTable = new TableView[RunHaplogroupCall](tableData) {
    prefHeight = 150
    maxHeight = 200
    columnResizePolicy = TableView.ConstrainedResizePolicy
    // Dark theme styling for table, headers, and rows
    style = """-fx-background-color: #333333; -fx-border-color: #444444;
              |-fx-control-inner-background: #333333;
              |-fx-control-inner-background-alt: #3a3a3a;
              |-fx-table-header-border-color: #444444;
              |-fx-table-cell-border-color: transparent;""".stripMargin
    // Add style class for dark theme table header text
    styleClass.add("dark-reconciliation-table")

    columns ++= Seq(
      new TableColumn[RunHaplogroupCall, String] {
        text = t("reconciliation.source")
        prefWidth = 130
        cellValueFactory = { p =>
          val source = p.value.sourceRef
          val shortSource = extractSourceName(source)
          StringProperty(shortSource)
        }
        cellFactory = { (_: TableColumn[RunHaplogroupCall, String]) =>
          new TableCell[RunHaplogroupCall, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                style = "-fx-text-fill: #ffffff; -fx-font-size: 11px;"
              } else {
                text = ""
              }
            }
          }
        }
      },
      new TableColumn[RunHaplogroupCall, String] {
        text = t("reconciliation.technology")
        prefWidth = 80
        cellValueFactory = { p =>
          val tech = p.value.technology.map(formatTechnology).getOrElse("-")
          StringProperty(tech)
        }
        cellFactory = { (_: TableColumn[RunHaplogroupCall, String]) =>
          new TableCell[RunHaplogroupCall, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                style = "-fx-text-fill: #b0b0b0; -fx-font-size: 11px;"
              } else {
                text = ""
              }
            }
          }
        }
      },
      new TableColumn[RunHaplogroupCall, String] {
        text = t("reconciliation.haplogroup")
        prefWidth = 120
        cellValueFactory = { p =>
          StringProperty(p.value.haplogroup)
        }
        cellFactory = { (_: TableColumn[RunHaplogroupCall, String]) =>
          new TableCell[RunHaplogroupCall, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                style = "-fx-text-fill: #4CAF50; -fx-font-weight: bold; -fx-font-family: monospace; -fx-font-size: 12px;"
              } else {
                text = ""
              }
            }
          }
        }
      },
      new TableColumn[RunHaplogroupCall, String] {
        text = t("reconciliation.snps")
        prefWidth = 80
        cellValueFactory = { p =>
          val snps = (p.value.supportingSnps, p.value.conflictingSnps) match {
            case (Some(s), Some(c)) => s"+$s/-$c"
            case (Some(s), None) => s"+$s"
            case (None, Some(c)) => s"-$c"
            case _ => "-"
          }
          StringProperty(snps)
        }
        cellFactory = { (_: TableColumn[RunHaplogroupCall, String]) =>
          new TableCell[RunHaplogroupCall, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                val color = if (newValue.contains("-") && !newValue.startsWith("+")) "#F44336" else "#b0b0b0"
                style = s"-fx-text-fill: $color; -fx-font-size: 11px; -fx-font-family: monospace;"
              } else {
                text = ""
              }
            }
          }
        }
      },
      new TableColumn[RunHaplogroupCall, String] {
        text = t("reconciliation.quality")
        prefWidth = 100
        cellValueFactory = { p =>
          val quality = formatQuality(p.value.confidence, p.value.meanCoverage)
          StringProperty(quality)
        }
        cellFactory = { (_: TableColumn[RunHaplogroupCall, String]) =>
          new TableCell[RunHaplogroupCall, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                val color = qualityColor(newValue)
                style = s"-fx-text-fill: $color; -fx-font-size: 11px;"
              } else {
                text = ""
              }
            }
          }
        }
      },
      new TableColumn[RunHaplogroupCall, String] {
        text = t("reconciliation.tree")
        prefWidth = 80
        cellValueFactory = { p =>
          val tree = p.value.treeProvider.map(_.capitalize).getOrElse("-")
          StringProperty(tree)
        }
        cellFactory = { (_: TableColumn[RunHaplogroupCall, String]) =>
          new TableCell[RunHaplogroupCall, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                style = "-fx-text-fill: #888888; -fx-font-size: 11px;"
              } else {
                text = ""
              }
            }
          }
        }
      }
    )
  }

  // Placeholder when no data
  private val noDataPlaceholder = new VBox(10) {
    alignment = Pos.Center
    padding = Insets(20)
    children = Seq(
      new Label(t("reconciliation.single_source")) {
        style = "-fx-text-fill: #888888; -fx-font-size: 12px;"
      }
    )
    visible = false
    managed = false
  }

  children = Seq(headerBox, runCallsTable, noDataPlaceholder)

  /**
   * Update the panel with reconciliation data.
   */
  def setReconciliation(reconciliation: Option[HaplogroupReconciliation]): Unit = {
    currentReconciliation = reconciliation

    reconciliation match {
      case Some(recon) if recon.runCalls.size > 1 =>
        // Multiple sources - show table
        tableData.clear()
        tableData ++= recon.runCalls.sortBy(c => -c.confidence)

        runCallsTable.visible = true
        runCallsTable.managed = true
        noDataPlaceholder.visible = false
        noDataPlaceholder.managed = false

        updateStatusIndicator(recon.status)
        this.visible = true
        this.managed = true

      case Some(recon) if recon.runCalls.size == 1 =>
        // Single source - show simplified view
        tableData.clear()
        tableData ++= recon.runCalls

        runCallsTable.visible = true
        runCallsTable.managed = true
        noDataPlaceholder.visible = false
        noDataPlaceholder.managed = false

        statusIndicator.text = s"[${t("reconciliation.single")}]"
        statusIndicator.style = "-fx-font-size: 12px; -fx-text-fill: #888888;"
        this.visible = true
        this.managed = true

      case _ =>
        // No reconciliation data - hide panel
        this.visible = false
        this.managed = false
    }
  }

  private def updateStatusIndicator(status: ReconciliationStatus): Unit = {
    val (color, symbol, label) = status.compatibilityLevel match {
      case CompatibilityLevel.COMPATIBLE =>
        ("#4CAF50", "OK", t("reconciliation.compatible"))
      case CompatibilityLevel.MINOR_DIVERGENCE =>
        ("#FF9800", "!", t("reconciliation.minor_divergence"))
      case CompatibilityLevel.MAJOR_DIVERGENCE =>
        ("#F44336", "!!", t("reconciliation.major_divergence"))
      case CompatibilityLevel.INCOMPATIBLE =>
        ("#9C27B0", "X", t("reconciliation.incompatible"))
    }

    statusIndicator.text = s"[$symbol] $label"
    statusIndicator.style = s"-fx-font-size: 12px; -fx-font-weight: bold; -fx-text-fill: $color;"
    statusIndicator.tooltip = Tooltip(s"${status.runCount} ${t("reconciliation.runs")}")
  }

  private def extractSourceName(sourceRef: String): String = {
    // Extract meaningful name from source reference
    // Format: "local:sequencerun:accession" or "local:chipprofile:accession"
    if (sourceRef.contains(":")) {
      val parts = sourceRef.split(":")
      val name = parts.lastOption.getOrElse(sourceRef)
      if (name.length > 18) name.take(18) + "..." else name
    } else if (sourceRef.length > 20) {
      sourceRef.take(20) + "..."
    } else {
      sourceRef
    }
  }

  private def formatTechnology(tech: HaplogroupTechnology): String = tech match {
    case HaplogroupTechnology.WGS => "WGS"
    case HaplogroupTechnology.WES => "WES"
    case HaplogroupTechnology.BIG_Y => "BigY"
    case HaplogroupTechnology.SNP_ARRAY => "Chip"
    case HaplogroupTechnology.AMPLICON => "Amplicon"
    case HaplogroupTechnology.STR_PANEL => "STR"
  }

  private def formatQuality(confidence: Double, coverage: Option[Double]): String = {
    val stars = (confidence * 5).toInt match {
      case n if n >= 5 => "*****"
      case n if n >= 4 => "****"
      case n if n >= 3 => "***"
      case n if n >= 2 => "**"
      case _ => "*"
    }
    val label = confidence match {
      case c if c >= 0.9 => t("analysis.quality.excellent")
      case c if c >= 0.7 => t("analysis.quality.good")
      case c if c >= 0.5 => t("analysis.quality.fair")
      case _ => t("analysis.quality.poor")
    }
    s"$stars $label"
  }

  private def qualityColor(qualityStr: String): String = {
    if (qualityStr.contains("*****") || qualityStr.toLowerCase.contains("excellent")) "#4CAF50"
    else if (qualityStr.contains("****") || qualityStr.toLowerCase.contains("good")) "#8BC34A"
    else if (qualityStr.contains("***") || qualityStr.toLowerCase.contains("fair")) "#FF9800"
    else "#F44336"
  }
}

object SourceReconciliationPanel {
  /**
   * Creates a panel for displaying Y-DNA or mtDNA reconciliation data.
   */
  def apply(): SourceReconciliationPanel = new SourceReconciliationPanel()
}
