package com.decodingus.ui.components

import com.decodingus.workspace.model.*
import scalafx.Includes.*
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.*

import java.time.format.DateTimeFormatter

/**
 * Dialog showing haplogroup reconciliation details for a biosample.
 * Displays all run calls, consensus result, and compatibility status.
 */
class ReconciliationDetailDialog(
                                  subject: Biosample,
                                  yDnaReconciliation: Option[HaplogroupReconciliation],
                                  mtDnaReconciliation: Option[HaplogroupReconciliation]
                                ) extends Dialog[Unit] {

  title = s"Haplogroup Reconciliation - ${subject.donorIdentifier}"
  headerText = "Multi-Run Haplogroup Analysis"
  dialogPane().buttonTypes = Seq(ButtonType.Close)
  resizable = true
  dialogPane().setPrefSize(700, 500)

  private val dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm")

  private def createReconciliationPane(
                                        dnaLabel: String,
                                        reconciliationOpt: Option[HaplogroupReconciliation]
                                      ): VBox = {
    new VBox(8) {
      padding = Insets(10)
      style = "-fx-border-color: #E0E0E0; -fx-border-radius: 4; -fx-background-color: #FAFAFA; -fx-background-radius: 4;"

      reconciliationOpt match {
        case None =>
          children = Seq(
            new Label(s"$dnaLabel Haplogroup") {
              style = "-fx-font-size: 14px; -fx-font-weight: bold;"
            },
            new Label("No analysis performed yet") {
              style = "-fx-text-fill: #757575;"
            }
          )

        case Some(reconciliation) =>
          val status = reconciliation.status
          val statusIndicator = createStatusIndicator(status.compatibilityLevel)

          // Header with consensus and status
          val headerBox = new HBox(10) {
            alignment = Pos.CenterLeft
            children = Seq(
              new Label(s"$dnaLabel Haplogroup") {
                style = "-fx-font-size: 14px; -fx-font-weight: bold;"
              },
              statusIndicator,
              new Region {
                HBox.setHgrow(this, Priority.Always)
              },
              new Label(s"${status.runCount} run${if (status.runCount != 1) "s" else ""}") {
                style = "-fx-text-fill: #757575; -fx-font-size: 12px;"
              }
            )
          }

          // Consensus result
          val consensusBox = new HBox(10) {
            alignment = Pos.CenterLeft
            children = Seq(
              new Label("Consensus:") {
                style = "-fx-font-weight: bold;"
              },
              new Label(status.consensusHaplogroup) {
                style = "-fx-font-family: monospace; -fx-font-size: 14px;"
              },
              new Label(f"(${status.confidence * 100}%.1f%% confidence)") {
                style = "-fx-text-fill: #757575;"
              }
            )
          }

          // Warnings if any
          val warningsBox = if (status.warnings.nonEmpty) {
            Some(new VBox(4) {
              children = status.warnings.map { warning =>
                new Label(s"! $warning") {
                  style = "-fx-text-fill: #FF9800; -fx-font-size: 12px;"
                }
              }
            })
          } else None

          // Run calls table
          val runCallsTable = createRunCallsTable(reconciliation.runCalls)

          children = Seq(
            headerBox,
            consensusBox
          ) ++ warningsBox.toSeq ++ Seq(
            new Label("Individual Run Calls:") {
              style = "-fx-font-weight: bold; -fx-padding: 10 0 5 0;"
            },
            runCallsTable
          )
      }
    }
  }

  private def createStatusIndicator(level: CompatibilityLevel): HBox = {
    val (color, symbol, tooltipText) = level match {
      case CompatibilityLevel.COMPATIBLE =>
        ("#4CAF50", "OK", "All runs are compatible - same haplogroup branch")
      case CompatibilityLevel.MINOR_DIVERGENCE =>
        ("#FF9800", "!", "Minor differences - tip-level variations between runs")
      case CompatibilityLevel.MAJOR_DIVERGENCE =>
        ("#F44336", "!!", "Major divergence - branch-level differences detected")
      case CompatibilityLevel.INCOMPATIBLE =>
        ("#9C27B0", "X", "Incompatible - results suggest different individuals")
    }

    new HBox(4) {
      alignment = Pos.CenterLeft
      val indicator = new Label(s"[$symbol]") {
        style = s"-fx-text-fill: $color; -fx-font-size: 12px; -fx-font-weight: bold;"
      }
      val statusLabel = new Label(level.toString.replace("_", " ").toLowerCase.capitalize) {
        style = s"-fx-text-fill: $color; -fx-font-size: 12px;"
        tooltip = Tooltip(tooltipText)
      }
      children = Seq(indicator, statusLabel)
    }
  }

  private def createRunCallsTable(runCalls: List[RunHaplogroupCall]): TableView[RunHaplogroupCall] = {
    val tableData = ObservableBuffer.from(runCalls)

    val table = new TableView[RunHaplogroupCall](tableData) {
      prefHeight = 150
      columnResizePolicy = TableView.ConstrainedResizePolicy

      columns ++= Seq(
        new TableColumn[RunHaplogroupCall, String] {
          text = "Source"
          prefWidth = 120
          cellValueFactory = { p =>
            val source = p.value.sourceRef
            val shortSource = if (source.contains(":")) {
              val parts = source.split(":")
              parts.lastOption.map(s => if (s.length > 12) s.take(12) + "..." else s).getOrElse(source)
            } else if (source.length > 15) source.take(15) + "..." else source
            StringProperty(shortSource)
          }
        },
        new TableColumn[RunHaplogroupCall, String] {
          text = "Technology"
          prefWidth = 80
          cellValueFactory = { p =>
            val tech = p.value.technology.map(_.toString).getOrElse("Unknown")
            StringProperty(tech)
          }
        },
        new TableColumn[RunHaplogroupCall, String] {
          text = "Haplogroup"
          prefWidth = 120
          cellValueFactory = { p =>
            StringProperty(p.value.haplogroup)
          }
        },
        new TableColumn[RunHaplogroupCall, String] {
          text = "SNPs"
          prefWidth = 80
          cellValueFactory = { p =>
            val snps = (p.value.supportingSnps, p.value.conflictingSnps) match {
              case (Some(s), Some(c)) => s"+$s/-$c"
              case (Some(s), None) => s"+$s"
              case _ => "-"
            }
            StringProperty(snps)
          }
        },
        new TableColumn[RunHaplogroupCall, String] {
          text = "Confidence"
          prefWidth = 80
          cellValueFactory = { p =>
            StringProperty(f"${p.value.confidence * 100}%.1f%%")
          }
        },
        new TableColumn[RunHaplogroupCall, String] {
          text = "Tree"
          prefWidth = 80
          cellValueFactory = { p =>
            val tree = p.value.treeProvider.getOrElse("Unknown")
            StringProperty(tree.capitalize)
          }
        }
      )
    }
    VBox.setVgrow(table, Priority.Always)
    table
  }

  // Build the dialog content
  private val dialogContent = new VBox(15) {
    padding = Insets(15)
    children = Seq(
      createReconciliationPane("Y-DNA", yDnaReconciliation),
      createReconciliationPane("mtDNA", mtDnaReconciliation)
    )
  }

  private val scrollPane = new ScrollPane {
    fitToWidth = true
    hbarPolicy = ScrollPane.ScrollBarPolicy.Never
  }
  scrollPane.content = dialogContent

  dialogPane().content = scrollPane

  resultConverter = { _ => () }
}

object ReconciliationDetailDialog {
  /**
   * Creates a compact status indicator widget for use in subject detail views.
   * Returns (indicator label, tooltip text).
   */
  def createCompactStatusIndicator(
                                    yDnaReconciliation: Option[HaplogroupReconciliation],
                                    mtDnaReconciliation: Option[HaplogroupReconciliation]
                                  ): Option[(String, String, String)] = {
    val reconciliations = List(yDnaReconciliation, mtDnaReconciliation).flatten

    if (reconciliations.isEmpty) {
      None
    } else {
      // Find worst compatibility level
      val worstLevel = reconciliations.map(_.status.compatibilityLevel).maxBy {
        case CompatibilityLevel.COMPATIBLE => 0
        case CompatibilityLevel.MINOR_DIVERGENCE => 1
        case CompatibilityLevel.MAJOR_DIVERGENCE => 2
        case CompatibilityLevel.INCOMPATIBLE => 3
      }

      val totalRuns = reconciliations.map(_.status.runCount).sum

      val (color, symbol, statusText) = worstLevel match {
        case CompatibilityLevel.COMPATIBLE =>
          ("#4CAF50", "OK", "Compatible")
        case CompatibilityLevel.MINOR_DIVERGENCE =>
          ("#FF9800", "!", "Minor differences")
        case CompatibilityLevel.MAJOR_DIVERGENCE =>
          ("#F44336", "!!", "Major divergence")
        case CompatibilityLevel.INCOMPATIBLE =>
          ("#9C27B0", "X", "Incompatible")
      }

      val tooltipText = s"$totalRuns analysis run${if (totalRuns != 1) "s" else ""} - $statusText. Click for details."
      Some((symbol, color, tooltipText))
    }
  }
}
