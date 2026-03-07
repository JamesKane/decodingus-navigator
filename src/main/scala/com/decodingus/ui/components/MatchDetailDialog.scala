package com.decodingus.ui.components

import com.decodingus.i18n.I18n.t
import com.decodingus.workspace.model.{IbdSegment, MatchResult, RelationshipEstimate}
import scalafx.Includes.*
import scalafx.beans.property.StringProperty
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.*

/**
 * Dialog showing full details for a confirmed IBD match,
 * including match summary, segment list, and chromosome browser.
 */
class MatchDetailDialog(matchResult: MatchResult) extends Dialog[Unit]:

  title = t("ibd.match_detail_title")
  headerText = matchResult.relationshipEstimate.map(_.label).getOrElse("IBD Match") +
    s" — ${f"${matchResult.totalSharedCm}%.1f"} cM shared"

  private val chromosomeBrowser = new ChromosomeBrowserPanel
  chromosomeBrowser.setMatch(matchResult)

  private val summaryGrid = new GridPane {
    hgap = 15
    vgap = 8
    padding = Insets(10, 0, 10, 0)
  }

  private def addSummaryRow(grid: GridPane, row: Int, label: String, value: String): Unit =
    grid.add(new Label(label) {
      style = "-fx-text-fill: #999999; -fx-font-size: 12px;"
    }, 0, row)
    grid.add(new Label(value) {
      style = "-fx-text-fill: #ffffff; -fx-font-size: 12px; -fx-font-weight: bold;"
    }, 1, row)

  addSummaryRow(summaryGrid, 0, t("ibd.match_name"),
    matchResult.matchedCitizenDid.getOrElse(matchResult.matchedBiosampleRef.takeRight(16)))
  addSummaryRow(summaryGrid, 1, t("ibd.shared_cm"), f"${matchResult.totalSharedCm}%.1f cM")
  addSummaryRow(summaryGrid, 2, t("ibd.segments"), matchResult.segmentCount.toString)
  addSummaryRow(summaryGrid, 3, t("ibd.longest"),
    matchResult.longestSegmentCm.map(cm => f"$cm%.1f cM").getOrElse("-"))
  addSummaryRow(summaryGrid, 4, t("ibd.relationship"),
    matchResult.relationshipEstimate.map(_.label).getOrElse("Unknown"))
  matchResult.attestationHash.foreach { hash =>
    addSummaryRow(summaryGrid, 5, t("ibd.attestation_hash"), hash.take(16) + "...")
  }

  // Segment table
  private val segmentTable = new TableView[IbdSegment] {
    prefHeight = 180
    style = "-fx-background-color: #333333;"
    columnResizePolicy = TableView.ConstrainedResizePolicy

    columns ++= Seq(
      new TableColumn[IbdSegment, String] {
        text = "Chr"
        prefWidth = 50
        cellValueFactory = p => StringProperty(p.value.chromosome)
      },
      new TableColumn[IbdSegment, String] {
        text = "Start"
        prefWidth = 100
        cellValueFactory = p => StringProperty(ChromosomeBrowserRenderer.formatBp(p.value.startPosition))
      },
      new TableColumn[IbdSegment, String] {
        text = "End"
        prefWidth = 100
        cellValueFactory = p => StringProperty(ChromosomeBrowserRenderer.formatBp(p.value.endPosition))
      },
      new TableColumn[IbdSegment, String] {
        text = "Length (cM)"
        prefWidth = 80
        cellValueFactory = p => StringProperty(f"${p.value.lengthCm}%.2f")
      },
      new TableColumn[IbdSegment, String] {
        text = "SNPs"
        prefWidth = 60
        cellValueFactory = p => StringProperty(p.value.snpCount.map(_.toString).getOrElse("-"))
      }
    )

    items = scalafx.collections.ObservableBuffer.from(
      matchResult.sharedSegments.sortBy(s => (chrSortKey(s.chromosome), s.startPosition))
    )
  }

  private val content = new VBox(15) {
    padding = Insets(15)
    prefWidth = 750
    children = Seq(
      summaryGrid,
      new Separator,
      chromosomeBrowser,
      new Label(t("ibd.shared_segments_list")) {
        style = "-fx-font-weight: bold; -fx-text-fill: #ffffff;"
      },
      segmentTable
    )
  }

  dialogPane().content = content
  dialogPane().buttonTypes = Seq(ButtonType.Close)
  dialogPane().setPrefWidth(800)

  // Dark theme for the dialog
  dialogPane().setStyle("-fx-background-color: #1e1e1e;")

  private def chrSortKey(chr: String): Int =
    chr.toIntOption.getOrElse(chr match
      case "X" => 23
      case "Y" => 24
      case _ => 99
    )
