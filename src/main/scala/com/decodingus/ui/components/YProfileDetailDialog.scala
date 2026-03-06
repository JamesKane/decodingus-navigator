package com.decodingus.ui.components

import com.decodingus.refgenome.YRegionAnnotator
import com.decodingus.yprofile.model.*
import scalafx.Includes.*
import scalafx.beans.property.{ObjectProperty, StringProperty}
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.{HBox, Priority, Region, VBox}
import scalafx.scene.web.WebView

import java.time.format.DateTimeFormatter
import java.util.UUID

/**
 * Comprehensive dialog showing Y Chromosome Profile details.
 * Displays unified profile data with tabs for Summary, Variants, Sources, Concordance, Audit Trail, and Ideogram.
 *
 * @param yRegionAnnotator Optional annotator for displaying chromosome ideogram visualization
 */
class YProfileDetailDialog(
                            profile: YChromosomeProfileEntity,
                            variants: List[YProfileVariantEntity],
                            sources: List[YProfileSourceEntity],
                            variantCalls: Map[UUID, List[YVariantSourceCallEntity]],
                            auditEntries: List[YVariantAuditEntity],
                            biosampleName: String,
                            yRegionAnnotator: Option[YRegionAnnotator] = None
                          ) extends Dialog[Unit] {

  private val dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm")

  title = "Y Chromosome Profile"
  headerText = s"Y Profile for $biosampleName"

  dialogPane().buttonTypes = Seq(ButtonType.Close)
  dialogPane().setPrefSize(1000, 750)

  // Summary Panel
  private val summaryPanel = createSummaryPanel()

  // Tabs
  private val summaryTab = createSummaryTab()
  private val variantsTab = createVariantsTab()
  private val sourcesTab = createSourcesTab()
  private val concordanceTab = createConcordanceTab()
  private val auditTab = createAuditTab()
  private val ideogramTab: Option[Tab] = yRegionAnnotator.map(createIdeogramTab)

  private val tabPane = new TabPane {
    tabs = Seq(summaryTab, variantsTab, sourcesTab, concordanceTab, auditTab) ++ ideogramTab.toSeq
  }
  VBox.setVgrow(tabPane, Priority.Always)

  private val dialogContent = new VBox(10) {
    padding = Insets(15)
    children = Seq(summaryPanel, tabPane)
  }
  VBox.setVgrow(dialogContent, Priority.Always)

  dialogPane().content = dialogContent

  // Make dialog resizable
  dialogPane().getScene.getWindow match {
    case stage: javafx.stage.Stage => stage.setResizable(true)
    case _ =>
  }

  // --- Panel Creation Methods ---

  private def createSummaryPanel(): VBox = {
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
          createStatBox("Total Variants", profile.totalVariants.toString),
          createStatBox("Confirmed", profile.confirmedCount.toString),
          createStatBox("Novel", profile.novelCount.toString),
          createStatBox("Conflict", profile.conflictCount.toString),
          createStatBox("Sources", profile.sourceCount.toString)
        )
      }

      // Callable region progress bar
      val callableBox = profile.callableRegionPct.map { pct =>
        new VBox(4) {
          val progressBar = new ProgressBar {
            progress = pct
            prefWidth = 200
          }
          children = Seq(
            new Label(f"Callable Region: ${pct * 100}%.1f%%") {
              style = "-fx-text-fill: #aaa; -fx-font-size: 11px;"
            },
            progressBar
          )
        }
      }

      val metadataBox = new HBox(20) {
        children = Seq(
          profile.haplogroupTreeProvider.map(p => new Label(s"Tree: $p") {
            style = "-fx-text-fill: #aaa; -fx-font-size: 11px;"
          }),
          profile.haplogroupTreeVersion.map(v => new Label(s"Version: $v") {
            style = "-fx-text-fill: #aaa; -fx-font-size: 11px;"
          }),
          Some(new Label(s"Updated: ${profile.meta.updatedAt.format(dateFormatter)}") {
            style = "-fx-text-fill: #aaa; -fx-font-size: 11px;"
          })
        ).flatten
      }

      children = Seq(haplogroupDisplay, statsBox) ++ callableBox.toSeq ++ Seq(metadataBox)
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

  // --- Tab Creation Methods ---

  private def createSummaryTab(): Tab = {
    val statusBreakdown = new VBox(10) {
      padding = Insets(15)

      // Group variants by status (YVariantStatus, not YConsensusState)
      val byStatus = variants.groupBy(_.status)
      val statusCounts = Seq(
        ("Confirmed", YVariantStatus.CONFIRMED, "#4CAF50"),
        ("Novel", YVariantStatus.NOVEL, "#2196F3"),
        ("Conflict", YVariantStatus.CONFLICT, "#F44336"),
        ("No Coverage", YVariantStatus.NO_COVERAGE, "#9E9E9E"),
        ("Pending", YVariantStatus.PENDING, "#FF9800")
      )

      val statusRows = statusCounts.map { case (name, status, color) =>
        val count = byStatus.get(status).map(_.size).getOrElse(0)
        new HBox(10) {
          alignment = Pos.CenterLeft
          children = Seq(
            new Label("●") {
              style = s"-fx-text-fill: $color; -fx-font-size: 14px;"
            },
            new Label(s"$name: $count") {
              style = "-fx-font-size: 14px;"
            }
          )
        }
      }

      // Variant type breakdown
      val byType = variants.groupBy(_.variantType)
      val typeRows = byType.toSeq.sortBy(-_._2.size).map { case (vType, vs) =>
        new Label(s"${vType.toString}: ${vs.size}") {
          style = "-fx-font-size: 12px; -fx-text-fill: #666;"
        }
      }

      children = Seq(
        new Label("Status Breakdown") {
          style = "-fx-font-size: 16px; -fx-font-weight: bold;"
        }
      ) ++ statusRows ++ Seq(
        new Region {
          prefHeight = 20
        },
        new Label("Variant Types") {
          style = "-fx-font-size: 16px; -fx-font-weight: bold;"
        }
      ) ++ typeRows
    }

    new Tab {
      text = "Summary"
      closable = false
      this.content = new ScrollPane {
        fitToWidth = true
        this.content = statusBreakdown
      }
    }
  }

  private def createVariantsTab(): Tab = {
    // Filter controls
    val statusFilter = new ComboBox[String] {
      items = ObservableBuffer("All", "Confirmed", "Novel", "Conflict", "No Coverage", "Pending")
      value = "All"
      prefWidth = 120
    }

    val buildFilter = new ComboBox[String] {
      items = ObservableBuffer("GRCh38", "GRCh37", "hs1 (CHM13v2)")
      value = "GRCh38"
      prefWidth = 140
    }

    val searchField = new TextField {
      promptText = "Search variants..."
      prefWidth = 200
    }

    val filterBox = new HBox(10) {
      padding = Insets(5)
      alignment = Pos.CenterLeft
      children = Seq(
        new Label("Build:"),
        buildFilter,
        new Label("Status:"),
        statusFilter,
        new Label("Search:"),
        searchField
      )
    }

    // Table data model with multi-reference support
    case class VariantRow(
                           name: String,                    // Canonical name or "Novel-N"
                           position: Option[Long],          // Position in selected build
                           variantType: String,
                           ancAllele: String,               // Ancestral allele
                           derAllele: String,               // Derived allele
                           status: String,
                           call: String,                    // "+" / "-" / "?"
                           isDerived: Option[Boolean],      // For color coding
                           sourceCount: Int,
                           confidence: String,
                           allCoordinates: Map[String, (Long, String, String)], // build -> (pos, anc, der)
                           variantEntity: YProfileVariantEntity // Keep reference for filtering
                         )

    // Helper to normalize build name for lookup
    def normalizeBuildName(build: String): String = build match {
      case "hs1 (CHM13v2)" => "hs1"
      case other => other
    }

    // Helper to get coordinates for a build from variant entity
    def getCoordinatesForBuild(v: YProfileVariantEntity, build: String): Option[(Long, String, String)] = {
      val normalizedBuild = normalizeBuildName(build)
      // First check novelCoordinates map
      v.novelCoordinates.flatMap(_.get(normalizedBuild)).map(nc => (nc.position, nc.ref, nc.alt))
        // Fallback to legacy fields for GRCh38 (default build)
        .orElse {
          if (normalizedBuild == "GRCh38") Some((v.position, v.refAllele, v.altAllele))
          else None
        }
    }

    // Helper to collect all available coordinates
    def getAllCoordinates(v: YProfileVariantEntity): Map[String, (Long, String, String)] = {
      val fromNovelCoords = v.novelCoordinates.getOrElse(Map.empty).map { case (build, nc) =>
        build -> (nc.position, nc.ref, nc.alt)
      }
      // Always include legacy GRCh38 coordinates if not already present
      if (!fromNovelCoords.contains("GRCh38")) {
        fromNovelCoords + ("GRCh38" -> (v.position, v.refAllele, v.altAllele))
      } else {
        fromNovelCoords
      }
    }

    // Helper to determine call symbol and isDerived based on consensus state
    def determineCall(v: YProfileVariantEntity): (String, Option[Boolean]) = v.consensusState match {
      case YConsensusState.DERIVED => ("+", Some(true))
      case YConsensusState.ANCESTRAL => ("-", Some(false))
      case YConsensusState.HETEROPLASMY => ("~", None) // Rare on Y
      case YConsensusState.NO_CALL => ("?", None)
    }

    // Generate Novel-N IDs for unnamed variants, ordered by position
    val sortedVariants = variants.sortBy(v => (v.canonicalName.isEmpty && v.variantName.isEmpty, v.position))
    var novelCounter = 0
    val variantToName: Map[UUID, String] = sortedVariants.map { v =>
      val name = v.canonicalName.orElse(v.variantName).getOrElse {
        novelCounter += 1
        s"Novel-$novelCounter"
      }
      v.id -> name
    }.toMap

    // Build rows for current build selection
    def buildRows(selectedBuild: String): List[VariantRow] = {
      sortedVariants.map { v =>
        val name = variantToName(v.id)
        val coords = getCoordinatesForBuild(v, selectedBuild)
        val allCoords = getAllCoordinates(v)
        val callCount = variantCalls.get(v.id).map(_.size).getOrElse(0)
        val (callSymbol, isDerived) = determineCall(v)

        VariantRow(
          name = name,
          position = coords.map(_._1),
          variantType = v.variantType.toString,
          ancAllele = coords.map(_._2).getOrElse("-"),
          derAllele = coords.map(_._3).getOrElse("-"),
          status = v.status.toString,
          call = callSymbol,
          isDerived = isDerived,
          sourceCount = callCount,
          confidence = f"${v.confidenceScore * 100}%.0f%%",
          allCoordinates = allCoords,
          variantEntity = v
        )
      }
    }

    val tableData = ObservableBuffer.from(buildRows(buildFilter.value.value))

    // Position column with dynamic header
    val positionColumn = new TableColumn[VariantRow, String] {
      text = s"Position (${buildFilter.value.value})"
      cellValueFactory = { r =>
        StringProperty(r.value.position.map(p => f"$p%,d").getOrElse("-"))
      }
      cellFactory = { (_: TableColumn[VariantRow, String]) =>
        new TableCell[VariantRow, String] {
          item.onChange { (_, _, newValue) =>
            text = newValue
            Option(tableRow.value).flatMap(tr => Option(tr.item.value)).foreach { row =>
              if (row.allCoordinates.size > 1) {
                val tooltipText = row.allCoordinates.toSeq.sortBy(_._1).map { case (build, (pos, anc, der)) =>
                  f"$build: chrY:$pos%,d ($anc>$der)"
                }.mkString("\n")
                tooltip = new Tooltip(tooltipText)
              } else {
                tooltip = null
              }
            }
          }
        }
      }
      prefWidth = 140
    }

    // Call column with color coding
    val callColumn = new TableColumn[VariantRow, String] {
      text = "Call"
      cellValueFactory = { r => StringProperty(r.value.call) }
      cellFactory = { (_: TableColumn[VariantRow, String]) =>
        new TableCell[VariantRow, String] {
          item.onChange { (_, _, newValue) =>
            text = newValue
            Option(tableRow.value).flatMap(tr => Option(tr.item.value)).foreach { row =>
              style = row.isDerived match {
                case Some(true) => "-fx-background-color: #22c55e33; -fx-text-fill: #22c55e; -fx-font-weight: bold;"
                case Some(false) => "-fx-background-color: #6b728033; -fx-text-fill: #9ca3af;"
                case None => ""
              }
            }
          }
        }
      }
      prefWidth = 50
    }

    val table = new TableView[VariantRow](tableData) {
      columnResizePolicy = TableView.ConstrainedResizePolicy

      columns ++= Seq(
        new TableColumn[VariantRow, String] {
          text = "Name"
          cellValueFactory = { r => StringProperty(r.value.name) }
          prefWidth = 100
        },
        positionColumn,
        new TableColumn[VariantRow, String] {
          text = "Type"
          cellValueFactory = { r => StringProperty(r.value.variantType) }
          prefWidth = 60
        },
        new TableColumn[VariantRow, String] {
          text = "Anc"
          cellValueFactory = { r => StringProperty(r.value.ancAllele) }
          prefWidth = 50
        },
        new TableColumn[VariantRow, String] {
          text = "Der"
          cellValueFactory = { r => StringProperty(r.value.derAllele) }
          prefWidth = 50
        },
        new TableColumn[VariantRow, String] {
          text = "Status"
          cellValueFactory = { r => StringProperty(r.value.status) }
          prefWidth = 90
        },
        callColumn,
        new TableColumn[VariantRow, Int] {
          text = "Sources"
          cellValueFactory = { r => ObjectProperty(r.value.sourceCount) }
          prefWidth = 60
        },
        new TableColumn[VariantRow, String] {
          text = "Confidence"
          cellValueFactory = { r => StringProperty(r.value.confidence) }
          prefWidth = 80
        }
      )
    }
    VBox.setVgrow(table, Priority.Always)

    // Filter logic
    def applyFilters(): Unit = {
      val selectedBuild = buildFilter.value.value
      val statusValue = statusFilter.value.value
      val searchText = searchField.text.value.toLowerCase

      val statusMapping = Map(
        "Confirmed" -> YVariantStatus.CONFIRMED,
        "Novel" -> YVariantStatus.NOVEL,
        "Conflict" -> YVariantStatus.CONFLICT,
        "No Coverage" -> YVariantStatus.NO_COVERAGE,
        "Pending" -> YVariantStatus.PENDING
      )

      // Rebuild rows with current build selection
      val allRows = buildRows(selectedBuild)

      val filtered = allRows.filter { row =>
        val statusMatch = statusValue == "All" || statusMapping.get(statusValue).contains(row.variantEntity.status)
        val searchMatch = searchText.isEmpty ||
          row.name.toLowerCase.contains(searchText) ||
          row.position.exists(_.toString.contains(searchText))
        statusMatch && searchMatch
      }

      tableData.clear()
      tableData.addAll(filtered: _*)

      // Update position column header
      positionColumn.text = s"Position ($selectedBuild)"
    }

    buildFilter.value.onChange { (_, _, _) => applyFilters() }
    statusFilter.value.onChange { (_, _, _) => applyFilters() }
    searchField.text.onChange { (_, _, _) => applyFilters() }

    new Tab {
      text = s"Variants (${variants.size})"
      closable = false
      this.content = new VBox(5) {
        children = Seq(filterBox, table)
        VBox.setVgrow(table, Priority.Always)
      }
    }
  }

  private def createSourcesTab(): Tab = {
    case class SourceRow(
                          vendor: String,
                          testName: String,
                          sourceType: String,
                          tier: String,
                          meanDepth: String,
                          coverage: String,
                          variantCount: Int,
                          importedAt: String
                        )

    val tableData = ObservableBuffer.from(sources.sortBy(-_.methodTier).map { s =>
      SourceRow(
        s.vendor.getOrElse("-"),
        s.testName.getOrElse("-"),
        s.sourceType.toString,
        s"Tier ${s.methodTier}",
        s.meanReadDepth.map(d => f"$d%.1fx").getOrElse("-"),
        s.coveragePct.map(p => f"${p * 100}%.1f%%").getOrElse("-"),
        s.variantCount,
        s.importedAt.format(dateFormatter)
      )
    })

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
        new TableColumn[SourceRow, String] {
          text = "Depth"
          cellValueFactory = { r => StringProperty(r.value.meanDepth) }
          prefWidth = 80
        },
        new TableColumn[SourceRow, String] {
          text = "Coverage"
          cellValueFactory = { r => StringProperty(r.value.coverage) }
          prefWidth = 80
        },
        new TableColumn[SourceRow, Int] {
          text = "Variants"
          cellValueFactory = { r => ObjectProperty(r.value.variantCount) }
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

    new Tab {
      text = s"Sources (${sources.size})"
      closable = false
      this.content = table
    }
  }

  private def createConcordanceTab(): Tab = {
    // Show variants with multiple sources or conflicts
    val multiSourceVariants = variants.filter { v =>
      variantCalls.get(v.id).exists(_.size > 1) || v.status == YVariantStatus.CONFLICT
    }

    case class ConcordanceRow(
                               variantName: String,
                               position: Long,
                               sourceName: String,
                               calledAllele: String,
                               callState: String,
                               weight: String,
                               depth: String,
                               isVariant: Boolean,
                               indent: Int
                             )

    val rows = multiSourceVariants.flatMap { v =>
      val calls = variantCalls.getOrElse(v.id, Nil)
      val variantRow = ConcordanceRow(
        v.variantName.getOrElse(s"pos:${v.position}"),
        v.position,
        "",
        "",
        v.status.toString,
        "",
        "",
        isVariant = true,
        indent = 0
      )
      val callRows = calls.map { c =>
        val sourceName = sources.find(_.id == c.sourceId).flatMap(_.vendor).getOrElse("Unknown")
        ConcordanceRow(
          "",
          v.position,
          sourceName,
          c.calledAllele,
          c.callState.toString,
          f"${c.concordanceWeight}%.2f",
          c.readDepth.map(_.toString).getOrElse("-"),
          isVariant = false,
          indent = 1
        )
      }
      Seq(variantRow) ++ callRows
    }

    val tableData = ObservableBuffer.from(rows)

    val table = new TableView[ConcordanceRow](tableData) {
      columnResizePolicy = TableView.ConstrainedResizePolicy

      columns ++= Seq(
        new TableColumn[ConcordanceRow, String] {
          text = "Variant / Source"
          cellValueFactory = { r =>
            val prefix = if (r.value.indent > 0) "    └ " else ""
            StringProperty(prefix + (if (r.value.isVariant) r.value.variantName else r.value.sourceName))
          }
          prefWidth = 200
        },
        new TableColumn[ConcordanceRow, String] {
          text = "Called Allele"
          cellValueFactory = { r => StringProperty(r.value.calledAllele) }
          prefWidth = 100
        },
        new TableColumn[ConcordanceRow, String] {
          text = "State"
          cellValueFactory = { r => StringProperty(r.value.callState) }
          prefWidth = 100
        },
        new TableColumn[ConcordanceRow, String] {
          text = "Weight"
          cellValueFactory = { r => StringProperty(r.value.weight) }
          prefWidth = 80
        },
        new TableColumn[ConcordanceRow, String] {
          text = "Depth"
          cellValueFactory = { r => StringProperty(r.value.depth) }
          prefWidth = 80
        }
      )
    }
    VBox.setVgrow(table, Priority.Always)

    val headerLabel = new Label(s"Showing ${multiSourceVariants.size} variants with multiple sources or conflicts") {
      style = "-fx-font-size: 12px; -fx-text-fill: #666;"
      padding = Insets(5)
    }

    new Tab {
      text = "Concordance"
      closable = false
      this.content = new VBox(5) {
        children = Seq(headerLabel, table)
        VBox.setVgrow(table, Priority.Always)
      }
    }
  }

  private def createAuditTab(): Tab = {
    case class AuditRow(
                         timestamp: String,
                         variantName: String,
                         action: String,
                         previousStatus: String,
                         newStatus: String,
                         userId: String,
                         reason: String
                       )

    val tableData = ObservableBuffer.from(auditEntries.sortBy(-_.createdAt.toEpochSecond(java.time.ZoneOffset.UTC)).map { a =>
      val variantName = variants.find(_.id == a.variantId).flatMap(_.variantName).getOrElse(s"Variant ${a.variantId.toString.take(8)}")
      AuditRow(
        a.createdAt.format(dateFormatter),
        variantName,
        a.action.toString,
        a.previousStatus.map(_.toString).getOrElse("-"),
        a.newStatus.map(_.toString).getOrElse("-"),
        a.userId.getOrElse("System"),
        a.reason
      )
    })

    val table = new TableView[AuditRow](tableData) {
      columnResizePolicy = TableView.ConstrainedResizePolicy

      columns ++= Seq(
        new TableColumn[AuditRow, String] {
          text = "Timestamp"
          cellValueFactory = { r => StringProperty(r.value.timestamp) }
          prefWidth = 130
        },
        new TableColumn[AuditRow, String] {
          text = "Variant"
          cellValueFactory = { r => StringProperty(r.value.variantName) }
          prefWidth = 100
        },
        new TableColumn[AuditRow, String] {
          text = "Action"
          cellValueFactory = { r => StringProperty(r.value.action) }
          prefWidth = 100
        },
        new TableColumn[AuditRow, String] {
          text = "Previous"
          cellValueFactory = { r => StringProperty(r.value.previousStatus) }
          prefWidth = 100
        },
        new TableColumn[AuditRow, String] {
          text = "New"
          cellValueFactory = { r => StringProperty(r.value.newStatus) }
          prefWidth = 100
        },
        new TableColumn[AuditRow, String] {
          text = "By"
          cellValueFactory = { r => StringProperty(r.value.userId) }
          prefWidth = 100
        },
        new TableColumn[AuditRow, String] {
          text = "Reason"
          cellValueFactory = { r => StringProperty(r.value.reason) }
          prefWidth = 200
        }
      )
    }
    VBox.setVgrow(table, Priority.Always)

    new Tab {
      text = s"Audit Trail (${auditEntries.size})"
      closable = false
      this.content = table
    }
  }

  /**
   * Create the ideogram tab showing Y chromosome regions and variant positions.
   */
  private def createIdeogramTab(annotator: YRegionAnnotator): Tab = {
    // Convert variants to markers (only show derived/novel/conflict)
    val variantMarkers = variants
      .filter(v => v.consensusState == YConsensusState.DERIVED ||
        v.status == YVariantStatus.NOVEL ||
        v.status == YVariantStatus.CONFLICT)
      .map(YChromosomeIdeogramRenderer.VariantMarker.fromVariantEntity)

    // Generate SVG
    val svgContent = YChromosomeIdeogramRenderer.render(annotator, variantMarkers)
    val statsHtml = YChromosomeIdeogramRenderer.renderStatsHtml(variants, annotator)

    // Wrap in HTML for WebView
    val html =
      s"""<!DOCTYPE html>
         |<html>
         |<head>
         |  <style>
         |    body { margin: 0; padding: 15px; background: #1a1a1a; font-family: system-ui, sans-serif; }
         |    svg { max-width: 100%; height: auto; }
         |  </style>
         |</head>
         |<body>
         |$svgContent
         |$statsHtml
         |</body>
         |</html>""".stripMargin

    val webView = new WebView {
      prefHeight = 350
    }
    webView.engine.loadContent(html)

    new Tab {
      text = "Ideogram"
      closable = false
      this.content = new VBox(10) {
        padding = Insets(10)
        children = Seq(
          new Label("Y Chromosome Region Map") {
            style = "-fx-font-size: 14px; -fx-font-weight: bold;"
          },
          webView
        )
        VBox.setVgrow(webView, Priority.Always)
      }
    }
  }
}
