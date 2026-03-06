package com.decodingus.ui.components

import com.decodingus.haplogroup.model.HaplogroupResult as AnalysisHaplogroupResult
import com.decodingus.haplogroup.tree.TreeType
import com.decodingus.workspace.model.{PrivateVariantData, VariantCall, HaplogroupResult as WorkspaceHaplogroupResult}
import javafx.collections.transformation.FilteredList
import scalafx.Includes.*
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.*

import java.nio.file.{Files, Path}
import scala.io.Source
import scala.util.Using

/**
 * Comprehensive dialog showing haplogroup analysis results in a GUI format.
 * Parses and displays cached report files with tabs for candidates, lineage, SNP details, and private variants.
 * Can display from either workspace model or cached report file.
 */
class HaplogroupReportDialog(
                              treeType: TreeType,
                              workspaceResult: Option[WorkspaceHaplogroupResult] = None,
                              analysisResults: Option[List[AnalysisHaplogroupResult]] = None,
                              artifactDir: Option[Path] = None,
                              sampleName: Option[String] = None
                            ) extends Dialog[Unit] {

  private val dnaType = if (treeType == TreeType.YDNA) "Y-DNA" else "mtDNA"

  title = s"$dnaType Haplogroup Results"
  headerText = s"$dnaType Haplogroup Analysis"

  dialogPane().buttonTypes = Seq(ButtonType.OK)
  dialogPane().setPrefSize(900, 700)

  // Try to parse cached report if available
  private val reportData: Option[ParsedReport] = artifactDir.flatMap { dir =>
    val prefix = if (treeType == TreeType.YDNA) "ydna" else "mtdna"
    val reportFile = dir.resolve(s"${prefix}_haplogroup_report.txt")
    if (Files.exists(reportFile)) {
      Some(parseReport(reportFile))
    } else None
  }

  // Format large numbers with commas
  private def formatNumber(n: Int): String = f"$n%,d"

  // Score to confidence display
  private def scoreToConfidence(score: Double, depth: Int): String = {
    val normalizedScore = score / math.max(depth, 1)
    if (normalizedScore > 25) "Very High"
    else if (normalizedScore > 20) "High"
    else if (normalizedScore > 15) "Moderate"
    else if (normalizedScore > 10) "Low"
    else "Very Low"
  }

  // Quality to stars
  private def qualityToStars(quality: Option[Double]): String = {
    quality match {
      case None => "-"
      case Some(q) if q < 10 => "☆☆☆☆☆"
      case Some(q) if q < 20 => "★☆☆☆☆"
      case Some(q) if q < 30 => "★★☆☆☆"
      case Some(q) if q < 40 => "★★★☆☆"
      case Some(q) if q < 50 => "★★★★☆"
      case Some(_) => "★★★★★"
    }
  }

  // Summary Panel
  private val summaryPanel = createSummaryPanel()

  // Candidates Table
  private val candidatesTab = createCandidatesTab()

  // Lineage Path Tab
  private val lineageTab = createLineageTab()

  // SNP Details Tab
  private val snpDetailsTab = createSnpDetailsTab()

  // Private Variant Tabs (3 separate tabs)
  private val novelSnpsTab = createNovelSnpsTab()
  private val strIndelsTab = createStrIndelsTab()
  private val otherIndelsTab = createOtherIndelsTab()

  private val tabPane = new TabPane {
    tabs = Seq(candidatesTab, lineageTab, snpDetailsTab, novelSnpsTab, strIndelsTab, otherIndelsTab).flatten
  }
  VBox.setVgrow(tabPane, Priority.Always)

  private val content = new VBox(10) {
    padding = Insets(15)
    children = Seq(summaryPanel, tabPane)
  }
  VBox.setVgrow(content, Priority.Always)

  dialogPane().content = content

  // Make dialog resizable
  dialogPane().getScene.getWindow match {
    case stage: javafx.stage.Stage => stage.setResizable(true)
    case _ =>
  }

  // --- Panel Creation Methods ---

  private def createSummaryPanel(): VBox = {
    val topResult = reportData.flatMap(_.topCandidate)
      .orElse(analysisResults.flatMap(_.headOption).map(r => CandidateRow(r.name, r.score, r.matchingSnps, r.ancestralMatches, r.noCalls, r.depth)))
      .orElse(workspaceResult.map(r => CandidateRow(r.haplogroupName, r.score, r.matchingSnps.getOrElse(0), r.ancestralMatches.getOrElse(0), 0, r.treeDepth.getOrElse(0))))

    val metadata = reportData.map(_.metadata).getOrElse(Map.empty)

    new VBox(8) {
      padding = Insets(15)
      style = "-fx-background-color: linear-gradient(to bottom, #2d5a2d, #1a3a1a); -fx-background-radius: 8;"

      children = topResult match {
        case Some(result) =>
          val confidence = scoreToConfidence(result.score, result.depth)
          val confidenceColor = confidence match {
            case "Very High" => "#00ff00"
            case "High" => "#88ff00"
            case "Moderate" => "#ffff00"
            case "Low" => "#ff8800"
            case _ => "#ff4400"
          }

          Seq(
            new HBox(20) {
              alignment = Pos.CenterLeft
              children = Seq(
                new Label(result.haplogroup) {
                  style = "-fx-font-size: 28px; -fx-font-weight: bold; -fx-text-fill: white;"
                },
                new Label(s"Confidence: $confidence") {
                  style = s"-fx-font-size: 14px; -fx-text-fill: $confidenceColor; -fx-font-weight: bold;"
                }
              )
            },
            new HBox(30) {
              children = Seq(
                createStatBox("Score", f"${result.score}%.0f"),
                createStatBox("Derived SNPs", formatNumber(result.derived)),
                createStatBox("Ancestral", formatNumber(result.ancestral)),
                createStatBox("No Calls", formatNumber(result.noCalls)),
                createStatBox("Tree Depth", result.depth.toString)
              )
            },
            new HBox(20) {
              children = Seq(
                metadata.get("treeProvider").map(p => new Label(s"Tree: $p") {
                  style = "-fx-text-fill: #aaa; -fx-font-size: 11px;"
                }),
                metadata.get("treeBuild").map(b => new Label(s"Build: $b") {
                  style = "-fx-text-fill: #aaa; -fx-font-size: 11px;"
                }),
                metadata.get("liftover").map(l => new Label(s"Liftover: $l") {
                  style = "-fx-text-fill: #aaa; -fx-font-size: 11px;"
                }),
                sampleName.map(n => new Label(s"Sample: $n") {
                  style = "-fx-text-fill: #aaa; -fx-font-size: 11px;"
                })
              ).flatten
            }
          )

        case None =>
          Seq(
            new Label("No haplogroup could be determined") {
              style = "-fx-font-size: 18px; -fx-text-fill: #ff6666;"
            }
          )
      }
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

  private def createCandidatesTab(): Option[Tab] = {
    val candidates = reportData.map(_.candidates)
      .orElse(analysisResults.map(_.take(20).map(r => CandidateRow(r.name, r.score, r.matchingSnps, r.ancestralMatches, r.noCalls, r.depth))))
      .getOrElse(List.empty)

    if (candidates.isEmpty) return None

    val tableData = ObservableBuffer.from(candidates)

    val table = new TableView[CandidateRow](tableData) {
      columnResizePolicy = TableView.ConstrainedResizePolicy

      columns ++= Seq(
        new TableColumn[CandidateRow, String] {
          text = "Rank"
          cellValueFactory = r => StringProperty((tableData.indexOf(r.value) + 1).toString)
          prefWidth = 50
        },
        new TableColumn[CandidateRow, String] {
          text = "Haplogroup"
          cellValueFactory = r => StringProperty(r.value.haplogroup)
          prefWidth = 180
        },
        new TableColumn[CandidateRow, String] {
          text = "Score"
          cellValueFactory = r => StringProperty(f"${r.value.score}%.1f")
          prefWidth = 80
        },
        new TableColumn[CandidateRow, String] {
          text = "Derived"
          cellValueFactory = r => StringProperty(formatNumber(r.value.derived))
          prefWidth = 80
        },
        new TableColumn[CandidateRow, String] {
          text = "Ancestral"
          cellValueFactory = r => StringProperty(formatNumber(r.value.ancestral))
          prefWidth = 80
        },
        new TableColumn[CandidateRow, String] {
          text = "No Calls"
          cellValueFactory = r => StringProperty(formatNumber(r.value.noCalls))
          prefWidth = 80
        },
        new TableColumn[CandidateRow, String] {
          text = "Depth"
          cellValueFactory = r => StringProperty(r.value.depth.toString)
          prefWidth = 60
        }
      )
    }
    VBox.setVgrow(table, Priority.Always)

    val tab = new Tab {
      text = "Top Candidates"
      closable = false
    }
    tab.content = new VBox(10) {
      padding = Insets(10)
      children = Seq(
        new Label(s"Top ${candidates.size} candidate haplogroups by score:") {
          style = "-fx-font-weight: bold;"
        },
        table
      )
    }
    Some(tab)
  }

  private def createLineageTab(): Option[Tab] = {
    val pathNodes = reportData.map(_.lineagePath)
      .orElse(workspaceResult.flatMap(_.lineagePath).map(_.map(name => LineageNode(name, 0, ""))))
      .getOrElse(List.empty)

    if (pathNodes.isEmpty) return None

    // Build tree structure
    def buildTreeItem(nodes: List[LineageNode], depth: Int): Option[TreeItem[String]] = {
      nodes.filter(_.depth == depth) match {
        case Nil => None
        case node :: _ =>
          val label = if (node.derivedInfo.nonEmpty) s"${node.name} ${node.derivedInfo}" else node.name
          val item = new TreeItem[String](label) {
            expanded = true
          }
          buildTreeItem(nodes, depth + 1).foreach(child => item.children.add(child))
          Some(item)
      }
    }

    val rootItem = buildTreeItem(pathNodes, 0).getOrElse(new TreeItem[String]("(empty)"))

    val treeView = new TreeView[String](rootItem) {
      showRoot = true
    }
    VBox.setVgrow(treeView, Priority.Always)

    val tab = new Tab {
      text = "Lineage Path"
      closable = false
    }
    tab.content = new VBox(10) {
      padding = Insets(10)
      children = Seq(
        new Label("Phylogenetic path from root to predicted haplogroup:") {
          style = "-fx-font-weight: bold;"
        },
        treeView
      )
    }
    Some(tab)
  }

  // Helper to get parsed variant data, separated into SNPs, STR indels, and other indels
  private lazy val parsedVariants: (List[PrivateVariantRow], List[PrivateVariantRow], List[PrivateVariantRow]) = {
    val (snps, allIndels) = reportData.map(r => (r.privateSnps, r.privateIndels))
      .orElse(workspaceResult.flatMap(_.privateVariants).map { pvd =>
        val (s, i) = pvd.variants.partition(v => v.referenceAllele.length == 1 && v.alternateAllele.length == 1)
        (s.map(v => PrivateVariantRow(v.contigAccession, v.position, v.referenceAllele, v.alternateAllele, qualityToStars(v.quality), None)),
          i.map(v => PrivateVariantRow(v.contigAccession, v.position, v.referenceAllele, v.alternateAllele, qualityToStars(v.quality), None)))
      })
      .getOrElse((List.empty, List.empty))

    // Separate STR indels (have strInfo) from other indels
    val (strIndels, otherIndels) = allIndels.partition(_.strInfo.isDefined)
    (snps, strIndels, otherIndels)
  }

  /** Creates a filterable variant table with contig and quality filters */
  private def createFilterableVariantTab(
                                          title: String,
                                          variants: List[PrivateVariantRow],
                                          showStrInfo: Boolean
                                        ): Option[Tab] = {
    if (variants.isEmpty) return None

    val sourceData = ObservableBuffer.from(variants)
    val filteredList = new FilteredList[PrivateVariantRow](sourceData)

    // Get unique contigs for filter dropdown
    val contigs = variants.map(_.contig).distinct.sorted
    val qualityLevels = List("All", "★★★★★ (5)", "★★★★☆+ (4+)", "★★★☆☆+ (3+)", "★★☆☆☆+ (2+)", "★☆☆☆☆+ (1+)")

    // Filter controls
    val contigFilter = new ComboBox[String](ObservableBuffer.from("All" +: contigs)) {
      value = "All"
      prefWidth = 100
    }

    val qualityFilter = new ComboBox[String](ObservableBuffer.from(qualityLevels)) {
      value = "All"
      prefWidth = 120
    }

    val positionFilter = new TextField {
      promptText = "Position contains..."
      prefWidth = 150
    }

    val refAltFilter = new TextField {
      promptText = "Ref/Alt contains..."
      prefWidth = 120
    }

    // Count label that updates with filter
    val countLabel = new Label(s"Showing ${variants.size} of ${variants.size}")

    // Apply filters
    def updateFilter(): Unit = {
      val contigValue = contigFilter.value.value
      val qualityValue = qualityFilter.value.value
      val posValue = positionFilter.text.value.trim
      val refAltValue = refAltFilter.text.value.trim.toUpperCase

      filteredList.setPredicate { row =>
        val contigMatch = contigValue == "All" || row.contig == contigValue

        val qualityMatch = qualityValue match {
          case "All" => true
          case s if s.contains("(5)") => row.quality.count(_ == '★') >= 5
          case s if s.contains("(4+)") => row.quality.count(_ == '★') >= 4
          case s if s.contains("(3+)") => row.quality.count(_ == '★') >= 3
          case s if s.contains("(2+)") => row.quality.count(_ == '★') >= 2
          case s if s.contains("(1+)") => row.quality.count(_ == '★') >= 1
          case _ => true
        }

        val posMatch = posValue.isEmpty || row.position.toString.contains(posValue)
        val refAltMatch = refAltValue.isEmpty ||
          row.ref.toUpperCase.contains(refAltValue) ||
          row.alt.toUpperCase.contains(refAltValue)

        contigMatch && qualityMatch && posMatch && refAltMatch
      }

      countLabel.text = s"Showing ${filteredList.size} of ${variants.size}"
    }

    // Bind filter updates
    contigFilter.value.onChange { (_, _, _) => updateFilter() }
    qualityFilter.value.onChange { (_, _, _) => updateFilter() }
    positionFilter.text.onChange { (_, _, _) => updateFilter() }
    refAltFilter.text.onChange { (_, _, _) => updateFilter() }

    val filterBar = new HBox(10) {
      alignment = Pos.CenterLeft
      children = Seq(
        new Label("Contig:"),
        contigFilter,
        new Label("Min Quality:"),
        qualityFilter,
        new Label("Position:"),
        positionFilter,
        new Label("Ref/Alt:"),
        refAltFilter,
        new Region {
          HBox.setHgrow(this, Priority.Always)
        },
        countLabel
      )
    }

    // Create table with filtered data
    val table = new TableView[PrivateVariantRow]() {
      columnResizePolicy = TableView.ConstrainedResizePolicy
      delegate.setItems(filteredList)

      columns ++= Seq(
        new TableColumn[PrivateVariantRow, String] {
          text = "Contig"
          cellValueFactory = r => StringProperty(r.value.contig)
          prefWidth = 70
        },
        new TableColumn[PrivateVariantRow, String] {
          text = "Position"
          cellValueFactory = r => StringProperty(formatNumber(r.value.position))
          prefWidth = 100
        },
        new TableColumn[PrivateVariantRow, String] {
          text = "Ref"
          cellValueFactory = r => StringProperty(if (r.value.ref.length > 10) r.value.ref.take(8) + ".." else r.value.ref)
          prefWidth = 90
        },
        new TableColumn[PrivateVariantRow, String] {
          text = "Alt"
          cellValueFactory = r => StringProperty(if (r.value.alt.length > 10) r.value.alt.take(8) + ".." else r.value.alt)
          prefWidth = 90
        },
        new TableColumn[PrivateVariantRow, String] {
          text = "Quality"
          cellValueFactory = r => StringProperty(r.value.quality)
          prefWidth = 85
        }
      )

      // Add optional depth and region columns if data is available
      val hasDepthData = variants.exists(_.readDepth.isDefined)
      val hasRegionData = variants.exists(_.region.isDefined)

      if (hasDepthData) {
        // Insert depth column before quality
        val qualityCol = columns.last
        columns.remove(columns.size - 1)
        columns += new TableColumn[PrivateVariantRow, String] {
          text = "Depth"
          cellValueFactory = r => StringProperty(r.value.readDepth.getOrElse("-"))
          prefWidth = 60
        }
        columns += qualityCol
      }

      if (hasRegionData) {
        // Insert region column before quality
        val qualityCol = columns.last
        columns.remove(columns.size - 1)
        columns += new TableColumn[PrivateVariantRow, String] {
          text = "Region"
          cellValueFactory = r => StringProperty(r.value.region.getOrElse("-"))
          prefWidth = 120
        }
        columns += qualityCol
      }

      if (showStrInfo) {
        columns += new TableColumn[PrivateVariantRow, String] {
          text = "STR Type"
          cellValueFactory = r => StringProperty(r.value.strInfo.getOrElse("-"))
          prefWidth = 200
        }
      }
    }
    VBox.setVgrow(table, Priority.Always)

    val tab = new Tab {
      text = title
      closable = false
    }
    tab.content = new VBox(10) {
      padding = Insets(10)
      children = Seq(filterBar, table)
    }
    Some(tab)
  }

  private def createNovelSnpsTab(): Option[Tab] = {
    val (snps, _, _) = parsedVariants
    createFilterableVariantTab(s"Novel SNPs (${snps.size})", snps, showStrInfo = false)
  }

  private def createStrIndelsTab(): Option[Tab] = {
    val (_, strIndels, _) = parsedVariants
    createFilterableVariantTab(s"STR Indels (${strIndels.size})", strIndels, showStrInfo = true)
  }

  private def createOtherIndelsTab(): Option[Tab] = {
    val (_, _, otherIndels) = parsedVariants
    createFilterableVariantTab(s"Other Indels (${otherIndels.size})", otherIndels, showStrInfo = false)
  }

  private def createSnpDetailsTab(): Option[Tab] = {
    val snpDetails = reportData.map(_.snpDetails).getOrElse(List.empty)

    if (snpDetails.isEmpty) return None

    val tableData = ObservableBuffer.from(snpDetails)

    val table = new TableView[SnpDetailRow](tableData) {
      columnResizePolicy = TableView.ConstrainedResizePolicy

      columns ++= Seq(
        new TableColumn[SnpDetailRow, String] {
          text = "Contig"
          cellValueFactory = r => StringProperty(r.value.contig)
          prefWidth = 70
        },
        new TableColumn[SnpDetailRow, String] {
          text = "Position"
          cellValueFactory = r => StringProperty(formatNumber(r.value.position))
          prefWidth = 100
        },
        new TableColumn[SnpDetailRow, String] {
          text = "SNP Name"
          cellValueFactory = r => StringProperty(r.value.snpName)
          prefWidth = 150
        },
        new TableColumn[SnpDetailRow, String] {
          text = "Anc"
          cellValueFactory = r => StringProperty(r.value.ancestral)
          prefWidth = 50
        },
        new TableColumn[SnpDetailRow, String] {
          text = "Der"
          cellValueFactory = r => StringProperty(r.value.derived)
          prefWidth = 50
        },
        new TableColumn[SnpDetailRow, String] {
          text = "Call"
          cellValueFactory = r => StringProperty(r.value.call)
          prefWidth = 50
        },
        new TableColumn[SnpDetailRow, String] {
          text = "State"
          cellValueFactory = r => StringProperty(r.value.state)
          prefWidth = 80
        }
      )

      // Add optional columns if data is available
      val hasDepthData = snpDetails.exists(_.readDepth.isDefined)
      val hasRegionData = snpDetails.exists(_.region.isDefined)
      val hasQualityData = snpDetails.exists(_.quality.isDefined)

      if (hasDepthData) {
        columns += new TableColumn[SnpDetailRow, String] {
          text = "Depth"
          cellValueFactory = r => StringProperty(r.value.readDepth.getOrElse("-"))
          prefWidth = 60
        }
      }
      if (hasRegionData) {
        columns += new TableColumn[SnpDetailRow, String] {
          text = "Region"
          cellValueFactory = r => StringProperty(r.value.region.getOrElse("-"))
          prefWidth = 120
        }
      }
      if (hasQualityData) {
        columns += new TableColumn[SnpDetailRow, String] {
          text = "Quality"
          cellValueFactory = r => StringProperty(r.value.quality.getOrElse("-"))
          prefWidth = 85
        }
      }
    }
    VBox.setVgrow(table, Priority.Always)

    val tab = new Tab {
      text = s"SNP Details (${snpDetails.size})"
      closable = false
    }
    tab.content = new VBox(10) {
      padding = Insets(10)
      children = Seq(
        new Label("SNPs along the predicted haplogroup path:") {
          style = "-fx-font-weight: bold;"
        },
        table
      )
    }
    Some(tab)
  }

  // --- Report Parsing ---

  private case class ParsedReport(
                                   metadata: Map[String, String],
                                   topCandidate: Option[CandidateRow],
                                   candidates: List[CandidateRow],
                                   lineagePath: List[LineageNode],
                                   snpDetails: List[SnpDetailRow],
                                   privateSnps: List[PrivateVariantRow],
                                   privateIndels: List[PrivateVariantRow]
                                 )

  private case class CandidateRow(haplogroup: String, score: Double, derived: Int, ancestral: Int, noCalls: Int, depth: Int)

  private case class LineageNode(name: String, depth: Int, derivedInfo: String)

  private case class SnpDetailRow(
                                   contig: String,
                                   position: Int,
                                   snpName: String,
                                   ancestral: String,
                                   derived: String,
                                   call: String,
                                   state: String,
                                   readDepth: Option[String] = None,
                                   region: Option[String] = None,
                                   quality: Option[String] = None
                                 )

  private case class PrivateVariantRow(
                                        contig: String,
                                        position: Int,
                                        ref: String,
                                        alt: String,
                                        quality: String,
                                        strInfo: Option[String],
                                        readDepth: Option[String] = None,
                                        region: Option[String] = None
                                      )

  private def parseReport(reportPath: Path): ParsedReport = {
    val lines = Using.resource(Source.fromFile(reportPath.toFile))(_.getLines().toList)

    var metadata = Map.empty[String, String]
    var candidates = List.empty[CandidateRow]
    var lineagePath = List.empty[LineageNode]
    var snpDetails = List.empty[SnpDetailRow]
    var privateSnps = List.empty[PrivateVariantRow]
    var privateIndels = List.empty[PrivateVariantRow]

    var currentSection = ""
    var inSnpSection = false
    var inNovelSnpSection = false
    var inNovelIndelSection = false

    lines.foreach { line =>
      val trimmed = line.trim

      // Section detection
      if (trimmed == "HAPLOGROUP PREDICTION") currentSection = "prediction"
      else if (trimmed == "TOP 10 CANDIDATES") currentSection = "candidates"
      else if (trimmed == "HAPLOGROUP PATH") currentSection = "path"
      else if (trimmed == "SNP DETAILS (along predicted path)") {
        currentSection = "snp_details"; inSnpSection = false
      }
      else if (trimmed == "NOVEL/UNPLACED SNPs") {
        currentSection = "novel_snps"; inNovelSnpSection = false
      }
      else if (trimmed == "NOVEL/UNPLACED INDELS") {
        currentSection = "novel_indels"; inNovelIndelSection = false
      }
      else if (trimmed == "SUMMARY STATISTICS") currentSection = "summary"

      // Metadata parsing
      if (trimmed.startsWith("Tree Provider:")) metadata += ("treeProvider" -> trimmed.stripPrefix("Tree Provider:").trim)
      if (trimmed.startsWith("Tree Build:")) metadata += ("treeBuild" -> trimmed.stripPrefix("Tree Build:").trim)
      if (trimmed.startsWith("Sample Build:")) metadata += ("sampleBuild" -> trimmed.stripPrefix("Sample Build:").trim)
      if (trimmed.startsWith("Liftover:")) metadata += ("liftover" -> trimmed.stripPrefix("Liftover:").trim)

      // Candidates parsing
      if (currentSection == "candidates" && trimmed.matches("^\\d+\\s+.*")) {
        val parts = trimmed.split("\\s+")
        if (parts.length >= 6) {
          try {
            candidates = candidates :+ CandidateRow(
              parts(1),
              parts(2).toDouble,
              parts(3).toInt,
              parts(4).toInt,
              0, // no calls not in this section
              parts(5).toInt
            )
          } catch {
            case _: Exception =>
          }
        }
      }

      // Path parsing - detect indentation
      if (currentSection == "path" && !trimmed.startsWith("-") && trimmed.nonEmpty && !trimmed.startsWith("HAPLOGROUP")) {
        val depth = (line.length - line.stripLeading.length) / 2
        val derivedMatch = "\\[([^\\]]+)\\]".r.findFirstMatchIn(trimmed)
        val name = trimmed.split("\\s+")(0)
        val derivedInfo = derivedMatch.map(_.group(0)).getOrElse("")
        lineagePath = lineagePath :+ LineageNode(name, depth, derivedInfo)
      }

      // SNP details parsing - handles both old (7 cols) and new (10 cols with depth/region/quality) formats
      if (currentSection == "snp_details") {
        if (trimmed.startsWith("Contig")) inSnpSection = true
        else if (inSnpSection && trimmed.matches("^[A-Za-z0-9]+\\s+\\d+.*")) {
          val parts = trimmed.split("\\s+")
          if (parts.length >= 10) {
            // New format with depth, region, quality: Contig Position SNP Anc Der Call State Depth Region Quality
            snpDetails = snpDetails :+ SnpDetailRow(
              parts(0), parts(1).toInt, parts(2), parts(3), parts(4), parts(5), parts(6),
              readDepth = Some(parts(7)),
              region = Some(parts(8)),
              quality = Some(parts.drop(9).mkString(" "))
            )
          } else if (parts.length >= 7) {
            // Old format: Contig Position SNP Anc Der Call State
            snpDetails = snpDetails :+ SnpDetailRow(parts(0), parts(1).toInt, parts(2), parts(3), parts(4), parts(5), parts(6))
          }
        }
      }

      // Novel SNPs parsing - handles both old (5 cols) and new (7 cols with depth/region) formats
      if (currentSection == "novel_snps") {
        if (trimmed.startsWith("Contig")) inNovelSnpSection = true
        else if (inNovelSnpSection && trimmed.matches("^[A-Za-z0-9]+\\s+\\d+.*")) {
          val parts = trimmed.split("\\s+")
          if (parts.length >= 7) {
            // New format with depth, region: Contig Position Ref Alt Depth Region Quality
            privateSnps = privateSnps :+ PrivateVariantRow(
              parts(0), parts(1).toInt, parts(2), parts(3),
              quality = parts.drop(6).mkString(" "),
              strInfo = None,
              readDepth = Some(parts(4)),
              region = Some(parts(5))
            )
          } else if (parts.length >= 5) {
            // Old format: Contig Position Ref Alt Quality
            privateSnps = privateSnps :+ PrivateVariantRow(parts(0), parts(1).toInt, parts(2), parts(3), parts.drop(4).mkString(" "), None)
          }
        }
      }

      // Novel indels parsing - handles both old (5+ cols) and new (7+ cols with depth/region) formats
      if (currentSection == "novel_indels") {
        if (trimmed.startsWith("Contig")) inNovelIndelSection = true
        else if (inNovelIndelSection && trimmed.matches("^[A-Za-z0-9]+\\s+\\d+.*")) {
          val parts = trimmed.split("\\s+")
          if (parts.length >= 8) {
            // New format with depth, region, quality, STR: Contig Position Ref Alt Depth Region Quality STR...
            privateIndels = privateIndels :+ PrivateVariantRow(
              parts(0), parts(1).toInt, parts(2), parts(3),
              quality = parts(6),
              strInfo = if (parts.length > 7) Some(parts.drop(7).mkString(" ")) else None,
              readDepth = Some(parts(4)),
              region = Some(parts(5))
            )
          } else if (parts.length >= 7) {
            // New format with depth, region, quality: Contig Position Ref Alt Depth Region Quality
            privateIndels = privateIndels :+ PrivateVariantRow(
              parts(0), parts(1).toInt, parts(2), parts(3),
              quality = parts.drop(6).mkString(" "),
              strInfo = None,
              readDepth = Some(parts(4)),
              region = Some(parts(5))
            )
          } else if (parts.length >= 5) {
            // Old format: Contig Position Ref Alt Quality [STR...]
            val strInfo = if (parts.length > 5) Some(parts.drop(5).mkString(" ")) else None
            privateIndels = privateIndels :+ PrivateVariantRow(parts(0), parts(1).toInt, parts(2), parts(3), parts(4), strInfo)
          }
        }
      }
    }

    ParsedReport(
      metadata = metadata,
      topCandidate = candidates.headOption,
      candidates = candidates,
      lineagePath = lineagePath,
      snpDetails = snpDetails,
      privateSnps = privateSnps,
      privateIndels = privateIndels
    )
  }
}
