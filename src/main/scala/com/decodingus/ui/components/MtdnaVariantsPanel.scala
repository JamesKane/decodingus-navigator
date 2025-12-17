package com.decodingus.ui.components

import com.decodingus.i18n.I18n.t
import com.decodingus.workspace.model.{HaplogroupResult, VariantCall}
import scalafx.Includes.*
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.*

/**
 * Panel showing mtDNA variants from rCRS (revised Cambridge Reference Sequence).
 * Displays a table of all mutations with position, reference, alternate, and region.
 */
class MtdnaVariantsPanel extends VBox {

  spacing = 10
  padding = Insets(15)
  style = "-fx-background-color: #2a2a2a; -fx-background-radius: 10;"

  private val tableData = ObservableBuffer.empty[MtdnaVariant]
  private var variantCount: Int = 0

  // Header with title and count
  private val titleLabel = new Label(t("mtdna.variants_from_rcrs")) {
    style = "-fx-font-size: 14px; -fx-font-weight: bold; -fx-text-fill: #ffffff;"
  }

  private val countLabel = new Label {
    style = "-fx-font-size: 12px; -fx-text-fill: #888888;"
  }

  private val headerBox = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(
      titleLabel,
      countLabel,
      new Region { hgrow = Priority.Always },
      new Button(t("mtdna.export_fasta")) {
        style = "-fx-font-size: 11px;"
        disable = true
        tooltip = Tooltip(t("mtdna.export_coming_soon"))
      }
    )
  }

  // Variants table
  private val variantsTable = new TableView[MtdnaVariant](tableData) {
    prefHeight = 200
    maxHeight = 300
    columnResizePolicy = TableView.ConstrainedResizePolicy
    style = "-fx-background-color: #333333; -fx-border-color: #444444;"

    columns ++= Seq(
      new TableColumn[MtdnaVariant, String] {
        text = t("mtdna.position")
        prefWidth = 80
        cellValueFactory = { p =>
          StringProperty(p.value.positionDisplay)
        }
        cellFactory = { (_: TableColumn[MtdnaVariant, String]) =>
          new TableCell[MtdnaVariant, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                style = "-fx-text-fill: #ffffff; -fx-font-family: monospace; -fx-font-size: 12px;"
              } else {
                text = ""
              }
            }
          }
        }
      },
      new TableColumn[MtdnaVariant, String] {
        text = t("mtdna.rcrs")
        prefWidth = 60
        cellValueFactory = { p =>
          StringProperty(p.value.reference)
        }
        cellFactory = { (_: TableColumn[MtdnaVariant, String]) =>
          new TableCell[MtdnaVariant, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                style = "-fx-text-fill: #888888; -fx-font-family: monospace; -fx-font-size: 12px;"
              } else {
                text = ""
              }
            }
          }
        }
      },
      new TableColumn[MtdnaVariant, String] {
        text = t("mtdna.sample")
        prefWidth = 60
        cellValueFactory = { p =>
          StringProperty(p.value.alternate)
        }
        cellFactory = { (_: TableColumn[MtdnaVariant, String]) =>
          new TableCell[MtdnaVariant, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                style = "-fx-text-fill: #60a5fa; -fx-font-weight: bold; -fx-font-family: monospace; -fx-font-size: 12px;"
              } else {
                text = ""
              }
            }
          }
        }
      },
      new TableColumn[MtdnaVariant, String] {
        text = t("mtdna.region")
        prefWidth = 120
        cellValueFactory = { p =>
          StringProperty(p.value.regionDisplay)
        }
        cellFactory = { (_: TableColumn[MtdnaVariant, String]) =>
          new TableCell[MtdnaVariant, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                val color = regionColor(newValue)
                style = s"-fx-text-fill: $color; -fx-font-size: 11px;"
              } else {
                text = ""
              }
            }
          }
        }
      },
      new TableColumn[MtdnaVariant, String] {
        text = t("mtdna.type")
        prefWidth = 80
        cellValueFactory = { p =>
          StringProperty(p.value.mutationType)
        }
        cellFactory = { (_: TableColumn[MtdnaVariant, String]) =>
          new TableCell[MtdnaVariant, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                text = newValue
                val color = if (newValue == "SNP") "#b0b0b0" else "#FF9800"
                style = s"-fx-text-fill: $color; -fx-font-size: 11px;"
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
      new Label(t("mtdna.no_variants")) {
        style = "-fx-text-fill: #888888; -fx-font-size: 12px;"
      }
    )
    visible = false
    managed = false
  }

  children = Seq(headerBox, variantsTable, noDataPlaceholder)

  /**
   * Update the panel with mtDNA haplogroup result containing variants.
   */
  def setMtdnaResult(result: Option[HaplogroupResult]): Unit = {
    result.flatMap(_.privateVariants) match {
      case Some(pvd) if pvd.variants.nonEmpty =>
        // Convert VariantCall to MtdnaVariant with region classification
        val variants = pvd.variants
          .filter(isMtdnaVariant)
          .map(toMtdnaVariant)
          .sortBy(_.position)

        variantCount = variants.size
        countLabel.text = s"(${t("mtdna.mutations_count", variantCount.toString)})"

        tableData.clear()
        tableData ++= variants

        variantsTable.visible = true
        variantsTable.managed = true
        noDataPlaceholder.visible = false
        noDataPlaceholder.managed = false
        this.visible = true
        this.managed = true

      case _ =>
        // No variants - hide panel
        this.visible = false
        this.managed = false
    }
  }

  private def isMtdnaVariant(v: VariantCall): Boolean = {
    // mtDNA contig is NC_012920 (rCRS) or chrM
    v.contigAccession.contains("NC_012920") ||
      v.contigAccession.equalsIgnoreCase("chrM") ||
      v.contigAccession.equalsIgnoreCase("MT")
  }

  private def toMtdnaVariant(v: VariantCall): MtdnaVariant = {
    val region = classifyMtdnaRegion(v.position)
    val mutType = if (v.referenceAllele.length == 1 && v.alternateAllele.length == 1) "SNP"
    else if (v.referenceAllele == "-" || v.alternateAllele.length > v.referenceAllele.length) "Insertion"
    else "Deletion"

    MtdnaVariant(
      position = v.position,
      reference = if (v.referenceAllele == "-") "-" else v.referenceAllele,
      alternate = if (v.alternateAllele == "-") "-" else v.alternateAllele,
      region = region,
      mutationType = mutType
    )
  }

  private def classifyMtdnaRegion(position: Int): MtdnaRegion = {
    // mtDNA regions based on rCRS coordinates
    position match {
      case p if p >= 16024 && p <= 16569 => MtdnaRegion.HVS1
      case p if p >= 1 && p <= 576 => MtdnaRegion.HVS2
      case p if p >= 438 && p <= 574 => MtdnaRegion.HVS3 // Overlaps with HVS2
      case _ => MtdnaRegion.Coding
    }
  }

  private def regionColor(region: String): String = region match {
    case s if s.contains("HVS1") => "#FF9800"
    case s if s.contains("HVS2") => "#2196F3"
    case s if s.contains("HVS3") => "#9C27B0"
    case _ => "#888888"
  }

  /**
   * Represents a single mtDNA variant for display.
   */
  private case class MtdnaVariant(
    position: Int,
    reference: String,
    alternate: String,
    region: MtdnaRegion,
    mutationType: String
  ) {
    def positionDisplay: String = position.toString

    def regionDisplay: String = region match {
      case MtdnaRegion.HVS1 => "HVS1 (Control)"
      case MtdnaRegion.HVS2 => "HVS2 (Control)"
      case MtdnaRegion.HVS3 => "HVS3 (Control)"
      case MtdnaRegion.Coding => "Coding"
    }
  }

  /**
   * mtDNA region classification.
   */
  private enum MtdnaRegion:
    case HVS1, HVS2, HVS3, Coding
}

object MtdnaVariantsPanel {
  def apply(): MtdnaVariantsPanel = new MtdnaVariantsPanel()
}
