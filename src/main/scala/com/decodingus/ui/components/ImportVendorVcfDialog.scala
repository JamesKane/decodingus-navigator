package com.decodingus.ui.components

import com.decodingus.analysis.VcfVendor
import scalafx.Includes.*
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.control.*
import scalafx.scene.layout.{GridPane, HBox, Priority, VBox}
import scalafx.stage.FileChooser

import java.io.File
import java.nio.file.Path

/**
 * Result from the import vendor VCF dialog.
 */
case class VendorVcfImportRequest(
                                   vcfPath: Path,
                                   bedPath: Option[Path],
                                   vendor: VcfVendor,
                                   referenceBuild: String,
                                   notes: Option[String]
                                 )

/**
 * Dialog for importing vendor-provided VCF files (e.g., FTDNA Big Y).
 * Allows selection of:
 * - VCF file (required)
 * - Target regions BED file (optional)
 * - Vendor type
 * - Reference genome build
 */
class ImportVendorVcfDialog extends Dialog[Option[VendorVcfImportRequest]] {

  title = "Import Vendor VCF"
  headerText = "Import vendor-provided VCF and target regions"

  val importButtonType = new ButtonType("Import", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(importButtonType, ButtonType.Cancel)

  // State
  private var selectedVcfFile: Option[File] = None
  private var selectedBedFile: Option[File] = None

  // VCF file selection
  private val vcfPathField = new TextField {
    editable = false
    promptText = "Select VCF file..."
    prefWidth = 350
  }

  private val vcfBrowseButton = new Button("Browse...") {
    onAction = _ => browseForVcf()
  }

  // BED file selection (optional)
  private val bedPathField = new TextField {
    editable = false
    promptText = "(Optional) Select target regions BED..."
    prefWidth = 350
  }

  private val bedBrowseButton = new Button("Browse...") {
    onAction = _ => browseForBed()
  }

  private val bedClearButton = new Button("Clear") {
    onAction = _ => {
      selectedBedFile = None
      bedPathField.text = ""
    }
  }

  // Vendor selection
  private val vendorCombo = new ComboBox[String] {
    items = ObservableBuffer(
      "FTDNA Big Y",
      "FTDNA mtFull Sequence",
      "YSEQ",
      "Nebula Genomics",
      "Dante Labs",
      "Full Genomes Corp",
      "Other"
    )
    value = "FTDNA Big Y"
    prefWidth = 200
  }

  // Reference build selection
  private val refBuildCombo = new ComboBox[String] {
    items = ObservableBuffer(
      "GRCh38",
      "GRCh37",
      "T2T-CHM13v2.0",
      "hs1"
    )
    value = "GRCh38"
    prefWidth = 200
  }

  // Notes field
  private val notesField = new TextField {
    promptText = "(Optional) Notes about this import..."
    prefWidth = 350
  }

  // Layout
  private val grid = new GridPane {
    hgap = 10
    vgap = 15
    padding = Insets(20)

    // Row 0: VCF file
    add(new Label("VCF File:") {
      style = "-fx-font-weight: bold;"
    }, 0, 0)
    add(new HBox(10) {
      children = Seq(vcfPathField, vcfBrowseButton)
      hgrow = Priority.Always
    }, 1, 0)

    // Row 1: BED file
    add(new Label("Target Regions:"), 0, 1)
    add(new HBox(10) {
      children = Seq(bedPathField, bedBrowseButton, bedClearButton)
      hgrow = Priority.Always
    }, 1, 1)

    // Row 2: Vendor
    add(new Label("Vendor:"), 0, 2)
    add(vendorCombo, 1, 2)

    // Row 3: Reference build
    add(new Label("Reference Build:"), 0, 3)
    add(refBuildCombo, 1, 3)

    // Row 4: Notes
    add(new Label("Notes:"), 0, 4)
    add(notesField, 1, 4)
  }

  // Help text
  private val helpLabel = new Label(
    "Import a VCF file from a vendor test (e.g., FTDNA Big Y).\n" +
      "The target regions BED file is optional and helps with coverage analysis."
  ) {
    style = "-fx-text-fill: #666666; -fx-font-size: 11px;"
    wrapText = true
  }

  dialogPane().content = new VBox(15) {
    padding = Insets(10)
    children = Seq(helpLabel, grid)
  }

  // Disable import button until VCF is selected
  private val importButton = dialogPane().lookupButton(importButtonType)
  importButton.disable = true

  // Result converter
  resultConverter = dialogButton => {
    if (dialogButton == importButtonType && selectedVcfFile.isDefined) {
      val vendor = vendorCombo.value.value match {
        case "FTDNA Big Y" => VcfVendor.FTDNA_BIGY
        case "FTDNA mtFull Sequence" => VcfVendor.FTDNA_MTFULL
        case "YSEQ" => VcfVendor.YSEQ
        case "Nebula Genomics" => VcfVendor.NEBULA
        case "Dante Labs" => VcfVendor.DANTE
        case "Full Genomes Corp" => VcfVendor.FULL_GENOMES
        case _ => VcfVendor.OTHER
      }

      Some(VendorVcfImportRequest(
        vcfPath = selectedVcfFile.get.toPath,
        bedPath = selectedBedFile.map(_.toPath),
        vendor = vendor,
        referenceBuild = refBuildCombo.value.value,
        notes = Option(notesField.text.value).filter(_.nonEmpty)
      ))
    } else {
      None
    }
  }

  private def browseForVcf(): Unit = {
    val fileChooser = new FileChooser {
      title = "Select VCF File"
      extensionFilters.addAll(
        new FileChooser.ExtensionFilter("VCF Files", Seq("*.vcf", "*.vcf.gz")),
        new FileChooser.ExtensionFilter("All Files", Seq("*.*"))
      )
    }

    Option(fileChooser.showOpenDialog(dialogPane().getScene.getWindow)).foreach { file =>
      selectedVcfFile = Some(file)
      vcfPathField.text = file.getName
      importButton.disable = false

      // Auto-detect reference build from filename if possible
      val filename = file.getName.toLowerCase
      if (filename.contains("grch37") || filename.contains("hg19") || filename.contains("b37")) {
        refBuildCombo.value = "GRCh37"
      } else if (filename.contains("chm13") || filename.contains("t2t")) {
        refBuildCombo.value = "T2T-CHM13v2.0"
      } else if (filename.contains("hs1")) {
        refBuildCombo.value = "hs1"
      }
    }
  }

  private def browseForBed(): Unit = {
    val fileChooser = new FileChooser {
      title = "Select Target Regions BED File"
      extensionFilters.addAll(
        new FileChooser.ExtensionFilter("BED Files", Seq("*.bed", "*.bed.gz")),
        new FileChooser.ExtensionFilter("All Files", Seq("*.*"))
      )
    }

    Option(fileChooser.showOpenDialog(dialogPane().getScene.getWindow)).foreach { file =>
      selectedBedFile = Some(file)
      bedPathField.text = file.getName
    }
  }
}
