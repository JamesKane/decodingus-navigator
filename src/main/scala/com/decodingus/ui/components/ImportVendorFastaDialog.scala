package com.decodingus.ui.components

import scalafx.Includes._
import scalafx.scene.control.{ButtonType, Dialog, Label, TextField, ButtonBar, Button, ComboBox, Alert}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.layout.{VBox, HBox, GridPane, Priority}
import scalafx.geometry.{Insets, Pos}
import scalafx.stage.FileChooser
import scalafx.collections.ObservableBuffer
import com.decodingus.analysis.VcfVendor

import java.io.File
import java.nio.file.Path

/**
 * Result from the import vendor FASTA dialog.
 */
case class VendorFastaImportRequest(
  fastaPath: Path,
  vendor: VcfVendor,
  notes: Option[String]
)

/**
 * Dialog for importing vendor-provided mtDNA FASTA files (e.g., FTDNA mtFull Sequence, YSEQ mtDNA).
 * These files contain the full mtDNA sequence aligned to rCRS.
 */
class ImportVendorFastaDialog extends Dialog[Option[VendorFastaImportRequest]] {

  title = "Import mtDNA FASTA"
  headerText = "Import vendor-provided mtDNA FASTA sequence"

  val importButtonType = new ButtonType("Import", ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(importButtonType, ButtonType.Cancel)

  // State
  private var selectedFastaFile: Option[File] = None

  // FASTA file selection
  private val fastaPathField = new TextField {
    editable = false
    promptText = "Select FASTA file..."
    prefWidth = 350
  }

  private val fastaBrowseButton = new Button("Browse...") {
    onAction = _ => browseForFasta()
  }

  // Vendor selection - only mtDNA-capable vendors
  private val vendorCombo = new ComboBox[String] {
    items = ObservableBuffer(
      "FTDNA mtFull Sequence",
      "YSEQ",
      "Nebula Genomics",
      "Dante Labs",
      "Full Genomes Corp",
      "Other"
    )
    value = "FTDNA mtFull Sequence"
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

    // Row 0: FASTA file
    add(new Label("FASTA File:") { style = "-fx-font-weight: bold;" }, 0, 0)
    add(new HBox(10) {
      children = Seq(fastaPathField, fastaBrowseButton)
      hgrow = Priority.Always
    }, 1, 0)

    // Row 1: Vendor
    add(new Label("Vendor:"), 0, 1)
    add(vendorCombo, 1, 1)

    // Row 2: Notes
    add(new Label("Notes:"), 0, 2)
    add(notesField, 1, 2)
  }

  // Help text
  private val helpLabel = new Label(
    "Import an mtDNA FASTA file from a vendor test (e.g., FTDNA mtFull Sequence).\n" +
    "The sequence will be compared against rCRS to identify variants for haplogroup analysis.\n\n" +
    "Supported formats:\n" +
    "- Standard FASTA (.fa, .fasta, .fna)\n" +
    "- Full mtDNA sequence (~16,569 bp)"
  ) {
    style = "-fx-text-fill: #666666; -fx-font-size: 11px;"
    wrapText = true
  }

  dialogPane().content = new VBox(15) {
    padding = Insets(10)
    children = Seq(helpLabel, grid)
  }

  // Disable import button until FASTA is selected
  private val importButton = dialogPane().lookupButton(importButtonType)
  importButton.disable = true

  // Result converter
  resultConverter = dialogButton => {
    if (dialogButton == importButtonType && selectedFastaFile.isDefined) {
      val vendor = vendorCombo.value.value match {
        case "FTDNA mtFull Sequence" => VcfVendor.FTDNA_MTFULL
        case "YSEQ" => VcfVendor.YSEQ
        case "Nebula Genomics" => VcfVendor.NEBULA
        case "Dante Labs" => VcfVendor.DANTE
        case "Full Genomes Corp" => VcfVendor.FULL_GENOMES
        case _ => VcfVendor.OTHER
      }

      Some(VendorFastaImportRequest(
        fastaPath = selectedFastaFile.get.toPath,
        vendor = vendor,
        notes = Option(notesField.text.value).filter(_.nonEmpty)
      ))
    } else {
      None
    }
  }

  private def browseForFasta(): Unit = {
    val fileChooser = new FileChooser {
      title = "Select mtDNA FASTA File"
      extensionFilters.addAll(
        new FileChooser.ExtensionFilter("FASTA Files", Seq("*.fa", "*.fasta", "*.fna", "*.fas")),
        new FileChooser.ExtensionFilter("All Files", Seq("*.*"))
      )
    }

    Option(fileChooser.showOpenDialog(dialogPane().getScene.getWindow)).foreach { file =>
      // Validate the file appears to be mtDNA FASTA
      validateFastaFile(file) match {
        case Right(_) =>
          selectedFastaFile = Some(file)
          fastaPathField.text = file.getName
          importButton.disable = false

        case Left(error) =>
          new Alert(AlertType.Warning) {
            initOwner(dialogPane().getScene.getWindow)
            title = "Invalid FASTA File"
            headerText = "The selected file may not be a valid mtDNA FASTA"
            contentText = error + "\n\nDo you want to import it anyway?"
            buttonTypes = Seq(ButtonType.Yes, ButtonType.No)
          }.showAndWait() match {
            case Some(ButtonType.Yes) =>
              selectedFastaFile = Some(file)
              fastaPathField.text = file.getName
              importButton.disable = false
            case _ =>
              // User declined, keep current state
          }
      }
    }
  }

  /**
   * Basic validation of FASTA file format and length.
   */
  private def validateFastaFile(file: File): Either[String, Unit] = {
    try {
      val source = scala.io.Source.fromFile(file)
      try {
        val lines = source.getLines().toList
        if (lines.isEmpty) {
          return Left("File is empty")
        }

        // Check for FASTA header
        val hasHeader = lines.head.startsWith(">")
        if (!hasHeader) {
          return Left("File does not start with a FASTA header (>)")
        }

        // Calculate sequence length
        val sequence = lines.tail.filterNot(_.startsWith(">")).map(_.trim.toUpperCase).mkString
        val seqLength = sequence.length

        // mtDNA should be approximately 16,569 bp
        if (seqLength < 16000) {
          return Left(s"Sequence is too short ($seqLength bp). Expected ~16,569 bp for mtDNA.")
        }
        if (seqLength > 17000) {
          return Left(s"Sequence is too long ($seqLength bp). Expected ~16,569 bp for mtDNA.")
        }

        // Check for valid nucleotides
        val validBases = Set('A', 'C', 'G', 'T', 'N')
        val invalidBases = sequence.filterNot(validBases.contains)
        if (invalidBases.nonEmpty) {
          val sample = invalidBases.take(5).mkString
          return Left(s"Sequence contains invalid characters: $sample...")
        }

        Right(())
      } finally {
        source.close()
      }
    } catch {
      case e: Exception =>
        Left(s"Error reading file: ${e.getMessage}")
    }
  }
}
