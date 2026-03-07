package com.decodingus.ui.components

import com.decodingus.i18n.I18n.t
import com.decodingus.util.{DetectedFileType, FileTypeDetector}
import scalafx.Includes.*
import scalafx.beans.property.StringProperty
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.{HBox, Priority, Region, VBox}
import scalafx.stage.FileChooser

import java.io.File

/**
 * Result from the batch import dialog — list of files with their detected types.
 */
case class BatchImportEntry(
  file: File,
  detectedType: DetectedFileType,
  subjectName: String
)

/**
 * Dialog for selecting and importing multiple data files at once.
 * Auto-detects file types and assigns subject names based on filenames.
 */
class BatchImportDialog extends Dialog[Option[List[BatchImportEntry]]] {

  title = t("batch.title")
  headerText = t("batch.select_files")
  resizable = true
  dialogPane().setPrefSize(700, 500)

  private val importButtonType = new ButtonType(t("batch.import_all"), ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(importButtonType, ButtonType.Cancel)

  private val fileEntries = ObservableBuffer.empty[BatchFileRow]

  private val tableView = new TableView[BatchFileRow](fileEntries) {
    columnResizePolicy = TableView.ConstrainedResizePolicy
    placeholder = new Label(t("batch.no_files")) {
      style = "-fx-text-fill: #888888;"
    }

    columns ++= Seq(
      new TableColumn[BatchFileRow, String] {
        text = t("batch.filename")
        prefWidth = 200
        cellValueFactory = r => StringProperty(r.value.fileName)
      },
      new TableColumn[BatchFileRow, String] {
        text = t("batch.type")
        prefWidth = 120
        cellValueFactory = r => StringProperty(r.value.detectedType)
      },
      new TableColumn[BatchFileRow, String] {
        text = t("batch.subject")
        prefWidth = 150
        cellValueFactory = r => StringProperty(r.value.subjectName)
        cellFactory = { (_: TableColumn[BatchFileRow, String]) =>
          new TableCell[BatchFileRow, String] {
            item.onChange { (_, _, newValue) =>
              if (newValue != null) {
                val tf = new TextField {
                  text = newValue
                  prefWidth = 140
                }
                tf.text.onChange { (_, _, updated) =>
                  Option(tableRow.value).flatMap(tr => Option(tr.item.value)).foreach { row =>
                    row.subjectName = updated
                  }
                }
                graphic = tf
              } else {
                graphic = null
              }
            }
          }
        }
      },
      new TableColumn[BatchFileRow, String] {
        text = t("batch.size")
        prefWidth = 80
        cellValueFactory = r => StringProperty(r.value.fileSize)
      }
    )
  }
  VBox.setVgrow(tableView, Priority.Always)

  private val addFilesButton = new Button(t("batch.add_files")) {
    onAction = _ => browseFiles()
  }

  private val removeButton = new Button(t("batch.remove")) {
    disable = true
    onAction = _ => {
      val selected = tableView.selectionModel.value.getSelectedItem
      if (selected != null) fileEntries -= selected
      updateImportButton()
    }
  }

  private val clearButton = new Button(t("batch.clear")) {
    onAction = _ => {
      fileEntries.clear()
      updateImportButton()
    }
  }

  private val countLabel = new Label {
    style = "-fx-text-fill: #888888; -fx-font-size: 12px;"
  }

  tableView.selectionModel.value.selectedItemProperty.onChange { (_, _, sel) =>
    removeButton.disable = sel == null
  }

  fileEntries.onChange { (_, _) =>
    countLabel.text = s"${fileEntries.size} file${if (fileEntries.size != 1) "s" else ""}"
    updateImportButton()
  }

  private val buttonBar = new HBox(10) {
    alignment = Pos.CenterLeft
    children = Seq(
      addFilesButton,
      removeButton,
      clearButton,
      new Region { HBox.setHgrow(this, Priority.Always) },
      countLabel
    )
  }

  dialogPane().content = new VBox(10) {
    padding = Insets(10)
    children = Seq(buttonBar, tableView)
  }

  updateImportButton()

  resultConverter = dialogButton => {
    if (dialogButton == importButtonType && fileEntries.nonEmpty) {
      Some(fileEntries.toList.map(row =>
        BatchImportEntry(row.file, row.fileType, row.subjectName)
      ))
    } else {
      None
    }
  }

  private def browseFiles(): Unit = {
    val fileChooser = new FileChooser {
      this.title = t("batch.select_files")
      extensionFilters.addAll(
        new FileChooser.ExtensionFilter("All Supported Files", Seq("*.bam", "*.cram", "*.vcf", "*.vcf.gz", "*.fasta", "*.fa", "*.fna", "*.csv", "*.tsv", "*.txt")),
        new FileChooser.ExtensionFilter("Alignment Files", Seq("*.bam", "*.cram")),
        new FileChooser.ExtensionFilter("VCF Files", Seq("*.vcf", "*.vcf.gz")),
        new FileChooser.ExtensionFilter("FASTA Files", Seq("*.fasta", "*.fa", "*.fna")),
        new FileChooser.ExtensionFilter("All Files", "*.*")
      )
    }
    val javaFiles = fileChooser.delegate.showOpenMultipleDialog(dialogPane().getScene.getWindow)
    if (javaFiles != null) {
      import scala.jdk.CollectionConverters.*
      val existingPaths = fileEntries.map(_.file.getAbsolutePath).toSet
      javaFiles.asScala.foreach { file =>
        if (!existingPaths.contains(file.getAbsolutePath)) {
          val detected = FileTypeDetector.detect(file)
          val subjectName = deriveSubjectName(file)
          fileEntries += BatchFileRow(file, detected, subjectName)
        }
      }
    }
  }

  private def deriveSubjectName(file: File): String = {
    val name = file.getName
    // Strip common extensions to get a clean subject name
    val stripped = name
      .replaceAll("\\.(bam|cram|vcf|vcf\\.gz|fasta|fa|fna|csv|tsv|txt)$", "")
      .replaceAll("[._]", " ")
      .trim
    if (stripped.nonEmpty) stripped else name
  }

  private def formatFileSize(file: File): String = {
    val size = file.length()
    if (size < 1024) s"${size}B"
    else if (size < 1024 * 1024) f"${size / 1024.0}%.0fKB"
    else if (size < 1024L * 1024 * 1024) f"${size / (1024.0 * 1024)}%.1fMB"
    else f"${size / (1024.0 * 1024 * 1024)}%.1fGB"
  }

  private def updateImportButton(): Unit = {
    dialogPane().lookupButton(importButtonType).disable = fileEntries.isEmpty
  }

  private case class BatchFileRow(
    file: File,
    fileType: DetectedFileType,
    var subjectName: String
  ) {
    def fileName: String = file.getName
    def detectedType: String = fileType.description
    def fileSize: String = formatFileSize(file)
  }
}
