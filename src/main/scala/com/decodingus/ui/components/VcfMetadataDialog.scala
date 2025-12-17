package com.decodingus.ui.components

import com.decodingus.genotype.model.{TestTypeDefinition, TestTypes}
import com.decodingus.i18n.I18n.t
import com.decodingus.workspace.model.FileInfo
import scalafx.Includes.*
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.{GridPane, VBox}

/**
 * Metadata for a VCF-only import where no BAM/CRAM is available.
 */
case class VcfMetadata(
  testType: TestTypeDefinition,
  platform: String,
  notes: Option[String]
)

/**
 * Dialog for collecting metadata when importing a VCF file without an existing sequence run.
 * Common for Big Y, Veritas WGS, and other services that deliver VCF as the primary format.
 */
class VcfMetadataDialog(fileInfo: FileInfo) extends Dialog[Option[VcfMetadata]] {

  title = t("data.vcf_metadata.title")
  headerText = t("data.vcf_metadata.header")

  val addButtonType = new ButtonType(t("action.import"), ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(addButtonType, ButtonType.Cancel)

  // Test types that commonly deliver VCF as primary format
  private val vcfTestTypes: Seq[TestTypeDefinition] = Seq(
    TestTypes.BIG_Y_700,
    TestTypes.BIG_Y_500,
    TestTypes.Y_PRIME,
    TestTypes.Y_ELITE,
    TestTypes.WGS,
    TestTypes.WGS_LOW_PASS,
    TestTypes.WGS_HIFI,
    TestTypes.WGS_NANOPORE,
    TestTypes.MT_FULL_SEQUENCE,
    TestTypes.MT_PLUS
  )

  // Test type selection
  private val testTypeCombo = new ComboBox[TestTypeDefinition] {
    items = ObservableBuffer.from(vcfTestTypes)
    cellFactory = (lv: scalafx.scene.control.ListView[TestTypeDefinition]) => new ListCell[TestTypeDefinition] {
      item.onChange { (_, _, newType) =>
        text = Option(newType).map(_.displayName).getOrElse("")
      }
    }
    buttonCell = new ListCell[TestTypeDefinition] {
      item.onChange { (_, _, newType) =>
        text = Option(newType).map(_.displayName).getOrElse("")
      }
    }
    prefWidth = 300

    // Try to auto-detect from filename
    val fileName = fileInfo.fileName.toLowerCase
    val autoDetected = if (fileName.contains("bigy") || fileName.contains("big_y") || fileName.contains("big-y")) {
      Some(TestTypes.BIG_Y_700)
    } else if (fileName.contains("yprime") || fileName.contains("y_prime") || fileName.contains("yseq")) {
      Some(TestTypes.Y_PRIME)
    } else if (fileName.contains("yelite") || fileName.contains("y_elite") || fileName.contains("fullgenomes")) {
      Some(TestTypes.Y_ELITE)
    } else if (fileName.contains("veritas") || fileName.contains("dante") || fileName.contains("nebula")) {
      Some(TestTypes.WGS)
    } else if (fileName.contains("mtdna") || fileName.contains("mt_")) {
      Some(TestTypes.MT_FULL_SEQUENCE)
    } else {
      None
    }

    autoDetected match {
      case Some(detected) => selectionModel.value.select(detected)
      case None => selectionModel.value.selectFirst()
    }
  }

  // Platform field (auto-fills based on test type)
  private val platformField = new TextField {
    prefWidth = 300
    promptText = t("data.vcf_metadata.platform_hint")
  }

  // Auto-fill platform when test type changes
  testTypeCombo.selectionModel.value.selectedItemProperty.onChange { (_, _, newType) =>
    if (newType != null) {
      newType.vendor match {
        case Some("FamilyTreeDNA") => platformField.text = "Illumina"
        case Some("YSEQ") => platformField.text = "Illumina"
        case Some("Full Genomes") => platformField.text = "Illumina"
        case Some("PacBio") => platformField.text = "PacBio"
        case Some("Oxford Nanopore") => platformField.text = "Nanopore"
        case _ =>
          if (platformField.text.value.isEmpty) platformField.text = "Illumina"
      }
    }
  }
  // Trigger initial auto-fill
  Option(testTypeCombo.selectionModel.value.getSelectedItem).foreach { t =>
    t.vendor match {
      case Some("FamilyTreeDNA") => platformField.text = "Illumina"
      case _ => platformField.text = "Illumina"
    }
  }

  // Notes field
  private val notesField = new TextArea {
    prefWidth = 300
    prefHeight = 60
    promptText = t("data.vcf_metadata.notes_hint")
    wrapText = true
  }

  // File info display
  private val fileInfoLabel = new Label(fileInfo.fileName) {
    style = "-fx-font-weight: bold;"
  }

  // Layout
  private val grid = new GridPane {
    hgap = 10
    vgap = 10
    padding = Insets(20)

    add(new Label(t("data.file") + ":"), 0, 0)
    add(fileInfoLabel, 1, 0)

    add(new Label(t("data.vcf_metadata.test_type") + ":"), 0, 1)
    add(testTypeCombo, 1, 1)

    add(new Label(t("data.vcf_metadata.platform") + ":"), 0, 2)
    add(platformField, 1, 2)

    add(new Label(t("data.vcf_metadata.notes") + ":"), 0, 3)
    add(notesField, 1, 3)
  }

  private val helpText = new Label(t("data.vcf_metadata.help")) {
    wrapText = true
    prefWidth = 380
    style = "-fx-text-fill: #666666; -fx-font-size: 11px;"
  }

  dialogPane().content = new VBox(15) {
    padding = Insets(10)
    children = Seq(grid, helpText)
  }

  dialogPane().setPrefWidth(450)

  // Validation
  private val addButton = dialogPane().lookupButton(addButtonType)

  private def updateAddButton(): Unit = {
    val isValid = testTypeCombo.selectionModel.value.getSelectedItem != null &&
                  platformField.text.value.nonEmpty
    addButton.disable = !isValid
  }

  platformField.text.onChange { (_, _, _) => updateAddButton() }
  testTypeCombo.selectionModel.value.selectedItemProperty.onChange { (_, _, _) => updateAddButton() }
  updateAddButton()

  // Result converter
  resultConverter = dialogButton => {
    if (dialogButton == addButtonType) {
      val selectedType = testTypeCombo.selectionModel.value.getSelectedItem
      val platform = platformField.text.value.trim
      val notes = Option(notesField.text.value).map(_.trim).filter(_.nonEmpty)

      if (selectedType != null && platform.nonEmpty) {
        Some(VcfMetadata(selectedType, platform, notes))
      } else {
        None
      }
    } else {
      None
    }
  }
}
