package com.decodingus.ui.components

import com.decodingus.config.LabsConfig
import com.decodingus.genotype.model.TestTypes
import com.decodingus.i18n.I18n.t
import com.decodingus.ui.theme.Theme
import com.decodingus.workspace.model.SequenceRun
import scalafx.Includes.*
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.*
import scalafx.scene.layout.{GridPane, HBox, Priority, VBox}

/**
 * Result of editing a sequence run - contains the updated values.
 */
case class SequenceRunEditResult(
  testType: String,
  sequencingFacility: Option[String]
)

/**
 * Dialog for editing editable fields of a SequenceRun.
 *
 * Allows editing:
 * - Test Type (WGS, Y_ELITE, BIG_Y_700, etc.)
 * - Lab (sequencing facility name)
 *
 * @param seqRun The sequence run to edit
 */
class EditSequenceRunDialog(seqRun: SequenceRun) extends Dialog[Option[SequenceRunEditResult]] {

  title = t("dialog.edit_sequence_run.title")
  headerText = t("dialog.edit_sequence_run.header")

  val saveButtonType = new ButtonType(t("common.save"), ButtonBar.ButtonData.OKDone)
  dialogPane().buttonTypes = Seq(saveButtonType, ButtonType.Cancel)

  // Get all known test types for the dropdown
  private val testTypeOptions = TestTypes.all.map { defn =>
    (defn.code, defn.displayName)
  }.sortBy(_._2)

  // Test type combo box with display name rendering
  private val testTypeCombo = new ComboBox[String]() {
    items = ObservableBuffer.from(testTypeOptions.map(_._1))
    value = seqRun.testType
    prefWidth = 280
    // Use converter for display
    converter = new javafx.util.StringConverter[String] {
      override def toString(code: String): String =
        Option(code).flatMap(c => testTypeOptions.find(_._1 == c).map(_._2)).getOrElse("")
      override def fromString(string: String): String =
        testTypeOptions.find(_._2 == string).map(_._1).getOrElse(string)
    }
  }

  // Auto-detected indicator
  private val autoDetectedLabel = new Label(t("dialog.edit_sequence_run.auto_detected")) {
    style = s"-fx-font-size: 11px; -fx-text-fill: ${Theme.current.textMuted};"
  }

  // Lab/sequencing facility combo box (editable for custom entry)
  private val labOptions = ObservableBuffer.from(
    "" +: LabsConfig.sequenceRunLabNames :+ "Other..."
  )
  private val labComboBox = new ComboBox[String]() {
    items = labOptions
    editable = true
    // Set current value - find matching lab name or use existing value
    val currentLab = seqRun.sequencingFacility.getOrElse("")
    value = if (currentLab.isEmpty) ""
            else LabsConfig.findLab(currentLab).map(_.displayName).getOrElse(currentLab)
    promptText = t("dialog.edit_sequence_run.lab_placeholder")
    prefWidth = 280
  }

  // Handle "Other..." selection to clear for custom entry
  labComboBox.selectionModel().selectedItem.onChange { (_, _, newVal) =>
    if (newVal == "Other...") {
      labComboBox.editor.value.clear()
      labComboBox.editor.value.requestFocus()
    }
  }

  // Build the form
  private val formGrid = new GridPane {
    hgap = 15
    vgap = 10
    padding = Insets(20)

    add(new Label(t("dialog.edit_sequence_run.test_type") + ":") {
      style = "-fx-font-weight: bold;"
    }, 0, 0)
    add(new HBox(10) {
      alignment = Pos.CenterLeft
      children = Seq(testTypeCombo, autoDetectedLabel)
    }, 1, 0)

    add(new Label(t("dialog.edit_sequence_run.lab") + ":") {
      style = "-fx-font-weight: bold;"
    }, 0, 1)
    add(labComboBox, 1, 1)

    // Info about current detection (read-only)
    add(new Label(t("dialog.edit_sequence_run.platform") + ":"), 0, 2)
    add(new Label(seqRun.platformName), 1, 2)

    add(new Label(t("dialog.edit_sequence_run.instrument") + ":"), 0, 3)
    add(new Label(seqRun.instrumentModel.getOrElse(t("common.unknown"))), 1, 3)

    if (seqRun.totalReads.isDefined) {
      add(new Label(t("dialog.edit_sequence_run.total_reads") + ":"), 0, 4)
      add(new Label(seqRun.totalReads.map(r => f"$r%,d").getOrElse(t("common.unknown"))), 1, 4)
    }
  }

  // Help text
  private val helpText = new Label(t("dialog.edit_sequence_run.help_text")) {
    style = s"-fx-font-size: 11px; -fx-text-fill: ${Theme.current.textMuted};"
    wrapText = true
    maxWidth = 350
  }

  dialogPane().content = new VBox(15) {
    padding = Insets(10)
    children = Seq(formGrid, helpText)
  }

  // Enable save button only when test type is selected
  val saveButton = dialogPane().lookupButton(saveButtonType)
  saveButton.disable <== testTypeCombo.value.isNull

  // Convert result
  resultConverter = dialogButton => {
    if (dialogButton == saveButtonType) {
      // Get value from editor (for editable combo) or selection
      val rawValue = Option(labComboBox.editor.value.getText).getOrElse(labComboBox.value.value)
      val labValue = Option(rawValue).map(_.trim).filter(v => v.nonEmpty && v != "Other...")
      Some(SequenceRunEditResult(
        testType = testTypeCombo.value.value,
        sequencingFacility = labValue
      ))
    } else {
      None
    }
  }
}
