package com.decodingus.ui.components

import com.decodingus.workspace.WorkbenchViewModel
import com.decodingus.workspace.model.*
import scalafx.Includes.*
import scalafx.application.Platform
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Pos}
import scalafx.scene.control.Alert.AlertType
import scalafx.scene.control.ControlIncludes.*
import scalafx.scene.control.*
import scalafx.scene.input.{ClipboardContent, DragEvent, MouseEvent, TransferMode}
import scalafx.scene.layout.*

import java.util.{Timer, TimerTask, UUID}

class WorkbenchView(val viewModel: WorkbenchViewModel) extends SplitPane {
  println(s"[DEBUG] WorkbenchView: Initializing WorkbenchView. ViewModel Projects: ${viewModel.projects.size}, ViewModel Samples: ${viewModel.samples.size}")

  // Track drag state to prevent click-on-drag from triggering navigation
  private var dragInProgress = false
  // Timer for delayed selection check
  private val selectionTimer = new Timer("SelectionTimer", true) // daemon thread
  private var pendingSelectionTask: Option[TimerTask] = None
  private val clickDelayMs = 150 // milliseconds to wait before applying selection

  // Use the shared DataFormat from ProjectDetailView companion object
  private val biosampleFormat = ProjectDetailView.biosampleFormat

  // Observable buffers for UI lists - now using filtered versions from ViewModel
  private val projectBuffer: ObservableBuffer[Project] = viewModel.filteredProjects
  private val sampleBuffer: ObservableBuffer[Biosample] = viewModel.filteredSamples

  println(s"[DEBUG] WorkbenchView: After binding buffers. projectBuffer size: ${projectBuffer.size}, sampleBuffer size: ${sampleBuffer.size}")

  // Detail view for right panel
  private val detailView = new VBox(10) {
    padding = Insets(10)
    children = Seq(
      new Label("Select an item to view details") {
        style = "-fx-font-size: 18px; -fx-font-weight: bold;"
      }
    )
  }
  VBox.setVgrow(detailView, Priority.Always)

  /** Renders the detail view for a selected subject with Edit/Delete actions */
  private def renderSubjectDetail(subject: Biosample): Unit = {
    detailView.children.clear()

    // === HEADER ROW: Title + Action Buttons ===
    val headerRow = new HBox(10) {
      alignment = Pos.CenterLeft
      children = Seq(
        new Label(s"Subject: ${subject.donorIdentifier}") {
          style = "-fx-font-size: 20px; -fx-font-weight: bold;"
        },
        new Region { hgrow = Priority.Always },
        new Button("Edit") {
          onAction = _ => handleEditSubject(subject)
        },
        new Button("Delete") {
          style = "-fx-text-fill: #D32F2F;"
          onAction = _ => handleDeleteSubject(subject)
        }
      )
    }

    // === INFO SECTION: Two-column grid for better horizontal space usage ===
    val infoGrid = new GridPane {
      hgap = 30
      vgap = 4
      padding = Insets(5, 0, 10, 0)
    }
    // Left column
    infoGrid.add(new Label(s"Accession: ${subject.sampleAccession}") { style = "-fx-font-size: 12px;" }, 0, 0)
    infoGrid.add(new Label(s"Sex: ${subject.sex.getOrElse("N/A")}") { style = "-fx-font-size: 12px;" }, 0, 1)
    infoGrid.add(new Label(s"Created: ${subject.meta.createdAt.toLocalDate}") { style = "-fx-font-size: 12px;" }, 0, 2)
    // Right column
    infoGrid.add(new Label(s"Center: ${subject.centerName.getOrElse("N/A")}") { style = "-fx-font-size: 12px;" }, 1, 0)
    infoGrid.add(new Label(s"Description: ${subject.description.getOrElse("N/A")}") { style = "-fx-font-size: 12px;" }, 1, 1)

    // === HAPLOGROUP SECTION: Combines Y Haplogroup with Profile management ===
    def formatHaplogroup(result: Option[com.decodingus.workspace.model.HaplogroupResult]): String = {
      result match {
        case Some(h) =>
          val name = h.haplogroupName
          val derived = h.matchingSnps.map(n => s"+$n").getOrElse("")
          val ancestral = h.ancestralMatches.map(n => s"-$n").getOrElse("")
          if (derived.nonEmpty || ancestral.nonEmpty) s"$name ($derived/$ancestral)" else name
        case None => "—"
      }
    }

    val yDnaReconciliation = viewModel.workspace.value.main.getYDnaReconciliation(subject)
    val mtDnaReconciliation = viewModel.workspace.value.main.getMtDnaReconciliation(subject)
    val biosampleId = viewModel.getBiosampleIdByAccession(subject.sampleAccession)

    // Y-DNA haplogroup box with integrated profile management
    val yHaplogroupBox = new HBox(8) {
      alignment = Pos.CenterLeft
      padding = Insets(8)
      style = "-fx-background-color: #2d3a2d; -fx-background-radius: 6;"

      val yHaplogroup = subject.haplogroups.flatMap(_.yDna)
      val yText = formatHaplogroup(yHaplogroup)

      children = Seq(
        new Label("Y:") { style = "-fx-font-size: 13px; -fx-font-weight: bold; -fx-text-fill: #888;" },
        new Label(yText) { style = "-fx-font-size: 14px; -fx-font-weight: bold; -fx-text-fill: #8f8;" },
        new Region { hgrow = Priority.Always }
      )

      // Add reconciliation indicator if available
      yDnaReconciliation.foreach { recon =>
        val (color, _) = recon.status.compatibilityLevel match {
          case CompatibilityLevel.COMPATIBLE => ("#4CAF50", "Compatible")
          case CompatibilityLevel.MINOR_DIVERGENCE => ("#FF9800", "Minor differences")
          case CompatibilityLevel.MAJOR_DIVERGENCE => ("#F44336", "Major divergence")
          case CompatibilityLevel.INCOMPATIBLE => ("#9C27B0", "Incompatible")
        }
        children.add(new Label(s"● ${recon.status.runCount}") {
          style = s"-fx-text-fill: $color; -fx-font-size: 11px;"
          tooltip = Tooltip("Runs contributing to haplogroup")
        })
      }

      // Y Profile management button (replaces separate Y Profile section)
      biosampleId.foreach { bsId =>
        val hasProfile = viewModel.getYProfileSummary(bsId).isDefined
        children.add(new Button(if (hasProfile) "Profile" else "+ Profile") {
          style = "-fx-font-size: 10px; -fx-padding: 2 6;"
          tooltip = Tooltip(if (hasProfile) "Manage Y Chromosome Profile" else "Create Y Chromosome Profile")
          onAction = _ => handleManageYProfile(subject, bsId)
        })
      }
    }
    HBox.setHgrow(yHaplogroupBox, Priority.Always)

    // MT-DNA haplogroup box
    val mtHaplogroupBox = new HBox(8) {
      alignment = Pos.CenterLeft
      padding = Insets(8)
      style = "-fx-background-color: #2d2d3a; -fx-background-radius: 6;"

      val mtHaplogroup = subject.haplogroups.flatMap(_.mtDna)
      val mtText = formatHaplogroup(mtHaplogroup)

      children = Seq(
        new Label("MT:") { style = "-fx-font-size: 13px; -fx-font-weight: bold; -fx-text-fill: #888;" },
        new Label(mtText) { style = "-fx-font-size: 14px; -fx-font-weight: bold; -fx-text-fill: #88f;" },
        new Region { hgrow = Priority.Always }
      )

      // Add reconciliation indicator if available
      mtDnaReconciliation.foreach { recon =>
        val (color, _) = recon.status.compatibilityLevel match {
          case CompatibilityLevel.COMPATIBLE => ("#4CAF50", "Compatible")
          case CompatibilityLevel.MINOR_DIVERGENCE => ("#FF9800", "Minor differences")
          case CompatibilityLevel.MAJOR_DIVERGENCE => ("#F44336", "Major divergence")
          case CompatibilityLevel.INCOMPATIBLE => ("#9C27B0", "Incompatible")
        }
        children.add(new Label(s"● ${recon.status.runCount}") {
          style = s"-fx-text-fill: $color; -fx-font-size: 11px;"
          tooltip = Tooltip("Runs contributing to haplogroup")
        })
      }
    }
    HBox.setHgrow(mtHaplogroupBox, Priority.Always)

    val haplogroupRow = new HBox(10) {
      padding = Insets(5, 0, 10, 0)
      children = Seq(yHaplogroupBox, mtHaplogroupBox)
    }

    // === DATA TABS: Sequencing, Chip/Array, STR ===
    val sequenceRuns = viewModel.workspace.value.main.getSequenceRunsForBiosample(subject)
    val allAlignments = viewModel.workspace.value.main.alignments
    val chipProfiles = viewModel.getChipProfilesForBiosample(subject.sampleAccession)
    val strProfiles = viewModel.getStrProfilesForBiosample(subject.sampleAccession)

    val sequenceTable = new SequenceDataTable(
      viewModel = viewModel,
      subject = subject,
      sequenceRuns = sequenceRuns,
      alignments = allAlignments,
      onAnalyze = (index: Int) => handleAnalyzeSequenceData(subject.sampleAccession, index),
      onRemove = (index: Int) => handleRemoveSequenceData(subject.sampleAccession, index)
    )

    val chipTable = new ChipDataTable(
      viewModel = viewModel,
      subject = subject,
      chipProfiles = chipProfiles,
      onRemove = (uri: String) => handleRemoveChipProfile(subject.sampleAccession, uri)
    )

    val strTable = new StrProfileTable(
      viewModel = viewModel,
      subject = subject,
      strProfiles = strProfiles,
      onRemove = (uri: String) => handleRemoveStrProfile(subject.sampleAccession, uri)
    )

    // Create tabs for data types
    val sequenceTab = new Tab {
      text = s"Sequencing (${sequenceRuns.size})"
      closable = false
      content = sequenceTable
    }

    val chipTab = new Tab {
      text = s"Chip/Array (${chipProfiles.size})"
      closable = false
      content = chipTable
    }

    val strTab = new Tab {
      text = s"STR Profiles (${strProfiles.size})"
      closable = false
      content = strTable
    }

    val dataTabPane = new TabPane {
      tabs = Seq(sequenceTab, chipTab, strTab)
      tabMinWidth = 100
    }
    VBox.setVgrow(dataTabPane, Priority.Always)

    // === ASSEMBLE THE VIEW ===
    detailView.children.addAll(
      headerRow,
      infoGrid,
      haplogroupRow,
      dataTabPane
    )
  }

  /** Shows the reconciliation detail dialog for a subject */
  private def showReconciliationDetails(
                                         subject: Biosample,
                                         yDnaReconciliation: Option[HaplogroupReconciliation],
                                         mtDnaReconciliation: Option[HaplogroupReconciliation]
                                       ): Unit = {
    val dialog = new ReconciliationDetailDialog(subject, yDnaReconciliation, mtDnaReconciliation)
    dialog.showAndWait()
  }

  /** Creates the Y Chromosome Profile summary section for a subject */
  private def createYProfileSection(subject: Biosample): Option[VBox] = {
    if (!viewModel.isYProfileAvailable) return None

    // Try to get biosample UUID from atUri
    val biosampleId = viewModel.getBiosampleIdByAccession(subject.sampleAccession)
    biosampleId.map { bsId =>
      val profileSummary = viewModel.getYProfileSummary(bsId)

      profileSummary match {
        case Some(summary) =>
          // Profile exists - show summary with View and Manage buttons
          new VBox(8) {
            padding = Insets(10, 0, 10, 0)
            style = "-fx-background-color: #f5f5f5; -fx-background-radius: 6; -fx-padding: 10;"

            val headerBox = new HBox(10) {
              alignment = Pos.CenterLeft
              children = Seq(
                new Label("Y Chromosome Profile") {
                  style = "-fx-font-size: 14px; -fx-font-weight: bold;"
                },
                new Region {
                  HBox.setHgrow(this, Priority.Always)
                },
                new Button("View Details") {
                  style = "-fx-font-size: 11px;"
                  onAction = _ => handleViewYProfile(subject, bsId)
                },
                new Button("Manage") {
                  style = "-fx-font-size: 11px;"
                  onAction = _ => handleManageYProfile(subject, bsId)
                }
              )
            }

            // Haplogroup display
            val haplogroupLabel = summary.consensusHaplogroup match {
              case Some(hg) =>
                val confidenceText = summary.haplogroupConfidence.map(c => f" (${c * 100}%.0f%%)").getOrElse("")
                new Label(s"$hg$confidenceText") {
                  style = "-fx-font-size: 16px; -fx-font-weight: bold; -fx-text-fill: #2d5a2d;"
                }
              case None =>
                new Label("Haplogroup pending") {
                  style = "-fx-font-size: 14px; -fx-text-fill: #666;"
                }
            }

            // Status badges
            val badgeBox = new HBox(8) {
              alignment = Pos.CenterLeft
              children = Seq(
                if (summary.confirmedCount > 0) Some(createBadge(s"${summary.confirmedCount} Confirmed", "#4CAF50")) else None,
                if (summary.novelCount > 0) Some(createBadge(s"${summary.novelCount} Novel", "#2196F3")) else None,
                if (summary.conflictCount > 0) Some(createBadge(s"${summary.conflictCount} Conflict", "#F44336")) else None,
                Some(createBadge(s"${summary.sourceCount} Source${if (summary.sourceCount != 1) "s" else ""}", "#9E9E9E"))
              ).flatten
            }

            // Callable region (if available)
            val callableLabel = summary.callableRegionPct.map { pct =>
              new Label(f"Callable: ${pct * 100}%.1f%%") {
                style = "-fx-font-size: 11px; -fx-text-fill: #666;"
              }
            }

            children = Seq(headerBox, haplogroupLabel, badgeBox) ++ callableLabel.toSeq
          }

        case None =>
          // No profile exists - show create button
          new VBox(8) {
            padding = Insets(10, 0, 10, 0)
            style = "-fx-background-color: #f0f0f0; -fx-background-radius: 6; -fx-padding: 10;"

            val headerBox = new HBox(10) {
              alignment = Pos.CenterLeft
              children = Seq(
                new Label("Y Chromosome Profile") {
                  style = "-fx-font-size: 14px; -fx-font-weight: bold;"
                },
                new Region {
                  HBox.setHgrow(this, Priority.Always)
                },
                new Button("Create Profile") {
                  style = "-fx-font-size: 11px;"
                  onAction = _ => handleManageYProfile(subject, bsId)
                }
              )
            }

            val descLabel = new Label("No Y-DNA profile exists. Create one to combine data from multiple Y-DNA tests.") {
              style = "-fx-font-size: 12px; -fx-text-fill: #666;"
              wrapText = true
            }

            children = Seq(headerBox, descLabel)
          }
      }
    }
  }

  /** Creates a colored badge label */
  private def createBadge(text: String, color: String): Label = {
    new Label(text) {
      style = s"-fx-background-color: $color; -fx-text-fill: white; -fx-padding: 2 6 2 6; -fx-background-radius: 3; -fx-font-size: 11px;"
    }
  }

  /** Handles opening the Y Profile detail dialog */
  private def handleViewYProfile(subject: Biosample, biosampleId: UUID): Unit = {
    // Show loading indicator
    val loadingAlert = new Alert(AlertType.Information) {
      title = "Loading"
      headerText = "Loading Y Profile..."
      contentText = "Please wait while the profile data is loaded."
      buttonTypes = Seq.empty // No buttons - auto-close when loaded
    }

    // Load data asynchronously
    viewModel.loadYProfileForBiosample(biosampleId, {
      case Right(data) =>
        Platform.runLater {
          loadingAlert.close()
          val dialog = new YProfileDetailDialog(
            data.profile,
            data.variants,
            data.sources,
            data.variantCalls,
            data.auditEntries,
            subject.donorIdentifier,
            data.yRegionAnnotator
          )
          dialog.showAndWait()
        }
      case Left(error) =>
        Platform.runLater {
          loadingAlert.close()
          new Alert(AlertType.Error) {
            title = "Error"
            headerText = "Could not load Y Profile"
            contentText = error
          }.showAndWait()
        }
    })

    loadingAlert.show()
  }

  /** Handles opening the Y Profile management dialog */
  private def handleManageYProfile(subject: Biosample, biosampleId: UUID): Unit = {
    // Load data for the management dialog
    viewModel.loadYProfileManagementData(biosampleId, {
      case Right(data) =>
        Platform.runLater {
          val dialog = new YProfileManagementDialog(
            biosampleId = biosampleId,
            biosampleName = subject.donorIdentifier,
            yProfileService = data.yProfileService,
            existingProfile = data.profile,
            sources = data.sources,
            variants = data.variants,
            snpPanels = data.snpPanels,
            onRefresh = () => {
              // Refresh the subject detail view after profile changes
              renderSubjectDetail(subject)
            }
          )
          dialog.showAndWait()
        }
      case Left(error) =>
        Platform.runLater {
          new Alert(AlertType.Error) {
            title = "Error"
            headerText = "Could not load Y Profile data"
            contentText = error
          }.showAndWait()
        }
    })
  }

  /** Handles triggering analysis for a sequence run */
  private def handleAnalyzeSequenceData(sampleAccession: String, index: Int): Unit = {
    viewModel.findSubject(sampleAccession) match {
      case Some(subject) =>
        val sequenceRuns = viewModel.workspace.value.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(index) match {
          case Some(sequenceRun) =>
            sequenceRun.files.headOption match {
              case Some(fileInfo) =>
                // Check if initial analysis has been run (has alignments)
                val runAlignments = viewModel.workspace.value.main.getAlignmentsForSequenceRun(sequenceRun)
                val hasAlignments = runAlignments.nonEmpty
                val hasMetrics = runAlignments.exists(_.metrics.isDefined)

                if (!hasAlignments) {
                  // Run initial analysis
                  showAnalysisChoiceDialog(sampleAccession, index, fileInfo.fileName, "initial")
                } else if (!hasMetrics) {
                  // Offer to run deep coverage analysis
                  showAnalysisChoiceDialog(sampleAccession, index, fileInfo.fileName, "wgs")
                } else {
                  // Both analyses complete - offer to re-run
                  showAnalysisChoiceDialog(sampleAccession, index, fileInfo.fileName, "both_complete")
                }
              case None =>
                new Alert(AlertType.Warning) {
                  title = "No File"
                  headerText = "No alignment file associated"
                  contentText = "Please add a BAM/CRAM file to this sequencing run."
                }.showAndWait()
            }
          case None =>
            println(s"[View] Sequence run not found at index $index")
        }
      case None =>
        println(s"[View] Subject not found: $sampleAccession")
    }
  }

  /** Shows a dialog to choose which analysis to run */
  private def showAnalysisChoiceDialog(sampleAccession: String, index: Int, fileName: String, state: String): Unit = {
    val (dialogHeader, dialogContent, options) = state match {
      case "initial" =>
        ("Run Initial Analysis",
          s"Analyze $fileName to detect platform, reference build, and collect library statistics.",
          Seq(("Run Initial Analysis", () => runInitialAnalysis(sampleAccession, index))))
      case "wgs" =>
        ("Run Deep Coverage Analysis",
          s"Initial analysis complete. Would you like to run WGS metrics analysis?\n\nThis will calculate detailed coverage statistics using GATK and may take several minutes for large genomes.",
          Seq(
            ("Run WGS Metrics", () => runWgsMetricsAnalysis(sampleAccession, index)),
            ("Re-run Initial Analysis", () => runInitialAnalysis(sampleAccession, index))
          ))
      case "both_complete" =>
        ("Analysis Complete",
          s"Both initial and WGS metrics analysis have been completed for $fileName.\n\nWould you like to re-run any analysis?",
          Seq(
            ("Re-run WGS Metrics", () => runWgsMetricsAnalysis(sampleAccession, index)),
            ("Re-run Initial Analysis", () => runInitialAnalysis(sampleAccession, index))
          ))
      case _ =>
        ("Analysis", "Choose an analysis to run.", Seq.empty)
    }

    if (options.size == 1) {
      // Single option - just confirm
      val confirm = new Alert(AlertType.Confirmation) {
        title = "Analysis"
        headerText = dialogHeader
        contentText = dialogContent
      }
      confirm.showAndWait() match {
        case Some(ButtonType.OK) => options.head._2()
        case _ =>
      }
    } else if (options.nonEmpty) {
      // Multiple options - use custom buttons
      val alert = new Alert(AlertType.Confirmation) {
        title = "Analysis Options"
        headerText = dialogHeader
        contentText = dialogContent
        buttonTypes = options.map(o => new ButtonType(o._1)) :+ ButtonType.Cancel
      }
      val result = alert.showAndWait()
      result.foreach { btn =>
        options.find(_._1 == btn.text).foreach(_._2())
      }
    }
  }

  /** Runs initial analysis with progress dialog */
  private def runInitialAnalysis(sampleAccession: String, index: Int): Unit = {
    val progressDialog = new AnalysisProgressDialog(
      "Initial Analysis",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runInitialAnalysis(sampleAccession, index, {
      case Right(libraryStats) =>
        Platform.runLater {
          // Calculate mean read length from distribution
          val meanReadLength = if (libraryStats.lengthDistribution.nonEmpty) {
            val totalReads = libraryStats.lengthDistribution.values.sum.toDouble
            val weightedSum = libraryStats.lengthDistribution.map { case (len, count) => len.toLong * count }.sum
            f"${weightedSum / totalReads}%.1f bp"
          } else "N/A"

          // Calculate mean insert size from distribution
          val meanInsertSize = if (libraryStats.insertSizeDistribution.nonEmpty) {
            val totalPairs = libraryStats.insertSizeDistribution.values.sum.toDouble
            val weightedSum = libraryStats.insertSizeDistribution.map { case (size, count) => size * count }.sum
            f"${weightedSum / totalPairs}%.1f bp"
          } else "N/A"

          new Alert(AlertType.Information) {
            title = "Analysis Complete"
            headerText = "Initial Analysis Results"
            contentText =
              s"""Sample: ${libraryStats.sampleName}
                 |Platform: ${libraryStats.inferredPlatform}
                 |Instrument: ${libraryStats.mostFrequentInstrument}
                 |Reference: ${libraryStats.referenceBuild}
                 |Aligner: ${libraryStats.aligner}
                 |Mean Read Length: $meanReadLength
                 |Mean Insert Size: $meanInsertSize""".stripMargin
          }.showAndWait()
          // Refresh the detail view
          viewModel.selectedSubject.value.foreach(renderSubjectDetail)
        }
      case Left(error) =>
        Platform.runLater {
          new Alert(AlertType.Error) {
            title = "Analysis Failed"
            headerText = "Initial analysis encountered an error"
            contentText = error
          }.showAndWait()
        }
    })

    progressDialog.show()
  }

  /** Runs WGS metrics analysis with progress dialog */
  private def runWgsMetricsAnalysis(sampleAccession: String, index: Int): Unit = {
    val progressDialog = new AnalysisProgressDialog(
      "WGS Metrics Analysis",
      viewModel.analysisProgress,
      viewModel.analysisProgressPercent,
      viewModel.analysisInProgress
    )

    viewModel.runWgsMetricsAnalysis(sampleAccession, index, {
      case Right(wgsMetrics) =>
        Platform.runLater {
          new Alert(AlertType.Information) {
            title = "Analysis Complete"
            headerText = "WGS Metrics Results"
            contentText =
              f"""Mean Coverage: ${wgsMetrics.meanCoverage}%.1fx
                 |Median Coverage: ${wgsMetrics.medianCoverage}%.1fx
                 |SD Coverage: ${wgsMetrics.sdCoverage}%.2f
                 |PCT 10x: ${wgsMetrics.pct10x * 100}%.1f%%
                 |PCT 20x: ${wgsMetrics.pct20x * 100}%.1f%%
                 |PCT 30x: ${wgsMetrics.pct30x * 100}%.1f%%
                 |Het SNP Sensitivity: ${wgsMetrics.hetSnpSensitivity}%.4f""".stripMargin
          }.showAndWait()
          // Refresh the detail view
          viewModel.selectedSubject.value.foreach(renderSubjectDetail)
        }
      case Left(error) =>
        Platform.runLater {
          new Alert(AlertType.Error) {
            title = "Analysis Failed"
            headerText = "WGS metrics analysis encountered an error"
            contentText = error
          }.showAndWait()
        }
    })

    progressDialog.show()
  }

  /** Handles removing a sequence data entry */
  private def handleRemoveSequenceData(sampleAccession: String, index: Int): Unit = {
    viewModel.removeSequenceData(sampleAccession, index)
  }

  /** Handles removing an STR profile */
  private def handleRemoveStrProfile(sampleAccession: String, profileUri: String): Unit = {
    viewModel.deleteStrProfile(sampleAccession, profileUri) match {
      case Right(_) =>
        // Refresh the detail view
        viewModel.selectedSubject.value.foreach(renderSubjectDetail)
      case Left(error) =>
        new Alert(AlertType.Error) {
          title = "Error"
          headerText = "Could not remove STR profile"
          contentText = error
        }.showAndWait()
    }
  }

  /** Handles removing a chip profile */
  private def handleRemoveChipProfile(sampleAccession: String, profileUri: String): Unit = {
    viewModel.deleteChipProfile(sampleAccession, profileUri) match {
      case Right(_) =>
        // Refresh the detail view
        viewModel.selectedSubject.value.foreach(renderSubjectDetail)
      case Left(error) =>
        new Alert(AlertType.Error) {
          title = "Error"
          headerText = "Could not remove chip profile"
          contentText = error
        }.showAndWait()
    }
  }

  /** Renders the detail view for a selected project */
  private def renderProjectDetail(project: Project): Unit = {
    detailView.children.clear()
    val projectDetailView = new ProjectDetailView(
      viewModel = viewModel,
      project = project,
      onEdit = handleEditProject,
      onDelete = handleDeleteProject
    )
    VBox.setVgrow(projectDetailView, Priority.Always)
    detailView.children.add(projectDetailView)
  }

  /** Handles the Edit Project action */
  private def handleEditProject(project: Project): Unit = {
    val dialog = new EditProjectDialog(project)
    val result = dialog.showAndWait().asInstanceOf[Option[Option[Project]]]

    result match {
      case Some(Some(updatedProject)) =>
        viewModel.updateProject(updatedProject)
        // Refresh the detail view with updated project
        viewModel.findProject(updatedProject.projectName).foreach(renderProjectDetail)
      case _ => // User cancelled
    }
  }

  /** Handles the Delete Project action with confirmation */
  private def handleDeleteProject(project: Project): Unit = {
    val confirmDialog = new Alert(AlertType.Confirmation) {
      title = "Delete Project"
      headerText = s"Delete ${project.projectName}?"
      contentText = "This action cannot be undone. The project will be removed but subjects will remain in the workspace."
    }

    val result = confirmDialog.showAndWait()
    result match {
      case Some(ButtonType.OK) =>
        viewModel.deleteProject(project.projectName)
        renderEmptyDetail("Select an item to view details")
      case _ => // User cancelled
    }
  }

  /** Renders the empty state when nothing is selected */
  private def renderEmptyDetail(message: String): Unit = {
    detailView.children.clear()
    detailView.children.add(
      new Label(message) {
        style = "-fx-font-size: 18px; -fx-font-weight: bold;"
      }
    )
  }

  /** Handles the Edit Subject action */
  private def handleEditSubject(subject: Biosample): Unit = {
    val dialog = new EditSubjectDialog(subject)
    val result = dialog.showAndWait().asInstanceOf[Option[Option[Biosample]]]

    result match {
      case Some(Some(updatedBiosample)) =>
        viewModel.updateSubject(updatedBiosample)
      case _ => // User cancelled
    }
  }

  /** Handles the Delete Subject action with confirmation */
  private def handleDeleteSubject(subject: Biosample): Unit = {
    val confirmDialog = new Alert(AlertType.Confirmation) {
      title = "Delete Subject"
      headerText = s"Delete ${subject.donorIdentifier}?"
      contentText = "This action cannot be undone. All associated sequence data and analysis results will be removed."
    }

    val result = confirmDialog.showAndWait()
    result match {
      case Some(ButtonType.OK) =>
        viewModel.deleteSubject(subject.sampleAccession)
      case _ => // User cancelled
    }
  }

  // Unified detail view rendering based on selection state
  // This prevents race conditions between project and subject selection
  private def updateDetailView(): Unit = {
    Platform.runLater {
      (viewModel.selectedProject.value, viewModel.selectedSubject.value) match {
        case (Some(project), _) =>
          // Project takes precedence when selected
          renderProjectDetail(project)
        case (None, Some(subject)) =>
          renderSubjectDetail(subject)
        case (None, None) =>
          renderEmptyDetail("Select an item to view details")
      }
    }
  }

  // Listen to ViewModel's selectedSubject changes to update detailView
  viewModel.selectedSubject.onChange { (_, _, _) =>
    updateDetailView()
  }

  // Listen to ViewModel's selectedProject changes to update detailView
  viewModel.selectedProject.onChange { (_, _, _) =>
    updateDetailView()
  }

  // Left Panel - Navigation
  private val projectList = new ListView[Project]() {
    items = projectBuffer // Explicitly set items
    vgrow = Priority.Always
    prefHeight = 200 // Initial height for projects
    cellFactory = { (v: ListView[Project]) =>
      new ListCell[Project] {
        item.onChange { (_, _, newProject) =>
          text = if (newProject != null) newProject.projectName else null
        }
      }
    }
  }
  // UI to ViewModel sync
  projectList.selectionModel().selectedItem.onChange { (_, _, newProject) =>
    if (newProject != null) {
      // Clear subject selection first (both UI and ViewModel) to avoid race conditions
      if (viewModel.selectedSubject.value.isDefined) {
        viewModel.selectedSubject.value = None
      }
      sampleList.selectionModel().clearSelection()
      // Then set project selection
      viewModel.selectedProject.value = Some(newProject)
    } else if (viewModel.selectedProject.value.isDefined && projectList.selectionModel().getSelectedItem == null) {
      // Clear ViewModel selection if UI selection is cleared manually
      viewModel.selectedProject.value = None
    }
  }
  // ViewModel to UI sync
  viewModel.selectedProject.onChange { (_, _, newViewModelProjectOpt) =>
    if (newViewModelProjectOpt.isDefined && projectList.selectionModel().getSelectedItem != newViewModelProjectOpt.getOrElse(null)) {
      projectList.selectionModel().select(newViewModelProjectOpt.get)
    } else if (newViewModelProjectOpt.isEmpty && projectList.selectionModel().getSelectedItem != null) {
      projectList.selectionModel().clearSelection()
    }
  }

  private val sampleList = new ListView[Biosample]() {
    items = sampleBuffer // Explicitly set items
    vgrow = Priority.Always
    cellFactory = { (v: ListView[Biosample]) =>
      new ListCell[Biosample] {
        item.onChange { (_, _, newBiosample) =>
          text = if (newBiosample != null) s"${newBiosample.donorIdentifier} (${newBiosample.sampleAccession.take(8)}...)" else null
        }

        // On mouse press, schedule delayed selection check
        onMousePressed = (_: MouseEvent) => {
          dragInProgress = false
          // Cancel any existing pending task
          pendingSelectionTask.foreach(_.cancel())

          Option(item.value).foreach { biosample =>
            val task = new TimerTask {
              override def run(): Unit = {
                // Check if drag started during the delay
                if (!dragInProgress) {
                  Platform.runLater {
                    // Clear project selection first
                    if (viewModel.selectedProject.value.isDefined) {
                      viewModel.selectedProject.value = None
                    }
                    projectList.selectionModel().clearSelection()
                    // Then set subject selection
                    viewModel.selectedSubject.value = Some(biosample)
                  }
                }
                pendingSelectionTask = None
              }
            }
            pendingSelectionTask = Some(task)
            selectionTimer.schedule(task, clickDelayMs)
          }
        }

        // Drag source - enable dragging subjects to project members lists
        onDragDetected = (event: MouseEvent) => {
          Option(item.value).foreach { biosample =>
            dragInProgress = true
            // Cancel pending selection
            pendingSelectionTask.foreach(_.cancel())
            pendingSelectionTask = None

            val db = startDragAndDrop(TransferMode.Move)
            val content = new ClipboardContent()
            content.put(biosampleFormat, biosample.sampleAccession)
            content.putString(biosample.sampleAccession)
            db.setContent(content)
            event.consume()
          }
        }

        // Reset drag flag when drag completes
        onDragDone = (_: DragEvent) => {
          dragInProgress = false
        }
      }
    }
  }
  // UI to ViewModel sync - only handle deselection
  sampleList.selectionModel().selectedItem.onChange { (_, _, newBiosample) =>
    if (newBiosample == null && viewModel.selectedSubject.value.isDefined) {
      // Clear ViewModel selection if UI selection is cleared
      viewModel.selectedSubject.value = None
    }
  }
  // ViewModel to UI sync
  viewModel.selectedSubject.onChange { (_, _, newViewModelSubjectOpt) =>
    if (newViewModelSubjectOpt.isDefined && sampleList.selectionModel().getSelectedItem != newViewModelSubjectOpt.getOrElse(null)) {
      sampleList.selectionModel().select(newViewModelSubjectOpt.get)
    } else if (newViewModelSubjectOpt.isEmpty && sampleList.selectionModel().getSelectedItem != null) {
      sampleList.selectionModel().clearSelection()
    }
  }

  private val newProjectButton = new Button("New Project") {
    onAction = _ => {
      val dialog = new AddProjectDialog()
      val result = dialog.showAndWait().asInstanceOf[Option[Option[Project]]]

      result match {
        case Some(Some(newProject)) =>
          viewModel.addProject(newProject)
        case _ => // User cancelled
      }
    }
  }

  // Filter controls
  private val projectFilterField = new TextField() {
    promptText = "Filter projects..."
    prefWidth = 150
  }
  projectFilterField.text.bindBidirectional(viewModel.projectFilter)

  private val subjectFilterField = new TextField() {
    promptText = "Filter subjects..."
    prefWidth = 150
  }
  subjectFilterField.text.bindBidirectional(viewModel.subjectFilter)

  private val addSampleButton = new Button("Add Subject") {
    onAction = _ => {
      val dialog = new AddSubjectDialog()
      val result = dialog.showAndWait().asInstanceOf[Option[Option[Biosample]]]

      result match {
        case Some(Some(newBiosample)) =>
          viewModel.addSubject(newBiosample) // Delegate to ViewModel
        case _ => // User cancelled or closed dialog
      }
    }
  }

  private val saveButton = new Button("Save Workspace") {
    styleClass.add("button-primary")
    onAction = _ => viewModel.saveWorkspace() // Delegate to ViewModel
  }

  // Note: Sync status is now displayed in the application's StatusBar (bottom)
  // This simplifies the left panel and provides a more standard UX

  private val leftPanel = new VBox(10) {
    padding = Insets(10)
    children = Seq(
      new HBox(10) {
        alignment = Pos.CenterLeft
        children = Seq(
          new Label("Projects:") {
            style = "-fx-font-weight: bold;"
          },
          new Region {
            HBox.setHgrow(this, Priority.Always)
          },
          projectFilterField
        )
      },
      projectList,
      newProjectButton,
      new HBox(10) {
        alignment = Pos.CenterLeft
        children = Seq(
          new Label("Subjects:") {
            style = "-fx-font-weight: bold;"
          },
          new Region {
            HBox.setHgrow(this, Priority.Always)
          },
          subjectFilterField
        )
      },
      sampleList,
      addSampleButton,
      saveButton
    )
  }
  SplitPane.setResizableWithParent(leftPanel, false) // Make left panel not resize with parent by default

  // Right Panel - Details/Content Area
  private val rightPanel = new VBox(10) {
    padding = Insets(10)
    children = Seq(detailView) // rightPanel now contains the dynamic detailView
  }
  VBox.setVgrow(rightPanel, Priority.Always) // Allow right panel to grow vertically

  // Set the items of the SplitPane
  items.addAll(leftPanel, rightPanel)
  dividerPositions = 0.25 // Initial divider position

  // Listen for reference download prompts from the ViewModel
  viewModel.pendingReferenceDownload.onChange { (_, _, request) =>
    request match {
      case viewModel.PendingDownload(build, url, sizeMB, onConfirm, onCancel) =>
        Platform.runLater {
          val dialog = new ReferenceDownloadPromptDialog(build, url, sizeMB)
          dialog.showAndWait() match {
            case Some(ReferenceDownloadPromptDialog.Result.Download) =>
              onConfirm()
            case Some(ReferenceDownloadPromptDialog.Result.Configure) =>
              // Open settings dialog
              val configDialog = new ReferenceConfigDialog()
              configDialog.showAndWait()
              onCancel() // Cancel the current operation - user can retry after configuring
            case _ =>
              onCancel()
          }
        }
      case viewModel.NoDownloadPending =>
      // Nothing to do
    }
  }
}