package com.decodingus.ui.v2

import com.decodingus.i18n.I18n.{t, bind}
import com.decodingus.i18n.Formatters
import com.decodingus.ui.v2.BiosampleExtensions.*
import com.decodingus.workspace.WorkbenchViewModel
import com.decodingus.workspace.model.{Biosample, Project}
import scalafx.Includes.*
import scalafx.beans.property.{ObjectProperty, StringProperty}
import scalafx.collections.ObservableBuffer
import scalafx.geometry.{Insets, Orientation, Pos}
import scalafx.scene.control.*
import scalafx.scene.input.{ClipboardContent, DataFormat, DragEvent, MouseEvent, TransferMode}
import scalafx.scene.layout.*

import java.time.format.DateTimeFormatter
import java.util.UUID

/**
 * Projects view with project list and member management.
 *
 * Features:
 * - Searchable project list
 * - Project detail with member list
 * - Drag-and-drop support for adding members
 */
class ProjectsView(viewModel: WorkbenchViewModel) extends SplitPane {

  orientation = Orientation.Horizontal
  dividerPositions = 0.35
  styleClass += "projects-view"

  // ============================================================================
  // State
  // ============================================================================

  private val searchText = StringProperty("")
  private val selectedProject: ObjectProperty[Option[Project]] = ObjectProperty(None)

  // ============================================================================
  // Left Panel: Project List
  // ============================================================================

  private val searchField = new TextField {
    promptText = t("projects.search")
    prefWidth = 200
    text.onChange { (_, _, newValue) =>
      searchText.value = newValue
      applyFilter()
    }
  }

  private val addProjectButton = new Button {
    text = t("projects.add")
    styleClass += "button-primary"
    onAction = _ => handleAddProject()
  }

  private val searchBar = new HBox(10) {
    alignment = Pos.CenterLeft
    padding = Insets(10)
    children = Seq(
      searchField,
      new Region { hgrow = Priority.Always },
      addProjectButton
    )
  }

  private val projectListView: ListView[Project] = {
    val lv = new ListView[Project] {
      vgrow = Priority.Always
      placeholder = new Label(t("projects.no_projects"))
      styleClass += "project-list"

      selectionModel.value.selectedItemProperty.onChange { (_, _, newProject) =>
        if (newProject != null) {
          selectedProject.value = Some(newProject)
        } else {
          selectedProject.value = None
        }
      }
    }

    lv.cellFactory = { (_: ListView[Project]) =>
      new ListCell[Project] {
        item.onChange { (_, _, project) =>
          if (project != null) {
            graphic = createProjectListItem(project)
            text = null
          } else {
            graphic = null
            text = null
          }
        }
      }
    }

    lv
  }

  private val leftPanel = new VBox {
    vgrow = Priority.Always
    children = Seq(searchBar, projectListView)
  }

  // ============================================================================
  // Right Panel: Project Detail
  // ============================================================================

  private val projectNameLabel = new Label {
    style = "-fx-font-size: 20px; -fx-font-weight: bold;"
  }

  private val projectMetaLabel = new Label {
    style = "-fx-font-size: 12px; -fx-text-fill: #888888;"
  }

  private val editProjectButton = new Button {
    text = t("action.edit")
    onAction = _ => handleEditProject()
  }

  private val deleteProjectButton = new Button {
    text = t("action.delete")
    styleClass += "button-danger"
    onAction = _ => handleDeleteProject()
  }

  private val projectHeader = new HBox(15) {
    alignment = Pos.CenterLeft
    padding = Insets(15)
    style = "-fx-background-color: #2a2a2a;"
    children = Seq(
      new VBox(5) {
        children = Seq(projectNameLabel, projectMetaLabel)
      },
      new Region { hgrow = Priority.Always },
      editProjectButton,
      deleteProjectButton
    )
  }

  // Member list
  private val memberTableView = new TableView[Biosample] {
    vgrow = Priority.Always
    placeholder = new Label(t("projects.no_members"))
    styleClass += "member-table"

    columns ++= Seq(
      createColumn[String]("column.name", 150, s => s.donorId.getOrElse(s.accession)),
      createColumn[String]("column.id", 100, _.accession),
      createColumn[String]("column.ydna", 100, s => s.yHaplogroup.getOrElse("-")),
      createColumn[String]("column.mtdna", 80, s => s.mtHaplogroup.getOrElse("-"))
    )
  }

  private val addMemberButton = new Button {
    text = t("projects.add_member")
    onAction = _ => handleAddMember()
  }

  private val removeMemberButton = new Button {
    text = t("projects.remove_member")
    disable = true
    onAction = _ => handleRemoveMember()
  }

  private val memberActionsBar = new HBox(10) {
    alignment = Pos.CenterLeft
    padding = Insets(10)
    children = Seq(
      addMemberButton,
      removeMemberButton,
      new Region { hgrow = Priority.Always }
    )
  }

  // Drag-drop zone
  private val dropZone = new VBox {
    alignment = Pos.Center
    padding = Insets(20)
    style = "-fx-border-style: dashed; -fx-border-color: #555555; -fx-border-radius: 10; -fx-background-color: #252525; -fx-background-radius: 10;"
    children = Seq(
      new Label(t("projects.drag_hint")) {
        style = "-fx-text-fill: #666666;"
      }
    )

    onDragOver = (event: DragEvent) => {
      if (event.dragboard.hasString) {
        event.acceptTransferModes(TransferMode.Copy)
        style = "-fx-border-style: dashed; -fx-border-color: #4a9eff; -fx-border-radius: 10; -fx-background-color: #2a3a4a; -fx-background-radius: 10;"
      }
      event.consume()
    }

    onDragExited = (_: DragEvent) => {
      style = "-fx-border-style: dashed; -fx-border-color: #555555; -fx-border-radius: 10; -fx-background-color: #252525; -fx-background-radius: 10;"
    }

    onDragDropped = (event: DragEvent) => {
      val db = event.dragboard
      if (db.hasString) {
        val accession = db.getString
        handleDroppedMember(accession)
        event.dropCompleted = true
      }
      event.consume()
    }
  }

  private val memberSection = new VBox(10) {
    vgrow = Priority.Always
    padding = Insets(15)
    children = Seq(
      new Label { text <== bind("projects.members"); style = "-fx-font-weight: bold;" },
      memberTableView,
      memberActionsBar,
      dropZone
    )
  }

  private val emptyProjectPane = new VBox {
    alignment = Pos.Center
    children = Seq(
      new Label(t("info.no_data")) {
        style = "-fx-font-size: 16px; -fx-text-fill: #666666;"
      }
    )
  }

  private val projectDetailPane = new VBox {
    vgrow = Priority.Always
    visible = false
    children = Seq(projectHeader, memberSection)
  }

  private val rightPanel = new StackPane {
    children = Seq(emptyProjectPane, projectDetailPane)
  }

  // ============================================================================
  // Split Pane Setup
  // ============================================================================

  items.addAll(leftPanel, rightPanel)

  // ============================================================================
  // Data Binding
  // ============================================================================

  viewModel.projects.onChange { (_, _) =>
    applyFilter()
  }

  selectedProject.onChange { (_, _, newProject) =>
    newProject match {
      case Some(project) =>
        updateProjectDetail(project)
        projectDetailPane.visible = true
        emptyProjectPane.visible = false
      case None =>
        projectDetailPane.visible = false
        emptyProjectPane.visible = true
    }
  }

  // Enable/disable remove button based on selection
  memberTableView.selectionModel.value.selectedItemProperty.onChange { (_, _, selected) =>
    removeMemberButton.disable = selected == null
  }

  // Initial load
  applyFilter()

  // ============================================================================
  // Helper Methods
  // ============================================================================

  private def createProjectListItem(project: Project): HBox = {
    val memberCount = project.memberAccessions.size

    new HBox(10) {
      alignment = Pos.CenterLeft
      padding = Insets(10)
      children = Seq(
        new VBox(3) {
          hgrow = Priority.Always
          children = Seq(
            new Label(project.name) {
              style = "-fx-font-weight: bold;"
            },
            new Label(t("projects.member_count", memberCount)) {
              style = "-fx-font-size: 11px; -fx-text-fill: #888888;"
            }
          )
        }
      )
    }
  }

  private def createColumn[T](
    headerKey: String,
    colWidth: Double,
    valueExtractor: Biosample => T
  ): TableColumn[Biosample, T] = {
    new TableColumn[Biosample, T] {
      text <== bind(headerKey)
      prefWidth = colWidth
      cellValueFactory = { cellData =>
        ObjectProperty(valueExtractor(cellData.value))
      }
    }
  }

  private def applyFilter(): Unit = {
    val query = searchText.value.toLowerCase.trim
    val filtered = if (query.isEmpty) {
      viewModel.projects.toSeq
    } else {
      viewModel.projects.filter { p =>
        p.name.toLowerCase.contains(query) ||
          p.description.exists(_.toLowerCase.contains(query))
      }.toSeq
    }

    projectListView.items = ObservableBuffer.from(filtered)
  }

  private def updateProjectDetail(project: Project): Unit = {
    projectNameLabel.text = project.name

    val createdStr = project.createdAt.map { dt =>
      Formatters.formatDate(dt.toLocalDate)
    }.getOrElse("-")

    projectMetaLabel.text = s"${t("projects.created")}: $createdStr • ${t("projects.members")}: ${project.memberAccessions.size}"

    // Load members
    val members = project.memberAccessions.flatMap { accession =>
      viewModel.samples.find(_.accession == accession)
    }

    memberTableView.items = ObservableBuffer.from(members)
  }

  // ============================================================================
  // Action Handlers
  // ============================================================================

  private def handleAddProject(): Unit = {
    // TODO: Open AddProjectDialog
    println("[ProjectsView] Add project")
  }

  private def handleEditProject(): Unit = {
    selectedProject.value.foreach { project =>
      // TODO: Open EditProjectDialog
      println(s"[ProjectsView] Edit project: ${project.name}")
    }
  }

  private def handleDeleteProject(): Unit = {
    selectedProject.value.foreach { project =>
      // TODO: Show confirmation and delete
      println(s"[ProjectsView] Delete project: ${project.name}")
    }
  }

  private def handleAddMember(): Unit = {
    selectedProject.value.foreach { project =>
      // TODO: Open member picker dialog
      println(s"[ProjectsView] Add member to project: ${project.name}")
    }
  }

  private def handleRemoveMember(): Unit = {
    val selected = memberTableView.selectionModel.value.getSelectedItem
    if (selected != null) {
      selectedProject.value.foreach { project =>
        // TODO: Remove member from project
        println(s"[ProjectsView] Remove ${selected.accession} from ${project.name}")
      }
    }
  }

  private def handleDroppedMember(accession: String): Unit = {
    selectedProject.value.foreach { project =>
      viewModel.samples.find(_.accession == accession).foreach { biosample =>
        // TODO: Add member to project via ViewModel
        println(s"[ProjectsView] Add ${accession} to ${project.name} via drag-drop")
      }
    }
  }

  // ============================================================================
  // Public API
  // ============================================================================

  /**
   * Refresh the projects list.
   */
  def refresh(): Unit = {
    applyFilter()
    // Refresh member list if a project is selected
    selectedProject.value.foreach(updateProjectDetail)
  }

  /**
   * Select a project by ID.
   */
  def selectProject(projectIdStr: String): Unit = {
    viewModel.projects.find(_.projectId == projectIdStr).foreach { project =>
      projectListView.selectionModel.value.select(project)
      selectedProject.value = Some(project)
    }
  }
}
