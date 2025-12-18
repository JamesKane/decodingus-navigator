package com.decodingus.ui.v2

import com.decodingus.i18n.I18n
import com.decodingus.i18n.I18n.{t, bind}
import com.decodingus.workspace.WorkbenchViewModel
import scalafx.Includes.*
import scalafx.geometry.{Insets, Pos, Side}
import scalafx.scene.control.*
import scalafx.scene.layout.*

/**
 * Main application shell with tab-based navigation.
 *
 * Provides three top-level tabs:
 * - Dashboard: Workspace overview and pending work
 * - Subjects: Subject management with grid and detail views
 * - Projects: Project management
 *
 * This replaces the previous SplitPane-based WorkbenchView.
 */
class MainShell(viewModel: WorkbenchViewModel) extends BorderPane {

  // Initialize i18n from preferences at shell creation
  I18n.initializeFromPreferences()

  // ============================================================================
  // Tab Content Views
  // ============================================================================

  private val dashboardView = new DashboardView(viewModel)
  private val subjectsView = new SubjectsView(viewModel)
  private val projectsView = new ProjectsView(viewModel)

  // ============================================================================
  // Tab Navigation
  // ============================================================================

  private val dashboardTab = createTab("nav.dashboard", "\uD83D\uDCCA", dashboardView)
  private val subjectsTab = createTab("nav.subjects", "\uD83D\uDC65", subjectsView)
  private val projectsTab = createTab("nav.projects", "\uD83D\uDCC1", projectsView)

  private val tabPane = new TabPane {
    tabClosingPolicy = TabPane.TabClosingPolicy.Unavailable
    side = Side.Top
    styleClass += "main-tab-pane"
    tabs = Seq(dashboardTab, subjectsTab, projectsTab)

    // Track tab selection for any needed side effects
    selectionModel.value.selectedItemProperty.onChange { (_, _, newTab) =>
      if (newTab != null) {
        onTabSelected(newTab)
      }
    }
  }

  // ============================================================================
  // Layout
  // ============================================================================

  center = tabPane
  styleClass += "main-shell"

  // ============================================================================
  // Helper Methods
  // ============================================================================

  /**
   * Creates a tab with an icon and i18n-bound text.
   */
  private def createTab(i18nKey: String, icon: String, tabContent: javafx.scene.Node): Tab = {
    val tab = new Tab {
      // Bind text to i18n key with icon prefix
      text <== bind(i18nKey).map(s => s"$icon $s")
      closable = false
      styleClass += "main-tab"
    }
    tab.content = tabContent
    tab
  }

  /**
   * Called when a tab is selected.
   * Can be used for lazy loading or refresh logic.
   */
  private def onTabSelected(tab: javafx.scene.control.Tab): Unit = {
    // Refresh data when switching tabs if needed
    tab match {
      case t if t == dashboardTab.delegate => dashboardView.refresh()
      case t if t == subjectsTab.delegate => subjectsView.refresh()
      case t if t == projectsTab.delegate => projectsView.refresh()
      case _ => // Unknown tab
    }
  }

  // ============================================================================
  // Public API
  // ============================================================================

  /**
   * Navigate to the Subjects tab and optionally select a subject.
   */
  def navigateToSubject(subjectId: String): Unit = {
    tabPane.selectionModel.value.select(subjectsTab)
    subjectsView.selectSubject(subjectId)
  }

  /**
   * Navigate to the Projects tab and optionally select a project.
   */
  def navigateToProject(projectId: String): Unit = {
    tabPane.selectionModel.value.select(projectsTab)
    projectsView.selectProject(projectId)
  }

  /**
   * Navigate to the Dashboard tab.
   */
  def navigateToDashboard(): Unit = {
    tabPane.selectionModel.value.select(dashboardTab)
  }
}
