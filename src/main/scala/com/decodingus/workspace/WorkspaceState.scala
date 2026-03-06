package com.decodingus.workspace

import com.decodingus.workspace.model.Workspace

/**
 * Immutable wrapper for workspace state, used by services for functional updates.
 *
 * Services receive this state, perform transformations, and return updated state.
 * The ViewModel is responsible for persisting changes and notifying observers.
 */
case class WorkspaceState(workspace: Workspace)

object WorkspaceState {
  def empty: WorkspaceState = WorkspaceState(Workspace.empty)
}
