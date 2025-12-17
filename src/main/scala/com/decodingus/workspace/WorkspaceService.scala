package com.decodingus.workspace

import com.decodingus.workspace.model.*

/**
 * Workspace persistence service trait.
 *
 * Implementation: H2WorkspaceAdapter (backed by H2 database)
 */
trait WorkspaceService {
  def load(): Either[String, Workspace]

  def save(workspace: Workspace): Either[String, Unit]
}
