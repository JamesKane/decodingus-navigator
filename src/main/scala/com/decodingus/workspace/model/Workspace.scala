package com.decodingus.workspace.model

case class WorkspaceContent(
  samples: List[Biosample],
  projects: List[Project]
)

case class Workspace(
  lexicon: Int,
  id: String,
  main: WorkspaceContent
)
