package com.decodingus.workspace

import com.decodingus.service.{H2WorkspaceService, DatabaseContext}
import com.decodingus.workspace.model.{Workspace, WorkspaceContent, RecordMeta}

/**
 * Adapter that implements the legacy WorkspaceService interface using H2 backend.
 *
 * This bridges the old load/save API with the new H2 database layer,
 * allowing gradual migration from JSON to H2 persistence.
 *
 * Key behaviors:
 * - `load()` reads from H2 database, returns Workspace
 * - `save()` writes granularly to H2 (individual entity operations)
 * - Maintains backwards compatibility with existing ViewModel
 */
class H2WorkspaceAdapter(h2Service: H2WorkspaceService) extends WorkspaceService:

  /**
   * Load workspace content from H2 database.
   *
   * Converts the WorkspaceContent from H2 to the legacy Workspace format.
   */
  override def load(): Either[String, Workspace] =
    h2Service.loadWorkspaceContent().map { content =>
      Workspace(
        lexicon = Workspace.CurrentLexiconVersion,
        id = Workspace.NamespaceId,
        main = content
      )
    }

  /**
   * Save workspace to H2 database.
   *
   * This performs incremental updates rather than wholesale replacement:
   * - New samples/projects are created
   * - Existing samples/projects are updated
   * - Removed samples/projects are deleted
   *
   * Note: Currently performs a simplified save that updates existing entities.
   * A full implementation would diff the old and new state.
   */
  override def save(workspace: Workspace): Either[String, Unit] =
    // For each entity type, sync with H2
    // This is a simplified implementation - a production version would:
    // 1. Load current H2 state
    // 2. Diff with incoming workspace
    // 3. Apply create/update/delete operations

    val result = for
      // Sync biosamples first (sequence runs depend on them)
      _ <- syncBiosamples(workspace.main.samples)
      // Sync sequence runs (alignments depend on them)
      _ <- syncSequenceRuns(workspace.main.sequenceRuns)
      // Sync alignments
      _ <- syncAlignments(workspace.main.alignments)
      // Sync projects
      _ <- syncProjects(workspace.main.projects)
    yield ()

    result

  /**
   * Sync biosamples from workspace to H2.
   * Creates new biosamples, updates existing ones.
   */
  private def syncBiosamples(samples: List[model.Biosample]): Either[String, Unit] =
    samples.foldLeft[Either[String, Unit]](Right(())) { (acc, sample) =>
      acc.flatMap { _ =>
        // Check if biosample exists by accession
        h2Service.getBiosampleByAccession(sample.sampleAccession).flatMap {
          case Some(_) =>
            // Update existing
            h2Service.updateBiosample(sample).map(_ => ())
          case None =>
            // Create new
            h2Service.createBiosample(sample).map(_ => ())
        }
      }
    }

  /**
   * Sync projects from workspace to H2.
   * Creates new projects, updates existing ones.
   */
  private def syncProjects(projects: List[model.Project]): Either[String, Unit] =
    projects.foldLeft[Either[String, Unit]](Right(())) { (acc, project) =>
      acc.flatMap { _ =>
        // Check if project exists by name
        h2Service.getProjectByName(project.projectName).flatMap {
          case Some(_) =>
            // Update existing
            h2Service.updateProject(project).map(_ => ())
          case None =>
            // Create new
            h2Service.createProject(project).map(_ => ())
        }
      }
    }

  /**
   * Sync sequence runs from workspace to H2.
   * Creates new sequence runs, updates existing ones.
   */
  private def syncSequenceRuns(sequenceRuns: List[model.SequenceRun]): Either[String, Unit] =
    sequenceRuns.foldLeft[Either[String, Unit]](Right(())) { (acc, seqRun) =>
      acc.flatMap { _ =>
        // Extract biosample ID from biosampleRef
        val biosampleId = parseIdFromRef(seqRun.biosampleRef)
        biosampleId match {
          case None =>
            // Skip if we can't determine the parent biosample
            println(s"[H2WorkspaceAdapter] Skipping sequence run - invalid biosampleRef: ${seqRun.biosampleRef}")
            Right(())
          case Some(bsId) =>
            // Check if sequence run exists by its ID
            val seqRunId = seqRun.atUri.flatMap(parseIdFromRef)
            seqRunId match {
              case Some(id) =>
                h2Service.getSequenceRun(id).flatMap {
                  case Some(_) =>
                    // Update existing
                    h2Service.updateSequenceRun(seqRun).map(_ => ())
                  case None =>
                    // Create new (atUri suggests it should exist, but it doesn't - create it)
                    h2Service.createSequenceRun(seqRun, bsId).map(_ => ())
                }
              case None =>
                // No ID yet, create new
                h2Service.createSequenceRun(seqRun, bsId).map(_ => ())
            }
        }
      }
    }

  /**
   * Sync alignments from workspace to H2.
   * Creates new alignments, updates existing ones.
   */
  private def syncAlignments(alignments: List[model.Alignment]): Either[String, Unit] =
    alignments.foldLeft[Either[String, Unit]](Right(())) { (acc, alignment) =>
      acc.flatMap { _ =>
        // Extract sequence run ID from sequenceRunRef
        val seqRunId = parseIdFromRef(alignment.sequenceRunRef)
        seqRunId match {
          case None =>
            // Skip if we can't determine the parent sequence run
            println(s"[H2WorkspaceAdapter] Skipping alignment - invalid sequenceRunRef: ${alignment.sequenceRunRef}")
            Right(())
          case Some(srId) =>
            // Check if alignment exists by its ID
            val alignmentId = alignment.atUri.flatMap(parseIdFromRef)
            alignmentId match {
              case Some(id) =>
                h2Service.getAlignment(id).flatMap {
                  case Some(_) =>
                    // Update existing
                    h2Service.updateAlignment(alignment).map(_ => ())
                  case None =>
                    // Create new (atUri suggests it should exist, but it doesn't - create it)
                    h2Service.createAlignment(alignment, srId).map(_ => ())
                }
              case None =>
                // No ID yet, create new
                h2Service.createAlignment(alignment, srId).map(_ => ())
            }
        }
      }
    }

  /**
   * Parse a UUID from a reference (AT URI or local URI).
   */
  private def parseIdFromRef(ref: String): Option[java.util.UUID] =
    import com.decodingus.service.EntityConversions.parseIdFromRef as parse
    parse(ref)

object H2WorkspaceAdapter:
  /**
   * Create an adapter from a DatabaseContext.
   */
  def apply(context: DatabaseContext): H2WorkspaceAdapter =
    new H2WorkspaceAdapter(context.workspaceService)

  /**
   * Create an adapter directly from an H2WorkspaceService.
   */
  def apply(h2Service: H2WorkspaceService): H2WorkspaceAdapter =
    new H2WorkspaceAdapter(h2Service)
