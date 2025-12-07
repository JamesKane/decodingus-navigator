package com.decodingus.workspace.model

/**
 * The root container for a Researcher's workspace, holding all first-class records.
 * This is the local cache representation of the Atmosphere Lexicon workspace.
 *
 * In the AT Protocol model, the Workspace contains reference arrays (sampleRefs, projectRefs)
 * pointing to individual records. For local storage, we denormalize and store all records
 * together for efficient access.
 *
 * @param meta          Record metadata for tracking changes and sync
 * @param sampleRefs    AT URIs of biosample records in this workspace (for PDS sync)
 * @param projectRefs   AT URIs of project records in this workspace (for PDS sync)
 * @param samples       Denormalized biosample records for local access
 * @param projects      Denormalized project records for local access
 * @param sequenceRuns  Denormalized sequence run records for local access
 * @param alignments    Denormalized alignment records for local access
 * @param strProfiles   Denormalized STR profile records for local access
 */
case class WorkspaceContent(
  meta: Option[RecordMeta] = None,
  sampleRefs: List[String] = List.empty,
  projectRefs: List[String] = List.empty,
  samples: List[Biosample] = List.empty,
  projects: List[Project] = List.empty,
  sequenceRuns: List[SequenceRun] = List.empty,
  alignments: List[Alignment] = List.empty,
  strProfiles: List[StrProfile] = List.empty
) {
  /**
   * Returns sequence runs for a given biosample.
   * Resolves sequenceRunRefs to actual SequenceRun records.
   */
  def getSequenceRunsForBiosample(biosample: Biosample): List[SequenceRun] = {
    biosample.sequenceRunRefs.flatMap { ref =>
      sequenceRuns.find(_.atUri.contains(ref))
    }
  }

  /**
   * Returns alignments for a given sequence run.
   * Resolves alignmentRefs to actual Alignment records.
   */
  def getAlignmentsForSequenceRun(sequenceRun: SequenceRun): List[Alignment] = {
    sequenceRun.alignmentRefs.flatMap { ref =>
      alignments.find(_.atUri.contains(ref))
    }
  }

  /**
   * Returns legacy SequenceData view for a biosample.
   * This provides backward compatibility for UI components that expect embedded data.
   *
   * @deprecated Use getSequenceRunsForBiosample with the new model instead
   */
  @deprecated("Use getSequenceRunsForBiosample instead", "2.0")
  def getLegacySequenceData(biosample: Biosample): List[SequenceData] = {
    getSequenceRunsForBiosample(biosample).map { run =>
      val runAlignments = getAlignmentsForSequenceRun(run)
      SequenceData.fromSequenceRun(run, runAlignments)
    }
  }

  /**
   * Finds a biosample by sample accession.
   */
  def findBiosample(sampleAccession: String): Option[Biosample] = {
    samples.find(_.sampleAccession == sampleAccession)
  }

  /**
   * Finds a project by name.
   */
  def findProject(projectName: String): Option[Project] = {
    projects.find(_.projectName == projectName)
  }
}

/**
 * Root workspace container following the Atmosphere Lexicon structure.
 *
 * @param lexicon Schema version number
 * @param id      Namespace identifier (com.decodingus.atmosphere.workspace)
 * @param main    The workspace content
 */
case class Workspace(
  lexicon: Int,
  id: String,
  main: WorkspaceContent
)

object Workspace {
  /** Current lexicon version - increment when schema changes */
  val CurrentLexiconVersion: Int = 2

  /** Namespace ID */
  val NamespaceId: String = "com.decodingus.atmosphere.workspace"

  /** Create an empty workspace with current schema version */
  def empty: Workspace = Workspace(
    lexicon = CurrentLexiconVersion,
    id = NamespaceId,
    main = WorkspaceContent()
  )
}
