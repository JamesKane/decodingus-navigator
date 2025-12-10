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
 * @param strProfiles              Denormalized STR profile records for local access
 * @param chipProfiles             Denormalized chip/array genotype profiles for local access
 * @param ySnpPanels               Denormalized Y-DNA SNP panel results for local access
 * @param haplogroupReconciliations Haplogroup reconciliation records (multi-run tracking)
 */
case class WorkspaceContent(
  meta: Option[RecordMeta] = None,
  sampleRefs: List[String] = List.empty,
  projectRefs: List[String] = List.empty,
  samples: List[Biosample] = List.empty,
  projects: List[Project] = List.empty,
  sequenceRuns: List[SequenceRun] = List.empty,
  alignments: List[Alignment] = List.empty,
  strProfiles: List[StrProfile] = List.empty,
  chipProfiles: List[ChipProfile] = List.empty,
  ySnpPanels: List[com.decodingus.genotype.model.YDnaSnpPanelResult] = List.empty,
  haplogroupReconciliations: List[HaplogroupReconciliation] = List.empty
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

  /**
   * Returns chip profiles for a given biosample.
   * Resolves genotypeRefs to actual ChipProfile records.
   */
  def getChipProfilesForBiosample(biosample: Biosample): List[ChipProfile] = {
    biosample.genotypeRefs.flatMap { ref =>
      chipProfiles.find(_.atUri.contains(ref))
    }
  }

  /**
   * Returns Y-DNA SNP panel results for a given biosample.
   * Resolves ySnpPanelRefs to actual YDnaSnpPanelResult records.
   */
  def getYSnpPanelsForBiosample(biosample: Biosample): List[com.decodingus.genotype.model.YDnaSnpPanelResult] = {
    biosample.ySnpPanelRefs.flatMap { ref =>
      ySnpPanels.find(_.atUri.contains(ref))
    }
  }

  /**
   * Returns haplogroup reconciliation records for a biosample.
   */
  def getReconciliationsForBiosample(biosample: Biosample): List[HaplogroupReconciliation] = {
    val biosampleRef = biosample.atUri.getOrElse(s"local:biosample:${biosample.sampleAccession}")
    haplogroupReconciliations.filter(_.biosampleRef == biosampleRef)
  }

  /**
   * Returns Y-DNA haplogroup reconciliation for a biosample, if any.
   */
  def getYDnaReconciliation(biosample: Biosample): Option[HaplogroupReconciliation] = {
    getReconciliationsForBiosample(biosample).find(_.dnaType == DnaType.Y_DNA)
  }

  /**
   * Returns mtDNA haplogroup reconciliation for a biosample, if any.
   */
  def getMtDnaReconciliation(biosample: Biosample): Option[HaplogroupReconciliation] = {
    getReconciliationsForBiosample(biosample).find(_.dnaType == DnaType.MT_DNA)
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
