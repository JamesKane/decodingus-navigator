package com.decodingus.service

import com.decodingus.workspace.model.*

import java.util.UUID

/**
 * Clean workspace service API for the H2-backed database layer.
 *
 * This replaces the JSON-based persistence with proper CRUD operations.
 * All operations return Either for explicit error handling.
 *
 * Design principles:
 * - Operations are granular (not load/save entire workspace)
 * - Returns domain models (Biosample, not BiosampleEntity)
 * - Sync status is managed internally
 * - Thread-safe for UI usage
 */
trait WorkspaceService:

  // ============================================
  // Biosample Operations
  // ============================================

  /**
   * Get all biosamples in the workspace.
   */
  def getAllBiosamples(): Either[String, List[Biosample]]

  /**
   * Find a biosample by its database ID.
   */
  def getBiosample(id: UUID): Either[String, Option[Biosample]]

  /**
   * Find a biosample by its sample accession.
   */
  def getBiosampleByAccession(accession: String): Either[String, Option[Biosample]]

  /**
   * Create a new biosample.
   * Returns the created biosample with assigned ID.
   */
  def createBiosample(biosample: Biosample): Either[String, Biosample]

  /**
   * Update an existing biosample.
   * The biosample must have a valid ID (stored in atUri or looked up by accession).
   */
  def updateBiosample(biosample: Biosample): Either[String, Biosample]

  /**
   * Delete a biosample by ID.
   * Also deletes related sequence runs and alignments (CASCADE).
   */
  def deleteBiosample(id: UUID): Either[String, Boolean]

  /**
   * Update haplogroup assignments for a biosample.
   */
  def updateBiosampleHaplogroups(id: UUID, haplogroups: HaplogroupAssignments): Either[String, Boolean]

  // ============================================
  // Project Operations
  // ============================================

  /**
   * Get all projects in the workspace.
   */
  def getAllProjects(): Either[String, List[Project]]

  /**
   * Find a project by its database ID.
   */
  def getProject(id: UUID): Either[String, Option[Project]]

  /**
   * Find a project by name.
   */
  def getProjectByName(name: String): Either[String, Option[Project]]

  /**
   * Create a new project.
   */
  def createProject(project: Project): Either[String, Project]

  /**
   * Update an existing project.
   */
  def updateProject(project: Project): Either[String, Project]

  /**
   * Delete a project by ID.
   */
  def deleteProject(id: UUID): Either[String, Boolean]

  /**
   * Get biosample members of a project.
   */
  def getProjectMembers(projectId: UUID): Either[String, List[Biosample]]

  /**
   * Add a biosample to a project.
   */
  def addProjectMember(projectId: UUID, biosampleId: UUID): Either[String, Boolean]

  /**
   * Remove a biosample from a project.
   */
  def removeProjectMember(projectId: UUID, biosampleId: UUID): Either[String, Boolean]

  // ============================================
  // SequenceRun Operations
  // ============================================

  /**
   * Get all sequence runs in the workspace.
   */
  def getAllSequenceRuns(): Either[String, List[SequenceRun]]

  /**
   * Get sequence runs for a biosample.
   */
  def getSequenceRunsForBiosample(biosampleId: UUID): Either[String, List[SequenceRun]]

  /**
   * Find a sequence run by ID.
   */
  def getSequenceRun(id: UUID): Either[String, Option[SequenceRun]]

  /**
   * Create a new sequence run.
   * The biosampleRef should contain the biosample ID or accession.
   */
  def createSequenceRun(sequenceRun: SequenceRun, biosampleId: UUID): Either[String, SequenceRun]

  /**
   * Update an existing sequence run.
   */
  def updateSequenceRun(sequenceRun: SequenceRun): Either[String, SequenceRun]

  /**
   * Delete a sequence run by ID.
   * Also deletes related alignments (CASCADE).
   */
  def deleteSequenceRun(id: UUID): Either[String, Boolean]

  // ============================================
  // Alignment Operations
  // ============================================

  /**
   * Get all alignments in the workspace.
   */
  def getAllAlignments(): Either[String, List[Alignment]]

  /**
   * Get alignments for a sequence run.
   */
  def getAlignmentsForSequenceRun(sequenceRunId: UUID): Either[String, List[Alignment]]

  /**
   * Get alignments for a biosample (via sequence runs).
   */
  def getAlignmentsForBiosample(biosampleId: UUID): Either[String, List[Alignment]]

  /**
   * Find an alignment by ID.
   */
  def getAlignment(id: UUID): Either[String, Option[Alignment]]

  /**
   * Create a new alignment.
   */
  def createAlignment(alignment: Alignment, sequenceRunId: UUID): Either[String, Alignment]

  /**
   * Update an existing alignment.
   */
  def updateAlignment(alignment: Alignment): Either[String, Alignment]

  /**
   * Update alignment metrics.
   */
  def updateAlignmentMetrics(id: UUID, metrics: AlignmentMetrics): Either[String, Boolean]

  /**
   * Delete an alignment by ID.
   */
  def deleteAlignment(id: UUID): Either[String, Boolean]

  // ============================================
  // STR Profile Operations
  // ============================================

  /**
   * Get all STR profiles in the workspace.
   */
  def getAllStrProfiles(): Either[String, List[StrProfile]]

  /**
   * Get STR profiles for a biosample.
   */
  def getStrProfilesForBiosample(biosampleId: UUID): Either[String, List[StrProfile]]

  /**
   * Find an STR profile by ID.
   */
  def getStrProfile(id: UUID): Either[String, Option[StrProfile]]

  /**
   * Create a new STR profile.
   */
  def createStrProfile(profile: StrProfile, biosampleId: UUID): Either[String, StrProfile]

  /**
   * Update an existing STR profile.
   */
  def updateStrProfile(profile: StrProfile): Either[String, StrProfile]

  /**
   * Delete an STR profile by ID.
   */
  def deleteStrProfile(id: UUID): Either[String, Boolean]

  // ============================================
  // Chip Profile Operations
  // ============================================

  /**
   * Get all chip profiles in the workspace.
   */
  def getAllChipProfiles(): Either[String, List[ChipProfile]]

  /**
   * Get chip profiles for a biosample.
   */
  def getChipProfilesForBiosample(biosampleId: UUID): Either[String, List[ChipProfile]]

  /**
   * Find a chip profile by ID.
   */
  def getChipProfile(id: UUID): Either[String, Option[ChipProfile]]

  /**
   * Create a new chip profile.
   */
  def createChipProfile(profile: ChipProfile, biosampleId: UUID): Either[String, ChipProfile]

  /**
   * Update an existing chip profile.
   */
  def updateChipProfile(profile: ChipProfile): Either[String, ChipProfile]

  /**
   * Delete a chip profile by ID.
   */
  def deleteChipProfile(id: UUID): Either[String, Boolean]

  /**
   * Find a chip profile by source file hash (for deduplication).
   */
  def getChipProfileBySourceHash(hash: String): Either[String, Option[ChipProfile]]

  // ============================================
  // Sync Status Operations
  // ============================================

  /**
   * Get counts of entities by sync status.
   */
  def getSyncStatusSummary(): Either[String, SyncStatusSummary]

  /**
   * Get all entities pending sync (Local or Modified status).
   */
  def getPendingSyncEntities(): Either[String, PendingSyncEntities]

  // ============================================
  // Bulk Operations
  // ============================================

  /**
   * Load the full workspace content (for compatibility/migration).
   * Prefer granular operations for normal use.
   */
  def loadWorkspaceContent(): Either[String, WorkspaceContent]

/**
 * Summary of sync status across all entity types.
 */
case class SyncStatusSummary(
                              localCount: Int,
                              syncedCount: Int,
                              modifiedCount: Int,
                              conflictCount: Int
                            )

/**
 * Entities pending synchronization to PDS.
 */
case class PendingSyncEntities(
                                biosamples: List[Biosample],
                                projects: List[Project],
                                sequenceRuns: List[SequenceRun],
                                alignments: List[Alignment]
                              )
