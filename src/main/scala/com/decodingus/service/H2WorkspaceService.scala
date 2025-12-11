package com.decodingus.service

import com.decodingus.db.Transactor
import com.decodingus.repository.{
  BiosampleRepository, ProjectRepository, SequenceRunRepository, AlignmentRepository,
  SyncStatus as RepoSyncStatus
}
import com.decodingus.service.EntityConversions.*
import com.decodingus.workspace.model.{
  Biosample, Project, SequenceRun, Alignment,
  WorkspaceContent, RecordMeta, HaplogroupAssignments, AlignmentMetrics
}
import java.util.UUID

/**
 * H2 database-backed implementation of WorkspaceService.
 *
 * All operations are transactional and thread-safe.
 * Uses the repository layer for data access.
 */
class H2WorkspaceService(
  transactor: Transactor,
  biosampleRepo: BiosampleRepository,
  projectRepo: ProjectRepository,
  sequenceRunRepo: SequenceRunRepository,
  alignmentRepo: AlignmentRepository
) extends WorkspaceService:

  // ============================================
  // Biosample Operations
  // ============================================

  override def getAllBiosamples(): Either[String, List[Biosample]] =
    transactor.readOnly {
      biosampleRepo.findAll().map(fromBiosampleEntity)
    }

  override def getBiosample(id: UUID): Either[String, Option[Biosample]] =
    transactor.readOnly {
      biosampleRepo.findById(id).map(fromBiosampleEntity)
    }

  override def getBiosampleByAccession(accession: String): Either[String, Option[Biosample]] =
    transactor.readOnly {
      biosampleRepo.findByAccession(accession).map(fromBiosampleEntity)
    }

  override def createBiosample(biosample: Biosample): Either[String, Biosample] =
    transactor.readWrite {
      // Check for duplicate accession
      biosampleRepo.findByAccession(biosample.sampleAccession) match
        case Some(_) =>
          throw new IllegalArgumentException(s"Biosample with accession ${biosample.sampleAccession} already exists")
        case None =>
          val entity = toBiosampleEntity(biosample)
          val saved = biosampleRepo.insert(entity)
          fromBiosampleEntity(saved)
    }

  override def updateBiosample(biosample: Biosample): Either[String, Biosample] =
    transactor.readWrite {
      // Find existing entity by accession or atUri
      val existingId = biosample.atUri.flatMap(parseIdFromRef)
        .orElse(biosampleRepo.findByAccession(biosample.sampleAccession).map(_.id))

      existingId match
        case Some(id) =>
          val entity = toBiosampleEntity(biosample, Some(id))
          val updated = biosampleRepo.update(entity)
          fromBiosampleEntity(updated)
        case None =>
          throw new IllegalArgumentException(s"Biosample not found: ${biosample.sampleAccession}")
    }

  override def deleteBiosample(id: UUID): Either[String, Boolean] =
    transactor.readWrite {
      biosampleRepo.delete(id)
    }

  override def updateBiosampleHaplogroups(id: UUID, haplogroups: HaplogroupAssignments): Either[String, Boolean] =
    transactor.readWrite {
      biosampleRepo.updateHaplogroups(id, haplogroups)
    }

  // ============================================
  // Project Operations
  // ============================================

  override def getAllProjects(): Either[String, List[Project]] =
    transactor.readOnly {
      projectRepo.findAll().map { entity =>
        val memberRefs = projectRepo.getMemberIds(entity.id).map(id => localUri("biosample", id))
        fromProjectEntity(entity, memberRefs)
      }
    }

  override def getProject(id: UUID): Either[String, Option[Project]] =
    transactor.readOnly {
      projectRepo.findById(id).map { entity =>
        val memberRefs = projectRepo.getMemberIds(entity.id).map(id => localUri("biosample", id))
        fromProjectEntity(entity, memberRefs)
      }
    }

  override def getProjectByName(name: String): Either[String, Option[Project]] =
    transactor.readOnly {
      projectRepo.findByName(name).map { entity =>
        val memberRefs = projectRepo.getMemberIds(entity.id).map(id => localUri("biosample", id))
        fromProjectEntity(entity, memberRefs)
      }
    }

  override def createProject(project: Project): Either[String, Project] =
    transactor.readWrite {
      // Check for duplicate name
      projectRepo.findByName(project.projectName) match
        case Some(_) =>
          throw new IllegalArgumentException(s"Project with name ${project.projectName} already exists")
        case None =>
          val entity = toProjectEntity(project)
          val saved = projectRepo.insert(entity)
          fromProjectEntity(saved)
    }

  override def updateProject(project: Project): Either[String, Project] =
    transactor.readWrite {
      val existingId = project.atUri.flatMap(parseIdFromRef)
        .orElse(projectRepo.findByName(project.projectName).map(_.id))

      existingId match
        case Some(id) =>
          val entity = toProjectEntity(project, Some(id))
          val updated = projectRepo.update(entity)
          val memberRefs = projectRepo.getMemberIds(id).map(mid => localUri("biosample", mid))
          fromProjectEntity(updated, memberRefs)
        case None =>
          throw new IllegalArgumentException(s"Project not found: ${project.projectName}")
    }

  override def deleteProject(id: UUID): Either[String, Boolean] =
    transactor.readWrite {
      projectRepo.delete(id)
    }

  override def getProjectMembers(projectId: UUID): Either[String, List[Biosample]] =
    transactor.readOnly {
      val memberIds = projectRepo.getMemberIds(projectId)
      memberIds.flatMap(id => biosampleRepo.findById(id).map(fromBiosampleEntity))
    }

  override def addProjectMember(projectId: UUID, biosampleId: UUID): Either[String, Boolean] =
    transactor.readWrite {
      projectRepo.addMember(projectId, biosampleId)
    }

  override def removeProjectMember(projectId: UUID, biosampleId: UUID): Either[String, Boolean] =
    transactor.readWrite {
      projectRepo.removeMember(projectId, biosampleId)
    }

  // ============================================
  // SequenceRun Operations
  // ============================================

  override def getAllSequenceRuns(): Either[String, List[SequenceRun]] =
    transactor.readOnly {
      sequenceRunRepo.findAll().map { entity =>
        val biosampleRef = biosampleRepo.findById(entity.biosampleId)
          .map(b => localUri("biosample", b.id))
          .getOrElse(s"unknown:biosample:${entity.biosampleId}")
        fromSequenceRunEntity(entity, biosampleRef)
      }
    }

  override def getSequenceRunsForBiosample(biosampleId: UUID): Either[String, List[SequenceRun]] =
    transactor.readOnly {
      val biosampleRef = localUri("biosample", biosampleId)
      sequenceRunRepo.findByBiosample(biosampleId).map { entity =>
        fromSequenceRunEntity(entity, biosampleRef)
      }
    }

  override def getSequenceRun(id: UUID): Either[String, Option[SequenceRun]] =
    transactor.readOnly {
      sequenceRunRepo.findById(id).map { entity =>
        val biosampleRef = biosampleRepo.findById(entity.biosampleId)
          .map(b => localUri("biosample", b.id))
          .getOrElse(s"unknown:biosample:${entity.biosampleId}")
        fromSequenceRunEntity(entity, biosampleRef)
      }
    }

  override def createSequenceRun(sequenceRun: SequenceRun, biosampleId: UUID): Either[String, SequenceRun] =
    transactor.readWrite {
      // Verify biosample exists
      biosampleRepo.findById(biosampleId) match
        case None =>
          throw new IllegalArgumentException(s"Biosample not found: $biosampleId")
        case Some(biosample) =>
          val entity = toSequenceRunEntity(sequenceRun, biosampleId)
          val saved = sequenceRunRepo.insert(entity)
          fromSequenceRunEntity(saved, localUri("biosample", biosampleId))
    }

  override def updateSequenceRun(sequenceRun: SequenceRun): Either[String, SequenceRun] =
    transactor.readWrite {
      val existingId = sequenceRun.atUri.flatMap(parseIdFromRef)

      existingId match
        case Some(id) =>
          sequenceRunRepo.findById(id) match
            case Some(existing) =>
              val entity = toSequenceRunEntity(sequenceRun, existing.biosampleId, Some(id))
              val updated = sequenceRunRepo.update(entity)
              val biosampleRef = localUri("biosample", updated.biosampleId)
              fromSequenceRunEntity(updated, biosampleRef)
            case None =>
              throw new IllegalArgumentException(s"SequenceRun not found: $id")
        case None =>
          throw new IllegalArgumentException("SequenceRun has no valid ID")
    }

  override def deleteSequenceRun(id: UUID): Either[String, Boolean] =
    transactor.readWrite {
      sequenceRunRepo.delete(id)
    }

  // ============================================
  // Alignment Operations
  // ============================================

  override def getAllAlignments(): Either[String, List[Alignment]] =
    transactor.readOnly {
      alignmentRepo.findAll().map { entity =>
        val sequenceRunRef = localUri("sequencerun", entity.sequenceRunId)
        fromAlignmentEntity(entity, sequenceRunRef)
      }
    }

  override def getAlignmentsForSequenceRun(sequenceRunId: UUID): Either[String, List[Alignment]] =
    transactor.readOnly {
      val sequenceRunRef = localUri("sequencerun", sequenceRunId)
      alignmentRepo.findBySequenceRun(sequenceRunId).map { entity =>
        fromAlignmentEntity(entity, sequenceRunRef)
      }
    }

  override def getAlignmentsForBiosample(biosampleId: UUID): Either[String, List[Alignment]] =
    transactor.readOnly {
      alignmentRepo.findByBiosample(biosampleId).map { entity =>
        val sequenceRunRef = localUri("sequencerun", entity.sequenceRunId)
        fromAlignmentEntity(entity, sequenceRunRef)
      }
    }

  override def getAlignment(id: UUID): Either[String, Option[Alignment]] =
    transactor.readOnly {
      alignmentRepo.findById(id).map { entity =>
        val sequenceRunRef = localUri("sequencerun", entity.sequenceRunId)
        fromAlignmentEntity(entity, sequenceRunRef)
      }
    }

  override def createAlignment(alignment: Alignment, sequenceRunId: UUID): Either[String, Alignment] =
    transactor.readWrite {
      // Verify sequence run exists
      sequenceRunRepo.findById(sequenceRunId) match
        case None =>
          throw new IllegalArgumentException(s"SequenceRun not found: $sequenceRunId")
        case Some(_) =>
          val entity = toAlignmentEntity(alignment, sequenceRunId)
          val saved = alignmentRepo.insert(entity)
          fromAlignmentEntity(saved, localUri("sequencerun", sequenceRunId))
    }

  override def updateAlignment(alignment: Alignment): Either[String, Alignment] =
    transactor.readWrite {
      val existingId = alignment.atUri.flatMap(parseIdFromRef)

      existingId match
        case Some(id) =>
          alignmentRepo.findById(id) match
            case Some(existing) =>
              val entity = toAlignmentEntity(alignment, existing.sequenceRunId, Some(id))
              val updated = alignmentRepo.update(entity)
              val sequenceRunRef = localUri("sequencerun", updated.sequenceRunId)
              fromAlignmentEntity(updated, sequenceRunRef)
            case None =>
              throw new IllegalArgumentException(s"Alignment not found: $id")
        case None =>
          throw new IllegalArgumentException("Alignment has no valid ID")
    }

  override def updateAlignmentMetrics(id: UUID, metrics: AlignmentMetrics): Either[String, Boolean] =
    transactor.readWrite {
      alignmentRepo.updateMetrics(id, metrics)
    }

  override def deleteAlignment(id: UUID): Either[String, Boolean] =
    transactor.readWrite {
      alignmentRepo.delete(id)
    }

  // ============================================
  // Sync Status Operations
  // ============================================

  override def getSyncStatusSummary(): Either[String, SyncStatusSummary] =
    transactor.readOnly {
      val local = biosampleRepo.findByStatus(RepoSyncStatus.Local).size +
        projectRepo.findByStatus(RepoSyncStatus.Local).size +
        sequenceRunRepo.findByStatus(RepoSyncStatus.Local).size +
        alignmentRepo.findByStatus(RepoSyncStatus.Local).size

      val synced = biosampleRepo.findByStatus(RepoSyncStatus.Synced).size +
        projectRepo.findByStatus(RepoSyncStatus.Synced).size +
        sequenceRunRepo.findByStatus(RepoSyncStatus.Synced).size +
        alignmentRepo.findByStatus(RepoSyncStatus.Synced).size

      val modified = biosampleRepo.findByStatus(RepoSyncStatus.Modified).size +
        projectRepo.findByStatus(RepoSyncStatus.Modified).size +
        sequenceRunRepo.findByStatus(RepoSyncStatus.Modified).size +
        alignmentRepo.findByStatus(RepoSyncStatus.Modified).size

      val conflict = biosampleRepo.findByStatus(RepoSyncStatus.Conflict).size +
        projectRepo.findByStatus(RepoSyncStatus.Conflict).size +
        sequenceRunRepo.findByStatus(RepoSyncStatus.Conflict).size +
        alignmentRepo.findByStatus(RepoSyncStatus.Conflict).size

      SyncStatusSummary(local, synced, modified, conflict)
    }

  override def getPendingSyncEntities(): Either[String, PendingSyncEntities] =
    transactor.readOnly {
      val biosamples = biosampleRepo.findPendingSync().map(fromBiosampleEntity)

      val projects = projectRepo.findPendingSync().map { entity =>
        val memberRefs = projectRepo.getMemberIds(entity.id).map(id => localUri("biosample", id))
        fromProjectEntity(entity, memberRefs)
      }

      val sequenceRuns = sequenceRunRepo.findPendingSync().map { entity =>
        val biosampleRef = localUri("biosample", entity.biosampleId)
        fromSequenceRunEntity(entity, biosampleRef)
      }

      val alignments = alignmentRepo.findPendingSync().map { entity =>
        val sequenceRunRef = localUri("sequencerun", entity.sequenceRunId)
        fromAlignmentEntity(entity, sequenceRunRef)
      }

      PendingSyncEntities(biosamples, projects, sequenceRuns, alignments)
    }

  // ============================================
  // Bulk Operations
  // ============================================

  override def loadWorkspaceContent(): Either[String, WorkspaceContent] =
    transactor.readOnly {
      val samples = biosampleRepo.findAll().map(fromBiosampleEntity)

      val projects = projectRepo.findAll().map { entity =>
        val memberRefs = projectRepo.getMemberIds(entity.id).map(id => localUri("biosample", id))
        fromProjectEntity(entity, memberRefs)
      }

      val sequenceRuns = sequenceRunRepo.findAll().map { entity =>
        val biosampleRef = localUri("biosample", entity.biosampleId)
        fromSequenceRunEntity(entity, biosampleRef)
      }

      val alignments = alignmentRepo.findAll().map { entity =>
        val sequenceRunRef = localUri("sequencerun", entity.sequenceRunId)
        fromAlignmentEntity(entity, sequenceRunRef)
      }

      WorkspaceContent(
        meta = Some(RecordMeta.initial),
        samples = samples,
        projects = projects,
        sequenceRuns = sequenceRuns,
        alignments = alignments
      )
    }
