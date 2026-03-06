package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*

import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Project entity for database persistence.
 *
 * Projects group biosamples for research purposes.
 * Membership is managed via a junction table, not stored directly on the entity.
 */
case class ProjectEntity(
                          id: UUID,
                          projectName: String,
                          description: Option[String],
                          administratorDid: String,
                          meta: EntityMeta
                        ) extends Entity[UUID]

object ProjectEntity:
  /**
   * Create a new ProjectEntity with generated ID and initial metadata.
   */
  def create(
              projectName: String,
              administratorDid: String,
              description: Option[String] = None
            ): ProjectEntity = ProjectEntity(
    id = UUID.randomUUID(),
    projectName = projectName,
    description = description,
    administratorDid = administratorDid,
    meta = EntityMeta.create()
  )

/**
 * A project membership record from the junction table.
 */
case class ProjectMembership(
                              projectId: UUID,
                              biosampleId: UUID,
                              addedAt: LocalDateTime
                            )

/**
 * Repository for project persistence operations.
 *
 * Includes methods for managing project membership via the junction table.
 */
class ProjectRepository extends SyncableRepositoryBase[ProjectEntity]:

  override protected def tableName: String = "project"

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[ProjectEntity] =
    queryOne(
      "SELECT * FROM project WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[ProjectEntity] =
    queryList("SELECT * FROM project ORDER BY project_name")(mapRow)

  override def insert(entity: ProjectEntity)(using conn: Connection): ProjectEntity =
    executeUpdate(
      """INSERT INTO project (
        |  id, project_name, description, administrator_did, sync_status,
        |  at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.projectName,
        entity.description,
        entity.administratorDid,
        entity.meta.syncStatus,
        entity.meta.atUri,
        entity.meta.atCid,
        entity.meta.version,
        entity.meta.createdAt,
        entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: ProjectEntity)(using conn: Connection): ProjectEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)

    executeUpdate(
      """UPDATE project SET
        |  project_name = ?, description = ?, administrator_did = ?,
        |  sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.projectName,
        entity.description,
        entity.administratorDid,
        updatedMeta.syncStatus,
        updatedMeta.atUri,
        updatedMeta.atCid,
        updatedMeta.version,
        updatedMeta.updatedAt,
        entity.id
      )
    )
    entity.copy(meta = updatedMeta)

  override def delete(id: UUID)(using conn: Connection): Boolean =
    // CASCADE delete will remove project_member entries
    executeUpdate("DELETE FROM project WHERE id = ?", Seq(id)) > 0

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM project WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Project-Specific Queries
  // ============================================

  /**
   * Find a project by name.
   */
  def findByName(name: String)(using conn: Connection): Option[ProjectEntity] =
    queryOne(
      "SELECT * FROM project WHERE project_name = ?",
      Seq(name)
    )(mapRow)

  /**
   * Find all projects administered by a user.
   */
  def findByAdministrator(did: String)(using conn: Connection): List[ProjectEntity] =
    queryList(
      "SELECT * FROM project WHERE administrator_did = ? ORDER BY project_name",
      Seq(did)
    )(mapRow)

  /**
   * Search projects by name prefix.
   */
  def searchByName(prefix: String)(using conn: Connection): List[ProjectEntity] =
    queryList(
      "SELECT * FROM project WHERE project_name LIKE ? ORDER BY project_name",
      Seq(s"$prefix%")
    )(mapRow)

  // ============================================
  // Membership Operations (Junction Table)
  // ============================================

  /**
   * Add a biosample to a project.
   * Idempotent - silently succeeds if already a member.
   */
  def addMember(projectId: UUID, biosampleId: UUID)(using conn: Connection): Boolean =
    try
      executeUpdate(
        """INSERT INTO project_member (project_id, biosample_id, added_at)
          |VALUES (?, ?, ?)
        """.stripMargin,
        Seq(projectId, biosampleId, LocalDateTime.now())
      )
      // Mark project as modified if it was synced
      executeUpdate(
        """UPDATE project SET
          |  sync_status = CASE WHEN sync_status = 'Synced' THEN 'Modified' ELSE sync_status END,
          |  updated_at = ?
          |WHERE id = ?
        """.stripMargin,
        Seq(LocalDateTime.now(), projectId)
      )
      true
    catch
      case _: java.sql.SQLIntegrityConstraintViolationException => true // Already a member

  /**
   * Remove a biosample from a project.
   */
  def removeMember(projectId: UUID, biosampleId: UUID)(using conn: Connection): Boolean =
    val removed = executeUpdate(
      "DELETE FROM project_member WHERE project_id = ? AND biosample_id = ?",
      Seq(projectId, biosampleId)
    ) > 0

    if removed then
      // Mark project as modified if it was synced
      executeUpdate(
        """UPDATE project SET
          |  sync_status = CASE WHEN sync_status = 'Synced' THEN 'Modified' ELSE sync_status END,
          |  updated_at = ?
          |WHERE id = ?
        """.stripMargin,
        Seq(LocalDateTime.now(), projectId)
      )

    removed

  /**
   * Check if a biosample is a member of a project.
   */
  def isMember(projectId: UUID, biosampleId: UUID)(using conn: Connection): Boolean =
    queryOne(
      "SELECT 1 FROM project_member WHERE project_id = ? AND biosample_id = ?",
      Seq(projectId, biosampleId)
    ) { _ => true }.isDefined

  /**
   * Get all member biosample IDs for a project.
   */
  def getMemberIds(projectId: UUID)(using conn: Connection): List[UUID] =
    queryList(
      "SELECT biosample_id FROM project_member WHERE project_id = ? ORDER BY added_at",
      Seq(projectId)
    ) { rs =>
      getUUID(rs, "biosample_id")
    }

  /**
   * Get all memberships for a project with timestamps.
   */
  def getMemberships(projectId: UUID)(using conn: Connection): List[ProjectMembership] =
    queryList(
      "SELECT * FROM project_member WHERE project_id = ? ORDER BY added_at",
      Seq(projectId)
    ) { rs =>
      ProjectMembership(
        projectId = getUUID(rs, "project_id"),
        biosampleId = getUUID(rs, "biosample_id"),
        addedAt = getDateTime(rs, "added_at")
      )
    }

  /**
   * Get all project IDs that a biosample belongs to.
   */
  def getProjectsForBiosample(biosampleId: UUID)(using conn: Connection): List[UUID] =
    queryList(
      "SELECT project_id FROM project_member WHERE biosample_id = ?",
      Seq(biosampleId)
    ) { rs =>
      getUUID(rs, "project_id")
    }

  /**
   * Count members in a project.
   */
  def countMembers(projectId: UUID)(using conn: Connection): Long =
    queryOne(
      "SELECT COUNT(*) FROM project_member WHERE project_id = ?",
      Seq(projectId)
    ) { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  /**
   * Replace all members of a project.
   * Useful for bulk updates.
   */
  def setMembers(projectId: UUID, biosampleIds: List[UUID])(using conn: Connection): Unit =
    // Remove existing members
    executeUpdate("DELETE FROM project_member WHERE project_id = ?", Seq(projectId))

    // Add new members
    val now = LocalDateTime.now()
    for biosampleId <- biosampleIds do
      executeUpdate(
        "INSERT INTO project_member (project_id, biosample_id, added_at) VALUES (?, ?, ?)",
        Seq(projectId, biosampleId, now)
      )

    // Mark project as modified if it was synced
    executeUpdate(
      """UPDATE project SET
        |  sync_status = CASE WHEN sync_status = 'Synced' THEN 'Modified' ELSE sync_status END,
        |  updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(now, projectId)
    )

  // ============================================
  // Result Set Mapping
  // ============================================

  override protected def mapRow(rs: ResultSet): ProjectEntity =
    ProjectEntity(
      id = getUUID(rs, "id"),
      projectName = rs.getString("project_name"),
      description = getOptString(rs, "description"),
      administratorDid = rs.getString("administrator_did"),
      meta = readEntityMeta(rs)
    )
