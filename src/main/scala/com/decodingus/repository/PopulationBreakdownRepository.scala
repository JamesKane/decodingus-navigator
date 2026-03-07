package com.decodingus.repository

import com.decodingus.ancestry.model.{ConfidenceInterval, PopulationComponent, SuperPopulationSummary}
import com.decodingus.repository.SqlHelpers.*
import io.circe.*
import io.circe.generic.semiauto.*
import io.circe.parser.*
import io.circe.syntax.*

import java.sql.{Connection, ResultSet}
import java.time.LocalDateTime
import java.util.UUID

/**
 * Population Breakdown entity for database persistence.
 * Maps to the population_breakdown table (V010 migration).
 */
case class PopulationBreakdownEntity(
                                      id: UUID,
                                      biosampleId: UUID,
                                      analysisMethod: String,
                                      panelType: String,
                                      referencePopulations: String,
                                      snpsAnalyzed: Int,
                                      snpsWithGenotype: Int,
                                      snpsMissing: Int,
                                      confidenceLevel: Double,
                                      components: List[PopulationComponent],
                                      superPopulationSummary: List[SuperPopulationSummary],
                                      pcaCoordinates: Option[List[Double]],
                                      analysisDate: Option[LocalDateTime],
                                      pipelineVersion: Option[String],
                                      referenceVersion: Option[String],
                                      meta: EntityMeta
                                    ) extends Entity[UUID]

object PopulationBreakdownEntity:

  def create(
              biosampleId: UUID,
              analysisMethod: String,
              panelType: String,
              referencePopulations: String,
              snpsAnalyzed: Int,
              snpsWithGenotype: Int,
              snpsMissing: Int,
              confidenceLevel: Double,
              components: List[PopulationComponent],
              superPopulationSummary: List[SuperPopulationSummary],
              pcaCoordinates: Option[List[Double]] = None,
              analysisDate: Option[LocalDateTime] = None,
              pipelineVersion: Option[String] = None,
              referenceVersion: Option[String] = None
            ): PopulationBreakdownEntity = PopulationBreakdownEntity(
    id = UUID.randomUUID(),
    biosampleId = biosampleId,
    analysisMethod = analysisMethod,
    panelType = panelType,
    referencePopulations = referencePopulations,
    snpsAnalyzed = snpsAnalyzed,
    snpsWithGenotype = snpsWithGenotype,
    snpsMissing = snpsMissing,
    confidenceLevel = confidenceLevel,
    components = components,
    superPopulationSummary = superPopulationSummary,
    pcaCoordinates = pcaCoordinates,
    analysisDate = analysisDate,
    pipelineVersion = pipelineVersion,
    referenceVersion = referenceVersion,
    meta = EntityMeta.create()
  )

/**
 * Circe codecs for PopulationBreakdown JSON fields.
 */
object PopulationBreakdownCodecs:
  given Encoder[ConfidenceInterval] = deriveEncoder
  given Decoder[ConfidenceInterval] = deriveDecoder
  given Encoder[PopulationComponent] = deriveEncoder
  given Decoder[PopulationComponent] = deriveDecoder
  given Encoder[SuperPopulationSummary] = deriveEncoder
  given Decoder[SuperPopulationSummary] = deriveDecoder

/**
 * Repository for Population Breakdown persistence operations.
 */
class PopulationBreakdownRepository extends SyncableRepository[PopulationBreakdownEntity, UUID]:

  import PopulationBreakdownCodecs.given

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[PopulationBreakdownEntity] =
    queryOne("SELECT * FROM population_breakdown WHERE id = ?", Seq(id))(mapRow)

  override def findAll()(using conn: Connection): List[PopulationBreakdownEntity] =
    queryList("SELECT * FROM population_breakdown ORDER BY created_at DESC")(mapRow)

  override def insert(entity: PopulationBreakdownEntity)(using conn: Connection): PopulationBreakdownEntity =
    val componentsJson = JsonValue(entity.components.asJson.noSpaces)
    val summaryJson = JsonValue(entity.superPopulationSummary.asJson.noSpaces)
    val pcaJson = entity.pcaCoordinates.map(c => JsonValue(c.asJson.noSpaces))

    executeUpdate(
      """INSERT INTO population_breakdown (
        |  id, biosample_id, analysis_method, panel_type, reference_populations,
        |  snps_analyzed, snps_with_genotype, snps_missing, confidence_level,
        |  components, super_population_summary, pca_coordinates,
        |  analysis_date, pipeline_version, reference_version,
        |  sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.biosampleId,
        entity.analysisMethod,
        entity.panelType,
        entity.referencePopulations,
        entity.snpsAnalyzed,
        entity.snpsWithGenotype,
        entity.snpsMissing,
        entity.confidenceLevel,
        componentsJson,
        summaryJson,
        pcaJson,
        entity.analysisDate,
        entity.pipelineVersion,
        entity.referenceVersion,
        entity.meta.syncStatus,
        entity.meta.atUri,
        entity.meta.atCid,
        entity.meta.version,
        entity.meta.createdAt,
        entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: PopulationBreakdownEntity)(using conn: Connection): PopulationBreakdownEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    val componentsJson = JsonValue(entity.components.asJson.noSpaces)
    val summaryJson = JsonValue(entity.superPopulationSummary.asJson.noSpaces)
    val pcaJson = entity.pcaCoordinates.map(c => JsonValue(c.asJson.noSpaces))

    executeUpdate(
      """UPDATE population_breakdown SET
        |  biosample_id = ?, analysis_method = ?, panel_type = ?, reference_populations = ?,
        |  snps_analyzed = ?, snps_with_genotype = ?, snps_missing = ?, confidence_level = ?,
        |  components = ?, super_population_summary = ?, pca_coordinates = ?,
        |  analysis_date = ?, pipeline_version = ?, reference_version = ?,
        |  sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.biosampleId,
        entity.analysisMethod,
        entity.panelType,
        entity.referencePopulations,
        entity.snpsAnalyzed,
        entity.snpsWithGenotype,
        entity.snpsMissing,
        entity.confidenceLevel,
        componentsJson,
        summaryJson,
        pcaJson,
        entity.analysisDate,
        entity.pipelineVersion,
        entity.referenceVersion,
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
    executeUpdate("DELETE FROM population_breakdown WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM population_breakdown") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM population_breakdown WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Syncable Repository Operations
  // ============================================

  override def findByStatus(status: SyncStatus)(using conn: Connection): List[PopulationBreakdownEntity] =
    queryList(
      "SELECT * FROM population_breakdown WHERE sync_status = ? ORDER BY updated_at DESC",
      Seq(status.toString)
    )(mapRow)

  override def updateStatus(id: UUID, status: SyncStatus)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE population_breakdown SET sync_status = ?, updated_at = ? WHERE id = ?",
      Seq(status.toString, LocalDateTime.now(), id)
    ) > 0

  override def markSynced(id: UUID, atUri: String, atCid: String)(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE population_breakdown SET
        |  sync_status = ?, at_uri = ?, at_cid = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(SyncStatus.Synced.toString, atUri, atCid, LocalDateTime.now(), id)
    ) > 0

  // ============================================
  // PopulationBreakdown-Specific Queries
  // ============================================

  /**
   * Find breakdown for a biosample and panel type.
   */
  def findByBiosampleAndPanel(biosampleId: UUID, panelType: String)(using conn: Connection): Option[PopulationBreakdownEntity] =
    queryOne(
      "SELECT * FROM population_breakdown WHERE biosample_id = ? AND panel_type = ?",
      Seq(biosampleId, panelType)
    )(mapRow)

  /**
   * Find all breakdowns for a biosample.
   */
  def findByBiosample(biosampleId: UUID)(using conn: Connection): List[PopulationBreakdownEntity] =
    queryList(
      "SELECT * FROM population_breakdown WHERE biosample_id = ? ORDER BY panel_type",
      Seq(biosampleId)
    )(mapRow)

  /**
   * Upsert breakdown (insert if not exists, update if exists).
   * Uses the unique constraint on (biosample_id, panel_type).
   */
  def upsert(entity: PopulationBreakdownEntity)(using conn: Connection): PopulationBreakdownEntity =
    findByBiosampleAndPanel(entity.biosampleId, entity.panelType) match
      case Some(existing) =>
        update(entity.copy(id = existing.id, meta = existing.meta))
      case None =>
        insert(entity)

  /**
   * Find breakdowns pending sync.
   */
  def findPendingSync()(using conn: Connection): List[PopulationBreakdownEntity] =
    queryList(
      """SELECT * FROM population_breakdown
        |WHERE sync_status IN ('Local', 'Modified')
        |ORDER BY updated_at ASC
      """.stripMargin
    )(mapRow)

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): PopulationBreakdownEntity =
    val componentsJson = getOptJsonString(rs, "components")
    val components = componentsJson.flatMap(json =>
      parse(json).flatMap(_.as[List[PopulationComponent]]).toOption
    ).getOrElse(List.empty)

    val summaryJson = getOptJsonString(rs, "super_population_summary")
    val superPopSummary = summaryJson.flatMap(json =>
      parse(json).flatMap(_.as[List[SuperPopulationSummary]]).toOption
    ).getOrElse(List.empty)

    val pcaJson = getOptJsonString(rs, "pca_coordinates")
    val pcaCoords = pcaJson.flatMap(json =>
      parse(json).flatMap(_.as[List[Double]]).toOption
    )

    PopulationBreakdownEntity(
      id = getUUID(rs, "id"),
      biosampleId = getUUID(rs, "biosample_id"),
      analysisMethod = rs.getString("analysis_method"),
      panelType = rs.getString("panel_type"),
      referencePopulations = rs.getString("reference_populations"),
      snpsAnalyzed = rs.getInt("snps_analyzed"),
      snpsWithGenotype = rs.getInt("snps_with_genotype"),
      snpsMissing = rs.getInt("snps_missing"),
      confidenceLevel = rs.getDouble("confidence_level"),
      components = components,
      superPopulationSummary = superPopSummary,
      pcaCoordinates = pcaCoords,
      analysisDate = getOptDateTime(rs, "analysis_date"),
      pipelineVersion = getOptString(rs, "pipeline_version"),
      referenceVersion = getOptString(rs, "reference_version"),
      meta = EntityMeta(
        syncStatus = SyncStatus.fromString(rs.getString("sync_status")),
        atUri = getOptString(rs, "at_uri"),
        atCid = getOptString(rs, "at_cid"),
        version = rs.getInt("version"),
        createdAt = getDateTime(rs, "created_at"),
        updatedAt = getDateTime(rs, "updated_at")
      )
    )
