package com.decodingus.repository

import com.decodingus.repository.SqlHelpers.*
import com.decodingus.workspace.model.FileInfo
import io.circe.*
import io.circe.generic.semiauto.*
import io.circe.parser.*
import io.circe.syntax.*
import java.sql.{Connection, ResultSet, Timestamp}
import java.time.LocalDateTime
import java.util.UUID

/**
 * SequenceRun entity for database persistence.
 *
 * Represents a single sequencing session with platform info, metrics, and file references.
 */
case class SequenceRunEntity(
  id: UUID,
  biosampleId: UUID,
  platform: String,
  instrumentModel: Option[String],
  instrumentId: Option[String],
  testType: String,
  libraryId: Option[String],
  platformUnit: Option[String],
  libraryLayout: Option[String],
  sampleName: Option[String],
  sequencingFacility: Option[String],
  runFingerprint: Option[String],
  totalReads: Option[Long],
  pfReads: Option[Long],
  pfReadsAligned: Option[Long],
  readLength: Option[Int],
  meanInsertSize: Option[Double],
  medianInsertSize: Option[Double],
  stdInsertSize: Option[Double],
  flowcellId: Option[String],
  runDate: Option[LocalDateTime],
  files: List[FileInfo],
  meta: EntityMeta
) extends Entity[UUID]

object SequenceRunEntity:
  // Circe codecs for FileInfo JSON
  given Encoder[FileInfo] = deriveEncoder
  given Decoder[FileInfo] = deriveDecoder

  /**
   * Create a new SequenceRunEntity with generated ID and initial metadata.
   */
  def create(
    biosampleId: UUID,
    platform: String,
    testType: String,
    instrumentModel: Option[String] = None,
    instrumentId: Option[String] = None,
    libraryId: Option[String] = None,
    platformUnit: Option[String] = None,
    libraryLayout: Option[String] = None,
    sampleName: Option[String] = None,
    sequencingFacility: Option[String] = None,
    runFingerprint: Option[String] = None
  ): SequenceRunEntity = SequenceRunEntity(
    id = UUID.randomUUID(),
    biosampleId = biosampleId,
    platform = platform,
    instrumentModel = instrumentModel,
    instrumentId = instrumentId,
    testType = testType,
    libraryId = libraryId,
    platformUnit = platformUnit,
    libraryLayout = libraryLayout,
    sampleName = sampleName,
    sequencingFacility = sequencingFacility,
    runFingerprint = runFingerprint,
    totalReads = None,
    pfReads = None,
    pfReadsAligned = None,
    readLength = None,
    meanInsertSize = None,
    medianInsertSize = None,
    stdInsertSize = None,
    flowcellId = None,
    runDate = None,
    files = List.empty,
    meta = EntityMeta.create()
  )

/**
 * Repository for sequence run persistence operations.
 */
class SequenceRunRepository extends SyncableRepository[SequenceRunEntity, UUID]:

  import SequenceRunEntity.given

  // ============================================
  // Core Repository Operations
  // ============================================

  override def findById(id: UUID)(using conn: Connection): Option[SequenceRunEntity] =
    queryOne(
      "SELECT * FROM sequence_run WHERE id = ?",
      Seq(id)
    )(mapRow)

  override def findAll()(using conn: Connection): List[SequenceRunEntity] =
    queryList("SELECT * FROM sequence_run ORDER BY created_at DESC")(mapRow)

  override def insert(entity: SequenceRunEntity)(using conn: Connection): SequenceRunEntity =
    val filesJson = JsonValue(entity.files.asJson.noSpaces)

    executeUpdate(
      """INSERT INTO sequence_run (
        |  id, biosample_id, platform, instrument_model, instrument_id, test_type,
        |  library_id, platform_unit, library_layout, sample_name, sequencing_facility,
        |  run_fingerprint, total_reads, pf_reads, pf_reads_aligned, read_length,
        |  mean_insert_size, median_insert_size, std_insert_size, flowcell_id, run_date,
        |  files, sync_status, at_uri, at_cid, version, created_at, updated_at
        |) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
      """.stripMargin,
      Seq(
        entity.id,
        entity.biosampleId,
        entity.platform,
        entity.instrumentModel,
        entity.instrumentId,
        entity.testType,
        entity.libraryId,
        entity.platformUnit,
        entity.libraryLayout,
        entity.sampleName,
        entity.sequencingFacility,
        entity.runFingerprint,
        entity.totalReads,
        entity.pfReads,
        entity.pfReadsAligned,
        entity.readLength,
        entity.meanInsertSize,
        entity.medianInsertSize,
        entity.stdInsertSize,
        entity.flowcellId,
        entity.runDate.map(Timestamp.valueOf).orNull,
        filesJson,
        entity.meta.syncStatus,
        entity.meta.atUri,
        entity.meta.atCid,
        entity.meta.version,
        entity.meta.createdAt,
        entity.meta.updatedAt
      )
    )
    entity

  override def update(entity: SequenceRunEntity)(using conn: Connection): SequenceRunEntity =
    val updatedMeta = EntityMeta.forUpdate(entity.meta)
    val filesJson = JsonValue(entity.files.asJson.noSpaces)

    executeUpdate(
      """UPDATE sequence_run SET
        |  biosample_id = ?, platform = ?, instrument_model = ?, instrument_id = ?,
        |  test_type = ?, library_id = ?, platform_unit = ?, library_layout = ?,
        |  sample_name = ?, sequencing_facility = ?, run_fingerprint = ?,
        |  total_reads = ?, pf_reads = ?, pf_reads_aligned = ?, read_length = ?,
        |  mean_insert_size = ?, median_insert_size = ?, std_insert_size = ?,
        |  flowcell_id = ?, run_date = ?, files = ?,
        |  sync_status = ?, at_uri = ?, at_cid = ?, version = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        entity.biosampleId,
        entity.platform,
        entity.instrumentModel,
        entity.instrumentId,
        entity.testType,
        entity.libraryId,
        entity.platformUnit,
        entity.libraryLayout,
        entity.sampleName,
        entity.sequencingFacility,
        entity.runFingerprint,
        entity.totalReads,
        entity.pfReads,
        entity.pfReadsAligned,
        entity.readLength,
        entity.meanInsertSize,
        entity.medianInsertSize,
        entity.stdInsertSize,
        entity.flowcellId,
        entity.runDate.map(Timestamp.valueOf).orNull,
        filesJson,
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
    // CASCADE delete will remove alignments
    executeUpdate("DELETE FROM sequence_run WHERE id = ?", Seq(id)) > 0

  override def count()(using conn: Connection): Long =
    queryOne("SELECT COUNT(*) FROM sequence_run") { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  override def exists(id: UUID)(using conn: Connection): Boolean =
    queryOne("SELECT 1 FROM sequence_run WHERE id = ?", Seq(id)) { _ => true }.isDefined

  // ============================================
  // Syncable Repository Operations
  // ============================================

  override def findByStatus(status: SyncStatus)(using conn: Connection): List[SequenceRunEntity] =
    queryList(
      "SELECT * FROM sequence_run WHERE sync_status = ? ORDER BY updated_at DESC",
      Seq(status.toString)
    )(mapRow)

  override def updateStatus(id: UUID, status: SyncStatus)(using conn: Connection): Boolean =
    executeUpdate(
      "UPDATE sequence_run SET sync_status = ?, updated_at = ? WHERE id = ?",
      Seq(status.toString, LocalDateTime.now(), id)
    ) > 0

  override def markSynced(id: UUID, atUri: String, atCid: String)(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE sequence_run SET
        |  sync_status = ?, at_uri = ?, at_cid = ?, updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(SyncStatus.Synced.toString, atUri, atCid, LocalDateTime.now(), id)
    ) > 0

  // ============================================
  // SequenceRun-Specific Queries
  // ============================================

  /**
   * Find all sequence runs for a biosample.
   */
  def findByBiosample(biosampleId: UUID)(using conn: Connection): List[SequenceRunEntity] =
    queryList(
      "SELECT * FROM sequence_run WHERE biosample_id = ? ORDER BY created_at DESC",
      Seq(biosampleId)
    )(mapRow)

  /**
   * Find sequence run by library ID.
   */
  def findByLibraryId(libraryId: String)(using conn: Connection): Option[SequenceRunEntity] =
    queryOne(
      "SELECT * FROM sequence_run WHERE library_id = ?",
      Seq(libraryId)
    )(mapRow)

  /**
   * Find sequence run by platform unit (flowcell.lane.barcode).
   */
  def findByPlatformUnit(platformUnit: String)(using conn: Connection): Option[SequenceRunEntity] =
    queryOne(
      "SELECT * FROM sequence_run WHERE platform_unit = ?",
      Seq(platformUnit)
    )(mapRow)

  /**
   * Find sequence run by fingerprint.
   * Used for matching same run across different reference alignments.
   */
  def findByFingerprint(fingerprint: String)(using conn: Connection): List[SequenceRunEntity] =
    queryList(
      "SELECT * FROM sequence_run WHERE run_fingerprint = ?",
      Seq(fingerprint)
    )(mapRow)

  /**
   * Find sequence runs by platform.
   */
  def findByPlatform(platform: String)(using conn: Connection): List[SequenceRunEntity] =
    queryList(
      "SELECT * FROM sequence_run WHERE platform = ? ORDER BY created_at DESC",
      Seq(platform)
    )(mapRow)

  /**
   * Find sequence runs by test type.
   */
  def findByTestType(testType: String)(using conn: Connection): List[SequenceRunEntity] =
    queryList(
      "SELECT * FROM sequence_run WHERE test_type = ? ORDER BY created_at DESC",
      Seq(testType)
    )(mapRow)

  /**
   * Find sequence runs pending sync.
   */
  def findPendingSync()(using conn: Connection): List[SequenceRunEntity] =
    queryList(
      """SELECT * FROM sequence_run
        |WHERE sync_status IN ('Local', 'Modified')
        |ORDER BY updated_at ASC
      """.stripMargin
    )(mapRow)

  /**
   * Update metrics for a sequence run.
   * Commonly called after analysis completes.
   */
  def updateMetrics(
    id: UUID,
    totalReads: Option[Long] = None,
    pfReads: Option[Long] = None,
    pfReadsAligned: Option[Long] = None,
    readLength: Option[Int] = None,
    meanInsertSize: Option[Double] = None,
    medianInsertSize: Option[Double] = None,
    stdInsertSize: Option[Double] = None
  )(using conn: Connection): Boolean =
    executeUpdate(
      """UPDATE sequence_run SET
        |  total_reads = COALESCE(?, total_reads),
        |  pf_reads = COALESCE(?, pf_reads),
        |  pf_reads_aligned = COALESCE(?, pf_reads_aligned),
        |  read_length = COALESCE(?, read_length),
        |  mean_insert_size = COALESCE(?, mean_insert_size),
        |  median_insert_size = COALESCE(?, median_insert_size),
        |  std_insert_size = COALESCE(?, std_insert_size),
        |  sync_status = CASE WHEN sync_status = 'Synced' THEN 'Modified' ELSE sync_status END,
        |  updated_at = ?
        |WHERE id = ?
      """.stripMargin,
      Seq(
        totalReads, pfReads, pfReadsAligned, readLength,
        meanInsertSize, medianInsertSize, stdInsertSize,
        LocalDateTime.now(), id
      )
    ) > 0

  /**
   * Add a file to the sequence run.
   */
  def addFile(id: UUID, file: FileInfo)(using conn: Connection): Boolean =
    // Get current files
    findById(id) match
      case Some(entity) =>
        val updatedFiles = entity.files :+ file
        val filesJson = JsonValue(updatedFiles.asJson.noSpaces)
        executeUpdate(
          """UPDATE sequence_run SET
            |  files = ?,
            |  sync_status = CASE WHEN sync_status = 'Synced' THEN 'Modified' ELSE sync_status END,
            |  updated_at = ?
            |WHERE id = ?
          """.stripMargin,
          Seq(filesJson, LocalDateTime.now(), id)
        ) > 0
      case None => false

  /**
   * Count sequence runs per biosample (for summary stats).
   */
  def countByBiosample(biosampleId: UUID)(using conn: Connection): Long =
    queryOne(
      "SELECT COUNT(*) FROM sequence_run WHERE biosample_id = ?",
      Seq(biosampleId)
    ) { rs =>
      rs.getLong(1)
    }.getOrElse(0L)

  // ============================================
  // Result Set Mapping
  // ============================================

  private def mapRow(rs: ResultSet): SequenceRunEntity =
    val filesJson = getOptJsonString(rs, "files").getOrElse("[]")
    val files = parse(filesJson).flatMap(_.as[List[FileInfo]]).getOrElse(List.empty)

    SequenceRunEntity(
      id = getUUID(rs, "id"),
      biosampleId = getUUID(rs, "biosample_id"),
      platform = rs.getString("platform"),
      instrumentModel = getOptString(rs, "instrument_model"),
      instrumentId = getOptString(rs, "instrument_id"),
      testType = rs.getString("test_type"),
      libraryId = getOptString(rs, "library_id"),
      platformUnit = getOptString(rs, "platform_unit"),
      libraryLayout = getOptString(rs, "library_layout"),
      sampleName = getOptString(rs, "sample_name"),
      sequencingFacility = getOptString(rs, "sequencing_facility"),
      runFingerprint = getOptString(rs, "run_fingerprint"),
      totalReads = getOptLong(rs, "total_reads"),
      pfReads = getOptLong(rs, "pf_reads"),
      pfReadsAligned = getOptLong(rs, "pf_reads_aligned"),
      readLength = getOptInt(rs, "read_length"),
      meanInsertSize = getOptDouble(rs, "mean_insert_size"),
      medianInsertSize = getOptDouble(rs, "median_insert_size"),
      stdInsertSize = getOptDouble(rs, "std_insert_size"),
      flowcellId = getOptString(rs, "flowcell_id"),
      runDate = getOptDateTime(rs, "run_date"),
      files = files,
      meta = EntityMeta(
        syncStatus = SyncStatus.fromString(rs.getString("sync_status")),
        atUri = getOptString(rs, "at_uri"),
        atCid = getOptString(rs, "at_cid"),
        version = rs.getInt("version"),
        createdAt = getDateTime(rs, "created_at"),
        updatedAt = getDateTime(rs, "updated_at")
      )
    )
