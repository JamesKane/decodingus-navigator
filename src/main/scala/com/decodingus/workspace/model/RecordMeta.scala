package com.decodingus.workspace.model

import java.time.LocalDateTime

/**
 * Metadata for tracking changes and enabling efficient sync.
 * Part of the Atmosphere Lexicon (com.decodingus.atmosphere.defs#recordMeta).
 *
 * @param version           Monotonically increasing version number for this record. Incremented on each update.
 * @param createdAt         Timestamp when this record was first created.
 * @param updatedAt         Timestamp of the most recent update.
 * @param lastModifiedField Hint about what field changed in the last update (e.g., 'haplogroups.yDna', 'description').
 */
case class RecordMeta(
                       version: Int,
                       createdAt: LocalDateTime,
                       updatedAt: Option[LocalDateTime] = None,
                       lastModifiedField: Option[String] = None
                     ) {
  /** Creates a new meta with incremented version and updated timestamp */
  def updated(field: String): RecordMeta = copy(
    version = version + 1,
    updatedAt = Some(LocalDateTime.now()),
    lastModifiedField = Some(field)
  )
}

object RecordMeta {
  /** Creates initial metadata for a new record */
  def initial: RecordMeta = RecordMeta(
    version = 1,
    createdAt = LocalDateTime.now()
  )
}
