package com.decodingus.workspace.model

import java.time.LocalDateTime

case class FileInfo(
  fileName: String,
  fileSizeBytes: Option[Long],
  fileFormat: String,
  checksum: Option[String],
  location: String
)
