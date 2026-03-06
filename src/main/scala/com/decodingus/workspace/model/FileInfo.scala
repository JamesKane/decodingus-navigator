package com.decodingus.workspace.model

/**
 * Metadata about a data file for provenance tracking.
 * Part of the Atmosphere Lexicon (com.decodingus.atmosphere.defs#fileInfo).
 *
 * NOTE: This is metadata ONLY - DecodingUs never accesses the actual file content.
 * Files remain on the user's local system or personal storage.
 *
 * @param fileName          Name of the file
 * @param fileSizeBytes     Size of the file in bytes
 * @param fileFormat        Format of the file (FASTQ, BAM, CRAM, VCF, GVCF, BED, 23ANDME, ANCESTRY, FTDNA)
 * @param checksum          SHA-256 or similar checksum for data integrity verification
 * @param checksumAlgorithm Algorithm used for checksum (SHA-256, MD5, CRC32)
 * @param location          User's personal reference to file location (local path, personal cloud).
 *                          DecodingUs does NOT access this URI.
 */
case class FileInfo(
                     fileName: String,
                     fileSizeBytes: Option[Long],
                     fileFormat: String,
                     checksum: Option[String],
                     checksumAlgorithm: Option[String] = None,
                     location: Option[String] = None
                   )
