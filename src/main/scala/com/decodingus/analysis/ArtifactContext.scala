package com.decodingus.analysis

import java.nio.file.Path

/**
 * Context for organizing analysis artifacts by subject/run/alignment.
 */
case class ArtifactContext(
                            sampleAccession: String,
                            sequenceRunUri: Option[String],
                            alignmentUri: Option[String]
                          ) {
  /** Gets the artifact directory for this context */
  def getArtifactDir: Path = SubjectArtifactCache.getAlignmentDirFromUris(sampleAccession, sequenceRunUri, alignmentUri)

  /** Gets a subdirectory for a specific artifact type */
  def getSubdir(name: String): Path = {
    val runId = sequenceRunUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-run")
    val alignId = alignmentUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-alignment")
    SubjectArtifactCache.getArtifactSubdir(sampleAccession, runId, alignId, name)
  }
}
