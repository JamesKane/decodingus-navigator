package com.decodingus.analysis

import java.io.File
import java.nio.file.{Files, Path, Paths}

/**
 * Manages artifact storage organized by subject, sequence run, and alignment.
 *
 * Directory structure:
 * ~/.decodingus/cache/subjects/{sampleAccession}/runs/{runId}/alignments/{alignmentId}/
 *   ├── wgs_metrics.txt
 *   ├── callable_loci/
 *   │   ├── chr1.callable.bed
 *   │   ├── chr1.table.txt
 *   │   └── chr1.callable.svg
 *   └── haplogroup/
 *       ├── alleles.vcf
 *       └── called.vcf
 */
object SubjectArtifactCache {

  private val CacheRoot = Paths.get(System.getProperty("user.home"), ".decodingus", "cache", "subjects")

  /**
   * Gets the base directory for a subject's artifacts.
   * Creates the directory if it doesn't exist.
   */
  def getSubjectDir(sampleAccession: String): Path = {
    val dir = CacheRoot.resolve(sanitizePath(sampleAccession))
    ensureDir(dir)
    dir
  }

  /**
   * Gets the directory for a specific sequence run's artifacts.
   * The runId should be derived from the SequenceRun's atUri (last segment).
   */
  def getSequenceRunDir(sampleAccession: String, runId: String): Path = {
    val dir = getSubjectDir(sampleAccession).resolve("runs").resolve(sanitizePath(runId))
    ensureDir(dir)
    dir
  }

  /**
   * Gets the directory for a specific alignment's artifacts.
   * The alignmentId should be derived from the Alignment's atUri (last segment).
   */
  def getAlignmentDir(sampleAccession: String, runId: String, alignmentId: String): Path = {
    val dir = getSequenceRunDir(sampleAccession, runId).resolve("alignments").resolve(sanitizePath(alignmentId))
    ensureDir(dir)
    dir
  }

  /**
   * Gets a specific artifact path within an alignment directory.
   * Creates parent directories as needed.
   */
  def getArtifactPath(sampleAccession: String, runId: String, alignmentId: String, artifactName: String): Path = {
    val alignDir = getAlignmentDir(sampleAccession, runId, alignmentId)
    alignDir.resolve(artifactName)
  }

  /**
   * Gets a subdirectory for an artifact type (e.g., "callable_loci", "haplogroup").
   * Creates the directory if it doesn't exist.
   */
  def getArtifactSubdir(sampleAccession: String, runId: String, alignmentId: String, subdirName: String): Path = {
    val dir = getAlignmentDir(sampleAccession, runId, alignmentId).resolve(subdirName)
    ensureDir(dir)
    dir
  }

  /**
   * Extracts the ID portion from an AT URI.
   * e.g., "local:alignment:SAMPLE001:abc12345" -> "abc12345"
   *       "at://did:plc:123/com.decodingus.atmosphere.alignment/xyz" -> "xyz"
   */
  def extractIdFromUri(uri: String): String = {
    if (uri.contains(":")) {
      uri.split(":").last.split("/").last
    } else {
      uri.split("/").last
    }
  }

  /**
   * Convenience method to get artifact directory using AT URIs directly.
   */
  def getAlignmentDirFromUris(
    sampleAccession: String,
    sequenceRunUri: Option[String],
    alignmentUri: Option[String]
  ): Path = {
    val runId = sequenceRunUri.map(extractIdFromUri).getOrElse("unknown-run")
    val alignId = alignmentUri.map(extractIdFromUri).getOrElse("unknown-alignment")
    getAlignmentDir(sampleAccession, runId, alignId)
  }

  /**
   * Sanitizes a path component to be filesystem-safe.
   * Replaces problematic characters with underscores.
   */
  private def sanitizePath(name: String): String = {
    name.replaceAll("[^a-zA-Z0-9._-]", "_")
  }

  /**
   * Ensures a directory exists, creating it if necessary.
   */
  private def ensureDir(path: Path): Unit = {
    if (!Files.exists(path)) {
      Files.createDirectories(path)
    }
  }

  /**
   * Checks if an artifact exists at the given path.
   */
  def artifactExists(sampleAccession: String, runId: String, alignmentId: String, artifactName: String): Boolean = {
    val path = getArtifactPath(sampleAccession, runId, alignmentId, artifactName)
    Files.exists(path)
  }

  /**
   * Lists all artifacts in an alignment directory.
   */
  def listArtifacts(sampleAccession: String, runId: String, alignmentId: String): List[Path] = {
    val dir = getAlignmentDir(sampleAccession, runId, alignmentId)
    if (Files.exists(dir)) {
      import scala.jdk.CollectionConverters._
      Files.walk(dir).iterator().asScala.filter(Files.isRegularFile(_)).toList
    } else {
      List.empty
    }
  }

  /**
   * Deletes all artifacts for an alignment.
   */
  def deleteAlignmentArtifacts(sampleAccession: String, runId: String, alignmentId: String): Unit = {
    val dir = getAlignmentDir(sampleAccession, runId, alignmentId)
    if (Files.exists(dir)) {
      import scala.jdk.CollectionConverters._
      Files.walk(dir)
        .sorted(java.util.Comparator.reverseOrder())
        .iterator()
        .asScala
        .foreach(Files.delete)
    }
  }

  /**
   * Deletes all artifacts for a sequence run (including all alignments).
   */
  def deleteSequenceRunArtifacts(sampleAccession: String, runId: String): Unit = {
    val dir = getSequenceRunDir(sampleAccession, runId)
    if (Files.exists(dir)) {
      import scala.jdk.CollectionConverters._
      Files.walk(dir)
        .sorted(java.util.Comparator.reverseOrder())
        .iterator()
        .asScala
        .foreach(Files.delete)
    }
  }

  /**
   * Deletes all artifacts for a subject.
   */
  def deleteSubjectArtifacts(sampleAccession: String): Unit = {
    val dir = getSubjectDir(sampleAccession)
    if (Files.exists(dir)) {
      import scala.jdk.CollectionConverters._
      Files.walk(dir)
        .sorted(java.util.Comparator.reverseOrder())
        .iterator()
        .asScala
        .foreach(Files.delete)
    }
  }
}
