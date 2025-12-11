package com.decodingus.analysis

import io.circe.{Codec, parser}
import io.circe.syntax._

import java.io.File
import java.nio.file.{Files, Path}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import scala.util.{Either, Left, Right, Using}

/**
 * Information about a cached whole-genome VCF.
 */
case class CachedVcfInfo(
  vcfPath: String,
  indexPath: String,
  referenceBuild: String,
  callerVersion: String,
  gatkVersion: String,
  createdAt: String,
  fileSizeBytes: Long,
  variantCount: Long,
  contigs: List[String],
  inferredSex: Option[String] = None
) derives Codec.AsObject {

  def vcfFile: File = new File(vcfPath)
  def indexFile: File = new File(indexPath)
  def createdAtDateTime: LocalDateTime = LocalDateTime.parse(createdAt, DateTimeFormatter.ISO_LOCAL_DATE_TIME)

  def isValid: Boolean = {
    vcfFile.exists() && indexFile.exists()
  }
}

/**
 * Manages cached VCF storage and retrieval for whole-genome variant calling.
 *
 * VCF files are stored in:
 * ~/.decodingus/cache/subjects/{sample}/runs/{run}/alignments/{align}/vcf/
 *   - whole_genome.vcf.gz
 *   - whole_genome.vcf.gz.tbi
 *   - vcf_metadata.json
 */
object VcfCache {

  private val VCF_SUBDIR = "vcf"
  private val VCF_FILENAME = "whole_genome.vcf.gz"
  private val INDEX_FILENAME = "whole_genome.vcf.gz.tbi"
  private val METADATA_FILENAME = "vcf_metadata.json"

  /**
   * Get the VCF directory for an alignment.
   */
  def getVcfDir(sampleAccession: String, runId: String, alignmentId: String): Path = {
    SubjectArtifactCache.getArtifactSubdir(sampleAccession, runId, alignmentId, VCF_SUBDIR)
  }

  /**
   * Get the VCF file path for an alignment.
   */
  def getVcfPath(sampleAccession: String, runId: String, alignmentId: String): Path = {
    getVcfDir(sampleAccession, runId, alignmentId).resolve(VCF_FILENAME)
  }

  /**
   * Get the VCF index file path for an alignment.
   */
  def getIndexPath(sampleAccession: String, runId: String, alignmentId: String): Path = {
    getVcfDir(sampleAccession, runId, alignmentId).resolve(INDEX_FILENAME)
  }

  /**
   * Get the metadata file path for an alignment.
   */
  def getMetadataPath(sampleAccession: String, runId: String, alignmentId: String): Path = {
    getVcfDir(sampleAccession, runId, alignmentId).resolve(METADATA_FILENAME)
  }

  /**
   * Check if a cached VCF exists for an alignment.
   */
  def exists(sampleAccession: String, runId: String, alignmentId: String): Boolean = {
    val vcfPath = getVcfPath(sampleAccession, runId, alignmentId)
    val indexPath = getIndexPath(sampleAccession, runId, alignmentId)
    Files.exists(vcfPath) && Files.exists(indexPath)
  }

  /**
   * Load cached VCF metadata for an alignment.
   */
  def loadMetadata(sampleAccession: String, runId: String, alignmentId: String): Either[String, CachedVcfInfo] = {
    val metadataPath = getMetadataPath(sampleAccession, runId, alignmentId)

    if (!Files.exists(metadataPath)) {
      return Left(s"VCF metadata not found: $metadataPath")
    }

    Using(scala.io.Source.fromFile(metadataPath.toFile)) { source =>
      val json = source.mkString
      parser.decode[CachedVcfInfo](json) match {
        case Right(info) =>
          if (info.isValid) Right(info)
          else Left(s"VCF files missing or invalid at ${info.vcfPath}")
        case Left(error) =>
          Left(s"Failed to parse VCF metadata: ${error.getMessage}")
      }
    }.fold(
      error => Left(s"Failed to read VCF metadata: ${error.getMessage}"),
      identity
    )
  }

  /**
   * Load cached VCF info using AT URIs.
   */
  def loadMetadataFromUris(
    sampleAccession: String,
    sequenceRunUri: Option[String],
    alignmentUri: Option[String]
  ): Either[String, CachedVcfInfo] = {
    val runId = sequenceRunUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-run")
    val alignId = alignmentUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-alignment")
    loadMetadata(sampleAccession, runId, alignId)
  }

  /**
   * Save VCF metadata after generation.
   */
  def saveMetadata(
    sampleAccession: String,
    runId: String,
    alignmentId: String,
    info: CachedVcfInfo
  ): Either[String, Unit] = {
    val metadataPath = getMetadataPath(sampleAccession, runId, alignmentId)

    try {
      Files.writeString(metadataPath, info.asJson.spaces2)
      Right(())
    } catch {
      case e: Exception => Left(s"Failed to save VCF metadata: ${e.getMessage}")
    }
  }

  /**
   * Create VCF metadata from a generated VCF file.
   */
  def createMetadata(
    vcfPath: Path,
    indexPath: Path,
    referenceBuild: String,
    callerVersion: String,
    gatkVersion: String,
    contigs: List[String],
    variantCount: Long,
    inferredSex: Option[String] = None
  ): CachedVcfInfo = {
    CachedVcfInfo(
      vcfPath = vcfPath.toAbsolutePath.toString,
      indexPath = indexPath.toAbsolutePath.toString,
      referenceBuild = referenceBuild,
      callerVersion = callerVersion,
      gatkVersion = gatkVersion,
      createdAt = LocalDateTime.now().format(DateTimeFormatter.ISO_LOCAL_DATE_TIME),
      fileSizeBytes = Files.size(vcfPath),
      variantCount = variantCount,
      contigs = contigs,
      inferredSex = inferredSex
    )
  }

  /**
   * Delete cached VCF for an alignment.
   */
  def delete(sampleAccession: String, runId: String, alignmentId: String): Unit = {
    val vcfDir = getVcfDir(sampleAccession, runId, alignmentId)

    if (Files.exists(vcfDir)) {
      import scala.jdk.CollectionConverters._
      Files.walk(vcfDir)
        .sorted(java.util.Comparator.reverseOrder())
        .iterator()
        .asScala
        .foreach(Files.delete)
    }
  }

  /**
   * Get a summary of the cached VCF status for display.
   */
  def getStatus(sampleAccession: String, runId: String, alignmentId: String): VcfStatus = {
    loadMetadata(sampleAccession, runId, alignmentId) match {
      case Right(info) =>
        VcfStatus.Available(info)
      case Left(_) =>
        if (exists(sampleAccession, runId, alignmentId)) {
          // VCF exists but metadata is missing/invalid
          VcfStatus.Incomplete
        } else {
          VcfStatus.NotGenerated
        }
    }
  }

  /**
   * Validate that a cached VCF matches the expected reference build.
   */
  def validateBuild(
    sampleAccession: String,
    runId: String,
    alignmentId: String,
    expectedBuild: String
  ): Either[String, CachedVcfInfo] = {
    loadMetadata(sampleAccession, runId, alignmentId).flatMap { info =>
      if (info.referenceBuild == expectedBuild) {
        Right(info)
      } else {
        Left(s"VCF reference build mismatch: expected $expectedBuild, found ${info.referenceBuild}")
      }
    }
  }
}

/**
 * Status of a cached VCF.
 */
enum VcfStatus:
  case Available(info: CachedVcfInfo)
  case InProgress(startedAt: LocalDateTime, progress: Double, currentContig: Option[String])
  case NotGenerated
  case Incomplete  // Files exist but metadata is invalid
  case Stale       // VCF exists but alignment has been modified since

  def isAvailable: Boolean = this match {
    case Available(_) => true
    case _ => false
  }

  def isInProgress: Boolean = this match {
    case InProgress(_, _, _) => true
    case _ => false
  }
