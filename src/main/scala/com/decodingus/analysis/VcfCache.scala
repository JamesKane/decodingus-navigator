package com.decodingus.analysis

import io.circe.{Codec, parser}
import io.circe.syntax._

import java.io.File
import java.nio.file.{Files, Path}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import scala.jdk.CollectionConverters._
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
 * Known vendor types for VCF imports.
 */
enum VcfVendor(val code: String, val displayName: String):
  case FTDNA_BIGY extends VcfVendor("ftdna_bigy", "FTDNA Big Y")
  case FTDNA_MTFULL extends VcfVendor("ftdna_mtfull", "FTDNA mtFull Sequence")
  case YSEQ extends VcfVendor("yseq", "YSEQ")
  case NEBULA extends VcfVendor("nebula", "Nebula Genomics")
  case DANTE extends VcfVendor("dante", "Dante Labs")
  case FULL_GENOMES extends VcfVendor("fullgenomes", "Full Genomes Corp")
  case OTHER extends VcfVendor("other", "Other")

object VcfVendor:
  def fromCode(code: String): Option[VcfVendor] =
    VcfVendor.values.find(_.code.equalsIgnoreCase(code))

  given Codec[VcfVendor] = Codec.from(
    io.circe.Decoder.decodeString.map(code => fromCode(code).getOrElse(OTHER)),
    io.circe.Encoder.encodeString.contramap(_.code)
  )

/**
 * Information about a vendor-provided VCF (e.g., FTDNA Big Y).
 *
 * Vendors like FTDNA provide:
 * - A VCF file with variant calls
 * - A BED file with target capture regions (where sequencing was performed)
 *
 * Note: The BED file defines the target/capture regions, NOT quality-filtered
 * callable regions. It indicates where the assay was designed to sequence,
 * not what regions achieved adequate coverage.
 *
 * @param vcfPath           Path to the vendor VCF file
 * @param indexPath         Path to the VCF index (.tbi)
 * @param targetBedPath     Optional path to target capture regions BED file
 * @param vendor            Vendor that provided the file
 * @param originalVcfName   Original filename of the VCF
 * @param originalBedName   Original filename of the BED (if provided)
 * @param referenceBuild    Reference genome build (GRCh37, GRCh38, etc.)
 * @param importedAt        When the file was imported
 * @param variantCount      Number of variants in the VCF
 * @param contigs           Contigs present in the VCF
 * @param notes             Optional notes about the import
 */
case class VendorVcfInfo(
  vcfPath: String,
  indexPath: String,
  targetBedPath: Option[String],
  vendor: VcfVendor,
  originalVcfName: String,
  originalBedName: Option[String],
  referenceBuild: String,
  importedAt: String,
  variantCount: Long,
  contigs: List[String],
  notes: Option[String] = None
) derives Codec.AsObject {

  def vcfFile: File = new File(vcfPath)
  def indexFile: File = new File(indexPath)
  def targetBedFile: Option[File] = targetBedPath.map(new File(_))
  def importedAtDateTime: LocalDateTime = LocalDateTime.parse(importedAt, DateTimeFormatter.ISO_LOCAL_DATE_TIME)

  def isValid: Boolean = {
    vcfFile.exists() && indexFile.exists() &&
      targetBedPath.forall(p => new File(p).exists())
  }

  /** Check if this is a Y-DNA vendor VCF */
  def isYDna: Boolean = contigs.exists(c =>
    c.equalsIgnoreCase("chrY") || c.equalsIgnoreCase("Y"))

  /** Check if this is an mtDNA vendor VCF */
  def isMtDna: Boolean = contigs.exists(c =>
    c.equalsIgnoreCase("chrM") || c.equalsIgnoreCase("MT") || c.equalsIgnoreCase("chrMT"))
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

  // --- Vendor VCF Support ---

  private val VENDOR_SUBDIR = "vendor"
  private val VENDOR_METADATA_SUFFIX = "_metadata.json"

  /**
   * Get the vendor VCF directory for an alignment.
   */
  def getVendorVcfDir(sampleAccession: String, runId: String, alignmentId: String): Path = {
    getVcfDir(sampleAccession, runId, alignmentId).resolve(VENDOR_SUBDIR)
  }

  /**
   * Import a vendor-provided VCF (and optional target regions BED) into the cache.
   *
   * The VCF will be indexed if not already indexed.
   *
   * @param sampleAccession Sample accession
   * @param runId Sequence run ID
   * @param alignmentId Alignment ID
   * @param vcfSourcePath Path to the source VCF file
   * @param bedSourcePath Optional path to the target regions BED file
   * @param vendor The vendor that provided the files
   * @param referenceBuild Reference genome build
   * @param notes Optional notes about this import
   * @return Either error message or the VendorVcfInfo for the imported files
   */
  def importVendorVcf(
    sampleAccession: String,
    runId: String,
    alignmentId: String,
    vcfSourcePath: Path,
    bedSourcePath: Option[Path],
    vendor: VcfVendor,
    referenceBuild: String,
    notes: Option[String] = None
  ): Either[String, VendorVcfInfo] = {
    try {
      val vendorDir = getVendorVcfDir(sampleAccession, runId, alignmentId)
      Files.createDirectories(vendorDir)

      val originalVcfName = vcfSourcePath.getFileName.toString
      val originalBedName = bedSourcePath.map(_.getFileName.toString)

      // Determine destination filename based on vendor
      val baseFilename = vendor.code
      val vcfDest = vendorDir.resolve(s"$baseFilename.vcf.gz")
      val indexDest = vendorDir.resolve(s"$baseFilename.vcf.gz.tbi")
      val bedDest = bedSourcePath.map(_ => vendorDir.resolve(s"${baseFilename}_targets.bed"))

      // Copy or compress the VCF
      val vcfIsGzipped = originalVcfName.endsWith(".gz")
      if (vcfIsGzipped) {
        Files.copy(vcfSourcePath, vcfDest, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
      } else {
        // Compress using GATK SortVcf which also handles bgzip
        GatkRunner.run(Array(
          "SortVcf",
          "-I", vcfSourcePath.toString,
          "-O", vcfDest.toString,
          "--CREATE_INDEX", "true"
        )) match {
          case Left(error) => return Left(s"Failed to compress VCF: $error")
          case Right(_) =>
        }
      }

      // Create index if needed
      if (vcfIsGzipped && !Files.exists(indexDest)) {
        // Check for existing index
        val sourceIndex = vcfSourcePath.resolveSibling(originalVcfName + ".tbi")
        if (Files.exists(sourceIndex)) {
          Files.copy(sourceIndex, indexDest, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
        } else {
          // Create index using GATK
          GatkRunner.run(Array(
            "IndexFeatureFile",
            "-I", vcfDest.toString
          )) match {
            case Left(error) => return Left(s"Failed to index VCF: $error")
            case Right(_) =>
          }
        }
      }

      // Copy BED file if provided
      bedSourcePath.zip(bedDest).foreach { case (src, dest) =>
        Files.copy(src, dest, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
      }

      // Extract contig list and count variants
      val (contigs, variantCount) = extractVcfInfo(vcfDest)

      // Create metadata
      val metadata = VendorVcfInfo(
        vcfPath = vcfDest.toAbsolutePath.toString,
        indexPath = indexDest.toAbsolutePath.toString,
        targetBedPath = bedDest.map(_.toAbsolutePath.toString),
        vendor = vendor,
        originalVcfName = originalVcfName,
        originalBedName = originalBedName,
        referenceBuild = referenceBuild,
        importedAt = LocalDateTime.now().format(DateTimeFormatter.ISO_LOCAL_DATE_TIME),
        variantCount = variantCount,
        contigs = contigs,
        notes = notes
      )

      // Save metadata
      val metadataPath = vendorDir.resolve(s"$baseFilename$VENDOR_METADATA_SUFFIX")
      Files.writeString(metadataPath, metadata.asJson.spaces2)

      Right(metadata)

    } catch {
      case e: Exception => Left(s"Failed to import vendor VCF: ${e.getMessage}")
    }
  }

  /**
   * Extract contig list and variant count from a VCF file.
   */
  private def extractVcfInfo(vcfPath: Path): (List[String], Long) = {
    import htsjdk.variant.vcf.VCFFileReader
    val reader = new VCFFileReader(vcfPath, false)
    try {
      val contigs = reader.getFileHeader.getContigLines.asScala.map(_.getID).toList
      var count = 0L
      val iter = reader.iterator()
      while (iter.hasNext) {
        iter.next()
        count += 1
      }
      (contigs, count)
    } finally {
      reader.close()
    }
  }

  /**
   * List all vendor VCFs imported for an alignment.
   */
  def listVendorVcfs(sampleAccession: String, runId: String, alignmentId: String): List[VendorVcfInfo] = {
    val vendorDir = getVendorVcfDir(sampleAccession, runId, alignmentId)

    if (!Files.exists(vendorDir)) {
      return List.empty
    }

    import scala.jdk.CollectionConverters._
    Files.list(vendorDir).iterator().asScala
      .filter(p => p.toString.endsWith(VENDOR_METADATA_SUFFIX))
      .flatMap { metadataPath =>
        Using(scala.io.Source.fromFile(metadataPath.toFile)) { source =>
          parser.decode[VendorVcfInfo](source.mkString).toOption
        }.toOption.flatten
      }
      .filter(_.isValid)
      .toList
  }

  /**
   * Load a specific vendor VCF by vendor type.
   */
  def loadVendorVcf(
    sampleAccession: String,
    runId: String,
    alignmentId: String,
    vendor: VcfVendor
  ): Option[VendorVcfInfo] = {
    val vendorDir = getVendorVcfDir(sampleAccession, runId, alignmentId)
    val metadataPath = vendorDir.resolve(s"${vendor.code}$VENDOR_METADATA_SUFFIX")

    if (!Files.exists(metadataPath)) {
      return None
    }

    Using(scala.io.Source.fromFile(metadataPath.toFile)) { source =>
      parser.decode[VendorVcfInfo](source.mkString).toOption.filter(_.isValid)
    }.toOption.flatten
  }

  /**
   * Find the best vendor VCF for Y-DNA haplogroup analysis.
   * Prefers FTDNA Big Y, then other Y-DNA specific vendors.
   */
  def findYDnaVendorVcf(sampleAccession: String, runId: String, alignmentId: String): Option[VendorVcfInfo] = {
    val vendorVcfs = listVendorVcfs(sampleAccession, runId, alignmentId)
    val yDnaVcfs = vendorVcfs.filter(_.isYDna)

    // Priority order for Y-DNA
    val priorityOrder = List(VcfVendor.FTDNA_BIGY, VcfVendor.YSEQ, VcfVendor.FULL_GENOMES)

    priorityOrder.flatMap(vendor => yDnaVcfs.find(_.vendor == vendor)).headOption
      .orElse(yDnaVcfs.headOption)
  }

  /**
   * Find the best vendor VCF for mtDNA haplogroup analysis.
   */
  def findMtDnaVendorVcf(sampleAccession: String, runId: String, alignmentId: String): Option[VendorVcfInfo] = {
    val vendorVcfs = listVendorVcfs(sampleAccession, runId, alignmentId)
    val mtDnaVcfs = vendorVcfs.filter(_.isMtDna)

    // Priority order for mtDNA
    val priorityOrder = List(VcfVendor.FTDNA_MTFULL, VcfVendor.NEBULA, VcfVendor.DANTE)

    priorityOrder.flatMap(vendor => mtDnaVcfs.find(_.vendor == vendor)).headOption
      .orElse(mtDnaVcfs.headOption)
  }

  /**
   * Delete a vendor VCF from the cache.
   */
  def deleteVendorVcf(
    sampleAccession: String,
    runId: String,
    alignmentId: String,
    vendor: VcfVendor
  ): Boolean = {
    val vendorDir = getVendorVcfDir(sampleAccession, runId, alignmentId)
    val baseFilename = vendor.code

    var deleted = false
    List(
      vendorDir.resolve(s"$baseFilename.vcf.gz"),
      vendorDir.resolve(s"$baseFilename.vcf.gz.tbi"),
      vendorDir.resolve(s"${baseFilename}_targets.bed"),
      vendorDir.resolve(s"$baseFilename$VENDOR_METADATA_SUFFIX")
    ).foreach { path =>
      if (Files.exists(path)) {
        Files.delete(path)
        deleted = true
      }
    }
    deleted
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
