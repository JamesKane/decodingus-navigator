package com.decodingus.analysis

import io.circe.syntax.*
import io.circe.{Codec, parser}

import java.io.File
import java.nio.file.{Files, Path}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import scala.jdk.CollectionConverters.*
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
 * @param vcfPath         Path to the vendor VCF file
 * @param indexPath       Path to the VCF index (.tbi)
 * @param targetBedPath   Optional path to target capture regions BED file
 * @param vendor          Vendor that provided the file
 * @param originalVcfName Original filename of the VCF
 * @param originalBedName Original filename of the BED (if provided)
 * @param referenceBuild  Reference genome build (GRCh37, GRCh38, etc.)
 * @param importedAt      When the file was imported
 * @param variantCount    Number of variants in the VCF
 * @param contigs         Contigs present in the VCF
 * @param notes           Optional notes about the import
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
 * Information about a vendor-provided mtDNA FASTA sequence.
 *
 * Vendors like FTDNA (mtFull Sequence) and YSEQ provide complete mitochondrial
 * genome sequences as FASTA files. These need to be compared against rCRS
 * (revised Cambridge Reference Sequence) to identify variants.
 *
 * @param fastaPath        Path to the FASTA file
 * @param vendor           Vendor that provided the file
 * @param originalFileName Original filename of the FASTA
 * @param importedAt       When the file was imported
 * @param sequenceLength   Length of the mtDNA sequence
 * @param notes            Optional notes about the import
 */
case class VendorFastaInfo(
                            fastaPath: String,
                            vendor: VcfVendor,
                            originalFileName: String,
                            importedAt: String,
                            sequenceLength: Int,
                            notes: Option[String] = None
                          ) derives Codec.AsObject {

  def fastaFile: File = new File(fastaPath)

  def importedAtDateTime: LocalDateTime = LocalDateTime.parse(importedAt, DateTimeFormatter.ISO_LOCAL_DATE_TIME)

  def isValid: Boolean = fastaFile.exists()
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
      import scala.jdk.CollectionConverters.*
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
   * @param runId           Sequence run ID
   * @param alignmentId     Alignment ID
   * @param vcfSourcePath   Path to the source VCF file
   * @param bedSourcePath   Optional path to the target regions BED file
   * @param vendor          The vendor that provided the files
   * @param referenceBuild  Reference genome build
   * @param notes           Optional notes about this import
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

    import scala.jdk.CollectionConverters.*
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

  // --- SequenceRun-Level Vendor VCF Support (for VCF-only imports, no alignment) ---

  /**
   * Get the vendor VCF directory for a sequence run (no alignment required).
   * Used for vendor VCFs like FTDNA Big Y that don't have associated BAM files.
   */
  def getRunVendorVcfDir(sampleAccession: String, runId: String): Path = {
    SubjectArtifactCache.getSequenceRunDir(sampleAccession, runId).resolve(VENDOR_SUBDIR)
  }

  /**
   * Import a vendor-provided VCF at the SequenceRun level (no alignment required).
   * This is for vendor deliverables like FTDNA Big Y that don't include BAM files.
   *
   * @param sampleAccession Sample accession
   * @param runId           Sequence run ID
   * @param vcfSourcePath   Path to the source VCF file
   * @param bedSourcePath   Optional path to the target regions BED file
   * @param vendor          The vendor that provided the files
   * @param referenceBuild  Reference genome build
   * @param notes           Optional notes about this import
   * @return Either error message or the VendorVcfInfo for the imported files
   */
  def importRunVendorVcf(
                          sampleAccession: String,
                          runId: String,
                          vcfSourcePath: Path,
                          bedSourcePath: Option[Path],
                          vendor: VcfVendor,
                          referenceBuild: String,
                          notes: Option[String] = None
                        ): Either[String, VendorVcfInfo] = {
    try {
      val vendorDir = getRunVendorVcfDir(sampleAccession, runId)
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
        val sourceIndex = vcfSourcePath.resolveSibling(originalVcfName + ".tbi")
        if (Files.exists(sourceIndex)) {
          Files.copy(sourceIndex, indexDest, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
        } else {
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

      // Write metadata
      val metadataPath = vendorDir.resolve(s"$baseFilename$VENDOR_METADATA_SUFFIX")
      Files.writeString(metadataPath, metadata.asJson.spaces2)

      Right(metadata)
    } catch {
      case e: Exception =>
        Left(s"Failed to import vendor VCF: ${e.getMessage}")
    }
  }

  /**
   * List all vendor VCFs imported at the SequenceRun level.
   */
  def listRunVendorVcfs(sampleAccession: String, runId: String): List[VendorVcfInfo] = {
    val vendorDir = getRunVendorVcfDir(sampleAccession, runId)

    if (!Files.exists(vendorDir)) {
      return List.empty
    }

    import scala.jdk.CollectionConverters.*
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
   * Load a specific vendor VCF by vendor type at the SequenceRun level.
   */
  def loadRunVendorVcf(
                        sampleAccession: String,
                        runId: String,
                        vendor: VcfVendor
                      ): Option[VendorVcfInfo] = {
    val vendorDir = getRunVendorVcfDir(sampleAccession, runId)
    val metadataPath = vendorDir.resolve(s"${vendor.code}$VENDOR_METADATA_SUFFIX")

    if (!Files.exists(metadataPath)) {
      return None
    }

    Using(scala.io.Source.fromFile(metadataPath.toFile)) { source =>
      parser.decode[VendorVcfInfo](source.mkString).toOption.filter(_.isValid)
    }.toOption.flatten
  }

  /**
   * Find the best vendor VCF for Y-DNA haplogroup analysis at the SequenceRun level.
   */
  def findYDnaRunVendorVcf(sampleAccession: String, runId: String): Option[VendorVcfInfo] = {
    val vendorVcfs = listRunVendorVcfs(sampleAccession, runId)
    val yDnaVcfs = vendorVcfs.filter(_.isYDna)

    val priorityOrder = List(VcfVendor.FTDNA_BIGY, VcfVendor.YSEQ, VcfVendor.FULL_GENOMES)
    priorityOrder.flatMap(vendor => yDnaVcfs.find(_.vendor == vendor)).headOption
      .orElse(yDnaVcfs.headOption)
  }

  /**
   * Find the best vendor VCF for mtDNA haplogroup analysis at the SequenceRun level.
   */
  def findMtDnaRunVendorVcf(sampleAccession: String, runId: String): Option[VendorVcfInfo] = {
    val vendorVcfs = listRunVendorVcfs(sampleAccession, runId)
    val mtDnaVcfs = vendorVcfs.filter(_.isMtDna)

    val priorityOrder = List(VcfVendor.FTDNA_MTFULL, VcfVendor.NEBULA, VcfVendor.DANTE)
    priorityOrder.flatMap(vendor => mtDnaVcfs.find(_.vendor == vendor)).headOption
      .orElse(mtDnaVcfs.headOption)
  }

  /**
   * Delete a vendor VCF from the run-level cache.
   */
  def deleteRunVendorVcf(
                          sampleAccession: String,
                          runId: String,
                          vendor: VcfVendor
                        ): Boolean = {
    val vendorDir = getRunVendorVcfDir(sampleAccession, runId)
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

  // --- Vendor mtDNA FASTA Support ---

  private val FASTA_SUBDIR = "fasta"
  private val FASTA_METADATA_SUFFIX = "_fasta_metadata.json"

  /**
   * Get the vendor FASTA directory for a sequence run.
   */
  def getRunFastaDir(sampleAccession: String, runId: String): Path = {
    SubjectArtifactCache.getSequenceRunDir(sampleAccession, runId).resolve(FASTA_SUBDIR)
  }

  /**
   * Import a vendor-provided mtDNA FASTA file at the SequenceRun level.
   *
   * @param sampleAccession Sample accession
   * @param runId           Sequence run ID
   * @param fastaSourcePath Path to the source FASTA file
   * @param vendor          The vendor that provided the file
   * @param notes           Optional notes about this import
   * @return Either error message or the VendorFastaInfo for the imported file
   */
  def importRunFasta(
                      sampleAccession: String,
                      runId: String,
                      fastaSourcePath: Path,
                      vendor: VcfVendor,
                      notes: Option[String] = None
                    ): Either[String, VendorFastaInfo] = {
    try {
      val fastaDir = getRunFastaDir(sampleAccession, runId)
      Files.createDirectories(fastaDir)

      val originalFileName = fastaSourcePath.getFileName.toString
      val baseFilename = vendor.code
      val fastaDest = fastaDir.resolve(s"$baseFilename.fasta")

      // Copy the FASTA file
      Files.copy(fastaSourcePath, fastaDest, java.nio.file.StandardCopyOption.REPLACE_EXISTING)

      // Read sequence length
      val sequenceLength = readFastaSequenceLength(fastaDest)

      // Create metadata
      val metadata = VendorFastaInfo(
        fastaPath = fastaDest.toAbsolutePath.toString,
        vendor = vendor,
        originalFileName = originalFileName,
        importedAt = LocalDateTime.now().format(DateTimeFormatter.ISO_LOCAL_DATE_TIME),
        sequenceLength = sequenceLength,
        notes = notes
      )

      // Write metadata
      val metadataPath = fastaDir.resolve(s"$baseFilename$FASTA_METADATA_SUFFIX")
      Files.writeString(metadataPath, metadata.asJson.spaces2)

      Right(metadata)
    } catch {
      case e: Exception =>
        Left(s"Failed to import vendor FASTA: ${e.getMessage}")
    }
  }

  /**
   * Read the sequence length from a FASTA file.
   */
  private def readFastaSequenceLength(fastaPath: Path): Int = {
    Using(scala.io.Source.fromFile(fastaPath.toFile)) { source =>
      source.getLines()
        .filterNot(_.startsWith(">")) // Skip header lines
        .map(_.trim.length)
        .sum
    }.getOrElse(0)
  }

  /**
   * List all vendor FASTA files imported at the SequenceRun level.
   */
  def listRunFastas(sampleAccession: String, runId: String): List[VendorFastaInfo] = {
    val fastaDir = getRunFastaDir(sampleAccession, runId)

    if (!Files.exists(fastaDir)) {
      return List.empty
    }

    import scala.jdk.CollectionConverters.*
    Files.list(fastaDir).iterator().asScala
      .filter(p => p.toString.endsWith(FASTA_METADATA_SUFFIX))
      .flatMap { metadataPath =>
        Using(scala.io.Source.fromFile(metadataPath.toFile)) { source =>
          parser.decode[VendorFastaInfo](source.mkString).toOption
        }.toOption.flatten
      }
      .filter(_.isValid)
      .toList
  }

  /**
   * Load a specific vendor FASTA by vendor type at the SequenceRun level.
   */
  def loadRunFasta(
                    sampleAccession: String,
                    runId: String,
                    vendor: VcfVendor
                  ): Option[VendorFastaInfo] = {
    val fastaDir = getRunFastaDir(sampleAccession, runId)
    val metadataPath = fastaDir.resolve(s"${vendor.code}$FASTA_METADATA_SUFFIX")

    if (!Files.exists(metadataPath)) {
      return None
    }

    Using(scala.io.Source.fromFile(metadataPath.toFile)) { source =>
      parser.decode[VendorFastaInfo](source.mkString).toOption.filter(_.isValid)
    }.toOption.flatten
  }

  /**
   * Find the best vendor FASTA for mtDNA haplogroup analysis at the SequenceRun level.
   */
  def findMtDnaRunFasta(sampleAccession: String, runId: String): Option[VendorFastaInfo] = {
    val fastas = listRunFastas(sampleAccession, runId)

    // Priority order for mtDNA FASTA
    val priorityOrder = List(VcfVendor.FTDNA_MTFULL, VcfVendor.YSEQ, VcfVendor.OTHER)
    priorityOrder.flatMap(vendor => fastas.find(_.vendor == vendor)).headOption
      .orElse(fastas.headOption)
  }

  /**
   * Delete a vendor FASTA from the run-level cache.
   */
  def deleteRunFasta(
                      sampleAccession: String,
                      runId: String,
                      vendor: VcfVendor
                    ): Boolean = {
    val fastaDir = getRunFastaDir(sampleAccession, runId)
    val baseFilename = vendor.code

    var deleted = false
    List(
      fastaDir.resolve(s"$baseFilename.fasta"),
      fastaDir.resolve(s"$baseFilename$FASTA_METADATA_SUFFIX")
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
  case Incomplete // Files exist but metadata is invalid
  case Stale // VCF exists but alignment has been modified since

  def isAvailable: Boolean = this match {
    case Available(_) => true
    case _ => false
  }

  def isInProgress: Boolean = this match {
    case InProgress(_, _, _) => true
    case _ => false
  }
