package com.decodingus.analysis

import com.decodingus.refgenome.{StrAnnotator, StrRegion}
import htsjdk.samtools.{SAMRecord, SamReaderFactory, ValidationStringency}

import java.io.File
import java.nio.file.Path
import scala.collection.mutable
import scala.jdk.CollectionConverters.*

/**
 * Result of calling an STR locus from BAM reads.
 *
 * @param chrom       Chromosome (e.g., "chrY")
 * @param start       Start position (0-based, from HipSTR reference)
 * @param end         End position (0-based exclusive)
 * @param period      Repeat unit length in bp
 * @param refRepeats  Reference repeat count
 * @param name        Locus name if known (e.g., "DYS393")
 * @param calledRepeats Called repeat count (None if no call)
 * @param confidence  Call confidence (0.0-1.0)
 * @param quality     Quality category: HIGH, MEDIUM, LOW, NO_CALL
 * @param readDepth   Number of spanning reads used for calling
 * @param alleleDistribution Distribution of observed repeat counts across reads
 * @param stutterFiltered Number of reads filtered as likely stutter
 */
case class StrCall(
  chrom: String,
  start: Long,
  end: Long,
  period: Int,
  refRepeats: Double,
  name: Option[String],
  calledRepeats: Option[Int],
  confidence: Double,
  quality: String,
  readDepth: Int,
  alleleDistribution: Map[Int, Int],
  stutterFiltered: Int
) {
  /** GRCh38 1-based start position for display/export */
  def start1Based: Long = start + 1

  /** Region span in base pairs */
  def regionSpan: Long = end - start

  /** Difference from reference repeat count */
  def deltaFromRef: Option[Double] = calledRepeats.map(_ - refRepeats)
}

/**
 * Configuration for STR calling.
 *
 * @param minReadDepth      Minimum spanning reads for a call
 * @param minMapQ           Minimum mapping quality for reads
 * @param minBaseQ          Minimum base quality in repeat region
 * @param stutterThreshold  Fraction of reads at which off-by-one is considered stutter
 * @param consensusThreshold Fraction of reads required for consensus call
 * @param flankingBases     Required flanking bases on each side of repeat
 */
case class StrCallerConfig(
  minReadDepth: Int = 5,
  minMapQ: Int = 20,
  minBaseQ: Int = 20,
  stutterThreshold: Double = 0.15,
  consensusThreshold: Double = 0.7,
  flankingBases: Int = 5
)

/**
 * STR caller for Y-chromosome loci using HTSJDK.
 *
 * Calls Short Tandem Repeat alleles from BAM/CRAM files by:
 * 1. Fetching reads spanning known STR loci from HipSTR reference
 * 2. Counting repeat units in each read
 * 3. Filtering stutter artifacts (off-by-one from PCR)
 * 4. Calling consensus allele with confidence scoring
 *
 * Designed for haploid Y-chromosome calling (simpler than diploid autosomal).
 */
class StrCaller(
  strAnnotator: StrAnnotator,
  config: StrCallerConfig = StrCallerConfig()
) {

  /**
   * Call STRs at all Y-chromosome loci from a BAM file.
   *
   * @param bamPath       Path to BAM/CRAM file
   * @param referencePath Path to reference FASTA (required for CRAM, recommended for all)
   * @param onProgress    Progress callback (message, current, total)
   * @return Either error message or list of STR calls
   */
  def callYChromosomeStrs(
    bamPath: String,
    referencePath: String,
    onProgress: (String, Long, Long) => Unit = (_, _, _) => ()
  ): Either[String, List[StrCall]] = {
    try {
      val samReaderFactory = SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)

      val samReader = if (bamPath.toLowerCase.endsWith(".cram")) {
        samReaderFactory.referenceSequence(new File(referencePath)).open(new File(bamPath))
      } else {
        samReaderFactory.open(new File(bamPath))
      }

      // Get Y chromosome contig name from BAM header
      val header = samReader.getFileHeader
      val yContig = findYContig(header.getSequenceDictionary.getSequences.asScala.map(_.getSequenceName).toList)

      if (yContig.isEmpty) {
        samReader.close()
        return Left("No Y chromosome contig found in BAM header")
      }

      val yChrom = yContig.get
      onProgress(s"Calling STRs on $yChrom...", 0, 1)

      // Get all STR regions for Y chromosome from annotator
      val yStrRegions = getYStrRegions(yChrom)

      if (yStrRegions.isEmpty) {
        samReader.close()
        return Left(s"No STR regions found for $yChrom in reference")
      }

      val totalLoci = yStrRegions.size
      var processed = 0
      val calls = mutable.ListBuffer[StrCall]()

      for (region <- yStrRegions) {
        processed += 1
        if (processed % 100 == 0) {
          onProgress(s"Processing locus $processed of $totalLoci...", processed, totalLoci)
        }

        // Query reads overlapping this STR region
        // Use 1-based coordinates for HTSJDK query
        val regionStart = region.start.toInt + 1
        val regionEnd = region.end.toInt

        val overlappingReads = samReader.query(yChrom, regionStart, regionEnd, false)
          .asScala
          .filter(isUsableRead)
          .toList

        val call = callLocus(region, overlappingReads, yChrom)
        calls += call
      }

      samReader.close()
      onProgress(s"Called ${calls.size} STR loci", totalLoci, totalLoci)

      Right(calls.toList)
    } catch {
      case e: Exception =>
        Left(s"STR calling failed: ${e.getMessage}")
    }
  }

  /**
   * Find the Y chromosome contig name (handles chrY, Y, NC_000024.10, etc.)
   */
  private def findYContig(contigs: List[String]): Option[String] = {
    contigs.find { c =>
      val normalized = c.toLowerCase
      normalized == "chry" || normalized == "y" || normalized.contains("nc_000024")
    }
  }

  /**
   * Get STR regions for Y chromosome from the annotator.
   * Tries multiple naming conventions.
   */
  private def getYStrRegions(yChrom: String): List[StrRegion] = {
    strAnnotator.getRegionsForChromosome(yChrom)
  }

  /**
   * Check if a read is usable for STR calling.
   */
  private def isUsableRead(record: SAMRecord): Boolean = {
    !record.getReadUnmappedFlag &&
    !record.isSecondaryOrSupplementary &&
    !record.getDuplicateReadFlag &&
    record.getMappingQuality >= config.minMapQ
  }

  /**
   * Call a single STR locus from overlapping reads.
   */
  private def callLocus(region: StrRegion, reads: List[SAMRecord], chrom: String): StrCall = {
    if (reads.size < config.minReadDepth) {
      return StrCall(
        chrom = chrom,
        start = region.start,
        end = region.end,
        period = region.period,
        refRepeats = region.numRepeats,
        name = region.name,
        calledRepeats = None,
        confidence = 0.0,
        quality = "NO_CALL",
        readDepth = reads.size,
        alleleDistribution = Map.empty,
        stutterFiltered = 0
      )
    }

    // Count repeat units in each read
    val repeatCounts = reads.flatMap { read =>
      countRepeatsInRead(read, region)
    }

    if (repeatCounts.isEmpty) {
      return StrCall(
        chrom = chrom,
        start = region.start,
        end = region.end,
        period = region.period,
        refRepeats = region.numRepeats,
        name = region.name,
        calledRepeats = None,
        confidence = 0.0,
        quality = "NO_CALL",
        readDepth = reads.size,
        alleleDistribution = Map.empty,
        stutterFiltered = 0
      )
    }

    // Build allele distribution
    val distribution = repeatCounts.groupBy(identity).view.mapValues(_.size).toMap

    // Find modal allele
    val (modalAllele, modalCount) = distribution.maxBy(_._2)
    val totalCounted = repeatCounts.size

    // Apply stutter filtering
    val (filteredDistribution, stutterCount) = filterStutter(distribution, modalAllele)

    // Calculate consensus
    val filteredTotal = filteredDistribution.values.sum
    val consensusFraction = if (filteredTotal > 0) {
      filteredDistribution.getOrElse(modalAllele, 0).toDouble / filteredTotal
    } else 0.0

    // Determine quality and confidence
    val (quality, confidence) = if (filteredTotal < config.minReadDepth) {
      ("LOW", consensusFraction * 0.5)
    } else if (consensusFraction >= config.consensusThreshold) {
      if (filteredTotal >= config.minReadDepth * 2) {
        ("HIGH", consensusFraction)
      } else {
        ("MEDIUM", consensusFraction * 0.9)
      }
    } else {
      ("LOW", consensusFraction * 0.7)
    }

    StrCall(
      chrom = chrom,
      start = region.start,
      end = region.end,
      period = region.period,
      refRepeats = region.numRepeats,
      name = region.name,
      calledRepeats = Some(modalAllele),
      confidence = confidence,
      quality = quality,
      readDepth = filteredTotal,
      alleleDistribution = filteredDistribution,
      stutterFiltered = stutterCount
    )
  }

  /**
   * Count repeat units in a read for a given STR region.
   *
   * This uses a simple approach:
   * 1. Find the read segment covering the STR region
   * 2. Calculate the observed length in the read (accounting for indels via CIGAR)
   * 3. Divide by repeat period to get repeat count
   */
  private def countRepeatsInRead(read: SAMRecord, region: StrRegion): Option[Int] = {
    val readStart = read.getAlignmentStart // 1-based
    val readEnd = read.getAlignmentEnd     // 1-based

    // Region coordinates (0-based in BED, convert to 1-based)
    val regionStart1 = region.start.toInt + 1
    val regionEnd1 = region.end.toInt

    // Check if read spans the entire region with flanking bases
    val requiredStart = regionStart1 - config.flankingBases
    val requiredEnd = regionEnd1 + config.flankingBases

    if (readStart > requiredStart || readEnd < requiredEnd) {
      return None // Read doesn't fully span the region
    }

    // Get the portion of the read aligned to the STR region
    // Use CIGAR to properly account for insertions/deletions
    val observedLength = getObservedLengthInRegion(read, regionStart1, regionEnd1)

    observedLength.map { len =>
      // Round to nearest integer repeat count
      Math.round(len.toDouble / region.period).toInt
    }
  }

  /**
   * Calculate the observed sequence length in a genomic region,
   * accounting for insertions and deletions via CIGAR parsing.
   */
  private def getObservedLengthInRegion(read: SAMRecord, regionStart: Int, regionEnd: Int): Option[Int] = {
    val cigar = read.getCigar
    if (cigar == null) return None

    var refPos = read.getAlignmentStart
    var observedBases = 0
    var inRegion = false
    var regionComplete = false

    for (element <- cigar.getCigarElements.asScala if !regionComplete) {
      val op = element.getOperator
      val length = element.getLength

      op.toString match {
        case "M" | "=" | "X" => // Match/mismatch - consumes both ref and read
          for (_ <- 0 until length if !regionComplete) {
            if (refPos >= regionStart && refPos <= regionEnd) {
              inRegion = true
              observedBases += 1
            } else if (inRegion) {
              regionComplete = true
            }
            refPos += 1
          }

        case "I" => // Insertion - consumes read only
          if (inRegion) {
            observedBases += length
          }

        case "D" => // Deletion - consumes ref only
          for (_ <- 0 until length if !regionComplete) {
            if (refPos >= regionStart && refPos <= regionEnd) {
              inRegion = true
              // Deletion - don't add to observed bases
            } else if (inRegion) {
              regionComplete = true
            }
            refPos += 1
          }

        case "N" => // Skipped region (intron) - consumes ref only
          refPos += length
          if (inRegion) {
            regionComplete = true // Can't span an intron
          }

        case "S" | "H" => // Soft/hard clip - doesn't affect alignment
          // S consumes read, H consumes neither

        case _ => // Other operators
      }
    }

    if (observedBases > 0) Some(observedBases) else None
  }

  /**
   * Filter stutter artifacts from allele distribution.
   *
   * Stutter typically appears as -1 repeat (occasionally -2 or +1)
   * from the true allele due to PCR slippage.
   */
  private def filterStutter(
    distribution: Map[Int, Int],
    modalAllele: Int
  ): (Map[Int, Int], Int) = {
    val totalReads = distribution.values.sum
    var stutterCount = 0

    val filtered = distribution.filter { case (allele, count) =>
      val diff = Math.abs(allele - modalAllele)
      val fraction = count.toDouble / totalReads

      // Consider it stutter if:
      // - Off by 1 or 2 repeats from modal
      // - Below stutter threshold
      val isStutter = diff > 0 && diff <= 2 && fraction < config.stutterThreshold

      if (isStutter) {
        stutterCount += count
        false
      } else {
        true
      }
    }

    (filtered, stutterCount)
  }
}

object StrCaller {
  /**
   * Create an STR caller for the given reference build.
   */
  def forBuild(referenceBuild: String, config: StrCallerConfig = StrCallerConfig()): Either[String, StrCaller] = {
    StrAnnotator.forBuild(referenceBuild).map { annotator =>
      new StrCaller(annotator, config)
    }
  }
}
