package com.decodingus.analysis.sv

import htsjdk.samtools.{SAMRecord, SamReaderFactory, ValidationStringency}

import java.io.File
import scala.collection.mutable
import scala.jdk.CollectionConverters.*

/**
 * Single-pass BAM walker that collects structural variant evidence.
 *
 * Follows the same pattern as UnifiedMetricsWalker but collects:
 * - Per-bin read depth (for CNV detection)
 * - Discordant read pairs (for SV breakpoint detection)
 * - Split reads (for precise breakpoint localization)
 *
 * Uses HTSJDK SamReader directly for efficient traversal.
 *
 * References:
 * - Discordant pair detection: Chen et al. "BreakDancer: an algorithm for high-resolution
 *   mapping of genomic structural variation." Nature Methods 6.9 (2009): 677-681.
 *   https://doi.org/10.1038/nmeth.1363
 *
 * - Split-read analysis: Ye et al. "Pindel: a pattern growth approach to detect break points
 *   of large deletions and medium sized insertions from paired-end short reads."
 *   Bioinformatics 25.21 (2009): 2865-2871.
 *   https://doi.org/10.1093/bioinformatics/btp394
 *
 * - Combined evidence approach: Layer et al. "LUMPY: a probabilistic framework for
 *   structural variant discovery." Genome Biology 15.6 (2014): R84.
 *   https://doi.org/10.1186/gb-2014-15-6-r84
 */
class SvEvidenceWalker(config: SvCallerConfig = SvCallerConfig.default) {

  private val BIN_SIZE = config.binSize

  /**
   * Collect SV evidence from a BAM/CRAM file in a single pass.
   *
   * @param bamPath           Path to BAM/CRAM file
   * @param referencePath     Path to reference genome (required for CRAM)
   * @param contigLengths     Map of contig names to lengths
   * @param expectedInsertSize Expected insert size (from library metrics)
   * @param insertSizeSd      Insert size standard deviation
   * @param onProgress        Progress callback (message, current, total)
   * @return Either error or SvEvidenceCollection
   */
  def collectEvidence(
    bamPath: String,
    referencePath: String,
    contigLengths: Map[String, Long],
    expectedInsertSize: Double,
    insertSizeSd: Double,
    onProgress: (String, Long, Long) => Unit
  ): Either[String, SvEvidenceCollection] = {
    try {
      val samReaderFactory = SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)

      val samReader = if (bamPath.toLowerCase.endsWith(".cram")) {
        samReaderFactory.referenceSequence(new File(referencePath)).open(new File(bamPath))
      } else {
        samReaderFactory.open(new File(bamPath))
      }

      onProgress("Collecting SV evidence from BAM/CRAM...", 0, 1)

      // Get sample name from BAM header
      val sampleName = samReader.getFileHeader.getReadGroups.asScala
        .headOption
        .map(_.getSample)
        .getOrElse("unknown")

      // Initialize depth bins for each contig
      val depthBins = mutable.Map[String, Array[Int]]()
      contigLengths.foreach { case (contig, length) =>
        val numBins = ((length + BIN_SIZE - 1) / BIN_SIZE).toInt
        depthBins(contig) = new Array[Int](numBins)
      }

      // Evidence accumulators
      val discordantPairs = mutable.ListBuffer[DiscordantPair]()
      val splitReads = mutable.ListBuffer[SplitRead]()

      // For tracking read pairs (to get mate info for discordant detection)
      val seenReads = mutable.Set[String]()

      val iterator = samReader.iterator()
      var recordCount = 0L
      val progressInterval = 1000000L

      // Calculate insert size thresholds
      val insertSizeMax = expectedInsertSize + (config.insertSizeZThreshold * insertSizeSd)
      val insertSizeMin = math.max(0, expectedInsertSize - (config.insertSizeZThreshold * insertSizeSd))

      while (iterator.hasNext) {
        val record = iterator.next()
        recordCount += 1

        if (recordCount % progressInterval == 0) {
          onProgress(s"Processed ${recordCount / 1000000}M reads...", recordCount, recordCount + progressInterval)
        }

        // Skip unmapped reads entirely
        if (!record.getReadUnmappedFlag) {
          val contig = record.getContig

          // 1. Depth tracking (primary and non-supplementary alignments only)
          if (!record.isSecondaryOrSupplementary && depthBins.contains(contig)) {
            val binIndex = (record.getAlignmentStart / BIN_SIZE).toInt
            if (binIndex >= 0 && binIndex < depthBins(contig).length) {
              depthBins(contig)(binIndex) += 1
            }
          }

          // 2. Discordant pair detection (primary alignments only)
          if (!record.isSecondaryOrSupplementary && record.getReadPairedFlag) {
            detectDiscordantPair(
              record, expectedInsertSize, insertSizeMin, insertSizeMax
            ).foreach { dp =>
              discordantPairs += dp
            }
          }

          // 3. Split read detection (look for SA tag)
          if (record.getMappingQuality >= config.minMapQ) {
            extractSplitRead(record).foreach { sr =>
              splitReads += sr
            }
          }
        }
      }

      samReader.close()
      onProgress("SV evidence collection complete.", recordCount, recordCount)

      Right(SvEvidenceCollection(
        discordantPairs = discordantPairs.toList,
        splitReads = splitReads.toList,
        depthBins = depthBins.map { case (k, v) => k -> v }.toMap,
        sampleName = sampleName,
        expectedInsertSize = expectedInsertSize,
        insertSizeSd = insertSizeSd
      ))

    } catch {
      case e: Exception =>
        Left(s"Failed to collect SV evidence: ${e.getMessage}")
    }
  }

  /**
   * Detect if a read pair is discordant.
   *
   * A pair is discordant if:
   * 1. Insert size is significantly larger/smaller than expected
   * 2. Pair orientation is unexpected (not FR for standard Illumina)
   * 3. Reads map to different chromosomes
   */
  private def detectDiscordantPair(
    record: SAMRecord,
    expectedInsertSize: Double,
    insertSizeMin: Double,
    insertSizeMax: Double
  ): Option[DiscordantPair] = {
    // Skip if mate is unmapped or mapping quality is too low
    if (record.getMateUnmappedFlag) return None
    if (record.getMappingQuality < config.minMapQ) return None

    val readStrand = if (record.getReadNegativeStrandFlag) '-' else '+'
    val mateStrand = if (record.getMateNegativeStrandFlag) '-' else '+'

    // Check for inter-chromosomal (different chromosomes)
    if (record.getReferenceIndex != record.getMateReferenceIndex) {
      return Some(DiscordantPair(
        readName = record.getReadName,
        chrom1 = record.getContig,
        pos1 = record.getAlignmentStart,
        strand1 = readStrand,
        chrom2 = record.getMateReferenceName,
        pos2 = record.getMateAlignmentStart,
        strand2 = mateStrand,
        insertSize = 0,
        mapQ = record.getMappingQuality,
        reason = DiscordantReason.InterChromosomal
      ))
    }

    val insertSize = math.abs(record.getInferredInsertSize)

    // Check for abnormal insert size
    if (insertSize > insertSizeMax || (insertSize > 0 && insertSize < insertSizeMin)) {
      return Some(DiscordantPair(
        readName = record.getReadName,
        chrom1 = record.getContig,
        pos1 = record.getAlignmentStart,
        strand1 = readStrand,
        chrom2 = record.getMateReferenceName,
        pos2 = record.getMateAlignmentStart,
        strand2 = mateStrand,
        insertSize = insertSize,
        mapQ = record.getMappingQuality,
        reason = DiscordantReason.InsertSizeOutlier
      ))
    }

    // Check for unexpected orientation
    // Standard Illumina is FR (forward-reverse)
    if (!isExpectedOrientation(record)) {
      return Some(DiscordantPair(
        readName = record.getReadName,
        chrom1 = record.getContig,
        pos1 = record.getAlignmentStart,
        strand1 = readStrand,
        chrom2 = record.getMateReferenceName,
        pos2 = record.getMateAlignmentStart,
        strand2 = mateStrand,
        insertSize = insertSize,
        mapQ = record.getMappingQuality,
        reason = DiscordantReason.WrongOrientation
      ))
    }

    None
  }

  /**
   * Check if read pair has expected FR orientation.
   * For standard Illumina paired-end:
   * - First read in pair should be on + strand, mate on - strand (when read is upstream)
   * - Or vice versa (when read is downstream)
   */
  private def isExpectedOrientation(record: SAMRecord): Boolean = {
    if (record.getReadNegativeStrandFlag == record.getMateNegativeStrandFlag) {
      // Same strand orientation (tandem) - not expected for FR
      false
    } else if (record.getAlignmentStart < record.getMateAlignmentStart) {
      // Read is upstream of mate - expect read on + strand
      !record.getReadNegativeStrandFlag && record.getMateNegativeStrandFlag
    } else {
      // Read is downstream of mate - expect read on - strand
      record.getReadNegativeStrandFlag && !record.getMateNegativeStrandFlag
    }
  }

  /**
   * Extract split read evidence from the SA (supplementary alignment) tag.
   *
   * The SA tag format is: rname,pos,strand,CIGAR,mapQ,NM;...
   */
  private def extractSplitRead(record: SAMRecord): Option[SplitRead] = {
    val saTag = record.getStringAttribute("SA")
    if (saTag == null || saTag.isEmpty) return None

    // Parse first supplementary alignment from SA tag
    val parts = saTag.split(";")(0).split(",")
    if (parts.length < 5) return None

    try {
      val suppChrom = parts(0)
      val suppPos = parts(1).toLong
      val suppStrand = parts(2).head
      val suppMapQ = parts(4).toInt

      // Calculate clip length from the CIGAR
      val cigar = record.getCigar
      val clipLength = if (cigar != null) {
        cigar.getCigarElements.asScala
          .filter(e => e.getOperator.name == "S" || e.getOperator.name == "H")
          .map(_.getLength)
          .sum
      } else {
        0
      }

      // Only include if both alignments have sufficient quality
      if (suppMapQ >= config.minMapQ && clipLength >= 10) {
        Some(SplitRead(
          readName = record.getReadName,
          primaryChrom = record.getContig,
          primaryPos = record.getAlignmentStart,
          primaryStrand = if (record.getReadNegativeStrandFlag) '-' else '+',
          suppChrom = suppChrom,
          suppPos = suppPos,
          suppStrand = suppStrand,
          clipLength = clipLength,
          mapQ = math.min(record.getMappingQuality, suppMapQ)
        ))
      } else {
        None
      }
    } catch {
      case _: Exception => None
    }
  }
}
