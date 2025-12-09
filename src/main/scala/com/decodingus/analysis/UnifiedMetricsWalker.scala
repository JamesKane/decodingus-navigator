package com.decodingus.analysis

import htsjdk.samtools.{SAMRecord, SamReaderFactory, ValidationStringency}

import java.io.File
import scala.collection.mutable
import scala.jdk.CollectionConverters.*

/**
 * Read-level metrics collected during BAM traversal.
 * Replaces CollectAlignmentSummaryMetrics + CollectInsertSizeMetrics (no R dependency).
 */
case class ReadMetrics(
  // Read counts
  totalReads: Long,
  pfReads: Long,                    // Pass-filter reads (non-vendor-failed)
  pfReadsAligned: Long,             // PF reads that aligned
  readsAlignedInPairs: Long,        // Reads aligned as part of a pair
  properPairs: Long,                // Reads aligned as proper pairs

  // Percentages
  pctPfReadsAligned: Double,
  pctReadsAlignedInPairs: Double,
  pctProperPairs: Double,

  // Read length statistics
  medianReadLength: Double,
  meanReadLength: Double,
  stdReadLength: Double,
  minReadLength: Int,
  maxReadLength: Int,
  readLengthHistogram: Map[Int, Long],

  // Insert size (from proper pairs, first-of-pair only)
  medianInsertSize: Double,
  meanInsertSize: Double,
  stdInsertSize: Double,
  minInsertSize: Int,
  maxInsertSize: Int,
  insertSizeHistogram: Map[Int, Long],

  // Pair orientation (FR, RF, TANDEM)
  pairOrientation: String,

  // Quality metrics
  pctChimeras: Double,
  meanMappingQuality: Double
)

/**
 * Single-pass BAM walker that collects read-level metrics.
 *
 * Replaces GATK CollectMultipleMetrics (CollectAlignmentSummaryMetrics + CollectInsertSizeMetrics)
 * without requiring R for histogram generation.
 *
 * Uses HTSJDK SamReader directly for efficient traversal.
 */
class UnifiedMetricsWalker {

  // Insert size histogram bounds
  private val MAX_INSERT_SIZE = 10000

  /**
   * Process a BAM/CRAM file and collect read-level metrics.
   *
   * @param bamPath Path to BAM/CRAM file
   * @param referencePath Path to reference genome (required for CRAM)
   * @param onProgress Progress callback (message, current, total)
   * @return ReadMetrics containing all collected statistics
   */
  def collectReadMetrics(
    bamPath: String,
    referencePath: String,
    onProgress: (String, Long, Long) => Unit
  ): Either[String, ReadMetrics] = {
    try {
      val samReaderFactory = SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT)

      val samReader = if (bamPath.toLowerCase.endsWith(".cram")) {
        samReaderFactory.referenceSequence(new File(referencePath)).open(new File(bamPath))
      } else {
        samReaderFactory.open(new File(bamPath))
      }

      onProgress("Counting reads in BAM/CRAM...", 0, 1)

      // Accumulators
      var totalReads = 0L
      var pfReads = 0L
      var pfReadsAligned = 0L
      var readsAlignedInPairs = 0L
      var properPairs = 0L
      var chimericReads = 0L

      // Read length tracking
      val readLengthHist = mutable.Map[Int, Long]().withDefaultValue(0L)
      var readLengthSum = 0L
      var readLengthSumSq = 0L
      var minReadLen = Int.MaxValue
      var maxReadLen = 0

      var totalMappingQuality = 0L
      var mappedReadsForMQ = 0L

      // Insert size tracking (only from proper pairs, first-of-pair)
      val insertSizeHist = mutable.Map[Int, Long]().withDefaultValue(0L)
      var insertSizeSum = 0L
      var insertSizeSumSq = 0L
      var insertSizeCount = 0L
      var minInsert = Int.MaxValue
      var maxInsert = 0

      // Pair orientation tracking
      var frCount = 0L
      var rfCount = 0L
      var tandemCount = 0L

      val iterator = samReader.iterator()
      var recordCount = 0L
      val progressInterval = 1000000L

      while (iterator.hasNext) {
        val record = iterator.next()
        recordCount += 1

        if (recordCount % progressInterval == 0) {
          onProgress(s"Processed ${recordCount / 1000000}M reads...", recordCount, recordCount + progressInterval)
        }

        // Skip secondary and supplementary alignments for primary metrics
        if (!record.isSecondaryOrSupplementary) {
          totalReads += 1

          // Pass-filter check (not vendor failed)
          if (!record.getReadFailsVendorQualityCheckFlag) {
            pfReads += 1

            val readLength = record.getReadLength
            readLengthHist(readLength) += 1
            readLengthSum += readLength
            readLengthSumSq += readLength.toLong * readLength
            if (readLength < minReadLen) minReadLen = readLength
            if (readLength > maxReadLen) maxReadLen = readLength

            // Alignment metrics
            if (!record.getReadUnmappedFlag) {
              pfReadsAligned += 1

              val mapQ = record.getMappingQuality
              if (mapQ != 255) { // 255 means unavailable
                totalMappingQuality += mapQ
                mappedReadsForMQ += 1
              }

              // Paired read metrics
              if (record.getReadPairedFlag) {
                if (!record.getMateUnmappedFlag) {
                  readsAlignedInPairs += 1
                }

                if (record.getProperPairFlag) {
                  properPairs += 1

                  // Collect insert size from first-of-pair only to avoid double counting
                  if (record.getFirstOfPairFlag) {
                    val insertSize = math.abs(record.getInferredInsertSize)
                    if (insertSize > 0 && insertSize < MAX_INSERT_SIZE) {
                      insertSizeHist(insertSize) += 1
                      insertSizeSum += insertSize
                      insertSizeSumSq += insertSize.toLong * insertSize
                      insertSizeCount += 1
                      if (insertSize < minInsert) minInsert = insertSize
                      if (insertSize > maxInsert) maxInsert = insertSize

                      // Determine pair orientation
                      val orientation = detectPairOrientation(record)
                      orientation match {
                        case "FR" => frCount += 1
                        case "RF" => rfCount += 1
                        case "TANDEM" => tandemCount += 1
                        case _ => // ignore
                      }
                    }
                  }
                }

                // Chimeric read detection (mapped to different chromosome than mate)
                if (!record.getReadUnmappedFlag && !record.getMateUnmappedFlag) {
                  if (record.getReferenceIndex != record.getMateReferenceIndex) {
                    chimericReads += 1
                  }
                }
              }
            }
          }
        }
      }

      samReader.close()
      onProgress("Calculating statistics...", recordCount, recordCount)

      // Calculate derived metrics
      val pctPfReadsAligned = if (pfReads > 0) pfReadsAligned.toDouble / pfReads else 0.0
      val pctReadsAlignedInPairs = if (pfReadsAligned > 0) readsAlignedInPairs.toDouble / pfReadsAligned else 0.0
      val pctProperPairs = if (pfReadsAligned > 0) properPairs.toDouble / pfReadsAligned else 0.0
      val pctChimeras = if (readsAlignedInPairs > 0) chimericReads.toDouble / readsAlignedInPairs else 0.0

      val meanMappingQuality = if (mappedReadsForMQ > 0) totalMappingQuality.toDouble / mappedReadsForMQ else 0.0

      // Read length statistics
      val (medianReadLen, meanReadLen, stdReadLen) = if (pfReads > 0) {
        val mean = readLengthSum.toDouble / pfReads
        val variance = (readLengthSumSq.toDouble / pfReads) - (mean * mean)
        val std = math.sqrt(math.max(0, variance))
        val median = calculateMedianFromHistogram(readLengthHist.toMap, pfReads)
        (median, mean, std)
      } else {
        (0.0, 0.0, 0.0)
      }

      // Insert size statistics
      val (medianInsert, meanInsert, stdInsert) = if (insertSizeCount > 0) {
        val mean = insertSizeSum.toDouble / insertSizeCount
        val variance = (insertSizeSumSq.toDouble / insertSizeCount) - (mean * mean)
        val std = math.sqrt(math.max(0, variance))
        val median = calculateMedianFromHistogram(insertSizeHist.toMap, insertSizeCount)
        (median, mean, std)
      } else {
        (0.0, 0.0, 0.0)
      }

      // Determine dominant pair orientation
      val pairOrientation = if (frCount >= rfCount && frCount >= tandemCount) "FR"
        else if (rfCount >= frCount && rfCount >= tandemCount) "RF"
        else "TANDEM"

      Right(ReadMetrics(
        totalReads = totalReads,
        pfReads = pfReads,
        pfReadsAligned = pfReadsAligned,
        readsAlignedInPairs = readsAlignedInPairs,
        properPairs = properPairs,
        pctPfReadsAligned = pctPfReadsAligned,
        pctReadsAlignedInPairs = pctReadsAlignedInPairs,
        pctProperPairs = pctProperPairs,
        medianReadLength = medianReadLen,
        meanReadLength = meanReadLen,
        stdReadLength = stdReadLen,
        minReadLength = if (minReadLen == Int.MaxValue) 0 else minReadLen,
        maxReadLength = maxReadLen,
        readLengthHistogram = readLengthHist.toMap,
        medianInsertSize = medianInsert,
        meanInsertSize = meanInsert,
        stdInsertSize = stdInsert,
        minInsertSize = if (minInsert == Int.MaxValue) 0 else minInsert,
        maxInsertSize = maxInsert,
        insertSizeHistogram = insertSizeHist.toMap,
        pairOrientation = pairOrientation,
        pctChimeras = pctChimeras,
        meanMappingQuality = meanMappingQuality
      ))

    } catch {
      case e: Exception =>
        Left(s"Failed to collect read metrics: ${e.getMessage}")
    }
  }

  /**
   * Detect pair orientation from a properly paired first-of-pair read.
   * FR = forward-reverse (standard Illumina)
   * RF = reverse-forward (mate-pair libraries)
   * TANDEM = same orientation
   */
  private def detectPairOrientation(record: SAMRecord): String = {
    val readNegStrand = record.getReadNegativeStrandFlag
    val mateNegStrand = record.getMateNegativeStrandFlag

    if (readNegStrand == mateNegStrand) {
      "TANDEM"
    } else if (record.getAlignmentStart < record.getMateAlignmentStart) {
      // Read is upstream of mate
      if (!readNegStrand && mateNegStrand) "FR" else "RF"
    } else {
      // Read is downstream of mate
      if (readNegStrand && !mateNegStrand) "FR" else "RF"
    }
  }

  /**
   * Calculate median from a histogram.
   */
  private def calculateMedianFromHistogram(histogram: Map[Int, Long], totalCount: Long): Double = {
    if (totalCount == 0) {
      0.0
    } else {
      val sortedKeys = histogram.keys.toSeq.sorted
      val medianPos = totalCount / 2
      var cumulative = 0L
      var result: Option[Double] = None

      val iter = sortedKeys.iterator
      while (iter.hasNext && result.isEmpty) {
        val key = iter.next()
        cumulative += histogram(key)
        if (cumulative >= medianPos) {
          result = Some(key.toDouble)
        }
      }

      result.getOrElse(sortedKeys.lastOption.map(_.toDouble).getOrElse(0.0))
    }
  }
}
