package com.decodingus.analysis

import com.decodingus.model.LibraryStats
import htsjdk.samtools.{SamReaderFactory, SAMRecord, SAMProgramRecord, ValidationStringency}
import java.io.File
import scala.jdk.CollectionConverters._
import scala.collection.mutable

class LibraryStatsProcessor {

  def process(bamPath: String, onProgress: (String, Long, Long) => Unit): LibraryStats = {
    val samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamPath))
    val header = samReader.getFileHeader

    val aligner = detectAligner(header.getProgramRecords.asScala.toList)
    val referenceBuild = header.getSequenceDictionary.getSequences.asScala.headOption.map(_.getSequenceName).getOrElse("Unknown")
    val genomeSize = header.getSequenceDictionary.getReferenceLength

    var readCount = 0
    var totalReadLength = 0L
    var pairedReads = 0
    var totalInsertSize = 0L
    var pairedCount = 0
    val lengthDistribution = mutable.Map[Int, Int]()
    val insertSizeDistribution = mutable.Map[Long, Int]()
    val platformCounts = mutable.Map[String, Int]()

    val recordIterator = samReader.iterator().asScala
    var processedRecords = 0L

    // A full BAM scan can be long, so we need a way to estimate progress.
    // We can't know the total number of records without iterating, so we'll use a large number as a proxy for progress updates.
    val totalRecordsProxy = 100000000L // A large number for progress reporting

    for (record <- recordIterator) {
      processedRecords += 1
      if (!record.isSecondaryOrSupplementary) {
        readCount += 1
        val seqLen = record.getReadLength
        totalReadLength += seqLen
        lengthDistribution(seqLen) = lengthDistribution.getOrElse(seqLen, 0) + 1

        val qname = record.getReadName
        val platform = detectPlatformFromQname(qname)
        platformCounts(platform) = platformCounts.getOrElse(platform, 0) + 1

        if (record.getReadPairedFlag) {
          pairedReads += 1
          if (record.getProperPairFlag && record.getFirstOfPairFlag) {
            val insertSize = record.getInferredInsertSize.abs
            if (insertSize > 0) {
              insertSizeDistribution(insertSize) = insertSizeDistribution.getOrElse(insertSize, 0) + 1
              totalInsertSize += insertSize
              pairedCount += 1
            }
          }
        }
      }

      if (processedRecords % 100000 == 0) {
        onProgress(s"Processed $processedRecords reads...", processedRecords, totalRecordsProxy)
      }
    }

    samReader.close()

    val averageDepth = if (genomeSize > 0) totalReadLength.toDouble / genomeSize else 0.0

    LibraryStats(
      readCount = readCount,
      totalReadLength = totalReadLength,
      pairedReads = pairedReads,
      totalInsertSize = totalInsertSize,
      pairedCount = pairedCount,
      lengthDistribution = lengthDistribution,
      insertSizeDistribution = insertSizeDistribution,
      aligner = aligner,
      referenceBuild = referenceBuild,
      platformCounts = platformCounts,
      genomeSize = genomeSize,
      averageDepth = averageDepth
    )
  }

  private def detectAligner(programRecords: List[SAMProgramRecord]): String = {
    programRecords.headOption.map(_.getProgramName).getOrElse("Unknown")
  }

  private def detectPlatformFromQname(qname: String): String = {
    if (qname.matches(".*:[0-9]+:[0-9]+:[0-9]+:[0-9]+#.*") || qname.matches(".*:[0-9]+:[A-Z0-9]+:[0-9]:[0-9]+:[0-9]+:[0-9]+.*")) {
      "Illumina"
    } else if (qname.matches("^[a-f0-9]{8}-([a-f0-9]{4}-){3}[a-f0-9]{12}.*")) {
      "Nanopore"
    } else if (qname.matches("^m[0-9]{5,}.*")) {
      "PacBio"
    } else {
      "Unknown"
    }
  }
}