package com.decodingus.analysis

import com.decodingus.model.LibraryStats
import htsjdk.samtools.{SAMProgramRecord, SamReaderFactory, ValidationStringency}

import java.io.File
import scala.collection.mutable
import scala.jdk.CollectionConverters.*

class LibraryStatsProcessor {

  def process(bamPath: String, onProgress: (String, Long, Long) => Unit): LibraryStats = {
    val samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamPath))
    val header = samReader.getFileHeader

    val aligner = detectAligner(header.getProgramRecords.asScala.toList)
    val referenceBuild = header.getSequenceDictionary.getSequences.asScala.headOption.map(_.getSequenceName).getOrElse("Unknown")
    val genomeSize = header.getSequenceDictionary.getReferenceLength
    val sampleName = header.getReadGroups.asScala.headOption.map(_.getSample).getOrElse("Unknown")

    var readCount = 0
    var totalReadLength = 0L
    var pairedReads = 0
    var totalInsertSize = 0L
    var pairedCount = 0
    val lengthDistribution = mutable.Map[Int, Int]()
    val insertSizeDistribution = mutable.Map[Long, Int]()
    val platformCounts = mutable.Map[String, Int]()
    val instruments = mutable.Map[String, Int]()
    val flowCells = mutable.Map[String, Int]()

    val recordIterator = samReader.iterator().asScala
    var processedRecords = 0L
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

        parseInstrumentAndFlowcell(qname, platform).foreach { case (instrument, flowcell) =>
          instrument.foreach(i => instruments(i) = instruments.getOrElse(i, 0) + 1)
          flowcell.foreach(f => flowCells(f) = flowCells.getOrElse(f, 0) + 1)
        }

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
    val mostFrequentInstrument = instruments.toSeq.sortBy(-_._2).headOption.map(_._1).getOrElse("Unknown")

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
      sampleName = sampleName,
      flowCells = flowCells,
      instruments = instruments,
      mostFrequentInstrument = mostFrequentInstrument,
      platformCounts = platformCounts,
      genomeSize = genomeSize,
      averageDepth = averageDepth
    )
  }

  private def detectAligner(programRecords: List[SAMProgramRecord]): String = {
    programRecords.headOption.map(_.getProgramName).getOrElse("Unknown")
  }

  private def detectPlatformFromQname(qname: String): String = {
    if (qname.length > 15) {
      val prefix = qname.substring(0, 5).toUpperCase
      if (prefix.startsWith("V300") || prefix.startsWith("E100") || prefix.startsWith("CL100") || prefix.startsWith("G400") || prefix.startsWith("G99")) {
        return "MGI"
      }
      if (qname.count(_ == ':') >= 6) {
        val parts = qname.split(':')
        if (parts(0).startsWith("V") || parts(0).startsWith("E") || parts(0).startsWith("CL") || parts(0).startsWith("G")) {
          if (parts.length >= 3 && parts(2).startsWith("L")) {
            return "MGI"
          }
        }
      }
    }
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

  private def parseInstrumentAndFlowcell(qname: String, platform: String): Option[(Option[String], Option[String])] = {
    platform match {
      case "Illumina" =>
        val parts = qname.split(":")
        if (parts.length >= 3) Some(Some(parts(0)), Some(parts(2))) else None
      case "PacBio" =>
        val parts = qname.split("/")
        if (parts.nonEmpty) Some(Some(parts(0).split("_")(0)), None) else None
      case "MGI" =>
        if (qname.count(_ == ':') >= 3) {
          val parts = qname.split(':')
          if (parts.length >= 3) {
            Some(Some(parts(0)), Some(parts(1)))
          } else {
            None
          }
        } else {
          if (qname.length > 10) {
            val lPos = qname.indexOf('L')
            if (lPos > 0) {
              val instrument = qname.substring(0, lPos)
              val cPos = qname.substring(lPos).indexOf('C')
              if (cPos > 0) {
                val rPos = qname.substring(lPos).indexOf('R')
                val endPos = if (rPos > 0) rPos else qname.substring(lPos).length
                val flowcell = qname.substring(lPos, lPos + endPos)
                Some(Some(instrument), Some(flowcell))
              } else {
                None
              }
            } else {
              None
            }
          } else {
            None
          }
        }
      case _ => None
    }
  }
}