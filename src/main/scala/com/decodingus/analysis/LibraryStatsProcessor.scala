package com.decodingus.analysis

import com.decodingus.model.LibraryStats
import htsjdk.samtools.{SAMProgramRecord, SamReaderFactory, ValidationStringency}

import java.io.File
import scala.collection.mutable
import scala.jdk.CollectionConverters.*
import scala.util.boundary
import scala.util.boundary.break

class LibraryStatsProcessor {

  private val MAX_SAMPLES = 10000

  def process(bamPath: String, onProgress: (String, Long, Long) => Unit): LibraryStats = {
    val samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamPath))
    val header = samReader.getFileHeader

    val aligner = detectAligner(header.getProgramRecords.asScala.toList)
    val referenceBuild = header.getSequenceDictionary.getSequences.asScala.headOption.map(_.getSequenceName).getOrElse("Unknown")
    val sampleName = header.getReadGroups.asScala.headOption.map(_.getSample).getOrElse("Unknown")

    var readCount = 0
    var pairedReads = 0
    val lengthDistribution = mutable.Map[Int, Int]()
    val insertSizeDistribution = mutable.Map[Long, Int]()
    val platformCounts = mutable.Map[String, Int]()
    val instruments = mutable.Map[String, Int]()
    val flowCells = mutable.Map[String, Int]()

    val recordIterator = samReader.iterator().asScala
    var processedRecords = 0

    boundary {
      for (record <- recordIterator) {
        if (processedRecords >= MAX_SAMPLES) {
          break(buildLibraryStats(readCount, pairedReads, lengthDistribution, insertSizeDistribution, aligner, referenceBuild, sampleName, flowCells, instruments, platformCounts))
        }

        processedRecords += 1
        if (!record.isSecondaryOrSupplementary) {
          readCount += 1
          val seqLen = record.getReadLength
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
              }
            }
          }
        }

        if (processedRecords % 1000 == 0) {
          onProgress(s"Scanned $processedRecords reads...", processedRecords, MAX_SAMPLES)
        }
      }
    }

    samReader.close()
    buildLibraryStats(readCount, pairedReads, lengthDistribution, insertSizeDistribution, aligner, referenceBuild, sampleName, flowCells, instruments, platformCounts)
  }

  private def buildLibraryStats(readCount: Int, pairedReads: Int, lengthDistribution: mutable.Map[Int, Int], insertSizeDistribution: mutable.Map[Long, Int], aligner: String, referenceBuild: String, sampleName: String, flowCells: mutable.Map[String, Int], instruments: mutable.Map[String, Int], platformCounts: mutable.Map[String, Int]): LibraryStats = {
    val mostFrequentInstrument = instruments.toSeq.sortBy(-_._2).headOption.map(_._1).getOrElse("Unknown")
    val primaryPlatform = platformCounts.toSeq.sortBy(-_._2).headOption.map(_._1).getOrElse("Unknown")
    val inferredPlatform = inferPlatform(primaryPlatform, mostFrequentInstrument)

    LibraryStats(
      readCount = readCount,
      pairedReads = pairedReads,
      lengthDistribution = lengthDistribution,
      insertSizeDistribution = insertSizeDistribution,
      aligner = aligner,
      referenceBuild = referenceBuild,
      sampleName = sampleName,
      flowCells = flowCells,
      instruments = instruments,
      mostFrequentInstrument = mostFrequentInstrument,
      inferredPlatform = inferredPlatform,
      platformCounts = platformCounts
    )
  }

  private def inferPlatform(platform: String, instrumentId: String): String = {
    platform match {
      case "Illumina" =>
        instrumentId.headOption match {
          case Some('A' | 'a') => "NovaSeq"
          case Some('D' | 'd') => "HiSeq 2500"
          case Some('J' | 'j') => "HiSeq 3000"
          case Some('K' | 'k') => "HiSeq 4000"
          case Some('E' | 'e') => "HiSeq X"
          case Some('N' | 'n') => "NextSeq"
          case Some('M' | 'm') => "MiSeq"
          case Some('V' | 'v') => "NovaSeq X"
          case Some('F' | 'f') => "iSeq"
          case _ => "Unknown Illumina"
        }
      case "PacBio" =>
        if (instrumentId.startsWith("m84")) "PacBio Revio"
        else if (instrumentId.startsWith("m64")) "PacBio Sequel II/IIe"
        else if (instrumentId.startsWith("m54")) "PacBio Sequel"
        else "PacBio"
      case "MGI" =>
        if (instrumentId.startsWith("V300")) "MGI DNBSEQ/MGISEQ-2000"
        else if (instrumentId.startsWith("E100")) "MGI MGISEQ-200"
        else if (instrumentId.startsWith("CL100")) "MGI MGISEQ-T7"
        else if (instrumentId.startsWith("G400")) "MGI DNBSEQ-G400"
        else if (instrumentId.startsWith("G99")) "MGI MGISEQ-T1"
        else "MGI DNBseq"
      case "Nanopore" => "Oxford Nanopore"
      case _ => "Unknown"
    }
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
