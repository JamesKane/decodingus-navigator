package com.decodingus.refgenome

import htsjdk.samtools.reference.{ReferenceSequenceFile, ReferenceSequenceFileFactory}

import java.io.File

class ReferenceQuerier(referencePath: String) extends AutoCloseable {
  private val referenceFile: ReferenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(referencePath))

  def getBase(contig: String, position: Long): Char = {
    val sequence = referenceFile.getSequence(contig)
    sequence.getBases()(position.toInt - 1).toChar
  }

  override def close(): Unit = {
    referenceFile.close()
  }
}
