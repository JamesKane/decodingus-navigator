package com.decodingus.refgenome

import htsjdk.samtools.reference.{ReferenceSequenceFile, ReferenceSequenceFileFactory}

import java.io.File

/**
 * Queries bases from a reference genome for a single contig.
 * Loads the entire contig into memory on first access for fast lookups.
 *
 * @param referencePath Path to the reference FASTA file
 * @param contig The contig to query (e.g., "chrY", "chrM")
 */
class ReferenceQuerier(referencePath: String, contig: String) extends AutoCloseable {
  private val referenceFile: ReferenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(referencePath))
  private val bases: Array[Byte] = referenceFile.getSequence(contig).getBases

  /** Length of the contig in bases */
  val length: Int = bases.length

  /**
   * Check if a position is within the contig bounds.
   * @param position 1-based position
   */
  def isValidPosition(position: Long): Boolean = {
    position >= 1 && position <= length
  }

  /**
   * Get the base at a 1-based position.
   * @return Some(base) if position is valid, None if out of bounds
   */
  def getBase(position: Long): Option[Char] = {
    if (isValidPosition(position)) {
      Some(bases(position.toInt - 1).toChar)
    } else {
      None
    }
  }

  /**
   * Get the base at a 1-based position, throwing if out of bounds.
   * @throws IndexOutOfBoundsException if position is outside contig
   */
  def getBaseUnsafe(position: Long): Char = {
    bases(position.toInt - 1).toChar
  }

  override def close(): Unit = {
    referenceFile.close()
  }
}

/**
 * Queries bases from a reference genome across multiple contigs.
 * Caches loaded contigs in memory for fast repeated lookups.
 *
 * @param referencePath Path to the reference FASTA file
 */
class MultiContigReferenceQuerier(referencePath: String) extends AutoCloseable {
  private val referenceFile: ReferenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(referencePath))
  private var cachedContig: String = _
  private var cachedBases: Array[Byte] = _

  private def ensureContig(contig: String): Unit = {
    if (cachedContig != contig) {
      cachedContig = contig
      cachedBases = referenceFile.getSequence(contig).getBases
    }
  }

  /** Get the length of a contig */
  def contigLength(contig: String): Int = {
    ensureContig(contig)
    cachedBases.length
  }

  /**
   * Check if a position is within the contig bounds.
   * @param position 1-based position
   */
  def isValidPosition(contig: String, position: Long): Boolean = {
    ensureContig(contig)
    position >= 1 && position <= cachedBases.length
  }

  /**
   * Get the base at a 1-based position.
   * @return Some(base) if position is valid, None if out of bounds
   */
  def getBase(contig: String, position: Long): Option[Char] = {
    ensureContig(contig)
    if (position >= 1 && position <= cachedBases.length) {
      Some(cachedBases(position.toInt - 1).toChar)
    } else {
      None
    }
  }

  /**
   * Get the base at a 1-based position, throwing if out of bounds.
   * @throws IndexOutOfBoundsException if position is outside contig
   */
  def getBaseUnsafe(contig: String, position: Long): Char = {
    ensureContig(contig)
    cachedBases(position.toInt - 1).toChar
  }

  override def close(): Unit = {
    cachedBases = null
    referenceFile.close()
  }
}
