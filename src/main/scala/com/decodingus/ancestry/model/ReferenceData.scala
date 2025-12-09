package com.decodingus.ancestry.model

import java.io.{DataInputStream, DataOutputStream, File, FileInputStream, FileOutputStream}
import java.nio.{ByteBuffer, ByteOrder}
import java.nio.channels.FileChannel
import java.nio.file.{Path, StandardOpenOption}
import scala.util.Using

/**
 * Binary format for allele frequencies across populations.
 *
 * File format:
 * - Header: magic (4 bytes) + version (2 bytes) + numPops (2 bytes) + numSnps (4 bytes)
 * - Population codes: numPops * 32 bytes (null-padded strings)
 * - SNP IDs: numSnps * 32 bytes (null-padded strings)
 * - Frequencies: numPops * numSnps * 4 bytes (floats, SNP-major order)
 */
case class AlleleFrequencyMatrix(
  populations: Array[String],
  snpIds: Array[String],
  frequencies: Array[Float]  // Flattened: frequencies[snpIdx * numPops + popIdx]
) {
  require(frequencies.length == populations.length * snpIds.length,
    s"Frequency array size mismatch: expected ${populations.length * snpIds.length}, got ${frequencies.length}")

  val numPopulations: Int = populations.length
  val numSnps: Int = snpIds.length

  /**
   * Get allele frequency for a population at a SNP.
   */
  def getFrequency(popIndex: Int, snpIndex: Int): Float = {
    frequencies(snpIndex * numPopulations + popIndex)
  }

  /**
   * Get allele frequency by population code and SNP ID.
   */
  def getFrequency(popCode: String, snpId: String): Option[Float] = {
    for {
      popIdx <- populations.indexOf(popCode) match { case -1 => None; case i => Some(i) }
      snpIdx <- snpIds.indexOf(snpId) match { case -1 => None; case i => Some(i) }
    } yield getFrequency(popIdx, snpIdx)
  }

  /**
   * Get all frequencies for a SNP across populations.
   */
  def getFrequenciesForSnp(snpIndex: Int): Array[Float] = {
    val result = new Array[Float](numPopulations)
    val offset = snpIndex * numPopulations
    System.arraycopy(frequencies, offset, result, 0, numPopulations)
    result
  }
}

object AlleleFrequencyMatrix {
  private val MAGIC = 0x41464D58  // "AFMX"
  private val VERSION: Short = 1
  private val STRING_SIZE = 32

  /**
   * Load allele frequency matrix from binary file.
   */
  def load(path: Path): Either[String, AlleleFrequencyMatrix] = {
    Using(new DataInputStream(new FileInputStream(path.toFile))) { dis =>
      // Read header
      val magic = dis.readInt()
      if (magic != MAGIC) {
        throw new IllegalArgumentException(s"Invalid magic number: expected $MAGIC, got $magic")
      }
      val version = dis.readShort()
      if (version != VERSION) {
        throw new IllegalArgumentException(s"Unsupported version: $version")
      }
      val numPops = dis.readShort().toInt
      val numSnps = dis.readInt()

      // Read population codes
      val populations = (0 until numPops).map { _ =>
        val bytes = new Array[Byte](STRING_SIZE)
        dis.readFully(bytes)
        new String(bytes).trim.stripSuffix("\u0000")
      }.toArray

      // Read SNP IDs
      val snpIds = (0 until numSnps).map { _ =>
        val bytes = new Array[Byte](STRING_SIZE)
        dis.readFully(bytes)
        new String(bytes).trim.stripSuffix("\u0000")
      }.toArray

      // Read frequencies
      val freqCount = numPops * numSnps
      val frequencies = new Array[Float](freqCount)
      for (i <- 0 until freqCount) {
        frequencies(i) = dis.readFloat()
      }

      AlleleFrequencyMatrix(populations, snpIds, frequencies)
    }.toEither.left.map(_.getMessage)
  }

  /**
   * Save allele frequency matrix to binary file.
   */
  def save(matrix: AlleleFrequencyMatrix, path: Path): Either[String, Unit] = {
    Using(new DataOutputStream(new FileOutputStream(path.toFile))) { dos =>
      // Write header
      dos.writeInt(MAGIC)
      dos.writeShort(VERSION)
      dos.writeShort(matrix.numPopulations.toShort)
      dos.writeInt(matrix.numSnps)

      // Write population codes
      matrix.populations.foreach { code =>
        val bytes = code.take(STRING_SIZE).padTo(STRING_SIZE, '\u0000').getBytes
        dos.write(bytes)
      }

      // Write SNP IDs
      matrix.snpIds.foreach { id =>
        val bytes = id.take(STRING_SIZE).padTo(STRING_SIZE, '\u0000').getBytes
        dos.write(bytes)
      }

      // Write frequencies
      matrix.frequencies.foreach(dos.writeFloat)
    }.toEither.left.map(_.getMessage)
  }
}

/**
 * PCA loadings for projecting samples onto reference population space.
 *
 * File format:
 * - Header: magic (4 bytes) + version (2 bytes) + numSnps (4 bytes) + numComponents (2 bytes) + numPops (2 bytes)
 * - SNP IDs: numSnps * 32 bytes
 * - SNP means: numSnps * 4 bytes (floats) - for centering
 * - Loadings: numSnps * numComponents * 4 bytes (floats)
 * - Population codes: numPops * 32 bytes
 * - Centroids: numPops * numComponents * 4 bytes
 * - Variances: numPops * numComponents * 4 bytes (diagonal covariance)
 */
case class PCALoadings(
  snpIds: Array[String],
  snpMeans: Array[Float],              // Mean genotype for each SNP (for centering)
  loadings: Array[Float],              // Flattened: loadings[snpIdx * numComponents + pcIdx]
  numComponents: Int,
  populations: Array[String],
  centroids: Array[Float],             // Flattened: centroids[popIdx * numComponents + pcIdx]
  variances: Array[Float]              // Flattened: variances[popIdx * numComponents + pcIdx]
) {
  val numSnps: Int = snpIds.length
  val numPopulations: Int = populations.length

  /**
   * Get PCA loading for a SNP and component.
   */
  def getLoading(snpIndex: Int, componentIndex: Int): Float = {
    loadings(snpIndex * numComponents + componentIndex)
  }

  /**
   * Get centroid for a population.
   */
  def getCentroid(popIndex: Int): Array[Float] = {
    val result = new Array[Float](numComponents)
    val offset = popIndex * numComponents
    System.arraycopy(centroids, offset, result, 0, numComponents)
    result
  }

  /**
   * Get variance (diagonal covariance) for a population.
   */
  def getVariance(popIndex: Int): Array[Float] = {
    val result = new Array[Float](numComponents)
    val offset = popIndex * numComponents
    System.arraycopy(variances, offset, result, 0, numComponents)
    result
  }
}

object PCALoadings {
  private val MAGIC = 0x50434C44  // "PCLD"
  private val VERSION: Short = 1
  private val STRING_SIZE = 32

  /**
   * Load PCA loadings from binary file.
   */
  def load(path: Path): Either[String, PCALoadings] = {
    Using(new DataInputStream(new FileInputStream(path.toFile))) { dis =>
      // Read header
      val magic = dis.readInt()
      if (magic != MAGIC) {
        throw new IllegalArgumentException(s"Invalid magic number: expected $MAGIC, got $magic")
      }
      val version = dis.readShort()
      if (version != VERSION) {
        throw new IllegalArgumentException(s"Unsupported version: $version")
      }
      val numSnps = dis.readInt()
      val numComponents = dis.readShort().toInt
      val numPops = dis.readShort().toInt

      // Read SNP IDs
      val snpIds = (0 until numSnps).map { _ =>
        val bytes = new Array[Byte](STRING_SIZE)
        dis.readFully(bytes)
        new String(bytes).trim.stripSuffix("\u0000")
      }.toArray

      // Read SNP means
      val snpMeans = (0 until numSnps).map(_ => dis.readFloat()).toArray

      // Read loadings
      val loadingsCount = numSnps * numComponents
      val loadings = (0 until loadingsCount).map(_ => dis.readFloat()).toArray

      // Read population codes
      val populations = (0 until numPops).map { _ =>
        val bytes = new Array[Byte](STRING_SIZE)
        dis.readFully(bytes)
        new String(bytes).trim.stripSuffix("\u0000")
      }.toArray

      // Read centroids
      val centroidsCount = numPops * numComponents
      val centroids = (0 until centroidsCount).map(_ => dis.readFloat()).toArray

      // Read variances
      val variances = (0 until centroidsCount).map(_ => dis.readFloat()).toArray

      PCALoadings(snpIds, snpMeans, loadings, numComponents, populations, centroids, variances)
    }.toEither.left.map(_.getMessage)
  }

  /**
   * Save PCA loadings to binary file.
   */
  def save(pca: PCALoadings, path: Path): Either[String, Unit] = {
    Using(new DataOutputStream(new FileOutputStream(path.toFile))) { dos =>
      // Write header
      dos.writeInt(MAGIC)
      dos.writeShort(VERSION)
      dos.writeInt(pca.numSnps)
      dos.writeShort(pca.numComponents.toShort)
      dos.writeShort(pca.numPopulations.toShort)

      // Write SNP IDs
      pca.snpIds.foreach { id =>
        val bytes = id.take(STRING_SIZE).padTo(STRING_SIZE, '\u0000').getBytes
        dos.write(bytes)
      }

      // Write SNP means
      pca.snpMeans.foreach(dos.writeFloat)

      // Write loadings
      pca.loadings.foreach(dos.writeFloat)

      // Write population codes
      pca.populations.foreach { code =>
        val bytes = code.take(STRING_SIZE).padTo(STRING_SIZE, '\u0000').getBytes
        dos.write(bytes)
      }

      // Write centroids
      pca.centroids.foreach(dos.writeFloat)

      // Write variances
      pca.variances.foreach(dos.writeFloat)
    }.toEither.left.map(_.getMessage)
  }
}

/**
 * Panel type for ancestry analysis.
 */
enum AncestryPanelType {
  case Aims        // ~5,000 ancestry-informative markers
  case GenomeWide  // ~500,000 common SNPs
}

object AncestryPanelType {
  def fromString(s: String): Option[AncestryPanelType] = s.toLowerCase match {
    case "aims" => Some(Aims)
    case "genome-wide" | "genomewide" | "full" => Some(GenomeWide)
    case _ => None
  }
}
