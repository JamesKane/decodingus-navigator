package com.decodingus.ibd.engine

import htsjdk.variant.vcf.VCFFileReader

import java.io.File
import scala.collection.mutable.ArrayBuffer
import scala.jdk.CollectionConverters.*

/**
 * Extracts autosomal biallelic SNP genotypes from VCF files
 * into position-sorted arrays suitable for pairwise IBD detection.
 */
object IbdVariantExtractor:

  private val autosomePattern = "^(chr)?([1-9]|1[0-9]|2[0-2])$".r

  /**
   * Extract autosomal genotypes from a VCF file.
   *
   * @param vcfFile    Path to VCF or VCF.gz file
   * @param sampleName Sample name in the VCF (uses first sample if None)
   * @param onProgress Optional progress callback (chromosome, fraction)
   * @return Either error or map of chromosome → genotypes
   */
  def extractFromVcf(
                      vcfFile: File,
                      sampleName: Option[String] = None,
                      onProgress: Option[(String, Double) => Unit] = None
                    ): Either[String, Map[String, ChromosomeGenotypes]] =
    try
      val reader = new VCFFileReader(vcfFile, false)
      val header = reader.getFileHeader
      val sample = sampleName.getOrElse(header.getGenotypeSamples.get(0))

      // Accumulate per chromosome
      val chrData = scala.collection.mutable.Map.empty[String,
        (ArrayBuffer[Int], ArrayBuffer[Byte])]

      var totalVariants = 0
      var skippedNonBiallelic = 0
      var skippedNonAutosomal = 0

      for vc <- reader.iterator().asScala do
        val contig = vc.getContig
        val chrNum = normalizeAutosome(contig)

        chrNum match
          case Some(chr) =>
            // Only biallelic SNPs
            if vc.isSNP && vc.isBiallelic then
              val genotype = vc.getGenotype(sample)
              val gt: Byte =
                if genotype.isNoCall then -1
                else if genotype.isHomRef then 0
                else if genotype.isHet then 1
                else if genotype.isHomVar then 2
                else -1

              val (positions, genotypes) = chrData.getOrElseUpdate(chr,
                (ArrayBuffer.empty[Int], ArrayBuffer.empty[Byte]))
              positions += vc.getStart
              genotypes += gt
              totalVariants += 1
            else
              skippedNonBiallelic += 1
          case None =>
            skippedNonAutosomal += 1

      reader.close()

      onProgress.foreach(_(s"Extracted $totalVariants autosomal SNPs", 1.0))

      val result = chrData.map { case (chr, (positions, genotypes)) =>
        chr -> ChromosomeGenotypes(chr, positions.toArray, genotypes.toArray)
      }.toMap

      Right(result)
    catch
      case e: Exception => Left(s"Failed to extract genotypes: ${e.getMessage}")

  /**
   * Extract genotypes from a pre-computed numeric array (for encrypted exchange).
   * The format is a compact binary: position (4 bytes) + genotype (1 byte) per entry.
   *
   * @param data       Serialized genotype data
   * @param chromosome Chromosome this data belongs to
   * @return ChromosomeGenotypes
   */
  def fromCompactBytes(data: Array[Byte], chromosome: String): ChromosomeGenotypes =
    val entrySize = 5 // 4 bytes position + 1 byte genotype
    val n = data.length / entrySize
    val positions = new Array[Int](n)
    val genotypes = new Array[Byte](n)

    var i = 0
    var offset = 0
    while i < n do
      positions(i) = ((data(offset) & 0xFF) << 24) |
        ((data(offset + 1) & 0xFF) << 16) |
        ((data(offset + 2) & 0xFF) << 8) |
        (data(offset + 3) & 0xFF)
      genotypes(i) = data(offset + 4)
      i += 1
      offset += entrySize

    ChromosomeGenotypes(chromosome, positions, genotypes)

  /**
   * Serialize genotypes to compact binary format for encrypted exchange.
   */
  def toCompactBytes(genotypes: ChromosomeGenotypes): Array[Byte] =
    val entrySize = 5
    val data = new Array[Byte](genotypes.size * entrySize)

    var i = 0
    var offset = 0
    while i < genotypes.size do
      val pos = genotypes.positions(i)
      data(offset) = ((pos >> 24) & 0xFF).toByte
      data(offset + 1) = ((pos >> 16) & 0xFF).toByte
      data(offset + 2) = ((pos >> 8) & 0xFF).toByte
      data(offset + 3) = (pos & 0xFF).toByte
      data(offset + 4) = genotypes.genotypes(i)
      i += 1
      offset += entrySize

    data

  /**
   * Normalize autosomal chromosome name.
   * Returns Some("1") through Some("22") for autosomes, None otherwise.
   */
  private def normalizeAutosome(contig: String): Option[String] =
    contig match
      case autosomePattern(_, num) => Some(num)
      case _ => None
