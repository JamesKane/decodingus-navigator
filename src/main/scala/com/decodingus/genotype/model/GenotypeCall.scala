package com.decodingus.genotype.model

import io.circe.Codec

/**
 * A single genotype call from a chip/array file.
 *
 * @param markerId rsID or vendor marker name (e.g., "rs12345")
 * @param chromosome Chromosome name (e.g., "1", "X", "Y", "MT")
 * @param position 1-based position on the chromosome
 * @param allele1 First allele (A, C, G, T, I, D, or - for no call)
 * @param allele2 Second allele (A, C, G, T, I, D, or - for no call)
 */
case class GenotypeCall(
  markerId: String,
  chromosome: String,
  position: Int,
  allele1: Char,
  allele2: Char
) derives Codec.AsObject {

  /**
   * Check if this is a no-call (missing data).
   */
  def isNoCall: Boolean = allele1 == '-' || allele2 == '-' || allele1 == '0' || allele2 == '0'

  /**
   * Check if this is homozygous.
   */
  def isHomozygous: Boolean = !isNoCall && allele1 == allele2

  /**
   * Check if this is heterozygous.
   */
  def isHeterozygous: Boolean = !isNoCall && allele1 != allele2

  /**
   * Get genotype as string (e.g., "AA", "AG", "--").
   */
  def genotype: String = s"$allele1$allele2"

  /**
   * Normalize the genotype to have alleles in alphabetical order.
   * This helps with comparison across different file formats.
   */
  def normalized: GenotypeCall =
    if (allele1 <= allele2) this
    else copy(allele1 = allele2, allele2 = allele1)

  /**
   * Check if this is on the Y chromosome.
   */
  def isYChromosome: Boolean =
    chromosome.equalsIgnoreCase("Y") ||
    chromosome.equalsIgnoreCase("chrY") ||
    chromosome == "24"

  /**
   * Check if this is on mitochondrial DNA.
   */
  def isMitochondrial: Boolean =
    chromosome.equalsIgnoreCase("MT") ||
    chromosome.equalsIgnoreCase("M") ||
    chromosome.equalsIgnoreCase("chrM") ||
    chromosome.equalsIgnoreCase("chrMT") ||
    chromosome == "26"

  /**
   * Check if this is on the X chromosome.
   */
  def isXChromosome: Boolean =
    chromosome.equalsIgnoreCase("X") ||
    chromosome.equalsIgnoreCase("chrX") ||
    chromosome == "23"

  /**
   * Check if this is autosomal (chromosomes 1-22).
   */
  def isAutosomal: Boolean = {
    val chr = chromosome.toLowerCase.stripPrefix("chr")
    chr.toIntOption.exists(n => n >= 1 && n <= 22)
  }

  /**
   * Get numeric genotype value (0 = hom ref, 1 = het, 2 = hom alt).
   * Returns -1 for no-call.
   *
   * @param referenceAllele The reference allele at this position
   */
  def numericGenotype(referenceAllele: Char): Int = {
    if (isNoCall) -1
    else {
      val refCount = Seq(allele1, allele2).count(_ == referenceAllele)
      2 - refCount // 0 if both ref, 1 if one ref, 2 if no ref
    }
  }
}

object GenotypeCall {
  /**
   * Create a no-call genotype.
   */
  def noCall(markerId: String, chromosome: String, position: Int): GenotypeCall =
    GenotypeCall(markerId, chromosome, position, '-', '-')

  /**
   * Parse genotype from a single string (e.g., "AA", "AG", "--").
   */
  def fromGenotype(markerId: String, chromosome: String, position: Int, genotype: String): GenotypeCall = {
    val cleaned = genotype.trim.toUpperCase
    if (cleaned.length == 2) {
      GenotypeCall(markerId, chromosome, position, cleaned(0), cleaned(1))
    } else if (cleaned.length == 1) {
      // Haploid call (Y, MT)
      GenotypeCall(markerId, chromosome, position, cleaned(0), cleaned(0))
    } else {
      noCall(markerId, chromosome, position)
    }
  }
}

/**
 * Summary statistics from a genotyping test.
 */
case class GenotypingTestSummary(
  testType: TestTypeDefinition,
  totalMarkersCalled: Int,
  totalMarkersPossible: Int,
  noCallCount: Int,
  noCallRate: Double,

  // Y-DNA marker coverage
  yMarkersCalled: Option[Int],
  yMarkersTotal: Option[Int],
  yCoverageRate: Option[Double],

  // mtDNA marker coverage
  mtMarkersCalled: Option[Int],
  mtMarkersTotal: Option[Int],
  mtCoverageRate: Option[Double],

  // Autosomal markers
  autosomalMarkersCalled: Int,

  // Quality indicators
  hetRate: Option[Double],
  chipVersion: Option[String],
  sourceFileHash: Option[String]
) derives Codec.AsObject {

  /**
   * Overall call rate.
   */
  def callRate: Double = if (totalMarkersPossible > 0) {
    totalMarkersCalled.toDouble / totalMarkersPossible
  } else 0.0

  /**
   * Check if quality is acceptable for ancestry analysis.
   */
  def isAcceptableForAncestry: Boolean =
    noCallRate < 0.05 && autosomalMarkersCalled >= 100000

  /**
   * Check if Y-DNA coverage is sufficient for haplogroup analysis.
   */
  def hasSufficientYCoverage: Boolean =
    yMarkersCalled.exists(_ >= 50) // At least 50 Y markers

  /**
   * Check if mtDNA coverage is sufficient for haplogroup analysis.
   */
  def hasSufficientMtCoverage: Boolean =
    mtMarkersCalled.exists(_ >= 20) // At least 20 mtDNA markers
}

object GenotypingTestSummary {
  /**
   * Compute summary from a list of genotype calls.
   */
  def fromCalls(
    calls: List[GenotypeCall],
    testType: TestTypeDefinition,
    chipVersion: Option[String] = None,
    sourceFileHash: Option[String] = None
  ): GenotypingTestSummary = {
    val totalCalls = calls.size
    val validCalls = calls.filterNot(_.isNoCall)
    val noCalls = calls.filter(_.isNoCall)

    val yCalls = calls.filter(_.isYChromosome)
    val yValid = yCalls.filterNot(_.isNoCall)

    val mtCalls = calls.filter(_.isMitochondrial)
    val mtValid = mtCalls.filterNot(_.isNoCall)

    val autosomalCalls = validCalls.filter(_.isAutosomal)
    val hetCalls = autosomalCalls.filter(_.isHeterozygous)

    GenotypingTestSummary(
      testType = testType,
      totalMarkersCalled = validCalls.size,
      totalMarkersPossible = totalCalls,
      noCallCount = noCalls.size,
      noCallRate = if (totalCalls > 0) noCalls.size.toDouble / totalCalls else 0.0,
      yMarkersCalled = if (yCalls.nonEmpty) Some(yValid.size) else None,
      yMarkersTotal = if (yCalls.nonEmpty) Some(yCalls.size) else None,
      yCoverageRate = if (yCalls.nonEmpty) Some(yValid.size.toDouble / yCalls.size) else None,
      mtMarkersCalled = if (mtCalls.nonEmpty) Some(mtValid.size) else None,
      mtMarkersTotal = if (mtCalls.nonEmpty) Some(mtCalls.size) else None,
      mtCoverageRate = if (mtCalls.nonEmpty) Some(mtValid.size.toDouble / mtCalls.size) else None,
      autosomalMarkersCalled = autosomalCalls.size,
      hetRate = if (autosomalCalls.nonEmpty) Some(hetCalls.size.toDouble / autosomalCalls.size) else None,
      chipVersion = chipVersion,
      sourceFileHash = sourceFileHash
    )
  }
}
