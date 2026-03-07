package com.decodingus.ibd.engine

import java.io.{BufferedReader, InputStreamReader}
import java.nio.file.{Files, Path}
import java.util.zip.GZIPInputStream
import scala.collection.mutable.ArrayBuffer
import scala.jdk.CollectionConverters.*
import scala.util.Using

/**
 * A genetic recombination map for converting physical positions (bp)
 * to genetic distances (cM) using linear interpolation.
 *
 * Uses HapMap-format genetic maps. Each chromosome has its own map.
 * The map stores sorted arrays of (position, cM) pairs for O(log n) lookup.
 */
class GeneticMap private (maps: Map[String, GeneticMap.ChromosomeMap]):

  /**
   * Convert a physical position to genetic distance in cM.
   * Uses linear interpolation between flanking map positions.
   *
   * @param chromosome Chromosome name (e.g., "1", "chr1")
   * @param position   Physical position in base pairs
   * @return Genetic position in centiMorgans, or None if chromosome not in map
   */
  def positionToCm(chromosome: String, position: Int): Option[Double] =
    val chr = GeneticMap.normalizeChromosome(chromosome)
    maps.get(chr).map(_.interpolate(position))

  /**
   * Convert a physical interval to genetic distance in cM.
   */
  def intervalCm(chromosome: String, startBp: Int, endBp: Int): Option[Double] =
    for
      startCm <- positionToCm(chromosome, startBp)
      endCm <- positionToCm(chromosome, endBp)
    yield math.abs(endCm - startCm)

  def hasChromosome(chromosome: String): Boolean =
    maps.contains(GeneticMap.normalizeChromosome(chromosome))

  def chromosomes: Set[String] = maps.keySet

object GeneticMap:

  /**
   * Sorted position/cM arrays for a single chromosome.
   * Uses binary search for O(log n) interpolation.
   */
  class ChromosomeMap(val positions: Array[Int], val cmValues: Array[Double]):
    require(positions.length == cmValues.length && positions.nonEmpty,
      "Positions and cM arrays must be non-empty and same length")

    def interpolate(position: Int): Double =
      val idx = java.util.Arrays.binarySearch(positions, position)
      if idx >= 0 then
        // Exact match
        cmValues(idx)
      else
        val insertionPoint = -(idx + 1)
        if insertionPoint == 0 then
          // Before first marker — extrapolate using first rate
          if positions.length >= 2 then
            val rate = (cmValues(1) - cmValues(0)) / (positions(1) - positions(0)).toDouble
            cmValues(0) + rate * (position - positions(0))
          else cmValues(0)
        else if insertionPoint >= positions.length then
          // After last marker — extrapolate using last rate
          if positions.length >= 2 then
            val n = positions.length
            val rate = (cmValues(n - 1) - cmValues(n - 2)) / (positions(n - 1) - positions(n - 2)).toDouble
            cmValues(n - 1) + rate * (position - positions(n - 1))
          else cmValues.last
        else
          // Between two markers — linear interpolation
          val lo = insertionPoint - 1
          val hi = insertionPoint
          val frac = (position - positions(lo)).toDouble / (positions(hi) - positions(lo)).toDouble
          cmValues(lo) + frac * (cmValues(hi) - cmValues(lo))

  /**
   * Normalize chromosome name to bare number (e.g., "chr1" → "1", "chrX" → "X").
   */
  def normalizeChromosome(chr: String): String =
    chr.toLowerCase.stripPrefix("chr") match
      case s if s.toIntOption.isDefined => s.toInt.toString
      case s => s.toUpperCase

  /**
   * Load a genetic map from HapMap-format files.
   *
   * Expected format (tab-separated, with header):
   * ```
   * Chromosome  Position(bp)  Rate(cM/Mb)  Map(cM)
   * chr1        55550         2.981822      0.000000
   * ```
   *
   * @param mapDir Directory containing per-chromosome map files
   * @param pattern File name pattern with %s for chromosome (e.g., "genetic_map_chr%s_b37.txt")
   * @param chromosomes Which chromosomes to load
   * @return Either error or loaded GeneticMap
   */
  def fromHapMapFiles(mapDir: Path, pattern: String,
                      chromosomes: Seq[String] = (1 to 22).map(_.toString)): Either[String, GeneticMap] =
    try
      val maps = chromosomes.flatMap { chr =>
        val fileName = pattern.format(chr)
        val filePath = mapDir.resolve(fileName)
        if Files.exists(filePath) then
          loadChromosomeMap(filePath).map(chr -> _).toOption
        else
          None
      }.toMap

      if maps.isEmpty then Left("No genetic map files found")
      else Right(new GeneticMap(maps))
    catch
      case e: Exception => Left(s"Failed to load genetic map: ${e.getMessage}")

  /**
   * Load a genetic map from gzipped HapMap-format files.
   */
  def fromHapMapGzFiles(mapDir: Path, pattern: String,
                        chromosomes: Seq[String] = (1 to 22).map(_.toString)): Either[String, GeneticMap] =
    try
      val maps = chromosomes.flatMap { chr =>
        val fileName = pattern.format(chr)
        val filePath = mapDir.resolve(fileName)
        if Files.exists(filePath) then
          loadChromosomeMapGz(filePath).map(chr -> _).toOption
        else
          None
      }.toMap

      if maps.isEmpty then Left("No genetic map files found")
      else Right(new GeneticMap(maps))
    catch
      case e: Exception => Left(s"Failed to load genetic map: ${e.getMessage}")

  /**
   * Create a uniform-rate genetic map for testing.
   * Assumes a constant recombination rate across all chromosomes.
   *
   * @param cmPerMb cM per megabase (default 1.0, roughly average human rate)
   * @param chromosomeLengths Map of chromosome → length in bp
   */
  def uniformRate(cmPerMb: Double = 1.0,
                  chromosomeLengths: Map[String, Int] = defaultChromosomeLengths): GeneticMap =
    val maps = chromosomeLengths.map { case (chr, length) =>
      val positions = Array(1, length)
      val cmValues = Array(0.0, length.toDouble / 1_000_000.0 * cmPerMb)
      chr -> ChromosomeMap(positions, cmValues)
    }
    new GeneticMap(maps)

  private def loadChromosomeMap(file: Path): Either[String, ChromosomeMap] =
    try
      val lines = Files.readAllLines(file).asScala
      parseMapLines(lines.iterator)
    catch
      case e: Exception => Left(s"Failed to read ${file.getFileName}: ${e.getMessage}")

  private def loadChromosomeMapGz(file: Path): Either[String, ChromosomeMap] =
    try
      Using(new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(file))))) { reader =>
        val lines = Iterator.continually(reader.readLine()).takeWhile(_ != null)
        parseMapLines(lines)
      }.get
    catch
      case e: Exception => Left(s"Failed to read ${file.getFileName}: ${e.getMessage}")

  private def parseMapLines(lines: Iterator[String]): Either[String, ChromosomeMap] =
    val positions = ArrayBuffer.empty[Int]
    val cmValues = ArrayBuffer.empty[Double]

    for line <- lines do
      val trimmed = line.trim
      if trimmed.nonEmpty && !trimmed.startsWith("#") && !trimmed.startsWith("Chromosome") && !trimmed.startsWith("chr") then
        val fields = trimmed.split("\\s+")
        if fields.length >= 4 then
          try
            val pos = fields(1).toInt
            val cm = fields(3).toDouble
            positions += pos
            cmValues += cm
          catch
            case _: NumberFormatException => // skip malformed lines
        else if fields.length >= 2 then
          // Some formats only have position and cM
          try
            val pos = fields(0).toInt
            val cm = fields(1).toDouble
            positions += pos
            cmValues += cm
          catch
            case _: NumberFormatException => // skip

    if positions.isEmpty then Left("No valid map entries found")
    else Right(ChromosomeMap(positions.toArray, cmValues.toArray))

  private def defaultChromosomeLengths: Map[String, Int] =
    Map(
      "1" -> 248956422, "2" -> 242193529, "3" -> 198295559,
      "4" -> 190214555, "5" -> 181538259, "6" -> 170805979,
      "7" -> 159345973, "8" -> 145138636, "9" -> 138394717,
      "10" -> 133797422, "11" -> 135086622, "12" -> 133275309,
      "13" -> 114364328, "14" -> 107043718, "15" -> 101991189,
      "16" -> 90338345, "17" -> 83257441, "18" -> 80373285,
      "19" -> 58617616, "20" -> 64444167, "21" -> 46709983,
      "22" -> 50818468
    )
