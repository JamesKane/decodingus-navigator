package com.decodingus.refgenome

import com.decodingus.refgenome.model.*

import scala.collection.mutable

/**
 * Types of genomic regions on the Y chromosome with their quality modifiers.
 *
 * Modifiers represent confidence in variant calls within each region type:
 * - 1.0 = reliable, high confidence
 * - <1.0 = reduced confidence due to structural characteristics
 *
 * Modifiers combine multiplicatively for overlapping regions.
 */
enum RegionType(val modifier: Double, val displayName: String, val description: String):
  case Cytoband extends RegionType(1.0, "Cytoband", "Chromosome band (display only)")
  case XDegenerate extends RegionType(1.0, "X-degenerate", "Stable, single-copy regions - gold standard")
  case Normal extends RegionType(1.0, "Normal", "Normal callable region")
  case PAR extends RegionType(0.5, "PAR", "Pseudoautosomal region - recombines with X")
  case Palindrome extends RegionType(0.4, "Palindrome", "Palindromic region - gene conversion risk")
  case XTR extends RegionType(0.3, "XTR", "X-transposed region - 99% X-identical, contamination risk")
  case Ampliconic extends RegionType(0.3, "Ampliconic", "Ampliconic region - high copy number, mapping artifacts")
  case STR extends RegionType(0.25, "STR", "Short tandem repeat - recLOH risk")
  case Centromere extends RegionType(0.1, "Centromere", "Centromeric region - nearly unmappable")
  case Heterochromatin extends RegionType(0.1, "Heterochromatin", "Heterochromatic region (Yq12) - unmappable")
  case NonCallable extends RegionType(0.5, "Non-callable", "Failed callable loci criteria")
  case LowDepth extends RegionType(0.7, "Low depth", "Read depth below threshold (<10x)")

/**
 * A genomic region with interval bounds and metadata.
 *
 * Uses 1-based, inclusive coordinates (matching VCF/GFF3 conventions).
 *
 * @param contig     Chromosome (e.g., "chrY")
 * @param start      Start position (1-based, inclusive)
 * @param end        End position (1-based, inclusive)
 * @param regionType Type of region affecting quality modifier
 * @param name       Optional name for display (e.g., "P8", "DYS389", "Yq11.223")
 */
case class GenomicRegion(
                          contig: String,
                          start: Long,
                          end: Long,
                          regionType: RegionType,
                          name: Option[String] = None
                        ) {
  def contains(position: Long): Boolean = position >= start && position <= end

  def displayDescription: String = name match {
    case Some(n) => s"$n ${regionType.displayName}"
    case None => regionType.displayName
  }
}

/**
 * Result of annotating a genomic position with region information.
 *
 * @param regions         All overlapping regions at the position
 * @param qualityModifier Combined modifier from all regions (multiplicative)
 * @param cytoband        The cytoband containing this position (for display)
 * @param primaryRegion   The most significant region affecting quality
 */
case class RegionAnnotation(
                             regions: List[GenomicRegion],
                             qualityModifier: Double,
                             cytoband: Option[GenomicRegion],
                             primaryRegion: Option[GenomicRegion]
                           ) {
  /**
   * Human-readable description of the annotation.
   * Prioritizes specific region names over generic types.
   */
  def description: String = {
    // Priority: specific named region > cytoband
    val regionDesc = primaryRegion.map(_.displayDescription)
    val cytobandDesc = cytoband.flatMap(_.name)

    (regionDesc, cytobandDesc) match {
      case (Some(r), Some(c)) => s"$r ($c)"
      case (Some(r), None) => r
      case (None, Some(c)) => c
      case (None, None) => "-"
    }
  }

  /**
   * Short display string for table columns.
   */
  def shortDisplay: String = {
    primaryRegion.map(r => r.name.getOrElse(r.regionType.displayName))
      .orElse(cytoband.flatMap(_.name))
      .getOrElse("-")
  }

  /**
   * Whether quality was adjusted from baseline.
   */
  def isAdjusted: Boolean = qualityModifier < 1.0

  /**
   * Get tooltip explaining the modifier.
   */
  def modifierTooltip: String = {
    if (!isAdjusted) "No quality adjustment"
    else {
      val modifiers = regions
        .filter(_.regionType.modifier < 1.0)
        .map(r => s"${r.regionType.displayName}: ${r.regionType.modifier}x")
      s"Modifiers: ${modifiers.mkString(", ")} → ${f"$qualityModifier%.2f"}x"
    }
  }
}

object RegionAnnotation {
  val empty: RegionAnnotation = RegionAnnotation(List.empty, 1.0, None, None)
}

/**
 * Annotator for Y chromosome genomic positions.
 *
 * Loads region files and provides fast lookup via binary search.
 * Regions are stored sorted by start position for efficient interval queries.
 *
 * @param cytobands         Cytoband regions (display only, no modifier)
 * @param palindromes       Palindromic regions P1-P8
 * @param strs              Short tandem repeat regions
 * @param pars              Pseudoautosomal regions PAR1/PAR2
 * @param xtrs              X-transposed regions
 * @param ampliconic        Ampliconic/multicopy regions
 * @param centromeres       Centromeric regions
 * @param heterochromatin   Heterochromatic regions (Yq12)
 * @param xdegenerate       X-degenerate regions
 * @param callablePositions Optional set of callable positions from callable_loci.bed
 */
class YRegionAnnotator(
                        cytobands: IndexedSeq[GenomicRegion] = IndexedSeq.empty,
                        palindromes: IndexedSeq[GenomicRegion] = IndexedSeq.empty,
                        strs: IndexedSeq[GenomicRegion] = IndexedSeq.empty,
                        pars: IndexedSeq[GenomicRegion] = IndexedSeq.empty,
                        xtrs: IndexedSeq[GenomicRegion] = IndexedSeq.empty,
                        ampliconic: IndexedSeq[GenomicRegion] = IndexedSeq.empty,
                        centromeres: IndexedSeq[GenomicRegion] = IndexedSeq.empty,
                        heterochromatin: IndexedSeq[GenomicRegion] = IndexedSeq.empty,
                        xdegenerate: IndexedSeq[GenomicRegion] = IndexedSeq.empty,
                        callablePositions: Option[Set[Long]] = None
                      ) {
  // All region sets for iteration (exclude cytobands from modifier calculation)
  private val modifierRegions: List[IndexedSeq[GenomicRegion]] = List(
    palindromes, strs, pars, xtrs, ampliconic, centromeres, heterochromatin
  )

  /**
   * Annotate a genomic position.
   *
   * @param contig   Chromosome (accepts "chrY", "Y", etc.)
   * @param position 1-based position
   * @param depth    Optional read depth for depth-based modifier
   * @return RegionAnnotation with all overlapping regions and combined modifier
   */
  def annotate(contig: String, position: Long, depth: Option[Int] = None): RegionAnnotation = {
    if (!isYChromosome(contig)) return RegionAnnotation.empty

    val overlapping = mutable.ListBuffer[GenomicRegion]()

    // Find cytoband (display only)
    val cytoband = findOverlapping(cytobands, position)
    cytoband.foreach(overlapping += _)

    // Find modifier regions
    val palindrome = findOverlapping(palindromes, position)
    palindrome.foreach(overlapping += _)

    val str = findOverlapping(strs, position)
    str.foreach(overlapping += _)

    val par = findOverlapping(pars, position)
    par.foreach(overlapping += _)

    val xtr = findOverlapping(xtrs, position)
    xtr.foreach(overlapping += _)

    val ampliconicRegion = findOverlapping(ampliconic, position)
    ampliconicRegion.foreach(overlapping += _)

    val centromere = findOverlapping(centromeres, position)
    centromere.foreach(overlapping += _)

    val heterochrom = findOverlapping(heterochromatin, position)
    heterochrom.foreach(overlapping += _)

    // Check callable status
    val nonCallable = callablePositions.exists(callable => !callable.contains(position))
    if (nonCallable) {
      overlapping += GenomicRegion("chrY", position, position, RegionType.NonCallable)
    }

    // Check depth
    val lowDepth = depth.exists(_ < 10)
    if (lowDepth) {
      overlapping += GenomicRegion("chrY", position, position, RegionType.LowDepth)
    }

    // Calculate combined modifier (multiplicative)
    val modifier = overlapping.toList
      .filter(_.regionType != RegionType.Cytoband) // Don't count cytobands in modifier
      .map(_.regionType.modifier)
      .foldLeft(1.0)(_ * _)

    // Determine primary region (lowest modifier = most impactful)
    val primaryRegion = overlapping.toList
      .filter(r => r.regionType != RegionType.Cytoband && r.regionType.modifier < 1.0)
      .sortBy(_.regionType.modifier)
      .headOption

    RegionAnnotation(overlapping.toList, modifier, cytoband, primaryRegion)
  }

  /**
   * Get quality modifier for a position (shortcut for just the modifier value).
   */
  def getQualityModifier(contig: String, position: Long, depth: Option[Int] = None): Double = {
    annotate(contig, position, depth).qualityModifier
  }

  /**
   * Binary search to find a region overlapping a position.
   *
   * Regions must be sorted by start position.
   */
  private def findOverlapping(regions: IndexedSeq[GenomicRegion], position: Long): Option[GenomicRegion] = {
    if (regions.isEmpty) return None

    // Binary search for the rightmost region with start <= position
    var lo = 0
    var hi = regions.length - 1
    var result: Option[GenomicRegion] = None

    while (lo <= hi) {
      val mid = lo + (hi - lo) / 2
      val region = regions(mid)

      if (region.start <= position) {
        // Check if this region contains the position
        if (position <= region.end) {
          result = Some(region)
        }
        lo = mid + 1
      } else {
        hi = mid - 1
      }
    }

    // Check the region at hi index if we haven't found a match
    if (result.isEmpty && hi >= 0 && hi < regions.length) {
      val region = regions(hi)
      if (region.start <= position && position <= region.end) {
        result = Some(region)
      }
    }

    result
  }

  /**
   * Find all overlapping regions at a position.
   */
  private def findAllOverlapping(regions: IndexedSeq[GenomicRegion], position: Long): List[GenomicRegion] = {
    // For simplicity, use linear scan after binary search finds first candidate
    // Could be optimized with interval tree if needed
    val result = mutable.ListBuffer[GenomicRegion]()

    var i = binarySearchLowerBound(regions, position)
    while (i < regions.length && regions(i).start <= position) {
      if (position <= regions(i).end) {
        result += regions(i)
      }
      i += 1
    }

    result.toList
  }

  /**
   * Binary search for lower bound (first region with start <= position).
   */
  private def binarySearchLowerBound(regions: IndexedSeq[GenomicRegion], position: Long): Int = {
    var lo = 0
    var hi = regions.length

    while (lo < hi) {
      val mid = lo + (hi - lo) / 2
      if (regions(mid).end < position) {
        lo = mid + 1
      } else {
        hi = mid
      }
    }

    lo
  }

  private def isYChromosome(contig: String): Boolean = {
    val normalized = contig.toLowerCase
    normalized == "chry" || normalized == "y" || normalized.startsWith("chry_") || normalized == "nc_000024"
  }

  // Statistics
  def cytobandCount: Int = cytobands.size

  def palindromeCount: Int = palindromes.size

  def strCount: Int = strs.size

  def parCount: Int = pars.size

  def xtrCount: Int = xtrs.size

  def ampliconicCount: Int = ampliconic.size

  def totalRegionCount: Int = cytobands.size + palindromes.size + strs.size + pars.size + xtrs.size + ampliconic.size + centromeres.size + heterochromatin.size

  /**
   * Get all regions grouped by type for visualization.
   * Useful for drawing chromosome ideograms with color-coded regions.
   */
  def getAllRegions: Map[RegionType, IndexedSeq[GenomicRegion]] = Map(
    RegionType.Cytoband -> cytobands,
    RegionType.Palindrome -> palindromes,
    RegionType.STR -> strs,
    RegionType.PAR -> pars,
    RegionType.XTR -> xtrs,
    RegionType.Ampliconic -> ampliconic,
    RegionType.Centromere -> centromeres,
    RegionType.Heterochromatin -> heterochromatin,
    RegionType.XDegenerate -> xdegenerate
  ).filter(_._2.nonEmpty)

  /**
   * Get chromosome length (max end position from all regions).
   * Falls back to CHM13v2 Y length if no regions loaded.
   */
  def getChromosomeLength: Long = {
    val allEnds = (cytobands ++ pars ++ heterochromatin ++ xdegenerate).map(_.end)
    if (allEnds.isEmpty) 62_460_029L else allEnds.max // Default to CHM13v2 Y length
  }
}

object YRegionAnnotator {
  /**
   * Create an annotator from parsed region files.
   */
  def fromRegions(
                   cytobands: List[GenomicRegion] = Nil,
                   palindromes: List[GenomicRegion] = Nil,
                   strs: List[GenomicRegion] = Nil,
                   pars: List[GenomicRegion] = Nil,
                   xtrs: List[GenomicRegion] = Nil,
                   ampliconic: List[GenomicRegion] = Nil,
                   centromeres: List[GenomicRegion] = Nil,
                   heterochromatin: List[GenomicRegion] = Nil,
                   xdegenerate: List[GenomicRegion] = Nil,
                   callablePositions: Option[Set[Long]] = None
                 ): YRegionAnnotator = {
    new YRegionAnnotator(
      cytobands = cytobands.sortBy(_.start).toIndexedSeq,
      palindromes = palindromes.sortBy(_.start).toIndexedSeq,
      strs = strs.sortBy(_.start).toIndexedSeq,
      pars = pars.sortBy(_.start).toIndexedSeq,
      xtrs = xtrs.sortBy(_.start).toIndexedSeq,
      ampliconic = ampliconic.sortBy(_.start).toIndexedSeq,
      centromeres = centromeres.sortBy(_.start).toIndexedSeq,
      heterochromatin = heterochromatin.sortBy(_.start).toIndexedSeq,
      xdegenerate = xdegenerate.sortBy(_.start).toIndexedSeq,
      callablePositions = callablePositions
    )
  }

  /**
   * Create an empty annotator (no region data loaded).
   */
  def empty: YRegionAnnotator = new YRegionAnnotator()

  /**
   * Convert GFF3 records to GenomicRegion list.
   */
  def gff3ToRegions(records: List[Gff3Record], regionType: RegionType): List[GenomicRegion] = {
    records.map { rec =>
      GenomicRegion(
        contig = rec.seqId,
        start = rec.start,
        end = rec.end,
        regionType = regionType,
        name = rec.name.orElse(rec.getAttribute("ID"))
      )
    }
  }

  /**
   * Convert BED records to GenomicRegion list.
   * BED uses 0-based, half-open coordinates; converts to 1-based inclusive.
   */
  def bedToRegions(records: List[BedRecord], regionType: RegionType): List[GenomicRegion] = {
    records.map { rec =>
      val (start1, end1) = RegionFileParser.bedToOneBased(rec.start, rec.end)
      GenomicRegion(
        contig = rec.chrom,
        start = start1,
        end = end1,
        regionType = regionType,
        name = rec.name
      )
    }
  }

  /**
   * GRCh38 Y chromosome heterochromatin (Yq12) hardcoded boundaries.
   * This region is ~30 Mbp of highly repetitive sequence, essentially unmappable.
   */
  val grch38Heterochromatin: List[GenomicRegion] = List(
    GenomicRegion("chrY", 26673237, 56887902, RegionType.Heterochromatin, Some("Yq12"))
  )

  /**
   * GRCh37 Y chromosome heterochromatin boundaries.
   */
  val grch37Heterochromatin: List[GenomicRegion] = List(
    GenomicRegion("chrY", 25294945, 59034049, RegionType.Heterochromatin, Some("Yq12"))
  )

  /**
   * CHM13v2.0 (T2T) Y chromosome heterochromatin boundaries.
   *
   * The T2T-Y assembly is 62,460,029 bp with complete heterochromatin sequence.
   * Yq12 contains ~30 Mbp of DYZ1, DYZ2, DYZ3 satellite arrays.
   *
   * Note: More precise boundaries available from chm13v2.0_censat_v2.1.bed
   * This is a conservative estimate based on T2T-Y publication.
   */
  val chm13v2Heterochromatin: List[GenomicRegion] = List(
    GenomicRegion("chrY", 26637971, 62122809, RegionType.Heterochromatin, Some("Yq12"))
  )

  /**
   * Create an annotator from centralized GenomeRegions API response.
   *
   * Converts the API's GenomeRegions model to YRegionAnnotator format,
   * extracting chrY-specific regions including cytobands, centromere,
   * heterochromatin, PAR, XTR, ampliconic, palindromes, and STR markers.
   *
   * @param regions           GenomeRegions from the API or bundled resource
   * @param callablePositions Optional set of callable positions from callable_loci.bed
   * @return YRegionAnnotator configured with all available region data
   */
  def fromGenomeRegions(regions: GenomeRegions, callablePositions: Option[Set[Long]] = None): YRegionAnnotator = {
    // Get chrY data if available
    val chrY = regions.chromosomes.get("chrY")
    if (chrY.isEmpty) {
      println(s"[YRegionAnnotator] No chrY data in GenomeRegions for ${regions.build}")
      return empty
    }

    val yData = chrY.get

    // Convert cytobands
    val cytobands = yData.cytobands.map { c =>
      GenomicRegion("chrY", c.start, c.end, RegionType.Cytoband, Some(c.name))
    }

    // Convert centromere
    val centromeres = yData.centromere.toList.map { r =>
      GenomicRegion("chrY", r.start, r.end, RegionType.Centromere, Some("centromere"))
    }

    // Convert Y-specific regions if available
    val (pars, xtrs, ampliconic, palindromes, heterochromatin, xdegenerate) = yData.regions match {
      case Some(yRegions) =>
        val parList = List(
          GenomicRegion("chrY", yRegions.par1.start, yRegions.par1.end, RegionType.PAR, Some("PAR1")),
          GenomicRegion("chrY", yRegions.par2.start, yRegions.par2.end, RegionType.PAR, Some("PAR2"))
        )

        val xtrList = yRegions.xtr.map { r =>
          GenomicRegion("chrY", r.start, r.end, RegionType.XTR, r.regionType)
        }

        val ampliconicList = yRegions.ampliconic.map { r =>
          GenomicRegion("chrY", r.start, r.end, RegionType.Ampliconic, r.regionType)
        }

        val palindromeList = yRegions.palindromes.map { nr =>
          GenomicRegion("chrY", nr.start, nr.end, RegionType.Palindrome, Some(nr.name))
        }

        val heterochromatinList = List(
          GenomicRegion("chrY", yRegions.heterochromatin.start, yRegions.heterochromatin.end,
            RegionType.Heterochromatin, Some("Yq12"))
        )

        val xdegList = yRegions.xDegenerate.map { r =>
          GenomicRegion("chrY", r.start, r.end, RegionType.XDegenerate, r.regionType)
        }

        (parList, xtrList, ampliconicList, palindromeList, heterochromatinList, xdegList)

      case None =>
        (Nil, Nil, Nil, Nil, Nil, Nil)
    }

    // Convert STR markers
    val strs = yData.strMarkers.getOrElse(Nil).map { marker =>
      GenomicRegion("chrY", marker.start, marker.end, RegionType.STR, Some(marker.name))
    }

    fromRegions(
      cytobands = cytobands,
      palindromes = palindromes,
      strs = strs,
      pars = pars,
      xtrs = xtrs,
      ampliconic = ampliconic,
      centromeres = centromeres,
      heterochromatin = heterochromatin,
      xdegenerate = xdegenerate,
      callablePositions = callablePositions
    )
  }
}
