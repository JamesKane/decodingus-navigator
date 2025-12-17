package com.decodingus.haplogroup.model

import com.decodingus.refgenome.{RegionAnnotation, RegionType, YRegionAnnotator}

/**
 * A variant call enriched with depth information and genomic region annotations.
 *
 * This data model adds display-only annotations to variant calls for reporting:
 * - Read depth at the call site
 * - Genomic region information (cytobands, palindromes, STRs, etc.)
 * - Adjusted quality based on region modifiers
 *
 * Note: Quality adjustments are for display only and don't affect haplogroup scoring.
 *
 * @param contig           Chromosome (e.g., "chrY")
 * @param position         1-based genomic position
 * @param ref              Reference allele
 * @param alt              Alternate allele
 * @param call             Called base at this position
 * @param quality          PHRED quality score from VCF
 * @param depth            Read depth at this position (from DP field)
 * @param regionAnnotation Region annotation for display and quality adjustment
 */
case class EnrichedVariantCall(
                                contig: String,
                                position: Long,
                                ref: String,
                                alt: String,
                                call: String,
                                quality: Option[Double],
                                depth: Option[Int],
                                regionAnnotation: Option[RegionAnnotation] = None
                              ) {
  /**
   * Get the adjusted quality based on region modifiers.
   * Returns None if no base quality available.
   */
  def adjustedQuality: Option[Double] = {
    quality.map { q =>
      val depthModifier = if (depth.exists(_ < 10)) RegionType.LowDepth.modifier else 1.0
      val regionModifier = regionAnnotation.map(_.qualityModifier).getOrElse(1.0)
      q * depthModifier * regionModifier
    }
  }

  /**
   * Convert PHRED quality to 5-star rating.
   * Uses adjusted quality if region annotation is available.
   */
  def qualityStars: String = qualityToStars(quality)

  /**
   * Get quality stars with adjustment applied.
   */
  def adjustedQualityStars: String = qualityToStars(adjustedQuality)

  /**
   * Whether quality was adjusted from baseline.
   */
  def isQualityAdjusted: Boolean = {
    val depthAdjusted = depth.exists(_ < 10)
    val regionAdjusted = regionAnnotation.exists(_.isAdjusted)
    depthAdjusted || regionAdjusted
  }

  /**
   * Display string for depth column (e.g., "32x" or "-").
   */
  def depthDisplay: String = depth.map(d => s"${d}x").getOrElse("-")

  /**
   * Short display string for region column.
   */
  def regionDisplay: String = regionAnnotation.map(_.shortDisplay).getOrElse("-")

  /**
   * Full region description for tooltips.
   */
  def regionDescription: String = regionAnnotation.map(_.description).getOrElse("-")

  /**
   * Tooltip explaining quality adjustment.
   */
  def qualityTooltip: String = {
    if (!isQualityAdjusted) {
      quality.map(q => f"Quality: Q=$q%.0f").getOrElse("No quality score")
    } else {
      val baseStars = qualityStars
      val adjStars = adjustedQualityStars
      val baseQ = quality.map(q => f"Q=$q%.0f").getOrElse("Q=?")

      val modifiers = scala.collection.mutable.ListBuffer[String]()
      if (depth.exists(_ < 10)) {
        modifiers += f"Low depth (${depth.get}x): ${RegionType.LowDepth.modifier}x"
      }
      regionAnnotation.foreach { ann =>
        ann.regions.filter(_.regionType.modifier < 1.0).foreach { r =>
          modifiers += s"${r.regionType.displayName}: ${r.regionType.modifier}x"
        }
      }

      val adjustedQ = adjustedQuality.map(q => f"Q=$q%.0f").getOrElse("Q=?")
      s"Base: $baseStars ($baseQ)\n${modifiers.mkString("\n")}\nAdjusted: $adjStars ($adjustedQ)"
    }
  }

  /**
   * Get the call state relative to ref/alt.
   */
  def callState: String = call match {
    case c if c.equalsIgnoreCase(alt) => "Derived"
    case c if c.equalsIgnoreCase(ref) => "Ancestral"
    case "-" | "" => "No Call"
    case _ => "Unknown"
  }

  /**
   * Convert PHRED quality to star rating.
   */
  private def qualityToStars(qual: Option[Double]): String = qual match {
    case None => "-"
    case Some(q) if q < 10 => "☆☆☆☆☆" // 0 stars
    case Some(q) if q < 20 => "★☆☆☆☆" // 1 star
    case Some(q) if q < 30 => "★★☆☆☆" // 2 stars
    case Some(q) if q < 40 => "★★★☆☆" // 3 stars
    case Some(q) if q < 50 => "★★★★☆" // 4 stars
    case Some(_) => "★★★★★" // 5 stars
  }
}

object EnrichedVariantCall {
  /**
   * Create an enriched variant call with region annotation.
   *
   * @param contig    Chromosome
   * @param position  1-based position
   * @param ref       Reference allele
   * @param alt       Alternate allele
   * @param call      Called base
   * @param quality   PHRED quality score
   * @param depth     Read depth
   * @param annotator Optional region annotator for region info
   */
  def create(
              contig: String,
              position: Long,
              ref: String,
              alt: String,
              call: String,
              quality: Option[Double] = None,
              depth: Option[Int] = None,
              annotator: Option[YRegionAnnotator] = None
            ): EnrichedVariantCall = {
    val regionAnnotation = annotator.map(_.annotate(contig, position, depth))
    EnrichedVariantCall(contig, position, ref, alt, call, quality, depth, regionAnnotation)
  }

  /**
   * Quality modifier legend for reports.
   */
  val qualityModifierLegend: String =
    """Region Quality Modifiers:
      |  X-degenerate: 1.0x (reliable)    PAR: 0.5x           Palindrome: 0.4x
      |  XTR: 0.3x                        Ampliconic: 0.3x    STR: 0.25x
      |  Centromere: 0.1x                 Heterochromatin: 0.1x
      |  Low depth (<10x): 0.7x           Non-callable: 0.5x""".stripMargin
}

/**
 * Container for enriched variant data for a haplogroup analysis.
 *
 * @param snpCalls       Map of position to enriched call for SNPs on the haplogroup path
 * @param novelSnps      Enriched novel/private SNP variants
 * @param novelIndels    Enriched novel/private indel variants
 * @param referenceBuild Reference build used for annotation
 */
case class EnrichedVariantData(
                                snpCalls: Map[Long, EnrichedVariantCall],
                                novelSnps: List[EnrichedVariantCall],
                                novelIndels: List[EnrichedVariantCall],
                                referenceBuild: String
                              ) {
  def hasRegionAnnotations: Boolean = {
    snpCalls.values.exists(_.regionAnnotation.isDefined) ||
      novelSnps.exists(_.regionAnnotation.isDefined) ||
      novelIndels.exists(_.regionAnnotation.isDefined)
  }

  def hasDepthData: Boolean = {
    snpCalls.values.exists(_.depth.isDefined) ||
      novelSnps.exists(_.depth.isDefined) ||
      novelIndels.exists(_.depth.isDefined)
  }
}

object EnrichedVariantData {
  val empty: EnrichedVariantData = EnrichedVariantData(Map.empty, Nil, Nil, "unknown")
}
