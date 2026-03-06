package com.decodingus.haplogroup.report

import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupResult, Locus, NamedVariant}
import com.decodingus.haplogroup.processor.{PrivateVariant, SnpCallInfo}
import com.decodingus.haplogroup.tree.{TreeProviderType, TreeType}
import com.decodingus.haplogroup.vendor.NamedVariantCache
import com.decodingus.refgenome.{StrAnnotator, YRegionAnnotator}

import java.io.{File, PrintWriter}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import scala.util.Using

/**
 * Generates Yleaf-style haplogroup reports.
 */
object HaplogroupReportWriter {

  /**
   * Convert a PHRED quality score to a 5-star rating.
   * 0 stars = shouldn't be included (Q < 10)
   * 1 star  = very low confidence (Q 10-19)
   * 2 stars = low confidence (Q 20-29)
   * 3 stars = moderate confidence (Q 30-39)
   * 4 stars = high confidence (Q 40-49)
   * 5 stars = very high confidence (Q >= 50)
   */
  private def qualityToStars(quality: Option[Double]): String = {
    quality match {
      case None => "-"
      case Some(q) if q < 10 => "☆☆☆☆☆" // 0 stars (empty)
      case Some(q) if q < 20 => "★☆☆☆☆" // 1 star
      case Some(q) if q < 30 => "★★☆☆☆" // 2 stars
      case Some(q) if q < 40 => "★★★☆☆" // 3 stars
      case Some(q) if q < 50 => "★★★★☆" // 4 stars
      case Some(_) => "★★★★★" // 5 stars
    }
  }

  /**
   * Format variant aliases for display, showing rsIds and common names.
   * Truncates if too long for the column width.
   */
  private def formatVariantAliases(variant: NamedVariant): String = {
    val parts = scala.collection.mutable.ListBuffer[String]()

    // Add first rsId if available
    variant.aliases.rsIds.headOption.foreach(parts += _)

    // Add common names (excluding the canonical name which is shown separately)
    val otherNames = variant.aliases.commonNames.filterNot(n =>
      variant.canonicalName.contains(n)
    ).take(2)
    parts ++= otherNames

    if (parts.isEmpty) {
      "-"
    } else {
      val result = parts.mkString(", ")
      if (result.length > 23) result.take(20) + "..." else result
    }
  }

  /**
   * Format depth for display.
   */
  private def formatDepth(depth: Option[Int]): String = {
    depth.map(d => s"${d}x").getOrElse("-")
  }

  /**
   * Calculate adjusted quality stars based on depth and region modifiers.
   */
  private def adjustedQualityToStars(
                                      quality: Option[Double],
                                      depth: Option[Int],
                                      regionModifier: Double
                                    ): String = {
    quality match {
      case None => "-"
      case Some(q) =>
        val depthModifier = if (depth.exists(_ < 10)) 0.7 else 1.0
        val adjustedQ = q * depthModifier * regionModifier
        qualityToStars(Some(adjustedQ))
    }
  }

  /**
   * Write a haplogroup analysis report to the specified directory.
   *
   * @param outputDir         Directory to write the report
   * @param treeType          Y-DNA or MT-DNA
   * @param results           Scored haplogroup results
   * @param tree              The haplogroup tree used for analysis
   * @param snpCalls          The SNP calls from the VCF
   * @param sampleName        Optional sample name
   * @param privateVariants   Optional list of private/novel variants
   * @param treeProvider      Optional tree provider used for analysis
   * @param strAnnotator      Optional STR annotator for indel classification
   * @param sampleBuild       Optional reference build of the sample BAM/CRAM
   * @param treeBuild         Optional reference build of the tree coordinates
   * @param namedVariantCache Optional cache for looking up variant aliases and rsIds
   * @param snpCallInfo        Optional map of position to full SNP call info (quality, depth)
   * @param yRegionAnnotator   Optional Y chromosome region annotator for region info
   * @param expectedYCoverage  Optional expected Y chromosome coverage for excessive depth detection
   */
  def writeReport(
                   outputDir: File,
                   treeType: TreeType,
                   results: List[HaplogroupResult],
                   tree: List[Haplogroup],
                   snpCalls: Map[Long, String],
                   sampleName: Option[String] = None,
                   privateVariants: Option[List[PrivateVariant]] = None,
                   treeProvider: Option[TreeProviderType] = None,
                   strAnnotator: Option[StrAnnotator] = None,
                   sampleBuild: Option[String] = None,
                   treeBuild: Option[String] = None,
                   namedVariantCache: Option[NamedVariantCache] = None,
                   snpCallInfo: Option[Map[Long, SnpCallInfo]] = None,
                   yRegionAnnotator: Option[YRegionAnnotator] = None,
                   expectedYCoverage: Option[Double] = None
                 ): File = {
    outputDir.mkdirs()

    val prefix = treeType match {
      case TreeType.YDNA => "ydna"
      case TreeType.MTDNA => "mtdna"
    }

    val reportFile = new File(outputDir, s"${prefix}_haplogroup_report.txt")

    Using.resource(new PrintWriter(reportFile)) { writer =>
      val timestamp = LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss"))
      val dnaType = if (treeType == TreeType.YDNA) "Y-DNA" else "MT-DNA"
      val treeProviderName = treeProvider match {
        case Some(TreeProviderType.FTDNA) => "FTDNA"
        case Some(TreeProviderType.DECODINGUS) => "Decoding Us"
        case None => "Unknown"
      }

      // Pre-compute result lookup map - O(n) instead of O(n²) for path lookups
      val resultsByName: Map[String, HaplogroupResult] = results.map(r => r.name -> r).toMap

      writer.println("=" * 80)
      writer.println(s"  $dnaType Haplogroup Analysis Report")
      writer.println("=" * 80)
      writer.println()
      writer.println(s"Generated: $timestamp")
      writer.println(s"Tree Provider: $treeProviderName")
      treeBuild.foreach(build => writer.println(s"Tree Build: $build"))
      sampleBuild.foreach(build => writer.println(s"Sample Build: $build"))
      // Show liftover status
      (sampleBuild, treeBuild) match {
        case (Some(sample), Some(tree)) if sample != tree =>
          writer.println(s"Liftover: Yes ($tree → $sample for calling, $sample → $tree for scoring)")
        case (Some(sample), Some(tree)) =>
          writer.println(s"Liftover: No (native $sample coordinates)")
        case _ =>
        // Don't show liftover status if builds are unknown
      }
      sampleName.foreach(name => writer.println(s"Sample: $name"))
      writer.println()

      // Top prediction
      val topResult = results.headOption
      writer.println("-" * 80)
      writer.println("HAPLOGROUP PREDICTION")
      writer.println("-" * 80)
      topResult match {
        case Some(result) =>
          writer.println(s"  Predicted Haplogroup: ${result.name}")
          writer.println(s"  Score: ${result.score}")
          writer.println(s"  Derived SNPs: ${result.matchingSnps}")
          writer.println(s"  Ancestral SNPs: ${result.ancestralMatches}")
          writer.println(s"  No Calls: ${result.noCalls}")
          writer.println(s"  Tree Depth: ${result.depth}")
        case None =>
          writer.println("  No haplogroup could be determined.")
      }
      writer.println()

      // Top 10 candidates
      writer.println("-" * 80)
      writer.println("TOP 10 CANDIDATES")
      writer.println("-" * 80)
      writer.println(f"${"Rank"}%5s  ${"Haplogroup"}%-25s  ${"Score"}%8s  ${"Derived"}%8s  ${"Ancestral"}%10s  ${"Depth"}%6s")
      writer.println("-" * 80)
      results.take(10).zipWithIndex.foreach { case (result, idx) =>
        writer.println(f"${idx + 1}%5d  ${result.name}%-25s  ${result.score}%8.1f  ${result.matchingSnps}%8d  ${result.ancestralMatches}%10d  ${result.depth}%6d")
      }
      writer.println()

      // Build path to predicted haplogroup (compute once, reuse)
      val pathOpt = topResult.map(top => findPathToHaplogroup(tree, top.name))

      pathOpt.foreach { path =>
        writer.println("-" * 80)
        writer.println("HAPLOGROUP PATH")
        writer.println("-" * 80)
        path.foreach { haplo =>
          val indent = "  " * haplo.depth
          val result = resultsByName.get(haplo.name)
          val parentDerived = haplo.parent.flatMap(resultsByName.get).map(_.matchingSnps).getOrElse(0)
          val scoreInfo = result.map(r => s"[+${r.matchingSnps - parentDerived} derived]").getOrElse("")
          writer.println(s"$indent${haplo.name} $scoreInfo")
        }
        writer.println()
      }

      // SNP details for top haplogroup path
      pathOpt.foreach { path =>
        writer.println("-" * 80)
        writer.println("SNP DETAILS (along predicted path)")
        writer.println("-" * 80)

        val allLociOnPath = path.flatMap(_.loci)

        // Check if we have named variant data and determine which build to use for lookups
        val variantLookupBuild = namedVariantCache.flatMap { cache =>
          // Try treeBuild first, then sampleBuild, then common builds
          treeBuild.orElse(sampleBuild).orElse(Some("GRCh38"))
        }

        // Determine what columns to show
        val showAliases = namedVariantCache.isDefined && namedVariantCache.exists(_.isLoaded)
        val showDepthRegion = snpCallInfo.isDefined || yRegionAnnotator.isDefined

        // Print header based on available columns
        if (showDepthRegion) {
          if (showAliases) {
            writer.println(f"${"Contig"}%6s  ${"Position"}%12s  ${"SNP Name"}%-12s  ${"Aliases"}%-18s  ${"Anc"}%4s  ${"Der"}%4s  ${"Call"}%4s  ${"State"}%-10s  ${"Depth"}%6s  ${"Region"}%-15s  ${"Quality"}%-10s")
          } else {
            writer.println(f"${"Contig"}%6s  ${"Position"}%12s  ${"SNP Name"}%-15s  ${"Anc"}%4s  ${"Der"}%4s  ${"Call"}%4s  ${"State"}%-10s  ${"Depth"}%6s  ${"Region"}%-15s  ${"Quality"}%-10s")
          }
        } else {
          if (showAliases) {
            writer.println(f"${"Contig"}%6s  ${"Position"}%12s  ${"SNP Name"}%-15s  ${"Aliases/rsIds"}%-25s  ${"Anc"}%5s  ${"Der"}%5s  ${"Call"}%5s  ${"State"}%10s")
          } else {
            writer.println(f"${"Contig"}%6s  ${"Position"}%12s  ${"SNP Name"}%-20s  ${"Anc"}%5s  ${"Der"}%5s  ${"Call"}%5s  ${"State"}%10s")
          }
        }
        writer.println("-" * (if (showDepthRegion) 120 else 80))

        allLociOnPath.sortBy(l => (l.contig, l.position)).foreach { locus =>
          val called = snpCalls.get(locus.position).getOrElse("-")
          val state = snpCalls.get(locus.position) match {
            case Some(base) if base.equalsIgnoreCase(locus.alt) => "Derived"
            case Some(base) if base.equalsIgnoreCase(locus.ref) => "Ancestral"
            case Some(_) => "Unknown"
            case None => "No Call"
          }

          // Get call info (quality, depth) if available
          val callInfo = snpCallInfo.flatMap(_.get(locus.position))
          val depth = callInfo.flatMap(_.depth)
          val quality = callInfo.flatMap(_.quality)

          // Get region annotation if available
          val regionAnn = yRegionAnnotator.map(_.annotate(locus.contig, locus.position, depth, expectedYCoverage))
          val regionDisplay = regionAnn.map(_.shortDisplay).getOrElse("-")
          val regionModifier = regionAnn.map(_.qualityModifier).getOrElse(1.0)

          // Format quality - show adjusted if there's a modifier
          val qualityDisplay = if (regionModifier < 1.0 || depth.exists(_ < 10)) {
            val adjStars = adjustedQualityToStars(quality, depth, regionModifier)
            s"$adjStars (adj)"
          } else {
            qualityToStars(quality)
          }

          if (showDepthRegion) {
            val depthStr = formatDepth(depth)
            if (showAliases) {
              val aliases = variantLookupBuild.flatMap { build =>
                namedVariantCache.flatMap(_.getByPosition(build, locus.position))
              }.map(formatVariantAliases).getOrElse("-")
              val aliasShort = if (aliases.length > 16) aliases.take(14) + ".." else aliases
              writer.println(f"${locus.contig}%6s  ${locus.position}%12d  ${locus.name}%-12s  $aliasShort%-18s  ${locus.ref}%4s  ${locus.alt}%4s  $called%4s  $state%-10s  $depthStr%6s  $regionDisplay%-15s  $qualityDisplay%-10s")
            } else {
              writer.println(f"${locus.contig}%6s  ${locus.position}%12d  ${locus.name}%-15s  ${locus.ref}%4s  ${locus.alt}%4s  $called%4s  $state%-10s  $depthStr%6s  $regionDisplay%-15s  $qualityDisplay%-10s")
            }
          } else {
            if (showAliases) {
              val aliases = variantLookupBuild.flatMap { build =>
                namedVariantCache.flatMap(_.getByPosition(build, locus.position))
              }.map(formatVariantAliases).getOrElse("-")
              writer.println(f"${locus.contig}%6s  ${locus.position}%12d  ${locus.name}%-15s  $aliases%-25s  ${locus.ref}%5s  ${locus.alt}%5s  $called%5s  $state%10s")
            } else {
              writer.println(f"${locus.contig}%6s  ${locus.position}%12d  ${locus.name}%-20s  ${locus.ref}%5s  ${locus.alt}%5s  $called%5s  $state%10s")
            }
          }
        }
        writer.println()
      }

      // Novel/Unplaced Variants section - separated into SNPs and Indels
      privateVariants.foreach { variants =>
        // Separate SNPs from indels (potential STRs)
        val (snpVariants, indelVariants) = variants.partition { v =>
          v.ref.length == 1 && v.alt.length == 1
        }

        // Check if we have region annotation available
        val showDepthRegionNovel = yRegionAnnotator.isDefined || snpVariants.exists(_.depth.isDefined)

        // SNPs section
        writer.println("-" * (if (showDepthRegionNovel) 120 else 80))
        writer.println("NOVEL/UNPLACED SNPs")
        writer.println("-" * (if (showDepthRegionNovel) 120 else 80))
        if (snpVariants.isEmpty) {
          writer.println("  No novel SNPs discovered.")
        } else {
          writer.println(s"  Total novel/unplaced SNPs: ${snpVariants.size}")
          writer.println()
          if (showDepthRegionNovel) {
            writer.println(f"${"Contig"}%6s  ${"Position"}%12s  ${"Ref"}%5s  ${"Alt"}%5s  ${"Depth"}%6s  ${"Region"}%-15s  ${"Quality"}%-10s")
          } else {
            writer.println(f"${"Contig"}%6s  ${"Position"}%12s  ${"Ref"}%5s  ${"Alt"}%5s  ${"Quality"}%-10s")
          }
          writer.println("-" * (if (showDepthRegionNovel) 120 else 80))
          snpVariants.sortBy(v => (v.contig, v.position)).foreach { v =>
            if (showDepthRegionNovel) {
              val regionAnn = yRegionAnnotator.map(_.annotate(v.contig, v.position, v.depth, expectedYCoverage))
              val regionDisplay = regionAnn.map(_.shortDisplay).getOrElse("-")
              val regionModifier = regionAnn.map(_.qualityModifier).getOrElse(1.0)
              val qualityDisplay = if (regionModifier < 1.0 || v.depth.exists(_ < 10)) {
                val adjStars = adjustedQualityToStars(v.quality, v.depth, regionModifier)
                s"$adjStars (adj)"
              } else {
                qualityToStars(v.quality)
              }
              writer.println(f"${v.contig}%6s  ${v.position}%12d  ${v.ref}%5s  ${v.alt}%5s  ${formatDepth(v.depth)}%6s  $regionDisplay%-15s  $qualityDisplay%-10s")
            } else {
              val stars = qualityToStars(v.quality)
              writer.println(f"${v.contig}%6s  ${v.position}%12d  ${v.ref}%5s  ${v.alt}%5s  $stars%-10s")
            }
          }
        }
        writer.println()

        // Indels section (potential STRs) - annotate with HipSTR if available
        val showDepthRegionIndel = yRegionAnnotator.isDefined || indelVariants.exists(_.depth.isDefined)
        writer.println("-" * (if (showDepthRegionIndel) 120 else 80))
        writer.println("NOVEL/UNPLACED INDELS")
        writer.println("-" * (if (showDepthRegionIndel) 120 else 80))
        if (indelVariants.isEmpty) {
          writer.println("  No novel indels discovered.")
        } else {
          // Separate into known STRs and other indels
          val (strIndels, otherIndels) = strAnnotator match {
            case Some(annotator) =>
              indelVariants.partition { v =>
                annotator.findOverlapping(v.contig, v.position).isDefined
              }
            case None =>
              (List.empty[PrivateVariant], indelVariants)
          }

          writer.println(s"  Total novel/unplaced indels: ${indelVariants.size}")
          if (strAnnotator.isDefined) {
            writer.println(s"    - In known STR regions: ${strIndels.size}")
            writer.println(s"    - Other indels: ${otherIndels.size}")
          }
          writer.println()

          // STR indels with annotation
          if (strIndels.nonEmpty) {
            writer.println("  Known STR Regions:")
            if (showDepthRegionIndel) {
              writer.println(f"${"Contig"}%6s  ${"Position"}%12s  ${"Ref"}%-10s  ${"Alt"}%-10s  ${"Depth"}%6s  ${"Region"}%-12s  ${"Quality"}%-10s  ${"STR Type"}%-20s")
            } else {
              writer.println(f"${"Contig"}%6s  ${"Position"}%12s  ${"Ref"}%-10s  ${"Alt"}%-10s  ${"Quality"}%-10s  ${"STR Type"}%-20s")
            }
            writer.println("-" * (if (showDepthRegionIndel) 120 else 80))
            strIndels.sortBy(v => (v.contig, v.position)).foreach { v =>
              val refDisplay = if (v.ref.length > 8) v.ref.take(6) + ".." else v.ref
              val altDisplay = if (v.alt.length > 8) v.alt.take(6) + ".." else v.alt
              val strInfo = strAnnotator.flatMap(_.findOverlapping(v.contig, v.position))
                .map(r => strAnnotator.get.formatAnnotation(r))
                .getOrElse("-")
              if (showDepthRegionIndel) {
                val regionAnn = yRegionAnnotator.map(_.annotate(v.contig, v.position, v.depth, expectedYCoverage))
                val regionDisplay = regionAnn.map(_.shortDisplay).getOrElse("-")
                val regionModifier = regionAnn.map(_.qualityModifier).getOrElse(1.0)
                val qualityDisplay = if (regionModifier < 1.0 || v.depth.exists(_ < 10)) {
                  val adjStars = adjustedQualityToStars(v.quality, v.depth, regionModifier)
                  s"$adjStars (adj)"
                } else {
                  qualityToStars(v.quality)
                }
                writer.println(f"${v.contig}%6s  ${v.position}%12d  $refDisplay%-10s  $altDisplay%-10s  ${formatDepth(v.depth)}%6s  $regionDisplay%-12s  $qualityDisplay%-10s  $strInfo%-20s")
              } else {
                val stars = qualityToStars(v.quality)
                writer.println(f"${v.contig}%6s  ${v.position}%12d  $refDisplay%-10s  $altDisplay%-10s  $stars%-10s  $strInfo%-20s")
              }
            }
            writer.println()
          }

          // Other indels
          if (otherIndels.nonEmpty) {
            if (strIndels.nonEmpty) {
              writer.println("  Other Indels:")
            }
            if (showDepthRegionIndel) {
              writer.println(f"${"Contig"}%6s  ${"Position"}%12s  ${"Ref"}%-12s  ${"Alt"}%-12s  ${"Depth"}%6s  ${"Region"}%-12s  ${"Quality"}%-10s")
            } else {
              writer.println(f"${"Contig"}%6s  ${"Position"}%12s  ${"Ref"}%-12s  ${"Alt"}%-12s  ${"Quality"}%-10s")
            }
            writer.println("-" * (if (showDepthRegionIndel) 120 else 80))
            otherIndels.sortBy(v => (v.contig, v.position)).foreach { v =>
              val refDisplay = if (v.ref.length > 10) v.ref.take(8) + ".." else v.ref
              val altDisplay = if (v.alt.length > 10) v.alt.take(8) + ".." else v.alt
              if (showDepthRegionIndel) {
                val regionAnn = yRegionAnnotator.map(_.annotate(v.contig, v.position, v.depth, expectedYCoverage))
                val regionDisplay = regionAnn.map(_.shortDisplay).getOrElse("-")
                val regionModifier = regionAnn.map(_.qualityModifier).getOrElse(1.0)
                val qualityDisplay = if (regionModifier < 1.0 || v.depth.exists(_ < 10)) {
                  val adjStars = adjustedQualityToStars(v.quality, v.depth, regionModifier)
                  s"$adjStars (adj)"
                } else {
                  qualityToStars(v.quality)
                }
                writer.println(f"${v.contig}%6s  ${v.position}%12d  $refDisplay%-12s  $altDisplay%-12s  ${formatDepth(v.depth)}%6s  $regionDisplay%-12s  $qualityDisplay%-10s")
              } else {
                val stars = qualityToStars(v.quality)
                writer.println(f"${v.contig}%6s  ${v.position}%12d  $refDisplay%-12s  $altDisplay%-12s  $stars%-10s")
              }
            }
          }
        }
        writer.println()
      }

      // Summary statistics
      writer.println("-" * 80)
      writer.println("SUMMARY STATISTICS")
      writer.println("-" * 80)
      // Use results.size as proxy for tree size - avoids full tree traversal
      writer.println(s"  Haplogroups evaluated: ${results.size}")
      writer.println(s"  SNPs with calls: ${snpCalls.size}")
      pathOpt.foreach { path =>
        val pathLociCount = path.map(_.loci.size).sum
        writer.println(s"  SNPs on predicted path: $pathLociCount")
      }
      privateVariants.foreach { variants =>
        val (snps, indels) = variants.partition(v => v.ref.length == 1 && v.alt.length == 1)
        writer.println(s"  Novel/unplaced SNPs: ${snps.size}")
        writer.println(s"  Novel/unplaced indels: ${indels.size}")
      }
      writer.println()

      // Add quality modifier legend if region annotations were used
      if (yRegionAnnotator.isDefined) {
        writer.println("-" * 80)
        writer.println("QUALITY MODIFIER LEGEND")
        writer.println("-" * 80)
        writer.println("  Quality stars marked '(adj)' have been adjusted for region or depth:")
        writer.println()
        writer.println("  Region Modifiers:")
        writer.println("    X-degenerate: 1.0x (reliable)    PAR: 0.5x           Palindrome: 0.4x")
        writer.println("    XTR: 0.3x                        Ampliconic: 0.3x    STR: 0.25x")
        writer.println("    Centromere: 0.1x                 Heterochromatin: 0.1x")
        writer.println()
        writer.println("  Depth Modifiers:")
        writer.println("    Low depth (<10x): 0.7x           Non-callable: 0.5x")
        writer.println("    Excessive depth (2-3x): 0.4x     Highly excessive (3-4x): 0.2x")
        writer.println("    Extreme depth (>4x): 0.1x")
        writer.println()
        writer.println("  Note: Modifiers combine multiplicatively. Adjustments affect display only,")
        writer.println("        not internal haplogroup scoring.")
        writer.println()
      }

      writer.println("=" * 80)
    }

    reportFile
  }

  private case class HaplogroupWithDepth(name: String, parent: Option[String], loci: List[Locus], depth: Int)

  private def findPathToHaplogroup(tree: List[Haplogroup], targetName: String): List[HaplogroupWithDepth] = {
    def search(haplogroup: Haplogroup, depth: Int): Option[List[HaplogroupWithDepth]] = {
      val current = HaplogroupWithDepth(haplogroup.name, haplogroup.parent, haplogroup.loci, depth)
      if (haplogroup.name == targetName) {
        Some(List(current))
      } else {
        haplogroup.children.flatMap(child => search(child, depth + 1)).headOption.map(path => current :: path)
      }
    }

    tree.flatMap(root => search(root, 0)).headOption.getOrElse(List.empty)
  }

  /**
   * Escape a string for CSV output (handles commas, quotes, newlines).
   */
  private def escapeCsv(s: String): String = {
    if (s.contains(",") || s.contains("\"") || s.contains("\n") || s.contains("\r")) {
      "\"" + s.replace("\"", "\"\"") + "\""
    } else {
      s
    }
  }

  /**
   * Write a CSV export of haplogroup analysis results with region annotations.
   *
   * Creates three CSV files:
   * - {prefix}_snp_details.csv - SNPs along the predicted haplogroup path
   * - {prefix}_novel_snps.csv - Novel/unplaced SNPs
   * - {prefix}_novel_indels.csv - Novel/unplaced indels
   *
   * @param outputDir        Directory to write the CSV files
   * @param treeType         Y-DNA or MT-DNA
   * @param results          Scored haplogroup results
   * @param tree             The haplogroup tree used for analysis
   * @param snpCalls         The SNP calls from the VCF
   * @param sampleName       Optional sample name
   * @param privateVariants  Optional list of private/novel variants
   * @param snpCallInfo       Optional map of position to full SNP call info (quality, depth)
   * @param yRegionAnnotator  Optional Y chromosome region annotator for region info
   * @param strAnnotator      Optional STR annotator for indel classification
   * @param expectedYCoverage Optional expected Y chromosome coverage for excessive depth detection
   * @return List of created CSV files
   */
  def writeCsvReport(
                      outputDir: File,
                      treeType: TreeType,
                      results: List[HaplogroupResult],
                      tree: List[Haplogroup],
                      snpCalls: Map[Long, String],
                      sampleName: Option[String] = None,
                      privateVariants: Option[List[PrivateVariant]] = None,
                      snpCallInfo: Option[Map[Long, SnpCallInfo]] = None,
                      yRegionAnnotator: Option[YRegionAnnotator] = None,
                      strAnnotator: Option[StrAnnotator] = None,
                      expectedYCoverage: Option[Double] = None
                    ): List[File] = {
    outputDir.mkdirs()

    val prefix = treeType match {
      case TreeType.YDNA => "ydna"
      case TreeType.MTDNA => "mtdna"
    }

    val createdFiles = scala.collection.mutable.ListBuffer[File]()

    // Build path to predicted haplogroup
    val topResult = results.headOption
    val pathOpt = topResult.map(top => findPathToHaplogroup(tree, top.name))

    // SNP Details CSV
    pathOpt.foreach { path =>
      val snpDetailsFile = new File(outputDir, s"${prefix}_snp_details.csv")
      Using.resource(new PrintWriter(snpDetailsFile)) { writer =>
        // Header
        writer.println("contig,position,snp_name,haplogroup,ref,alt,call,state,depth,quality,region_type,region_name,quality_modifier,adjusted_quality")

        val allLociOnPath = path.flatMap { hp =>
          hp.loci.map(l => (l, hp.name))
        }

        allLociOnPath.sortBy { case (l, _) => (l.contig, l.position) }.foreach { case (locus, haplogroup) =>
          val called = snpCalls.get(locus.position).getOrElse("")
          val state = snpCalls.get(locus.position) match {
            case Some(base) if base.equalsIgnoreCase(locus.alt) => "Derived"
            case Some(base) if base.equalsIgnoreCase(locus.ref) => "Ancestral"
            case Some(_) => "Unknown"
            case None => "NoCall"
          }

          val callInfo = snpCallInfo.flatMap(_.get(locus.position))
          val depth = callInfo.flatMap(_.depth)
          val quality = callInfo.flatMap(_.quality)

          val regionAnn = yRegionAnnotator.map(_.annotate(locus.contig, locus.position, depth, expectedYCoverage))
          val regionType = regionAnn.map(_.primaryRegion.map(_.regionType.toString).getOrElse("")).getOrElse("")
          val regionName = regionAnn.map(_.shortDisplay).getOrElse("")
          val regionModifier = regionAnn.map(_.qualityModifier).getOrElse(1.0)

          val depthModifier = if (depth.exists(_ < 10)) 0.7 else 1.0
          val adjustedQuality = quality.map(q => q * depthModifier * regionModifier)

          val row = Seq(
            escapeCsv(locus.contig),
            locus.position.toString,
            escapeCsv(locus.name),
            escapeCsv(haplogroup),
            escapeCsv(locus.ref),
            escapeCsv(locus.alt),
            escapeCsv(called),
            state,
            depth.map(_.toString).getOrElse(""),
            quality.map(q => f"$q%.1f").getOrElse(""),
            regionType,
            escapeCsv(regionName),
            f"$regionModifier%.2f",
            adjustedQuality.map(q => f"$q%.1f").getOrElse("")
          ).mkString(",")

          writer.println(row)
        }
      }
      createdFiles += snpDetailsFile
    }

    // Novel SNPs CSV
    privateVariants.foreach { variants =>
      val snpVariants = variants.filter(v => v.ref.length == 1 && v.alt.length == 1)

      if (snpVariants.nonEmpty) {
        val novelSnpsFile = new File(outputDir, s"${prefix}_novel_snps.csv")
        Using.resource(new PrintWriter(novelSnpsFile)) { writer =>
          writer.println("contig,position,ref,alt,depth,quality,region_type,region_name,quality_modifier,adjusted_quality")

          snpVariants.sortBy(v => (v.contig, v.position)).foreach { v =>
            val regionAnn = yRegionAnnotator.map(_.annotate(v.contig, v.position, v.depth, expectedYCoverage))
            val regionType = regionAnn.map(_.primaryRegion.map(_.regionType.toString).getOrElse("")).getOrElse("")
            val regionName = regionAnn.map(_.shortDisplay).getOrElse("")
            val regionModifier = regionAnn.map(_.qualityModifier).getOrElse(1.0)

            val depthModifier = if (v.depth.exists(_ < 10)) 0.7 else 1.0
            val adjustedQuality = v.quality.map(q => q * depthModifier * regionModifier)

            val row = Seq(
              escapeCsv(v.contig),
              v.position.toString,
              escapeCsv(v.ref),
              escapeCsv(v.alt),
              v.depth.map(_.toString).getOrElse(""),
              v.quality.map(q => f"$q%.1f").getOrElse(""),
              regionType,
              escapeCsv(regionName),
              f"$regionModifier%.2f",
              adjustedQuality.map(q => f"$q%.1f").getOrElse("")
            ).mkString(",")

            writer.println(row)
          }
        }
        createdFiles += novelSnpsFile
      }
    }

    // Novel Indels CSV
    privateVariants.foreach { variants =>
      val indelVariants = variants.filterNot(v => v.ref.length == 1 && v.alt.length == 1)

      if (indelVariants.nonEmpty) {
        val novelIndelsFile = new File(outputDir, s"${prefix}_novel_indels.csv")
        Using.resource(new PrintWriter(novelIndelsFile)) { writer =>
          writer.println("contig,position,ref,alt,depth,quality,region_type,region_name,quality_modifier,adjusted_quality,str_marker,str_period,str_repeats")

          indelVariants.sortBy(v => (v.contig, v.position)).foreach { v =>
            val regionAnn = yRegionAnnotator.map(_.annotate(v.contig, v.position, v.depth, expectedYCoverage))
            val regionType = regionAnn.map(_.primaryRegion.map(_.regionType.toString).getOrElse("")).getOrElse("")
            val regionName = regionAnn.map(_.shortDisplay).getOrElse("")
            val regionModifier = regionAnn.map(_.qualityModifier).getOrElse(1.0)

            val depthModifier = if (v.depth.exists(_ < 10)) 0.7 else 1.0
            val adjustedQuality = v.quality.map(q => q * depthModifier * regionModifier)

            // STR annotation if available
            val strInfo = strAnnotator.flatMap(_.findOverlapping(v.contig, v.position))
            val strMarker = strInfo.flatMap(_.name).getOrElse("")
            val strPeriod = strInfo.map(_.period.toString).getOrElse("")
            val strRepeats = strInfo.map(r => f"${r.numRepeats}%.1f").getOrElse("")

            val row = Seq(
              escapeCsv(v.contig),
              v.position.toString,
              escapeCsv(v.ref),
              escapeCsv(v.alt),
              v.depth.map(_.toString).getOrElse(""),
              v.quality.map(q => f"$q%.1f").getOrElse(""),
              regionType,
              escapeCsv(regionName),
              f"$regionModifier%.2f",
              adjustedQuality.map(q => f"$q%.1f").getOrElse(""),
              escapeCsv(strMarker),
              strPeriod,
              strRepeats
            ).mkString(",")

            writer.println(row)
          }
        }
        createdFiles += novelIndelsFile
      }
    }

    createdFiles.toList
  }
}
