package com.decodingus.haplogroup.report

import com.decodingus.analysis.PrivateVariant
import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupResult, Locus}
import com.decodingus.haplogroup.tree.{TreeProviderType, TreeType}

import java.io.{File, PrintWriter}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import scala.util.Using

/**
 * Generates Yleaf-style haplogroup reports.
 */
object HaplogroupReportWriter {

  /**
   * Write a haplogroup analysis report to the specified directory.
   *
   * @param outputDir Directory to write the report
   * @param treeType Y-DNA or MT-DNA
   * @param results Scored haplogroup results
   * @param tree The haplogroup tree used for analysis
   * @param snpCalls The SNP calls from the VCF
   * @param sampleName Optional sample name
   * @param privateVariants Optional list of private/novel variants
   * @param treeProvider Optional tree provider used for analysis
   */
  def writeReport(
                   outputDir: File,
                   treeType: TreeType,
                   results: List[HaplogroupResult],
                   tree: List[Haplogroup],
                   snpCalls: Map[Long, String],
                   sampleName: Option[String] = None,
                   privateVariants: Option[List[PrivateVariant]] = None,
                   treeProvider: Option[TreeProviderType] = None
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
        writer.println(f"${"Position"}%12s  ${"SNP Name"}%-20s  ${"Ancestral"}%10s  ${"Derived"}%10s  ${"Called"}%10s  ${"State"}%10s")
        writer.println("-" * 80)

        val allLociOnPath = path.flatMap(_.loci)

        allLociOnPath.sortBy(_.position).foreach { locus =>
          val called = snpCalls.get(locus.position).getOrElse("-")
          val state = snpCalls.get(locus.position) match {
            case Some(base) if base.equalsIgnoreCase(locus.alt) => "Derived"
            case Some(base) if base.equalsIgnoreCase(locus.ref) => "Ancestral"
            case Some(_) => "Unknown"
            case None => "No Call"
          }
          writer.println(f"${locus.position}%12d  ${locus.name}%-20s  ${locus.ref}%10s  ${locus.alt}%10s  $called%10s  $state%10s")
        }
        writer.println()
      }

      // Private/Novel Variants section
      privateVariants.foreach { variants =>
        writer.println("-" * 80)
        writer.println("PRIVATE/NOVEL VARIANTS")
        writer.println("-" * 80)
        if (variants.isEmpty) {
          writer.println("  No private variants discovered.")
        } else {
          writer.println(s"  Total private variants: ${variants.size}")
          writer.println()
          writer.println(f"${"Position"}%12s  ${"Ref"}%6s  ${"Alt"}%6s  ${"Quality"}%10s")
          writer.println("-" * 80)
          variants.sortBy(_.position).foreach { v =>
            val qualStr = v.quality.map(q => f"$q%.1f").getOrElse("-")
            writer.println(f"${v.position}%12d  ${v.ref}%6s  ${v.alt}%6s  $qualStr%10s")
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
        writer.println(s"  Private variants discovered: ${variants.size}")
      }
      writer.println()
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
}
