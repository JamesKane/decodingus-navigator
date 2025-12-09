package com.decodingus.analysis

import com.decodingus.haplogroup.caller.{GatkHaplotypeCallerProcessor, TwoPassCallerResult}
import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupResult, Locus}
import com.decodingus.haplogroup.report.HaplogroupReportWriter
import com.decodingus.haplogroup.scoring.HaplogroupScorer
import com.decodingus.haplogroup.tree.{TreeCache, TreeProvider, TreeProviderType, TreeType}
import com.decodingus.haplogroup.vendor.{DecodingUsTreeProvider, FtdnaTreeProvider}
import com.decodingus.liftover.LiftoverProcessor
import com.decodingus.model.LibraryStats
import com.decodingus.refgenome.{LiftoverGateway, MultiContigReferenceQuerier, ReferenceGateway, ReferenceQuerier}
import htsjdk.variant.vcf.VCFFileReader

import java.io.{File, PrintWriter}
import java.nio.file.{Files, Path}
import scala.jdk.CollectionConverters.*
import scala.util.Using

/**
 * A private/novel variant not found in the haplogroup tree.
 */
case class PrivateVariant(
  contig: String,
  position: Long,
  ref: String,
  alt: String,
  quality: Option[Double]
)

class HaplogroupProcessor {

  private val standardContigOrder: Map[String, Int] = (1 to 22).map(i => s"chr$i" -> i).toMap ++
    Map("chrX" -> 23, "chrY" -> 24, "chrM" -> 25)

  private val ARTIFACT_SUBDIR_NAME = "haplogroup"

  /**
   * Analyze a BAM/CRAM file for haplogroup assignment.
   *
   * @param bamPath Path to the BAM/CRAM file
   * @param libraryStats Library statistics from initial analysis
   * @param treeType Y-DNA or MT-DNA tree type
   * @param treeProviderType Tree data provider (FTDNA or DecodingUs)
   * @param onProgress Progress callback
   * @param artifactContext Optional context for organizing output artifacts by subject/run/alignment
   */
  def analyze(
               bamPath: String,
               libraryStats: LibraryStats,
               treeType: TreeType,
               treeProviderType: TreeProviderType,
               onProgress: (String, Double, Double) => Unit,
               artifactContext: Option[ArtifactContext] = None
             ): Either[String, List[HaplogroupResult]] = {

    onProgress("Loading haplogroup tree...", 0.0, 1.0)
    val treeProvider: TreeProvider = treeProviderType match {
      case TreeProviderType.FTDNA => new FtdnaTreeProvider(treeType)
      case TreeProviderType.DECODINGUS => new DecodingUsTreeProvider(treeType)
    }
    val treeCache = new TreeCache()

    treeProvider.loadTree(libraryStats.referenceBuild).flatMap { tree =>
      val allLoci = collectAllLoci(tree).distinct
      // The 'contig' variable is no longer directly used here for createVcfAllelesFile,
      // as Locus objects now carry their own contig information.
      // However, it's still used for liftover and other contig-specific operations.
      val primaryContig = if (treeType == TreeType.YDNA) "chrY" else "chrM"
      val outputPrefix = if (treeType == TreeType.YDNA) "ydna" else "mtdna"

      val referenceBuild = libraryStats.referenceBuild
      val treeSourceBuild = if (treeProvider.supportedBuilds.contains(referenceBuild)) {
        referenceBuild
      } else {
        treeProvider.sourceBuild
      }

      val referenceGateway = new ReferenceGateway((_, _) => {})

      referenceGateway.resolve(treeSourceBuild).flatMap { treeRefPath =>
        // Check if cached sites VCF exists and is valid
        val cachedSitesVcf = treeCache.getSitesVcfPath(treeProvider.cachePrefix, treeSourceBuild)
        val initialAllelesVcf = if (treeCache.isSitesVcfValid(treeProvider.cachePrefix, treeSourceBuild)) {
          onProgress("Using cached sites VCF...", 0.05, 1.0)
          cachedSitesVcf
        } else {
          onProgress("Creating sites VCF...", 0.05, 1.0)
          createVcfAllelesFile(allLoci, treeRefPath.toString, treeType, Some(cachedSitesVcf))
        }

        val (allelesForCalling, performReverseLiftover) = if (referenceBuild == treeSourceBuild) {
          onProgress("Reference builds match.", 0.1, 1.0)
          (Right(initialAllelesVcf), false)
        } else {
          onProgress(s"Reference mismatch: tree is $treeSourceBuild, BAM/CRAM is $referenceBuild. Performing liftover...", 0.1, 1.0)
          // Note: The contig parameter here for liftover still refers to the primary contig for the tree type.
          val lifted = performLiftover(initialAllelesVcf, primaryContig, treeSourceBuild, referenceBuild, onProgress)
          (lifted, true)
        }

        allelesForCalling.flatMap { vcf =>
          referenceGateway.resolve(referenceBuild).flatMap { referencePath =>
            val caller = new GatkHaplotypeCallerProcessor()
            // Pass artifact directory for called VCF and logs
            val artifactDir = artifactContext.map(_.getSubdir(ARTIFACT_SUBDIR_NAME))

            // Two-pass calling: tree sites first, then private variants
            caller.callTwoPass(
              bamPath,
              referencePath.toString,
              vcf,
              (msg, done, total) => onProgress(msg, 0.2 + (done * 0.5), 1.0),
              artifactDir,
              Some(outputPrefix)
            ).flatMap { twoPassResult =>
              // Handle reverse liftover for tree sites VCF if needed
              val finalTreeVcf = if (performReverseLiftover) {
                onProgress("Performing reverse liftover on tree sites...", 0.72, 1.0)
                performLiftover(twoPassResult.treeSitesVcf, primaryContig, referenceBuild, treeSourceBuild, onProgress)
              } else {
                Right(twoPassResult.treeSitesVcf)
              }

              finalTreeVcf.flatMap { scoredVcf =>
                onProgress("Scoring haplogroups...", 0.8, 1.0)
                val snpCalls = parseVcf(scoredVcf)
                val scorer = new HaplogroupScorer()
                val results = scorer.score(tree, snpCalls)

                // Identify private variants - only exclude positions on path to terminal haplogroup
                // Positions on other branches could be legitimate private variants for undiscovered sub-clades
                onProgress("Identifying private variants...", 0.85, 1.0)
                val terminalHaplogroup = results.headOption.map(_.name).getOrElse("")
                val pathPositions = collectPathPositions(tree, terminalHaplogroup)
                val privateVariants = parsePrivateVariants(twoPassResult.privateVariantsVcf, pathPositions)

                // Write report to artifact directory if available
                artifactDir.foreach { dir =>
                  onProgress("Writing haplogroup report...", 0.9, 1.0)
                  HaplogroupReportWriter.writeReport(
                    outputDir = dir.toFile,
                    treeType = treeType,
                    results = results,
                    tree = tree,
                    snpCalls = snpCalls,
                    sampleName = None,
                    privateVariants = Some(privateVariants)
                  )
                }

                onProgress("Analysis complete.", 1.0, 1.0)
                Right(results)
              }
            }
          }
        }
      }
    }
  }

  /**
   * Parse private variants from VCF, excluding known tree positions.
   */
  private def parsePrivateVariants(vcfFile: File, treePositions: Set[Long]): List[PrivateVariant] = {
    val reader = new VCFFileReader(vcfFile, false)
    val variants = reader.iterator().asScala.flatMap { vc =>
      val pos = vc.getStart.toLong
      if (!treePositions.contains(pos)) {
        val genotype = vc.getGenotypes.get(0)
        val allele = genotype.getAlleles.get(0).getBaseString
        val ref = vc.getReference.getBaseString
        val qual = if (vc.hasLog10PError) Some(vc.getPhredScaledQual) else None
        Some(PrivateVariant(
          contig = vc.getContig,
          position = pos,
          ref = ref,
          alt = allele,
          quality = qual
        ))
      } else {
        None
      }
    }.toList
    reader.close()
    variants
  }

  private def performLiftover(
                               vcfFile: File,
                               contig: String, // This contig parameter is still relevant for liftover operations
                               fromBuild: String,
                               toBuild: String,
                               onProgress: (String, Double, Double) => Unit
                             ): Either[String, File] = {
    val liftoverGateway = new LiftoverGateway((_, _) => {})
    val referenceGateway = new ReferenceGateway((_, _) => {})

    for {
      chainFile <- liftoverGateway.resolve(fromBuild, toBuild)
      targetRef <- referenceGateway.resolve(toBuild)
      liftedVcf <- new LiftoverProcessor().liftoverVcf(vcfFile, chainFile, targetRef, (msg, done, total) => onProgress(msg, 0.2 + (done * 0.2), 1.0))
    } yield liftedVcf
  }

  private def createVcfAllelesFile(
    loci: List[Locus],
    referencePath: String,
    treeType: TreeType,
    outputFile: Option[File]
  ): File = {
    // Use provided output file or create temp file
    val vcfFile = outputFile match {
      case Some(file) =>
        // Ensure parent directory exists
        Option(file.getParentFile).foreach(_.mkdirs())
        file
      case None =>
        val tempFile = File.createTempFile("alleles", ".vcf")
        tempFile.deleteOnExit()
        tempFile
    }

    // Group loci by contig first, then by position within each contig
    val lociByContig = loci.groupBy(_.contig)
    val sortedContigs = lociByContig.keys.toList.sortBy(c => standardContigOrder.getOrElse(c, 999))

    Using.resource(new PrintWriter(vcfFile)) { writer =>
      writer.println("##fileformat=VCFv4.2")
      sortedContigs.foreach { c =>
        writer.println(s"##contig=<ID=$c>")
      }
      writer.println("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">")
      writer.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

      // Process one contig at a time - load reference once per contig
      sortedContigs.foreach { contig =>
        Using.resource(new ReferenceQuerier(referencePath, contig)) { refQuerier =>
          val contigLoci = lociByContig(contig)
          // Group by position to combine alternates at the same site
          val groupedByPosition = contigLoci.groupBy(_.position)
          // Filter out positions beyond contig bounds, then sort
          val sortedPositions = groupedByPosition.keys.toList
            .filter(refQuerier.isValidPosition)
            .sorted

          sortedPositions.foreach { position =>
            val lociAtPosition = groupedByPosition(position)
            // Safe to use get since we filtered valid positions above
            val refBase = refQuerier.getBase(position).get

            // Filter to valid SNPs only:
            // - Single base ref and alt (no indels)
            // - Only valid nucleotides A, C, G, T (no dashes, dots, or other characters)
            val validBases = Set("A", "C", "G", "T")
            val snpLoci = lociAtPosition.filter { l =>
              l.ref.length == 1 && l.alt.length == 1 &&
                validBases.contains(l.ref.toUpperCase) &&
                validBases.contains(l.alt.toUpperCase)
            }

            // Collect all valid alternates at this position
            val refAndAlts = snpLoci.flatMap { locus =>
              if (refBase.toString.equalsIgnoreCase(locus.ref)) {
                Some((locus.ref.toUpperCase, locus.alt.toUpperCase))
              } else if (refBase.toString.equalsIgnoreCase(locus.alt)) {
                Some((locus.alt.toUpperCase, locus.ref.toUpperCase))
              } else {
                // This locus is problematic, the reference doesn't match ANC or DER
                None
              }
            }

            if (refAndAlts.nonEmpty) {
              // All refs should be the same (the actual reference base)
              val ref = refAndAlts.head._1
              // Collect unique alternates
              val alts = refAndAlts.map(_._2).distinct.filterNot(_ == ref)

              if (alts.nonEmpty) {
                writer.println(s"$contig\t$position\t.\t$ref\t${alts.mkString(",")}\t.\t.\t.")
              }
            }
          }
        }
      }
    }
    vcfFile
  }

  /**
   * Recursively collect all loci from the haplogroup tree.
   */
  private def collectAllLoci(tree: List[Haplogroup]): List[Locus] = {
    tree.flatMap(collectLociFromHaplogroup)
  }

  private def collectLociFromHaplogroup(haplogroup: Haplogroup): List[Locus] = {
    haplogroup.loci ++ haplogroup.children.flatMap(collectLociFromHaplogroup)
  }

  /**
   * Collect all positions along the path from root to the specified terminal haplogroup.
   * Only these positions should be excluded from private variant detection.
   */
  private def collectPathPositions(tree: List[Haplogroup], terminalName: String): Set[Long] = {
    def findPath(haplogroup: Haplogroup): Option[List[Haplogroup]] = {
      if (haplogroup.name == terminalName) {
        Some(List(haplogroup))
      } else {
        haplogroup.children.flatMap(findPath).headOption.map(path => haplogroup :: path)
      }
    }

    val path = tree.flatMap(findPath).headOption.getOrElse(List.empty)
    path.flatMap(_.loci.map(_.position)).toSet
  }

  private def parseVcf(vcfFile: File): Map[Long, String] = {
    val reader = new VCFFileReader(vcfFile, false)
    val snpCalls = reader.iterator().asScala.map {
      vc =>
        val pos = vc.getStart.toLong
        val genotype = vc.getGenotypes.get(0) // Assuming single sample VCF
        val allele = genotype.getAlleles.get(0).getBaseString
        pos -> allele
    }.toMap
    reader.close()
    snpCalls
  }
}