package com.decodingus.analysis

import com.decodingus.haplogroup.caller.GatkHaplotypeCallerProcessor
import com.decodingus.haplogroup.model.{HaplogroupResult, Locus}
import com.decodingus.haplogroup.scoring.HaplogroupScorer
import com.decodingus.haplogroup.tree.{TreeProvider, TreeProviderType, TreeType}
import com.decodingus.haplogroup.vendor.{DecodingUsTreeProvider, FtdnaTreeProvider}
import com.decodingus.liftover.LiftoverProcessor
import com.decodingus.model.LibraryStats
import com.decodingus.refgenome.{LiftoverGateway, ReferenceGateway, ReferenceQuerier}
import htsjdk.variant.vcf.VCFFileReader

import java.io.{File, PrintWriter}
import scala.jdk.CollectionConverters.*
import scala.util.Using

class HaplogroupProcessor {

  private val standardContigOrder: Map[String, Int] = (1 to 22).map(i => s"chr$i" -> i).toMap ++
    Map("chrX" -> 23, "chrY" -> 24, "chrM" -> 25)

  def analyze(
               bamPath: String,
               libraryStats: LibraryStats,
               treeType: TreeType,
               treeProviderType: TreeProviderType,
               onProgress: (String, Double, Double) => Unit
             ): Either[String, List[HaplogroupResult]] = {

    onProgress("Loading haplogroup tree...", 0.0, 1.0)
    val treeProvider: TreeProvider = treeProviderType match {
      case TreeProviderType.FTDNA => new FtdnaTreeProvider(treeType)
      case TreeProviderType.DECODINGUS => new DecodingUsTreeProvider(treeType)
    }

    treeProvider.loadTree(libraryStats.referenceBuild).flatMap { tree =>
      val allLoci = tree.flatMap(h => h.loci).distinct
      // The 'contig' variable is no longer directly used here for createVcfAllelesFile,
      // as Locus objects now carry their own contig information.
      // However, it's still used for liftover and other contig-specific operations.
      val primaryContig = if (treeType == TreeType.YDNA) "chrY" else "chrM"

      val referenceBuild = libraryStats.referenceBuild
      val treeSourceBuild = if (treeProvider.supportedBuilds.contains(referenceBuild)) {
        referenceBuild
      } else {
        treeProvider.sourceBuild
      }

      val referenceGateway = new ReferenceGateway((_, _) => {})

      referenceGateway.resolve(treeSourceBuild).flatMap { treeRefPath =>
        val initialAllelesVcf = createVcfAllelesFile(allLoci, treeRefPath.toString)

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
            caller.callSnps(bamPath, referencePath.toString, vcf, (msg, done, total) => onProgress(msg, 0.4 + (done * 0.4), 1.0)).flatMap { calledVcf =>
              val finalVcf = if (performReverseLiftover) {
                onProgress("Performing reverse liftover on results...", 0.8, 1.0)
                // Note: The contig parameter here for liftover still refers to the primary contig for the tree type.
                performLiftover(calledVcf, primaryContig, referenceBuild, treeSourceBuild, onProgress)
              } else {
                Right(calledVcf)
              }

              finalVcf.flatMap { scoredVcf =>
                onProgress("Scoring haplogroups...", 0.9, 1.0)
                val snpCalls = parseVcf(scoredVcf)
                val scorer = new HaplogroupScorer()
                val results = scorer.score(tree, snpCalls)
                onProgress("Analysis complete.", 1.0, 1.0)
                Right(results)
              }
            }
          }
        }
      }
    }
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

  private def createVcfAllelesFile(loci: List[Locus], referencePath: String): File = {
    val vcfFile = File.createTempFile("alleles", ".vcf")
    vcfFile.deleteOnExit()

    Using.resource(new PrintWriter(vcfFile)) { writer =>
      writer.println("##fileformat=VCFv4.2")
      // Dynamically add contig headers based on unique contigs in loci
      loci.map(_.contig).distinct.foreach { c =>
        writer.println(s"##contig=<ID=$c>")
      }
      writer.println("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">")
      writer.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

      Using.resource(new ReferenceQuerier(referencePath)) { refQuerier =>
        loci.sortBy(locus => (standardContigOrder.getOrElse(locus.contig, 999), locus.position)).foreach {
          case locus =>
            val refBase = refQuerier.getBase(locus.contig, locus.position)
            val (ref, alt) = if (refBase.toString.equalsIgnoreCase(locus.ref)) {
              (locus.ref, locus.alt)
            } else if (refBase.toString.equalsIgnoreCase(locus.alt)) {
              (locus.alt, locus.ref)
            } else {
              // This locus is problematic, the reference doesn't match ANC or DER
              // We'll skip it for now
              ("", "")
            }

            if (ref.nonEmpty && alt.nonEmpty) {
              writer.println(s"${locus.contig}\t${locus.position}\t.\t$ref\t$alt\t.\t.\t.")
            }
        }
      }
    }
    vcfFile
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