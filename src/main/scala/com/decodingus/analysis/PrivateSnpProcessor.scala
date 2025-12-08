package com.decodingus.analysis

import com.decodingus.haplogroup.caller.GatkHaplotypeCallerProcessor
import com.decodingus.haplogroup.model.Locus
import com.decodingus.liftover.LiftoverProcessor
import com.decodingus.refgenome.{LiftoverGateway, ReferenceGateway, ReferenceQuerier}
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.variantcontext.VariantContext

import java.io.{File, PrintWriter}
import scala.jdk.CollectionConverters._
import scala.util.Using

class PrivateSnpProcessor {

  def findPrivateSnps(
    bamPath: String,
    referencePath: String, // Path to the BAM/CRAM's reference
    contig: String,
    knownLoci: Set[Locus], // Loci from the tree (in treeSourceBuild)
    treeSourceBuild: String, // The build of the knownLoci
    referenceBuild: String, // The build of the BAM/CRAM
    onProgress: (String, Double, Double) => Unit
  ): List[VariantContext] = {

    val knownPositionsInReferenceBuild: Set[Long] = if (treeSourceBuild != referenceBuild) {
      onProgress(s"Lifting over known loci from $treeSourceBuild to $referenceBuild...", 0.0, 1.0)
      
      val tempKnownLociVcf = createVcfFromLoci(knownLoci.toList, contig, treeSourceBuild)

      val liftoverGateway = new LiftoverGateway((_, _) => {})
      val refGateway = new ReferenceGateway((_, _) => {})

      val liftedKnownLociVcfEither = for {
        chainFile <- liftoverGateway.resolve(treeSourceBuild, referenceBuild)
        targetRef <- refGateway.resolve(referenceBuild)
        liftedVcf <- new LiftoverProcessor().liftoverVcf(tempKnownLociVcf, chainFile, targetRef, (msg, done, total) => onProgress(msg, 0.1 + (done * 0.1), 1.0))
      } yield liftedVcf

      liftedKnownLociVcfEither match {
        case Right(liftedVcf) =>
          val reader = new VCFFileReader(liftedVcf, false)
          val positions = reader.iterator().asScala.map(_.getStart.toLong).toSet
          reader.close()
          positions
        case Left(error) =>
          throw new RuntimeException(s"Failed to liftover known loci: $error")
      }
    } else {
      knownLoci.map(_.position)
    }

    val caller = new GatkHaplotypeCallerProcessor()
    caller.callAllVariantsInContig(bamPath, referencePath, contig, onProgress) match {
      case Right(callerResult) =>
        onProgress("Filtering for private SNPs...", 0.8, 1.0)
        val reader = new VCFFileReader(callerResult.vcfFile, false)

        val privateVariants = reader.iterator().asScala.filterNot {
          vc =>
            knownPositionsInReferenceBuild.contains(vc.getStart.toLong)
        }.toList

        reader.close()
        onProgress("Private SNP analysis complete.", 1.0, 1.0)
        privateVariants
      case Left(error) =>
        throw new RuntimeException(s"Failed to call variants in $contig: $error")
    }
  }

  private def createVcfFromLoci(loci: List[Locus], contig: String, referenceBuild: String): File = {
    val vcfFile = File.createTempFile("known_alleles", ".vcf")
    vcfFile.deleteOnExit()

    val refGateway = new ReferenceGateway((_, _) => {})
    val referencePathEither = refGateway.resolve(referenceBuild)

    referencePathEither match {
      case Right(referencePath) =>
        Using.resource(new PrintWriter(vcfFile)) {
          writer =>
          writer.println("##fileformat=VCFv4.2")
          writer.println(s"##contig=<ID=$contig>")
          writer.println("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">")
          writer.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

          Using.resource(new ReferenceQuerier(referencePath.toString)) {
            refQuerier =>
            loci.foreach {
              locus =>
              val refBase = refQuerier.getBase(contig, locus.position)
              val (ref, alt) = if (refBase.toString.equalsIgnoreCase(locus.ref)) {
                (locus.ref, locus.alt)
              } else if (refBase.toString.equalsIgnoreCase(locus.alt)) {
                (locus.alt, locus.ref)
              } else {
                // This locus is problematic, the reference doesn't match ANC or DER
                // We'll skip it for now
                ("", "")
              }

              if (ref.nonEmpty) {
                writer.println(s"$contig\t${locus.position}\t.\t$ref\t$alt\t.\t.\t.")
              }
            }
          }
        }
        vcfFile
      case Left(error) =>
        throw new RuntimeException(s"Failed to resolve reference for VCF creation: $error")
    }
  }

  def writeReport(privateVariants: List[VariantContext], reportFile: File): Unit = {
    val writer = new PrintWriter(reportFile)
    try {
      writer.println("Private SNPs Report")
      writer.println("=" * 20)
      privateVariants.foreach {
        vc =>
        writer.println(s"Position: ${vc.getContig}:${vc.getStart}")
        writer.println(s"  REF: ${vc.getReference.getBaseString}")
        writer.println(s"  ALT: ${vc.getAlternateAlleles.asScala.map(_.getBaseString).mkString(", ")}")
        writer.println(s"  Genotype: ${vc.getGenotypes.get(0).getGenotypeString}")
        writer.println("-" * 20)
      }
    } finally {
      writer.close()
    }
  }
}