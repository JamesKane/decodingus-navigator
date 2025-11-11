package com.decodingus.analysis

import com.decodingus.haplogroup.caller.GatkHaplotypeCallerProcessor
import com.decodingus.haplogroup.model.Locus
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.variantcontext.VariantContext

import java.io.{File, PrintWriter}
import scala.jdk.CollectionConverters._

class PrivateSnpProcessor {

  def findPrivateSnps(
    bamPath: String,
    referencePath: String,
    contig: String,
    knownLoci: Set[Locus],
    onProgress: (String, Double, Double) => Unit
  ): List[VariantContext] = {
    val caller = new GatkHaplotypeCallerProcessor()
    val allVariantsVcf = caller.callAllVariantsInContig(bamPath, referencePath, contig, onProgress)

    onProgress("Filtering for private SNPs...", 0.8, 1.0)
    val reader = new VCFFileReader(allVariantsVcf, false)
    val knownPositions = knownLoci.map(_.position)

    val privateVariants = reader.iterator().asScala.filterNot { vc =>
      knownPositions.contains(vc.getStart.toLong)
    }.toList

    reader.close()
    onProgress("Private SNP analysis complete.", 1.0, 1.0)
    privateVariants
  }

  def writeReport(privateVariants: List[VariantContext], reportFile: File): Unit = {
    val writer = new PrintWriter(reportFile)
    try {
      writer.println("Private SNPs Report")
      writer.println("=" * 20)
      privateVariants.foreach { vc =>
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
