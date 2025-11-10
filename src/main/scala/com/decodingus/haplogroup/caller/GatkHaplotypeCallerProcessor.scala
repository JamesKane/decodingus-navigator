package com.decodingus.haplogroup.caller

import com.decodingus.haplogroup.model.Locus
import htsjdk.variant.vcf.VCFFileReader
import org.broadinstitute.hellbender.Main

import java.io.{File, PrintWriter}
import scala.jdk.CollectionConverters._

class GatkHaplotypeCallerProcessor {

  def callSnps(
    bamPath: String,
    referencePath: String,
    loci: List[Locus],
    onProgress: (String, Double, Double) => Unit
  ): Map[Long, String] = {
    onProgress("Calling SNPs with GATK HaplotypeCaller...", 0.0, 1.0)

    val allelesFile = createVcfAllelesFile(loci)
    val vcfFile = File.createTempFile("haplotypes", ".vcf")
    vcfFile.deleteOnExit()
    allelesFile.deleteOnExit()

    val args = Array(
      "HaplotypeCaller",
      "-I", bamPath,
      "-R", referencePath,
      "-O", vcfFile.getAbsolutePath,
      "-L", allelesFile.getAbsolutePath,
      "--alleles", allelesFile.getAbsolutePath,
      "--genotyping-mode", "GENOTYPE_GIVEN_ALLELES"
    )
    Main.main(args)

    onProgress("Parsing VCF output...", 0.9, 1.0)
    val snpCalls = parseVcf(vcfFile)
    onProgress("SNP calling complete.", 1.0, 1.0)
    snpCalls
  }

  private def createVcfAllelesFile(loci: List[Locus]): File = {
    val file = File.createTempFile("alleles", ".vcf")
    val writer = new PrintWriter(file)
    try {
      writer.println("##fileformat=VCFv4.2")
      writer.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
      loci.foreach {
        locus =>
          locus.coordinates.get("GRCh38").foreach { coord => // Assuming GRCh38 for now
            writer.println(s"${coord.chromosome}\t${coord.position}\t${locus.name}\t${coord.ancestral}\t${coord.derived}\t.\t.\t.")
          }
      }
    } finally {
      writer.close()
    }
    file
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