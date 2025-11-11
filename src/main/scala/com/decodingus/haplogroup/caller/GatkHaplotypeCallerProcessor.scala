package com.decodingus.haplogroup.caller

import org.broadinstitute.hellbender.Main
import java.io.File

class GatkHaplotypeCallerProcessor {

  def callSnps(
    bamPath: String,
    referencePath: String,
    allelesVcf: File,
    onProgress: (String, Double, Double) => Unit
  ): File = {
    onProgress("Calling SNPs with GATK HaplotypeCaller...", 0.0, 1.0)

    val vcfFile = File.createTempFile("haplotypes", ".vcf")
    vcfFile.deleteOnExit()

    val args = Array(
      "HaplotypeCaller",
      "-I", bamPath,
      "-R", referencePath,
      "-O", vcfFile.getAbsolutePath,
      "-L", allelesVcf.getAbsolutePath,
      "--alleles", allelesVcf.getAbsolutePath,
      "--genotyping-mode", "GENOTYPE_GIVEN_ALLELES"
    )
    Main.main(args)

    onProgress("SNP calling complete.", 1.0, 1.0)
    vcfFile
  }

  def callAllVariantsInContig(
    bamPath: String,
    referencePath: String,
    contig: String,
    onProgress: (String, Double, Double) => Unit
  ): File = {
    onProgress(s"Calling all variants in $contig with GATK HaplotypeCaller...", 0.0, 1.0)

    val vcfFile = File.createTempFile(s"variants-$contig", ".vcf")
    vcfFile.deleteOnExit()

    val args = Array(
      "HaplotypeCaller",
      "-I", bamPath,
      "-R", referencePath,
      "-O", vcfFile.getAbsolutePath,
      "-L", contig
    )
    Main.main(args)

    onProgress(s"Variant calling for $contig complete.", 1.0, 1.0)
    vcfFile
  }
}
