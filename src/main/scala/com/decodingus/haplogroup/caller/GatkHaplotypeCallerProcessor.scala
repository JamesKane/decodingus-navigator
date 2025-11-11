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
}