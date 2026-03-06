package com.decodingus.analysis.sv

import htsjdk.variant.variantcontext.{Allele, GenotypeBuilder, VariantContextBuilder}
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriterBuilder}
import htsjdk.variant.vcf.*

import java.io.File
import java.nio.file.Path
import java.time.LocalDate
import java.time.format.DateTimeFormatter
import scala.jdk.CollectionConverters.*

/**
 * Writes structural variant calls to VCF format.
 *
 * Follows VCF 4.3 specification for structural variants:
 * - SVTYPE, SVLEN, END in INFO field
 * - CIPOS, CIEND for confidence intervals
 * - PE, SR for paired-end and split-read support
 * - RD for relative depth (CNVs)
 *
 * References:
 * - VCF 4.3 Specification: https://samtools.github.io/hts-specs/VCFv4.3.pdf
 * - SV representation in VCF: https://samtools.github.io/hts-specs/VCFv4.3.pdf#page=12
 */
class SvVcfWriter(config: SvCallerConfig = SvCallerConfig.default) {

  /**
   * Write SV calls to a gzipped VCF file with tabix index.
   *
   * @param calls          List of SV calls to write
   * @param outputPath     Output path (will create .vcf.gz and .vcf.gz.tbi)
   * @param sampleName     Sample name for the VCF
   * @param referenceBuild Reference build name
   */
  def write(
    calls: List[SvCall],
    outputPath: Path,
    sampleName: String,
    referenceBuild: String
  ): Unit = {
    val header = createHeader(sampleName, referenceBuild)

    // Create VCF writer with block compression
    val writer = new VariantContextWriterBuilder()
      .setOutputFile(outputPath.toFile)
      .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
      .setOption(Options.INDEX_ON_THE_FLY)
      .build()

    writer.writeHeader(header)

    // Sort calls by position and write
    val sortedCalls = calls.sortBy(c => (c.chrom, c.start))
    sortedCalls.foreach { call =>
      val vc = createVariantContext(call, sampleName, header)
      writer.add(vc)
    }

    writer.close()
  }

  /**
   * Create VCF header with SV-specific fields.
   */
  private def createHeader(sampleName: String, referenceBuild: String): VCFHeader = {
    val headerLines = new java.util.LinkedHashSet[VCFHeaderLine]()

    // File format
    headerLines.add(new VCFHeaderLine("fileformat", "VCFv4.3"))
    headerLines.add(new VCFHeaderLine("fileDate", LocalDate.now().format(DateTimeFormatter.BASIC_ISO_DATE)))
    headerLines.add(new VCFHeaderLine("source", "DUNavigator_SvCaller"))
    headerLines.add(new VCFHeaderLine("reference", referenceBuild))

    // INFO fields
    headerLines.add(new VCFInfoHeaderLine(
      "SVTYPE", 1, VCFHeaderLineType.String,
      "Type of structural variant"
    ))
    headerLines.add(new VCFInfoHeaderLine(
      "SVLEN", 1, VCFHeaderLineType.Integer,
      "Difference in length between REF and ALT alleles"
    ))
    headerLines.add(new VCFInfoHeaderLine(
      "END", 1, VCFHeaderLineType.Integer,
      "End position of the variant"
    ))
    headerLines.add(new VCFInfoHeaderLine(
      "CIPOS", 2, VCFHeaderLineType.Integer,
      "Confidence interval around POS"
    ))
    headerLines.add(new VCFInfoHeaderLine(
      "CIEND", 2, VCFHeaderLineType.Integer,
      "Confidence interval around END"
    ))
    headerLines.add(new VCFInfoHeaderLine(
      "PE", 1, VCFHeaderLineType.Integer,
      "Number of paired-end reads supporting the variant"
    ))
    headerLines.add(new VCFInfoHeaderLine(
      "SR", 1, VCFHeaderLineType.Integer,
      "Number of split reads supporting the variant"
    ))
    headerLines.add(new VCFInfoHeaderLine(
      "RD", 1, VCFHeaderLineType.Float,
      "Relative read depth (observed/expected)"
    ))
    headerLines.add(new VCFInfoHeaderLine(
      "MATEID", 1, VCFHeaderLineType.String,
      "ID of mate breakend for translocations"
    ))

    // ALT alleles for SV types
    headerLines.add(new VCFSimpleHeaderLine(
      "ALT", "<DEL>", "Deletion"
    ))
    headerLines.add(new VCFSimpleHeaderLine(
      "ALT", "<DUP>", "Duplication"
    ))
    headerLines.add(new VCFSimpleHeaderLine(
      "ALT", "<INV>", "Inversion"
    ))
    headerLines.add(new VCFSimpleHeaderLine(
      "ALT", "<INS>", "Insertion"
    ))

    // FORMAT fields
    headerLines.add(new VCFFormatHeaderLine(
      "GT", 1, VCFHeaderLineType.String,
      "Genotype"
    ))
    headerLines.add(new VCFFormatHeaderLine(
      "GQ", 1, VCFHeaderLineType.Integer,
      "Genotype Quality"
    ))

    // FILTER fields
    headerLines.add(new VCFFilterHeaderLine(
      "LowQual", "Quality score below threshold"
    ))
    headerLines.add(new VCFFilterHeaderLine(
      "LowSupport", "Insufficient read support"
    ))

    new VCFHeader(headerLines, java.util.Collections.singleton(sampleName))
  }

  /**
   * Create a VariantContext from an SV call.
   */
  private def createVariantContext(
    call: SvCall,
    sampleName: String,
    header: VCFHeader
  ): htsjdk.variant.variantcontext.VariantContext = {
    val refAllele = Allele.create("N", true)

    val altAllele = call.svType match {
      case SvType.BND =>
        // BND format: N]chr:pos] or [chr:pos[N
        val mateChrom = call.mateChrom.getOrElse(call.chrom)
        val matePos = call.matePos.getOrElse(call.end)
        Allele.create(s"N]$mateChrom:$matePos]", false)
      case SvType.DEL => Allele.create("<DEL>", false)
      case SvType.DUP => Allele.create("<DUP>", false)
      case SvType.INV => Allele.create("<INV>", false)
      case SvType.INS => Allele.create("<INS>", false)
    }

    val alleles = java.util.Arrays.asList(refAllele, altAllele)

    val builder = new VariantContextBuilder()
      .chr(call.chrom)
      .start(call.start)
      .stop(if (call.svType == SvType.BND) call.start else call.end)
      .id(call.id)
      .alleles(alleles)
      .log10PError(-call.quality / 10.0) // Convert Phred to log10

    // Add INFO fields
    builder.attribute("SVTYPE", call.svType.toString)
    builder.attribute("SVLEN", call.svLen.toInt)
    builder.attribute("END", call.end.toInt)
    builder.attribute("CIPOS", java.util.Arrays.asList(call.ciPos._1: Integer, call.ciPos._2: Integer))
    builder.attribute("CIEND", java.util.Arrays.asList(call.ciEnd._1: Integer, call.ciEnd._2: Integer))

    if (call.pairedEndSupport > 0) {
      builder.attribute("PE", call.pairedEndSupport)
    }
    if (call.splitReadSupport > 0) {
      builder.attribute("SR", call.splitReadSupport)
    }
    call.relativeDepth.foreach { rd =>
      builder.attribute("RD", rd.toFloat)
    }

    // Set filter
    if (call.filter == "PASS") {
      builder.passFilters()
    } else {
      builder.filter(call.filter)
    }

    // Add genotype
    val genotypeAlleles = call.genotype match {
      case "0/1" => java.util.Arrays.asList(refAllele, altAllele)
      case "1/1" => java.util.Arrays.asList(altAllele, altAllele)
      case _ => java.util.Arrays.asList(refAllele, altAllele)
    }

    val genotype = new GenotypeBuilder(sampleName)
      .alleles(genotypeAlleles)
      .GQ(call.quality.toInt)
      .make()

    builder.genotypes(genotype)

    builder.make()
  }
}
