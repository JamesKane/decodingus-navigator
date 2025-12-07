package com.decodingus.workspace.model

/**
 * A single variant call representing a mutation.
 * Part of the Atmosphere Lexicon (com.decodingus.atmosphere.defs#variantCall).
 *
 * @param contigAccession  GenBank accession for the contig (e.g., 'NC_000024.10' for chrY)
 * @param position         1-based position on the contig
 * @param referenceAllele  Reference allele
 * @param alternateAllele  Alternate (mutant) allele
 * @param rsId             dbSNP rsID if known (e.g., 'rs123456')
 * @param variantName      Common name if known (e.g., 'M269', 'L21')
 * @param genotype         Called genotype (e.g., 'A', 'T', 'het')
 * @param quality          Variant call quality score
 * @param depth            Read depth at this position
 */
case class VariantCall(
  contigAccession: String,
  position: Int,
  referenceAllele: String,
  alternateAllele: String,
  rsId: Option[String] = None,
  variantName: Option[String] = None,
  genotype: Option[String] = None,
  quality: Option[Double] = None,
  depth: Option[Int] = None
)

/**
 * Detailed private variant calls that extend beyond the terminal haplogroup.
 * Part of the Atmosphere Lexicon (com.decodingus.atmosphere.defs#privateVariantData).
 *
 * @param variants        List of private (novel) variant calls
 * @param analysisVersion Version of the haplogroup analysis pipeline that identified these variants
 * @param referenceTree   The haplogroup tree version used (e.g., 'ISOGG-2024', 'PhyloTree-17')
 */
case class PrivateVariantData(
  variants: List[VariantCall] = List.empty,
  analysisVersion: Option[String] = None,
  referenceTree: Option[String] = None
)
