package com.decodingus.ui.v2

import com.decodingus.workspace.model.{Biosample, Project}

import java.time.LocalDateTime

/**
 * Extension methods for Biosample and Project to provide convenient accessors.
 * This simplifies UI code by providing shorter names for commonly accessed fields.
 */
object BiosampleExtensions {

  // ============================================================================
  // Biosample Extensions
  // ============================================================================

  extension (b: Biosample) {
    /** Alias for sampleAccession */
    def accession: String = b.sampleAccession

    /** Alias for donorIdentifier, as Option for UI compatibility */
    def donorId: Option[String] = Option(b.donorIdentifier).filter(_.nonEmpty)

    /** Alias for centerName */
    def center: Option[String] = b.centerName

    /** Extract Y-DNA haplogroup name if present */
    def yHaplogroup: Option[String] =
      b.haplogroups.flatMap(_.yDna).map(_.haplogroupName)

    /** Extract mtDNA haplogroup name if present */
    def mtHaplogroup: Option[String] =
      b.haplogroups.flatMap(_.mtDna).map(_.haplogroupName)

    /** Get Y-DNA haplogroup score if present */
    def yHaplogroupScore: Option[Double] =
      b.haplogroups.flatMap(_.yDna).map(_.score)

    /** Get mtDNA haplogroup score if present */
    def mtHaplogroupScore: Option[Double] =
      b.haplogroups.flatMap(_.mtDna).map(_.score)

    /** Get Y-DNA matching SNP count if present */
    def yMatchingSnps: Option[Int] =
      b.haplogroups.flatMap(_.yDna).flatMap(_.matchingSnps)

    /** Get Y-DNA ancestral match count if present */
    def yAncestralMatches: Option[Int] =
      b.haplogroups.flatMap(_.yDna).flatMap(_.ancestralMatches)

    /** Get atUri for identification (use sampleAccession as fallback ID) */
    def id: String = b.atUri.getOrElse(b.sampleAccession)

    /** Check if biosample has any sequencing data references */
    def hasSequenceData: Boolean = b.sequenceRunRefs.nonEmpty

    /** Check if biosample has any genotype data references */
    def hasGenotypeData: Boolean = b.genotypeRefs.nonEmpty
  }

  // ============================================================================
  // Project Extensions
  // ============================================================================

  extension (p: Project) {
    /** Alias for projectName */
    def name: String = p.projectName

    /** Alias for memberRefs (member AT URIs or accessions) */
    def memberAccessions: List[String] = p.memberRefs

    /** Get created timestamp from meta */
    def createdAt: Option[LocalDateTime] = Some(p.meta.createdAt)

    /** Get atUri for identification (use projectName as fallback) */
    def projectId: String = p.atUri.getOrElse(p.projectName)
  }
}
