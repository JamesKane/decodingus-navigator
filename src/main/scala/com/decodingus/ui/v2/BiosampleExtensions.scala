package com.decodingus.ui.v2

import com.decodingus.workspace.model.{Biosample, HaplogroupResult, Project}

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

    /** Get the full Y-DNA HaplogroupResult if present */
    def yHaplogroupResult: Option[HaplogroupResult] =
      b.haplogroups.flatMap(_.yDna)

    /** Get the full mtDNA HaplogroupResult if present */
    def mtHaplogroupResult: Option[HaplogroupResult] =
      b.haplogroups.flatMap(_.mtDna)

    /** Get Y-DNA haplogroup score if present */
    def yHaplogroupScore: Option[Double] =
      b.haplogroups.flatMap(_.yDna).map(_.score)

    /** Get mtDNA haplogroup score if present */
    def mtHaplogroupScore: Option[Double] =
      b.haplogroups.flatMap(_.mtDna).map(_.score)

    /** Get Y-DNA matching SNP count if present */
    def yMatchingSnps: Option[Int] =
      b.haplogroups.flatMap(_.yDna).flatMap(_.matchingSnps)

    /** Get Y-DNA mismatching SNP count if present */
    def yMismatchingSnps: Option[Int] =
      b.haplogroups.flatMap(_.yDna).flatMap(_.mismatchingSnps)

    /** Get Y-DNA ancestral match count if present */
    def yAncestralMatches: Option[Int] =
      b.haplogroups.flatMap(_.yDna).flatMap(_.ancestralMatches)

    /** Get Y-DNA lineage path if present */
    def yLineagePath: Option[List[String]] =
      b.haplogroups.flatMap(_.yDna).flatMap(_.lineagePath)

    /** Get mtDNA lineage path if present */
    def mtLineagePath: Option[List[String]] =
      b.haplogroups.flatMap(_.mtDna).flatMap(_.lineagePath)

    /** Get atUri for identification (use sampleAccession as fallback ID) */
    def id: String = b.atUri.getOrElse(b.sampleAccession)

    /** Check if biosample has any sequencing data references */
    def hasSequenceData: Boolean = b.sequenceRunRefs.nonEmpty

    /** Check if biosample has any genotype data references */
    def hasGenotypeData: Boolean = b.genotypeRefs.nonEmpty

    /** Count of sequencing run references */
    def sequenceRunCount: Int = b.sequenceRunRefs.size

    /** Count of genotype (chip) profile references */
    def genotypeCount: Int = b.genotypeRefs.size

    /** Count of STR profile references */
    def strProfileCount: Int = b.strProfileRefs.size
  }

  // ============================================================================
  // HaplogroupResult Extensions
  // ============================================================================

  extension (h: HaplogroupResult) {
    /** Format the lineage path as a displayable string */
    def formattedPath: String =
      h.lineagePath.map(_.mkString(" → ")).getOrElse("")

    /** Confidence level as a human-readable string */
    def confidenceLevel: String = h.score match {
      case s if s >= 0.95 => "HIGH"
      case s if s >= 0.80 => "MEDIUM"
      case _ => "LOW"
    }

    /** Confidence as a percentage string */
    def confidencePercent: String = f"${h.score * 100}%.1f%%"

    /** Derived SNP count (matching SNPs) */
    def derivedCount: Int = h.matchingSnps.getOrElse(0)

    /** Ancestral count */
    def ancestralCount: Int = h.ancestralMatches.getOrElse(0)

    /** Source display name */
    def sourceDisplay: String = h.source match {
      case Some("wgs") => "WGS"
      case Some("bigy") => "Big Y"
      case Some("chip") => "SNP Array"
      case Some(other) => other.toUpperCase
      case None => "Unknown"
    }

    /** Quality rating based on source and score */
    def qualityRating: String = {
      val sourceQuality = h.source match {
        case Some("wgs") => 2.0
        case Some("bigy") => 1.5
        case Some("chip") => 1.0
        case _ => 0.5
      }
      val combined = (h.score + sourceQuality) / 2.0
      combined match {
        case q if q >= 1.4 => "Excellent"
        case q if q >= 1.1 => "Good"
        case q if q >= 0.8 => "Fair"
        case _ => "Poor"
      }
    }
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
