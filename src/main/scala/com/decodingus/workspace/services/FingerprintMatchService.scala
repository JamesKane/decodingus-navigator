package com.decodingus.workspace.services

import com.decodingus.model.LibraryStats
import com.decodingus.workspace.model.SequenceRun

/**
 * Result of fingerprint matching against existing sequence runs.
 * Used to detect when the same sequencing run has been aligned to different references.
 */
sealed trait FingerprintMatchResult

object FingerprintMatchResult {
  /**
   * A matching sequence run was found.
   *
   * @param sequenceRun The matching sequence run
   * @param index       The index of the run within the biosample's runs
   * @param confidence  Match confidence level: "HIGH", "MEDIUM", or "LOW"
   */
  case class MatchFound(sequenceRun: SequenceRun, index: Int, confidence: String) extends FingerprintMatchResult

  /**
   * No matching sequence run was found - this is a new run.
   */
  case object NoMatch extends FingerprintMatchResult
}

/**
 * Service for matching sequence runs based on fingerprint data.
 *
 * Fingerprints are computed from BAM/CRAM metadata:
 * - Platform Unit (PU): Unique run identifier from read groups
 * - Library ID (LB): Library preparation identifier
 * - Sample Name (SM): Sample identifier in the BAM
 *
 * This service determines if a newly added alignment file belongs to
 * an existing sequence run (just aligned to a different reference) or
 * represents a completely new sequencing run.
 */
class FingerprintMatchService {

  /**
   * Find an existing sequence run that matches the given fingerprint.
   *
   * Matching is simple: SM (Sample Name) + Platform must match.
   * Same sample name on same platform = same sequencing run aligned to different references.
   *
   * If SM or Platform is unknown, no automatic matching occurs - files are treated as separate runs.
   *
   * @param candidateRuns Sequence runs to search within (typically all runs for a biosample)
   * @param fingerprint   The computed run fingerprint from LibraryStats (kept for API compatibility)
   * @param libraryStats  Full stats for matching criteria
   * @return Match result with confidence level
   */
  def findMatch(
                 candidateRuns: List[(SequenceRun, Int)],
                 fingerprint: String,
                 libraryStats: LibraryStats
               ): FingerprintMatchResult = {
    val sm = libraryStats.sampleName
    val platform = libraryStats.inferredPlatform

    // SM + Platform match (HIGH confidence)
    // Same sample name on same platform = same sequencing run
    if (sm != "Unknown" && platform != "Unknown") {
      candidateRuns.find { case (run, _) =>
        run.sampleName.contains(sm) && run.platformName.equalsIgnoreCase(platform)
      }.map { case (run, idx) =>
        FingerprintMatchResult.MatchFound(run, idx, "HIGH")
      }.getOrElse(FingerprintMatchResult.NoMatch)
    } else {
      // Can't reliably match without SM + Platform - treat as new run
      FingerprintMatchResult.NoMatch
    }
  }

  /**
   * Checks if a fingerprint is likely LOW confidence and should prompt user.
   * LOW confidence typically means we're matching on less reliable criteria.
   *
   * @param confidence The confidence level from findMatch
   * @return true if user confirmation should be requested
   */
  def requiresUserConfirmation(confidence: String): Boolean = {
    confidence == "LOW"
  }

  /**
   * Checks if we should auto-group (no user interaction needed).
   * HIGH and MEDIUM confidence matches are auto-grouped.
   *
   * @param confidence The confidence level from findMatch
   * @return true if the match can be auto-grouped
   */
  def canAutoGroup(confidence: String): Boolean = {
    confidence == "HIGH" || confidence == "MEDIUM"
  }
}
