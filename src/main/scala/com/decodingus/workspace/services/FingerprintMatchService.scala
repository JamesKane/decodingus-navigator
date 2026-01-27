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
   * Matching logic (in order of precedence):
   * 1. Exact Fingerprint Match (HIGH): The computed fingerprint string is identical.
   * 2. Platform Unit Match (HIGH): The PU tag matches (and is not Unknown).
   * 3. Library + Sample Match (MEDIUM): LB and SM tags match (and are not Unknown).
   *
   * If none of these match, it is treated as a new sequencing run.
   *
   * @param candidateRuns Sequence runs to search within (typically all runs for a biosample)
   * @param fingerprint   The computed run fingerprint from LibraryStats
   * @param libraryStats  Full stats for matching criteria
   * @return Match result with confidence level
   */
  def findMatch(
                 candidateRuns: List[(SequenceRun, Int)],
                 fingerprint: String,
                 libraryStats: LibraryStats
               ): FingerprintMatchResult = {

    // Helper to check for valid (known) values
    def isValid(value: Option[String]): Boolean = value.exists(v => v.nonEmpty && v != "Unknown")
    def isValidStr(value: String): Boolean = value.nonEmpty && value != "Unknown"

    // 1. Exact Fingerprint Match (HIGH)
    val fingerprintMatch = candidateRuns.find { case (run, _) =>
      run.runFingerprint.contains(fingerprint)
    }

    if (fingerprintMatch.isDefined) {
      val (run, idx) = fingerprintMatch.get
      return FingerprintMatchResult.MatchFound(run, idx, "HIGH")
    }

    // 2. Platform Unit Match (HIGH)
    // PU is usually globally unique for a lane/run
    if (isValid(libraryStats.platformUnit)) {
      val puMatch = candidateRuns.find { case (run, _) =>
        isValid(run.platformUnit) && run.platformUnit == libraryStats.platformUnit
      }

      if (puMatch.isDefined) {
        val (run, idx) = puMatch.get
        return FingerprintMatchResult.MatchFound(run, idx, "HIGH")
      }
    }

    // 3. Library + Sample Match (MEDIUM)
    // Same library and sample name usually implies same physical sequencing run,
    // but without PU we can't be 100% certain it's not a re-run of the same library.
    if (isValidStr(libraryStats.libraryId) && isValidStr(libraryStats.sampleName)) {
      val lbSmMatch = candidateRuns.find { case (run, _) =>
        isValid(run.libraryId) &&
          isValid(run.sampleName) &&
          run.libraryId.contains(libraryStats.libraryId) &&
          run.sampleName.contains(libraryStats.sampleName)
      }

      if (lbSmMatch.isDefined) {
        val (run, idx) = lbSmMatch.get
        return FingerprintMatchResult.MatchFound(run, idx, "MEDIUM")
      }
    }

    // No reliable match found
    FingerprintMatchResult.NoMatch
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