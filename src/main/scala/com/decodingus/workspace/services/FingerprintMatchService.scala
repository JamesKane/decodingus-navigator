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
   * Uses a tiered matching approach with different confidence levels.
   *
   * Matching tiers:
   * 1. Exact fingerprint hash match → HIGH confidence
   * 2. Platform Unit (PU) match → HIGH confidence
   * 3. Library ID + Sample Name match → MEDIUM confidence (may require user confirmation)
   *
   * @param candidateRuns     Sequence runs to search within (typically all runs for a biosample)
   * @param fingerprint       The computed run fingerprint from LibraryStats
   * @param libraryStats      Full stats for additional matching criteria
   * @return Match result with confidence level
   */
  def findMatch(
    candidateRuns: List[(SequenceRun, Int)],
    fingerprint: String,
    libraryStats: LibraryStats
  ): FingerprintMatchResult = {
    // Tier 1: Exact fingerprint match (HIGH confidence)
    candidateRuns.find { case (run, _) =>
      run.runFingerprint.contains(fingerprint)
    }.map { case (run, idx) =>
      FingerprintMatchResult.MatchFound(run, idx, "HIGH")
    }.getOrElse {
      // Tier 2: PU match (HIGH confidence)
      libraryStats.platformUnit.flatMap { pu =>
        candidateRuns.find { case (run, _) =>
          run.platformUnit.contains(pu)
        }.map { case (run, idx) =>
          FingerprintMatchResult.MatchFound(run, idx, "HIGH")
        }
      }.getOrElse {
        // Tier 3: LB + SM match (MEDIUM confidence)
        if (libraryStats.libraryId != "Unknown" && libraryStats.sampleName != "Unknown") {
          candidateRuns.find { case (run, _) =>
            run.libraryId.contains(libraryStats.libraryId) &&
            run.sampleName.contains(libraryStats.sampleName)
          }.map { case (run, idx) =>
            FingerprintMatchResult.MatchFound(run, idx, "MEDIUM")
          }.getOrElse(FingerprintMatchResult.NoMatch)
        } else {
          FingerprintMatchResult.NoMatch
        }
      }
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
