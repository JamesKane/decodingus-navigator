package com.decodingus.str

import com.decodingus.workspace.model.*

/**
 * Result of comparing STR markers across multiple profiles.
 *
 * @param allMarkers       All markers found across profiles: marker -> List[(provider, value)]
 * @param conflicts        Markers with differing values between providers
 * @param uniqueToProvider Markers unique to each provider
 * @param agreementCount   Number of markers where all providers agree
 */
case class MarkerComparisonResult(
  allMarkers: Map[String, List[(String, StrValue)]],
  conflicts: List[MarkerConflict],
  uniqueToProvider: Map[String, Set[String]],
  agreementCount: Int
) {
  /** Total markers across all providers */
  def totalMarkers: Int = allMarkers.size

  /** Number of markers with conflicts */
  def conflictCount: Int = conflicts.size

  /** Whether any conflicts exist */
  def hasConflicts: Boolean = conflicts.nonEmpty
}

/**
 * A marker value conflict between providers.
 *
 * @param markerName Normalized marker name
 * @param values     Map of provider -> value for this marker
 */
case class MarkerConflict(
  markerName: String,
  values: Map[String, StrValue]
) {
  /** Format for display in tooltip */
  def formatForDisplay: String = {
    val lines = values.map { case (provider, value) =>
      s"  $provider: ${formatValue(value)}"
    }.mkString("\n")
    s"$markerName:\n$lines"
  }

  private def formatValue(value: StrValue): String = value match {
    case SimpleStrValue(repeats) => repeats.toString
    case MultiCopyStrValue(copies) => copies.mkString("-")
    case ComplexStrValue(alleles, rawNotation) =>
      rawNotation.getOrElse(alleles.map(a => s"${a.repeats}${a.designation.getOrElse("")}").mkString("-"))
  }
}

/**
 * Service for comparing STR markers across profiles from different providers.
 * Used to detect discrepancies in marker values between FTDNA, YSEQ, and other labs.
 */
object StrMarkerComparator {

  /**
   * Compare STR values for equality.
   * Handles different value types: simple, multi-copy, and complex.
   */
  def valuesMatch(v1: StrValue, v2: StrValue): Boolean = {
    (v1, v2) match {
      case (SimpleStrValue(r1), SimpleStrValue(r2)) =>
        r1 == r2

      case (MultiCopyStrValue(c1), MultiCopyStrValue(c2)) =>
        // Compare sorted lists to handle order differences
        c1.sorted == c2.sorted

      case (ComplexStrValue(a1, _), ComplexStrValue(a2, _)) =>
        // Compare alleles by repeat count and designation
        val sorted1 = a1.sortBy(a => (a.repeats, a.designation.getOrElse("")))
        val sorted2 = a2.sortBy(a => (a.repeats, a.designation.getOrElse("")))
        sorted1 == sorted2

      case _ =>
        // Different types - try to compare by string representation
        formatValue(v1) == formatValue(v2)
    }
  }

  /**
   * Compare markers from multiple STR profiles.
   *
   * @param profiles List of STR profiles to compare
   * @return Comparison result with conflicts and statistics
   */
  def compare(profiles: List[StrProfile]): MarkerComparisonResult = {
    if (profiles.isEmpty) {
      return MarkerComparisonResult(Map.empty, Nil, Map.empty, 0)
    }

    // Build map: markerName -> List[(provider, value)]
    val allMarkers: Map[String, List[(String, StrValue)]] = profiles
      .flatMap { profile =>
        val provider = profile.importedFrom.getOrElse("UNKNOWN")
        profile.markers.map { m =>
          (m.marker, provider, m.value)
        }
      }
      .groupBy(_._1)
      .view
      .mapValues(_.map(e => (e._2, e._3)))
      .toMap

    // Find conflicts (same marker, different values across providers)
    val conflicts = allMarkers.flatMap { case (marker, entries) =>
      val byProvider = entries.toMap // provider -> value
      if (byProvider.size > 1) {
        // Multiple providers have this marker - check for value differences
        val distinctValues = byProvider.values.toList
          .foldLeft(List.empty[StrValue]) { (acc, v) =>
            if (acc.exists(existing => valuesMatch(existing, v))) acc
            else acc :+ v
          }

        if (distinctValues.size > 1) {
          Some(MarkerConflict(marker, byProvider))
        } else {
          None
        }
      } else {
        None
      }
    }.toList.sortBy(_.markerName)

    // Find markers unique to each provider
    val providers = profiles.flatMap(_.importedFrom).distinct
    val uniqueToProvider = providers.map { provider =>
      val providerMarkers = allMarkers.filter { case (_, entries) =>
        entries.size == 1 && entries.head._1 == provider
      }.keySet
      provider -> providerMarkers
    }.toMap

    // Count agreements (markers where all providers agree)
    val agreementCount = allMarkers.count { case (_, entries) =>
      val byProvider = entries.toMap
      if (byProvider.size <= 1) false // Need multiple providers to "agree"
      else {
        val values = byProvider.values.toList
        values.tail.forall(v => valuesMatch(values.head, v))
      }
    }

    MarkerComparisonResult(allMarkers, conflicts, uniqueToProvider, agreementCount)
  }

  /**
   * Compare markers from two specific providers.
   *
   * @param profile1 First profile
   * @param profile2 Second profile
   * @return Comparison result
   */
  def compareTwoProfiles(profile1: StrProfile, profile2: StrProfile): MarkerComparisonResult = {
    compare(List(profile1, profile2))
  }

  /**
   * Format a value for display.
   */
  private def formatValue(value: StrValue): String = value match {
    case SimpleStrValue(repeats) => repeats.toString
    case MultiCopyStrValue(copies) => copies.mkString("-")
    case ComplexStrValue(alleles, rawNotation) =>
      rawNotation.getOrElse(alleles.map(a => s"${a.repeats}${a.designation.getOrElse("")}").mkString("-"))
  }
}
