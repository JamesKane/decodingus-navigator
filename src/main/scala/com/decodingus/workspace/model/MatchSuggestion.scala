package com.decodingus.workspace.model

enum SuggestionReason:
  case PopulationOverlap, HaplogroupProximity, GeographicCluster, ProjectMember, MutualInterest

object SuggestionReason:
  def fromString(s: String): SuggestionReason = s match
    case "POPULATION_OVERLAP" => PopulationOverlap
    case "HAPLOGROUP_PROXIMITY" => HaplogroupProximity
    case "GEOGRAPHIC_CLUSTER" => GeographicCluster
    case "PROJECT_MEMBER" => ProjectMember
    case "MUTUAL_INTEREST" => MutualInterest
    case other => throw new IllegalArgumentException(s"Unknown suggestion reason: $other")

  extension (sr: SuggestionReason)
    def toDbString: String = sr match
      case PopulationOverlap => "POPULATION_OVERLAP"
      case HaplogroupProximity => "HAPLOGROUP_PROXIMITY"
      case GeographicCluster => "GEOGRAPHIC_CLUSTER"
      case ProjectMember => "PROJECT_MEMBER"
      case MutualInterest => "MUTUAL_INTEREST"

/**
 * A match suggestion from the AppView discovery engine.
 * Not persisted locally — fetched on demand from the AppView API.
 *
 * @param suggestionId   Unique ID assigned by AppView
 * @param biosampleRef   AT URI of the local biosample this suggestion is for
 * @param matchedBiosampleRef AT URI of the suggested match biosample
 * @param matchedDid     DID of the suggested match citizen (if consent allows)
 * @param matchedLabel   Display label for the match (anonymized or real name depending on consent)
 * @param score          Discovery score (0.0 to 1.0, higher = stronger signal)
 * @param reasonType     Why this match was suggested
 * @param reasonDetail   Human-readable explanation of the suggestion reason
 * @param populationOverlap Estimated population overlap percentage (if available)
 * @param dismissed      Whether the user has dismissed this suggestion
 */
case class MatchSuggestion(
                            suggestionId: String,
                            biosampleRef: String,
                            matchedBiosampleRef: String,
                            matchedDid: Option[String] = None,
                            matchedLabel: String,
                            score: Double,
                            reasonType: SuggestionReason,
                            reasonDetail: Option[String] = None,
                            populationOverlap: Option[Double] = None,
                            dismissed: Boolean = false
                          )
