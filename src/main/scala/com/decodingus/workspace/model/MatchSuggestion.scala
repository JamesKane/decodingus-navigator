package com.decodingus.workspace.model

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
                            reasonType: String,
                            reasonDetail: Option[String] = None,
                            populationOverlap: Option[Double] = None,
                            dismissed: Boolean = false
                          )

object MatchSuggestion:
  val ReasonTypes: Set[String] = Set(
    "POPULATION_OVERLAP",
    "HAPLOGROUP_PROXIMITY",
    "GEOGRAPHIC_CLUSTER",
    "PROJECT_MEMBER",
    "MUTUAL_INTEREST"
  )
