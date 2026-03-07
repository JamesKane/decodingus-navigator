package com.decodingus.workspace.model

import munit.FunSuite

class MatchSuggestionSpec extends FunSuite:

  test("MatchSuggestion stores discovery data") {
    val suggestion = MatchSuggestion(
      suggestionId = "sug-001",
      biosampleRef = "at://did:plc:a/bio/1",
      matchedBiosampleRef = "at://did:plc:b/bio/1",
      matchedDid = Some("did:plc:b"),
      matchedLabel = "Anonymous Match #42",
      score = 0.85,
      reasonType = SuggestionReason.PopulationOverlap,
      reasonDetail = Some("Shared 92% population overlap in Northern European reference panel"),
      populationOverlap = Some(92.0)
    )
    assertEquals(suggestion.suggestionId, "sug-001")
    assertEquals(suggestion.score, 0.85)
    assertEquals(suggestion.reasonType, SuggestionReason.PopulationOverlap)
    assert(suggestion.populationOverlap.contains(92.0))
    assert(!suggestion.dismissed)
  }

  test("MatchSuggestion defaults") {
    val suggestion = MatchSuggestion(
      suggestionId = "sug-002",
      biosampleRef = "at://did:plc:a/bio/1",
      matchedBiosampleRef = "at://did:plc:c/bio/1",
      matchedLabel = "Anonymous",
      score = 0.5,
      reasonType = SuggestionReason.HaplogroupProximity
    )
    assert(suggestion.matchedDid.isEmpty)
    assert(suggestion.reasonDetail.isEmpty)
    assert(suggestion.populationOverlap.isEmpty)
    assert(!suggestion.dismissed)
  }

  test("SuggestionReason has all expected values") {
    assertEquals(SuggestionReason.values.length, 5)
    assertEquals(SuggestionReason.fromString("POPULATION_OVERLAP"), SuggestionReason.PopulationOverlap)
    assertEquals(SuggestionReason.fromString("HAPLOGROUP_PROXIMITY"), SuggestionReason.HaplogroupProximity)
    assertEquals(SuggestionReason.fromString("GEOGRAPHIC_CLUSTER"), SuggestionReason.GeographicCluster)
    assertEquals(SuggestionReason.fromString("PROJECT_MEMBER"), SuggestionReason.ProjectMember)
    assertEquals(SuggestionReason.fromString("MUTUAL_INTEREST"), SuggestionReason.MutualInterest)
  }

  test("MatchSuggestion dismissed flag") {
    val suggestion = MatchSuggestion(
      suggestionId = "sug-003",
      biosampleRef = "at://did:plc:a/bio/1",
      matchedBiosampleRef = "at://did:plc:d/bio/1",
      matchedLabel = "Anonymous",
      score = 0.3,
      reasonType = SuggestionReason.ProjectMember,
      dismissed = true
    )
    assert(suggestion.dismissed)
  }
