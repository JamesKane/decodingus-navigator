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
      reasonType = "POPULATION_OVERLAP",
      reasonDetail = Some("Shared 92% population overlap in Northern European reference panel"),
      populationOverlap = Some(92.0)
    )
    assertEquals(suggestion.suggestionId, "sug-001")
    assertEquals(suggestion.score, 0.85)
    assertEquals(suggestion.reasonType, "POPULATION_OVERLAP")
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
      reasonType = "HAPLOGROUP_PROXIMITY"
    )
    assert(suggestion.matchedDid.isEmpty)
    assert(suggestion.reasonDetail.isEmpty)
    assert(suggestion.populationOverlap.isEmpty)
    assert(!suggestion.dismissed)
  }

  test("MatchSuggestion.ReasonTypes contains expected values") {
    assert(MatchSuggestion.ReasonTypes.contains("POPULATION_OVERLAP"))
    assert(MatchSuggestion.ReasonTypes.contains("HAPLOGROUP_PROXIMITY"))
    assert(MatchSuggestion.ReasonTypes.contains("GEOGRAPHIC_CLUSTER"))
    assert(MatchSuggestion.ReasonTypes.contains("PROJECT_MEMBER"))
    assert(MatchSuggestion.ReasonTypes.contains("MUTUAL_INTEREST"))
  }

  test("MatchSuggestion dismissed flag") {
    val suggestion = MatchSuggestion(
      suggestionId = "sug-003",
      biosampleRef = "at://did:plc:a/bio/1",
      matchedBiosampleRef = "at://did:plc:d/bio/1",
      matchedLabel = "Anonymous",
      score = 0.3,
      reasonType = "PROJECT_MEMBER",
      dismissed = true
    )
    assert(suggestion.dismissed)
  }
