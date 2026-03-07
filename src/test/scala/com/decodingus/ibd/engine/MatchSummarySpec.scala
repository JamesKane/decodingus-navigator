package com.decodingus.ibd.engine

import com.decodingus.workspace.model.{IbdSegment, RelationshipEstimate}
import munit.FunSuite

class MatchSummarySpec extends FunSuite:

  test("fromSegments computes correct totals") {
    val segments = List(
      IbdSegment("1", 1000000, 5000000, 12.5, Some(500), Some(true)),
      IbdSegment("2", 2000000, 8000000, 25.3, Some(800), Some(true)),
      IbdSegment("5", 10000000, 15000000, 8.2, Some(300), Some(true))
    )

    val summary = MatchSummary.fromSegments(segments)
    assertEquals(summary.segmentCount, 3)
    assertEquals(summary.totalSharedCm, 46.0) // 12.5 + 25.3 + 8.2 = 46.0
    assertEquals(summary.longestSegmentCm, 25.3)
    assertEquals(summary.relationshipEstimate, RelationshipEstimate.FourthCousin)
  }

  test("fromSegments handles empty segments") {
    val summary = MatchSummary.fromSegments(Nil)
    assertEquals(summary.segmentCount, 0)
    assertEquals(summary.totalSharedCm, 0.0)
    assertEquals(summary.longestSegmentCm, 0.0)
    assertEquals(summary.relationshipEstimate, RelationshipEstimate.Unknown)
    assert(summary.summaryHash.nonEmpty)
  }

  test("hash is deterministic for same segments") {
    val segments = List(
      IbdSegment("1", 1000000, 5000000, 12.5, Some(500)),
      IbdSegment("2", 2000000, 8000000, 25.3, Some(800))
    )

    val hash1 = MatchSummary.computeHash(segments)
    val hash2 = MatchSummary.computeHash(segments)
    assertEquals(hash1, hash2)
  }

  test("hash is order-independent (canonical sorting)") {
    val seg1 = IbdSegment("1", 1000000, 5000000, 12.5)
    val seg2 = IbdSegment("2", 2000000, 8000000, 25.3)

    val hash1 = MatchSummary.computeHash(List(seg1, seg2))
    val hash2 = MatchSummary.computeHash(List(seg2, seg1))
    assertEquals(hash1, hash2, "Hash should be order-independent")
  }

  test("hash differs for different segments") {
    val seg1 = List(IbdSegment("1", 1000000, 5000000, 12.5))
    val seg2 = List(IbdSegment("1", 1000000, 6000000, 15.0))

    val hash1 = MatchSummary.computeHash(seg1)
    val hash2 = MatchSummary.computeHash(seg2)
    assert(hash1 != hash2, "Different segments should produce different hashes")
  }

  test("hash is consistent across chromosome sort order") {
    val segments = List(
      IbdSegment("10", 1000000, 5000000, 12.5),
      IbdSegment("2", 2000000, 8000000, 25.3),
      IbdSegment("1", 3000000, 7000000, 10.0)
    )

    // Canonical sort: 1, 2, 10 (numeric)
    val hash1 = MatchSummary.computeHash(segments)
    val hash2 = MatchSummary.computeHash(segments.reverse)
    assertEquals(hash1, hash2)
  }

  test("summary hash matches for two parties with same computation") {
    // Simulate both parties computing the same IBD segments
    val aliceSegments = List(
      IbdSegment("1", 5000000, 10000000, 15.23),
      IbdSegment("3", 20000000, 30000000, 28.76)
    )
    val bobSegments = List(
      IbdSegment("3", 20000000, 30000000, 28.76),
      IbdSegment("1", 5000000, 10000000, 15.23)
    )

    val aliceSummary = MatchSummary.fromSegments(aliceSegments)
    val bobSummary = MatchSummary.fromSegments(bobSegments)

    assertEquals(aliceSummary.summaryHash, bobSummary.summaryHash,
      "Both parties should produce identical hashes")
    assertEquals(aliceSummary.totalSharedCm, bobSummary.totalSharedCm)
    assertEquals(aliceSummary.segmentCount, bobSummary.segmentCount)
  }
