package com.decodingus.ibd.engine

import com.decodingus.workspace.model.RelationshipEstimate
import munit.FunSuite

class RelationshipEstimatorSpec extends FunSuite:

  test("parent-child at ~3400 cM") {
    assertEquals(RelationshipEstimator.estimate(3500), RelationshipEstimate.ParentChild)
  }

  test("full sibling at ~2550 cM") {
    assertEquals(RelationshipEstimator.estimate(2600), RelationshipEstimate.FullSibling)
  }

  test("grandparent at ~1700 cM") {
    assertEquals(RelationshipEstimator.estimate(1800), RelationshipEstimate.Grandparent)
  }

  test("aunt/uncle at ~1200 cM") {
    assertEquals(RelationshipEstimator.estimate(1300), RelationshipEstimate.AuntUncle)
  }

  test("1st cousin at ~880 cM") {
    assertEquals(RelationshipEstimator.estimate(880), RelationshipEstimate.FirstCousin)
  }

  test("1st cousin once removed at ~440 cM") {
    assertEquals(RelationshipEstimator.estimate(440), RelationshipEstimate.FirstCousinOnceRemoved)
  }

  test("2nd cousin at ~230 cM") {
    assertEquals(RelationshipEstimator.estimate(230), RelationshipEstimate.SecondCousin)
  }

  test("3rd cousin at ~73 cM") {
    assertEquals(RelationshipEstimator.estimate(73), RelationshipEstimate.ThirdCousin)
  }

  test("4th cousin at ~35 cM") {
    assertEquals(RelationshipEstimator.estimate(35), RelationshipEstimate.FourthCousin)
  }

  test("5th cousin at ~18 cM") {
    assertEquals(RelationshipEstimator.estimate(18), RelationshipEstimate.FifthCousin)
  }

  test("distant relative at low cM") {
    assertEquals(RelationshipEstimator.estimate(8), RelationshipEstimate.Distant)
  }

  test("unknown below threshold") {
    assertEquals(RelationshipEstimator.estimate(5), RelationshipEstimate.Unknown)
  }

  test("labels are human-readable") {
    assertEquals(RelationshipEstimate.ParentChild.label, "Parent/Child")
    assertEquals(RelationshipEstimate.SecondCousin.label, "2nd Cousin")
    assertEquals(RelationshipEstimate.Unknown.label, "Unknown")
  }
