package com.decodingus.ibd.engine

import munit.FunSuite

class GeneticMapSpec extends FunSuite:

  test("uniformRate map converts positions correctly") {
    val gmap = GeneticMap.uniformRate(cmPerMb = 1.0)

    // 1 Mb = 1.0 cM at uniform rate
    val cm = gmap.positionToCm("1", 1_000_001)
    assert(cm.isDefined)
    assert(math.abs(cm.get - 1.0) < 0.01, s"Expected ~1.0 cM, got ${cm.get}")
  }

  test("uniformRate map interval calculation") {
    val gmap = GeneticMap.uniformRate(cmPerMb = 1.0)

    val cm = gmap.intervalCm("1", 10_000_000, 20_000_000)
    assert(cm.isDefined)
    assert(math.abs(cm.get - 10.0) < 0.1, s"Expected ~10.0 cM, got ${cm.get}")
  }

  test("uniformRate map handles all autosomes") {
    val gmap = GeneticMap.uniformRate()
    for chr <- 1 to 22 do
      assert(gmap.hasChromosome(chr.toString), s"Missing chromosome $chr")
  }

  test("normalizeChromosome handles chr prefix") {
    assertEquals(GeneticMap.normalizeChromosome("chr1"), "1")
    assertEquals(GeneticMap.normalizeChromosome("chr22"), "22")
    assertEquals(GeneticMap.normalizeChromosome("chrX"), "X")
    assertEquals(GeneticMap.normalizeChromosome("1"), "1")
    assertEquals(GeneticMap.normalizeChromosome("X"), "X")
  }

  test("ChromosomeMap interpolation exact match") {
    val cmap = GeneticMap.ChromosomeMap(
      Array(1000, 2000, 3000),
      Array(0.0, 1.0, 3.0)
    )
    assertEquals(cmap.interpolate(2000), 1.0)
  }

  test("ChromosomeMap interpolation between points") {
    val cmap = GeneticMap.ChromosomeMap(
      Array(1000, 3000),
      Array(0.0, 2.0)
    )
    // Midpoint should be 1.0
    assertEquals(cmap.interpolate(2000), 1.0)
  }

  test("ChromosomeMap extrapolation before first point") {
    val cmap = GeneticMap.ChromosomeMap(
      Array(1000, 2000),
      Array(0.0, 1.0)
    )
    // Rate = 1.0 cM / 1000 bp = 0.001 cM/bp
    // At position 500: 0.0 + 0.001 * (500 - 1000) = -0.5
    val result = cmap.interpolate(500)
    assert(math.abs(result - (-0.5)) < 0.001)
  }

  test("ChromosomeMap extrapolation after last point") {
    val cmap = GeneticMap.ChromosomeMap(
      Array(1000, 2000),
      Array(0.0, 1.0)
    )
    // Rate = 1.0 cM / 1000 bp = 0.001 cM/bp
    // At position 3000: 1.0 + 0.001 * (3000 - 2000) = 2.0
    val result = cmap.interpolate(3000)
    assert(math.abs(result - 2.0) < 0.001)
  }

  test("positionToCm returns None for unknown chromosome") {
    val gmap = GeneticMap.uniformRate()
    assert(gmap.positionToCm("99", 1000).isEmpty)
  }

  test("intervalCm returns None for unknown chromosome") {
    val gmap = GeneticMap.uniformRate()
    assert(gmap.intervalCm("99", 1000, 2000).isEmpty)
  }
