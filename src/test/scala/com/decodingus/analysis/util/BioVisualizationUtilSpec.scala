package com.decodingus.analysis.util

import munit.FunSuite
import java.nio.file.{Files, Path}

class BioVisualizationUtilSpec extends FunSuite {

  test("binIntervalsFromBed correctly bins intervals") {
    val tempBed = Files.createTempFile("test_callable", ".bed")
    try {
      // Create a test BED file
      // contig, start, stop, state
      // Stride is default 10000 (10kb)
      // Bin 0: 0-10000
      // Bin 1: 10000-20000
      // Bin 2: 20000-30000
      val content =
        """chr1 0 5000 CALLABLE
          |chr1 5000 10000 POOR_MAPPING_QUALITY
          |chr1 10000 15000 CALLABLE
          |chr1 15000 20000 NO_COVERAGE
          |chr1 20000 25000 CALLABLE
          |chr1 25000 30000 CALLABLE
          |""".stripMargin
      Files.writeString(tempBed, content)

      // Total length 30000 -> 3 bins
      val bins = BioVisualizationUtil.binIntervalsFromBed(tempBed, "chr1", 30000, 10000)

      assertEquals(bins.length, 3)

      // Bin 0: 5000 CALLABLE, 5000 POOR_MAPPING_QUALITY
      assertEquals(bins(0)(0), 5000) // CALLABLE
      assertEquals(bins(0)(1), 5000) // POOR_MAPPING_QUALITY
      assertEquals(bins(0)(2), 0)    // Other

      // Bin 1: 5000 CALLABLE, 5000 NO_COVERAGE (Other)
      assertEquals(bins(1)(0), 5000) // CALLABLE
      assertEquals(bins(1)(1), 0)    // POOR_MAPPING_QUALITY
      assertEquals(bins(1)(2), 5000) // Other

      // Bin 2: 10000 CALLABLE
      assertEquals(bins(2)(0), 10000) // CALLABLE
      assertEquals(bins(2)(1), 0)     // POOR_MAPPING_QUALITY
      assertEquals(bins(2)(2), 0)     // Other

    } finally {
      Files.deleteIfExists(tempBed)
    }
  }

  test("generateSvgForContig produces valid SVG structure") {
    // 3 bins, stride 10000
    val bins = Array(
      Array(10000, 0, 0), // Full callable
      Array(5000, 5000, 0), // Mixed callable/poor
      Array(0, 0, 10000) // Full other
    )

    val svg = BioVisualizationUtil.generateSvgForContig("chr1", 30000, 30000, bins)

    assert(svg.contains("<svg width="))
    assert(svg.contains("chr1 (Stride: 10kb)"))
    assert(svg.contains(BioVisualizationUtil.COLOR_GREEN)) // Should have green bar
    assert(svg.contains(BioVisualizationUtil.COLOR_RED))   // Should have red bar
    assert(svg.contains(BioVisualizationUtil.COLOR_GREY))  // Should have grey bar
    assert(svg.contains("</svg>"))
  }

  test("generateSvgForContig handles empty bins gracefully") {
    val bins = Array(Array(0, 0, 0))
    val svg = BioVisualizationUtil.generateSvgForContig("chr1", 10000, 10000, bins)
    
    assert(svg.contains("<svg"))
    assert(svg.contains("</svg>"))
    // Should not crash even if bars are height 0
  }
}
