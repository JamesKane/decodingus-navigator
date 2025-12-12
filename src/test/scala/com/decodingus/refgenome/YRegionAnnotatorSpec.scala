package com.decodingus.refgenome

import munit.FunSuite

class YRegionAnnotatorSpec extends FunSuite {

  test("RegionType has correct quality modifiers") {
    assertEquals(RegionType.XDegenerate.modifier, 1.0)
    assertEquals(RegionType.Normal.modifier, 1.0)
    assertEquals(RegionType.PAR.modifier, 0.5)
    assertEquals(RegionType.Palindrome.modifier, 0.4)
    assertEquals(RegionType.XTR.modifier, 0.3)
    assertEquals(RegionType.Ampliconic.modifier, 0.3)
    assertEquals(RegionType.STR.modifier, 0.25)
    assertEquals(RegionType.Centromere.modifier, 0.1)
    assertEquals(RegionType.Heterochromatin.modifier, 0.1)
    assertEquals(RegionType.NonCallable.modifier, 0.5)
    assertEquals(RegionType.LowDepth.modifier, 0.7)
  }

  test("GenomicRegion.contains works correctly") {
    val region = GenomicRegion("chrY", 1000, 2000, RegionType.Palindrome, Some("P1"))

    assert(region.contains(1000))
    assert(region.contains(1500))
    assert(region.contains(2000))
    assert(!region.contains(999))
    assert(!region.contains(2001))
  }

  test("YRegionAnnotator.annotate returns empty for non-Y chromosome") {
    val annotator = YRegionAnnotator.empty

    val result = annotator.annotate("chrX", 1000)
    assertEquals(result, RegionAnnotation.empty)
  }

  test("YRegionAnnotator handles Y chromosome variants") {
    val annotator = YRegionAnnotator.empty

    // Should return empty annotation but not fail
    val result1 = annotator.annotate("chrY", 1000)
    assertEquals(result1.qualityModifier, 1.0)

    val result2 = annotator.annotate("Y", 1000)
    assertEquals(result2.qualityModifier, 1.0)
  }

  test("YRegionAnnotator finds overlapping regions") {
    val palindromes = List(
      GenomicRegion("chrY", 1000, 2000, RegionType.Palindrome, Some("P1")),
      GenomicRegion("chrY", 5000, 6000, RegionType.Palindrome, Some("P2"))
    )

    val annotator = YRegionAnnotator.fromRegions(palindromes = palindromes)

    val result1 = annotator.annotate("chrY", 1500)
    assertEquals(result1.qualityModifier, 0.4)
    assertEquals(result1.primaryRegion.map(_.name), Some(Some("P1")))

    val result2 = annotator.annotate("chrY", 3000)
    assertEquals(result2.qualityModifier, 1.0)
    assert(result2.primaryRegion.isEmpty)
  }

  test("YRegionAnnotator combines modifiers multiplicatively") {
    val palindromes = List(GenomicRegion("chrY", 1000, 2000, RegionType.Palindrome, Some("P1")))
    val strs = List(GenomicRegion("chrY", 1500, 2500, RegionType.STR, Some("DYS389")))

    val annotator = YRegionAnnotator.fromRegions(palindromes = palindromes, strs = strs)

    // Position 1500 overlaps both palindrome and STR
    val result = annotator.annotate("chrY", 1500)

    // Modifiers combine: 0.4 * 0.25 = 0.1
    assertEqualsDouble(result.qualityModifier, 0.1, 0.001)
  }

  test("YRegionAnnotator includes cytoband for display") {
    val cytobands = List(
      GenomicRegion("chrY", 1, 2781479, RegionType.Cytoband, Some("Yp11.32")),
      GenomicRegion("chrY", 2781480, 3000000, RegionType.Cytoband, Some("Yp11.31"))
    )

    val annotator = YRegionAnnotator.fromRegions(cytobands = cytobands)

    val result = annotator.annotate("chrY", 1000000)
    assertEquals(result.cytoband.flatMap(_.name), Some("Yp11.32"))

    // Cytobands don't affect quality modifier
    assertEquals(result.qualityModifier, 1.0)
  }

  test("YRegionAnnotator applies low depth modifier") {
    val annotator = YRegionAnnotator.empty

    // Depth >= 10 - no modifier
    val result1 = annotator.annotate("chrY", 1000, Some(10))
    assertEquals(result1.qualityModifier, 1.0)

    // Depth < 10 - applies 0.7x modifier
    val result2 = annotator.annotate("chrY", 1000, Some(5))
    assertEquals(result2.qualityModifier, 0.7)
  }

  test("YRegionAnnotator applies non-callable modifier") {
    val callablePositions = Set(1000L, 2000L, 3000L)
    val annotator = YRegionAnnotator.fromRegions(callablePositions = Some(callablePositions))

    // Callable position - no modifier
    val result1 = annotator.annotate("chrY", 1000)
    assertEquals(result1.qualityModifier, 1.0)

    // Non-callable position - applies 0.5x modifier
    val result2 = annotator.annotate("chrY", 1500)
    assertEquals(result2.qualityModifier, 0.5)
  }

  test("RegionAnnotation.description formats correctly") {
    val palindrome = GenomicRegion("chrY", 1000, 2000, RegionType.Palindrome, Some("P8"))
    val cytoband = GenomicRegion("chrY", 1, 5000000, RegionType.Cytoband, Some("Yq11.223"))

    val ann1 = RegionAnnotation(List(palindrome, cytoband), 0.4, Some(cytoband), Some(palindrome))
    assertEquals(ann1.description, "P8 Palindrome (Yq11.223)")

    val ann2 = RegionAnnotation(List(cytoband), 1.0, Some(cytoband), None)
    assertEquals(ann2.description, "Yq11.223")

    val ann3 = RegionAnnotation.empty
    assertEquals(ann3.description, "-")
  }

  test("GRCh38 heterochromatin boundaries are defined") {
    val heterochromatin = YRegionAnnotator.grch38Heterochromatin
    assertEquals(heterochromatin.size, 1)
    assertEquals(heterochromatin.head.name, Some("Yq12"))
    assertEquals(heterochromatin.head.start, 26673237L)
    assertEquals(heterochromatin.head.end, 56887902L)
  }

  test("YRegionAnnotator.gff3ToRegions converts records correctly") {
    val records = List(
      Gff3Record("chrY", "ybrowse", "palindrome", 1000, 2000, None, ".", None, Map("Name" -> "P8"))
    )

    val regions = YRegionAnnotator.gff3ToRegions(records, RegionType.Palindrome)
    assertEquals(regions.size, 1)
    assertEquals(regions.head.start, 1000L)
    assertEquals(regions.head.end, 2000L)
    assertEquals(regions.head.name, Some("P8"))
    assertEquals(regions.head.regionType, RegionType.Palindrome)
  }

  test("YRegionAnnotator.bedToRegions converts coordinates correctly") {
    val records = List(
      BedRecord("chrY", 999, 2000, Some("PAR1"))  // BED is 0-based
    )

    val regions = YRegionAnnotator.bedToRegions(records, RegionType.PAR)
    assertEquals(regions.size, 1)
    assertEquals(regions.head.start, 1000L)  // Converted to 1-based
    assertEquals(regions.head.end, 2000L)
    assertEquals(regions.head.name, Some("PAR1"))
    assertEquals(regions.head.regionType, RegionType.PAR)
  }
}
