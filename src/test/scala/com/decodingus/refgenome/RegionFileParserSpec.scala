package com.decodingus.refgenome

import munit.FunSuite
import java.nio.file.{Files, Path}

class RegionFileParserSpec extends FunSuite {

  def withTempFile[T](content: String, suffix: String)(f: Path => T): T = {
    val file = Files.createTempFile("test-", suffix)
    try {
      Files.writeString(file, content)
      f(file)
    } finally {
      Files.deleteIfExists(file)
    }
  }

  test("parseGff3 parses valid GFF3 content") {
    val gff3Content =
      """##gff-version 3
        |chrY	ybrowse	palindrome	14969754	15077740	.	.	.	Name=P6;Note=Palindrome P6
        |chrY	ybrowse	palindrome	22234328	22418595	.	.	.	Name=P5;Note=Palindrome P5
        |""".stripMargin

    withTempFile(gff3Content, ".gff3") { file =>
      val result = RegionFileParser.parseGff3(file)
      assert(result.isRight)
      val records = result.toOption.get
      assertEquals(records.size, 2)

      val p6 = records.head
      assertEquals(p6.seqId, "chrY")
      assertEquals(p6.featureType, "palindrome")
      assertEquals(p6.start, 14969754L)
      assertEquals(p6.end, 15077740L)
      assertEquals(p6.name, Some("P6"))
      assertEquals(p6.note, Some("Palindrome P6"))
    }
  }

  test("parseGff3 skips comment lines") {
    val gff3Content =
      """##gff-version 3
        |# This is a comment
        |chrY	ybrowse	cytoband	1	2781479	.	.	.	Name=Yp11.32
        |""".stripMargin

    withTempFile(gff3Content, ".gff3") { file =>
      val result = RegionFileParser.parseGff3(file)
      assert(result.isRight)
      assertEquals(result.toOption.get.size, 1)
    }
  }

  test("parseGff3 handles URL-encoded attributes") {
    val gff3Content =
      """##gff-version 3
        |chrY	ybrowse	str	1000	2000	.	.	.	Name=DYS%20389;Note=STR%20Region
        |""".stripMargin

    withTempFile(gff3Content, ".gff3") { file =>
      val result = RegionFileParser.parseGff3(file)
      assert(result.isRight)
      val records = result.toOption.get
      assertEquals(records.head.name, Some("DYS 389"))
      assertEquals(records.head.note, Some("STR Region"))
    }
  }

  test("parseBed parses valid BED content") {
    val bedContent =
      """chrY	10000	2781479	PAR1	0	.
        |chrY	56887902	57217415	PAR2	0	.
        |""".stripMargin

    withTempFile(bedContent, ".bed") { file =>
      val result = RegionFileParser.parseBed(file)
      assert(result.isRight)
      val records = result.toOption.get
      assertEquals(records.size, 2)

      val par1 = records.head
      assertEquals(par1.chrom, "chrY")
      assertEquals(par1.start, 10000L)  // 0-based
      assertEquals(par1.end, 2781479L)  // exclusive
      assertEquals(par1.name, Some("PAR1"))
    }
  }

  test("parseBed skips header lines") {
    val bedContent =
      """track name="test"
        |browser position chrY:1-1000000
        |#comment
        |chrY	1000	2000	region1
        |""".stripMargin

    withTempFile(bedContent, ".bed") { file =>
      val result = RegionFileParser.parseBed(file)
      assert(result.isRight)
      assertEquals(result.toOption.get.size, 1)
    }
  }

  test("parseBed handles minimal 3-column format") {
    val bedContent =
      """chrY	1000	2000
        |chrY	3000	4000
        |""".stripMargin

    withTempFile(bedContent, ".bed") { file =>
      val result = RegionFileParser.parseBed(file)
      assert(result.isRight)
      val records = result.toOption.get
      assertEquals(records.size, 2)
      assert(records.head.name.isEmpty)
    }
  }

  test("bedToOneBased converts coordinates correctly") {
    // BED: 0-based, half-open [start, end)
    // GFF3/VCF: 1-based, closed [start, end]

    val (start, end) = RegionFileParser.bedToOneBased(0, 100)
    assertEquals(start, 1L)
    assertEquals(end, 100L)

    val (start2, end2) = RegionFileParser.bedToOneBased(999, 2000)
    assertEquals(start2, 1000L)
    assertEquals(end2, 2000L)
  }

  test("oneBasedToBed converts coordinates correctly") {
    val (start, end) = RegionFileParser.oneBasedToBed(1, 100)
    assertEquals(start, 0L)
    assertEquals(end, 100L)

    val (start2, end2) = RegionFileParser.oneBasedToBed(1000, 2000)
    assertEquals(start2, 999L)
    assertEquals(end2, 2000L)
  }

  test("filterYChromosome filters correctly") {
    val records = List(
      BedRecord("chrY", 1000, 2000),
      BedRecord("chrX", 1000, 2000),
      BedRecord("Y", 3000, 4000),
      BedRecord("chr1", 1000, 2000)
    )

    val filtered = RegionFileParser.filterYChromosome(records, _.chrom)
    assertEquals(filtered.size, 2)
    assertEquals(filtered.map(_.chrom).toSet, Set("chrY", "Y"))
  }

  test("Gff3Record getAttribute works for different keys") {
    val record = Gff3Record(
      "chrY", "source", "feature", 1000, 2000, None, ".", None,
      Map("Name" -> "test", "ID" -> "id1", "Note" -> "note text")
    )

    assertEquals(record.getAttribute("Name"), Some("test"))
    assertEquals(record.getAttribute("ID"), Some("id1"))
    assertEquals(record.getAttribute("Note"), Some("note text"))
    assertEquals(record.getAttribute("Missing"), None)
  }

  test("parseGff3 returns error for invalid file") {
    val result = RegionFileParser.parseGff3(Path.of("/nonexistent/file.gff3"))
    assert(result.isLeft)
  }

  test("parseBed returns error for invalid file") {
    val result = RegionFileParser.parseBed(Path.of("/nonexistent/file.bed"))
    assert(result.isLeft)
  }
}
