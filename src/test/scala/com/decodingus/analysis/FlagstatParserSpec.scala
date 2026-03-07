package com.decodingus.analysis

import munit.FunSuite

import java.io.File
import java.nio.file.Files

class FlagstatParserSpec extends FunSuite {

  private def writeTempFlagstat(content: String): File = {
    val file = Files.createTempFile("flagstat_test_", ".flagstat").toFile
    file.deleteOnExit()
    Files.writeString(file.toPath, content)
    file
  }

  private val sampleFlagstat =
    """49876543 + 0 in total (QC-passed reads + QC-failed reads)
      |123456 + 0 secondary
      |78901 + 0 supplementary
      |2345678 + 0 duplicates
      |48500000 + 0 mapped (97.24% : N/A)
      |49674186 + 0 paired in sequencing
      |24837093 + 0 read1
      |24837093 + 0 read2
      |46000000 + 0 properly paired (92.60% : N/A)
      |47500000 + 0 with itself and mate mapped
      |500000 + 0 singletons (1.01% : N/A)
      |100000 + 0 with mate mapped to a different chr
      |50000 + 0 with mate mapped to a different chr (mapQ>=5)""".stripMargin

  test("parse standard samtools flagstat output") {
    val file = writeTempFlagstat(sampleFlagstat)
    val result = FlagstatParser.parse(file)

    assert(result.isRight)
    val fs = result.toOption.get

    assertEquals(fs.totalReads, 49876543L)
    assertEquals(fs.secondary, 123456L)
    assertEquals(fs.supplementary, 78901L)
    assertEquals(fs.duplicates, 2345678L)
    assertEquals(fs.mapped, 48500000L)
    assertEquals(fs.mappedPercent, Some(97.24))
    assertEquals(fs.paired, 49674186L)
    assertEquals(fs.read1, 24837093L)
    assertEquals(fs.read2, 24837093L)
    assertEquals(fs.properlyPaired, 46000000L)
    assertEquals(fs.properlyPairedPercent, Some(92.60))
    assertEquals(fs.withItselfAndMateMapped, 47500000L)
    assertEquals(fs.singletons, 500000L)
    assertEquals(fs.singletonsPercent, Some(1.01))
    assertEquals(fs.mateMappedToDiffChr, 100000L)
    assertEquals(fs.mateMappedToDiffChrMapQ5, 50000L)
  }

  test("isPairedEnd returns true for paired data") {
    val file = writeTempFlagstat(sampleFlagstat)
    val fs = FlagstatParser.parse(file).toOption.get
    assert(fs.isPairedEnd)
  }

  test("primaryReads excludes secondary and supplementary") {
    val file = writeTempFlagstat(sampleFlagstat)
    val fs = FlagstatParser.parse(file).toOption.get
    assertEquals(fs.primaryReads, 49876543L - 123456L - 78901L)
  }

  test("duplicationRate computes correctly") {
    val file = writeTempFlagstat(sampleFlagstat)
    val fs = FlagstatParser.parse(file).toOption.get
    val rate = fs.duplicationRate.get
    assertEqualsDouble(rate, 2345678.0 / 49876543.0, 0.0001)
  }

  test("parse empty file returns default result") {
    val file = writeTempFlagstat("")
    val result = FlagstatParser.parse(file)
    assert(result.isRight)
    assertEquals(result.toOption.get.totalReads, 0L)
  }

  test("parse nonexistent file returns Left") {
    val result = FlagstatParser.parse(new File("/nonexistent/flagstat.txt"))
    assert(result.isLeft)
  }
}
