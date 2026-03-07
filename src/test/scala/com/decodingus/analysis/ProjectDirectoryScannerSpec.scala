package com.decodingus.analysis

import munit.FunSuite

import java.io.File
import java.nio.file.Files

class ProjectDirectoryScannerSpec extends FunSuite {

  private def createTempProjectDir(
    samples: Map[String, List[String]]
  ): File = {
    val projectDir = Files.createTempDirectory("PRJEB_test_").toFile
    projectDir.deleteOnExit()

    samples.foreach { case (sampleId, fileNames) =>
      val sampleDir = new File(projectDir, sampleId)
      sampleDir.mkdirs()
      sampleDir.deleteOnExit()

      fileNames.foreach { name =>
        val f = new File(sampleDir, name)
        f.createNewFile()
        f.deleteOnExit()
      }
    }

    projectDir
  }

  test("scan discovers samples with BAM files") {
    val dir = createTempProjectDir(Map(
      "HG02759" -> List("HG02759.GRCh38.cram", "HG02759.GRCh38.cram.crai"),
      "HG03456" -> List("HG03456.GRCh38.bam", "HG03456.GRCh38.bam.bai")
    ))

    val result = ProjectDirectoryScanner.scan(dir)
    assert(result.isRight)

    val project = result.toOption.get
    assertEquals(project.samples.size, 2)
    assertEquals(project.totalAlignmentFiles, 2)
  }

  test("scan classifies file types correctly") {
    val dir = createTempProjectDir(Map(
      "SAMPLE1" -> List(
        "SAMPLE1.cram",
        "SAMPLE1.cram.crai",
        "SAMPLE1.g.vcf.gz",
        "SAMPLE1.flagstat",
        "SAMPLE1.wgs_metrics.txt",
        "readme.txt"
      )
    ))

    val result = ProjectDirectoryScanner.scan(dir)
    val sample = result.toOption.get.samples.head

    assertEquals(sample.alignmentFiles.size, 1)
    assertEquals(sample.variantFiles.size, 1)
    assertEquals(sample.indexFiles.size, 1)
    assertEquals(sample.flagstatFiles.size, 1)
    assertEquals(sample.wgsMetricsFiles.size, 1)
    assert(sample.hasAlignments)
    assert(sample.hasVariants)
    assert(sample.hasPrecomputedMetrics)
  }

  test("scan skips directories without alignment or variant files") {
    val dir = createTempProjectDir(Map(
      "HG02759" -> List("HG02759.cram"),
      "logs" -> List("pipeline.log", "errors.txt")
    ))

    val result = ProjectDirectoryScanner.scan(dir)
    val project = result.toOption.get
    assertEquals(project.samples.size, 1)
    assertEquals(project.samples.head.sampleId, "HG02759")
  }

  test("scan returns error for nonexistent directory") {
    val result = ProjectDirectoryScanner.scan(new File("/nonexistent/dir"))
    assert(result.isLeft)
  }

  test("scan returns error for empty directory") {
    val dir = Files.createTempDirectory("empty_proj_").toFile
    dir.deleteOnExit()

    val result = ProjectDirectoryScanner.scan(dir)
    assert(result.isLeft)
  }

  test("scan returns error for directory with no data files") {
    val dir = createTempProjectDir(Map(
      "logs" -> List("pipeline.log")
    ))

    val result = ProjectDirectoryScanner.scan(dir)
    assert(result.isLeft)
  }

  test("isProjectAccession recognizes ENA accessions") {
    assert(ProjectDirectoryScanner.isProjectAccession("PRJEB31736"))
    assert(ProjectDirectoryScanner.isProjectAccession("PRJNA12345"))
    assert(ProjectDirectoryScanner.isProjectAccession("SRP123456"))
    assert(!ProjectDirectoryScanner.isProjectAccession("HG02759"))
    assert(!ProjectDirectoryScanner.isProjectAccession("random_folder"))
  }

  test("scan uses directory name as projectId") {
    val dir = createTempProjectDir(Map(
      "HG02759" -> List("HG02759.bam")
    ))

    val result = ProjectDirectoryScanner.scan(dir)
    assert(result.toOption.get.projectId.startsWith("PRJEB_test_"))
  }

  test("samplesWithMetrics counts correctly") {
    val dir = createTempProjectDir(Map(
      "S1" -> List("S1.cram", "S1.flagstat"),
      "S2" -> List("S2.cram"),
      "S3" -> List("S3.bam", "S3.wgs_metrics.txt")
    ))

    val project = ProjectDirectoryScanner.scan(dir).toOption.get
    assertEquals(project.samplesWithMetrics, 2)
  }

  test("VCF variants detected with multiple extensions") {
    val dir = createTempProjectDir(Map(
      "S1" -> List("S1.cram", "S1.vcf", "S1.vcf.gz", "S1.g.vcf.gz")
    ))

    val sample = ProjectDirectoryScanner.scan(dir).toOption.get.samples.head
    assertEquals(sample.variantFiles.size, 3)
  }
}
