package com.decodingus.haplogroup.scoring

import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupTree, HaplogroupNode, Locus}
import com.decodingus.haplogroup.vendor.FtdnaTreeProvider
import com.decodingus.haplogroup.tree.TreeType
import htsjdk.variant.vcf.VCFFileReader

import java.io.File
import scala.io.Source
import scala.jdk.CollectionConverters.*
import scala.util.Using

class HaplogroupScorerSpec extends munit.FunSuite {

  /**
   * Parse VCF file with genotypes - same as HaplogroupProcessor.parseVcf
   * Extracts the called allele from the GT field for each position.
   */
  private def parseVcf(vcfPath: String): Map[Long, String] = {
    val vcfFile = new File(vcfPath)
    val reader = new VCFFileReader(vcfFile, false)
    val snpCalls = reader.iterator().asScala.map { vc =>
      val pos = vc.getStart.toLong
      val genotype = vc.getGenotypes.get(0)
      val allele = genotype.getAlleles.get(0).getBaseString
      pos -> allele
    }.toMap
    reader.close()
    snpCalls
  }

  // Load FTDNA tree from JSON
  private def loadFtdnaTree(jsonPath: String): List[Haplogroup] = {
    val data = Using(Source.fromFile(jsonPath))(_.mkString).get
    val treeProvider = new FtdnaTreeProvider(TreeType.MTDNA)
    treeProvider.parseTree(data, "GRCh38") match {
      case Right(tree) => treeProvider.buildTree(tree)
      case Left(error) => throw new RuntimeException(s"Failed to parse tree: $error")
    }
  }

  test("HaplogroupScorer correctly scores MT-DNA haplogroup U5a1b1g sample") {
    // Load test data - use the real GATK-called VCF with genotypes
    val vcfPath = getClass.getResource("/haplogroup/mtdna_calls.vcf").getPath
    val treePath = getClass.getResource("/haplogroup/ftdna-mttree.json").getPath

    val snpCalls = parseVcf(vcfPath)
    val tree = loadFtdnaTree(treePath)

    // Score the haplogroups
    val scorer = new HaplogroupScorer()
    val results = scorer.score(tree, snpCalls)

    // Debug output
    println(s"Total SNP calls in VCF: ${snpCalls.size}")
    println(s"Sample of derived (ALT) calls: ${snpCalls.filter(_._2 != snpCalls.values.head).take(10)}")

    // Check N-defining positions specifically
    val nPositions = List(8701L, 9540L, 10398L, 10873L, 15301L)
    println(s"\nN-defining positions (tree says ancestral->derived):")
    println(s"  8701: G->A, sample has: ${snpCalls.get(8701L)}")
    println(s"  9540: C->T, sample has: ${snpCalls.get(9540L)}")
    println(s"  10398: G->A, sample has: ${snpCalls.get(10398L)}")
    println(s"  10873: C->T, sample has: ${snpCalls.get(10873L)}")
    println(s"  15301: A->G, sample has: ${snpCalls.get(15301L)}")

    // Check U-defining positions
    println(s"\nU-defining positions (tree says ancestral->derived):")
    println(s"  11467: A->G, sample has: ${snpCalls.get(11467L)}")
    println(s"  12308: A->G, sample has: ${snpCalls.get(12308L)}")
    println(s"  12372: G->A, sample has: ${snpCalls.get(12372L)}")
    println(s"\nTop 10 haplogroup results:")
    results.take(10).foreach { r =>
      println(s"  ${r.name}: score=${r.score}, matches=${r.matchingSnps}, ancestral=${r.ancestralMatches}, noCalls=${r.noCalls}, totalSnps=${r.totalSnps}")
    }

    // Find U5a1b1g in results
    val u5a1b1g = results.find(_.name == "U5a1b1g")
    println(s"\nU5a1b1g result: $u5a1b1g")

    // The expected top haplogroup should be U5a1b1g or a parent/child of it
    val topResult = results.head
    println(s"\nTop result: ${topResult.name}")

    // Verify U5a1b1g is in results
    assert(u5a1b1g.isDefined, "U5a1b1g should be in results")

    // Top result should NOT be just "N" - that indicates scoring failure
    assertNotEquals(topResult.name, "N", "Top result should not be 'N' - that indicates scoring failure")

    // The top scoring haplogroup should be in the U5 branch
    assert(
      topResult.name.startsWith("U5") || topResult.name == "U",
      s"Top result should be in U5 branch, got: ${topResult.name}"
    )
  }

  test("HaplogroupScorer handles case-insensitive base comparison") {
    // Create a simple test case with mixed case
    val loci = List(
      Locus("test1", "chrM", 100L, "A", "G"),  // uppercase
      Locus("test2", "chrM", 200L, "C", "t"),  // lowercase derived
      Locus("test3", "chrM", 300L, "g", "A")   // lowercase ref
    )
    val testHaplogroup = Haplogroup("TestHaplo", None, loci, List.empty)

    // SNP calls with opposite case than loci
    val snpCalls = Map(
      100L -> "g",  // lowercase, should match uppercase G in tree
      200L -> "T",  // uppercase, should match lowercase t in tree
      300L -> "a"   // lowercase, should match uppercase A in tree
    )

    val scorer = new HaplogroupScorer()
    val results = scorer.score(List(testHaplogroup), snpCalls)

    val result = results.find(_.name == "TestHaplo").get
    println(s"Case sensitivity test: matches=${result.matchingSnps}, total=${result.totalSnps}")

    // All 3 should match if comparison is case-insensitive
    assertEquals(result.matchingSnps, 3, "All 3 SNPs should match with case-insensitive comparison")
  }
}
