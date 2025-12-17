package com.decodingus.haplogroup.scoring

import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupTree, HaplogroupNode, Locus}
import com.decodingus.haplogroup.tree.TreeType

import scala.jdk.CollectionConverters.*

class HaplogroupScorerSpec extends munit.FunSuite {

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