package com.decodingus.haplogroup.vendor

import com.decodingus.haplogroup.tree.TreeType
import munit.FunSuite

class DecodingUsTreeProviderSpec extends FunSuite {

  test("correctly select coordinates for the target build (GRCh38)") {
    val provider = new DecodingUsTreeProvider(TreeType.YDNA)
    
    val json =
      """
        |[
        |  {
        |    "name": "R-CTS4466",
        |    "parentName": "R-M269",
        |    "variants": [
        |      {
        |        "name": "CTS4466",
        |        "coordinates": {
        |          "CM000686.1": { "start": 100, "stop": 101, "anc": "A", "der": "G" },
        |          "CM000686.2": { "start": 200, "stop": 201, "anc": "A", "der": "G" }
        |        },
        |        "variantType": "SNP"
        |      }
        |    ],
        |    "lastUpdated": "2023-01-01",
        |    "isBackbone": false
        |  }
        |]
        |""".stripMargin

    // Target Build: GRCh38 (maps to CM000686.2)
    val result = provider.parseTree(json, "GRCh38")

    assert(result.isRight)
    val tree = result.toOption.get
    
    // Find the node
    val node = tree.allNodes.values.find(_.name == "R-CTS4466").get
    
    // Check loci
    assertEquals(node.loci.size, 1)
    val locus = node.loci.head
    
    // Should pick the coordinate for GRCh38 (200), not GRCh37 (100)
    assertEquals(locus.position, 200L)
    assertEquals(locus.contig, "chrY")
  }
  
  test("correctly select coordinates when target build is GRCh37") {
    val provider = new DecodingUsTreeProvider(TreeType.YDNA)

    val json =
      """
        |[
        |  {
        |    "name": "R-CTS4466",
        |    "parentName": "R-M269",
        |    "variants": [
        |      {
        |        "name": "CTS4466",
        |        "coordinates": {
        |          "CM000686.1": { "start": 100, "stop": 101, "anc": "A", "der": "G" },
        |          "CM000686.2": { "start": 200, "stop": 201, "anc": "A", "der": "G" }
        |        },
        |        "variantType": "SNP"
        |      }
        |    ],
        |    "lastUpdated": "2023-01-01",
        |    "isBackbone": false
        |  }
        |]
        |""".stripMargin

    // Target Build: GRCh37 (maps to CM000686.1)
    val result = provider.parseTree(json, "GRCh37")

    assert(result.isRight)
    val tree = result.toOption.get

    val node = tree.allNodes.values.find(_.name == "R-CTS4466").get
    assertEquals(node.loci.size, 1)
    val locus = node.loci.head

    // Should pick the coordinate for GRCh37 (100)
    assertEquals(locus.position, 100L)
    // DecodingUsTreeProvider returns "Y" for GRCh37 Y-DNA
    assertEquals(locus.contig, "Y") 
  }

  test("return no loci if target build is not present in coordinates") {
    val provider = new DecodingUsTreeProvider(TreeType.YDNA)

    val json =
      """
        |[
        |  {
        |    "name": "R-CTS4466",
        |    "parentName": "R-M269",
        |    "variants": [
        |      {
        |        "name": "CTS4466",
        |        "coordinates": {
        |          "CM000686.1": { "start": 100, "stop": 101, "anc": "A", "der": "G" }
        |        },
        |        "variantType": "SNP"
        |      }
        |    ],
        |    "lastUpdated": "2023-01-01",
        |    "isBackbone": false
        |  }
        |]
        |""".stripMargin

    // Target Build: GRCh38 (maps to CM000686.2 - missing)
    val result = provider.parseTree(json, "GRCh38")

    assert(result.isRight)
    val tree = result.toOption.get

    val node = tree.allNodes.values.find(_.name == "R-CTS4466").get
    assertEquals(node.loci.size, 0)
  }
}