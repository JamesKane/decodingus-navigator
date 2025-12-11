package com.decodingus.haplogroup.vendor

import com.decodingus.haplogroup.model.{NamedVariant, VariantAliases, VariantCoordinate, DefiningHaplogroup}
import io.circe.parser.decode

class NamedVariantCacheSpec extends munit.FunSuite {

  test("NamedVariant can be decoded from JSON") {
    val json = """{
      "variantId": 12345,
      "canonicalName": "M269",
      "variantType": "SNP",
      "namingStatus": "named",
      "coordinates": {
        "GRCh38": {
          "contig": "chrY",
          "position": 22127403,
          "ref": "A",
          "alt": "C"
        },
        "GRCh37": {
          "contig": "Y",
          "position": 22755414,
          "ref": "A",
          "alt": "C"
        }
      },
      "aliases": {
        "commonNames": ["M269", "S257"],
        "rsIds": ["rs9786076"],
        "sources": {
          "ISOGG": ["M269"],
          "YFull": ["M269"]
        }
      },
      "definingHaplogroup": {
        "haplogroupId": 1234,
        "haplogroupName": "R-M269"
      }
    }"""

    val result = decode[NamedVariant](json)
    assert(result.isRight, s"Failed to decode: ${result.left.getOrElse("unknown")}")

    val variant = result.toOption.get
    assertEquals(variant.variantId, 12345)
    assertEquals(variant.canonicalName, Some("M269"))
    assertEquals(variant.variantType, "SNP")
    assertEquals(variant.namingStatus, "named")

    // Check coordinates
    assert(variant.coordinates.contains("GRCh38"))
    val grch38Coord = variant.coordinates("GRCh38")
    assertEquals(grch38Coord.contig, "chrY")
    assertEquals(grch38Coord.position, 22127403)
    assertEquals(grch38Coord.ref, "A")
    assertEquals(grch38Coord.alt, "C")

    // Check aliases
    assertEquals(variant.aliases.commonNames, List("M269", "S257"))
    assertEquals(variant.aliases.rsIds, List("rs9786076"))
    assert(variant.aliases.sources.contains("ISOGG"))

    // Check defining haplogroup
    assert(variant.definingHaplogroup.isDefined)
    assertEquals(variant.definingHaplogroup.get.haplogroupName, "R-M269")
  }

  test("NamedVariant can be decoded without optional fields") {
    val json = """{
      "variantId": 99999,
      "variantType": "SNP",
      "namingStatus": "unnamed",
      "coordinates": {
        "GRCh38": {
          "contig": "chrY",
          "position": 12345678,
          "ref": "G",
          "alt": "T"
        }
      },
      "aliases": {
        "commonNames": [],
        "rsIds": [],
        "sources": {}
      }
    }"""

    val result = decode[NamedVariant](json)
    assert(result.isRight, s"Failed to decode: ${result.left.getOrElse("unknown")}")

    val variant = result.toOption.get
    assertEquals(variant.variantId, 99999)
    assertEquals(variant.canonicalName, None)
    assertEquals(variant.definingHaplogroup, None)
    assertEquals(variant.displayName, s"var_99999")
  }

  test("NamedVariant.displayName returns canonical name when present") {
    val variant = NamedVariant(
      variantId = 1,
      canonicalName = Some("M269"),
      variantType = "SNP",
      namingStatus = "named",
      coordinates = Map.empty,
      aliases = VariantAliases(commonNames = List("S257", "U106"))
    )
    assertEquals(variant.displayName, "M269")
  }

  test("NamedVariant.displayName falls back to common name") {
    val variant = NamedVariant(
      variantId = 1,
      canonicalName = None,
      variantType = "SNP",
      namingStatus = "unnamed",
      coordinates = Map.empty,
      aliases = VariantAliases(commonNames = List("S257", "U106"))
    )
    assertEquals(variant.displayName, "S257")
  }

  test("NamedVariant.allNames combines canonical and aliases") {
    val variant = NamedVariant(
      variantId = 1,
      canonicalName = Some("M269"),
      variantType = "SNP",
      namingStatus = "named",
      coordinates = Map.empty,
      aliases = VariantAliases(commonNames = List("S257", "U106"))
    )
    assertEquals(variant.allNames, List("M269", "S257", "U106"))
  }

  test("NamedVariant.coordinateFor returns correct coordinate for build") {
    val variant = NamedVariant(
      variantId = 1,
      canonicalName = Some("M269"),
      variantType = "SNP",
      namingStatus = "named",
      coordinates = Map(
        "GRCh38" -> VariantCoordinate("chrY", 22127403, "A", "C"),
        "GRCh37" -> VariantCoordinate("Y", 22755414, "A", "C")
      ),
      aliases = VariantAliases()
    )

    val grch38 = variant.coordinateFor("GRCh38")
    assert(grch38.isDefined)
    assertEquals(grch38.get.position, 22127403)

    val grch37 = variant.coordinateFor("GRCh37")
    assert(grch37.isDefined)
    assertEquals(grch37.get.position, 22755414)

    val chm13 = variant.coordinateFor("CHM13v2")
    assert(chm13.isEmpty)
  }

  test("NamedVariantCache is singleton") {
    val cache1 = NamedVariantCache()
    val cache2 = NamedVariantCache()
    assert(cache1 eq cache2, "NamedVariantCache should return the same instance")
  }

  test("NamedVariantCache reports needsRefresh when cache file missing") {
    // Create a new cache instance with a unique cache dir
    val cache = new NamedVariantCache()
    // A fresh cache with no file should need refresh
    assert(cache.needsRefresh, "Cache should need refresh when file doesn't exist")
  }
}
