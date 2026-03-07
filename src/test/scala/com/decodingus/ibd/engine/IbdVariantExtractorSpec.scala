package com.decodingus.ibd.engine

import munit.FunSuite

class IbdVariantExtractorSpec extends FunSuite:

  test("compact bytes round trip") {
    val positions = Array(1000000, 2000000, 3000000, 4000000, 5000000)
    val genotypes = Array[Byte](0, 1, 2, 1, 0)
    val original = ChromosomeGenotypes("1", positions, genotypes)

    val bytes = IbdVariantExtractor.toCompactBytes(original)
    assertEquals(bytes.length, 5 * 5) // 5 entries × 5 bytes each

    val decoded = IbdVariantExtractor.fromCompactBytes(bytes, "1")
    assertEquals(decoded.chromosome, "1")
    assertEquals(decoded.size, 5)

    for i <- 0 until 5 do
      assertEquals(decoded.positions(i), original.positions(i),
        s"Position $i mismatch")
      assertEquals(decoded.genotypes(i), original.genotypes(i),
        s"Genotype $i mismatch")
  }

  test("compact bytes handles large positions") {
    // Test with position > 2^24 (requires all 4 bytes)
    val positions = Array(248_956_422) // chr1 length
    val genotypes = Array[Byte](2)
    val original = ChromosomeGenotypes("1", positions, genotypes)

    val bytes = IbdVariantExtractor.toCompactBytes(original)
    val decoded = IbdVariantExtractor.fromCompactBytes(bytes, "1")
    assertEquals(decoded.positions(0), 248_956_422)
  }

  test("compact bytes handles no-call genotypes") {
    val positions = Array(1000, 2000)
    val genotypes = Array[Byte](-1, 1)
    val original = ChromosomeGenotypes("1", positions, genotypes)

    val bytes = IbdVariantExtractor.toCompactBytes(original)
    val decoded = IbdVariantExtractor.fromCompactBytes(bytes, "1")
    assertEquals(decoded.genotypes(0), (-1).toByte)
    assertEquals(decoded.genotypes(1), 1.toByte)
  }

  test("compact bytes empty data") {
    val original = ChromosomeGenotypes("1", Array.empty[Int], Array.empty[Byte])
    val bytes = IbdVariantExtractor.toCompactBytes(original)
    assertEquals(bytes.length, 0)

    val decoded = IbdVariantExtractor.fromCompactBytes(bytes, "1")
    assertEquals(decoded.size, 0)
  }

  test("ChromosomeGenotypes requires matching array lengths") {
    intercept[IllegalArgumentException] {
      ChromosomeGenotypes("1", Array(1000, 2000), Array[Byte](0))
    }
  }
