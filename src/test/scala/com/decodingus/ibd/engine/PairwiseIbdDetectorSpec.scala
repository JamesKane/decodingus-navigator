package com.decodingus.ibd.engine

import com.decodingus.workspace.model.IbdSegment
import munit.FunSuite

class PairwiseIbdDetectorSpec extends FunSuite:

  // Uniform rate: 1 cM per Mb
  private val geneticMap = GeneticMap.uniformRate(cmPerMb = 1.0)

  /**
   * Create synthetic genotype data where two individuals share an IBD segment.
   * Outside the segment, genotypes are randomized. Inside, they match (IBS-2).
   */
  private def createIbdPair(
                              chromosome: String,
                              totalSnps: Int,
                              spacing: Int, // bp between SNPs
                              ibdStartSnp: Int,
                              ibdEndSnp: Int,
                              seed: Long = 42L
                            ): (ChromosomeGenotypes, ChromosomeGenotypes) =
    val rng = new java.util.Random(seed)
    val positions = Array.tabulate(totalSnps)(i => (i + 1) * spacing)
    val g1 = Array.tabulate(totalSnps)(i => (rng.nextInt(3)).toByte) // 0, 1, or 2
    val g2 = new Array[Byte](totalSnps)

    // Outside IBD region: independent random genotypes
    val rng2 = new java.util.Random(seed + 1)
    for i <- 0 until totalSnps do
      if i >= ibdStartSnp && i < ibdEndSnp then
        g2(i) = g1(i) // IBD: identical genotypes
      else
        g2(i) = rng2.nextInt(3).toByte // Independent

    (ChromosomeGenotypes(chromosome, positions, g1),
      ChromosomeGenotypes(chromosome, positions, g2))

  test("detect a clear IBD segment") {
    // 2000 SNPs spaced 10kb apart = 20Mb total
    // IBD from SNP 500-1000 = 5Mb = ~5 cM
    // Use relaxed config for this small segment
    val config = IbdDetectorConfig(
      minSegmentCm = 3.0,
      minSnpCount = 50,
      windowSize = 50,
      ibsThreshold = 0.65,
      errorTolerance = 0.02
    )
    val detector = PairwiseIbdDetector(config)

    val (g1, g2) = createIbdPair("1", 2000, 10000, 500, 1000)
    val segments = detector.detectChromosomeSegments(g1, g2, geneticMap)

    assert(segments.nonEmpty, "Should detect at least one segment")
    val seg = segments.maxBy(_.lengthCm)
    assert(seg.lengthCm >= 3.0, s"Segment should be >= 3 cM, got ${seg.lengthCm}")
    // The segment should roughly cover the IBD region (500*10kb to 1000*10kb = 5M-10M bp)
    assert(seg.startPosition >= 3_000_000 && seg.startPosition <= 7_000_000,
      s"Start should be near 5M, got ${seg.startPosition}")
  }

  test("no segments detected for unrelated individuals") {
    val rng1 = new java.util.Random(100L)
    val rng2 = new java.util.Random(200L)
    val n = 1000
    val spacing = 10000
    val positions = Array.tabulate(n)(i => (i + 1) * spacing)
    val g1 = Array.tabulate(n)(_ => rng1.nextInt(3).toByte)
    val g2 = Array.tabulate(n)(_ => rng2.nextInt(3).toByte)

    val config = IbdDetectorConfig(minSegmentCm = 5.0, minSnpCount = 50, windowSize = 50)
    val detector = PairwiseIbdDetector(config)

    val geno1 = ChromosomeGenotypes("1", positions, g1)
    val geno2 = ChromosomeGenotypes("1", positions, g2)

    val segments = detector.detectChromosomeSegments(geno1, geno2, geneticMap)
    // Unrelated individuals should have few or no large IBD segments
    val largeSeg = segments.filter(_.lengthCm >= 5.0)
    assert(largeSeg.isEmpty, s"Unrelated should have no large segments, got ${largeSeg.size}")
  }

  test("identical twins share entire chromosome") {
    val rng = new java.util.Random(42L)
    val n = 2000
    val spacing = 10000 // 10kb spacing → 20Mb total
    val positions = Array.tabulate(n)(i => (i + 1) * spacing)
    val g1 = Array.tabulate(n)(_ => rng.nextInt(3).toByte)
    val g2 = g1.clone() // Identical

    val config = IbdDetectorConfig(minSegmentCm = 5.0, minSnpCount = 50, windowSize = 50)
    val detector = PairwiseIbdDetector(config)

    val geno1 = ChromosomeGenotypes("1", positions, g1)
    val geno2 = ChromosomeGenotypes("1", positions, g2)

    val segments = detector.detectChromosomeSegments(geno1, geno2, geneticMap)
    assert(segments.nonEmpty, "Identical genotypes should produce segments")
    val totalCm = segments.map(_.lengthCm).sum
    // 20Mb at 1cM/Mb = ~20 cM total
    assert(totalCm >= 15.0, s"Expected ~20 cM total, got $totalCm")
  }

  test("detectSegments works across multiple chromosomes") {
    val config = IbdDetectorConfig(
      minSegmentCm = 3.0,
      minSnpCount = 50,
      windowSize = 50,
      ibsThreshold = 0.65,
      errorTolerance = 0.02
    )
    val detector = PairwiseIbdDetector(config)

    val (g1_chr1, g2_chr1) = createIbdPair("1", 1500, 10000, 300, 700)
    val (g1_chr2, g2_chr2) = createIbdPair("2", 1500, 10000, 400, 800, seed = 99L)

    val sample1 = Map("1" -> g1_chr1, "2" -> g1_chr2)
    val sample2 = Map("1" -> g2_chr1, "2" -> g2_chr2)

    val segments = detector.detectSegments(sample1, sample2, geneticMap)
    val chrs = segments.map(_.chromosome).toSet
    assert(chrs.contains("1"), "Should detect segments on chr1")
    assert(chrs.contains("2"), "Should detect segments on chr2")
  }

  test("too few SNPs produces no segments") {
    val config = IbdDetectorConfig(minSnpCount = 100)
    val detector = PairwiseIbdDetector(config)

    val positions = Array(1000, 2000, 3000)
    val g1 = Array[Byte](0, 1, 2)
    val g2 = Array[Byte](0, 1, 2)

    val geno1 = ChromosomeGenotypes("1", positions, g1)
    val geno2 = ChromosomeGenotypes("1", positions, g2)

    val segments = detector.detectChromosomeSegments(geno1, geno2, geneticMap)
    assert(segments.isEmpty, "Should not detect segments with too few SNPs")
  }

  test("no-call genotypes are excluded from comparison") {
    val positions = Array.tabulate(500)(i => (i + 1) * 10000)
    val g1 = Array.tabulate(500)(_ => 1.toByte) // All het
    val g2 = Array.tabulate(500)(i =>
      if i < 250 then (-1).toByte // First half no-call
      else 1.toByte // Second half matches
    )

    val config = IbdDetectorConfig(minSegmentCm = 1.0, minSnpCount = 20, windowSize = 20)
    val detector = PairwiseIbdDetector(config)

    val geno1 = ChromosomeGenotypes("1", positions, g1)
    val geno2 = ChromosomeGenotypes("1", positions, g2)

    val segments = detector.detectChromosomeSegments(geno1, geno2, geneticMap)
    // Only the second half (250 SNPs from ~2.5M to 5M) should contribute
    segments.foreach { seg =>
      assert(seg.startPosition >= 2_000_000,
        s"Segment should start after no-call region, got ${seg.startPosition}")
    }
  }
