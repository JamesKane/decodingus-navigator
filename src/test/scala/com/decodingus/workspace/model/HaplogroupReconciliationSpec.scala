package com.decodingus.workspace.model

import munit.FunSuite

class HaplogroupReconciliationSpec extends FunSuite:

  // ============================================
  // LCA / Branch Compatibility
  // ============================================

  test("longestCommonPrefix with identical paths") {
    val paths = List(
      List("R", "R-M269", "R-L21", "R-DF13"),
      List("R", "R-M269", "R-L21", "R-DF13")
    )
    val lca = HaplogroupReconciliation.longestCommonPrefix(paths)
    assertEquals(lca, List("R", "R-M269", "R-L21", "R-DF13"))
  }

  test("longestCommonPrefix with ancestor/descendant paths") {
    val paths = List(
      List("R", "R-M269", "R-L21", "R-DF13", "R-FGC11134"),
      List("R", "R-M269", "R-L21")
    )
    val lca = HaplogroupReconciliation.longestCommonPrefix(paths)
    assertEquals(lca, List("R", "R-M269", "R-L21"))
  }

  test("longestCommonPrefix with divergent paths") {
    val paths = List(
      List("R", "R-M269", "R-L21", "R-DF13", "R-FGC11134"),
      List("R", "R-M269", "R-L21", "R-L1065", "R-S5668")
    )
    val lca = HaplogroupReconciliation.longestCommonPrefix(paths)
    assertEquals(lca, List("R", "R-M269", "R-L21"))
  }

  test("longestCommonPrefix with completely different paths") {
    val paths = List(
      List("R", "R-M269"),
      List("I", "I-M253")
    )
    val lca = HaplogroupReconciliation.longestCommonPrefix(paths)
    assertEquals(lca, List.empty)
  }

  test("longestCommonPrefix with three paths") {
    val paths = List(
      List("R", "R-M269", "R-L21", "R-DF13", "R-FGC11134"),
      List("R", "R-M269", "R-L21", "R-DF13"),
      List("R", "R-M269")
    )
    val lca = HaplogroupReconciliation.longestCommonPrefix(paths)
    assertEquals(lca, List("R", "R-M269"))
  }

  test("scoreToCompatibilityLevel thresholds") {
    assertEquals(HaplogroupReconciliation.scoreToCompatibilityLevel(1.0), CompatibilityLevel.COMPATIBLE)
    assertEquals(HaplogroupReconciliation.scoreToCompatibilityLevel(0.8), CompatibilityLevel.COMPATIBLE)
    assertEquals(HaplogroupReconciliation.scoreToCompatibilityLevel(0.79), CompatibilityLevel.MINOR_DIVERGENCE)
    assertEquals(HaplogroupReconciliation.scoreToCompatibilityLevel(0.5), CompatibilityLevel.MINOR_DIVERGENCE)
    assertEquals(HaplogroupReconciliation.scoreToCompatibilityLevel(0.49), CompatibilityLevel.MAJOR_DIVERGENCE)
    assertEquals(HaplogroupReconciliation.scoreToCompatibilityLevel(0.3), CompatibilityLevel.MAJOR_DIVERGENCE)
    assertEquals(HaplogroupReconciliation.scoreToCompatibilityLevel(0.29), CompatibilityLevel.INCOMPATIBLE)
    assertEquals(HaplogroupReconciliation.scoreToCompatibilityLevel(0.0), CompatibilityLevel.INCOMPATIBLE)
  }

  // ============================================
  // assessBranchCompatibility
  // ============================================

  test("single call is always compatible") {
    val call = makeCall("R-FGC11134", lineage = Some(List("R", "R-M269", "R-L21", "R-DF13", "R-FGC11134")))
    val (level, score, divergence, warnings) =
      HaplogroupReconciliation.assessBranchCompatibility(List(call), DnaType.Y_DNA)
    assertEquals(level, CompatibilityLevel.COMPATIBLE)
    assertEquals(score, Some(1.0))
    assertEquals(divergence, None)
    assert(warnings.isEmpty)
  }

  test("compatible calls on same branch at different depths") {
    val calls = List(
      makeCall("R-FGC11134", lineage = Some(List("R", "R-M269", "R-L21", "R-DF13", "R-FGC11134"))),
      makeCall("R-M269", lineage = Some(List("R", "R-M269")))
    )
    val (level, score, divergence, warnings) =
      HaplogroupReconciliation.assessBranchCompatibility(calls, DnaType.Y_DNA)
    // One path is a prefix of the other — ancestor/descendant, fully compatible
    assertEquals(level, CompatibilityLevel.COMPATIBLE)
    assertEquals(score, Some(1.0))
    assertEquals(divergence, None)
    assert(warnings.isEmpty)
  }

  test("divergent sibling branches detected") {
    val calls = List(
      makeCall("R-FGC11134", lineage = Some(List("R", "R-M269", "R-L21", "R-DF13", "R-FGC11134"))),
      makeCall("R-S5668", lineage = Some(List("R", "R-M269", "R-L21", "R-L1065", "R-S5668")))
    )
    val (level, score, divergence, warnings) =
      HaplogroupReconciliation.assessBranchCompatibility(calls, DnaType.Y_DNA)
    // LCA = ["R", "R-M269", "R-L21"], depth 3, max 5 => 3/5 = 0.6
    assertEquals(level, CompatibilityLevel.MINOR_DIVERGENCE)
    assertEquals(score, Some(0.6))
    assertEquals(divergence, Some("R-L21"))
    assert(warnings.exists(_.contains("diverge")))
  }

  test("incompatible haplogroups from different major branches") {
    val calls = List(
      makeCall("R-M269", lineage = Some(List("R", "R-M269"))),
      makeCall("I-M253", lineage = Some(List("I", "I-M253")))
    )
    val (level, score, divergence, warnings) =
      HaplogroupReconciliation.assessBranchCompatibility(calls, DnaType.Y_DNA)
    // LCA = [], depth 0, max 2 => 0/2 = 0.0
    assertEquals(level, CompatibilityLevel.INCOMPATIBLE)
    assertEquals(score, Some(0.0))
    assert(warnings.exists(_.contains("verification")))
  }

  test("missing lineage data falls back gracefully") {
    val calls = List(
      makeCall("R-M269", lineage = None),
      makeCall("R-L21", lineage = None)
    )
    val (level, _, _, warnings) =
      HaplogroupReconciliation.assessBranchCompatibility(calls, DnaType.Y_DNA)
    assertEquals(level, CompatibilityLevel.COMPATIBLE)
    assert(warnings.exists(_.contains("Insufficient lineage")))
  }

  // ============================================
  // SNP Conflict Detection
  // ============================================

  test("SNP concordance calculated from supporting/conflicting counts") {
    val calls = List(
      makeCall("R-FGC11134", supporting = Some(480), conflicting = Some(5)),
      makeCall("R-FGC11134", supporting = Some(475), conflicting = Some(8))
    )
    val (concordance, _) = HaplogroupReconciliation.detectSnpConflicts(calls, DnaType.Y_DNA)
    assert(concordance.isDefined)
    // (480+475) / (480+475+5+8) = 955/968 ≈ 0.9866
    assert(concordance.get > 0.98)
    assert(concordance.get < 0.99)
  }

  test("no SNP concordance when data missing") {
    val calls = List(
      makeCall("R-FGC11134"),
      makeCall("R-M269")
    )
    val (concordance, conflicts) = HaplogroupReconciliation.detectSnpConflicts(calls, DnaType.Y_DNA)
    assertEquals(concordance, None)
    assert(conflicts.isEmpty)
  }

  // ============================================
  // Full recalculate()
  // ============================================

  test("recalculate with empty calls") {
    val recon = HaplogroupReconciliation(
      meta = RecordMeta.initial,
      biosampleRef = "local:biosample:TEST001",
      dnaType = DnaType.Y_DNA,
      status = ReconciliationStatus(CompatibilityLevel.COMPATIBLE, "R-M269", 0.9, runCount = 1),
      runCalls = List.empty
    )
    val result = recon.recalculate()
    assertEquals(result.status.consensusHaplogroup, "")
    assertEquals(result.status.confidence, 0.0)
    assertEquals(result.status.runCount, 0)
    assert(result.lastReconciliationAt.isDefined)
  }

  test("recalculate picks best call by quality tier") {
    val calls = List(
      makeCall("R-M269", confidence = 0.95, tech = Some(HaplogroupTechnology.SNP_ARRAY)),
      makeCall("R-FGC11134", confidence = 0.85, tech = Some(HaplogroupTechnology.WGS))
    )
    val recon = HaplogroupReconciliation(
      meta = RecordMeta.initial,
      biosampleRef = "local:biosample:TEST001",
      dnaType = DnaType.Y_DNA,
      status = ReconciliationStatus(CompatibilityLevel.COMPATIBLE, "", 0.0, runCount = 0),
      runCalls = calls
    )
    val result = recon.recalculate()
    // WGS (tier 3) should beat SNP_ARRAY (tier 1) even with lower confidence
    assertEquals(result.status.consensusHaplogroup, "R-FGC11134")
  }

  test("recalculate detects compatible ancestor/descendant") {
    val calls = List(
      makeCall("R-FGC11134", confidence = 0.9, tech = Some(HaplogroupTechnology.WGS),
        lineage = Some(List("R", "R-M269", "R-L21", "R-DF13", "R-FGC11134"))),
      makeCall("R-DF13", confidence = 0.8, tech = Some(HaplogroupTechnology.SNP_ARRAY),
        lineage = Some(List("R", "R-M269", "R-L21", "R-DF13")))
    )
    val recon = HaplogroupReconciliation(
      meta = RecordMeta.initial,
      biosampleRef = "local:biosample:TEST001",
      dnaType = DnaType.Y_DNA,
      status = ReconciliationStatus(CompatibilityLevel.COMPATIBLE, "", 0.0, runCount = 0),
      runCalls = calls
    )
    val result = recon.recalculate()
    assertEquals(result.status.consensusHaplogroup, "R-FGC11134")
    // LCA depth 4 / max 5 = 0.8 => COMPATIBLE
    assertEquals(result.status.compatibilityLevel, CompatibilityLevel.COMPATIBLE)
    assert(result.status.branchCompatibilityScore.isDefined)
  }

  test("recalculate detects divergent branches") {
    val calls = List(
      makeCall("R-FGC11134", confidence = 0.9, tech = Some(HaplogroupTechnology.WGS),
        lineage = Some(List("R", "R-M269", "R-L21", "R-DF13", "R-FGC11134"))),
      makeCall("R-S5668", confidence = 0.85, tech = Some(HaplogroupTechnology.BIG_Y),
        lineage = Some(List("R", "R-M269", "R-L21", "R-L1065", "R-S5668")))
    )
    val recon = HaplogroupReconciliation(
      meta = RecordMeta.initial,
      biosampleRef = "local:biosample:TEST001",
      dnaType = DnaType.Y_DNA,
      status = ReconciliationStatus(CompatibilityLevel.COMPATIBLE, "", 0.0, runCount = 0),
      runCalls = calls
    )
    val result = recon.recalculate()
    assertEquals(result.status.compatibilityLevel, CompatibilityLevel.MINOR_DIVERGENCE)
    assertEquals(result.status.divergencePoint, Some("R-L21"))
    assert(result.status.warnings.nonEmpty)
  }

  test("recalculate detects incompatible haplogroups") {
    val calls = List(
      makeCall("R-M269", confidence = 0.9, tech = Some(HaplogroupTechnology.WGS),
        lineage = Some(List("R", "R-M269"))),
      makeCall("I-M253", confidence = 0.85, tech = Some(HaplogroupTechnology.WGS),
        lineage = Some(List("I", "I-M253")))
    )
    val recon = HaplogroupReconciliation(
      meta = RecordMeta.initial,
      biosampleRef = "local:biosample:TEST001",
      dnaType = DnaType.Y_DNA,
      status = ReconciliationStatus(CompatibilityLevel.COMPATIBLE, "", 0.0, runCount = 0),
      runCalls = calls
    )
    val result = recon.recalculate()
    assertEquals(result.status.compatibilityLevel, CompatibilityLevel.INCOMPATIBLE)
    assertEquals(result.status.branchCompatibilityScore, Some(0.0))
    assert(result.status.warnings.exists(_.contains("verification")))
  }

  test("recalculate adds SNP concordance warning when low") {
    val calls = List(
      makeCall("R-FGC11134", supporting = Some(400), conflicting = Some(50),
        lineage = Some(List("R", "R-M269"))),
      makeCall("R-FGC11134", supporting = Some(380), conflicting = Some(60),
        lineage = Some(List("R", "R-M269")))
    )
    val recon = HaplogroupReconciliation(
      meta = RecordMeta.initial,
      biosampleRef = "local:biosample:TEST001",
      dnaType = DnaType.Y_DNA,
      status = ReconciliationStatus(CompatibilityLevel.COMPATIBLE, "", 0.0, runCount = 0),
      runCalls = calls
    )
    val result = recon.recalculate()
    // (400+380) / (400+380+50+60) = 780/890 ≈ 0.876 — below 0.95
    assert(result.status.warnings.exists(_.contains("sample mix-up")))
  }

  // ============================================
  // withRunCall / removeRunCall
  // ============================================

  test("withRunCall replaces existing call from same source") {
    val call1 = makeCall("R-M269", sourceRef = "local:run:1")
    val call2 = makeCall("R-FGC11134", sourceRef = "local:run:1")
    val recon = HaplogroupReconciliation.fromSingleRun("local:biosample:TEST", DnaType.Y_DNA, call1)
    val updated = recon.withRunCall(call2)
    assertEquals(updated.runCalls.size, 1)
    assertEquals(updated.runCalls.head.haplogroup, "R-FGC11134")
  }

  test("removeRunCall removes and empty recalculate clears consensus") {
    val call = makeCall("R-M269", sourceRef = "local:run:1")
    val recon = HaplogroupReconciliation.fromSingleRun("local:biosample:TEST", DnaType.Y_DNA, call)
    val updated = recon.removeRunCall("local:run:1").recalculate()
    assertEquals(updated.runCalls.size, 0)
    assertEquals(updated.status.consensusHaplogroup, "")
    assertEquals(updated.status.runCount, 0)
  }

  // ============================================
  // Helpers
  // ============================================

  private def makeCall(
                        haplogroup: String,
                        confidence: Double = 0.9,
                        tech: Option[HaplogroupTechnology] = None,
                        lineage: Option[List[String]] = None,
                        supporting: Option[Int] = None,
                        conflicting: Option[Int] = None,
                        sourceRef: String = s"local:run:${java.util.UUID.randomUUID()}"
                      ): RunHaplogroupCall =
    RunHaplogroupCall(
      sourceRef = sourceRef,
      haplogroup = haplogroup,
      confidence = confidence,
      callMethod = CallMethod.SNP_PHYLOGENETIC,
      technology = tech,
      lineagePath = lineage,
      supportingSnps = supporting,
      conflictingSnps = conflicting
    )
