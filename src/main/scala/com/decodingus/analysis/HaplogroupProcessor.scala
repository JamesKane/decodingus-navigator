package com.decodingus.analysis

import com.decodingus.haplogroup.caller.GatkHaplotypeCallerProcessor
import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupResult}
import com.decodingus.haplogroup.scoring.HaplogroupScorer
import com.decodingus.haplogroup.tree.{TreeCache, TreeProviderType, TreeType}

class HaplogroupProcessor {

  def analyze(
    bamPath: String,
    referencePath: String,
    treeType: TreeType,
    providerType: TreeProviderType.Value,
    onProgress: (String, Double, Double) => Unit
  ): Either[String, List[HaplogroupResult]] = {

    onProgress("Loading haplogroup tree...", 0.0, 1.0)
    val treeCache = new TreeCache(treeType, providerType)
    val treeResult = treeCache.getTree

    treeResult.flatMap { haplogroupTree =>
      val rootNodeId = haplogroupTree.allNodes.values.find(_.is_root).map(_.haplogroup_id).getOrElse(0L)
      treeCache.provider.buildTree(haplogroupTree, rootNodeId, treeType) match {
        case None => Left("Failed to build tree structure")
        case Some(rootHaplogroup) =>
          val allLoci = collectAllLoci(rootHaplogroup)

          val caller = new GatkHaplotypeCallerProcessor()
          val snpCalls = caller.callSnps(bamPath, referencePath, allLoci, onProgress)

          val scorer = new HaplogroupScorer()
          val scores = scorer.calculateScores(rootHaplogroup, snpCalls, "GRCh38") // Assuming GRCh38

          // Post-process scores (filter, sort)
          val finalResults = postProcessScores(scores)

          Right(finalResults)
      }
    }
  }

  private def collectAllLoci(haplogroup: Haplogroup): List[com.decodingus.haplogroup.model.Locus] = {
    haplogroup.loci ++ haplogroup.children.flatMap(collectAllLoci)
  }

  private def postProcessScores(scores: List[HaplogroupResult]): List[HaplogroupResult] = {
    // Simple sort by score for now. The logic from the Rust code is more complex.
    scores.sortBy(-_.score)
  }
}