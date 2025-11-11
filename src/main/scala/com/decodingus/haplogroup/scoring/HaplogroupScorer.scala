package com.decodingus.haplogroup.scoring

import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupResult, HaplogroupScore, Locus}

import scala.collection.mutable

class HaplogroupScorer {

  def score(tree: List[Haplogroup], snpCalls: Map[Long, String]): List[HaplogroupResult] = {
    val scores = mutable.ListBuffer[HaplogroupResult]()
    tree.foreach(rootNode => calculateHaplogroupScore(rootNode, snpCalls, scores, None, 0))
    scores.toList.sortBy(-_.score)
  }

  private def calculateHaplogroupScore(
                                        haplogroup: Haplogroup,
                                        snpCalls: Map[Long, String],
                                        scores: mutable.ListBuffer[HaplogroupResult],
                                        parentScore: Option[HaplogroupScore],
                                        depth: Int
                                      ): HaplogroupScore = {
    var currentScore = parentScore.getOrElse(HaplogroupScore())

    var branchDerived = 0
    var branchAncestral = 0
    var branchNoCalls = 0

    for (locus <- haplogroup.loci) {
      snpCalls.get(locus.position) match {
        case Some(calledBase) =>
          if (calledBase == locus.alt) {
            branchDerived += 1
          } else if (calledBase == locus.ref) {
            branchAncestral += 1
          }
        case None =>
          branchNoCalls += 1
      }
    }

    currentScore = currentScore.copy(
      matches = currentScore.matches + branchDerived,
      ancestralMatches = currentScore.ancestralMatches + branchAncestral,
      noCalls = currentScore.noCalls + branchNoCalls,
      totalSnps = currentScore.totalSnps + haplogroup.loci.length
    )

    val scoreValue = (branchDerived + 1).toDouble / (branchAncestral + 1).toDouble
    currentScore = currentScore.copy(score = scoreValue)

    scores += HaplogroupResult(
      name = haplogroup.name,
      score = currentScore.score,
      matchingSnps = currentScore.matches,
      mismatchingSnps = 0, // Mismatch logic not in scoring.rs, seems to be part of low quality
      ancestralMatches = currentScore.ancestralMatches,
      noCalls = currentScore.noCalls,
      totalSnps = haplogroup.loci.length,
      cumulativeSnps = currentScore.totalSnps,
      depth = depth
    )

    for (child <- haplogroup.children) {
      calculateHaplogroupScore(child, snpCalls, scores, Some(currentScore), depth + 1)
    }

    currentScore
  }
}