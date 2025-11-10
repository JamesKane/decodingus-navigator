package com.decodingus.haplogroup.scoring

import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupResult, HaplogroupScore}

import scala.collection.mutable

class HaplogroupScorer {

  def calculateScores(
    haplogroup: Haplogroup,
    snpCalls: Map[Long, String],
    buildId: String
  ): List[HaplogroupResult] = {
    val scores = mutable.ListBuffer[HaplogroupResult]()
    calculateHaplogroupScore(haplogroup, snpCalls, scores, None, 0, buildId)
    scores.toList
  }

  private def calculateHaplogroupScore(
    haplogroup: Haplogroup,
    snpCalls: Map[Long, String],
    scores: mutable.ListBuffer[HaplogroupResult],
    parentScore: Option[HaplogroupScore],
    depth: Int,
    buildId: String
  ): HaplogroupScore = {
    var currentScore = parentScore.getOrElse(HaplogroupScore())

    val definingLoci = haplogroup.loci.filter(_.coordinates.contains(buildId))

    var branchDerived = 0
    var branchAncestral = 0
    var branchNoCalls = 0

    for (locus <- definingLoci) {
      val coord = locus.coordinates(buildId)
      snpCalls.get(coord.position) match {
        case Some(calledBase) =>
          if (calledBase == coord.derived) {
            branchDerived += 1
          } else if (calledBase == coord.ancestral) {
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
      totalSnps = currentScore.totalSnps + definingLoci.length
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
      totalSnps = currentScore.totalSnps,
      cumulativeSnps = currentScore.totalSnps, // Placeholder
      depth = depth
    )

    for (child <- haplogroup.children) {
      calculateHaplogroupScore(child, snpCalls, scores, Some(currentScore), depth + 1, buildId)
    }

    currentScore
  }
}
