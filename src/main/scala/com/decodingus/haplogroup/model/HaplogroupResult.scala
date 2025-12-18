package com.decodingus.haplogroup.model

case class HaplogroupScore(
                            matches: Int = 0,
                            ancestralMatches: Int = 0,
                            noCalls: Int = 0,
                            totalSnps: Int = 0,
                            score: Double = 0.0,
                            depth: Int = 0
                          )

case class HaplogroupResult(
                             name: String,
                             score: Double,
                             matchingSnps: Int,
                             mismatchingSnps: Int,
                             ancestralMatches: Int,
                             noCalls: Int,
                             totalSnps: Int,
                             cumulativeSnps: Int,
                             depth: Int,
                             lineagePath: List[String] = List.empty
                           )
