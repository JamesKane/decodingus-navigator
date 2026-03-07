package com.decodingus.ibd.engine

import com.decodingus.workspace.model.RelationshipEstimate

/**
 * Estimates relationship category from total shared centiMorgans.
 *
 * Thresholds based on the Shared cM Project (Blaine Bettinger)
 * and ISOGG wiki reference values. These are population averages;
 * individual values vary significantly.
 */
object RelationshipEstimator:

  /**
   * Estimate relationship category from total shared cM.
   *
   * @param totalSharedCm Total cM shared across all segments
   * @return Relationship estimate
   */
  def estimate(totalSharedCm: Double): RelationshipEstimate =
    if totalSharedCm >= 3400 then RelationshipEstimate.ParentChild
    else if totalSharedCm >= 2550 then RelationshipEstimate.FullSibling
    else if totalSharedCm >= 1700 then RelationshipEstimate.Grandparent
    else if totalSharedCm >= 1200 then RelationshipEstimate.AuntUncle
    else if totalSharedCm >= 680 then RelationshipEstimate.FirstCousin
    else if totalSharedCm >= 400 then RelationshipEstimate.FirstCousinOnceRemoved
    else if totalSharedCm >= 200 then RelationshipEstimate.SecondCousin
    else if totalSharedCm >= 90 then RelationshipEstimate.SecondCousinOnceRemoved
    else if totalSharedCm >= 50 then RelationshipEstimate.ThirdCousin
    else if totalSharedCm >= 25 then RelationshipEstimate.FourthCousin
    else if totalSharedCm >= 10 then RelationshipEstimate.FifthCousin
    else if totalSharedCm >= 7 then RelationshipEstimate.Distant
    else RelationshipEstimate.Unknown
