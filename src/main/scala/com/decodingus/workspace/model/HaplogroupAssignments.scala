package com.decodingus.workspace.model

/**
 * Container for paternal (Y-DNA) and maternal (mtDNA) haplogroup classifications.
 *
 * This stores the consensus/best haplogroup result for the biosample.
 * Multi-run tracking and reconciliation details are stored in HaplogroupReconciliation.
 *
 * Matches global Atmosphere schema: com.decodingus.atmosphere.defs#haplogroupAssignments
 */
case class HaplogroupAssignments(
  yDna: Option[HaplogroupResult] = None,
  mtDna: Option[HaplogroupResult] = None
)
