package com.decodingus.pds

import com.decodingus.workspace.model.*

/**
 * Validates records against Atmosphere Lexicon required fields before PDS submission.
 *
 * Each record type has its own validation method that returns either a list of validation
 * errors or Unit. Records must pass validation before being synced to PDS.
 */
object PdsSyncValidation:

  /**
   * Validate a Biosample for PDS sync.
   * Required: citizenDid, centerName, atUri
   */
  def validateBiosample(b: Biosample): Either[List[String], Unit] =
    val errors = List.newBuilder[String]
    if b.atUri.isEmpty then errors += "Biosample missing atUri"
    if b.citizenDid.isEmpty then errors += "Biosample missing citizenDid"
    if b.centerName.isEmpty then errors += "Biosample missing centerName"
    buildResult(errors)

  /**
   * Validate a SequenceRun for PDS sync.
   * Required: atUri, biosampleRef, platformName, testType, files (non-empty)
   */
  def validateSequenceRun(sr: SequenceRun): Either[List[String], Unit] =
    val errors = List.newBuilder[String]
    if sr.atUri.isEmpty then errors += "SequenceRun missing atUri"
    if sr.biosampleRef.isEmpty then errors += "SequenceRun missing biosampleRef"
    if sr.platformName.isEmpty then errors += "SequenceRun missing platformName"
    if sr.testType.isEmpty then errors += "SequenceRun missing testType"
    if sr.files.isEmpty then errors += "SequenceRun has no files"
    buildResult(errors)

  /**
   * Validate an Alignment for PDS sync.
   * Required: atUri, sequenceRunRef, referenceBuild, aligner
   */
  def validateAlignment(a: Alignment): Either[List[String], Unit] =
    val errors = List.newBuilder[String]
    if a.atUri.isEmpty then errors += "Alignment missing atUri"
    if a.sequenceRunRef.isEmpty then errors += "Alignment missing sequenceRunRef"
    if a.referenceBuild.isEmpty then errors += "Alignment missing referenceBuild"
    if a.aligner.isEmpty then errors += "Alignment missing aligner"
    buildResult(errors)

  /**
   * Validate a ChipProfile (genotype) for PDS sync.
   * Required: atUri, biosampleRef, testTypeCode, provider
   */
  def validateChipProfile(cp: ChipProfile): Either[List[String], Unit] =
    val errors = List.newBuilder[String]
    if cp.atUri.isEmpty then errors += "ChipProfile missing atUri"
    if cp.biosampleRef.isEmpty then errors += "ChipProfile missing biosampleRef"
    if cp.testTypeCode.isEmpty then errors += "ChipProfile missing testTypeCode"
    if cp.provider.isEmpty then errors += "ChipProfile missing provider"
    buildResult(errors)

  /**
   * Validate a PopulationBreakdown for PDS sync.
   * Required: biosampleRef, analysisMethod, panelType, components (non-empty)
   */
  def validatePopulationBreakdown(pb: PopulationBreakdown): Either[List[String], Unit] =
    val errors = List.newBuilder[String]
    if pb.biosampleRef.isEmpty then errors += "PopulationBreakdown missing biosampleRef"
    if pb.analysisMethod.isEmpty then errors += "PopulationBreakdown missing analysisMethod"
    if pb.panelType.isEmpty then errors += "PopulationBreakdown missing panelType"
    if pb.components.isEmpty then errors += "PopulationBreakdown has no components"
    buildResult(errors)

  /**
   * Validate a HaplogroupReconciliation for PDS sync.
   * Required: biosampleRef, dnaType, status, runCalls (non-empty)
   */
  def validateHaplogroupReconciliation(hr: HaplogroupReconciliation): Either[List[String], Unit] =
    val errors = List.newBuilder[String]
    if hr.biosampleRef.isEmpty then errors += "HaplogroupReconciliation missing biosampleRef"
    if hr.runCalls.isEmpty then errors += "HaplogroupReconciliation has no runCalls"
    buildResult(errors)

  private def buildResult(errors: scala.collection.mutable.Builder[String, List[String]]): Either[List[String], Unit] =
    val result = errors.result()
    if result.isEmpty then Right(()) else Left(result)
