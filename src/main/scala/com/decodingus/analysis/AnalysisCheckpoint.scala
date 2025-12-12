package com.decodingus.analysis

import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser._
import io.circe.syntax._

import java.nio.file.{Files, Path}
import java.time.Instant
import scala.util.{Try, Using}

/**
 * Tracks the completion status of analysis pipeline steps.
 * Allows resuming from the last successful step on re-run.
 *
 * @param wgsMetricsCompleted WGS metrics step completed successfully
 * @param callableLociCompleted Callable loci step completed successfully
 * @param sexInferenceCompleted Sex inference step completed successfully
 * @param mtDnaHaplogroupCompleted mtDNA haplogroup step completed successfully
 * @param yDnaHaplogroupCompleted Y-DNA haplogroup step completed successfully (or skipped if female)
 * @param yDnaSkipped Y-DNA was skipped (female sample)
 * @param ancestryCompleted Ancestry analysis step completed successfully
 * @param lastUpdated Timestamp of last checkpoint update
 * @param bamPath Path to the BAM file being analyzed (for validation)
 * @param bamModified BAM file modification time (for invalidation)
 */
case class AnalysisCheckpoint(
  wgsMetricsCompleted: Boolean = false,
  callableLociCompleted: Boolean = false,
  sexInferenceCompleted: Boolean = false,
  mtDnaHaplogroupCompleted: Boolean = false,
  yDnaHaplogroupCompleted: Boolean = false,
  yDnaSkipped: Boolean = false,
  ancestryCompleted: Boolean = false,
  lastUpdated: Instant = Instant.now(),
  bamPath: Option[String] = None,
  bamModified: Option[Long] = None
) {
  /** Check if all steps are complete */
  def isComplete: Boolean =
    wgsMetricsCompleted &&
    callableLociCompleted &&
    sexInferenceCompleted &&
    mtDnaHaplogroupCompleted &&
    (yDnaHaplogroupCompleted || yDnaSkipped) &&
    ancestryCompleted

  /** Get the next incomplete step number (1-6) or None if all complete */
  def nextStep: Option[Int] = {
    if (!wgsMetricsCompleted) Some(1)
    else if (!callableLociCompleted) Some(2)
    else if (!sexInferenceCompleted) Some(3)
    else if (!mtDnaHaplogroupCompleted) Some(4)
    else if (!yDnaHaplogroupCompleted && !yDnaSkipped) Some(5)
    else if (!ancestryCompleted) Some(6)
    else None
  }

  /** Get count of completed steps */
  def completedSteps: Int = {
    var count = 0
    if (wgsMetricsCompleted) count += 1
    if (callableLociCompleted) count += 1
    if (sexInferenceCompleted) count += 1
    if (mtDnaHaplogroupCompleted) count += 1
    if (yDnaHaplogroupCompleted || yDnaSkipped) count += 1
    if (ancestryCompleted) count += 1
    count
  }

  /** Check if checkpoint is valid for the given BAM file */
  def isValidFor(bamFilePath: String): Boolean = {
    bamPath.contains(bamFilePath) && {
      Try {
        val currentModified = Files.getLastModifiedTime(Path.of(bamFilePath)).toMillis
        bamModified.contains(currentModified)
      }.getOrElse(false)
    }
  }
}

object AnalysisCheckpoint {
  private val CHECKPOINT_FILENAME = "analysis_checkpoint.json"

  // Circe codecs
  implicit val instantEncoder: Encoder[Instant] = Encoder.encodeString.contramap(_.toString)
  implicit val instantDecoder: Decoder[Instant] = Decoder.decodeString.emap { str =>
    Try(Instant.parse(str)).toEither.left.map(_.getMessage)
  }
  implicit val checkpointEncoder: Encoder[AnalysisCheckpoint] = deriveEncoder
  implicit val checkpointDecoder: Decoder[AnalysisCheckpoint] = deriveDecoder

  /**
   * Load checkpoint from artifact directory.
   * Returns empty checkpoint if file doesn't exist or is invalid.
   */
  def load(artifactDir: Path): AnalysisCheckpoint = {
    val checkpointFile = artifactDir.resolve(CHECKPOINT_FILENAME)
    if (Files.exists(checkpointFile)) {
      Try {
        val content = Files.readString(checkpointFile)
        decode[AnalysisCheckpoint](content).getOrElse(AnalysisCheckpoint())
      }.getOrElse(AnalysisCheckpoint())
    } else {
      AnalysisCheckpoint()
    }
  }

  /**
   * Load and validate checkpoint for a specific BAM file.
   * Returns empty checkpoint if invalid or BAM has changed.
   */
  def loadAndValidate(artifactDir: Path, bamPath: String): AnalysisCheckpoint = {
    val checkpoint = load(artifactDir)
    if (checkpoint.isValidFor(bamPath)) {
      println(s"[AnalysisCheckpoint] Resuming from step ${checkpoint.nextStep.getOrElse("complete")} (${checkpoint.completedSteps}/6 steps done)")
      checkpoint
    } else {
      if (checkpoint.bamPath.isDefined) {
        println(s"[AnalysisCheckpoint] BAM file changed, starting fresh analysis")
      }
      // Create new checkpoint with BAM info
      val bamModified = Try(Files.getLastModifiedTime(Path.of(bamPath)).toMillis).toOption
      AnalysisCheckpoint(bamPath = Some(bamPath), bamModified = bamModified)
    }
  }

  /**
   * Save checkpoint to artifact directory.
   */
  def save(artifactDir: Path, checkpoint: AnalysisCheckpoint): Unit = {
    Files.createDirectories(artifactDir)
    val checkpointFile = artifactDir.resolve(CHECKPOINT_FILENAME)
    val updated = checkpoint.copy(lastUpdated = Instant.now())
    Files.writeString(checkpointFile, updated.asJson.spaces2)
  }

  /**
   * Update a specific step as completed and save.
   */
  def markStepComplete(artifactDir: Path, checkpoint: AnalysisCheckpoint, step: Int, skipped: Boolean = false): AnalysisCheckpoint = {
    val updated = step match {
      case 1 => checkpoint.copy(wgsMetricsCompleted = true)
      case 2 => checkpoint.copy(callableLociCompleted = true)
      case 3 => checkpoint.copy(sexInferenceCompleted = true)
      case 4 => checkpoint.copy(mtDnaHaplogroupCompleted = true)
      case 5 => if (skipped) checkpoint.copy(yDnaSkipped = true) else checkpoint.copy(yDnaHaplogroupCompleted = true)
      case 6 => checkpoint.copy(ancestryCompleted = true)
      case _ => checkpoint
    }
    save(artifactDir, updated)
    updated
  }

  /**
   * Clear checkpoint (start fresh).
   */
  def clear(artifactDir: Path): Unit = {
    val checkpointFile = artifactDir.resolve(CHECKPOINT_FILENAME)
    Files.deleteIfExists(checkpointFile)
  }
}
