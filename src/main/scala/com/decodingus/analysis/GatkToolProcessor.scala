package com.decodingus.analysis

import java.io.File
import java.nio.file.Path

/**
 * Abstract base class for GATK tool processors to reduce code duplication.
 * Encapsulates the common workflow:
 * 1. Ensure BAM index exists
 * 2. Resolve output path (context-aware or temp)
 * 3. Construct arguments
 * 4. Run GATK with progress mapping
 * 5. Parse output
 *
 * @tparam T The type of the result object produced by parsing the output
 */
abstract class GatkToolProcessor[T] {

  /**
   * The name of the GATK tool (for logging/progress messages).
   */
  protected def getToolName: String

  /**
   * Core logic to run a GATK tool.
   * Subclasses should call this from their public process() method.
   *
   * @param bamPath           Path to the BAM/CRAM file
   * @param referencePath     Path to the reference genome
   * @param onProgress        Progress callback
   * @param artifactContext   Optional context for organizing output artifacts
   * @param totalReads        Optional total read count for progress estimation
   * @param buildArgs         Function to build GATK arguments: (bamPath, refPath, outputPath) => args
   * @param parseOutput       Function to parse the output file(s): (outputPath) => Result
   * @param resolveOutputPath Function to resolve output path: (Option[ArtifactContext]) => (Option[OutputDir], OutputPathString)
   * @return Either an error or the result object
   */
  protected def executeGatkTool(
                                 bamPath: String,
                                 referencePath: String,
                                 onProgress: (String, Double, Double) => Unit,
                                 artifactContext: Option[ArtifactContext],
                                 totalReads: Option[Long],
                                 buildArgs: (String, String, String) => Array[String],
                                 parseOutput: String => T,
                                 resolveOutputPath: Option[ArtifactContext] => (Option[Path], String)
                               ): Either[Throwable, T] = {

    onProgress("Checking BAM index...", 0.0, 1.0)
    GatkRunner.ensureIndex(bamPath) match {
      case Left(error) => return Left(new RuntimeException(error))
      case Right(_) => // index exists or was created
    }

    onProgress(s"Running GATK $getToolName...", 0.05, 1.0)

    val (outputDir, outputPath) = resolveOutputPath(artifactContext)

    // Ensure output directory exists if applicable
    outputDir.foreach(dir => dir.toFile.mkdirs())

    val args = buildArgs(bamPath, referencePath, outputPath)

    // Progress callback that maps GATK progress (0-1) to our range (0.05-0.95)
    val gatkProgress: (String, Double) => Unit = (msg, fraction) => {
      val mappedProgress = 0.05 + (fraction * 0.9)
      onProgress(msg, mappedProgress, 1.0)
    }

    GatkRunner.runWithProgress(args, Some(gatkProgress), totalReads, None) match {
      case Right(_) =>
        onProgress(s"Parsing GATK $getToolName output...", 0.95, 1.0)
        try {
          val result = parseOutput(outputPath)
          onProgress(s"GATK $getToolName complete.", 1.0, 1.0)
          Right(result)
        } catch {
          case e: Throwable => Left(e)
        }
      case Left(error) =>
        Left(new RuntimeException(error))
    }
  }
}
