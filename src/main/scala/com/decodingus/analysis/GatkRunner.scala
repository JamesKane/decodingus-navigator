package com.decodingus.analysis

import org.broadinstitute.hellbender.Main

import java.io.{ByteArrayOutputStream, PrintStream}

/**
 * Safely executes GATK tools and captures stdout/stderr.
 * Uses GATK's instanceMain() method which returns exit codes instead of calling System.exit().
 */
object GatkRunner {

  case class GatkResult(exitCode: Int, stdout: String, stderr: String)

  /**
   * Runs a GATK tool with the given arguments.
   * Uses instanceMain() to avoid System.exit() calls.
   *
   * @param args Command line arguments for GATK (tool name first, then options)
   * @return Either an error message (Left) or success (Right with exit code 0)
   */
  def run(args: Array[String]): Either[String, GatkResult] = {
    val originalOut = System.out
    val originalErr = System.err

    val stdoutCapture = new ByteArrayOutputStream()
    val stderrCapture = new ByteArrayOutputStream()

    try {
      // Capture stdout/stderr
      System.setOut(new PrintStream(stdoutCapture))
      System.setErr(new PrintStream(stderrCapture))

      // Use instanceMain which returns exit code instead of calling System.exit()
      val gatkMain = new Main()
      val exitCodeObj = gatkMain.instanceMain(args)
      val exitCode = exitCodeObj match {
        case i: java.lang.Integer => i.intValue()
        case _ => 0
      }

      val stdout = stdoutCapture.toString
      val stderr = stderrCapture.toString

      if (exitCode == 0) {
        Right(GatkResult(exitCode, stdout, stderr))
      } else {
        val toolName = args.headOption.getOrElse("Unknown")
        Left(s"$toolName failed with exit code $exitCode.\n$stderr")
      }

    } catch {
      case e: Exception =>
        val stdout = stdoutCapture.toString
        val stderr = stderrCapture.toString
        val toolName = args.headOption.getOrElse("Unknown")
        Left(s"$toolName threw an exception: ${e.getMessage}\n$stderr")
    } finally {
      // Restore original streams
      System.setOut(originalOut)
      System.setErr(originalErr)
    }
  }
}
