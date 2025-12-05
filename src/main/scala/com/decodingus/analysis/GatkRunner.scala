package com.decodingus.analysis

import org.broadinstitute.hellbender.Main

import java.io.{ByteArrayOutputStream, PrintStream}
import java.security.Permission
import scala.util.{Try, Success, Failure}

/**
 * Safely executes GATK tools without allowing System.exit() to crash the application.
 * Captures stdout/stderr and converts exit codes to Either results.
 */
object GatkRunner {

  case class GatkResult(exitCode: Int, stdout: String, stderr: String)

  /**
   * Runs a GATK tool with the given arguments.
   * Prevents System.exit() from terminating the JVM and captures output.
   *
   * @param args Command line arguments for GATK (tool name first, then options)
   * @return Either an error message (Left) or success (Right with exit code 0)
   */
  def run(args: Array[String]): Either[String, GatkResult] = {
    val originalOut = System.out
    val originalErr = System.err
    val originalSecurityManager = System.getSecurityManager

    val stdoutCapture = new ByteArrayOutputStream()
    val stderrCapture = new ByteArrayOutputStream()

    try {
      // Install security manager to catch System.exit()
      System.setSecurityManager(new NoExitSecurityManager())

      // Capture stdout/stderr
      System.setOut(new PrintStream(stdoutCapture))
      System.setErr(new PrintStream(stderrCapture))

      val result = Try {
        Main.main(args)
        0 // If we get here, exit code is 0
      }

      val exitCode = result match {
        case Success(code) => code
        case Failure(e: ExitException) => e.exitCode
        case Failure(e) => throw e // Re-throw unexpected exceptions
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
      case e: ExitException =>
        val stdout = stdoutCapture.toString
        val stderr = stderrCapture.toString
        if (e.exitCode == 0) {
          Right(GatkResult(e.exitCode, stdout, stderr))
        } else {
          val toolName = args.headOption.getOrElse("Unknown")
          Left(s"$toolName failed with exit code ${e.exitCode}.\n$stderr")
        }
      case e: Exception =>
        val toolName = args.headOption.getOrElse("Unknown")
        Left(s"$toolName threw an exception: ${e.getMessage}")
    } finally {
      // Restore original streams and security manager
      System.setOut(originalOut)
      System.setErr(originalErr)
      System.setSecurityManager(originalSecurityManager)
    }
  }

  /** Exception thrown when System.exit() is called */
  private class ExitException(val exitCode: Int) extends SecurityException(s"System.exit($exitCode) blocked")

  /** Security manager that prevents System.exit() */
  private class NoExitSecurityManager extends SecurityManager {
    override def checkPermission(perm: Permission): Unit = {
      // Allow everything except exit
    }

    override def checkPermission(perm: Permission, context: Any): Unit = {
      // Allow everything except exit
    }

    override def checkExit(status: Int): Unit = {
      throw new ExitException(status)
    }
  }
}
