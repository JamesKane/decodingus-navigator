package com.decodingus.workspace

import com.decodingus.workspace.model._
import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser._
import io.circe.syntax._

import java.io.File
import java.nio.file.{Files, Path, Paths, StandardOpenOption}
import scala.util.{Try, Success, Failure}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

trait WorkspaceService {
  def load(): Either[String, Workspace]
  def save(workspace: Workspace): Either[String, Unit]
}

object LiveWorkspaceService extends WorkspaceService {

  private val CONFIG_DIR: Path = Paths.get(System.getProperty("user.home"), ".config", "decodingus-tools")
  private val WORKSPACE_FILE: Path = CONFIG_DIR.resolve("workspace.json")

  // Ensure config directory exists
  private def ensureConfigDir(): Unit = {
    if (!Files.exists(CONFIG_DIR)) {
      Files.createDirectories(CONFIG_DIR)
    }
  }

  // --- Circe Codecs for WorkspaceModels ---
  // Custom LocalDateTime encoder/decoder
  implicit val encodeLocalDateTime: Encoder[LocalDateTime] = Encoder.encodeString.contramap[LocalDateTime](_.format(DateTimeFormatter.ISO_LOCAL_DATE_TIME))
  implicit val decodeLocalDateTime: Decoder[LocalDateTime] = Decoder.decodeString.emap { str =>
    Try(LocalDateTime.parse(str, DateTimeFormatter.ISO_LOCAL_DATE_TIME)).toEither.left.map(t => s"LocalDateTime: $t")
  }

  implicit val fileInfoCodec: Codec[FileInfo] = deriveCodec
  implicit val contigMetricsCodec: Codec[ContigMetrics] = deriveCodec
  implicit val alignmentMetricsCodec: Codec[AlignmentMetrics] = deriveCodec
  implicit val alignmentDataCodec: Codec[AlignmentData] = deriveCodec
  implicit val sequenceDataCodec: Codec[SequenceData] = deriveCodec
  implicit val haplogroupResultCodec: Codec[HaplogroupResult] = deriveCodec
  implicit val haplogroupAssignmentsCodec: Codec[HaplogroupAssignments] = deriveCodec
  implicit val biosampleCodec: Codec[Biosample] = deriveCodec
  implicit val projectCodec: Codec[Project] = deriveCodec
  implicit val workspaceContentCodec: Codec[WorkspaceContent] = deriveCodec // New codec for WorkspaceContent
  implicit val workspaceCodec: Codec[Workspace] = deriveCodec
  // --- End Circe Codecs ---

  /**
   * Loads the Workspace from ~/.config/decodingus-tools/workspace.json.
   * If the file does not exist, returns an empty Workspace.
   *
   * @return Either an error message or the loaded Workspace.
   */
  override def load(): Either[String, Workspace] = {
    println(s"[DEBUG] Attempting to load workspace from: $WORKSPACE_FILE")
    if (!Files.exists(WORKSPACE_FILE)) {
      println(s"[DEBUG] Workspace file not found at $WORKSPACE_FILE. Initializing with empty workspace.")
      Right(Workspace(lexicon = 1, id = "com.decodingus.atmosphere.workspace", main = WorkspaceContent(samples = List.empty, projects = List.empty)))
    } else {
      Try(Files.readString(WORKSPACE_FILE)) match {
        case Success(jsonString) =>
          println(s"[DEBUG] File content read. Length: ${jsonString.length}. Attempting to parse JSON.")
          parse(jsonString).flatMap(_.as[Workspace]) match {
            case Right(workspace) =>
              println(s"[DEBUG] Successfully parsed workspace: ${workspace.main.samples.size} samples, ${workspace.main.projects.size} projects.")
              Right(workspace)
            case Left(error) =>
              println(s"[DEBUG] Failed to parse workspace JSON: ${error.getMessage()}. Content: ${jsonString.take(200)}...")
              Left(s"Failed to parse workspace JSON: ${error.getMessage()}")
          }
        case Failure(exception) =>
          println(s"[DEBUG] Failed to read workspace file: ${exception.getMessage}")
          Left(s"Failed to read workspace file: ${exception.getMessage}")
      }
    }
  }

  /**
   * Saves the given Workspace to ~/.config/decodingus-tools/workspace.json.
   *
   * @param workspace The Workspace object to save.
   * @return Either an error message or Unit on success.
   */
  override def save(workspace: Workspace): Either[String, Unit] = {
    ensureConfigDir()
    val jsonString = workspace.asJson.spaces2
    Try(Files.writeString(WORKSPACE_FILE, jsonString, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)) match {
      case Success(_) =>
        println(s"[DEBUG] Workspace saved to $WORKSPACE_FILE")
        Right(())
      case Failure(exception) => Left(s"Failed to write workspace to file: ${exception.getMessage}")
    }
  }
}