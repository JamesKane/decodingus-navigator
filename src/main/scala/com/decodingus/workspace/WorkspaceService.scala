package com.decodingus.workspace

import com.decodingus.workspace.model._
import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser._
import io.circe.syntax._

import java.io.File
import java.nio.file.{Files, Paths, StandardOpenOption}
import scala.util.{Try, Success, Failure}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

trait WorkspaceService {
  def load(filePath: String = "workspace.json"): Either[String, Workspace]
  def save(workspace: Workspace, filePath: String = "workspace.json"): Either[String, Unit]
}

object LiveWorkspaceService extends WorkspaceService {

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
  implicit val workspaceCodec: Codec[Workspace] = deriveCodec
  // --- End Circe Codecs ---

  private val WORKSPACE_FILE_NAME = "workspace.json"

  /**
   * Loads the Workspace from a local JSON file.
   * If the file does not exist, returns an empty Workspace.
   *
   * @param filePath The path to the workspace JSON file.
   * @return Either an error message or the loaded Workspace.
   */
  override def load(filePath: String = WORKSPACE_FILE_NAME): Either[String, Workspace] = {
    val path = Paths.get(filePath)
    if (!Files.exists(path)) {
      println(s"Workspace file not found at $filePath. Initializing with empty workspace.")
      Right(Workspace(samples = List.empty, projects = List.empty))
    } else {
      Try(Files.readString(path)) match {
        case Success(jsonString) =>
          parse(jsonString).flatMap(_.as[Workspace]) match {
            case Right(workspace) => Right(workspace)
            case Left(error) => Left(s"Failed to parse workspace JSON: ${error.getMessage()}")
          }
        case Failure(exception) => Left(s"Failed to read workspace file: ${exception.getMessage}")
      }
    }
  }

  /**
   * Saves the given Workspace to a local JSON file.
   *
   * @param workspace The Workspace object to save.
   * @param filePath The path to save the workspace JSON file.
   * @return Either an error message or Unit on success.
   */
  override def save(workspace: Workspace, filePath: String = WORKSPACE_FILE_NAME): Either[String, Unit] = {
    val jsonString = workspace.asJson.spaces2 // Use spaces2 for pretty printing
    Try(Files.writeString(Paths.get(filePath), jsonString, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)) match {
      case Success(_) => Right(())
      case Failure(exception) => Left(s"Failed to write workspace to file: ${exception.getMessage}")
    }
  }
}