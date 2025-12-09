package com.decodingus.workspace

import com.decodingus.workspace.model._
import io.circe._
import io.circe.generic.semiauto._
import io.circe.parser._
import io.circe.syntax._

import java.nio.file.{Files, Path, Paths, StandardOpenOption}
import scala.util.{Try, Success, Failure}
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

trait WorkspaceService {
  def load(): Either[String, Workspace]
  def save(workspace: Workspace): Either[String, Unit]
}

object LiveWorkspaceService extends WorkspaceService {

  private val CONFIG_DIR: Path = Paths.get(System.getProperty("user.home"), ".decodingus", "config")
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

  // Basic model codecs
  implicit val fileInfoCodec: Codec[FileInfo] = deriveCodec
  implicit val contigMetricsCodec: Codec[ContigMetrics] = deriveCodec
  implicit val alignmentMetricsCodec: Codec[AlignmentMetrics] = deriveCodec
  implicit val recordMetaCodec: Codec[RecordMeta] = deriveCodec

  // Variant call codecs for haplogroup discovery
  implicit val variantCallCodec: Codec[VariantCall] = deriveCodec
  implicit val privateVariantDataCodec: Codec[PrivateVariantData] = deriveCodec

  // Haplogroup codecs
  implicit val haplogroupResultCodec: Codec[HaplogroupResult] = deriveCodec
  implicit val haplogroupAssignmentsCodec: Codec[HaplogroupAssignments] = deriveCodec

  // STR value codecs - discriminated union using "type" field
  implicit val strAlleleCodec: Codec[StrAllele] = deriveCodec

  implicit val encodeStrValue: Encoder[StrValue] = Encoder.instance {
    case SimpleStrValue(repeats) =>
      Json.obj("type" -> Json.fromString("simple"), "repeats" -> Json.fromInt(repeats))
    case MultiCopyStrValue(copies) =>
      Json.obj("type" -> Json.fromString("multiCopy"), "copies" -> Json.arr(copies.map(Json.fromInt): _*))
    case ComplexStrValue(alleles, rawNotation) =>
      val base = Json.obj(
        "type" -> Json.fromString("complex"),
        "alleles" -> alleles.asJson
      )
      rawNotation match {
        case Some(rn) => base.deepMerge(Json.obj("rawNotation" -> Json.fromString(rn)))
        case None => base
      }
  }

  implicit val decodeStrValue: Decoder[StrValue] = Decoder.instance { cursor =>
    cursor.downField("type").as[String].flatMap {
      case "simple" =>
        cursor.downField("repeats").as[Int].map(SimpleStrValue.apply)
      case "multiCopy" =>
        cursor.downField("copies").as[List[Int]].map(MultiCopyStrValue.apply)
      case "complex" =>
        for {
          alleles <- cursor.downField("alleles").as[List[StrAllele]]
          rawNotation <- cursor.downField("rawNotation").as[Option[String]]
        } yield ComplexStrValue(alleles, rawNotation)
      case other =>
        Left(DecodingFailure(s"Unknown StrValue type: $other", cursor.history))
    }
  }

  implicit val strMarkerValueCodec: Codec[StrMarkerValue] = deriveCodec
  implicit val strPanelCodec: Codec[StrPanel] = deriveCodec
  implicit val strProfileCodec: Codec[StrProfile] = deriveCodec

  // Chip profile codec
  implicit val chipProfileCodec: Codec[ChipProfile] = deriveCodec

  // Y-DNA SNP panel codecs
  import com.decodingus.genotype.model.{YSnpResult, YDnaSnpCall, YDnaSnpPanelResult}
  implicit val ySnpResultCodec: Codec[YSnpResult] = Codec.from(
    Decoder.decodeString.emap {
      case "Positive" => Right(YSnpResult.Positive)
      case "Negative" => Right(YSnpResult.Negative)
      case "NoCall" => Right(YSnpResult.NoCall)
      case other => Left(s"Unknown YSnpResult: $other")
    },
    Encoder.encodeString.contramap(_.toString)
  )
  implicit val yDnaSnpCallCodec: Codec[YDnaSnpCall] = deriveCodec
  implicit val yDnaSnpPanelResultCodec: Codec[YDnaSnpPanelResult] = deriveCodec

  // First-class record codecs
  implicit val sequenceRunCodec: Codec[SequenceRun] = deriveCodec
  implicit val alignmentCodec: Codec[Alignment] = deriveCodec
  implicit val biosampleCodec: Codec[Biosample] = deriveCodec
  implicit val projectCodec: Codec[Project] = deriveCodec

  // Workspace codecs
  implicit val workspaceContentCodec: Codec[WorkspaceContent] = deriveCodec
  implicit val workspaceCodec: Codec[Workspace] = deriveCodec
  // --- End Circe Codecs ---

  /**
   * Loads the Workspace from ~/.decodingus/config/workspace.json.
   * If the file does not exist, returns an empty Workspace.
   *
   * @return Either an error message or the loaded Workspace.
   */
  override def load(): Either[String, Workspace] = {
    println(s"[DEBUG] Attempting to load workspace from: $WORKSPACE_FILE")
    if (!Files.exists(WORKSPACE_FILE)) {
      println(s"[DEBUG] Workspace file not found at $WORKSPACE_FILE. Initializing with empty workspace.")
      Right(Workspace.empty)
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
   * Saves the given Workspace to ~/.decodingus/config/workspace.json.
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
