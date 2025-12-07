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

  // Legacy embedded data codecs (needed for deprecated fields in Biosample)
  implicit val alignmentDataCodec: Codec[AlignmentData] = deriveCodec
  implicit val sequenceDataCodec: Codec[SequenceData] = deriveCodec

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
   * Loads the Workspace from ~/.config/decodingus-tools/workspace.json.
   * If the file does not exist, returns an empty Workspace.
   * Handles migration from older lexicon versions.
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
          // First try to parse with current schema
          parse(jsonString).flatMap(_.as[Workspace]) match {
            case Right(workspace) =>
              println(s"[DEBUG] Successfully parsed workspace: ${workspace.main.samples.size} samples, ${workspace.main.projects.size} projects.")
              // Check if migration is needed
              if (workspace.lexicon < Workspace.CurrentLexiconVersion) {
                println(s"[DEBUG] Migrating workspace from lexicon ${workspace.lexicon} to ${Workspace.CurrentLexiconVersion}")
                Right(migrateWorkspace(workspace))
              } else {
                Right(workspace)
              }
            case Left(error) =>
              // Try to migrate from legacy format
              println(s"[DEBUG] Failed to parse with current schema, attempting legacy migration: ${error.getMessage()}")
              migrateLegacyWorkspace(jsonString) match {
                case Right(migrated) =>
                  println(s"[DEBUG] Successfully migrated legacy workspace: ${migrated.main.samples.size} samples")
                  Right(migrated)
                case Left(migrationError) =>
                  println(s"[DEBUG] Failed to parse workspace JSON: ${error.getMessage()}. Content: ${jsonString.take(200)}...")
                  Left(s"Failed to parse workspace JSON: ${error.getMessage()}")
              }
          }
        case Failure(exception) =>
          println(s"[DEBUG] Failed to read workspace file: ${exception.getMessage}")
          Left(s"Failed to read workspace file: ${exception.getMessage}")
      }
    }
  }

  /**
   * Migrates a workspace from an older lexicon version to the current version.
   */
  private def migrateWorkspace(workspace: Workspace): Workspace = {
    workspace.copy(lexicon = Workspace.CurrentLexiconVersion)
  }

  /**
   * Attempts to migrate from legacy workspace format (lexicon 1 with embedded SequenceData/AlignmentData).
   */
  private def migrateLegacyWorkspace(jsonString: String): Either[String, Workspace] = {
    // Define legacy types and codecs for parsing old format
    case class LegacyFileInfo(
      fileName: String,
      fileSizeBytes: Option[Long],
      fileFormat: String,
      checksum: Option[String],
      location: String // Legacy format had required String
    )

    case class LegacyAlignmentData(
      referenceBuild: String,
      aligner: String,
      files: List[LegacyFileInfo],
      metrics: Option[AlignmentMetrics]
    )

    case class LegacySequenceData(
      platformName: String,
      instrumentModel: Option[String],
      testType: String,
      libraryLayout: Option[String],
      totalReads: Option[Long],
      readLength: Option[Int],
      meanInsertSize: Option[Double],
      files: List[LegacyFileInfo],
      alignments: List[LegacyAlignmentData]
    )

    case class LegacyBiosample(
      sampleAccession: String,
      donorIdentifier: String,
      atUri: Option[String],
      description: Option[String],
      centerName: Option[String],
      sex: Option[String],
      sequenceData: List[LegacySequenceData],
      haplogroups: Option[HaplogroupAssignments],
      createdAt: Option[LocalDateTime]
    )

    case class LegacyProject(
      projectName: String,
      atUri: Option[String],
      description: Option[String],
      administrator: String,
      members: List[String]
    )

    def convertFileInfo(legacy: LegacyFileInfo): FileInfo = FileInfo(
      fileName = legacy.fileName,
      fileSizeBytes = legacy.fileSizeBytes,
      fileFormat = legacy.fileFormat,
      checksum = legacy.checksum,
      checksumAlgorithm = legacy.checksum.map(_ => "SHA-256"),
      location = Some(legacy.location)
    )

    case class LegacyWorkspaceContent(
      samples: List[LegacyBiosample],
      projects: List[LegacyProject]
    )

    case class LegacyWorkspace(
      lexicon: Int,
      id: String,
      main: LegacyWorkspaceContent
    )

    // Legacy codecs
    implicit val legacyFileInfoCodec: Codec[LegacyFileInfo] = deriveCodec
    implicit val legacyAlignmentDataCodec: Codec[LegacyAlignmentData] = deriveCodec
    implicit val legacySequenceDataCodec: Codec[LegacySequenceData] = deriveCodec
    implicit val legacyBiosampleCodec: Codec[LegacyBiosample] = deriveCodec
    implicit val legacyProjectCodec: Codec[LegacyProject] = deriveCodec
    implicit val legacyWorkspaceContentCodec: Codec[LegacyWorkspaceContent] = deriveCodec
    implicit val legacyWorkspaceCodec: Codec[LegacyWorkspace] = deriveCodec

    parse(jsonString).flatMap(_.as[LegacyWorkspace]) match {
      case Right(legacy) =>
        // Convert legacy format to new format
        var sequenceRuns: List[SequenceRun] = List.empty
        var alignments: List[Alignment] = List.empty

        val samples = legacy.main.samples.zipWithIndex.map { case (legacySample, sampleIdx) =>
          val biosampleUri = legacySample.atUri.getOrElse(s"local:biosample:${legacySample.sampleAccession}")
          val meta = RecordMeta(
            version = 1,
            createdAt = legacySample.createdAt.getOrElse(LocalDateTime.now())
          )

          // Convert embedded sequence data to first-class records
          val sequenceRunRefs = legacySample.sequenceData.zipWithIndex.map { case (legacySeq, seqIdx) =>
            val seqRunUri = s"local:sequencerun:${legacySample.sampleAccession}:$seqIdx"

            // Convert embedded alignments to first-class records
            val alignmentRefs = legacySeq.alignments.zipWithIndex.map { case (legacyAlign, alignIdx) =>
              val alignUri = s"local:alignment:${legacySample.sampleAccession}:$seqIdx:$alignIdx"

              val alignment = Alignment(
                atUri = Some(alignUri),
                meta = meta,
                sequenceRunRef = seqRunUri,
                biosampleRef = Some(biosampleUri),
                referenceBuild = legacyAlign.referenceBuild,
                aligner = legacyAlign.aligner,
                files = legacyAlign.files.map(convertFileInfo),
                metrics = legacyAlign.metrics
              )
              alignments = alignments :+ alignment
              alignUri
            }

            val sequenceRun = SequenceRun(
              atUri = Some(seqRunUri),
              meta = meta,
              biosampleRef = biosampleUri,
              platformName = legacySeq.platformName,
              instrumentModel = legacySeq.instrumentModel,
              testType = legacySeq.testType,
              libraryLayout = legacySeq.libraryLayout,
              totalReads = legacySeq.totalReads,
              readLength = legacySeq.readLength,
              meanInsertSize = legacySeq.meanInsertSize,
              files = legacySeq.files.map(convertFileInfo),
              alignmentRefs = alignmentRefs
            )
            sequenceRuns = sequenceRuns :+ sequenceRun
            seqRunUri
          }

          Biosample(
            atUri = Some(biosampleUri),
            meta = meta,
            sampleAccession = legacySample.sampleAccession,
            donorIdentifier = legacySample.donorIdentifier,
            description = legacySample.description,
            centerName = legacySample.centerName,
            sex = legacySample.sex,
            haplogroups = legacySample.haplogroups,
            sequenceRunRefs = sequenceRunRefs
          )
        }

        val projects = legacy.main.projects.map { legacyProject =>
          Project(
            atUri = legacyProject.atUri,
            meta = RecordMeta.initial,
            projectName = legacyProject.projectName,
            description = legacyProject.description,
            administrator = legacyProject.administrator,
            memberRefs = legacyProject.members
          )
        }

        Right(Workspace(
          lexicon = Workspace.CurrentLexiconVersion,
          id = Workspace.NamespaceId,
          main = WorkspaceContent(
            meta = Some(RecordMeta.initial),
            sampleRefs = samples.flatMap(_.atUri),
            projectRefs = projects.flatMap(_.atUri),
            samples = samples,
            projects = projects,
            sequenceRuns = sequenceRuns,
            alignments = alignments
          )
        ))

      case Left(error) =>
        Left(s"Failed to parse legacy workspace: ${error.getMessage()}")
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