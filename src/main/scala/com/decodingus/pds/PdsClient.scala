package com.decodingus.pds

import com.decodingus.auth.User
import com.decodingus.model.{ContigSummary, CoverageSummary, LibraryStats, WgsMetrics}
import com.decodingus.workspace.model._
import sttp.client3._
import sttp.client3.circe._
import io.circe.{Decoder, Encoder, Json}
import io.circe.generic.semiauto.{deriveDecoder, deriveEncoder}
import io.circe.syntax._
import io.circe.parser.decode

import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import scala.concurrent.{ExecutionContext, Future}
import scala.util.Try

object PdsClient {

  implicit val libraryStatsEncoder: Encoder[LibraryStats] = deriveEncoder
  implicit val wgsMetricsEncoder: Encoder[WgsMetrics] = deriveEncoder
  implicit val contigSummaryEncoder: Encoder[ContigSummary] = deriveEncoder
  implicit val coverageSummaryEncoder: Encoder[CoverageSummary] = deriveEncoder

  // --- Workspace Codecs (mirrored from LiveWorkspaceService for PDS serialization) ---
  implicit val encodeLocalDateTime: Encoder[LocalDateTime] = Encoder.encodeString.contramap[LocalDateTime](_.format(DateTimeFormatter.ISO_LOCAL_DATE_TIME))
  implicit val decodeLocalDateTime: Decoder[LocalDateTime] = Decoder.decodeString.emap { str =>
    Try(LocalDateTime.parse(str, DateTimeFormatter.ISO_LOCAL_DATE_TIME)).toEither.left.map(t => s"LocalDateTime: $t")
  }

  // Basic type codecs (order matters for derivation - leaf types first)
  implicit val fileInfoEncoder: Encoder[FileInfo] = deriveEncoder
  implicit val fileInfoDecoder: Decoder[FileInfo] = deriveDecoder
  implicit val recordMetaEncoder: Encoder[RecordMeta] = deriveEncoder
  implicit val recordMetaDecoder: Decoder[RecordMeta] = deriveDecoder
  implicit val contigMetricsEncoder: Encoder[ContigMetrics] = deriveEncoder
  implicit val contigMetricsDecoder: Decoder[ContigMetrics] = deriveDecoder
  implicit val alignmentMetricsEncoder: Encoder[AlignmentMetrics] = deriveEncoder
  implicit val alignmentMetricsDecoder: Decoder[AlignmentMetrics] = deriveDecoder

  // Variant call codecs (needed before HaplogroupResult)
  implicit val variantCallEncoder: Encoder[VariantCall] = deriveEncoder
  implicit val variantCallDecoder: Decoder[VariantCall] = deriveDecoder
  implicit val privateVariantDataEncoder: Encoder[PrivateVariantData] = deriveEncoder
  implicit val privateVariantDataDecoder: Decoder[PrivateVariantData] = deriveDecoder

  // Haplogroup codecs (after variant call types)
  implicit val haplogroupResultEncoder: Encoder[HaplogroupResult] = deriveEncoder
  implicit val haplogroupResultDecoder: Decoder[HaplogroupResult] = deriveDecoder
  implicit val haplogroupAssignmentsEncoder: Encoder[HaplogroupAssignments] = deriveEncoder
  implicit val haplogroupAssignmentsDecoder: Decoder[HaplogroupAssignments] = deriveDecoder

  // First-class record codecs
  implicit val sequenceRunEncoder: Encoder[SequenceRun] = deriveEncoder
  implicit val sequenceRunDecoder: Decoder[SequenceRun] = deriveDecoder
  implicit val alignmentEncoder: Encoder[Alignment] = deriveEncoder
  implicit val alignmentDecoder: Decoder[Alignment] = deriveDecoder

  // STR profile codecs
  implicit val strAlleleEncoder: Encoder[StrAllele] = deriveEncoder
  implicit val strAlleleDecoder: Decoder[StrAllele] = deriveDecoder

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
        Left(io.circe.DecodingFailure(s"Unknown StrValue type: $other", cursor.history))
    }
  }

  implicit val strMarkerValueEncoder: Encoder[StrMarkerValue] = deriveEncoder
  implicit val strMarkerValueDecoder: Decoder[StrMarkerValue] = deriveDecoder
  implicit val strPanelEncoder: Encoder[StrPanel] = deriveEncoder
  implicit val strPanelDecoder: Decoder[StrPanel] = deriveDecoder
  implicit val strProfileEncoder: Encoder[StrProfile] = deriveEncoder
  implicit val strProfileDecoder: Decoder[StrProfile] = deriveDecoder

  // Chip profile codecs
  implicit val chipProfileEncoder: Encoder[ChipProfile] = deriveEncoder
  implicit val chipProfileDecoder: Decoder[ChipProfile] = deriveDecoder

  // Y-DNA SNP panel codecs
  import com.decodingus.genotype.model.{YSnpResult, YDnaSnpCall, YDnaSnpPanelResult}
  implicit val ySnpResultEncoder: Encoder[YSnpResult] = Encoder.encodeString.contramap(_.toString)
  implicit val ySnpResultDecoder: Decoder[YSnpResult] = Decoder.decodeString.emap {
    case "Positive" => Right(YSnpResult.Positive)
    case "Negative" => Right(YSnpResult.Negative)
    case "NoCall" => Right(YSnpResult.NoCall)
    case other => Left(s"Unknown YSnpResult: $other")
  }
  implicit val yDnaSnpCallEncoder: Encoder[YDnaSnpCall] = deriveEncoder
  implicit val yDnaSnpCallDecoder: Decoder[YDnaSnpCall] = deriveDecoder
  implicit val yDnaSnpPanelResultEncoder: Encoder[YDnaSnpPanelResult] = deriveEncoder
  implicit val yDnaSnpPanelResultDecoder: Decoder[YDnaSnpPanelResult] = deriveDecoder

  // HaplogroupReconciliation enum codecs
  import com.decodingus.workspace.model.{DnaType, CompatibilityLevel, HaplogroupTechnology, CallMethod, ConflictResolution}
  import com.decodingus.workspace.model.{RunHaplogroupCall, SnpCallFromRun, SnpConflict, ReconciliationStatus, HaplogroupReconciliation}

  implicit val dnaTypeEncoder: Encoder[DnaType] = Encoder.encodeString.contramap(_.toString)
  implicit val dnaTypeDecoder: Decoder[DnaType] = Decoder.decodeString.emap {
    case "Y_DNA" => Right(DnaType.Y_DNA)
    case "MT_DNA" => Right(DnaType.MT_DNA)
    case other => Left(s"Unknown DnaType: $other")
  }

  implicit val compatibilityLevelEncoder: Encoder[CompatibilityLevel] = Encoder.encodeString.contramap(_.toString)
  implicit val compatibilityLevelDecoder: Decoder[CompatibilityLevel] = Decoder.decodeString.emap {
    case "COMPATIBLE" => Right(CompatibilityLevel.COMPATIBLE)
    case "MINOR_DIVERGENCE" => Right(CompatibilityLevel.MINOR_DIVERGENCE)
    case "MAJOR_DIVERGENCE" => Right(CompatibilityLevel.MAJOR_DIVERGENCE)
    case "INCOMPATIBLE" => Right(CompatibilityLevel.INCOMPATIBLE)
    case other => Left(s"Unknown CompatibilityLevel: $other")
  }

  implicit val haplogroupTechnologyEncoder: Encoder[HaplogroupTechnology] = Encoder.encodeString.contramap(_.toString)
  implicit val haplogroupTechnologyDecoder: Decoder[HaplogroupTechnology] = Decoder.decodeString.emap {
    case "WGS" => Right(HaplogroupTechnology.WGS)
    case "WES" => Right(HaplogroupTechnology.WES)
    case "BIG_Y" => Right(HaplogroupTechnology.BIG_Y)
    case "SNP_ARRAY" => Right(HaplogroupTechnology.SNP_ARRAY)
    case "AMPLICON" => Right(HaplogroupTechnology.AMPLICON)
    case "STR_PANEL" => Right(HaplogroupTechnology.STR_PANEL)
    case other => Left(s"Unknown HaplogroupTechnology: $other")
  }

  implicit val callMethodEncoder: Encoder[CallMethod] = Encoder.encodeString.contramap(_.toString)
  implicit val callMethodDecoder: Decoder[CallMethod] = Decoder.decodeString.emap {
    case "SNP_PHYLOGENETIC" => Right(CallMethod.SNP_PHYLOGENETIC)
    case "STR_PREDICTION" => Right(CallMethod.STR_PREDICTION)
    case "VENDOR_REPORTED" => Right(CallMethod.VENDOR_REPORTED)
    case other => Left(s"Unknown CallMethod: $other")
  }

  implicit val conflictResolutionEncoder: Encoder[ConflictResolution] = Encoder.encodeString.contramap(_.toString)
  implicit val conflictResolutionDecoder: Decoder[ConflictResolution] = Decoder.decodeString.emap {
    case "ACCEPT_MAJORITY" => Right(ConflictResolution.ACCEPT_MAJORITY)
    case "ACCEPT_HIGHER_QUALITY" => Right(ConflictResolution.ACCEPT_HIGHER_QUALITY)
    case "ACCEPT_HIGHER_COVERAGE" => Right(ConflictResolution.ACCEPT_HIGHER_COVERAGE)
    case "UNRESOLVED" => Right(ConflictResolution.UNRESOLVED)
    case "HETEROPLASMY" => Right(ConflictResolution.HETEROPLASMY)
    case other => Left(s"Unknown ConflictResolution: $other")
  }

  implicit val runHaplogroupCallEncoder: Encoder[RunHaplogroupCall] = deriveEncoder
  implicit val runHaplogroupCallDecoder: Decoder[RunHaplogroupCall] = deriveDecoder
  implicit val snpCallFromRunEncoder: Encoder[SnpCallFromRun] = deriveEncoder
  implicit val snpCallFromRunDecoder: Decoder[SnpCallFromRun] = deriveDecoder
  implicit val snpConflictEncoder: Encoder[SnpConflict] = deriveEncoder
  implicit val snpConflictDecoder: Decoder[SnpConflict] = deriveDecoder
  implicit val reconciliationStatusEncoder: Encoder[ReconciliationStatus] = deriveEncoder
  implicit val reconciliationStatusDecoder: Decoder[ReconciliationStatus] = deriveDecoder
  implicit val haplogroupReconciliationEncoder: Encoder[HaplogroupReconciliation] = deriveEncoder
  implicit val haplogroupReconciliationDecoder: Decoder[HaplogroupReconciliation] = deriveDecoder

  // Main entity codecs
  implicit val biosampleEncoder: Encoder[Biosample] = deriveEncoder
  implicit val biosampleDecoder: Decoder[Biosample] = deriveDecoder
  implicit val projectEncoder: Encoder[Project] = deriveEncoder
  implicit val projectDecoder: Decoder[Project] = deriveDecoder
  implicit val workspaceContentEncoder: Encoder[WorkspaceContent] = deriveEncoder
  implicit val workspaceContentDecoder: Decoder[WorkspaceContent] = deriveDecoder
  implicit val workspaceEncoder: Encoder[Workspace] = deriveEncoder
  implicit val workspaceDecoder: Decoder[Workspace] = deriveDecoder

  // Using the Future backend
  private val backend = HttpClientFutureBackend()

  // Placeholder endpoint - would normally be in config
  // In a real scenario, this URL would handle the ingestion of the JSON payload.
  private val PdsEndpoint = uri"https://decoding.us.com/api/v1/pds/ingest"

  /**
   * Legacy Upload: Transmits the summary data to the user's PDS data vault via the DecodingUs ingestion API.
   *
   * @param summary The CoverageSummary to upload.
   * @param ec      The execution context for the future.
   * @return A Future that completes when the upload is finished.
   */
  def uploadSummary(summary: CoverageSummary)(implicit ec: ExecutionContext): Future[Unit] = {
    val request = basicRequest
      .post(PdsEndpoint)
      .body(summary)
      .response(asString)

    println(s"Initiating PDS upload for user ${summary.pdsUserId} to $PdsEndpoint...")

    request.send(backend).flatMap { response =>
      response.body match {
        case Right(body) =>
          println(s"Successfully uploaded summary. Server response: $body")
          Future.successful(())
        case Left(error) =>
          println(s"Failed to upload summary. Status: ${response.code}, Error: $error")
          Future.failed(new RuntimeException(s"PDS Upload failed with status ${response.code}: $error"))
      }
    }.recoverWith { case e: Exception =>
      println(s"Exception during PDS upload: ${e.getMessage}")
      Future.failed(e)
    }
  }

  /**
   * AT Protocol Upload: Transmits the summary data to the user's PDS data vault via the AT Protocol.
   * Uses com.atproto.repo.createRecord to store the data.
   *
   * @param user    The authenticated User.
   * @param summary The CoverageSummary to upload.
   * @param ec      The execution context for the future.
   * @return        A Future that completes when the upload is finished.
   */
  def uploadSummaryAtProto(user: User, summary: CoverageSummary)(implicit ec: ExecutionContext): Future[Unit] = {
    // Check if PDS URL is valid (it should be if user is logged in)
    val pdsUrl = if (user.pdsUrl.endsWith("/")) user.pdsUrl.dropRight(1) else user.pdsUrl
    val endpoint = uri"$pdsUrl/xrpc/com.atproto.repo.createRecord"
    
    val collection = "com.decodingus.genome.summary"
    // Generate a random Record Key (rkey)
    val rkey = java.util.UUID.randomUUID().toString

    val payload = Json.obj(
      "repo" -> Json.fromString(user.did),
      "collection" -> Json.fromString(collection),
      "rkey" -> Json.fromString(rkey),
      "record" -> summary.asJson
    )

    val request = basicRequest
      .post(endpoint)
      .header("Authorization", s"Bearer ${user.token}")
      .body(payload)
      .response(asString)

    println(s"Initiating PDS upload for user ${user.did} to $endpoint...")

    request.send(backend).flatMap { response =>
      response.body match {
        case Right(body) =>
          println(s"Successfully uploaded summary to PDS. Response: $body")
          Future.successful(())
        case Left(error) =>
          println(s"Failed to upload summary. Status: ${response.code}, Error: $error")
          Future.failed(new RuntimeException(s"PDS Upload failed with status ${response.code}: $error"))
      }
    }.recoverWith { case e: Exception =>
      println(s"Exception during PDS upload: ${e.getMessage}")
      Future.failed(e)
    }
  }

  // --- Workspace AT Protocol Methods ---

  private val WorkspaceCollection = "com.decodingus.atmosphere.workspace"
  // Use a fixed rkey for the singleton workspace record
  private val WorkspaceRkey = "self"

  /**
   * Saves the workspace to the user's PDS via AT Protocol.
   * Uses putRecord to create or update the singleton workspace record.
   *
   * @param user      The authenticated User.
   * @param workspace The Workspace to save.
   * @param ec        The execution context for the future.
   * @return          A Future that completes when the save is finished.
   */
  def saveWorkspace(user: User, workspace: Workspace)(implicit ec: ExecutionContext): Future[Unit] = {
    val pdsUrl = if (user.pdsUrl.endsWith("/")) user.pdsUrl.dropRight(1) else user.pdsUrl
    val endpoint = uri"$pdsUrl/xrpc/com.atproto.repo.putRecord"

    val payload = Json.obj(
      "repo" -> Json.fromString(user.did),
      "collection" -> Json.fromString(WorkspaceCollection),
      "rkey" -> Json.fromString(WorkspaceRkey),
      "record" -> workspace.asJson
    )

    val request = basicRequest
      .post(endpoint)
      .header("Authorization", s"Bearer ${user.token}")
      .body(payload)
      .response(asString)

    println(s"[PDS] Saving workspace for user ${user.did} to $endpoint...")

    request.send(backend).flatMap { response =>
      response.body match {
        case Right(body) =>
          println(s"[PDS] Successfully saved workspace. Response: $body")
          Future.successful(())
        case Left(error) =>
          println(s"[PDS] Failed to save workspace. Status: ${response.code}, Error: $error")
          Future.failed(new RuntimeException(s"PDS save failed with status ${response.code}: $error"))
      }
    }.recoverWith { case e: Exception =>
      println(s"[PDS] Exception during workspace save: ${e.getMessage}")
      Future.failed(e)
    }
  }

  /**
   * Loads the workspace from the user's PDS via AT Protocol.
   *
   * @param user The authenticated User.
   * @param ec   The execution context for the future.
   * @return     A Future containing the Workspace, or a failed future if not found/error.
   */
  def loadWorkspace(user: User)(implicit ec: ExecutionContext): Future[Workspace] = {
    val pdsUrl = if (user.pdsUrl.endsWith("/")) user.pdsUrl.dropRight(1) else user.pdsUrl
    val endpoint = uri"$pdsUrl/xrpc/com.atproto.repo.getRecord?repo=${user.did}&collection=$WorkspaceCollection&rkey=$WorkspaceRkey"

    val request = basicRequest
      .get(endpoint)
      .header("Authorization", s"Bearer ${user.token}")
      .response(asString)

    println(s"[PDS] Loading workspace for user ${user.did} from $endpoint...")

    request.send(backend).flatMap { response =>
      response.body match {
        case Right(body) =>
          // The response contains { uri, cid, value } - we need to extract "value"
          io.circe.parser.parse(body).flatMap(_.hcursor.downField("value").as[Workspace]) match {
            case Right(workspace) =>
              println(s"[PDS] Successfully loaded workspace: ${workspace.main.samples.size} samples, ${workspace.main.projects.size} projects")
              Future.successful(workspace)
            case Left(parseError) =>
              println(s"[PDS] Failed to parse workspace response: $parseError")
              Future.failed(new RuntimeException(s"Failed to parse workspace: $parseError"))
          }
        case Left(error) =>
          // 400 with RecordNotFound is expected for new users
          if (response.code.code == 400 && error.contains("RecordNotFound")) {
            println(s"[PDS] No workspace found for user (new user), returning empty workspace")
            Future.successful(Workspace.empty)
          } else {
            println(s"[PDS] Failed to load workspace. Status: ${response.code}, Error: $error")
            Future.failed(new RuntimeException(s"PDS load failed with status ${response.code}: $error"))
          }
      }
    }.recoverWith { case e: Exception =>
      println(s"[PDS] Exception during workspace load: ${e.getMessage}")
      Future.failed(e)
    }
  }
}