package com.decodingus.yprofile.model

import io.circe.*
import io.circe.generic.semiauto.*
import io.circe.syntax.*

/**
 * Circe JSON codecs for Y chromosome profile types.
 */
object YProfileCodecs:

  // ============================================
  // Enum Codecs
  // ============================================

  given Encoder[YProfileSourceType] = Encoder.encodeString.contramap(_.toString)
  given Decoder[YProfileSourceType] = Decoder.decodeString.emap { s =>
    try Right(YProfileSourceType.fromString(s))
    catch case e: IllegalArgumentException => Left(e.getMessage)
  }

  given Encoder[YVariantType] = Encoder.encodeString.contramap(_.toString)
  given Decoder[YVariantType] = Decoder.decodeString.emap { s =>
    try Right(YVariantType.fromString(s))
    catch case e: IllegalArgumentException => Left(e.getMessage)
  }

  given Encoder[YConsensusState] = Encoder.encodeString.contramap(_.toString)
  given Decoder[YConsensusState] = Decoder.decodeString.emap { s =>
    try Right(YConsensusState.fromString(s))
    catch case e: IllegalArgumentException => Left(e.getMessage)
  }

  given Encoder[YVariantStatus] = Encoder.encodeString.contramap(_.toString)
  given Decoder[YVariantStatus] = Decoder.decodeString.emap { s =>
    try Right(YVariantStatus.fromString(s))
    catch case e: IllegalArgumentException => Left(e.getMessage)
  }

  given Encoder[YCallableState] = Encoder.encodeString.contramap(_.toString)
  given Decoder[YCallableState] = Decoder.decodeString.emap { s =>
    try Right(YCallableState.fromString(s))
    catch case e: IllegalArgumentException => Left(e.getMessage)
  }

  given Encoder[YAuditAction] = Encoder.encodeString.contramap(_.toString)
  given Decoder[YAuditAction] = Decoder.decodeString.emap { s =>
    try Right(YAuditAction.fromString(s))
    catch case e: IllegalArgumentException => Left(e.getMessage)
  }

  // ============================================
  // StrMetadata Codec (stored as JSON in database)
  // ============================================

  given Encoder[StrMetadata] = Encoder.instance { meta =>
    Json.obj(
      "repeatMotif" -> meta.repeatMotif.fold(Json.Null)(Json.fromString),
      "repeatUnit" -> meta.repeatUnit.fold(Json.Null)(Json.fromInt),
      "copies" -> meta.copies.fold(Json.Null)(c => Json.arr(c.map(Json.fromInt)*)),
      "rawNotation" -> meta.rawNotation.fold(Json.Null)(Json.fromString)
    ).dropNullValues
  }

  given Decoder[StrMetadata] = Decoder.instance { cursor =>
    for
      repeatMotif <- cursor.get[Option[String]]("repeatMotif")
      repeatUnit <- cursor.get[Option[Int]]("repeatUnit")
      copies <- cursor.get[Option[List[Int]]]("copies")
      rawNotation <- cursor.get[Option[String]]("rawNotation")
    yield StrMetadata(repeatMotif, repeatUnit, copies, rawNotation)
  }
