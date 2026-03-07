package com.decodingus.ibd.crypto

import io.circe.*
import java.util.Base64

enum PayloadDataType:
  case Genotypes, Attestation, Summary, KeyExchange

object PayloadDataType:
  def fromString(s: String): PayloadDataType = s match
    case "GENOTYPES" => Genotypes
    case "ATTESTATION" => Attestation
    case "SUMMARY" => Summary
    case "KEY_EXCHANGE" => KeyExchange
    case other => throw IllegalArgumentException(s"Unknown payload data type: $other")

  extension (pdt: PayloadDataType)
    def toDbString: String = pdt match
      case Genotypes => "GENOTYPES"
      case Attestation => "ATTESTATION"
      case Summary => "SUMMARY"
      case KeyExchange => "KEY_EXCHANGE"

/**
 * An encrypted data payload for IBD relay transport.
 *
 * All fields are Base64-encoded for JSON/AT Protocol transport.
 * The relay sees only this opaque structure — it cannot decrypt the content.
 *
 * @param sessionId     Unique session identifier for this comparison
 * @param encryptedData Base64-encoded ciphertext (AES-256-GCM output)
 * @param iv            Base64-encoded initialization vector (12 bytes for GCM)
 * @param dataType      What kind of data is encrypted
 * @param senderKeyId   Identifier for the sender's public key (for key lookup)
 */
case class EncryptedPayload(
                             sessionId: String,
                             encryptedData: String,
                             iv: String,
                             dataType: PayloadDataType,
                             senderKeyId: Option[String] = None
                           )

object EncryptedPayload:

  given Encoder[EncryptedPayload] = Encoder.instance { p =>
    Json.obj(
      "sessionId" -> Json.fromString(p.sessionId),
      "encryptedData" -> Json.fromString(p.encryptedData),
      "iv" -> Json.fromString(p.iv),
      "dataType" -> Json.fromString(p.dataType.toDbString),
      "senderKeyId" -> p.senderKeyId.fold(Json.Null)(Json.fromString)
    )
  }

  given Decoder[EncryptedPayload] = Decoder.instance { c =>
    for
      sessionId <- c.get[String]("sessionId")
      encryptedData <- c.get[String]("encryptedData")
      iv <- c.get[String]("iv")
      dataType <- c.get[String]("dataType").map(PayloadDataType.fromString)
      senderKeyId <- c.get[Option[String]]("senderKeyId")
    yield EncryptedPayload(sessionId, encryptedData, iv, dataType, senderKeyId)
  }
