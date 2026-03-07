package com.decodingus.ibd.crypto

import io.circe.parser.decode
import io.circe.syntax.*
import munit.FunSuite

class EncryptedPayloadSpec extends FunSuite:

  test("EncryptedPayload JSON round trip") {
    val payload = EncryptedPayload(
      sessionId = "session-001",
      encryptedData = "dGVzdCBkYXRh",
      iv = "aXYxMjM0NTY3ODkw",
      dataType = PayloadDataType.Genotypes,
      senderKeyId = Some("key-123")
    )

    val json = payload.asJson.noSpaces
    val decoded = decode[EncryptedPayload](json)

    assert(decoded.isRight)
    assertEquals(decoded.toOption.get, payload)
  }

  test("EncryptedPayload JSON round trip without senderKeyId") {
    val payload = EncryptedPayload(
      sessionId = "session-002",
      encryptedData = "Y2lwaGVydGV4dA==",
      iv = "aXYxMjM0NTY=",
      dataType = PayloadDataType.Attestation
    )

    val json = payload.asJson.noSpaces
    val decoded = decode[EncryptedPayload](json)

    assert(decoded.isRight)
    decoded.foreach { p =>
      assertEquals(p.sessionId, "session-002")
      assertEquals(p.dataType, PayloadDataType.Attestation)
      assertEquals(p.senderKeyId, None)
    }
  }

  test("PayloadDataType has all expected values") {
    assertEquals(PayloadDataType.values.length, 4)
    assertEquals(PayloadDataType.fromString("GENOTYPES"), PayloadDataType.Genotypes)
    assertEquals(PayloadDataType.fromString("ATTESTATION"), PayloadDataType.Attestation)
    assertEquals(PayloadDataType.fromString("SUMMARY"), PayloadDataType.Summary)
    assertEquals(PayloadDataType.fromString("KEY_EXCHANGE"), PayloadDataType.KeyExchange)
  }
