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
      dataType = "GENOTYPES",
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
      dataType = "ATTESTATION"
    )

    val json = payload.asJson.noSpaces
    val decoded = decode[EncryptedPayload](json)

    assert(decoded.isRight)
    decoded.foreach { p =>
      assertEquals(p.sessionId, "session-002")
      assertEquals(p.dataType, "ATTESTATION")
      assertEquals(p.senderKeyId, None)
    }
  }

  test("EncryptedPayload.DataTypes contains expected values") {
    assert(EncryptedPayload.DataTypes.contains("GENOTYPES"))
    assert(EncryptedPayload.DataTypes.contains("ATTESTATION"))
    assert(EncryptedPayload.DataTypes.contains("SUMMARY"))
    assert(EncryptedPayload.DataTypes.contains("KEY_EXCHANGE"))
  }
