package com.decodingus.ibd.crypto

import munit.FunSuite

import java.nio.charset.StandardCharsets
import java.util.Base64

class IbdCryptoServiceSpec extends FunSuite:

  // ============================================
  // X25519 Key Exchange
  // ============================================

  test("generateX25519KeyPair produces valid key pair") {
    val kp = IbdCryptoService.generateX25519KeyPair()
    assert(kp.getPublic != null)
    assert(kp.getPrivate != null)
    assertEquals(kp.getPublic.getAlgorithm, "XDH")
  }

  test("two parties derive the same shared secret") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()

    val aliceSecret = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)
    val bobSecret = IbdCryptoService.deriveSharedSecret(bob.getPrivate, alice.getPublic)

    // Both should derive identical AES keys
    assert(java.util.Arrays.equals(aliceSecret.getEncoded, bobSecret.getEncoded))
    assertEquals(aliceSecret.getAlgorithm, "AES")
    assertEquals(aliceSecret.getEncoded.length, 32) // 256 bits
  }

  test("different key pairs produce different shared secrets") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()
    val carol = IbdCryptoService.generateX25519KeyPair()

    val aliceBob = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)
    val aliceCarol = IbdCryptoService.deriveSharedSecret(alice.getPrivate, carol.getPublic)

    assert(!java.util.Arrays.equals(aliceBob.getEncoded, aliceCarol.getEncoded))
  }

  // ============================================
  // AES-256-GCM Encryption
  // ============================================

  test("encrypt and decrypt round trip") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()
    val key = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)

    val plaintext = "IBD segment data for chromosome 1".getBytes(StandardCharsets.UTF_8)
    val (ciphertext, iv) = IbdCryptoService.encrypt(plaintext, key)
    val decrypted = IbdCryptoService.decrypt(ciphertext, iv, key)

    assert(java.util.Arrays.equals(plaintext, decrypted))
  }

  test("encrypt produces different ciphertext each time (random IV)") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()
    val key = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)

    val plaintext = "same data".getBytes(StandardCharsets.UTF_8)
    val (ct1, iv1) = IbdCryptoService.encrypt(plaintext, key)
    val (ct2, iv2) = IbdCryptoService.encrypt(plaintext, key)

    // IVs should be different (random)
    assert(!java.util.Arrays.equals(iv1, iv2))
    // Ciphertext should be different
    assert(!java.util.Arrays.equals(ct1, ct2))
  }

  test("decrypt with wrong key fails") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()
    val carol = IbdCryptoService.generateX25519KeyPair()

    val correctKey = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)
    val wrongKey = IbdCryptoService.deriveSharedSecret(alice.getPrivate, carol.getPublic)

    val plaintext = "sensitive data".getBytes(StandardCharsets.UTF_8)
    val (ciphertext, iv) = IbdCryptoService.encrypt(plaintext, correctKey)

    intercept[javax.crypto.AEADBadTagException] {
      IbdCryptoService.decrypt(ciphertext, iv, wrongKey)
    }
  }

  test("decrypt with tampered ciphertext fails") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()
    val key = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)

    val plaintext = "sensitive data".getBytes(StandardCharsets.UTF_8)
    val (ciphertext, iv) = IbdCryptoService.encrypt(plaintext, key)

    // Tamper with the ciphertext
    val tampered = ciphertext.clone()
    tampered(0) = (tampered(0) ^ 0xFF).toByte

    intercept[javax.crypto.AEADBadTagException] {
      IbdCryptoService.decrypt(tampered, iv, key)
    }
  }

  test("encrypt handles empty data") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()
    val key = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)

    val plaintext = Array.emptyByteArray
    val (ciphertext, iv) = IbdCryptoService.encrypt(plaintext, key)
    val decrypted = IbdCryptoService.decrypt(ciphertext, iv, key)

    assert(java.util.Arrays.equals(plaintext, decrypted))
  }

  test("encrypt handles large data") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()
    val key = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)

    // 1MB of data (simulating genotype payload)
    val plaintext = new Array[Byte](1024 * 1024)
    new java.security.SecureRandom().nextBytes(plaintext)

    val (ciphertext, iv) = IbdCryptoService.encrypt(plaintext, key)
    val decrypted = IbdCryptoService.decrypt(ciphertext, iv, key)

    assert(java.util.Arrays.equals(plaintext, decrypted))
  }

  // ============================================
  // EncryptedPayload Round Trip
  // ============================================

  test("encryptToPayload and decryptPayload round trip") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()
    val key = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)

    val plaintext = """{"chromosome":"1","start":1000000,"end":5000000}""".getBytes(StandardCharsets.UTF_8)
    val payload = IbdCryptoService.encryptToPayload(
      plaintext, key,
      sessionId = "session-001",
      dataType = PayloadDataType.Genotypes,
      senderKeyId = Some("alice-key-1")
    )

    assertEquals(payload.sessionId, "session-001")
    assertEquals(payload.dataType, PayloadDataType.Genotypes)
    assertEquals(payload.senderKeyId, Some("alice-key-1"))

    // Bob decrypts with the same shared secret
    val bobKey = IbdCryptoService.deriveSharedSecret(bob.getPrivate, alice.getPublic)
    val decrypted = IbdCryptoService.decryptPayload(payload, bobKey)
    assert(java.util.Arrays.equals(plaintext, decrypted))
  }

  test("EncryptedPayload fields are valid Base64") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()
    val key = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)

    val plaintext = "test data".getBytes(StandardCharsets.UTF_8)
    val payload = IbdCryptoService.encryptToPayload(plaintext, key, "s1", PayloadDataType.Genotypes)

    // All fields should be valid Base64
    val decoder = Base64.getDecoder
    assert(decoder.decode(payload.encryptedData).nonEmpty)
    assertEquals(decoder.decode(payload.iv).length, 12) // GCM IV is 12 bytes
  }

  // ============================================
  // Ed25519 Signatures
  // ============================================

  test("generateEd25519KeyPair produces valid key pair") {
    val kp = IbdCryptoService.generateEd25519KeyPair()
    assert(kp.getPublic != null)
    assert(kp.getPrivate != null)
    assertEquals(kp.getPublic.getAlgorithm, "EdDSA")
  }

  test("sign and verify attestation round trip") {
    val kp = IbdCryptoService.generateEd25519KeyPair()

    val attestation = "match:session-001:150.5cM:5segments:sha256hash".getBytes(StandardCharsets.UTF_8)
    val signature = IbdCryptoService.signAttestation(attestation, kp.getPrivate)

    assert(signature.nonEmpty)
    assert(IbdCryptoService.verifyAttestation(attestation, signature, kp.getPublic))
  }

  test("verify fails with wrong public key") {
    val alice = IbdCryptoService.generateEd25519KeyPair()
    val bob = IbdCryptoService.generateEd25519KeyPair()

    val data = "attestation data".getBytes(StandardCharsets.UTF_8)
    val signature = IbdCryptoService.signAttestation(data, alice.getPrivate)

    // Verify with wrong key should fail
    assert(!IbdCryptoService.verifyAttestation(data, signature, bob.getPublic))
  }

  test("verify fails with tampered data") {
    val kp = IbdCryptoService.generateEd25519KeyPair()

    val data = "attestation data".getBytes(StandardCharsets.UTF_8)
    val signature = IbdCryptoService.signAttestation(data, kp.getPrivate)

    val tampered = "tampered data".getBytes(StandardCharsets.UTF_8)
    assert(!IbdCryptoService.verifyAttestation(tampered, signature, kp.getPublic))
  }

  test("verify fails with tampered signature") {
    val kp = IbdCryptoService.generateEd25519KeyPair()

    val data = "attestation data".getBytes(StandardCharsets.UTF_8)
    val signature = IbdCryptoService.signAttestation(data, kp.getPrivate)

    val tamperedSig = signature.clone()
    tamperedSig(0) = (tamperedSig(0) ^ 0xFF).toByte
    assert(!IbdCryptoService.verifyAttestation(data, tamperedSig, kp.getPublic))
  }

  // ============================================
  // Key Serialization
  // ============================================

  test("X25519 public key encode/decode round trip") {
    val kp = IbdCryptoService.generateX25519KeyPair()
    val encoded = IbdCryptoService.encodePublicKey(kp.getPublic)

    assert(encoded.nonEmpty)
    // Should be valid Base64
    val decoded = IbdCryptoService.decodeX25519PublicKey(encoded)
    assert(java.util.Arrays.equals(kp.getPublic.getEncoded, decoded.getEncoded))
  }

  test("Ed25519 public key encode/decode round trip") {
    val kp = IbdCryptoService.generateEd25519KeyPair()
    val encoded = IbdCryptoService.encodePublicKey(kp.getPublic)

    assert(encoded.nonEmpty)
    val decoded = IbdCryptoService.decodeEd25519PublicKey(encoded)
    assert(java.util.Arrays.equals(kp.getPublic.getEncoded, decoded.getEncoded))
  }

  test("decoded X25519 key works for ECDH") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()

    // Simulate Alice sending her public key as Base64 over the network
    val alicePubEncoded = IbdCryptoService.encodePublicKey(alice.getPublic)
    val alicePubDecoded = IbdCryptoService.decodeX25519PublicKey(alicePubEncoded)

    // Bob uses the decoded key for key agreement
    val bobSecret = IbdCryptoService.deriveSharedSecret(bob.getPrivate, alicePubDecoded)
    val aliceSecret = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)

    assert(java.util.Arrays.equals(aliceSecret.getEncoded, bobSecret.getEncoded))
  }

  test("decoded Ed25519 key works for verification") {
    val kp = IbdCryptoService.generateEd25519KeyPair()

    // Simulate sending the public key over the network
    val pubEncoded = IbdCryptoService.encodePublicKey(kp.getPublic)
    val pubDecoded = IbdCryptoService.decodeEd25519PublicKey(pubEncoded)

    val data = "attestation data".getBytes(StandardCharsets.UTF_8)
    val signature = IbdCryptoService.signAttestation(data, kp.getPrivate)

    assert(IbdCryptoService.verifyAttestation(data, signature, pubDecoded))
  }

  // ============================================
  // SHA-256 Hashing
  // ============================================

  test("sha256Hex produces consistent hex output") {
    val data = "match:session-001:150.5cM".getBytes(StandardCharsets.UTF_8)
    val hash1 = IbdCryptoService.sha256Hex(data)
    val hash2 = IbdCryptoService.sha256Hex(data)

    assertEquals(hash1, hash2)
    assertEquals(hash1.length, 64) // SHA-256 = 32 bytes = 64 hex chars
    assert(hash1.matches("[0-9a-f]+"))
  }

  test("sha256Hex produces different output for different input") {
    val hash1 = IbdCryptoService.sha256Hex("data1".getBytes(StandardCharsets.UTF_8))
    val hash2 = IbdCryptoService.sha256Hex("data2".getBytes(StandardCharsets.UTF_8))
    assert(hash1 != hash2)
  }

  // ============================================
  // Full E2E Crypto Flow
  // ============================================

  test("full E2E: key exchange, encrypt, transport, decrypt, sign, verify") {
    // 1. Both parties generate ephemeral X25519 keys
    val aliceX = IbdCryptoService.generateX25519KeyPair()
    val bobX = IbdCryptoService.generateX25519KeyPair()

    // 2. Both derive shared secret
    val aliceKey = IbdCryptoService.deriveSharedSecret(aliceX.getPrivate, bobX.getPublic)
    val bobKey = IbdCryptoService.deriveSharedSecret(bobX.getPrivate, aliceX.getPublic)

    // 3. Alice encrypts genotype data and sends as payload
    val genotypeData = """{"chr1":[0,1,2,0,1,1,2,0]}""".getBytes(StandardCharsets.UTF_8)
    val payload = IbdCryptoService.encryptToPayload(
      genotypeData, aliceKey, "session-001", PayloadDataType.Genotypes, Some("alice-x-key")
    )

    // 4. Bob decrypts (simulating receiving via relay)
    val decrypted = IbdCryptoService.decryptPayload(payload, bobKey)
    assert(java.util.Arrays.equals(genotypeData, decrypted))

    // 5. Both compute IBD result hash
    val ibdResult = "150.5cM:5segments:45.2longest"
    val aliceHash = IbdCryptoService.sha256Hex(ibdResult.getBytes(StandardCharsets.UTF_8))
    val bobHash = IbdCryptoService.sha256Hex(ibdResult.getBytes(StandardCharsets.UTF_8))
    assertEquals(aliceHash, bobHash) // Both computed the same result

    // 6. Both sign attestations with their Ed25519 keys
    val aliceEd = IbdCryptoService.generateEd25519KeyPair()
    val bobEd = IbdCryptoService.generateEd25519KeyPair()

    val attestationData = s"match:session-001:$aliceHash".getBytes(StandardCharsets.UTF_8)
    val aliceSig = IbdCryptoService.signAttestation(attestationData, aliceEd.getPrivate)
    val bobSig = IbdCryptoService.signAttestation(attestationData, bobEd.getPrivate)

    // 7. Both can verify each other's signatures
    assert(IbdCryptoService.verifyAttestation(attestationData, aliceSig, aliceEd.getPublic))
    assert(IbdCryptoService.verifyAttestation(attestationData, bobSig, bobEd.getPublic))
  }
