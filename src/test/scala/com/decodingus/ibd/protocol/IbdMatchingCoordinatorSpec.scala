package com.decodingus.ibd.protocol

import com.decodingus.ibd.crypto.{EncryptedPayload, IbdCryptoService, PayloadDataType}
import com.decodingus.ibd.engine.*
import com.decodingus.workspace.model.{IbdSegment, RelationshipEstimate}
import munit.FunSuite

import java.nio.charset.StandardCharsets
import java.security.KeyPair
import java.util.Base64
import java.util.concurrent.{ConcurrentLinkedQueue, CountDownLatch, TimeUnit}
import javax.crypto.SecretKey

/**
 * Tests for the IbdMatchingCoordinator protocol logic.
 *
 * These tests verify the coordinator's internal protocol steps using
 * direct method calls rather than a live WebSocket relay. The relay
 * transport is tested separately via IbdRelayClient tests.
 */
class IbdMatchingCoordinatorSpec extends FunSuite:

  test("state machine starts at Idle") {
    val config = MatchingConfig(relayUrl = "wss://test/relay")
    val kp = IbdCryptoService.generateEd25519KeyPair()
    val coordinator = IbdMatchingCoordinator(
      config = config,
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "test-session",
      authToken = "test-token",
      localDid = "did:plc:a",
      localSampleRef = "at://did:plc:a/bio/1",
      partnerSampleRef = "at://did:plc:b/bio/1",
      vcfFile = java.io.File("nonexistent.vcf"),
      signingKeyPair = kp
    )
    assertEquals(coordinator.currentState, MatchingProtocolState.Idle)
  }

  test("two-party key exchange produces matching shared secrets") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()

    val aliceSecret = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)
    val bobSecret = IbdCryptoService.deriveSharedSecret(bob.getPrivate, alice.getPublic)

    // Encrypt with Alice's key, decrypt with Bob's
    val plaintext = "test data for key exchange".getBytes(StandardCharsets.UTF_8)
    val (ciphertext, iv) = IbdCryptoService.encrypt(plaintext, aliceSecret)
    val decrypted = IbdCryptoService.decrypt(ciphertext, iv, bobSecret)

    assert(java.util.Arrays.equals(plaintext, decrypted),
      "Shared secrets should produce interoperable encryption")
  }

  test("variant compact encoding round-trips through encryption") {
    val alice = IbdCryptoService.generateX25519KeyPair()
    val bob = IbdCryptoService.generateX25519KeyPair()
    val sharedKey = IbdCryptoService.deriveSharedSecret(alice.getPrivate, bob.getPublic)

    // Create synthetic genotype data
    val positions = Array.tabulate(100)(i => (i + 1) * 10000)
    val genotypes = Array.tabulate(100)(i => (i % 3).toByte)
    val original = ChromosomeGenotypes("1", positions, genotypes)

    // Compact encode → encrypt → decrypt → compact decode
    val compactBytes = IbdVariantExtractor.toCompactBytes(original)
    val header = s"1:${compactBytes.length}"
    val headerBytes = header.getBytes(StandardCharsets.UTF_8)
    val payload = IbdCryptoService.encryptToPayload(
      headerBytes ++ compactBytes,
      sharedKey, "session-1", PayloadDataType.Genotypes
    )

    val bobKey = IbdCryptoService.deriveSharedSecret(bob.getPrivate, alice.getPublic)
    val decrypted = IbdCryptoService.decryptPayload(payload, bobKey)

    // Parse header from raw bytes — find the colon separator in header region
    // Header format: "chr:length" as UTF-8, followed by raw binary data
    val headerLen = headerBytes.length
    val headerStr = new String(decrypted, 0, headerLen, StandardCharsets.UTF_8)
    val colonIdx = headerStr.indexOf(':')
    val chr = headerStr.substring(0, colonIdx)
    val dataLength = headerStr.substring(colonIdx + 1).toInt
    val extractedBytes = java.util.Arrays.copyOfRange(decrypted, headerLen, headerLen + dataLength)

    assertEquals(chr, "1")
    assertEquals(dataLength, compactBytes.length)
    val decoded = IbdVariantExtractor.fromCompactBytes(extractedBytes, "1")
    assertEquals(decoded.size, 100)
    for i <- 0 until 100 do
      assertEquals(decoded.positions(i), original.positions(i), s"Position $i mismatch")
      assertEquals(decoded.genotypes(i), original.genotypes(i), s"Genotype $i mismatch")
  }

  test("hash comparison detects agreement") {
    val segments = List(
      IbdSegment("1", 5000000, 10000000, 15.23),
      IbdSegment("3", 20000000, 30000000, 28.76)
    )

    // Both parties compute independently from same segments
    val hash1 = MatchSummary.computeHash(segments)
    val hash2 = MatchSummary.computeHash(segments.reverse) // Different order

    assertEquals(hash1, hash2, "Canonical hashing should be order-independent")
  }

  test("hash comparison detects disagreement") {
    val segments1 = List(IbdSegment("1", 5000000, 10000000, 15.23))
    val segments2 = List(IbdSegment("1", 5000000, 10000000, 15.24)) // Slightly different cM

    val hash1 = MatchSummary.computeHash(segments1)
    val hash2 = MatchSummary.computeHash(segments2)

    assert(hash1 != hash2, "Different segments should produce different hashes")
  }

  test("full protocol simulation with two in-process parties") {
    // This test simulates the protocol without a real WebSocket relay.
    // Instead, we use a direct message queue between two "navigators".

    val geneticMap = GeneticMap.uniformRate(cmPerMb = 1.0)
    val ibdConfig = IbdDetectorConfig(
      minSegmentCm = 3.0,
      minSnpCount = 50,
      windowSize = 50,
      ibsThreshold = 0.65,
      errorTolerance = 0.02
    )

    // Create synthetic genotype data where the two share an IBD segment
    val rng = new java.util.Random(42L)
    val n = 2000
    val spacing = 10000 // 10kb spacing → 20Mb total
    val positions = Array.tabulate(n)(i => (i + 1) * spacing)

    // Sample 1: random genotypes
    val g1 = Array.tabulate(n)(_ => rng.nextInt(3).toByte)

    // Sample 2: IBD from SNP 500-1000 (5M-10M bp), random elsewhere
    val rng2 = new java.util.Random(99L)
    val g2 = Array.tabulate(n) { i =>
      if i >= 500 && i < 1000 then g1(i) // IBD region
      else rng2.nextInt(3).toByte
    }

    val aliceGenotypes = Map("1" -> ChromosomeGenotypes("1", positions, g1))
    val bobGenotypes = Map("1" -> ChromosomeGenotypes("1", positions, g2))

    // Both parties run IBD detection independently
    val detector = PairwiseIbdDetector(ibdConfig)
    val aliceSegments = detector.detectSegments(aliceGenotypes, bobGenotypes, geneticMap)
    val bobSegments = detector.detectSegments(bobGenotypes, aliceGenotypes, geneticMap)

    // Both should find at least one segment
    assert(aliceSegments.nonEmpty, "Alice should detect segments")
    assert(bobSegments.nonEmpty, "Bob should detect segments")

    // Both compute summaries
    val aliceSummary = MatchSummary.fromSegments(aliceSegments)
    val bobSummary = MatchSummary.fromSegments(bobSegments)

    // Hashes should match (mutual agreement)
    assertEquals(aliceSummary.summaryHash, bobSummary.summaryHash,
      "Both parties should compute matching hashes")
    assertEquals(aliceSummary.totalSharedCm, bobSummary.totalSharedCm)
    assertEquals(aliceSummary.segmentCount, bobSummary.segmentCount)

    // Both sign attestations
    val aliceSigningKey = IbdCryptoService.generateEd25519KeyPair()
    val bobSigningKey = IbdCryptoService.generateEd25519KeyPair()

    val aliceAttestation = IbdAttestation.create(
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "session-001",
      attestingDid = "did:plc:a",
      attestingSampleRef = "at://did:plc:a/bio/1",
      partnerSampleRef = "at://did:plc:b/bio/1",
      summary = aliceSummary,
      partnerSummaryHash = bobSummary.summaryHash,
      signingKeyPair = aliceSigningKey
    )

    val bobAttestation = IbdAttestation.create(
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "session-001",
      attestingDid = "did:plc:b",
      attestingSampleRef = "at://did:plc:b/bio/1",
      partnerSampleRef = "at://did:plc:a/bio/1",
      summary = bobSummary,
      partnerSummaryHash = aliceSummary.summaryHash,
      signingKeyPair = bobSigningKey
    )

    // Both attestations should verify
    assert(IbdAttestation.verify(aliceAttestation), "Alice's attestation should verify")
    assert(IbdAttestation.verify(bobAttestation), "Bob's attestation should verify")

    // Both should have matching summary hashes
    assertEquals(aliceAttestation.summaryHash, bobAttestation.summaryHash)

    // Shared cM should be reasonable for the IBD region (500 SNPs * 10kb = 5Mb ≈ 5 cM)
    assert(aliceSummary.totalSharedCm >= 3.0,
      s"Expected >= 3.0 cM shared, got ${aliceSummary.totalSharedCm}")
  }

  test("encrypted key exchange between two parties") {
    val sessionId = "test-session-key-exchange"

    // Derive same session key (used before ECDH)
    val digest = java.security.MessageDigest.getInstance("SHA-256")
    val sessionKeyBytes = digest.digest(s"ibd-session-key:$sessionId".getBytes(StandardCharsets.UTF_8))
    val sessionKey = javax.crypto.spec.SecretKeySpec(sessionKeyBytes, "AES")

    // Alice generates X25519 key pair and encrypts public key
    val aliceX25519 = IbdCryptoService.generateX25519KeyPair()
    val aliceKeyBase64 = IbdCryptoService.encodePublicKey(aliceX25519.getPublic)
    val alicePayload = IbdCryptoService.encryptToPayload(
      aliceKeyBase64.getBytes(StandardCharsets.UTF_8),
      sessionKey, sessionId, PayloadDataType.KeyExchange
    )

    // Bob decrypts Alice's key
    val decryptedBytes = IbdCryptoService.decryptPayload(alicePayload, sessionKey)
    val recoveredKeyBase64 = new String(decryptedBytes, StandardCharsets.UTF_8)
    val aliceRecoveredKey = IbdCryptoService.decodeX25519PublicKey(recoveredKeyBase64)

    // Bob generates his own key pair
    val bobX25519 = IbdCryptoService.generateX25519KeyPair()

    // Both derive shared secret
    val aliceShared = IbdCryptoService.deriveSharedSecret(aliceX25519.getPrivate, bobX25519.getPublic)
    val bobShared = IbdCryptoService.deriveSharedSecret(bobX25519.getPrivate, aliceRecoveredKey)

    // Verify they can encrypt/decrypt
    val testData = "mutual authentication test".getBytes(StandardCharsets.UTF_8)
    val encrypted = IbdCryptoService.encryptToPayload(testData, aliceShared, sessionId, PayloadDataType.Summary)
    val decrypted = IbdCryptoService.decryptPayload(encrypted, bobShared)
    assert(java.util.Arrays.equals(testData, decrypted))
  }

  test("attestation exchange and cross-verification") {
    val segments = List(
      IbdSegment("1", 5000000, 10000000, 15.23),
      IbdSegment("3", 20000000, 30000000, 28.76)
    )
    val summary = MatchSummary.fromSegments(segments)

    // Alice creates and signs attestation
    val aliceKey = IbdCryptoService.generateEd25519KeyPair()
    val aliceAttestation = IbdAttestation.create(
      matchRequestUri = "at://did:plc:a/matchrequest/1",
      sessionId = "s1",
      attestingDid = "did:plc:a",
      attestingSampleRef = "at://did:plc:a/bio/1",
      partnerSampleRef = "at://did:plc:b/bio/1",
      summary = summary,
      partnerSummaryHash = summary.summaryHash,
      signingKeyPair = aliceKey
    )

    // Encrypt attestation for transport
    val x25519Alice = IbdCryptoService.generateX25519KeyPair()
    val x25519Bob = IbdCryptoService.generateX25519KeyPair()
    val sharedKey = IbdCryptoService.deriveSharedSecret(x25519Alice.getPrivate, x25519Bob.getPublic)

    import io.circe.syntax.*
    val attestationJson = aliceAttestation.asJson.noSpaces
    val encPayload = IbdCryptoService.encryptToPayload(
      attestationJson.getBytes(StandardCharsets.UTF_8),
      sharedKey, "s1", PayloadDataType.Attestation
    )

    // Bob decrypts and verifies
    val bobKey = IbdCryptoService.deriveSharedSecret(x25519Bob.getPrivate, x25519Alice.getPublic)
    val decryptedBytes = IbdCryptoService.decryptPayload(encPayload, bobKey)
    val decryptedJson = new String(decryptedBytes, StandardCharsets.UTF_8)

    import io.circe.parser.decode
    val Right(received) = decode[IbdAttestation](decryptedJson): @unchecked

    assert(IbdAttestation.verify(received), "Bob should verify Alice's attestation")
    assertEquals(received.attestingDid, "did:plc:a")
    assertEquals(received.summaryHash, summary.summaryHash)
  }

  test("hash mismatch detection between parties with different results") {
    val geneticMap = GeneticMap.uniformRate(cmPerMb = 1.0)
    val ibdConfig = IbdDetectorConfig(
      minSegmentCm = 1.0,
      minSnpCount = 20,
      windowSize = 20
    )

    // Create two completely unrelated individuals
    val rng1 = new java.util.Random(100L)
    val rng2 = new java.util.Random(200L)
    val n = 500
    val positions = Array.tabulate(n)(i => (i + 1) * 10000)
    val g1 = Array.tabulate(n)(_ => rng1.nextInt(3).toByte)
    val g2 = Array.tabulate(n)(_ => rng2.nextInt(3).toByte)
    val g3 = Array.tabulate(n)(_ => rng1.nextInt(3).toByte) // Different from g2

    val genotypes1 = Map("1" -> ChromosomeGenotypes("1", positions, g1))
    val genotypes2 = Map("1" -> ChromosomeGenotypes("1", positions, g2))
    val genotypes3 = Map("1" -> ChromosomeGenotypes("1", positions, g3))

    val detector = PairwiseIbdDetector(ibdConfig)

    // Alice compares with Bob
    val aliceSegments = detector.detectSegments(genotypes1, genotypes2, geneticMap)
    // Alice accidentally compares with Charlie (wrong partner!)
    val wrongSegments = detector.detectSegments(genotypes1, genotypes3, geneticMap)

    val aliceHash = MatchSummary.fromSegments(aliceSegments).summaryHash
    val wrongHash = MatchSummary.fromSegments(wrongSegments).summaryHash

    // If the segments differ at all, the hashes should differ
    // (They might both be empty if no segments found, which would match — that's OK)
    if aliceSegments != wrongSegments then
      assert(aliceHash != wrongHash, "Different segment sets should produce different hashes")
  }
