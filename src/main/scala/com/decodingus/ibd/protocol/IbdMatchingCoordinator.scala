package com.decodingus.ibd.protocol

import com.decodingus.analysis.ProgressReporter.ProgressCallback
import com.decodingus.ibd.crypto.{EncryptedPayload, IbdCryptoService, PayloadDataType}
import com.decodingus.ibd.engine.{GeneticMap, IbdDetectorConfig, IbdVariantExtractor, MatchSummary, PairwiseIbdDetector}
import com.decodingus.ibd.relay.IbdRelayClient
import com.decodingus.util.Logger
import com.decodingus.workspace.model.{IbdSegment, RelationshipEstimate}

import java.io.File
import java.nio.charset.StandardCharsets
import java.security.{KeyPair, PublicKey}
import java.util.Base64
import java.util.concurrent.{CountDownLatch, TimeUnit}
import java.util.concurrent.atomic.AtomicReference
import javax.crypto.SecretKey
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Failure, Success, Try}

/**
 * Configuration for the IBD matching protocol.
 *
 * @param relayUrl           WebSocket relay base URL
 * @param timeoutSeconds     Maximum time to wait for partner in any phase
 * @param ibdConfig          IBD detection algorithm configuration
 * @param geneticMap         Genetic map for cM conversion
 */
case class MatchingConfig(
                           relayUrl: String,
                           timeoutSeconds: Int = 300,
                           ibdConfig: IbdDetectorConfig = IbdDetectorConfig(),
                           geneticMap: GeneticMap = GeneticMap.uniformRate()
                         )

/**
 * Result of a completed matching protocol run.
 */
case class MatchingProtocolResult(
                                   summary: MatchSummary,
                                   attestation: IbdAttestation,
                                   partnerAttestation: Option[IbdAttestation]
                                 )

/**
 * Orchestrates the full IBD matching protocol between two Navigator instances.
 *
 * The protocol proceeds through these steps:
 *   1. Connect to relay
 *   2. Exchange X25519 public keys
 *   3. Derive shared AES-256 secret
 *   4. Extract local variant genotypes from VCF
 *   5. Encrypt and send variants through relay
 *   6. Receive and decrypt partner's variants
 *   7. Run pairwise IBD detection
 *   8. Exchange summary hashes for mutual verification
 *   9. If hashes agree, sign attestation
 *  10. Submit attestation to AppView
 *  11. Persist match result locally
 *
 * Error handling:
 *   - Partner disconnect → TimedOut state
 *   - Hash mismatch → Failed state (logged for investigation)
 *   - Relay failure → Retry with backoff (handled by IbdRelayClient)
 *
 * @param config           Protocol configuration
 * @param matchRequestUri  AT URI of the match request being fulfilled
 * @param sessionId        Unique session ID for this comparison
 * @param authToken        Auth token for relay access
 * @param localDid         DID of the local citizen
 * @param localSampleRef   AT URI of the local biosample
 * @param partnerSampleRef AT URI of the partner biosample
 * @param vcfFile          Local VCF/GVCF file with variant calls
 * @param signingKeyPair   Ed25519 key pair for attestation signing
 * @param onProgress       Optional progress callback
 */
class IbdMatchingCoordinator(
                              config: MatchingConfig,
                              matchRequestUri: String,
                              sessionId: String,
                              authToken: String,
                              localDid: String,
                              localSampleRef: String,
                              partnerSampleRef: String,
                              vcfFile: File,
                              signingKeyPair: KeyPair,
                              onProgress: ProgressCallback = (_, _, _) => ()
                            ):
  private val log = Logger[IbdMatchingCoordinator]
  private val state = AtomicReference[MatchingProtocolState](MatchingProtocolState.Idle)

  // Shared mutable state protected by latches for synchronization
  private val partnerKeyReceived = CountDownLatch(1)
  private val partnerVariantsReceived = CountDownLatch(1)
  private val partnerHashReceived = CountDownLatch(1)
  private val partnerAttestationReceived = CountDownLatch(1)

  private val partnerPublicKey = AtomicReference[Option[PublicKey]](None)
  private val partnerVariantData = AtomicReference[Option[Map[String, Array[Byte]]]](None)
  private val partnerHash = AtomicReference[Option[String]](None)
  private val partnerAttestationRef = AtomicReference[Option[IbdAttestation]](None)

  private var relayClient: Option[IbdRelayClient] = None

  def currentState: MatchingProtocolState = state.get()

  /**
   * Execute the full matching protocol.
   *
   * This method blocks the calling thread until completion, timeout, or failure.
   * Run it on a background thread.
   *
   * @return Either an error message or the protocol result
   */
  def execute()(implicit ec: ExecutionContext): Either[String, MatchingProtocolResult] =
    try
      executeSteps()
    catch
      case e: ProtocolAbortException =>
        cleanup()
        Left(e.getMessage)
      case e: Exception =>
        val msg = s"Protocol error: ${e.getMessage}"
        log.error(msg, e)
        transition(MatchingProtocolState.Failed(msg))
        cleanup()
        Left(msg)

  private def executeSteps()(implicit ec: ExecutionContext): Either[String, MatchingProtocolResult] =
    transition(MatchingProtocolState.Connecting)

    // Step 1: Connect to relay
    val relay = IbdRelayClient(config.relayUrl, sessionId, authToken)
    relay.onMessage(handleIncomingMessage)
    relay.onError(e => log.error(s"Relay error: ${e.getMessage}"))
    relayClient = Some(relay)

    awaitFuture(relay.connect(), "relay connection")

    // Step 2: Exchange keys
    transition(MatchingProtocolState.ExchangingKeys)
    val x25519KeyPair = IbdCryptoService.generateX25519KeyPair()
    val keyPayload = IbdCryptoService.encryptToPayload(
      IbdCryptoService.encodePublicKey(x25519KeyPair.getPublic).getBytes(StandardCharsets.UTF_8),
      deriveSessionKey(sessionId),
      sessionId,
      PayloadDataType.KeyExchange
    )
    awaitFuture(relay.send(keyPayload), "key exchange send")

    // Wait for partner's key
    if !partnerKeyReceived.await(config.timeoutSeconds, TimeUnit.SECONDS) then
      transition(MatchingProtocolState.TimedOut)
      throw ProtocolAbortException("Timed out waiting for partner's public key")

    val theirPublicKey = partnerPublicKey.get().getOrElse(
      throw ProtocolAbortException("Partner public key not received")
    )

    // Step 3: Derive shared secret
    val sharedSecret = IbdCryptoService.deriveSharedSecret(x25519KeyPair.getPrivate, theirPublicKey)
    setSharedSecretForReceive(sharedSecret)

    // Step 4: Extract local variants
    transition(MatchingProtocolState.ExtractingVariants)
    val localVariants = IbdVariantExtractor.extractFromVcf(vcfFile) match
      case Right(variants) => variants
      case Left(err) =>
        transition(MatchingProtocolState.Failed(s"Variant extraction failed: $err"))
        throw ProtocolAbortException(s"Variant extraction failed: $err")

    // Step 5: Encrypt and send variants
    transition(MatchingProtocolState.SendingVariants)
    for (chr, genotypes) <- localVariants do
      val compactBytes = IbdVariantExtractor.toCompactBytes(genotypes)
      val headerBytes = s"$chr:${compactBytes.length}".getBytes(StandardCharsets.UTF_8)
      val payload = IbdCryptoService.encryptToPayload(
        headerBytes ++ compactBytes,
        sharedSecret, sessionId, PayloadDataType.Genotypes
      )
      awaitFuture(relay.send(payload), s"send chr $chr variants")

    // Send "done" marker
    val donePayload = IbdCryptoService.encryptToPayload(
      "VARIANTS_COMPLETE".getBytes(StandardCharsets.UTF_8),
      sharedSecret, sessionId, PayloadDataType.Genotypes
    )
    awaitFuture(relay.send(donePayload), "send variants complete marker")

    // Step 6: Wait for partner's variants
    transition(MatchingProtocolState.ReceivingVariants)
    if !partnerVariantsReceived.await(config.timeoutSeconds, TimeUnit.SECONDS) then
      transition(MatchingProtocolState.TimedOut)
      throw ProtocolAbortException("Timed out waiting for partner's variants")

    val receivedVariantBytes = partnerVariantData.get().getOrElse(
      throw ProtocolAbortException("Partner variant data not received")
    )

    // Decode partner variants
    val partnerGenotypes = receivedVariantBytes.map { (chr, bytes) =>
      chr -> IbdVariantExtractor.fromCompactBytes(bytes, chr)
    }

    // Step 7: Compute IBD
    transition(MatchingProtocolState.ComputingIbd)
    val detector = PairwiseIbdDetector(config.ibdConfig)
    val segments = detector.detectSegments(localVariants, partnerGenotypes, config.geneticMap)
    val summary = MatchSummary.fromSegments(segments)

    log.info(s"IBD detection complete: ${summary.segmentCount} segments, ${summary.totalSharedCm} cM total")

    // Step 8: Exchange hashes
    transition(MatchingProtocolState.ExchangingHashes)
    val hashPayload = IbdCryptoService.encryptToPayload(
      summary.summaryHash.getBytes(StandardCharsets.UTF_8),
      sharedSecret, sessionId, PayloadDataType.Summary
    )
    awaitFuture(relay.send(hashPayload), "send summary hash")

    if !partnerHashReceived.await(config.timeoutSeconds, TimeUnit.SECONDS) then
      transition(MatchingProtocolState.TimedOut)
      throw ProtocolAbortException("Timed out waiting for partner's hash")

    val receivedHash = partnerHash.get().getOrElse(
      throw ProtocolAbortException("Partner hash not received")
    )

    // Step 9: Verify hashes agree
    transition(MatchingProtocolState.VerifyingHashes)
    if summary.summaryHash != receivedHash then
      val msg = s"Hash mismatch: local=${summary.summaryHash.take(16)}... partner=${receivedHash.take(16)}..."
      log.warn(msg)
      transition(MatchingProtocolState.Failed(msg))
      throw ProtocolAbortException(msg)

    log.info("Summary hashes match — mutual agreement confirmed")

    // Step 10: Sign attestation
    transition(MatchingProtocolState.SigningAttestation)
    val attestation = IbdAttestation.create(
      matchRequestUri = matchRequestUri,
      sessionId = sessionId,
      attestingDid = localDid,
      attestingSampleRef = localSampleRef,
      partnerSampleRef = partnerSampleRef,
      summary = summary,
      partnerSummaryHash = receivedHash,
      signingKeyPair = signingKeyPair
    )

    // Send attestation to partner
    import io.circe.syntax.*
    val attestationJson = attestation.asJson.noSpaces
    val attestationPayload = IbdCryptoService.encryptToPayload(
      attestationJson.getBytes(StandardCharsets.UTF_8),
      sharedSecret, sessionId, PayloadDataType.Attestation
    )
    awaitFuture(relay.send(attestationPayload), "send attestation")

    // Wait briefly for partner's attestation (optional — match is valid without it)
    partnerAttestationReceived.await(
      math.min(config.timeoutSeconds, 30).toLong,
      TimeUnit.SECONDS
    )
    val receivedAttestation = partnerAttestationRef.get()

    // Step 11: Submit to AppView
    transition(MatchingProtocolState.SubmittingAttestation)
    // The actual AppView submission is handled by the caller (IbdMatchService)
    // since it requires the DecodingUsClient. We return the attestation for submission.

    // Step 12: Persist locally
    transition(MatchingProtocolState.PersistingResult)
    // Also handled by the caller, since it requires Transactor/Repository access.

    transition(MatchingProtocolState.Completed)
    cleanup()

    Right(MatchingProtocolResult(summary, attestation, receivedAttestation))

  private class ProtocolAbortException(msg: String) extends RuntimeException(msg)

  /**
   * Abort the protocol and disconnect.
   */
  def abort()(implicit ec: ExecutionContext): Unit =
    transition(MatchingProtocolState.Failed("Aborted by user"))
    cleanup()

  private def transition(newState: MatchingProtocolState): Unit =
    val old = state.getAndSet(newState)
    log.debug(s"Protocol state: $old → $newState")
    onProgress(newState.description, newState.progressFraction, 1.0)

  private def handleIncomingMessage(payload: EncryptedPayload): Unit =
    try
      payload.dataType match
        case PayloadDataType.KeyExchange =>
          handleKeyExchange(payload)
        case PayloadDataType.Genotypes =>
          handleGenotypes(payload)
        case PayloadDataType.Summary =>
          handleSummaryHash(payload)
        case PayloadDataType.Attestation =>
          handleAttestation(payload)
    catch
      case e: Exception =>
        log.error(s"Error handling incoming message (${payload.dataType}): ${e.getMessage}")

  private def handleKeyExchange(payload: EncryptedPayload): Unit =
    val sessionKey = deriveSessionKey(sessionId)
    val keyBytes = IbdCryptoService.decryptPayload(payload, sessionKey)
    val keyBase64 = new String(keyBytes, StandardCharsets.UTF_8)
    val publicKey = IbdCryptoService.decodeX25519PublicKey(keyBase64)
    partnerPublicKey.set(Some(publicKey))
    partnerKeyReceived.countDown()
    log.debug("Received partner's X25519 public key")

  private def handleGenotypes(payload: EncryptedPayload): Unit =
    // We need the shared secret to decrypt, but we can't decrypt genotypes
    // until key exchange is complete. The relay guarantees ordering per session,
    // so key exchange will have completed by the time genotypes arrive.
    // Store raw payloads and decrypt them when we have the shared secret.
    // For simplicity, we accumulate into the variant map.
    // The actual decryption is deferred — see the note below.
    accumulateVariantPayload(payload)

  private val pendingVariantPayloads = java.util.concurrent.ConcurrentLinkedQueue[EncryptedPayload]()

  private def accumulateVariantPayload(payload: EncryptedPayload): Unit =
    pendingVariantPayloads.add(payload)
    // We process them when we have the shared secret — triggered by the coordinator thread
    // after key exchange completes. For the "VARIANTS_COMPLETE" marker check, we
    // need to try decrypting to see if it's the done marker. Since key exchange should
    // be complete before any genotype data arrives (protocol ordering), we can
    // process eagerly in a background thread or lazily when needed.

    // Try to process if we already have the shared key
    tryProcessPendingVariants()

  @volatile private var sharedSecretForReceive: Option[SecretKey] = None
  private val receivedChromosomes = java.util.concurrent.ConcurrentHashMap[String, Array[Byte]]()

  private[protocol] def setSharedSecretForReceive(key: SecretKey): Unit =
    sharedSecretForReceive = Some(key)
    tryProcessPendingVariants()

  private val variantsCompleteMarker = "VARIANTS_COMPLETE".getBytes(StandardCharsets.UTF_8)

  private def tryProcessPendingVariants(): Unit =
    sharedSecretForReceive.foreach { key =>
      var payload = pendingVariantPayloads.poll()
      while payload != null do
        val decrypted = IbdCryptoService.decryptPayload(payload, key)
        if java.util.Arrays.equals(decrypted, variantsCompleteMarker) then
          partnerVariantData.set(Some(
            scala.jdk.CollectionConverters.MapHasAsScala(receivedChromosomes).asScala.toMap
          ))
          partnerVariantsReceived.countDown()
          log.debug(s"Received all partner variants (${receivedChromosomes.size()} chromosomes)")
        else
          // Parse "chr:length" header followed by compact binary data.
          // Find the colon in the first few bytes (chromosome names are short).
          val colonIdx = decrypted.indexOf(':'.toByte)
          if colonIdx > 0 then
            val chr = new String(decrypted, 0, colonIdx, StandardCharsets.UTF_8)
            // Find the end of the length digits — scan until a non-digit byte
            var lengthEnd = colonIdx + 1
            while lengthEnd < decrypted.length && decrypted(lengthEnd) >= '0' && decrypted(lengthEnd) <= '9' do
              lengthEnd += 1
            val dataLength = new String(decrypted, colonIdx + 1, lengthEnd - colonIdx - 1, StandardCharsets.UTF_8).toInt
            val headerLength = lengthEnd
            val compactBytes = java.util.Arrays.copyOfRange(decrypted, headerLength, headerLength + dataLength)
            receivedChromosomes.put(chr, compactBytes)
            log.debug(s"Received partner chr $chr variants ($dataLength bytes)")
        payload = pendingVariantPayloads.poll()
    }

  private def handleSummaryHash(payload: EncryptedPayload): Unit =
    sharedSecretForReceive.foreach { key =>
      val decrypted = IbdCryptoService.decryptPayload(payload, key)
      val hash = new String(decrypted, StandardCharsets.UTF_8)
      partnerHash.set(Some(hash))
      partnerHashReceived.countDown()
      log.debug(s"Received partner summary hash: ${hash.take(16)}...")
    }

  private def handleAttestation(payload: EncryptedPayload): Unit =
    sharedSecretForReceive.foreach { key =>
      val decrypted = IbdCryptoService.decryptPayload(payload, key)
      val json = new String(decrypted, StandardCharsets.UTF_8)
      import io.circe.parser.decode as jsonDecode
      jsonDecode[IbdAttestation](json) match
        case Right(attestation) =>
          if IbdAttestation.verify(attestation) then
            partnerAttestationRef.set(Some(attestation))
            partnerAttestationReceived.countDown()
            log.info("Received and verified partner's attestation")
          else
            log.warn("Partner attestation signature verification failed")
        case Left(err) =>
          log.warn(s"Failed to decode partner attestation: ${err.getMessage}")
    }

  /**
   * Derive a deterministic session key from the session ID.
   * Used only for the initial key exchange phase (before ECDH).
   * This is NOT secret — it just provides structure for the relay message format.
   */
  private def deriveSessionKey(sessionId: String): SecretKey =
    val digest = java.security.MessageDigest.getInstance("SHA-256")
    val keyBytes = digest.digest(s"ibd-session-key:$sessionId".getBytes(StandardCharsets.UTF_8))
    javax.crypto.spec.SecretKeySpec(keyBytes, "AES")

  private def awaitFuture[T](f: Future[T], description: String)(implicit ec: ExecutionContext): T =
    import scala.concurrent.Await
    import scala.concurrent.duration.*
    Try(Await.result(f, config.timeoutSeconds.seconds)) match
      case Success(result) => result
      case Failure(e) =>
        throw RuntimeException(s"Failed to $description: ${e.getMessage}", e)

  private def cleanup()(implicit ec: ExecutionContext = ExecutionContext.global): Unit =
    relayClient.foreach { relay =>
      Try(relay.disconnect())
      relayClient = None
    }
