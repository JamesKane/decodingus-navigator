package com.decodingus.ibd.crypto

import java.security.*
import java.security.spec.{NamedParameterSpec, X509EncodedKeySpec}
import java.util.Base64
import javax.crypto.spec.{GCMParameterSpec, SecretKeySpec}
import javax.crypto.{Cipher, KeyAgreement, SecretKey}

/**
 * Cryptographic operations for IBD matching protocol.
 *
 * Uses JDK 17 built-in APIs exclusively — no external crypto dependencies:
 * - X25519 ECDH for ephemeral key exchange (JEP 324, JDK 11+)
 * - AES-256-GCM for authenticated encryption
 * - Ed25519 for attestation signatures (JEP 339, JDK 15+)
 */
object IbdCryptoService:

  private val GCM_IV_LENGTH = 12
  private val GCM_TAG_LENGTH = 128 // bits

  // ============================================
  // X25519 Key Exchange
  // ============================================

  /**
   * Generate an ephemeral X25519 key pair for ECDH key exchange.
   * Each IBD comparison session should use a fresh key pair.
   */
  def generateX25519KeyPair(): KeyPair =
    val kpg = KeyPairGenerator.getInstance("X25519")
    kpg.generateKeyPair()

  /**
   * Derive a shared AES-256 secret from local private key and remote public key.
   * Uses X25519 ECDH followed by SHA-256 key derivation.
   */
  def deriveSharedSecret(myPrivateKey: PrivateKey, theirPublicKey: PublicKey): SecretKey =
    val ka = KeyAgreement.getInstance("XDH")
    ka.init(myPrivateKey)
    ka.doPhase(theirPublicKey, true)
    val sharedSecret = ka.generateSecret()

    // Derive AES-256 key via SHA-256 hash of the raw shared secret
    val digest = MessageDigest.getInstance("SHA-256")
    val keyBytes = digest.digest(sharedSecret)
    SecretKeySpec(keyBytes, "AES")

  // ============================================
  // AES-256-GCM Encryption
  // ============================================

  /**
   * Encrypt data using AES-256-GCM with a random IV.
   * Returns the ciphertext (which includes the GCM auth tag appended by the JDK implementation).
   */
  def encrypt(data: Array[Byte], key: SecretKey): (Array[Byte], Array[Byte]) =
    val iv = new Array[Byte](GCM_IV_LENGTH)
    SecureRandom().nextBytes(iv)

    val cipher = Cipher.getInstance("AES/GCM/NoPadding")
    cipher.init(Cipher.ENCRYPT_MODE, key, GCMParameterSpec(GCM_TAG_LENGTH, iv))
    val ciphertext = cipher.doFinal(data)
    (ciphertext, iv)

  /**
   * Decrypt AES-256-GCM ciphertext.
   * The ciphertext must include the appended GCM auth tag (JDK default behavior).
   */
  def decrypt(ciphertext: Array[Byte], iv: Array[Byte], key: SecretKey): Array[Byte] =
    val cipher = Cipher.getInstance("AES/GCM/NoPadding")
    cipher.init(Cipher.DECRYPT_MODE, key, GCMParameterSpec(GCM_TAG_LENGTH, iv))
    cipher.doFinal(ciphertext)

  /**
   * Encrypt data and package as an EncryptedPayload for relay transport.
   */
  def encryptToPayload(data: Array[Byte], key: SecretKey, sessionId: String,
                       dataType: String, senderKeyId: Option[String] = None): EncryptedPayload =
    val (ciphertext, iv) = encrypt(data, key)
    val encoder = Base64.getEncoder
    EncryptedPayload(
      sessionId = sessionId,
      encryptedData = encoder.encodeToString(ciphertext),
      iv = encoder.encodeToString(iv),
      dataType = dataType,
      senderKeyId = senderKeyId
    )

  /**
   * Decrypt an EncryptedPayload received from the relay.
   */
  def decryptPayload(payload: EncryptedPayload, key: SecretKey): Array[Byte] =
    val decoder = Base64.getDecoder
    val ciphertext = decoder.decode(payload.encryptedData)
    val iv = decoder.decode(payload.iv)
    decrypt(ciphertext, iv, key)

  // ============================================
  // Ed25519 Signatures
  // ============================================

  /**
   * Generate an Ed25519 signing key pair for attestation signatures.
   * This key pair should be long-lived (tied to the user's identity).
   */
  def generateEd25519KeyPair(): KeyPair =
    val kpg = KeyPairGenerator.getInstance("Ed25519")
    kpg.generateKeyPair()

  /**
   * Sign data with an Ed25519 private key.
   * Used for attestation signing — proving that this Navigator computed the match.
   */
  def signAttestation(data: Array[Byte], signingKey: PrivateKey): Array[Byte] =
    val sig = Signature.getInstance("Ed25519")
    sig.initSign(signingKey)
    sig.update(data)
    sig.sign()

  /**
   * Verify an Ed25519 signature.
   * Used to verify the other party's attestation.
   */
  def verifyAttestation(data: Array[Byte], signature: Array[Byte], publicKey: PublicKey): Boolean =
    try
      val sig = Signature.getInstance("Ed25519")
      sig.initVerify(publicKey)
      sig.update(data)
      sig.verify(signature)
    catch
      case _: SignatureException => false

  // ============================================
  // Key Serialization (for AT Protocol transport)
  // ============================================

  /**
   * Encode a public key to Base64 (X.509 SubjectPublicKeyInfo format).
   * Compatible with standard key exchange in DID documents and AT Protocol.
   */
  def encodePublicKey(key: PublicKey): String =
    Base64.getEncoder.encodeToString(key.getEncoded)

  /**
   * Decode an X25519 public key from Base64 (X.509 SubjectPublicKeyInfo format).
   */
  def decodeX25519PublicKey(base64: String): PublicKey =
    val keyBytes = Base64.getDecoder.decode(base64)
    val keySpec = X509EncodedKeySpec(keyBytes)
    val kf = KeyFactory.getInstance("X25519")
    kf.generatePublic(keySpec)

  /**
   * Decode an Ed25519 public key from Base64 (X.509 SubjectPublicKeyInfo format).
   */
  def decodeEd25519PublicKey(base64: String): PublicKey =
    val keyBytes = Base64.getDecoder.decode(base64)
    val keySpec = X509EncodedKeySpec(keyBytes)
    val kf = KeyFactory.getInstance("Ed25519")
    kf.generatePublic(keySpec)

  /**
   * Compute SHA-256 hash of data, returned as hex string.
   * Used for attestation hashes to compare IBD computation results.
   */
  def sha256Hex(data: Array[Byte]): String =
    val digest = MessageDigest.getInstance("SHA-256")
    digest.digest(data).map("%02x".format(_)).mkString
