package com.decodingus.ibd.protocol

/**
 * States of the IBD matching protocol state machine.
 *
 * The protocol proceeds linearly through these phases.
 * Any error transitions to Failed; partner disconnect
 * or timeout transitions to TimedOut.
 */
enum MatchingProtocolState:
  case Idle
  case Connecting
  case ExchangingKeys
  case ExtractingVariants
  case SendingVariants
  case ReceivingVariants
  case ComputingIbd
  case ExchangingHashes
  case VerifyingHashes
  case SigningAttestation
  case SubmittingAttestation
  case PersistingResult
  case Completed
  case Failed(reason: String)
  case TimedOut

object MatchingProtocolState:
  extension (s: MatchingProtocolState)
    def isTerminal: Boolean = s match
      case Completed | Failed(_) | TimedOut => true
      case _ => false

    def description: String = s match
      case Idle => "Ready"
      case Connecting => "Connecting to relay..."
      case ExchangingKeys => "Exchanging encryption keys..."
      case ExtractingVariants => "Extracting variant data..."
      case SendingVariants => "Sending encrypted variants..."
      case ReceivingVariants => "Waiting for partner's variants..."
      case ComputingIbd => "Computing IBD segments..."
      case ExchangingHashes => "Exchanging match hashes..."
      case VerifyingHashes => "Verifying match agreement..."
      case SigningAttestation => "Signing attestation..."
      case SubmittingAttestation => "Submitting attestation to AppView..."
      case PersistingResult => "Saving match result..."
      case Completed => "Match completed"
      case Failed(reason) => s"Failed: $reason"
      case TimedOut => "Timed out waiting for partner"

    def progressFraction: Double = s match
      case Idle => 0.0
      case Connecting => 0.05
      case ExchangingKeys => 0.10
      case ExtractingVariants => 0.20
      case SendingVariants => 0.35
      case ReceivingVariants => 0.45
      case ComputingIbd => 0.60
      case ExchangingHashes => 0.75
      case VerifyingHashes => 0.80
      case SigningAttestation => 0.85
      case SubmittingAttestation => 0.90
      case PersistingResult => 0.95
      case Completed => 1.0
      case Failed(_) => 1.0
      case TimedOut => 1.0
