package com.decodingus.ibd.protocol

import munit.FunSuite

class MatchingProtocolStateSpec extends FunSuite:

  test("terminal states are correctly identified") {
    assert(MatchingProtocolState.Completed.isTerminal)
    assert(MatchingProtocolState.Failed("error").isTerminal)
    assert(MatchingProtocolState.TimedOut.isTerminal)

    assert(!MatchingProtocolState.Idle.isTerminal)
    assert(!MatchingProtocolState.Connecting.isTerminal)
    assert(!MatchingProtocolState.ExchangingKeys.isTerminal)
    assert(!MatchingProtocolState.ExtractingVariants.isTerminal)
    assert(!MatchingProtocolState.ComputingIbd.isTerminal)
    assert(!MatchingProtocolState.VerifyingHashes.isTerminal)
  }

  test("all states have descriptions") {
    val allStates = List(
      MatchingProtocolState.Idle,
      MatchingProtocolState.Connecting,
      MatchingProtocolState.ExchangingKeys,
      MatchingProtocolState.ExtractingVariants,
      MatchingProtocolState.SendingVariants,
      MatchingProtocolState.ReceivingVariants,
      MatchingProtocolState.ComputingIbd,
      MatchingProtocolState.ExchangingHashes,
      MatchingProtocolState.VerifyingHashes,
      MatchingProtocolState.SigningAttestation,
      MatchingProtocolState.SubmittingAttestation,
      MatchingProtocolState.PersistingResult,
      MatchingProtocolState.Completed,
      MatchingProtocolState.TimedOut,
      MatchingProtocolState.Failed("test error")
    )
    for state <- allStates do
      assert(state.description.nonEmpty, s"State $state has empty description")

    assertEquals(MatchingProtocolState.Failed("test error").description, "Failed: test error")
  }

  test("progress fractions are monotonically increasing for linear states") {
    val linearStates = List(
      MatchingProtocolState.Idle,
      MatchingProtocolState.Connecting,
      MatchingProtocolState.ExchangingKeys,
      MatchingProtocolState.ExtractingVariants,
      MatchingProtocolState.SendingVariants,
      MatchingProtocolState.ReceivingVariants,
      MatchingProtocolState.ComputingIbd,
      MatchingProtocolState.ExchangingHashes,
      MatchingProtocolState.VerifyingHashes,
      MatchingProtocolState.SigningAttestation,
      MatchingProtocolState.SubmittingAttestation,
      MatchingProtocolState.PersistingResult,
      MatchingProtocolState.Completed
    )

    for i <- 0 until linearStates.length - 1 do
      assert(
        linearStates(i).progressFraction < linearStates(i + 1).progressFraction,
        s"${linearStates(i)} (${linearStates(i).progressFraction}) should be < ${linearStates(i + 1)} (${linearStates(i + 1).progressFraction})"
      )
  }

  test("progress starts at 0 and ends at 1") {
    assertEquals(MatchingProtocolState.Idle.progressFraction, 0.0)
    assertEquals(MatchingProtocolState.Completed.progressFraction, 1.0)
  }

  test("Failed state preserves reason") {
    val failed = MatchingProtocolState.Failed("hash mismatch")
    assertEquals(failed.description, "Failed: hash mismatch")
    assert(failed.isTerminal)
  }
