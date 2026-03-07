# IBD Matching System — Navigator Implementation Plan

Last updated: 2026-03-07

## Overview

Implements backlog item 2.1 (IBD Matching System) in Navigator, coordinating with the
DecodingUs AppView for discovery, consent tracking, and attestation indexing.

**Architecture:** Dual-computation mutual-attestation model.
- Navigator A and Navigator B each compute IBD segments independently
- Results are hashed and compared; matching hashes confirm agreement
- Both sign attestations; AppView "stamps" the match when both attestations arrive

**Key constraint:** Raw variant data never leaves the local machine. Only match summaries
(total cM, segment count, largest segment) and cryptographic attestations are shared.

---

## Library Decisions

### Cryptography — Java 17 Built-in (zero new dependencies)

All three crypto needs are covered by JDK 17's standard library:

| Need | JDK API | Since |
|------|---------|-------|
| X25519 ECDH key exchange | `KeyPairGenerator.getInstance("X25519")`, `KeyAgreement.getInstance("XDH")` | JDK 11 (JEP 324) |
| AES-256-GCM encryption | `Cipher.getInstance("AES/GCM/NoPadding")`, `GCMParameterSpec` | JDK 11 |
| Ed25519 signatures | `KeyPairGenerator.getInstance("Ed25519")`, `Signature.getInstance("Ed25519")` | JDK 15 (JEP 339) |

No Bouncy Castle, Tink, or lazysodium needed. This avoids dependency bloat, native
library deployment issues (lazysodium), and opinionated keyset formats (Tink) that
would complicate AT Protocol interoperability.

### P2P Communication — STTP WebSocket (already in project)

STTP client3 3.9.7 (already a dependency) supports WebSockets via `HttpClientFutureBackend`
(Java 11+'s `java.net.http.HttpClient`). No new client dependency needed.

The communication model uses a **WebSocket relay** hosted by the AppView, not direct
peer-to-peer. This avoids NAT traversal complexity for a desktop application. The relay
sees only encrypted payloads — it cannot read the content.

### IBD Segment Detection — Custom using HTSJDK (already in project)

No IBD detection library exists on Maven Central. hap-IBD (Browning Lab) is Apache 2.0
pure Java but distributed only as source — not published to any Maven repository.
GATK, Picard, Hail, KING, PLINK, and BioJava were all evaluated; none provide
IBD segment detection as a Maven dependency.

**Approach:** Build a custom pairwise IBS/IBD detector using HTSJDK 3.0.5 (already
in `build.sbt`). For comparing two individuals, the algorithm is O(n) in variant count:

1. Read VCF genotypes via HTSJDK `VCFFileReader` + `VariantContext.getGenotype()`
2. Classify each position as IBS-0, IBS-1, or IBS-2 (count shared alleles)
3. Sliding window scan along each chromosome, tracking IBS-2 fraction
4. Mark runs of high-IBS sites exceeding minimum length as candidate segments
5. Convert physical positions to genetic distance (cM) using bundled HapMap
   recombination maps (~30MB for GRCh37+GRCh38)
6. Filter segments below minimum cM threshold
7. Apply LOD scoring with error tolerance for genotyping noise

The core IBS comparison is ~50 lines. The complexity is in segment boundary detection
(error tolerance, gap handling) and genetic map interpolation. Both are well-understood
algorithms with published reference implementations.

### Summary — New Dependencies

| Dependency | Purpose | Required? |
|-----------|---------|-----------|
| None (JDK crypto) | ECDH, AES-GCM, Ed25519 | Already available |
| None (STTP) | WebSocket client | Already in build.sbt |
| None (HTSJDK) | VCF reading for IBD detector | Already in build.sbt |

---

## Consent Model Reconciliation

The Atmosphere Lexicon (`06-IBD-Matching-Records.md`) and the planning doc
(`ibd-matching-system.md`) define `matchConsent` differently:

| Spec | Meaning | Purpose |
|------|---------|---------|
| Lexicon | Global opt-in toggle (presence = willing to be matched) | Gate for discovery |
| Planning | Per-request vote (ACCEPT/REJECT/DEFER) | Gate for IBD analysis |

**Resolution:** Both are needed. Implement as two separate concepts:

1. **`MatchConsent` (global)** — Atmosphere Lexicon record synced to PDS.
   Presence enables this biosample for match discovery. Deletion revokes.
   Maps to `match_consent` local table in Navigator.

2. **`MatchComparisonApproval` (per-request)** — Local-only record tracking
   the user's approval for a specific comparison. Not synced to PDS directly;
   instead, triggers the Edge App coordination protocol when mutual.

The AppView tracks consent records via Firehose and per-request approvals
via its own `match_consent_tracking` table.

---

## Implementation Phases

### Phase 1: Models, Persistence, and Consent Management

**Goal:** Local data model, consent opt-in/out, feature toggle.

**Navigator work:**

1. **Feature toggle** — Add `ibd-matching` toggle to `feature_toggles.conf` (default off)

2. **Domain models** in `com.decodingus.workspace.model`:
   - `MatchConsent` — global opt-in with consentLevel, allowedMatchTypes, minimumSegmentCm
   - `MatchEntry` — a single confirmed match (matchedBiosampleRef, totalSharedCm, segmentCount, longestSegmentCm, relationshipEstimate, sharedSegments)
   - `MatchRequest` — incoming/outgoing request (fromBiosampleRef, toBiosampleRef, status, message, discoveryReason)
   - `IbdSegment` — shared segment (chromosome, startPosition, endPosition, lengthCm, snpCount, isHalfIdentical)

3. **Database migration** — `V012__ibd_matching_tables.sql`:
   ```sql
   CREATE TABLE match_consent (
     id UUID PRIMARY KEY,
     biosample_id UUID NOT NULL,
     consent_level VARCHAR(20) NOT NULL,
     allowed_match_types JSON DEFAULT '["IBD"]',
     minimum_segment_cm DOUBLE DEFAULT 7.0,
     share_contact_info BOOLEAN DEFAULT FALSE,
     consented_at TIMESTAMP NOT NULL,
     expires_at TIMESTAMP,
     at_uri VARCHAR(500),
     at_cid VARCHAR(100),
     at_version INTEGER DEFAULT 1,
     at_status VARCHAR(20) DEFAULT 'active',
     meta JSON NOT NULL
   );

   CREATE TABLE match_request (
     id UUID PRIMARY KEY,
     from_biosample_ref VARCHAR(500) NOT NULL,
     to_biosample_ref VARCHAR(500) NOT NULL,
     status VARCHAR(20) NOT NULL DEFAULT 'PENDING',
     request_type VARCHAR(50) DEFAULT 'AUTOSOMAL',
     message TEXT,
     discovery_reason JSON,
     created_at TIMESTAMP NOT NULL,
     expires_at TIMESTAMP,
     responded_at TIMESTAMP,
     at_uri VARCHAR(500),
     meta JSON NOT NULL
   );

   CREATE TABLE match_result (
     id UUID PRIMARY KEY,
     biosample_id UUID NOT NULL,
     matched_biosample_ref VARCHAR(500) NOT NULL,
     matched_citizen_did VARCHAR(255),
     relationship_estimate VARCHAR(50),
     total_shared_cm DOUBLE NOT NULL,
     longest_segment_cm DOUBLE,
     segment_count INTEGER NOT NULL,
     shared_segments JSON DEFAULT '[]',
     x_match_shared_cm DOUBLE,
     matched_at TIMESTAMP NOT NULL,
     attestation_hash VARCHAR(64),
     at_uri VARCHAR(500),
     meta JSON NOT NULL
   );
   ```

4. **Repositories** — `MatchConsentRepository`, `MatchRequestRepository`, `MatchResultRepository`
   following existing `SyncableRepositoryBase` pattern

5. **Entity conversions** — Add to `EntityConversions.scala`

6. **SyncEntityType** — Add `MatchConsent`, `MatchRequest` cases to enum

7. **PDS validation** — Add `validateMatchConsent` to `PdsSyncValidation`

8. **PDS codecs** — Add Circe encoders/decoders to `PdsClient.scala` for match types

9. **Consent UI** — Add consent toggle to Subject detail view (opt-in/out with consent level picker)

**AppView work:** None for Phase 1 (local-only).

**Tests:** Consent CRUD, validation, entity conversions.

---

### Phase 2: Match Discovery Integration

**Goal:** Receive and display match suggestions from the AppView discovery engine.

**Navigator work:**

1. **DecodingUsClient extensions** — Add methods to call AppView discovery API:
   - `getSuggestions(biosampleGuid, limit)` → `Seq[MatchSuggestion]`
   - `dismissSuggestion(suggestionId)` → `Boolean`
   - `getPopulationOverlap(sampleGuid1, sampleGuid2)` → `Double`

2. **Match suggestion model** — `MatchSuggestion` with reason type, score, matched biosample info

3. **Match request flow** — UI to:
   - Browse suggestions (ranked by score)
   - Send match request (creates `matchRequest` in PDS)
   - View incoming requests
   - Accept/decline incoming requests

4. **SubjectDetailView IBD tab** — Replace placeholder with live suggestion list and
   request management (the placeholder UI structure at lines 744-877 already has the
   right layout; wire it to real data)

**AppView work:** See supplementary backlog (IBD-AV-1, IBD-AV-2).

**Tests:** Client API integration, suggestion display.

---

### Phase 3: Cryptographic Infrastructure

**Goal:** Key exchange and encrypted data exchange capabilities.

**Navigator work:**

1. **`IbdCryptoService`** in `com.decodingus.ibd.crypto`:
   ```
   generateX25519KeyPair(): (PublicKey, PrivateKey)
   deriveSharedSecret(myPrivateKey, theirPublicKey): SecretKey
   encrypt(data: Array[Byte], key: SecretKey): EncryptedPayload
   decrypt(payload: EncryptedPayload, key: SecretKey): Array[Byte]
   signAttestation(data: Array[Byte], signingKey: PrivateKey): Array[Byte]
   verifyAttestation(data: Array[Byte], signature: Array[Byte], publicKey: PublicKey): Boolean
   ```
   All using JDK 17 built-in APIs (X25519, AES-256-GCM, Ed25519).

2. **`EncryptedPayload`** model — sessionId, encryptedData (Base64), iv, authTag, dataType

3. **Key serialization** — Methods to encode/decode public keys for AT Protocol
   transport (Base64 or multibase format, compatible with DID document keys)

4. **`IbdRelayClient`** — STTP WebSocket client connecting to AppView relay:
   - `connect(sessionId, authToken)` → WebSocket session
   - `send(encryptedPayload)` → Unit
   - `onReceive(handler: EncryptedPayload => Unit)` → Unit
   - Reconnection with exponential backoff

**AppView work:** See supplementary backlog (IBD-AV-3: WebSocket relay endpoint).

**Tests:** Round-trip encrypt/decrypt, key exchange between two in-process parties,
signature verify, relay client mock.

---

### Phase 4: IBD Computation Engine

**Goal:** Local IBD segment detection between two individuals.

**Navigator work:**

1. **Variant extraction** — `IbdVariantExtractor` in `com.decodingus.ibd`:
   - Extract biallelic SNP genotypes from VCF/GVCF for autosomal chromosomes
   - Filter to common variant positions (from reference panel)
   - Output: position-sorted genotype array per chromosome

2. **Pairwise IBD detector** — `PairwiseIbdDetector`:
   - Input: two genotype arrays (one local, one received encrypted)
   - Algorithm: sliding window IBS (identity by state) scan with
     genetic map interpolation for cM conversion
   - Output: `List[IbdSegment]` above minimum cM threshold
   - IBS scan with LOD scoring and error-tolerant segment boundaries
   - Uses HTSJDK VCFFileReader for genotype access (already in build.sbt)

3. **Genetic map** — Bundle HapMap GRCh37/GRCh38 recombination maps for
   physical→genetic position conversion (small files, ~30MB total)

4. **Match summary computation**:
   - `totalSharedCm`, `segmentCount`, `longestSegmentCm`, `relationshipEstimate`
   - `matchSummaryHash` — SHA-256 of canonical JSON summary (deterministic)

5. **Relationship estimator** — Map total shared cM to relationship category
   using standard genetic genealogy thresholds (Shared cM Project data)

**AppView work:** None (computation is entirely local).

**Tests:** Known IBD pairs from reference data, hash determinism, relationship estimation.

---

### Phase 5: End-to-End Matching Protocol

**Goal:** Complete the matching flow from mutual consent to confirmed match.

**Navigator work:**

1. **`IbdMatchingCoordinator`** — Orchestrates the full protocol:
   ```
   Step 1: Receive "match ready" notification from AppView
   Step 2: Generate X25519 session keypair
   Step 3: Exchange public keys via relay
   Step 4: Derive shared AES-256 key
   Step 5: Extract and encrypt local variant positions
   Step 6: Send encrypted variants via relay
   Step 7: Receive partner's encrypted variants
   Step 8: Decrypt and run pairwise IBD detection
   Step 9: Compute match summary + hash
   Step 10: Exchange hashes for verification
   Step 11: If hashes match → sign attestation
   Step 12: Submit attestation to AppView
   Step 13: Persist match result locally
   ```

2. **Attestation model** — `IbdAttestation` with matchRequestUri, attestingDid,
   attestingSampleGuid, matchSummary, matchSummaryHash, signature, partnerSummaryHash

3. **Background execution** — Run matching in a background thread with progress
   reporting to the UI (similar to existing analysis processors)

4. **Error handling** — Timeout, partner disconnect, hash mismatch (→ abort + log),
   relay failure (→ retry with backoff)

**AppView work:** See supplementary backlog (IBD-AV-4, IBD-AV-5).

**Tests:** Full protocol mock (two in-process "navigators"), timeout handling,
hash mismatch rejection.

---

### Phase 6: UI — Match List and Chromosome Browser

**Goal:** Display confirmed matches and shared segment visualization.

**Navigator work:**

1. **Match list table** — Wire `ibdMatchesTable` in SubjectDetailView to real
   `MatchResult` data from repository. Columns: name, shared cM, segments,
   longest, relationship estimate.

2. **Chromosome browser** — Replace placeholder with actual visualization:
   - 22 autosomal chromosome ideograms + X
   - Highlighted regions for shared segments
   - Color-coded by segment length (longer = more saturated)
   - Click segment → show detail (start/end positions, cM, SNP count)

3. **Match detail dialog** — On row click, show:
   - Full match summary
   - Shared segments list
   - Chromosome browser focused on this match
   - Contact request button (if consent allows)

4. **Consent management panel** — Settings section for:
   - Enable/disable matching
   - Consent level (FULL, ANONYMOUS, PROJECT_ONLY)
   - Minimum segment cM threshold
   - Allowed match types (IBD, Y_STR, MT_SEQUENCE, AUTOSOMAL)

**Tests:** UI integration tests where feasible.

---

## Phase Dependencies

```
Phase 1 (Models/Consent) ─────────────────────────┐
                                                    │
Phase 2 (Discovery) ──────────────────────────┐    │
   depends on: Phase 1, AppView IBD-AV-1/2    │    │
                                               │    │
Phase 3 (Crypto) ─────────────────────────┐   │    │
   depends on: Phase 1                     │   │    │
                                           │   │    │
Phase 4 (IBD Engine) ────────────────┐    │   │    │
   depends on: nothing (parallel)     │    │   │    │
                                      │    │   │    │
Phase 5 (E2E Protocol) ──────────────┤    │   │    │
   depends on: Phase 3, Phase 4,     │    │   │    │
   AppView IBD-AV-3/4/5              │    │   │    │
                                      │    │   │    │
Phase 6 (UI) ─────────────────────────┘    │   │    │
   depends on: Phase 2, Phase 5            │   │    │
                                           │   │    │
                                     [All complete]
```

Phases 3 and 4 can be developed in parallel. Phase 2 can proceed once the AppView
implements its discovery engine (IBD-AV-1/2).

---

## Navigator Package Structure

```
com.decodingus.ibd/
  crypto/
    IbdCryptoService.scala        — X25519, AES-GCM, Ed25519 operations
    EncryptedPayload.scala        — Encrypted data transport model
  engine/
    IbdVariantExtractor.scala     — Extract genotypes from VCF for comparison
    PairwiseIbdDetector.scala     — IBS/IBD segment detection algorithm
    RelationshipEstimator.scala   — Shared cM → relationship mapping
    GeneticMap.scala              — Physical ↔ genetic position conversion
  protocol/
    IbdMatchingCoordinator.scala  — Orchestrates full matching flow
    IbdRelayClient.scala          — WebSocket client for relay communication
    IbdAttestation.scala          — Attestation model + signing
  service/
    MatchConsentService.scala     — Consent management + PDS sync
    MatchDiscoveryService.scala   — Integration with AppView discovery API
    MatchResultService.scala      — Local match result persistence
```

---

## Existing Infrastructure to Leverage

| What | Where | Use |
|------|-------|-----|
| SubjectDetailView IBD tab | `ui/v2/SubjectDetailView.scala:744-877` | Placeholder UI to wire up |
| 19 i18n keys | `messages.properties:263-282` | Already defined for IBD UI |
| `supportsAutosomalIbd` | `TestTypeDefinition` | Filter eligible biosamples |
| SHA-256 hashing | `AnalysisCache`, `LibraryStats` | Pattern for match summary hashing |
| STTP + Circe | `build.sbt` | HTTP + WebSocket client + JSON |
| SyncableRepositoryBase | `repository/` | Pattern for IBD repositories |
| PdsSyncValidation | `pds/PdsSyncValidation.scala` | Add match consent validation |
| PdsClient codecs | `pds/PdsClient.scala` | Add IBD record codecs |
| Analysis processor pattern | `analysis/` | Background IBD computation |
| Feature toggles | `feature_toggles.conf` | Gate IBD features |
| RecordMeta + AT URI | `workspace/model/` | Standard sync metadata |
