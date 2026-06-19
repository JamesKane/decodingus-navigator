# Authentication & Permissions (OAuth)

**Status:** Design — supersedes app-password auth and the REST/Kafka relay path
**Spec reference:** https://atproto.com/specs/permission

This document records how the landing of AT Protocol **OAuth + granular permissions**
changes the Atmosphere integration design, and what existing scaffolding it lets us
remove.

---

## 1. What the permission spec does — and does not — cover

The AT Protocol permission spec governs **write-side authorization** only. It defines
five resource types a client can request scoped access to:

| Resource | Grants |
|:---|:---|
| `repo` | Writing records to the user's repo, scoped by collection (NSID) and operation (create/update/delete) |
| `rpc` | Authenticated calls to remote services (service tokens / proxied requests) |
| `blob` | Uploading media, constrained by MIME type |
| `account` | Hosting details (email, repo import) |
| `identity` | DID document and handle |

**It explicitly does NOT cover:**

- **Read permissions** for repository records
- **Subscription / firehose access** (`com.atproto.sync.subscribeRepos`)

> Consequence: OAuth does **not** replace the firehose. The AppView's read/ingest
> path stays. See §4.

---

## 2. Why this matters: it removes our credential-holding intermediary

Our phased plan (REST → Kafka → direct-PDS) existed almost entirely to work around the
fact that the only auth mechanism — **app passwords** — grants *full account access*.
Handing a full-access credential to an Edge client (Navigator) was unacceptable, so
writes were funneled through a trusted backend (Nexus BGS Node) that held the
credentials and pushed via REST, with Kafka planned as the scaled version of the same
intermediary.

OAuth permission sets eliminate that reason to exist:

- The Navigator client requests a **narrow scope** (e.g. create/update/delete on
  `com.decodingus.atmosphere.*` collections) and writes **directly to the user's PDS**.
- No backend ever holds the user's account credentials.
- The user sees exactly which record types the app can touch, and can revoke.

**Net effect:** Phase 1 (REST relay) and Phase 2 (Kafka) collapse into a single
direct-to-PDS write path gated by OAuth scopes. We jump straight to the Phase 3
architecture.

---

## 3. Permission set for Navigator

Rather than presenting the user a wall of individual scopes, we publish a
**permission set** Lexicon bundling the writes Navigator needs:

```
com.decodingus.atmosphere.navigatorCore   (permission set)
  repo:com.decodingus.atmosphere.biosample?action=create,update,delete
  repo:com.decodingus.atmosphere.sequencerun?action=create,update,delete
  repo:com.decodingus.atmosphere.alignment?action=create,update,delete
  repo:com.decodingus.atmosphere.project?action=create,update,delete
  repo:com.decodingus.atmosphere.workspace?action=create,update,delete
  repo:com.decodingus.atmosphere.genotype?action=create,update,delete
  repo:com.decodingus.atmosphere.populationBreakdown?action=create,update,delete
  repo:com.decodingus.atmosphere.haplogroupReconciliation?action=create,update,delete
```

Notes:
- **Namespace authority:** a set under `com.decodingus.*` may only grant permissions
  within `com.decodingus.*` and children — never sibling/parent namespaces. Our records
  live entirely under `com.decodingus.atmosphere.*`, so this fits cleanly.
- **Temporal dynamism:** new record types added to the set later (e.g. `strProfile`,
  `matchConsent`) are picked up on re-resolution without forcing re-authentication.
  Caches expire on the spec's recommended 24h-stale / 90d-hard schedule.
- Wildcards (`*`) are allowed in raw scope strings but **not inside a permission set**,
  so each collection is enumerated explicitly.

---

## 4. The firehose stays (read/ingest path)

Because read and subscription access are out of the permission spec's scope, the AppView
still ingests records via `com.atproto.sync.subscribeRepos`. See
[08-AppView-Lifecycle.md](./08-AppView-Lifecycle.md).

Two clarifications this re-evaluation produces:

1. **Consume the standard relay, don't build custom firehose infrastructure.** The
   AppView subscribes to a standard relay (or Jetstream for filtered JSON) for our
   collection NSIDs. There is no "custom firehose" to maintain — that label applied to
   the REST/Kafka relay (§2), which we are removing.
2. **There is no per-user read gate.** Since the spec offers no read-consent mechanism,
   records the AppView indexes are effectively public-by-design once written to the PDS.
   This is consistent with our core principle (only computed summaries/metadata leave
   the device; raw genomic data never does). Anything that must remain private must
   simply **not be written as a public PDS record** — privacy is enforced at write time,
   not by read permissions.

---

## 5. Backfeed writes (AppView → user PDS) now need an explicit grant

The "backfeed" records (`haplogroupUpdate`, `branchDiscovery`, `matchList`,
`haplogroupAncestralStr`, etc.) involve **DecodingUs writing into a user's repo**. Under
OAuth this is no longer "the AppView holds the token." It requires either:

- A user-granted `rpc` / service scope authorizing the AppView to act on their behalf, or
- The AppView writing to **its own** repo and the user's client subscribing/reading, so
  the user's repo is never mutated by a third party.

Recommendation: prefer the second model where possible (AppView-owned records the client
reads), and reserve user-repo writes for cases that genuinely require them, behind an
explicit grant.

---

## 6. Client type — Navigator is a public/native client (NOT the web confidential client)

The decodingus web app registers as a **confidential client** (`private_key_jwt` ES256,
hosted `client-metadata.json` at the site URL, `/oauth/callback`, signed-cookie session) —
see `decodingus/rust/docs/atproto-oauth-findings.md`. **Navigator is a different client and
must register separately.** As a desktop app it cannot safely hold a signing key, so:

| Aspect | decodingus web (confidential) | Navigator desktop (public/native) |
|:---|:---|:---|
| Client auth | `private_key_jwt` (ES256) | **PKCE only** (no client assertion) |
| Token binding | DPoP | **DPoP** (same) |
| `client_id` / metadata | hosted at site URL | its **own** native `client-metadata.json` |
| Redirect URI | `https://…/oauth/callback` | **loopback** `http://127.0.0.1:<port>/callback` (or claimed/custom scheme) |
| Token storage | signed cookie | **OS keychain** |
| Refresh lifetime | longer | typically **shorter** (public clients) |
| Scope | **read** the genomic collections | **write** (`navigatorCore`, §3) |

> Note the **two-clients/two-scopes** split: the AppView requests *read* on the collections
> (for private notify-fetch; public records arrive via firehose regardless), while Navigator
> requests *write*. They are distinct grants the user approves separately.

**`du-atproto` reuse is partial:** the crypto primitives — PKCE (S256), DPoP proofs,
DID/handle resolution, PDS discovery — reuse directly. The confidential-client pieces
(`private_key_jwt` client assertion, served client-metadata/JWKS, cookie session) do **not**
apply; confirm `du-atproto`'s token-exchange builder supports PKCE-without-client-assertion.

---

## 7. Open joint points (with the decodingus/Edge team)

From `atproto-oauth-findings.md`, the items that are Navigator's to settle:

- **Permission set (their #3):** confirm `navigatorCore` (§3) as the *write* scope, and the
  separate *read* scope the AppView requests. Edge team enumerates the collections.
- **DID method (their #6):** confirm edge accounts use `did:plc` (via `plc.directory`) vs
  `did:web`; `du-atproto` resolves both.
- **Notify-fetch for private data (their #7):** Navigator is the *producer* and may be
  offline → favor **push-on-write from Navigator to an AppView notification endpoint**, with
  the AppView then fetching the record from the PDS under its read scope. **Design-deferred**
  — the upstream group-private spec is still maturing.

---

## 8. Migration checklist

- [ ] Replace `AuthenticationService.loginAtProto` (app-password `createSession`) with the
      OAuth authorization-code flow.
- [ ] Implement **DPoP**-bound access tokens and **refresh-token rotation** (today the
      `refreshJwt` is received and discarded).
- [ ] Register Navigator as a **public/native client** with its own `client-metadata.json`
      and a loopback redirect; store tokens in the OS keychain (§6).
- [ ] Publish the `com.decodingus.atmosphere.navigatorCore` permission set Lexicon.
- [ ] Point `AsyncSyncService.pushCreate/Update/Delete` (currently stubbed) at direct
      PDS writes using the OAuth session.
- [ ] Remove the Kafka integration plan; demote the REST relay to a legacy/bootstrap path.
- [ ] Rework the backfeed design per §5.
- [ ] Keep the firehose ingest path; switch any "custom relay" language to standard
      relay/Jetstream consumption.
