# AT Protocol federation — Edge/Navigator reply

**Re:** `decodingus/rust/docs/atproto-oauth-findings.md`
**From:** Navigator (Edge) team
**Date:** 2026-06-01
**Related:** [11-Auth-and-Permissions.md](./11-Auth-and-Permissions.md), `DUNavigator/documents/design/RustRewrite_Plan.md`

Thanks — the design pivot (drop the custom firehose; OAuth-scoped PDS access +
notify/fetch) matches our side exactly. One correction up front that shapes several
answers below, then point-by-point responses and our crate-reuse plan.

## Up front: Navigator is a *public/native* client, not a confidential web client

The findings doc describes decodingus as a confidential web client (`private_key_jwt`,
hosted `client-metadata.json`, cookie session) and asks (your #1) whether the PDS requires
a public client. For **Navigator the answer is yes** — it's a desktop app, so it registers
**separately** as a public/native client. This is additive, not a conflict: decodingus
(server/AppView) stays confidential; Navigator (desktop) is public. **Two clients, two
scopes** (see #3).

---

## Point-by-point

**#1 — Client auth method.** decodingus web: keep `private_key_jwt` (ES256) + DPoP.
Navigator desktop: **PKCE (S256) only, no client assertion**, DPoP-bound tokens, **loopback
redirect** (`http://127.0.0.1:<port>/callback`), tokens in the **OS keychain**, and we
expect a shorter public-client refresh lifetime. We need the auth server to support a
public client with PKCE + DPoP + loopback redirect.

**#2 — Hosting / registration.** Your `client_id` at the production base URL is your call;
we'll confirm it once set. Navigator hosts its **own** native `client-metadata.json`
(separate `client_id`) — likely a static file under the same domain (e.g.
`https://decoding-us.com/navigator/client-metadata.json`). Can you host a Navigator client
metadata file alongside the web one? If so we'll supply its contents.

**#3 — Scopes / permission sets.** Simplified by the AppView **scope reduction** (the
AppView no longer mirrors the network — see 08-AppView-Lifecycle.md):
- **Navigator → write.** The permission set `com.decodingus.atmosphere.navigatorCore`:
  create/update/delete on `com.decodingus.atmosphere.{biosample, sequencerun, alignment,
  project, workspace, genotype, populationBreakdown, haplogroupReconciliation}` (full
  enumeration in 11-Auth-and-Permissions.md §3). This now **also covers Navigator publishing
  per-sample coverage summaries** as public records for the AppView to aggregate.
- **AppView → read: likely none for now.** Public coverage summaries need no read scope (no
  read-consent gate in the spec); the variant catalog is fed by **direct submission** to an
  AppView curation API (Navigator-authenticated), not a PDS read. A broad AppView read scope
  only resurfaces if/when **private** match data uses notify-fetch — deferred (#7).

We (Edge) own defining/publishing the `navigatorCore` set lexicon. The AppView likely needs
no PDS read scope at this stage — please confirm.

**#4 — Signing key lifecycle.** Applies to your confidential web client (persist
`OAUTH_EC_KEY`, rotate via JWKS `kid`). Navigator as a public client has **no long-lived
signing key** to manage — PKCE is per-flow. So nothing for us here beyond keychain token
handling.

**#5 — DPoP nonce.** Your single-retry-on-`DPoP-Nonce` behavior is what we'll implement too
(reusing `du-atproto`'s DPoP proof builder). We'll match whatever the auth server's nonce
behavior turns out to be at PAR + token endpoints — please share once known.

**#6 — Identity resolution.** Confirm edge accounts: `did:plc` (resolved via
`plc.directory` — self-hosted PLC?) or `did:web`? `du-atproto` resolves both and does
handle→DID via HTTPS well-known. We're fine either way; just need to know for test fixtures.

**#7 — Private data / notify-fetch.** Largely moot for the two surviving AppView flows:
- **Variants** arrive by **direct proposal submission** to an AppView curation API
  (Navigator-authenticated) — no firehose, no notify-fetch.
- **Coverage** is **aggregated on demand** from public summary records. Open mechanism: how
  the AppView *discovers* which public summaries exist for a cohort — a lightweight
  firehose-derived **URI index** (pointers only, not a mirror) or a query against a
  relay/AppView-of-record. This is the only residual firehose use, and it's discovery-only.

Genuine notify-fetch only returns for **private** match data (matchConsent/matchList/IBD).
For that, since Navigator is the producer and may be offline, we'd favor **push-on-write
from Navigator** over AppView polling — but it's **design-deferred** until the upstream
group-private spec settles.

---

## Crate-reuse plan

Navigator is being rewritten in Rust (egui desktop; see `RustRewrite_Plan.md`) and intends
to **build on your crates rather than fork**. Proposed approach: **extract the genuinely
shared crates** out of the decodingus repo into a shared location (own repo or path/git
deps) so both apps depend on them and fixes flow both ways.

| Crate | Navigator reuse | Notes |
|:---|:---|:---|
| `du-domain` | **Full** | Pure types/enums/IDs, zero IO. Also fixes our type-duplication debt. |
| `du-atproto` | **Partial** | Reuse PKCE, DPoP proofs, DID/handle resolution, PDS discovery. *Not* the confidential-client pieces (`private_key_jwt` assertion, served metadata/JWKS, cookie session). **Ask:** confirm the token-exchange builder can run PKCE-without-client-assertion (public client) — small add if not. |
| `du-bio` | **Full + extend** | noodles I/O + liftover + callable. Navigator will add a pure-Rust **haploid variant caller** (force-call + de-novo Y/MT for private variants/branch creation) — candidate to live here so the AppView can reuse it too. |
| `du-db` | **Pattern only** | We use **SQLite** (desktop), not Postgres; reuse the query-module + JSONB→domain mapping pattern, not the crate. |
| `du-external`, `du-web`, `du-jobs` | No | Server/web-specific. |

**Asks of decodingus team:**
1. Agree to extract `du-domain`, `du-atproto`, `du-bio` to a shared crate location, and
   where it lives (new `decodingus-shared` repo / submodule / published).
2. Confirm `du-atproto` token exchange supports the public-client (PKCE-only) path.
3. Decide whether the Navigator haploid caller lands in `du-bio` (shared) or a
   Navigator-only crate.

---

## What we need from you to test end-to-end

Same as your closing ask, plus the public-client bits:
- A **test PDS + account** (handle + DID) and its auth-server endpoints.
- Confirmation the auth server accepts a **public client (PKCE + DPoP + loopback
  redirect)** — so `navigator login --handle <test>` can complete a real PAR → redirect →
  token flow.
- The **read scope** string the AppView will request (#3) and the **DID method** edge
  accounts use (#6).
