# Navigator OAuth client (native/public)

Navigator authenticates to a user's PDS with AT Protocol OAuth as a **native/public** client:
it catches the redirect on a `http://127.0.0.1:{random-port}/callback` loopback server and uses
**PKCE only** — no client assertion, because a desktop app can't safely hold a signing key or
receive a server-side redirect.

This is a *different* client from the decoding-us **web** app. The web client
(`/oauth/client-metadata.json`) is confidential: `token_endpoint_auth_method: private_key_jwt`,
`application_type: web`, and a single hosted redirect (`/oauth/callback`). Pointing the desktop
app at it breaks login — the auth server would demand a client assertion Navigator can't produce
and would reject the loopback redirect URI.

## What must be published

Serve [`navigator-client-metadata.json`](./navigator-client-metadata.json) verbatim at:

    https://decoding-us.org/oauth/navigator-client-metadata.json

Key fields (must not drift from what the app sends):

- `client_id` — the document's own URL. Navigator sends exactly this as `client_id`.
- `token_endpoint_auth_method: "none"` — public client (PKCE, no assertion).
- `application_type: "native"`.
- `redirect_uris: ["http://127.0.0.1/callback"]` — loopback with **no port**. Per RFC 8252 §7.3
  the authorization server matches loopback (`127.0.0.1`) redirects ignoring the runtime port, so
  the app's `http://127.0.0.1:{port}/callback` is accepted. The path (`/callback`) must match.
- `scope: "atproto transition:generic"` — identity **and** transitional write access. Navigator
  publishes federated records to the user's PDS, which needs write scope; identity-only `atproto`
  is not enough. The app requests this exact string ([`OAUTH_SCOPE`](../../crates/navigator-app/src/lib.rs)).
- `dpop_bound_access_tokens: true`.

No `jwks_uri` — a public client has no signing key.

## How the app resolves the client

`crates/navigator-app/src/lib.rs`:

- `DEFAULT_OAUTH_CLIENT_ID` = the URL above (the default when nothing is overridden).
- `OAUTH_SCOPE` = `atproto transition:generic`.
- `resolve_oauth_config()` / env `DECODINGUS_OAUTH_CLIENT_ID`:
  - unset → hosted default (production).
  - `loopback` → the atproto dev loopback client (for logging in against a **local/test PDS** that
    hasn't registered the production document).
  - any other value → treated as a hosted client-metadata URL (e.g. a local dev document).

## Verifying the cutover

Once the document is live:

1. `curl -s https://decoding-us.org/oauth/navigator-client-metadata.json | jq .` — confirm it
   matches this file (especially `token_endpoint_auth_method`, `redirect_uris`, `scope`).
2. In Navigator, sign in with a real handle → browser authorize → loopback callback should
   complete and persist a DPoP-bound session.
3. Confirm a write works (e.g. publish a coverage/ancestry record) — this exercises
   `transition:generic`; identity-only scope would fail here.

If the loopback-port match is rejected by the entryway, fall back to a custom-scheme redirect
(`com.decoding-us.navigator:/callback`) in both the document and `oauth::redirect_uri`.
