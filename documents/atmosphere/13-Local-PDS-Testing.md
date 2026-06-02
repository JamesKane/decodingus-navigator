# Local PDS for Navigator OAuth dev/test

**From:** decodingus (AppView) team
**Date:** 2026-06-02
**Re:** A self-contained local AT Protocol PDS for developing/testing Navigator's
OAuth client and record writes — no public network, no real account needed.

We stood up the official PDS in a container and validated the live OAuth
**discovery + PAR + DPoP** path against it from `du-atproto`. This doc is the
reproducible setup + the gotchas, so Navigator can drive the same loop (and the
full browser login, with the TLS note in §6).

---

## 1. What you get

- A real AT Protocol authorization server (`@atproto/oauth-provider`) you can hit
  locally: PAR, authorize, token, DPoP, `use_dpop_nonce`.
- Real account creation → a `did:plc` you own.
- Confirmed server capabilities (from `/.well-known/oauth-authorization-server`):
  `code` flow, `S256` PKCE, `token_endpoint_auth_methods = [none, private_key_jwt]`
  (so **public/loopback clients are accepted** — Navigator's path), and
  `client_id_metadata_document_supported = true`. PDS service DID is `did:web:pds.test`.

## 2. Boot the PDS

Image: `ghcr.io/bluesky-social/pds:latest` (multi-arch; arm64 works on Apple
Silicon). Pick one runtime:

**Docker (Navigator devs):**
```sh
mkdir -p ./pdsdata/blocks
docker run -d --name pds -p 3000:3000 -v "$PWD/pdsdata:/pds" \
  -e PDS_HOSTNAME=pds.test -e PDS_PORT=3000 \
  -e PDS_JWT_SECRET=$(openssl rand --hex 16) \
  -e PDS_ADMIN_PASSWORD=$(openssl rand --hex 16) \
  -e PDS_PLC_ROTATION_KEY_K256_PRIVATE_KEY_HEX=$(openssl rand --hex 32) \
  -e PDS_DATA_DIRECTORY=/pds -e PDS_BLOBSTORE_DISK_LOCATION=/pds/blocks \
  -e PDS_DID_PLC_URL=https://plc.directory \
  -e PDS_INVITE_REQUIRED=false -e PDS_DEV_MODE=true \
  ghcr.io/bluesky-social/pds:latest
# reachable at http://localhost:3000
```

**Apple `container` (what we used; each container gets its own IP, no port-map):**
```sh
container image pull ghcr.io/bluesky-social/pds:latest
mkdir -p /tmp/pdsdata/blocks
container run -d --name pds -v /tmp/pdsdata:/pds \
  -e PDS_HOSTNAME=pds.test -e PDS_PORT=3000 \
  -e PDS_JWT_SECRET=$(openssl rand --hex 16) \
  -e PDS_ADMIN_PASSWORD=$(openssl rand --hex 16) \
  -e PDS_PLC_ROTATION_KEY_K256_PRIVATE_KEY_HEX=$(openssl rand --hex 32) \
  -e PDS_DATA_DIRECTORY=/pds -e PDS_BLOBSTORE_DISK_LOCATION=/pds/blocks \
  -e PDS_DID_PLC_URL=https://plc.directory \
  -e PDS_INVITE_REQUIRED=false -e PDS_DEV_MODE=true \
  ghcr.io/bluesky-social/pds:latest
IP=$(container ls | awk '$1=="pds"{print $6}' | cut -d/ -f1)   # e.g. 192.168.64.5
```

Save the generated `PDS_ADMIN_PASSWORD` (echo it / capture from the env) for
`pdsadmin` operations. `container logs pds` / `docker logs pds` for boot output;
`http://<host>/xrpc/_health` returns `{"version":"…"}`.

**Boot gotchas (learned the hard way):**
- `PDS_HOSTNAME` **must not end in `.local`** (PDS rejects it). Use a `.test`
  domain. It also can't be a bare IP.
- The data dir (`/pds`) **must exist** → bind-mount it (above), or the SQLite open
  fails on boot.
- `PDS_INVITE_REQUIRED=false` lets you create accounts without an invite code.

## 3. Create a test account

With invites off, call `createAccount` directly (handle must be under the PDS host):
```sh
curl -s -X POST "http://<host>:3000/xrpc/com.atproto.server.createAccount" \
  -H "Content-Type: application/json" \
  -d '{"handle":"alice.pds.test","email":"alice@example.test","password":"alice-pw-12345"}'
# → { "handle":"alice.pds.test", "did":"did:plc:…", "accessJwt":"…", "refreshJwt":"…" }
```
You now own that `did:plc` and can log in at the PDS's `/oauth/authorize` screen
with the handle + password. (`pdsadmin account create` inside the container also
works if you prefer the admin path.)

## 4. OAuth endpoints

From `http://<host>:3000/.well-known/oauth-authorization-server` (note the
endpoints are advertised under the **canonical** `https://pds.test` host — see §5):
- `pushed_authorization_request_endpoint` → `https://pds.test/oauth/par`
- `authorization_endpoint` → `https://pds.test/oauth/authorize`
- `token_endpoint` → `https://pds.test/oauth/token`

Protected-resource metadata: `http://<host>:3000/.well-known/oauth-protected-resource`.

## 5. Navigator client specifics (public / loopback)

Navigator is a **public native client**: PKCE only (no `client_assertion`), DPoP,
loopback redirect. For local dev use the atproto **loopback client** — no hosted
client-metadata document needed:

- `client_id` = `http://localhost?redirect_uri=<urlenc 127.0.0.1 cb>&scope=atproto`
- `redirect_uri` = `http://127.0.0.1:<port>/callback` (loopback only)
- PKCE `S256`; DPoP on PAR + token.

**Two gotchas that will bite you:**

1. **DPoP `htu` = the server's CANONICAL endpoint**, not the transport URL you
   connect over. Over a local container you POST to `http://<ip>:3000/oauth/par`,
   but the DPoP proof's `htu` must be `https://pds.test/oauth/par` (from metadata).
   Signing the transport URL → `401 invalid_dpop_proof: DPoP "htu" mismatch`.

2. **`use_dpop_nonce` single-retry is mandatory.** The first PAR returns
   `400 { "error":"use_dpop_nonce" }` with a `DPoP-Nonce` response header; re-sign
   the DPoP proof with that nonce and re-POST → `201 Created` + `request_uri`.
   (Same dance on the token endpoint.) We verified this end-to-end.

**Reference harness:** `decodingus-shared/crates/du-atproto/tests/live_pds.rs` does
discovery → public-client PAR → DPoP → nonce-retry against this PDS
(`PDS_TEST_URL=http://<ip>:3000 cargo test -p du-atproto --test live_pds -- --nocapture`).
Reusable `du-atproto` builders: `Pkce`, `dpop_proof`, `par_form_public`,
`token_form_public`.

## 6. Full browser login loop over TLS (verified setup)

Discovery + PAR work over plain `http://<ip>:3000` **if** you sign the canonical
`htu` (§5). The redirect → consent → `code` → token loop additionally needs the
auth server reachable over **HTTPS at its canonical host** (`https://pds.test`) —
the issuer and DPoP `htu` are https-canonical and the browser is redirected to the
advertised `https://pds.test/oauth/authorize`. We stood this up and verified the
handshake end-to-end **up to the consent screen** (the consent click itself is
browser-gated — §6.3). Recipe:

### 6.1 TLS proxy (Caddy internal CA)

```sh
PDS_IP=$(container ls | awk '$1=="pds"{print $6}' | cut -d/ -f1)   # docker: 127.0.0.1
printf '{\n  auto_https disable_redirects\n}\npds.test {\n  tls internal\n  reverse_proxy %s:3000\n}\n' "$PDS_IP" > /tmp/Caddyfile
container run -d --name caddy -v /tmp/Caddyfile:/etc/caddy/Caddyfile docker.io/library/caddy:2
CADDY_IP=$(container ls | awk '$1=="caddy"{print $6}' | cut -d/ -f1)
# Export Caddy's internal root CA so clients can trust https://pds.test:
container exec caddy cat /data/caddy/pki/authorities/local/root.crt > /tmp/caddy_ca.crt
# Verify (no /etc/hosts needed for curl — use --resolve):
curl --resolve pds.test:443:$CADDY_IP --cacert /tmp/caddy_ca.crt \
  https://pds.test/.well-known/oauth-authorization-server
```

(Docker users: `-p 443:443` on Caddy and `127.0.0.1` as `PDS_IP`, then a hosts
entry `127.0.0.1 pds.test`.)

### 6.2 Reaching `https://pds.test` from your client without `/etc/hosts`

For **programmatic** HTTP clients, pin the name → IP rather than editing
`/etc/hosts` (decodingus' Rust client does exactly this via env — useful pattern
for Navigator too):
- trust the CA: load `/tmp/caddy_ca.crt` as an extra root;
- resolve override: map `pds.test` → `$CADDY_IP:443`.

decodingus' reference: `du-web/src/oauth.rs` reads `DU_OAUTH_DEV_CA` (PEM path) and
`DU_OAUTH_DEV_RESOLVE` (`pds.test:<caddy-ip>`) to build a reqwest client that trusts
the CA + pins the host; `DU_OAUTH_DEV_PDS` makes it use the PDS as a fixed auth
server (skipping handle→DID→PDS resolution); the loopback redirect comes from
`DU_OAUTH_LOOPBACK`. Over this, discovery + PAR + DPoP + `use_dpop_nonce` complete
over canonical `https://pds.test` and the authorize page renders with the **loopback
client accepted** (`__authorizeData.clientMetadata.redirect_uris = [http://127.0.0.1:<port>/...]`).

### 6.3 The consent step (real browser)

The **browser** still needs to reach `https://pds.test`, so for the interactive
part add a hosts entry `pds.test → $CADDY_IP` and trust `/tmp/caddy_ca.crt` in the
browser (or click through the warning). Then drive `…/oauth/authorize`, sign in, and
approve; the PDS redirects to your loopback `…/oauth/callback?code=…&state=…` and
the client exchanges the code. The authorize endpoint enforces browser **`Sec-Fetch-*`**
headers + a CSRF SPA, so this step is a genuine browser action, not easily scripted.

### 6.4 Identity-resolution caveat

A handle like `alice.pds.test` won't resolve without wildcard DNS (handle→DID is
HTTPS well-known on the *handle's* host). For local dev, point the client directly
at the known PDS as the authorization server (skip handle→DID→PDS resolution) and
pass the handle as `login_hint`, or add a DNS/hosts entry. `did:plc` resolution also
depends on the genesis op reaching the configured `PDS_DID_PLC_URL`.

## 7. Teardown / notes

`docker rm -f pds caddy` / `container rm -f pds caddy`; the data dir is recreatable
from §2. Everything here is throwaway local state — generate fresh secrets per boot;
don't reuse these creds anywhere real. Coordinate with decodingus if you want a
shared long-lived instance or a scripted `compose`/`Makefile` target.
