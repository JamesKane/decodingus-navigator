//! AT Proto OAuth for a **public/native client** (plan §7): PKCE only (no client
//! assertion), a loopback redirect, and tokens in the OS keychain. Reuses the shared
//! `du-atproto` primitives — the same handshake as the decodingus web client, swapping
//! the confidential `par_form`/`token_form` for the `_public` (PKCE-only) builders and
//! catching the redirect on a local `127.0.0.1` server instead of a web route.
//!
//! The interactive [`login`] needs a live PDS (an env-gated integration test drives it
//! against the test PDS); the local pieces — callback parsing, the loopback server,
//! redirect URIs — are unit-tested here.

use std::io::{BufRead, BufReader, Write};
use std::net::TcpListener;
use std::time::{SystemTime, UNIX_EPOCH};

use du_atproto::did::Did;
use du_atproto::oauth::{
    authorize_url, discover_auth_server, dpop_proof, par_form_public, random_token,
    token_form_public, EcKey, Pkce,
};
use du_atproto::Resolver;

use crate::error::SyncError;
use crate::tokens::Session;

/// Public/native client configuration. `client_id` is the (hosted) client-metadata.json
/// URL the authorization server fetches; `scope` is e.g. `"atproto navigatorCore"`.
#[derive(Debug, Clone)]
pub struct OAuthConfig {
    pub client_id: String,
    pub scope: String,
}

/// The loopback redirect URI for a chosen port.
pub fn redirect_uri(port: u16) -> String {
    format!("http://127.0.0.1:{port}/callback")
}

fn now() -> i64 {
    SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs() as i64).unwrap_or(0)
}

/// Parameters parsed from the authorization redirect.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct CallbackParams {
    pub code: Option<String>,
    pub state: String,
    pub error: Option<String>,
}

/// Parse the `code`/`state`/`error` from a redirect request target (`/callback?...`).
pub fn parse_callback_query(target: &str) -> CallbackParams {
    let mut params = CallbackParams::default();
    let query = target.split_once('?').map(|(_, q)| q).unwrap_or("");
    for pair in query.split('&') {
        if let Some((k, v)) = pair.split_once('=') {
            let v = percent_decode(v);
            match k {
                "code" => params.code = Some(v),
                "state" => params.state = v,
                "error" => params.error = Some(v),
                _ => {}
            }
        }
    }
    params
}

/// Minimal percent-decode (`+` -> space, `%XX` -> byte) for redirect query values.
fn percent_decode(s: &str) -> String {
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'+' => out.push(b' '),
            b'%' if i + 2 < bytes.len() => {
                let hi = (bytes[i + 1] as char).to_digit(16);
                let lo = (bytes[i + 2] as char).to_digit(16);
                if let (Some(h), Some(l)) = (hi, lo) {
                    out.push((h * 16 + l) as u8);
                    i += 3;
                    continue;
                }
                out.push(b'%');
            }
            b => out.push(b),
        }
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// A one-shot loopback HTTP server that catches the OAuth redirect.
pub struct LoopbackServer {
    listener: TcpListener,
}

impl LoopbackServer {
    /// Bind an ephemeral `127.0.0.1` port.
    pub fn bind() -> Result<Self, SyncError> {
        Ok(LoopbackServer { listener: TcpListener::bind("127.0.0.1:0")? })
    }

    pub fn port(&self) -> u16 {
        self.listener.local_addr().map(|a| a.port()).unwrap_or(0)
    }

    /// Block until the redirect arrives, parse its query, and reply with a closing page.
    pub fn wait(self) -> Result<CallbackParams, SyncError> {
        let (mut stream, _) = self.listener.accept()?;
        let mut request_line = String::new();
        BufReader::new(&stream).read_line(&mut request_line)?;
        // "GET /callback?code=...&state=... HTTP/1.1"
        let target = request_line.split_whitespace().nth(1).unwrap_or("");
        let params = parse_callback_query(target);

        let body = "<html><body><h3>DUNavigator</h3><p>Authentication complete — you can close this window.</p></body></html>";
        let response = format!(
            "HTTP/1.1 200 OK\r\nContent-Type: text/html\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
            body.len(),
            body
        );
        let _ = stream.write_all(response.as_bytes());
        Ok(params)
    }
}

/// POST a form with a DPoP proof, retrying once with the server-supplied nonce (the
/// AT Proto OAuth nonce dance). Returns the JSON body on success.
async fn post_with_dpop(
    http: &reqwest::Client,
    key: &EcKey,
    url: &str,
    form: &[(String, String)],
) -> Result<serde_json::Value, SyncError> {
    let proof = dpop_proof(key, "POST", url, now(), None, None);
    let resp = http.post(url).header("DPoP", proof).form(form).send().await?;
    if resp.status().is_success() {
        return Ok(resp.json().await?);
    }
    let nonce = resp.headers().get("DPoP-Nonce").and_then(|v| v.to_str().ok()).map(String::from);
    let Some(nonce) = nonce else {
        return Err(SyncError::Oauth(format!("{url}: {}", resp.status())));
    };
    let proof = dpop_proof(key, "POST", url, now(), Some(&nonce), None);
    let retry = http.post(url).header("DPoP", proof).form(form).send().await?;
    if retry.status().is_success() {
        Ok(retry.json().await?)
    } else {
        Err(SyncError::Oauth(format!("{url}: {}", retry.status())))
    }
}

/// Best-effort: open `url` in the user's browser (also printed so a headless user can
/// copy it).
fn open_browser(url: &str) {
    println!("Open this URL to authorize DUNavigator:\n  {url}");
    let cmd = if cfg!(target_os = "macos") {
        Some("open")
    } else if cfg!(target_os = "windows") {
        Some("explorer")
    } else if cfg!(target_os = "linux") {
        Some("xdg-open")
    } else {
        None
    };
    if let Some(cmd) = cmd {
        let _ = std::process::Command::new(cmd).arg(url).spawn();
    }
}

/// Run the full public-client OAuth login for `handle` (a handle or DID): resolve →
/// discover → PAR → browser authorize → loopback callback → token exchange. Returns the
/// authenticated [`Session`] (the caller persists it via `TokenStore`).
pub async fn login(
    http: &reqwest::Client,
    resolver: &Resolver,
    config: &OAuthConfig,
    handle: &str,
) -> Result<Session, SyncError> {
    let did = if handle.starts_with("did:") {
        Did::parse(handle)?
    } else {
        resolver.resolve_handle(handle).await?
    };
    let pds = resolver.resolve_pds(&did).await?;
    let meta = discover_auth_server(http, &pds).await?;
    let par_endpoint = meta
        .pushed_authorization_request_endpoint
        .clone()
        .ok_or_else(|| SyncError::Oauth("authorization server has no PAR endpoint".into()))?;

    let ec = EcKey::generate();
    let pkce = Pkce::generate();
    let state = random_token();
    let server = LoopbackServer::bind()?;
    let redirect = redirect_uri(server.port());

    // Pushed authorization request (PKCE, no client assertion).
    let par_form = par_form_public(&config.client_id, &redirect, &config.scope, &state, &pkce.challenge, Some(handle));
    let par = post_with_dpop(http, &ec, &par_endpoint, &par_form).await?;
    let request_uri = par
        .get("request_uri")
        .and_then(|v| v.as_str())
        .ok_or_else(|| SyncError::Oauth("PAR response missing request_uri".into()))?;

    open_browser(&authorize_url(&meta.authorization_endpoint, &config.client_id, request_uri));

    // Wait for the redirect on the loopback server (blocking accept off the runtime).
    let params = tokio::task::spawn_blocking(move || server.wait())
        .await
        .map_err(|e| SyncError::Oauth(format!("callback task failed: {e}")))??;
    if params.state != state {
        return Err(SyncError::Oauth("state mismatch (possible CSRF)".into()));
    }
    if let Some(err) = params.error {
        return Err(SyncError::Oauth(format!("authorization denied: {err}")));
    }
    let code = params.code.ok_or_else(|| SyncError::Oauth("callback missing code".into()))?;

    // Exchange the code for tokens (PKCE verifier, DPoP-bound).
    let token_form = token_form_public(&config.client_id, &redirect, &code, &pkce.verifier);
    let tok = post_with_dpop(http, &ec, &meta.token_endpoint, &token_form).await?;
    let access_token = tok
        .get("access_token")
        .and_then(|v| v.as_str())
        .ok_or_else(|| SyncError::Oauth("token response missing access_token".into()))?
        .to_string();
    let refresh_token = tok.get("refresh_token").and_then(|v| v.as_str()).unwrap_or("").to_string();
    let scope = tok.get("scope").and_then(|v| v.as_str()).unwrap_or(&config.scope).to_string();

    Ok(Session {
        did: did.as_str().to_string(),
        pds,
        access_token,
        refresh_token,
        dpop_key_b64: ec.to_base64(),
        scope,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_callback_code_and_state() {
        let p = parse_callback_query("/callback?code=abc123&state=xyz");
        assert_eq!(p.code.as_deref(), Some("abc123"));
        assert_eq!(p.state, "xyz");
        assert!(p.error.is_none());
    }

    #[test]
    fn parses_callback_error_and_decodes() {
        let p = parse_callback_query("/callback?error=access_denied&state=a%20b");
        assert_eq!(p.error.as_deref(), Some("access_denied"));
        assert_eq!(p.state, "a b");
        assert!(p.code.is_none());
    }

    #[test]
    fn redirect_uri_is_loopback() {
        assert_eq!(redirect_uri(49152), "http://127.0.0.1:49152/callback");
    }

    #[test]
    fn loopback_server_binds_an_ephemeral_port() {
        let s = LoopbackServer::bind().unwrap();
        assert!(s.port() >= 1024);
    }

    /// The loopback server captures a real redirect over a TCP connection.
    #[test]
    fn loopback_server_captures_a_redirect() {
        use std::io::Write;
        use std::net::TcpStream;

        let server = LoopbackServer::bind().unwrap();
        let port = server.port();
        let client = std::thread::spawn(move || {
            let mut stream = TcpStream::connect(("127.0.0.1", port)).unwrap();
            stream
                .write_all(b"GET /callback?code=THECODE&state=THESTATE HTTP/1.1\r\nHost: localhost\r\n\r\n")
                .unwrap();
        });
        let params = server.wait().unwrap();
        client.join().unwrap();
        assert_eq!(params.code.as_deref(), Some("THECODE"));
        assert_eq!(params.state, "THESTATE");
    }
}
