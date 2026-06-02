//! Authenticated PDS writes (`com.atproto.repo`) for a logged-in account.
//!
//! Resource requests use the DPoP-bound access token: `Authorization: DPoP <token>` plus
//! a DPoP proof whose `ath` is the access-token hash, with the same `use_dpop_nonce`
//! retry as the OAuth endpoints. A `Bearer` mode (a plain `createAccount` session token,
//! no DPoP) exists so repo CRUD can be exercised live without a browser-minted token.

use std::time::{SystemTime, UNIX_EPOCH};

use du_atproto::oauth::{dpop_proof, EcKey};

use crate::error::SyncError;
use crate::tokens::Session;

fn now() -> i64 {
    SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs() as i64).unwrap_or(0)
}

/// Build a reqwest client, optionally trusting a dev CA (`NAVIGATOR_DEV_CA`, a PEM path)
/// and pinning a host→IP (`NAVIGATOR_DEV_RESOLVE`, `host:ip`) so a TLS-proxied local PDS
/// at its canonical `https://` name is reachable without `/etc/hosts` (PDS doc §6.2).
pub fn dev_http_client() -> reqwest::Client {
    let mut builder = reqwest::Client::builder();
    if let Some(ca) = std::env::var("NAVIGATOR_DEV_CA").ok().filter(|s| !s.is_empty()) {
        if let Ok(pem) = std::fs::read(&ca) {
            if let Ok(cert) = reqwest::Certificate::from_pem(&pem) {
                builder = builder.add_root_certificate(cert);
            }
        }
    }
    if let Some(spec) = std::env::var("NAVIGATOR_DEV_RESOLVE").ok().filter(|s| !s.is_empty()) {
        if let Some((host, ip)) = spec.rsplit_once(':') {
            if let Ok(addr) = format!("{ip}:443").parse::<std::net::SocketAddr>() {
                builder = builder.resolve(host, addr);
            }
        }
    }
    builder.build().unwrap_or_default()
}

/// A reference to a written record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordRef {
    pub uri: String,
    pub cid: String,
}

impl RecordRef {
    /// The record key (last `at://…/<collection>/<rkey>` segment).
    pub fn rkey(&self) -> &str {
        self.uri.rsplit('/').next().unwrap_or("")
    }
}

/// How a PDS request authenticates.
enum Auth {
    /// OAuth DPoP-bound access token (the production path).
    Dpop { token: String, key: EcKey },
    /// A plain Bearer session token (e.g. `createAccount`'s `accessJwt`).
    Bearer(String),
}

/// Authenticated client for one repo (account) on a PDS.
pub struct PdsClient {
    http: reqwest::Client,
    pds_base: String,
    did: String,
    auth: Auth,
}

impl PdsClient {
    /// Build from a logged-in OAuth [`Session`] (DPoP-bound writes).
    pub fn from_session(http: reqwest::Client, session: &Session) -> Result<Self, SyncError> {
        let key = EcKey::from_base64(&session.dpop_key_b64)?;
        Ok(PdsClient {
            http,
            pds_base: session.pds.trim_end_matches('/').to_string(),
            did: session.did.clone(),
            auth: Auth::Dpop { token: session.access_token.clone(), key },
        })
    }

    /// Bearer-auth client (no DPoP) — for repo CRUD against a `createAccount` session.
    pub fn bearer(http: reqwest::Client, pds_base: &str, did: &str, token: &str) -> Self {
        PdsClient {
            http,
            pds_base: pds_base.trim_end_matches('/').to_string(),
            did: did.to_string(),
            auth: Auth::Bearer(token.to_string()),
        }
    }

    /// `com.atproto.repo.createRecord` — put `record` into `collection` (optional `rkey`).
    pub async fn create_record(
        &self,
        collection: &str,
        record: serde_json::Value,
        rkey: Option<&str>,
    ) -> Result<RecordRef, SyncError> {
        let body = create_record_body(&self.did, collection, record, rkey);
        let v = self.post("com.atproto.repo.createRecord", &body).await?;
        let uri = v
            .get("uri")
            .and_then(|x| x.as_str())
            .ok_or_else(|| SyncError::Oauth("createRecord response missing uri".into()))?
            .to_string();
        let cid = v.get("cid").and_then(|x| x.as_str()).unwrap_or_default().to_string();
        Ok(RecordRef { uri, cid })
    }

    /// `com.atproto.repo.getRecord` — read a record's value (public read).
    pub async fn get_record(&self, collection: &str, rkey: &str) -> Result<serde_json::Value, SyncError> {
        let url = format!("{}/xrpc/com.atproto.repo.getRecord", self.pds_base);
        let resp = self
            .http
            .get(&url)
            .query(&[("repo", self.did.as_str()), ("collection", collection), ("rkey", rkey)])
            .send()
            .await?;
        if !resp.status().is_success() {
            return Err(SyncError::Oauth(format!("getRecord: {}", resp.status())));
        }
        Ok(resp.json().await?)
    }

    /// POST an XRPC procedure with the configured auth (DPoP nonce-retry for DPoP).
    async fn post(&self, nsid: &str, body: &serde_json::Value) -> Result<serde_json::Value, SyncError> {
        let url = format!("{}/xrpc/{nsid}", self.pds_base);
        match &self.auth {
            Auth::Bearer(token) => {
                let resp = self.http.post(&url).bearer_auth(token).json(body).send().await?;
                if resp.status().is_success() {
                    Ok(resp.json().await?)
                } else {
                    Err(xrpc_error(nsid, resp).await)
                }
            }
            Auth::Dpop { token, key } => {
                let proof = dpop_proof(key, "POST", &url, now(), None, Some(token));
                let resp = self
                    .http
                    .post(&url)
                    .header("Authorization", format!("DPoP {token}"))
                    .header("DPoP", proof)
                    .json(body)
                    .send()
                    .await?;
                if resp.status().is_success() {
                    return Ok(resp.json().await?);
                }
                let nonce = resp.headers().get("DPoP-Nonce").and_then(|v| v.to_str().ok()).map(String::from);
                let Some(nonce) = nonce else {
                    return Err(SyncError::Oauth(format!("{nsid}: {}", resp.status())));
                };
                let proof = dpop_proof(key, "POST", &url, now(), Some(&nonce), Some(token));
                let retry = self
                    .http
                    .post(&url)
                    .header("Authorization", format!("DPoP {token}"))
                    .header("DPoP", proof)
                    .json(body)
                    .send()
                    .await?;
                if retry.status().is_success() {
                    Ok(retry.json().await?)
                } else {
                    Err(xrpc_error(nsid, retry).await)
                }
            }
        }
    }
}

/// Turn a non-2xx XRPC response into an error that includes the PDS's reason body.
async fn xrpc_error(nsid: &str, resp: reqwest::Response) -> SyncError {
    let status = resp.status();
    let body = resp.text().await.unwrap_or_default();
    SyncError::Oauth(format!("{nsid}: {status}: {body}"))
}

/// The `createRecord` request body.
fn create_record_body(did: &str, collection: &str, record: serde_json::Value, rkey: Option<&str>) -> serde_json::Value {
    let mut body = serde_json::json!({ "repo": did, "collection": collection, "record": record });
    if let Some(rkey) = rkey {
        body["rkey"] = serde_json::Value::String(rkey.to_string());
    }
    body
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn create_record_body_shape() {
        let body = create_record_body(
            "did:plc:abc",
            "com.decodingus.test.rec",
            json!({"$type": "com.decodingus.test.rec", "n": 1}),
            None,
        );
        assert_eq!(body["repo"], "did:plc:abc");
        assert_eq!(body["collection"], "com.decodingus.test.rec");
        assert_eq!(body["record"]["n"], 1);
        assert!(body.get("rkey").is_none());

        let with_rkey = create_record_body("did:plc:abc", "c", json!({}), Some("self"));
        assert_eq!(with_rkey["rkey"], "self");
    }

    #[test]
    fn record_ref_rkey() {
        let r = RecordRef { uri: "at://did:plc:x/com.decodingus.test.rec/3kabc".into(), cid: "bafy".into() };
        assert_eq!(r.rkey(), "3kabc");
    }

    /// Live repo CRUD against a local PDS: create a throwaway account, write a record via
    /// Bearer auth, read it back. Proves the XRPC write/read path against a real PDS
    /// (the DPoP-bound variant uses the same code with a DPoP proof; the OAuth DPoP path
    /// is proven by the PAR test + dpop_proof's `ath`). Set PDS_TEST_URL to run.
    #[tokio::test]
    #[ignore = "requires PDS_TEST_URL (local atproto PDS container)"]
    async fn create_and_read_record_against_live_pds() {
        let Ok(pds) = std::env::var("PDS_TEST_URL") else {
            eprintln!("PDS_TEST_URL unset — skipping live record CRUD test");
            return;
        };
        let pds = pds.trim_end_matches('/').to_string();
        let http = reqwest::Client::new();

        // Throwaway account (unique handle under the PDS host).
        // Short unique suffix (low digits of the nanos clock) — handle labels are length-limited.
        let n = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos() % 1_000_000_000;
        let handle = format!("nav{n}.pds.test");
        let acct: serde_json::Value = http
            .post(format!("{pds}/xrpc/com.atproto.server.createAccount"))
            .json(&json!({ "handle": handle, "email": format!("nav{n}@example.test"), "password": "nav-pw-123456" }))
            .send()
            .await
            .expect("createAccount request")
            .json()
            .await
            .expect("createAccount json");
        let did = acct["did"].as_str().expect("did");
        let jwt = acct["accessJwt"].as_str().expect("accessJwt");

        let client = PdsClient::bearer(http, &pds, did, jwt);
        // Publish the real typed coverage-summary record (floats as strings — atproto
        // DAG-CBOR rejects floats).
        let rec = crate::records::CoverageSummaryRecord::new(
            "chm13v2.0", 178.81, 182.0, 28.9, 1.0, 1.0, 1.0, 16569, 16292, "2026-06-02T00:00:00Z",
        );
        let record = serde_json::to_value(&rec).unwrap();
        let r = client
            .create_record(crate::records::COVERAGE_SUMMARY_COLLECTION, record, None)
            .await
            .expect("createRecord");
        assert!(r.uri.starts_with("at://"), "uri: {}", r.uri);

        let got = client.get_record(crate::records::COVERAGE_SUMMARY_COLLECTION, r.rkey()).await.expect("getRecord");
        assert_eq!(got["value"]["meanCoverage"], "178.81");
        assert_eq!(got["value"]["callableBases"], 16292);
        assert_eq!(got["value"]["referenceBuild"], "chm13v2.0");
        eprintln!("✓ wrote + read {}", r.uri);
    }
}
