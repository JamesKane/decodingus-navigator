//! The one way Navigator talks to the AppView's `/api/v1/*` Edge API.
//!
//! Three clients grew here independently — IBD exchange, social, recruitment — and each arrived at
//! the same two shapes: an unauthenticated-looking POST whose body carries the device-key
//! signature, and a replay-guarded signed GET whose `did`/`ts`/`sig` ride on the query string. The
//! IBD and social versions were byte-for-byte identical, and the remaining one-off calls in
//! `sync.rs` / `matching.rs` open-coded the same thing a fourth and fifth time. They are all this
//! module now, so the error mapping, the signing-query layout, and the non-2xx classification are
//! decided once.
//!
//! What travels: a DID, a timestamp, a signature, and whatever the caller chose to send. Never
//! genotypes, never coordinates.

use super::*;

/// A transport failure (connection refused, timeout, TLS) on an AppView call.
///
/// The AppView is reached with a bare `reqwest` client rather than through the sync engine, but a
/// network failure means the same thing either way, so it lands in the same error variant the PDS
/// paths use and the offline indicator already understands.
pub(crate) fn transport(e: reqwest::Error) -> AppError {
    AppError::Sync(navigator_sync::SyncError::from(e))
}

/// Classify a non-2xx AppView response into a user-facing [`AppError::AppView`]. Consumes `resp` to
/// read the body (so capture the status first at the call site if it is also needed).
pub(crate) async fn status_error(api: &str, resp: reqwest::Response) -> AppError {
    let status = resp.status();
    let body = resp.text().await.unwrap_or_default();
    match status.as_u16() {
        403 => AppError::AppView(format!(
            "{api}: device key not yet registered or verified by the AppView (403)"
        )),
        422 => AppError::AppView(format!(
            "{api}: request rejected, likely clock skew (422) — check the system clock"
        )),
        _ => AppError::AppView(format!("{api}: {status}: {body}")),
    }
}

impl App {
    /// The absolute URL of an `/api/v1/<path>` endpoint on the configured AppView.
    pub(crate) fn appview_url(&self, path: &str) -> String {
        format!("{}/api/v1/{path}", decodingus_appview_url())
    }

    /// POST a JSON body to an `/api/v1/<path>` endpoint and return the decoded response.
    ///
    /// The signature (and the DID it is over) belongs in `body` — these endpoints authenticate the
    /// device key per call, not the HTTP request — so this deliberately takes an already-signed
    /// body rather than signing on the caller's behalf: the canonical string differs per endpoint
    /// and only the caller knows it.
    pub(crate) async fn appview_post(
        &self,
        path: &str,
        body: serde_json::Value,
    ) -> Result<serde_json::Value, AppError> {
        let resp = self
            .auth
            .http
            .post(self.appview_url(path))
            .json(&body)
            .send()
            .await
            .map_err(transport)?;
        if !resp.status().is_success() {
            return Err(status_error(path, resp).await);
        }
        resp.json().await.map_err(transport)
    }

    /// Device-key-signed GET to an `/api/v1/<path>` endpoint, decoded into `T`.
    ///
    /// `build_msg(did, ts)` produces the canonical string to sign — the one thing that varies
    /// between a poll, a thread read, and an exchange pull. `did`/`ts`/`sig` plus `extra` go on the
    /// query; the timestamp is what makes the signature replay-guarded.
    pub(crate) async fn appview_get_signed<T, F>(
        &self,
        path: &str,
        build_msg: F,
        extra: &[(&str, &str)],
    ) -> Result<T, AppError>
    where
        T: serde::de::DeserializeOwned,
        F: Fn(&str, i64) -> String,
    {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = Utc::now().timestamp();
        let sig = dev.sign(&build_msg(&did, ts));
        let ts_s = ts.to_string();
        let mut query: Vec<(&str, &str)> = vec![("did", did.as_str()), ("ts", ts_s.as_str()), ("sig", sig.as_str())];
        query.extend_from_slice(extra);
        let resp = self
            .auth
            .http
            .get(self.appview_url(path))
            .query(&query)
            .send()
            .await
            .map_err(transport)?;
        if !resp.status().is_success() {
            return Err(status_error(path, resp).await);
        }
        resp.json().await.map_err(transport)
    }
}
