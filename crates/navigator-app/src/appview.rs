//! The one way that Navigator speaks to the `/api/v1/*` Edge API of the AppView.
//!
//! Three clients grew here separately: IBD exchange, social, and recruitment. Each client made the
//! same two request shapes.
//!
//! The first shape is a POST. The body of the POST holds the device-key signature, so the request
//! looks unauthenticated. The second shape is a signed GET with a replay guard. Its `did`, `ts`,
//! and `sig` values go on the query string.
//!
//! The IBD version and the social version were the same code. The single calls in `sync.rs` and
//! `matching.rs` wrote the same request a fourth time and a fifth time. All of this code is now in
//! this module. So the error map, the layout of the signature query, and the class of a non-2xx
//! response have one definition.
//!
//! These values cross the network: a DID, a timestamp, a signature, and the content that the caller
//! chose to send. A genotype never crosses. A coordinate never crosses.

use super::*;

/// A transport failure on a call to the AppView. Examples are a refused connection, a timeout, and
/// a TLS fault.
///
/// A plain `reqwest` client makes these calls. The calls do not go through the sync engine. But a
/// network failure has the same result on both paths. So this function returns the error variant
/// that the PDS paths use, and the offline indicator already knows that variant.
pub(crate) fn transport(e: reqwest::Error) -> AppError {
    AppError::Sync(navigator_sync::SyncError::from(e))
}

/// Change a non-2xx response from the AppView into an [`AppError::AppView`] for the user. The
/// function consumes `resp` to read the body. So the caller must keep the status first, if the
/// caller also needs it.
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

    /// Send a JSON body to an `/api/v1/<path>` endpoint with POST, and return the decoded
    /// response.
    ///
    /// The signature and the DID for that signature belong in `body`. These endpoints authenticate
    /// the device key for each call. They do not authenticate the HTTP request. For this reason the
    /// function takes a body that the caller signed. It does not sign the body, because the
    /// canonical string is different at each endpoint and only the caller knows it.
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

    /// Send a GET to an `/api/v1/<path>` endpoint with a device-key signature, and decode the
    /// response into `T`.
    ///
    /// `build_msg(did, ts)` makes the canonical string for the signature. That string is the only
    /// difference between a poll, a thread read, and an exchange pull. The `did`, `ts`, and `sig`
    /// values go on the query, together with `extra`. The timestamp gives the signature its replay
    /// guard.
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
