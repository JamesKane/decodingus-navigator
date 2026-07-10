//! `impl App` methods for the AppView's signed recruitment Edge API (`/api/v1/recruitment/*`,
//! social roadmap 3c) — the **response** side: list the caller's open invitations and accept/decline
//! them. Campaign creation stays on the AppView web flow (it's gated to a group-project admin, which
//! the Navigator can't yet act as). Device-key-signed like the social/exchange clients; reuses the
//! shared [`social_post`](App::social_post) / [`social_get`](App::social_get) transport. Invitations
//! also arrive as SYSTEM notifications, so this pairs with the Community → Notifications surface.
use super::*;

use navigator_sync::recruitment::messages;

/// One open recruitment invitation, as the AppView's `/recruitment/invitations` returns it.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct RecruitmentInvitation {
    pub campaign_id: i64,
    pub title: String,
    pub message: String,
    pub project_name: String,
}

impl App {
    /// The signed-in account's open (INVITED) recruitment invitations.
    pub async fn recruitment_invitations(&self) -> Result<Vec<RecruitmentInvitation>, AppError> {
        #[derive(serde::Deserialize)]
        struct Resp {
            #[serde(default)]
            items: Vec<RecruitmentInvitation>,
        }
        let r: Resp = self.social_get("recruitment/invitations", messages::poll, &[]).await?;
        Ok(r.items)
    }

    /// Accept (`true`) or decline (`false`) a recruitment invitation. Returns whether it changed
    /// (a no-op if already responded). On acceptance the AppView notifies the researcher.
    pub async fn recruitment_respond(&self, campaign_id: i64, accept: bool) -> Result<bool, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = chrono::Utc::now().timestamp();
        let sig = dev.sign_fresh(ts, &messages::respond(&did, campaign_id, accept));
        let body = serde_json::json!({
            "did": did,
            "campaign_id": campaign_id,
            "accept": accept,
            "ts": ts,
            "signature": sig,
        });
        let v = self.social_post("recruitment/respond", body).await?;
        Ok(v.get("changed").and_then(|x| x.as_bool()).unwrap_or(false))
    }
}
