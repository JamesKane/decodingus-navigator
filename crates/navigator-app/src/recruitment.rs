//! `impl App` methods for the signed recruitment Edge API of the AppView
//! (`/api/v1/recruitment/*`, social roadmap 3c).
//!
//! This module is the **response** side. It lists the open invitations of the caller. It also
//! accepts or declines them. To make a campaign, the user must use the web flow of the AppView.
//! Only an administrator of a group project can make a campaign, and Navigator can not yet act as
//! one.
//!
//! The device key signs each call, as it does for the social client and the exchange client. This
//! module uses the shared [`appview_post`](App::appview_post) and
//! [`appview_get_signed`](App::appview_get_signed) transport. An invitation also arrives as a
//! SYSTEM notification. So this module works together with the Community → Notifications view.
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
        let r: Resp = self
            .appview_get_signed("recruitment/invitations", messages::poll, &[])
            .await?;
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
        let v = self.appview_post("recruitment/respond", body).await?;
        Ok(v.get("changed").and_then(|x| x.as_bool()).unwrap_or(false))
    }
}
