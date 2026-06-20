//! `impl App` methods for the AppView's signed social Edge API (`/api/v1/social/*`) — the
//! communication core the alpha/beta testers use to reach the team (support threads), read the
//! community feed (+ federated posts), and receive notifications.
//!
//! Every call is **device-key-signed** (no per-call OAuth), exactly like the IBD `exchange` client:
//! reads are a replay-guarded signed GET (`did`/`ts`/`sig` query); writes carry `did` + `signature`
//! (and `ts` where the canonical string includes it) in the JSON body. Canonical signing strings live
//! in [`navigator_sync::social::messages`] and mirror the AppView byte-for-byte. PII-free: only a DID,
//! a signature, and content the user chose to send crosses the wire.

use super::*;

use navigator_sync::social::messages;

/// One of the caller's support threads (team↔tester), as listed by `GET /social/threads`.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct SocialThreadSummary {
    /// Conversation id (UUID string) — the key for reading/replying.
    pub conversation_id: String,
    #[serde(default)]
    pub subject: Option<String>,
    /// `open` | `replied` | `closed`.
    pub status: String,
    #[serde(default)]
    pub last_message_at: Option<String>,
    /// The team has posted since the user last read.
    #[serde(default)]
    pub unread: bool,
}

/// One message within a support thread (`GET /social/thread/:id`).
#[derive(Debug, Clone, serde::Deserialize)]
pub struct SocialMessage {
    /// Posted by a Curator/Admin (vs. the tester).
    #[serde(default)]
    pub from_team: bool,
    #[serde(default)]
    pub author: Option<String>,
    pub body: String,
    #[serde(default)]
    pub at: Option<String>,
}

/// The community feed (`GET /social/feed`): team announcements, AppView-native community posts, and
/// read-only PDS-federated posts.
#[derive(Debug, Clone, Default, serde::Deserialize)]
pub struct FeedView {
    #[serde(default)]
    pub announcements: Vec<FeedItem>,
    #[serde(default)]
    pub community: Vec<FeedItem>,
    #[serde(default)]
    pub federated: Vec<FederatedItem>,
}

/// An announcement or community post.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct FeedItem {
    pub id: String,
    /// `ANNOUNCEMENT` | `COMMUNITY`.
    #[serde(default)]
    pub kind: Option<String>,
    #[serde(default)]
    pub author: Option<String>,
    #[serde(default)]
    pub topic: Option<String>,
    pub content: String,
    #[serde(default)]
    pub pinned: bool,
    #[serde(default)]
    pub at: Option<String>,
    #[serde(default)]
    pub reply_count: i64,
    #[serde(default)]
    pub score: i64,
    #[serde(default)]
    pub parent_post_id: Option<String>,
}

/// A PDS-federated community post mirrored into the feed (read-only — voting/reply/block stay native).
#[derive(Debug, Clone, serde::Deserialize)]
pub struct FederatedItem {
    #[serde(default)]
    pub did: Option<String>,
    #[serde(default)]
    pub author: Option<String>,
    pub text: String,
    #[serde(default)]
    pub topic: Option<String>,
    #[serde(default)]
    pub uri: Option<String>,
    #[serde(default)]
    pub at: Option<String>,
}

/// One notification (`GET /social/notifications`).
#[derive(Debug, Clone, serde::Deserialize)]
pub struct SocialNotification {
    pub id: String,
    /// `THREAD_REPLY` | `FEED_REPLY` | `MATCH` | `SYSTEM` | …
    pub kind: String,
    pub title: String,
    #[serde(default)]
    pub body: Option<String>,
    #[serde(default)]
    pub link: Option<String>,
    #[serde(default)]
    pub actor: Option<String>,
    #[serde(default)]
    pub at: Option<String>,
    #[serde(default)]
    pub unread: bool,
}

/// The notifications response: the list + the server's unread count (for the app-bar bell).
#[derive(Debug, Clone, Default, serde::Deserialize)]
pub struct NotificationList {
    #[serde(default)]
    pub items: Vec<SocialNotification>,
    #[serde(default)]
    pub unread: i64,
}

impl App {
    // ---- support threads ---------------------------------------------------

    /// List the signed-in account's support threads (newest first).
    pub async fn support_threads(&self) -> Result<Vec<SocialThreadSummary>, AppError> {
        #[derive(serde::Deserialize)]
        struct Resp {
            #[serde(default)]
            items: Vec<SocialThreadSummary>,
        }
        let r: Resp = self
            .social_get("social/threads", messages::poll, &[])
            .await?;
        Ok(r.items)
    }

    /// Read one thread's messages (oldest first). Marks the thread read on the user side (server-side).
    pub async fn support_thread(&self, conversation_id: &str) -> Result<Vec<SocialMessage>, AppError> {
        #[derive(serde::Deserialize)]
        struct Resp {
            #[serde(default)]
            items: Vec<SocialMessage>,
        }
        let path = format!("social/thread/{conversation_id}");
        let r: Resp = self
            .social_get(&path, |d, ts| messages::thread_read(d, conversation_id, ts), &[])
            .await?;
        Ok(r.items)
    }

    /// Open a new support thread; returns the new conversation id.
    pub async fn open_support_thread(&self, subject: &str, body: &str) -> Result<String, AppError> {
        self.write_thread(None, Some(subject), body).await
    }

    /// Reply to an existing thread the caller owns; returns its conversation id.
    pub async fn reply_support_thread(&self, conversation_id: &str, body: &str) -> Result<String, AppError> {
        self.write_thread(Some(conversation_id), None, body).await
    }

    async fn write_thread(
        &self,
        conversation_id: Option<&str>,
        subject: Option<&str>,
        body: &str,
    ) -> Result<String, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let sig = dev.sign(&messages::thread(&did, conversation_id));
        let mut b = serde_json::json!({ "did": did, "body": body, "signature": sig });
        if let Some(c) = conversation_id {
            b["conversation_id"] = serde_json::json!(c);
        }
        if let Some(s) = subject {
            b["subject"] = serde_json::json!(s);
        }
        let v = self.social_post("social/thread", b).await?;
        Ok(v.get("conversation_id").and_then(|x| x.as_str()).unwrap_or_default().to_string())
    }

    // ---- community feed ----------------------------------------------------

    /// Read the community feed: announcements + community posts + federated mirror.
    pub async fn community_feed(&self) -> Result<FeedView, AppError> {
        self.social_get("social/feed", messages::poll, &[]).await
    }

    /// Post to the community feed (optionally tagged with a `topic`, or as a reply to `parent`);
    /// returns the new post id. A reputation gate maps to [`AppError`] (HTTP 403) — surface it as a
    /// "not enough reputation yet" hint in the UI.
    pub async fn post_community(
        &self,
        content: &str,
        topic: Option<&str>,
        parent: Option<&str>,
    ) -> Result<String, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let sig = dev.sign(&messages::post(&did, parent));
        let mut b = serde_json::json!({ "did": did, "content": content, "signature": sig });
        if let Some(t) = topic {
            b["topic"] = serde_json::json!(t);
        }
        if let Some(p) = parent {
            b["parent_post_id"] = serde_json::json!(p);
        }
        let v = self.social_post("social/post", b).await?;
        Ok(v.get("id").and_then(|x| x.as_str()).unwrap_or_default().to_string())
    }

    // ---- notifications -----------------------------------------------------

    /// The signed-in account's notifications + unread count.
    pub async fn notifications(&self) -> Result<NotificationList, AppError> {
        self.social_get("social/notifications", messages::poll, &[]).await
    }

    /// Mark one notification read (`id = Some`) or all (`id = None`); returns how many were marked.
    pub async fn mark_notification_read(&self, id: Option<&str>) -> Result<i64, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = Utc::now().timestamp();
        let sig = dev.sign(&messages::notif_read(&did, id, ts));
        let mut b = serde_json::json!({ "did": did, "ts": ts, "signature": sig });
        if let Some(i) = id {
            b["id"] = serde_json::json!(i);
        }
        let v = self.social_post("social/notifications/read", b).await?;
        Ok(v.get("marked").and_then(|x| x.as_i64()).unwrap_or(0))
    }

    // ---- transport helpers (mirror the IBD exchange client) ----------------

    /// POST a JSON body to a `/api/v1/social/<…>` endpoint, mapping non-2xx to an AppView error.
    async fn social_post(&self, path: &str, body: serde_json::Value) -> Result<serde_json::Value, AppError> {
        let url = format!("{}/api/v1/{path}", decodingus_appview_url());
        let resp = self
            .auth
            .http
            .post(&url)
            .json(&body)
            .send()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        if !resp.status().is_success() {
            return Err(appview_status_error(path, resp).await);
        }
        resp.json()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))
    }

    /// Device-key-signed GET to a `/api/v1/social/<…>` endpoint. `build_msg(did, ts)` produces the
    /// canonical string to sign (poll or thread-read); `did`/`ts`/`sig` + `extra` go on the query.
    async fn social_get<T, F>(&self, path: &str, build_msg: F, extra: &[(&str, &str)]) -> Result<T, AppError>
    where
        T: serde::de::DeserializeOwned,
        F: Fn(&str, i64) -> String,
    {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let url = format!("{}/api/v1/{path}", decodingus_appview_url());
        let ts = Utc::now().timestamp();
        let sig = dev.sign(&build_msg(&did, ts));
        let ts_s = ts.to_string();
        let mut query: Vec<(&str, &str)> = vec![("did", did.as_str()), ("ts", ts_s.as_str()), ("sig", sig.as_str())];
        query.extend_from_slice(extra);
        let resp = self
            .auth
            .http
            .get(&url)
            .query(&query)
            .send()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))?;
        if !resp.status().is_success() {
            return Err(appview_status_error(path, resp).await);
        }
        resp.json()
            .await
            .map_err(|e| AppError::Sync(navigator_sync::SyncError::from(e)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The DTOs deserialize from the AppView's exact wire shapes (`social_edge.rs` json! bodies).
    #[test]
    fn dtos_match_appview_wire_shapes() {
        let threads: Vec<SocialThreadSummary> = serde_json::from_value(serde_json::json!([
            { "conversation_id": "c1", "subject": "hi", "status": "replied", "last_message_at": "2026-06-20T00:00:00Z", "unread": true },
            { "conversation_id": "c2", "subject": null, "status": "open", "last_message_at": null, "unread": false }
        ]))
        .unwrap();
        assert_eq!(threads.len(), 2);
        assert!(threads[0].unread);
        assert_eq!(threads[1].subject, None);

        let msgs: Vec<SocialMessage> = serde_json::from_value(serde_json::json!([
            { "from_team": false, "author": "Tester", "body": "first", "at": "2026-06-20T00:00:00Z" },
            { "from_team": true, "author": null, "body": "reply", "at": "2026-06-20T00:01:00Z" }
        ]))
        .unwrap();
        assert!(msgs[1].from_team);

        let feed: FeedView = serde_json::from_value(serde_json::json!({
            "announcements": [
                { "id": "a1", "kind": "ANNOUNCEMENT", "author": "Team", "topic": null, "content": "v2 out",
                  "pinned": true, "at": "2026-06-20T00:00:00Z", "reply_count": 5, "score": 12, "parent_post_id": null }
            ],
            "community": [
                { "id": "p1", "kind": "COMMUNITY", "author": "User", "topic": "haplogroup:R-M269", "content": "neat",
                  "pinned": false, "at": "2026-06-20T00:00:00Z", "reply_count": 0, "score": 1, "parent_post_id": null }
            ],
            "federated": [
                { "did": "did:plc:x", "author": "Remote", "text": "hello", "topic": "general",
                  "uri": "at://did:plc:x/com.decodingus.atmosphere.feed.post/3k", "at": "2026-06-20T00:00:00Z" }
            ]
        }))
        .unwrap();
        assert_eq!(feed.announcements.len(), 1);
        assert!(feed.announcements[0].pinned);
        assert_eq!(feed.community[0].topic.as_deref(), Some("haplogroup:R-M269"));
        assert_eq!(feed.federated[0].text, "hello");

        let notifs: NotificationList = serde_json::from_value(serde_json::json!({
            "items": [
                { "id": "n1", "kind": "THREAD_REPLY", "title": "The team replied", "body": null,
                  "link": "/messages/c1", "actor": "Team", "at": "2026-06-20T00:00:00Z", "unread": true }
            ],
            "unread": 1
        }))
        .unwrap();
        assert_eq!(notifs.unread, 1);
        assert_eq!(notifs.items[0].kind, "THREAD_REPLY");

        // An empty feed (all sections defaulted) is valid.
        let empty: FeedView = serde_json::from_value(serde_json::json!({})).unwrap();
        assert!(empty.community.is_empty());
    }
}
