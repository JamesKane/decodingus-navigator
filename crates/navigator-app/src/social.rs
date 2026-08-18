//! `impl App` methods for the signed social Edge API of the AppView (`/api/v1/social/*`).
//!
//! This module is the communication core. A tester uses it for three tasks. The tester speaks to
//! the team in a support thread, reads the community feed with its federated posts, and receives a
//! notification.
//!
//! The device key signs **each call**, and no call uses OAuth. The IBD `exchange` client works in
//! the same way. A read is a signed GET with a replay guard, and its `did`, `ts`, and `sig` values
//! go on the query. A write puts `did` and `signature` in the JSON body. The write also puts `ts`
//! there when the canonical string holds a timestamp.
//!
//! [`navigator_sync::social::messages`] holds the canonical strings for a signature. These strings
//! are the same as the strings of the AppView. The module sends no personal data. Only a DID, a
//! signature, and the content that the user chose to send cross the network.

use super::*;

use navigator_sync::social::messages;

/// One of the caller's support threads (team↔tester), as listed by `GET /social/threads`.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct SocialThreadSummary {
    /// The conversation id, as a UUID string. It is the key to read the thread and to reply.
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

/// A federated community post from a PDS, copied into the feed. The user can only read it. A vote,
/// a reply, and a block stay in the native AppView records.
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
        let r: Resp = self.appview_get_signed("social/threads", messages::poll, &[]).await?;
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
            .appview_get_signed(&path, |d, ts| messages::thread_read(d, conversation_id, ts), &[])
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
        let ts = Utc::now().timestamp();
        let sig = dev.sign_fresh(ts, &messages::thread(&did, conversation_id));
        let mut b = serde_json::json!({ "did": did, "body": body, "ts": ts, "signature": sig });
        if let Some(c) = conversation_id {
            b["conversation_id"] = serde_json::json!(c);
        }
        if let Some(s) = subject {
            b["subject"] = serde_json::json!(s);
        }
        let v = self.appview_post("social/thread", b).await?;
        Ok(v.get("conversation_id")
            .and_then(|x| x.as_str())
            .unwrap_or_default()
            .to_string())
    }

    // ---- community feed ----------------------------------------------------

    /// Read the community feed: announcements + community posts + federated mirror.
    pub async fn community_feed(&self) -> Result<FeedView, AppError> {
        self.appview_get_signed("social/feed", messages::poll, &[]).await
    }

    /// Send a post to the community feed and return the id of the new post. The caller can add a
    /// `topic` tag, or make the post a reply to `parent`.
    ///
    /// A reputation gate gives HTTP 403, which becomes an [`AppError`]. Show it in the UI as a hint
    /// that the user does not have enough reputation.
    pub async fn post_community(
        &self,
        content: &str,
        topic: Option<&str>,
        parent: Option<&str>,
    ) -> Result<String, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = Utc::now().timestamp();
        let sig = dev.sign_fresh(ts, &messages::post(&did, parent));
        let mut b = serde_json::json!({ "did": did, "content": content, "ts": ts, "signature": sig });
        if let Some(t) = topic {
            b["topic"] = serde_json::json!(t);
        }
        if let Some(p) = parent {
            b["parent_post_id"] = serde_json::json!(p);
        }
        let v = self.appview_post("social/post", b).await?;
        Ok(v.get("id").and_then(|x| x.as_str()).unwrap_or_default().to_string())
    }

    /// Publish a community post to the PDS of the active account. The record is a federated
    /// `com.decodingus.atmosphere.feed.post` record (roadmap 3b).
    ///
    /// The Jetstream consumer of the AppView copies the record into the community feed as a
    /// read-only "Atmosphere" entry. So this method is the portable, federated form of
    /// [`post_community`](Self::post_community), which writes a native AppView record.
    ///
    /// The record goes through the sync **outbox**, so the publish is durable. It continues after a
    /// restart, and it tries again with a longer delay after a temporary failure or an offline
    /// failure.
    ///
    /// Each post is a **separate** record. The outbox `entity_ref` is a new id, because the app only
    /// appends a post and never joins two posts. A summary record for one entity behaves in a
    /// different way. `rkey: None` lets the PDS choose the TID.
    ///
    /// A federated post is **not** in `PUBLISHED_COLLECTIONS`, by design. A PULL reconcile must never
    /// return a post that the user deleted on their PDS.
    ///
    /// The method fails when no account is active. It also fails for a local `did:key` identity,
    /// because that identity certifies itself and has no PDS repository to write to. The federated
    /// feed needs an OAuth account with a PDS. The UI gates the option and shows the error as a
    /// hint.
    pub async fn publish_feed_post(&self, content: &str, topic: Option<&str>) -> Result<(), AppError> {
        let did = self.require_account()?;
        if did.starts_with("did:key:") {
            return Err(AppError::Import(
                "publishing to the federated feed needs a signed-in PDS account — the local identity has no repo"
                    .into(),
            ));
        }
        let value = feed_post_record(content, topic, None);
        let entity_ref = format!("feed_post:{}", Uuid::new_v4());
        self.enqueue_publish("feed_post", &entity_ref, NS_FEED_POST, None, value)
            .await
    }

    // ---- notifications -----------------------------------------------------

    /// The signed-in account's notifications + unread count.
    pub async fn notifications(&self) -> Result<NotificationList, AppError> {
        self.appview_get_signed("social/notifications", messages::poll, &[])
            .await
    }

    /// Mark one notification as read with `id = Some`, or mark all of them with `id = None`. The
    /// method returns the count of the notifications that it marked.
    pub async fn mark_notification_read(&self, id: Option<&str>) -> Result<i64, AppError> {
        let did = self.current_account().ok_or(AppError::NotAuthenticated)?;
        let dev = self.ensure_device_key().await?;
        let ts = Utc::now().timestamp();
        let sig = dev.sign(&messages::notif_read(&did, id, ts));
        let mut b = serde_json::json!({ "did": did, "ts": ts, "signature": sig });
        if let Some(i) = id {
            b["id"] = serde_json::json!(i);
        }
        let v = self.appview_post("social/notifications/read", b).await?;
        Ok(v.get("marked").and_then(|x| x.as_i64()).unwrap_or(0))
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

    /// The feed-post record we publish matches the `com.decodingus.atmosphere.feed.post` contract
    /// the AppView's Jetstream consumer reads: `$type`, top-level `text` + `createdAt`, optional
    /// `topic` dropped when blank.
    #[test]
    fn feed_post_record_wire_shape() {
        let v = crate::feed_post_record("hello community", Some("haplogroup:R-M269"), None);
        assert_eq!(v["$type"], crate::NS_FEED_POST);
        assert_eq!(v["text"], "hello community");
        assert_eq!(v["topic"], "haplogroup:R-M269");
        assert!(v.get("createdAt").and_then(|c| c.as_str()).is_some());
        assert!(v.get("meta").is_none() && v.get("reply").is_none());

        // The code removes a blank topic.
        let v2 = crate::feed_post_record("no topic", Some("   "), None);
        assert!(v2.get("topic").is_none());
    }
}
