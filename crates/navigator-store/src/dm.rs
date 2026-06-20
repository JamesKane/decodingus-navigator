//! Peer direct-message persistence (social roadmap 3a) — the LOCAL, decrypted side of an encrypted
//! conversation carried over the D1 relay. One [`DmConversation`] per broker session (it holds the
//! persisted session key + the outgoing seq counter so messaging is async and restart-safe), and one
//! [`DmMessage`] per relayed line. Bodies are plaintext and citizen-private: never federated, never
//! sent to the AppView (which only ever relays ciphertext). Mirrors the `ibd_exchange` store style.

use sqlx::SqlitePool;

use crate::StoreError;

/// A persisted DM conversation (keyed by the broker session id).
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct DmConversation {
    pub session_id: String,
    pub request_uri: String,
    pub my_did: String,
    pub partner_did: String,
    pub purpose: String,
    /// Base64 of the 32-byte derived AES session key.
    pub session_key: String,
    pub next_send_seq: i64,
    pub last_read_seq: i64,
    pub created_at: String,
    pub updated_at: String,
}

/// One stored message line (decrypted). `outgoing` is `from_did == my_did`, resolved by the caller.
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct DmMessage {
    pub id: i64,
    pub session_id: String,
    pub from_did: String,
    pub seq: i64,
    pub body: String,
    pub created_at: String,
}

/// A conversation row as the list view needs it: the conversation plus its newest body and the
/// count of unread (partner) messages above `last_read_seq`.
#[derive(Debug, Clone, PartialEq, sqlx::FromRow)]
pub struct DmConversationSummary {
    pub session_id: String,
    pub partner_did: String,
    pub purpose: String,
    pub last_body: Option<String>,
    pub last_at: Option<String>,
    pub unread: i64,
    pub updated_at: String,
}

/// Insert (or update the mutable fields of) a conversation. Called when a session is established:
/// the session key + partner are set; the seq counters keep their existing values on re-connect.
#[allow(clippy::too_many_arguments)]
pub async fn upsert_conversation(
    pool: &SqlitePool,
    session_id: &str,
    request_uri: &str,
    my_did: &str,
    partner_did: &str,
    purpose: &str,
    session_key: &str,
    now: &str,
) -> Result<(), StoreError> {
    sqlx::query(
        "INSERT INTO dm_conversation \
           (session_id, request_uri, my_did, partner_did, purpose, session_key, next_send_seq, last_read_seq, created_at, updated_at) \
         VALUES (?, ?, ?, ?, ?, ?, 1, 0, ?, ?) \
         ON CONFLICT(session_id) DO UPDATE SET \
           request_uri = excluded.request_uri, my_did = excluded.my_did, partner_did = excluded.partner_did, \
           purpose = excluded.purpose, session_key = excluded.session_key, updated_at = excluded.updated_at",
    )
    .bind(session_id)
    .bind(request_uri)
    .bind(my_did)
    .bind(partner_did)
    .bind(purpose)
    .bind(session_key)
    .bind(now)
    .bind(now)
    .execute(pool)
    .await?;
    Ok(())
}

/// Fetch one conversation by session id.
pub async fn get_conversation(pool: &SqlitePool, session_id: &str) -> Result<Option<DmConversation>, StoreError> {
    Ok(
        sqlx::query_as::<_, DmConversation>("SELECT * FROM dm_conversation WHERE session_id = ?")
            .bind(session_id)
            .fetch_optional(pool)
            .await?,
    )
}

/// All conversations for an account, newest activity first, with last-message + unread counts.
pub async fn list_conversations(pool: &SqlitePool, my_did: &str) -> Result<Vec<DmConversationSummary>, StoreError> {
    Ok(sqlx::query_as::<_, DmConversationSummary>(
        "SELECT c.session_id, c.partner_did, c.purpose, \
                m.body AS last_body, m.created_at AS last_at, \
                (SELECT COUNT(*) FROM dm_message u \
                   WHERE u.session_id = c.session_id AND u.from_did = c.partner_did AND u.seq > c.last_read_seq) AS unread, \
                c.updated_at \
         FROM dm_conversation c \
         LEFT JOIN dm_message m ON m.id = ( \
           SELECT id FROM dm_message mm WHERE mm.session_id = c.session_id ORDER BY mm.created_at DESC, mm.id DESC LIMIT 1) \
         WHERE c.my_did = ? \
         ORDER BY COALESCE(m.created_at, c.updated_at) DESC",
    )
    .bind(my_did)
    .fetch_all(pool)
    .await?)
}

/// Insert a message line, ignoring a duplicate `(session_id, from_did, seq)` (idempotent re-delivery).
/// Returns `true` if a new row was inserted.
pub async fn insert_message(
    pool: &SqlitePool,
    session_id: &str,
    from_did: &str,
    seq: i64,
    body: &str,
    created_at: &str,
) -> Result<bool, StoreError> {
    let res = sqlx::query(
        "INSERT OR IGNORE INTO dm_message (session_id, from_did, seq, body, created_at) VALUES (?, ?, ?, ?, ?)",
    )
    .bind(session_id)
    .bind(from_did)
    .bind(seq)
    .bind(body)
    .bind(created_at)
    .execute(pool)
    .await?;
    Ok(res.rows_affected() > 0)
}

/// The full transcript for a conversation, oldest first.
pub async fn messages(pool: &SqlitePool, session_id: &str) -> Result<Vec<DmMessage>, StoreError> {
    Ok(
        sqlx::query_as::<_, DmMessage>("SELECT * FROM dm_message WHERE session_id = ? ORDER BY created_at ASC, id ASC")
            .bind(session_id)
            .fetch_all(pool)
            .await?,
    )
}

/// Advance the outgoing seq counter by one and return the seq the caller should use for THIS send
/// (the value before the bump). `None` if the conversation is unknown.
pub async fn take_send_seq(pool: &SqlitePool, session_id: &str, now: &str) -> Result<Option<i64>, StoreError> {
    let mut tx = pool.begin().await?;
    let seq: Option<i64> = sqlx::query_scalar("SELECT next_send_seq FROM dm_conversation WHERE session_id = ?")
        .bind(session_id)
        .fetch_optional(&mut *tx)
        .await?;
    let Some(seq) = seq else {
        return Ok(None);
    };
    sqlx::query("UPDATE dm_conversation SET next_send_seq = next_send_seq + 1, updated_at = ? WHERE session_id = ?")
        .bind(now)
        .bind(session_id)
        .execute(&mut *tx)
        .await?;
    tx.commit().await?;
    Ok(Some(seq))
}

/// Mark a conversation read up to the partner's highest stored seq (clears the unread count).
pub async fn set_last_read(
    pool: &SqlitePool,
    session_id: &str,
    partner_did: &str,
    now: &str,
) -> Result<(), StoreError> {
    sqlx::query(
        "UPDATE dm_conversation SET last_read_seq = COALESCE( \
           (SELECT MAX(seq) FROM dm_message WHERE session_id = ? AND from_did = ?), last_read_seq), \
           updated_at = ? \
         WHERE session_id = ?",
    )
    .bind(session_id)
    .bind(partner_did)
    .bind(now)
    .bind(session_id)
    .execute(pool)
    .await?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const ME: &str = "did:plc:me";
    const PARTNER: &str = "did:plc:partner";

    async fn convo(pool: &SqlitePool) {
        upsert_conversation(
            pool,
            "sess-1",
            "exchange:req-1",
            ME,
            PARTNER,
            "GENEALOGY_PII",
            "a2V5",
            "2026-06-20T00:00:00Z",
        )
        .await
        .unwrap();
    }

    #[tokio::test]
    async fn upsert_keeps_seq_counters_on_reconnect() {
        let store = crate::Store::open_in_memory().await.unwrap();
        let pool = store.pool();
        convo(pool).await;
        // Advance the send seq, then re-connect (re-upsert): the counter must not reset to 1.
        assert_eq!(take_send_seq(pool, "sess-1", "t").await.unwrap(), Some(1));
        assert_eq!(take_send_seq(pool, "sess-1", "t").await.unwrap(), Some(2));
        upsert_conversation(
            pool,
            "sess-1",
            "exchange:req-1",
            ME,
            PARTNER,
            "GENEALOGY_PII",
            "bmV3a2V5",
            "t2",
        )
        .await
        .unwrap();
        let c = get_conversation(pool, "sess-1").await.unwrap().unwrap();
        assert_eq!(c.next_send_seq, 3); // preserved
        assert_eq!(c.session_key, "bmV3a2V5"); // key refreshed
    }

    #[tokio::test]
    async fn messages_dedupe_and_unread_then_read() {
        let store = crate::Store::open_in_memory().await.unwrap();
        let pool = store.pool();
        convo(pool).await;
        // Outgoing + two incoming.
        assert!(insert_message(pool, "sess-1", ME, 1, "hi", "2026-06-20T00:01:00Z")
            .await
            .unwrap());
        assert!(
            insert_message(pool, "sess-1", PARTNER, 1, "hello", "2026-06-20T00:02:00Z")
                .await
                .unwrap()
        );
        assert!(
            insert_message(pool, "sess-1", PARTNER, 2, "again", "2026-06-20T00:03:00Z")
                .await
                .unwrap()
        );
        // Re-delivery of partner seq 1 is ignored.
        assert!(
            !insert_message(pool, "sess-1", PARTNER, 1, "hello", "2026-06-20T00:02:00Z")
                .await
                .unwrap()
        );
        assert_eq!(messages(pool, "sess-1").await.unwrap().len(), 3);

        let before = list_conversations(pool, ME).await.unwrap();
        assert_eq!(before.len(), 1);
        assert_eq!(before[0].unread, 2); // two partner messages above last_read_seq=0
        assert_eq!(before[0].last_body.as_deref(), Some("again"));

        set_last_read(pool, "sess-1", PARTNER, "t").await.unwrap();
        let after = list_conversations(pool, ME).await.unwrap();
        assert_eq!(after[0].unread, 0);
    }
}
