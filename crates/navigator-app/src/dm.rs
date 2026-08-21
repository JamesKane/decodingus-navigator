//! `impl App` methods for **peer direct messages** (social roadmap 3a).
//!
//! This module is a thin DM layer on the generic D1 encrypted exchange, which is `ibd_exchange.rs`
//! and `navigator_sync::exchange`. It uses the cryptography, the discovery, the consent, and the
//! relay of that exchange without a change.
//!
//! This module adds three parts: the DM `purpose`, storage for the session key, and a store for
//! each message. Storage of the session key makes a conversation asynchronous, and the conversation
//! continues after a restart. The AppView relays only ciphertext. The plaintext of a message never
//! leaves the device.
use super::*;

use base64::engine::general_purpose::STANDARD;
use base64::Engine;

/// The exchange `purpose` tag for a peer DM. The AppView makes the title of its consent
/// notification from this tag. An IBD request uses a tag that starts with `IBD_`. So a filter on
/// the tag keeps the DM inbox separate from the IBD tab.
pub const DM_PURPOSE: &str = "GENEALOGY_PII";

impl App {
    /// Open a DM request to a partner DID (paste-a-DID / "Message this match"). Returns the request
    /// URI; the partner sees it in their DM inbox ([`dm_incoming`](Self::dm_incoming)) and consents.
    pub async fn dm_initiate(&self, partner_did: &str) -> Result<String, AppError> {
        self.exchange_request(partner_did, DM_PURPOSE, None).await
    }

    /// The DM requests that arrived and that need our consent. The view is symmetric-blind, and it
    /// holds only requests with the DM purpose.
    pub async fn dm_incoming(&self) -> Result<Vec<IncomingRequest>, AppError> {
        Ok(self
            .exchange_incoming()
            .await?
            .into_iter()
            .filter(|r| r.purpose == DM_PURPOSE)
            .collect())
    }

    /// The DM sessions that both parties agreed to, but that we did not connect. These sessions
    /// have no conversation and no key in the store. The list holds only sessions with the DM
    /// purpose.
    pub async fn dm_ready(&self) -> Result<Vec<ExchangeSessionInfo>, AppError> {
        let mut out = Vec::new();
        for info in self.exchange_pending().await? {
            if info.purpose != DM_PURPOSE {
                continue;
            }
            if navigator_store::dm::get_conversation(self.store.pool(), &info.session_id)
                .await?
                .is_none()
            {
                out.push(info);
            }
        }
        Ok(out)
    }

    /// Consent to (or decline) a DM request.
    pub async fn dm_consent(&self, request_uri: &str, given: bool) -> Result<ConsentOutcome, AppError> {
        self.exchange_consent(request_uri, given).await
    }

    /// Connect a DM session that both parties agreed to.
    ///
    /// The function does the X3DH-lite handshake. Both peers must be online for this one step. It
    /// then writes the derived session key and the partner to the store. So each later send and
    /// receive is asynchronous, and it continues after a restart.
    ///
    /// A second call is safe. It makes a new key, and the seq counters of the conversation keep
    /// their values.
    pub async fn dm_connect(&self, info: &ExchangeSessionInfo) -> Result<(), AppError> {
        let did = self.require_account()?;
        let session = self.open_exchange_session(info).await?;
        let key_b64 = STANDARD.encode(session.key);
        navigator_store::dm::upsert_conversation(
            self.store.pool(),
            &info.session_id,
            &info.request_uri,
            &did,
            &info.partner_did,
            DM_PURPOSE,
            &key_b64,
            &Utc::now().to_rfc3339(),
        )
        .await?;
        Ok(())
    }

    /// All persisted DM conversations for the signed-in account, newest activity first.
    pub async fn dm_conversations(&self) -> Result<Vec<DmConversationSummary>, AppError> {
        let did = self.require_account()?;
        Ok(navigator_store::dm::list_conversations(self.store.pool(), &did).await?)
    }

    /// The transcript for a conversation (oldest first); also marks it read.
    pub async fn dm_messages(&self, session_id: &str) -> Result<Vec<DmMessage>, AppError> {
        let convo = self.dm_conversation_or_err(session_id).await?;
        navigator_store::dm::set_last_read(
            self.store.pool(),
            session_id,
            &convo.partner_did,
            &Utc::now().to_rfc3339(),
        )
        .await?;
        Ok(navigator_store::dm::messages(self.store.pool(), session_id).await?)
    }

    /// Encrypt a message, relay it on an open conversation, and write it to the local store. The
    /// function returns the seq of the message.
    pub async fn dm_send(&self, session_id: &str, text: &str) -> Result<i64, AppError> {
        let convo = self.dm_conversation_or_err(session_id).await?;
        let session = self.rebuild_session(&convo)?;
        let now = Utc::now().to_rfc3339();
        let seq = navigator_store::dm::take_send_seq(self.store.pool(), session_id, &now)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("dm conversation {session_id}"))))?;
        self.exchange_send(&session, seq as i32, text.as_bytes()).await?;
        navigator_store::dm::insert_message(self.store.pool(), session_id, &convo.my_did, seq, text, &now).await?;
        Ok(seq)
    }

    /// Read, decrypt, store, and acknowledge each message that a conversation holds. The function
    /// returns the count of new messages. It does not count a duplicate message.
    pub async fn dm_sync(&self, session_id: &str) -> Result<usize, AppError> {
        let did = self.require_account()?;
        let convo = self.dm_conversation_or_err(session_id).await?;
        let session = self.rebuild_session(&convo)?;
        let mut new = 0usize;
        for env in self.exchange_relay_pull(session_id).await? {
            let Ok(parsed) = exchange::Envelope::from_blob(&env.blob) else {
                continue;
            };
            // An old handshake has seq 0, and the code can not decrypt it as data. Acknowledge
            // the handshake and remove it.
            if matches!(parsed, exchange::Envelope::Handshake { .. }) {
                let _ = self.exchange_relay_ack(env.id).await;
                continue;
            }
            // The AAD binds the route of the sender. `from` is the partner, and `to` is us.
            let aad = exchange::relay_aad(session_id, &env.from_did, &did, env.seq);
            let Ok(pt) = exchange::open(&session.key, &aad, &parsed) else {
                continue; // not for this session key / tampered — leave un-acked
            };
            let body = String::from_utf8_lossy(&pt);
            let inserted = navigator_store::dm::insert_message(
                self.store.pool(),
                session_id,
                &env.from_did,
                env.seq as i64,
                &body,
                &Utc::now().to_rfc3339(),
            )
            .await?;
            if inserted {
                new += 1;
            }
            let _ = self.exchange_relay_ack(env.id).await;
        }
        Ok(new)
    }

    /// Load a conversation or map it to a not-found store error.
    async fn dm_conversation_or_err(&self, session_id: &str) -> Result<navigator_store::dm::DmConversation, AppError> {
        navigator_store::dm::get_conversation(self.store.pool(), session_id)
            .await?
            .ok_or_else(|| AppError::Store(StoreError::NotFound(format!("dm conversation {session_id}"))))
    }

    /// Rebuild an [`EstablishedSession`] from a persisted conversation's stored key.
    fn rebuild_session(&self, convo: &navigator_store::dm::DmConversation) -> Result<EstablishedSession, AppError> {
        let bytes = STANDARD
            .decode(convo.session_key.trim())
            .map_err(|e| AppError::Import(format!("corrupt DM session key: {e}")))?;
        let key: [u8; 32] = bytes
            .try_into()
            .map_err(|_| AppError::Import("DM session key is not 32 bytes".into()))?;
        Ok(EstablishedSession {
            session_id: convo.session_id.clone(),
            partner_did: convo.partner_did.clone(),
            key,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The session key survives a base64 round-trip through the store key column, so a rebuilt
    /// session decrypts what the original sealed (no live AppView needed).
    #[test]
    fn session_key_roundtrips_through_storage() {
        let key: [u8; 32] = std::array::from_fn(|i| i as u8);
        let stored = STANDARD.encode(key);
        let back: [u8; 32] = STANDARD.decode(&stored).unwrap().try_into().unwrap();
        assert_eq!(key, back);

        // The value also works as an AES session key. Seal the data here, then open the data with
        // the key that the code makes again.
        let aad = exchange::relay_aad("sess", "did:plc:a", "did:plc:b", 1);
        let blob = exchange::seal(&key, &aad, b"hello").and_then(|e| e.to_blob()).unwrap();
        let parsed = exchange::Envelope::from_blob(&blob).unwrap();
        let pt = exchange::open(&back, &aad, &parsed).unwrap();
        assert_eq!(pt, b"hello");
    }
}
