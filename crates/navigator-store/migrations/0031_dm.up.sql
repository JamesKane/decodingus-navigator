-- Peer direct messages over the D1 encrypted relay (social roadmap 3a). The crypto + relay live in
-- navigator-sync::exchange / the AppView broker (which only ever sees ciphertext); this stores the
-- LOCAL, decrypted side of a conversation so it survives restart and stays async.
--
-- Privacy: like the MDKA table, the message bodies are plaintext and citizen-private — they are
-- NEVER federated/published and NEVER reach the AppView (which is a blind ciphertext relay). The
-- per-session symmetric key persists here so a once-established session can send/receive without a
-- fresh handshake; this matches the locally-stored-plaintext posture of the rest of the workspace DB.
CREATE TABLE dm_conversation (
    session_id    TEXT PRIMARY KEY,        -- the broker session UUID (one conversation per session)
    request_uri   TEXT NOT NULL,           -- the exchange request that opened it
    my_did        TEXT NOT NULL,           -- the local account DID (sender of outgoing messages)
    partner_did   TEXT NOT NULL,
    purpose       TEXT NOT NULL,           -- 'GENEALOGY_PII' for peer DMs
    session_key   TEXT NOT NULL,           -- base64 of the 32-byte derived AES key
    next_send_seq INTEGER NOT NULL DEFAULT 1,  -- monotonic outgoing seq (handshake was seq 0)
    last_read_seq INTEGER NOT NULL DEFAULT 0,  -- highest partner seq the user has seen (unread = above this)
    created_at    TEXT NOT NULL,           -- ISO-8601
    updated_at    TEXT NOT NULL
);

CREATE TABLE dm_message (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    session_id TEXT NOT NULL REFERENCES dm_conversation(session_id) ON DELETE CASCADE,
    from_did   TEXT NOT NULL,              -- == my_did for outgoing, partner_did for incoming
    seq        INTEGER NOT NULL,           -- the relay seq the body rode under
    body       TEXT NOT NULL,              -- decrypted plaintext (local-only)
    created_at TEXT NOT NULL,
    UNIQUE (session_id, from_did, seq)     -- idempotent re-delivery + restart-safe dedupe
);
CREATE INDEX idx_dm_message_session ON dm_message(session_id, seq);
