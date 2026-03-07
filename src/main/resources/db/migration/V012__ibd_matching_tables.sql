-- V012: IBD Matching System tables
-- Supports match consent, match requests, and match results

-- ============================================================================
-- Match Consent — global opt-in for IBD matching per biosample
-- ============================================================================

CREATE TABLE match_consent (
    id UUID PRIMARY KEY,
    biosample_id UUID NOT NULL,
    consent_level VARCHAR(20) NOT NULL,
    allowed_match_types JSON DEFAULT '["IBD"]',
    minimum_segment_cm DOUBLE DEFAULT 7.0,
    share_contact_info BOOLEAN DEFAULT FALSE,
    consented_at TIMESTAMP NOT NULL,
    expires_at TIMESTAMP,
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(500),
    at_cid VARCHAR(100),
    version INTEGER DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_mc_biosample FOREIGN KEY (biosample_id) REFERENCES biosample(id),
    CONSTRAINT uq_mc_biosample UNIQUE (biosample_id)
);

CREATE INDEX idx_mc_biosample ON match_consent(biosample_id);

-- ============================================================================
-- Match Request — incoming/outgoing IBD comparison requests
-- ============================================================================

CREATE TABLE match_request (
    id UUID PRIMARY KEY,
    from_biosample_ref VARCHAR(500) NOT NULL,
    to_biosample_ref VARCHAR(500) NOT NULL,
    status VARCHAR(20) NOT NULL DEFAULT 'PENDING',
    request_type VARCHAR(50) NOT NULL DEFAULT 'AUTOSOMAL',
    message TEXT,
    shared_ancestor_hint VARCHAR(500),
    discovery_reason TEXT,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    expires_at TIMESTAMP,
    responded_at TIMESTAMP,
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(500),
    at_cid VARCHAR(100),
    version INTEGER DEFAULT 1,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_mr_from ON match_request(from_biosample_ref);
CREATE INDEX idx_mr_to ON match_request(to_biosample_ref);
CREATE INDEX idx_mr_status ON match_request(status);

-- ============================================================================
-- Match Result — confirmed IBD matches
-- ============================================================================

CREATE TABLE match_result (
    id UUID PRIMARY KEY,
    biosample_id UUID NOT NULL,
    matched_biosample_ref VARCHAR(500) NOT NULL,
    matched_citizen_did VARCHAR(255),
    relationship_estimate VARCHAR(50),
    total_shared_cm DOUBLE NOT NULL,
    longest_segment_cm DOUBLE,
    segment_count INTEGER NOT NULL,
    shared_segments JSON DEFAULT '[]',
    x_match_shared_cm DOUBLE,
    matched_at TIMESTAMP NOT NULL,
    attestation_hash VARCHAR(64),
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(500),
    at_cid VARCHAR(100),
    version INTEGER DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT fk_mres_biosample FOREIGN KEY (biosample_id) REFERENCES biosample(id)
);

CREATE INDEX idx_mres_biosample ON match_result(biosample_id);
CREATE INDEX idx_mres_matched ON match_result(matched_biosample_ref);
CREATE INDEX idx_mres_cm ON match_result(total_shared_cm DESC);
