-- V003__sync_queue_tables.sql
-- PDS sync queue and history tables for async synchronization

-- ============================================
-- SYNC_QUEUE: Pending sync operations
-- ============================================
-- Supports indefinite offline operation - no max_attempts limit
-- Queue persists forever until sync succeeds or user cancels

CREATE TABLE sync_queue (
    id UUID PRIMARY KEY,

    -- Target entity
    entity_type VARCHAR(50) NOT NULL,
    entity_id UUID NOT NULL,

    -- Operation
    operation VARCHAR(20) NOT NULL,

    -- Queue state
    status VARCHAR(20) NOT NULL DEFAULT 'PENDING',
    priority INT NOT NULL DEFAULT 5,

    -- Timing
    queued_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP,
    completed_at TIMESTAMP,

    -- Retry handling (no max_attempts - queue indefinitely for offline support)
    attempt_count INT NOT NULL DEFAULT 0,
    next_retry_at TIMESTAMP,
    last_error TEXT,

    -- Payload snapshot (JSON of entity state at queue time)
    payload_snapshot JSON,

    -- Metadata
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT chk_sq_entity_type CHECK (entity_type IN (
        'BIOSAMPLE', 'PROJECT', 'SEQUENCE_RUN', 'ALIGNMENT',
        'STR_PROFILE', 'CHIP_PROFILE', 'Y_SNP_PANEL', 'HAPLOGROUP_RECONCILIATION'
    )),
    CONSTRAINT chk_sq_operation CHECK (operation IN ('CREATE', 'UPDATE', 'DELETE')),
    CONSTRAINT chk_sq_status CHECK (status IN ('PENDING', 'IN_PROGRESS', 'COMPLETED', 'FAILED', 'CANCELLED'))
);

CREATE INDEX idx_sync_queue_status ON sync_queue(status, priority, queued_at);
CREATE INDEX idx_sync_queue_entity ON sync_queue(entity_type, entity_id);
CREATE INDEX idx_sync_queue_next_retry ON sync_queue(next_retry_at) WHERE status = 'PENDING';

-- ============================================
-- SYNC_HISTORY: Audit trail of sync operations
-- ============================================
CREATE TABLE sync_history (
    id UUID PRIMARY KEY,

    -- Target entity
    entity_type VARCHAR(50) NOT NULL,
    entity_id UUID NOT NULL,
    at_uri VARCHAR(512),

    -- Operation details
    operation VARCHAR(20) NOT NULL,
    direction VARCHAR(10) NOT NULL,

    -- Result
    status VARCHAR(20) NOT NULL,
    error_message TEXT,

    -- Timing
    started_at TIMESTAMP NOT NULL,
    completed_at TIMESTAMP NOT NULL,
    duration_ms BIGINT,

    -- Version tracking
    local_version_before INT,
    local_version_after INT,
    remote_version_before INT,
    remote_version_after INT,

    -- CID tracking for AT Protocol
    local_cid_before VARCHAR(128),
    local_cid_after VARCHAR(128),
    remote_cid VARCHAR(128),

    -- Metadata
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT chk_sh_entity_type CHECK (entity_type IN (
        'BIOSAMPLE', 'PROJECT', 'SEQUENCE_RUN', 'ALIGNMENT',
        'STR_PROFILE', 'CHIP_PROFILE', 'Y_SNP_PANEL', 'HAPLOGROUP_RECONCILIATION'
    )),
    CONSTRAINT chk_sh_operation CHECK (operation IN ('CREATE', 'UPDATE', 'DELETE')),
    CONSTRAINT chk_sh_direction CHECK (direction IN ('PUSH', 'PULL')),
    CONSTRAINT chk_sh_status CHECK (status IN ('SUCCESS', 'FAILED', 'CONFLICT', 'SKIPPED'))
);

CREATE INDEX idx_sync_history_entity ON sync_history(entity_type, entity_id);
CREATE INDEX idx_sync_history_time ON sync_history(completed_at DESC);
CREATE INDEX idx_sync_history_status ON sync_history(status);

-- ============================================
-- SYNC_CONFLICT: Unresolved conflicts
-- ============================================
CREATE TABLE sync_conflict (
    id UUID PRIMARY KEY,

    -- Target entity
    entity_type VARCHAR(50) NOT NULL,
    entity_id UUID NOT NULL,
    at_uri VARCHAR(512),

    -- Conflict details
    detected_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    local_version INT NOT NULL,
    remote_version INT NOT NULL,

    -- Changed fields (JSON arrays)
    local_changes JSON,
    remote_changes JSON,
    overlapping_fields JSON,

    -- Suggested resolution
    suggested_resolution VARCHAR(30),
    resolution_reason TEXT,

    -- Resolution state
    status VARCHAR(20) NOT NULL DEFAULT 'UNRESOLVED',
    resolved_at TIMESTAMP,
    resolved_by VARCHAR(50),
    resolution_action VARCHAR(30),

    -- Snapshots for comparison
    local_snapshot JSON,
    remote_snapshot JSON,

    -- Metadata
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT chk_sc_entity_type CHECK (entity_type IN (
        'BIOSAMPLE', 'PROJECT', 'SEQUENCE_RUN', 'ALIGNMENT',
        'STR_PROFILE', 'CHIP_PROFILE', 'Y_SNP_PANEL', 'HAPLOGROUP_RECONCILIATION'
    )),
    CONSTRAINT chk_sc_suggested CHECK (suggested_resolution IS NULL OR suggested_resolution IN (
        'KEEP_LOCAL', 'ACCEPT_REMOTE', 'MERGE', 'MANUAL'
    )),
    CONSTRAINT chk_sc_status CHECK (status IN ('UNRESOLVED', 'RESOLVED', 'DISMISSED')),
    CONSTRAINT chk_sc_action CHECK (resolution_action IS NULL OR resolution_action IN (
        'KEPT_LOCAL', 'ACCEPTED_REMOTE', 'MERGED', 'MANUAL_EDIT'
    ))
);

CREATE INDEX idx_sync_conflict_entity ON sync_conflict(entity_type, entity_id);
CREATE INDEX idx_sync_conflict_status ON sync_conflict(status);
CREATE INDEX idx_sync_conflict_unresolved ON sync_conflict(detected_at) WHERE status = 'UNRESOLVED';

-- ============================================
-- SYNC_SETTINGS: User sync preferences
-- ============================================
CREATE TABLE sync_settings (
    id UUID PRIMARY KEY,

    -- Sync behavior
    incoming_sync_enabled BOOLEAN NOT NULL DEFAULT TRUE,
    incoming_sync_interval_minutes INT NOT NULL DEFAULT 60,
    auto_resolve_non_overlapping BOOLEAN NOT NULL DEFAULT TRUE,
    auto_accept_appview_updates BOOLEAN NOT NULL DEFAULT TRUE,

    -- Notification preferences
    notify_on_conflict BOOLEAN NOT NULL DEFAULT TRUE,
    notify_on_sync_error BOOLEAN NOT NULL DEFAULT TRUE,

    -- Last sync timestamps
    last_outgoing_sync_at TIMESTAMP,
    last_incoming_sync_at TIMESTAMP,
    last_successful_sync_at TIMESTAMP,

    -- Metadata
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

-- Insert default settings (singleton row)
INSERT INTO sync_settings (id, incoming_sync_enabled, incoming_sync_interval_minutes)
VALUES ('00000000-0000-0000-0000-000000000001', TRUE, 60);
