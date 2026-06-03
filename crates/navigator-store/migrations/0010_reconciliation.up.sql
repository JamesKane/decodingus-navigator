-- Donor-level reconciliation: manual override of the consensus + an audit log.
CREATE TABLE reconciliation_override (
    biosample_guid TEXT NOT NULL,
    dna_type       TEXT NOT NULL,   -- 'Y' | 'Mt'
    haplogroup     TEXT NOT NULL,
    reason         TEXT,
    PRIMARY KEY (biosample_guid, dna_type)
);
CREATE TABLE reconciliation_audit (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid TEXT NOT NULL,
    dna_type       TEXT NOT NULL,
    ts             TEXT NOT NULL,   -- RFC3339
    action         TEXT NOT NULL,
    note           TEXT NOT NULL
);
CREATE INDEX idx_recon_audit ON reconciliation_audit(biosample_guid, dna_type);
