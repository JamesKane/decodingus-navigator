-- Vendor-neutral Subject identity (FTDNA project-import design §4.2). The platform is NOT locked to
-- FTDNA — a kit number is one identifier among many. `(source, external_id)` is the global
-- cross-project dedup anchor the matching engine keys on.
--
-- PII / never-federated: `external_id` and `ftdna_member` hold member-linkage and FTDNA-reported
-- labels. They must NEVER be derived into a public PDS `fed` record nor put in an AppView-bound
-- payload; they may only enter the encrypted Edge-to-Edge tier. Keep distinct from our own computed
-- haplogroup calls (`haplogroup_call`) — different provenance.
CREATE TABLE external_id (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    biosample_guid TEXT NOT NULL REFERENCES biosample(guid),
    source         TEXT NOT NULL,        -- 'FTDNA' | 'YSEQ' | 'NEBULA' | 'WGS' | 'MANUAL' | ...
    external_id    TEXT NOT NULL,        -- kit number / vendor id
    UNIQUE (source, external_id)         -- one Subject per (vendor, id)
);
CREATE INDEX idx_external_id_lookup ON external_id(source, external_id);
CREATE INDEX idx_external_id_biosample ON external_id(biosample_guid);

-- FTDNA-reported member labels only (ancestry → mdka, §4.3). Computed haplogroups stay in
-- haplogroup_call.
CREATE TABLE ftdna_member (
    biosample_guid      TEXT PRIMARY KEY REFERENCES biosample(guid),
    member_name         TEXT,            -- member/contact name as in GAP (PII)
    y_haplogroup_ftdna  TEXT,            -- as reported by FTDNA (label only)
    mt_haplogroup_ftdna TEXT,
    haplo_status        TEXT,            -- 'predicted' | 'confirmed'
    access_granted      TEXT,            -- 'Advanced' | 'Limited' | 'None' — pose-as gate + Big Y data tier
    publicly_shares     INTEGER          -- consent flag (col 11): 0/1; gates federation
);
