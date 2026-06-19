-- Generalize the per-biosample Y-variant snapshot (0020) into a DNA-type-agnostic consensus profile.
-- One row per (biosample, dna_type) — 'Y' today; 'Mt' / autosomal adapters reuse the same engine and
-- table. The full reconciled profile (variants + per-source provenance + summary) lives in `payload`
-- as JSON; the scalar columns mirror the header for quick listing without decoding. `consensus_label`
-- holds the haplogroup for Y/mt (NULL where a DNA type has no lineage label).
CREATE TABLE consensus_profile (
    biosample_guid     TEXT NOT NULL REFERENCES biosample(guid),
    dna_type           TEXT NOT NULL,
    consensus_label    TEXT,
    overall_confidence REAL NOT NULL DEFAULT 0,
    source_count       INTEGER NOT NULL DEFAULT 0,
    total              INTEGER NOT NULL DEFAULT 0,
    confirmed          INTEGER NOT NULL DEFAULT 0,
    novel              INTEGER NOT NULL DEFAULT 0,
    conflict           INTEGER NOT NULL DEFAULT 0,
    single_source      INTEGER NOT NULL DEFAULT 0,
    tree_provider      TEXT,
    payload            TEXT NOT NULL,
    last_reconciled_at TEXT NOT NULL,
    PRIMARY KEY (biosample_guid, dna_type)
);

-- Carry forward any existing Y snapshots as dna_type='Y'.
INSERT INTO consensus_profile (biosample_guid, dna_type, consensus_label, overall_confidence, source_count,
    total, confirmed, novel, conflict, single_source, tree_provider, payload, last_reconciled_at)
SELECT biosample_guid, 'Y', consensus_haplogroup, overall_confidence, source_count,
    total, confirmed, novel, conflict, single_source, tree_provider, payload, last_reconciled_at
FROM y_profile;

DROP TABLE y_profile;
