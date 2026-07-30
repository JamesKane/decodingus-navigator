-- Cached Tier B archaic SEGMENT calls for a subject, keyed to the alignment they were called from
-- rather than to the autosomal consensus: segments come from genome-wide de-novo diploid calls on
-- one alignment, not from the pooled consensus (which only carries the 1240k panel loci).
-- `source_sig` is the alignment id plus the caller's genotype version, so re-calling with a newer
-- caller invalidates the cache. `segments` is the full ArchaicSegmentResult as JSON.
CREATE TABLE consensus_archaic_segments (
    biosample_guid TEXT PRIMARY KEY REFERENCES biosample(guid),
    source_sig     TEXT NOT NULL,
    segments       TEXT NOT NULL,
    computed_at    TEXT NOT NULL
);
