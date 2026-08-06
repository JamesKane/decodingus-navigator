-- Private-Y buckets derived from a variant set — the VCF counterpart of the per-alignment
-- `private_y` analysis artifact.
--
-- Private-Y has been alignment-keyed since it was written: it walks a BAM/CRAM (or its GVCF sidecar)
-- and caches against `alignment_id`. A subject whose Y data arrived as an externally processed VCF
-- has no alignment at all, so the option was never offered — on R1b-CTS4466Plus that is ~1,600 of
-- 1,881 members. Everything the engine needs is in the call set once the import keeps it (migration
-- 0042's per-call evidence, and 0043's tree-position genotypes), so the bucket is cached here,
-- keyed the same way: the set, plus a hash of the tree sites it was classified against.
CREATE TABLE variant_set_private_y (
    variant_set_id INTEGER NOT NULL REFERENCES variant_set(id),
    -- Contig + tree site-set hash, as the alignment path's `algorithm_version` does. A changed tree
    -- re-classifies rather than serving a bucket sorted against sites that moved.
    cache_key      TEXT NOT NULL,
    -- JSON `PrivateBucket` — the same shape the alignment path stores.
    bucket         TEXT NOT NULL,
    computed_at    TEXT NOT NULL,
    PRIMARY KEY (variant_set_id, cache_key)
);

CREATE INDEX ix_variant_set_private_y_set ON variant_set_private_y(variant_set_id);
