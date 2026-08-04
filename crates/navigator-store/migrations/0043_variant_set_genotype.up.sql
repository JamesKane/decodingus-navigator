-- Tree-position genotypes for a variant set — the VCF counterpart of the `tree-genotype` artifact
-- the BAM/CRAM path caches per alignment.
--
-- Placement needs to know the donor's state at *every* tree position: derived, ancestral, or not
-- covered. The CRAM path gets that by genotyping the alignment at all of them and storing the
-- observed base whatever it is. A VCF import stored only the non-reference rows, so the workspace
-- held the derived calls and no record of where the donor was confidently ancestral — placement then
-- cannot tell "ancestral" from "uncovered", every backbone node reads as no-call, and scoring runs on
-- a few dozen sites instead of hundreds (observed: matched 32 of 1152 on a Big Y).
--
-- The source VCF carries exactly what is missing (an aengine Big Y export is ~218k PASS records,
-- mostly 0/0), so `source_path` records where that file is — the same role `alignment.bam_path`
-- plays for the CRAM path — and the genotypes are cached here rather than re-parsed per placement.

-- Where the set was imported from, so it can be re-read to genotype at tree positions. NULL for
-- sets imported before this, and for sources that never had a file (hand entry).
ALTER TABLE variant_set ADD COLUMN source_path TEXT;

CREATE TABLE variant_set_genotype (
    variant_set_id INTEGER NOT NULL REFERENCES variant_set(id),
    -- Same shape as the artifact `algorithm_version` the CRAM path uses: contig + source build +
    -- a hash of the sorted target positions. A changed tree changes the hash, so a stale genotype is
    -- a cache miss rather than a silently wrong placement.
    cache_key      TEXT NOT NULL,
    -- JSON `[[position, "base"], …]` — identical to the `tree-genotype` payload, so the placement
    -- path consumes a VCF-sourced genotype exactly as it does an alignment-sourced one.
    calls          TEXT NOT NULL,
    computed_at    TEXT NOT NULL,
    PRIMARY KEY (variant_set_id, cache_key)
);

CREATE INDEX ix_variant_set_genotype_set ON variant_set_genotype(variant_set_id);
