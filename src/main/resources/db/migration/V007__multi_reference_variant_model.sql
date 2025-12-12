-- V007: Multi-Reference Variant Model
--
-- Fixes the evidence grouping problem where the same sequencing data aligned
-- to multiple references (GRCh37, GRCh38, hs1) would be counted as multiple
-- pieces of evidence instead of one.
--
-- Key insight: Variant identity is (canonical_name, defining_haplogroup), not position.
-- Evidence is per-SOURCE, not per-alignment. Coordinates are just representations.
--
-- Hierarchy:
--   y_profile_variant      = Variant identity (what mutation)
--   y_variant_source_call  = Evidence from one source (one vote in concordance)
--   y_source_call_alignment = Coordinate representation (doesn't affect voting)

-- ============================================
-- Step 1: Add new columns to y_profile_variant
-- ============================================

-- Add variant identity columns
ALTER TABLE y_profile_variant ADD COLUMN IF NOT EXISTS canonical_name VARCHAR(100);
ALTER TABLE y_profile_variant ADD COLUMN IF NOT EXISTS naming_status VARCHAR(50) DEFAULT 'UNNAMED';

-- For novel variants, store normalized coordinates (GRCh38 as canonical)
ALTER TABLE y_profile_variant ADD COLUMN IF NOT EXISTS novel_coordinates JSONB;

-- Migrate existing data: copy variant_name to canonical_name
UPDATE y_profile_variant SET canonical_name = variant_name WHERE variant_name IS NOT NULL;

-- ============================================
-- Step 2: Create alignment table for coordinates
-- ============================================

CREATE TABLE IF NOT EXISTS y_source_call_alignment (
    id UUID PRIMARY KEY,
    source_call_id UUID NOT NULL REFERENCES y_variant_source_call(id) ON DELETE CASCADE,

    -- Reference build this alignment used
    reference_build VARCHAR(50) NOT NULL,  -- "GRCh38", "GRCh37", "hs1", "T2T-CHM13"

    -- Coordinates in this reference
    contig VARCHAR(50) NOT NULL DEFAULT 'chrY',
    position BIGINT NOT NULL,
    ref_allele VARCHAR(255) NOT NULL,
    alt_allele VARCHAR(255) NOT NULL,

    -- Call details for this alignment
    called_allele VARCHAR(255) NOT NULL,

    -- Quality metrics specific to this alignment
    read_depth INT,
    mapping_quality DOUBLE PRECISION,
    base_quality DOUBLE PRECISION,
    variant_allele_frequency DOUBLE PRECISION,

    -- For graph references (future)
    graph_node VARCHAR(100),      -- Node ID for pangenome
    graph_offset INT,             -- Offset within node

    -- Metadata
    alignment_id UUID,            -- Optional link to alignment entity
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    -- One coordinate representation per (source_call, reference_build)
    UNIQUE(source_call_id, reference_build)
);

CREATE INDEX idx_alignment_source_call ON y_source_call_alignment(source_call_id);
CREATE INDEX idx_alignment_reference ON y_source_call_alignment(reference_build);
CREATE INDEX idx_alignment_position ON y_source_call_alignment(reference_build, contig, position);

-- ============================================
-- Step 3: Migrate existing coordinate data
-- ============================================

-- Move position/ref/alt from source_call to alignment table
-- Assume existing data is GRCh38 (most common) - can be corrected later
INSERT INTO y_source_call_alignment (
    id, source_call_id, reference_build, contig, position,
    ref_allele, alt_allele, called_allele,
    read_depth, mapping_quality, variant_allele_frequency
)
SELECT
    gen_random_uuid(),
    sc.id,
    COALESCE(s.reference_build, 'GRCh38'),  -- Use source's reference_build if known
    'chrY',
    v.position,
    v.ref_allele,
    v.alt_allele,
    sc.called_allele,
    sc.read_depth,
    sc.mapping_quality,
    sc.variant_allele_frequency
FROM y_variant_source_call sc
JOIN y_profile_variant v ON sc.variant_id = v.id
JOIN y_profile_source s ON sc.source_id = s.id
WHERE v.position IS NOT NULL;

-- ============================================
-- Step 4: Update unique constraints
-- ============================================

-- Drop old position-based unique constraint
ALTER TABLE y_profile_variant DROP CONSTRAINT IF EXISTS y_profile_variant_y_profile_id_position_key;

-- Add new identity-based unique constraint for named variants
-- (canonical_name + defining_haplogroup within a profile)
-- Note: H2 doesn't support functional indexes with COALESCE, so we use a simpler index
-- and enforce uniqueness at application level for variants with NULL defining_haplogroup
CREATE INDEX IF NOT EXISTS idx_variant_identity ON y_profile_variant(
    y_profile_id,
    canonical_name,
    defining_haplogroup
);

-- For unnamed/novel variants, we use position-based index as fallback
-- Full JSONB functional indexes and partial indexes are PostgreSQL-specific
CREATE INDEX IF NOT EXISTS idx_variant_novel ON y_profile_variant(
    y_profile_id,
    position,
    canonical_name
);

-- ============================================
-- Step 5: Add helpful views
-- ============================================

-- View that shows variant with best alignment data
CREATE OR REPLACE VIEW v_variant_with_alignments AS
SELECT
    v.id AS variant_id,
    v.y_profile_id,
    v.canonical_name,
    v.defining_haplogroup,
    v.naming_status,
    v.consensus_allele,
    v.consensus_state,
    v.status,
    v.confidence_score,
    v.source_count,
    v.concordant_count,
    v.discordant_count,
    sc.id AS source_call_id,
    s.source_type,
    s.vendor,
    a.reference_build,
    a.position,
    a.ref_allele,
    a.alt_allele,
    a.called_allele,
    a.read_depth,
    a.mapping_quality
FROM y_profile_variant v
JOIN y_variant_source_call sc ON sc.variant_id = v.id
JOIN y_profile_source s ON sc.source_id = s.id
LEFT JOIN y_source_call_alignment a ON a.source_call_id = sc.id;

-- View for concordance calculation (one row per source, not per alignment)
CREATE OR REPLACE VIEW v_variant_evidence AS
SELECT
    v.id AS variant_id,
    v.y_profile_id,
    v.canonical_name,
    v.defining_haplogroup,
    v.variant_type,
    sc.id AS source_call_id,
    sc.source_id,
    s.source_type,
    s.method_tier,
    sc.called_allele,
    sc.call_state,
    sc.concordance_weight,
    sc.called_repeat_count,
    -- Aggregate best quality metrics across all alignments for this source
    MAX(a.read_depth) AS max_read_depth,
    MAX(a.mapping_quality) AS max_mapping_quality,
    COUNT(a.id) AS alignment_count
FROM y_profile_variant v
JOIN y_variant_source_call sc ON sc.variant_id = v.id
JOIN y_profile_source s ON sc.source_id = s.id
LEFT JOIN y_source_call_alignment a ON a.source_call_id = sc.id
GROUP BY v.id, v.y_profile_id, v.canonical_name, v.defining_haplogroup, v.variant_type,
         sc.id, sc.source_id, s.source_type, s.method_tier,
         sc.called_allele, sc.call_state, sc.concordance_weight, sc.called_repeat_count;

-- ============================================
-- Comments
-- ============================================

COMMENT ON TABLE y_source_call_alignment IS
'Coordinate representations of a source call. Multiple alignments of the same
sequencing data to different references are ONE piece of evidence, stored as
multiple rows here but counted once in concordance.';

COMMENT ON COLUMN y_profile_variant.canonical_name IS
'Primary variant name (e.g., M269, L21). NULL for novel/unnamed variants.
Combined with defining_haplogroup forms the variant identity.';

COMMENT ON COLUMN y_profile_variant.naming_status IS
'UNNAMED = novel variant, PENDING = submitted for naming, NAMED = has canonical name';

COMMENT ON COLUMN y_profile_variant.novel_coordinates IS
'For unnamed variants: normalized coordinates in GRCh38 for deduplication.
Format: {"GRCh38": {"position": 12345, "ref": "G", "alt": "A"}}';

COMMENT ON VIEW v_variant_evidence IS
'Flattened view for concordance calculation. Groups by SOURCE not alignment,
so multi-reference alignments count as one vote.';
