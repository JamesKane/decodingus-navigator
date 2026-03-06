-- V006: Y Chromosome Unified Profile System
-- Replaces HaplogroupReconciliation for Y-DNA
-- Combines SNPs, INDELs, MNPs, and STRs with quality-weighted concordance

-- ============================================
-- Y_CHROMOSOME_PROFILE: Main entity per biosample
-- ============================================

CREATE TABLE y_chromosome_profile (
    id UUID PRIMARY KEY,
    biosample_id UUID NOT NULL REFERENCES biosample(id) ON DELETE CASCADE,

    -- Consensus haplogroup derived from unified variant calls
    consensus_haplogroup VARCHAR(100),
    haplogroup_confidence DOUBLE PRECISION,
    haplogroup_tree_provider VARCHAR(50),
    haplogroup_tree_version VARCHAR(50),

    -- Summary statistics - SNPs/INDELs
    total_variants INT NOT NULL DEFAULT 0,
    confirmed_count INT NOT NULL DEFAULT 0,
    novel_count INT NOT NULL DEFAULT 0,
    conflict_count INT NOT NULL DEFAULT 0,
    no_coverage_count INT NOT NULL DEFAULT 0,

    -- Summary statistics - STRs
    str_marker_count INT NOT NULL DEFAULT 0,
    str_confirmed_count INT NOT NULL DEFAULT 0,

    -- Quality metrics
    overall_confidence DOUBLE PRECISION,
    callable_region_pct DOUBLE PRECISION,
    mean_coverage DOUBLE PRECISION,

    -- Source tracking
    source_count INT NOT NULL DEFAULT 0,
    primary_source_type VARCHAR(30),

    -- Timestamps
    last_reconciled_at TIMESTAMP,

    -- Sync tracking
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(512) UNIQUE,
    at_cid VARCHAR(128),

    -- Versioning
    version INT NOT NULL DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    -- Ensure one profile per biosample
    CONSTRAINT uq_y_profile_biosample UNIQUE (biosample_id),
    CONSTRAINT chk_y_profile_sync CHECK (sync_status IN ('Local', 'Synced', 'Modified', 'Conflict'))
);

CREATE INDEX idx_y_profile_biosample ON y_chromosome_profile(biosample_id);
CREATE INDEX idx_y_profile_haplogroup ON y_chromosome_profile(consensus_haplogroup);
CREATE INDEX idx_y_profile_sync ON y_chromosome_profile(sync_status);

-- ============================================
-- Y_PROFILE_SOURCE: Contributing test sources
-- ============================================

CREATE TABLE y_profile_source (
    id UUID PRIMARY KEY,
    y_profile_id UUID NOT NULL REFERENCES y_chromosome_profile(id) ON DELETE CASCADE,

    -- Source identification
    source_type VARCHAR(30) NOT NULL,
    source_ref VARCHAR(512),

    -- Provider metadata
    vendor VARCHAR(100),
    test_name VARCHAR(100),
    test_date TIMESTAMP,

    -- Method-based quality tier
    -- SNPs: SANGER=5, WGS_LONG_READ=4, WGS_SHORT_READ=3, TARGETED_NGS=2, CHIP=1, MANUAL=0
    -- STRs: CAPILLARY_ELECTROPHORESIS=5, SANGER=4, WGS_LONG_READ=3, WGS_SHORT_READ=2, CHIP=1, MANUAL=0
    method_tier INT NOT NULL DEFAULT 1,

    -- Quality metrics
    mean_read_depth DOUBLE PRECISION,
    mean_mapping_quality DOUBLE PRECISION,
    coverage_pct DOUBLE PRECISION,

    -- Contribution tracking
    variant_count INT NOT NULL DEFAULT 0,
    str_marker_count INT NOT NULL DEFAULT 0,
    novel_variant_count INT NOT NULL DEFAULT 0,

    -- Link to alignment (for sequencing sources)
    alignment_id UUID REFERENCES alignment(id) ON DELETE SET NULL,
    reference_build VARCHAR(50),

    -- Timestamps
    imported_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT chk_source_type CHECK (source_type IN (
        'SANGER', 'CAPILLARY_ELECTROPHORESIS', 'WGS_LONG_READ', 'WGS_SHORT_READ',
        'TARGETED_NGS', 'CHIP', 'MANUAL'
    ))
);

CREATE INDEX idx_y_source_profile ON y_profile_source(y_profile_id);
CREATE INDEX idx_y_source_type ON y_profile_source(source_type);
CREATE INDEX idx_y_source_alignment ON y_profile_source(alignment_id);

-- ============================================
-- Y_PROFILE_REGION: Callable regions per source
-- ============================================

CREATE TABLE y_profile_region (
    id UUID PRIMARY KEY,
    y_profile_id UUID NOT NULL REFERENCES y_chromosome_profile(id) ON DELETE CASCADE,
    source_id UUID NOT NULL REFERENCES y_profile_source(id) ON DELETE CASCADE,

    -- Region definition (GRCh38)
    contig VARCHAR(50) NOT NULL DEFAULT 'chrY',
    start_position BIGINT NOT NULL,
    end_position BIGINT NOT NULL,

    -- Coverage quality
    callable_state VARCHAR(30) NOT NULL,
    mean_coverage DOUBLE PRECISION,
    mean_mapping_quality DOUBLE PRECISION,

    -- Link to callable_loci_cache if available
    callable_loci_cache_id UUID REFERENCES callable_loci_cache(id) ON DELETE SET NULL,

    CONSTRAINT chk_region_state CHECK (callable_state IN (
        'CALLABLE', 'NO_COVERAGE', 'LOW_COVERAGE',
        'EXCESSIVE_COVERAGE', 'POOR_MAPPING_QUALITY', 'REF_N', 'SUMMARY'
    ))
);

CREATE INDEX idx_y_region_profile ON y_profile_region(y_profile_id);
CREATE INDEX idx_y_region_source ON y_profile_region(source_id);
CREATE INDEX idx_y_region_position ON y_profile_region(start_position, end_position);

-- ============================================
-- Y_PROFILE_VARIANT: Individual variant calls with concordance
-- ============================================

CREATE TABLE y_profile_variant (
    id UUID PRIMARY KEY,
    y_profile_id UUID NOT NULL REFERENCES y_chromosome_profile(id) ON DELETE CASCADE,

    -- Variant identification (GRCh38)
    contig VARCHAR(50) NOT NULL DEFAULT 'chrY',
    position BIGINT NOT NULL,
    end_position BIGINT,
    ref_allele VARCHAR(500) NOT NULL,
    alt_allele VARCHAR(500) NOT NULL,

    -- Variant metadata
    variant_type VARCHAR(20) NOT NULL DEFAULT 'SNP',
    variant_name VARCHAR(100),
    rs_id VARCHAR(50),

    -- STR-specific fields
    marker_name VARCHAR(100),
    repeat_count INT,
    str_metadata JSON,

    -- Consensus call
    consensus_allele VARCHAR(500),
    consensus_state VARCHAR(30) NOT NULL DEFAULT 'NO_CALL',

    -- Variant status
    status VARCHAR(30) NOT NULL DEFAULT 'PENDING',

    -- Concordance metrics
    source_count INT NOT NULL DEFAULT 0,
    concordant_count INT NOT NULL DEFAULT 0,
    discordant_count INT NOT NULL DEFAULT 0,
    confidence_score DOUBLE PRECISION NOT NULL DEFAULT 0.0,

    -- Best quality metrics across sources
    max_read_depth INT,
    max_quality_score DOUBLE PRECISION,

    -- Haplogroup association
    defining_haplogroup VARCHAR(100),
    haplogroup_branch_depth INT,

    -- Timestamps
    last_updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    -- Unique per profile+position+alleles
    CONSTRAINT uq_variant_position UNIQUE (y_profile_id, contig, position, ref_allele, alt_allele),
    CONSTRAINT chk_variant_type CHECK (variant_type IN ('SNP', 'INDEL', 'MNP', 'STR')),
    CONSTRAINT chk_variant_status CHECK (status IN ('CONFIRMED', 'NOVEL', 'CONFLICT', 'NO_COVERAGE', 'PENDING')),
    CONSTRAINT chk_consensus_state CHECK (consensus_state IN ('DERIVED', 'ANCESTRAL', 'HETEROPLASMY', 'NO_CALL'))
);

CREATE INDEX idx_y_variant_profile ON y_profile_variant(y_profile_id);
CREATE INDEX idx_y_variant_position ON y_profile_variant(position);
CREATE INDEX idx_y_variant_name ON y_profile_variant(variant_name);
CREATE INDEX idx_y_variant_marker ON y_profile_variant(marker_name);
CREATE INDEX idx_y_variant_status ON y_profile_variant(status);
CREATE INDEX idx_y_variant_type ON y_profile_variant(variant_type);
CREATE INDEX idx_y_variant_haplogroup ON y_profile_variant(defining_haplogroup);

-- ============================================
-- Y_VARIANT_SOURCE_CALL: Per-source variant calls
-- ============================================

CREATE TABLE y_variant_source_call (
    id UUID PRIMARY KEY,
    variant_id UUID NOT NULL REFERENCES y_profile_variant(id) ON DELETE CASCADE,
    source_id UUID NOT NULL REFERENCES y_profile_source(id) ON DELETE CASCADE,

    -- Call data
    called_allele VARCHAR(500) NOT NULL,
    call_state VARCHAR(30) NOT NULL,

    -- STR-specific
    called_repeat_count INT,

    -- Quality metrics from this source
    read_depth INT,
    quality_score DOUBLE PRECISION,
    mapping_quality DOUBLE PRECISION,
    variant_allele_frequency DOUBLE PRECISION,

    -- Callable state at this position
    callable_state VARCHAR(30),

    -- Weight in concordance calculation (computed)
    concordance_weight DOUBLE PRECISION NOT NULL DEFAULT 1.0,

    -- Timestamps
    called_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    -- Unique per variant+source
    CONSTRAINT uq_variant_source_call UNIQUE (variant_id, source_id),
    CONSTRAINT chk_call_state CHECK (call_state IN ('DERIVED', 'ANCESTRAL', 'NO_CALL', 'HETEROPLASMY'))
);

CREATE INDEX idx_variant_call_variant ON y_variant_source_call(variant_id);
CREATE INDEX idx_variant_call_source ON y_variant_source_call(source_id);
CREATE INDEX idx_variant_call_state ON y_variant_source_call(call_state);

-- ============================================
-- Y_VARIANT_AUDIT: Manual override audit trail
-- ============================================

CREATE TABLE y_variant_audit (
    id UUID PRIMARY KEY,
    variant_id UUID NOT NULL REFERENCES y_profile_variant(id) ON DELETE CASCADE,

    -- Change tracking
    action VARCHAR(30) NOT NULL,

    -- Previous values
    previous_consensus_allele VARCHAR(500),
    previous_consensus_state VARCHAR(30),
    previous_status VARCHAR(30),
    previous_confidence DOUBLE PRECISION,

    -- New values
    new_consensus_allele VARCHAR(500),
    new_consensus_state VARCHAR(30),
    new_status VARCHAR(30),
    new_confidence DOUBLE PRECISION,

    -- User and reason
    user_id VARCHAR(255),
    reason TEXT NOT NULL,
    supporting_evidence TEXT,

    -- Timestamps
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT chk_audit_action CHECK (action IN ('OVERRIDE', 'CONFIRM', 'REJECT', 'ANNOTATE', 'REVERT'))
);

CREATE INDEX idx_variant_audit_variant ON y_variant_audit(variant_id);
CREATE INDEX idx_variant_audit_action ON y_variant_audit(action);
CREATE INDEX idx_variant_audit_created ON y_variant_audit(created_at);
CREATE INDEX idx_variant_audit_user ON y_variant_audit(user_id);

-- ============================================
-- Add deprecation comment to haplogroup_reconciliation
-- ============================================

COMMENT ON TABLE haplogroup_reconciliation IS
    'DEPRECATED for Y-DNA: Use y_chromosome_profile instead. Retained for MT-DNA only.';
