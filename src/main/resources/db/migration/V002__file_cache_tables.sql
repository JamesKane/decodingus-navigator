-- V002__file_cache_tables.sql
-- File cache and analysis artifact tracking tables

-- ============================================
-- SOURCE_FILE: Tracks user's BAM/CRAM files
-- ============================================
CREATE TABLE source_file (
    id UUID PRIMARY KEY,
    alignment_id UUID REFERENCES alignment(id) ON DELETE SET NULL,

    -- File identity (checksum is stable; path may change)
    file_path VARCHAR(2000),
    file_checksum VARCHAR(128) NOT NULL,
    file_size BIGINT,
    file_format VARCHAR(20),

    -- Last known state
    last_verified_at TIMESTAMP,
    is_accessible BOOLEAN DEFAULT TRUE,

    -- Analysis state
    has_been_analyzed BOOLEAN DEFAULT FALSE,
    analysis_completed_at TIMESTAMP,

    -- Metadata
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT uq_source_file_checksum UNIQUE (file_checksum),
    CONSTRAINT chk_sf_format CHECK (file_format IS NULL OR file_format IN ('BAM', 'CRAM', 'FASTQ', 'VCF', 'GVCF'))
);

CREATE INDEX idx_source_file_checksum ON source_file(file_checksum);
CREATE INDEX idx_source_file_alignment ON source_file(alignment_id);
CREATE INDEX idx_source_file_accessible ON source_file(is_accessible);

-- ============================================
-- ANALYSIS_ARTIFACT: Tracks cached analysis outputs
-- ============================================
CREATE TABLE analysis_artifact (
    id UUID PRIMARY KEY,
    alignment_id UUID NOT NULL REFERENCES alignment(id) ON DELETE CASCADE,

    -- Artifact type
    artifact_type VARCHAR(50) NOT NULL,

    -- Cache location (relative to ~/.decodingus/cache/)
    cache_path VARCHAR(1000) NOT NULL,

    -- File metadata
    file_size BIGINT,
    file_checksum VARCHAR(128),
    file_format VARCHAR(20),

    -- Generation tracking
    generated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    generator_version VARCHAR(100),
    generation_params JSON,

    -- Status
    status VARCHAR(20) NOT NULL DEFAULT 'AVAILABLE',
    stale_reason VARCHAR(255),

    -- Dependencies (for invalidation)
    depends_on_source_checksum VARCHAR(128),
    depends_on_reference_build VARCHAR(50),

    -- Metadata
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT chk_aa_type CHECK (artifact_type IN (
        'WGS_METRICS', 'CALLABLE_LOCI', 'HAPLOGROUP_VCF', 'WHOLE_GENOME_VCF',
        'PRIVATE_VARIANTS', 'INSERT_SIZE_METRICS', 'ALIGNMENT_SUMMARY',
        'COVERAGE_SUMMARY', 'SEX_INFERENCE', 'DUPLICATE_METRICS'
    )),
    CONSTRAINT chk_aa_status CHECK (status IN ('AVAILABLE', 'IN_PROGRESS', 'STALE', 'DELETED', 'ERROR'))
);

CREATE INDEX idx_artifact_alignment ON analysis_artifact(alignment_id);
CREATE INDEX idx_artifact_type ON analysis_artifact(artifact_type);
CREATE INDEX idx_artifact_status ON analysis_artifact(status);

-- ============================================
-- VCF_CACHE: VCF-specific metadata
-- ============================================
CREATE TABLE vcf_cache (
    id UUID PRIMARY KEY,
    artifact_id UUID NOT NULL REFERENCES analysis_artifact(id) ON DELETE CASCADE,

    -- VCF-specific metadata
    reference_build VARCHAR(50) NOT NULL,
    variant_caller VARCHAR(100),
    variant_caller_version VARCHAR(50),
    gatk_version VARCHAR(50),

    -- Statistics
    total_variants INT,
    snp_count INT,
    indel_count INT,
    ti_tv_ratio DOUBLE PRECISION,

    -- Contig statistics (JSON for flexibility)
    contig_stats JSON,

    -- Inferred metadata
    inferred_sex VARCHAR(20),
    sex_confidence DOUBLE PRECISION,

    -- Index file tracking
    index_path VARCHAR(1000),
    index_checksum VARCHAR(128),

    CONSTRAINT uq_vcf_artifact UNIQUE (artifact_id),
    CONSTRAINT chk_vcf_ref CHECK (reference_build IN ('GRCh38', 'GRCh37', 'T2T-CHM13', 'hg19', 'hg38'))
);

-- ============================================
-- CALLABLE_LOCI_CACHE: Callable loci summary
-- ============================================
CREATE TABLE callable_loci_cache (
    id UUID PRIMARY KEY,
    artifact_id UUID NOT NULL REFERENCES analysis_artifact(id) ON DELETE CASCADE,

    -- Summary statistics
    total_callable_bases BIGINT,
    total_ref_n_bases BIGINT,
    total_no_coverage_bases BIGINT,
    total_low_coverage_bases BIGINT,
    total_excessive_coverage_bases BIGINT,
    total_poor_mapping_bases BIGINT,

    -- Completion status
    is_complete BOOLEAN DEFAULT FALSE,
    completion_pct DOUBLE PRECISION,

    -- Output files
    bed_file_path VARCHAR(1000),
    summary_table_path VARCHAR(1000),
    visualization_path VARCHAR(1000),

    CONSTRAINT uq_callable_artifact UNIQUE (artifact_id)
);

-- ============================================
-- CALLABLE_LOCI_CONTIG: Per-contig breakdown
-- ============================================
CREATE TABLE callable_loci_contig (
    id UUID PRIMARY KEY,
    callable_loci_cache_id UUID NOT NULL REFERENCES callable_loci_cache(id) ON DELETE CASCADE,

    contig_name VARCHAR(50) NOT NULL,
    callable_bases BIGINT NOT NULL DEFAULT 0,
    ref_n_bases BIGINT NOT NULL DEFAULT 0,
    no_coverage_bases BIGINT NOT NULL DEFAULT 0,
    low_coverage_bases BIGINT NOT NULL DEFAULT 0,
    excessive_coverage_bases BIGINT NOT NULL DEFAULT 0,
    poor_mapping_bases BIGINT NOT NULL DEFAULT 0,
    mean_coverage DOUBLE PRECISION,
    callable_pct DOUBLE PRECISION,

    CONSTRAINT uq_contig_per_cache UNIQUE (callable_loci_cache_id, contig_name)
);

CREATE INDEX idx_callable_contig_cache ON callable_loci_contig(callable_loci_cache_id);

-- ============================================
-- COVERAGE_SUMMARY_CACHE: WGS metrics summary
-- ============================================
CREATE TABLE coverage_summary_cache (
    id UUID PRIMARY KEY,
    artifact_id UUID NOT NULL REFERENCES analysis_artifact(id) ON DELETE CASCADE,

    -- WGS Metrics
    genome_territory BIGINT,
    mean_coverage DOUBLE PRECISION,
    median_coverage DOUBLE PRECISION,
    sd_coverage DOUBLE PRECISION,
    mad_coverage DOUBLE PRECISION,
    pct_exc_adapter DOUBLE PRECISION,
    pct_exc_mapq DOUBLE PRECISION,
    pct_exc_dupe DOUBLE PRECISION,
    pct_exc_unpaired DOUBLE PRECISION,
    pct_exc_baseq DOUBLE PRECISION,
    pct_exc_overlap DOUBLE PRECISION,
    pct_exc_capped DOUBLE PRECISION,
    pct_exc_total DOUBLE PRECISION,
    pct_1x DOUBLE PRECISION,
    pct_5x DOUBLE PRECISION,
    pct_10x DOUBLE PRECISION,
    pct_15x DOUBLE PRECISION,
    pct_20x DOUBLE PRECISION,
    pct_25x DOUBLE PRECISION,
    pct_30x DOUBLE PRECISION,
    pct_40x DOUBLE PRECISION,
    pct_50x DOUBLE PRECISION,
    pct_60x DOUBLE PRECISION,
    pct_70x DOUBLE PRECISION,
    pct_80x DOUBLE PRECISION,
    pct_90x DOUBLE PRECISION,
    pct_100x DOUBLE PRECISION,
    het_snp_sensitivity DOUBLE PRECISION,
    het_snp_q DOUBLE PRECISION,

    CONSTRAINT uq_coverage_artifact UNIQUE (artifact_id)
);
