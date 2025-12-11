-- V001__initial_schema.sql
-- Clean schema for Phase 1 core entities
-- H2 database with PostgreSQL compatibility mode

-- ============================================
-- BIOSAMPLE: Primary research subject
-- ============================================
CREATE TABLE biosample (
    id UUID PRIMARY KEY,
    sample_accession VARCHAR(255) NOT NULL,
    donor_identifier VARCHAR(255) NOT NULL,
    description TEXT,
    center_name VARCHAR(255),
    sex VARCHAR(20),
    citizen_did VARCHAR(512),

    -- Haplogroup assignments as JSON (complex nested object)
    haplogroups JSON,

    -- Sync tracking for PDS integration
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(512),
    at_cid VARCHAR(128),

    -- Record versioning
    version INT NOT NULL DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT uq_biosample_accession UNIQUE (sample_accession),
    CONSTRAINT uq_biosample_at_uri UNIQUE (at_uri),
    CONSTRAINT chk_biosample_sex CHECK (sex IS NULL OR sex IN ('Male', 'Female', 'Other', 'Unknown')),
    CONSTRAINT chk_biosample_sync CHECK (sync_status IN ('Local', 'Synced', 'Modified', 'Conflict'))
);

CREATE INDEX idx_biosample_accession ON biosample(sample_accession);
CREATE INDEX idx_biosample_donor ON biosample(donor_identifier);
CREATE INDEX idx_biosample_sync ON biosample(sync_status);
CREATE INDEX idx_biosample_citizen ON biosample(citizen_did);

-- ============================================
-- PROJECT: Groups biosamples for research
-- ============================================
CREATE TABLE project (
    id UUID PRIMARY KEY,
    project_name VARCHAR(255) NOT NULL,
    description TEXT,
    administrator_did VARCHAR(512) NOT NULL,

    -- Sync tracking
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(512),
    at_cid VARCHAR(128),

    -- Record versioning
    version INT NOT NULL DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT uq_project_name UNIQUE (project_name),
    CONSTRAINT uq_project_at_uri UNIQUE (at_uri),
    CONSTRAINT chk_project_sync CHECK (sync_status IN ('Local', 'Synced', 'Modified', 'Conflict'))
);

CREATE INDEX idx_project_name ON project(project_name);
CREATE INDEX idx_project_admin ON project(administrator_did);
CREATE INDEX idx_project_sync ON project(sync_status);

-- ============================================
-- PROJECT_MEMBER: Junction table for project membership
-- Proper relational design instead of JSON arrays
-- ============================================
CREATE TABLE project_member (
    project_id UUID NOT NULL,
    biosample_id UUID NOT NULL,
    added_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    PRIMARY KEY (project_id, biosample_id),
    CONSTRAINT fk_pm_project FOREIGN KEY (project_id) REFERENCES project(id) ON DELETE CASCADE,
    CONSTRAINT fk_pm_biosample FOREIGN KEY (biosample_id) REFERENCES biosample(id) ON DELETE CASCADE
);

CREATE INDEX idx_pm_project ON project_member(project_id);
CREATE INDEX idx_pm_biosample ON project_member(biosample_id);

-- ============================================
-- SEQUENCE_RUN: Single sequencing session
-- ============================================
CREATE TABLE sequence_run (
    id UUID PRIMARY KEY,
    biosample_id UUID NOT NULL,

    -- Platform information
    platform VARCHAR(50) NOT NULL,
    instrument_model VARCHAR(255),
    instrument_id VARCHAR(255),
    test_type VARCHAR(50) NOT NULL,

    -- Library information
    library_id VARCHAR(255),
    platform_unit VARCHAR(255),
    library_layout VARCHAR(20),
    sample_name VARCHAR(255),
    sequencing_facility VARCHAR(255),
    run_fingerprint VARCHAR(128),

    -- Metrics (nullable, populated after analysis)
    total_reads BIGINT,
    pf_reads BIGINT,
    pf_reads_aligned BIGINT,
    read_length INT,
    mean_insert_size DOUBLE PRECISION,
    median_insert_size DOUBLE PRECISION,
    std_insert_size DOUBLE PRECISION,

    -- Run metadata
    flowcell_id VARCHAR(255),
    run_date TIMESTAMP,

    -- Files as JSON array (FileInfo objects)
    files JSON DEFAULT '[]',

    -- Sync tracking
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(512),
    at_cid VARCHAR(128),

    -- Record versioning
    version INT NOT NULL DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT fk_sr_biosample FOREIGN KEY (biosample_id) REFERENCES biosample(id) ON DELETE CASCADE,
    CONSTRAINT uq_sequence_run_at_uri UNIQUE (at_uri),
    CONSTRAINT chk_sr_platform CHECK (platform IN ('ILLUMINA', 'PACBIO', 'NANOPORE', 'ION_TORRENT', 'BGI', 'ELEMENT', 'ULTIMA', 'Unknown')),
    CONSTRAINT chk_sr_layout CHECK (library_layout IS NULL OR library_layout IN ('PAIRED', 'SINGLE')),
    CONSTRAINT chk_sr_sync CHECK (sync_status IN ('Local', 'Synced', 'Modified', 'Conflict'))
);

CREATE INDEX idx_sr_biosample ON sequence_run(biosample_id);
CREATE INDEX idx_sr_platform ON sequence_run(platform);
CREATE INDEX idx_sr_test_type ON sequence_run(test_type);
CREATE INDEX idx_sr_library_id ON sequence_run(library_id);
CREATE INDEX idx_sr_platform_unit ON sequence_run(platform_unit);
CREATE INDEX idx_sr_sync ON sequence_run(sync_status);

-- ============================================
-- ALIGNMENT: Mapping to reference genome
-- ============================================
CREATE TABLE alignment (
    id UUID PRIMARY KEY,
    sequence_run_id UUID NOT NULL,

    -- Reference information
    reference_build VARCHAR(50) NOT NULL,
    aligner VARCHAR(255) NOT NULL,
    variant_caller VARCHAR(255),

    -- Metrics as JSON (AlignmentMetrics - complex nested structure)
    metrics JSON,

    -- Files as JSON array (FileInfo objects)
    files JSON DEFAULT '[]',

    -- Sync tracking
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(512),
    at_cid VARCHAR(128),

    -- Record versioning
    version INT NOT NULL DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT fk_align_seqrun FOREIGN KEY (sequence_run_id) REFERENCES sequence_run(id) ON DELETE CASCADE,
    CONSTRAINT uq_alignment_at_uri UNIQUE (at_uri),
    CONSTRAINT chk_align_ref CHECK (reference_build IN ('GRCh38', 'GRCh37', 'T2T-CHM13', 'hg19', 'hg38')),
    CONSTRAINT chk_align_sync CHECK (sync_status IN ('Local', 'Synced', 'Modified', 'Conflict'))
);

CREATE INDEX idx_align_seqrun ON alignment(sequence_run_id);
CREATE INDEX idx_align_ref ON alignment(reference_build);
CREATE INDEX idx_align_sync ON alignment(sync_status);
