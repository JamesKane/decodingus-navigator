-- V005: Phase 2 Entity Tables
-- STR profiles, Chip profiles, and Haplogroup reconciliation

-- ============================================
-- STR Profile Table
-- ============================================

CREATE TABLE str_profile (
    id UUID PRIMARY KEY,
    biosample_id UUID NOT NULL REFERENCES biosample(id) ON DELETE CASCADE,
    sequence_run_id UUID REFERENCES sequence_run(id) ON DELETE SET NULL,

    -- Panel information (JSON array)
    panels JSON DEFAULT '[]',

    -- Marker values (JSON array of marker objects)
    markers JSON DEFAULT '[]',

    -- Summary
    total_markers INT,
    source VARCHAR(50),                   -- DIRECT_TEST, WGS_DERIVED, BIG_Y_DERIVED, IMPORTED, MANUAL_ENTRY
    imported_from VARCHAR(100),           -- FTDNA, YSEQ, YFULL, etc.
    derivation_method VARCHAR(100),       -- HIPSTR, GANGSTR, LOBSTR, etc.

    -- Files (JSON array)
    files JSON DEFAULT '[]',

    -- Sync tracking
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(512) UNIQUE,
    at_cid VARCHAR(128),

    -- Versioning
    version INT NOT NULL DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_str_profile_biosample ON str_profile(biosample_id);
CREATE INDEX idx_str_profile_sequence_run ON str_profile(sequence_run_id);
CREATE INDEX idx_str_profile_sync ON str_profile(sync_status);

-- ============================================
-- Chip Profile Table
-- ============================================

CREATE TABLE chip_profile (
    id UUID PRIMARY KEY,
    biosample_id UUID NOT NULL REFERENCES biosample(id) ON DELETE CASCADE,

    -- Vendor and test info
    vendor VARCHAR(100) NOT NULL,         -- 23andMe, AncestryDNA, FTDNA, etc.
    test_type_code VARCHAR(50) NOT NULL,
    chip_version VARCHAR(50),

    -- Marker statistics
    total_markers_called INT NOT NULL,
    total_markers_possible INT NOT NULL,
    no_call_rate DOUBLE NOT NULL,
    y_markers_called INT,
    mt_markers_called INT,
    autosomal_markers_called INT NOT NULL,
    het_rate DOUBLE,

    -- Import info
    import_date TIMESTAMP NOT NULL,
    source_file_hash VARCHAR(128),
    source_file_name VARCHAR(512),

    -- Files (JSON array)
    files JSON DEFAULT '[]',

    -- Sync tracking
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(512) UNIQUE,
    at_cid VARCHAR(128),

    -- Versioning
    version INT NOT NULL DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_chip_profile_biosample ON chip_profile(biosample_id);
CREATE INDEX idx_chip_profile_vendor ON chip_profile(vendor);
CREATE INDEX idx_chip_profile_sync ON chip_profile(sync_status);

-- ============================================
-- Haplogroup Reconciliation Table
-- ============================================

CREATE TABLE haplogroup_reconciliation (
    id UUID PRIMARY KEY,
    biosample_id UUID NOT NULL REFERENCES biosample(id) ON DELETE CASCADE,

    -- DNA type (Y_DNA or MT_DNA)
    dna_type VARCHAR(10) NOT NULL,

    -- Reconciliation status (JSON object with consensus, compatibility level, etc.)
    status JSON NOT NULL,

    -- Individual run calls (JSON array)
    run_calls JSON DEFAULT '[]',

    -- SNP-level conflicts (JSON array)
    snp_conflicts JSON DEFAULT '[]',

    -- Last reconciliation timestamp
    last_reconciliation_at TIMESTAMP,

    -- Sync tracking
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(512) UNIQUE,
    at_cid VARCHAR(128),

    -- Versioning
    version INT NOT NULL DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    -- Ensure one reconciliation per biosample per DNA type
    CONSTRAINT uq_reconciliation_biosample_dna UNIQUE (biosample_id, dna_type),
    CONSTRAINT chk_dna_type CHECK (dna_type IN ('Y_DNA', 'MT_DNA'))
);

CREATE INDEX idx_haplogroup_reconciliation_biosample ON haplogroup_reconciliation(biosample_id);
CREATE INDEX idx_haplogroup_reconciliation_dna_type ON haplogroup_reconciliation(dna_type);
CREATE INDEX idx_haplogroup_reconciliation_sync ON haplogroup_reconciliation(sync_status);

-- ============================================
-- Y-SNP Panel Table (called SNPs from testing)
-- ============================================

CREATE TABLE y_snp_panel (
    id UUID PRIMARY KEY,
    biosample_id UUID NOT NULL REFERENCES biosample(id) ON DELETE CASCADE,
    alignment_id UUID REFERENCES alignment(id) ON DELETE SET NULL,

    -- Panel info
    panel_name VARCHAR(100),              -- Big Y-700, WGS, etc.
    provider VARCHAR(100),                -- FTDNA, YSEQ, etc.
    test_date TIMESTAMP,

    -- SNP counts
    total_snps_tested INT,
    derived_count INT,
    ancestral_count INT,
    no_call_count INT,

    -- Terminal haplogroup from this panel
    terminal_haplogroup VARCHAR(100),
    confidence DOUBLE,

    -- SNP calls (JSON array of {name, position, allele, derived, quality})
    snp_calls JSON DEFAULT '[]',

    -- Novel/private variants (JSON array)
    private_variants JSON DEFAULT '[]',

    -- Files (JSON array)
    files JSON DEFAULT '[]',

    -- Sync tracking
    sync_status VARCHAR(20) NOT NULL DEFAULT 'Local',
    at_uri VARCHAR(512) UNIQUE,
    at_cid VARCHAR(128),

    -- Versioning
    version INT NOT NULL DEFAULT 1,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_y_snp_panel_biosample ON y_snp_panel(biosample_id);
CREATE INDEX idx_y_snp_panel_alignment ON y_snp_panel(alignment_id);
CREATE INDEX idx_y_snp_panel_haplogroup ON y_snp_panel(terminal_haplogroup);
CREATE INDEX idx_y_snp_panel_sync ON y_snp_panel(sync_status);
