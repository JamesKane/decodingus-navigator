-- V008__add_sequence_run_metrics_columns.sql
-- Add missing metrics columns to sequence_run table
-- These fields are populated during initial analysis but were not being persisted

-- ============================================
-- CRITICAL FIX: Update platform constraint to match normalized values
-- The original constraint was case-sensitive and too restrictive
-- ============================================
ALTER TABLE sequence_run DROP CONSTRAINT IF EXISTS chk_sr_platform;
ALTER TABLE sequence_run ADD CONSTRAINT chk_sr_platform
    CHECK (platform IN ('ILLUMINA', 'PACBIO', 'NANOPORE', 'ION_TORRENT', 'BGI', 'ELEMENT', 'ULTIMA', 'Unknown', 'Other'));

-- ============================================
-- CRITICAL FIX: Update library_layout constraint to match normalized values
-- Analysis was setting "Paired-End"/"Single-End" but constraint required "PAIRED"/"SINGLE"
-- ============================================
ALTER TABLE sequence_run DROP CONSTRAINT IF EXISTS chk_sr_layout;
ALTER TABLE sequence_run ADD CONSTRAINT chk_sr_layout
    CHECK (library_layout IS NULL OR library_layout IN ('PAIRED', 'SINGLE', 'Paired-End', 'Single-End'));

-- ============================================
-- Add missing alignment percentage metrics
-- ============================================
ALTER TABLE sequence_run ADD COLUMN IF NOT EXISTS pct_pf_reads_aligned DOUBLE PRECISION;
ALTER TABLE sequence_run ADD COLUMN IF NOT EXISTS reads_paired BIGINT;
ALTER TABLE sequence_run ADD COLUMN IF NOT EXISTS pct_reads_paired DOUBLE PRECISION;
ALTER TABLE sequence_run ADD COLUMN IF NOT EXISTS pct_proper_pairs DOUBLE PRECISION;

-- ============================================
-- Add missing read length metrics
-- ============================================
ALTER TABLE sequence_run ADD COLUMN IF NOT EXISTS max_read_length INT;

-- ============================================
-- Add pair orientation (FR, RF, TANDEM)
-- ============================================
ALTER TABLE sequence_run ADD COLUMN IF NOT EXISTS pair_orientation VARCHAR(20);

-- Add constraint for pair_orientation values
ALTER TABLE sequence_run ADD CONSTRAINT chk_sr_pair_orientation
    CHECK (pair_orientation IS NULL OR pair_orientation IN ('FR', 'RF', 'TANDEM'));
