-- Absolute path of the imported chip raw-data file, so ancestry-from-chip can re-read the full
-- ~600k genotypes on demand (we only persist the QC summary) — mirrors how `alignment.bam_path`
-- lets analyses re-read the BAM. NULL for pre-existing rows / when the path is unknown.
ALTER TABLE chip_profile ADD COLUMN source_path TEXT;
