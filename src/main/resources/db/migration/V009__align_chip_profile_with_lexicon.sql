-- V009: Align chip_profile with Atmosphere Lexicon genotype record
-- Renames vendor -> provider and adds new fields from the web app schema.

-- Rename vendor to provider (Atmosphere Lexicon field name)
ALTER TABLE chip_profile ALTER COLUMN vendor RENAME TO provider;

-- Drop and recreate vendor index with new column name
DROP INDEX IF EXISTS idx_chip_profile_vendor;
CREATE INDEX idx_chip_profile_provider ON chip_profile(provider);

-- Add new fields from com.decodingus.atmosphere.genotype schema
ALTER TABLE chip_profile ADD COLUMN y_markers_total INT;
ALTER TABLE chip_profile ADD COLUMN mt_markers_total INT;
ALTER TABLE chip_profile ADD COLUMN test_date TIMESTAMP;
ALTER TABLE chip_profile ADD COLUMN processed_at TIMESTAMP;
ALTER TABLE chip_profile ADD COLUMN build_version VARCHAR(20);
ALTER TABLE chip_profile ADD COLUMN derived_haplogroups JSON;
ALTER TABLE chip_profile ADD COLUMN population_breakdown_ref VARCHAR(500);
ALTER TABLE chip_profile ADD COLUMN imputation_ref VARCHAR(500);
