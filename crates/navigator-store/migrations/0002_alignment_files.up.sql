-- An alignment knows its BAM/CRAM and the reference it was aligned to, so analysis
-- can be run directly from a stored alignment.
ALTER TABLE alignment ADD COLUMN bam_path TEXT;
ALTER TABLE alignment ADD COLUMN reference_path TEXT;
