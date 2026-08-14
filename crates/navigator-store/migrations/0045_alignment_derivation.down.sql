DROP INDEX IF EXISTS ix_alignment_derived_from;
ALTER TABLE alignment DROP COLUMN derivation;
ALTER TABLE alignment DROP COLUMN derived_from_alignment_id;
