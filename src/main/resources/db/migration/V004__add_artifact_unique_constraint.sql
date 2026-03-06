-- V004__add_artifact_unique_constraint.sql
-- Add unique constraint on alignment_id + artifact_type to prevent duplicate artifacts

ALTER TABLE analysis_artifact ADD CONSTRAINT uq_artifact_alignment_type UNIQUE (alignment_id, artifact_type);
