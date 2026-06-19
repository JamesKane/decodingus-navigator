-- Bind each analysis artifact to its source file's signature (mtime:size) at compute time, so a
-- changed BAM/CRAM (re-aligned at the same path) invalidates stale cached results. NULL for rows
-- written before this column existed → treated as fresh (back-compat; no forced recompute).
ALTER TABLE analysis_artifact ADD COLUMN source_sig TEXT;
