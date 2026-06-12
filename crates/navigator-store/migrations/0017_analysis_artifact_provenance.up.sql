-- Provenance for the analysis-artifact cache so the manual deep pass can be additive:
-- distinguish fast-path sidecar results from full CRAM walks, and partial (upgradeable)
-- results from complete ones. Both nullable → existing rows read as the defaults
-- (navigator-walk / full) the app applies on read.
ALTER TABLE analysis_artifact ADD COLUMN source TEXT;        -- 'navigator-walk' | 'pipeline-sidecar'
ALTER TABLE analysis_artifact ADD COLUMN completeness TEXT;  -- 'full' | 'partial'
