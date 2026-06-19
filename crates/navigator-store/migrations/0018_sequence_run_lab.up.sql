-- The lab/instrument identity block lost in the Rust rewrite (Scala SequenceRun parity):
-- inferred from the alignment at import (read-name scan + @RG tags). `sequencing_facility` is
-- the lab; `instrument_id` is the crowd-source key resolved against the AppView (roadmap D8).
ALTER TABLE sequence_run ADD COLUMN sequencing_facility TEXT;
ALTER TABLE sequence_run ADD COLUMN instrument_id TEXT;
ALTER TABLE sequence_run ADD COLUMN sample_name TEXT;
ALTER TABLE sequence_run ADD COLUMN library_id TEXT;
ALTER TABLE sequence_run ADD COLUMN platform_unit TEXT;
ALTER TABLE sequence_run ADD COLUMN flowcell_id TEXT;
