-- Read-profile fields backing the standardized DTC test label (design doc
-- `documents/atmosphere/11-Standardized-Test-Profiles.md`). `total_bases` is the exact sequenced
-- yield (Σ read_length_histogram) — the "Gbases" figure; `read_type` is the read chemistry/mode
-- (SHORT/HIFI/CLR/ONT_*), the only signal that tells HiFi from CLR for the long-read label.
ALTER TABLE sequence_run ADD COLUMN total_bases INTEGER;
ALTER TABLE sequence_run ADD COLUMN read_type TEXT;
