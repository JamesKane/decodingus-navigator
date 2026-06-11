-- Fingerprint of the inputs that produced this call (file content hash + tree content hash).
-- Lets scoring skip re-running when neither the alignment file nor the tree has changed.
ALTER TABLE haplogroup_call ADD COLUMN source_fingerprint TEXT;
