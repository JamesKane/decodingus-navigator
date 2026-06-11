-- Content SHA-256 (hex) of the alignment file, computed at import. Identity of the file's
-- content; lets cached analyses (haplogroup scoring) invalidate only when the file changes.
-- NULL until computed (batch-imported files are hashed lazily on first analysis).
ALTER TABLE alignment ADD COLUMN content_sha256 TEXT;
