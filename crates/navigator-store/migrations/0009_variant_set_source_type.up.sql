-- Source type per variant set (concordance weighting; Sanger = gold standard).
ALTER TABLE variant_set ADD COLUMN source_type TEXT NOT NULL DEFAULT 'IMPORTED';
