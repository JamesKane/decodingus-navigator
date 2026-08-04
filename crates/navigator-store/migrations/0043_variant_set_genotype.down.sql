DROP INDEX IF EXISTS ix_variant_set_genotype_set;
DROP TABLE IF EXISTS variant_set_genotype;
ALTER TABLE variant_set DROP COLUMN source_path;
