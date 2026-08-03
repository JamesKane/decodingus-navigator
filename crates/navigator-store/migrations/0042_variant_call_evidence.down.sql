-- SQLite has supported ALTER TABLE ... DROP COLUMN since 3.35; these are plain columns with no
-- index or generated-column dependency, so dropping them is safe. Re-applying the up migration
-- restores the columns but not the evidence, which only a re-import can supply.
ALTER TABLE variant_set  DROP COLUMN call_schema;
ALTER TABLE variant_call DROP COLUMN ad_alt;
ALTER TABLE variant_call DROP COLUMN ad_ref;
ALTER TABLE variant_call DROP COLUMN gq;
ALTER TABLE variant_call DROP COLUMN dp;
ALTER TABLE variant_call DROP COLUMN filter;
ALTER TABLE variant_call DROP COLUMN qual;
