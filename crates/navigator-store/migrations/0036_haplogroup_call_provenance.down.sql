-- Reverse the external re-key (mt first so ':ext:mt' is undone before the ':ext' strip), then drop
-- the column.
UPDATE haplogroup_call
   SET source_key = replace(source_key, ':ext:mt', ':mt')
 WHERE provenance = 'external' AND source_key LIKE '%:ext:mt';

UPDATE haplogroup_call
   SET source_key = substr(source_key, 1, length(source_key) - 4)
 WHERE provenance = 'external' AND source_key LIKE '%:ext';

ALTER TABLE haplogroup_call DROP COLUMN provenance;
