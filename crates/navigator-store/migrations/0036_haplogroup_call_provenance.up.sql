-- Provenance tier for each per-source Y/mtDNA haplogroup call.
--
-- Before this, the sidecar fast path (external GATK4 GVCF) and Navigator's internal CRAM walk
-- wrote the *same* source_key (`aln:{id}` / `aln:{id}:mt`), so the upsert let an internal
-- re-analysis silently overwrite an externally-imported call — the PRJEB37976 ancient-DNA
-- "haplogroups changed" bug. Provenance makes external calls a distinct, preferable tier; the
-- external rows are also re-keyed (`:ext`) so a walk can never collide with them again.
--
-- Legacy rows default to the internal 'navigator-walk' tier. Rows the fast path wrote carry a
-- 'gv:'-prefixed source_fingerprint and are the external ones.
ALTER TABLE haplogroup_call ADD COLUMN provenance TEXT NOT NULL DEFAULT 'navigator-walk';

UPDATE haplogroup_call SET provenance = 'external' WHERE source_fingerprint LIKE 'gv:%';

-- Re-key surviving external rows so future fast-path upserts target them (never a walk's key).
-- Order matters: mt first (its key ends ':mt'), then Y (everything else). ids are numeric, so
-- ':mt' cannot appear elsewhere in the key.
UPDATE haplogroup_call
   SET source_key = replace(source_key, ':mt', ':ext:mt')
 WHERE provenance = 'external' AND source_key LIKE '%:mt';

UPDATE haplogroup_call
   SET source_key = source_key || ':ext'
 WHERE provenance = 'external' AND source_key NOT LIKE '%:mt';
