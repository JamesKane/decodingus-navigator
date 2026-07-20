-- Purge the pre-GRCh37/38-support IBD-panel genotype cache.
--
-- Before build-aware genotyping landed (app commit 3cf4956, 2026-07-13), a GRCh37/GRCh38 alignment
-- was genotyped at the panel's CHM13 coordinates on its non-CHM13 BAM: wrong positions, so the
-- dosages were near-random and never heterozygous. The code fix did not change the cache key
-- (panel-manifest salt + genotype version), so those corrupt genotypes kept feeding the autosomal
-- consensus (IBD / ancestry / identity). The app now writes this cache under the kind stem
-- `ibd_panel_genotypes2`; drop every row under the old stem so the stale calls can't be read and the
-- (often hundreds of MB per alignment) dead payloads are reclaimed. Re-genotyping is automatic and
-- build-aware on next use.
DELETE FROM analysis_artifact
WHERE kind = 'ibd_panel_genotypes'
   OR kind LIKE 'ibd_panel_genotypes:%';
