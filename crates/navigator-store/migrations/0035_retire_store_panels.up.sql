-- Retire the legacy, build-naive store-panel genotyping/IBD/identity path. The `panel`/`panel_site`
-- tables carried no reference-build column, so genotyping an alignment's BAM directly at a panel's
-- raw (contig, position) silently produced wrong genotypes whenever the panel and the alignment were
-- on different builds. The build-aware IBD-panel (asset-backed) and consensus paths supersede it.
-- Drop the child table first to respect the panel_site → panel foreign key.
DROP TABLE IF EXISTS panel_site;
DROP TABLE IF EXISTS panel;
