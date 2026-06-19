-- Genotyping panels: named sets of biallelic SNP sites (ancestry-informative markers /
-- IBD sites) that samples are genotyped against for population + IBD analysis.
CREATE TABLE panel (
    id   INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL UNIQUE
);

CREATE TABLE panel_site (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    panel_id        INTEGER NOT NULL REFERENCES panel(id),
    chrom           TEXT NOT NULL,
    position        INTEGER NOT NULL,
    reference_allele TEXT NOT NULL,
    alternate_allele TEXT NOT NULL,
    name            TEXT NOT NULL
);

CREATE INDEX idx_panel_site_panel ON panel_site(panel_id);
