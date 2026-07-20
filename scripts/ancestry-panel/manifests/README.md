# Ancestry-panel manifests

## `consumer_array_1240k_rsids.txt.gz` — deep-ancestry ascertainment floor (Option A′)

One rsID per line (gzipped): the sites that mainstream consumer DNA arrays assay, restricted to the
AADR 1240k universe. Fed to `panelbuild ancient-panel --ascertain-sites` (via `$CHIP_MANIFEST` in
`config.sh`) so the deep-ancestry panel is built only on array-ascertained sites.

**Why this exists.** Allele-frequency admixture is only valid when the sample and the reference share
ascertainment. The AADR/1240k universe includes ancient-capture sites consumer arrays don't assay, and
on those the deep estimate is unstable across data sources — the *same* subject reads ~90% Steppe from
WGS but ~58% from his own chip. Intersecting the panel with the sites arrays actually assay removes
that split. See `docs/design/ancient-ancestry-rebuild.md` §4.

**Provenance.** The union of the probe rsIDs of two GSA-based consumer arrays — 23andMe v5 and
AncestryDNA v2 (the two most common consumer tests) — intersected with the AADR 1240k `.snp` (v66.p1).
This is the arrays' *design* (which sites they type), not anyone's genotypes: no personal data.
- array probe rsIDs (23andMe v5 ∪ AncestryDNA v2): 1,401,714
- AADR 1240k rs* universe: 1,156,490
- **manifest = array ∩ 1240k: 649,478**

**Validation.** The ancient panel built through this floor is 9,971 of 19,727 sites and clears every
§3.4 gate: stability (WGS-on-array-sites 57.9% Steppe vs chips 56–58%), recovery (20/30/50 →
20.2/29.2/50.6), sanity band (GBR 50.0 Steppe, all Europeans in band, PJL/CHB/JPT/YRI/LWK rejected),
and density (GBR 50.0 → 51.8 at half the sites).

**To refresh / broaden** (e.g. new AADR release, more arrays):
```sh
# array probe rsIDs (col 1 of any 23andMe/AncestryDNA/FTDNA/etc. raw export — array design, public)
awk 'BEGIN{FS="[ \t,]"} $0!~/^#/ && $1~/^rs/ {print $1}' chip1.txt chip2.txt ... | sort -u > array.txt
# AADR 1240k rsIDs
awk '$1~/^rs/ {print $1}' "$RAW/${AADR_FILE_PREFIX}.snp" | sort -u > k1240.txt
comm -12 array.txt k1240.txt | gzip > consumer_array_1240k_rsids.txt.gz
```
A canonical alternative is the AADR **Human Origins** panel (a published ~600k array ascertainment,
same Dataverse as 1240k); it was not used here only because it wasn't downloaded — validating it is
the natural next step if a broader, vendor-neutral ascertainment is wanted.
