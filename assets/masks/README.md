# chrY private-Y filtering masks (CHM13v2 / hs1)

Bundled reference data for the private-Y variant filter (see
`docs/design/private-y-variant-filtering.md`). Small enough to live in the repo (gzipped BEDs, ~1.8 MB
total), unlike the ancestry `.bin` panels. They are staged into the installer by
`packaging/stage-assets.sh` (resource target `masks`) and seeded to `~/.decodingus/masks/` on first run
by `navigator_app::seed_bundled_masks()`. `RegionMask::from_bed` reads the `.gz` transparently.

Both are in **CHM13v2.0 (hs1) chrY** coordinates (0-based, half-open BED) and are applied **only** to
CHM13 alignments; other builds fall back to DecodingUs-tree classification + self-callable + region-class.

| File | Meaning | Filter |
|------|---------|--------|
| `chrY_callable_mask.chm13v2.bed.gz` | Poznik-style callable mask — positions CALLABLE (depth ≥4, MQ ≥20) in ≥90 % of a ~3,000-male cohort (14.96 Mbp = ~25 % of non-PAR chrY) | L2 (keep only these) |
| `chrY_cohort_shared_sites.chm13v2.bed.gz` | Every cohort position that varies with ≥2 carriers (joint-VCF `AC≥2`) ∪ homoplasy hotspots (recur on ≥3 tree branches). A true private has exactly one carrier (`AC=1`) and survives; a shared variant absent from the DecodingUs tree is a suspect artifact | L3 (drop these) |

## Regeneration (offline, from the ytree de-novo pipeline)

Source: `/Users/jkane/Genomics/ytree` (the CHM13 chrY/mtDNA de-novo tree workflow).

```bash
# L2 — callable mask (already produced by stage 3 of the pipeline)
gzip -c ytree/results/chrY.callable_mask.chm13v2.bed > chrY_callable_mask.chm13v2.bed.gz

# L3 — cohort-shared sites = joint-VCF AC>=2 positions  ∪  homoplasy hotspots
bcftools query -f '%CHROM\t%POS\t%INFO/AC\n' ytree/results/chrY.joint.vcf.gz \
  | awk 'BEGIN{OFS="\t"}{n=split($3,a,","); m=0; for(i=1;i<=n;i++) if(a[i]+0>m) m=a[i]+0; if(m>=2) print $1,$2-1,$2}' \
  > shared.bed
awk 'BEGIN{OFS="\t"} $1=="chrY"{print $1,$2,$3}' ytree/refs/branch_recurrent_exclude.chm13v2.bed >> shared.bed
sort -k1,1 -k2,2n -u shared.bed | gzip -c > chrY_cohort_shared_sites.chm13v2.bed.gz
```

Refreshing the cohort means re-running the pipeline's joint-genotyping and re-deriving both files. A
future step should add these to the sha256 `AssetManifest` alongside the ancestry assets.
