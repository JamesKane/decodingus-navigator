# CHM13 / T2T Reference Resources

Scraped from [marbl/CHM13](https://github.com/marbl/CHM13) (2026-06-03) for the Rust
rewrite. These resources did **not** exist when the Scala version was written.

All assemblies below are **T2T-CHM13v2.0** (the v2.0 release that includes the Y
chromosome). Everything on the `human-pangenomics` S3 bucket is public domain (CC0)
and fetchable anonymously — no AWS credentials needed.

## Access note for Rust code

Everything under `s3-us-west-2.amazonaws.com/human-pangenomics/...` is also reachable
as `s3://human-pangenomics/T2T/CHM13/...` (anonymous). For the `ReferenceGateway`
equivalent in Rust, just HTTP-GET the HTTPS URLs directly.

---

## Liftover chain files

UCSC liftOver chain format, but these are **1:1 chains** (single alignment block set),
not the multi-level chains UCSC normally ships. Pair them with the `unique_to_*` BEDs
below to detect coordinates that won't lift.

### GRCh38 ↔ CHM13v2.0
- `grch38-chm13v2.chain` — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain
- `chm13v2-grch38.chain` — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain
- PAF alignment — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.paf

### hg19/GRCh37 ↔ CHM13v2.0
- `hg19-chm13v2.chain` — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/hg19-chm13v2.chain
- `chm13v2-hg19.chain` — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-hg19.chain
- PAF alignment — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/hg19-chm13v2.paf

### Non-syntenic / unique regions (no clean 1:1 mapping — flag as unliftable)
- Unique to GRCh38 — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed
- Unique to GRCh37 — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg19.bed

## Reference genome FASTA

- Full assembly (with Y) — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
- No Y — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_noY.fa.gz
- **Y PAR-masked** — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY.fa.gz
- **Y PAR-masked + rCRS mito** — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz

> Each FASTA above has a companion `.fa.gz.gzi` bgzip index alongside it. The bucket
> also carries pre-built BWA indexes under `analysis_set/masked_DJ_rDNA_PHR*/` if needed.

README guidance: for short-read variant calling use a **Y-PAR-masked** FASTA. The
`.rCRS.fa.gz` variant carries the rCRS mitochondrion and is likely the most relevant
for this app (matches the haplogroup/mtDNA work).

## Variant catalogs / reference panels

- **dbSNP build 155** (lifted) — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_dbSNPv155.vcf.gz
- **ClinVar 20220313** (lifted) — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_ClinVar20220313.vcf.gz
- **GWAS catalog v1.0** (lifted, with rsids) — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_GWASv1.0rsids_e100_r2022-03-08.vcf.gz
- **gnomAD v3.1.2** (Ensembl, lifted) — https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_10/vcf/Homo_sapiens-GCA_009914755.4-2022_10-gnomad.vcf.gz
- **1000 Genomes** on CHM13 (dir listing) — https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/
  - Per-population allele freqs — `.../unrelated_samples_2504/allele_freq/`
  - Phased (SHAPEIT5) — https://github.com/JosephLalli/phasing_T2T (imputation-panel material)
- **Simons Genome Diversity Project** on CHM13 — https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/variants/SGDP/chm13v2.0/

## Annotation / masks (callable-loci & QC)

- RefSeqv110 + Liftoff gene annotation — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz
- RepeatMasker — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed
- Centromere/satellite (CenSat v2.1) — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed
- Segmental duplications — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed
- Telomere coords — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_telomere.bed
- **Short-read accessibility mask** — https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/accessibility/combined_mask.bed.gz

## S3 browse

- Web interface — https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13
- CLI — `s3://human-pangenomics/T2T/CHM13/` (anonymous, e.g. `aws s3 ls --no-sign-request`)

---

*Caveat: URLs were pulled via a page summary. Verify exact filenames/paths before
wiring into code — the README occasionally reorganizes paths.*
