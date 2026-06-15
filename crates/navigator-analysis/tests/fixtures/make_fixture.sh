#!/usr/bin/env bash
# Regenerates the synthetic walker fixtures. Requires samtools on PATH.
#
# Produces three small, hand-computable BAMs (+ indexes) and a reference:
#   coverage.bam  reference chrM (50 bp, N at 25) — coverage + callable + caller
#   paired.bam    two FR proper pairs on chrM     — read_metrics (insert size/orient)
#   sex.bam       chr1 + chrX with biased counts  — sex inference (BAI metadata)
#
# See the per-walker tests for the expected values derived from these layouts.
set -euo pipefail
cd "$(dirname "$0")"

seq10="AAAAAAAAAA"; qual10="IIIIIIIIII"
seq5="AAAAA";       qual5="IIIII"

# ---- coverage.bam + ref.fa ------------------------------------------------
# 50 bp reference: 24 bp + N + 25 bp.
#   pos 1-10   depth 4  MAPQ 60  -> CALLABLE
#   pos 11-20  depth 2  MAPQ 60  -> LOW_COVERAGE
#   pos 21-24  depth 0           -> NO_COVERAGE
#   pos 25     ref N             -> REF_N
#   pos 26-30  depth 5  MAPQ 0   -> POOR_MAPPING_QUALITY
#   pos 31-50  depth 0           -> NO_COVERAGE
{
  echo ">chrM"
  echo "ACGTACGTACGTACGTACGTACGTNACGTACGTACGTACGTACGTACGTA"
} > ref.fa
samtools faidx ref.fa
[ "$(cut -f2 ref.fa.fai)" = "50" ] || { echo "ref length != 50" >&2; exit 1; }
{
  echo -e "@HD\tVN:1.6\tSO:coordinate"
  echo -e "@SQ\tSN:chrM\tLN:50"
  emit() { echo -e "$1\t0\tchrM\t$2\t$3\t$4\t*\t0\t0\t$5\t$6"; }
  for i in 1 2 3 4;   do emit "d4_$i"  1  60 10M "$seq10" "$qual10"; done
  for i in 1 2;       do emit "d2_$i"  11 60 10M "$seq10" "$qual10"; done
  for i in 1 2 3 4 5; do emit "mq0_$i" 26 0  5M  "$seq5"  "$qual5";  done
} > coverage.sam
samtools sort -o coverage.bam coverage.sam && samtools index coverage.bam
# CRAM counterpart (same reads, ref-compressed) — proves the reader unification is parity-clean.
samtools view -T ref.fa -C -o coverage.cram coverage.bam && samtools index coverage.cram

# ---- ychr.bam -------------------------------------------------------------
# A chrY analogue of coverage.bam (chrY, 50 bp; 4 'A' reads at pos 1 → 'A' callable at 1-10).
# Used by the haplogroup-liftover test: a seeded GRCh38→chm13 chain maps GRCh38 tree
# positions onto these chrY positions, exercising lift → query → map-back → score.
{
  echo -e "@HD\tVN:1.6\tSO:coordinate"
  echo -e "@SQ\tSN:chrY\tLN:50"
  for i in 1 2 3 4; do echo -e "y$i\t0\tchrY\t1\t60\t10M\t*\t0\t0\t$seq10\t$qual10"; done
} > ychr.sam
samtools sort -o ychr.bam ychr.sam && samtools index ychr.bam

# ---- paired.bam -----------------------------------------------------------
# Two FR proper pairs on chrM. flags: /1 = 99 (paired+proper+mate_rev+first),
# /2 = 147 (paired+proper+rev+last). Insert sizes 40 and 30 (first-of-pair TLEN).
{
  echo -e "@HD\tVN:1.6\tSO:coordinate"
  echo -e "@SQ\tSN:chrM\tLN:50"
  pe() { # name flag pos pnext tlen
    echo -e "$1\t$2\tchrM\t$3\t60\t10M\t=\t$4\t$5\t$seq10\t$qual10"
  }
  pe "pairA" 99  1  31  40
  pe "pairB" 99  5  25  30
  pe "pairB" 147 25 5  -30
  pe "pairA" 147 31 1  -40
} > paired.sam
samtools sort -o paired.bam paired.sam && samtools index paired.bam
# CRAM counterpart for the read-metrics CRAM parity test.
samtools view -T ref.fa -C -o paired.cram paired.bam && samtools index paired.cram

# ---- sex.bam --------------------------------------------------------------
# chr1 (100 bp) gets 10 reads, chrX (100 bp) gets 2 -> autosome cov 10x, X cov 2x,
# ratio 0.2 -> Male, high confidence. Single mapped reads (flag 0).
{
  echo -e "@HD\tVN:1.6\tSO:coordinate"
  echo -e "@SQ\tSN:chr1\tLN:100"
  echo -e "@SQ\tSN:chrX\tLN:100"
  sr() { echo -e "$1\t0\t$2\t$3\t60\t10M\t*\t0\t0\t$seq10\t$qual10"; }
  for i in $(seq 1 10); do sr "a_$i" chr1 $(( (i-1)*9 + 1 )); done
  sr "x_1" chrX 1
  sr "x_2" chrX 40
} > sex.sam
samtools sort -o sex.bam sex.sam && samtools index sex.bam
# CRAM counterpart: reads are all "A", so a chr1+chrX reference of 100 A's each matches.
{ echo ">chr1"; printf 'A%.0s' $(seq 1 100); echo; echo ">chrX"; printf 'A%.0s' $(seq 1 100); echo; } > sexref.fa
samtools faidx sexref.fa
samtools view -T sexref.fa -C -o sex.cram sex.bam && samtools index sex.cram

# ---- sv.bam ---------------------------------------------------------------
# chr1 + chr2 (5000 bp each, bin size 1000 -> 5 bins). Evidence:
#   inter pair  (chr1:100 <-> chr2:200)         -> 2 InterChromosomal (one per mate)
#   big-insert  (chr1:300 <-> chr1:4000, TLEN 5000) -> 2 InsertSizeOutlier
#   split read  (chr1:1000, 20S30M + SA tag)    -> 1 SplitRead (clip 20)
# Depth bins: chr1 = [2,1,0,0,1], chr2 = [1,0,0,0,0].
seq50="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
qual50="IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
{
  echo -e "@HD\tVN:1.6\tSO:coordinate"
  echo -e "@SQ\tSN:chr1\tLN:5000"
  echo -e "@SQ\tSN:chr2\tLN:5000"
  # name flag rname pos mapq cigar rnext pnext tlen seq qual [extra]
  row() { echo -e "$1\t$2\t$3\t$4\t60\t$5\t$6\t$7\t$8\t$seq50\t$qual50${9:+\t$9}"; }
  row r_inter  65  chr1 100  50M  chr2 200  0
  row r_big    97  chr1 300  50M  =    4000 5000
  row r_split  0   chr1 1000 20S30M '*' 0  0    "SA:Z:chr1,2000,+,30M20S,60,0"
  row r_big2   145 chr1 4000 50M  =    300  -5000
  row r_inter2 129 chr2 200  50M  chr1 100  0
} > sv.sam
samtools sort -o sv.bam sv.sam && samtools index sv.bam

# ---- diploid.bam ------------------------------------------------------------
# chr1 (10 bp), two haplotypes at depth 20 (10 reads each):
#   H1 = ACGTACGAAC, H2 = AGGTTCGAAC  (conceptual ref ACGTACGTAC)
#   pos2 C/G het, pos5 A/T het, pos8 T->A hom-alt, pos1 hom-ref. For the GL genotyper.
{
  echo -e "@HD\tVN:1.6\tSO:coordinate"
  echo -e "@SQ\tSN:chr1\tLN:10"
  dr() { echo -e "$1\t0\tchr1\t1\t60\t10M\t*\t0\t0\t$2\tIIIIIIIIII"; }
  for i in $(seq 1 10); do dr "h1_$i" ACGTACGAAC; done
  for i in $(seq 1 10); do dr "h2_$i" AGGTTCGAAC; done
} > diploid.sam
samtools sort -o diploid.bam diploid.sam && samtools index diploid.bam

# ---- indel.bam --------------------------------------------------------------
# chrM: a heterozygous 2 bp deletion of ref pos 6-7 (C,G). 10 ref reads (50M) + 10 deletion reads
# (5M2D43M), depth 20 -> the de-novo diploid caller emits REF=ACG / ALT=A, GT 0/1 at pos 5.
REF="ACGTACGTACGTACGTACGTACGTNACGTACGTACGTACGTACGTACGTA"
DEL="${REF:0:5}${REF:7:43}"
QREF=$(printf 'I%.0s' $(seq 1 50)); QDEL=$(printf 'I%.0s' $(seq 1 48))
{
  echo -e "@HD\tVN:1.6\tSO:coordinate"
  echo -e "@SQ\tSN:chrM\tLN:50"
  for i in $(seq 1 10); do echo -e "ref$i\t0\tchrM\t1\t60\t50M\t*\t0\t0\t$REF\t$QREF"; done
  for i in $(seq 1 10); do echo -e "del$i\t0\tchrM\t1\t60\t5M2D43M\t*\t0\t0\t$DEL\t$QDEL"; done
} > indel.sam
samtools sort -o indel.bam indel.sam && samtools index indel.bam

rm -f coverage.sam paired.sam sex.sam sv.sam diploid.sam indel.sam
echo "wrote ref.fa(.fai), coverage.bam, paired.bam, sex.bam, sv.bam, diploid.bam, indel.bam (+ .bai)"
