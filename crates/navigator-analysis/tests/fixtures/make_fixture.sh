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

rm -f coverage.sam paired.sam sex.sam
echo "wrote ref.fa(.fai), coverage.bam, paired.bam, sex.bam (+ .bai)"
