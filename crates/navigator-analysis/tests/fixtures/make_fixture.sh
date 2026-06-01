#!/usr/bin/env bash
# Regenerates the synthetic coverage-walker fixture. Requires samtools on PATH.
#
# Reference `chrM`, 50 bp, with an N at position 25. Reads are placed to make every
# coverage/callable output hand-computable and to exercise all callable states:
#
#   pos 1-10   depth 4  MAPQ 60  -> CALLABLE
#   pos 11-20  depth 2  MAPQ 60  -> LOW_COVERAGE        (qc-pass 2 < min_depth 4)
#   pos 21-24  depth 0           -> NO_COVERAGE
#   pos 25     depth 0  ref N    -> REF_N
#   pos 26-30  depth 5  MAPQ 0   -> POOR_MAPPING_QUALITY (100% reads <= max_low_mapq)
#   pos 31-50  depth 0           -> NO_COVERAGE
#
# Base quality is 'I' (Phred 40) throughout, so the base-quality filter never trips;
# states are driven by depth and MAPQ alone. See coverage.rs tests for expected values.
set -euo pipefail
cd "$(dirname "$0")"

# 50 bp reference: 24 bp + N + 25 bp. faidx confirms LN below.
{
  echo ">chrM"
  echo "ACGTACGTACGTACGTACGTACGTNACGTACGTACGTACGTACGTACGTA"
} > ref.fa
samtools faidx ref.fa
ln_actual=$(cut -f2 ref.fa.fai)
[ "$ln_actual" = "50" ] || { echo "ref length is $ln_actual, expected 50" >&2; exit 1; }

# Build the SAM. 10-base reads use CIGAR 10M / SEQ+QUAL of length 10; the 5-base
# MAPQ-0 reads use 5M / length 5. Flag 0 = mapped, forward, single.
seq10="AAAAAAAAAA"; qual10="IIIIIIIIII"
seq5="AAAAA";       qual5="IIIII"
{
  echo -e "@HD\tVN:1.6\tSO:coordinate"
  echo -e "@SQ\tSN:chrM\tLN:50"
  emit() { # name pos mapq cigar seq qual
    echo -e "$1\t0\tchrM\t$2\t$3\t$4\t*\t0\t0\t$5\t$6"
  }
  for i in 1 2 3 4;     do emit "d4_$i"  1  60 10M "$seq10" "$qual10"; done
  for i in 1 2;         do emit "d2_$i"  11 60 10M "$seq10" "$qual10"; done
  for i in 1 2 3 4 5;   do emit "mq0_$i" 26 0  5M  "$seq5"  "$qual5";  done
} > reads.sam

samtools sort -o coverage.bam reads.sam
samtools index coverage.bam
rm -f reads.sam
echo "wrote ref.fa(.fai) and coverage.bam(.bai)"
