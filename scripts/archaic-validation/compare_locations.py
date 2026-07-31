#!/usr/bin/env python3
"""Compare Navigator's archaic segment calls against the lifted hmmix truth for one sample.

Usage: compare_archaic.py <ours.json> <truth.bed> [contig ...]

`ours.json` is `navigator archaic-segments --json` output; `truth.bed` is the hmmix callset for the
same individual, lifted hg38 -> CHM13 and unioned across haplotypes (hmmix reports per haplotype,
our caller is unphased, so union -- summing would double their figure).

Reports three things, in increasing order of how hard they are to fake:
  1. total extent            -- a caller that only matches the cohort mean passes this
  2. base-level overlap      -- did we find the SAME sequence, not just the same amount
  3. per-segment recovery    -- how many of their tracts we hit at all
"""
import json
import sys
from collections import defaultdict


def load_bed(path, keep=None):
    iv = defaultdict(list)
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3 or p[0].startswith("#"):
            continue
        if keep and p[0] not in keep:
            continue
        iv[p[0]].append((int(p[1]), int(p[2])))
    return {c: union(v) for c, v in iv.items()}


def union(iv):
    iv = sorted(iv)
    out, cs, ce = [], *iv[0]
    for s, e in iv[1:]:
        if s > ce:
            out.append((cs, ce))
            cs, ce = s, e
        else:
            ce = max(ce, e)
    out.append((cs, ce))
    return out


def total(d):
    return sum(e - s for v in d.values() for s, e in v)


def intersect(a, b):
    """Base-level intersection of two contig->intervals dicts."""
    out = 0
    for c, av in a.items():
        bv = b.get(c)
        if not bv:
            continue
        i = j = 0
        while i < len(av) and j < len(bv):
            lo, hi = max(av[i][0], bv[j][0]), min(av[i][1], bv[j][1])
            if hi > lo:
                out += hi - lo
            if av[i][1] < bv[j][1]:
                i += 1
            else:
                j += 1
    return out


def main():
    ours_path, truth_path = sys.argv[1], sys.argv[2]
    contigs = set(sys.argv[3:]) or None

    doc = json.load(open(ours_path))
    segs = doc.get("segments", doc if isinstance(doc, list) else [])
    ours_iv = defaultdict(list)
    for s in segs:
        c = s.get("contig") or s.get("chrom")
        if contigs and c not in contigs:
            continue
        ours_iv[c].append((int(s["start"]), int(s["end"])))
    ours = {c: union(v) for c, v in ours_iv.items()} if ours_iv else {}
    truth = load_bed(truth_path, contigs)

    o_mb, t_mb = total(ours) / 1e6, total(truth) / 1e6
    inter = intersect(ours, truth) / 1e6
    union_mb = o_mb + t_mb - inter

    print(f"contigs compared      : {sorted(set(ours) | set(truth))}")
    print()
    print("1. EXTENT")
    print(f"   ours                : {o_mb:7.3f} Mb  in {sum(len(v) for v in ours.values()):5d} segments")
    print(f"   hmmix (truth)       : {t_mb:7.3f} Mb  in {sum(len(v) for v in truth.values()):5d} merged tracts")
    print(f"   ratio ours/theirs   : {o_mb / t_mb:7.3f}" if t_mb else "   ratio: n/a")
    print()
    print("2. BASE-LEVEL AGREEMENT  (the test extent alone cannot pass)")
    print(f"   overlap             : {inter:7.3f} Mb")
    print(f"   sensitivity         : {inter / t_mb * 100:6.1f}%  of their archaic bases we also call")
    print(f"   precision           : {inter / o_mb * 100:6.1f}%  of our archaic bases they also call" if o_mb else "")
    print(f"   Jaccard             : {inter / union_mb:7.3f}" if union_mb else "")
    print()
    print("3. PER-TRACT RECOVERY")
    hit = miss = 0
    for c, tv in truth.items():
        ov = ours.get(c, [])
        for s, e in tv:
            if any(min(e, oe) > max(s, os_) for os_, oe in ov):
                hit += 1
            else:
                miss += 1
    print(f"   their tracts hit    : {hit}/{hit + miss}  ({hit / (hit + miss) * 100:.1f}%)" if hit + miss else "   n/a")
    print()
    print("   A random caller with our extent would score sensitivity ~= our_Mb / callable_Mb,")
    print("   i.e. a few percent. Sensitivity near that floor means we match the AMOUNT but not")
    print("   the LOCATION -- which would mean the 1.01x cohort agreement was luck.")


if __name__ == "__main__":
    main()
