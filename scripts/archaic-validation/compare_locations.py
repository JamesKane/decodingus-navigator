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


def load_bed(path, keep=None, tol=0):
    iv = defaultdict(list)
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3 or p[0].startswith("#"):
            continue
        if keep and p[0] not in keep:
            continue
        iv[p[0]].append((int(p[1]), int(p[2])))
    return {c: union(v, tol) for c, v in iv.items()}


def union(iv, tol=0):
    iv = sorted(iv)
    out, cs, ce = [], *iv[0]
    for s, e in iv[1:]:
        if s > ce + tol:
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


def null_sensitivity(ours, truth, truth_bp, draws=400):
    """Sensitivity a caller of OUR extent would score by placing its segments at random.

    Sensitivity rises with how much sequence you call, so the raw number is meaningless without
    this. The density caller scored 2.1% against a 5.0% null -- below chance -- and reporting only
    the 2.1% would have looked like weak performance rather than none.
    """
    import random
    rnd = random.Random(0)
    lens = {c: [e - s for s, e in v] for c, v in ours.items()}
    vals = []
    for _ in range(draws):
        r = {}
        for c, L in lens.items():
            if c not in truth:
                continue
            lo = min(s for s, _ in truth[c])
            hi = max(e for _, e in truth[c])
            r[c] = union([(p, p + l) for l in L for p in [rnd.randint(lo, max(lo, hi - l))]])
        vals.append(intersect(r, truth) / truth_bp * 100)  # both in bp
    vals.sort()
    return sum(vals) / len(vals), vals[int(0.95 * len(vals)) - 1], vals[-1]


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
    # tol=1000: the lift splits tracts at median 2 bp gaps; a strict union reports 423
    # shards where there are 48 real tracts, which makes per-tract recovery meaningless.
    truth = load_bed(truth_path, contigs, tol=1000)

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
    if t_mb and ours:
        mean, p95, mx = null_sensitivity(ours, truth, total(truth))
        sens = inter / t_mb * 100
        print()
        print("4. AGAINST THE NULL  (sensitivity alone rises with how much you call)")
        print(f"   random placement    : mean {mean:5.1f}%   p95 {p95:5.1f}%   max {mx:5.1f}%")
        verdict = ("ABOVE the null's full range" if sens > mx
                   else "above p95" if sens > p95
                   else "AT OR BELOW CHANCE")
        print(f"   observed {sens:5.1f}%      : {verdict}")


if __name__ == "__main__":
    main()
