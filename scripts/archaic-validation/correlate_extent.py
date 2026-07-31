#!/usr/bin/env python3
"""Does Tier B's archaic extent track the truth PER PERSON, or only on average?

The shipped validation showed one genome landing at 1.01x the cohort MEAN. A caller whose output
is pure noise, scaled to the right average, passes that test. This asks the question that separates
the two: across individuals, does our extent rise and fall with hmmix's?

Reports Pearson r (linear agreement) and Spearman rho (rank agreement, robust to a scale error),
each with a permutation p-value, plus the same for segment COUNT. A wide spread in the truth is what
gives the test power, so the observed spread is printed too.
"""
import collections
import glob
import json
import math
import os
import random

SC = os.environ.get("ARCHAIC_RUNS", os.path.dirname(os.path.abspath(__file__)))


def union_bp(iv):
    iv = sorted(iv)
    tot, cs, ce = 0, *iv[0]
    for s, e in iv[1:]:
        if s > ce:
            tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + ce - cs


def truth_extent():
    """hmmix extent per individual on chr21+22, haplotypes UNIONED (they report per haplotype)."""
    src = os.path.expanduser('~/.decodingus/ancestry-build/tmp/hmmix_segments_chr21_22.tsv')
    segs = collections.defaultdict(list)
    nseg = collections.Counter()
    with open(src) as f:
        next(f)
        for line in f:
            p = line.rstrip('\n').split('\t')
            segs[(p[0], p[4])].append((int(p[5]), int(p[6])))
    per, cnt = collections.defaultdict(int), collections.Counter()
    for (n, _c), iv in segs.items():
        per[n] += union_bp(iv)
        cnt[n] += len(iv)
    return per, cnt


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    dx = math.sqrt(sum((a - mx) ** 2 for a in x))
    dy = math.sqrt(sum((b - my) ** 2 for b in y))
    return num / (dx * dy) if dx and dy else 0.0


def rank(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def perm_p(x, y, stat, draws=20000):
    obs = abs(stat(x, y))
    rnd = random.Random(0)
    yy = list(y)
    hits = 0
    for _ in range(draws):
        rnd.shuffle(yy)
        if abs(stat(x, yy)) >= obs:
            hits += 1
    return (hits + 1) / (draws + 1)


def main():
    per, cnt = truth_extent()
    rows = []
    for f in sorted(glob.glob(f'{SC}/runs/*.json')):
        name = os.path.basename(f)[:-5]
        try:
            d = json.load(open(f))
        except Exception:
            continue
        s = d.get('summary') or {}
        if name not in per or not s:
            continue
        rows.append((name, per[name] / 1e6, s['total_mb'], cnt[name], s['n_segments']))

    if len(rows) < 5:
        print(f'only {len(rows)} samples — not enough')
        return
    rows.sort(key=lambda r: r[1])
    print(f'n = {len(rows)} European individuals, chr21+22\n')
    print(f"{'sample':10s} {'hmmix Mb':>9s} {'ours Mb':>8s} {'ratio':>6s} {'hmmix segs':>11s} {'our segs':>9s}")
    for n, t, o, tc, oc in rows:
        print(f'{n:10s} {t:9.3f} {o:8.3f} {o / t:6.2f} {tc:11d} {oc:9d}')

    t = [r[1] for r in rows]
    o = [r[2] for r in rows]
    print(f'\ntruth spread: {min(t):.2f} - {max(t):.2f} Mb   ours: {min(o):.2f} - {max(o):.2f} Mb')
    print(f'mean ratio ours/theirs: {sum(o) / sum(t):.3f}  <- the statistic the shipped validation used')

    r = pearson(t, o)
    rho = pearson(rank(t), rank(o))
    print(f'\nPER-INDIVIDUAL AGREEMENT (the statistic that was never measured)')
    print(f'  Pearson  r   = {r:+.3f}   permutation p = {perm_p(t, o, pearson):.4f}')
    print(f'  Spearman rho = {rho:+.3f}   permutation p = {perm_p(rank(t), rank(o), pearson):.4f}')
    tc = [float(r_[3]) for r_ in rows]
    oc = [float(r_[4]) for r_ in rows]
    print(f'  segment count: Pearson r = {pearson(tc, oc):+.3f}')

    # Spread carries the same message without needing a correlation: a measurement of a varying
    # quantity should vary about as much as the quantity does. A calibrated constant does not.
    def sd(v):
        m = sum(v) / len(v)
        return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))

    st, so = sd(t), sd(o)
    print(f'\nSPREAD  truth SD = {st:.3f} Mb   ours SD = {so:.3f} Mb   (ours/truth = {so / st:.2f})')
    print(f'  truth CV = {st / (sum(t) / len(t)):.3f}   ours CV = {so / (sum(o) / len(o)):.3f}')
    print('\n  r ~ 0 with this spread => the total is a calibrated constant, not a measurement of')
    print('  the individual. r > 0 => the headline number carries real per-person signal even')
    print('  though the LOCATIONS do not.')


if __name__ == '__main__':
    main()
