"""Calibrate the reference-based archaic caller on a TRAIN split, report on a held-out TEST split.

The density caller looked validated because it was tuned until a cohort statistic matched, and the
statistic was then reported as evidence. The split exists so that cannot happen again: every number
quoted as performance comes from individuals whose data never touched the fit.

Objective is base-level F1 against the external callset. F1 rather than sensitivity because
sensitivity alone is bought by calling more sequence -- the current caller over-calls 2.2x and still
scores 45% sensitivity. F1 makes over-calling cost something.

Reading the output: TRAIN F1 is the fit, TEST F1 is the claim. A large gap between them means the
grid found the split rather than the signal.
"""

import collections
import glob
import json
import math
import os
import random
import subprocess
import sys

SC = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.join(SC, 'v2')
CONTIGS = ('chr21', 'chr22')


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


def load_bed(path, tol=0):
    d = collections.defaultdict(list)
    for line in open(path):
        p = line.split()
        if len(p) >= 3:
            d[p[0]].append((int(p[1]), int(p[2])))
    return {c: union(v, tol) for c, v in d.items()}


def bp(d):
    return sum(e - s for v in d.values() for s, e in v)


def intersect(a, b):
    tot = 0
    for c, av in a.items():
        for s, e in av:
            for s2, e2 in b.get(c, []):
                lo, hi = max(s, s2), min(e, e2)
                if hi > lo:
                    tot += hi - lo
    return tot


def pearson(x, y):
    n = len(x)
    if n < 3:
        return 0.0
    mx, my = sum(x) / n, sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    dx = math.sqrt(sum((a - mx) ** 2 for a in x))
    dy = math.sqrt(sum((b - my) ** 2 for b in y))
    return num / (dx * dy) if dx and dy else 0.0


def segments_from(doc, min_post, min_sites, min_bp):
    """Re-filter a caller run at new thresholds.

    The HMM posterior is fixed per run, so the three post-hoc thresholds can be swept without
    re-running the caller. `archaic_ratio` cannot -- it changes the emissions -- so it is swept by
    re-running the probe, outside this function.
    """
    oi = collections.defaultdict(list)
    for s in doc['segments']:
        if s['contig'] not in CONTIGS:
            continue
        if s['posterior'] < min_post or s['n_private'] < min_sites:
            continue
        if s['end'] - s['start'] < min_bp:
            continue
        oi[s['contig']].append((s['start'], s['end']))
    return {c: union(v) for c, v in oi.items()} if oi else {}


def score(samples, runs, truths, min_post, min_sites, min_bp):
    sens_n = sens_d = prec_d = 0
    t_mb, o_mb = [], []
    for s in samples:
        ours = segments_from(runs[s], min_post, min_sites, min_bp)
        truth = truths[s]
        tb, ob = bp(truth), bp(ours)
        ov = intersect(ours, truth) if ours else 0
        sens_n += ov
        sens_d += tb
        prec_d += ob
        t_mb.append(tb / 1e6)
        o_mb.append(ob / 1e6)
    sens = sens_n / sens_d if sens_d else 0.0
    prec = sens_n / prec_d if prec_d else 0.0
    f1 = 2 * sens * prec / (sens + prec) if (sens + prec) else 0.0
    return {
        'sens': sens * 100, 'prec': prec * 100, 'f1': f1 * 100,
        'ratio': (sum(o_mb) / sum(t_mb)) if sum(t_mb) else 0.0,
        'r': pearson(t_mb, o_mb),
    }


def main():
    runs, truths = {}, {}
    for f in sorted(glob.glob(f'{V2}/*.match.json')):
        s = os.path.basename(f)[:-len('.match.json')]
        t = f'{V2}/truth_{s}.bed'
        if os.path.exists(t):
            runs[s] = json.load(open(f))
            truths[s] = load_bed(t, tol=1000)
    samples = sorted(runs)
    if len(samples) < 20:
        print(f'only {len(samples)} scored runs — need the cohort first')
        return

    rnd = random.Random(20260731)
    shuffled = samples[:]
    rnd.shuffle(shuffled)
    half = len(shuffled) // 2
    train, test = sorted(shuffled[:half]), sorted(shuffled[half:])
    print(f'n = {len(samples)}   train {len(train)}   test {len(test)}   (seed fixed)\n')

    grid = []
    for mp in (0.80, 0.85, 0.90, 0.95, 0.98):
        for ms in (8, 12, 16, 24, 32):
            for mb in (5_000, 10_000, 20_000, 40_000):
                grid.append((mp, ms, mb))

    scored = [(score(train, runs, truths, *g), g) for g in grid]
    scored.sort(key=lambda x: -x[0]['f1'])

    print('TOP 5 ON TRAIN (fitted)')
    print(f"  {'post':>5s} {'sites':>5s} {'minbp':>7s}   {'F1':>6s} {'sens':>6s} {'prec':>6s} "
          f"{'ratio':>6s} {'r':>6s}")
    for sc, g in scored[:5]:
        print(f'  {g[0]:5.2f} {g[1]:5d} {g[2]:7d}   {sc["f1"]:5.1f}% {sc["sens"]:5.1f}% '
              f'{sc["prec"]:5.1f}% {sc["ratio"]:6.2f} {sc["r"]:+6.3f}')

    best = scored[0][1]
    base = (0.80, 8, 5_000)
    print(f'\nCHOSEN (best TRAIN F1): min_posterior {best[0]}, min_sites {best[1]}, '
          f'min_segment_bp {best[2]}')
    print('\n                        F1     sens    prec   ratio       r')
    for label, samp in (('TRAIN (fitted)', train), ('TEST  (held out)', test)):
        b = score(samp, runs, truths, *base)
        c = score(samp, runs, truths, *best)
        print(f'  {label:22s}')
        print(f'    before (defaults)  {b["f1"]:5.1f}% {b["sens"]:6.1f}% {b["prec"]:6.1f}% '
              f'{b["ratio"]:6.2f} {b["r"]:+7.3f}')
        print(f'    after  (calibrated){c["f1"]:5.1f}% {c["sens"]:6.1f}% {c["prec"]:6.1f}% '
              f'{c["ratio"]:6.2f} {c["r"]:+7.3f}')
    tr = score(train, runs, truths, *best)['f1']
    te = score(test, runs, truths, *best)['f1']
    gap = tr - te
    # Overfitting is train >> test. Test scoring HIGHER is split-to-split variation, not the grid
    # having found the split -- reading any large gap as overfitting would be its own error.
    if gap > 5:
        verdict = 'OVERFIT — the grid found the split, not the signal'
    elif gap < -5:
        verdict = 'test split is easier; not overfitting, but the gap is split variation at n=30'
    else:
        verdict = 'generalises'
    print(f'\n  train F1 {tr:.1f}%  test F1 {te:.1f}%  gap {gap:+.1f} points  ({verdict})')

    # Is the held-out extent correlation real, or n=30 noise?
    t_mb, o_mb = [], []
    for s in test:
        ours = segments_from(runs[s], *best)
        t_mb.append(bp(truths[s]) / 1e6)
        o_mb.append(bp(ours) / 1e6 if ours else 0.0)
    r = pearson(t_mb, o_mb)
    rnd = random.Random(0)
    yy = list(o_mb)
    hits = 0
    for _ in range(20000):
        rnd.shuffle(yy)
        if abs(pearson(t_mb, yy)) >= abs(r):
            hits += 1
    print(f'  held-out extent correlation r = {r:+.3f}, permutation p = {(hits + 1) / 20001:.4f}'
          f'   (density caller: r = -0.018, p = 0.94)')


if __name__ == '__main__':
    main()
