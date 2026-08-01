"""Full calibration: sweep the emission ratio AND the three post-hoc thresholds, train/test split.

`archaic_ratio` is the one parameter that cannot be re-filtered from a finished run -- it changes
the emissions, so the HMM has to be re-decoded. The probe sweeps it in-process (one reference read
per sample, not per value) and writes `{ratio: result}`; this searches the joint grid.

Fitted on TRAIN only. Every number reported as performance comes from TEST.
"""

import collections
import glob
import json
import math
import os
import random

SC = os.path.dirname(os.path.abspath(__file__))
SWEEP = os.path.join(SC, 'v2', 'sweep')
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


def segs(doc, mp, ms, mb):
    oi = collections.defaultdict(list)
    for s in doc['segments']:
        if s['contig'] not in CONTIGS or s['posterior'] < mp or s['n_private'] < ms:
            continue
        if s['end'] - s['start'] < mb:
            continue
        oi[s['contig']].append((s['start'], s['end']))
    return {c: union(v) for c, v in oi.items()} if oi else {}


def score(samples, sweeps, truths, ratio, mp, ms, mb):
    n = d_t = d_o = 0
    t_mb, o_mb = [], []
    for s in samples:
        doc = sweeps[s].get(ratio)
        if doc is None:
            continue
        ours = segs(doc, mp, ms, mb)
        tb, ob = bp(truths[s]), (bp(ours) if ours else 0)
        n += intersect(ours, truths[s]) if ours else 0
        d_t += tb
        d_o += ob
        t_mb.append(tb / 1e6)
        o_mb.append(ob / 1e6)
    sens = n / d_t if d_t else 0.0
    prec = n / d_o if d_o else 0.0
    f1 = 2 * sens * prec / (sens + prec) if (sens + prec) else 0.0
    return {'f1': f1 * 100, 'sens': sens * 100, 'prec': prec * 100,
            'ratio_mb': (sum(o_mb) / sum(t_mb)) if sum(t_mb) else 0.0,
            'r': pearson(t_mb, o_mb), 't': t_mb, 'o': o_mb}


def perm_p(x, y, draws=20000):
    obs = abs(pearson(x, y))
    rnd = random.Random(0)
    yy = list(y)
    hits = sum(1 for _ in range(draws) if (rnd.shuffle(yy), abs(pearson(x, yy)) >= obs)[1])
    return (hits + 1) / (draws + 1)


def main():
    sweeps, truths = {}, {}
    for f in sorted(glob.glob(f'{SWEEP}/*.json')):
        s = os.path.basename(f)[:-5]
        t = f'{V2}/truth_{s}.bed'
        if os.path.exists(t):
            sweeps[s] = json.load(open(f))
            truths[s] = load_bed(t, tol=1000)
    samples = sorted(sweeps)
    if len(samples) < 40:
        print(f'only {len(samples)} swept runs so far')
        return
    ratios = sorted(next(iter(sweeps.values())).keys(), key=float)

    rnd = random.Random(20260731)
    sh = samples[:]
    rnd.shuffle(sh)
    h = len(sh) // 2
    train, test = sorted(sh[:h]), sorted(sh[h:])
    print(f'n = {len(samples)}   train {len(train)}   test {len(test)}   ratios {ratios}\n')

    grid = [(r, mp, ms, mb)
            for r in ratios
            for mp in (0.90, 0.95, 0.98)
            for ms in (8, 16, 24, 32)
            for mb in (5_000, 10_000)]
    scored = sorted(((score(train, sweeps, truths, *g), g) for g in grid), key=lambda x: -x[0]['f1'])

    print('TOP 8 ON TRAIN (fitted)')
    print(f"  {'ratio':>5s} {'post':>5s} {'sites':>5s} {'minbp':>6s}   {'F1':>6s} {'sens':>6s} "
          f"{'prec':>6s} {'ext':>5s} {'r':>7s}")
    for sc, g in scored[:8]:
        print(f'  {g[0]:>5s} {g[1]:5.2f} {g[2]:5d} {g[3]:6d}   {sc["f1"]:5.1f}% {sc["sens"]:5.1f}% '
              f'{sc["prec"]:5.1f}% {sc["ratio_mb"]:5.2f} {sc["r"]:+7.3f}')

    best = scored[0][1]
    print(f'\nCHOSEN (best TRAIN F1): archaic_ratio {best[0]}, min_posterior {best[1]}, '
          f'min_sites {best[2]}, min_segment_bp {best[3]}')
    prev = ('3.04', 0.95, 24, 5_000)
    print('\n                              F1     sens    prec    ext        r')
    for label, samp in (('TRAIN (fitted)', train), ('TEST  (held out)', test)):
        print(f'  {label}')
        for tag, g in (('previous (ratio 3.04)', prev), ('new      (swept)     ', best)):
            sc = score(samp, sweeps, truths, *g)
            print(f'    {tag} {sc["f1"]:5.1f}% {sc["sens"]:6.1f}% {sc["prec"]:6.1f}% '
                  f'{sc["ratio_mb"]:6.2f} {sc["r"]:+8.3f}')
    tr = score(train, sweeps, truths, *best)['f1']
    te_s = score(test, sweeps, truths, *best)
    gap = tr - te_s['f1']
    verdict = ('OVERFIT — the grid found the split' if gap > 5
               else 'test split is easier; split variation, not overfitting' if gap < -5
               else 'generalises')
    print(f'\n  train F1 {tr:.1f}%  test F1 {te_s["f1"]:.1f}%  gap {gap:+.1f}  ({verdict})')
    print(f'  held-out extent r = {te_s["r"]:+.3f}, permutation p = '
          f'{perm_p(te_s["t"], te_s["o"]):.4f}')


if __name__ == '__main__':
    main()
