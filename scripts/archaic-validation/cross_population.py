"""Does the caller transfer to East Asians, using parameters fitted only on Europeans?

This is the sharp test, and it asks something the calibration could not have bought. Two things
have to hold:

  1. TRANSFER -- per-individual locations and extent hold up on a population the thresholds were
     never fitted to. If performance collapses, the calibration learned European structure.

  2. THE POPULATION PREDICTION -- hmmix's own data puts East Asian archaic extent at 2.45 Mb
     against Europe's 2.09, a ratio of ~1.18. Reproducing that ORDERING is a prediction the caller
     was never shown: nothing in the fit knows which population a sample comes from. A caller that
     merely reproduces whatever it was tuned on would return the same number for both.

Parameters are frozen at the European-fitted values; nothing here is refitted.
"""

import collections
import glob
import json
import math
import os
import random

SC = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.join(SC, 'v2')
CONTIGS = ('chr21', 'chr22')
# European-fitted, frozen.
RATIO, MIN_POST, MIN_SITES, MIN_BP = '4.5', 0.98, 16, 5_000


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


def segs(doc):
    oi = collections.defaultdict(list)
    for s in doc['segments']:
        if s['contig'] not in CONTIGS or s['posterior'] < MIN_POST or s['n_private'] < MIN_SITES:
            continue
        if s['end'] - s['start'] < MIN_BP:
            continue
        oi[s['contig']].append((s['start'], s['end']))
    return {c: union(v) for c, v in oi.items()} if oi else {}


def pearson(x, y):
    n = len(x)
    if n < 3:
        return 0.0
    mx, my = sum(x) / n, sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    dx = math.sqrt(sum((a - mx) ** 2 for a in x))
    dy = math.sqrt(sum((b - my) ** 2 for b in y))
    return num / (dx * dy) if dx and dy else 0.0


def null_p95(ours, truth, tbp, draws=200, seed=0):
    rnd = random.Random(seed)
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
        vals.append(intersect(r, truth) / tbp * 100)
    vals.sort()
    return sum(vals) / len(vals), vals[int(0.95 * len(vals)) - 1]


def group(samples, source):
    rows = []
    for s in samples:
        f = f'{source}/{s}.json'
        if not os.path.exists(f):
            continue
        doc = json.load(open(f))
        doc = doc[RATIO] if RATIO in doc else doc
        ours = segs(doc)
        truth = load_bed(f'{V2}/truth_{s}.bed', tol=1000)
        if not truth:
            continue
        tb, ob = bp(truth), (bp(ours) if ours else 0)
        ov = intersect(ours, truth) if ours else 0
        mean, p95 = null_p95(ours, truth, tb) if ours else (0, 0)
        rows.append({'s': s, 't': tb / 1e6, 'o': ob / 1e6, 'ov': ov,
                     'sens': ov / tb * 100, 'prec': (ov / ob * 100) if ob else 0,
                     'null': mean, 'p95': p95})
    return rows


def summarize(rows, label):
    if not rows:
        print(f'{label}: no data')
        return None
    sens = sum(r['ov'] for r in rows) / sum(r['t'] * 1e6 for r in rows) * 100
    prec = sum(r['ov'] for r in rows) / max(sum(r['o'] * 1e6 for r in rows), 1) * 100
    f1 = 2 * sens * prec / (sens + prec) if sens + prec else 0
    beat = sum(1 for r in rows if r['sens'] > r['p95'])
    t = [r['t'] for r in rows]
    o = [r['o'] for r in rows]
    print(f'{label:16s} n={len(rows):3d}  F1 {f1:5.1f}%  sens {sens:5.1f}%  prec {prec:5.1f}%  '
          f'ext {sum(o) / sum(t):5.2f}  r {pearson(t, o):+6.3f}  above-null {beat}/{len(rows)}')
    return {'truth_mean': sum(t) / len(t), 'ours_mean': sum(o) / len(o), 'r': pearson(t, o),
            'n': len(rows), 't': t, 'o': o}


def main():
    eur = [l.strip() for l in open(f'{SC}/eur60.txt') if l.strip()]
    eas = [l.strip() for l in open(f'{SC}/eas30.txt') if l.strip()]
    print(f'parameters FROZEN at the European fit: ratio {RATIO}, posterior {MIN_POST}, '
          f'sites {MIN_SITES}, min_bp {MIN_BP}\n')

    e = summarize(group(eur, f'{V2}/sweep'), 'EUROPE (fitted)')
    a = summarize(group(eas, f'{V2}/sweep_eas'), 'EAST ASIA (new)')
    if not (e and a):
        return

    print('\nTHE POPULATION PREDICTION (nothing in the fit knows a sample\'s population)')
    print(f'  hmmix truth  EAS {a["truth_mean"]:.3f} Mb / EUR {e["truth_mean"]:.3f} Mb'
          f'  = {a["truth_mean"] / e["truth_mean"]:.3f}x')
    print(f'  our calls    EAS {a["ours_mean"]:.3f} Mb / EUR {e["ours_mean"]:.3f} Mb'
          f'  = {a["ours_mean"] / e["ours_mean"]:.3f}x')
    print('  published expectation: East Asians carry ~1.2x the Neanderthal ancestry of Europeans.')

    # Is the elevation real, or within noise? Permutation over the pooled samples.
    pooled = e['o'] + a['o']
    labels = [0] * len(e['o']) + [1] * len(a['o'])
    obs = (sum(a['o']) / len(a['o'])) - (sum(e['o']) / len(e['o']))
    rnd = random.Random(0)
    hits = 0
    for _ in range(20000):
        rnd.shuffle(labels)
        ea = [v for v, l in zip(pooled, labels) if l == 1]
        eu = [v for v, l in zip(pooled, labels) if l == 0]
        if (sum(ea) / len(ea)) - (sum(eu) / len(eu)) >= obs:
            hits += 1
    print(f'  our EAS-EUR difference {obs:+.3f} Mb, permutation p = {(hits + 1) / 20001:.4f}')


if __name__ == '__main__':
    main()
