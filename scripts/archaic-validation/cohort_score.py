"""Score the reference-based caller across the cohort, on the tests the density caller failed.

Two questions, both of which the density caller answered badly:

  LOCATIONS -- is base-level overlap above the random-placement null, per individual?
               (density: 2.1% against a 5.0% null, i.e. below chance)
  AMOUNTS   -- does extent track the individual across people?
               (density: Pearson r = -0.018, p = 0.94, against a 2.5x range of true values)

Reports both, plus precision, because over-calling inflates sensitivity and the pair has to be
read together.
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
        if len(p) < 3:
            continue
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


def null_mean_p95(ours, truth, tbp, draws=200, seed=0):
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
    for pos, i in enumerate(order):
        r[i] = pos + 1
    return r


def perm_p(x, y, draws=20000):
    obs = abs(pearson(x, y))
    rnd = random.Random(0)
    yy = list(y)
    hits = 0
    for _ in range(draws):
        rnd.shuffle(yy)
        if abs(pearson(x, yy)) >= obs:
            hits += 1
    return (hits + 1) / (draws + 1)


def main():
    rows = []
    for f in sorted(glob.glob(f'{V2}/*.match.json')):
        s = os.path.basename(f)[:-len('.match.json')]
        tpath = f'{V2}/truth_{s}.bed'
        if not os.path.exists(tpath):
            continue
        doc = json.load(open(f))
        oi = collections.defaultdict(list)
        for seg in doc['segments']:
            if seg['contig'] in CONTIGS:
                oi[seg['contig']].append((seg['start'], seg['end']))
        if not oi:
            continue
        ours = {c: union(v) for c, v in oi.items()}
        truth = load_bed(tpath, tol=1000)
        tbp, obp = bp(truth), bp(ours)
        ov = intersect(ours, truth)
        mean, p95 = null_mean_p95(ours, truth, tbp)
        rows.append({
            'sample': s, 'truth_mb': tbp / 1e6, 'ours_mb': obp / 1e6,
            'sens': ov / tbp * 100, 'prec': ov / obp * 100,
            'null': mean, 'p95': p95, 'nseg': sum(len(v) for v in ours.values()),
        })

    rows.sort(key=lambda r: r['truth_mb'])
    print(f'n = {len(rows)} Europeans, chr21+22, reference-based caller\n')
    print(f"{'sample':10s} {'truth':>7s} {'ours':>7s} {'segs':>5s} {'sens':>7s} {'prec':>6s} "
          f"{'null':>6s} {'p95':>6s}  verdict")
    beat = 0
    for r in rows:
        ok = r['sens'] > r['p95']
        beat += ok
        print(f"{r['sample']:10s} {r['truth_mb']:7.3f} {r['ours_mb']:7.3f} {r['nseg']:5d} "
              f"{r['sens']:6.1f}% {r['prec']:5.1f}% {r['null']:5.1f}% {r['p95']:5.1f}%  "
              f"{'above null' if ok else 'AT/BELOW'}")

    print(f'\nLOCATIONS: {beat}/{len(rows)} individuals score above their own p95 null')
    print(f'  mean sensitivity {sum(r["sens"] for r in rows) / len(rows):.1f}%   '
          f'mean null {sum(r["null"] for r in rows) / len(rows):.1f}%   '
          f'mean precision {sum(r["prec"] for r in rows) / len(rows):.1f}%')
    print('  (density caller: 2.1% against a 5.0% null -- below chance)')

    t = [r['truth_mb'] for r in rows]
    o = [r['ours_mb'] for r in rows]
    r_ = pearson(t, o)
    rho = pearson(rank(t), rank(o))
    print(f'\nAMOUNTS: per-individual extent correlation')
    print(f'  Pearson  r   = {r_:+.3f}   permutation p = {perm_p(t, o):.4f}')
    print(f'  Spearman rho = {rho:+.3f}   permutation p = {perm_p(rank(t), rank(o)):.4f}')
    print(f'  truth range {min(t):.2f}-{max(t):.2f} Mb   ours {min(o):.2f}-{max(o):.2f} Mb')
    print(f'  mean ratio ours/theirs = {sum(o) / sum(t):.3f}')
    print('  (density caller: r = -0.018, p = 0.94)')


if __name__ == '__main__':
    main()
