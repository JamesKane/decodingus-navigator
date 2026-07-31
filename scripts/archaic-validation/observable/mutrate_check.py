"""Does African-outgroup site density explain the background variation the emission model ignores?

If it does, it is a usable mutation-rate normalizer and the fix is an asset build. If it does not,
a different proxy is needed and building this one would waste the effort.

Two questions, in order:
  1. Correlation between outgroup density and private-variant density in BACKGROUND regions
     (callable, non-archaic). This is the normalizer's whole job.
  2. Whether normalizing by it actually flattens the background — the overdispersion should fall
     from the measured 14.6x toward 1.0x. Correlation alone does not guarantee that.
"""

import bisect
import collections
import math

BIN = 100_000
CONTIGS = ('chr21', 'chr22')


def merge(iv, tol=0):
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
        c, s, e = line.split()[:3]
        d[c].append((int(s), int(e)))
    return {c: merge(v, tol) for c, v in d.items()}


def covered(regions, s, e):
    tot = 0
    for s2, e2 in regions:
        lo, hi = max(s, s2), min(e, e2)
        if hi > lo:
            tot += hi - lo
    return tot


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    dx = math.sqrt(sum((a - mx) ** 2 for a in x))
    dy = math.sqrt(sum((b - my) ** 2 for b in y))
    return num / (dx * dy) if dx and dy else 0.0


def main():
    priv = collections.defaultdict(list)
    for c in CONTIGS:
        with open(f'private.{c}.tsv') as f:
            next(f)
            for line in f:
                p = line.split('\t')
                priv[p[0]].append(int(p[1]))
    for c in priv:
        priv[c].sort()

    og = collections.defaultdict(collections.Counter)
    with open('og_density.tsv') as f:
        next(f)
        for line in f:
            c, w, k = line.split('\t')
            og[c][int(w)] = int(k)

    callable_r = load_bed('callable.bed')
    truth = load_bed('truth_HG00096.chm13.bed', 1000)

    xs, ys, weights = [], [], []
    for c in CONTIGS:
        a = priv[c]
        lo = min(s for s, _ in callable_r[c])
        hi = max(e for _, e in callable_r[c])
        for b in range(lo, hi, BIN):
            e = b + BIN
            cal = covered(callable_r[c], b, e)
            if cal < BIN * 0.5 or covered(truth.get(c, []), b, e) > 0:
                continue
            n = bisect.bisect_left(a, e) - bisect.bisect_left(a, b)
            ogn = sum(og[c].get(w, 0) for w in range(b, e, 1000))
            if ogn == 0:
                continue
            xs.append(ogn / (cal / 1e6))     # outgroup sites per callable Mb
            ys.append(n / (cal / 1e6))       # our private variants per callable Mb
            weights.append(cal)

    r = pearson(xs, ys)
    print(f'BACKGROUND bins ({BIN // 1000} kb, callable, non-archaic): n = {len(xs)}')
    print(f'  outgroup-site density vs private-variant density:  Pearson r = {r:+.3f}'
          f'   (r^2 = {r * r:.2f})')
    print(f'  -> outgroup density explains {r * r * 100:.0f}% of the background variance')
    print()

    # Does dividing by the proxy actually flatten it?
    mean_x = sum(xs) / len(xs)
    raw = ys
    norm = [y / (x / mean_x) for x, y in zip(xs, ys)]

    def stats(v, label):
        v = sorted(v)
        n = len(v)
        m = sum(v) / n
        q = lambda f: v[int(f * (n - 1))]
        # overdispersion of the implied per-bin count
        cnt = [d * BIN / 1e6 for d in v]
        cm = sum(cnt) / len(cnt)
        cv = sum((z - cm) ** 2 for z in cnt) / (len(cnt) - 1)
        print(f'  {label:22s} p10 {q(.1):6.0f}  median {q(.5):6.0f}  p90 {q(.9):6.0f}'
              f'   p90/p10 {q(.9) / max(q(.1), 1e-9):5.1f}x   overdispersion {cv / cm:5.1f}x')

    print('BACKGROUND FLATNESS (the emission model assumes this is flat):')
    stats(raw, 'raw density')
    stats(norm, 'normalized by proxy')
    print()
    print('  The archaic signal is 2.89x. For the model to separate it, the background spread')
    print('  must be well below that. Overdispersion near 1.0x is what a Poisson emission needs.')


if __name__ == '__main__':
    main()
