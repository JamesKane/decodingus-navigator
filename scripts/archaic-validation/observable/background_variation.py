"""Is the private-variant background uniform, as the caller's emission model assumes?

The model says: background windows emit Poisson(lambda) with ONE genome-wide lambda, so a window
with several private variants is evidence of an archaic tract. That inference only holds if the
non-archaic background is actually flat. hmmix does not assume this — it requires a mutation-rate
map and scales the expected density per window by it. We have no such asset.

This measures the spread of background density directly, in regions that are callable and NOT in
hmmix's archaic tracts, and asks how a Poisson model would fare against it.
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

    callable_r = load_bed('callable.bed')
    truth = load_bed('truth_HG00096.chm13.bed', 1000)

    def covered(regions, s, e):
        tot = 0
        for s2, e2 in regions:
            lo, hi = max(s, s2), min(e, e2)
            if hi > lo:
                tot += hi - lo
        return tot

    dens = []
    for c in CONTIGS:
        a = priv[c]
        lo = min(s for s, _ in callable_r[c])
        hi = max(e for _, e in callable_r[c])
        for b in range(lo, hi, BIN):
            e = b + BIN
            cal = covered(callable_r[c], b, e)
            if cal < BIN * 0.5:
                continue                        # mostly uncallable: not a background sample
            if covered(truth.get(c, []), b, e) > 0:
                continue                        # overlaps a real archaic tract: not background
            n = bisect.bisect_left(a, e) - bisect.bisect_left(a, b)
            dens.append(n / (cal / 1e6))
    dens.sort()
    n = len(dens)
    mean = sum(dens) / n
    var = sum((d - mean) ** 2 for d in dens) / (n - 1)
    q = lambda f: dens[int(f * (n - 1))]

    print(f'BACKGROUND private-variant density, {BIN // 1000} kb bins, callable and NOT archaic')
    print(f'  n bins {n}')
    print(f'  mean {mean:.0f}/Mb   median {q(.5):.0f}   p10 {q(.1):.0f}   p90 {q(.9):.0f}'
          f'   max {dens[-1]:.0f}')
    print(f'  p90/p10 spread = {q(.9) / max(q(.1), 1e-9):.1f}x')
    print()
    print(f'  enrichment inside real archaic tracts, for comparison: 2.89x')
    print(f'  fraction of BACKGROUND bins already above 2.89x the median: '
          f'{sum(1 for d in dens if d > 2.89 * q(.5)) / n * 100:.1f}%')
    print()
    # Poisson would predict variance == mean for the per-bin COUNT; overdispersion is the degree to
    # which a single-lambda model is simply the wrong distribution.
    counts = [d * BIN / 1e6 for d in dens]
    cm = sum(counts) / len(counts)
    cv = sum((x - cm) ** 2 for x in counts) / (len(counts) - 1)
    print(f'  per-bin counts: mean {cm:.1f}, variance {cv:.1f} -> overdispersion {cv / cm:.1f}x')
    print('  (Poisson assumes variance = mean, i.e. 1.0x. Anything well above that means a')
    print('   single-lambda background will call its own upper tail archaic.)')
    print()
    print(f'  sd/mean of background density = {math.sqrt(var) / mean:.2f}')


if __name__ == '__main__':
    main()
