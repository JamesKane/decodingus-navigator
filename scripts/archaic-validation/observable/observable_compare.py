"""Compare candidate observables for the Tier B HMM on signal-to-noise.

The current observable is "any private variant". Measured, that gives 2.89x enrichment inside real
archaic tracts against a background that varies 5.3x (p10-p90) and is 14.6x overdispersed — the
noise is bigger than the signal, which is why no parameter setting works.

A more specific observable should trade count for contrast: fewer observations, but a much higher
ratio inside tracts. What matters is whether the enrichment clears the background spread, because
that ratio is what any two-state model has to separate.
"""

import bisect
import collections

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


def assess(pos_by_contig, callable_r, truth, label):
    in_n = in_bp = out_n = out_bp = 0
    dens = []
    for c in CONTIGS:
        a = sorted(pos_by_contig.get(c, []))
        cnt = lambda s, e: bisect.bisect_left(a, e) - bisect.bisect_left(a, s)
        in_bp += sum(covered(callable_r[c], s, e) for s, e in truth[c])
        in_n += sum(cnt(s, e) for s, e in truth[c])
        cb = sum(e - s for s, e in callable_r[c])
        out_bp += cb - sum(covered(callable_r[c], s, e) for s, e in truth[c])
        out_n += sum(cnt(s, e) for s, e in callable_r[c]) - sum(cnt(s, e) for s, e in truth[c])
        lo = min(s for s, _ in callable_r[c])
        hi = max(e for _, e in callable_r[c])
        for b in range(lo, hi, BIN):
            e = b + BIN
            cal = covered(callable_r[c], b, e)
            if cal < BIN * 0.5 or covered(truth.get(c, []), b, e) > 0:
                continue
            dens.append(cnt(b, e) / (cal / 1e6))
    if not in_bp or not out_n:
        print(f'{label:38s} (insufficient data)')
        return
    enrich = (in_n / in_bp) / (out_n / out_bp)
    dens.sort()
    n = len(dens)
    q = lambda f: dens[int(f * (n - 1))]
    spread = q(.9) / q(.1) if q(.1) > 0 else float('inf')
    print(f'{label:38s} n={in_n + out_n:6d}  in-tract {in_n / (in_bp / 1e6):6.1f}/Mb  '
          f'bg {out_n / (out_bp / 1e6):6.1f}/Mb  ENRICH {enrich:5.2f}x  bg p90/p10 {spread:6.1f}x')


def main():
    priv = collections.defaultdict(list)
    for c in CONTIGS:
        with open(f'private.{c}.tsv') as f:
            next(f)
            for line in f:
                p = line.rstrip('\n').split('\t')
                priv[p[0]].append((int(p[1]), int(p[2]), int(p[3]), int(p[4])))

    diag = collections.defaultdict(dict)
    with open('classify.tsv') as f:
        next(f)
        for line in f:
            c, p, d, k = line.rstrip('\n').split('\t')
            diag[c][int(p)] = (d, int(k))

    callable_r = load_bed('callable.bed')
    truth = load_bed('truth_HG00096.chm13.bed', 1000)

    print('observable                             count     density in/out         contrast')
    assess({c: [r[0] for r in v] for c, v in priv.items()}, callable_r, truth,
           'all private variants (current)')

    # Private variants that also sit at a known archaic-diagnostic site.
    d_all = {c: [r[0] for r in v if r[0] in diag.get(c, {})] for c, v in priv.items()}
    assess(d_all, callable_r, truth, 'private AND archaic-diagnostic')

    for k, name in ((0, 'Neanderthal-diagnostic'), (2, 'shared-archaic')):
        f = {c: [r[0] for r in v if diag.get(c, {}).get(r[0], (None, -1))[1] == k]
             for c, v in priv.items()}
        assess(f, callable_r, truth, f'  ...of which {name}')

    # Control: diagnostic sites the subject does NOT carry should show no enrichment.
    carried = {c: {r[0] for r in v} for c, v in priv.items()}
    notc = {c: [p for p in diag.get(c, {}) if p not in carried.get(c, set())] for c in CONTIGS}
    assess(notc, callable_r, truth, 'diagnostic sites NOT carried (control)')

    print()
    print('  Contrast is what matters: the enrichment has to clear the background spread for a')
    print('  two-state model to separate the states. "all private" fails that (2.89x vs 5.3x).')


if __name__ == '__main__':
    main()
