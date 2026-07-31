"""Is the background's excess variance real, or is it this caller's error rate varying by region?

The emission model needs a flat background. Measured raw, it is 14.6x overdispersed; normalizing by
African-outgroup density gets it to 7.4x, still above the 2.89x archaic signal. If the remainder is
low-confidence calls -- artifacts clustering in hard regions -- then filtering should flatten the
background AND raise the archaic enrichment, because artifacts dilute real tracts too.

Reports both numbers for each filter, since a filter that flattens the background by discarding the
signal has bought nothing.
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


def load_private():
    rows = collections.defaultdict(list)
    for c in CONTIGS:
        with open(f'private.{c}.tsv') as f:
            next(f)
            for line in f:
                p = line.rstrip('\n').split('\t')
                rows[p[0]].append((int(p[1]), int(p[2]), int(p[3]), int(p[4])))
    return rows


def assess(rows, callable_r, truth, label):
    pos = {c: sorted(r[0] for r in v) for c, v in rows.items()}
    in_n = in_bp = out_n = out_bp = 0
    dens = []
    for c in CONTIGS:
        a = pos.get(c, [])
        cnt = lambda s, e: bisect.bisect_left(a, e) - bisect.bisect_left(a, s)
        tb = sum(covered(callable_r[c], s, e) for s, e in truth[c])
        tn = sum(cnt(s, e) for s, e in truth[c])
        cb = sum(e - s for s, e in callable_r[c])
        cn = sum(cnt(s, e) for s, e in callable_r[c])
        in_n += tn
        in_bp += tb
        out_n += cn - tn
        out_bp += cb - tb
        lo = min(s for s, _ in callable_r[c])
        hi = max(e for _, e in callable_r[c])
        for b in range(lo, hi, BIN):
            e = b + BIN
            cal = covered(callable_r[c], b, e)
            if cal < BIN * 0.5 or covered(truth.get(c, []), b, e) > 0:
                continue
            dens.append(cnt(b, e) / (cal / 1e6))
    enrich = (in_n / in_bp) / (out_n / out_bp) if out_n and in_bp else 0
    dens.sort()
    n = len(dens)
    q = lambda f: dens[int(f * (n - 1))]
    counts = [d * BIN / 1e6 for d in dens]
    cm = sum(counts) / len(counts)
    cv = sum((z - cm) ** 2 for z in counts) / (len(counts) - 1)
    kept = sum(len(v) for v in rows.values())
    print(f'{label:26s} kept {kept:6d}  enrich {enrich:4.2f}x  '
          f'bg p90/p10 {q(.9) / max(q(.1), 1e-9):5.1f}x  overdisp {cv / max(cm, 1e-9):5.1f}x')


def main():
    raw = load_private()
    callable_r = load_bed('callable.bed')
    truth = load_bed('truth_HG00096.chm13.bed', 1000)

    print('filter                      variants   archaic signal   background flatness')
    assess(raw, callable_r, truth, 'none (current)')
    for gq in (20, 30, 50):
        f = {c: [r for r in v if r[2] >= gq] for c, v in raw.items()}
        assess(f, callable_r, truth, f'GQ >= {gq}')
    for dp in (10, 15, 20):
        f = {c: [r for r in v if r[3] >= dp] for c, v in raw.items()}
        assess(f, callable_r, truth, f'depth >= {dp}')
    # Homozygous-derived only: an introgressed tract is usually heterozygous, so this should HURT
    # if the signal is real -- a useful control that the enrichment is not an artifact of genotype.
    f = {c: [r for r in v if r[1] == 1] for c, v in raw.items()}
    assess(f, callable_r, truth, 'het only (control)')
    f = {c: [r for r in v if r[1] == 2] for c, v in raw.items()}
    assess(f, callable_r, truth, 'hom-alt only (control)')
    for gq, dp in ((30, 15), (50, 20)):
        f = {c: [r for r in v if r[2] >= gq and r[3] >= dp] for c, v in raw.items()}
        assess(f, callable_r, truth, f'GQ>={gq} & depth>={dp}')
    print()
    print('  Wanted: enrichment UP and overdispersion DOWN together. Overdispersion must fall well')
    print('  below the enrichment for a Poisson emission to separate the two states.')


if __name__ == '__main__':
    main()
