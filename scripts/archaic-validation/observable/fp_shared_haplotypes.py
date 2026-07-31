"""Do our "false positives" land where OTHER people have archaic tracts?

The carrying-rate version of this question is circular: the caller selects regions of high carrying
rate, so every segment it emits has one. This test does not use our caller's own evidence at all.

Introgressed haplotypes are SHARED -- the same archaic haplotype segregates across a population. So
if a segment we call and hmmix did not is really archaic, it should coincide with tracts hmmix calls
in OTHER individuals. If instead it is our noise, it should fall where hmmix never calls anything,
i.e. no better than a random region of the same size.

The individual's own truth is excluded from the union, so a segment cannot vindicate itself.
"""

import collections
import json
import os
import random

SC = os.path.dirname(os.path.abspath(__file__))
V2 = f'{SC}/v2'
CONTIGS = ('chr21', 'chr22')
RATIO, MIN_POST, MIN_SITES, MIN_BP = '4.5', 0.98, 16, 5_000


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
        p = line.split()
        if len(p) >= 3:
            d[p[0]].append((int(p[1]), int(p[2])))
    return {c: merge(v, tol) for c, v in d.items()}


def overlaps(regs, s, e):
    for a, b in regs:
        if min(e, b) > max(s, a):
            return True
    return False


def main():
    eur = [x.strip() for x in open(f'{SC}/eur60.txt') if x.strip()]
    eas = [x.strip() for x in open(f'{SC}/eas30.txt') if x.strip()]
    truths = {}
    for s in eur + eas:
        p = f'{V2}/truth_{s}.bed'
        if os.path.exists(p):
            truths[s] = load_bed(p, tol=1000)

    # The population-wide map of where archaic tracts occur at all, per contig.
    pool = collections.defaultdict(list)
    for s, t in truths.items():
        for c, v in t.items():
            pool[c].extend(v)
    print(f'population archaic map from {len(truths)} individuals: ' + ', '.join(
        f'{c} {sum(e - s for s, e in merge(v)) / 1e6:.1f} Mb' for c, v in sorted(pool.items())))

    callable_r = load_bed(f'{SC}/callable.bed')
    rnd = random.Random(0)
    print(f"\n{'pop':4s} {'class':16s} {'n':>5s} {'in others truth':>16s} {'random null':>12s}")
    for grp, names, source in (('EUR', eur[:25], f'{V2}/sweep'), ('EAS', eas[:25], f'{V2}/sweep_eas')):
        counts = collections.Counter()
        tot = collections.Counter()
        null_hit = null_tot = 0
        for s in names:
            f = f'{source}/{s}.json'
            if not (os.path.exists(f) and s in truths):
                continue
            doc = json.load(open(f))
            doc = doc[RATIO] if RATIO in doc else doc
            # union of everyone ELSE's tracts
            others = {}
            for c in CONTIGS:
                iv = [x for o, t in truths.items() if o != s for x in t.get(c, [])]
                others[c] = merge(iv) if iv else []
            for seg in doc['segments']:
                c = seg['contig']
                if c not in CONTIGS or seg['posterior'] < MIN_POST or seg['n_private'] < MIN_SITES:
                    continue
                if seg['end'] - seg['start'] < MIN_BP:
                    continue
                cls = 'true positive' if overlaps(truths[s].get(c, []), seg['start'], seg['end']) \
                    else 'FALSE positive'
                tot[cls] += 1
                if overlaps(others[c], seg['start'], seg['end']):
                    counts[cls] += 1
                if cls == 'FALSE positive':
                    # a same-size region placed at random in callable territory
                    L = seg['end'] - seg['start']
                    cal = callable_r.get(c, [])
                    if cal:
                        lo = min(a for a, _ in cal)
                        hi = max(b for _, b in cal)
                        p = rnd.randint(lo, max(lo, hi - L))
                        null_tot += 1
                        null_hit += overlaps(others[c], p, p + L)
        for cls in ('true positive', 'FALSE positive'):
            if tot[cls]:
                extra = (f'{null_hit / null_tot * 100:11.1f}%'
                         if cls == 'FALSE positive' and null_tot else ' ' * 12)
                print(f'{grp:4s} {cls:16s} {tot[cls]:5d} {counts[cls] / tot[cls] * 100:15.1f}%{extra}')
    print('\n  If FALSE positives hit other people\'s tracts far above the random null, they are')
    print('  real archaic haplotypes hmmix missed in THIS individual -- and precision against')
    print('  hmmix is then measuring their sensitivity, not our specificity.')


if __name__ == '__main__':
    main()
