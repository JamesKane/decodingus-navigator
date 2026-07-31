"""Is my lifted truth actually where hmmix says it is?

Everything downstream rests on it: I concluded the segment caller's locations are below chance, and
a shipped feature was gated off on that basis. But measuring private-variant density inside those
tracts gives only ~2x enrichment even in 1000G's OWN calls -- and hmmix found these tracts BY that
density, so it should be far higher. Either their method is weaker than advertised, or my lift put
the tracts in the wrong place.

This checks the truth in native hg38, with no lifting anywhere, using only hmmix's own two files:
their segment calls and their DAV (derived-in-archaic) SNP list. If their segments are enriched for
their own archaic SNPs in hg38, the callset is self-consistent and any failure is in my lift.
"""

import bisect
import collections
import os
import random

SEG = os.path.expanduser('~/.decodingus/ancestry-build/tmp/hmmix_segments_chr21_22.tsv')
SNPS = os.path.expanduser('~/.decodingus/ancestry-build/raw/hmmix/hg38_1000g_SNPS.txt')
SAMPLE = 'HG00096'
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


def main():
    # hmmix's own segments for this sample, in hg38, unioned across haplotypes.
    seg = collections.defaultdict(list)
    seg_all = collections.defaultdict(list)
    with open(SEG) as f:
        next(f)
        for line in f:
            p = line.rstrip('\n').split('\t')
            if p[4] not in CONTIGS:
                continue
            seg_all[p[4]].append((int(p[5]), int(p[6])))
            if p[0] == SAMPLE:
                seg[p[4]].append((int(p[5]), int(p[6])))
    S = {c: merge(v) for c, v in seg.items()}

    # hmmix's own archaic SNP list, in hg38.
    dav = collections.defaultdict(list)
    with open(SNPS) as f:
        next(f)
        for line in f:
            p = line.split('\t')
            if p[0] in CONTIGS:
                dav[p[0]].append(int(p[1]))
    for c in dav:
        dav[c].sort()

    print(f'hmmix segments for {SAMPLE} (hg38, unioned): '
          f'{sum(len(v) for v in S.values())} tracts, '
          f'{sum(e - s for v in S.values() for s, e in v) / 1e6:.3f} Mb')
    print(f'hmmix DAV archaic SNPs on chr21+22: {sum(len(v) for v in dav.values())}')
    print()

    # Span over which everyone's tracts fall -- the region the sampling actually covers.
    bounds = {c: (min(s for s, _ in v), max(e for _, e in v)) for c, v in seg_all.items()}

    tot_in = tot_bp = 0
    tot_out = tot_out_bp = 0
    for c in CONTIGS:
        a = dav[c]
        cnt = lambda s, e: bisect.bisect_left(a, e) - bisect.bisect_left(a, s)
        lo, hi = bounds[c]
        in_n = sum(cnt(s, e) for s, e in S[c])
        in_bp = sum(e - s for s, e in S[c])
        all_n = cnt(lo, hi)
        all_bp = hi - lo
        tot_in += in_n
        tot_bp += in_bp
        tot_out += all_n - in_n
        tot_out_bp += all_bp - in_bp

    d_in = tot_in / (tot_bp / 1e6)
    d_out = tot_out / (tot_out_bp / 1e6)
    print('IN NATIVE hg38, NO LIFTING:')
    print(f'  archaic SNPs inside {SAMPLE}\'s own hmmix tracts : {d_in:7.1f} /Mb  ({tot_in} in '
          f'{tot_bp / 1e6:.2f} Mb)')
    print(f'  archaic SNPs elsewhere                        : {d_out:7.1f} /Mb')
    print(f'  ENRICHMENT                                    : {d_in / d_out:6.2f}x')
    print()
    print('  hmmix called these tracts from archaic-variant density, so a strong enrichment here')
    print('  means their callset is internally consistent and my LIFT is what to distrust.')
    print('  A weak enrichment here means the truth was never as sharp as assumed.')

    # Null: same tract lengths placed at random in the same span.
    rnd = random.Random(0)
    vals = []
    for _ in range(200):
        n = bp = 0
        for c in CONTIGS:
            a = dav[c]
            lo, hi = bounds[c]
            for s, e in S[c]:
                L = e - s
                p = rnd.randint(lo, max(lo, hi - L))
                n += bisect.bisect_left(a, p + L) - bisect.bisect_left(a, p)
                bp += L
        vals.append(n / (bp / 1e6) / d_out)
    vals.sort()
    print(f'  null (same tracts placed at random): mean {sum(vals) / len(vals):.2f}x  '
          f'p95 {vals[189]:.2f}x')


if __name__ == '__main__':
    main()
