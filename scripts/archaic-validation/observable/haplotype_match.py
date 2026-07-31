"""Should the observable be archaic-ALLELE MATCHING rather than private-mutation density?

The current caller counts private variants per window: measured contrast inside real archaic tracts
is ~2x, and at ~20-120 private variants/Mb a median 31 kb tract carries about ONE informative
variant. That is why no parameter setting works.

But an introgressed tract is a haplotype inherited intact from an archaic ancestor, and we HOLD the
archaic genomes -- 2,031,406 diagnostic sites where archaics carry a derived allele. A real tract
should carry the archaic allele at a large fraction of the diagnostic sites it spans, while a
non-introgressed region carries it only at the background rate (~4%, per the design's own
measurement).

Diagnostic sites are the denominator here, not megabases, which controls for their uneven density
for free -- the thing a mutation-rate map was going to have to correct.
"""

import bisect
import collections
import json
import random

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
    # Diagnostic sites: position -> (derived base, class)
    diag = collections.defaultdict(dict)
    with open('classify.tsv') as f:
        next(f)
        for line in f:
            c, p, d, k = line.rstrip('\n').split('\t')
            diag[c][int(p)] = (d, int(k))

    # The subject's calls, so we can ask whether he carries the archaic allele.
    carries = collections.defaultdict(dict)
    for c in CONTIGS:
        for rec in json.load(open(f'HG00096.{c}.calls.json')):
            ref = rec['reference_allele']
            alt = rec['alternate_allele']
            dos = rec['dosage']
            if dos is None or dos < 0:
                continue
            alleles = set()
            if dos < 2:
                alleles.add(ref)
            if dos > 0:
                alleles.add(alt)
            carries[rec['contig']][rec['position']] = alleles

    callable_r = load_bed('callable.bed')
    truth = load_bed('truth_HG00096.chm13.bed', 1000)

    def in_regions(regions, p):
        i = bisect.bisect_right([s for s, _ in regions], p) - 1
        return i >= 0 and regions[i][1] > p

    stats = {'in': [0, 0], 'out': [0, 0]}
    for c in CONTIGS:
        cal = callable_r[c]
        tr = truth.get(c, [])
        for p, (derived, _k) in diag[c].items():
            if not in_regions(cal, p):
                continue
            # EVERY diagnostic site in callable territory is in the denominator. The caller emits
            # only variant records, so no record means hom-reference — i.e. NOT carrying the
            # archaic allele, given the panel orients derived as ALT. Conditioning on "has a call"
            # instead samples only sites where he already has a variant, which is why that version
            # reported an impossible ~80% carrying rate against a known 4.3% background.
            called = carries[c].get(p)
            bucket = 'in' if in_regions(tr, p) else 'out'
            stats[bucket][1] += 1
            if called is not None and derived in called:
                stats[bucket][0] += 1

    print('ARCHAIC-ALLELE CARRYING RATE at diagnostic sites the subject has a call for')
    for b, label in (('in', "inside hmmix's tracts"), ('out', 'elsewhere')):
        hit, tot = stats[b]
        rate = hit / tot * 100 if tot else 0
        print(f'  {label:26s} {hit:6d} / {tot:6d} = {rate:5.1f}%')
    ri = stats['in'][0] / max(stats['in'][1], 1)
    ro = stats['out'][0] / max(stats['out'][1], 1)
    print(f'  ENRICHMENT {ri / ro if ro else 0:.2f}x')
    print()
    print('  Compare: private-variant DENSITY gives 2.89x with ~1 informative variant per tract.')
    print('  Sites here are the denominator, so uneven diagnostic-site density cancels out.')

    # How much evidence does a typical tract actually carry under this observable?
    spans = [e - s for c in CONTIGS for s, e in truth.get(c, [])]
    spans.sort()
    med = spans[len(spans) // 2]
    dens = sum(len(diag[c]) for c in CONTIGS) / sum(
        e - s for c in CONTIGS for s, e in merge([(min(x for x, _ in callable_r[c]),
                                                   max(y for _, y in callable_r[c]))]))
    print(f'\n  median tract {med / 1000:.0f} kb; diagnostic sites ~{dens * 1e6:.0f}/Mb'
          f' -> ~{med * dens:.0f} informative sites per tract')
    print('  (the density model gets ~1 private variant for the same tract)')


if __name__ == '__main__':
    main()
