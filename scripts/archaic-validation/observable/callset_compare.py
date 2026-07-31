"""Same individual, same outgroup, same build: 1000G's call set versus ours.

The Tier B observable is "variants this person carries that no African carries". Its contrast
inside real archaic tracts is 2.89x for our calls -- too weak to separate from a background that
varies 5.3x. hmmix's documented example emission rates imply a background near 40/Mb and a ~10x
contrast, and ours is 124/Mb at 2.89x.

Our outgroup track was just verified complete against its source, so if the observable is diluted
the dilution must come from the variant calls. HG00096 is in the 1000G callset, so the same
quantity can be computed from their calls and ours and compared directly.
"""

import bisect
import collections


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


def stats(positions, callable_r, truth, label):
    a = sorted(positions)
    cnt = lambda s, e: bisect.bisect_left(a, e) - bisect.bisect_left(a, s)
    in_bp = sum(covered(callable_r, s, e) for s, e in truth)
    in_n = sum(cnt(s, e) for s, e in truth)
    cb = sum(e - s for s, e in callable_r)
    cn = sum(cnt(s, e) for s, e in callable_r)
    out_n, out_bp = cn - in_n, cb - in_bp
    enrich = (in_n / in_bp) / (out_n / out_bp)
    print(f'{label:26s} {len(a):7d} private   in-tract {in_n / (in_bp / 1e6):6.1f}/Mb   '
          f'background {out_n / (out_bp / 1e6):6.1f}/Mb   CONTRAST {enrich:5.2f}x')
    return enrich


def main():
    afr = set(int(x) for x in open('afr_sites.chr21.txt'))
    kgp_carried = [int(x) for x in open('kgp_carried.chr21.txt')]
    ours_carried, ours_private = [], []
    with open('private.chr21.tsv') as f:
        next(f)
        for line in f:
            p = line.split('\t')
            ours_private.append(int(p[1]))

    callable_r = load_bed('callable.bed')['chr21']
    truth = load_bed('truth_HG00096.chm13.bed', 1000)['chr21']

    kgp_private = [p for p in kgp_carried if p not in afr]

    print(f'HG00096 chr21, carried SNVs:  1000G {len(kgp_carried):7d}   ours {93643:7d}'
          f'   ratio {93643 / len(kgp_carried):.2f}x')
    print(f'African outgroup: {len(afr)} segregating sites (652 unrelated individuals)\n')

    print('observable                    count      density in tract / background      contrast')
    e_kgp = stats(kgp_private, callable_r, truth, "1000G's calls")
    e_ours = stats(ours_private, callable_r, truth, 'our calls')
    print()
    print(f'  contrast ratio 1000G/ours = {e_kgp / e_ours:.2f}x')
    print()
    print('  If 1000G\'s calls give a much higher contrast on the SAME person, the observable is')
    print('  diluted by our variant calling, not by the outgroup or the model -- and the fix is')
    print('  upstream of the HMM entirely.')

    # Where do the extra calls sit? If ours are mostly at positions 1000G did not call at all,
    # they are ours to explain.
    kgp_set = set(kgp_carried)
    extra = [p for p in ours_private if p not in kgp_set]
    shared = len(ours_private) - len(extra)
    print(f'\n  of our {len(ours_private)} private calls: {shared} also called by 1000G, '
          f'{len(extra)} ours alone ({len(extra) / max(len(ours_private), 1) * 100:.0f}%)')


if __name__ == '__main__':
    main()
