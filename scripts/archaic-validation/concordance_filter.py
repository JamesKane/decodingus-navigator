"""Filter called segments by archaic-genome concordance — fitted on the three Neanderthals,
validated on Denisova, which the filter never sees.

Why: the reported extent orders the populations backwards, and the decomposition says why. The
true-positive component reproduces the truth ordering almost exactly (1.219 against 1.217); the
false positives, which are about twice the true positives, run the other way (0.877). The headline
number is dominated by noise whose population ordering is inverted. So the fix is precision.

The arbiter discriminates -- false positives score 81.3 %/72.9 % against true positives' ~93.6 % --
and it uses information the caller does not: which archaic genome carries what. Using it as a filter
is therefore a real gain rather than re-tuning something already fitted.

The cost is that filtering on the arbiter spends it as an independent referee. Holding out Denisova
keeps one: the filter sees only Altai, Vindija and Chagyrskaya, so Denisova concordance remains an
untouched check on whether the kept segments are genuinely archaic rather than merely
Neanderthal-shaped.
"""

import bisect
import collections
import json
import os
import subprocess

SC = os.path.dirname(os.path.abspath(__file__))
V2 = f'{SC}/v2'
DB = os.path.expanduser('~/.decodingus/navigator-rs.db')
CONTIGS = ('chr21', 'chr22')
RATIO, MIN_POST, MIN_SITES, MIN_BP = '4.5', 0.98, 16, 5_000
NEANDERTHALS = ['AltaiNeanderthal', 'Vindija33.19', 'Chagyrskaya8']
HELD_OUT = 'Denisova3'


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


def sql(q):
    return subprocess.run(['sqlite3', DB, q], capture_output=True, text=True).stdout.strip()


def carries_for(sample):
    aln = sql(f"SELECT a.id FROM biosample b JOIN sequence_run r ON r.biosample_guid=b.guid "
              f"JOIN alignment a ON a.sequence_run_id=r.id WHERE b.donor_identifier='{sample}' "
              f"AND a.reference_build='chm13v2.0' LIMIT 1;")
    if not aln:
        return None
    out = {}
    for c in CONTIGS:
        payload = sql(f"SELECT payload FROM analysis_artifact WHERE alignment_id={aln} "
                      f"AND kind='diploid_denovo:{c}';")
        if not payload:
            return None
        for rec in json.loads(payload):
            d = rec['dosage']
            if d is None or d < 0:
                continue
            al = set()
            if d < 2:
                al.add(rec['reference_allele'])
            if d > 0:
                al.add(rec['alternate_allele'])
            out[(rec['contig'], rec['position'])] = al
    return out


def load_panel():
    by_contig = collections.defaultdict(list)
    with open(f'{SC}/panel.tsv') as f:
        hdr = next(f).rstrip('\n').split('\t')
        ni = [hdr.index(g) for g in NEANDERTHALS]
        di = hdr.index(HELD_OUT)
        for line in f:
            p = line.rstrip('\n').split('\t')
            by_contig[p[0]].append((int(p[1]), p[2], [p[i] for i in ni], p[di]))
    return {c: sorted(v) for c, v in by_contig.items()}


def overlaps(regs, s, e):
    for a, b in regs:
        if min(e, b) > max(s, a):
            return True
    return False


def scores(panel, contig, s, e, carries):
    """(neanderthal concordance, denisova concordance) over the segment; None when too few sites."""
    v = panel.get(contig, [])
    keys = [x[0] for x in v]
    lo, hi = bisect.bisect_left(keys, s), bisect.bisect_left(keys, e)
    nh = [0] * len(NEANDERTHALS)
    nd = [0] * len(NEANDERTHALS)
    dh = dd = 0
    for pos, derived, ncalls, dcall in v[lo:hi]:
        got = carries.get((contig, pos))
        has = bool(got) and derived in got
        for i, c in enumerate(ncalls):
            if c == 'D':
                nd[i] += 1
                nh[i] += has
        if dcall == 'D':
            dd += 1
            dh += has
    best = max((h / d for h, d in zip(nh, nd) if d >= 3), default=None)
    den = (dh / dd) if dd >= 3 else None
    return best, den


def main():
    panel = load_panel()
    rows = []
    for grp, listfile, source, n in (('EUR', 'eur60.txt', f'{V2}/sweep', 20),
                                     ('EAS', 'eas30.txt', f'{V2}/sweep_eas', 20)):
        for s in [x.strip() for x in open(f'{SC}/{listfile}') if x.strip()][:n]:
            f_ = f'{source}/{s}.json'
            tp = f'{V2}/truth_{s}.bed'
            if not (os.path.exists(f_) and os.path.exists(tp)):
                continue
            carries = carries_for(s)
            if carries is None:
                continue
            truth = load_bed(tp, tol=1000)
            doc = json.load(open(f_))
            doc = doc[RATIO] if RATIO in doc else doc
            for seg in doc['segments']:
                c = seg['contig']
                if c not in CONTIGS or seg['posterior'] < MIN_POST or seg['n_private'] < MIN_SITES:
                    continue
                if seg['end'] - seg['start'] < MIN_BP:
                    continue
                nea, den = scores(panel, c, seg['start'], seg['end'], carries)
                rows.append({
                    'grp': grp, 'sample': s, 'contig': c, 's': seg['start'], 'e': seg['end'],
                    'mb': (seg['end'] - seg['start']) / 1e6,
                    'tp': overlaps(truth.get(c, []), seg['start'], seg['end']),
                    'nea': nea, 'den': den,
                })
    print(f'segments scored: {len(rows)}  '
          f'(with a Neanderthal score: {sum(1 for r in rows if r["nea"] is not None)})\n')

    truth_mb = {}
    for grp, listfile in (('EUR', 'eur60.txt'), ('EAS', 'eas30.txt')):
        tot = n = 0
        for s in [x.strip() for x in open(f'{SC}/{listfile}') if x.strip()][:20]:
            p = f'{V2}/truth_{s}.bed'
            if os.path.exists(p):
                t = load_bed(p, tol=1000)
                tot += sum(e - a for v in t.values() for a, e in v) / 1e6
                n += 1
        truth_mb[grp] = tot / max(n, 1)
    print(f'truth Mb/person: EUR {truth_mb["EUR"]:.3f}  EAS {truth_mb["EAS"]:.3f}  '
          f'-> ordering to reproduce {truth_mb["EAS"] / truth_mb["EUR"]:.3f}\n')

    npeople = {g: len({r['sample'] for r in rows if r['grp'] == g}) for g in ('EUR', 'EAS')}
    print(f"{'threshold':>9s} {'kept':>6s} {'prec':>6s} {'EUR Mb':>7s} {'EAS Mb':>7s} "
          f"{'ordering':>9s} {'Denisova kept/dropped':>22s}")
    for thr in (0.0, 0.70, 0.80, 0.85, 0.90, 0.95, 1.0):
        kept = [r for r in rows if r['nea'] is not None and r['nea'] >= thr]
        dropped = [r for r in rows if r['nea'] is not None and r['nea'] < thr]
        if not kept:
            continue
        prec = sum(1 for r in kept if r['tp']) / len(kept) * 100
        eur = sum(r['mb'] for r in kept if r['grp'] == 'EUR') / npeople['EUR']
        eas = sum(r['mb'] for r in kept if r['grp'] == 'EAS') / npeople['EAS']
        dk = [r['den'] for r in kept if r['den'] is not None]
        dd = [r['den'] for r in dropped if r['den'] is not None]
        dtxt = (f'{sum(dk) / len(dk) * 100:8.1f}% / '
                f'{sum(dd) / len(dd) * 100:6.1f}%') if dk and dd else f'{sum(dk) / len(dk) * 100:8.1f}% /      -'
        print(f'{thr:9.2f} {len(kept):6d} {prec:5.1f}% {eur:7.3f} {eas:7.3f} '
              f'{eas / eur if eur else 0:9.3f} {dtxt:>22s}')
    print('\n  ordering: 1.217 is the target. Denisova is the HELD-OUT check -- it never enters the')
    print('  filter, so kept segments scoring higher on it than dropped ones means the filter is')
    print('  selecting genuinely archaic sequence, not just Neanderthal-shaped noise.')


if __name__ == '__main__':
    main()
