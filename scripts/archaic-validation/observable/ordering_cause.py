"""Does the population difference in BACKGROUND archaic carriage explain the ordering inversion?

The blocker: reported extent orders the populations backwards (we call 0.937x where the truth says
1.217x), because the false-positive load differs (precision 32.2% EUR vs 41.9% EAS) at identical
thresholds and identical in-tract contrast.

The arbiter turned up a candidate: outside tracts, Europeans match archaic genomes at 59.0% against
East Asians' 45.5%. If Europeans carry archaic-derived alleles more often in non-introgressed
sequence, then a fixed threshold turns more European background into calls -- more false positives,
more spurious extent, and an inflated European total.

The model's own background parameter should absorb this: p_background is estimated per individual
from the genome-wide carrying rate. So the question is whether that estimate actually tracks the
difference, or whether the two rates measure different things.

Measured here, per population:
  * the carrying rate the MODEL estimates (all diagnostic sites, what p_background sees)
  * the carrying rate the ARBITER sees outside tracts (per-genome concordance)
  * the resulting false-positive Mb per person
If the model's estimate does NOT differ while the arbiter's does, p_background is blind to exactly
the thing driving the inversion, and that is the defect to fix.
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
GENOMES = ['AltaiNeanderthal', 'Vindija33.19', 'Chagyrskaya8', 'Denisova3']


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


def overlaps(regs, s, e):
    for a, b in regs:
        if min(e, b) > max(s, a):
            return True
    return False


def in_reg(regs, p):
    i = bisect.bisect_right([s for s, _ in regs], p) - 1
    return i >= 0 and regs[i][1] > p


def main():
    # every diagnostic site (what the MODEL's p_background is estimated over)
    diag = collections.defaultdict(dict)
    with open(f'{SC}/classify.tsv') as f:
        next(f)
        for line in f:
            c, p, d, _ = line.rstrip('\n').split('\t')
            diag[c][int(p)] = d
    callable_r = load_bed(f'{SC}/callable.bed')

    print(f"{'pop':4s} {'model p_bg (all sites)':>24s} {'model p_bg OUTSIDE tracts':>27s} "
          f"{'FP Mb/person':>13s}")
    out = {}
    for grp, listfile, source, n in (('EUR', 'eur60.txt', f'{V2}/sweep', 15),
                                     ('EAS', 'eas30.txt', f'{V2}/sweep_eas', 15)):
        all_hit = all_tot = out_hit = out_tot = 0
        fp_mb = 0.0
        people = 0
        for s in [x.strip() for x in open(f'{SC}/{listfile}') if x.strip()][:n]:
            f_ = f'{source}/{s}.json'
            tp = f'{V2}/truth_{s}.bed'
            if not (os.path.exists(f_) and os.path.exists(tp)):
                continue
            carries = carries_for(s)
            if carries is None:
                continue
            people += 1
            truth = load_bed(tp, tol=1000)
            for c in CONTIGS:
                cal = callable_r.get(c, [])
                tr = truth.get(c, [])
                for p, derived in diag[c].items():
                    if not cal or not in_reg(cal, p):
                        continue
                    got = carries.get((c, p))
                    hit = bool(got) and derived in got
                    all_tot += 1
                    all_hit += hit
                    if not (tr and in_reg(tr, p)):
                        out_tot += 1
                        out_hit += hit
            doc = json.load(open(f_))
            doc = doc[RATIO] if RATIO in doc else doc
            for seg in doc['segments']:
                c = seg['contig']
                if c not in CONTIGS or seg['posterior'] < MIN_POST or seg['n_private'] < MIN_SITES:
                    continue
                if seg['end'] - seg['start'] < MIN_BP:
                    continue
                if not overlaps(truth.get(c, []), seg['start'], seg['end']):
                    fp_mb += (seg['end'] - seg['start']) / 1e6
        a = all_hit / all_tot * 100
        o = out_hit / out_tot * 100
        print(f'{grp:4s} {a:23.2f}% {o:26.2f}% {fp_mb / people:13.3f}')
        out[grp] = (a, o, fp_mb / people)

    print()
    e, a = out['EUR'], out['EAS']
    print(f'  model p_background   EUR {e[0]:.2f}% vs EAS {a[0]:.2f}%   ratio {a[0] / e[0]:.3f}')
    print(f'  outside tracts       EUR {e[1]:.2f}% vs EAS {a[1]:.2f}%   ratio {a[1] / e[1]:.3f}')
    print(f'  false-positive Mb    EUR {e[2]:.3f}  vs EAS {a[2]:.3f}    ratio {a[2] / e[2]:.3f}')
    print()
    print('  If the model estimate is flat across populations while the FP load is not, then')
    print('  p_background is blind to whatever drives the extra European calls, and no amount of')
    print('  per-individual estimation of it will fix the ordering.')


if __name__ == '__main__':
    main()
