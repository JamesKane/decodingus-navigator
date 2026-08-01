"""Is the diagnostic-site panel less informative inside EAST ASIAN archaic tracts?

The caller transfers per-individual (30/30 above null, precision 41.9%) but inverts the population
ordering: the truth puts East Asian archaic extent at 1.22x Europe's, and we call 0.94x. Ruled out
already: background contamination (carrying rates 11.9% vs 12.2%, and both states scale together)
and tract length (median 29 kb in both; EAS simply have MORE tracts).

What remains is the observable itself. The model detects a tract by an elevated rate of carrying the
archaic allele at panel sites. If East Asian tracts carry those sites at a LOWER rate -- because the
panel's sites were ascertained on data that better represents the haplotypes introgressed into
Europeans -- then the same tract is less visible in an East Asian, and the caller under-calls exactly
where the truth says there is more.

This measures the in-tract and background carrying rates per population. The contrast, not the
absolute rate, is what the model separates on.
"""

import bisect
import collections
import json
import os
import subprocess
import sys

SC = os.path.dirname(os.path.abspath(__file__))
V2 = f'{SC}/v2'
DB = os.path.expanduser('~/.decodingus/navigator-rs.db')
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
        p = line.split()
        if len(p) >= 3:
            d[p[0]].append((int(p[1]), int(p[2])))
    return {c: merge(v, tol) for c, v in d.items()}


def sql(q):
    return subprocess.run(['sqlite3', DB, q], capture_output=True, text=True).stdout.strip()


def calls_for(sample):
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
            alleles = set()
            if d < 2:
                alleles.add(rec['reference_allele'])
            if d > 0:
                alleles.add(rec['alternate_allele'])
            out[(rec['contig'], rec['position'])] = alleles
    return out


def main():
    diag = collections.defaultdict(dict)
    with open(f'{SC}/classify.tsv') as f:
        next(f)
        for line in f:
            c, p, d, k = line.rstrip('\n').split('\t')
            diag[c][int(p)] = d

    callable_r = load_bed(f'{SC}/callable.bed')

    def in_reg(regs, p):
        i = bisect.bisect_right([s for s, _ in regs], p) - 1
        return i >= 0 and regs[i][1] > p

    print(f"{'pop':4s} {'sample':10s} {'in-tract':>9s} {'background':>11s} {'contrast':>9s}")
    agg = collections.defaultdict(lambda: [0, 0, 0, 0])
    for grp, listfile, n in (('EUR', 'eur60.txt', 8), ('EAS', 'eas30.txt', 8)):
        for s in [x.strip() for x in open(f'{SC}/{listfile}') if x.strip()][:n]:
            tpath = f'{V2}/truth_{s}.bed'
            if not os.path.exists(tpath):
                continue
            carries = calls_for(s)
            if carries is None:
                continue
            truth = load_bed(tpath, tol=1000)
            hits = [0, 0, 0, 0]  # in_carried, in_total, bg_carried, bg_total
            for c in CONTIGS:
                cal = callable_r.get(c, [])
                tr = truth.get(c, [])
                for p, derived in diag[c].items():
                    if not cal or not in_reg(cal, p):
                        continue
                    got = carries.get((c, p))
                    carried = got is not None and derived in got
                    if tr and in_reg(tr, p):
                        hits[1] += 1
                        hits[0] += carried
                    else:
                        hits[3] += 1
                        hits[2] += carried
            if hits[1] == 0 or hits[3] == 0:
                continue
            r_in = hits[0] / hits[1] * 100
            r_bg = hits[2] / hits[3] * 100
            print(f'{grp:4s} {s:10s} {r_in:8.1f}% {r_bg:10.1f}% {r_in / r_bg:8.2f}x')
            for i in range(4):
                agg[grp][i] += hits[i]

    print()
    for grp in ('EUR', 'EAS'):
        a = agg[grp]
        if a[1] and a[3]:
            r_in, r_bg = a[0] / a[1] * 100, a[2] / a[3] * 100
            print(f'{grp} POOLED  in-tract {r_in:.1f}%   background {r_bg:.1f}%   '
                  f'CONTRAST {r_in / r_bg:.2f}x')
    e, a = agg['EUR'], agg['EAS']
    if e[1] and a[1]:
        ce = (e[0] / e[1]) / (e[2] / e[3])
        ca = (a[0] / a[1]) / (a[2] / a[3])
        print(f'\ncontrast EAS/EUR = {ca / ce:.3f}')
        print('  < 1 means the panel is less informative inside East Asian tracts, which would')
        print('  explain under-calling exactly where the truth says there is MORE archaic DNA.')


if __name__ == '__main__':
    main()
