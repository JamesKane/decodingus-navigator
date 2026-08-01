"""Are our "false positives" archaic-looking, or noise? And does that differ by population?

The blocker is that reported extent orders the populations backwards (0.937x where the truth says
1.217x), driven by precision differing between them (32.2% EUR vs 41.9% EAS) at identical
thresholds and identical in-tract contrast. Two very different explanations fit that:

  (a) OUR FAULT -- we emit more spurious calls in Europeans. Those segments should look like
      background: carrying rate near 13%.

  (b) THE TRUTH'S FAULT -- hmmix is less sensitive in Europeans, so real tracts we find are
      scored as false positives. Those segments should look like true positives: carrying rate
      near 40%.

The archaic-allele carrying rate inside each segment class separates the two, and it is measured
independently of hmmix -- it uses only our calls and the diagnostic panel, so it does not assume
the reference callset is right.
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


def main():
    diag = collections.defaultdict(dict)
    with open(f'{SC}/classify.tsv') as f:
        next(f)
        for line in f:
            c, p, d, _k = line.rstrip('\n').split('\t')
            diag[c][int(p)] = d
    diag_sorted = {c: sorted(v) for c, v in diag.items()}
    callable_r = load_bed(f'{SC}/callable.bed')

    def rate(contig, s, e, carries):
        keys = diag_sorted.get(contig, [])
        lo = bisect.bisect_left(keys, s)
        hi = bisect.bisect_left(keys, e)
        tot = hit = 0
        for p in keys[lo:hi]:
            tot += 1
            got = carries.get((contig, p))
            if got and diag[contig][p] in got:
                hit += 1
        return hit, tot

    print(f"{'pop':4s} {'class':18s} {'segments':>9s} {'Mb':>7s} {'carrying rate':>14s}")
    agg = collections.defaultdict(lambda: [0, 0, 0, 0.0])
    for grp, listfile, source, n in (('EUR', 'eur60.txt', f'{V2}/sweep', 12),
                                     ('EAS', 'eas30.txt', f'{V2}/sweep_eas', 12)):
        for s in [x.strip() for x in open(f'{SC}/{listfile}') if x.strip()][:n]:
            f = f'{source}/{s}.json'
            tp = f'{V2}/truth_{s}.bed'
            if not (os.path.exists(f) and os.path.exists(tp)):
                continue
            doc = json.load(open(f))
            doc = doc[RATIO] if RATIO in doc else doc
            truth = load_bed(tp, tol=1000)
            carries = carries_for(s)
            if carries is None:
                continue
            for seg in doc['segments']:
                c = seg['contig']
                if c not in CONTIGS or seg['posterior'] < MIN_POST or seg['n_private'] < MIN_SITES:
                    continue
                if seg['end'] - seg['start'] < MIN_BP:
                    continue
                cls = 'true positive' if overlaps(truth.get(c, []), seg['start'], seg['end']) \
                    else 'FALSE positive'
                hit, tot = rate(c, seg['start'], seg['end'], carries)
                a = agg[(grp, cls)]
                a[0] += 1
                a[1] += hit
                a[2] += tot
                a[3] += (seg['end'] - seg['start']) / 1e6
            # background: the callable territory outside both our calls and the truth
            ours = merge([(x['start'], x['end']) for x in doc['segments']
                          if x['contig'] == 'chr21' and x['posterior'] >= MIN_POST
                          and x['n_private'] >= MIN_SITES] or [(0, 0)])
            for a0, b0 in callable_r.get('chr21', [])[:400]:
                if overlaps(ours, a0, b0) or overlaps(truth.get('chr21', []), a0, b0):
                    continue
                hit, tot = rate('chr21', a0, b0, carries)
                a = agg[(grp, 'background')]
                a[1] += hit
                a[2] += tot

    for grp in ('EUR', 'EAS'):
        for cls in ('true positive', 'FALSE positive', 'background'):
            a = agg[(grp, cls)]
            if not a[2]:
                continue
            print(f'{grp:4s} {cls:18s} {a[0]:9d} {a[3]:7.2f} {a[1] / a[2] * 100:13.1f}%')
    print()
    for grp in ('EUR', 'EAS'):
        tp = agg[(grp, 'true positive')]
        fp = agg[(grp, 'FALSE positive')]
        bg = agg[(grp, 'background')]
        if tp[2] and fp[2] and bg[2]:
            r_tp, r_fp, r_bg = tp[1] / tp[2], fp[1] / fp[2], bg[1] / bg[2]
            # 1.0 = our false positives look exactly like real tracts; 0.0 = like background.
            pos = (r_fp - r_bg) / (r_tp - r_bg) if r_tp > r_bg else 0
            print(f'{grp}: false positives sit {pos * 100:.0f}% of the way from background to '
                  f'true positive')
    print('\n  Near 0% => they are noise and the precision gap is ours.')
    print('  Near 100% => they are archaic-looking and hmmix simply did not call them, which makes')
    print('  "precision" against this reference a measure of the reference as much as of us.')


if __name__ == '__main__':
    main()
