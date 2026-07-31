"""An arbiter for Tier B calls that does not ask another caller's opinion.

Precision has been measured against hmmix, but a call absent from hmmix is not necessarily wrong:
their callset is incomplete by an unknown amount (their own tracts are enriched just 1.84x for their
own archaic SNPs). Two attempts to settle whether our extra calls are real both failed --
carrying rate is circular (the caller selects on it), and shared-haplotype overlap is saturated
(hmmix tracts cover 67% of callable territory, so the null is ~60%).

This asks the archaic genomes instead. A real introgressed tract was inherited from ONE archaic
individual, so the derived alleles the subject carries inside it should concentrate on the genomes
that share that haplotype -- and in particular should be present at sites where ONE genome is
derived and another is positively called ancestral. Those discordant sites are the informative ones:
a site where all four archaics are derived says nothing about which haplotype was inherited.

Crucially, the caller never sees per-genome calls -- it reads only a derived base and a lineage
class from `ArchaicClassify` -- so this is evidence it cannot have fitted to.

Reports, for true positives, false positives and background:
  * DISCORDANT-SITE concordance -- of the sites where the subject carries the derived allele and the
    archaic genomes DISAGREE, what fraction match the best-matching genome. Random carriage gives
    the base rate; an inherited haplotype gives much more.
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


def load_panel():
    """Sites where the archaic genomes DISAGREE — the only ones that identify a haplotype."""
    by_contig = collections.defaultdict(list)
    with open(f'{SC}/panel.tsv') as f:
        hdr = next(f).rstrip('\n').split('\t')
        gi = [hdr.index(g) for g in GENOMES]
        for line in f:
            p = line.rstrip('\n').split('\t')
            calls = [p[i] for i in gi]
            if 'D' in calls and 'A' in calls:  # informative: some derived, some POSITIVELY ancestral
                by_contig[p[0]].append((int(p[1]), p[2], calls))
    return {c: sorted(v) for c, v in by_contig.items()}


def overlaps(regs, s, e):
    for a, b in regs:
        if min(e, b) > max(s, a):
            return True
    return False


def concordance(panel, contig, s, e, carries):
    """For each archaic genome g: of the sites where g is DERIVED, how many does the subject carry?

    Conditioning on the GENOME, not on the subject. An earlier version conditioned on the subject
    already carrying the derived allele, which is vacuous: at a discordant site at least one genome
    is derived by construction, so the best-matching genome scored ~100 % everywhere including
    background, and the statistic separated nothing.

    Read this way the quantity is the subject's sensitivity to a specific archaic haplotype:
    background sits at the genome-wide carrying rate (~13 %), while a tract inherited from that
    lineage should be far higher.
    """
    v = panel.get(contig, [])
    keys = [x[0] for x in v]
    lo, hi = bisect.bisect_left(keys, s), bisect.bisect_left(keys, e)
    hits = [0] * len(GENOMES)
    dens = [0] * len(GENOMES)
    for pos, derived, calls in v[lo:hi]:
        got = carries.get((contig, pos))
        subject_has = bool(got) and derived in got
        for i, c in enumerate(calls):
            if c != 'D':
                continue
            dens[i] += 1
            if subject_has:
                hits[i] += 1
    # Best genome by rate, requiring a few sites so a 1/1 does not win.
    best_rate, best_hits, best_den = 0.0, 0, 0
    for h, d in zip(hits, dens):
        if d >= 3 and (h / d) > best_rate:
            best_rate, best_hits, best_den = h / d, h, d
    return best_hits, best_den


def main():
    panel = load_panel()
    n_inf = sum(len(v) for v in panel.values())
    print(f'informative (discordant) panel sites on chr21+22: {n_inf}\n')

    callable_r = load_bed(f'{SC}/callable.bed')
    agg = collections.defaultdict(lambda: [0, 0, 0])  # best_hits, total, n_segments

    for grp, listfile, source, n in (('EUR', 'eur60.txt', f'{V2}/sweep', 15),
                                     ('EAS', 'eas30.txt', f'{V2}/sweep_eas', 15)):
        for s in [x.strip() for x in open(f'{SC}/{listfile}') if x.strip()][:n]:
            f = f'{source}/{s}.json'
            tp = f'{V2}/truth_{s}.bed'
            if not (os.path.exists(f) and os.path.exists(tp)):
                continue
            carries = carries_for(s)
            if carries is None:
                continue
            doc = json.load(open(f))
            doc = doc[RATIO] if RATIO in doc else doc
            truth = load_bed(tp, tol=1000)
            called = collections.defaultdict(list)
            for seg in doc['segments']:
                c = seg['contig']
                if c not in CONTIGS or seg['posterior'] < MIN_POST or seg['n_private'] < MIN_SITES:
                    continue
                if seg['end'] - seg['start'] < MIN_BP:
                    continue
                called[c].append((seg['start'], seg['end']))
                cls = 'true positive' if overlaps(truth.get(c, []), seg['start'], seg['end']) \
                    else 'FALSE positive'
                h, d = concordance(panel, c, seg['start'], seg['end'], carries)
                if d >= 3:
                    a = agg[(grp, cls)]
                    a[0] += h
                    a[1] += d
                    a[2] += 1
            # background: callable windows we did NOT call and hmmix did not either
            for c in CONTIGS:
                cal = callable_r.get(c, [])
                mine = merge(called[c]) if called[c] else []
                for a0, b0 in cal[:600]:
                    if overlaps(mine, a0, b0) or overlaps(truth.get(c, []), a0, b0):
                        continue
                    h, d = concordance(panel, c, a0, b0, carries)
                    if d >= 3:
                        a = agg[(grp, 'background')]
                        a[0] += h
                        a[1] += d
                        a[2] += 1

    print(f"{'pop':4s} {'class':16s} {'regions':>8s} {'sites':>7s} {'best-genome concordance':>24s}")
    for grp in ('EUR', 'EAS'):
        for cls in ('true positive', 'FALSE positive', 'background'):
            a = agg[(grp, cls)]
            if a[1]:
                print(f'{grp:4s} {cls:16s} {a[2]:8d} {a[1]:7d} {a[0] / a[1] * 100:23.1f}%')
    print()
    for grp in ('EUR', 'EAS'):
        tp, fp, bg = agg[(grp, 'true positive')], agg[(grp, 'FALSE positive')], agg[(grp, 'background')]
        if tp[1] and fp[1] and bg[1]:
            r = lambda a: a[0] / a[1]
            span = r(tp) - r(bg)
            pos = (r(fp) - r(bg)) / span if span > 0 else 0
            print(f'{grp}: false positives sit {pos * 100:.0f}% of the way from background to true '
                  f'positive on ARCHAIC-GENOME concordance')
    print('\n  This statistic is not one the caller optimises: it reads per-genome calls, which the')
    print('  caller never sees. Near 100% means our extra calls carry the same haplotype signature')
    print('  as confirmed tracts, and precision against hmmix understates us. Near 0% means they')
    print('  are noise and the precision figure is fair.')


if __name__ == '__main__':
    main()
