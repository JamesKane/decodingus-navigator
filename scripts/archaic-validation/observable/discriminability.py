"""How detectable is a single tract under each observable?

Contrast alone does not decide this. What decides it is how much evidence one tract carries, and the
two observables differ 30-fold there:

  private-variant density : ~1 informative variant per 36 kb tract (Poisson 0.6 vs 1.3)
  archaic-allele matching : ~30 diagnostic sites  (Binomial 30 x 13% vs 30 x 39.5%)

Same ~3x contrast, wildly different decidability. This computes, for each, the sensitivity
achievable at a fixed false-positive rate -- i.e. what fraction of real tracts a perfect classifier
could call while keeping background calls rare.
"""

import math


def binom_pmf(k, n, p):
    return math.comb(n, k) * p ** k * (1 - p) ** (n - k)


def pois_pmf(k, lam):
    return math.exp(-lam) * lam ** k / math.factorial(k)


def roc(bg, sig, kmax):
    """Sensitivity at the threshold where background false-positive rate first drops below 5%/1%."""
    out = {}
    for target in (0.05, 0.01):
        for t in range(kmax + 1):
            fp = sum(bg(k) for k in range(t, kmax + 1))
            if fp <= target:
                out[target] = (t, sum(sig(k) for k in range(t, kmax + 1)), fp)
                break
        else:
            out[target] = (None, 0.0, 0.0)
    return out


print('DETECTABILITY OF ONE TRACT (36 kb, the measured median)\n')

# --- private-variant density ---------------------------------------------------------------
lam_bg, lam_sig = 0.64, 1.30
r = roc(lambda k: pois_pmf(k, lam_bg), lambda k: pois_pmf(k, lam_sig), 30)
print(f'private-variant density   background Poisson({lam_bg})  tract Poisson({lam_sig})')
for tgt, (t, sens, fp) in r.items():
    print(f'   at <= {tgt:.0%} false positives: threshold k >= {t}, sensitivity {sens:6.1%}')

# --- archaic-allele matching ---------------------------------------------------------------
n, p_bg, p_sig = 30, 0.130, 0.395
r = roc(lambda k: binom_pmf(k, n, p_bg), lambda k: binom_pmf(k, n, p_sig), n)
print(f'\narchaic-allele matching   background Binom({n}, {p_bg})  tract Binom({n}, {p_sig})')
for tgt, (t, sens, fp) in r.items():
    print(f'   at <= {tgt:.0%} false positives: threshold k >= {t}, sensitivity {sens:6.1%}')

# What tract size does each need to become usable?
print('\nTRACT SIZE NEEDED FOR 80% SENSITIVITY AT 5% FALSE POSITIVES')
for label, kind in (('private density', 'pois'), ('allele matching', 'binom')):
    for kb in (10, 20, 36, 50, 100, 200, 500):
        scale = kb / 36
        if kind == 'pois':
            bg, sig = lam_bg * scale, lam_sig * scale
            rr = roc(lambda k: pois_pmf(k, bg), lambda k: pois_pmf(k, sig), 60)
        else:
            nn = max(int(round(n * scale)), 1)
            rr = roc(lambda k: binom_pmf(k, nn, p_bg), lambda k: binom_pmf(k, nn, p_sig), nn)
        if rr[0.05][1] >= 0.80:
            print(f'  {label:18s} {kb:4d} kb  -> sensitivity {rr[0.05][1]:.0%}')
            break
    else:
        print(f'  {label:18s} not reached by 500 kb')

print('\n  hmmix p10 tract is 7 kb and the median 31-36 kb, so a method needing hundreds of kb')
print('  cannot report tracts at the resolution the feature claims.')
