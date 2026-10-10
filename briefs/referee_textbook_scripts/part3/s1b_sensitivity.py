"""Sensitivity: with a mixed (non-global) relabelling or a real basis of a degenerate pair, Dorey's rule / area rule fail."""
import itertools, numpy as np
from lie import Algebra
exec(open('s1_masses_dorey.py').read().split('def check')[0].split('# ---- e8 masses')[0])
exec('def bicolour_coxeter' + open('s1_masses_dorey.py').read().split('def bicolour_coxeter')[1].split('def check')[0])
Al = Algebra('e', 6); h = Al.h; nodes = sorted(Al.u)
pairs = sorted({tuple(sorted((a, Al.conj[a]))) for a in nodes if Al.conj[a] != a})
print('e6 conjugate pairs (Bourbaki nodes):', pairs)
w, orb = bicolour_coxeter(Al)
def test(u):
    bad = 0; consts = []
    for a, b, c in itertools.combinations_with_replacement(nodes, 3):
        val = sum(Al.n[j]*(Al.alpha[j] @ u[a])*(Al.alpha[j] @ u[b])*(Al.alpha[j] @ u[c]) for j in range(7))
        rule = any(np.allclose(x+y+z, 0, atol=1e-8) for x in orb[a] for y in orb[b] for z in orb[c])
        bad += (abs(val) > 1e-8) != rule
    return bad
u = dict(Al.u); a, b = pairs[0]; u[a], u[b] = Al.u[b], Al.u[a]
print('swap only the pair', pairs[0], '-> rule mismatches:', test(u))
u = dict(Al.u); a, b = pairs[0]
u[a] = (Al.u[a] + Al.u[b])/np.sqrt(2); u[b] = 1j*(Al.u[a] - Al.u[b])/np.sqrt(2)
print('real basis for pair', pairs[0], '-> rule mismatches:', test(u))
# area rule in a real orthonormal basis of the degenerate pairs (as literally prescribed at 09-lagrangian.qmd:120)
u = dict(Al.u)
for (a, b) in pairs:
    va, vb = Al.u[a], Al.u[b]
    u[a] = ((va + vb)/np.sqrt(2)).real; u[b] = ((1j*(va - vb))/np.sqrt(2)).real
    assert np.allclose(u[a] @ u[a], 1) and np.allclose(u[b] @ u[b], 1) and abs(u[a] @ u[b]) < 1e-12
ratios = []
for a, b, c in itertools.combinations_with_replacement(nodes, 3):
    val = sum(Al.n[j]*(Al.alpha[j] @ u[a])*(Al.alpha[j] @ u[b])*(Al.alpha[j] @ u[c]) for j in range(7))
    ma, mb, mc = Al.mass[a], Al.mass[b], Al.mass[c]; s = (ma+mb+mc)/2; ar2 = s*(s-ma)*(s-mb)*(s-mc)
    if abs(val) > 1e-8 and ar2 > 1e-12: ratios.append(abs(val)/np.sqrt(ar2)/(4/np.sqrt(h)))
print('real orthonormal basis of e6: |C_abc| / ((4/sqrt h) area) takes values', sorted({round(r, 4) for r in ratios}))
