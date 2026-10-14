"""RSOS kinks of a_2^(1) at q and -q (book: prp-rsos-kinks).

At fixed spectral parameter x, q -> -q leaves every quantum dimension and every diagonal face weight
unchanged and reverses the off-diagonal ones (checked below against gam -> gam + pi with gam*u fixed).
The restricted transfer matrices then differ, but their spectra agree.
For contrast: gam -> gam + pi at fixed u = 3 i theta/2 pi changes x(theta), i.e. the coupling, and does
break the p' = p + 2 sectors; that is not a change of the sign of q.
"""
import numpy as np
from irf_a2 import A2RSOS, Fast

_W33 = A2RSOS.W33
def flipped(self, a, b, c, d, u, sqrt_branch=None):
    w = _W33(self, a, b, c, d, u, sqrt_branch); return -w if d != b else w

# the flip is q -> -q at fixed x
p, pp = 7, 9; m1 = A2RSOS(p, pp); lam = (pp - p)/p; m2 = A2RSOS(p, lam=lam + 1)
u = 0.31 + 0.27j; ratios = {}
for a in m1.H:
    for b in m1.H:
        for c in m1.H:
            for d in m1.H:
                w1 = m1.W33(a, b, c, d, u); w2 = m2.W33(a, b, c, d, u*lam/(lam + 1))
                if abs(w1) > 1e-12: ratios.setdefault('diagonal' if d == b else 'off-diagonal', set()).add(round((w2/w1).real, 8))
print("weights at gam+pi (same x) / weights at gam:", {k: sorted(v) for k, v in ratios.items()},
      "; quantum dimensions equal:", max(abs(m1.qdim(h) - m2.qdim(h)) for h in m1.H) < 1e-12)

rng = np.random.default_rng(11)
for (p, pp) in [(7, 9), (9, 11), (10, 13), (8, 11)]:
    m = A2RSOS(p, pp)
    for types in [(0, 1), (0, 0, 0), (0, 1, 0, 1), (0, 0, 1, 1)]:
        th = rng.normal(size=len(types))*1.5
        A2RSOS.W33 = _W33; T1 = Fast(m).transfer(types, th)
        A2RSOS.W33 = flipped; T2 = Fast(m).transfer(types, th)
        A2RSOS.W33 = _W33
        e1, e2 = np.linalg.eigvals(T1), np.linalg.eigvals(T2)
        d = max(np.min(np.abs(e2 - z)) for z in e1)
        lab = ''.join('3' if t == 0 else 'b' for t in types)
        print(f"W3({p},{pp}) {lab}: |T(q)-T(-q)| = {np.abs(T1 - T2).max():.1e}, spectral distance {d:.0e}, "
              f"max||s|-1| {np.abs(np.abs(e1) - 1).max():.1e} (q), {np.abs(np.abs(e2) - 1).max():.1e} (-q)")

print("contrast, gam -> gam + pi at fixed u (a change of coupling):")
for (p, pp) in [(7, 9), (9, 11)]:
    for lab, m in (("gam = pi(p'-p)/p", A2RSOS(p, pp)), ("gam = pi p'/p", A2RSOS(p, lam=pp/p))):
        f = Fast(m); out = []
        for types in [(0, 1), (0, 0, 0), (0, 1, 0, 1)]:
            w = max(np.abs(np.abs(np.linalg.eigvals(f.transfer(types, rng.normal(size=len(types))*1.5))) - 1).max() for _ in range(8))
            out.append(f"{''.join('3' if t == 0 else 'b' for t in types)}: {w:.1e}")
        print(f"  W3({p},{pp}) {lab}: " + "; ".join(out))
