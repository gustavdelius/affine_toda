"""Charge-conserving intertwiners at roots of unity: the loop system on one sector (numerical SVD, complex parameters)."""
import numpy as np, itertools, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import integration_uniqueness as cu
def loop_sv(h, q, sign, conj, m=0.43-0.21j, y=0.8+0.3j, t=None):
    t = t or tuple([1.0]*h)
    W = cu.Weyl(h, q, m, t, sign); P = cu.particle(h, q, y, conj)
    Q = tuple([1] + [0]*(h - 1))
    LA = [np.array(L, complex) for L in cu.loops(P, W, Q, True)]; LB = [np.array(L, complex) for L in cu.loops(P, W, Q, False)]
    I = np.eye(h)
    M = np.vstack([np.kron(B, I) - np.kron(I, A.T) for A, B in zip(LA, LB)])     # vec(T) row-major: (B T - T A)
    s = np.linalg.svd(M, compute_uv=False)
    return s[-1]/s[0], s[-2]/s[0]
for name, q in (('q^5=1 (omega=0.6)', np.exp(2j*np.pi/5)), ('q=i (omega=1/2)', 1j), ('omega=1.37 (q^200=1)', -np.exp(-1j*np.pi*1.37)),
                ('q^3=1', np.exp(2j*np.pi/3)), ('q^4=1', 1j*np.exp(0j)), ('q^6=1', np.exp(1j*np.pi/3))):
    res = []
    for h in (3, 4):
        for sign in (-1, 1):
            for conj in (False, True):
                a, b = loop_sv(h, q, sign, conj)
                res.append((a, b))
    print(f'{name:22s}: smallest relative sv {max(r[0] for r in res):.1e}, next {min(r[1] for r in res):.3f} (h = 3, 4, both Borels, V and Vbar)')
