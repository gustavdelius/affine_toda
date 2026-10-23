"""Integration check, ch. 43 item 1: higher sine-Gordon breathers depend on the coupling; physical threshold of the
lowest a_r breather.  Own code; uses only the book helpers sg_tmatrix, invariant_vector (ch. 37 matrices)."""
import sys, numpy as np
sys.path.insert(0, 'code')
from boundary_common import q_phys, prop_residual
from defect_common import sg_tmatrix, invariant_vector

TD = 0.25
def Tcl(D):  # eq-shg-linear-transmission
    return -1j*np.sinh(D/2 - 1j*np.pi/4)/np.sinh(D/2 + 1j*np.pi/4)

def breather_n(th, lam, n, n0=0):
    """B_n = singlet of V(th - i u_n/2) (x) V(th + i u_n/2), u_n = pi - n pi/lam, crossing eq-sg-transmission-matrix,
    right constituent first."""
    qq = q_phys(lam); u = np.pi - n*np.pi/lam
    t1, t2 = th - 1j*u/2, th + 1j*u/2
    rel, w = invariant_vector(np.exp(lam*t1), np.exp(lam*t2), qq)
    T1, V, W = sg_tmatrix(t1, lam, qq, TD, 6)
    T2 = sg_tmatrix(t2, lam, qq, TD, 6)[0]
    L = np.kron(T1, np.eye(2)) @ np.kron(np.eye(2), T2)
    e = np.zeros(W['dim']); e[list(W['n']).index(n0)] = 1
    res, c = prop_residual(L @ np.kron(w, e), np.kron(e, w))
    return c, max(rel, res)

th = TD + 0.13
out = {}
for lam in (2.6, 3.1):
    for n in (1, 2):
        c, fit = breather_n(th, lam, n)
        prod = np.prod([Tcl(th - TD + 1j*np.pi*(n + 1 - 2*k)/(2*lam)) for k in range(1, n + 1)])
        out[(lam, n)] = c
        print(f'lam={lam} B_{n}: bootstrap {c:.9f}  product formula {prod:.9f}  |diff| {abs(c-prod):.1e}  fit {fit:.0e}')
print('B_1 coupling dependence', abs(out[(2.6, 1)] - out[(3.1, 1)]))
print('B_2 coupling dependence', abs(out[(2.6, 2)] - out[(3.1, 2)]))
# physical threshold: lowest soliton-antisoliton pole at relative rapidity i pi - 2 pi i/(h omega): 0 < Im < pi iff h omega > 2
for h, om in ((3, 0.6), (3, 1.37), (4, 0.6), (5, 0.6), (2, 1.0)):
    ub = np.pi/2 - np.pi/(h*om)
    print(f'h={h} omega={om}: h*omega={h*om:.2f}, 2 ub = {2*ub:.3f}, physical: {0 < 2*ub < np.pi}')
