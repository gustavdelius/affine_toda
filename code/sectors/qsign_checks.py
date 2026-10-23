"""Direct comparison of q and -q (book: eq-q-physical, sec-crossing-sign, sec-sector-status).

1. Spinor R-matrices of c_n^(1) (U_q(d_{n+1}^(2))) and a_{2n-1}^(2) (U_q(b_n^(1))), n = 2, 3:
   R(x; -q) = L R(x; q) L with L diagonal, entries +-1, and L is not swap-symmetric,
   so it is not a product A (x) A of one-soliton basis changes.
2. a_2^(2) kinks in the Takacs-Watts conventions: two-kink eigenvalues (closed form) and
   N = 3, 4 kink transfer matrices have the same eigenvalue moduli at q and -q.
"""
import numpy as np, sys

# ---- 1. spinor R-matrices
for N_, alg_ in ((2, 'c'), (3, 'c'), (2, 'a'), (3, 'a')):
    sys.argv = ['x', str(N_), alg_]
    exec(open('spinor_scan.py').read().split('if __name__ == "__main__":')[0])
    q = np.exp(-1j*np.pi*2.37)
    Pp, Pm = setup(q), setup(-q)
    ok = True; Lv = None
    for x in (1.7 + 0.4j, 0.3 - 1.1j):
        Rp, Rm = Rmat(x, q, Pp), Rmat(x, -q, Pm)
        nz = np.abs(Rp) > 1e-10; D2 = Rp.shape[0]
        ok &= np.array_equal(nz, np.abs(Rm) > 1e-10)
        r = np.zeros((D2, D2)); r[nz] = (Rm[nz]/Rp[nz]).real
        ok &= np.allclose(np.abs(r[nz]), 1, atol=1e-10)
        sgn = np.round(r).astype(int)
        L = {}
        for seed in range(D2):                       # R conserves weight: one seed per connected block
            if seed in L: continue
            L[seed] = 1; stack = [seed]
            while stack:
                i = stack.pop()
                for j in np.nonzero(nz[i] | nz[:, i])[0]:
                    s = sgn[i, j] if nz[i, j] else sgn[j, i]
                    if j not in L: L[j] = s*L[i]; stack.append(j)
                    elif L[i]*L[j] != s: ok = False
        Lv = np.array([L[i] for i in range(D2)]).reshape(D, D)
    asym = int(np.sum(Lv*Lv.T < 0))//2
    print(f"{alg_}_{N_} spinor: R(-q) = L R(q) L with diagonal L = +-1: {bool(ok)};  "
          f"pairs of one-soliton states with L_ij != L_ji: {asym} of {D*(D-1)//2} (so L is not A(x)A)")

# ---- 2. a_2^(2) kinks
exec(open('a22cubic.py').read().split('lxs = ')[0])      # cubic() and allev() only, without the module's table
lxs = np.concatenate([np.linspace(1e-3, 6, 120), np.linspace(6, 40, 20)])
for phi in (0.13, 0.37, 0.61, 0.83):
    d = [max(np.abs(np.abs(allev(np.exp(l), s*np.exp(1j*np.pi*phi))) - 1).max() for l in lxs) for s in (1, -1)]
    print(f"a_2^(2) two kinks, q = +-e^(i pi {phi}): max||s|-1| {d[0]:.2e} (+q), {d[1]:.2e} (-q)")
from lib import scanN
exec(open('a22three.py').read().split('q = np.exp(-1j*0.77)')[0])
for xi_pi in (0.30, 0.45, 2.0, 3.0):
    xi = (xi_pi + 0.0013)*np.pi; qTW = np.conj(1j*np.exp(1j*np.pi**2/(3*xi)))
    row = []
    for lab, qq in (("q", qTW), ("-q", -qTW)):
        Rpl = setup(qq)
        (d3, _), _ = scanN(Rpl, 3, 3, L=6, n=150)
        (d4, _), _ = scanN(Rpl, 3, 4, L=6, n=60)
        row.append(f"{lab}: N=3 {d3:.1e}, N=4 {d4:.1e}")
    print(f"a_2^(2) kinks, xi = {xi_pi} pi: " + " | ".join(row), flush=True)
