"""Part (b), matrix level.  Full 2x2 boundary bootstrap for the soliton multiplet (s, sbar) with the book's coideal conventions.
From the reflection equation (braid form, code/boundary_common.re_preserving)
   R(y1/y2)(1 x K(y1))R(y1 y2)(1 x K(y2)) = (1 x K(y2))R(y1 y2)(1 x K(y1))R(y1/y2)
the residue at a pole of K(y2) (y2 = e^{i lam p}, K(y2) -> rho_out rho_in^T) gives, with A: V(y1) -> V(1/y1) the reflection
matrix on the new boundary state,
   (1 x rho_in^T) R(y1 y2)(1 x K(y1)) R(y1/y2) (v x w) = A v * <rho_in^T, w>        for all v, w.
We check (i) that the left side factorizes in this way (independence of w: the 'rank-one' consistency), and (ii) A is proportional to
the closed-form matrix claimed in BPTT (3.6)-(3.8), i.e. X K(etabar) X on |n>, K(eta) on |n,m>, X K(etabar) X on |n,m,k>, ...
(all labels n = 0..5, 2- and 3-label states), at several (lambda, eta, vartheta) including eta beyond the strip conditions
(the identities are analytic in eta).  The scalar factor is checked in p12b_an_bootstrap.py.
Runtime: seconds.
"""
import sys, numpy as np
sys.path.insert(0, '/home/gustav/Git/affine_toda/code')
from boundary_common import q_phys, sl2_spin, rcheck, prop_residual
def report(name, ok): print(('PASS ' if ok else 'FAIL ') + name, flush=True)
X = np.array([[0, 1], [1, 0]])

def setup(lam):
    q = q_phys(lam)
    Rsg = lambda z: rcheck(sl2_spin(1, q, z), sl2_spin(1, q, 1.0))[0]
    def Kc(y, eta, vt):
        e1 = -np.cos(eta - 1j*vt)/np.sin(np.pi*lam); e0 = -np.cos(eta + 1j*vt)/np.sin(np.pi*lam)
        d = y**2 - y**-2
        return np.array([[(q - 1/q)*(e1*y + e0/y), d], [d, (q - 1/q)*(e0*y + e1/y)]])/d
    return q, Rsg, Kc

def bootstrap(Rsg, Kfun, y1, ypole):
    """Return A(y1) and the w-independence defect.  Kfun(y) = matrix of K (map V(y)->V(1/y), matrix part only)."""
    Kp = Kfun(ypole)
    U, s, Vh = np.linalg.svd(Kp)
    rank_ratio = s[1]/s[0]
    rho_out, rho_in = U[:, 0], Vh[0, :]            # K(ypole) = s0 rho_out rho_in^T (rho_in as row)
    K1 = Kfun(y1)
    D = np.kron(np.eye(2), rho_in.reshape(1, 2)) @ Rsg(y1*ypole) @ np.kron(np.eye(2), K1) @ Rsg(y1/ypole)   # (2 x 4)
    # D (v x w) = A v <rho_in, w>
    A_list = []
    for w in (np.array([1, 0.]), np.array([0, 1.]), np.array([0.3, -0.8])):
        cols = [D @ np.kron(e, w)/ (rho_in @ w if abs(rho_in @ w) > 1e-8 else np.nan) for e in np.eye(2)]
        A_list.append(np.array(cols).T.reshape(2, 2))
    ws = [np.array([1, 0.]), np.array([0, 1.]), np.array([0.3, -0.8])]
    A_ref = A_list[2]
    # independence of w, via the full bilinear structure: D(v x w) - A v <rho_in,w> = 0 for w in ker(rho_in)
    w_ker = np.array([rho_in[1], -rho_in[0]])
    dker = max(np.linalg.norm(D @ np.kron(e, w_ker)) for e in np.eye(2))/np.linalg.norm(D)
    return A_ref, rank_ratio, dker

lam_list = [(3.3, 5.5, 0.3), (7.1, 4.4, 0.9), (4.6, 6.0, 0.0), (6.3, 5.0, 1.7)]
worst = {'rank': 0, 'ker': 0, 'prop': 0}; ncases = 0
for lam, eta, vt in lam_list:
    q, Rsg, Kc = setup(lam)
    etab = np.pi*(lam + 1) - eta
    nu = lambda e, n: e/lam - (2*n + 1)*np.pi/(2*lam)
    y1 = np.exp(lam*(0.37 + 0.11j))
    # label sequences up to depth 3, n = 0..5; label at odd positions: nu_n(eta); even positions: w_m = nu_m(etab)
    seqs = [(n,) for n in range(6)] + [(n, m) for n in (0, 1, 3) for m in (0, 2, 4)] + [(1, 2, 3), (4, 0, 5)]
    for S in seqs:
        # build the matrix of the state step by step: Kcur(y)
        # ground state: Kc(y, eta); after one label: X Kc(etab) X ; after two: Kc(eta) ...
        def Kstate(k):
            ee = eta if k % 2 == 0 else etab
            if k % 2 == 0: return lambda y: Kc(y, ee, vt)
            return lambda y: X @ Kc(y, ee, vt) @ X
        for k, n in enumerate(S):
            ee = eta if k % 2 == 0 else etab
            ypole = np.exp(1j*lam*nu(ee, n))        # y = e^{lam theta}, theta = i p
            A, rk, dk = bootstrap(Rsg, Kstate(k), y1, ypole)
            target = Kstate(k + 1)(y1)
            res = prop_residual(A, target)[0]
            worst['rank'] = max(worst['rank'], rk); worst['ker'] = max(worst['ker'], dk); worst['prop'] = max(worst['prop'], res)
            ncases += 1
print(worst)
report(f'residue of the K matrix at every pole nu_n / w_m is rank one (second singular value / first) [{worst["rank"]:.0e}]', worst['rank'] < 1e-9)
report(f'bootstrap factorization (1 x rho_in^T) R K R (v x w) = A v <rho_in,w>: vanishes on ker rho_in [{worst["ker"]:.0e}]', worst['ker'] < 1e-9)
report(f'full-matrix bootstrap A = X K(etabar) X (ground -> |n>, |n,m,k>), = K(eta) (|n> -> |n,m>), {ncases} steps, n <= 5, depth <= 3 [{worst["prop"]:.0e}]', worst['prop'] < 1e-9)

# negative controls: the same A is NOT proportional to K(etabar) without the s <-> sbar swap, nor to K(eta) on |n>
lam, eta, vt = 3.3, 5.5, 0.3
q, Rsg, Kc = setup(lam); etab = np.pi*(lam + 1) - eta
y1 = np.exp(lam*(0.37 + 0.11j)); ypole = np.exp(1j*lam*(eta/lam - np.pi/(2*lam)))
A, _, _ = bootstrap(Rsg, lambda y: Kc(y, eta, vt), y1, ypole)
r_swap, r_noswap, r_same = prop_residual(A, X @ Kc(y1, etab, vt) @ X)[0], prop_residual(A, Kc(y1, etab, vt))[0], prop_residual(A, Kc(y1, eta, vt))[0]
report(f'negative controls: |0> matrix = X K(etabar) X [{r_swap:.0e}] but not K(etabar) [{r_noswap:.2f}] nor K(eta) [{r_same:.2f}]', r_swap < 1e-9 and r_noswap > 1e-2 and r_same > 1e-2)
