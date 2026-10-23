"""d_n, n >= 4, C = +1: the vacuum as the unique static solution in the sector invariant under the affine diagram
automorphisms (unique by convexity of the energy for C_j >= 0), found by shooting from the decaying tail at -L, and the
linearized reflection about it vs prod_j S_aj(theta, B = 1). Own code; book potentials (eq-boundary-potential)."""
import sys
import itertools
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from cl_lin import droots

TC = '/home/gustav/Git/affine_toda/code/'
ns = {}
src = open(TC + 'smatrix_simply_laced.py').read().split('\nB = 0.37')[0]
src = src.replace("open('toda_classical.py')", "open('" + TC + "toda_classical.py')")
exec(src, ns)
setup, S = ns['setup'], ns['S']

n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
A, nj = droots(n)
N1 = n + 1
# orbits of the affine D_n diagram automorphism group on the nodes 0..n (0,1 at one end, n-1, n at the other)
orb = {}
for j in range(N1):
    jj = min(j, n - j)            # flip j <-> n - j
    if jj <= 1:
        key = 'end'
    else:
        key = jj
    orb.setdefault(key, []).append(j)
# invariant subspace: psi with a_j.psi equal on each orbit
rows = []
for o in orb.values():
    for j in o[1:]:
        rows.append(A[o[0]] - A[j])
Rm = np.array(rows)
_, sv, Vt = np.linalg.svd(Rm)
rank = int(np.sum(sv > 1e-10))
P = Vt[rank:].T                    # columns span the invariant subspace
print(f'd_{n}: orbits {list(orb.values())}; invariant subspace dim {P.shape[1]}')

M0 = (A.T*nj) @ A
ev, U = np.linalg.eigh(M0)
ma = np.sqrt(ev)
# invariant mass eigenvectors
inv = [k for k in range(n) if np.linalg.norm(P @ (P.T @ U[:, k]) - U[:, k]) < 1e-9]
print('   masses', np.round(ma, 6), ' invariant modes:', [round(ma[k], 6) for k in inv])


AP = A @ P                       # alpha_j . (P z)
d = P.shape[1]


def rhs(x, y):
    z = y[:d]
    return np.concatenate([y[d:], P.T @ ((A.T*nj) @ np.expm1(AP @ z))])


L = 26/min(ma[k] for k in inv)


def blow(x, y):
    return 20 - np.abs(y[:d]).max()


blow.terminal = True


def shoot(c, dense=False, rtol=1e-13):
    x0 = -L
    psi0 = sum(c[i]*U[:, k]*np.exp(ma[k]*x0) for i, k in enumerate(inv))
    dpsi0 = sum(c[i]*ma[k]*U[:, k]*np.exp(ma[k]*x0) for i, k in enumerate(inv))
    return solve_ivp(rhs, [x0, 0.0], np.concatenate([P.T @ psi0, P.T @ dpsi0]), method='DOP853', rtol=rtol,
                     atol=1e-14*np.exp(-min(ma[k] for k in inv)*L), dense_output=dense, events=blow)


def bc(c):
    s = shoot(c)
    if s.status == 1:
        return 1e3*np.ones(P.shape[1])
    psi, dpsi = P @ s.y[:d, -1], P @ s.y[d:, -1]
    r = dpsi + (A.T*np.sqrt(nj)) @ np.exp(A @ psi/2)
    return P.T @ r


best = None
for guess in ([0.2488, -0.085],) + tuple(itertools.product([-0.25, 0.25, -1, 1, -4, 4], repeat=len(inv))):
    try:
        c, info, ier, msg = fsolve(bc, np.array(guess, float), full_output=True, xtol=1e-14)
    except Exception:
        continue
    if np.abs(bc(c)).max() < 1e-9:
        if best is None or np.abs(c - best).max() > 1e-6:
            print('   solution c =', c, ' residual', np.abs(bc(c)).max())
        best = c
        break
sol = shoot(best, dense=True)
psi = lambda x: P @ sol.sol(x)[:d]
dpsi = lambda x: P @ sol.sol(x)[d:]
print('   psi(0) =', np.round(psi(0.0), 10), ' alpha_j.psi(0) =', np.round(A @ psi(0.0), 10))
print('   e^{alpha_j.psi(0)} =', np.round(np.exp(A @ psi(0.0)), 12))
# energy
xs = np.linspace(-L, 0, 20001)
dens = [0.5*dpsi(x) @ dpsi(x) + nj @ (np.exp(A @ psi(x)) - 1) for x in xs]
E = np.trapz(dens, xs) + 2*np.sqrt(nj) @ (np.exp(A @ psi(0.0)/2) - 1)
print(f'   energy {E:.10f} m/beta^2 (trapezoid)')


# linearized reflection
def M2(x):
    return (A.T*(nj*np.exp(A @ psi(x)))) @ A


Rb = 0.5*(A.T*(np.sqrt(nj)*np.exp(A @ psi(0.0)/2))) @ A


def kmatrix(Om, Lk=None):
    Lk = Lk if Lk is not None else min(L, 34/min(ma[k] for k in inv))
    op = [a for a in range(n) if ma[a] < Om]
    cl = [a for a in range(n) if ma[a] >= Om]
    starts = []
    for a in op:
        p = np.sqrt(Om**2 - ma[a]**2)
        for sg in (1, -1):
            starts.append(('o', sg, np.concatenate([U[:, a], 1j*sg*p*U[:, a]])*np.exp(-1j*sg*p*Lk)))
    for a in cl:
        k = np.sqrt(ma[a]**2 - Om**2)
        starts.append(('c', 0, np.exp(-k*Lk)*np.concatenate([U[:, a], k*U[:, a]]).astype(complex)))
    Y0 = np.array([s[2] for s in starts]).T
    m = Y0.shape[1]
    def f(x, y):
        Y = y.reshape(2*n, m)
        return np.vstack([Y[n:], (M2(x) - Om**2*np.eye(n)) @ Y[:n]]).ravel()
    s = solve_ivp(f, [-Lk, 0.0], Y0.ravel(), method='DOP853', rtol=1e-12, atol=1e-14)
    Y = s.y[:, -1].reshape(2*n, m)
    Bm = Y[n:] + Rb @ Y[:n]
    iin = [i for i, t in enumerate(starts) if t[0] == 'o' and t[1] == 1]
    iout = [i for i, t in enumerate(starts) if t[0] == 'o' and t[1] == -1]
    icl = [i for i, t in enumerate(starts) if t[0] == 'c']
    Km = np.zeros((len(op), len(op)), complex)
    for b, ib in enumerate(iin):
        Km[:, b] = np.linalg.solve(Bm[:, iout + icl], -Bm[:, ib])[:len(op)]
    return op, Km


h, blocks, mass, _ = setup('d', n)
scale = ma.min()/mass.min()
worst, woff = 0, 0
for Om in (1.2*ma.min(), 0.6*(ma.min() + ma.max()), 1.4*ma.max()):
    op, Km = kmatrix(Om)
    woff = max(woff, np.abs(Km - np.diag(np.diag(Km))).max())
    for i, a in enumerate(op):
        th = np.arccosh(Om/ma[a])
        cands = [b for b in range(n) if abs(mass[b]*scale - ma[a]) < 1e-8]
        pred = np.prod([S(blocks[(cands[0], j)], h, 1.0, th) for j in range(n)])
        worst = max(worst, abs(Km[i, i]/pred - 1))
print(f'   linearized reflection: max |K_aa/prod_j S_aj(theta,1) - 1| = {worst:.1e}, max |offdiag| = {woff:.1e}')
