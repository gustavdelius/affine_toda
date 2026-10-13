"""Coordinator check of A6 task 1: C=+1 vacuum of d_n by a full-space BVP (no symmetry ansatz), energy,
and linearized reflection vs prod_b S_ab(theta, B=1). Book root data and S-matrix; m = beta = 1."""
import sys, numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.linalg import sqrtm
exec(open('smatrix_simply_laced.py').read().split('\nB = 0.37')[0])
n = int(sys.argv[1]); L = 14.0
C_, al, roots, h, A, nj, M2 = data('d', n)
sq = np.sqrt(2*nj/np.einsum('ij,ij->i', A, A))       # sqrt(2 n_j/alpha_j^2)
def rhs(x, y):
    phi, dphi = y[:n], y[n:]
    e = np.exp(A @ phi)                               # (n+1, npts)
    return np.vstack([dphi, A.T @ (nj[:, None]*e)])
M0 = np.real(sqrtm(M2))
def bc(ya, yb):
    phi0, dphi0 = yb[:n], yb[n:]
    return np.concatenate([ya[n:] - M0 @ ya[:n],       # decaying at -L
                           dphi0 + A.T @ (sq*np.exp(A @ phi0/2))])
x = -L*(1 - np.linspace(0, 1, 400))**1.0
x = np.unique(np.concatenate([x, -np.geomspace(1e-3, L, 300)[::-1], [0.0]]))
y0 = np.zeros((2*n, x.size))
sol = solve_bvp(rhs, bc, x, y0, tol=1e-10, max_nodes=200000)
print('bvp', sol.status, sol.message)
phi0 = sol.sol(0.0)[:n]
print('alpha_j.phi(0) =', np.round(A @ phi0, 10), ' exp:', np.round(np.exp(A @ phi0), 12))
# energy
xs = np.linspace(-L, 0, 200001); Y = sol.sol(xs)
dens = 0.5*np.sum(Y[n:]**2, 0) + np.sum(nj[:, None]*(np.exp(A @ Y[:n]) - 1), 0)
from scipy.integrate import simpson
E = simpson(dens, x=xs) + 2*np.sum(sq*(np.exp(A @ phi0/2) - 1))
ref = {4: 16*np.sqrt(2) - 8 - 6*np.sqrt(6), 5: 32*(1 + np.sqrt(2))*np.sin(np.pi/8) - 24 - 4*np.sqrt(2)}
print('energy', E, ' closed form', ref.get(n))
# ---- linearized reflection
R = 0.5*sum(sq[j]*np.exp(A[j] @ phi0/2)*np.outer(A[j], A[j]) for j in range(n + 1))
def M2x(x):
    p = sol.sol(x)[:n] if x > -L else np.zeros(n)
    return sum(nj[j]*np.exp(A[j] @ p)*np.outer(A[j], A[j]) for j in range(n + 1))
ev, U = np.linalg.eigh(M2); ms = np.sqrt(ev)
hh, blocks, pf, orbits = setup('d', n)
scale = ms.min()
lab = [[a for a in range(n) if abs(pf[a]*scale - m) < 1e-8] for m in ms]
def Kmat(Om):
    cols = []
    def integ(v, dv):
        f = lambda x, y: np.concatenate([y[n:], (M2x(x) - Om**2*np.eye(n)) @ y[:n]])
        s = solve_ivp(f, [-L, 0], np.concatenate([v, dv]).astype(complex), rtol=1e-12, atol=1e-14, method='DOP853')
        return s.y[:, -1]
    basis = []; kinds = []
    for k in range(n):
        u = U[:, k]
        if Om > ms[k]:
            p = np.sqrt(Om**2 - ms[k]**2)
            for sgn, kind in ((1, 'in'), (-1, 'out')):
                basis.append(integ(u*np.exp(1j*sgn*p*(-L)), 1j*sgn*p*u*np.exp(1j*sgn*p*(-L)))); kinds.append((kind, k))
        else:
            kap = np.sqrt(ms[k]**2 - Om**2)
            basis.append(integ(u*np.exp(-kap*L)*np.exp(kap*L/2), kap*u*np.exp(-kap*L)*np.exp(kap*L/2))); kinds.append(('dec', k))
    G = lambda y: y[n:] + R @ y[:n]                    # boundary functional
    ins = [i for i, kd in enumerate(kinds) if kd[0] == 'in']
    oth = [i for i, kd in enumerate(kinds) if kd[0] != 'in']
    Mo = np.array([G(basis[i]) for i in oth]).T
    out = {}
    for i in ins:
        c = np.linalg.solve(Mo, -G(basis[i]))
        out[kinds[i][1]] = {kinds[o][1]: c[j] for j, o in enumerate(oth) if kinds[o][0] == 'out'}
    return out
worst = 0; offd = 0
for Om in (1.05*ms.max(), 0.5*(ms[0] + ms[1]) if ms[1] - ms[0] > 1e-6 else 1.3*ms[0]):
    K = Kmat(Om)
    for k, row in K.items():
        th = np.arccosh(Om/ms[k])
        for k2, v in row.items():
            if abs(ms[k2] - ms[k]) > 1e-8 or k2 != k:
                offd = max(offd, abs(v)) if k2 != k else offd
        pred = [np.prod([S(blocks[(a, b)], h, 1.0, th) for b in range(n)]) for a in lab[k]]
        worst = max(worst, min(abs(row[k] - p) for p in pred))
        print(f'Om={Om:.4f} mass {ms[k]:.4f} labels {lab[k]} K={row[k]:.12f} prodS={pred[0]:.12f}')
print('max |K - prod S| =', worst, ' max off-diagonal =', offd)
# K_1 vs B -> 0 limit of eq-dn-uniform-amplitude
pr = lambda x, th: np.sinh(th/2 + 1j*np.pi*x/(2*h))/np.sinh(th/2 - 1j*np.pi*x/(2*h))
k1 = [k for k in range(n) if 0 in lab[k]][0]
for th in (0.4, 1.1):
    Om = ms[k1]*np.cosh(th)
    if Om < ms.max():
        lim = -pr(h/2 - 1, th)*pr(h/2 + 1, th)/(pr(1, th)*pr(h - 1, th)*pr(h/2, th)**2)
        print('K_1 vs eq-dn-uniform-amplitude at B=0:', abs(Kmat(Om)[k1][k1] - lim))
