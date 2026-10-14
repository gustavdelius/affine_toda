"""B4 spot check: linearized reflection about a regular continued vacuum of a_n (prp-solitonic-linear, ch. 31),
by direct ODE integration (not Bowcock's Hirota dressing used in solitonic_linear.py).
Units m = beta = 1.  alpha_j = e_j - e_{j+1 mod h} in R^h, phi = -sum_j alpha_j ln tau_j, tau_j of eq-vacuum-pair
with s = sin(pi a/h).  Pattern C found by testing all sign patterns against eq-boundary-condition (real).
Fluctuations: eps'' = (M^2(x) - Omega^2) eps, M^2 = sum_j e^{alpha_j.phi} alpha_j alpha_j^T;
eps'(0) + H eps(0) = 0, H = (1/2) sum_j C_j e^{alpha_j.phi(0)/2} alpha_j alpha_j^T.
K[c', c]: eps -> rho_c e^{ipx} + sum K[c', c] rho_c' e^{-ipx}, rho_c = (w^{jc})_j, w = e^{2 pi i/h}.
Compare with DG at B = 0: K_c = -1/((c)(h-c)) (c+b-h/2)(-c-b-h/2)(c-b+h/2)(b-c+h/2) (prp-solitonic-amplitude (iii)),
b = a (class of the pair of species a).
"""
import itertools
import numpy as np
from scipy.integrate import solve_ivp


def blk(x, th, h):
    return np.sinh((th + 1j*np.pi*x/h)/2)/np.sinh((th - 1j*np.pi*x/h)/2)


def run(n, a, chi, Om, L=22.0):
    h = n + 1
    al = np.zeros((h, h))
    for j in range(h):
        al[j, j] += 1; al[j, (j+1) % h] -= 1
    s = np.sin(np.pi*a/h); ma = 2*s
    def taus(x):
        return np.array([1 + 2*np.cos(chi + 2*np.pi*j*a/h)/(1 + s)*np.exp(ma*x) + (1 - s)/(1 + s)*np.exp(2*ma*x)
                         for j in range(h)])
    def dtaus(x):
        return np.array([2*np.cos(chi + 2*np.pi*j*a/h)/(1 + s)*ma*np.exp(ma*x) + (1 - s)/(1 + s)*2*ma*np.exp(2*ma*x)
                         for j in range(h)])
    xs = np.linspace(-L, 0, 400)
    assert all((taus(x) > 0).all() for x in xs), 'not regular'
    phi = lambda x: -al.T @ np.log(taus(x))
    dphi0 = -al.T @ (dtaus(0.0)/taus(0.0))
    E = lambda x: np.exp(al @ phi(x))
    # pattern
    best = None
    for C in itertools.product((1, -1), repeat=h):
        r = np.abs(dphi0 + al.T @ (np.array(C)*np.exp(al @ phi(0.0)/2))).max()
        if best is None or r < best[0]:
            best = (r, C)
    resC, C = best
    M2 = lambda x: (al.T*E(x)) @ al
    H = 0.5*(al.T*(np.array(C)*np.exp(al @ phi(0.0)/2))) @ al
    # basis of the root hyperplane
    Q = np.linalg.qr(np.vstack([np.ones(h), np.eye(h)[:-1]]).T)[0][:, 1:]
    sols = []
    for k in range(n):
        v = Q[:, k]
        y0 = np.concatenate([v, -H @ v]).astype(complex)
        rhs = lambda x, y: np.concatenate([y[h:], (M2(x) - Om**2*np.eye(h)) @ y[:h]])
        sol = solve_ivp(rhs, [0.0, -L], y0, method='DOP853', rtol=1e-12, atol=1e-14)
        sols.append(sol.y[:, -1])
    w = np.exp(2j*np.pi/h)
    rho = np.array([[w**(j*c) for j in range(h)] for c in range(1, h)]).T     # columns c = 1..n
    mc = 2*np.sin(np.pi*np.arange(1, h)/h)
    p = np.sqrt(Om**2 - mc**2)
    x0 = -L
    # eps = sum_c (A_c e^{ipx} + B_c e^{-ipx}) rho_c ; eps' = sum_c i p (A e^{ipx} - B e^{-ipx}) rho_c
    Mat = np.block([[rho*np.exp(1j*p*x0), rho*np.exp(-1j*p*x0)],
                    [rho*(1j*p)*np.exp(1j*p*x0), rho*(-1j*p)*np.exp(-1j*p*x0)]])
    # equations live in the root hyperplane: project with Q
    P = np.block([[Q.T, np.zeros((n, h))], [np.zeros((n, h)), Q.T]])
    coeffs = np.array([np.linalg.solve(P @ Mat, P @ y) for y in sols]).T       # (2n, n)
    A, B = coeffs[:n], coeffs[n:]
    K = B @ np.linalg.inv(A)
    th = np.arccosh(Om/mc)
    dg = []
    for c in range(1, h):
        t = th[c-1]
        Ku = -1/(blk(c, t, h)*blk(h - c, t, h))
        b = a
        M = blk(c + b - h/2, t, h)*blk(-c - b - h/2, t, h)*blk(c - b + h/2, t, h)*blk(b - c + h/2, t, h)
        dg.append(Ku*M)
    dg = np.array(dg)
    off = np.abs(K - np.diag(np.diag(K))).max()
    dev = np.abs(np.diag(K)/dg - 1).max()
    return C, resC, off, dev


for n, a, chi in ((2, 1, 0.4), (3, 1, 0.3), (4, 2, 0.5), (4, 1, 2.0)):
    for Om in (2.3, 3.1):
        C, resC, off, dev = run(n, a, chi, Om)
        print(f'a_{n} species {a} chi = {chi} Omega = {Om}: pattern {C} (BC residual {resC:.1e}); off-diagonal {off:.1e}; '
              f'|K_cc/DG - 1| = {dev:.1e}')
