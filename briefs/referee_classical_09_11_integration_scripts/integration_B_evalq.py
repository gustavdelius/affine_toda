"""Evaluate E, P and the a_2^(1) spin-2/4 charges on Hirota solutions of a_n^(1), from the book's formulas.

tau_j = sum_S prod_{k in S} omega^{j a_k} E_k prod_{k<l in S} A_{a_k a_l}(theta_k - theta_l),
E_k = exp(m_a (x cosh th - t sinh th) + xi), m_a = 2 sin(pi a/h)   (m = beta = 1).
psi = i beta phi = -sum_j alpha_j ln tau_j (alpha_j = e_j - e_{j+1}), components c_c = ln tau_{c-1} - ln tau_c.
"""
import itertools
import numpy as np
import sympy as sp
import integration_B_currents as C


def A_int(h, a, b, th):
    c = np.cosh(th)
    return (c - np.cos(np.pi*(a-b)/h))/(c - np.cos(np.pi*(a+b)/h))


def tau_terms(h, sols):
    om = np.exp(2j*np.pi/h); terms = []
    N = len(sols)
    for r in range(N+1):
        for S in itertools.combinations(range(N), r):
            k = sum(2*np.sin(np.pi*sols[i][0]/h)*np.cosh(sols[i][1]) for i in S)
            w = -sum(2*np.sin(np.pi*sols[i][0]/h)*np.sinh(sols[i][1]) for i in S)
            c0 = np.prod([np.exp(sols[i][2]) for i in S]) if S else 1.0
            for i1, i2 in itertools.combinations(S, 2):
                c0 = c0*A_int(h, sols[i1][0], sols[i2][0], sols[i1][1]-sols[i2][1])
            cj = np.array([c0*om**(j*sum(sols[i][0] for i in S)) for j in range(h)])
            terms.append((cj, k, w))
    return terms


def tau_derivs(h, sols, x, t, nmax=6):
    """D^n tau_j (D = d_-) for n <= nmax, d_+ tau_j and d_+ d_- tau_j; arrays shape (h, len(x))."""
    terms = tau_terms(h, sols)
    Dn = np.zeros((nmax+1, h, len(x)), complex); P = np.zeros((h, len(x)), complex); PM = np.zeros_like(P)
    for cj, k, w in terms:
        ex = np.exp(k*x + w*t)
        for n in range(nmax+1):
            Dn[n] += np.outer(cj, (w-k)**n*ex)
        P += np.outer(cj, (w+k)*ex); PM += np.outer(cj, (w+k)*(w-k)*ex)
    return Dn, P, PM


def log_derivs(Dn):
    """Cumulants: D^n ln tau from r_n = D^n tau / tau (n = 1..5)."""
    r = [None] + [Dn[n]/Dn[0] for n in range(1, Dn.shape[0])]
    L = [None]*6
    L[1] = r[1]
    L[2] = r[2] - r[1]**2
    L[3] = r[3] - 3*r[1]*r[2] + 2*r[1]**3
    L[4] = r[4] - 4*r[1]*r[3] - 3*r[2]**2 + 12*r[1]**2*r[2] - 6*r[1]**4
    L[5] = (r[5] - 5*r[1]*r[4] - 10*r[2]*r[3] + 20*r[1]**2*r[3] + 30*r[1]*r[2]**2
            - 60*r[1]**3*r[2] + 24*r[1]**5)
    return L


def fields(h, sols, x, t):
    Dn, P, PM = tau_derivs(h, sols, x, t)
    tau = Dn[0]
    L = log_derivs(Dn)
    lp = P/tau                      # d_+ ln tau
    lpm = PM/tau - P*Dn[1]/tau**2   # d_+ d_- ln tau
    Ej = np.array([tau[(j-1) % h]*tau[(j+1) % h]/tau[j]**2 for j in range(h)])
    comp = lambda arr: np.array([arr[(c-1) % h] - arr[c] for c in range(h)])   # c_c from ln tau terms
    return dict(tau=tau, Ej=Ej, q={n: comp(L[n]) for n in range(1, 6)}, qp=comp(lp), qpm=comp(lpm))


def trap(f, x):
    return np.trapz(f, x)


def energy_momentum(h, sols, t=0.0, L=40.0, N=40001):
    x = np.linspace(-L, L, N); F = fields(h, sols, x, t)
    dm, dp = F['q'][1], F['qp']                  # d_- psi, d_+ psi (all h components)
    sq = lambda v: np.sum(v*v, axis=0)
    eps = -0.25*(sq(dp) + sq(dm)) + np.sum(1 - F['Ej'], axis=0)
    mom = 0.25*(sq(dp) - sq(dm))
    eom = np.max(np.abs(F['qpm'] + np.array([F['Ej'][c] - F['Ej'][(c-1) % h] for c in range(h)])))
    return trap(eps, x), trap(mom, x), eom, np.min(np.abs(F['tau']))


_cache = {}


def charge(s, sols, t=0.0, L=40.0, N=40001):
    """a_2^(1) spin-s charge int (T - X) dx with the current of currents.py (raw normalization)."""
    if s not in _cache:
        n, rt, rall, T, X = C.currents(s)
        syms = [C.q[k][c] for k in range(1, 6) for c in range(2)] + list(C.E)
        _cache[s] = sp.lambdify(syms, T - X, 'numpy')
    x = np.linspace(-L, L, N); F = fields(3, sols, x, t)
    args = [F['q'][k][c] for k in range(1, 6) for c in range(2)] + [F['Ej'][j] for j in range(3)]
    return trap(_cache[s](*args)*np.ones_like(x), x), np.min(np.abs(F['tau']))
