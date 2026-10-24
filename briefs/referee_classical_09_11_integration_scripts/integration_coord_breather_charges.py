"""Coordinator's independent check of prp-breather-charges (a_2^(1), spins 2 and 4).

Written from the book's definitions only: e^{i beta alpha_k.phi} = prod_j tau_j^{-a_jk}
(eq-hirota-ansatz), the a_n solitons (eq-an-soliton, eq-n-soliton) with m_a = 2 m sin(pi a/h),
A_ab of eq-interaction-coefficient, d_pm = d_t pm d_x. Fields enter through y_k = alpha_k . psi,
psi = i beta phi (k = 0, 1, 2; y_0 = -y_1 - y_2), which obey d+d- y_k = -m^2 sum_j C_kj e^{y_j}.
Currents: d+ T = d- X, T polynomial in d- jets of y_1, y_2, X linear in the e^{y_j};
Q = int (T - X) dx.  m = 1.
"""
import itertools, numpy as np, sympy as sp

C = sp.Matrix([[2, -1, -1], [-1, 2, -1], [-1, -1, 2]])   # affine Cartan matrix of a_2^(1)

def monomials(weight, nmax):
    """Monomials in jets J[k][n] (k=1,2; n>=1) of total weight `weight`, as tuples of (k,n)."""
    gens = [(k, n) for n in range(1, nmax + 1) for k in (1, 2)]
    out = []
    def rec(start, w, cur):
        if w == 0:
            out.append(tuple(cur)); return
        for i in range(start, len(gens)):
            if gens[i][1] <= w:
                rec(i, w - gens[i][1], cur + [gens[i]])
    rec(0, weight, [])
    return out

def build_current(s):
    N = s + 2
    J = {(k, n): sp.Symbol(f'J{k}_{n}') for k in (1, 2) for n in range(1, N + 1)}
    e = sp.symbols('e0 e1 e2')
    Y = {0: -J[1, 1] - J[2, 1], 1: J[1, 1], 2: J[2, 1]}
    def Dm(f):   # total d- derivative on the jet algebra
        r = sum(sp.diff(f, J[k, n]) * J[k, n + 1] for (k, n) in J if n < N)
        r += sum(sp.diff(f, e[j]) * e[j] * Y[j] for j in range(3))
        return sp.expand(r)
    # d+ of the jets: d+ d-^n y_k = d-^(n-1) (-sum_j C_kj e_j)
    dplus = {}
    for k in (1, 2):
        f = -sum(C[k, j] * e[j] for j in range(3))
        for n in range(1, N):
            dplus[k, n] = f
            f = Dm(f)
    def Dp(f):   # f polynomial in jets only (no exponentials)
        return sp.expand(sum(sp.diff(f, J[k, n]) * dplus[k, n] for (k, n) in J if sp.diff(f, J[k, n]) != 0))
    Tm = [sp.Mul(*[J[g] for g in mon]) for mon in monomials(s + 1, s + 1)]
    Xm = [e[j] * sp.Mul(*[J[g] for g in mon]) for j in range(3) for mon in monomials(s - 1, s - 1)]
    # image of each basis element under T -> d+T, X -> -d-X, as coefficient dictionaries
    cols = [sp.expand(Dp(m)).as_coefficients_dict() for m in Tm] + \
           [sp.expand(-Dm(m)).as_coefficients_dict() for m in Xm]
    keys = sorted({k for c in cols for k in c}, key=sp.default_sort_key)
    row = {k: i for i, k in enumerate(keys)}
    M = sp.zeros(len(keys), len(cols))
    for jc, c in enumerate(cols):
        for k, v in c.items():
            M[row[k], jc] += v
    ns = M.nullspace()
    return J, e, Tm, Xm, ns

def tau_terms(data):
    """tau_j as lists of (coefficient, exponent vector k) with tau_j = sum c prod_i E_i^{k_i}."""
    h = 3; om = np.exp(2j * np.pi / h); n = len(data)
    def A(a, b, th):
        return (np.cosh(th) - np.cos(np.pi * (a - b) / h)) / (np.cosh(th) - np.cos(np.pi * (a + b) / h))
    taus = []
    for j in range(3):
        terms = []
        for sub in itertools.product((0, 1), repeat=n):
            c = 1 + 0j; idx = [i for i in range(n) if sub[i]]
            for i in idx:
                c *= om ** (j * data[i][0])
            for p, q in itertools.combinations(idx, 2):
                c *= A(data[p][0], data[q][0], data[p][1] - data[q][1])
            terms.append((c, np.array(sub)))
        taus.append(terms)
    return taus

def tau_derivs(terms, Ex, lam, nmax, lam2=None):
    """d^n tau for n = 0..nmax along a derivation with d E_i = lam_i E_i (numerically exact).
    With lam2, also return d2 d1 tau (mixed, first order in each)."""
    out = [0] * (nmax + 1); mixed = 0
    for c, k in terms:
        mon = c * np.prod([Ex[i] ** k[i] for i in range(len(k))], axis=0)
        rate = sum(lam[i] * k[i] for i in range(len(k)))
        for n in range(nmax + 1):
            out[n] = out[n] + rate ** n * mon
        if lam2 is not None:
            mixed = mixed + rate * sum(lam2[i] * k[i] for i in range(len(k))) * mon
    return out, mixed

def log_derivs(td):
    """derivatives g^(n), n = 1..N, of g = ln tau from tau^(n) by the Leibniz recursion."""
    from math import comb
    N = len(td) - 1; g = [None] * (N + 1)
    for n in range(1, N + 1):
        acc = td[n]
        for k in range(0, n - 1):
            acc = acc - comb(n - 1, k) * g[k + 1] * td[n - 1 - k]
        g[n] = acc / td[0]
    return g

def charges(s, J, e, Tm, Xm, vec, data, t0):
    m = [2 * np.sin(np.pi * a / 3) for (a, th, xi) in data]
    kappa = min(m[i] * np.real(np.cosh(data[i][1])) for i in range(len(data)))   # slowest decay rate
    L = 36.0 / kappa; npts = int(800 * L) + 1                                      # tails below e^-36
    lam_m = [-m[i] * np.exp(data[i][1]) for i in range(len(data))]   # d- E_i = lam E_i
    lam_p = [m[i] * np.exp(-data[i][1]) for i in range(len(data))]   # d+ E_i = lam E_i
    xs = np.linspace(-L, L, npts)
    Ex = [np.exp(m[i] * np.cosh(data[i][1]) * xs - m[i] * np.sinh(data[i][1]) * t0 + data[i][2]) for i in range(len(data))]
    taus = tau_terms(data)
    tds, gs, mixed, d1p, d1m = [], [], [], [], []
    for terms in taus:
        td, mx = tau_derivs(terms, Ex, lam_m, s + 1, lam2=lam_p)
        tp, _ = tau_derivs(terms, Ex, lam_p, 1)
        tds.append(td); gs.append(log_derivs(td))
        # d+ d- ln tau = (tau d+d- tau - d+tau d-tau)/tau^2
        mixed.append((td[0] * mx - tp[1] * td[1]) / td[0] ** 2)
    ev = [tds[(j - 1) % 3][0] * tds[(j + 1) % 3][0] / tds[j][0] ** 2 for j in range(3)]
    def ycomb(k, arr): return -(2 * arr[k] - arr[(k - 1) % 3] - arr[(k + 1) % 3])
    res = max(np.max(np.abs(ycomb(k, mixed) + sum(float(C[k, j]) * ev[j] for j in range(3)))) for k in range(3))
    vals = {}
    for k in (1, 2):
        for n in range(1, s + 2):
            vals[J[k, n]] = ycomb(k, [gs[j][n] for j in range(3)])
    for j in range(3):
        vals[e[j]] = ev[j]
    nT = len(Tm)
    dens = sum(vec[i] * Tm[i] for i in range(nT)) - sum(vec[nT + i] * Xm[i] for i in range(len(Xm)))
    syms = list(vals.keys())
    f = sp.lambdify(syms, dens, 'numpy')
    Q = np.trapz(f(*[vals[k] for k in syms]) * np.ones_like(xs), xs)
    return Q, res, min(np.min(np.abs(td[0])) for td in tds)

if __name__ == '__main__':
    ok = True
    for s in (2, 3, 4):
        J, e, Tm, Xm, ns = build_current(s)
        if s == 3:
            zs = [abs(charges(s, J, e, Tm, Xm, [complex(sp.N(c)) for c in v], [(1, 0.2, 0.3)], 0.0)[0]) for v in ns]
            print(f"spin 3: nullspace dimension {len(ns)}; largest charge on a soliton {max(zs):.1e} (all total derivatives)")
            ok &= max(zs) < 1e-9; continue
        print(f"spin {s}: nullspace dimension {len(ns)}")
        rng = np.random.default_rng(s)
        vec = sum((complex(rng.normal()) * v for v in ns), sp.zeros(len(ns[0]), 1))
        vec = [complex(sp.N(c)) for c in vec]
        q1, r1, _ = charges(s, J, e, Tm, Xm, vec, [(1, 0, 0.3)], 0.0)
        q1b, _, _ = charges(s, J, e, Tm, Xm, vec, [(1, 0.4 + 0.3j, 0.2 - 0.5j)], 0.7)
        q2, _, _ = charges(s, J, e, Tm, Xm, vec, [(2, 0, -0.1)], 0.0)
        print(f"  q2/q1 = {q2 / q1:.12f}   soliton at theta=0.4+0.3i: Q/(q1 e^(s theta)) = {q1b / (q1 * np.exp(s * (0.4 + 0.3j))):.12f}")
        ok &= abs(q2 / q1 + 1) < 1e-9 and abs(q1b / (q1 * np.exp(s * (0.4 + 0.3j))) - 1) < 1e-9
        for (th0, u, xi1, t0) in [(0.3, 0.8, 0.4 + 0.2j, 0.0), (-0.5, 1.1, -0.3 + 1.0j, 1.3), (0.1, 1.4, 0.0, -0.6)]:
            data = [(1, th0 + 1j * u, xi1), (2, th0 - 1j * u, np.conj(xi1) + 1j * np.pi)]
            A = 1 - np.sin(np.pi / 3) ** 2 / np.cos(u) ** 2
            Q, res, tmin = charges(s, J, e, Tm, Xm, vec, data, t0)
            pred = 2j * np.exp(s * th0) * np.sin(s * u)
            print(f"  breather th0={th0}, u={u} (A={A:.3f}), t={t0}: Q/q1 = {Q / q1:.12f}, predicted {pred:.12f}; "
                  f"EOM residual {res:.1e}, min|tau| {tmin:.2f}")
            ok &= abs(Q / q1 - pred) < 1e-9 and res < 1e-8 and tmin > 1e-3
    print("PASS" if ok else "FAIL", "breather charges Q_s = 2i e^{s theta0} sin(s u) q_1, s = 2, 4, q_2 = -q_1")
