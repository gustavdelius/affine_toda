"""Referee check of prp-dn-uniform-classical (ch. 35): linearized reflection about the C = +1 vacuum.

Own code: vacuum from the closed forms of prp-dn-uniform-vacuum (d_4, d_5) and from a full-field BVP
(d_6, no symmetry assumed); multichannel scattering integrated FORWARD from x = -L, starting from
the exact asymptotic modes (incoming/outgoing waves in open channels, decaying exponentials in
closed channels), boundary condition eps'(0) = -R eps(0), R = Hess B(psi(0)).
Comparison: prod_b S_ab(theta, B=1) with the block content from the book's smatrix_simply_laced
(the brief allows importing it), particles matched by mass.
Also: no bound state (scan of the bound-state determinant for 0 < Omega < m_min and for
Omega^2 < 0), and the positions of the physical-strip poles of prod_b S_ab(theta, 1) in p = m_a sinh(theta).
Units m = beta = 1, |alpha|^2 = 2; roots alpha_0 = -e1-e2, alpha_i = e_i - e_{i+1}, alpha_n = e_{n-1}+e_n.
"""
import os
import sys
import io
import contextlib
import cmath
import numpy as np
from scipy.integrate import solve_ivp, solve_bvp
from scipy.linalg import sqrtm

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.join(HERE, '..', '..', 'code')


def load_smatrix():
    cwd = os.getcwd()
    os.chdir(TC)
    sys.path.insert(0, TC)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            import smatrix_simply_laced as sm
    finally:
        os.chdir(cwd)
    return sm


def dn_roots(n):
    al = np.zeros((n + 1, n))
    al[0, 0] = -1; al[0, 1] = -1
    for i in range(1, n):
        al[i, i - 1] = 1; al[i, i] = -1
    al[n, n - 2] = 1; al[n, n - 1] = 1
    mk = np.array([1, 1] + [2] * (n - 3) + [1, 1], float)
    return al, mk


def vacuum(n):
    """returns psi(x) (vectorized over x <= 0) and psi(0)"""
    al, mk = dn_roots(n)
    if n == 4:
        v = np.zeros(4); v[1] = 1; v[2] = -1
        E0 = (3 * np.sqrt(3) - 5) / (3 * np.sqrt(3) + 5)

        def f(x):
            E = E0 * np.exp(np.sqrt(6) * x)
            s = np.sqrt(3) * (1 - E) / (1 + E)
            return np.log(2 / (s ** 2 - 1))
        return (lambda x: np.outer(f(np.atleast_1d(x)), v)), v * f(0.0)
    if n == 5:
        v = np.zeros(5); v[1] = 1; v[3] = -1
        f0 = np.log(4 - 2 * np.sqrt(2))
        x0 = -0.5 * np.log(np.tanh(f0 / 4))      # tanh(f/4) = e^{2(x - x0)}

        def f(x):
            return 4 * np.arctanh(np.exp(2 * (x - x0)))
        return (lambda x: np.outer(f(np.atleast_1d(x)), v)), v * f(0.0)
    # general n: full-field BVP
    H = (al.T * mk) @ al
    Msq = np.real(sqrtm(H))
    L = 16.0

    def fun(x, y):
        p = y[:n]
        return np.vstack([y[n:], al.T @ (mk[:, None] * np.exp(al @ p))])

    def bc(ya, yb):
        return np.concatenate([ya[n:] - Msq @ ya[:n], yb[n:] + al.T @ (np.sqrt(mk) * np.exp(al @ yb[:n] / 2))])
    x = -L + L * (1 - np.linspace(1, 0, 400) ** 2)
    rng = np.random.default_rng(0)
    kmin = np.sqrt(np.min(np.linalg.eigvalsh(H)))
    y0 = np.zeros((2 * n, x.size))
    y0[:n] = np.outer(0.3 * rng.normal(size=n), np.exp(kmin * x))
    y0[n:] = kmin * y0[:n]
    sol = solve_bvp(fun, bc, x, y0, tol=1e-10, max_nodes=300000)
    assert sol.status == 0, sol.message

    def psi(xx):
        xx = np.atleast_1d(xx)
        out = []
        for t in xx:
            if t >= -L:
                out.append(sol.sol(t)[:n])
            else:                                   # linear tail
                p0 = sol.sol(-L)[:n]
                w, U = np.linalg.eigh(H)
                out.append(U @ (np.exp(np.sqrt(w) * (t + L)) * (U.T @ p0)))
        return np.array(out)
    return psi, sol.sol(0.0)[:n]


def check_vacuum(n, psi, p0):
    al, mk = dn_roots(n)
    xs = np.linspace(-6, -0.01, 50)
    h = 1e-4
    res = 0
    for x in xs:
        d2 = (psi(x + h)[0] - 2 * psi(x)[0] + psi(x - h)[0]) / h ** 2
        res = max(res, np.abs(d2 - al.T @ (mk * np.exp(al @ psi(x)[0]))).max())
    d1 = (psi(-1e-6)[0] - psi(-1e-6 - 1e-5)[0]) / 1e-5
    bcres = np.abs(d1 + al.T @ (np.sqrt(mk) * np.exp(al @ p0 / 2))).max()
    return res, bcres


def kmatrix(n, psi, p0, Om, L=18.0):
    al, mk = dn_roots(n)
    H0 = (al.T * mk) @ al
    w, U = np.linalg.eigh(H0)
    masses = np.sqrt(w)
    R = 0.5 * (al.T * (np.sqrt(mk) * np.exp(al @ p0 / 2))) @ al
    opn = [a for a in range(n) if Om > masses[a]]
    cls = [a for a in range(n) if Om <= masses[a]]
    k = {a: np.sqrt(Om ** 2 - masses[a] ** 2) for a in opn}
    kap = {a: np.sqrt(masses[a] ** 2 - Om ** 2) for a in cls}

    def rhs(x, y):
        e = y[:n]
        M2 = (al.T * (mk * np.exp(al @ psi(x)[0]))) @ al
        return np.concatenate([y[n:], (M2 - Om ** 2 * np.eye(n)) @ e])
    starts = []          # (label, initial y at -L)
    for a in opn:
        for sgn, lab in ((+1, 'in'), (-1, 'out')):
            amp = np.exp(1j * sgn * k[a] * (-L))
            starts.append(((lab, a), np.concatenate([U[:, a] * amp, 1j * sgn * k[a] * U[:, a] * amp])))
    for b in cls:
        amp = np.exp(kap[b] * (-L))
        starts.append((('dec', b), np.concatenate([U[:, b] * amp, kap[b] * U[:, b] * amp])))
    finals = {}
    for lab, y0 in starts:
        sol = solve_ivp(rhs, (-L, 0.0), y0.astype(complex), method='DOP853', rtol=1e-12, atol=1e-14)
        yf = sol.y[:, -1]
        finals[lab] = yf[n:] + R @ yf[:n]          # boundary-condition residual vector
    # for each incoming channel a: solve sum_b r_b out_b + sum_c d_c dec_c = -in_a
    unknown = [('out', a) for a in opn] + [('dec', b) for b in cls]
    Mx = np.array([finals[u] for u in unknown]).T
    K = np.zeros((len(opn), len(opn)), complex)
    for i, a in enumerate(opn):
        sol = np.linalg.solve(Mx, -finals[('in', a)])
        K[:, i] = sol[:len(opn)]
    th = {a: np.arccosh(Om / masses[a]) for a in opn}
    return K, opn, masses, th, U


def bound_det(n, psi, p0, Om2, L=18.0):
    al, mk = dn_roots(n)
    H0 = (al.T * mk) @ al
    w, U = np.linalg.eigh(H0)
    R = 0.5 * (al.T * (np.sqrt(mk) * np.exp(al @ p0 / 2))) @ al

    def rhs(x, y):
        M2 = (al.T * (mk * np.exp(al @ psi(x)[0]))) @ al
        return np.concatenate([y[n:], (M2 - Om2 * np.eye(n)) @ y[:n]])
    cols = []
    for b in range(n):
        kap = np.sqrt(w[b] - Om2)
        y0 = np.concatenate([U[:, b], kap * U[:, b]]) * np.exp(-kap * L)
        sol = solve_ivp(rhs, (-L, 0.0), y0, method='DOP853', rtol=1e-11, atol=1e-30)
        yf = sol.y[:, -1]
        v = yf[n:] + R @ yf[:n]
        cols.append(v / np.linalg.norm(yf))
    return np.linalg.det(np.array(cols).T)


def main():
    sm = load_smatrix()
    ok = True
    for n in (4, 5, 6):
        al, mk = dn_roots(n)
        h, bl, pfm, _ = sm.setup('d', n)
        H0 = (al.T * mk) @ al
        masses = np.sqrt(np.linalg.eigvalsh(H0))
        mmin = masses.min()
        node_mass = {a: pfm[a] * mmin for a in range(n)}
        psi, p0 = vacuum(n)
        vres, bres = check_vacuum(n, psi, p0)
        print(f"d_{n}: vacuum field-equation residual {vres:.1e}, boundary condition {bres:.1e}")

        def prodS(m, th):
            nodes = [a for a in range(n) if abs(node_mass[a] - m) < 1e-8]
            vals = []
            for a in nodes:
                p = 1
                for b in range(n):
                    p *= sm.S(bl[(a, b)], h, 1.0, th)
                vals.append(p)
            assert max(abs(v - vals[0]) for v in vals) < 1e-10      # same for degenerate nodes
            return vals[0]
        thresholds = sorted(set(np.round(masses, 10)))
        Oms = [0.5 * (thresholds[0] + thresholds[1])] + [thresholds[-1] * 1.3, thresholds[-1] * 2.1]
        if len(thresholds) > 2:
            Oms.append(0.5 * (thresholds[-2] + thresholds[-1]))
        worst_d, worst_off = 0, 0
        for Om in Oms:
            K, opn, ms, th, U = kmatrix(n, psi, p0, Om)
            for i, a in enumerate(opn):
                target = prodS(ms[a], th[a])
                worst_d = max(worst_d, abs(K[i, i] - target))
                for j, b in enumerate(opn):
                    if j != i:
                        worst_off = max(worst_off, abs(K[i, j]))
        print(f"   K = diag(prod_b S_ab(theta_a, 1)) at Omega = {np.round(Oms, 3)}: diag {worst_d:.1e}, off-diagonal {worst_off:.1e}")
        ok &= worst_d < 1e-8 and worst_off < 1e-8
        if n in (4, 5):
            # vector particle: B -> 0 limit of eq-dn-uniform-amplitude
            Om = 1.3 * thresholds[-1]
            K, opn, ms, th, U = kmatrix(n, psi, p0, Om)
            a1 = [i for i, a in enumerate(opn) if abs(ms[a] - 2 * np.sqrt(2) * np.sin(np.pi / h)) < 1e-9]
            t = th[opn[a1[0]]]
            b = lambda x: cmath.sinh(t / 2 + 1j * np.pi * x / (2 * h)) / cmath.sinh(t / 2 - 1j * np.pi * x / (2 * h))
            K1cl = -b(h / 2 - 1) * b(h / 2 + 1) / (b(1) * b(h - 1) * b(h / 2) ** 2)
            print(f"   vector particle: K_11 vs B->0 limit of eq-dn-uniform-amplitude: {abs(K[a1[0], a1[0]] - K1cl):.1e}")
            ok &= abs(K[a1[0], a1[0]] - K1cl) < 1e-8
        # bound states: scan Omega^2 from -4 to m_min^2
        grid = np.linspace(-4.0, mmin ** 2 - 1e-4, 300 if n < 6 else 120)
        dets = np.array([bound_det(n, psi, p0, g) for g in grid])
        changes = int(np.sum(np.sign(dets[1:]) != np.sign(dets[:-1])))
        print(f"   bound-state determinant: {changes} sign changes for -4 < Omega^2 < m_min^2 (min |det| {np.abs(dets).min():.2e})")
        ok &= changes == 0
        # poles of prod_b S_ab(theta, 1) in the physical strip, exactly from the block content:
        # {x}|_{B=1} = (x-1)(x+1)/(x)^2; (y) has a pole at theta = i pi y/h and a zero at -i pi y/h (period 2h)
        poles = {}
        for a in range(n):
            e = {}
            for c in range(n):
                for x, k in bl[(a, c)].items():
                    for y, m in ((x - 1, k), (x + 1, k), (x, -2 * k)):
                        e[y % (2 * h)] = e.get(y % (2 * h), 0) + m
            for z in range(1, h):
                order = e.get(z % (2 * h), 0) - e.get((-z) % (2 * h), 0)
                if order > 0 and 2 * z < h:        # pole with 0 < Im theta < pi/2
                    rate = round(2 * node_mass[a] * np.sin(np.pi * z / h), 6)
                    poles.setdefault(rate, set()).add((a, z, order))
        ms = sorted(set(np.round(masses, 6)))
        print(f"   physical-strip poles of prod_b S_ab(theta,1): 2|p| in {sorted(poles)}; distinct masses {ms}")
        if n == 6:
            sums = sorted(set(round(a + b, 6) for a in ms for b in ms))
            print(f"   (d_6 sums of two masses: {sums})")
    print("PASS" if ok else "FAIL", "linearized reflection about the d_n C = +1 vacuum (n = 4, 5, 6)")


if __name__ == '__main__':
    main()
