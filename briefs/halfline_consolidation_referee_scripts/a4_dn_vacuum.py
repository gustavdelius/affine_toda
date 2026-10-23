"""Referee check of prp-dn-uniform-vacuum (ch. 35), independent of code.

Real coupling, m = beta = 1, |alpha|^2 = 2,
  V = sum_j n_j (e^{alpha_j.phi} - 1),  B = 2 sum_j C_j sqrt(n_j) (e^{alpha_j.phi/2} - 1),
  E = int_{-inf}^0 (phi'^2/2 + V) dx + B(phi(0)),  phi'(0) = -grad B.
Roots: alpha_0 = -e1-e2, alpha_i = e_i - e_{i+1} (1<=i<n), alpha_n = e_{n-1}+e_n.

1. closed forms: e^{f(0)} and the energies by mpmath quadrature against the book's
   (16 sqrt2 - 8 - 6 sqrt6) and (32(1+sqrt2) sin(pi/8) - 24 - 4 sqrt2);
2. full-field boundary-value problem (no symmetry assumed, random start) for
   d_4, d_5, d_6 with scipy solve_bvp, energy by quadrature; d_6 against -0.1032515;
3. convexity: Hessian of V positive definite and of B positive semidefinite on random
   points (C_j = +1), and a counterexample to convexity of B for one C_j = -1.
"""
import numpy as np
import mpmath as mp
from scipy.integrate import solve_bvp, quad
from scipy.linalg import sqrtm

mp.mp.dps = 30


def dn_roots(n):
    al = np.zeros((n + 1, n))
    al[0, 0] = -1; al[0, 1] = -1
    for i in range(1, n):
        al[i, i - 1] = 1; al[i, i] = -1
    al[n, n - 2] = 1; al[n, n - 1] = 1
    marks = np.array([1, 1] + [2] * (n - 3) + [1, 1], float)
    assert np.allclose(marks @ al, 0)
    return al, marks


def closed_forms():
    s2, s6 = mp.sqrt(2), mp.sqrt(6)
    # d4: f'' = 2(e^{2f} - e^{-f}); first integral f'^2 = 2e^{2f} + 4e^{-f} - 6
    u0 = mp.findroot(lambda u: (2 / mp.sqrt(u) - s2 * u) ** 2 - (2 * u ** 2 + 4 / u - 6), 1.1)
    f0 = mp.log(u0)
    bulk = 2 * mp.quad(lambda f: mp.sqrt(2 * mp.e ** (2 * f) + 4 * mp.e ** (-f) - 6), [0, f0])
    bdy = 2 * (4 * (mp.e ** (-f0 / 2) - 1) + s2 * (mp.e ** f0 - 1))
    E4 = bulk + bdy
    book4 = 16 * s2 - 8 - 6 * s6
    # d5: f'' = 4 sinh f; f'^2 = 8(cosh f - 1)
    u0b = mp.findroot(lambda u: (2 / mp.sqrt(u) - s2 * mp.sqrt(u)) ** 2 - (4 * u + 4 / u - 8), 1.17)
    f0b = mp.log(u0b)
    bulkb = 2 * mp.quad(lambda f: mp.sqrt(8 * (mp.cosh(f) - 1)), [0, f0b])
    bdyb = 2 * (4 * (mp.e ** (-f0b / 2) - 1) + 2 * s2 * (mp.e ** (f0b / 2) - 1))
    E5 = bulkb + bdyb
    book5 = 32 * (1 + s2) * mp.sin(mp.pi / 8) - 24 - 4 * s2
    return (u0, mp.mpf(9) / 8, E4, book4), (u0b, 4 - 2 * s2, E5, book5)


def bvp(n, C=None, L=16.0, seed=0):
    al, mk = dn_roots(n)
    if C is None:
        C = np.ones(n + 1)
    H = (al.T * mk) @ al            # Hessian of V at 0
    M = np.real(sqrtm(H))

    def Vgrad(p):
        return al.T @ (mk[:, None] * np.exp(al @ p))

    def Bgrad(p):
        return al.T @ (C * np.sqrt(mk) * np.exp(al @ p / 2))

    def fun(x, y):
        p, q = y[:n], y[n:]
        return np.vstack([q, Vgrad(p)])

    def bc(ya, yb):
        return np.concatenate([ya[n:] - M @ ya[:n], yb[n:] + Bgrad(yb[:n])])

    x = -L + L * (1 - np.linspace(1, 0, 400) ** 2)       # dense near the wall
    rng = np.random.default_rng(seed)
    y0 = np.zeros((2 * n, x.size))
    c = 0.3 * rng.normal(size=n)
    y0[:n] = np.outer(c, np.exp(np.sqrt(np.min(np.linalg.eigvalsh(H))) * x))
    y0[n:] = np.sqrt(np.min(np.linalg.eigvalsh(H))) * y0[:n]
    sol = solve_bvp(fun, bc, x, y0, tol=1e-10, max_nodes=200000)
    assert sol.status == 0, sol.message

    def edens(xx):
        y = sol.sol(xx)
        p, q = y[:n], y[n:]
        return 0.5 * q @ q + mk @ (np.exp(al @ p) - 1)
    bulk = quad(edens, -L, 0, limit=400, epsabs=1e-13, epsrel=1e-13)[0]
    # tail beyond -L, linear: int (1/2 q.q + 1/2 p.H.p) = p.M.p/2 at -L (exact for linear decay)
    p0 = sol.sol(-L)[:n]
    tail = 0.5 * p0 @ M @ p0
    pw = sol.sol(0.0)[:n]
    Bv = 2 * np.sum(C * np.sqrt(mk) * (np.exp(al @ pw / 2) - 1))
    return bulk + tail + Bv, pw, al, mk


def convexity(n=4):
    al, mk = dn_roots(n)
    rng = np.random.default_rng(3)
    mineig = np.inf
    for _ in range(200):
        p = rng.normal(size=n)
        HV = (al.T * (mk * np.exp(al @ p))) @ al
        HB = (al.T * (np.sqrt(mk) * np.exp(al @ p / 2) / 2)) @ al
        mineig = min(mineig, np.min(np.linalg.eigvalsh(HV)))
        assert np.min(np.linalg.eigvalsh(HB)) > -1e-12
    # C_0 = -1: Hessian of B has a negative direction
    C = np.ones(n + 1); C[0] = -1
    p = np.zeros(n); p += 3 * al[0] / 2
    HB = (al.T * (C * np.sqrt(mk) * np.exp(al @ p / 2) / 2)) @ al
    return mineig, np.min(np.linalg.eigvalsh(HB))


if __name__ == '__main__':
    ok = True
    (u4, b4, E4, B4), (u5, b5, E5, B5) = closed_forms()
    print(f"d4: e^f(0) = {mp.nstr(u4, 20)} (9/8 = {mp.nstr(b4, 6)}), E = {mp.nstr(E4, 20)}, book {mp.nstr(B4, 20)}")
    print(f"d5: e^f(0) = {mp.nstr(u5, 20)} (4-2sqrt2 = {mp.nstr(b5, 20)}), E = {mp.nstr(E5, 20)}, book {mp.nstr(B5, 20)}")
    ok &= abs(u4 - b4) < 1e-25 and abs(E4 - B4) < 1e-25 and abs(u5 - b5) < 1e-25 and abs(E5 - B5) < 1e-25
    for n, ref in [(4, float(B4)), (5, float(B5)), (6, -0.1032515)]:
        Es = []
        for seed in range(3):
            E, pw, al, mk = bvp(n, seed=seed)
            Es.append(E)
        e_al = al @ pw
        print(f"d{n} full-field BVP: E = {Es[0]:.9f} (spread over 3 random starts {max(Es)-min(Es):.1e}), ref {ref:.9f};"
              f" alpha_j.phi(0) = {np.round(e_al, 6)}")
        ok &= abs(Es[0] - ref) < 2e-7 and max(Es) - min(Es) < 1e-8
        if n == 4:
            ok &= abs(np.exp(e_al[2]) - 81 / 64) < 1e-7      # alpha_2.phi = 2f, e^{2f} = (9/8)^2
        if n == 5:
            ok &= abs(np.exp(e_al[2]) - (4 - 2 * np.sqrt(2))) < 1e-7
    mv, mb = convexity()
    print(f"convexity: min eigenvalue of Hess V on samples {mv:.3e} > 0; Hess B with C_0=-1 has eigenvalue {mb:.3f} < 0")
    ok &= mv > 0 and mb < 0
    print("PASS" if ok else "FAIL", "d_n uniform vacuum: closed forms, energies, full-field BVP")
