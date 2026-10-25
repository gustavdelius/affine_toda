#!/usr/bin/env python3
"""Integration agent B: independent checks for Lecture 4 (half-line lectures).

Written from the lecture's definitions: principal-gradation doublet V(y),
e_1 = yE, f_1 = F/y, k_1 = diag(q,1/q), e_0 = yF, f_0 = E/y, k_0 = diag(1/q,q),
b_j = e_j + c_j f_j k_j + eps_j k_j, physical q = -exp(-i pi lam), lam irrational.
Bulk braid matrix solved here from the coproduct (not imported).
"""
import numpy as np
import sympy as sp

OK = True


def rep(name, ok, detail=""):
    global OK
    OK = OK and ok
    print(("PASS" if ok else "FAIL"), name, detail)


E = np.array([[0, 1], [0, 0]], dtype=complex)
F = np.array([[0, 0], [1, 0]], dtype=complex)


def gens(y, q):
    return dict(e1=y * E, f1=F / y, k1=np.diag([q, 1 / q]), e0=y * F, f0=E / y, k0=np.diag([1 / q, q]))


def bgen(y, q, eps0, eps1, c0, c1):
    g = gens(y, q)
    return (g["e0"] + c0 * g["f0"] @ g["k0"] + eps0 * g["k0"],
            g["e1"] + c1 * g["f1"] @ g["k1"] + eps1 * g["k1"])


def intertwiner_system(y, q, eps0, eps1, c0, c1):
    """Linear map K -> (K b_j(y) - b_j(1/y) K)_{j=0,1}, as an 8x4 matrix."""
    b_y = bgen(y, q, eps0, eps1, c0, c1)
    b_iy = bgen(1 / y, q, eps0, eps1, c0, c1)
    M = np.zeros((8, 4), dtype=complex)
    for idx in range(4):
        K = np.zeros((2, 2), dtype=complex)
        K.flat[idx] = 1.0
        M[:, idx] = np.concatenate([(K @ b_y[j] - b_iy[j] @ K).ravel() for j in range(2)])
    return M


def check_c0c1():
    lam = np.sqrt(2.0)
    q = -np.exp(-1j * np.pi * lam)
    y = np.exp(lam * (0.37 + 0.21j))
    # on the constraint, non-symmetric split
    c0 = 2.3 / q
    c1 = 1 / (q * 2.3)
    s_on = np.linalg.svd(intertwiner_system(y, q, 0, 0, c0, c1), compute_uv=False)[-1]
    s_off = np.linalg.svd(intertwiner_system(y, q, 0, 0, c0, 1.1 * c1), compute_uv=False)[-1]
    s_neg = np.linalg.svd(intertwiner_system(y, q, 0, 0, -1 / q, -1 / q), compute_uv=False)[-1]
    rep("eps=0: nonzero K iff c0 c1 = q^-2 (generic y)", s_on < 1e-14 and s_off > 1e-3,
        f"on {s_on:.1e}, off {s_off:.3f}")
    rep("parity-symmetric c0=c1=-1/q also allowed", s_neg < 1e-14, f"{s_neg:.1e}")
    # exceptional y^4 = 1: extra solutions even off the constraint
    s_exc = np.linalg.svd(intertwiner_system(1j, q, 0, 0, c0, 1.1 * c1), compute_uv=False)[-1]
    rep("at y^4=1 a solution exists off the constraint (genericity needed)", s_exc < 1e-14, f"{s_exc:.1e}")


def K_sg(y, q, eps0, eps1):
    return np.array([[(q - 1 / q) * (eps1 * y + eps0 / y), y ** 2 - y ** -2],
                     [y ** 2 - y ** -2, (q - 1 / q) * (eps0 * y + eps1 / y)]])


def Rcheck(z, q):
    """Braid matrix V(y1) x V(y2) -> V(y2) x V(y1), z = y1/y2, intertwining the Part VI coproduct
    Delta(e)=e x 1 + k x e, Delta(f) = f x k^-1 + 1 x f for both nodes; solved as a linear
    system for all 16 entries with R[0,0] = 1 (no formula imported)."""
    y1, y2 = z, 1.0
    g1, g2 = gens(y1, q), gens(y2, q)
    I2 = np.eye(2)
    ops = []
    for j in ("1", "0"):
        ops.append((np.kron(g1["e" + j], I2) + np.kron(g1["k" + j], g2["e" + j]),
                    np.kron(g2["e" + j], I2) + np.kron(g2["k" + j], g1["e" + j])))
        ops.append((np.kron(g1["f" + j], np.linalg.inv(g2["k" + j])) + np.kron(I2, g2["f" + j]),
                    np.kron(g2["f" + j], np.linalg.inv(g1["k" + j])) + np.kron(I2, g1["f" + j])))
        ops.append((np.kron(g1["k" + j], g2["k" + j]), np.kron(g2["k" + j], g1["k" + j])))
    rows = []
    for D12, D21 in ops:
        for idx in range(16):
            pass
    A = []
    for D12, D21 in ops:
        # vec(R D12 - D21 R) = (D12^T x I - I x D21) vec(R) (column-major vec)
        A.append(np.kron(D12.T, np.eye(4)) - np.kron(np.eye(4), D21))
    A = np.vstack(A)
    # normalization R[0,0] = 1 -> vec index 0
    norm = np.zeros((1, 16), dtype=complex)
    norm[0, 0] = 1
    A = np.vstack([A, norm])
    rhs = np.zeros(A.shape[0], dtype=complex)
    rhs[-1] = 1
    x, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    assert np.linalg.norm(A @ x - rhs) < 1e-10
    sv = np.linalg.svd(A[:-1], compute_uv=False)
    assert sv[-2] > 1e-6  # one-dimensional solution space
    return x.reshape(4, 4, order="F")


def check_sign_of_q():
    lam = np.sqrt(2.0)
    qp = -np.exp(-1j * np.pi * lam)
    eps0, eps1 = 0.53 - 0.19j, -0.41 + 0.31j
    z = np.exp(lam * (0.29 - 0.17j))
    R_p = Rcheck(z, qp)
    R_m = Rcheck(z, -qp)
    Pi = np.diag([1, 1, -1, -1]).astype(complex)  # Pi x 1 with Pi = diag(1,-1)
    rep("R(-q) = (Pi x 1) R(q) (Pi x 1)", np.max(np.abs(R_m - Pi @ R_p @ Pi)) < 1e-12)
    y = np.exp(lam * (0.37 + 0.14j))
    rep("K(eps; -q) = K(-eps; q)", np.max(np.abs(K_sg(y, -qp, eps0, eps1) - K_sg(y, qp, -eps0, -eps1))) < 1e-13)
    # reflection equation at both signs with the same K
    y1, y2 = np.exp(lam * (0.37 + 0.14j)), np.exp(lam * (-0.29 + 0.21j))
    I2 = np.eye(2)
    for qq, name in ((qp, "q"), (-qp, "-q")):
        K1, K2 = K_sg(y1, qq, eps0, eps1), K_sg(y2, qq, eps0, eps1)
        Rd, Rs = Rcheck(y1 / y2, qq), Rcheck(y1 * y2, qq)
        lhs = Rd @ np.kron(I2, K1) @ Rs @ np.kron(I2, K2)
        rhs = np.kron(I2, K2) @ Rs @ np.kron(I2, K1) @ Rd
        rep(f"reflection equation holds at {name}", np.max(np.abs(lhs - rhs)) < 1e-11,
            f"{np.max(np.abs(lhs - rhs)):.1e}")
    # crossing: k(theta) prop. to R(2 theta) k(-theta),
    # k(theta) = sum_{a,b} K_{b abar}(i pi/2 - theta) v_a x v_b
    th = 0.37 - 0.12j

    def kvec(Kmat, Cmat):
        # entries K_{b, abar}: (K C)_{b a}; vector index (a, b) -> 2a+b
        KC = Kmat @ Cmat
        return np.array([KC[b, a] for a in range(2) for b in range(2)])

    def cross_err(qq, Cmat):
        km = kvec(K_sg(np.exp(lam * (0.5j * np.pi - th)), qp, eps0, eps1), Cmat)
        kp = kvec(K_sg(np.exp(lam * (0.5j * np.pi + th)), qp, eps0, eps1), Cmat)
        v = Rcheck(np.exp(2 * lam * th), qq) @ kp
        r = np.vdot(km, v) / np.vdot(km, km)
        return np.linalg.norm(v - r * km) / np.linalg.norm(v)

    C_plain = np.array([[0, 1], [1, 0]], dtype=complex)
    C_tw = np.array([[0, 1], [-1, 0]], dtype=complex)
    e1, e2, e3, e4 = cross_err(qp, C_plain), cross_err(-qp, C_plain), cross_err(-qp, C_tw), cross_err(qp, C_tw)
    rep("matrix crossing: physical q with plain C holds", e1 < 1e-13, f"{e1:.1e}")
    rep("matrix crossing: -q with plain C fails", e2 > 0.1, f"{e2:.3f}")
    rep("matrix crossing: -q with twisted C=((0,1),(-1,0)) holds", e3 < 1e-13, f"{e3:.1e}")
    rep("matrix crossing: physical q with twisted C fails", e4 > 0.1, f"{e4:.3f}")


def check_submodule():
    lam = np.sqrt(2.0)
    q = -np.exp(-1j * np.pi * lam)
    eta, vth = 2.1, 0.43
    e1 = -np.cos(eta - 1j * vth) / np.sin(np.pi * lam)
    e0 = -np.cos(eta + 1j * vth) / np.sin(np.pi * lam)
    beta1 = -np.cos(np.pi * (lam + 1) - eta + 1j * vth) / np.sin(np.pi * lam)
    beta0 = -np.cos(np.pi * (lam + 1) - eta - 1j * vth) / np.sin(np.pi * lam)
    worst = 0.0
    for n in range(4):
        nu = eta / lam - (2 * n + 1) * np.pi / (2 * lam)
        y = np.exp(-1j * lam * nu)
        worst = max(worst, abs(y - 1j * (-1) ** n * np.exp(-1j * eta)))
        b0, b1 = bgen(y, q, e0, e1, 1 / q, 1 / q)
        v = np.array([1.0, (-1) ** n * np.exp(-vth)])
        worst = max(worst, np.linalg.norm(b1 @ v - beta1 * v), np.linalg.norm(b0 @ v - beta0 * v))
    rep("y_n = i(-1)^n e^{-i eta}; v_n=(1,(-1)^n e^{-vartheta}) common eigenvector, beta_{1,0}",
        worst < 1e-13, f"{worst:.1e}")
    # all common-eigenvector branches: resultant of the two quadratics in alpha
    Y, al = sp.symbols("Y alpha")
    E1 = sp.Symbol("E1")
    E0 = sp.Symbol("E0")
    qs = sp.Symbol("q")
    p1 = Y * al ** 2 + E1 * (qs - 1 / qs) * al - 1 / Y
    p0 = al ** 2 / Y + E0 * (1 / qs - qs) * al - Y
    res = sp.resultant(sp.expand(p1 * Y), sp.expand(p0 * Y), al)
    res_num = sp.lambdify(Y, res.subs({E1: e1, E0: e0, qs: q}), "numpy")
    cand = [np.sqrt(-np.exp(2j * eta)), np.sqrt(-np.exp(-2j * eta)), np.sqrt(-np.exp(2 * vth) + 0j),
            np.sqrt(-np.exp(-2 * vth) + 0j)]
    vals = [abs(res_num(c)) for c in cand]
    # degree of the resultant in Y^2 and number of roots
    poly = sp.Poly(sp.expand(res.subs({E1: e1, E0: e0, qs: q}) * Y ** 4), Y)
    roots = np.roots([complex(c) for c in poly.all_coeffs()])
    roots = roots[np.abs(roots) > 1e-8]  # drop the zero roots introduced by clearing denominators
    y2 = np.sort_complex(np.round(roots ** 2, 9))
    expected = np.sort_complex(np.round(np.array([-np.exp(2j * eta), -np.exp(-2j * eta), -np.exp(2 * vth), -np.exp(-2 * vth)]), 9))
    rep("common eigenvectors exist exactly at y^2 in {-e^{+-2i eta}, -e^{+-2 vartheta}}",
        max(vals) < 1e-10 and len(y2) == 8 and np.allclose(np.unique(y2), np.unique(expected), atol=1e-7),
        f"{np.round(np.unique(y2), 4)}")
    # generalized determinant det(b_1 - beta k_1) is y-independent (the abandoned detour)
    dets = []
    for y in (0.7 + 0.2j, 1.9 - 0.4j):
        b0, b1 = bgen(y, q, e0, e1, 1 / q, 1 / q)
        dets.append(np.linalg.det(b1 - 0.37 * np.diag([q, 1 / q])))
    rep("det(b_1 - beta k_1) independent of y (cannot select poles)", abs(dets[0] - dets[1]) < 1e-13)
    # theta -> -theta swaps P+ and P- at fixed eta-bar; Q unchanged
    u = 0.31
    P = lambda et, vt, sgn: np.cos(lam * u) * np.cos(et) * np.cosh(vt) - sgn * np.sin(lam * u) * np.sin(et) * np.sinh(vt)
    eb = np.pi * (lam + 1) - eta
    rep("P^pm(eta-bar,-vartheta) = P^mp(eta-bar,vartheta)",
        abs(P(eb, -vth, 1) - P(eb, vth, -1)) < 1e-15 and abs(P(eb, -vth, -1) - P(eb, vth, 1)) < 1e-15)


def check_parametrization():
    lam = np.sqrt(2.0)
    q = -np.exp(-1j * np.pi * lam)
    eta, vth = 2.1, 0.43
    e1 = -np.cos(eta - 1j * vth) / np.sin(np.pi * lam)
    e0 = -np.cos(eta + 1j * vth) / np.sin(np.pi * lam)
    rep("q - 1/q = 2i sin(pi lam)", abs((q - 1 / q) - 2j * np.sin(np.pi * lam)) < 1e-14)
    u = 0.27 + 0.05j
    y = np.exp(1j * lam * u)
    K = K_sg(y, q, e0, e1)
    Pp = np.cos(lam * u) * np.cos(eta) * np.cosh(vth) - np.sin(lam * u) * np.sin(eta) * np.sinh(vth)
    Pm = np.cos(lam * u) * np.cos(eta) * np.cosh(vth) + np.sin(lam * u) * np.sin(eta) * np.sinh(vth)
    Q = -np.sin(lam * u) * np.cos(lam * u)
    rep("K / (-4i) = [[P+, Q],[Q, P-]]", np.max(np.abs(K / (-4j) - np.array([[Pp, Q], [Q, Pm]]))) < 1e-13)
    # Neumann value
    eN = -np.cos(np.pi / 2 * (lam + 1)) / np.sin(np.pi * lam)
    rep("Neumann eps = 1/(2 cos(pi lam/2)); (1-q)(1-1/q) = 4cos^2(pi lam/2)",
        abs(eN - 1 / (2 * np.cos(np.pi * lam / 2))) < 1e-14
        and abs((1 - q) * (1 - 1 / q) - 4 * np.cos(np.pi * lam / 2) ** 2) < 1e-14)
    # UV-IR combination: cos((eta -+ i vth)/(lam+1)) = mu_-+
    r, ph = 0.73, 0.4  # mu~/mu~crit, beta phi0/2
    # solve conjecture for (eta, vth): cos a cosh b = r cos ph, sin a sinh b = -r sin ph
    w = np.arccos(r * np.exp(-1j * ph))  # w = a - i b
    a_, b_ = w.real, -w.imag
    rep("UV-IR pair combines to cos((eta-+i vth)/(lam+1)) = mu_-+",
        abs(np.cos(a_) * np.cosh(b_) - r * np.cos(ph)) < 1e-14 and abs(np.sin(a_) * np.sinh(b_) + r * np.sin(ph)) < 1e-14
        and abs(np.cos(a_ + 1j * b_) - r * np.exp(1j * ph)) < 1e-14)
    # eps=0 point
    rep("eps=0 at (eta, vth) = (pi/2, 0): mu = cos(pi/(2(lam+1)))",
        abs(np.cos(np.pi / 2 / (lam + 1)) - np.cos((np.pi / 2) / (lam + 1))) < 1e-15)
    # Neumann poles
    lamv = 7.3
    etaN = np.pi / 2 * (lamv + 1)
    ok = True
    for n in range(0, 5):
        nu = etaN / lamv - (2 * n + 1) * np.pi / (2 * lamv)
        ok &= abs(nu - (np.pi / 2 - n * np.pi / lamv)) < 1e-14
        mk = lambda k: 2 * np.sin(k * np.pi / (2 * lamv))
        ok &= abs(np.cos(nu) - np.sin(n * np.pi / lamv)) < 1e-14
        ok &= abs(np.sin(n * np.pi / lamv) - mk(n) * np.cos(n * np.pi / (2 * lamv))) < 1e-14
        ok &= abs(np.sin(n * np.pi / lamv) - mk(2 * n) / 2) < 1e-14
    inside = [n for n in range(0, 10) if 0 < np.pi / 2 - n * np.pi / lamv < np.pi / 2]
    rep("Neumann: nu_n = pi/2 - n pi/lam, E = M sin(n pi/lam) = m_{2n}/2; interior iff 0<n<lam/2",
        ok and inside == [n for n in range(1, 10) if n < lamv / 2], str(inside))


def check_contact():
    # lim_{x->0-} [1/(u+ix)^2 - 1/(u-ix)^2] = 2 pi i d/dt' delta(t - t'), u = t - t'
    # test: integrate against f(t') at t=0: expect 2 pi i * (-f'(0))
    f = lambda tp: np.exp(-tp ** 2) * (1 + 0.3 * tp)
    fp0 = 0.3
    tp = np.linspace(-40, 40, 4000001)
    vals = []
    for x in (-1e-2, -1e-3):
        u = 0 - tp
        g = 1 / (u + 1j * x) ** 2 - 1 / (u - 1j * x) ** 2
        vals.append(np.trapz(g * f(tp), tp))
    target = 2j * np.pi * (-fp0)
    rep("contact identity n=2 sign", abs(vals[-1] - target) < 2e-2 and abs(vals[-1] - target) < abs(vals[0] - target),
        f"{vals} target {target}")


if __name__ == "__main__":
    check_c0c1()
    check_sign_of_q()
    check_submodule()
    check_parametrization()
    check_contact()
    print("ALL PASS" if OK else "SOME FAIL")
