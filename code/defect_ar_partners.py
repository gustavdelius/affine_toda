"""Type-I defects of affine Toda theory: the two halves of the defect potential, the partner and neighbour argument,
the a_r defect in components, unfreezing, and the reflection relation between the two types (Part VIII, ch. 38, 39).

Book conventions (ch. 38): fields phi on x < 0, psi on x > 0; V = (m^2/b^2) sum_j n_j (e^{b alpha_j.phi} - 1) with
b = beta (real coupling) or b = i beta (imaginary coupling); the defect conditions (eq-defect-conditions)
    phi_x = E phi_t + D psi_t - grad_phi B_d,   psi_x = D^T phi_t - F psi_t + grad_psi B_d   at x = 0,
the type-I_- defect D = 2 (1 - pi)^{-1}, E = F = 1 - D (eq-d-matrix), pi alpha_j = alpha_{j-1}, and
    B_d = (m/b^2) sum_j [sigma (f_j - 1) + sigma^{-1} (g_j - 1)],  f_j = e^{b alpha_j.(D^T phi + D psi)/2},
    g_j = e^{b alpha_j.D (phi - psi)/2}                                                         (eq-defect-potential);
type I_+ has (D, E, F) -> (-D, -E, -E) and sigma <-> sigma^{-1} in B_d with the same f_j, g_j. Solitons at imaginary
coupling: phi = (i/beta) sum_j alpha_j ln tau_j, tau_j = 1 + omega^{ja} E, omega = e^{2 pi i/h} (eq-an-soliton).
Written from these definitions, without defect_common.py. One PASS or FAIL line per check; exit status 1 on any FAIL.

Checks (book ids they support):
 1. Condition 2 of prp-defect-momentum in the coordinates xi = D^T phi + D psi, eta = phi - psi (inverse
    phi = (xi + D eta)/2, psi = (xi - D^T eta)/2): for B = f(xi) + g(eta) with arbitrary f, g it holds with P_d = g - f
    (sympy, two fields, E = eps J); a mixed term xi_1 eta_2 gives a curl.                          [sec-defect-momentum, eq-dc-halves]
 2. With D = 1 - E, F = E and B = f(xi) + g(eta), the defect conditions are d_+(phi - psi) = -2 grad f(xi) and
    d_-(D^T phi + D psi) = 2 grad g(eta), and condition 3 is 2 grad f . grad g = V(phi) - V(psi) (D is normal);
    for the a_r defect f and g are its sigma and sigma^{-1} halves, at both couplings.                [eq-dc-lightcone]
 3. prp-ar-only, the neighbour property: 23 affine diagrams of rank <= 4 and their 78 automorphisms (64 fix a node,
    5 have 1 - P singular, 2 fail by cross terms, 7 pass); for every fixed-point-free p with 1 - P invertible, the
    cross-term condition alpha_m.D alpha_j = 0 (j not in {m, p(m)}) holds iff every neighbour of a node i is p(i) or
    p^{-1}(i), iff the diagram is a cycle and p a rotation by one step; automorphisms keep the Kac labels; condition 3
    holds for the potential built from the partners exactly in the passing cases, and fails if the coefficients of
    g are not constant on an orbit.                                                                    [prp-ar-only]
 4. The a_r defect in components (h = 2..5): with mu_k the weights of the defining representation and phi_k =
    mu_k.phi, D alpha_j = -2 mu_{j+1}, D^T alpha_j = 2 mu_j, so f_j = e^{b(psi_j - phi_{j+1})}, g_j = e^{b(phi_j - psi_j)}
    at both couplings; the light-cone conditions in components, d_+(phi_k - psi_k) = -(m sigma/b)(f_k - f_{k-1}) and
    d_-(psi_k - phi_{k+1}) = (m/(b sigma))(g_k - g_{k+1}), and the field equation as the periodic 2D Toda lattice
    d_+ d_- phi_k = -(m^2/b)(e^{b(phi_k - phi_{k+1})} - e^{b(phi_{k-1} - phi_k)}).     [eq-dc-components, sec-defect-backlund-an]
 5. Unfreezing: if the light-cone conditions of check 2 hold at every x, then d_+ d_- phi = -grad V(phi) and
    d_+ d_- psi = -grad V(psi), from the second derivatives of f and g and the functional equation (a_2, a_3, a_4,
    both couplings, random points); an a_2 soliton pair satisfies the conditions at three values of x.   [sec-defect-backlund-an]
 6. The reflection rho alpha_j = alpha_{-j}: rho pi rho = pi^{-1}, rho D rho = D^T; the mirror defect 2(1 - pi^{-1})^{-1}
    is the type-I_- defect seen through rho on both sides; -D = (1 + E) pi^{-1}; and
    L^{I_+}(phi, psi; sigma) = L^{I_-}(rho phi, rho pi^{-1} psi; sigma), kinetic terms and potential, h = 2..5, both
    couplings (for h = 2, rho = 1 and pi^{-1} = -1); a type-I_- soliton pair of species h - a mapped by (rho, pi rho)
    solves the type-I_+ conditions and not the type-I_- ones; the species-a soliton crosses type I_+ with
    z_a(theta_d - theta) = omega^a z_{h-a}(theta - theta_d).                              [prp-dc-types, prp-an-defect-delay]
 7. The cyclic rotation conjugates both bulk Lax matrices by the cyclic shift, preserving the local hierarchy
    (h = 2..5, both couplings). Symmetries and the spin-2 charges: t -> -t maps the type-I_-(sigma) conditions to the type-I_+(1/sigma) ones;
    the composite with check 6 maps a type-I_-(sigma) soliton pair to a type-I_-(1/sigma) solution; the cubic
    invariant A(u) = 2 sum_k (mu_k.u)^3 of the spin-+-2 densities satisfies A(pi u) = A(u) and A(rho u) = -A(u).
                                                                                                         [prp-spin-two]
Runtime: a few seconds.
"""
import cmath
import itertools
import math

import numpy as np
import sympy as sp

OK = True


def report(name, ok):
    global OK
    OK = OK and bool(ok)
    print(("PASS " if ok else "FAIL ") + name)


rng = np.random.default_rng(38)
m, beta = 1.0, 0.7
COUPLINGS = {"real": beta, "imaginary": 1j * beta}


def ar_setup(h):
    """a_r^(1) roots alpha_j = e_j - e_{j+1}, weights mu_k = e_k - (1/h) sum e, pi e_k = e_{k-1}, in an orthonormal
    basis of the hyperplane sum = 0 of R^h (U maps coordinates to components)."""
    U = np.linalg.svd(np.eye(h) - np.ones((h, h)) / h)[0][:, :h - 1]
    e = np.eye(h)
    al = [U.T @ (e[j] - e[(j + 1) % h]) for j in range(h)]
    mu = [U.T @ e[k] for k in range(h)]
    pih = np.zeros((h, h))
    for k in range(h):
        pih[(k - 1) % h, k] = 1.0
    pi = U.T @ pih @ U
    D = 2 * np.linalg.inv(np.eye(h - 1) - pi)
    rho = np.array([al[(-j) % h] for j in range(1, h)]).T @ np.linalg.inv(np.array(al[1:]).T)
    return dict(h=h, r=h - 1, U=U, al=al, mu=mu, pi=pi, D=D, E=np.eye(h - 1) - D, rho=rho)


def Bd(S, phi, psi, sig, b, typ="-"):
    """The book's defect potential (eq-defect-potential) of type I_- or I_+ (D inside f_j, g_j unchanged)."""
    D = S["D"]
    f = sum(np.exp(b * a @ (D.T @ phi + D @ psi) / 2) - 1 for a in S["al"])
    g = sum(np.exp(b * a @ D @ (phi - psi) / 2) - 1 for a in S["al"])
    cf, cg = (sig, 1 / sig) if typ == "-" else (1 / sig, sig)
    return m / b**2 * (cf * f + cg * g)


def halves_grad(S, phi, psi, sig, b, typ="-"):
    """grad f(xi), grad g(eta) and Hessians, for B_d = f(xi) + g(eta)."""
    D, al = S["D"], S["al"]
    xi, eta = D.T @ phi + D @ psi, phi - psi
    cf, cg = (sig, 1 / sig) if typ == "-" else (1 / sig, sig)
    ys = [D.T @ a for a in al]
    gf = m / b**2 * cf * sum(b * a / 2 * np.exp(b * a @ xi / 2) for a in al)
    Hf = m / b**2 * cf * sum((b / 2)**2 * np.outer(a, a) * np.exp(b * a @ xi / 2) for a in al)
    gg = m / b**2 * cg * sum(b * y / 2 * np.exp(b * y @ eta / 2) for y in ys)
    Hg = m / b**2 * cg * sum((b / 2)**2 * np.outer(y, y) * np.exp(b * y @ eta / 2) for y in ys)
    return gf, gg, Hf, Hg


def V(S, phi, b):
    return m**2 / b**2 * sum(np.exp(b * a @ phi) - 1 for a in S["al"])


def gradV(S, phi, b):
    return m**2 / b * sum(a * np.exp(b * a @ phi) for a in S["al"])


# ---------------------------------------------------------------------------------------------------------
# 1. condition 2 in the coordinates xi, eta
eps = sp.Symbol("epsilon")
J = sp.Matrix([[0, 1], [-1, 0]])
Ec = eps * J
Dc = sp.eye(2) - Ec
p1, p2, s1, s2 = sp.symbols("p1 p2 s1 s2")
phv, psv = sp.Matrix([p1, p2]), sp.Matrix([s1, s2])
xi, eta = Dc.T * phv + Dc * psv, phv - psv
ok1 = sp.simplify((xi + Dc * eta) / 2 - phv) == sp.zeros(2, 1) and sp.simplify((xi - Dc.T * eta) / 2 - psv) == sp.zeros(2, 1)
fB = sp.exp(xi[0] - 2 * xi[1]) + sp.cos(xi[1]) * xi[0]**2
gB = sp.sin(eta[0]) * sp.exp(eta[1]) + eta[0]**3
B = fB + gB
X = Ec * sp.Matrix([B.diff(v) for v in (p1, p2)]) - Dc * sp.Matrix([B.diff(v) for v in (s1, s2)])
Y = -Dc.T * sp.Matrix([B.diff(v) for v in (p1, p2)]) - Ec * sp.Matrix([B.diff(v) for v in (s1, s2)])
Pd = gB - fB
ok1 &= sp.simplify(X - sp.Matrix([Pd.diff(v) for v in (p1, p2)])) == sp.zeros(2, 1)
ok1 &= sp.simplify(Y - sp.Matrix([Pd.diff(v) for v in (s1, s2)])) == sp.zeros(2, 1)
Bm = xi[0] * eta[1]
W = list(Ec * sp.Matrix([Bm.diff(v) for v in (p1, p2)]) - Dc * sp.Matrix([Bm.diff(v) for v in (s1, s2)])) + \
    list(-Dc.T * sp.Matrix([Bm.diff(v) for v in (p1, p2)]) - Ec * sp.Matrix([Bm.diff(v) for v in (s1, s2)]))
vs = (p1, p2, s1, s2)
ok1 &= any(sp.simplify(W[i].diff(vs[j]) - W[j].diff(vs[i])) != 0 for i in range(4) for j in range(i))
report("1. condition 2 in xi = D^T phi + D psi, eta = phi - psi: B = f(xi) + g(eta) satisfies it with P_d = g - f for "
       "arbitrary f, g; a mixed term xi_1 eta_2 does not", ok1)

# ---------------------------------------------------------------------------------------------------------
# 2. light-cone form and the functional equation
ok2 = True
for _ in range(4):
    A = rng.normal(size=(3, 3))
    E = 0.6 * (A - A.T)
    D = np.eye(3) - E
    gf, gg = rng.normal(size=3), rng.normal(size=3)
    Gp, Gq = D @ gf + gg, D.T @ gf - gg
    pt, qt = rng.normal(size=3), rng.normal(size=3)
    px, qx = E @ pt + D @ qt - Gp, D.T @ pt - E @ qt + Gq
    ok2 &= np.allclose((pt + px) - (qt + qx), -2 * gf) and np.allclose(D.T @ (pt - px) + D @ (qt - qx), 2 * gg)
    ok2 &= np.isclose(0.5 * Gp @ Gp - 0.5 * Gq @ Gq, 2 * gf @ gg)
for h in (3, 4):
    S = ar_setup(h)
    for b in COUPLINGS.values():
        phi, psi = rng.normal(size=h - 1) * 0.5, rng.normal(size=h - 1) * 0.5
        gf, gg, _, _ = halves_grad(S, phi, psi, 0.6, b)
        D = S["D"]
        num = np.array([(Bd(S, phi + 1e-6 * e, psi, 0.6, b) - Bd(S, phi - 1e-6 * e, psi, 0.6, b)) / 2e-6
                        for e in np.eye(h - 1)])
        ok2 &= np.allclose(num, D @ gf + gg, atol=1e-6)
        ok2 &= np.isclose(2 * gf @ gg, V(S, phi, b) - V(S, psi, b))
report("2. the defect conditions take the light-cone form d_+(phi - psi) = -2 grad f, d_-(D^T phi + D psi) = 2 grad g; "
       "condition 3 is 2 grad f.grad g = V(phi) - V(psi); the a_r halves at both couplings", ok2)


# ---------------------------------------------------------------------------------------------------------
# 3. the scan and the neighbour property
def gram(lengths, bonds, special=None):
    G = np.diag(np.array(lengths, float))
    for i, j in bonds:
        G[i, j] = G[j, i] = -max(lengths[i], lengths[j]) / 2
    return np.array(special, float) if special is not None else G


def chain(n):
    return [(i, i + 1) for i in range(n - 1)]


DIAGRAMS = {
    "a1^(1)": gram([2, 2], [], special=[[2, -2], [-2, 2]]),
    "a2^(1)": gram([2] * 3, [(0, 1), (1, 2), (2, 0)]),
    "a3^(1)": gram([2] * 4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    "a4^(1)": gram([2] * 5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
    "b3^(1)": gram([2, 2, 2, 1], [(0, 2), (1, 2), (2, 3)]),
    "b4^(1)": gram([2, 2, 2, 2, 1], [(0, 2), (1, 2), (2, 3), (3, 4)]),
    "c2^(1)": gram([2, 1, 2], chain(3)),
    "c3^(1)": gram([2, 1, 1, 2], chain(4)),
    "c4^(1)": gram([2, 1, 1, 1, 2], chain(5)),
    "d4^(1)": gram([2] * 5, [(0, 2), (1, 2), (3, 2), (4, 2)]),
    "f4^(1)": gram([2, 2, 2, 1, 1], chain(5)),
    "g2^(1)": gram([2, 2, 2 / 3], chain(3)),
    "a2^(2)": gram([2, 0.5], chain(2)),
    "a4^(2)": gram([0.5, 1, 2], chain(3)),
    "a6^(2)": gram([0.5, 1, 1, 2], chain(4)),
    "a8^(2)": gram([0.5, 1, 1, 1, 2], chain(5)),
    "d3^(2)": gram([1, 2, 1], chain(3)),
    "a5^(2)": gram([1, 1, 1, 2], [(0, 2), (1, 2), (2, 3)]),
    "a7^(2)": gram([1, 1, 1, 1, 2], [(0, 2), (1, 2), (2, 3), (3, 4)]),
    "d4^(2)": gram([1, 2, 2, 1], chain(4)),
    "d5^(2)": gram([1, 2, 2, 2, 1], chain(5)),
    "e6^(2)": gram([1, 1, 1, 2, 2], chain(5)),
    "d4^(3)": gram([2 / 3, 2 / 3, 2], chain(3)),
}


def roots_from_gram(G):
    h = G.shape[0]
    w, Vv = np.linalg.eigh(G)
    n = Vv[:, np.argmin(np.abs(w))]
    n = n / n[np.argmin(np.abs(n))]
    i0 = int(np.argmin(n))
    rest = [i for i in range(h) if i != i0]
    L = np.linalg.cholesky(G[np.ix_(rest, rest)])
    al = [None] * h
    for k, i in enumerate(rest):
        al[i] = L[k]
    al[i0] = -sum(n[i] * al[i] for i in rest) / n[i0]
    return al, n, rest


def cond3_residual(al, n, D, dcoef):
    """f = sum_i c_i e^{alpha_i.xi/2}, g = sum_i d_i e^{D^T alpha_i.eta/2}, c_i fixed by the e^{alpha_i.phi} terms."""
    h, r = len(al), len(al[0])
    res = 0.0
    for _ in range(3):
        phi, psi = rng.normal(size=r) * 0.6, rng.normal(size=r) * 0.6
        xi_, eta_ = D.T @ phi + D @ psi, phi - psi
        c = [2 * m**2 * n[i] / ((al[i] @ al[i]) * dcoef[i]) for i in range(h)]
        gf = sum(c[i] * al[i] / 2 * np.exp(al[i] @ xi_ / 2) for i in range(h))
        gg = sum(dcoef[i] * (D.T @ al[i]) / 2 * np.exp((D.T @ al[i]) @ eta_ / 2) for i in range(h))
        rhs = m**2 * sum(n[j] * (np.exp(al[j] @ phi) - np.exp(al[j] @ psi)) for j in range(h))
        res = max(res, abs(2 * gf @ gg - rhs))
    return res


counts = dict(fixes=0, singular=0, crossterm=0, passes=0)
ok3, total, passing = True, 0, set()
for name, G in DIAGRAMS.items():
    h = G.shape[0]
    al, n, rest = roots_from_gram(G)
    ok3 &= all(x > 0 for x in n) and np.allclose(n, np.round(n)) and np.linalg.matrix_rank(G, tol=1e-9) == h - 1
    nbrs = {i: [j for j in range(h) if j != i and abs(G[i, j]) > 1e-12] for i in range(h)}
    A = np.array([al[i] for i in rest]).T
    auts = [p for p in itertools.permutations(range(h))
            if all(abs(G[p[i], p[j]] - G[i, j]) < 1e-12 for i in range(h) for j in range(h))]
    total += len(auts)
    for p in auts:
        ok3 &= all(abs(n[p[i]] - n[i]) < 1e-9 for i in range(h))
        P = np.array([al[p[i]] for i in rest]).T @ np.linalg.inv(A)
        if any(p[i] == i for i in range(h)):
            counts["fixes"] += 1
            continue
        if abs(np.linalg.det(np.eye(h - 1) - P)) < 1e-9:
            counts["singular"] += 1
            continue
        D = 2 * np.linalg.inv(np.eye(h - 1) - P)
        pinv = {p[i]: i for i in range(h)}
        cross = all(abs(al[i] @ D @ al[j]) < 1e-9 for i in range(h) for j in range(h) if j not in (i, p[i]))
        nbr = all(set(nbrs[i]) <= {p[i], pinv[i]} for i in range(h))
        cyc = all(p[i] in nbrs[i] for i in range(h)) and all(len(nbrs[i]) == min(2, h - 1) for i in range(h))
        res = cond3_residual(al, n, D, [1.0] * h)
        ok3 &= (cross == nbr == cyc) and ((res < 1e-9) == cross)
        if cross:
            counts["passes"] += 1
            passing.add(name)
            # coefficients of g must be constant on the orbit
            ok3 &= cond3_residual(al, n, D, [1.0 + 0.3 * (i == 0) for i in range(h)]) > 1e-3
        else:
            counts["crossterm"] += 1
ok3 &= total == 78 and counts == dict(fixes=64, singular=5, crossterm=2, passes=7)
ok3 &= passing == {"a1^(1)", "a2^(1)", "a3^(1)", "a4^(1)"}
report(f"3. scan of 23 diagrams and {total} automorphisms {counts}: the cross-term condition holds iff every neighbour of "
       "i is p(i) or p^{-1}(i) iff the diagram is a cycle with p a one-step rotation; condition 3 holds exactly then, "
       "and needs the coefficients of g constant on the orbit; automorphisms keep the Kac labels", ok3)

# ---------------------------------------------------------------------------------------------------------
# 4. components
ok4 = True
for h in (2, 3, 4, 5):
    S = ar_setup(h)
    D, al, mu, U = S["D"], S["al"], S["mu"], S["U"]
    ok4 &= all(np.allclose(D @ al[j], -2 * mu[(j + 1) % h]) and np.allclose(D.T @ al[j], 2 * mu[j]) for j in range(h))
    for b in COUPLINGS.values():
        phi, psi = rng.normal(size=h - 1), rng.normal(size=h - 1)
        pc, qc = U @ phi, U @ psi
        ok4 &= all(np.isclose(np.exp(b * al[j] @ (D.T @ phi + D @ psi) / 2), np.exp(b * (qc[j] - pc[(j + 1) % h])))
                   for j in range(h))
        ok4 &= all(np.isclose(np.exp(b * al[j] @ D @ (phi - psi) / 2), np.exp(b * (pc[j] - qc[j]))) for j in range(h))
        # the light-cone conditions in components, and the two-dimensional Toda lattice
        sig = 0.65
        gf, gg, _, _ = halves_grad(S, phi, psi, sig, b)
        f = [np.exp(b * (qc[j] - pc[(j + 1) % h])) for j in range(h)]
        g = [np.exp(b * (pc[j] - qc[j])) for j in range(h)]
        lhs1 = U @ (-2 * gf)                                  # components of d_+(phi - psi)
        ok4 &= all(np.isclose(lhs1[k], -m * sig / b * (f[k] - f[(k - 1) % h])) for k in range(h))
        lhs2 = [al[k] @ (2 * gg) / 2 for k in range(h)]       # d_-(psi_k - phi_{k+1})
        ok4 &= all(np.isclose(lhs2[k], m / (b * sig) * (g[k] - g[(k + 1) % h])) for k in range(h))
        toda = U @ (-gradV(S, phi, b))
        ok4 &= all(np.isclose(toda[k], -m**2 / b * (np.exp(b * (pc[k] - pc[(k + 1) % h]))
                                                    - np.exp(b * (pc[(k - 1) % h] - pc[k])))) for k in range(h))
report("4. components: D alpha_j = -2 mu_{j+1}, D^T alpha_j = 2 mu_j, f_j = e^{b(psi_j - phi_{j+1})}, "
       "g_j = e^{b(phi_j - psi_j)}; the light-cone conditions d_+(phi_k - psi_k) = -(m sigma/b)(f_k - f_{k-1}), "
       "d_-(psi_k - phi_{k+1}) = (m/(b sigma))(g_k - g_{k+1}) and the 2D Toda lattice, h = 2..5, both couplings", ok4)


# ---------------------------------------------------------------------------------------------------------
# 5. unfreezing
def soliton(S, a, th, c0, fac, x, t):
    """Book soliton at imaginary coupling, phi = (i/beta) sum_j alpha_j ln tau_j, and d_+ phi, d_- phi."""
    h = S["h"]
    om = cmath.exp(2j * math.pi / h)
    ma = 2 * m * math.sin(math.pi * a / h)
    E = fac * cmath.exp(ma * (x * cmath.cosh(th) - t * cmath.sinh(th)) + c0)
    tau = [1 + om**(j * a) * E for j in range(h)]
    w = [om**(j * a) * E / tau[j] for j in range(h)]
    al = S["al"]
    phi = 1j / beta * sum(al[j] * cmath.log(tau[j]) for j in range(h))
    dp = 1j / beta * sum(al[j] * w[j] * ma * cmath.exp(-th) for j in range(h))
    dm = 1j / beta * sum(al[j] * w[j] * (-ma * cmath.exp(th)) for j in range(h))
    return phi, dp, dm


def residual(S, ph, phd, phm, ps, psd, psm, sig, b, typ="-"):
    """Residual of the book's defect conditions (matrix form) for complex fields."""
    D, E = S["D"], S["E"]
    Dk, Ek = (D, E) if typ == "-" else (-D, -E)
    pt, px = (phd + phm) / 2, (phd - phm) / 2
    qt, qx = (psd + psm) / 2, (psd - psm) / 2
    gf, gg, _, _ = halves_grad(S, ph, ps, sig, b, typ)
    Gp, Gq = D @ gf + gg, D.T @ gf - gg
    return max(np.max(np.abs(px - (Ek @ pt + Dk @ qt - Gp))), np.max(np.abs(qx - (Dk.T @ pt - Ek @ qt + Gq))))


def z_a(u, a, h):
    g = math.pi * a / h
    return (cmath.exp(-u) - 1j * cmath.exp(1j * g)) / (cmath.exp(-u) - 1j * cmath.exp(-1j * g))


ok5 = True
for h in (3, 4, 5):
    S = ar_setup(h)
    D = S["D"]
    for b in COUPLINGS.values():
        for _ in range(2):
            phi, psi = rng.normal(size=h - 1) * 0.5, rng.normal(size=h - 1) * 0.5
            gf, gg, Hf, Hg = halves_grad(S, phi, psi, 0.8, b)
            box_eta, box_xi = -4 * Hf @ gg, -4 * Hg @ gf
            ok5 &= np.allclose((box_xi + D @ box_eta) / 2, -gradV(S, phi, b))
            ok5 &= np.allclose((box_xi - D.T @ box_eta) / 2, -gradV(S, psi, b))
S = ar_setup(3)
th, thd = 0.5, 0.2
res = max(residual(S, *soliton(S, 1, th, 0.1 + 0.2j, 1.0, x, t), *soliton(S, 1, th, 0.1 + 0.2j, z_a(th - thd, 1, 3), x, t),
                   math.exp(-thd), 1j * beta) for x in (-0.7, 0.0, 1.3) for t in (-1.0, 0.5))
ok5 &= res < 1e-11
report("5. unfreezing: the light-cone conditions imposed everywhere give d_+d_- phi = -grad V(phi) and "
       "d_+d_- psi = -grad V(psi) (a_2, a_3, a_4, both couplings); an a_2 soliton pair satisfies them at x = -0.7, 0, 1.3",
       ok5)

# ---------------------------------------------------------------------------------------------------------
# 6. the reflection relation between the two types
ok6 = True
for h in (2, 3, 4, 5):
    S = ar_setup(h)
    D, E, pi, rho, rr = S["D"], S["E"], S["pi"], S["rho"], S["r"]
    pinv = np.linalg.inv(pi)
    ok6 &= np.allclose(rho.T @ rho, np.eye(rr)) and np.allclose(rho @ pi @ rho, pinv) and np.allclose(rho @ D @ rho, D.T)
    ok6 &= np.allclose(2 * np.linalg.inv(np.eye(rr) - pinv), D.T) and np.allclose(-D, (np.eye(rr) + E) @ pinv)
    A1, B1 = rho, rho @ pinv
    ok6 &= np.allclose(A1.T @ E @ A1, -E) and np.allclose(A1.T @ D @ B1, -D) and np.allclose(B1.T @ E @ B1, -E)
    if h == 2:
        ok6 &= np.allclose(rho, np.eye(1)) and np.allclose(B1, -np.eye(1))
    for b in COUPLINGS.values():
        for _ in range(2):
            phi, psi = rng.normal(size=rr) * 0.6, rng.normal(size=rr) * 0.6
            ok6 &= np.isclose(Bd(S, A1 @ phi, B1 @ psi, 0.55, b, "-"), Bd(S, phi, psi, 0.55, b, "+"))
            # mirror: the defect with D^T, built by the same formula, is the type-I_- defect seen through rho
            Dm = D.T
            Bm_ = m / b**2 * sum(0.55 * (np.exp(b * a @ (Dm.T @ phi + Dm @ psi) / 2) - 1)
                                 + (np.exp(b * a @ Dm @ (phi - psi) / 2) - 1) / 0.55 for a in S["al"])
            ok6 &= np.isclose(Bd(S, rho @ phi, rho @ psi, 0.55, b, "-"), Bm_)
res_p, res_ctrl, res_book = 0.0, 1.0, 0.0
for h in (3, 4):
    S = ar_setup(h)
    U, rho, pi = S["U"], S["rho"], S["pi"]
    for a in range(1, h):
        th, thd = 0.6, 0.15
        sig = math.exp(-thd)
        worst_ctrl = 0.0
        for t in (-0.8, 0.3, 1.5):
            Ph = soliton(S, h - a, th, 0.2 + 0.1j, 1.0, 0.0, t)
            Ps = soliton(S, h - a, th, 0.2 + 0.1j, z_a(th - thd, h - a, h), 0.0, t)
            mph = [rho @ v for v in Ph]
            mps = [pi @ rho @ v for v in Ps]
            res_p = max(res_p, residual(S, *mph, *mps, sig, 1j * beta, "+"))
            worst_ctrl = max(worst_ctrl, residual(S, *mph, *mps, sig, 1j * beta, "-"))
            Qh = soliton(S, a, th, 0.2 + 0.1j, 1.0, 0.0, t)
            Qs = soliton(S, a, th, 0.2 + 0.1j, z_a(thd - th, a, h), 0.0, t)
            res_book = max(res_book, residual(S, *Qh, *Qs, sig, 1j * beta, "+"))
        res_ctrl = min(res_ctrl, worst_ctrl)
ok6 &= res_p < 1e-11 and res_book < 1e-11 and res_ctrl > 1e-3
ok6 &= all(cmath.isclose(z_a(-u, a, h), cmath.exp(2j * math.pi * a / h) * z_a(u, h - a, h))
           for h in (2, 3, 4, 5) for a in range(1, h) for u in (-1.2, 0.4, 0.8 + 0.3j))
report("6. reflection: rho pi rho = pi^{-1}, rho D rho = D^T (the mirror defect), -D = (1 + E) pi^{-1}; "
       "L^{I+}(phi, psi; sigma) = L^{I-}(rho phi, rho pi^{-1} psi; sigma), kinetic and potential, h = 2..5, both "
       "couplings; mapped type-I_- soliton pairs solve type I_+ (not I_-); z_a(theta_d - theta) = omega^a z_{h-a}", ok6)

# ---------------------------------------------------------------------------------------------------------
# 7. cyclic invariance of the bulk hierarchy, time reversal, and the cubic invariant
def bulk_lax_pair(S, field, plus, minus, b, lam):
    """eq-lax-pair in the defining representation, with c_j = 1."""
    h, U = S["h"], S["U"]
    Ap = np.diag(b * (U @ plus) / 2).astype(complex)
    Am = np.diag(-b * (U @ minus) / 2).astype(complex)
    for j, alpha in enumerate(S["al"]):
        weight = m * np.exp(b * (alpha @ field) / 2)
        Ap[j, (j + 1) % h] += lam * weight
        Am[(j + 1) % h, j] -= weight / lam
    return Ap, Am


ok7 = True
for h in (2, 3, 4, 5):
    S = ar_setup(h)
    shift = np.roll(np.eye(h), -1, axis=0)
    jets = rng.normal(size=(3, h - 1)) * 0.5
    rotated = [S["pi"] @ v for v in jets]
    for b in COUPLINGS.values():
        original = bulk_lax_pair(S, *jets, b, 0.8 + 0.3j)
        transformed = bulk_lax_pair(S, *rotated, b, 0.8 + 0.3j)
        ok7 &= all(np.allclose(B, shift @ A @ shift.T) for A, B in zip(original, transformed))
for h in (3, 4):
    S = ar_setup(h)
    D, E = S["D"], S["E"]
    for b in COUPLINGS.values():
        phi, psi = rng.normal(size=h - 1) * 0.5, rng.normal(size=h - 1) * 0.5
        sig = 0.7
        pt, qt = rng.normal(size=h - 1), rng.normal(size=h - 1)
        gf, gg, _, _ = halves_grad(S, phi, psi, sig, b, "-")
        px, qx = E @ pt + D @ qt - (D @ gf + gg), D.T @ pt - E @ qt + (D.T @ gf - gg)
        # t -> -t: velocities change sign; the same x-derivatives must satisfy type I_+ with 1/sigma
        gf2, gg2, _, _ = halves_grad(S, phi, psi, 1 / sig, b, "+")
        ok7 &= np.allclose(px, -E @ (-pt) - D @ (-qt) - (D @ gf2 + gg2))
        ok7 &= np.allclose(qx, -D.T @ (-pt) + E @ (-qt) + (D.T @ gf2 - gg2))
    # composite: (Phi, Psi)(x, t) of type I_-(sigma) -> (rho Phi, rho pi^{-1} Psi)(x, -t) of type I_-(1/sigma)
    rho, pi = S["rho"], S["pi"]
    th, thd = 0.5, 0.3
    sig = math.exp(-thd)
    worst = 0.0
    for t in (-0.6, 0.4, 1.1):
        Ph = soliton(S, 1, th, 0.3, 1.0, 0.0, -t)
        Ps = soliton(S, 1, th, 0.3, z_a(th - thd, 1, h), 0.0, -t)
        # under t -> -t, d_+ -> -d_-, d_- -> -d_+
        nph = (rho @ Ph[0], -rho @ Ph[2], -rho @ Ph[1])
        nps = (rho @ np.linalg.inv(pi) @ Ps[0], -rho @ np.linalg.inv(pi) @ Ps[2], -rho @ np.linalg.inv(pi) @ Ps[1])
        worst = max(worst, residual(S, *nph, *nps, 1 / sig, 1j * beta, "-"))
    ok7 &= worst < 1e-11
    for _ in range(3):
        u = rng.normal(size=h - 1)
        Acub = lambda v: 2 * sum((mu @ v)**3 for mu in S["mu"])
        ok7 &= np.isclose(Acub(pi @ u), Acub(u)) and np.isclose(Acub(rho @ u), -Acub(u))
report("7. pi conjugates both bulk Lax matrices by the cyclic shift (h = 2..5, both couplings); t -> -t maps the type-I_-(sigma) conditions to type I_+(1/sigma); with the reflection it maps a type-I_-(sigma) "
       "soliton pair to a type-I_-(1/sigma) solution; A(pi u) = A(u), A(rho u) = -A(u) for A(u) = 2 sum (mu_k.u)^3", ok7)

raise SystemExit(0 if OK else 1)
