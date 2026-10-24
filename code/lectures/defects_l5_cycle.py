"""Checks for Lecture 5 of the defect lecture notes, "Why the cycle?"
(lecture-notes/defects/05-why-the-cycle.qmd).

Conventions of the lectures: dimensionless fields phi in R^r on x < 0 and psi in R^r on x > 0, action
(1/beta^2) int {int (1/2 d phi . d phi - V(phi)) + int (same for psi) + L_d}, with the real-coupling potential
    V(phi) = m^2 sum_{j=0}^r n_j (e^{alpha_j . phi} - 1),
alpha_j the affine simple roots (longest roots of squared length 2), n_j the Kac labels (n_0 = 1). The defect term is
    L_d = 1/2 phi . E phi_t + phi . D psi_t + 1/2 psi . F psi_t - B(phi, psi),     E, F antisymmetric,
with the defect conditions phi_x = E phi_t + D psi_t - grad_phi B and psi_x = D^T phi_t - F psi_t + grad_psi B at x = 0.
After condition 1, D = 1 - E and F = E; then xi = D^T phi + D psi, eta = phi - psi, B = f(xi) + g(eta), and the
light-cone form is d_+ eta = -2 grad f(xi), d_- xi = 2 grad g(eta), d_pm = d_t pm d_x.
For a_r^(1) (h = r + 1): alpha_j = e_j - e_{j+1} in the hyperplane sum = 0 of R^h, mu_k = e_k - (1/h) sum e,
pi e_k = e_{k-1} (so pi alpha_j = alpha_{j-1}), D = 2 (1 - pi)^{-1}, and
    B = m sum_j [sigma (f_j - 1) + sigma^{-1} (g_j - 1)],   f_j = e^{psi_j - phi_{j+1}},   g_j = e^{phi_j - psi_j},
with phi_k the components in R^h and sigma = e^{-theta_d}. Imaginary coupling is phi -> i phi; its solitons are
complex solutions of the same equations, phi_k = ln(tau_{k-1}/tau_k), tau_k = 1 + omega^{k a} E,
E = exp(m_a (x cosh theta - t sinh theta) + xi_0), omega = e^{2 pi i/h}, m_a = 2 m sin(pi a/h).
Written from these definitions, independently of code/defect_classical_an.py and defect_solitons_an.py.
One PASS or FAIL line per check; exit status 1 on any FAIL. Units m = 1 in the numerics.

  (a) Variation: the Euler-Lagrange expressions of L_d with constant matrices give the stated defect conditions;
      E + B/beta^2 is conserved for every B (symbolic, r = 2, generic D, antisymmetric E and F).
  (b) The momentum balance after substitution is a polynomial in the velocities; its quadratic part vanishes for
      D = (1 - E) O, F = O^T E O (O a rotation), its linear part is -(phi_t . X + psi_t . Y) with
      X = E grad_phi B - D grad_psi B, Y = -D^T grad_phi B - F grad_psi B, and its velocity-free part is
      1/2 |grad_psi B|^2 - 1/2 |grad_phi B|^2 + V(phi) - V(psi) (symbolic, r = 2).
  (c) Condition 1: O = (1 - E)^{-1} D is orthogonal when D D^T = 1 - E^2; then E D = D F forces F = O^T E O and
      D^T D = 1 - F^2 follows; psi' = O psi gives D' = 1 - E, F' = E; a non-orthogonal O fails (r = 3, numeric).
  (d) Condition 2 in the coordinates xi, eta: phi = (xi + D eta)/2, psi = (xi - D^T eta)/2; for B = f(xi) + g(eta)
      the vector field (X, Y) is the gradient of P_d = g - f; a mixed term xi_1 eta_2 gives a curl (symbolic, r = 2).
  (e) The light-cone form follows from the conditions: d_+ eta = -2 grad f, d_- xi = 2 grad g (r = 3, numeric).
  (f) Condition 3 becomes 2 grad f(xi) . grad g(eta) = V(phi) - V(psi), because D is normal (numeric).
  (g) Exponentials: (x, y) -> (D x + y, D^T x - y) is invertible; a product is e^{alpha . phi} only for
      (alpha, D^T alpha) and e^{alpha . psi} only for (alpha, -D alpha); alpha . D alpha = alpha . D^T alpha = |alpha|^2.
  (h) Partners: P = -D^{-1} D^T is orthogonal (the Cayley transform of E), D (1 - P) = 2, and for a_r
      D^T alpha_i = -D alpha_{i-1}.
  (i) The scan: 23 affine diagrams of rank <= 4 (untwisted and twisted) and their 78 automorphisms; 64 fix a
      node, 5 have 1 - P singular, 2 (pi^{+-2} of a_4^(1)) fail the cross-term condition, 7 pass (pi^{+-1} of
      a_2, a_3, a_4 and the exchange of a_1). For every fixed-point-free automorphism with 1 - P invertible:
      the cross-term condition holds exactly when every neighbour of i is p(i) or p^{-1}(i), exactly when the
      diagram is a cycle and p a rotation by one step; condition 3 for the potential built from the partners
      holds exactly in the passing cases; automorphisms preserve the Kac labels.
  (j) The a_r defect (h = 2..5): D alpha_j = -2 mu_{j+1}, D^T alpha_j = 2 mu_j, alpha_i . D alpha_j =
      2 delta_ij - 2 delta_{i,j+1}, D = 1 + i cot(pi k/h) on the eigenvectors of pi, D = the book's
      2 sum (lambda_i - lambda_{i+1}) lambda_i^T; f_j g_j = e^{alpha_j . phi}, f_j g_{j+1} = e^{alpha_j . psi};
      conditions 1-3 with P_d = m sum [sigma^{-1}(g_j - 1) - sigma (f_j - 1)], and condition 3 fails for a wrong
      coefficient; the matrix conditions are equivalent to the component light-cone equations; the light-cone
      split B +- P_d.
  (k) a_1^(1): the defect is Lecture 3's sinh-Gordon defect with Phi = 2 phi_1, m0 = 2 m, B_L3 = 2 B, and its
      continuation is the sine-Gordon defect.
  (l) Unfreezing: the light-cone system with any f, g obeying the functional equation gives box phi =
      -grad V(phi), box psi = -grad V(psi); for a_r these are the periodic two-dimensional Toda lattice; an a_2
      soliton pair solves the light-cone system at several x (a Baecklund pair).
  (m) Type I_+: (D, E, F) -> (-D, -E, -E), sigma <-> 1/sigma in B, satisfies condition 1 with O = pi^{-1};
      L_+(phi, psi; sigma) = L_-(rho phi, rho pi^{-1} psi; sigma) for the reflection rho: alpha_j -> alpha_{-j},
      kinetic and potential terms (h = 2..5); the mirror image D^T = 2 (1 - pi^{-1})^{-1} is L_d(rho phi, rho psi); a type-I_- soliton pair mapped by (rho, pi rho) solves the type-I_+
      conditions; z_a(-u) = omega^a z_{h-a}(u).
  (n) Linear waves: per mode, R_k = 0 and T_k = (e^{-theta} + sigma e^{i gamma_k})/(e^{-theta} + sigma e^{-i gamma_k})
      (sympy); the full vector linearized conditions for a_3, a_4 give the same, with Hessians of B by finite
      differences; |T_k| = 1, T_k(theta_d) = e^{i gamma_k}, T_k -> 1 and -> omega^k, the h = 2 reduction to
      Lecture 4's T; the large-sigma limit of the conditions is psi = pi phi on the vacuum branch.
  (o) Solitons: the species-a soliton with E -> z_a E solves the component light-cone equations for h = 2..5, all
      species, real and complex theta, at x = 0, -0.7, 1.3, with z_a(u) = (e^{-u} - i e^{i gamma_a})/
      (e^{-u} - i e^{-i gamma_a}), u = theta - theta_d, and fails for its conjugate; z_a(u) = T_a(u + theta_d - i pi/2);
      |z_a| > 1 on the real line; arg z_a = arg(e^{-2u} + e^{2 i gamma_a}); limits 1 and omega^a, and
      E -> omega^a E is psi = pi phi; the pole at u = -i (pi/2 - gamma_a) is real only for a = h/2.
  (p) Exercises: for a_2, E = -J/sqrt 3 (J the rotation by pi/2 in the orthonormal basis used), and P is the
      rotation by 2 pi/3; the a_3 reflection (0 1)(2 3) sends every node to a neighbour but has a neighbour outside
      {p(i), p^{-1}(i)}; a wave from the right has bar R_k = 0 and bar T_k(theta) = 1/T_k(-theta); the book's soliton
      -sum_j alpha_j ln tau_j has the components ln tau_{k-1} - ln tau_k, and solves the Hirota form
      tau_k d+d- tau_k - d+tau_k d-tau_k = m^2 (tau_{k-1} tau_{k+1} - tau_k^2); T_k has its zero at
      theta_d + i (pi - gamma_k) and its pole at theta_d - i (pi - gamma_k).

Runtime: under a minute.
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


rng = np.random.default_rng(5)
m = 1.0


def rand_antisym(r):
    A = rng.normal(size=(r, r))
    return A - A.T


def rand_orth(r):
    Q, R = np.linalg.qr(rng.normal(size=(r, r)))
    return Q @ np.diag(np.sign(np.diag(R)))


# ---------------------------------------------------------------------------------------------------------
# (a) variation and energy, symbolic, r = 2
t = sp.Symbol('t')
r2 = 2
ph = sp.Matrix([sp.Function(f'phi{a}')(t) for a in range(r2)])
ps = sp.Matrix([sp.Function(f'psi{a}')(t) for a in range(r2)])
e12, f12 = sp.symbols('e12 f12')
Es = sp.Matrix([[0, e12], [-e12, 0]])
Fs = sp.Matrix([[0, f12], [-f12, 0]])
Ds = sp.Matrix(2, 2, sp.symbols('d11 d12 d21 d22'))
b1, b2, b3 = sp.symbols('b1 b2 b3')
Bf = b1 * sp.cos(ph[0] + 2 * ps[1]) + b2 * sp.exp(ph[1] * ps[0]) + b3 * ph[0]**2 * ps[1] + sp.sin(ps[0] - ph[1])
Ld = (ph.T * Es * ph.diff(t))[0] / 2 + (ph.T * Ds * ps.diff(t))[0] + (ps.T * Fs * ps.diff(t))[0] / 2 - Bf


def EL(L, q):
    return sp.simplify(sp.diff(L, q) - sp.diff(sp.diff(L, q.diff(t)), t))


ELphi = sp.Matrix([EL(Ld, ph[a]) for a in range(r2)])
ELpsi = sp.Matrix([EL(Ld, ps[a]) for a in range(r2)])
gphi = sp.Matrix([sp.diff(Bf, ph[a]) for a in range(r2)])
gpsi = sp.Matrix([sp.diff(Bf, ps[a]) for a in range(r2)])
phx_claim = Es * ph.diff(t) + Ds * ps.diff(t) - gphi
psx_claim = Ds.T * ph.diff(t) - Fs * ps.diff(t) + gpsi
# coefficient of delta phi(0): -phi_x + EL_phi = 0; of delta psi(0): +psi_x + EL_psi = 0
ok_a = sp.simplify(ELphi - phx_claim) == sp.zeros(2, 1) and sp.simplify(-ELpsi - psx_claim) == sp.zeros(2, 1)
energy = (ph.diff(t).T * phx_claim)[0] - (ps.diff(t).T * psx_claim)[0] + sp.diff(Bf, t)
ok_a = ok_a and sp.simplify(sp.expand(energy)) == 0
report("(a) the defect term with constant matrices gives phi_x = E phi_t + D psi_t - grad_phi B, psi_x = D^T phi_t "
       "- F psi_t + grad_psi B, and E + B/beta^2 is conserved for every B", ok_a)

# ---------------------------------------------------------------------------------------------------------
# (b) the momentum polynomial, symbolic, r = 2
u1, u2, w1, w2, G1a, G1b, G2a, G2b, Vp, Vq = sp.symbols('u1 u2 w1 w2 G1a G1b G2a G2b Vp Vq')
ut = sp.Matrix([u1, u2])     # phi_t
wt = sp.Matrix([w1, w2])     # psi_t
G1 = sp.Matrix([G1a, G1b])   # grad_phi B
G2 = sp.Matrix([G2a, G2b])   # grad_psi B
phx = Es * ut + Ds * wt - G1
psx = Ds.T * ut - Fs * wt + G2
flux = (wt.dot(wt) / 2 + psx.dot(psx) / 2 - Vq) - (ut.dot(ut) / 2 + phx.dot(phx) / 2 - Vp)
flux = sp.expand(flux)
poly = sp.Poly(flux, u1, u2, w1, w2)
quad = sum(c * sp.prod([v**k for v, k in zip((u1, u2, w1, w2), mon)]) for mon, c in poly.terms() if sum(mon) == 2)
lin = sum(c * sp.prod([v**k for v, k in zip((u1, u2, w1, w2), mon)]) for mon, c in poly.terms() if sum(mon) == 1)
const = sum(c for mon, c in poly.terms() if sum(mon) == 0)
X = Es * G1 - Ds * G2
Yv = -Ds.T * G1 - Fs * G2
ok_b = sp.simplify(lin + ut.dot(X) + wt.dot(Yv)) == 0
ok_b &= sp.simplify(const - (G2.dot(G2) / 2 - G1.dot(G1) / 2 + Vp - Vq)) == 0
# quadratic part: vanishes for D = (1 - E) O, F = O^T E O with O a rotation or a reflection
ang = sp.Symbol('a')
for O in (sp.Matrix([[sp.cos(ang), -sp.sin(ang)], [sp.sin(ang), sp.cos(ang)]]),
          sp.Matrix([[sp.cos(ang), sp.sin(ang)], [sp.sin(ang), -sp.cos(ang)]])):
    sub = dict(zip(sp.symbols('d11 d12 d21 d22'), list((sp.eye(2) - Es) * O)))
    Fsub = sp.simplify(O.T * Es * O)
    q = quad.subs(sub).subs(f12, Fsub[0, 1])
    ok_b &= sp.simplify(q) == 0
# and does not vanish for a generic D
ok_b &= sp.simplify(quad.subs({'d11': 1, 'd12': 0, 'd21': 0, 'd22': 1, e12: sp.Rational(1, 3), f12: 0})) != 0
report("(b) the momentum balance is a polynomial in the velocities: quadratic part zero for D = (1-E)O, F = O^T E O, "
       "linear part -(phi_t.X + psi_t.Y), velocity-free part |grad_psi B|^2/2 - |grad_phi B|^2/2 + V(phi) - V(psi)",
       ok_b)

# ---------------------------------------------------------------------------------------------------------
# (c) condition 1, numeric, r = 3
r = 3
ok_c = True
for _ in range(5):
    E = rand_antisym(r)
    O = rand_orth(r)
    D = (np.eye(r) - E) @ O
    F = O.T @ E @ O
    I = np.eye(r)
    ok_c &= np.allclose(D @ D.T, I - E @ E) and np.allclose(D.T @ D, I - F @ F) and np.allclose(E @ D, D @ F)
    O2 = np.linalg.solve(I - E, D)
    ok_c &= np.allclose(O2.T @ O2, I)
    # F is forced by E D = D F
    ok_c &= np.allclose(np.linalg.solve(D, E @ D), F)
    # after psi' = O psi: phi . D psi_t = phi . (1 - E) psi'_t, 1/2 psi . F psi_t = 1/2 psi' . E psi'_t
    ok_c &= np.allclose(D @ O.T, I - E) and np.allclose(O @ F @ O.T, E)
    # a non-orthogonal O fails the first equation
    Ob = O @ np.diag([1.2, 1.0, 1.0])
    Db = (I - E) @ Ob
    ok_c &= not np.allclose(Db @ Db.T, I - E @ E)
report("(c) condition 1: D D^T = 1 - E^2 makes (1-E)^{-1} D orthogonal; E D = D F forces F = O^T E O; psi' = O psi "
       "gives D = 1 - E, F = E; a non-orthogonal O fails", ok_c)

# ---------------------------------------------------------------------------------------------------------
# (d) condition 2 in the coordinates xi, eta, symbolic, r = 2 with E = eps J
eps = sp.Symbol('epsilon')
J = sp.Matrix([[0, 1], [-1, 0]])
Ec = eps * J
Dc = sp.eye(2) - Ec
q1, q2, q3, q4 = sp.symbols('p1 p2 s1 s2')
phv = sp.Matrix([q1, q2])
psv = sp.Matrix([q3, q4])
xi = Dc.T * phv + Dc * psv
eta = phv - psv
ok_d = sp.simplify((xi + Dc * eta) / 2 - phv) == sp.zeros(2, 1)
ok_d &= sp.simplify((xi - Dc.T * eta) / 2 - psv) == sp.zeros(2, 1)
fB = sp.cos(xi[0]) * sp.exp(xi[1]) + xi[0]**3
gB = sp.sin(eta[0] * eta[1]) + eta[1]**2
B = fB + gB
gph = sp.Matrix([sp.diff(B, v) for v in (q1, q2)])
gps = sp.Matrix([sp.diff(B, v) for v in (q3, q4)])
Xd = Ec * gph - Dc * gps
Yd = -Dc.T * gph - Ec * gps
Pd = gB - fB
ok_d &= sp.simplify(Xd - sp.Matrix([sp.diff(Pd, v) for v in (q1, q2)])) == sp.zeros(2, 1)
ok_d &= sp.simplify(Yd - sp.Matrix([sp.diff(Pd, v) for v in (q3, q4)])) == sp.zeros(2, 1)
# a mixed term is not allowed: the one-form (X, Y) has a nonzero curl
Bm = xi[0] * eta[1]
gph = sp.Matrix([sp.diff(Bm, v) for v in (q1, q2)])
gps = sp.Matrix([sp.diff(Bm, v) for v in (q3, q4)])
W = list(Ec * gph - Dc * gps) + list(-Dc.T * gph - Ec * gps)
vars4 = (q1, q2, q3, q4)
curl = [sp.simplify(sp.diff(W[i], vars4[j]) - sp.diff(W[j], vars4[i])) for i in range(4) for j in range(i)]
ok_d &= any(c != 0 for c in curl)
report("(d) in xi = D^T phi + D psi, eta = phi - psi: phi = (xi + D eta)/2, psi = (xi - D^T eta)/2; B = f(xi) + g(eta) "
       "satisfies condition 2 with P_d = g - f, a mixed term xi_1 eta_2 does not", ok_d)

# ---------------------------------------------------------------------------------------------------------
# (e) light-cone form and (f) the functional equation, numeric, r = 3
ok_e, ok_f = True, True
for _ in range(5):
    E = rand_antisym(r) * 0.7
    D = np.eye(r) - E
    gf, gg = rng.normal(size=r), rng.normal(size=r)        # grad f(xi), grad g(eta)
    Gphi = D @ gf + gg
    Gpsi = D.T @ gf - gg
    pt, qt = rng.normal(size=r), rng.normal(size=r)
    px = E @ pt + D @ qt - Gphi
    qx = D.T @ pt - E @ qt + Gpsi
    dplus_eta = (pt + px) - (qt + qx)
    dminus_xi = D.T @ (pt - px) + D @ (qt - qx)
    ok_e &= np.allclose(dplus_eta, -2 * gf) and np.allclose(dminus_xi, 2 * gg)
    ok_f &= np.isclose(0.5 * Gphi @ Gphi - 0.5 * Gpsi @ Gpsi, 2 * gf @ gg)
    ok_f &= np.allclose(D @ D.T, D.T @ D)
report("(e) the defect conditions with D = 1 - E, F = E and B = f(xi) + g(eta) are d_+ eta = -2 grad f(xi), "
       "d_- xi = 2 grad g(eta)", ok_e)
report("(f) condition 3 becomes 2 grad f(xi).grad g(eta) = V(phi) - V(psi), because D is normal", ok_f)

# ---------------------------------------------------------------------------------------------------------
# (g) exponentials and (h) partners
ok_g, ok_h = True, True
for _ in range(5):
    E = rand_antisym(r) * 0.8
    D = np.eye(r) - E
    Mmap = np.block([[D, np.eye(r)], [D.T, -np.eye(r)]])      # (x, y) -> (D x + y, D^T x - y)
    ok_g &= abs(np.linalg.det(Mmap)) > 1e-8
    al = rng.normal(size=r)
    # solving D x + y = 2 alpha, D^T x - y = 0, and D x + y = 0, D^T x - y = 2 alpha
    s1 = np.linalg.solve(Mmap, np.concatenate([2 * al, np.zeros(r)]))
    s2 = np.linalg.solve(Mmap, np.concatenate([np.zeros(r), 2 * al]))
    ok_g &= np.allclose(s1, np.concatenate([al, D.T @ al])) and np.allclose(s2, np.concatenate([al, -D @ al]))
    ok_g &= np.isclose(al @ D @ al, al @ al) and np.isclose(al @ D.T @ al, al @ al)
    P = -np.linalg.solve(D, D.T)
    ok_h &= np.allclose(P.T @ P, np.eye(r)) and np.allclose(D @ (np.eye(r) - P), 2 * np.eye(r))
    ok_h &= np.allclose(P, -np.linalg.inv(np.eye(r) - E) @ (np.eye(r) + E))
report("(g) the exponent map (x, y) -> (Dx + y, D^T x - y) is invertible; e^{alpha.phi} only from (alpha, D^T alpha), "
       "e^{alpha.psi} only from (alpha, -D alpha); alpha.D alpha = |alpha|^2", ok_g)


# ---------------------------------------------------------------------------------------------------------
# a_r^(1) in R^h, in an orthonormal basis of the hyperplane
def ar_setup(h):
    Mh = np.eye(h) - np.ones((h, h)) / h
    U = np.linalg.svd(Mh)[0][:, :h - 1]
    e = np.eye(h)
    al = [U.T @ (e[j] - e[(j + 1) % h]) for j in range(h)]
    mu = [U.T @ e[k] for k in range(h)]
    pih = np.zeros((h, h))
    for k in range(h):
        pih[(k - 1) % h, k] = 1.0                                # pi e_k = e_{k-1}
    pi = U.T @ pih @ U
    D = 2 * np.linalg.inv(np.eye(h - 1) - pi)
    return dict(h=h, r=h - 1, U=U, al=al, mu=mu, pi=pi, D=D, E=np.eye(h - 1) - D)


for h in (3, 4, 5):
    S = ar_setup(h)
    D, al = S['D'], S['al']
    ok_h &= all(np.allclose(D.T @ al[i], -D @ al[(i - 1) % h]) for i in range(h))
    P = -np.linalg.solve(D, D.T)
    ok_h &= np.allclose(P, S['pi'])
report("(h) partners: P = -D^{-1} D^T is orthogonal, the Cayley transform -(1-E)^{-1}(1+E), D(1-P) = 2; for a_r, "
       "D^T alpha_i = -D alpha_{i-1} and P = pi", ok_h)


# ---------------------------------------------------------------------------------------------------------
# (i) the scan over affine Dynkin diagrams
def gram(lengths, bonds, special=None):
    h = len(lengths)
    G = np.diag(np.array(lengths, float))
    for i, j in bonds:
        G[i, j] = G[j, i] = -max(lengths[i], lengths[j]) / 2
    if special is not None:
        G = np.array(special, float)
    return G


def chain(n):
    return [(i, i + 1) for i in range(n - 1)]


DIAGRAMS = {
    'a1^(1)': gram([2, 2], [], special=[[2, -2], [-2, 2]]),
    'a2^(1)': gram([2] * 3, [(0, 1), (1, 2), (2, 0)]),
    'a3^(1)': gram([2] * 4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    'a4^(1)': gram([2] * 5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
    'b3^(1)': gram([2, 2, 2, 1], [(0, 2), (1, 2), (2, 3)]),
    'b4^(1)': gram([2, 2, 2, 2, 1], [(0, 2), (1, 2), (2, 3), (3, 4)]),
    'c2^(1)': gram([2, 1, 2], chain(3)),
    'c3^(1)': gram([2, 1, 1, 2], chain(4)),
    'c4^(1)': gram([2, 1, 1, 1, 2], chain(5)),
    'd4^(1)': gram([2] * 5, [(0, 2), (1, 2), (3, 2), (4, 2)]),
    'f4^(1)': gram([2, 2, 2, 1, 1], chain(5)),
    'g2^(1)': gram([2, 2, 2 / 3], chain(3)),
    'a2^(2)': gram([2, 0.5], chain(2)),
    'a4^(2)': gram([0.5, 1, 2], chain(3)),
    'a6^(2)': gram([0.5, 1, 1, 2], chain(4)),
    'a8^(2)': gram([0.5, 1, 1, 1, 2], chain(5)),
    'd3^(2)': gram([1, 2, 1], chain(3)),
    'a5^(2)': gram([1, 1, 1, 2], [(0, 2), (1, 2), (2, 3)]),
    'a7^(2)': gram([1, 1, 1, 1, 2], [(0, 2), (1, 2), (2, 3), (3, 4)]),
    'd4^(2)': gram([1, 2, 2, 1], chain(4)),
    'd5^(2)': gram([1, 2, 2, 2, 1], chain(5)),
    'e6^(2)': gram([1, 1, 1, 2, 2], chain(5)),
    'd4^(3)': gram([2 / 3, 2 / 3, 2], chain(3)),
}
KAC = {'b4^(1)': [1, 1, 2, 2, 2], 'c4^(1)': [1, 2, 2, 2, 1], 'd4^(1)': [1, 1, 2, 1, 1], 'f4^(1)': [1, 2, 3, 4, 2],
       'g2^(1)': [1, 2, 3], 'a4^(1)': [1] * 5}


def roots_from_gram(G):
    """Root vectors in R^r with Gram matrix G (h x h, rank r), the Kac labels (null vector, smallest label 1),
    and the index i0 of a node with label 1, whose root is -sum of the others weighted by their labels."""
    h = G.shape[0]
    w, V = np.linalg.eigh(G)
    n = V[:, np.argmin(np.abs(w))]
    n = n / n[np.argmin(np.abs(n))]
    i0 = int(np.argmin(n))
    rest = [i for i in range(h) if i != i0]
    L = np.linalg.cholesky(G[np.ix_(rest, rest)])
    al = [None] * h
    for k, i in enumerate(rest):
        al[i] = L[k]
    al[i0] = -sum(n[i] * al[i] for i in rest) / n[i0]
    return al, n, i0


def automorphisms(G):
    h = G.shape[0]
    out = []
    for p in itertools.permutations(range(h)):
        if all(abs(G[p[i], p[j]] - G[i, j]) < 1e-12 for i in range(h) for j in range(h)):
            out.append(p)
    return out


def cond3_residual(al, n, D, p):
    """Residual of condition 3 for f = sum c_i e^{alpha_i.xi/2}, g = sum d e^{D^T alpha_i . eta/2}, d = 1."""
    h = len(al)
    rr = len(al[0])
    res = 0.0
    for _ in range(4):
        phi, psi = rng.normal(size=rr) * 0.6, rng.normal(size=rr) * 0.6
        xi, eta = D.T @ phi + D @ psi, phi - psi
        c = [2 * m**2 * n[i] / (al[i] @ al[i]) for i in range(h)]
        gf = sum(c[i] * al[i] / 2 * np.exp(al[i] @ xi / 2) for i in range(h))
        gg = sum((D.T @ al[i]) / 2 * np.exp((D.T @ al[i]) @ eta / 2) for i in range(h))
        lhs = 2 * gf @ gg
        rhs = m**2 * sum(n[j] * (np.exp(al[j] @ phi) - np.exp(al[j] @ psi)) for j in range(h))
        res = max(res, abs(lhs - rhs))
    return res


counts = dict(fixes=0, singular=0, crossterm=0, passes=0)
passing = []
ok_i = True
total = 0
for name, G in DIAGRAMS.items():
    h = G.shape[0]
    al, n, i0 = roots_from_gram(G)
    rest = [i for i in range(h) if i != i0]
    ok_i &= all(abs(al[i] @ al[j] - G[i, j]) < 1e-9 for i in range(h) for j in range(h))
    ok_i &= abs(np.linalg.det(G)) < 1e-9 and np.linalg.matrix_rank(G, tol=1e-9) == h - 1
    ok_i &= all(x > 0 for x in n) and np.allclose(n, np.round(n))
    if name in KAC:
        ok_i &= np.allclose(n, KAC[name])
    A = np.array([al[i] for i in rest]).T
    nbrs = {i: [j for j in range(h) if j != i and abs(G[i, j]) > 1e-12] for i in range(h)}
    auts = automorphisms(G)
    total += len(auts)
    for p in auts:
        ok_i &= all(abs(n[p[i]] - n[i]) < 1e-9 for i in range(h))       # automorphisms keep the Kac labels
        Ap = np.array([al[p[i]] for i in rest]).T
        P = Ap @ np.linalg.inv(A)
        ok_i &= np.allclose(P @ al[i0], al[p[i0]]) and np.allclose(P.T @ P, np.eye(h - 1))
        if any(p[i] == i for i in range(h)):
            counts['fixes'] += 1
            continue
        if abs(np.linalg.det(np.eye(h - 1) - P)) < 1e-9:
            counts['singular'] += 1
            continue
        D = 2 * np.linalg.inv(np.eye(h - 1) - P)
        pinv = {p[i]: i for i in range(h)}
        cross_ok = all(abs(al[i] @ D @ al[j]) < 1e-9 for i in range(h) for j in range(h) if j not in (i, p[i]))
        nbr_ok = all(set(nbrs[i]) <= {p[i], pinv[i]} for i in range(h))
        # a cycle with p a rotation by one step: i ~ p(i) for all i, and h = 2 or every node has two neighbours
        cyc = all(p[i] in nbrs[i] for i in range(h)) and all(len(nbrs[i]) == min(2, h - 1) for i in range(h))
        res3 = cond3_residual(al, n, D, p)
        ok_i &= (cross_ok == nbr_ok == cyc) and ((res3 < 1e-9) == cross_ok)
        ok_i &= np.allclose(D @ (np.eye(h - 1) - P), 2 * np.eye(h - 1))
        if cross_ok:
            counts['passes'] += 1
            passing.append((name, p))
        else:
            counts['crossterm'] += 1
ok_i &= len(DIAGRAMS) == 23 and total == 78
ok_i &= counts == dict(fixes=64, singular=5, crossterm=2, passes=7)
ok_i &= sorted(set(nm for nm, _ in passing)) == ['a1^(1)', 'a2^(1)', 'a3^(1)', 'a4^(1)']
report(f"(i) scan of 23 affine diagrams, {total} automorphisms: {counts}; for fixed-point-free p with 1-P invertible, "
       "cross terms vanish iff every neighbour of i is p(i) or p^{-1}(i) iff the diagram is a cycle with p a "
       "one-step rotation, iff condition 3 holds; automorphisms keep the Kac labels", ok_i)


# ---------------------------------------------------------------------------------------------------------
# (j) the a_r defect
def comps(S, v):
    """Components of v (orthonormal hyperplane coordinates) in R^h."""
    return S['U'] @ v


def B_ar(S, phi, psi, sig, typ='-'):
    h = S['h']
    pc, qc = comps(S, phi), comps(S, psi)
    f = sum(np.exp(qc[j] - pc[(j + 1) % h]) - 1 for j in range(h))
    g = sum(np.exp(pc[j] - qc[j]) - 1 for j in range(h))
    return m * (sig * f + g / sig) if typ == '-' else m * (f / sig + sig * g)


def grads_ar(S, phi, psi, sig, typ='-'):
    """Analytic gradients of B_ar, written from B = f(xi) + g(eta)."""
    h, D, al, mu = S['h'], S['D'], S['al'], S['mu']
    xi, eta = D.T @ phi + D @ psi, phi - psi
    cf, cg = (m * sig, m / sig) if typ == '-' else (m / sig, m * sig)
    gf = cf * sum(al[j] / 2 * np.exp(al[j] @ xi / 2) for j in range(h))
    gg = cg * sum(mu[j] * np.exp(mu[j] @ eta) for j in range(h))
    return D @ gf + gg, D.T @ gf - gg, gf, gg


ok_j = True
for h in (2, 3, 4, 5):
    S = ar_setup(h)
    D, al, mu, rr = S['D'], S['al'], S['mu'], S['r']
    ok_j &= all(np.allclose(D @ al[j], -2 * mu[(j + 1) % h]) and np.allclose(D.T @ al[j], 2 * mu[j]) for j in range(h))
    ok_j &= all(np.isclose(al[i] @ D @ al[j], 2 * (i == j) - 2 * (i == (j + 1) % h)) for i in range(h) for j in range(h))
    ev = np.linalg.eigvals(D)
    ok_j &= all(np.min(np.abs(ev - (1 + 1j / math.tan(math.pi * k / h)))) < 1e-10 for k in range(1, h))
    lam = [sum(mu[1:i + 1], np.zeros(rr)) for i in range(h)] + [np.zeros(rr)]   # lambda_0 = lambda_h = 0
    Dbook = 2 * sum(np.outer(lam[i] - lam[i + 1], lam[i]) for i in range(1, h))
    ok_j &= np.allclose(D, Dbook)
    sig = 0.63
    for _ in range(3):
        phi, psi = rng.normal(size=rr) * 0.7, rng.normal(size=rr) * 0.7
        pc, qc = comps(S, phi), comps(S, psi)
        f = [np.exp(qc[j] - pc[(j + 1) % h]) for j in range(h)]
        g = [np.exp(pc[j] - qc[j]) for j in range(h)]
        xi, eta = D.T @ phi + D @ psi, phi - psi
        ok_j &= all(np.isclose(f[j], np.exp(al[j] @ xi / 2)) and np.isclose(g[j], np.exp((D.T @ al[j]) @ eta / 2))
                    for j in range(h))
        ok_j &= all(np.isclose(f[j] * g[j], np.exp(al[j] @ phi)) and np.isclose(f[j] * g[(j + 1) % h],
                                                                              np.exp(al[j] @ psi)) for j in range(h))
        Gp, Gq, gf, gg = grads_ar(S, phi, psi, sig)
        # numerical gradient of B_ar against the analytic one
        num = np.array([(B_ar(S, phi + 1e-6 * e, psi, sig) - B_ar(S, phi - 1e-6 * e, psi, sig)) / 2e-6
                        for e in np.eye(rr)])
        ok_j &= np.allclose(num, Gp, atol=1e-7)
        V = lambda v: m**2 * sum(np.exp(a @ v) - 1 for a in al)
        ok_j &= np.isclose(0.5 * Gp @ Gp - 0.5 * Gq @ Gq, V(phi) - V(psi))
        ok_j &= not np.isclose(0.5 * (D @ (1.1 * gf) + gg) @ (D @ (1.1 * gf) + gg)
                               - 0.5 * (D.T @ (1.1 * gf) - gg) @ (D.T @ (1.1 * gf) - gg), V(phi) - V(psi))
        E = S['E']
        Pd = lambda a, b: m * sum(np.exp(comps(S, a)[j] - comps(S, b)[j]) - 1 for j in range(h)) / sig - m * sig * \
            sum(np.exp(comps(S, b)[j] - comps(S, a)[(j + 1) % h]) - 1 for j in range(h))
        numPp = np.array([(Pd(phi + 1e-6 * e, psi) - Pd(phi - 1e-6 * e, psi)) / 2e-6 for e in np.eye(rr)])
        numPq = np.array([(Pd(phi, psi + 1e-6 * e) - Pd(phi, psi - 1e-6 * e)) / 2e-6 for e in np.eye(rr)])
        ok_j &= np.allclose(numPp, E @ Gp - D @ Gq, atol=1e-6) and np.allclose(numPq, -D.T @ Gp - E @ Gq, atol=1e-6)
        # light-cone split: B + P_d is the g half (times 2), B - P_d the f half
        ok_j &= np.isclose(B_ar(S, phi, psi, sig) + Pd(phi, psi), 2 * m / sig * sum(gj - 1 for gj in g))
        # matrix conditions <-> component light-cone equations
        pt, qt = rng.normal(size=rr), rng.normal(size=rr)
        px = E @ pt + D @ qt - Gp
        qx = D.T @ pt - E @ qt + Gq
        a_p = comps(S, (pt + px) - (qt + qx))
        ok_j &= all(np.isclose(a_p[k], -m * sig * (f[k] - f[(k - 1) % h])) for k in range(h))
        lc = comps(S, pt - px), comps(S, qt - qx)          # d_- phi, d_- psi in components
        ok_j &= all(np.isclose(lc[1][k] - lc[0][(k + 1) % h], m / sig * (g[k] - g[(k + 1) % h])) for k in range(h))
report("(j) a_r defect, h = 2..5: D alpha_j = -2 mu_{j+1}, D^T alpha_j = 2 mu_j, alpha_i.D alpha_j = 2 delta_ij - "
       "2 delta_{i,j+1}, eigenvalues 1 + i cot(pi k/h), the book's D; f_j = e^{psi_j - phi_{j+1}}, g_j = e^{phi_j - "
       "psi_j}, f_j g_j and f_j g_{j+1}; conditions 2 and 3 (and a wrong coefficient fails); the component form", ok_j)

# ---------------------------------------------------------------------------------------------------------
# (k) a_1^(1): Lecture 3's sinh-Gordon defect
S = ar_setup(2)
ok_k = abs(abs(S['U'][0, 0]) - 1 / math.sqrt(2)) < 1e-12
sig = 0.8
m0 = 2 * m
for _ in range(4):
    c, cp = rng.normal(), rng.normal()     # orthonormal coordinates: phi_1 = c U[0,0]
    phi, psi = np.array([c]), np.array([cp])
    Phi, Psi = 2 * comps(S, phi)[0], 2 * comps(S, psi)[0]
    BL3 = 2 * m0 * (sig * (math.cosh((Phi + Psi) / 2) - 1) + (math.cosh((Phi - Psi) / 2) - 1) / sig)
    ok_k &= math.isclose(2 * B_ar(S, phi, psi, sig), BL3, rel_tol=1e-12)
    # conditions: c_x = c'_t - dB/dc  <->  Phi_x = Psi_t - dB_L3/dPhi  with Phi = 2 U00 c
    Gp, Gq, _, _ = grads_ar(S, phi, psi, sig)
    fac = 2 * S['U'][0, 0]
    dBL3_dPhi = m0 * (sig * math.sinh((Phi + Psi) / 2) + math.sinh((Phi - Psi) / 2) / sig)
    dBL3_dPsi = m0 * (sig * math.sinh((Phi + Psi) / 2) - math.sinh((Phi - Psi) / 2) / sig)
    ok_k &= math.isclose(fac * Gp[0], dBL3_dPhi, rel_tol=1e-10) and math.isclose(fac * Gq[0], dBL3_dPsi, rel_tol=1e-10)
    ok_k &= abs(S['D'][0, 0] - 1) < 1e-12 and abs(S['E'][0, 0]) < 1e-12
# the bulk: Phi_tt - Phi_xx = -m0^2 sinh Phi, from box c = -dV/dc
cc = 0.37
Phi = 2 * comps(S, np.array([cc]))[0]
dV = m**2 * sum(a[0] * math.exp(a[0] * cc) for a in S['al'])
ok_k &= math.isclose(-2 * S['U'][0, 0] * dV, -m0**2 * math.sinh(Phi), rel_tol=1e-12)
# continuation phi -> i phi, beta^2 -> -beta^2 gives Lecture 3's sine-Gordon B = -2 m0 [sigma(cos-1) + (cos-1)/sigma]
Pr, Qr = 0.9, -0.4
ok_k &= cmath.isclose(-(2 * m0 * (sig * (cmath.cosh(1j * (Pr + Qr) / 2) - 1) + (cmath.cosh(1j * (Pr - Qr) / 2) - 1) / sig)),
                      -2 * m0 * (sig * (math.cos((Pr + Qr) / 2) - 1) + (math.cos((Pr - Qr) / 2) - 1) / sig))
report("(k) for a_1^(1) the defect is Lecture 3's sinh-Gordon defect with Phi = 2 phi_1, m0 = 2m, B_L3 = 2B, D = 1, "
       "E = 0, and its continuation is the sine-Gordon defect", ok_k)

# ---------------------------------------------------------------------------------------------------------
# (l) unfreezing: second derivatives from the light-cone system
ok_l = True
for h in (3, 4):
    S = ar_setup(h)
    D, al, mu, rr = S['D'], S['al'], S['mu'], S['r']
    sig = 0.7
    for _ in range(3):
        phi, psi = rng.normal(size=rr) * 0.6, rng.normal(size=rr) * 0.6
        xi, eta = D.T @ phi + D @ psi, phi - psi
        gf = m * sig * sum(al[j] / 2 * np.exp(al[j] @ xi / 2) for j in range(h))
        Hf = m * sig * sum(np.outer(al[j], al[j]) / 4 * np.exp(al[j] @ xi / 2) for j in range(h))
        gg = m / sig * sum(mu[j] * np.exp(mu[j] @ eta) for j in range(h))
        Hg = m / sig * sum(np.outer(mu[j], mu[j]) * np.exp(mu[j] @ eta) for j in range(h))
        box_eta = -2 * Hf @ (2 * gg)          # d_- d_+ eta = -2 Hess f d_- xi
        box_xi = 2 * Hg @ (-2 * gf)           # d_+ d_- xi = 2 Hess g d_+ eta
        box_phi = (box_xi + D @ box_eta) / 2
        box_psi = (box_xi - D.T @ box_eta) / 2
        gradV = lambda v: m**2 * sum(a * np.exp(a @ v) for a in al)
        ok_l &= np.allclose(box_phi, -gradV(phi)) and np.allclose(box_psi, -gradV(psi))
        bc = comps(S, box_phi)
        pc = comps(S, phi)
        ok_l &= all(np.isclose(bc[k], -m**2 * (np.exp(pc[k] - pc[(k + 1) % h]) - np.exp(pc[(k - 1) % h] - pc[k])))
                    for k in range(h))


def soliton_fields(h, a, th, xi0, fac, x, t):
    """phi_k = ln tau_{k-1} - ln tau_k with tau_k = 1 + omega^{ka} fac E, and d_+ phi_k, d_- phi_k (analytic).
    The principal logarithms telescope, so sum_k phi_k = 0 exactly and phi lies in the hyperplane."""
    om = cmath.exp(2j * math.pi / h)
    ma = 2 * m * math.sin(math.pi * a / h)
    E = fac * cmath.exp(ma * (x * cmath.cosh(th) - t * cmath.sinh(th)) + xi0)
    tau = [1 + om**(k * a) * E for k in range(h)]
    dlog_p = [om**(k * a) * E * ma * cmath.exp(-th) / tau[k] for k in range(h)]       # d_+ ln tau_k
    dlog_m = [-om**(k * a) * E * ma * cmath.exp(th) / tau[k] for k in range(h)]       # d_- ln tau_k
    f = np.array([cmath.log(tau[(k - 1) % h]) - cmath.log(tau[k]) for k in range(h)])
    fp = np.array([dlog_p[(k - 1) % h] - dlog_p[k] for k in range(h)])
    fm = np.array([dlog_m[(k - 1) % h] - dlog_m[k] for k in range(h)])
    return f, fp, fm


def lightcone_residual(h, ph, php, phm, ps, psp, psm, sig):
    f = [cmath.exp(ps[k] - ph[(k + 1) % h]) for k in range(h)]
    g = [cmath.exp(ph[k] - ps[k]) for k in range(h)]
    r1 = max(abs(php[k] - psp[k] + m * sig * (f[k] - f[(k - 1) % h])) for k in range(h))
    r2 = max(abs(psm[k] - phm[(k + 1) % h] - m / sig * (g[k] - g[(k + 1) % h])) for k in range(h))
    return max(r1, r2)


def z_a(u, a, h):
    ga = math.pi * a / h
    return (cmath.exp(-u) - 1j * cmath.exp(1j * ga)) / (cmath.exp(-u) - 1j * cmath.exp(-1j * ga))


res = 0.0
h, a, th, thd = 3, 1, 0.5, 0.2
sig = math.exp(-thd)
for x in (-0.7, 0.0, 1.3):
    for t in (-1.0, 0.4, 2.0):
        ph = soliton_fields(h, a, th, 0.1 + 0.3j, 1.0, x, t)
        ps = soliton_fields(h, a, th, 0.1 + 0.3j, z_a(th - thd, a, h), x, t)
        res = max(res, lightcone_residual(h, *ph, *ps, sig))
ok_l &= res < 1e-12
report("(l) unfreezing: the light-cone system gives box phi = -grad V(phi), box psi = -grad V(psi), for a_r the "
       "periodic 2D Toda lattice; an a_2 soliton pair solves the light-cone system at x = -0.7, 0, 1.3 (Baecklund pair)",
       ok_l)

# ---------------------------------------------------------------------------------------------------------
# (m) type I_+ and the reflection rho
ok_m = True
for h in (2, 3, 4, 5):
    S = ar_setup(h)
    D, E, al, pi, rr = S['D'], S['E'], S['al'], S['pi'], S['r']
    # condition 1 for (-D, -E, -E) with O = pi^{-1}: -D = (1 + E) O
    O = np.linalg.inv(pi)
    ok_m &= np.allclose(-D, (np.eye(rr) + E) @ O) and np.allclose(O.T @ (-E) @ O, -E)
    A = np.array(al[1:]).T
    rho = np.array([al[(-j) % h] for j in range(1, h)]).T @ np.linalg.inv(A)
    ok_m &= np.allclose(rho.T @ rho, np.eye(rr)) and np.allclose(rho @ pi @ rho, np.linalg.inv(pi))
    ok_m &= all(np.allclose(rho @ al[j], al[(-j) % h]) for j in range(h))
    A1, B1 = rho, rho @ np.linalg.inv(pi)
    ok_m &= np.allclose(A1.T @ E @ A1, -E) and np.allclose(A1.T @ D @ B1, -D) and np.allclose(B1.T @ E @ B1, -E)
    # the mirror image: D' = 2 (1 - pi^{-1})^{-1} = D^T, E' = 1 - D' = -E; L_d(rho phi, rho psi) is the defect built
    # from D' by the general recipe, and D' passes the cross-term condition for p(j) = j + 1
    Dm = 2 * np.linalg.inv(np.eye(rr) - np.linalg.inv(pi))
    ok_m &= np.allclose(Dm, D.T) and np.allclose(rho @ D @ rho, Dm) and np.allclose(rho @ E @ rho, np.eye(rr) - Dm)
    ok_m &= all(abs(al[i] @ Dm @ al[j]) < 1e-12 for i in range(h) for j in range(h) if j not in (i, (i + 1) % h))
    for _ in range(3):
        phi, psi = rng.normal(size=rr), rng.normal(size=rr)
        xi_m, eta_m = Dm.T @ phi + Dm @ psi, phi - psi
        Bm = m * sum(0.55 * (np.exp(a_ @ xi_m / 2) - 1) + (np.exp((Dm.T @ a_) @ eta_m / 2) - 1) / 0.55 for a_ in al)
        ok_m &= np.isclose(B_ar(S, rho @ phi, rho @ psi, 0.55, '-'), Bm)
    for _ in range(3):
        phi, psi = rng.normal(size=rr), rng.normal(size=rr)
        ok_m &= np.isclose(B_ar(S, A1 @ phi, B1 @ psi, 0.55, '-'), B_ar(S, phi, psi, 0.55, '+'))
    if h == 2:
        ok_m &= np.allclose(rho, np.eye(1)) and np.allclose(B1, -np.eye(1))


# a type-I_- soliton pair of species h - a, mapped by (rho, pi rho), solves the type-I_+ conditions
def matrix_residual(S, ph, php, phm, ps, psp, psm, sig, typ):
    """Residual of the matrix defect conditions for complex fields given in components."""
    U = S['U']
    D, E = S['D'], S['E']
    if typ == '+':
        D, E = -D, -E
    p, q = U.T @ ph, U.T @ ps
    pt, px = U.T @ (php + phm) / 2, U.T @ (php - phm) / 2
    qt, qx = U.T @ (psp + psm) / 2, U.T @ (psp - psm) / 2
    h, al, mu, Dm = S['h'], S['al'], S['mu'], S['D']
    xi, eta = Dm.T @ p + Dm @ q, p - q
    cf, cg = (m * sig, m / sig) if typ == '-' else (m / sig, m * sig)
    gf = cf * sum(al[j] / 2 * np.exp(al[j] @ xi / 2) for j in range(h))
    gg = cg * sum(mu[j] * np.exp(mu[j] @ eta) for j in range(h))
    Gp, Gq = Dm @ gf + gg, Dm.T @ gf - gg
    return max(np.max(np.abs(px - (E @ pt + D @ qt - Gp))), np.max(np.abs(qx - (D.T @ pt - E @ qt + Gq))))


res_minus, res_plus, res_ctrl = 0.0, 0.0, 0.0
for h in (3, 4):
    S = ar_setup(h)
    U = S['U']
    rho = np.array([S['al'][(-j) % h] for j in range(1, h)]).T @ np.linalg.inv(np.array(S['al'][1:]).T)
    for a in range(1, h):
        th, thd = 0.6, 0.15
        sig = math.exp(-thd)
        zz = z_a(th - thd, h - a, h)
        for t in (-0.8, 0.3, 1.5):
            Ph = soliton_fields(h, h - a, th, 0.2 + 0.1j, 1.0, 0.0, t)
            Ps = soliton_fields(h, h - a, th, 0.2 + 0.1j, zz, 0.0, t)
            res_minus = max(res_minus, matrix_residual(S, *Ph, *Ps, sig, '-'))
            mapped_ph = [U @ (rho @ (U.T @ v)) for v in Ph]
            mapped_ps = [U @ (S['pi'] @ rho @ (U.T @ v)) for v in Ps]
            res_plus = max(res_plus, matrix_residual(S, *mapped_ph, *mapped_ps, sig, '+'))
            res_ctrl = max(res_ctrl, matrix_residual(S, *mapped_ph, *mapped_ps, sig, '-'))
            # the mapped pair is the species-a soliton with E -> z_a(theta_d - theta) E (the book's type-I_+ delay)
            Qh = soliton_fields(h, a, th, 0.2 + 0.1j, 1.0, 0.0, t)
            Qs = soliton_fields(h, a, th, 0.2 + 0.1j, z_a(thd - th, a, h), 0.0, t)
            res_plus = max(res_plus, matrix_residual(S, *Qh, *Qs, sig, '+'))
ok_m &= res_minus < 1e-12 and res_plus < 1e-12 and res_ctrl > 1e-3
ok_m &= all(cmath.isclose(z_a(-u, a, h), cmath.exp(2j * math.pi * a / h) * z_a(u, h - a, h))
            for h in (3, 4, 5) for a in range(1, h) for u in (-1.1, 0.3, 0.7 + 0.2j))
report("(m) type I_+ is condition 1 with O = pi^{-1}; L_+(phi, psi; sigma) = L_-(rho phi, rho pi^{-1} psi; sigma), "
       "kinetic and potential terms; a type-I_- soliton pair mapped by (rho, pi rho) solves the type-I_+ conditions "
       "(and not the type-I_- ones); z_a(-u) = omega^a z_{h-a}(u)", ok_m)

# ---------------------------------------------------------------------------------------------------------
# (n) linear waves
ths, sgs, gms = sp.symbols('theta sigma gamma', positive=True)
Rk, Tk = sp.symbols('R T')
mk = 2 * sp.sin(gms)                                        # m = 1
d = 1 + sp.I * sp.cot(gms)
db = 1 - sp.I * sp.cot(gms)
a_in, a_out = mk * sp.exp(-ths), mk * sp.exp(ths)           # d_+ -> -i a on the incident wave, -i a_out on R
# incident e^{i(px - nu t)}: d_+ -> -i m_k e^{-theta}, d_- -> -i m_k e^{theta}; reflected: exchanged
eq1 = sp.Eq(-sp.I * (a_in * 1 + a_out * Rk - a_in * Tk), -2 * sgs * sp.sin(gms)**2 * (db * (1 + Rk) + d * Tk))
eq2 = sp.Eq(-sp.I * (a_out * db * 1 + a_in * db * Rk + a_out * d * Tk), 2 / sgs * (1 + Rk - Tk))
sol = sp.solve([eq1, eq2], [Rk, Tk], dict=True)[0]
Tclaim = (sp.exp(-ths) + sgs * sp.exp(sp.I * gms)) / (sp.exp(-ths) + sgs * sp.exp(-sp.I * gms))
ok_n = sp.simplify(sol[Rk]) == 0
Tf = sp.lambdify((ths, sgs, gms), sol[Tk], 'numpy')
Tc = sp.lambdify((ths, sgs, gms), Tclaim, 'numpy')
ok_n &= all(abs(Tf(a, b, c) - Tc(a, b, c)) < 1e-12 for a, b, c in rng.uniform(0.1, 1.5, size=(6, 3)))


def T_k(theta, thd, k, h):
    g = math.pi * k / h
    s = math.exp(-thd)
    return (cmath.exp(-theta) + s * cmath.exp(1j * g)) / (cmath.exp(-theta) + s * cmath.exp(-1j * g))


# the full vector linearized conditions, with Hessians of B by finite differences
for h in (3, 4):
    S = ar_setup(h)
    rr, pi = S['r'], S['pi']
    thd = 0.3
    sig = math.exp(-thd)
    eps_h = 1e-4

    def grad_num(phi, psi):
        out = []
        for e in np.eye(2 * rr):
            zp = np.concatenate([phi, psi]) + eps_h * e
            zm = np.concatenate([phi, psi]) - eps_h * e
            out.append((B_ar(S, zp[:rr], zp[rr:], sig) - B_ar(S, zm[:rr], zm[rr:], sig)) / (2 * eps_h))
        return np.array(out)

    Hs = np.array([(grad_num(*np.split(eps_h * e, 2)) - grad_num(*np.split(-eps_h * e, 2))) / (2 * eps_h)
                   for e in np.eye(2 * rr)])
    D, E = S['D'], S['E']
    w, V = np.linalg.eig(pi)
    for k in range(1, h):
        rho_k = V[:, np.argmin(np.abs(w - cmath.exp(2j * math.pi * k / h)))]
        mk_ = 2 * m * math.sin(math.pi * k / h)
        theta = 0.8
        nu, p = mk_ * math.cosh(theta), mk_ * math.sinh(theta)
        # unknowns R, T: phi = rho (e^{ipx} + R e^{-ipx}) e^{-i nu t}, psi = T rho e^{ipx} e^{-i nu t}, at x = 0
        cols, rhs = [], []
        # phi_x - E phi_t - D psi_t + H_phi z = 0, psi_x - D^T phi_t + E psi_t - H_psi z = 0
        def lin(phi_val, phi_t, phi_x, psi_val, psi_t, psi_x):
            zv = np.concatenate([phi_val, psi_val])
            Hz = Hs @ zv
            e1 = phi_x - E @ phi_t - D @ psi_t + Hz[:rr]
            e2 = psi_x - D.T @ phi_t + E @ psi_t - Hz[rr:]
            return np.concatenate([e1, e2])
        inc = lin(rho_k, -1j * nu * rho_k, 1j * p * rho_k, 0 * rho_k, 0 * rho_k, 0 * rho_k)
        cR = lin(rho_k, -1j * nu * rho_k, -1j * p * rho_k, 0 * rho_k, 0 * rho_k, 0 * rho_k)
        cT = lin(0 * rho_k, 0 * rho_k, 0 * rho_k, rho_k, -1j * nu * rho_k, 1j * p * rho_k)
        Amat = np.column_stack([cR, cT])
        solv, *_ = np.linalg.lstsq(Amat, -inc, rcond=None)
        resid = np.max(np.abs(Amat @ solv + inc))
        ok_n &= resid < 1e-6 and abs(solv[0]) < 1e-6 and abs(solv[1] - T_k(theta, thd, k, h)) < 1e-6
# properties
for h in (2, 3, 4, 5):
    for k in range(1, h):
        g = math.pi * k / h
        ok_n &= all(abs(abs(T_k(th, 0.4, k, h)) - 1) < 1e-13 for th in (-2.0, 0.1, 1.7))
        ok_n &= cmath.isclose(T_k(0.4, 0.4, k, h), cmath.exp(1j * g))
        ok_n &= abs(T_k(-30, 0.4, k, h) - 1) < 1e-10 and abs(T_k(30, 0.4, k, h) - cmath.exp(2j * g)) < 1e-10
ok_n &= all(cmath.isclose(T_k(th, 0.4, 1, 2), (cmath.exp(-th) + 1j * math.exp(-0.4)) /
                          (cmath.exp(-th) - 1j * math.exp(-0.4))) for th in (-1.0, 0.3, 2.2))
# large sigma: d_+ eta = -m sigma sum alpha_j f_j forces sum_j alpha_j f_j = 0, hence all f_j equal
for h in (3, 4, 5):
    S = ar_setup(h)
    Amat = np.array(S['al']).T            # r x h
    ns = np.linalg.svd(Amat)[2][-1]
    ok_n &= np.allclose(np.abs(ns), np.ones(h) / math.sqrt(h)) and np.linalg.matrix_rank(Amat) == h - 1
    # f_j = 1 for all j means psi_j = phi_{j+1}, i.e. psi = pi phi
    phi = rng.normal(size=h - 1)
    psi = S['pi'] @ phi
    pc, qc = comps(S, phi), comps(S, psi)
    ok_n &= all(abs(qc[j] - pc[(j + 1) % h]) < 1e-12 for j in range(h))
report("(n) linear waves: R_k = 0 and T_k = (e^{-theta} + sigma e^{i gamma_k})/(e^{-theta} + sigma e^{-i gamma_k}) "
       "(sympy, and the full vector conditions for a_3, a_4 with numerical Hessians); |T_k| = 1, T_k(theta_d) = "
       "e^{i gamma_k}, limits 1 and omega^k; h = 2 gives Lecture 4's T; large sigma gives psi = pi phi", ok_n)

# ---------------------------------------------------------------------------------------------------------
# (o) solitons through the a_r defect
ok_o = True
res, res_bad = 0.0, 1.0
for h in (2, 3, 4, 5):
    for a in range(1, h):
        for th, thd in ((0.5, 0.2), (1.4, 0.3), (-0.3 + 0.2j, 0.4)):
            sig = math.exp(-thd)
            zz = z_a(th - thd, a, h)
            worst_bad = 0.0
            for x in (-0.7, 0.0, 1.3):
                for t in (-1.2, 0.0, 0.9):
                    ph = soliton_fields(h, a, th, 0.3 - 0.2j, 1.0, x, t)
                    ps = soliton_fields(h, a, th, 0.3 - 0.2j, zz, x, t)
                    res = max(res, lightcone_residual(h, *ph, *ps, sig))
                    if h > 2 and 2 * a != h:
                        ps_b = soliton_fields(h, a, th, 0.3 - 0.2j, zz.conjugate(), x, t)
                        worst_bad = max(worst_bad, lightcone_residual(h, *ph, *ps_b, sig))
            if h > 2 and 2 * a != h:
                res_bad = min(res_bad, worst_bad)
ok_o &= res < 1e-11 and res_bad > 1e-3
for h in (2, 3, 4, 5):
    om = cmath.exp(2j * math.pi / h)
    for a in range(1, h):
        ga = math.pi * a / h
        thd = 0.35
        ok_o &= all(cmath.isclose(z_a(u, a, h), T_k(u + thd - 0.5j * math.pi, thd, a, h)) for u in (-0.8, 0.2, 1.1))
        ok_o &= all(abs(z_a(u, a, h)) > 1 for u in np.linspace(-4, 4, 41) if abs(u) > 1e-9)
        ok_o &= all(abs(cmath.phase(z_a(u, a, h)) - cmath.phase(cmath.exp(-2 * u) + cmath.exp(2j * ga))) < 1e-12
                    for u in (-1.0, 0.3, 0.9) if 2 * a != h)
        ok_o &= abs(z_a(-30, a, h) - 1) < 1e-10 and abs(z_a(30, a, h) - om**a) < 1e-10
        pole = -1j * (math.pi / 2 - ga)
        ok_o &= abs(1 / z_a(pole + 1e-9, a, h)) < 1e-6
        ok_o &= (abs(pole.imag) < 1e-12) == (2 * a == h)
    # E -> omega^a E turns tau_k into tau_{k+1}, so psi_k = phi_{k+1}: psi = pi phi
    S = ar_setup(h)
    for a in range(1, h):
        ph = soliton_fields(h, a, 0.4, 0.1, 1.0, 0.3, 0.2)[0]
        ps = soliton_fields(h, a, 0.4, 0.1, om**a, 0.3, 0.2)[0]
        pi_ph = comps(S, S['pi'] @ (S['U'].T @ ph))
        ok_o &= np.allclose(ps, pi_ph)
report("(o) solitons: E -> z_a E solves the component light-cone equations (h = 2..5, all species, real and complex "
       "theta, x = -0.7, 0, 1.3), the conjugate fails; z_a(u) = T_a(u + theta_d - i pi/2); |z_a| > 1; arg z_a = "
       "arg(e^{-2u} + e^{2 i gamma_a}); limits 1 and omega^a = psi = pi phi; the pole is real only for a = h/2", ok_o)

# ---------------------------------------------------------------------------------------------------------
# (p) exercises
ok_p = True
S = ar_setup(3)
E, pi = S['E'], S['pi']
Jm = np.array([[0, 1], [-1, 0]])
s3 = 1 / math.sqrt(3)
ok_p &= np.allclose(E, -s3 * Jm) or np.allclose(E, s3 * Jm)
ang = math.atan2(pi[1, 0], pi[0, 0])
ok_p &= np.isclose(abs(ang), 2 * math.pi / 3) and np.isclose(np.linalg.det(pi), 1)
# in the basis used, E = eps J gives P = -(1 - E)^{-1}(1 + E) = the rotation by pi + 2 arctan(eps)... check directly
P = -np.linalg.inv(np.eye(2) - E) @ (np.eye(2) + E)
ok_p &= np.allclose(P, pi)
# a_3 reflection (0 1)(2 3): every node goes to a neighbour, but node 0 has the neighbour 3 outside {p(0), p^{-1}(0)}
G = DIAGRAMS['a3^(1)']
p = (1, 0, 3, 2)
ok_p &= all(abs(G[p[i], p[j]] - G[i, j]) < 1e-12 for i in range(4) for j in range(4))
ok_p &= all(abs(G[i, p[i]]) > 0 for i in range(4)) and abs(G[0, 3]) > 0 and 3 not in (p[0], p.index(0))
# a wave from the right: Psi = e^{-i(px + nu t)} + barR e^{i(px - nu t)}, Phi = barT e^{-i(px + nu t)}
bR, bT = sp.symbols('bR bT')
# at x = 0: incident from the right (moving left): d_+ -> -i(nu - p) = -i a_out? for e^{-i(px+nu t)}: d_t -> -i nu,
# d_x -> -i p, so d_+ -> -i (nu + p) = -i m e^{theta}, d_- -> -i (nu - p) = -i m e^{-theta}
c_l_p, c_l_m = -sp.I * a_out, -sp.I * a_in        # left-moving waves
c_r_p, c_r_m = -sp.I * a_in, -sp.I * a_out        # right-moving waves
eq1 = sp.Eq(c_l_p * bT - (c_l_p * 1 + c_r_p * bR), -2 * sgs * sp.sin(gms)**2 * (db * bT + d * (1 + bR)))
eq2 = sp.Eq(c_l_m * db * bT + d * (c_l_m * 1 + c_r_m * bR), 2 / sgs * (bT - 1 - bR))
solR = sp.solve([eq1, eq2], [bR, bT], dict=True)[0]
ok_p &= sp.simplify(solR[bR]) == 0
bTf = sp.lambdify((ths, sgs, gms), solR[bT], 'numpy')
ok_p &= all(abs(bTf(a, b, c) * Tc(-a, b, c) - 1) < 1e-12 for a, b, c in rng.uniform(0.1, 1.4, size=(6, 3)))
# the book's soliton phi_c = -sum_j alpha_j ln tau_j has components ln tau_{k-1} - ln tau_k; Hirota form
for h in (3, 4):
    S = ar_setup(h)
    for a in range(1, h):
        om = cmath.exp(2j * math.pi / h)
        ma = 2 * m * math.sin(math.pi * a / h)
        th, x, t, c0 = 0.45, 0.3, -0.2, 0.1 + 0.4j
        Y = cmath.exp(ma * (x * math.cosh(th) - t * math.sinh(th)) + c0)
        tau = [1 + om**(k * a) * Y for k in range(h)]
        phic = -sum(S['al'][j] * cmath.log(tau[j]) for j in range(h))
        ok_p &= np.allclose(comps(S, phic), [cmath.log(tau[(k - 1) % h]) - cmath.log(tau[k]) for k in range(h)])
        # tau d+d- tau - d+tau d-tau = m^2 (tau_{k-1} tau_{k+1} - tau_k^2), with d+ Y = ma e^{-th} Y, d- Y = -ma e^{th} Y
        for k in range(h):
            w = om**(k * a) * Y
            lhs = tau[k] * (-ma**2 * w) - (ma * math.exp(-th) * w) * (-ma * math.exp(th) * w)
            rhs = m**2 * (tau[(k - 1) % h] * tau[(k + 1) % h] - tau[k]**2)
            ok_p &= abs(lhs - rhs) < 1e-12
# zero and pole of T_k: theta_d + i (pi - gamma_k) and theta_d - i (pi - gamma_k)
for h in (3, 4, 5):
    for k in range(1, h):
        g = math.pi * k / h
        ok_p &= abs(T_k(0.4 + 1j * (math.pi - g), 0.4, k, h)) < 1e-12
        ok_p &= abs(1 / T_k(0.4 - 1j * (math.pi - g) + 1e-9, 0.4, k, h)) < 1e-6
report("(p) exercises: for a_2, E = +-J/sqrt 3 and P = pi is the rotation by 2 pi/3; the a_3 reflection (0 1)(2 3) "
       "sends each node to a neighbour but violates the neighbour lemma; a wave from the right has bar R = 0 and "
       "bar T_k(theta) = 1/T_k(-theta); the book's soliton in components and its Hirota form; the zero and the pole "
       "of T_k at theta_d +- i (pi - gamma_k)", ok_p)

raise SystemExit(0 if OK else 1)
