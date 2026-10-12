"""Classical a_r^(1) affine Toda theory with a jump defect, book normalization (Part VIII, ch. 36, 38).

Conventions: defect_common.py (fields phi on x < 0, psi on x > 0, a_t of eq-boundary-lax-at, longest roots alpha^2 = 2).
BCZ = Bowcock, Corrigan, Zambon hep-th/0401020.

Checks (book ids they support):
 1. D identities, r = 1..5: D + D^T = 2, E = 1 - D antisymmetric, D alpha_j = 2(lambda_j - lambda_{j+1}), D alpha_0 = -2 lambda_1,
    D = 2(1 - pi)^{-1} with pi alpha_j = alpha_{j-1} orthogonal, D^T alpha_j = -D alpha_{j-1}, D normal,
    alpha_i.D alpha_j = 2 delta_ij - 2 delta_{i,j+1}, E = sum_i (lambda_{i+1} lambda_i^T - lambda_i lambda_{i+1}^T) (BCZ (3.30)).   [eq-d-matrix]
 2. Defect conditions from the action by variation (sympy, generic D and antisymmetric E, r = 2).                 [eq-defect-conditions]
 3. Integrability, r = 1..4, real and imaginary coupling, both types, random complex points:
    1/2|grad_phi B_d|^2 - 1/2|grad_psi B_d|^2 = V(phi) - V(psi) with the prefactor m/beta^2; along the defect conditions
    dE/dt = -dB_d/dt and dP/dt = -dP_d/dt with P = -int phi_t.phi_x (energy and momentum conservation).              [eq-defect-potential]
 4. Lax: the space of field-independent Ltilde with d_t L_d = L_d a_t(psi) - a_t(phi) L_d,
    L_d = exp(-b H.(E phi + D psi)/2) Ltilde exp(b H.(D^T phi - E psi)/2), on the defect conditions, defining rep,
    a_1, a_2, a_3, both couplings: one-dimensional; type I_- Ltilde = 1 + sigma lambda^{-1} sum_j E_{-alpha_j} (c = 1, i.e.
    BCZ (3.39) with d = sigma and their lambda = the book's lambda); type I_+ (B_d^{I+}(sigma) = B_d^{I-}(1/sigma))
    Ltilde^{-1} ~ 1 - sigma^{-1} lambda sum_j E_{alpha_j} (Ltilde = sum_k (lambda/sigma)^k P^k); a non-integrable D has none.
    For a soliton of rapidity theta the vacuum Lax solutions give lambda = -e^{-theta}.                                [eq-defect-lax-matrix, prp-defect-lax]
 5. Spin +-2 charges of a_2^(1) (sympy + Euler operator in t, and closed-loop integrals): with Q = int (T - Theta) dx
    built from d_eps = d_t + eps d_x (T = A(u,u,u) + w.Bu, u = d_eps phi; bulk current checked on shell; density not a total
    x-derivative), the defect term eps[(T + Theta)_phi - (T + Theta)_psi]_{x=0} is a total time derivative of a local
    function of phi(0,t), psi(0,t) and their t-derivatives for BOTH eps = +1 (spin -2) and eps = -1 (spin +2), for both
    types.  This contradicts BCZ hep-th/0401020 sec. 5 (and Penati-Polvara 1902.10690 App. B): their obstruction (5.11),
    A(Dv)^3 + A(D^T v)^3, vanishes identically because cot a cot b + cot b cot c + cot c cot a = 1 for a + b + c = 0 mod pi
    (D = 1 + i cot(pi k/h) on the Fourier mode k).  Controls: spin +-1; a non-integrable D fails.  (The same holds for
    a_3, checked separately at about 6 minutes per case.)  Book spin convention: spin s = e^{s theta} on a soliton of
    rapidity theta (sec-qg-gradation), so E + P has spin +1 and spin +s densities are built from d_- = d_t - d_x.  [prp-spin-two]
Chapter 38 additions (sec-defect-conditions):
 6. prp-ar-only: every affine Dynkin diagram of rank <= 4 (a_1..a_4, b_3, b_4, c_2..c_4, d_4, f_4, g_2; a_2, a_4, a_6, a_8,
    a_3 = d_3, a_5, a_7, d_4, d_5, e_6 twisted, d_4^(3)) and every diagram automorphism p: the type-I exponential ansatz forces
    D = 2(1 - P)^{-1} (P alpha_j = alpha_p(j)); p must be fixed-point free with 1 - P invertible, and alpha_i.D alpha_j = 0
    for j not in {i, p(i)}; the momentum identity is then tested at random points.  Only a_r^(1) with p = pi^{+-1} passes
    (pi^{+-2} of a_4 fail by cross terms).  Also cz2009a (2.7), (2.9) in book form (grad P_d).        [prp-ar-only, prp-defect-momentum]
 7. BCZ sec. 5.1: the defect conditions imposed for all x are a Baecklund transformation (cross-differentiation gives both
    field equations, r = 1..4, both types and couplings; a non-integrable D fails); BCZ (5.43) as printed is the transpose
    of their (3.28).                                                                                    [sec-defect-backlund-an]
 8. a_r solitons through the type-I_- defect in book conventions: delay z_a = (e^-(th-th_d) - i e^(i pi a/h))/(e^-(th-th_d)
    - i e^(-i pi a/h)) (by formula and by least squares), the pair (soliton, delayed soliton) satisfies the conditions at every
    x (Baecklund), = cz2009a (4.10) = BCZ (6.2) at theta -> -theta = cz2007 (2.20) at a -> h - a; twisted Bianchi form
    z_a = (A + eps_a A_d)/(A - A_d/eps_a), A = -e^-theta, eps_a = e^(i pi(2a-h)/2h); z_a z_abar = (cosh u + sin(pi a/h))/
    (cosh u - sin(pi a/h)); Theta conjugates the delay; self-conjugate solitons: -coth (cz2007 (2.21) has the wrong sign).
                                                                                                         [sec-defect-backlund-an]
 9. Theta: conj B_d(phi, psi; sigma) = B_d(-conj phi, -conj psi; conj sigma), so Theta maps the sigma-defect to the
    conj(sigma)-defect and is a symmetry iff sigma is real; checked on a_2 soliton solutions.               [prp-defect-theta]
10. Type II (cz2009b units): sine-Gordon (3.2) and Tzitzeica (3.7) conserve E + D and P + Omega (Omega = f - g); solitons pass
    with the delays (4.3) [z1 = tanh tanh in the sector e^lam -> -1/(sqrt2 cosh tau), z2 = coth coth for +1/(sqrt2 cosh tau);
    corrected (4.4)] and (4.9) [(4.11) as printed is z3].                                                 [prp-type-ii-momentum]
11. prp-spin-two for all r: BCZ (5.5) with B = 1 - D is solved by A = (2/3) sum_k (mu_k.u)^3, and (5.11) vanishes for
    r = 2..8 (cot identity, h = 3..11); explicit B_{+2}, B_{-2} for a_2 type I_- (homotopy formula; linear in velocities plus
    exponentials), dB/dt = defect flux.                                                                   [prp-spin-two]
Runtime: about 55 seconds (checks 5 and 11 dominate).
"""
import itertools
import numpy as np
import sympy as sp
from defect_common import an_data, dmatrix, defect_potential, bulk_potential, defect_dx, a_t, an_lax_rep


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


rng = np.random.default_rng(11)

# ---------------------------------------------------------------- 1. D identities
ok = True
for r in range(1, 6):
    d = an_data(r); D, E = dmatrix(r); al, lam, pi = d['alpha'], d['lam'], d['pi']
    I = np.eye(r)
    ok &= np.allclose(D + D.T, 2*I) and np.allclose(E, -E.T) and np.allclose(E, I - D)
    ok &= np.allclose(al @ al.T, [[2 if i == j else (-1 if (i - j) % (r + 1) in (1, r) else 0) for j in range(r + 1)]
                                  for i in range(r + 1)]) if r > 1 else np.allclose(al @ al.T, [[2, -2], [-2, 2]])
    ok &= all(np.allclose(D @ al[j], 2*(lam[j] - lam[j + 1])) for j in range(1, r + 1))
    ok &= np.allclose(D @ al[0], -2*lam[1])
    ok &= np.allclose(pi @ pi.T, I) and all(np.allclose(pi @ al[j], al[(j - 1) % (r + 1)]) for j in range(r + 1))
    ok &= np.allclose(D, 2*np.linalg.inv(I - pi))
    ok &= all(np.allclose(D.T @ al[j], -D @ al[(j - 1) % (r + 1)]) for j in range(r + 1))
    ok &= np.allclose(D @ D.T, D.T @ D)
    G = np.array([[al[i] @ D @ al[j] for j in range(r + 1)] for i in range(r + 1)])
    want = np.array([[2*(i == j) - 2*(i == (j + 1) % (r + 1)) for j in range(r + 1)] for i in range(r + 1)])
    ok &= np.allclose(G, want) if r > 1 else np.allclose(G, [[2, -2], [-2, 2]])
    ok &= np.allclose(E, sum(np.outer(lam[i + 1], lam[i]) - np.outer(lam[i], lam[i + 1]) for i in range(1, r + 1)))
    ok &= np.allclose(D, 2*sum(np.outer(lam[i] - lam[i + 1], lam[i]) for i in range(1, r + 1)))
report('D = 2 sum (lambda_i - lambda_{i+1}) lambda_i^T = 2(1 - pi)^-1 (pi alpha_j = alpha_{j-1}): D + D^T = 2, '
       'D^T alpha_j = -D alpha_{j-1}, alpha_i.D alpha_j table, BCZ (3.30), r = 1..5', ok)
D1, E1 = dmatrix(1)
report('a_1: D = 1, E = 0, alpha_1 = sqrt2 = -alpha_0', np.allclose(D1, 1) and np.allclose(E1, 0) and np.allclose(an_data(1)['alpha'][:, 0], [-np.sqrt(2), np.sqrt(2)]))

# ---------------------------------------------------------------- 2. variation (generic D, antisymmetric E)
t = sp.symbols('t')
r = 2
ph = sp.Matrix([sp.Function(f'phi{k}')(t) for k in range(r)])
ps = sp.Matrix([sp.Function(f'psi{k}')(t) for k in range(r)])
Dg = sp.Matrix(r, r, sp.symbols('d0:4'))
eps_ = sp.symbols('e12')
Eg = sp.Matrix([[0, eps_], [-eps_, 0]])
Bf = sp.Function('B')(*ph, *ps)
Ld = (ph.T*Eg*ph.diff(t)/2 + ph.T*Dg*ps.diff(t) + ps.T*Eg*ps.diff(t)/2)[0] - Bf
phx = sp.Matrix(sp.symbols('phix0:2')); psx = sp.Matrix(sp.symbols('psix0:2'))
ok = True
for k in range(r):
    el_phi = sp.diff(Ld, ph[k]) - sp.diff(sp.diff(Ld, ph[k].diff(t)), t) - phx[k]
    el_psi = sp.diff(Ld, ps[k]) - sp.diff(sp.diff(Ld, ps[k].diff(t)), t) + psx[k]
    want_phi = (Eg*ph.diff(t) + Dg*ps.diff(t))[k] - sp.diff(Bf, ph[k])
    want_psi = (Dg.T*ph.diff(t) - Eg*ps.diff(t))[k] + sp.diff(Bf, ps[k])
    ok &= sp.expand(sp.solve(el_phi, phx[k])[0] - want_phi) == 0
    ok &= sp.expand(sp.solve(el_psi, psx[k])[0] - want_psi) == 0
report('variation of S: phi_x = E phi_t + D psi_t - grad_phi B, psi_x = D^T phi_t - E psi_t + grad_psi B (generic D, E^T = -E)', ok)

# ---------------------------------------------------------------- 3. integrability, energy and momentum
def grad_Pd(phi, psi, sigma, r, beta, m, coupling, typ):
    h = 1e-6
    g1 = np.array([(defect_potential(phi + h*e, psi, sigma, r, beta, m, coupling, typ)[3]
                    - defect_potential(phi - h*e, psi, sigma, r, beta, m, coupling, typ)[3])/(2*h) for e in np.eye(r)])
    g2 = np.array([(defect_potential(phi, psi + h*e, sigma, r, beta, m, coupling, typ)[3]
                    - defect_potential(phi, psi - h*e, sigma, r, beta, m, coupling, typ)[3])/(2*h) for e in np.eye(r)])
    return g1, g2


ok_int = True; ok_em = True; worst = 0
for r in range(1, 5):
    for coupling in ('real', 'imag'):
        for typ in (1, -1):
            for _ in range(3):
                phi = rng.normal(size=r) + 0.2j*rng.normal(size=r); psi = rng.normal(size=r) + 0.2j*rng.normal(size=r)
                phit = rng.normal(size=r) + 0.2j*rng.normal(size=r); psit = rng.normal(size=r) + 0.2j*rng.normal(size=r)
                sig, be, m = 0.6 + 0.3*rng.random(), 0.7, 1.3
                B, gp, gq, Pd = defect_potential(phi, psi, sig, r, be, m, coupling, typ)
                Vp, gVp = bulk_potential(phi, r, be, m, coupling); Vq, _ = bulk_potential(psi, r, be, m, coupling)
                res = abs((gp @ gp - gq @ gq)/2 - (Vp - Vq))/(abs(Vp) + abs(Vq) + 1)
                ok_int &= res < 1e-12
                phix, psix = defect_dx(phit, psit, gp, gq, r, typ)
                # dE/dt from the bulk = phi_t.phi_x - psi_t.psi_x at x = 0;  dB/dt by the chain rule
                dE = phit @ phix - psit @ psix + gp @ phit + gq @ psit
                # dP/dt from the bulk (P = -int phi_t.phi_x) at x = 0
                dP = -(phix @ phix)/2 - (phit @ phit)/2 + Vp + (psix @ psix)/2 + (psit @ psit)/2 - Vq
                g1, g2 = grad_Pd(phi, psi, sig, r, be, m, coupling, typ)
                dP += g1 @ phit + g2 @ psit
                worst = max(worst, abs(dE), abs(dP)/(1 + abs(Vp)))
ok_em = worst < 1e-7
report('1/2|grad_phi B_d|^2 - 1/2|grad_psi B_d|^2 = V(phi) - V(psi), prefactor m/beta^2 (r = 1..4, both couplings, types I_-, I_+)', ok_int)
report(f'on the defect conditions d(E + B_d)/dt = 0 and d(P + P_d)/dt = 0, P_d = -(m/b^2) sum[sigma(f-1) - (g-1)/sigma] (type I_-)  [res {worst:.1e}]', ok_em)

# ---------------------------------------------------------------- 4. Lax compatibility
def ltilde_space(r, lam, sigma, coupling, typ, beta=0.8, m=1.1, nconf=4):
    """Null space of X Lt - Lt Y = 0 (Lt = Ltilde), stacked over random configurations on the defect conditions."""
    d = an_data(r); U = d['U']; D, E = dmatrix(r, typ)
    b = beta if coupling == 'real' else 1j*beta
    h = r + 1
    rows = []
    for _ in range(nconf):
        phi = rng.normal(size=r) + 0.2j*rng.normal(size=r); psi = rng.normal(size=r) + 0.2j*rng.normal(size=r)
        phit = rng.normal(size=r); psit = rng.normal(size=r)
        _, gp, gq, _ = defect_potential(phi, psi, sigma, r, beta, m, coupling, typ)
        phix, psix = defect_dx(phit, psit, gp, gq, r, typ)
        u, w = E @ phi + D @ psi, D.T @ phi - E @ psi
        ud, wd = E @ phit + D @ psit, D.T @ phit - E @ psit
        Gu = np.diag(np.exp(-b*(U @ u)/2)); Gw = np.diag(np.exp(b*(U @ w)/2))
        X = np.linalg.inv(Gu) @ a_t(lam, phi, phix, r, beta, m, coupling) @ Gu - (b/2)*np.diag(U @ ud)
        Y = Gw @ a_t(lam, psi, psix, r, beta, m, coupling) @ np.linalg.inv(Gw) - (b/2)*np.diag(U @ wd)
        rows.append(np.kron(np.eye(h), X) - np.kron(Y.T, np.eye(h)))      # vec(X L - L Y), column-major
    M = np.vstack(rows)
    _, s, Vh = np.linalg.svd(M)
    null = [Vh[i].conj().reshape(h, h).T for i in range(len(s)) if s[i] < 1e-9*s[0]]
    return null


ok_dim = True; ok_form = True; cvals = []
for r in (1, 2, 3):
    Es = an_lax_rep(r); Sm = sum(E.T for E in Es); Sp = sum(Es); h = r + 1
    for coupling in ('real', 'imag'):
        for typ in (1, -1):
            for lam in (0.7 + 0.3j, -1.9 + 0.4j):
                sig = 0.83
                null = ltilde_space(r, lam, sig, coupling, typ)
                ok_dim &= len(null) == 1
                if len(null) != 1:
                    continue
                L = null[0]/null[0][0, 0]
                if typ == 1:                   # Ltilde = 1 + c (sigma/lambda) sum E_-alpha
                    Lc, gen = L - np.eye(h), Sm*sig/lam
                else:                          # Ltilde^-1 ~ 1 - c (lambda/sigma) sum E_alpha  (normalized on the diagonal)
                    Li = np.linalg.inv(L); Li = Li/Li[0, 0]
                    Lc, gen = Li - np.eye(h), -Sp*lam/sig
                c = np.vdot(gen.ravel(), Lc.ravel())/np.vdot(gen.ravel(), gen.ravel())
                cvals.append(c)
                ok_form &= np.allclose(Lc, c*gen, atol=1e-8)
report('Ltilde unique up to scale (a_1, a_2, a_3; real and imaginary coupling; types I_-, I_+; two lambda)', ok_dim)
report(f'type I_-: Ltilde = 1 + c sigma lambda^-1 sum_j E_-alpha_j;  type I_+ (B_d^{{I+}}(sigma) = B_d^{{I-}}(1/sigma)): '
       f'Ltilde^-1 ~ 1 - c sigma^-1 lambda sum_j E_alpha_j;  c = 1  [max |c - 1| = {max(abs(np.array(cvals) - 1)):.1e}]',
       ok_form and np.allclose(cvals, 1, atol=1e-8))
# control: a non-integrable D (E = 1 - D antisymmetric, but not 2(1 - pi)^-1) admits no field-independent Ltilde
import defect_common as dc
_orig = dc.dmatrix
def _bad(rr, typ=1):
    D, E = _orig(rr, typ)
    return (np.eye(rr) - 0.6*E, 0.6*E) if rr > 1 else (D, E)
dc.dmatrix = _bad
globals()['dmatrix'] = _bad
nbad = [len(ltilde_space(r, 0.7 + 0.3j, 0.83, 'real', 1)) for r in (2, 3)]
dc.dmatrix = _orig; globals()['dmatrix'] = _orig
report('control: with E -> 0.6 E (D = 1 - E, still D + D^T = 2) there is no Ltilde (a_2, a_3)', all(n == 0 for n in nbad))

# lambda of a soliton: the vacuum Lax pair gives E = Psi_1/Psi_0 = exp(m0(x cosh theta - t sinh theta)) at lambda = -e^{-theta}
th0 = 0.37; lam_s = -np.exp(-th0); mm = 1.1
P = np.array([[0, 1], [1, 0]])
okl = True
for mu in (1, -1):                      # eigenvalue of Lambda_+ = Lambda_- = sigma_x on (1, mu)/sqrt2
    a = (mm/2)*mu*(1/lam_s - lam_s); bx = -(mm/2)*mu*(lam_s + 1/lam_s)         # Psi = exp(a t + bx x)
    # check: d_+ Psi = -A_+ Psi = -m lam mu Psi, d_- Psi = -A_- Psi = (m/lam) mu Psi
    okl &= np.isclose(a + bx, -mm*lam_s*mu) and np.isclose(a - bx, mm/lam_s*mu)
ratio_x = -(mm/2)*(lam_s + 1/lam_s)*2; ratio_t = (mm/2)*(1/lam_s - lam_s)*2
okl &= np.isclose(ratio_x, 2*mm*np.cosh(th0)) and np.isclose(ratio_t, -2*mm*np.sinh(th0))
report('a_1 vacuum Lax solutions: the soliton exponent m0(x cosh theta - t sinh theta) arises at lambda = -e^{-theta}', okl)

# ---------------------------------------------------------------- 5. spin +-2 (and +-1) charges of a_2 with a defect
r = 2
s2, s3 = sp.sqrt(2), sp.sqrt(3)
al = [sp.Matrix([-1/s2, -s3/s2]), sp.Matrix([s2, 0]), sp.Matrix([-1/s2, s3/s2])]        # alpha_0, alpha_1, alpha_2
ok_roots = np.allclose(np.array([[float(v) for v in a] for a in al]), an_data(2)['alpha'])
S = sp.Matrix.hstack(al[1], al[2])
lamv = (S.T).inv()                                     # columns lambda_1, lambda_2 (alpha_i.lambda_j = delta_ij)
lams = [sp.zeros(2, 1), lamv[:, 0], lamv[:, 1], sp.zeros(2, 1)]
Dm = 2*sum(((lams[i] - lams[i + 1])*lams[i].T for i in (1, 2)), sp.zeros(2, 2))
Dm = sp.simplify(Dm); Em = sp.eye(2) - Dm
ok_roots &= np.allclose(np.array(Dm, dtype=float), dmatrix(2)[0])
# bulk spin-2 densities: T = A_abc u_a u_b u_c + B_ab w_a u_b, Theta = sum_j (g_j.u) e^{alpha_j.phi}, g_j = -B alpha_j
Asym = {}
for idx in itertools.combinations_with_replacement(range(2), 3):
    Asym[idx] = sp.Symbol('A' + ''.join(map(str, idx)))
A = lambda a, b, c: Asym[tuple(sorted((a, b, c)))]
b12 = sp.Symbol('b12'); Bm = sp.Matrix([[0, b12], [-b12, 0]])
uu = sp.Matrix(sp.symbols('u0:2'))
eqs = []
for aj in al:
    lhs = 3*sum(A(a, b, c)*uu[a]*uu[b]*aj[c] for a in range(2) for b in range(2) for c in range(2))
    rhs = -2*(aj.T*Bm*uu)[0]*(aj.T*uu)[0]
    eqs += sp.Poly(sp.expand(lhs - rhs), *uu).coeffs()
solA = sp.solve(eqs, list(Asym.values()), dict=True)[0]
Avals = {k: sp.simplify(v.subs(solA).subs(b12, 1)) for k, v in Asym.items()}
Bnum = Bm.subs(b12, 1)
report('a_2 roots and D in the sympy check agree with defect_common; spin-2 bulk densities exist (A determined by B)',
       ok_roots and any(v != 0 for v in Avals.values()))

q0 = sp.symbols('f0 f1 g0 g1'); q1 = sp.symbols('fd0 fd1 gd0 gd1'); q2 = sp.symbols('fdd0 fdd1 gdd0 gdd1')
q3 = sp.symbols('f3_0 f3_1 g3_0 g3_1'); q4 = sp.symbols('f4_0 f4_1 g4_0 g4_1')
Q = [q0, q1, q2, q3, q4]
SIG = sp.Rational(7, 10)


def Dt(e):
    return sum(Q[k + 1][a]*sp.diff(e, Q[k][a]) for k in range(4) for a in range(4))


def boundary_G(eps, typ, spin, Dbase=None):
    Dm_ = Dm if Dbase is None else Dbase
    Em_ = sp.eye(2) - Dm_
    D_, E_ = (Dm_, Em_) if typ == 1 else (-Dm_, -Em_)
    s = SIG if typ == 1 else 1/SIG
    phi, psi = sp.Matrix(q0[:2]), sp.Matrix(q0[2:]); phid, psid = sp.Matrix(q1[:2]), sp.Matrix(q1[2:])
    phidd, psidd = sp.Matrix(q2[:2]), sp.Matrix(q2[2:])
    B = sum(s*(sp.exp((aj.T*(Dm_.T*phi + Dm_*psi))[0]/2) - 1) + (sp.exp((aj.T*Dm_*(phi - psi))[0]/2) - 1)/s for aj in al)
    gp = sp.Matrix([sp.diff(B, v) for v in phi]); gq = sp.Matrix([sp.diff(B, v) for v in psi])
    Fp = E_*phid + D_*psid - gp
    Fq = D_.T*phid - E_*psid + gq
    gV = lambda f: sum((aj*sp.exp((aj.T*f)[0]) for aj in al), sp.zeros(2, 1))
    V = lambda f: sum(sp.exp((aj.T*f)[0]) - 1 for aj in al)
    out = 0
    for sgn, f, fd, fdd, F in ((1, phi, phid, phidd, Fp), (-1, psi, psid, psidd, Fq)):
        u = fd + eps*F
        if spin == 1:
            TT = (u.T*u)[0]/2 - V(f)                       # T + Theta for the spin-1 charge (Theta = -V)
        else:
            Fd = sp.Matrix([Dt(c) for c in F])
            w = fdd + 2*eps*Fd + fdd + gV(f)               # d_eps^2 f = f_tt + 2 eps f_tx + f_xx, f_xx = f_tt + grad V
            T3 = sum(Avals[tuple(sorted((a, b, c)))]*u[a]*u[b]*u[c] for a in range(2) for b in range(2) for c in range(2))
            T3 += (w.T*Bnum*u)[0]
            Th = sum(((-Bnum*aj).T*u)[0]*sp.exp((aj.T*f)[0]) for aj in al)
            TT = T3 + Th
        out += sgn*TT
    return out


def euler_residual(G, npts=2):
    """max |Euler operator in t of G| / max |dG/dq| at random points (G = d/dt F, F local, iff the Euler operator vanishes)."""
    E = [sp.diff(G, q0[a]) - Dt(sp.diff(G, q1[a])) + Dt(Dt(sp.diff(G, q2[a]))) for a in range(4)]
    allv = [v for qq in Q for v in qq]
    f = sp.lambdify(allv, E, 'numpy', cse=True)
    fs = sp.lambdify(allv, [sp.diff(G, q0[a]) for a in range(4)], 'numpy', cse=True)
    worst = 0; scale = 0
    for _ in range(npts):
        vals = rng.normal(size=len(allv))*0.5
        worst = max(worst, np.abs(np.array(f(*vals), dtype=float)).max())
        scale = max(scale, np.abs(np.array(fs(*vals), dtype=float)).max())
    return worst/scale


# bulk: d_{-eps} T = d_eps Theta on shell (u = d_eps phi, w = d_eps^2 phi, d_{-eps} u = -grad V, d_{-eps} w = -Hess V u)
fv = sp.Matrix(sp.symbols('x0:2')); ww = sp.Matrix(sp.symbols('w0:2')); w3 = sp.Matrix(sp.symbols('w3_0:2')); vv = sp.Matrix(sp.symbols('v0:2'))
Vb = sum(sp.exp((aj.T*fv)[0]) - 1 for aj in al); gVb = sp.Matrix([sp.diff(Vb, c) for c in fv]); HVb = sp.hessian(Vb, fv)
Tb = sum(Avals[tuple(sorted((a, b, c)))]*uu[a]*uu[b]*uu[c] for a in range(2) for b in range(2) for c in range(2)) + (ww.T*Bnum*uu)[0]
Thb = sum(((-Bnum*aj).T*uu)[0]*sp.exp((aj.T*fv)[0]) for aj in al)
d_same = lambda e: sum(sp.diff(e, fv[a])*uu[a] + sp.diff(e, uu[a])*ww[a] + sp.diff(e, ww[a])*w3[a] for a in range(2))
d_other = lambda e: sum(sp.diff(e, fv[a])*vv[a] - sp.diff(e, uu[a])*gVb[a] - sp.diff(e, ww[a])*(HVb*uu)[a] for a in range(2))
fb = sp.lambdify(list(fv) + list(uu) + list(ww) + list(w3) + list(vv), d_other(Tb) - d_same(Thb))
bulk_res = max(abs(fb(*rng.normal(size=10))) for _ in range(3))
# non-triviality: the static density (u = eps phi_x, w = phi_xx) is not a total x-derivative (Euler operator in x != 0)
X0, X1, X2, X3, X4 = [sp.Matrix(sp.symbols(f'y{k}_0:2')) for k in range(5)]
Xs = [X0, X1, X2, X3, X4]
Dx = lambda e: sum(Xs[k + 1][a]*sp.diff(e, Xs[k][a]) for k in range(4) for a in range(2))
dens = (Tb - Thb).subs({**{uu[a]: X1[a] for a in range(2)}, **{ww[a]: X2[a] for a in range(2)}, **{fv[a]: X0[a] for a in range(2)}})
eulx = [sp.diff(dens, X0[a]) - Dx(sp.diff(dens, X1[a])) + Dx(Dx(sp.diff(dens, X2[a]))) for a in range(2)]
fx = sp.lambdify([c for X in Xs for c in X], eulx)
nontriv = max(np.abs(fx(*rng.normal(size=10))).max() for _ in range(3))
# BCZ (5.11): the cubic term A(Dv)^3 + A(D^T v)^3 of the 'wrong' chirality vanishes identically (cot identity)
Dn = np.array(Dm, dtype=float)
An = np.zeros((2, 2, 2))
for idx in itertools.product(range(2), repeat=3):
    An[idx] = float(Avals[tuple(sorted(idx))])
cub = lambda w: np.einsum('abc,a,b,c->', An, w, w, w)
vr = rng.normal(size=2)
report(f'a_2 spin-2 bulk current conserved on shell [res {bulk_res:.1e}]; its density is not a total x-derivative '
       f'[Euler_x {nontriv:.1e}]; BCZ (5.11) cubic A(Dv)^3 + A(D^T v)^3 = {cub(Dn @ vr) + cub(Dn.T @ vr):.1e} (vanishes)',
       bulk_res < 1e-10 and nontriv > 1e-2 and abs(cub(Dn @ vr) + cub(Dn.T @ vr)) < 1e-12)


def loop_integral(G, ncoef=3, npts=400):
    """|int_0^2pi G dt| / int |G| dt along a random periodic path q(t) (a total time derivative integrates to 0)."""
    f = sp.lambdify([c for qq in Q[:3] for c in qq], G)
    tt = np.linspace(0, 2*np.pi, npts, endpoint=False)
    c0 = rng.normal(size=4)*0.4; cs = rng.normal(size=(ncoef, 4))*0.3; sn = rng.normal(size=(ncoef, 4))*0.3
    k = np.arange(1, ncoef + 1)[:, None, None]
    q_ = c0[:, None] + (cs[:, :, None]*np.cos(k*tt) + sn[:, :, None]*np.sin(k*tt)).sum(0)
    qd = (-cs[:, :, None]*k*np.sin(k*tt) + sn[:, :, None]*k*np.cos(k*tt)).sum(0)
    qdd = (-cs[:, :, None]*k**2*np.cos(k*tt) - sn[:, :, None]*k**2*np.sin(k*tt)).sum(0)
    g = np.array(f(*q_, *qd, *qdd), dtype=float)
    return abs(g.mean())/np.abs(g).mean()


res = {}; loops = {}
for typ in (1, -1):
    for eps in (1, -1):
        G2 = boundary_G(eps, typ, 2)
        loops[(typ, eps)] = max(loop_integral(G2) for _ in range(3))
        if typ == 1:                                   # full Euler test for type I_- (about 30 s each)
            res[(typ, eps, 2)] = euler_residual(G2)
res[(1, 1, 1)] = euler_residual(boundary_G(1, 1, 1)); res[(1, -1, 1)] = euler_residual(boundary_G(-1, 1, 1))
Dbad = sp.Matrix([[1, sp.Rational(3, 10)], [-sp.Rational(3, 10), 1]])
bad_l = max(loop_integral(boundary_G(1, 1, 2, Dbad)) for _ in range(2))
bad_l2 = max(loop_integral(boundary_G(-1, 1, 2, Dbad)) for _ in range(2))
tol = 1e-9
report(f'control: spin +-1 (E +- P) defect contributions are total time derivatives (type I_-)  '
       f'[Euler residuals {res[(1, 1, 1)]:.1e}, {res[(1, -1, 1)]:.1e}]', res[(1, 1, 1)] < tol and res[(1, -1, 1)] < tol)
report(f'control: non-integrable D (E = [[0, -3/10], [3/10, 0]]) spoils both spin-2 charges  [closed-loop {bad_l:.1e}, {bad_l2:.1e}]',
       bad_l > 1e-3 and bad_l2 > 1e-3)
okall = all(res[(1, e_, 2)] < tol for e_ in (1, -1)) and all(loops[(t_, e_)] < 1e-10 for t_ in (1, -1) for e_ in (1, -1))
report('BOTH spin +2 (d_- densities) and spin -2 (d_+ densities) charges acquire local defect corrections, for type I_- AND '
       f'type I_+ (type I_- Euler residuals {res[(1, 1, 2)]:.0e}, {res[(1, -1, 2)]:.0e}; closed-loop integrals ' +
       ', '.join(f'{loops[(t_, e_)]:.0e}' for t_ in (1, -1) for e_ in (1, -1)) + ') -- contrary to BCZ hep-th/0401020 sec. 5', okall)


# ================================================================ ch. 38 additions (sec-defect-conditions)
print('--- chapter 38 additions')
from scipy.optimize import least_squares


# ---------------------------------------------------------------- 6. the scan of prp-ar-only
def affine_gram(name):
    """Gram matrix alpha_i.alpha_j (i, j = 0..r) of the affine simple roots, longest root^2 = 2 (book normalization,
    sec-lie-algebraic-data).  name = 'a3(1)', 'd4(2)', ... (Kac's notation; rank r = number of nodes - 1)."""
    typ, rest = name.split('(')
    k = int(rest[0]); X = typ[0]; n = int(typ[1:])

    def chain(L, edges):
        G = np.diag(np.array(L, float))
        for i, j, g in edges:
            G[i, j] = G[j, i] = g
        return G
    if (X, k) == ('a', 1):
        return chain([2, 2], [(0, 1, -2)]) if n == 1 else chain([2]*(n + 1), [(j, (j + 1) % (n + 1), -1) for j in range(n + 1)])
    if (X, k) == ('b', 1):                       # 0, 1 attached to 2; chain; (n-1) => n short
        return chain([2]*n + [1], [(0, 2, -1), (1, 2, -1)] + [(j, j + 1, -1) for j in range(2, n - 1)] + [(n - 1, n, -1)])
    if (X, k) == ('c', 1):                       # long ends 0, n; short chain
        return chain([2] + [1]*(n - 1) + [2], [(0, 1, -1)] + [(j, j + 1, -0.5) for j in range(1, n - 1)] + [(n - 1, n, -1)])
    if (X, k) == ('d', 1):                       # d_4^(1)
        return chain([2]*5, [(0, 2, -1), (1, 2, -1), (2, 3, -1), (2, 4, -1)])
    if (X, k) == ('f', 1):
        return chain([2, 2, 2, 1, 1], [(0, 1, -1), (1, 2, -1), (2, 3, -1), (3, 4, -0.5)])
    if (X, k) == ('g', 1):
        return chain([2, 2, 2/3], [(0, 1, -1), (1, 2, -1)])
    if (X, k) == ('a', 2) and n % 2 == 0:        # a_{2m}^(2): squared lengths 2, 1, ..., 1, 1/2 (sec-affine-folding)
        mm = n//2
        if mm == 1:
            return chain([2, 0.5], [(0, 1, -1)])
        return chain([2] + [1]*(mm - 1) + [0.5], [(0, 1, -1)] + [(j, j + 1, -0.5) for j in range(1, mm - 1)] + [(mm - 1, mm, -0.5)])
    if (X, k) == ('a', 2):                       # a_{2m-1}^(2), m >= 3: fork of short nodes, long end
        mm = (n + 1)//2
        return chain([1]*mm + [2], [(0, 2, -0.5), (1, 2, -0.5)] + [(j, j + 1, -0.5) for j in range(2, mm - 1)] + [(mm - 1, mm, -1)])
    if (X, k) == ('d', 2):                       # d_{m+1}^(2): short ends, long chain
        mm = n - 1
        return chain([1] + [2]*(mm - 1) + [1], [(0, 1, -1)] + [(j, j + 1, -1) for j in range(1, mm - 1)] + [(mm - 1, mm, -1)])
    if (X, k) == ('e', 2):
        return chain([1, 1, 1, 2, 2], [(0, 1, -0.5), (1, 2, -0.5), (2, 3, -1), (3, 4, -1)])
    if (X, k) == ('d', 3):
        return chain([2/3, 2/3, 2], [(0, 1, -1/3), (1, 2, -1)])
    raise ValueError(name)


def realize(G):
    """Root vectors (rows, in R^r) with Gram matrix G, and the coprime Kac labels (sum_j n_j alpha_j = 0)."""
    from math import gcd
    from functools import reduce
    N = len(G); r_ = N - 1
    w, V = np.linalg.eigh(G)
    nv = np.abs(V[:, np.argmin(abs(w))]); nv = nv/nv.min()
    for s in range(1, 13):
        if np.allclose(s*nv, np.round(s*nv), atol=1e-8):
            nv = np.round(s*nv).astype(int); break
    nv = nv//reduce(gcd, nv)
    i0 = int(np.argmin(nv)); others = [j for j in range(N) if j != i0]
    A = np.zeros((N, r_)); A[others] = np.linalg.cholesky(G[np.ix_(others, others)])
    A[i0] = -sum(nv[j]*A[j] for j in others)/nv[i0]
    assert np.allclose(A @ A.T, G, atol=1e-10)
    return A, nv


def scan_momentum(name, b=0.7, mm=1.3):
    """For every automorphism p of the affine diagram: the type-I exponential ansatz (cz2009a, BCZ hep-th/0401020
    (3.19)-(3.24)) forces D = 2(1 - P)^{-1}, P alpha_j = alpha_p(j), with p fixed-point free, P orthogonal and 1 - P
    invertible; then alpha_i.D alpha_j must vanish for j not in {i, p(i)}, and the momentum identity
    1/2|grad_phi B|^2 - 1/2|grad_psi B|^2 = V(phi) - V(psi) is tested with B = (m/b^2) sum_j c_j[s f_j + g_j/s]
    (c_j^2 = 2 n_j/alpha_j^2, one s per p-orbit) at random points.  Returns a list of (perm, reason or None)."""
    G = affine_gram(name); A, nv = realize(G); N, r_ = A.shape
    out = []
    for p in itertools.permutations(range(N)):
        p = list(p)
        if not np.allclose(G[np.ix_(p, p)], G):
            continue
        if any(p[j] == j for j in range(N)):
            out.append((tuple(p), 'fixed node', None)); continue
        P = np.linalg.lstsq(A, A[p], rcond=None)[0].T
        if not (np.allclose(A @ P.T, A[p]) and np.allclose(P @ P.T, np.eye(r_))):
            out.append((tuple(p), 'P not orthogonal', None)); continue
        if abs(np.linalg.det(np.eye(r_) - P)) < 1e-9:
            out.append((tuple(p), '1 - P singular', None)); continue
        D = 2*np.linalg.inv(np.eye(r_) - P)
        Gd = A @ D @ A.T
        nbad = sum(1 for i in range(N) for j in range(N) if j not in (i, p[i]) and abs(Gd[i, j]) > 1e-9)
        orb = {}; lab = 0
        for j in range(N):
            if j in orb:
                continue
            k = j
            while k not in orb:
                orb[k] = lab; k = p[k]
            lab += 1
        sig = 0.5 + rng.random(lab); worst = 0
        for _ in range(4):
            ph = rng.normal(size=r_)*0.4; ps = rng.normal(size=r_)*0.4
            gp = np.zeros(r_); gq = np.zeros(r_); Vp = Vq = 0
            for i in range(N):
                a = A[i]; c = np.sqrt(2*nv[i]/(a @ a)); s = sig[orb[i]]
                f = np.exp(b*(D @ a) @ ph/2 + b*(D.T @ a) @ ps/2); g = np.exp(b*(D.T @ a) @ (ph - ps)/2)
                gp += (mm/b)*c*(s*f*(D @ a) + g*(D.T @ a)/s)/2
                gq += (mm/b)*c*(s*f*(D.T @ a) - g*(D.T @ a)/s)/2
                Vp += (mm**2/b**2)*nv[i]*(np.exp(b*a @ ph) - 1); Vq += (mm**2/b**2)*nv[i]*(np.exp(b*a @ ps) - 1)
            worst = max(worst, abs((gp @ gp - gq @ gq)/2 - (Vp - Vq))/(1 + abs(Vp) + abs(Vq)))
        out.append((tuple(p), None if (nbad == 0 and worst < 1e-10) else f'{nbad} cross terms', worst))
    return out


names = ['a1(1)', 'a2(1)', 'a3(1)', 'a4(1)', 'b3(1)', 'b4(1)', 'c2(1)', 'c3(1)', 'c4(1)', 'd4(1)', 'f4(1)', 'g2(1)',
         'a2(2)', 'a4(2)', 'a6(2)', 'a8(2)', 'd3(2)', 'a5(2)', 'a7(2)', 'd4(2)', 'd5(2)', 'e6(2)', 'd4(3)']
ok_scan = True; npass = 0; naut = 0; reasons = {}; worst_pass = 0; best_fail = np.inf
for nm in names:
    res_ = scan_momentum(nm)
    naut += len(res_)
    N_ = len(affine_gram(nm))
    rot = {tuple((j + 1) % N_ for j in range(N_)), tuple((j - 1) % N_ for j in range(N_))}   # p(j) = j +- 1: pi^{-+1}
    for perm, why, w in res_:
        reasons[why if why is None else why.split(' ', 1)[-1] if 'cross' in why else why] = reasons.get(
            why if why is None else why.split(' ', 1)[-1] if 'cross' in why else why, 0) + 1
        expected = nm.endswith('(1)') and nm[0] == 'a' and perm in rot
        ok_scan &= (why is None) == expected
        if why is None:
            npass += 1; worst_pass = max(worst_pass, w)
        elif w is not None:
            best_fail = min(best_fail, w)
report(f'prp-ar-only scan: {len(names)} affine diagrams of rank <= 4 (untwisted and twisted), {naut} automorphisms; only '
       f'a_r^(1) with p = pi^(+-1) passes ({npass} cases, momentum residual {worst_pass:.0e}); failures: '
       + ', '.join(f'{k}: {v}' for k, v in reasons.items() if k is not None) + f' (smallest residual among them {best_fail:.2f})',
       ok_scan and worst_pass < 1e-10 and best_fail > 0.1)

# cz2009a (2.7)-(2.9) in book form: with D = 1 - E (E antisymmetric) the momentum flux has no quadratic velocity terms, and
# the linear terms are -dP_d/dt exactly when grad_phi P_d = E grad_phi B - D grad_psi B, grad_psi P_d = -D^T grad_phi B - E grad_psi B
ok = True
for r in (2, 3, 4):
    Ar = rng.normal(size=(r, r)); Ar = Ar - Ar.T; O, _ = np.linalg.qr(rng.normal(size=(r, r))); Bm_ = (np.eye(r) - Ar) @ O
    ok &= np.allclose(np.eye(r) - Ar @ Ar, Bm_ @ Bm_.T)                 # cz2009a (2.7)/(2.10): B = (1 - A) O
    D, E = dmatrix(r)
    ok &= np.allclose(np.eye(r) - E @ E, D @ D.T) and np.allclose(D @ D.T, D.T @ D) and np.allclose(E @ D, D @ E)
    for typ in (1, -1):
        D, E = dmatrix(r, typ)
        phi = rng.normal(size=r)*0.4; psi = rng.normal(size=r)*0.4
        _, gp, gq, _ = defect_potential(phi, psi, 0.8, r, 0.7, 1.3, 'real', typ)
        g1, g2 = grad_Pd(phi, psi, 0.8, r, 0.7, 1.3, 'real', typ)
        ok &= np.allclose(g1, E @ gp - D @ gq, atol=1e-7) and np.allclose(g2, -D.T @ gp - E @ gq, atol=1e-7)
report('cz2009a (2.7), (2.9) in book form: B = (1 - A)O solves 1 - A^2 = B B^T; D = 1 - E normal with 1 - E^2 = D D^T; '
       'grad_phi P_d = E grad_phi B - D grad_psi B, grad_psi P_d = -D^T grad_phi B - E grad_psi B (r = 2..4, both types)', ok)


# ---------------------------------------------------------------- 7. Baecklund: defect conditions for all x (BCZ sec. 5.1)
def backlund_residual(r, typ, coupling, Dpair=None, beta=0.7, m=1.3, sigma=0.8):
    """Cross-differentiate phi_x = E phi_t + D psi_t - grad_phi B, psi_x = D^T phi_t - E psi_t + grad_psi B (all x) and
    compare phi_xx - phi_tt, psi_xx - psi_tt with grad V(phi), grad V(psi); random point, random phi_tt, psi_tt."""
    d = an_data(r); al = d['alpha']
    D, E = dmatrix(r, typ) if Dpair is None else Dpair
    D1, _ = dmatrix(r, 1) if Dpair is None else (Dpair[0]*typ, None)
    b = beta if coupling == 'real' else 1j*beta
    s = sigma if typ == 1 else 1/sigma
    z = rng.normal(size=2*r)*0.4 + 0.2j*rng.normal(size=2*r); zt = rng.normal(size=2*r); ztt = rng.normal(size=2*r)
    grad = np.zeros(2*r, complex); hess = np.zeros((2*r, 2*r), complex)
    for a in al:
        for coef, vec in ((s, np.concatenate([D1 @ a, D1.T @ a])), (1/s, np.concatenate([D1.T @ a, -D1.T @ a]))):
            e = np.exp(b*(vec @ z)/2)
            grad += (m/b**2)*coef*(b/2)*vec*e; hess += (m/b**2)*coef*(b/2)**2*np.outer(vec, vec)*e
    # x-derivatives as functions of (z, z_t):  z_x = M z_t + J grad B, M = [[E, D], [D^T, -E]], J = diag(-1, +1)
    M = np.block([[E, D], [D.T, -E]]); J = np.diag(np.r_[-np.ones(r), np.ones(r)])
    zx = M @ zt + J @ grad
    zxt = M @ ztt + J @ hess @ zt
    zxx = M @ zxt + J @ hess @ zx
    gV = lambda f: (m**2/b)*sum(a*np.exp(b*(a @ f)) for a in al)
    want = np.concatenate([gV(z[:r]), gV(z[r:])])
    return np.abs(zxx - ztt - want).max()/(1 + np.abs(want).max())


ok = True; worst = 0
for r in (1, 2, 3, 4):
    for typ in (1, -1):
        for coupling in ('real', 'imag'):
            w = backlund_residual(r, typ, coupling); worst = max(worst, w)
Dg, Eg = dmatrix(3); Ebad = 0.6*Eg
bad = backlund_residual(3, 1, 'real', Dpair=(np.eye(3) - Ebad, Ebad))
lam_ = an_data(3)['lam']
E542 = sum(np.outer(lam_[i + 1], lam_[i]) - np.outer(lam_[i], lam_[i + 1]) for i in range(1, 4))
D543 = 2*sum(np.outer(lam_[i], lam_[i] - lam_[i + 1]) for i in range(1, 4))
mis = np.allclose(np.eye(3) - E542, Dg) and not np.allclose(D543, Dg) and np.allclose(D543, Dg.T)
report(f'BCZ sec. 5.1: the defect conditions imposed for all x are a Baecklund transformation (cross-differentiation gives '
       f'both field equations; r = 1..4, both types and couplings) [res {worst:.0e}; non-integrable D: {bad:.2f}]; BCZ (5.43) '
       f'as printed is D^T, i.e. 1 - E of (5.42) is the D of (3.28)', worst < 1e-12 and bad > 1e-2 and mis)


# ---------------------------------------------------------------- 8. delay of an a_r soliton = Baecklund (Bianchi) coefficient
BETA_, M_ = 0.8, 1.1


def an_soliton(r, a, theta, xi, z=1.0):
    """Book soliton (eq-an-soliton, eq-hirota-ansatz, imaginary coupling) with E -> zE: (phi, phi_x, phi_t)(x, t)."""
    d = an_data(r); h = r + 1; al = d['alpha']; om = np.exp(2j*np.pi/h)
    ma = 2*M_*np.sin(np.pi*a/h)

    def f(x, t):
        Ex = z*np.exp(ma*(x*np.cosh(theta) - t*np.sinh(theta)) + xi)
        tau = np.array([1 + om**(j*a)*Ex for j in range(h)]); dt_ = np.array([om**(j*a)*Ex for j in range(h)])
        ph = (1j/BETA_)*sum(al[j]*np.log(tau[j]) for j in range(h))
        phx = (1j/BETA_)*sum(al[j]*dt_[j]*ma*np.cosh(theta)/tau[j] for j in range(h))
        pht = (1j/BETA_)*sum(al[j]*dt_[j]*(-ma*np.sinh(theta))/tau[j] for j in range(h))
        return ph, phx, pht
    return f


def delay_vec(r, a, theta, xi, z, sigma, typ=1, ts=np.linspace(-2, 2, 7), xs=(0.0,), theta_map=False):
    """Residuals of the defect conditions for (phi, psi) = (soliton, soliton with E -> zE) at the points (x, t)."""
    F = an_soliton(r, a, theta, xi); Gs = an_soliton(r, a, theta, xi, z); out = []
    for x in xs:
        for t in ts:
            ph, phx, pht = F(x, t); ps, psx, pst = Gs(x, t)
            if theta_map:                                   # Theta: (phi, psi) -> (-conj phi, -conj psi)
                ph, phx, pht, ps, psx, pst = [-np.conj(v) for v in (ph, phx, pht, ps, psx, pst)]
            _, gp, gq, _ = defect_potential(ph, ps, sigma, r, BETA_, M_, 'imag', typ)
            px, sx = defect_dx(pht, pst, gp, gq, r, typ)
            out += list(phx - px) + list(psx - sx)
    return np.array(out)


def delay_residual(*args, **kw):
    return np.abs(delay_vec(*args, **kw)).max()


def z_book(a, h, theta, thd):          # = cz2009a (4.10) = BCZ hep-th/0401020 (6.2) at theta -> -theta
    x = np.exp(-(theta - thd)); g = np.pi*a/h
    return (x - 1j*np.exp(1j*g))/(x - 1j*np.exp(-1j*g))


def z_cz2007(a, h, theta, eta):         # cz2007 (2.20), 'clockwise' D = book D^T
    x = np.exp(-(theta - eta)); g = np.pi*a/h
    return (x + 1j*np.exp(-1j*g))/(x + 1j*np.exp(1j*g))


def z_bcz62(a, h, theta_bcz, d):        # BCZ hep-th/0401020 (6.2), their E = exp(m_a(x cosh th + t sinh th))
    g = np.pi*a/h
    return (1j*np.exp(theta_bcz) + d*np.exp(1j*g))/(1j*np.exp(theta_bcz) + d*np.exp(-1j*g))


thd = 0.3; sig_ = np.exp(-thd); worst_f = worst_fit = worst_dict = 0; okf = True
for r, a in [(2, 1), (2, 2), (3, 2), (1, 1)]:
    h = r + 1
    for theta in (-0.4, 0.9, 1.7):
        xi = 0.2 + 0.6j; zf = z_book(a, h, theta, thd)
        worst_f = max(worst_f, delay_residual(r, a, theta, xi, zf, sig_, xs=(0.0, -0.7, 1.3)))
        if theta == 0.9 and r == 2:                         # independent determination of z by least squares
            fun = lambda v, r=r, a=a, theta=theta, xi=xi: (lambda w: np.r_[w.real, w.imag])(
                delay_vec(r, a, theta, xi, v[0] + 1j*v[1], sig_, ts=np.linspace(-2, 2, 4)))
            sol = least_squares(fun, [1.0, 0.1], xtol=1e-15, ftol=1e-15, gtol=1e-15)
            worst_fit = max(worst_fit, abs(sol.x[0] + 1j*sol.x[1] - zf))
        worst_dict = max(worst_dict, abs(zf - z_cz2007(h - a, h, theta, thd)), abs(zf - z_bcz62(a, h, -theta, sig_)))
report(f'a_r soliton through the type-I_- defect (book conventions): psi = soliton with E -> z_a E, '
       f'z_a = (e^-(theta-theta_d) - i e^(i pi a/h))/(e^-(theta-theta_d) - i e^(-i pi a/h)) = cz2009a (4.10) = BCZ (6.2) at '
       f'theta -> -theta = cz2007 (2.20) at a -> h - a; the pair solves the conditions at x = 0, -0.7, 1.3 (Baecklund at all x) '
       f'(a_2 a = 1, 2; a_3 a = 2; a_1) [res {worst_f:.0e}; least-squares z for a_2 differs by {worst_fit:.0e}; dictionary {worst_dict:.0e}]',
       worst_f < 1e-12 and worst_dict < 1e-12 and worst_fit < 1e-8)

# the delay as a (twisted) Bianchi coefficient, self-conjugate solitons, soliton-antisoliton pairs, Theta on the delay
ok = True
for h in (2, 3, 4, 5):
    for a in range(1, h):
        for theta, td in ((0.4, -0.2), (1.3 + 0.2j, 0.5), (-0.7, 0.9 + 0.3j)):
            A_, Ad = -np.exp(-theta), -np.exp(-td); eps_a = np.exp(1j*np.pi*(2*a - h)/(2*h)); g = np.pi*a/h; u = theta - td
            z = z_book(a, h, theta, td)
            ok &= np.isclose(z, (A_ + eps_a*Ad)/(A_ - Ad/eps_a))                                   # twisted Bianchi form
            ok &= np.isclose(z, np.exp(1j*g)*np.cosh((u + 1j*g - 1j*np.pi/2)/2)/np.cosh((u - 1j*g - 1j*np.pi/2)/2))
            ok &= np.isclose(z*z_book(h - a, h, theta, td), (np.cosh(u) + np.sin(g))/(np.cosh(u) - np.sin(g)))
            ok &= np.isclose(z_book(h - a, h, np.conj(theta), np.conj(td)), np.conj(z))           # Theta: (a, theta) -> (abar, conj theta)
        if 2*a == h:
            th = 0.9; zc = z_book(a, h, th, thd)
            ok &= np.isclose(zc, -1/np.tanh((th - thd)/2)) and np.isclose(z_cz2007(a, h, th, thd), zc)
            ok &= not np.isclose((1 + np.exp(-(th - thd)))/(1 - np.exp(-(th - thd))), zc)        # cz2007 (2.21) as printed
report('delay = twisted Bianchi coefficient z_a = (A(theta) + eps_a A(theta_d))/(A(theta) - A(theta_d)/eps_a), A = -e^-theta, '
       'eps_a = e^(i pi(2a-h)/2h) (sine-Gordon: eps = 1); z_a z_abar = (cosh u + sin(pi a/h))/(cosh u - sin(pi a/h)); '
       'z_abar(conj theta; conj theta_d) = conj z_a; self-conjugate a = h/2: z = -coth((theta - theta_d)/2), so cz2007 (2.21) '
       'as printed has the wrong sign (h = 2..5, complex rapidities)', ok)

# ---------------------------------------------------------------- 9. Theta (prp-defect-theta)
ok_alg = True; dev = 0
for r in (2, 3):
    for typ in (1, -1):
        for coupling in ('imag', 'real'):
            for sig in (0.7, -1.3, 0.8*np.exp(0.6j)):
                z1 = rng.normal(size=r) + 0.4j*rng.normal(size=r); z2 = rng.normal(size=r) + 0.4j*rng.normal(size=r)
                B0 = defect_potential(z1, z2, sig, r, 0.7, 1.3, coupling, typ)[0]
                img = (-np.conj(z1), -np.conj(z2)) if coupling == 'imag' else (np.conj(z1), np.conj(z2))
                Bc = defect_potential(*img, np.conj(sig), r, 0.7, 1.3, coupling, typ)[0]
                Bs = defect_potential(*img, sig, r, 0.7, 1.3, coupling, typ)[0]
                ok_alg &= np.isclose(Bc, np.conj(B0))
                if np.isreal(sig):
                    ok_alg &= np.isclose(Bs, np.conj(B0))
                else:
                    dev = max(dev, abs(Bs - np.conj(B0)))
# Theta-image of a soliton passing the defect: solves the conditions with sigma-bar
res_real = delay_residual(2, 1, 0.9, 0.2 + 0.6j, z_book(1, 3, 0.9, thd), sig_, theta_map=True)
tdc = 0.3 + 0.4j; sc = np.exp(-tdc); zc = z_book(1, 3, 0.9, tdc)
res_orig = delay_residual(2, 1, 0.9, 0.2 + 0.6j, zc, sc)
res_bar = delay_residual(2, 1, 0.9, 0.2 + 0.6j, zc, np.conj(sc), theta_map=True)
res_same = delay_residual(2, 1, 0.9, 0.2 + 0.6j, zc, sc, theta_map=True)
report(f'Theta: conj B_d(phi, psi; sigma) = B_d(-conj phi, -conj psi; conj sigma) (imaginary coupling; plain conjugation at real '
       f'coupling), a_2, a_3, both types; so Theta maps the sigma-theory to the conj(sigma)-theory and is a symmetry iff sigma is '
       f'real [complex sigma: deviation {dev:.2f}]; a_2 soliton through the defect: Theta-image solves the conditions '
       f'(real sigma {res_real:.0e}; complex sigma: with conj sigma {res_bar:.0e}, with sigma {res_same:.2f}; original {res_orig:.0e})',
       ok_alg and dev > 1e-2 and max(res_real, res_bar, res_orig) < 1e-12 and res_same > 1e-3)


# ---------------------------------------------------------------- 10. type-II defects (cz2009b units; prp-type-ii-momentum)
# u on x < 0, v on x > 0, lam the defect field, p = (u+v)/2, q = (u-v)/2;  L_D = delta(x)[ (u v_t - v u_t)/2 + lam (u-v)_t
# - lam_t (u-v) - D ];  conditions (2.2)-(2.4): u_x = v_t - 2 lam_t - D_u, v_x = u_t - 2 lam_t + D_v, u_t = v_t + D_lam/2,
# D = f + g, P + Omega conserved with P = int u_x u_t (cz convention, = -P_book) and Omega = f - g.
def t2_sg(sig, tau):
    gam = np.exp(2*tau) + np.exp(-2*tau)
    return (lambda u: np.exp(u) + np.exp(-u),
            lambda p, q, l: (2*np.exp(p - l) + np.exp(-p + l)*(np.exp(q) + np.exp(-q) + gam))/sig,
            lambda p, q, l: sig*(np.exp(l)*(np.exp(q) + np.exp(-q) + gam) + 2*np.exp(-l)))


def t2_tz(sig):                         # cz2009b (3.7) with the whole of f multiplied by 1/sigma (as (4.6) requires)
    return (lambda u: np.exp(2*u) + 2*np.exp(-u),
            lambda p, q, l: (np.exp(2*p - 2*l) + np.exp(-p + l)*(np.exp(q) + np.exp(-q)))/sig,
            lambda p, q, l: (sig/2)*(8*np.exp(-l) + np.exp(2*l)*(np.exp(q) + np.exp(-q))**2))


def _nd(F, X, i, h=1e-6):
    e = np.zeros(len(X), complex); e[i] = h
    return (F(*(X + e)) - F(*(X - e)))/(2*h)


def t2_conditions(model, u, v, l, ut, vt, lt):
    U, f, g = model
    D = lambda u, v, l: f((u + v)/2, (u - v)/2, l) + g((u + v)/2, (u - v)/2, l)
    X = np.array([u, v, l], complex)
    return vt - 2*lt - _nd(D, X, 0), ut - 2*lt + _nd(D, X, 1), ut - vt - _nd(D, X, 2)/2


def t2_em(model):
    U, f, g = model
    D = lambda u, v, l: f((u + v)/2, (u - v)/2, l) + g((u + v)/2, (u - v)/2, l)
    Om = lambda u, v, l: f((u + v)/2, (u - v)/2, l) - g((u + v)/2, (u - v)/2, l)
    worst = 0; scale = 0
    for _ in range(5):
        u, v, l = rng.normal(size=3)*0.5 + 0.3j*rng.normal(size=3); vt, lt = rng.normal(size=2)
        X = np.array([u, v, l], complex); ut = vt + _nd(D, X, 2)/2
        ux, vx, _ = t2_conditions(model, u, v, l, ut, vt, lt)
        dE = ut*ux - vt*vx + _nd(D, X, 0)*ut + _nd(D, X, 1)*vt + _nd(D, X, 2)*lt
        dP = (ut**2 + ux**2)/2 - U(u) - (vt**2 + vx**2)/2 + U(v)
        dP += _nd(Om, X, 0)*ut + _nd(Om, X, 1)*vt + _nd(Om, X, 2)*lt
        worst = max(worst, abs(dE), abs(dP)); scale = max(scale, abs(U(u)) + abs(ux)**2)
    return worst/scale


em_sg, em_tz = t2_em(t2_sg(1.3, 0.4)), t2_em(t2_tz(1.3))


def sg_t2_soliton(theta, eta, tau, which, c=0.3 + 1j*np.pi/2):
    """SG soliton e^{u/2} = (1+E)/(1-E), E = exp(sqrt2(x cosh th - t sinh th) + c), through the type-II defect (sigma = e^eta);
    which = 1: z1 = tanh tanh (cz2009b (4.3)), 2: z2 = 1/z1.  lam from the two lam_t-free conditions (cz (4.1)); returns
    (|e^lam e^-lam - 1|, residual of the lam_t condition, e^lam at t -> -+inf, deviation from the corrected (4.4), from (4.4))."""
    sig = np.exp(eta); gam = np.exp(2*tau) + np.exp(-2*tau)
    rho, rhot = np.tanh((theta - eta + tau)/2), np.tanh((theta - eta - tau)/2)
    z = rho*rhot if which == 1 else 1/(rho*rhot)
    a_, b_ = np.sqrt(2)*np.cosh(theta), -np.sqrt(2)*np.sinh(theta)

    def fld(t, zz):
        E = zz*np.exp(b_*t + c)
        return 2*np.log((1 + E)/(1 - E)), 2*a_*(E/(1 + E) + E/(1 - E)), 2*b_*(E/(1 + E) + E/(1 - E))

    def lam_pair(t):
        u, ux, ut = fld(t, 1.0); v, vx, vt = fld(t, z)
        p, q = (u + v)/2, (u - v)/2; qx, qt = (ux - vx)/2, (ut - vt)/2
        el = -(2/sig)*(sig**2*(qx + qt) + np.exp(p)*(qx - qt))/((np.exp(p) - np.exp(-p))*(np.exp(q) + np.exp(-q) + gam))
        eml = -(1/sig)*(sig**2*(qx + qt) + np.exp(-p)*(qx - qt))/(np.exp(p) - np.exp(-p))
        return el, eml
    ts = np.linspace(-3, 3, 9)/max(abs(np.sinh(theta)), 0.3); hh = 1e-5; cons = res3 = dev = dev_pr = 0
    for t in ts:
        el, eml = lam_pair(t); cons = max(cons, abs(el*eml - 1))
        lt = (np.log(lam_pair(t + hh)[0]) - np.log(lam_pair(t - hh)[0]))/(2*hh)
        u, ux, ut = fld(t, 1.0); v, vx, vt = fld(t, z)
        px_, vx_, c3 = t2_conditions(t2_sg(sig, tau), u, v, np.log(el), ut, vt, lt)
        res3 = max(res3, abs(px_ - ux), abs(vx_ - vx), abs(c3))
        E0 = np.exp(b_*t + c)
        if which == 1:
            corr = -(1/(np.sqrt(2)*np.cosh(tau)))*(1 + E0)*(1 + z*E0)/((1 - rho*E0)*(1 - rhot*E0))
            prnt = (1/(np.sqrt(2)*np.cosh(tau)))*(1 + E0)*(1 + z*E0)/((1 + rho*E0)*(1 + rhot*E0))
            dev = max(dev, abs(el - corr)); dev_pr = max(dev_pr, abs(el - prnt))
    big = 13/abs(b_)                                    # |E0| ~ e^(+-13): asymptotic sectors
    return cons, res3, (lam_pair(-big)[0], lam_pair(big)[0]), dev, dev_pr


ok_sg = True; worst3 = 0; info = []
v0 = 1/(np.sqrt(2)*np.cosh(0.5))
for th in (-0.3, 0.9, 1.9):
    for w in (1, 2):
        cons, res3, asym, dev, dev_pr = sg_t2_soliton(th, 0.7, 0.5, w)
        worst3 = max(worst3, res3, cons)
        sgn = -1 if w == 1 else 1
        ok_sg &= all(abs(x - sgn*v0) < 1e-3 for x in asym)
        if w == 1:
            ok_sg &= dev < 1e-10 and dev_pr > 0.1
report(f'type II, sine-Gordon (cz2009b (3.2), sigma = e^eta): E + D and P + Omega conserved [rel. res {em_sg:.0e}]; a soliton '
       f'passes with z1 = tanh((th-eta+tau)/2) tanh((th-eta-tau)/2) when e^lam -> -1/(sqrt2 cosh tau) and with z2 = 1/z1 '
       f'(coth coth) when e^lam -> +1/(sqrt2 cosh tau) [conditions {worst3:.0e}]; for z1, e^lam = -(1/sqrt2 cosh tau)'
       f'(1+E0)(1+z1 E0)/((1-rho E0)(1-rhot E0)), not cz2009b (4.4) as printed', ok_sg and em_sg < 1e-7 and worst3 < 1e-6)


def tz_t2_soliton(theta, eta, which, c=0.1 + 0.9j):
    """Tzitzeica soliton e^u = (1+E)^2/(1-4E+E^2), E = exp(sqrt6(x cosh th - t sinh th) + c), through the type-II defect,
    sigma = e^eta/sqrt2, delays cz2009b (4.9); e^lam = common root of the two lam_t-free conditions (cubics (4.7));
    returns (z, max root mismatch, residual of the lam_t condition, e^lam at t -> -+inf)."""
    sig = np.exp(eta)/np.sqrt(2); X = np.exp(-theta + eta); w6 = np.exp(1j*np.pi/6)
    z1 = (X + w6)*(X + 1/w6)/((X - w6)*(X - 1/w6))
    z2 = (X - w6)*(X + np.exp(-1j*np.pi/2))/((X + w6)*(X - np.exp(-1j*np.pi/2)))
    z = {1: z1, 2: z2, 3: np.conj(z2)}[which]
    a_, b_ = np.sqrt(6)*np.cosh(theta), -np.sqrt(6)*np.sinh(theta)

    def fld(t, zz):
        E = zz*np.exp(b_*t + c); Q = 1 - 4*E + E**2; du = E*(2/(1 + E) - (-4 + 2*E)/Q)
        return 2*np.log(1 + E) - np.log(Q), a_*du, b_*du

    def roots(t):
        u, ux, ut = fld(t, 1.0); v, vx, vt = fld(t, z)
        p, q = (u + v)/2, (u - v)/2; qx, qt = (ux - vx)/2, (ut - vt)/2; C = np.exp(q) + np.exp(-q)
        r1 = np.roots([-(sig/2)*C**2, 0, -(qx - qt), 2*sig])
        r2 = np.roots([np.exp(-p)*C/(2*sig), -(qx + qt), 0, -np.exp(2*p)/sig])
        dd = np.abs(r1[:, None] - r2[None, :])
        return r1, dd
    ts = np.linspace(-2.5, 2.5, 11)/max(abs(np.sinh(theta)), 0.3); hh = 1e-5; mism = res3 = 0
    for t in ts:
        r1, dd = roots(t); i, j = np.unravel_index(np.argmin(dd), dd.shape); L0 = r1[i]; mism = max(mism, dd[i, j])
        Lp = min(roots(t + hh)[0], key=lambda y: abs(y - L0)); Lm = min(roots(t - hh)[0], key=lambda y: abs(y - L0))
        lt = (np.log(Lp) - np.log(Lm))/(2*hh)
        u, ux, ut = fld(t, 1.0); v, vx, vt = fld(t, z)
        px_, vx_, c3 = t2_conditions(t2_tz(sig), u, v, np.log(L0), ut, vt, lt)
        res3 = max(res3, abs(px_ - ux), abs(vx_ - vx), abs(c3))
    big = 13/abs(b_)
    asym = []
    for t in (-big, big):
        r1, dd = roots(t); i, j = np.unravel_index(np.argmin(dd), dd.shape); asym.append(r1[i])
    return z, mism, res3, asym


ok_tz = True; worst_tz = 0; zs = {}
for th in (0.4, 1.3):
    for w in (1, 2, 3):
        z, mism, res3, asym = tz_t2_soliton(th, 0.6, w)
        zs[(th, w)] = z; worst_tz = max(worst_tz, mism, res3)
        target = {1: 1, 2: np.exp(-2j*np.pi/3), 3: np.exp(2j*np.pi/3)}[w]
        ok_tz &= all(abs(x - target) < 1e-3 for x in asym)
    u = th - 0.6
    ok_tz &= np.isclose(zs[(th, 1)], 1/(np.tanh(u/2 - 1j*np.pi/12)*np.tanh(u/2 + 1j*np.pi/12)))       # (4.10)
    ok_tz &= np.isclose(zs[(th, 3)], 1/np.tanh(u/2 + 1j*np.pi/4)*np.tanh(u/2 - 1j*np.pi/12))          # (4.11) is z3 of (4.9)
    ok_tz &= not np.isclose(zs[(th, 2)], 1/np.tanh(u/2 + 1j*np.pi/4)*np.tanh(u/2 - 1j*np.pi/12))
    ok_tz &= np.isclose(zs[(th, 1)]*zs[(th, 2)]*zs[(th, 3)], 1)
report(f'type II, Tzitzeica a_2^(2) (cz2009b (3.7), e^eta = sqrt2 sigma): E + D and P + Omega conserved [rel. res {em_tz:.0e}]; '
       f'a soliton passes with each of the delays z1, z2, z3 = conj z2 of (4.9), z1 z2 z3 = 1, with e^lam -> 1, e^(-2 pi i/3), '
       f'e^(2 pi i/3) respectively [conditions {worst_tz:.0e}]; (4.10) = z1, but (4.11) as printed equals z3 = conj z2, not z2', ok_tz and em_tz < 1e-7 and worst_tz < 1e-6)


# ---------------------------------------------------------------- 11. prp-spin-two: the cubic obstruction for all r, and B_{+-2} for a_2
# BCZ (5.5) with B = 1 - D = E (their (5.6)) is solved by A(u,u,u) = (2/3) sum_k (mu_k.u)^3, mu_k the weights of the defining
# representation (H_k = e_kk); their obstruction (5.11), A(Dv)^3 + A(D^T v)^3, vanishes identically although each term does not.
ok55 = True; obs = 0; single = np.inf
for r in range(2, 9):
    d = an_data(r); mu = d['U']; al_r = d['alpha']; D, E = dmatrix(r)
    for w in rng.normal(size=(3, r)):
        ok55 &= all(np.isclose(2*np.sum((mu @ w)**2*(mu @ a)), -2*(a @ E @ w)*(a @ w)) for a in al_r)
    v = rng.normal(size=r); cubD = (2/3)*np.sum((mu @ (D @ v))**3); cubDT = (2/3)*np.sum((mu @ (D.T @ v))**3)
    obs = max(obs, abs(cubD + cubDT)/abs(cubD)); single = min(single, abs(cubD))
cot_ok = True
for h in range(3, 12):
    for k1, k2 in itertools.product(range(1, h), repeat=2):
        k3 = (-k1 - k2) % h
        if k3 == 0:
            continue
        c1, c2, c3 = [1/np.tan(np.pi*k/h) for k in (k1, k2, k3)]
        d1, d2, d3 = 1 + 1j*c1, 1 + 1j*c2, 1 + 1j*c3                  # D on the Fourier modes
        cot_ok &= np.isclose(c1*c2 + c2*c3 + c3*c1, 1) and np.isclose(d1*d2*d3 + np.conj(d1*d2*d3), 0)
report(f'BCZ (5.5) with B = 1 - D is solved by A = (2/3) sum_k (mu_k.u)^3 and the obstruction (5.11) A(Dv)^3 + A(D^T v)^3 '
       f'vanishes for r = 2..8 [relative {obs:.0e}; single terms >= {single:.2f}]; on Fourier modes D = 1 + i cot(pi k/h) and '
       f'Re(d1 d2 d3) = 1 - (c1c2 + c2c3 + c3c1) = 0 for k1 + k2 + k3 = 0 mod h (h = 3..11)', ok55 and obs < 1e-12 and cot_ok)


def construct_B2(G):
    """Local F(q, qdot) with dF/dt = G (G affine in qddot): F1 = int_0^1 A(q, s qdot).qdot ds with A = dG/dqddot (homotopy),
    then F0(q) from the remainder sum_a R_a(q) qdot_a by matching each exponential (R_a = c_a kappa e^{c.q})."""
    s_ = sp.Symbol('s')
    A_ = [sp.diff(G, q2[a]) for a in range(4)]
    F1 = sp.integrate(sp.expand(sum(A_[a].subs({q1[b]: s_*q1[b] for b in range(4)}, simultaneous=True)*q1[a]
                                    for a in range(4))), (s_, 0, 1))
    R_ = sp.expand(G - Dt(F1))
    Ra = [sp.diff(R_, q1[a]).subs({v: 0 for v in q1}) for a in range(4)]
    terms = {}
    for a in range(4):
        for term in sp.Add.make_args(sp.powsimp(sp.expand(Ra[a]), combine='exp')):
            coeff, rest = term.as_independent(*q0)
            lin = sp.expand(rest.args[0]) if isinstance(rest, sp.exp) else None
            terms.setdefault(lin, [0]*4)[a] += coeff
    F0 = 0; mism = 0
    for lin, co in terms.items():
        cvec = [sp.diff(lin, q0[a]) for a in range(4)]
        a0 = next(a for a in range(4) if cvec[a] != 0)
        kap = co[a0]/cvec[a0]
        mism += sum(abs(complex(sp.N(co[a] - kap*cvec[a]))) for a in range(4))
        F0 += kap*sp.exp(lin)
    return F1, F0, mism


ok_b2 = True; worst_b2 = 0; shapes = []
allv = [v for qq in Q for v in qq]
for eps in (-1, 1):
    G2 = boundary_G(eps, 1, 2)
    F1, F0, mism = construct_B2(G2)
    f = sp.lambdify(allv, [G2, Dt(F1 + F0)], 'numpy', cse=True)
    for _ in range(4):
        vals = rng.normal(size=len(allv))*0.5
        g_, df_ = f(*vals); worst_b2 = max(worst_b2, abs(g_ - df_)/(1 + abs(g_)))
    deg = max(sum(mo) for mo in sp.Poly(sp.expand(F1), *q1).monoms())
    shapes.append((len(sp.Add.make_args(sp.expand(F1))), deg, len(sp.Add.make_args(sp.powsimp(sp.expand(F0), combine='exp')))))
    ok_b2 &= mism < 1e-10 and deg == 1
report(f'explicit defect corrections B_(+2) (eps = -1) and B_(-2) (eps = +1) of a_2, type I_-, built by the homotopy formula: '
       f'B = (terms linear in phi_t, psi_t) + (exponentials) with dB/dt equal to the defect flux [rel. res {worst_b2:.0e}; '
       f'(velocity terms, degree, exponentials) = {shapes}]', ok_b2 and worst_b2 < 1e-12)
