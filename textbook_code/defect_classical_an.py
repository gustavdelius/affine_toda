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
Runtime: 35-100 seconds.
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
