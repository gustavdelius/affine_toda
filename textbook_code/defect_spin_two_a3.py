"""Spin +-2 charges of a_3^(1) with a defect, reproducing prp-spin-two item (3) for a_3 (Part VIII, ch. 38).

Companion to defect_classical_an.py's checks 5 and 11, which establish the a_2 case exactly (sympy) and the
cubic obstruction (5.11) for r = 2..8 numerically.  Those checks hard-code the a_2 root system and solve a
small linear system for the bulk spin-2 tensor A_abc.  Here that construction is generalized to any r, using
the closed-form A_abc = 2 sum_k mu_ka mu_kb mu_kc of check 11 (mu = weights of the defining representation,
an_data(r)['U']), and run for r = 3 (h = 4).  BCZ = Bowcock, Corrigan, Zambon hep-th/0401020.

Why this script is separate: the symbolic construction of B_{+-2} (check 11's homotopy formula, with sigma
symbolic) is slow at r = 3 (the non-integrable-D control alone takes about 2 minutes), so it is kept out of
defect_classical_an.py's main run and given its own script, as flagged in briefs/part8_defects_status.md.

Checks (book ids they support):
 1. Exact (sympy) Gram-Schmidt construction of the a_r^(1) root data (simple roots e_j - e_{j+1} in R^h,
    Gram-Schmidt on alpha_1..alpha_r, same order and sign convention as defect_common.an_data): agrees with
    an_data(r)/dmatrix(r) for r = 1..4.  This replaces a2's hand-written root vectors with a general-r
    construction.                                                                      [eq-d-matrix]
 2. The closed-form tensor A_abc = 2 sum_k mu_ka mu_kb mu_kc of check 11 solves, EXACTLY (sympy, symbolic
    vector), the defect_classical_an.py check-5 matching equation 3 A(w,w,alpha) = -2(alpha.B.w)(alpha.w)
    (for every simple and affine root alpha) with B = 3(1 - D) -- the same relation as check 11's direct
    form A(w,w,alpha) = -2(alpha.E.w)(alpha.w) with B = E, scaled by the factor of 3 that check 5's
    EOM-derived matching equation carries relative to BCZ's direct form.  Checked for a_2 (cross-check
    against check 5's solve) and a_3.
 3. BCZ (5.11)'s cubic obstruction A(Dv)^3 + A(D^T v)^3 vanishes identically for a_3, EXACTLY (sympy,
    symbolic vector; check 11 showed this numerically for r = 2..8).                    [prp-spin-two]
 4. The a_3 bulk spin-2 current is conserved on shell, EXACTLY (symbolic on-shell identity, not a numerical
    residual): d_{-eps} T = d_eps Theta with T = A(u,u,u) + w.Bu, Theta = sum_j (g_j.u) e^{alpha_j.phi},
    g_j = -B alpha_j, u = d_eps phi, w = d_eps^2 phi.
 5. Main result: for a_3^(1), both types (I_-, I_+) and both chiralities (eps = +-1, i.e. spin -+2), with
    sigma symbolic, the defect flux of the spin +-2 charge is exactly a total time derivative of a local
    function B_{+-2}(phi(0,t), psi(0,t), their t-derivatives): defect flux - dB_{+-2}/dt simplifies to 0
    EXACTLY (sympy; same homotopy construction as check 11, generalized to r = 3).  As for a_2, B_{+2}
    (eps = -1) contains the velocities only multiplied by sigma^{-1} times an exponential, and B_{-2}
    (eps = +1) only by sigma times an exponential.                                      [prp-spin-two]
 6. Controls: the spin +-1 defect contributions are likewise exact total time derivatives (same homotopy
    construction, spin = 1) for both types; a non-integrable D (one off-diagonal entry of E changed) makes
    the construction fail (nonzero mismatch and nonzero residual) for a_3.               [prp-spin-two]
Runtime: about 4-5 minutes (dominated by the sigma-symbolic homotopy construction for spin 2 and, most of
all, by the non-integrable-D control, whose messier exponential structure is slower to simplify).
"""
import itertools
import time
import numpy as np
import sympy as sp
from defect_common import an_data, dmatrix

T0 = time.time()


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + f'[{time.time() - T0:6.1f}s] ' + name)


# ---------------------------------------------------------------- 1. exact root data for general r
def exact_an_root_data(r):
    """Exact (sympy) version of defect_common.an_data(r)/dmatrix(r): simple roots alpha_h[j] = e_j - e_{j+1}
    (j = 0..r, mod h = r+1) in R^h, Gram-Schmidt on alpha_h[1]..alpha_h[r] in that order (reproducing
    numpy.linalg.qr's basis and sign, since standard Gram-Schmidt already gives a positive projection of
    each new basis vector onto the vector it was built from).  Returns (alpha, D, E, U): alpha a list of
    r + 1 exact r-vectors (rows alpha_0..alpha_r), D, E the r x r matrices of eq-d-matrix, U the h x r matrix
    of ON-basis coordinates (rows = weights mu_k of the defining representation, as in an_data(r)['U'])."""
    h = r + 1
    Ev = [sp.zeros(h, 1) for _ in range(h)]
    for i in range(h):
        Ev[i][i] = 1
    alpha_h = [Ev[j] - Ev[(j + 1) % h] for j in range(h)]
    simple = [alpha_h[j] for j in range(1, h)]
    basis = []
    for v in simple:
        w = v
        for b in basis:
            w = w - (v.T * b)[0] * b
        w = sp.simplify(w / sp.sqrt((w.T * w)[0]))
        basis.append(w)
    U = sp.Matrix.hstack(*basis)
    alpha = [sp.simplify(sp.Matrix([(alpha_h[j].T * basis[i])[0] for i in range(r)])) for j in range(h)]
    simpleM = sp.Matrix.hstack(*alpha[1:]).T
    lam_s = sp.simplify(simpleM.inv()).T
    lam = [sp.zeros(r, 1)] + [lam_s.row(i).T for i in range(r)] + [sp.zeros(r, 1)]
    D = sp.zeros(r, r)
    for i in range(1, r + 1):
        D += 2 * (lam[i] - lam[i + 1]) * lam[i].T
    D = sp.simplify(D)
    E = sp.eye(r) - D
    return alpha, D, E, U


ok_roots = True
for rr in (1, 2, 3, 4):
    al_, D_, E_, U_ = exact_an_root_data(rr)
    Dref, _ = dmatrix(rr)
    ok_roots &= np.allclose(np.array(D_.evalf(), dtype=float), Dref, atol=1e-10)
    ok_roots &= np.allclose(np.array([[float(x) for x in a] for a in al_]), an_data(rr)['alpha'], atol=1e-10)
report('exact Gram-Schmidt root data (general r, closed-form coordinates) agrees with '
       'defect_common.an_data/dmatrix for r = 1..4', ok_roots)

r = 3
al, Dm, Em, U = exact_an_root_data(r)
h = r + 1

# ---------------------------------------------------------------- 2. A_abc formula solves the matching equation
Avals = {idx: sp.simplify(2 * sum(U[k, idx[0]] * U[k, idx[1]] * U[k, idx[2]] for k in range(h)))
         for idx in itertools.combinations_with_replacement(range(r), 3)}


def Afun(a, b, c):
    return Avals[tuple(sorted((a, b, c)))]


Bnum = 3 * Em                     # matches check 5's "3 A(u,u,alpha) = -2(alpha.B.u)(alpha.u)" with A_abc = 2 sum mu mu mu

ok_match = True
for rr, (al_, _, E_, U_) in ((2, exact_an_root_data(2)), (3, (al, Dm, Em, U))):
    Av_ = {idx: sp.simplify(2 * sum(U_[k, idx[0]] * U_[k, idx[1]] * U_[k, idx[2]] for k in range(rr + 1)))
           for idx in itertools.combinations_with_replacement(range(rr), 3)}
    wv = sp.Matrix(sp.symbols(f'w0:{rr}'))
    Bn_ = 3 * E_
    for a in al_:
        lhs = 3 * sum(Av_[tuple(sorted((x, y, z)))] * wv[x] * wv[y] * a[z]
                      for x in range(rr) for y in range(rr) for z in range(rr))
        rhs = -2 * (a.T * Bn_ * wv)[0] * (a.T * wv)[0]
        ok_match &= sp.expand(lhs - rhs) == 0
report('A_abc = 2 sum_k mu_ka mu_kb mu_kc solves 3 A(w,w,alpha) = -2(alpha.B.w)(alpha.w) with B = 3(1 - D) '
       'EXACTLY (sympy, symbolic w), for a_2 and a_3 (the same relation as check 11, with B = E, scaled by 3)',
       ok_match)

# ---------------------------------------------------------------- 3. cubic obstruction (5.11) vanishes for a_3
v = sp.Matrix(sp.symbols(f'v0:{r}'))


def cubic(vec):
    return sum(Afun(x, y, z) * vec[x] * vec[y] * vec[z] for x in range(r) for y in range(r) for z in range(r))


obstruction = sp.expand(cubic(Dm * v) + cubic(Dm.T * v))
report('a_3: BCZ (5.11) cubic obstruction A(Dv)^3 + A(D^T v)^3 vanishes identically EXACTLY (sympy, symbolic v)',
       obstruction == 0)

# ---------------------------------------------------------------- 4. a_3 bulk spin-2 current conserved on shell
fv = sp.Matrix(sp.symbols(f'x0:{r}'))
uu = sp.Matrix(sp.symbols(f'u0:{r}'))
ww = sp.Matrix(sp.symbols(f'w0:{r}'))
w3 = sp.Matrix(sp.symbols(f'w3_0:{r}'))
vv = sp.Matrix(sp.symbols(f'vv0:{r}'))
Vb = sum(sp.exp((aj.T * fv)[0]) - 1 for aj in al)
gVb = sp.Matrix([sp.diff(Vb, c) for c in fv])
HVb = sp.hessian(Vb, fv)
Tb = sum(Afun(a, b, c) * uu[a] * uu[b] * uu[c] for a in range(r) for b in range(r) for c in range(r)) + (ww.T * Bnum * uu)[0]
Thb = sum(((-Bnum * aj).T * uu)[0] * sp.exp((aj.T * fv)[0]) for aj in al)
d_same = lambda e: sum(sp.diff(e, fv[a]) * uu[a] + sp.diff(e, uu[a]) * ww[a] + sp.diff(e, ww[a]) * w3[a] for a in range(r))
d_other = lambda e: sum(sp.diff(e, fv[a]) * vv[a] - sp.diff(e, uu[a]) * gVb[a] - sp.diff(e, ww[a]) * (HVb * uu)[a] for a in range(r))
bulk_expr = sp.expand(d_other(Tb) - d_same(Thb))
report('a_3 spin-2 bulk current conserved on shell: d_(-eps)T = d_eps Theta vanishes EXACTLY (sympy identity, '
       'not a numerical residual)', bulk_expr == 0)

# ---------------------------------------------------------------- 5-6. B_{+-2} for a_3, sigma symbolic, and controls
nf = 2 * r


def make_level(tag):
    return tuple(sp.symbols([f'{tag}f{i}' for i in range(r)] + [f'{tag}g{i}' for i in range(r)]))


q0, q1, q2, q3, q4 = [make_level(t) for t in ('', 'd1_', 'd2_', 'd3_', 'd4_')]
Q = [q0, q1, q2, q3, q4]
SIG = sp.Rational(7, 10)
sigma_s = sp.Symbol('sigma', positive=True)


def Dt(e):
    return sum(Q[k + 1][a] * sp.diff(e, Q[k][a]) for k in range(4) for a in range(nf))


def boundary_G(eps, typ, spin, Dbase=None, sig=None):
    """T(phi) - T(psi) - [Theta(phi) - Theta(psi)], the defect flux of the spin eps*(-2) [spin=2] or
    eps*(-1) [spin=1] charge, as a function of phi(0,t), psi(0,t) and their t-derivatives (q0, q1, q2).
    Same construction as defect_classical_an.py's boundary_G, generalized from r = 2 to general r."""
    Dm_ = Dm if Dbase is None else Dbase
    Em_ = sp.eye(r) - Dm_
    D_, E_ = (Dm_, Em_) if typ == 1 else (-Dm_, -Em_)
    sg = SIG if sig is None else sig
    s = sg if typ == 1 else 1 / sg
    phi, psi = sp.Matrix(q0[:r]), sp.Matrix(q0[r:])
    phid, psid = sp.Matrix(q1[:r]), sp.Matrix(q1[r:])
    phidd, psidd = sp.Matrix(q2[:r]), sp.Matrix(q2[r:])
    B = sum(s * (sp.exp((aj.T * (Dm_.T * phi + Dm_ * psi))[0] / 2) - 1) +
            (sp.exp((aj.T * Dm_ * (phi - psi))[0] / 2) - 1) / s for aj in al)
    gp = sp.Matrix([sp.diff(B, vx) for vx in phi])
    gq = sp.Matrix([sp.diff(B, vx) for vx in psi])
    Fp = E_ * phid + D_ * psid - gp
    Fq = D_.T * phid - E_ * psid + gq
    gV = lambda f: sum((aj * sp.exp((aj.T * f)[0]) for aj in al), sp.zeros(r, 1))
    V = lambda f: sum(sp.exp((aj.T * f)[0]) - 1 for aj in al)
    out = 0
    for sgn, f, fd, fdd, F in ((1, phi, phid, phidd, Fp), (-1, psi, psid, psidd, Fq)):
        u = fd + eps * F
        if spin == 1:
            TT = (u.T * u)[0] / 2 - V(f)
        else:
            Fd = sp.Matrix([Dt(c) for c in F])
            w = fdd + 2 * eps * Fd + fdd + gV(f)
            T3 = sum(Afun(a, b, c) * u[a] * u[b] * u[c] for a in range(r) for b in range(r) for c in range(r))
            T3 += (w.T * Bnum * u)[0]
            Th = sum(((-Bnum * aj).T * u)[0] * sp.exp((aj.T * f)[0]) for aj in al)
            TT = T3 + Th
        out += sgn * TT
    return out


def construct_B2(G):
    """Local F(q, qdot) with dF/dt = G (G affine in qddot): F1 = int_0^1 A(q, s qdot).qdot ds with
    A = dG/dqddot (homotopy formula), then F0(q) from the remainder sum_a R_a(q) qdot_a by matching each
    exponential (R_a = c_a kappa e^{c.q}).  Same construction as defect_classical_an.py check 11."""
    s_ = sp.Symbol('s')
    A_ = [sp.diff(G, q2[a]) for a in range(nf)]
    F1 = sp.integrate(sp.expand(sum(A_[a].subs({q1[b]: s_ * q1[b] for b in range(nf)}, simultaneous=True) * q1[a]
                                     for a in range(nf))), (s_, 0, 1))
    R_ = sp.expand(G - Dt(F1))
    Ra = [sp.diff(R_, q1[a]).subs({vx: 0 for vx in q1}) for a in range(nf)]
    terms = {}
    for a in range(nf):
        for term in sp.Add.make_args(sp.powsimp(sp.expand(Ra[a]), combine='exp')):
            coeff, rest = term.as_independent(*q0)
            lin = sp.expand(rest.args[0]) if isinstance(rest, sp.exp) else None
            terms.setdefault(lin, [0] * nf)[a] += coeff
    F0 = 0
    mism = 0
    for lin, co in terms.items():
        cvec = [sp.diff(lin, q0[a]) for a in range(nf)]
        a0 = next(a for a in range(nf) if cvec[a] != 0)
        kap = co[a0] / cvec[a0]
        mism += sum(1 for a in range(nf) if sp.simplify(co[a] - kap * cvec[a]) != 0)
        F0 += kap * sp.exp(lin)
    return F1, F0, mism


res2 = {}
shape_ok = True
for typ in (1, -1):
    for eps in (-1, 1):
        G2 = boundary_G(eps, typ, 2, sig=sigma_s)
        F1, F0, mism = construct_B2(G2)
        resid = sp.expand(sp.powsimp(sp.expand(G2 - Dt(F1 + F0)), combine='exp'))
        res2[(typ, eps)] = (mism == 0 and resid == 0)
        if typ == 1:   # velocity-structure check (type I_-): spin+2 (eps=-1) only 1/sigma, spin-2 (eps=1) only sigma

            def sigma_power(term):
                p = 0
                for factor in sp.Mul.make_args(term):
                    if factor.has(sigma_s):
                        base, ex = factor.as_base_exp()
                        if base == sigma_s:
                            p += ex
                return p

            terms = [t for t in sp.Add.make_args(sp.expand(F1)) if t != 0]
            want = -1 if eps == -1 else 1
            shape_ok &= all(sigma_power(t) == want for t in terms)
report(f'a_3 spin +-2 charges, sigma symbolic: defect flux - dB_(+-2)/dt = 0 EXACTLY, for type I_- AND type '
       f'I_+ (results {res2})', all(res2.values()))
report('a_3, type I_-: B_(+2) (eps = -1) contains the velocities only multiplied by sigma^-1 times an '
       'exponential, B_(-2) (eps = +1) only by sigma times an exponential (as for a_2)', shape_ok)

res1 = {}
for typ in (1, -1):
    for eps in (-1, 1):
        G1 = boundary_G(eps, typ, 1, sig=sigma_s)
        F1, F0, mism = construct_B2(G1)
        resid = sp.expand(sp.powsimp(sp.expand(G1 - Dt(F1 + F0)), combine='exp'))
        res1[(typ, eps)] = (mism == 0 and resid == 0)
report(f'control: a_3 spin +-1 (E +- P) defect contributions are likewise exact total time derivatives, '
       f'sigma symbolic, both types (results {res1})', all(res1.values()))

Ebad = sp.Matrix(Em)
Ebad[0, 1] = sp.Rational(3, 10)
Ebad[1, 0] = -sp.Rational(3, 10)
Dbad = sp.eye(r) - Ebad
Gbad = boundary_G(1, 1, 2, Dbase=Dbad, sig=sigma_s)
F1b, F0b, mismb = construct_B2(Gbad)
residb = sp.expand(sp.powsimp(sp.expand(Gbad - Dt(F1b + F0b)), combine='exp'))
report(f'control: non-integrable D (one off-diagonal entry of E changed) spoils the a_3 spin-2 construction '
       f'[mismatch count {mismb}, residual {"= 0" if residb == 0 else "!= 0"}]',
       mismb != 0 or residb != 0)

report(f'DONE in {time.time() - T0:.1f}s total', True)
