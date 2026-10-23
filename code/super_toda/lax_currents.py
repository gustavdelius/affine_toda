"""Zero curvature and conserved supercurrents of the super sinh-Gordon theory
(sec-super-lax, sec-super-shg).  Conventions of superspace.py.

Checks:
 1. A realization of the C(2)^(2) relations [h, e_1] = c e_1, [h, e_0] = -c e_0,
    {e_i, f_j} = delta_ij alpha_j.H in the loop algebra of sl(2|1) (3x3
    supermatrices of grading (0,0,1)), with c = alpha^2 = 1/2:
    h = c diag(1,-1,0), e_1 = E13 + E32, f_1 = c(E31 - E23),
    e_0 = -c(E23 + E31), f_0 = E13 - E32  (eq-ssg-realization).
 2. In the loop algebra of osp(1|2) the same assignment fails: there the
    generators of the two odd simple roots are built from the same pair of odd
    generators F_+-, and {F_-, F_-} = 2F_-^2 is nonzero.
 3. The superfield equation (eq-super-eom) follows from the superspace action:
    its theta components are sigma = ..., the fermion equations and the boson
    equation, and they agree with the Euler-Lagrange equations of the component
    Lagrangian eq-ssg-components (here derived again from the superfield).
 4. Zero curvature (eq-super-lax, eq-super-zero-curvature): with
    A_+ = beta D_+Phi.H + lambda(e_1 + e_0), A_- = lambda^{-1}(eta_1 e^{beta alpha.Phi} f_1
    + eta_0 e^{-beta alpha.Phi} f_0), the curvature D_+A_- + D_-A_+ + {A_+, A_-}
    equals (i beta/|alpha|) x (field equation) x h, for every lambda, when
    eta_0 = eta_1 = -i mu.
 5. Conserved supercurrents for L = (i/2) D_+Phi D_-Phi - m cosh Phi (the super
    sinh-Gordon Lagrangian (1.5) of Evans-Madsen 1999), eq-ssg-currents:
    D_-J = D_+Jbar on shell with
    J = T = D_+Phi D_+^2Phi - D_+^3Phi, Jbar = -i m D_+ e^{-Phi}  (spin 3/2), and
    J = T D_+T, Jbar = i m (D_+Phi (D_+^2Phi)^2 - D_+^2Phi D_+^3Phi) e^{-Phi}  (spin 7/2).
    The component equations (5.48)-(5.51) of Evans-Madsen are those of D_-J = D_+Jbar.
    Their printed (5.44) has D_+J = D_-Jbar, and (5.45), (5.46) have Jbar with the
    factor m/2 in place of -m and m; the script confirms that the printed forms
    do not hold (sec-super-notes, sec-app-corrections).
"""
import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from superspace import G, Jets, SMat, unit, sexp, berezin, d_theta, I, TP, TM   # noqa: E402

FAILS = []


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        FAILS.append(name)


def gb(X, Y, odd_pair):
    """[X, Y} of envelope elements: anticommutator if both are odd overall."""
    return X * Y + Y * X if odd_pair else X * Y - Y * X


# ============================================================ 1. realization
c = sp.Rational(1, 2)
par = [0, 0, 1]
E = lambda r, s: unit(par, r, s)
h = (E(0, 0).scale(c) + E(1, 1).scale(-c))
e1 = E(0, 2) + E(2, 1)
f1 = E(2, 0).scale(c) + E(1, 2).scale(-c)
e0 = E(1, 2).scale(-c) + E(2, 0).scale(-c)
f0 = E(0, 2) + E(2, 1).scale(-1)
ok = all([
    (gb(h, e1, False) - e1.scale(c)).is_zero(),
    (gb(h, e0, False) + e0.scale(c)).is_zero(),
    (gb(h, f1, False) + f1.scale(c)).is_zero(),
    (gb(h, f0, False) - f0.scale(c)).is_zero(),
    (gb(e1, f1, True) - h).is_zero(),
    (gb(e0, f0, True) + h).is_zero(),
    gb(e1, f0, True).is_zero(),
    gb(e0, f1, True).is_zero(),
])
str_h = c - c - 0          # supertrace of h: rows (0,0) even, (1,1) even, (2,2) odd
report('sl(2|1) realization of the C(2)^(2) relations with alpha^2 = 1/2 (eq-ssg-realization); '
       'h has zero supertrace', ok and str_h == 0)

# ==================================================== 2. the osp(1|2) obstruction
par2 = [0, 1, 1]
U = lambda r, s: unit(par2, r, s)
Fp = U(0, 1) + U(2, 0).scale(-1)
Fm = U(0, 2) + U(1, 0)
Hs = gb(Fp, Fm, True)
report('osp(1|2): {F+, F-} is diagonal, but {F-, F-} = 2F-^2 and {F+, F+} are nonzero, so '
       '{e_0, f_1} cannot vanish when e_0 and f_1 are both built from F-',
       all(k[0] == k[1] for k in Hs.e) and not gb(Fm, Fm, True).is_zero()
       and not gb(Fp, Fp, True).is_zero())

# ========================================== 3. field equations and components
beta, mu, lam = sp.symbols('beta mu lambda', positive=True)
J = Jets(1, 'q')
Phi = J.superfield(0)
a = sp.sqrt(c)                 # alpha.Phi = sqrt(c) Phi for the canonical scalar superfield
L = (J.D(Phi, +1) * J.D(Phi, -1)).map(lambda x: I / 2 * x) + \
    (sexp(Phi.map(lambda x: beta * a * x)) + sexp(Phi.map(lambda x: -beta * a * x))).map(
        lambda x: mu / beta**2 * x)
# superfield equation (Evans-Madsen (3.1) with alpha_j -> beta alpha_j, mu_j -> mu/beta^2):
#   i D_+D_-Phi = (mu/beta) a (e^{beta a Phi} - e^{-beta a Phi})
Eq = (J.D(J.D(Phi, -1), +1)).map(lambda x: I * x) - \
    (sexp(Phi.map(lambda x: beta * a * x)) - sexp(Phi.map(lambda x: -beta * a * x))).map(
        lambda x: mu / beta * a * x)
E0, Ep, Em, E2 = (J.component(Eq, w) for w in ('0', '+', '-', '2'))
phi, sig = J.phi(0), J.sig(0)
sig_sol = sp.solve(E0.body(), sig)[0]
# Euler-Lagrange equations of the component Lagrangian
Lc = berezin(L)
dL_dsig = sp.diff(Lc.body(), sig)
report('superfield equation, theta^0 component: sigma = (mu/beta) a (e^{beta a phi} - e^{-beta a phi}), '
       'the same as the Euler-Lagrange equation of sigma',
       sp.simplify(sig_sol - sp.solve(dL_dsig, sig)[0]) == 0)


def solve_linear_gen(Eodd, target):
    """Solve an odd equation, linear in the generator `target`, for that generator."""
    lab = list(target.t)[0][0]
    coef = Eodd.coeff((lab,))
    rest = Eodd - target.map(lambda x: coef * x)
    return rest.map(lambda x: -x / coef)


# theta^+ component contains d_+psi_-, theta^- component contains d_-psi_+
psim10 = J.psi('m', 0, 1, 0)
psip01 = J.psi('p', 0, 0, 1)
rule_psim = solve_linear_gen(Ep.map(lambda x: x.subs(sig, sig_sol)), psim10)
rule_psip = solve_linear_gen(Em.map(lambda x: x.subs(sig, sig_sol)), psip01)
phi11 = J.phi(0, 1, 1)
E2s = E2.map(lambda x: x.subs(sig, sig_sol))
k = sp.diff(E2s.body(), phi11)
rule_phi = (E2s - G.scalar(k * phi11)).map(lambda x: sp.expand(-x / k))
# Euler-Lagrange of the component Lagrangian after eliminating sigma
Lc_s = Lc.map(lambda x: sp.expand(x.subs(sig, sig_sol)))
ppl, pmi = J.psi('p', 0), J.psi('m', 0)
yuk = Lc_s.coeff((list(ppl.t)[0][0], list(pmi.t)[0][0]))
Ubody = -Lc_s.body().subs({J.phi(0, 1, 0): 0, J.phi(0, 0, 1): 0})
EL_phi = G.scalar(-sp.diff(Ubody, phi)) + (ppl * pmi).map(lambda x: sp.diff(yuk, phi) * x)
report('superfield equation, theta^+ theta^- component: d_+d_-phi equals the Euler-Lagrange '
       'equation of the component Lagrangian (Evans-Madsen (3.5) in these conventions)',
       (rule_phi - EL_phi).is_zero())
# fermion Euler-Lagrange: (i/2) psi+ d_- psi+ + yuk psi+ psi- -> i d_-psi+ + yuk psi- = 0
report('superfield equation, theta^- component: d_-psi_+ = i yuk psi_-, the Euler-Lagrange equation '
       'of the component Lagrangian',
       (rule_psip - pmi.map(lambda x: I * yuk * x)).is_zero())
report('superfield equation, theta^+ component: d_+psi_- = -i yuk psi_+',
       (rule_psim - ppl.map(lambda x: -I * yuk * x)).is_zero())


# ------------------------------------------------- on-shell reduction in jets
class OnShell:
    """Reduce jet expressions to the free data phi_{k,0}, phi_{0,l}, psi+_{k,0},
    psi-_{0,l} using sigma = sig_sol and the three rules above."""

    def __init__(self, J, sig_sol, rule_phi, rule_psip, rule_psim):
        self.J = J
        self.sig_sol, self.rule_phi = sig_sol, rule_phi
        self.rule_psip, self.rule_psim = rule_psip, rule_psim
        self.memo = {}

    def _jet(self, name):
        body = name[len(self.J.name):]
        if body.startswith('phi') or body.startswith('sig'):
            kind = body[:3]
            _, k, l = body[3:].split('_')
            return kind, int(k), int(l)
        s = body[3]
        _, k, l = body[4:].split('_')
        return 'psi' + s, int(k), int(l)

    def _dpow(self, x, k, l):
        for _ in range(k):
            x = self.J.d(x, +1)
        for _ in range(l):
            x = self.J.d(x, -1)
        return x

    def jet(self, name):
        """Reduced Grassmann element for a jet variable."""
        if name in self.memo:
            return self.memo[name]
        kind, k, l = self._jet(name)
        if kind == 'phi' and (k == 0 or l == 0):
            out = G.scalar(sp.Symbol(name))
        elif kind == 'phi':
            out = self.reduce(self._dpow(self.rule_phi, k - 1, l - 1))
        elif kind == 'sig':
            out = self.reduce(self._dpow(G.scalar(self.sig_sol), k, l))
        elif kind == 'psip' and l == 0:
            out = G({(name,): sp.Integer(1)})
        elif kind == 'psip':
            out = self.reduce(self._dpow(self.rule_psip, k, l - 1))
        elif kind == 'psim' and k == 0:
            out = G({(name,): sp.Integer(1)})
        else:
            out = self.reduce(self._dpow(self.rule_psim, k - 1, l))
        self.memo[name] = out
        return out

    def _free_symbol(self, s):
        if not s.name.startswith(self.J.name):
            return True
        kind, k, l = self._jet(s.name)
        return kind == 'phi' and (k == 0 or l == 0)

    def reduce_coef(self, c):
        c = sp.expand(c)
        bad = [s for s in c.free_symbols if not self._free_symbol(s)]
        if not bad:
            return G.scalar(c)
        poly = sp.Poly(c, *bad)
        out = G()
        for mon, cf in poly.terms():
            term = G.scalar(cf)
            for s, p in zip(bad, mon):
                for _ in range(p):
                    term = term * self.jet(s.name)
            out = out + term
        return out

    def reduce(self, x):
        out = G()
        for mono, cf in x.t.items():
            term = self.reduce_coef(cf)
            for lab in mono:
                if lab in ('thp', 'thm'):
                    term = term * G({(lab,): sp.Integer(1)})
                else:
                    term = term * self.jet(lab)
            out = out + term
        return out.map(lambda q: sp.expand(q))


# ============================================================ 4. zero curvature
eta = -I * mu
Hc = h.scale(1 / a)                              # H along alpha: alpha.H = h, Phi.H = Phi h / a
Ap = Hc.scale(J.D(Phi, +1).map(lambda x: beta * x)) + (e1 + e0).scale(lam)
Am = f1.scale(sexp(Phi.map(lambda x: beta * a * x)).map(lambda x: eta / lam * x)) + \
    f0.scale(sexp(Phi.map(lambda x: -beta * a * x)).map(lambda x: eta / lam * x))
DpAm = Am.map(lambda x: J.D(x, +1))
DmAp = Ap.map(lambda x: J.D(x, -1))
F = DpAm + DmAp + Ap * Am + Am * Ap
# expected: F = [(beta/a) D_-D_+Phi + eta (e^{beta a Phi} - e^{-beta a Phi})] h
#             = (i beta/a) Eq h   since D_-D_+ = -D_+D_- and eta = -i mu
target = h.scale(Eq.map(lambda x: I * beta / a * x))
report('zero curvature: D_+A_- + D_-A_+ + {A_+, A_-} = (i beta/|alpha|) (field equation) h for all '
       'lambda (eq-super-zero-curvature), with eta_0 = eta_1 = -i mu', (F - target).is_zero())

# ================================================= 5. conserved supercurrents
mm = sp.symbols('m', positive=True)
K = Jets(1, 'k')
P = K.superfield(0)
# L = (i/2) D+P D-P - m cosh P: superfield equation i D_+D_- P = -m sinh P
EqK = (K.D(K.D(P, -1), +1)).map(lambda x: I * x) + \
    (sexp(P) - sexp(P.map(lambda x: -x))).map(lambda x: mm / 2 * x)
E0, Ep, Em, E2 = (K.component(EqK, w) for w in ('0', '+', '-', '2'))
sK = K.sig(0)
sol_s = sp.solve(E0.body(), sK)[0]
r_psim = solve_linear_gen(Ep.map(lambda x: x.subs(sK, sol_s)), K.psi('m', 0, 1, 0))
r_psip = solve_linear_gen(Em.map(lambda x: x.subs(sK, sol_s)), K.psi('p', 0, 0, 1))
E2s = E2.map(lambda x: x.subs(sK, sol_s))
kk = sp.diff(E2s.body(), K.phi(0, 1, 1))
r_phi = (E2s - G.scalar(kk * K.phi(0, 1, 1))).map(lambda x: sp.expand(-x / kk))
OS = OnShell(K, sol_s, r_phi, r_psip, r_psim)
D1 = K.D(P, +1)
D2 = K.D(D1, +1)
D3 = K.D(D2, +1)
T = D1 * D2 - D3
em = sexp(P.map(lambda x: -x))
Tbar = K.D(em, +1).map(lambda x: -I * mm * x)
report('super sinh-Gordon: D_-T = D_+Tbar on shell with T = D_+Phi D_+^2Phi - D_+^3Phi and '
       'Tbar = -i m D_+ e^{-Phi} (eq-ssg-currents, spin 3/2)',
       OS.reduce(K.D(T, -1) - K.D(Tbar, +1)).is_zero())
J7 = T * K.D(T, +1)
X7 = (D1 * D2 * D2 - D2 * D3) * em
J7bar = X7.map(lambda x: I * mm * x)
report('super sinh-Gordon: D_-J = D_+Jbar on shell for J = T D_+T and '
       'Jbar = i m (D_+Phi (D_+^2Phi)^2 - D_+^2Phi D_+^3Phi) e^{-Phi} (eq-ssg-currents, spin 7/2)',
       OS.reduce(K.D(J7, -1) - K.D(J7bar, +1)).is_zero())
# the printed forms of Evans-Madsen 1999 (5.44)-(5.46): D_+J = D_-Jbar, Jbar = (m/2) i ...
pr45 = K.D(em, +1).map(lambda x: I * mm / 2 * x)
pr46 = X7.map(lambda x: I * mm / 2 * x)
bad1 = not OS.reduce(K.D(T, +1) - K.D(pr45, -1)).is_zero()
bad2 = not OS.reduce(K.D(T, -1) - K.D(pr45, +1)).is_zero()
bad3 = not OS.reduce(K.D(J7, -1) - K.D(pr46, +1)).is_zero()
report('the printed Evans-Madsen 1999 (5.44)-(5.46) do not hold, either as printed or with the '
       'D_+- exchanged (a misprint; sec-app-corrections)', bad1 and bad2 and bad3)
# spin-4 content: the theta^+ component of T D_+T at psi = 0
J7p = K.component(J7, '+')
dens = sp.factor(J7p.body())
report('at psi = 0 the theta^+ component of T D_+T is ((d_+phi)^2 - d_+^2 phi)^2, which differs from '
       '(d_+^2 phi)^2 + (d_+phi)^4 by a total d_+ derivative (the sinh-Gordon density of eq-spin3-current)',
       sp.expand(dens - ((K.phi(0, 1, 0)**2 - K.phi(0, 2, 0))**2)) == 0)

sys.exit(1 if FAILS else 0)
