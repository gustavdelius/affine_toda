"""Integration check (ch. 42, item 3): Corrigan-Zambon 1202.6016, direct fusion (2.9) of the general transmission
matrix (2.14) with the Dirichlet matrix (2.16), compared with their (2.17), (2.18) and (2.19).

CZ conventions (read from the PDF, printed pp. 3-6): T^{b beta}_{a alpha}, rows = incoming soliton a (and initial
charge alpha), columns = outgoing (b, beta); charge conservation a + alpha = b + beta; (2.9) is the row-incoming
product T(theta) R(theta) T^(theta) with T^(theta) = T(-theta)^{-1}.  x = e^{gamma theta}; Q^alpha = q^{-alpha/2},
so Q^2 = 1/q and q^alpha = Q^{-2 alpha}.  Symbol A stands for Q^alpha.  The conjugates bbar_+- are independent
symbols (|b|^2 = b bbar), which is more general than complex conjugation.

On the sector {(+, alpha), (-, alpha + 2)} of fixed a + alpha the fused matrix is the 2x2 matrix
  [[r_+(alpha), s_+(alpha)], [s_-(alpha + 2), r_-(alpha + 2)]]   (up to a common factor),
and T restricted to it is [[a(alpha), b(alpha)], [c(alpha + 2), d(alpha + 2)]].
"""
import sympy as sp

x, r, Q, A = sp.symbols('x r Q A', nonzero=True)
q = Q**-2
ap, am, bp, bm, cp, cm, dp, dm = sp.symbols('a_p a_m b_p b_m c_p c_m d_p d_m')
Bp, Bm, Bbp, Bbm = sp.symbols('b_+ b_- bb_+ bb_-')   # bb = bbar


def Tsec(xx, A_, pars):
    """(2.14) on the sector of (+, alpha) with Q^alpha = A_ (so Q^(alpha+2) = A_ Q^2)."""
    ap_, am_, bp_, bm_, cp_, cm_, dp_, dm_ = pars
    A2 = A_*Q**2
    a = ap_*A_/xx + am_*xx/A_
    b = bp_*A_ + bm_/A_
    c = cp_*A2 + cm_/A2
    d = dp_*A2*xx + dm_/(A2*xx)
    return sp.Matrix([[a, b], [c, d]])


RD = sp.diag(r*x + 1/(r*x), r/x + x/r)                  # (2.16), row-incoming = column-action (diagonal)


def fused(pars):
    return sp.simplify(Tsec(x, A, pars)*RD*Tsec(1/x, A, pars).inv())


def cz217(pars, A_):
    """(2.17) entries at Q^alpha = A_: r_+, r_-, s_+, s_- (q^alpha = A_^-2)."""
    ap_, am_, bp_, bm_, cp_, cm_, dp_, dm_ = pars
    qa = A_**-2
    g = am_*dp_/q + ap_*dm_*q
    h = bm_*cp_/q + bp_*cm_*q
    X = x**2 - x**-2
    rp = X*(am_*dm_*qa*q*r*x - ap_*dp_/(qa*q)/(r*x)) + x*(g*r - h/r) + (g/r - h*r)/x
    rm = X*(ap_*dp_*q/qa/r*x - am_*dm_*qa/q*r/x) + x*(g/r - h*r) + (g*r - h/r)/x
    sp_ = X*(ap_*bm_/r - am_*bp_*r - am_*bm_*qa*r + ap_*bp_/(qa*r))
    sm_ = X*(dm_*cp_*r - dp_*cm_/r + dm_*cm_*qa*r - dp_*cp_/(qa*r))
    return rp, rm, sp_, sm_


def compare(F, ent, label):
    """F = fused sector matrix; ent = (r_+(alpha), s_+(alpha), s_-(alpha+2), r_-(alpha+2)).  Ratio test."""
    M = sp.Matrix([[ent[0], ent[1]], [ent[2], ent[3]]])
    lam = sp.simplify(F[0, 0]/M[0, 0])
    diff = sp.simplify(F - lam*M)
    ok = diff == sp.zeros(2, 2)
    print(f'{label}: fused = c * formula entry by entry: {ok};  common factor c = {sp.factor(lam)}')
    return ok, lam, diff


# ---- 1. (2.17) for the general (2.14) with a_+- d_+- = b_+- c_+-
gen = (ap, am, bp, bm, cp, cm, dp, dm)
subs_gen = {dp: bp*cp/ap, dm: bm*cm/am}
genS = tuple(sp.sympify(p).subs(subs_gen) for p in gen)
F = fused(genS)
rp, rm, spl, smn = cz217(genS, A)
rp2, rm2, spl2, smn2 = cz217(genS, A*Q**2)            # entries at alpha + 2
ok217, lam217, _ = compare(F, (rp, spl, smn2, rm2), '(2.17), general type (2.14)')
print('   det of the sector T (CZ Delta) is alpha-independent:',
      sp.simplify(Tsec(x, A, genS).det() - Tsec(x, 3*A, genS).det()) == 0)

# ---- 2. type II substitution of CZ p. 5
typeII = (1, -Bbp*Bm, Bp, Bm, -Bbm*q**2, -Bbp, -Bp*Bbm*q**2, 1)
print('type II satisfies a d = b c:', sp.simplify(typeII[0]*typeII[6] - typeII[2]*typeII[4]) == 0,
      sp.simplify(typeII[1]*typeII[7] - typeII[3]*typeII[5]) == 0)
FII = fused(typeII)
# (2.17) with the type II parameters substituted
e = cz217(typeII, A); e2 = cz217(typeII, A*Q**2)
ok217II, _, _ = compare(FII, (e[0], e[2], e2[3], e2[1]), '(2.17) with the type-II parameters')
g_sub = sp.simplify(typeII[1]*typeII[6]/q + typeII[0]*typeII[7]*q)
h_sub = sp.simplify(typeII[3]*typeII[4]/q + typeII[2]*typeII[5]*q)
print('   g from (2.17) =', sp.factor(g_sub), ';   h from (2.17) =', sp.factor(h_sub))


def cz219(A_, hsel='printed', ssign=1):
    X = x**2 - x**-2
    g = (Bp*Bbp*Bm*Bbm + 1)*Q**-2
    h = -(Bm*Bbm*Q**-2 + Bp*Bbp)*Q**-2 if hsel == 'printed' else -(Bm*Bbm + Bp*Bbp)*Q**-2
    rp = X*(Bp*Bbm*A**2*Q**-2/(r*x) - Bbp*Bm*A**-2*Q**-2*r*x) + x*(g*r - h/r) + (g/r - h*r)/x
    rp = rp.subs(A, A_)
    rm = (X*(-Bp*Bbm*A**2*Q**-6/r*x + Bbp*Bm*A**-2*Q**2*r/x) + x*(g/r - h*r) + (g*r - h/r)/x).subs(A, A_)
    spl = (X*(Bm/r + Bp*Bbp*Bm*r + Bp*A**2/r + Bm**2*Bbp*A**-2*r)).subs(A, A_)
    smn = (ssign*X*(Bbm*Q**-4*r + Bp*Bbp*Bbm*Q**-4/r + Bbm**2*Bp*A**2*Q**-8/r + Bbp*A**-2*r)).subs(A, A_)
    return rp, rm, spl, smn


res = {}
for hsel in ('printed', 'corrected'):
    for ssign in (1, -1):
        e = cz219(A, hsel, ssign); e2 = cz219(A*Q**2, hsel, ssign)
        ok, _, _ = compare(FII, (e[0], e[2], e2[3], e2[1]), f'(2.19) with h {hsel}, s_- sign {ssign:+d}')
        res[(hsel, ssign)] = ok
# numerical size of the printed-form mismatch (gauge invariants at one point)
vals = {Q: sp.exp(-0.355j), Bp: 0.37 + 0.19j, Bm: -0.28 + 0.13j, r: sp.exp(0.23j), x: sp.exp(0.41)}
vals.update({Bbp: sp.conjugate(vals[Bp]), Bbm: sp.conjugate(vals[Bm])})
for hsel, ssign in (('printed', 1), ('corrected', 1), ('printed', -1), ('corrected', -1)):
    worst = 0
    for Aval in (1, sp.exp(-0.71j*sp.Rational(1, 2)*2), sp.exp(-0.71j*2)):
        e = cz219(Aval, hsel, ssign); e2 = cz219(Aval*Q**2, hsel, ssign)
        M = sp.Matrix([[e[0], e[2]], [e2[3], e2[1]]]).subs(vals)
        Fv = FII.subs(A, Aval).subs(vals)
        inv_f = [complex(Fv[0, 0]/Fv[1, 1]), complex(Fv[0, 1]*Fv[1, 0]/(Fv[0, 0]*Fv[1, 1]))]
        inv_m = [complex(M[0, 0]/M[1, 1]), complex(M[0, 1]*M[1, 0]/(M[0, 0]*M[1, 1]))]
        worst = max(worst, max(abs(a - b)/abs(a) for a, b in zip(inv_f, inv_m)))
    print(f'   gauge invariants (diag ratio, off-diag product), h {hsel}, s_- sign {ssign:+d}: max rel. dev. {worst:.3g}')

# charge dependence that no similarity removes (type II): diagonal ratio of the fused block at alpha = 0, 2, 4
print('type II fused diag ratio r_+(alpha)/r_-(alpha+2) at Q^alpha = 1, Q^2, Q^4:',
      [complex(sp.N((FII[0, 0]/FII[1, 1]).subs(A, Av).subs(vals))) for Av in (1, Q**2, Q**4)])

# ---- 3. type I (CZ p. 5: a_- = d_+ = c_- = b_+ = 0), (2.18)
nu, eps = sp.symbols('nu epsilon', positive=True)
aI, dI = sp.symbols('aI dI', nonzero=True)
# nu^{1/2} = (d_-/a_+)^{1/2}, b_- = c_+, eps = b_-/(a_+ d_-)^{1/2}; take a_+ = 1/sqrt(nu)*k, d_- = sqrt(nu)*k, k = 1
typeI = (1/sp.sqrt(nu), 0, 0, eps, eps, 0, 0, sp.sqrt(nu))
FI = fused(typeI)
X = x**2 - x**-2
rpI = x*(r*q - eps**2/(r*q)) + (q/r - eps**2*r/q)/x
rmI = x*(q/r - eps**2*r/q) + (r*q - eps**2/(r*q))/x
spI = X*eps/(r*sp.sqrt(nu))
smI = X*eps*r*sp.sqrt(nu)
okI, lamI, dI_ = compare(FI, (rpI, spI, smI, rmI), '(2.18), type I (multiplied by x as CZ do)')
if not okI:
    # (2.18) is stated after a similarity removing the alpha dependence of the off-diagonal entries: compare invariants
    inv_f = sp.simplify(FI[0, 0]/FI[1, 1] - rpI/rmI), sp.simplify(FI[0, 1]*FI[1, 0]/(FI[0, 0]*FI[1, 1]) - spI*smI/(rpI*rmI))
    print('   gauge invariants of (2.18) vs fusion:', inv_f)
