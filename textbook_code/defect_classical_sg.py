"""Classical sinh-Gordon and sine-Gordon with a jump defect, in the book's normalization (Part VIII, ch. 36-37).

Conventions: defect_common.py.  Field phi on x < 0, psi on x > 0; a_1: D = 1, E = 0 (type I_-), D = -1 (type I_+);
sinh-Gordon V = (m^2/beta^2)(e^{sqrt2 beta phi} + e^{-sqrt2 beta phi} - 2) (eq-imaginary-potential, real coupling);
sine-Gordon beta_SG = sqrt2 beta, m0 = 2m (sec-sg-dictionary); soliton eq-an-soliton with h = 2:
phi = (2i/beta_SG) ln((1 - E)/(1 + E)), E = exp(m0(x cosh theta - t sinh theta) + xi).

Checks (book ids they support):
 1. Defect conditions from the action by variation (sympy), a_1, both types.                    [eq-defect-action, eq-sg-defect-conditions]
 2. Integrability/momentum conditions (sympy): 1/2 B_phi^2 - 1/2 B_psi^2 = V(phi) - V(psi), B_phiphi = B_psipsi,
    and B_phi = -dP_d/dpsi, B_psi = -dP_d/dphi, so that P + P_d is conserved; sinh- and sine-Gordon.  [eq-sg-defect-potential, prp-sg-defect-momentum]
 3. Dictionary (sympy): BCZ hep-th/0305022 (2.13) (+ signs), BCZ hep-th/0506169 (2.3), (2.4) and their U;
    BCZ 2003 sigma = e^p has p = -theta_d; BCZ 2005 sigma = e^{-eta} has eta = theta_d.        [tbl-defect-dictionary]
 4. A sine-Gordon soliton (complex xi) on x < 0 and the shifted soliton E -> z E on x > 0 satisfy the defect conditions,
    with z found by root-finding = z_{I-} = (e^{-theta} + sigma)/(e^{-theta} - sigma) = -coth((theta - theta_d)/2),
    z_{I+} = +coth((theta - theta_d)/2); at theta = theta_d the soliton is absorbed (psi = 2 pi/beta_SG).  [eq-sg-delay, prp-sg-delay]
 5. E + B_d and P + P_d are conserved along these solutions (real kink, theta < theta_d and theta > theta_d, both
    types); totals M cosh theta, M sinh theta, M = 8 m0/beta_SG^2; the absorbed soliton is stored as
    B_d = M cosh theta_d, P_d = M sinh theta_d.                                                  [prp-sg-defect-momentum, prp-half-soliton]
 6. Linearization (sympy): R = 0 and T_cl(theta) = -i sinh((theta-theta_d)/2 - i pi/4)/sinh((theta-theta_d)/2 + i pi/4)
    (ratio of transmitted to incident amplitude, plane waves e^{-i(Et - px)}); about the odd vacuum pair
    (phi, psi) = (0, 2 pi/beta_SG) the sine-Gordon factor is T_odd = conj(T_even) (BCZ 2005 (2.9)).  [eq-shg-linear-transmission]
Chapter 37 (light-cone d_+- = d_t +- d_x; eq-sg-backlund: d_+(phi - psi) = (2 a m0/beta) sin(beta(phi+psi)/2),
d_-(phi + psi) = -(2 m0/(a beta)) sin(beta(phi-psi)/2), beta = beta_SG):
 7. (b) Backlund: eq-sg-backlund implies sine-Gordon for both fields (sympy); vacuum -> one soliton with a = -e^{-theta}
    (any complex xi); Bianchi permutability tan(beta(u12-u0)/4) = ((a1+a2)/(a1-a2)) tan(beta(u1-u2)/4) gives the
    two-soliton with A_11 = ((a1-a2)/(a1+a2))^2 = tanh^2 (eq-interaction-coefficient); the type I_- defect conditions are
    eq-sg-backlund with a = -sigma at x = 0, and z_{I-} is the Bianchi coefficient.  [eq-sg-backlund, prp-sg-backlund, prp-bianchi]
 8. a_1: type I_+ = type I_- with psi -> -psi (B_d, P_d, kinetic term); linear waves: I_+ gives -T_cl; a wave from the
    right has Tbar_cl(theta) = T_cl(i pi - theta) (depends on theta + theta_d); T_cl T_cl(+i pi) = 1, T_cl(theta) T_cl(-theta) != 1.
                                                                    [eq-sg-defect-conditions, eq-shg-linear-transmission]
 9. theta = theta_d: absorption and emission configurations for theta_d > 0 and < 0; two solitons through the defect,
    each E_k -> z(theta_k) E_k (both types).                                                    [prp-sg-delay]
10. (e) Half a soliton (pre-registered): odd-sector delay z_odd = 1/z_even = -tanh((theta-theta_d)/2); z_odd^2 = A_11(theta-theta_d),
    z_even^2 = A_11(theta - theta_d - i pi); with BCZ's constant the defect's (E, P) lies on the mass shell of mass M/2 at
    rapidity theta_d + i pi (even) or theta_d (odd); T_odd^2 = f_{1/2}(theta - theta_d).        [prp-half-soliton]
11. (g) Theta: conj B_d(phi, psi) = B_d(-conj phi, -conj psi) iff sigma real; Theta maps soliton solutions through the
    defect to solutions iff sigma real.                                                         [sec-sg-defect-lagrangian]
12. psi = phi_D turns B_d into the GZ boundary potential with C_{0,1} = cosh(theta_d +- i beta phi_D/2).  [sec-sg-breather-transmission]
13. The sinh-Gordon B_d satisfies the spin-3 conditions BCZ 2003 (2.6)-(2.9).                 [prp-sg-defect-classification]
Runtime: about 20 seconds.
"""
import numpy as np
import sympy as sp
import mpmath as mp
from defect_common import defect_potential, bulk_potential, defect_condition_residual


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


# ---------------------------------------------------------------- 1. variation of the action
t = sp.symbols('t', real=True)
ph, ps = sp.Function('phi')(t), sp.Function('psi')(t)
phx, psx = sp.symbols('phi_x psi_x')
Bf = sp.Function('B')
ok = True
for typ, Dv in [(1, 1), (-1, -1)]:
    Ld = ph*Dv*sp.diff(ps, t) - Bf(ph, ps)              # E = 0 for a_1
    # bulk boundary terms of the variation: -phi_x delta phi (x<0 region), +psi_x delta psi (x>0 region)
    el_phi = sp.diff(Ld, ph) - sp.diff(sp.diff(Ld, sp.diff(ph, t)), t) - phx
    el_psi = sp.diff(Ld, ps) - sp.diff(sp.diff(Ld, sp.diff(ps, t)), t) + psx
    want_phx = Dv*sp.diff(ps, t) - sp.diff(Bf(ph, ps), ph)
    want_psx = Dv*sp.diff(ph, t) + sp.diff(Bf(ph, ps), ps)
    ok &= sp.simplify(sp.solve(el_phi, phx)[0] - want_phx) == 0
    ok &= sp.simplify(sp.solve(el_psi, psx)[0] - want_psx) == 0
report('variation: phi_x = D psi_t - dB/dphi, psi_x = D phi_t + dB/dpsi at x = 0 (a_1, D = +-1)', ok)

# ---------------------------------------------------------------- 2. integrability conditions (sympy)
x_, y_, s_, b_, m_ = sp.symbols('phi psi sigma beta m', positive=True)
r2 = sp.sqrt(2)


def B_sym(bb, typ=1):
    s = s_ if typ == 1 else 1/s_
    f = sp.exp(bb*r2*(x_ + y_)/2) + sp.exp(-bb*r2*(x_ + y_)/2) - 2       # alpha_{0,1} = -+sqrt2, D = 1
    g = sp.exp(bb*r2*(x_ - y_)/2) + sp.exp(-bb*r2*(x_ - y_)/2) - 2
    return (m_/bb**2)*(s*f + g/s), -(m_/bb**2)*(s*f - g/s)*(1 if typ == 1 else -1)


def V_sym(bb, f):
    return (m_**2/bb**2)*(sp.exp(r2*bb*f) + sp.exp(-r2*bb*f) - 2)


ok = True
for bb, name in [(b_, 'sinh'), (sp.I*b_, 'sine')]:
    for typ in (1, -1):
        B, Pd = B_sym(bb, typ)
        Bx, By = sp.diff(B, x_), sp.diff(B, y_)
        ok &= sp.simplify((Bx**2 - By**2)/2 - (V_sym(bb, x_) - V_sym(bb, y_))) == 0
        ok &= sp.simplify(sp.diff(B, x_, 2) - sp.diff(B, y_, 2)) == 0
        # linear terms of dP/dt: D (psi_t B_phi + phi_t B_psi) = -dP_d/dt
        Dv = 1 if typ == 1 else -1
        ok &= sp.simplify(Dv*Bx + sp.diff(Pd, y_)) == 0 and sp.simplify(Dv*By + sp.diff(Pd, x_)) == 0
report('1/2 B_phi^2 - 1/2 B_psi^2 = V(phi) - V(psi), B_phiphi = B_psipsi, P_d with dP/dt = -dP_d/dt (sinh-, sine-Gordon; types I_-, I_+)', ok)

# numerical cross-check of the helper (real and imaginary coupling)
rng = np.random.default_rng(1)
ok = True
for coupling in ('real', 'imag'):
    for typ in (1, -1):
        p, q = rng.normal(size=1) + 0.3j*rng.normal(size=1), rng.normal(size=1) + 0.3j*rng.normal(size=1)
        B, gp, gq, Pd = defect_potential(p, q, 0.63, 1, beta=0.8, m=1.3, coupling=coupling, typ=typ)
        bb = 0.8 if coupling == 'real' else 0.8j
        Bs, Ps = B_sym(bb, typ)
        sub = {x_: complex(p[0]), y_: complex(q[0]), s_: 0.63, m_: 1.3}
        ok &= abs(complex(Bs.subs(sub)) - B) < 1e-12 and abs(complex(Ps.subs(sub)) - Pd) < 1e-12
        ok &= abs(complex(sp.diff(Bs, x_).subs(sub)) - gp[0]) < 1e-12 and abs(complex(sp.diff(Bs, y_).subs(sub)) - gq[0]) < 1e-12
report('defect_common.defect_potential agrees with the closed form (a_1, both couplings and types)', ok)

# ---------------------------------------------------------------- 3. dictionary to BCZ
# BCZ 2003 (2.13), '+' signs: (m, beta, sigma, phi1, phi2) = book (m, beta, sigma, phi, psi), up to constants
B03 = (m_*s_/b_**2)*(sp.exp(b_*(x_ + y_)/r2) + sp.exp(-b_*(x_ + y_)/r2)) + (m_/(b_**2*s_))*(sp.exp(b_*(x_ - y_)/r2) + sp.exp(-b_*(x_ - y_)/r2))
V03 = (m_**2/b_**2)*(sp.exp(r2*b_*x_) + sp.exp(-r2*b_*x_))
ok03 = sp.simplify(B03 - B_sym(b_)[0] - 2*m_*(s_ + 1/s_)/b_**2) == 0 and sp.simplify(V03 - V_sym(b_, x_) - 2*m_**2/b_**2) == 0
# BCZ 2005 (2.1), (2.3): (m, beta, u, v, sigma) = book (m0, beta_SG, phi, psi, sigma), m0 = 2m, beta_SG = sqrt2 beta
m0, bSG = 2*m_, r2*b_
B05 = -(2*m0*s_/bSG**2)*sp.cos(bSG*(x_ + y_)/2) - (2*m0/(s_*bSG**2))*sp.cos(bSG*(x_ - y_)/2)
Bimag = B_sym(sp.I*b_)[0]
ok05 = sp.simplify((B05 + 2*m0*(s_ + 1/s_)/bSG**2 - Bimag).rewrite(sp.exp)) == 0
# their scaled conditions (2.4) (m = beta = 1 in (2.1)): u_x - v_t = -sigma sin((u+v)/2) - sin((u-v)/2)/sigma, ...
Bsc = B05.subs(b_, 1/r2).subs(m_, sp.Rational(1, 2))              # m0 = beta_SG = 1
ok05 &= sp.simplify(-sp.diff(Bsc, x_) - (-s_*sp.sin((x_ + y_)/2) - sp.sin((x_ - y_)/2)/s_)) == 0
ok05 &= sp.simplify(sp.diff(Bsc, y_) - (s_*sp.sin((x_ + y_)/2) - sp.sin((x_ - y_)/2)/s_)) == 0
# their U = -2 sigma cos((u+v)/2) + (2/sigma) cos((u-v)/2) (scaled) and P_BCZ = int u_t u_x = -P_book: U = -P_d + const
Pd_sc = B_sym(sp.I*b_)[1].subs(b_, 1/r2).subs(m_, sp.Rational(1, 2))
U05 = -2*s_*sp.cos((x_ + y_)/2) + (2/s_)*sp.cos((x_ - y_)/2)
ok05 &= sp.simplify((U05 + Pd_sc).rewrite(sp.exp) - (-2*s_ + 2/s_)) == 0
report('BCZ 2003 (2.13): B = book B_d + 2m(sigma + 1/sigma)/beta^2, V = book V + 2m^2/beta^2 (same m, beta, sigma, fields)', ok03)
report('BCZ 2005 (2.3): (m, beta, u, v) = (m0, beta_SG, phi, psi), same sigma; (2.4) reproduced; U = -P_d + const', ok05)

# ---------------------------------------------------------------- 4. soliton through the defect
mp.mp.dps = 30
BETA = mp.mpf('0.9')                    # book beta; beta_SG = sqrt2 beta
M_ = mp.mpf('0.7')                      # book m; m0 = 2m
bsg, m0n = mp.sqrt(2)*BETA, 2*M_


def sol(x, t, th, xi, z=1):
    E = z*mp.e**(m0n*(x*mp.cosh(th) - t*mp.sinh(th)) + xi)
    return 2j/bsg*mp.log((1 - E)/(1 + E))


def resid(th, xi, sig, z, t0, typ):
    p = sol(0, t0, th, xi); q = sol(0, t0, th, xi, z)
    px = mp.diff(lambda x: sol(x, t0, th, xi), 0); qx = mp.diff(lambda x: sol(x, t0, th, xi, z), 0)
    pt = mp.diff(lambda t: sol(0, t, th, xi), t0); qt = mp.diff(lambda t: sol(0, t, th, xi, z), t0)
    c = lambda v: np.array([complex(v)])
    return defect_condition_residual(c(p), c(q), c(px), c(qx), c(pt), c(qt), float(sig), 1, beta=float(BETA),
                                     m=float(M_), coupling='imag', typ=typ)


def zI(th, thd):
    return -mp.coth((th - thd)/2)


zclosed = {1: zI, -1: lambda th, thd: -zI(th, thd)}
ok = True; worst = 0; okfit = True
thd = mp.mpf('0.4'); sig = mp.e**(-thd)
for typ in (1, -1):
    for th in (mp.mpf('-0.3'), mp.mpf('0.25'), mp.mpf('1.1')):
        xi = mp.mpc('0.2', '0.9')
        # derive z: solve the first defect condition at one point for z, then test everywhere
        def f(z):
            p = sol(0, 0.1, th, xi); q = sol(0, 0.1, th, xi, z)
            px = mp.diff(lambda x: sol(x, 0.1, th, xi), 0); qt = mp.diff(lambda t: sol(0, t, th, xi, z), 0.1)
            _, gp, _, _ = defect_potential(np.array([complex(p)]), np.array([complex(q)]), float(sig), 1,
                                           beta=float(BETA), m=float(M_), coupling='imag', typ=typ)
            return complex(px) - (typ*complex(qt) - gp[0])
        zfit = mp.findroot(lambda z: f(complex(z)), complex(zclosed[typ](th, thd)) * 1.05, tol=1e-20)
        okfit &= abs(zfit - zclosed[typ](th, thd)) < 1e-8
        for t0 in (-0.7, 0.0, 0.35, 1.3):
            worst = max(worst, resid(th, xi, sig, zclosed[typ](th, thd), t0, typ))
report(f'root-finding gives z_{{I-}} = (e^-theta + sigma)/(e^-theta - sigma) = -coth((theta-theta_d)/2), z_{{I+}} = coth((theta-theta_d)/2)', okfit)
report(f'soliton (complex xi) + shifted soliton E -> zE satisfy both defect conditions at 4 times, 3 rapidities, types I_-, I_+  [res {worst:.1e}]', worst < 1e-10)
# absorption at theta = theta_d: type I_- psi = 2 pi/beta_SG, type I_+ psi = 0
worst = 0
for typ, psi0 in ((1, 2*mp.pi/bsg), (-1, 2*mp.pi/bsg)):
    for t0 in (-0.5, 0.3, 1.0):
        xi = mp.mpc('0.1', '0.7')
        p = sol(0, t0, thd, xi); px = mp.diff(lambda x: sol(x, t0, thd, xi), 0); pt = mp.diff(lambda t: sol(0, t, thd, xi), t0)
        c = lambda v: np.array([complex(v)])
        worst = max(worst, defect_condition_residual(c(p), c(psi0), c(px), c(0), c(pt), c(0), float(sig), 1,
                                                     beta=float(BETA), m=float(M_), coupling='imag', typ=typ))
report(f'theta = theta_d (z_{{I-}} = z_{{I+}} = inf): soliton absorbed, psi = 2 pi/beta_SG constant on x > 0, both types  [res {worst:.1e}]', worst < 1e-12)

# ---------------------------------------------------------------- 5. conserved E + B_d, P + P_d
mp.mp.dps = 20
Mcl = 8*m0n/bsg**2


def kink(x, t, th, xr, z=1):
    Y = z*mp.e**(m0n*(x*mp.cosh(th) - t*mp.sinh(th)) + xr)
    phi = 4/bsg*mp.atan(Y)
    d = 4/bsg/(1 + Y**2)*Y*m0n
    return phi, d*mp.cosh(th), -d*mp.sinh(th)          # phi, phi_x, phi_t


def totals(t, th, xr, z, typ, sig):
    def dens(x, zz):
        f, fx, ft = kink(x, t, th, xr, zz)
        V = (m0n**2/bsg**2)*(1 - mp.cos(bsg*f))
        return ft**2/2 + fx**2/2 + V, -ft*fx
    L = 40/m0n
    pts_l = [-L, -L/2, -L/4, -L/8, 0]; pts_r = [0, L/8, L/4, L/2, L]
    El = mp.quad(lambda x: dens(x, 1)[0], pts_l); Er = mp.quad(lambda x: dens(x, z)[0], pts_r)
    Pl = mp.quad(lambda x: dens(x, 1)[1], pts_l); Pr = mp.quad(lambda x: dens(x, z)[1], pts_r)
    p0, q0 = kink(0, t, th, xr)[0], kink(0, t, th, xr, z)[0]
    B, _, _, Pd = defect_potential(np.array([float(p0)]), np.array([float(q0)]), float(sig), 1, beta=float(BETA),
                                   m=float(M_), coupling='imag', typ=typ)
    return El + Er + B.real, Pl + Pr + Pd.real


ok = True; dev = 0
for typ in (1, -1):
    for th in (mp.mpf('0.15'), mp.mpf('0.9')):
        z = zclosed[typ](th, thd)
        vals = [totals(tt, th, mp.mpf('0.3'), z, typ, sig) for tt in (-6, -0.5, 0.0, 0.4, 6)]
        Es = [float(v[0]) for v in vals]; Ps = [float(v[1]) for v in vals]
        dev = max(dev, np.ptp(Es)/float(Mcl), np.ptp(Ps)/float(Mcl))
        ok &= abs(Es[0] - float(Mcl*mp.cosh(th))) < 1e-8*float(Mcl) and abs(Ps[0] - float(Mcl*mp.sinh(th))) < 1e-8*float(Mcl)
report(f'E + B_d and P + P_d (P = -int phi_t phi_x) constant through the defect, = M cosh theta, M sinh theta '
       f'(theta < theta_d and > theta_d; types I_-, I_+)  [rel. spread {dev:.1e}]', ok and dev < 1e-9)
B, _, _, Pd = defect_potential(np.array([0.0]), np.array([float(2*mp.pi/bsg)]), float(sig), 1, beta=float(BETA), m=float(M_), coupling='imag')
report('a 2 pi/beta_SG jump stores B_d = M cosh theta_d, P_d = M sinh theta_d (sigma = e^{-theta_d})',
       abs(B - float(Mcl*mp.cosh(thd))) < 1e-12 and abs(Pd - float(Mcl*mp.sinh(thd))) < 1e-12)

# ---------------------------------------------------------------- 6. linearization
th, thd_s = sp.symbols('theta theta_d', real=True)
R, T = sp.symbols('R T')
mm = sp.symbols('m', positive=True)
sig_s = sp.exp(-thd_s)
En, pn = 2*mm*sp.cosh(th), 2*mm*sp.sinh(th)


def lin_T(shift_sign):
    """shift_sign = +1: sinh-Gordon / even sine-Gordon vacua; -1: sine-Gordon about (0, 2 pi/beta_SG)."""
    # quadratic part of B_d: (m/2)[sigma (phi+psi)^2 + sigma^-1 (phi-psi)^2], with sigma -> -sigma (odd vacuum pair)
    s = shift_sign*sig_s
    phi0, phix, phit = 1 + R, sp.I*pn*(1 - R), -sp.I*En*(1 + R)
    psi0, psix, psit = T, sp.I*pn*T, -sp.I*En*T
    Bphi = mm*(s*(phi0 + psi0) + (phi0 - psi0)/s)
    Bpsi = mm*(s*(phi0 + psi0) - (phi0 - psi0)/s)
    eqs = [phix - (psit - Bphi), psix - (phit + Bpsi)]
    return sp.solve(eqs, [R, T], dict=True)[0]


sol_e = lin_T(1)
Tcl = -sp.I*sp.sinh((th - thd_s)/2 - sp.I*sp.pi/4)/sp.sinh((th - thd_s)/2 + sp.I*sp.pi/4)
ok_e = sp.simplify(sol_e[R]) == 0
chk = [complex((sol_e[T] - Tcl).subs({th: a, thd_s: b, mm: 0.7}).evalf()) for a, b in [(0.3, -0.4), (1.2, 0.5), (-0.7, 0.2)]]
ok_e &= max(abs(c) for c in chk) < 1e-12
report('linearized defect: R = 0, T_cl = -i sinh((theta-theta_d)/2 - i pi/4)/sinh((theta-theta_d)/2 + i pi/4) (= BCZ 2003 (3.1), p = -theta_d)', ok_e)
sol_o = lin_T(-1)
chk = [complex((sol_o[T] - sp.conjugate(Tcl)).subs({th: a, thd_s: b, mm: 0.7}).evalf()) for a, b in [(0.3, -0.4), (1.2, 0.5)]]
report('sine-Gordon about (phi, psi) = (0, 2 pi/beta_SG): R = 0, T_odd = conj(T_even) (BCZ 2005 (2.9))',
       sp.simplify(sol_o[R]) == 0 and max(abs(c) for c in chk) < 1e-12)
chk = [abs(complex(Tcl.subs({th: a, thd_s: b}).evalf())) for a, b in [(0.3, -0.4), (1.2, 0.5)]]
v0 = complex(Tcl.subs({th: 0.2, thd_s: 0.2}).evalf())
report(f'|T_cl| = 1 for real rapidities; T_cl(theta = theta_d) = i; T_cl -> 1 as theta_d -> +inf',
       max(abs(c - 1) for c in chk) < 1e-12 and abs(v0 - 1j) < 1e-12 and abs(complex(Tcl.subs({th: 0.3, thd_s: 40}).evalf()) - 1) < 1e-8)

# ======================================================================== ch. 37 additions
# Light-cone derivatives d_+- = d_t +- d_x (sec-lagrangians), so d_+ d_- = d_t^2 - d_x^2 and the sine-Gordon equation is
# d_+ d_- phi = -(m0^2/beta) sin(beta phi), beta = beta_SG.  Backlund transformation (eq-sg-backlund), phi -> psi:
#   d_+(phi - psi) = (2 a m0/beta) sin(beta(phi + psi)/2),   d_-(phi + psi) = -(2 m0/(a beta)) sin(beta(phi - psi)/2).

# ---------------------------------------------------------------- 7. Backlund transformations
bS, mS, aS = sp.symbols('beta m0 a', positive=True)
Ssum, Sdif = sp.symbols('s d')                      # s = phi + psi, d = phi - psi at a point
F = (2*aS*mS/bS)*sp.sin(bS*Ssum/2)                  # = d_+(phi - psi)
G = -(2*mS/(aS*bS))*sp.sin(bS*Sdif/2)               # = d_-(phi + psi)
pmd = sp.diff(F, Ssum)*G                            # d_- d_+ (phi - psi) = F_s s_-
pms = sp.diff(G, Sdif)*F                            # d_+ d_- (phi + psi) = G_d d_+
phi_pm = (pms + pmd)/2; psi_pm = (pms - pmd)/2
ok = sp.simplify(phi_pm + (mS**2/bS)*sp.sin(bS*(Ssum + Sdif)/2)) == 0
ok &= sp.simplify(psi_pm + (mS**2/bS)*sp.sin(bS*(Ssum - Sdif)/2)) == 0
report('eq-sg-backlund is compatible: it implies d_+d_- phi = -(m0^2/beta) sin beta phi and the same for psi', ok)

mp.mp.dps = 25


def bt_res(f, g, a, x0, t0):
    """max residual of eq-sg-backlund for phi = f(x,t) -> psi = g(x,t) at (x0, t0)."""
    dpl = lambda h: mp.diff(lambda s: h(x0 + s, t0 + s), 0)             # d_+ = d_t + d_x
    dmi = lambda h: mp.diff(lambda s: h(x0 - s, t0 + s), 0)             # d_- = d_t - d_x
    f0, g0 = f(x0, t0), g(x0, t0)
    r1 = dpl(lambda x, t: f(x, t) - g(x, t)) - (2*a*m0n/bsg)*mp.sin(bsg*(f0 + g0)/2)
    r2 = dmi(lambda x, t: f(x, t) + g(x, t)) + (2*m0n/(a*bsg))*mp.sin(bsg*(f0 - g0)/2)
    return max(abs(r1), abs(r2))


def hirota1(th, xi):
    return lambda x, t: sol(x, t, th, xi)


pts = [(-0.4, 0.2), (0.3, -0.5), (1.1, 0.7)]
worst = 0
for th_v in (mp.mpf('-0.6'), mp.mpf('0.35'), mp.mpf('1.2')):
    for xi in (mp.mpc('0.3', '1.5707963267948966'), mp.mpc('-0.2', '-1.5707963267948966'), mp.mpc('0.1', '0.4')):
        for x0, t0 in pts:
            worst = max(worst, bt_res(lambda x, t: 0, hirota1(th_v, xi), -mp.e**(-th_v), x0, t0))
report(f'prp-sg-backlund: vacuum -> one soliton (eq-an-soliton, h = 2, any complex xi: kink, antikink, complex) with a = -e^(-theta) '
       f'[res {float(worst):.0e}]', worst < 1e-15)

# Bianchi permutability
th1, th2 = mp.mpf('0.9'), mp.mpf('-0.4')
a1, a2 = -mp.e**(-th1), -mp.e**(-th2)
c1, c2 = mp.mpf('0.3'), mp.mpf('-0.5')
Yf = lambda th, c: (lambda x, t: mp.e**(m0n*(x*mp.cosh(th) - t*mp.sinh(th)) + c))
u1 = lambda x, t: 4/bsg*mp.atan(Yf(th1, c1)(x, t))
u2 = lambda x, t: 4/bsg*mp.atan(Yf(th2, c2)(x, t))
kappa = (a1 + a2)/(a1 - a2)
u12 = lambda x, t: 4/bsg*mp.atan(kappa*mp.tan(bsg*(u1(x, t) - u2(x, t))/4))
worst = 0
for x0, t0 in pts:
    worst = max(worst, bt_res(lambda x, t: 0, u1, a1, x0, t0), bt_res(lambda x, t: 0, u2, a2, x0, t0),
                bt_res(u1, u12, a2, x0, t0), bt_res(u2, u12, a1, x0, t0))
# u12 is the Hirota two-soliton (eq-n-soliton) with A_11(theta1 - theta2) = tanh^2: tan(beta u/4) = (X1+X2)/(1 - A X1 X2),
# X1 = -Y1/t, X2 = Y2/t, t = tanh((theta1-theta2)/2)
Dl = th1 - th2; tt = mp.tanh(Dl/2)
A11 = (mp.cosh(Dl) - mp.cos(0))/(mp.cosh(Dl) - mp.cos(mp.pi))        # eq-interaction-coefficient, h = 2, a = b = 1
devH = 0
for x0, t0 in pts:
    X1, X2 = -Yf(th1, c1)(x0, t0)/tt, Yf(th2, c2)(x0, t0)/tt
    devH = max(devH, abs(mp.tan(bsg*u12(x0, t0)/4) - (X1 + X2)/(1 - A11*X1*X2)))
okA = abs(A11 - ((a1 - a2)/(a1 + a2))**2) < 1e-20 and abs(A11 - mp.tanh(Dl/2)**2) < 1e-20
report(f'prp-bianchi: tan(beta(u12-u0)/4) = ((a1+a2)/(a1-a2)) tan(beta(u1-u2)/4) gives u12 = B_a2 u1 = B_a1 u2 [res {float(worst):.0e}], '
       f'the two-soliton with A_11 = ((a1-a2)/(a1+a2))^2 = tanh^2((theta1-theta2)/2) of eq-interaction-coefficient [{float(devH):.0e}]',
       worst < 1e-15 and devH < 1e-15 and okA)

# the defect conditions are the Backlund transformation with a = -sigma, frozen at x = 0 (type I_-)
ph_, ps_, sg_ = sp.symbols('phi psi sigma', positive=True)
Bd = -(2*mS/bS**2)*(sg_*(sp.cos(bS*(ph_ + ps_)/2) - 1) + (sp.cos(bS*(ph_ - ps_)/2) - 1)/sg_)
pt_, st_ = sp.symbols('phi_t psi_t')
phx_ = st_ - sp.diff(Bd, ph_); psx_ = pt_ + sp.diff(Bd, ps_)          # eq-sg-defect-conditions, D = 1
dplus = (pt_ + phx_) - (st_ + psx_)                                   # d_+(phi - psi)
dminus = (pt_ - phx_) + (st_ - psx_)                                  # d_-(phi + psi)
ok = sp.simplify(dplus - (2*(-sg_)*mS/bS)*sp.sin(bS*(ph_ + ps_)/2)) == 0
ok &= sp.simplify(dminus + (2*mS/((-sg_)*bS))*sp.sin(bS*(ph_ - ps_)/2)) == 0
th_s, thd_ = sp.symbols('theta theta_d', real=True)
aa = lambda x: -sp.exp(-x)
ok &= sp.simplify(((aa(th_s) + aa(thd_))/(aa(th_s) - aa(thd_)) - (sp.exp(-th_s) + sp.exp(-thd_))/(sp.exp(-th_s) - sp.exp(-thd_)))) == 0
report('type I_- defect conditions = eq-sg-backlund with a = -sigma = a(theta_d) at x = 0; z_{I-}(theta) = (a(theta) + a(theta_d))/'
       '(a(theta) - a(theta_d)) is the Bianchi coefficient', ok)

# ---------------------------------------------------------------- 8. type I_+ = type I_- with psi -> -psi (a_1); linear waves
ok = True
for bb in (b_, sp.I*b_):
    Bp, Pp = B_sym(bb, -1); Bm, Pm = B_sym(bb, 1)
    ok &= sp.simplify((Bp - Bm.subs(y_, -y_)).rewrite(sp.exp)) == 0
    ok &= sp.simplify((Pp - Pm.subs(y_, -y_)).rewrite(sp.exp)) == 0
report('a_1: B_d^{I+}(phi, psi) = B_d^{I-}(phi, -psi), P_d^{I+}(phi, psi) = P_d^{I-}(phi, -psi), and the kinetic term phi D psi_t '
       'changes sign: type I_+ is type I_- followed by psi -> -psi (so z_{I+} = -z_{I-})', ok)


def lin_general(Dv, sigsgn, from_right=False):
    """linearized defect conditions about (0,0), types D = +-1; from_right: a wave of rapidity -theta from x > 0."""
    s = sigsgn*sig_s
    sB = s if Dv == 1 else 1/s                            # I_+: sigma -> 1/sigma
    if not from_right:
        phi0, phix, phit = 1 + R, sp.I*pn*(1 - R), -sp.I*En*(1 + R)
        psi0, psix, psit = T, sp.I*pn*T, -sp.I*En*T
    else:                                                 # psi = e^{-i(Et + px)} + R e^{-i(Et - px)}, phi = T e^{-i(Et + px)}
        psi0, psix, psit = 1 + R, -sp.I*pn*(1 - R), -sp.I*En*(1 + R)
        phi0, phix, phit = T, -sp.I*pn*T, -sp.I*En*T
    Bphi = mm*(sB*(phi0 + psi0) + (phi0 - psi0)/sB)
    Bpsi = mm*(sB*(phi0 + psi0) - (phi0 - psi0)/sB)
    eqs = [phix - (Dv*psit - Bphi), psix - (Dv*phit + Bpsi)]
    return sp.solve(eqs, [R, T], dict=True)[0]


solp = lin_general(-1, 1)
chk = [complex((solp[T] + Tcl).subs({th: a, thd_s: b, mm: 0.7}).evalf()) for a, b in [(0.3, -0.4), (1.2, 0.5)]]
ok1 = sp.simplify(solp[R]) == 0 and max(abs(c) for c in chk) < 1e-12
solr = lin_general(1, 1, from_right=True)
Tbar = Tcl.subs(th, sp.I*sp.pi - th)
chk = [complex((solr[T] - Tbar).subs({th: a, thd_s: b, mm: 0.7}).evalf()) for a, b in [(0.3, -0.4), (1.2, 0.5), (-0.6, 0.1)]]
ok2 = sp.simplify(solr[R]) == 0 and max(abs(c) for c in chk) < 1e-12
uni = [complex((Tcl*Tcl.subs(th, th + sp.I*sp.pi)).subs({th: a, thd_s: b}).evalf()) for a, b in [(0.3, -0.4), (0.2 + 0.5j, 1.1)]]
noP = complex((Tcl*Tcl.subs(th, -th)).subs({th: 0.7, thd_s: 0.3}).evalf())
report('type I_+: R = 0, T = -T_cl; a wave from the right: R = 0, Tbar_cl(theta) = T_cl(i pi - theta), a function of theta + theta_d; '
       f'T_cl(theta) T_cl(theta + i pi) = 1, but |T_cl(theta) T_cl(-theta) - 1| = {abs(noP - 1):.2f} (no parity)',
       ok1 and ok2 and max(abs(u - 1) for u in uni) < 1e-12 and abs(noP - 1) > 0.1)

# ---------------------------------------------------------------- 9. absorption, emission, several solitons
mp.mp.dps = 30
c_ = lambda v: np.array([complex(v)])


def cond_res(fphi, fpsi, t0, sgm, typ=1):
    p = fphi(0, t0); s_ = fpsi(0, t0)
    px = mp.diff(lambda x: fphi(x, t0), 0); sx = mp.diff(lambda x: fpsi(x, t0), 0)
    pt = mp.diff(lambda t: fphi(0, t), t0); st = mp.diff(lambda t: fpsi(0, t), t0)
    return defect_condition_residual(c_(p), c_(s_), c_(px), c_(sx), c_(pt), c_(st), float(sgm), 1, beta=float(BETA),
                                     m=float(M_), coupling='imag', typ=typ)


worst = 0
for thd_v in (mp.mpf('0.4'), mp.mpf('-0.5')):
    sgm = mp.e**(-thd_v)
    kink_ = lambda x, t: 4/bsg*mp.atan(mp.e**(m0n*(x*mp.cosh(thd_v) - t*mp.sinh(thd_v)) + mp.mpf('0.2')))
    for t0 in (-0.8, 0.1, 0.9):
        # (phi, psi) = (soliton, 2 pi/beta): absorbed from the left if theta_d > 0, emitted to the left if theta_d < 0
        worst = max(worst, cond_res(kink_, lambda x, t: 2*mp.pi/bsg, t0, sgm))
        # (phi, psi) = (0, soliton): emitted to the right if theta_d > 0, absorbed from the right if theta_d < 0
        worst = max(worst, cond_res(lambda x, t: 0, kink_, t0, sgm))
report(f'theta = theta_d: (soliton, 2 pi/beta) and (0, soliton) both satisfy the defect conditions for theta_d = 0.4 and -0.5: '
       f'absorption and emission, direction fixed by sign theta_d [res {float(worst):.0e}]', worst < 1e-12)


def two_sol(x, t, ths, xis, zs=(1, 1)):
    E = [z*mp.e**(m0n*(x*mp.cosh(a) - t*mp.sinh(a)) + c) for a, c, z in zip(ths, xis, zs)]
    A = (mp.cosh(ths[0] - ths[1]) - 1)/(mp.cosh(ths[0] - ths[1]) + 1)
    t0_ = 1 + E[0] + E[1] + A*E[0]*E[1]; t1_ = 1 - E[0] - E[1] + A*E[0]*E[1]     # tau_j, omega = -1
    return 2j/bsg*mp.log(t1_/t0_)


ths, xis = (mp.mpf('0.9'), mp.mpf('-0.3')), (mp.mpc('0.2', '1.5707963267948966'), mp.mpc('-0.4', '1.5707963267948966'))
thd = mp.mpf('0.4'); sig = mp.e**(-thd)
worst = 0
for typ in (1, -1):
    zs = tuple(zclosed[typ](a, thd) for a in ths)
    for t0 in (-0.6, 0.0, 0.5):
        worst = max(worst, cond_res(lambda x, t: two_sol(x, t, ths, xis), lambda x, t: two_sol(x, t, ths, xis, zs), t0, sig, typ))
report(f'two solitons (theta = 0.9 > theta_d, -0.3 < theta_d) through the defect: each E_k -> z(theta_k) E_k, types I_-, I_+ '
       f'[res {float(worst):.0e}]', worst < 1e-12)

# ---------------------------------------------------------------- 10. half a soliton (BCZ 2005 sec. 2), pre-registered test
worst = 0
for typ in (1, -1):
    for th_ in (mp.mpf('-0.3'), mp.mpf('1.1')):
        zodd = 1/zclosed[typ](th_, thd)                       # sigma -> -sigma: z_odd = 1/z_even
        xi = mp.mpc('0.2', '0.9')
        for t0 in (-0.5, 0.4):
            worst = max(worst, cond_res(lambda x, t: sol(x, t, th_, xi), lambda x, t: 2*mp.pi/bsg + sol(x, t, th_, xi, zodd), t0, sig, typ))
okodd = worst < 1e-12 and abs(1/zclosed[1](mp.mpf('0.7'), thd) + mp.tanh((mp.mpf('0.7') - thd)/2)) < 1e-25
A_ = lambda D: (mp.cosh(D) - 1)/(mp.cosh(D) + 1)
devA = max(abs((1/zclosed[1](a, thd))**2 - A_(a - thd)) + abs(zclosed[1](a, thd)**2 - A_(a - thd - 1j*mp.pi))
           for a in (mp.mpf('-0.6'), mp.mpf('0.1'), mp.mpf('1.3')))
report(f'odd sector (jump 2 pi/beta): z_odd = 1/z_even = -tanh((theta-theta_d)/2) [res {float(worst):.0e}]; z_odd^2 = A_11(theta - theta_d), '
       f'z_even^2 = A_11(theta - theta_d - i pi) = 1/A_11 [{float(devA):.0e}]: two odd defects shift a soliton as a soliton of rapidity '
       'theta_d does, two even ones as one of rapidity theta_d + i pi', okodd and devA < 1e-20)
Mh = float(Mcl)/2
ok = True
for jump, rap in ((0.0, float(thd) + 1j*np.pi), (float(2*mp.pi/bsg), float(thd))):
    B, _, _, Pd = defect_potential(np.array([0.0]), np.array([jump]), float(sig), 1, beta=float(BETA), m=float(M_), coupling='imag')
    E_ = B.real - Mh*np.cosh(float(thd)); P_ = Pd.real - Mh*np.sinh(float(thd))       # BCZ's additive constants
    ok &= np.isclose(E_, (Mh*np.cosh(rap)).real) and np.isclose(P_, (Mh*np.sinh(rap)).real)
report('with BCZ\'s additive constant the defect carries (E, P) = (M/2)(cosh, sinh)(theta_d + i pi) without a jump and '
       '(M/2)(cosh, sinh)(theta_d) with a 2 pi/beta jump: the mass shell of half a soliton', ok)
f12 = lambda D: (mp.sinh(D) + 1j)/(mp.sinh(D) - 1j)                 # f_{1/2}: classical limit of S_{s B_1}
Tc_ = lambda D: -1j*mp.sinh(D/2 - 1j*mp.pi/4)/mp.sinh(D/2 + 1j*mp.pi/4)
devB = max(abs(Tc_(D - 1j*mp.pi)**2 - f12(D)) + abs(Tc_(D)**2 - f12(D - 1j*mp.pi)) for D in (mp.mpf('0.4'), mp.mpf('-1.2')))
devB += abs(Tc_(mp.mpf('0.4') - 1j*mp.pi) - mp.conj(Tc_(mp.mpf('0.4'))))
report(f'linear waves: T_odd = T_cl(theta_d -> theta_d + i pi) = conj T_cl and T_odd^2 = f_(1/2)(theta - theta_d), the classical '
       f'soliton-B_1 amplitude [{float(devB):.0e}]', devB < 1e-20)

# ---------------------------------------------------------------- 11. reality: Theta iff sigma real
rng = np.random.default_rng(7)
devs = {}
for sname, sgm in (('real +', 0.63), ('real -', -1.7), ('complex', 0.6*np.exp(0.4j))):
    d = 0
    for _ in range(4):
        p, sv = rng.normal(size=1) + 1j*rng.normal(size=1), rng.normal(size=1) + 1j*rng.normal(size=1)
        B1_, gp1, gq1, _ = defect_potential(p, sv, sgm, 1, beta=0.8, m=1.3, coupling='imag')
        B2_, gp2, gq2, _ = defect_potential(-np.conj(p), -np.conj(sv), sgm, 1, beta=0.8, m=1.3, coupling='imag')
        d = max(d, abs(np.conj(B1_) - B2_), abs(np.conj(gp1[0]) + gp2[0]), abs(np.conj(gq1[0]) + gq2[0]))
    devs[sname] = d
# Theta maps the soliton through the defect to a solution: residual of (-conj phi, -conj psi) with sigma
th_, xi = mp.mpf('0.7'), mp.mpc('0.3', '0.8')
res_t = {}
for sname, sgm in (('real', mp.mpf('0.63')), ('complex', 0.63*mp.e**(0.4j))):
    zz = (mp.e**(-th_) + sgm)/(mp.e**(-th_) - sgm)
    fphi = lambda x, t: -mp.conj(sol(x, t, th_, xi)); fpsi = lambda x, t: -mp.conj(sol(x, t, th_, xi, zz))
    p = fphi(0, 0.2); sv = fpsi(0, 0.2)
    px = mp.diff(lambda x: fphi(x, 0.2), 0); sx = mp.diff(lambda x: fpsi(x, 0.2), 0)
    pt = mp.diff(lambda t: fphi(0, t), 0.2); st = mp.diff(lambda t: fpsi(0, t), 0.2)
    res_t[sname] = defect_condition_residual(c_(p), c_(sv), c_(px), c_(sx), c_(pt), c_(st), complex(sgm), 1,
                                             beta=float(BETA), m=float(M_), coupling='imag')
report(f'Theta (phi, psi) -> (-conj phi, -conj psi): conj B_d(phi, psi) = B_d(Theta(phi, psi)) for sigma = 0.63, -1.7 '
       f'[{max(devs["real +"], devs["real -"]):.0e}], not for complex sigma [{devs["complex"]:.2f}]; Theta maps the soliton solution '
       f'to a solution iff sigma is real [res {res_t["real"]:.0e} vs {res_t["complex"]:.2f}]',
       max(devs['real +'], devs['real -']) < 1e-12 and devs['complex'] > 1e-2 and res_t['real'] < 1e-12 and res_t['complex'] > 1e-2)

# ---------------------------------------------------------------- 12. fusing the defect onto a Dirichlet boundary
phD = sp.symbols('phi_D', real=True)
C0 = sp.cosh(thd_ + sp.I*bS*phD/2); C1 = sp.cosh(thd_ - sp.I*bS*phD/2)
Bdd = Bd.subs({ps_: phD, sg_: sp.exp(-thd_)})
gz = -(2*mS/bS**2)*(C0*sp.exp(-sp.I*bS*ph_/2) + C1*sp.exp(sp.I*bS*ph_/2))
ok = sp.simplify(sp.diff((Bdd - gz).rewrite(sp.exp), ph_)) == 0
report('psi = phi_D at the defect turns B_d into the Ghoshal-Zamolodchikov boundary potential with C_0 = cosh(theta_d + i beta phi_D/2), '
       'C_1 = cosh(theta_d - i beta phi_D/2) (eq-bsg-c-dictionary; Bajnok-Simon (4); Habibullin via CZ 1202.6016)', ok)

# ---------------------------------------------------------------- 13. BCZ 2003 spin-3 conditions (2.6)-(2.9)
lam3 = b_*r2/2                                                    # BCZ lambda = beta/sqrt2 (book beta) = beta_SG/2
B, _ = B_sym(b_)
B1d, B2d = sp.diff(B, x_), sp.diff(B, y_)
ok = all(sp.simplify(e) == 0 for e in (
    2*lam3**2*(B1d**2 - B2d**2) - sp.diff(V_sym(b_, x_), x_, 2) + sp.diff(V_sym(b_, y_), y_, 2),
    sp.diff(B, y_, 2, x_, 1) - lam3**2*B1d, sp.diff(B, x_, 2, y_, 1) - lam3**2*B2d,
    sp.diff(B, y_, 3) - lam3**2*B2d, sp.diff(B, x_, 3) - lam3**2*B1d))
report('the sinh-Gordon B_d satisfies BCZ 2003 (2.6)-(2.9) (modified spin-3 charge) with lambda_BCZ = beta_SG/2', ok)
