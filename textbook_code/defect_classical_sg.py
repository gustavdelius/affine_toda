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
Runtime: about 15 seconds.
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
