"""Independent Rule-1 checks of briefs/referee_bsg_dirichlet_sign_report.md (integration pass).

Written from the book's definitions (ch. 30) and the report's formulas taken as claims under test.
Imports neither the book's scripts nor the referee's.  Units m0 = 1, u = beta phi, u_tt = u_xx - sin u.
"""
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.integrate import quad

def report(name, ok, extra=''):
    print(('PASS' if ok else 'FAIL'), name, extra)

# ---------------------------------------------------------------- item 1: static states, charges, energies
ok = True
for a in (0.6, 1.4, 2.7):
    x = np.linspace(-60, 0, 600001)
    ug = 4*np.arctan(np.tan(a/4)*np.exp(x))
    ue = 2*np.pi - 4*np.arctan(1/np.tan(a/4)*np.exp(x))
    for w in (ug, ue):
        ok &= abs(w[-1] - a) < 1e-12
    # E = int u_x^2 dx for static BPS-type solutions (u_x^2/2 = 1 - cos u)
    E = lambda w: np.trapz(0.5*np.gradient(w, x, edge_order=2)**2 + 1 - np.cos(w), x)
    ok &= abs(E(ug) - 4*(1 - np.cos(a/2))) < 1e-6 and abs(E(ue) - 4*(1 + np.cos(a/2))) < 1e-6
    # charges (u(0) - u(-oo))/2pi
    ok &= abs((ug[-1] - ug[0])/(2*np.pi) - a/(2*np.pi)) < 1e-9 and abs((ue[-1] - ue[0])/(2*np.pi) - (a/(2*np.pi) - 1)) < 1e-9
report('item 1: static states, charges a/2pi and a/2pi-1, gap 8cos(a/2) (units m0/beta^2)', ok)

# direct-channel kinematics of the pole at theta = i pi/2: energy 0, momentum iM
th = 1j*np.pi/2
report('item 1: at theta = i pi/2, cosh theta = 0 (zero energy), sinh theta = i (momentum iM)',
       abs(np.cosh(th)) < 1e-15 and abs(np.sinh(th) - 1j) < 1e-15)

# ---------------------------------------------------------------- item 2: exact three-soliton Dirichlet scattering
def tau_parts(x, t, theta, a, s):
    g, h, v, w, r = np.cosh(theta), np.sinh(theta), np.tanh(theta), np.tanh(theta/2)**2, np.tan(a/4)
    D = np.sqrt((w + r*r)/(v*v*w*(1 + r*r*w)))
    # exponentials and their x,t derivatives: z1 = sD e^{gx-ht}, z2 = sD e^{gx+ht}, z3 = (r/w) e^x
    z1, z2, z3 = s*D*np.exp(g*x - h*t), s*D*np.exp(g*x + h*t), (r/w)*np.exp(x)
    # tau as sum of monomials c * z1^i z2^j z3^k; derivatives act as multiplication by (k1 i + k2 j + k3 k)
    mons = [(1, 0, 0, 0), (-v*v, 1, 1, 0), (-w, 1, 0, 1), (-w, 0, 1, 1),
            (1j, 1, 0, 0), (1j, 0, 1, 0), (1j, 0, 0, 1), (-1j*v*v*w*w, 1, 1, 1)]
    kx, kt = (g, g, 1.0), (-h, h, 0.0)
    T = Tx = Txx = Tt = Ttt = 0
    for c, i, j, k in mons:
        m = c * z1**i * z2**j * z3**k
        px = kx[0]*i + kx[1]*j + kx[2]*k
        pt = kt[0]*i + kt[1]*j + kt[2]*k
        T, Tx, Txx, Tt, Ttt = T + m, Tx + px*m, Txx + px*px*m, Tt + pt*m, Ttt + pt*pt*m
    return T, Tx, Txx, Tt, Ttt

maxres = maxbv = maxslope = maxI = 0.0
ok_charge = True
for theta in (0.4, 0.8, 1.3):
    for a in (0.6, 1.4, -0.9):
        h = np.sinh(theta)
        k = np.sin(a/2)/np.sqrt(h*h + np.sin(a/2)**2)
        Is = {}
        for s in (1, -1):
            # PDE residual: u = 4 Im log tau (+const)
            X, Tm = np.meshgrid(np.linspace(-8, 0, 41), np.linspace(-6, 6, 41))
            T, Tx, Txx, Tt, Ttt = tau_parts(X, Tm, theta, a, s)
            uxx = 4*np.imag(Txx/T - (Tx/T)**2)
            utt = 4*np.imag(Ttt/T - (Tt/T)**2)
            sinu = np.imag(T**4)/np.abs(T)**4
            maxres = max(maxres, np.max(np.abs(utt - uxx + sinu)))
            # boundary value mod 2pi
            tt = np.linspace(-10, 10, 201)
            T0 = tau_parts(0.0, tt, theta, a, s)
            u0 = 4*np.angle(T0[0]) - 2*np.pi*s
            maxbv = max(maxbv, np.max(np.abs(np.angle(np.exp(1j*(u0 - a))))))
            # closed-form slope
            ux0 = 4*np.imag(T0[1]/T0[0])
            slope = 4*h*np.sqrt(1 - k*k)/(s*np.cosh(h*tt) + k)
            maxslope = max(maxslope, np.max(np.abs(ux0 - 2*np.sin(a/2) - slope)))
            # flux integral by quadrature of the exact field (not of the closed-form slope)
            f = lambda t_: float(4*np.imag((lambda P: P[1]/P[0])(tau_parts(0.0, t_, theta, a, s)))) - 2*np.sin(a/2)
            Tq = 70/h; I = quad(lambda t_: f(float(t_)), -Tq, Tq, points=[0.0], limit=400, epsabs=1e-13, epsrel=1e-13)[0]
            Is[s] = float(I)
            maxI = max(maxI, abs(Is[s] - (4*np.pi*s - 8*np.arctan(np.sin(a/2)/h))))
            # charge of incoming particle: field far left vs field between particle and wall, at early time
            if a > 0:
                xs = np.linspace(-80, 0, 80001)
                P = tau_parts(xs, -40.0, theta, a, s)[0]
                uu = 4*np.unwrap(np.angle(P))
                uu += a - uu[-1]                      # fix the branch by the boundary value
                # wall region (x ~ -3) should be the ground state profile, far left = 2pi*(-s) relative
                ug_m3 = 4*np.arctan(np.tan(a/4)*np.exp(-3.0))
                i3 = np.argmin(abs(xs + 3))
                ok_charge &= abs(uu[i3] - ug_m3) < 1e-4 and abs(uu[0] - (-2*np.pi*s)) < 1e-6
        maxI = max(maxI, abs(Is[1] - Is[-1] - 8*np.pi))
report('item 2: exact solution solves sine-Gordon, u(0,t) = a mod 2pi', maxres < 1e-9 and maxbv < 1e-12,
       f'[res {maxres:.1e}, bv {maxbv:.1e}]')
report('item 2: incoming particle has charge s and the wall carries the ground state before the collision', ok_charge)
report('item 2: u_x(0,t) - 2 sin(a/2) = 4h sqrt(1-k^2)/(s cosh ht + k)', maxslope < 1e-10, f'[{maxslope:.1e}]')
report('item 2: I_s = 4 pi s - 8 arctan(sin(a/2)/sinh theta), I_+ - I_- = 8 pi', maxI < 1e-8, f'[{maxI:.1e}]')

# e^{+iS} convention: d_eta arg sigma(eta, -i theta) at eta = (lambda+1) a/2 -> (2/pi) arctan(sin(a/2)/sinh theta)
mp.mp.dps = 25
def log_sigma(x, u, lam):
    """sigma(x,u) of ch. 30 (BPTT form), written from the book's formula."""
    x, u, L = mp.mpmathify(x), mp.mpmathify(u), mp.mpf(lam)
    def term(l, uu):
        v, c = L*uu/mp.pi, x/mp.pi
        return (mp.loggamma(0.5 + c + (2*l - 1)*L - v) + mp.loggamma(0.5 - c + (2*l - 1)*L - v)
                - mp.loggamma(0.5 - c + (2*l - 2)*L - v) - mp.loggamma(0.5 + c + 2*l*L - v))
    return mp.log(mp.cos(x)/mp.cos(x + L*u)) + mp.nsum(lambda l: term(l, u) - term(l, -u), [1, mp.inf])

theta, a = 0.8, 1.4
target = 2/np.pi*np.arctan(np.sin(a/2)/np.sinh(theta))
vals = []
for lam in (200, 400, 800):
    eta, hh = (lam + 1)*a/2, mp.mpf('1e-5')
    d = (mp.im(log_sigma(eta + hh, -1j*theta, lam)) - mp.im(log_sigma(eta - hh, -1j*theta, lam)))/(2*hh)
    vals.append(float(d))
extrap = 2*vals[2] - vals[1]          # 1/lambda Richardson
report('item 2: d_eta arg sigma -> (2/pi) arctan(sin(a/2)/sinh theta) = -(I_+ + I_-)/8pi (sign fails for e^{-iS})',
       abs(extrap - target) < 2e-4, f'[lam=200,400,800: {vals}; extrap {extrap:.6f}; target {target:.6f}]')

# ---------------------------------------------------------------- item 3: Taylor coefficients and branch
ok = True
for lam in (np.sqrt(2), np.e, 2 + np.sqrt(2)):
    f = lambda z: -mp.cos((lam + 1)*mp.acos(z))/mp.sin(mp.pi*lam)
    ok &= abs(f(0) - 1/(2*np.cos(np.pi*lam/2))) < 1e-14
    ok &= abs(mp.diff(f, 0) + (lam + 1)/(2*np.sin(np.pi*lam/2))) < 1e-12
    # another branch of arccos (shift by 2 pi) gives a different value for non-integer lambda+1
    ok &= abs(-mp.cos((lam + 1)*(mp.acos(0.1) + 2*mp.pi))/mp.sin(mp.pi*lam) - f(0.1)) > 1e-3
report('item 3: f(0) = eps_*, f\'(0) = -(lambda+1)/(2 sin(pi lambda/2)); value depends on the arccos branch', ok)

# ---------------------------------------------------------------- item 4: Neumann representative, Dirichlet limit, C-conjugation
ok = True
lam = 2 + np.sqrt(2)
for a in (0.5, 1.9):
    # solve the corrected UV-IR relation for rho small and large with vartheta >= 0, continuing from Dirichlet
    for rho in (1e-4, 1e4):
        z = mp.acos(rho*mp.exp(-1j*a/2))           # (eta - i vth)/(lam+1) = +-z
        for sgn in (1, -1):
            A, B = float(mp.re(sgn*z)), float(-mp.im(sgn*z))
            if B >= 0:
                break
        ok &= abs(np.cos(A)*np.cosh(B) - rho*np.cos(a/2)) < 1e-8*max(1, rho) and abs(np.sin(A)*np.sinh(B) + rho*np.sin(a/2)) < 1e-8*max(1, rho)
        if rho < 1:
            ok &= abs(A + np.pi/2) < 1e-3 and B < 1e-3        # negative Neumann representative
        else:
            ok &= abs(A + a/2) < 1e-6                        # eta = -(lam+1) a/2
report('item 4: for phi_0 > 0, vartheta >= 0 the corrected relation tends to eta = -(pi/2)(lambda+1) at Neumann '
       'and to eta = -(lambda+1) beta phi_0/2 at Dirichlet', ok)
# K(-eta, vth) = C K(eta, vth) C from the P, Q formulas (sigma even in eta, cos even)
ok = True
for (eta, vth, uu) in ((1.3, 0.7, 0.4 + 0.2j), (2.9, 1.1, -0.3 + 0.9j)):
    P = lambda e, sg: np.cos(lam*uu)*np.cos(e)*np.cosh(vth) - sg*np.sin(lam*uu)*np.sin(e)*np.sinh(vth)
    ok &= abs(P(-eta, 1) - P(eta, -1)) < 1e-14 and abs(P(-eta, -1) - P(eta, 1)) < 1e-14
s_even = max(abs(mp.exp(log_sigma(-x_, u_, lam) - log_sigma(x_, u_, lam)) - 1)
             for x_ in (1.3, 2.9) for u_ in (mp.mpc(0.4, 0.2), mp.mpc(-0.3, 0.9)))
report('item 4: K(-eta, vartheta) = C K(eta, vartheta) C (P+ <-> P-, Q and sigma even)', ok and s_even < 1e-15,
       f'[sigma evenness {float(s_even):.0e}]')

# ---------------------------------------------------------------- item 5: Mattsson-Dorey arithmetic
b, ph0 = sp.symbols('beta varphi0', positive=True)
q_exc = b/(2*sp.pi)*(ph0 - 2*sp.pi/b)
report('item 5: MD charge (beta/2pi)[phi(0) - phi(-oo)] with phi(-oo) = 2pi/beta is beta phi0/2pi - 1',
       sp.simplify(q_exc - (b*ph0/(2*sp.pi) - 1)) == 0)

# ---------------------------------------------------------------- item 7: ordering of C0, C1 in defects_l9_fusion (e)
m0, sig, psi0 = sp.symbols('m0 sigma psi0', positive=True)
Fs, thd = sp.symbols('Phi theta_d', real=True)
Bdef = -2*m0*(sig*(sp.cos((Fs + psi0)/2) - 1) + (sp.cos((Fs - psi0)/2) - 1)/sig)
Bexp = sp.expand(Bdef.subs(sig, sp.exp(-thd)).rewrite(sp.exp))
cp = sp.simplify(Bexp.coeff(sp.exp(sp.I*Fs/2)))      # coefficient of e^{+i Phi/2}
cm = sp.simplify(Bexp.coeff(sp.exp(-sp.I*Fs/2)))     # coefficient of e^{-i Phi/2}
# book, eq-bsg-boundary and eq-bsg-c-dictionary at beta = 1: -M0 cos((Phi - Phi0)/2) has e^{+i Phi/2}-coefficient
# -(M0/2) e^{-i Phi0/2} = -2 m0 C1 and e^{-i Phi/2}-coefficient -2 m0 C0.  Verify that normalisation first.
M0, P0 = sp.symbols('M0 Phi0', positive=True)
Bgz = sp.expand((-M0*sp.cos((Fs - P0)/2)).rewrite(sp.exp))
okN = (sp.simplify(Bgz.coeff(sp.exp(sp.I*Fs/2)) + 2*m0*(M0/(4*m0))*sp.exp(-sp.I*P0/2)) == 0 and
       sp.simplify(Bgz.coeff(sp.exp(-sp.I*Fs/2)) + 2*m0*(M0/(4*m0))*sp.exp(sp.I*P0/2)) == 0)
C1 = sp.simplify(-cp/(2*m0)); C0 = sp.simplify(-cm/(2*m0))
okC = (sp.simplify((C0 - sp.cosh(thd + sp.I*psi0/2)).rewrite(sp.exp)) == 0 and
       sp.simplify((C1 - sp.cosh(thd - sp.I*psi0/2)).rewrite(sp.exp)) == 0)
print('   C0 =', C0, '  C1 =', C1)
report('item 7: frozen defect potential has C0 = cosh(theta_d + i psi0/2), C1 = cosh(theta_d - i psi0/2): "+-" ordering',
       okN and okC)
