"""The sign of phi_0 in the Dirichlet relation of boundary sine-Gordon theory (ch. 30, eq-bsg-dirichlet, sec-bsg-parameters).

Book conventions: bulk on x <= 0, phi(0,t) = phi_0, topological charge T = (beta/2pi)(phi(+oo) - phi(-oo)), soliton s with
T = +1, reflection matrix eq-bsg-k-full in the basis (s, sbar), so P^+ = K_ss.  On the half-line the Dirichlet boundary
conserves (beta/2pi)(phi(0) - phi(-oo)).  Classical units: u = beta phi, m_0 = 1, u_tt = u_xx - sin u, soliton mass
M = 8 (in units of m_0/beta^2).  GZ = Ghoshal-Zamolodchikov hep-th/9306002, whose (5.27) reads xi = 4 pi phi_0/beta,
i.e. eta = +(lambda+1) beta phi_0/2.

Checks:
 1. Static Dirichlet states for 0 < u_0 < pi: u = 4 arctan e^{x-a} (phi(-oo) = 0) and u = 2pi - 4 arctan e^{x-abar}
    (phi(-oo) = 2pi/beta) solve the static equation with u(0) = u_0; by quadrature their energies differ by
    M cos(u_0/2), and their charges are u_0/2pi and u_0/2pi - 1.  So the excited state has one unit of charge less:
    the boundary absorbs an antisoliton to reach it.
 2. In eq-bsg-k-full at vartheta -> oo, with eta > 0, the pole nu_0 = eta/lambda - pi/(2 lambda) of sigma(eta,u) is a
    simple pole of P^+ = K_ss and is cancelled in P^- by the zero of cos(lambda u - eta).  sigma(-x,u) = sigma(x,u), so
    K(-eta, vartheta) = C K(eta, vartheta) C: for eta < 0 the pole sits in P^- = K_sbar sbar.
 3. As lambda -> oo with |eta| = (lambda+1) u_0/2, the energy M cos nu_0 of the bound state tends to M cos(u_0/2), the
    excited static state of check 1.  With checks 1 and 2: the bound state is reached by an antisoliton, so for
    phi_0 > 0 the pole must be in P^-, i.e. eta < 0: eta = -(lambda+1) beta phi_0/2, not GZ's sign.
 4. Independent check of the same sign, by classical scattering.  For real theta, P^+/P^- = cos(eta - i lambda theta)/
    cos(eta + i lambda theta) -> e^{2 i eta} for lambda theta >> 1 (the scalar factor is common), so semiclassically
    d(delta_s - delta_sbar)/d phi_0 = 2 d eta/d phi_0.  The amplitude is e^{iS_sc}, with S_sc the action with the static
    ground state phi_g and the free asymptotic motion subtracted, and dS_sc/dphi_0 = -int dt [phi_x(0,t) - phi_x^g(0)],
    since dL/dphi_0 = (beta M_0/2) sin(beta(phi(0)-phi_0)/2) = -phi_x(0) by eq-bsg-boundary; u_x^g(0) = 2 sin(u_0/2).
    A leapfrog simulation of a soliton and an antisoliton reflecting off the Dirichlet boundary (u_0 = 0.9, 2.0;
    v = 0.5, 0.8) gives int dt [u_x^s(0,t) - u_x^sbar(0,t)] = +8 pi, independent of u_0 and v.  So
    d(delta_s - delta_sbar)/d phi_0 = -8 pi/beta = -(lambda+1) beta: eta = -(lambda+1) beta phi_0/2.  GZ's sign would
    need -8 pi.  The convention e^{+iS} is checked on the same runs, blind to the sign of eta: the common phase
    2 arg sigma(eta, -i theta) depends on |eta| only, and its eta-derivative at lambda = 400 equals -(I_s + I_sbar)/8pi.
 5. eq-bsg-dirichlet-flux.  The exact three-soliton solution u = 4 arg tau - 2 pi s, with
    tau = 1 - v^2 z1 z2 - w z3 (z1 + z2) + i [z1 + z2 + z3 - v^2 w^2 z1 z2 z3],
    z1,2 = s D e^{x cosh theta -+ t sinh theta}, z3 = (r/w) e^x, v = tanh theta, w = tanh^2(theta/2), r = tan(u_0/4),
    D^2 = (w + r^2)/(v^2 w (1 + r^2 w)), solves the sine-Gordon equation, has u(0,t) = u_0 mod 2pi, an incoming
    particle of charge s and the ground state at the wall before the collision.  Its flux is
    I_s = int dt [u_x(0,t) - 2 sin(u_0/2)] = 4 pi s - 8 arctan(sin(u_0/2)/sinh theta), by quadrature; the simulation of
    check 4 reproduces it, and -(I_s + I_sbar)/8pi = (2/pi) arctan(sin(u_0/2)/sinh theta) is the lambda -> oo limit of
    d_eta arg sigma(eta, -i theta) at eta = (lambda+1) u_0/2 (Richardson extrapolation from lambda = 400, 800).
Runtime: about 15 seconds.
"""
import numpy as np
import mpmath as mp
from boundary_common import report

mp.mp.dps = 20


# ---------------------------------------------------------------- 1. static Dirichlet states
def static_states(u0, L=40.0, n=400001):
    x = np.linspace(-L, 0, n)
    a = -np.log(np.tan(u0/4))                       # 4 arctan e^{-a} = u_0
    abar = -np.log(np.tan((2*np.pi - u0)/4))        # 2 pi - 4 arctan e^{-abar} = u_0
    return x, 4*np.arctan(np.exp(x - a)), 2*np.pi - 4*np.arctan(np.exp(x - abar))


def energy(x, u):
    ux = np.gradient(u, x)
    f = 0.5*ux**2 + 1 - np.cos(u)
    return np.sum((f[1:] + f[:-1])*np.diff(x))/2


ok = True
for u0 in (0.9, 2.0):
    x, ug, ue = static_states(u0)
    dx = x[1] - x[0]
    for w in (ug, ue):
        res = (w[2:] - 2*w[1:-1] + w[:-2])/dx**2 - np.sin(w[1:-1])
        ok &= np.max(np.abs(res)) < 1e-4 and abs(w[-1] - u0) < 1e-12
    dE = energy(x, ue) - energy(x, ug)
    ok &= abs(dE - 8*np.cos(u0/2)) < 1e-6
    Qg, Qe = (ug[-1] - ug[0])/(2*np.pi), (ue[-1] - ue[0])/(2*np.pi)
    ok &= abs(Qg - u0/(2*np.pi)) < 1e-12 and abs(Qe - (u0/(2*np.pi) - 1)) < 1e-12
report('static Dirichlet states, 0 < u_0 < pi (u_0 = 0.9, 2.0): phi(-oo) = 0 and 2pi/beta, energy difference M cos(u_0/2), '
       'charges u_0/2pi and u_0/2pi - 1: the excited state is reached by absorbing an antisoliton', ok)


# ---------------------------------------------------------------- 2. which amplitude has the pole nu_0
def lg(z):
    return mp.loggamma(z)


def log_sigma(x, u, lam):
    """BPTT form of GZ's sigma (as in boundary_sine_gordon.py)."""
    x, u, L = mp.mpmathify(x), mp.mpmathify(u), mp.mpf(lam)
    def t(l, uu):
        v = L*uu/mp.pi; a = x/mp.pi
        return (lg(0.5 + a + (2*l-1)*L - v) + lg(0.5 - a + (2*l-1)*L - v)
                - lg(0.5 - a + (2*l-2)*L - v) - lg(0.5 + a + 2*l*L - v))
    return mp.log(mp.cos(x)/mp.cos(x + L*u)) + mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])


def P_dirichlet(u, eta, lam, sign):
    """P^{+-} at vartheta -> oo without R_0 (common to both): sigma(eta,u) cos(lambda u +- eta)/cos(eta)."""
    return mp.exp(log_sigma(eta, u, lam))*mp.cos(lam*u + sign*eta)/mp.cos(eta)


def order(f, u0, d=mp.mpf('1e-7')):
    return float(mp.log(abs(f(u0 + d))/abs(f(u0 + 2*d)))/mp.log(2))


lam, eta = mp.sqrt(11), mp.mpf('2.1')              # generic lambda (not rational)
nu0 = eta/lam - mp.pi/(2*lam)
o_p = order(lambda u: P_dirichlet(u, eta, lam, +1), nu0)
o_m = order(lambda u: P_dirichlet(u, eta, lam, -1), nu0)
dev_even = max(abs(mp.exp(log_sigma(-x_, u_, lam) - log_sigma(x_, u_, lam)) - 1)
               for x_ in (eta, mp.mpf('0.7'), 1j*mp.mpf('0.9')) for u_ in (mp.mpc('0.3', '0.2'), mp.mpc('-0.4', '1.1')))
dev_C = max(abs(P_dirichlet(u_, -eta, lam, +1) - P_dirichlet(u_, eta, lam, -1)) for u_ in (mp.mpc('0.3', '0.2'), mp.mpf('0.21')))
report(f'eq-bsg-k-full at vartheta -> oo, eta = 2.1 > 0, lambda = sqrt(11): at nu_0 P^+ has a simple pole [order {o_p:.3f}], '
       f'P^- none [{o_m:.3f}]; sigma(-x,u) = sigma(x,u) [{float(dev_even):.0e}], so P^+(-eta) = P^-(eta) [{float(dev_C):.0e}]',
       abs(o_p - 1) < 1e-3 and abs(o_m) < 1e-3 and dev_even < 1e-15 and dev_C < 1e-15)

# ---------------------------------------------------------------- 3. semiclassical energy of the bound state
ok = True
for u0 in (0.9, 2.0):
    for lam_ in (1e3, 1e6):
        nu0_ = (lam_ + 1)*u0/2/lam_ - np.pi/(2*lam_)
        ok &= abs(np.cos(nu0_) - np.cos(u0/2)) < 2/lam_
report('lambda -> oo at |eta| = (lambda+1) u_0/2: M cos nu_0 -> M cos(u_0/2), the excited static state of check 1, '
       'so the pole belongs to the antisoliton and eta = -(lambda+1) beta phi_0/2 for phi_0 > 0', ok)


# ---------------------------------------------------------------- 4. classical scattering: phi_0-derivative of the phase
def flux_integral(u0, v, kind, L=45.0, dx=0.02, x0=-16.0):
    """int dt [u_x(0,t) - u_x^bg(0)] for a soliton ('s') or antisoliton ('a') of speed v reflecting at x = 0."""
    N = int(round(L/dx)); x = -L + dx*np.arange(N + 1)
    dt = 0.5*dx
    ubg = 4*np.arctan(np.exp(x + np.log(np.tan(u0/4))))
    g = 1/np.sqrt(1 - v*v); z = g*(x - x0)
    if kind == 's':      # phi(-oo) = -2pi/beta, then the ground state: T = +1
        u = ubg + 4*np.arctan(np.exp(z)) - 2*np.pi; ut = -2*v*g/np.cosh(z)
    else:                # phi(-oo) = +2pi/beta: T = -1
        u = ubg + 2*np.pi - 4*np.arctan(np.exp(z)); ut = 2*v*g/np.cosh(z)
    def acc(w):
        r = np.zeros_like(w)
        r[1:-1] = (w[2:] - 2*w[1:-1] + w[:-2])/dx**2 - np.sin(w[1:-1])
        return r
    uleft = u[0]
    uprev = u - dt*ut + 0.5*dt*dt*acc(u); uprev[0] = uleft; uprev[-1] = u0
    uxbg = 2*np.sin(u0/2)
    I = 0.0
    for _ in range(int((2*abs(x0)/v + 40)/dt)):
        I += ((3*u[-1] - 4*u[-2] + u[-3])/(2*dx) - uxbg)*dt
        unew = 2*u - uprev + dt*dt*acc(u); unew[0] = uleft; unew[-1] = u0
        uprev, u = u, unew
    return I


pairs = [(u0, v) for u0 in (0.9, 2.0) for v in (0.5, 0.8)]
flux = {p: (flux_integral(*p, 's'), flux_integral(*p, 'a')) for p in pairs}
# convention check, blind to the sign of eta: the common phase 2 arg sigma(eta, -i theta) depends on |eta| only, and
# d/d phi_0 of it must equal -(1/beta)(I_s + I_sbar); in eta at large lambda: d_eta arg sigma = -(I_s + I_sbar)/(8 pi)
dev_conv = 0
for u0, v in ((0.9, 0.5), (2.0, 0.8)):
    lamc, h = 400, mp.mpf('1e-4')
    etac, thc = (lamc + 1)*u0/2, np.arctanh(v)
    d = (mp.im(log_sigma(etac + h, -1j*thc, lamc)) - mp.im(log_sigma(etac - h, -1j*thc, lamc)))/(2*h)
    dev_conv = max(dev_conv, abs(float(d) + sum(flux[(u0, v)])/(8*np.pi)))
report(f'semiclassical amplitude e^(+iS): d/d eta of arg sigma(eta, -i theta) at lambda = 400 equals '
       f'-(I_s + I_sbar)/(8 pi) from the simulation [{dev_conv:.1e}] (u_0, v) = (0.9, 0.5), (2.0, 0.8)', dev_conv < 3e-3)
lam_, th_, eta_ = 40.0, 0.5, 1.3
ratio = np.cos(eta_ - 1j*lam_*th_)/np.cos(eta_ + 1j*lam_*th_)
diffs = [flux[p][0] - flux[p][1] for p in pairs]
report(f'P^+/P^- -> e^(2 i eta) at lambda theta = 20 [{abs(ratio - np.exp(2j*eta_)):.0e}]; classical reflection off '
       f'phi(0) = phi_0: int dt [u_x^s(0,t) - u_x^sbar(0,t)] = {", ".join(f"{d:.4f}" for d in diffs)} = +8 pi '
       f'({8*np.pi:.4f}) for u_0 = 0.9, 2.0, v = 0.5, 0.8; so d(delta_s - delta_sbar)/d phi_0 = -8pi/beta and '
       f'eta = -(lambda+1) beta phi_0/2 (GZ sign would need -8 pi)',
       abs(ratio - np.exp(2j*eta_)) < 1e-12 and all(abs(d - 8*np.pi) < 1e-2 for d in diffs))


# ---------------------------------------------------------------- 5. the exact flux, eq-bsg-dirichlet-flux
def tau_derivs(x, t, theta, u0, s):
    """tau of the three-soliton Dirichlet reflection and its x, t derivatives up to second order."""
    g, hh, v, w, r = np.cosh(theta), np.sinh(theta), np.tanh(theta), np.tanh(theta/2)**2, np.tan(u0/4)
    D = np.sqrt((w + r*r)/(v*v*w*(1 + r*r*w)))
    z = (s*D*np.exp(g*x - hh*t), s*D*np.exp(g*x + hh*t), (r/w)*np.exp(x))
    kx, kt = (g, g, 1.0), (-hh, hh, 0.0)
    out = [0, 0, 0, 0, 0]                          # tau, tau_x, tau_xx, tau_t, tau_tt
    for c, n in ((1, (0, 0, 0)), (-v*v, (1, 1, 0)), (-w, (1, 0, 1)), (-w, (0, 1, 1)),
                 (1j, (1, 0, 0)), (1j, (0, 1, 0)), (1j, (0, 0, 1)), (-1j*v*v*w*w, (1, 1, 1))):
        m = c*z[0]**n[0]*z[1]**n[1]*z[2]**n[2]
        px, pt = np.dot(kx, n), np.dot(kt, n)
        for i, f in enumerate((1, px, px*px, pt, pt*pt)):
            out[i] = out[i] + f*m
    return out


def exact_flux(theta, u0, s):
    """int dt [u_x(0,t) - 2 sin(u0/2)] for u = 4 arg tau - 2 pi s, by quadrature."""
    from scipy.integrate import quad
    def f(t):
        T = tau_derivs(0.0, t, theta, u0, s)
        return 4*np.imag(T[1]/T[0]) - 2*np.sin(u0/2)
    T = 70/np.sinh(theta)
    return quad(f, -T, T, points=[0.0], limit=400, epsabs=1e-13, epsrel=1e-13)[0]


res = bv = dI = 0.0
ok_bg = True
for theta in (0.4, 0.8, 1.3):
    for u0 in (0.6, 1.4, -0.9):
        for s in (1, -1):
            X, T = np.meshgrid(np.linspace(-8, 0, 41), np.linspace(-6, 6, 41))
            t0, tx, txx, tt, ttt = tau_derivs(X, T, theta, u0, s)
            # u = 4 Im log tau: u_tt - u_xx + sin u, with sin u = sin(4 arg tau) = Im tau^4/|tau|^4
            res = max(res, np.max(np.abs(4*np.imag(ttt/t0 - (tt/t0)**2) - 4*np.imag(txx/t0 - (tx/t0)**2)
                                         + np.imag(t0**4)/np.abs(t0)**4)))
            tb = tau_derivs(0.0, np.linspace(-10, 10, 201), theta, u0, s)[0]
            bv = max(bv, np.max(np.abs(np.angle(np.exp(1j*(4*np.angle(tb) - u0))))))
            closed = 4*np.pi*s - 8*np.arctan(np.sin(u0/2)/np.sinh(theta))
            dI = max(dI, abs(exact_flux(theta, u0, s) - closed))
            if u0 > 0:      # before the collision: phi(-oo) = -2 pi s/beta, ground state at the wall
                xs = np.linspace(-80, 0, 80001)
                uu = 4*np.unwrap(np.angle(tau_derivs(xs, -40.0, theta, u0, s)[0]))
                uu += u0 - uu[-1]
                i3 = np.argmin(abs(xs + 3))
                ok_bg &= (abs(uu[0] + 2*np.pi*s) < 1e-6
                          and abs(uu[i3] - 4*np.arctan(np.tan(u0/4)*np.exp(-3.0))) < 1e-4)
# the simulation of check 4 against the closed form (v = tanh theta)
dsim = max(abs(flux[(u0, v)][i] - (4*np.pi*s - 8*np.arctan(np.sin(u0/2)/np.sinh(np.arctanh(v)))))
           for (u0, v) in pairs for i, s in ((0, 1), (1, -1)))
report(f'exact three-soliton reflection: sine-Gordon residual {res:.0e}, u(0,t) = u_0 mod 2pi [{bv:.0e}], incoming '
       f'charge s with the ground state at the wall; int dt [u_x(0,t) - 2 sin(u_0/2)] = 4 pi s - 8 arctan(sin(u_0/2)/'
       f'sinh theta) [{dI:.0e}] for theta = 0.4, 0.8, 1.3, u_0 = 0.6, 1.4, -0.9; the simulation agrees [{dsim:.1e}]',
       res < 1e-9 and bv < 1e-12 and ok_bg and dI < 1e-8 and dsim < 2e-2)

# the convention e^{+iS}: lambda -> oo limit of d_eta arg sigma against the closed-form sum of the fluxes
th5, u05 = 0.8, 1.4
target = 2/np.pi*np.arctan(np.sin(u05/2)/np.sinh(th5))        # = -(I_s + I_sbar)/8pi
dvals = []
for lam5 in (400, 800):
    eta5, h5 = (lam5 + 1)*u05/2, mp.mpf('1e-5')
    dvals.append(float((mp.im(log_sigma(eta5 + h5, -1j*th5, lam5))
                        - mp.im(log_sigma(eta5 - h5, -1j*th5, lam5)))/(2*h5)))
extrap = 2*dvals[1] - dvals[0]
report(f'e^(+iS) in closed form: d_eta arg sigma(eta, -i theta) at eta = (lambda+1) u_0/2, extrapolated to lambda -> oo, '
       f'= {extrap:.6f}; -(I_s + I_sbar)/8pi = (2/pi) arctan(sin(u_0/2)/sinh theta) = {target:.6f} (theta = 0.8, '
       f'u_0 = 1.4); e^(-iS) would need {-target:.6f}', abs(extrap - target) < 2e-4)
