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
    d(delta_s - delta_sbar)/d phi_0 = 2 d eta/d phi_0.  The amplitude is e^{iS}, and dS/dphi_0 = int dt dL/dphi_0 =
    -int dt phi_x(0,t), since dL/dphi_0 = (beta M_0/2) sin(beta(phi(0)-phi_0)/2) = -phi_x(0) by eq-bsg-boundary.
    A leapfrog simulation of a soliton and an antisoliton reflecting off the Dirichlet boundary (u_0 = 0.9, 2.0;
    v = 0.5, 0.8) gives int dt [u_x^s(0,t) - u_x^sbar(0,t)] = +8 pi, independent of u_0 and v.  So
    d(delta_s - delta_sbar)/d phi_0 = -8 pi/beta = -(lambda+1) beta: eta = -(lambda+1) beta phi_0/2.  GZ's sign would
    need -8 pi.  The convention e^{+iS} is checked on the same runs, blind to the sign of eta: the common phase
    2 arg sigma(eta, -i theta) depends on |eta| only, and its eta-derivative at lambda = 400 equals -(I_s + I_sbar)/8pi.
Runtime: about 10 seconds.
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
