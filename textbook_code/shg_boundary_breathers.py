"""Sinh-Gordon boundary breathers (Part VII, ch. 29, sec-shg-boundary-breathers; Corrigan-Delius hep-th/9909145 = CD).

Book normalization (eq-boundary-condition with C_0 = C_1 = eps, real coupling, m = 1, particle mass 2m = 2):
  phi_tt - phi_xx + (sqrt8/beta) sinh(sqrt2 beta phi) = 0 on x <= 0,
  phi_x(0) = (sqrt2/beta)(eps e^{-beta phi/sqrt2} - eps e^{beta phi/sqrt2}),
  energy  = int [phi_t^2/2 + phi_x^2/2 + (2/beta^2)(cosh sqrt2 beta phi - 1)] + (2/beta^2) eps (e^{-..} + e^{..} - 2)|_0.
Checks:
 1. CD (3.1)-(3.2) tau-functions solve the field equation and the boundary condition.
 2. Regularity: no zero of tau_0 tau_1 on x <= 0 iff -1 < eps < 0 and cos rho < -eps (CD (3.4)); a point outside
    is singular.
 3. Energy by quadrature equals CD (3.7), E = (8/beta^2)(-cos rho - eps).
 4. Bohr-Sommerfeld action: int_0^T dt int dx phi_t^2 = (8 pi/beta^2)(rho - pi(1-a)), eps = cos(pi a) (CD (3.9)).
 5. Bootstrap levels (CD (4.6)-(4.7)) with E(eps, beta) = 2a(1 - B/2) + (1 - 2N_0)B/2 (CD (5.29)) and
    m(beta) = (8/(pi B)) sin(pi B/4) (CD (5.30)) reproduce the WKB level differences (CD (5.27)) for all n,
    for any N_0; the Neumann limit a -> 1/2 gives |E| = 1 - B/2 (Ghoshal-Zamolodchikov) only for N_0 = 1/2.
 6. m(beta) is the continuation beta_SG^2 -> -beta_SG^2 of the sine-Gordon m_1 = 2M sin(xi/2), M = m_0/xi.
 7. CD (4.4)-(4.6): the excited reflection factors K_1, K_2 follow from K_0 by the boundary bootstrap.
Runtime: about 10 seconds.
"""
import numpy as np
import sympy as sp
from boundary_common import report

x, t = sp.symbols('x t', real=True)


def phi_expr(rho, eps, beta):
    c = sp.cos(rho)
    r = (eps + c)/(eps - c)
    tau = [1 + (-1)**j*2*sp.cos(2*t*sp.sin(rho))*sp.exp(2*x*c)*sp.sqrt(r)/sp.tan(rho) - sp.exp(4*x*c)*r
           for j in (0, 1)]
    return sp.sqrt(2)/beta*sp.log(tau[0]/tau[1]), tau


# ---------------------------------------------------------------- 1. field equation and boundary condition
rho, eps, beta = sp.Rational(23, 20), sp.Rational(-7, 10), sp.Rational(9, 10)
ph, tau = phi_expr(rho, eps, beta)
eq = sp.diff(ph, t, 2) - sp.diff(ph, x, 2) + sp.sqrt(8)/beta*sp.sinh(sp.sqrt(2)*beta*ph)
bc = sp.diff(ph, x) - sp.sqrt(2)/beta*(eps*sp.exp(-beta*ph/sp.sqrt(2)) - eps*sp.exp(beta*ph/sp.sqrt(2)))
fe = sp.lambdify((x, t), eq, 'mpmath'); fb = sp.lambdify(t, bc.subs(x, 0), 'mpmath')
import mpmath as mp
mp.mp.dps = 30
r1 = max(abs(fe(mp.mpf(xx), mp.mpf(tt))) for xx in ('-0.3', '-1.7') for tt in ('0.2', '1.9'))
r2 = max(abs(fb(mp.mpf(tt))) for tt in ('0.2', '1.1', '2.6'))
report(f'CD (3.2) tau-functions: field equation [{float(r1):.0e}] and boundary condition with C_0 = C_1 = eps '
       f'[{float(r2):.0e}] (book normalization)', r1 < 1e-20 and r2 < 1e-20)


# ---------------------------------------------------------------- 2. regularity region
def min_tau(rho, eps):
    c = np.cos(rho); r = (eps + c)/(eps - c)
    if r < 0:
        return 0.0
    xs = np.linspace(-12, 0, 2401)[:, None]; ts = np.linspace(0, np.pi/np.sin(rho), 401)[None, :]
    vals = [1 + s*2*np.cos(2*ts*np.sin(rho))*np.exp(2*xs*c)*np.sqrt(r)/np.tan(rho) - np.exp(4*xs*c)*r for s in (1, -1)]
    return min(np.abs(v).min() for v in vals)


inside = [(1.2, -0.6), (1.4, -0.3), (1.0, -0.9)]
outside = [(1.2, -0.2), (0.9, -0.5)]               # cos rho > -eps
report('regular on x <= 0 for -1 < eps < 0, cos rho < -eps (three points); singular for cos rho > -eps (two points)',
       all(min_tau(*p) > 1e-2 for p in inside) and all(min_tau(*p) < 1e-2 for p in outside))

# ---------------------------------------------------------------- 3./4. energy and Bohr-Sommerfeld action
rho_n, eps_n, beta_n = 1.15, -0.7, 0.9
phx = sp.lambdify((x, t), sp.diff(ph, x), 'numpy'); pht = sp.lambdify((x, t), sp.diff(ph, t), 'numpy')
phf = sp.lambdify((x, t), ph, 'numpy')
from scipy.integrate import quad
dens1 = lambda xx, tt: (0.5*pht(xx, tt)**2 + 0.5*phx(xx, tt)**2
                        + 2/beta_n**2*(np.cosh(np.sqrt(2)*beta_n*phf(xx, tt)) - 1))
Ens = []
for tt in (0.0, 0.7, 1.9):
    p0 = phf(0.0, tt)
    bnd = 2/beta_n**2*eps_n*(np.exp(-beta_n*p0/np.sqrt(2)) + np.exp(beta_n*p0/np.sqrt(2)) - 2)
    Ens.append(quad(lambda xx: float(dens1(xx, tt)), -60, 0, limit=400, epsabs=1e-13, epsrel=1e-13)[0] + bnd)
E_cd = 8/beta_n**2*(-np.cos(rho_n) - eps_n)
report(f'energy by quadrature {Ens[0]:.10f} (conserved: spread {max(Ens)-min(Ens):.0e}) = CD (3.7) '
       f'(8/beta^2)(-cos rho - eps) = {E_cd:.10f}', abs(Ens[0] - E_cd) < 1e-8 and max(Ens) - min(Ens) < 1e-8)
T = np.pi/np.sin(rho_n)
ts = np.linspace(0, T, 129)[:-1]
act = sum(quad(lambda xx: float(pht(xx, tt)**2), -60, 0, limit=400, epsabs=1e-13, epsrel=1e-13)[0] for tt in ts)*T/len(ts)
a_par = np.arccos(eps_n)/np.pi
act_cd = 8*np.pi/beta_n**2*(rho_n - np.pi*(1 - a_par))
report(f'Bohr-Sommerfeld action int dt dx phi_t^2 = {act:.10f} = CD (3.9) (8pi/beta^2)(rho - pi(1-a)) = {act_cd:.10f}',
       abs(act - act_cd) < 1e-8)

# ---------------------------------------------------------------- 5. bootstrap levels vs WKB levels
ok, okN = True, True
for beta2 in (0.8, 3.1):
    B = 2*beta2/(4*np.pi + beta2)                    # eq-b-of-beta
    m = 8/(np.pi*B)*np.sin(np.pi*B/4)
    for a in (0.6, 0.83):
        for N0 in (0.5, 0.3):
            E = 2*a*(1 - B/2) + (1 - 2*N0)*B/2
            for n in range(4):
                boot = m*np.cos(np.pi/2*(n*B - E + 1))                                        # CD (4.7)
                wkb = 8/(np.pi*B)*np.sin(np.pi*B/4)*np.cos(np.pi/2*(2*np.pi*B/beta2*(2*a - 1) - (n + N0)*B))  # (5.27)
                ok &= abs(boot - wkb) < 1e-12
    EN = 2*0.5*(1 - B/2)
    okN &= abs(EN - (1 - B/2)) < 1e-14 and abs(2*0.5*(1 - B/2) + (1 - 2*0.3)*B/2 - (1 - B/2)) > 1e-3
report('bootstrap level spacings m cos(pi(nB - E + 1)/2) = WKB spacings (5.27) for n = 0..3, with E = 2a(1-B/2) + '
       '(1-2N_0)B/2 and m = (8/pi B) sin(pi B/4), for any N_0', ok)
report('Neumann limit a -> 1/2: E = 1 - B/2 (Ghoshal-Zamolodchikov, up to E -> -E) requires N_0 = 1/2', okN)

# ---------------------------------------------------------------- 6. mass formula vs sine-Gordon breather
bSG2 = 2.3                                           # sine-Gordon beta_SG^2; sinh-Gordon B -> -2 xi/pi
xi = bSG2/8/(1 - bSG2/(8*np.pi))
B = -2*xi/np.pi
m_shg = 8/(np.pi*B)*np.sin(np.pi*B/4)
m0 = 2.0
report(f'm(beta) continued to beta_SG^2 -> -beta_SG^2 equals m_1 = 2M sin(xi/2) with M = m_0/xi (DHN) '
       f'[{abs(m_shg - 2*(m0/xi)*np.sin(xi/2)):.0e}]', abs(m_shg - 2*(m0/xi)*np.sin(xi/2)) < 1e-12)

# ---------------------------------------------------------------- 7. excited reflection factors by the boundary bootstrap
import mpmath as mp
cdb = lambda x, th: mp.sinh(th/2 + 1j*mp.pi*x/4)/mp.sinh(th/2 - 1j*mp.pi*x/4)       # CD (2.3) = eq-blocks, h = 2
Bv, Ev = mp.mpf('0.3'), mp.mpf('1.75')                                               # 2 > E > 1 + B: two excited states
S = lambda th: -1/(cdb(Bv, th)*cdb(2 - Bv, th))
KD = lambda th: cdb(2 - Bv/2, th)*cdb(1 + Bv/2, th)/cdb(1, th)
Kn = lambda n: (lambda th: KD(th)/(cdb(1 - Ev, th)*cdb(1 + Ev, th)) if n == 0 else
                KD(th)*cdb(1 + Ev + Bv, th)*cdb(1 - Ev - Bv, th)
                / (cdb(1 - Ev + (n-1)*Bv, th)*cdb(1 + Ev - (n-1)*Bv, th)*cdb(1 - Ev + n*Bv, th)*cdb(1 + Ev - n*Bv, th)))
dev = 0
for n in (0, 1):
    psi = mp.pi*(Ev - 1 - n*Bv)/2                                                   # pole of K_n at i psi
    for th in (mp.mpc('0.3', '0.1'), mp.mpc('-0.7', '0.2')):
        dev = max(dev, abs(Kn(n)(th)*S(th - 1j*psi)*S(th + 1j*psi)/Kn(n+1)(th) - 1))
    dev = max(dev, 1/abs(Kn(n)(1j*psi + mp.mpf('1e-15')))*1e-5)
report(f'CD (4.4)-(4.6): K_{{n+1}}(theta) = K_n(theta) S(theta - i psi_n) S(theta + i psi_n), psi_n = pi(E - 1 - nB)/2 '
       f'a pole of K_n, for n = 0, 1 [{float(dev):.0e}]', dev < 1e-20)
