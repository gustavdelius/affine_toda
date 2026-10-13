"""Which Lagrangian boundary carries the uniform sinh-Gordon amplitude (Part VII, ch. 35, sec-uniform-boundary, "Rank one";
half-line open problem 12(a)).

Conventions (Bajnok-Rim-Zamolodchikov arXiv:0710.4789 = BRZ, Fateev-Zamolodchikov-Zamolodchikov hep-th/0001012 = FZZ):
  action (1/4pi)(d phi)^2 + 2 mu cosh(2 b phi) + M (e^{b phi} + e^{-b phi})|_bdy,  b^2 = beta^2/4pi = B/(2-B),
  M = M_0 cosh(pi b s),  M_0 = sqrt(mu/sin(pi b^2)),  p = b^2/(1+b^2) = B/2.
  The reflection factor is the book's (eq-shg-reflection) at F = 0: its E enters only through X = cos(pi E/2).
Checks:
 1. Ahn-Kim-Rim (akr2002) eq. (10), mu_B^2 = (mu/2) cot(pi b^2/2), is M/M_0 = cos(pi b^2/2), i.e. s = i b/2; the
    continuation cos(pi E/(2-B)) = M/M_0 of cnj-bsg-uv-ir then gives E = B/2 there and E = 0 at M = M_0.
 2. Fateev's conjectured boundary energy of the (+) boundary (hep-th/0103014, h = 2) equals the sinh-Gordon boundary
    energy at E = B/2, not at E = 0; the two differ by -m pi p/16 + O(p^2).
 3. The boundary TBA on an interval with both edges equal reproduces BRZ table 1 (b^2 = 0.8086, s = 0.80 i).
 4. With E = B/2 the TBA matches the boundary Liouville reflection phase of the Lagrangian boundary at s = i b/2
    (M/M_0 = cos(pi b^2/2)) and misses s = 0 (M = M_0); with E = 0 it is the other way round. b^2 = 0.8086, 0.4286.
Runtime: about 1 minute.
"""
import numpy as np
import mpmath as mp
from boundary_common import report

mp.mp.dps = 25

# ---------------------------------------------------------------- 1. AKR (10) and the UV-IR relation
dev1 = 0
for B in (mp.mpf('0.3'), mp.mpf('0.6'), mp.mpf('0.9')):      # 0 < b^2 < 1
    b2 = B/(2 - B)
    MM0 = mp.sqrt((mp.mpf(1)/2*mp.cot(mp.pi*b2/2))*mp.sin(mp.pi*b2))     # mu = 1
    E_of = lambda r: (2 - B)*mp.acos(r)/mp.pi                            # cos(pi E/(2-B)) = M/M_0
    dev1 = max(dev1, abs(MM0 - mp.cos(mp.pi*b2/2)), abs(E_of(MM0) - B/2), abs(E_of(1)))
report(f'AKR (10) mu_B^2 = (mu/2)cot(pi b^2/2) is M/M_0 = cos(pi b^2/2) (s = ib/2); cos(pi E/(2-B)) = M/M_0 puts '
       f'E = B/2 there and E = 0 at M = M_0 [{float(dev1):.0e}]', dev1 < 1e-20)


# ---------------------------------------------------------------- 2. boundary energies
def f_shg(X, p):
    """BRZ boundary energy per edge, symmetric boundary, X = cos(pi E/2), m = 1."""
    return (2*X + 1 - mp.sin(mp.pi*p/2) - mp.cos(mp.pi*p/2))/(4*mp.sin(mp.pi*p))


dev2, sep2 = 0, mp.inf
for p in (mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf('0.4471'), mp.mpf('0.7')):
    fateev = mp.sin(mp.pi/4)/(8*mp.sin(mp.pi*p/4)*mp.cos(mp.pi*(1 - p)/4))
    dev2 = max(dev2, abs(fateev - f_shg(mp.cos(mp.pi*p/2), p)))
    sep2 = min(sep2, abs(fateev - f_shg(1, p)))
p = mp.mpf('1e-3')
lead = (f_shg(mp.cos(mp.pi*p/2), p) - f_shg(1, p))/(-mp.pi*p/16) - 1
report(f'Fateev (+) boundary energy (h = 2) = sinh-Gordon boundary energy at E = B/2 [{float(dev2):.0e}], not at E = 0 '
       f'[min gap {float(sep2):.1e}]; difference -m pi p/16 (1 + O(p)) [{float(lead):.0e} at p = 1e-3]',
       dev2 < 1e-15 and sep2 > 1e-3 and abs(lead) < 1e-2)


# ---------------------------------------------------------------- 3./4. boundary TBA versus the Liouville phase
def kernel(th, p):
    return 4*np.sin(np.pi*p)*np.cosh(th)/(np.cosh(2*th) - np.cos(2*np.pi*p))/(2*np.pi)


def fugacity(th, p, X):
    c = np.cosh(th)
    cp, sp_ = np.cos(np.pi*p/2), np.sin(np.pi*p/2)
    return np.tanh(th/2)**2*(c + cp)/(c - cp)*(c + sp_)/(c - sp_)*((c - X)/(c + X))**2


def tba_energy(mR, p, X, N=4097, Lmax=25.0, tol=1e-13):
    th = np.linspace(-Lmax, Lmax, N)
    d = th[1] - th[0]
    lam, drive = fugacity(th, p, X), 2*mR*np.cosh(th)
    nfft = 2**int(np.ceil(np.log2(3*N)))
    Kf = np.fft.rfft(kernel(np.linspace(-2*Lmax, 2*Lmax, 2*N - 1), p)*d, nfft)
    eps = drive.copy()
    for _ in range(5000):
        conv = np.fft.irfft(Kf*np.fft.rfft(np.log1p(lam*np.exp(-eps)), nfft), nfft)[N - 1:2*N - 1]
        new = drive - conv
        err = np.max(np.abs(new - eps))
        eps = (eps + new)/2
        if err < tol:
            break
    return -np.sum(np.cosh(th)*np.log1p(lam*np.exp(-eps)))*d/(4*np.pi)


def delta_tba(mR, b2, X):
    """Delta^TBA(P) of BRZ: pi + 2 P Q log(R/pi) - Delta_L(P)/2, with P from c_eff."""
    p = float(b2/(1 + b2))
    E12 = tba_energy(mR, p, X) + mR/(8*np.sin(np.pi*p)) + 2*float(f_shg(mp.mpf(X), mp.mpf(p)))
    P = mp.sqrt((1 + 24*mR*E12/np.pi)/24)
    bb = mp.sqrt(b2)
    pm = mp.mpf(p)
    cm = (pm**pm*(1 - pm)**(1 - pm)*mp.gamma(pm/2)*mp.gamma((1 - pm)/2)/(8*mp.sqrt(mp.pi)))**(2 + 2*b2)
    dL = -(2*P/bb)*mp.log(cm) + 2*mp.im(mp.loggamma(1 + 2j*bb*P)) + 2*mp.im(mp.loggamma(1 + 2j*P/bb))
    return mp.pi + 2*P*(bb + 1/bb)*mp.log(mR/mp.pi) - dL/2, P


def delta_liouville(P, s, b2):
    """Boundary Liouville reflection phase Delta(P|s,s) (FZZ; BRZ eq. (Dint))."""
    b = mp.sqrt(b2)
    Q = 1/b - b
    f = lambda t: 2*mp.sin(2*P*t)/t*(mp.cos(s*t)**2 - mp.cosh(b*t/2)*mp.cosh(t/(2*b))*mp.cosh(Q*t/2)) \
        / (mp.sinh(b*t)*mp.sinh(t/b))
    return mp.re(mp.quad(f, [0, mp.inf]))


b2 = mp.mpf('0.8086')
sB = mp.mpc(0, '0.80')
dt, P = delta_tba(1e-6, b2, float(mp.re(mp.cosh(b2/(1 + b2)*mp.pi*sB/mp.sqrt(b2)))))
dev3 = abs(dt - delta_liouville(P, sB, b2))
report(f'boundary TBA reproduces BRZ table 1 (b^2 = 0.8086, s = 0.80i, mR = 1e-6) [{float(dev3):.0e}]', dev3 < 1e-10)

rows, ok = [], True
for b2 in (mp.mpf('0.8086'), mp.mpf('0.4286')):
    p, sa = b2/(1 + b2), 1j*mp.sqrt(b2)/2
    for E, s_match, s_other in ((p, sa, 0), (0, 0, sa)):
        dt, P = delta_tba(1e-6, b2, float(mp.cos(mp.pi*E/2)))
        hit, miss = abs(dt - delta_liouville(P, s_match, b2)), abs(dt - delta_liouville(P, s_other, b2))
        rows.append((float(hit), float(miss)))
        ok &= hit < 1e-10 and miss > 1e-2
report('TBA with E = B/2 matches the Lagrangian boundary at s = ib/2 (M/M_0 = cos(pi b^2/2)) and misses s = 0; '
       'E = 0 matches s = 0 and misses s = ib/2 (b^2 = 0.8086, 0.4286; mR = 1e-6) [hit/miss: '
       + ', '.join(f'{h:.0e}/{m:.2f}' for h, m in rows) + ']', ok)
