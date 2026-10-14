"""B2: the extra d_5 pole and the absence of classical boundary bound states (prp-dn-uniform-classical, ch. 35).
Units m = beta = 1.  Book S-matrix: smatrix_simply_laced.setup (Dorey's block content), blocks
(x) = sinh(theta/2 + i pi x/2h)/sinh(theta/2 - i pi x/2h), {x} = (x-1)(x+1)/((x-1+B)(x+1-B)).
1. Pole orders of prod_b S_ab(theta, B=1) at theta = i pi z/h, 1 <= z <= h/2, for every particle of d_4, d_5,
   counted from the block content and confirmed by epsilon scaling; momenta p = m_a sinh(theta) and frequencies.
2. d_5 vacuum phi = f (e2 - e4): bulk Hessian M^2(x) = sum_j n_j e^{alpha_j.phi} alpha_j alpha_j^T and boundary Hessian
   R = (1/2) sum_j sqrt(n_j) e^{alpha_j.phi(0)/2} alpha_j alpha_j^T: the e1, e5 (spinor) directions decouple.
3. Matrix Riccati scan: Y = F' F^-1 for the decaying solutions, Y' = M^2 - z - Y^2, Y(-L) = U sqrt(m^2 - z) U^T;
   a bound state at frequency^2 = z needs a zero eigenvalue of Y(0) + R.  d_4 and d_5, grid of z in (0, m_min^2).
"""
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp

TC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'textbook_code')
import os
cwd = os.getcwd(); os.chdir(TC)
ns = {}
exec(compile(open('smatrix_simply_laced.py').read().split('\nB = 0.37')[0], 'smatrix_simply_laced.py', 'exec'), ns)
os.chdir(cwd)
mp.mp.dps = 30


def blk(x, h, th):
    return mp.sinh(th/2 + 1j*mp.pi*x/(2*h))/mp.sinh(th/2 - 1j*mp.pi*x/(2*h))


def prodS(a, n, th, B=1):
    h, bl, mass, _ = ns['setup']('d', n)
    r = mp.mpc(1)
    for b in range(n):
        for x, k in bl[(a, b)].items():
            r *= (blk(x-1, h, th)*blk(x+1, h, th)/(blk(x-1+B, h, th)*blk(x+1-B, h, th)))**k
    return r


for n in (4, 5):
    h, bl, mass, _ = ns['setup']('d', n)
    # absolute masses: m_1 = 2 sqrt2 sin(pi/h) m for the vector; mass[] is normalized to the smallest
    m_vec = 2*np.sqrt(2)*np.sin(np.pi/h)
    scale = m_vec/mass[0]          # index 0 is the vector in this labelling (checked: masses of d_5)
    m2 = 2*np.sqrt(2)*np.sin(2*np.pi/h)
    print(f'd_{n}: h = {h}, masses (units m) {np.round(mass*scale, 6)}, m_2 = {m2:.6f}')
    for a in range(n):
        orders = {}
        for z in range(1, h//2 + 1):
            t0 = 1j*mp.pi*z/h
            vals = [abs(prodS(a, n, t0 + e)) for e in (mp.mpf('1e-6'), mp.mpf('1e-8'))]
            order = int(round(float(mp.log(vals[1]/vals[0])/mp.log(100))))     # |K| ~ eps^-order
            if order:
                orders[z] = order
        poles = {z: o for z, o in orders.items()}
        info = []
        for z, o in poles.items():
            th = np.pi*z/h
            ma = mass[a]*scale
            p = ma*np.sin(th); Om = ma*np.cos(th)
            info.append(f'z={z}: order {o:+d}, p = {p:.6f} i, Omega = {Om:.6f}')
        print(f'   particle {a+1} (m = {mass[a]*scale:.6f}): ' + '; '.join(info))
    # residue of particle 3 at 3 i pi/8 for d_5
    if n == 5:
        t0 = 3j*mp.pi/8
        res = [e*prodS(2, 5, t0 + e) for e in (mp.mpf('1e-8'), mp.mpf('1e-10'), mp.mpf('1e-12'))]
        print('   d_5 particle 3 at 3 i pi/8: eps*K =', [mp.nstr(r, 12) for r in res])

# ---------------------------------------------------------------- 2. Hessians, d_5
def dn_roots(n):
    e = np.eye(n)
    return np.array([-e[0]-e[1]] + [e[i]-e[i+1] for i in range(n-1)] + [e[n-2]+e[n-1]]), np.array([1, 1] + [2]*(n-3) + [1, 1])


A5, n5 = dn_roots(5)
f0 = np.log(4 - 2*np.sqrt(2)); c5 = np.tanh(f0/4)
v5 = np.array([0, 1, 0, -1, 0.])
phi5 = lambda x: 4*np.arctanh(c5*np.exp(2*x))*v5
M2 = lambda A, nj, phi, x: (A.T*(nj*np.exp(A@phi(x))))@A
Rb = lambda A, nj, phi: 0.5*(A.T*(np.sqrt(nj)*np.exp(A@phi(0.0)/2)))@A
sp_idx, mid = [0, 4], [1, 2, 3]
coup = max(np.abs(M2(A5, n5, phi5, x)[np.ix_(sp_idx, mid)]).max() for x in np.linspace(-6, 0, 13))
coupR = np.abs(Rb(A5, n5, phi5)[np.ix_(sp_idx, mid)]).max()
ev0 = np.linalg.eigvalsh(M2(A5, n5, lambda x: 0*v5, 0.0)[np.ix_(sp_idx, sp_idx)])
print(f'd_5 spinor block: masses^2 at -inf {ev0}; max coupling to e2,e3,e4: bulk {coup:.1e}, boundary {coupR:.1e}')

# ---------------------------------------------------------------- 3. Riccati scan
from scipy.integrate import solve_ivp as _ivp
A4, n4 = dn_roots(4)
v4 = np.array([0, 1, -1, 0.])
# u = e^f, u' = (u-1) sqrt(2u(u+2)) (prp-dn-uniform-vacuum); with w = u - 1: w' = w sqrt(2(1+w)(3+w)), in s = -x
sol4 = _ivp(lambda s, w: [-w[0]*np.sqrt(2*(1 + w[0])*(3 + w[0]))], [0, 40], [1/8], method='DOP853',
            rtol=1e-13, atol=1e-300, dense_output=True)
phi4 = lambda x: np.log1p(sol4.sol(-x)[0])*v4


def riccati_min_eig(A, nj, phi, z, L, rtol=1e-11):
    n = A.shape[1]
    M0 = M2(A, nj, lambda x: 0*phi(0.0), 0.0)
    ev, U = np.linalg.eigh(M0)
    Y0 = (U*np.sqrt(ev - z))@U.T
    rhs = lambda x, y: (M2(A, nj, phi, x) - z*np.eye(n) - y.reshape(n, n)@y.reshape(n, n)).ravel()
    s = solve_ivp(rhs, [-L, 0.0], Y0.ravel(), method='DOP853', rtol=rtol, atol=1e-13)
    Y = s.y[:, -1].reshape(n, n)
    J = Y + Rb(A, nj, phi)
    return s.success, np.abs(J - J.T).max(), np.linalg.eigvalsh((J + J.T)/2)[0]


for name, A, nj, phi in (('d_4', A4, n4, phi4), ('d_5', A5, n5, phi5)):
    mmin2 = np.linalg.eigvalsh(M2(A, nj, lambda x: 0*phi(0.0), 0.0)).min()
    zs = list(np.linspace(0, mmin2*(1 - 1e-3), 60)) + [mmin2*(1 - 10.0**-k) for k in range(4, 11)] + [1.0]
    out = [riccati_min_eig(A, nj, phi, z, 20.0) for z in zs]
    ok = all(o[0] for o in out)
    asym = max(o[1] for o in out)
    mins = [o[2] for o in out]
    i = int(np.argmin(mins))
    alt = riccati_min_eig(A, nj, phi, zs[i], 30.0, rtol=1e-12)
    print(f'{name}: m_min^2 = {mmin2:.9f}; {len(zs)} frequencies; success {ok}; asymmetry {asym:.1e}; '
          f'min eigenvalue of Y(0)+R = {mins[i]:.9f} at z = {zs[i]:.10f} (L=30, rtol 1e-12: {alt[2]:.9f}); '
          f'at z = 1: {mins[-1]:.9f}; monotone decreasing in z: {all(np.diff(mins[:60]) < 0)}')
