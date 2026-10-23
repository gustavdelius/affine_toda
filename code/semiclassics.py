"""Semiclassics checks for Part V (chapters 16-17).

1. Transmission-factor mass formula (eq-dhn-psi, = prp-one-loop-mass for one channel):
   Delta M = (m/2) sum_w ord_w(X) psi(arcsin w), psi(t) = (t cos t - sin t)/pi,
   reproduces the sine-Gordon kink (-m/pi) and the phi^4 kink m(1/(4 sqrt 3) - 3/(2 pi)),
   and equals the Cahill-Comtet-Glauber form -(m/pi) sum_bound (sin t - t cos t).
2. a_n^(1) solitons: summing over channels with X_ab(z) = (z - cos(A-B))/(z - cos(A+B)) gives
   Hollowood's Delta M_a = -m_a [h/(2 pi) - cot(pi/h)/4]  (eq-hollowood-mass).
3. Zero-dimensional i phi^3 integral (sec-picard-lefschetz): the saddle at x = i has lower action
   (-1/6) but does not contribute: |I(lambda)| <= sqrt(2 pi/lambda), and I matches the Gaussian
   expansion around x = 0.
4. Gaussian thimble integral counts algebraic multiplicity (prp-thimble-gaussian, exr-jordan-determinant).
5. H_g = -d^2 + g e^{ix} on the circle: spectrum k^2 with 2x2 Jordan blocks at n^2, n >= 1.
"""
import numpy as np
from mpmath import mp, quad, exp, sqrt, pi, inf, mpf

mp.dps = 30


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


psi = lambda t: (t*np.cos(t) - np.sin(t))/np.pi


def dM(m, zeros, poles):
    return 0.5*m*(sum(psi(np.arcsin(w)) for w in zeros) - sum(psi(np.arcsin(w)) for w in poles))


m = 1.0
sg = dM(m, [1.0], [-1.0])
report(f'sine-Gordon kink: Delta M = {sg:.12f} = -m/pi', abs(sg + 1/np.pi) < 1e-14)
phi4 = dM(m, [1.0, 0.5], [-1.0, -0.5])        # bound states omega = 0, sqrt(3)/2 m
report(f'phi^4 kink: Delta M = {phi4:.12f} = m(1/(4 sqrt3) - 3/(2pi))', abs(phi4 - (1/(4*np.sqrt(3)) - 3/(2*np.pi))) < 1e-14)
ccg = -(m/np.pi)*sum(np.sin(t) - t*np.cos(t) for t in (np.pi/2, np.pi/6))
report('psi-form = Cahill-Comtet-Glauber form for phi^4', abs(ccg - phi4) < 1e-14)

ok = True
for n in range(1, 8):
    h = n + 1
    mb = lambda b: 2*np.sin(np.pi*b/h)
    for a in range(1, n+1):
        A = np.pi*a/h
        tot = sum(dM(mb(b), [np.cos(A - np.pi*b/h)], [np.cos(A + np.pi*b/h)]) for b in range(1, n+1))
        holl = -mb(a)*(h/(2*np.pi) - 0.25/np.tan(np.pi/h))
        ok &= abs(tot - holl) < 1e-12
report("a_n^(1), n<=7: channel sum reproduces Hollowood's Delta M_a = -m_a[h/2pi - cot(pi/h)/4]", ok)

# 3. i phi^3 in zero dimensions
def I(lam):
    f = lambda x: exp(-lam*(x**2/2))*mp.cos(lam*x**3/3)   # imaginary part odd, cancels
    return 2*quad(f, [0, 1, 2, 4, inf])


ok3 = True
for lam in (mpf(5), mpf(10), mpf(20)):
    val = I(lam); gauss = sqrt(2*pi/lam)
    # asymptotic series around x=0: 1 - 5/(6 lam) + O(lam^-2), from -(lam^2/18) <x^6>, <x^6> = 15/lam^3
    ok3 &= abs(val) <= gauss and abs(val/gauss - (1 - mpf(5)/(6*lam))) < 6/lam**2   # next coefficient is +10395/1944 = 5.35
report('i phi^3 integral: bounded by the Gaussian, matches 1 - 5/(6 lambda) around x=0; x=i (S=-1/6) absent', ok3)

# 4. Gaussian over a Jordan block: A = [[2+i,1],[1,2-i]] has eigenvalue 2 (alg. mult. 2), det 4
A = np.array([[2+1j, 1], [1, 2-1j]])
ev = np.linalg.eigvals(A)
from scipy import integrate
re = lambda x, y: np.real(np.exp(-0.5*np.array([x, y]) @ A @ np.array([x, y])))
val = integrate.dblquad(re, -12, 12, -12, 12, epsabs=1e-11)[0]
report(f'Jordan-block Gaussian: eigenvalues {np.round(ev, 6)}, integral {val:.10f} = 2pi/sqrt(det A) = {2*np.pi/2:.10f}',
       abs(val - np.pi) < 1e-8 and np.allclose(ev, [2, 2], atol=1e-6) and abs(np.linalg.det(A) - 4) < 1e-12)

# 5. H_g on the circle in the Fourier basis e^{ikx}, |k| <= K: (H_g)_{k,k} = k^2, (H_g)_{k+1,k} = g
K, g = 8, 0.7
ks = np.arange(-K, K+1); N = len(ks)
H = np.diag(ks.astype(float)**2).astype(complex)
for i in range(N-1):
    H[i+1, i] = g
ok5 = True
for nn in range(1, 5):
    ranks = np.linalg.matrix_rank(H - nn**2*np.eye(N), tol=1e-9)
    ok5 &= (N - ranks == 1)                         # geometric multiplicity 1, algebraic 2
    ok5 &= np.linalg.matrix_rank(np.linalg.matrix_power(H - nn**2*np.eye(N), 2), tol=1e-9) == N - 2
report('H_g = -d^2 + g e^{ix}: eigenvalue n^2 (1<=n<=4) has algebraic multiplicity 2, geometric 1', ok5)
