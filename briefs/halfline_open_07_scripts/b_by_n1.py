"""Problem 7(b), one soliton of a_n on an interval, uniform (left) / s=+1 (right) and s=+1 / s=+1.
T = A_min K(theta) (UN) or (A_min K(theta))^2 (NN); K = function of the cyclic shift Z, so the Z-eigenvectors
u_zeta diagonalize it.  Checks |Lambda(zeta)| = rho(zeta)^p, rho = |e^{om th} zeta - w| / |e^{om th} zeta - wbar|,
w = e^{-i pi om/2}, then solves the Bethe-Yang equation 2 m L sinh(theta) - i log Lambda = 2 pi k for complex theta
(scalar phases omitted; they are real-analytic and unimodular and do not change Im theta at leading order)."""
import numpy as np
from p07_common import k_book, A_min
n = 2; N = 3
def lam(theta, zeta, om, p):
    u = zeta**(-np.arange(N))/np.sqrt(N)
    K = A_min(n, theta, 1, om)*k_book(n, theta, 1, om)
    return (u.conj() @ K @ u)**p
for om in (0.37, 1.43):
    w = np.exp(-1j*np.pi*om/2)
    for th in (0.0, 0.5, 2.0):
        out = []
        for zeta in np.exp(2j*np.pi*np.arange(3)/3):
            r = abs(np.exp(om*th)*zeta - w)/abs(np.exp(om*th)*zeta - np.conj(w))
            out.append(f'zeta={np.angle(zeta)/(2*np.pi/3):+.0f}: |L|={abs(lam(th, zeta, om, 1)):.4f} rho={r:.4f}')
        print(f'om={om} th={th}: ' + '; '.join(out))
# Newton for 2 mL sinh th - i log Lambda(th) = 2 pi k
def solve(zeta, om, p, mL, k, th0):
    th = complex(th0)
    for _ in range(60):
        f = lambda t: 2*mL*np.sinh(t) - 1j*np.log(lam(t, zeta, om, p)) - 2*np.pi*k
        h = 1e-6
        d = (f(th + h) - f(th - h))/(2*h)
        dth = f(th)/d; th -= dth
        if abs(dth) < 1e-13:
            break
    return th, abs(f(th))
print('\nBethe-Yang, one soliton, m L = 10 (energy E = m cosh theta in units of m):')
for cfg, p in (('UN', 1), ('NN', 2)):
    for om in (0.37, 1.43):
        for zeta in np.exp(2j*np.pi*np.arange(3)/3):
            rows = []
            for k in (1, 2, 4, 8):
                th, res = solve(zeta, om, p, 10.0, k, np.arcsinh(np.pi*k/10))
                rows.append(f'k={k}: E={np.cosh(th).real:.5f}{np.cosh(th).imag:+.5f}i')
            print(f'{cfg} om={om} zeta=e^(2 pi i {np.angle(zeta)/(2*np.pi/3):+.0f}/3): ' + '  '.join(rows))
