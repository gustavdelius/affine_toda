"""Independent test of prp-net-count part 3 and prp-one-loop-mass on the phi^4 kink (referee pass, Chapters 18-19).

Run from this directory:  python3 phi4_check.py     (about 20 s; output phi4_check.log)

A = -d^2 + 1 - (3/2) sech^2(x/2), meson mass m = 1, X(z) = (z-1)/(z+1) * (z-1/2)/(z+1/2).
(a) Fourier-grid Tr(e^{-tA} - e^{-tA0}) versus part 3: erf(sqrt t) + e^{-3t/4} erf(sqrt t / 2)
(b) Jost function from the ODE (Gel'fand-Yaglom) versus the product X at complex lambda
(c) Born constant c = int V dx versus -4 mu sum_w w (step 2 of the proof of prp-one-loop-mass)
(d) Mellin integral (eq-one-loop-mellin) with the directly computed counterterm versus DHN, m(1/(4 sqrt 3) - 3/(2 pi))
"""
import warnings
import numpy as np
import mpmath as mp
from scipy.special import erf
from scipy.integrate import solve_ivp, quad

warnings.filterwarnings('ignore')
V = lambda x: -1.5/np.cosh(x/2)**2


def H_part3(t):
    return erf(np.sqrt(t)) + np.exp(-0.75*t)*erf(0.5*np.sqrt(t))


def grid_spectra(N, L):
    x = L*(np.arange(N) - N//2)/N
    kk = 2*np.pi*np.fft.fftfreq(N, d=L/N)
    D2 = np.real(np.fft.ifft(-(kk**2)[:, None]*np.fft.fft(np.eye(N), axis=0), axis=0))
    A, A0 = -D2 + np.diag(1 + V(x)), -D2 + np.eye(N)
    return np.linalg.eigvalsh((A + A.T)/2), np.linalg.eigvalsh((A0 + A0.T)/2)


def X_ode(lam, Lx=40.0):
    """psi = e^{kappa x} phi, phi(-Lx) = 1; X = phi(+Lx)"""
    kap = np.sqrt(complex(1 - lam))
    kap = kap if kap.real >= 0 else -kap
    sol = solve_ivp(lambda x, y: [y[1], V(x)*y[0] - 2*kap*y[1]], [-Lx, Lx], [1+0j, 0j], rtol=1e-12, atol=1e-14)
    return sol.y[0, -1]


def X_prod(lam):
    z = np.sqrt(complex(1 - lam))
    z = z if z.real >= 0 else -z
    return (z-1)/(z+1)*(z-0.5)/(z+0.5)


if __name__ == "__main__":
    for N, L in [(512, 60), (1024, 80)]:
        ev, ev0 = grid_spectra(N, L)
        d = max(abs(np.sum(np.exp(-t*ev)) - np.sum(np.exp(-t*ev0)) - H_part3(t)) for t in (0.05, 0.3, 0.6, 1, 1.5, 5))
        print(f"(a) grid N={N}, L={L}: lowest eigenvalues {np.round(ev[:3], 9)}; "
              f"max |grid - part 3| at t = 0.05, 0.3, 0.6, 1, 1.5, 5: {d:.1e}")
    lams = [-2+0j, 0.3+0.5j, 0.75+0.01j, 3-1j, 10+2j]
    print(f"(b) Jost function vs product at 5 complex lambda: max diff {max(abs(X_ode(l) - X_prod(l)) for l in lams):.1e}")
    c = quad(V, -200, 200, limit=400)[0]
    print(f"(c) c = int V dx = {c:.12f}; from the roots -4(1 + 1/2) = -6")
    mp.mp.dps = 25
    F = lambda t: t**-1.5*(mp.erf(mp.sqrt(t)) + mp.exp(-0.75*t)*mp.erf(mp.sqrt(t)/2) + t*c*mp.exp(-t)/mp.sqrt(4*mp.pi*t))
    dM = -mp.quad(F, [0, 1, 10, 100, mp.inf])/(4*mp.sqrt(mp.pi))
    dhn = 1/(4*mp.sqrt(3)) - 3/(2*mp.pi)
    print(f"(d) Mellin Delta M = {mp.nstr(dM, 15)}; DHN = {mp.nstr(dhn, 15)}; diff {float(abs(dM - dhn)):.1e}")
