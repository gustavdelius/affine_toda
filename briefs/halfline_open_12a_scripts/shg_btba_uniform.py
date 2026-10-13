"""Sinh-Gordon boundary TBA versus the boundary Liouville reflection amplitude (BLRA)
at the eps_hat = 0 point (E = B/2 = p) and at E = 0.

Formulas from Bajnok-Rim-Zamolodchikov, arXiv:0710.4789 (BRZ):
  action (1/4pi)(d phi)^2 + 2 mu cosh(2 b phi) + M (e^{b phi}+e^{-b phi})|_bdy,
  M = M0 cosh(pi b s), M0 = sqrt(mu/sin(pi b^2)),  p = b^2/(1+b^2).
  BTBA:  eps = 2mR cosh(th) - phi * log(1+lam e^{-eps}),  E(R) = -(m/4pi) int cosh(th) L.
  E_12 = E(R) + calE R + 2 f,  c_eff = -24 R E_12/pi,  P = sqrt((1-c_eff)/24).
  Delta^TBA = pi + 2 P Q log(R/pi) - Delta_L(P)/2   (both edges symmetric, s1=s2=s)
  compared with Delta(P|s,s) of BRZ eq. (Dint).
The reflection factor is the book's (eq-shg-reflection) with F=0; X = cos(pi E/2).
Validation: BRZ table 1 (b^2=0.8086, s=0.80i).

Run: python3 shg_btba_uniform.py
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 20


def fugacity(th, p, X):
    c = np.cosh(th)
    cp, sp = np.cos(np.pi * p / 2), np.sin(np.pi * p / 2)
    return (np.tanh(th / 2) ** 2 * (c + cp) / (c - cp) * (c + sp) / (c - sp)
            * ((c - X) / (c + X)) ** 2)


def kernel(th, p):
    return 4 * np.sin(np.pi * p) * np.cosh(th) / (np.cosh(2 * th) - np.cos(2 * np.pi * p)) / (2 * np.pi)


def solve_btba(mR, p, X, h=0.004, tol=1e-12, maxit=5000):
    Lam = np.log(80.0 / mR)
    N = int(2 * Lam / h) | 1
    th = np.linspace(-Lam, Lam, N)
    dth = th[1] - th[0]
    lam = fugacity(th, p, X)
    # kernel on the difference grid for linear convolution
    M = 2 * N - 1
    kth = (np.arange(M) - (N - 1)) * dth
    K = kernel(kth, p) * dth
    nfft = 1 << int(np.ceil(np.log2(M + N)))
    Kf = np.fft.rfft(K, nfft)
    drive = 2 * mR * np.cosh(th)
    eps = drive.copy()
    for it in range(maxit):
        L = np.log1p(lam * np.exp(-eps))
        conv = np.fft.irfft(Kf * np.fft.rfft(L, nfft), nfft)[N - 1:2 * N - 1]
        new = drive - conv
        err = np.max(np.abs(new - eps))
        eps = 0.5 * eps + 0.5 * new
        if err < tol:
            break
    L = np.log1p(lam * np.exp(-eps))
    E = -1 / (4 * np.pi) * np.sum(np.cosh(th) * L) * dth
    return E, it, err


def liouville_phase(P, b2):
    b = mp.sqrt(b2)
    p = b2 / (1 + b2)
    pimugamma = (1 / (8 * mp.sqrt(mp.pi)) * p ** p * (1 - p) ** (1 - p)
                 * mp.gamma(p / 2) * mp.gamma((1 - p) / 2)) ** (2 + 2 * b2)
    return (-(2 * P / b) * mp.log(pimugamma) + 2 * mp.im(mp.loggamma(1 + 2j * b * P))
            + 2 * mp.im(mp.loggamma(1 + 2j * P / b)))


def Delta_s(P, s, b2):
    b = mp.sqrt(b2)
    q = 1 / b - b
    def f(t):
        num = mp.cos(s * t) ** 2 - mp.cosh(b * t / 2) * mp.cosh(t / (2 * b)) * mp.cosh(q * t / 2)
        return mp.re(mp.sin(2 * P * t) / t * num / (mp.sinh(b * t) * mp.sinh(t / b)))
    return 2 * mp.quad(f, [mp.mpf('1e-12'), 1, 5, 20, 60, 200])


def compare(b2, X, s_list, mR_list, label):
    b = np.sqrt(b2)
    p = b2 / (1 + b2)
    Q = b + 1 / b
    calE = 1 / (8 * np.sin(np.pi * p))
    f = 1 / (4 * np.sin(np.pi * p)) * (2 * X + 2 - np.sin(np.pi * p / 2) - np.cos(np.pi * p / 2) - 1)
    print(f"\n{label}: b^2={b2}, p={p:.6f}, X=cos(pi E/2)={X:.10f}, E={2/np.pi*np.arccos(X):.6f}")
    print("   mR        P_TBA          Delta^TBA          " + "  ".join(f"Delta(P|{s},{s})" for s in s_list))
    for mR in mR_list:
        E, it, err = solve_btba(mR, p, X)
        R = mR
        E12 = E + calE * R + 2 * f
        ceff = -24 * R * E12 / np.pi
        P = np.sqrt((1 - ceff) / 24)
        DT = np.pi + 2 * P * Q * np.log(R / np.pi) - float(liouville_phase(mp.mpf(P), mp.mpf(b2))) / 2
        vals = [float(Delta_s(mp.mpf(P), s, mp.mpf(b2))) for s in s_list]
        print(f"  {mR:.1e}  {P:.12f}  {DT:.12f}  " + "  ".join(f"{v:.12f}" for v in vals)
              + f"   (it={it}, err={err:.1e})")


if __name__ == "__main__":
    # validation: BRZ table 1, b^2=0.8086, s=0.80i
    b2 = 0.8086
    b = np.sqrt(b2); p = b2 / (1 + b2)
    s = 0.80j
    X = np.cosh(p * np.pi * s / b).real          # cos(pi E/2) = cosh(p eta), b eta = pi s
    compare(b2, X, [mp.mpc(0, 0.80)], [1e-8, 1e-6, 1e-4, 1e-2], "Validation (BRZ table 1)")
    # the test: eps_hat = 0 point E = p, against s = ib/2 (AKR) and s = 0 (M = M0)
    for b2 in [0.8086, 0.4286]:
        b = np.sqrt(b2); p = b2 / (1 + b2)
        sA, s0 = mp.mpc(0, b / 2), mp.mpc(0, 0)
        compare(b2, np.cos(np.pi * p / 2), [sA, s0], [1e-8, 1e-6, 1e-4, 1e-2], "eps_hat=0, E=p")
        compare(b2, 1.0, [sA, s0], [1e-8, 1e-6, 1e-4, 1e-2], "E=0")
