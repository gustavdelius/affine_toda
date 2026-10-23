"""Continuum Hirota solitons of imaginary-coupling a_n^(1) Toda, in u = beta*phi variables.

Static energy  E[u] = int (1/2) u'^2 + kappa W(u),  W(u) = -sum_j (e^{i alpha_j.u} - 1).
Field equation u'' = kappa grad W(u),  grad W = -i sum_j alpha_j e^{i alpha_j.u}.
Hirota form u = s * i * sum_j alpha_j log tau_j  (sign s fixed numerically below).
General (time-dependent) multi-soliton:
  tau_j = sum_{subsets} prod_{k in sub} omega^{j a_k} E_k * prod_{k<l in sub} A_{a_k a_l}(theta_k - theta_l)
  E_k = exp(m_{a_k} (x cosh th_k - t sinh th_k) + xi_k), m_a = 2 sqrt(kappa) sin(pi a/h).
"""
import itertools
import numpy as np
from lie import An


def A_int(n, a, b, th):
    """Two-soliton interaction coefficient for a_n^(1) (Hollowood/Fring-Liao-Olive form)."""
    h = n + 1
    A, B = np.pi * a / h, np.pi * b / h
    c = np.cosh(th)
    num = (c - np.cos(A - B)) * 1.0
    den = (c - np.cos(A + B))
    # full form: [cosh th - cos(A-B)] [..] / [cosh th - cos(A+B)] ; Hollowood's A_ab:
    # A_ab = (cosh th - cos((a-b)pi/h)) (cosh th - cos((a+b)pi/h))^{-1}  times similar? keep simplest
    return -num / den if False else num / den


def tau(g, species, thetas, xis, x, t=0.0, kappa=1.0, Afun=None):
    n = g.n
    h = g.h
    om = np.exp(2j * np.pi / h)
    m = 2 * np.sqrt(kappa) * np.sin(np.pi * np.array(species) / h)
    Ek = [np.exp(m[k] * (x * np.cosh(thetas[k]) - t * np.sinh(thetas[k])) + xis[k]) for k in range(len(species))]
    if Afun is None:
        Afun = lambda a, b, th: A_int(n, a, b, th)
    taus = []
    for j in range(h):
        tj = np.ones_like(x, dtype=complex)
        for size in range(1, len(species) + 1):
            for sub in itertools.combinations(range(len(species)), size):
                term = np.ones_like(x, dtype=complex)
                for k in sub:
                    term = term * om ** (j * species[k]) * Ek[k]
                for k, l in itertools.combinations(sub, 2):
                    term = term * Afun(species[k], species[l], thetas[k] - thetas[l])
                tj = tj + term
        taus.append(tj)
    return np.array(taus)  # (h, len(x))


def field(g, taus, sign=1):
    # u = sign * i * sum_j alpha_j log tau_j with continuous branch along x
    logs = np.log(np.abs(taus)) + 1j * np.unwrap(np.angle(taus), axis=1)
    return sign * 1j * (g.alpha.T @ logs)  # (r, len(x))


def check_eq(g, u, x, kappa=1.0):
    dx = x[1] - x[0]
    upp = (u[:, 2:] - 2 * u[:, 1:-1] + u[:, :-2]) / dx ** 2
    uu = u[:, 1:-1]
    gradW = -1j * (g.alpha.T @ np.exp(1j * (g.alpha @ uu)))
    return np.max(np.abs(upp - kappa * gradW))


def charge(g, species, thetas, xis, kappa=1.0, X=40.0, N=40001, sign=1):
    x = np.linspace(-X, X, N)
    taus = tau(g, species, thetas, xis, x, kappa=kappa)
    if np.min(np.abs(taus)) < 1e-9:
        return None
    u = field(g, taus, sign)
    q = (u[:, -1] - u[:, 0]) / (2 * np.pi)
    return q


if __name__ == "__main__":
    for n in (2, 3):
        g = An(n)
        x = np.linspace(-10, 10, 4001)
        for a in range(1, n + 1):
            for sign in (1, -1):
                taus = tau(g, [a], [0.0], [0.3 + 0.4j], x)
                u = field(g, taus, sign)
                print(n, a, sign, "residual", check_eq(g, u, x))
