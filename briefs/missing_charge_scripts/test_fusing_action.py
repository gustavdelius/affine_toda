"""Compute exact integrals of (u_x)^2/2, (u_t)^2/2, W(u), E_Mink, and L_Eucl for the fusing-angle solution tau_j^+(x,t)."""
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../textbook_code/thimbles')))

import numpy as np
from lie import An

g = An(3)  # a_3^(1)

x = np.linspace(-25.0, 25.0, 80001)
dx = x[1] - x[0]
j = np.arange(4)[:, None]

print("=== 1. General (C1, C4) chambers ===")
for label, C1, C4 in [
    ("static C1=0 (1,3)", 0.0, 1.0 + 0.5j),
    ("small C1 (2,3)", 0.5 + 0.3j, 1.0 + 0.5j),
    ("large C1 missing", 3.0 + 2.0j, 1.0 + 0.2j),
    ("very large C1 missing", 8.0 - 5.0j, 1.0 - 0.3j),
]:
    A = C1 * np.exp(x)[None, :]
    B = C4 * np.exp(2.0 * x)[None, :]
    tau = 1.0 + (1j ** j) * A + ((-1.0) ** j) * B
    tau_x = (1j ** j) * A + 2.0 * ((-1.0) ** j) * B
    tau_t = -1j * (1j ** j) * A

    ph = np.unwrap(np.angle(tau), axis=1)
    u = 1j * (g.alpha.T @ (np.log(np.abs(tau)) + 1j * ph))
    ux = 1j * (g.alpha.T @ (tau_x / tau))
    ut = 1j * (g.alpha.T @ (tau_t / tau))
    W = np.sum(1.0 - np.exp(1j * (g.alpha @ u)), axis=0)

    q = (u[:, -1] - u[:, 0]) / (2 * np.pi)
    lab = g.label(q.real, tol=1e-3)
    name = ''.join(str(i + 1) for i, c in enumerate(lab) for _ in range(c)) if lab else str(q.real)

    I_kin_x = np.sum(0.5 * np.sum(ux ** 2, axis=0)) * dx
    I_kin_t = np.sum(0.5 * np.sum(ut ** 2, axis=0)) * dx
    I_pot = np.sum(W) * dx
    E_mink = I_kin_t + I_kin_x + I_pot
    L_eucl = -I_kin_t + I_kin_x + I_pot
    print(f"{label:24s} [charge {name}]: E_Mink={E_mink.real:10.6f}{E_mink.imag:+10.6f}j | L_Eucl={L_eucl.real:10.6f}{L_eucl.imag:+10.6f}j")

print("\n=== 2. Testing C4 = -1.0, C1 = i * eps (missing charge 12, strictly real L_Eucl!) ===")
for eps in (1.9, 1.5, 1.0, 0.5, 0.2, 0.1, 0.05):
    C1 = 1j * eps
    C4 = -1.0
    A = C1 * np.exp(x)[None, :]
    B = C4 * np.exp(2.0 * x)[None, :]
    tau = 1.0 + (1j ** j) * A + ((-1.0) ** j) * B
    tau_x = (1j ** j) * A + 2.0 * ((-1.0) ** j) * B
    tau_t = -1j * (1j ** j) * A

    ph = np.unwrap(np.angle(tau), axis=1)
    u = 1j * (g.alpha.T @ (np.log(np.abs(tau)) + 1j * ph))
    ux = 1j * (g.alpha.T @ (tau_x / tau))
    ut = 1j * (g.alpha.T @ (tau_t / tau))
    W = np.sum(1.0 - np.exp(1j * (g.alpha @ u)), axis=0)

    q = (u[:, -1] - u[:, 0]) / (2 * np.pi)
    lab = g.label(q.real, tol=1e-3)
    name = ''.join(str(i + 1) for i, c in enumerate(lab) for _ in range(c)) if lab else str(q.real)

    I_kin_x = np.sum(0.5 * np.sum(ux ** 2, axis=0)) * dx
    I_kin_t = np.sum(0.5 * np.sum(ut ** 2, axis=0)) * dx
    I_pot = np.sum(W) * dx
    E_mink = I_kin_t + I_kin_x + I_pot
    L_eucl = -I_kin_t + I_kin_x + I_pot
    print(f"eps={eps:5.2f} [charge {name}]: int(u_x^2/2)={I_kin_x.real:9.5f} | int(u_t^2/2)={I_kin_t.real:9.5f} | int(W)={I_pot.real:9.5f} | E_Mink={E_mink.real:9.5f} | L_Eucl={L_eucl.real:9.5f} (Im={L_eucl.imag:.1e})")
