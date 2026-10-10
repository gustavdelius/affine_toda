import numpy as np
from scipy.special import loggamma
from common import *
import os; omega = float(os.environ.get('OMEGA', '2.37')); q = float(__import__("os").environ.get("QSIGN", "1"))*np.exp(-1j*np.pi*omega); T = 3*omega + 1; mu = 2*T
S = Spinor33(q); C = np.load('C.npy'); Ci = np.linalg.inv(C); tc = np.load('t.npy')[0]
A = [omega, 2*omega + 0.5, 3*omega]                         # zeros of the crossing factor c
def c_of(t): return np.prod([np.sin(np.pi*(t - a)) for a in A], axis=0)
def lf1(z):
    z = np.asarray(z, dtype=complex)
    return sum(loggamma(z - a) + loggamma(1 + z + a) for a in A) - len(A)*np.log(np.pi)
J = 6000
js = np.arange(1, J+1); Aj = 2*T*(js-1)
def logf_rel(t):
    d = lambda z, z0: lf1(z) - lf1(z0)
    g = (d(t + Aj, Aj) - d(t + Aj + T, Aj + T) + d(-t + Aj + T, Aj + T) - d(-t + Aj + 2*T, Aj + 2*T))
    return np.sum(g)
def F(theta, sign=1):
    t = mu*theta/(2j*np.pi)
    lf = logf_rel(t) - np.log(c_of(0.0) + 0j)
    return sign*c_of(t)*np.exp(lf)
def Sbraid(theta): return F(theta)*S.R(np.exp(mu*theta))
# --- unitarity and crossing at generic complex rapidities
for th in (0.31+0.17j, -0.8+0.9j):
    U = Sbraid(th) @ Sbraid(-th)
    print(f"unitarity  S(th)S(-th) = 1 : max dev {np.abs(U - np.eye(64)).max():.1e}   (th={th})")
def pt1(Mm): return Mm.reshape(8, 8, 8, 8).transpose(2, 1, 0, 3).reshape(64, 64)
for th in (0.31+0.17j, 1.1+0.6j):
    S12 = lambda z: P @ Sbraid(z)                   # R-form
    S21 = lambda z: S12(z) @ np.eye(64)             # placeholder replaced below
    lhs = P @ Sbraid(1j*np.pi - th)
    R21 = P @ (P @ Sbraid(th)) @ P                 # S_21 in R-form
    rhs = np.kron(C, np.eye(8)) @ pt1(R21) @ np.kron(Ci, np.eye(8))
    print(f"crossing   S_12(i pi - th) = (C x 1) S_21(th)^t1 (C^-1 x 1) : rel dev {np.abs(lhs - rhs).max()/np.abs(lhs).max():.1e}   (th={th})")
np.save('par.npy', np.array([omega]))

# ---------------- sign-fixing CDD factor and full checks ----------------
def g_cdd(theta): return np.sinh(theta/2 - 1j*np.pi/4)/np.sinh(theta/2 + 1j*np.pi/4)
def Sfull(theta): return g_cdd(theta)*Sbraid(theta)
th = 0.31+0.17j
lhs = P @ Sfull(1j*np.pi - th); R21 = P @ (P @ Sfull(th)) @ P
rhs = np.kron(C, np.eye(8)) @ pt1(R21) @ np.kron(Ci, np.eye(8))
print(f"with CDD factor: crossing rel dev {np.abs(lhs - rhs).max()/np.abs(lhs).max():.1e};  unitarity dev {np.abs(Sfull(th) @ Sfull(-th) - np.eye(64)).max():.1e}")
# ---------------- pole scan of the physical strip theta = i u ----------------
us = np.linspace(2e-4, np.pi - 2e-4, 8001)
norms = np.array([np.linalg.norm(Sfull(1j*u)) for u in us])
med = np.median(norms)
peaks = [us[i] for i in range(1, len(us)-1) if norms[i] > norms[i-1] and norms[i] > norms[i+1] and norms[i] > 20*med]
tt = lambda u: T*u/np.pi
pred = []
for a, lab, rk in ((omega, "soliton 2 (29)", 29), (2*omega+0.5, "soliton 1 (8)", 8), (3*omega, "breather (singlet)", 1)):
    k = 0
    while a - k > 0: pred.append((a - k, f"{lab}, k={k}")); k += 1
    k = 0
    while T - a + k < T: pred.append((T - a + k, f"crossed {lab}, k={k}")); k += 1
pred.sort()
def rank_at(t0):
    u0 = np.pi*t0/T; out = []
    for d in (1e-6, 1e-8):
        M = d*Sfull(1j*(u0 + d)); sv = np.linalg.svd(M, compute_uv=False); out.append(int(np.sum(sv > 1e-3*sv[0])))
    return out
print(f"\nS_33 at omega={omega}: {len(peaks)} norm peaks found on the grid; predicted poles: {len(pred)}")
matched = set()
for t0, lab in pred:
    near = [p for p in peaks if abs(tt(p) - t0) < 3*T/8000*3]
    rk = rank_at(t0)
    print(f"   t = {t0:7.4f} (u = {np.pi*t0/T:.5f})  {lab:32s}  grid peak: {'yes' if near else 'NO '}   residue rank {rk}")
    for p in near: matched.add(round(p, 6))
extra = [p for p in peaks if round(p, 6) not in matched]
print("unexpected peaks (t):", [round(tt(p), 4) for p in extra])
