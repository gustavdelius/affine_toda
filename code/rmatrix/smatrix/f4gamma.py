import numpy as np, time, sys
from fractions import Fraction as F
from scipy.special import loggamma
exec(open('f4halftower.py').read().split("for name, SB, sgn, Tf, c in")[0])
def blocks_val(ys, sgn, T, th):
    return sgn*np.prod([np.sinh(th/2 + 1j*np.pi*y/(2*T))/np.sinh(th/2 - 1j*np.pi*y/(2*T)) for y in ys])
def Fmaker(A, T, J=1500):
    js = np.arange(1, J+1); Aj = 2*T*(js-1)
    def lf1(z):
        z = np.asarray(z, dtype=complex); return sum(loggamma(z - a) + loggamma(1 + z + a) for a in A)
    def logf(t):
        d = lambda z, z0: lf1(z) - lf1(z0); r = 0.371
        return np.sum(d(t + Aj, Aj + r) - d(t + Aj + T, Aj + T + r) + d(-t + Aj + T, Aj + T + r) - d(-t + Aj + 2*T, Aj + 2*T + r))
    c_of = lambda t: np.prod([np.sin(np.pi*(t - a)) for a in A])
    K = c_of(0.123)*c_of(-0.123)*np.exp(logf(0.123) + logf(-0.123)); nrm = np.sqrt(K + 0j)
    Fr = lambda t: c_of(t)*np.exp(logf(t))/nrm; sg = np.sign(Fr(1e-6).real)
    return lambda theta: sg*Fr(T*theta/(1j*np.pi))
ths = [0.31+0.17j, -0.45+0.62j, 0.8-0.25j, 0.12+1.1j]
def scan(label, Tf, Aexact, ours_ex, cds_ex):
    T = val(Tf, w0); A = [val(a, w0) for a in Aexact]
    ours_y = [val(y, w0) for y in ours_ex[0]]; cds_y = [val(y, w0) for y in cds_ex[0]]
    need = np.array([blocks_val(cds_y, cds_ex[1], T, th)/blocks_val(ours_y, ours_ex[1], T, th) for th in ths])
    bet = np.pi*(T - 1)/T
    FA = Fmaker(A, T)
    base = {th: [FA(th + 1j*bet), FA(th), FA(th - 1j*bet)] for th in ths}
    res = []
    t0 = time.time()
    for al in range(0, 7):
        for b2 in range(-6, 13):
            e = al*w0 + b2/2
            if not (0 < e < T/2 + 1e-9): continue
            if any(abs(((e - a) % 1)) < 1e-9 or abs(((e - a) % 1) - 1) < 1e-9 for a in A + [T - a for a in A]): pass
            FE = Fmaker(A + [e, T - e], T)
            W = []
            for th in ths:
                z = [FE(th + 1j*bet)/base[th][0], FE(th)/base[th][1], FE(th - 1j*bet)/base[th][2]]
                W.append(z[0]*z[1]**2*z[2])
            err = np.abs(np.array(W)/need - 1).max()
            res.append((err, al, F(b2, 2)))
    res.sort()
    print(f"{label}: T = {Tf[0]}w + {Tf[1]};  scanned {len(res)} extra zero pairs {{e, T-e}}  [{time.time()-t0:.0f}s]")
    for err, al, b in res[:5]:
        print(f"   e = {al}w + {b}:  max |W_E / W_needed - 1| = {err:.3e}")
    return res
Ta = (6, F(9, 2)); s = (F(3), (F(9, 2) - 1)/2)
oursA = normalize([(y[0] + s[0], y[1] + s[1]) for y in SBa] + [(y[0] - s[0], y[1] - s[1]) for y in SBa], 1, Ta)
cdsA = cds_exact(Ta, 1)
scan("(a)", Ta, [(1, 1), (3, F(5, 2)), (4, 3), (6, F(9, 2))], oursA, cdsA)
Tb = (6, F(-3, 2))
oursB = normalize(BBb, int(sgn_b), Tb); cdsB = cds_exact(Tb, 0)
scan("(b)", Tb, [(1, 0), (3, F(-1, 2)), (4, -1), (6, F(-3, 2))], oursB, cdsB)
