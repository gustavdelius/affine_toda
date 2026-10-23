import numpy as np, os
from scipy.special import loggamma
from allS import *                     # c_3 fusion machinery (S_33 = c f R, omega = 2.37)
def Fmin_maker(A, T, J=2000):
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
blk = lambda y, th: np.sinh(th/2 + 1j*np.pi*y/(2*T))/np.sinh(th/2 - 1j*np.pi*y/(2*T))
ths = [0.37+0.21j, -0.8+1.3j, 1.1-0.4j, 0.2+2.2j]
Phi11 = lambda th: K_top(1, 1, th)[0]
for lifts in ([0, 0.5, 2.5, 3.0], [0, 0.5, 3.5, 3.0], [0, -0.5, 2.5, 3.0], [1, 0.5, 2.5, 3.0]):
    A = [omega + lifts[0], omega + lifts[1], 3*omega + lifts[2] - 2, 3*omega + lifts[3] - 2]   # zeros at w, w+1/2 (+/-q^-2), 3w+1/2, T=3w+1 (+/-q^-6)
    A = [omega + lifts[0], omega + lifts[1], 3*omega + lifts[2] - 2.0, 3*omega + 1.0]
    Fm = Fmin_maker(A, T)
    X = np.array([Phi11(th)/Fm(th) for th in ths])
    print(f"zeros {np.round(A, 3)}:  fused Phi_11 / minimal F = {np.round(X, 4)}")
