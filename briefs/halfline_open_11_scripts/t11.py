import numpy as np
from scal import *
n, om, l, eta = 2, 0.37, 1, 0.3+0.4j
T, eps, epsb, lam = params(n, om, l, eta)
# test F11 shift identity
x = 0.17+0.08j
print('F11 shift', F11(x-2*T, n+1, om)/F11(x, n+1, om)/(np.sin(np.pi*x)*np.sin(np.pi*(x-2*T))/(np.sin(np.pi*(x-2*T+om))*np.sin(np.pi*(x-om)))))
print('ST crossing unitarity', ST(T-x, n, om)*ST(T+x, n, om))
for nu in (0.11+0.05j, -0.2+0.13j):
    w = np.exp(logRt(nu+2*T, n, om) - logRt(nu, n, om))
    wt = Y(nu+T/2, n, om)/Y(nu+1.5*T, n, om)
    print('R shift ratio / W', w/wt)
    ga = np.exp(logg_a(nu+2*T, n, om, l, eta) - logg_a(nu, n, om, l, eta))
    Ub = lam**2*U(nu+T, eta + 1j*np.pi*(epsb-eps))
    print('g_a shift / (U/Ub(nu+T))', ga/(U(nu, eta)/Ub))
