import numpy as np
from scal import *
def ratfit(f, xs, d):
    M = []; b = []
    for x, y in zip(xs, f):
        M.append([x**k for k in range(d+1)] + [-y*x**k for k in range(1, d+1)]); b.append(y)
    sol = np.linalg.lstsq(np.array(M), np.array(b), rcond=None)[0]
    return sol[:d+1], np.concatenate([[1], sol[d+1:]])
if __name__ == "__main__":
 for n in (2, 3, 4):
     om = 0.37; T = (n+1)*om/2
     mus = [0.05*k - 0.3 + 0.07j*((k*k) % 5 - 2) for k in range(1, 14)]
     Q = [np.exp(logR(T/2-m, n, om) - logR(T/2+m, n, om))/Y(m, n, om) for m in mus]
     ws = [np.exp(1j*np.pi*m/T) for m in mus]
     for d in range(1, 6):
         num, den = ratfit(Q[:2*d+2], ws[:2*d+2], d)
         res = max(abs(np.polyval(num[::-1], w)/np.polyval(den[::-1], w) - q) for w, q in zip(ws, Q))
         if res < 1e-6:
             zs = np.roots(num[::-1]); ps = np.roots(den[::-1])
             # w0 = exp(i pi mu0/T) -> mu0 = T log(w0)/(i pi)
             f = lambda w: np.round(T*np.log(w+0j)/(1j*np.pi), 5)
             print(n, 'deg', d, 'res %.1e' % res, 'zeros mu:', [f(z) for z in zs], 'poles mu:', [f(p) for p in ps], 'T', T, 'om', om)
             break
     else:
         print(n, 'no fit')
