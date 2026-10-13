import numpy as np
from breather import KB
from scal import params
def lsfit(f, ws, d):
    M = np.array([[w**k for k in range(d+1)] + [-y*w**k for k in range(1, d+1)] for w, y in zip(ws, f)])
    sol = np.linalg.lstsq(M, np.array(f), rcond=None)[0]
    num, den = sol[:d+1], np.concatenate([[1], sol[d+1:]])
    res = max(abs(np.polyval(num[::-1], w)/np.polyval(den[::-1], w) - y)/abs(y) for w, y in zip(ws, f))
    return num, den, res
def show(f, T, mus, dmax=10, tol=1e-5, label=''):
    ws = [np.exp(1j*np.pi*m/T) for m in mus]
    for d in range(0, dmax+1):
        num, den, res = lsfit(f, ws, d)
        if res < tol:
            g = lambda w: np.round(T*np.log(w+0j)/(1j*np.pi), 4)
            print(label, 'deg', d, 'res %.1e' % res, 'zeros', sorted([g(z) for z in np.roots(num[::-1])], key=lambda z: z.real),
                  'poles', sorted([g(p) for p in np.roots(den[::-1])], key=lambda z: z.real), 'lead', np.round(num[-1]/den[-1], 5))
            return
    print(label, 'no fit')
om = 0.37
mus = [0.035*k - 0.4 + 0.04j*((k*k) % 7 - 3) for k in range(1, 30)]
for (n, l) in [(1, 1), (2, 1), (2, 0), (3, 1), (3, 2)]:
    T = (n+1)*om/2
    e1, e2 = 0.3+0.4j, -0.1+0.7j
    for conj in (False, True):
        f = [KB(m, n, om, l, e1, conj=conj)/KB(m, n, om, l, e2, conj=conj) for m in mus]
        show(f, T, mus, label=f'{n} {l} conj={conj} eta-ratio T={T:.3f} eps1={params(n,om,l,e1)[1]:.3f} eps2={params(n,om,l,e2)[1]:.3f} epsb1={params(n,om,l,e1)[2]:.3f}')
