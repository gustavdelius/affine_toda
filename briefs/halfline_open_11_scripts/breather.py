import numpy as np
from scal import amps, ST, params, blk
from bproj import cproj_sp
from t14 import ratfit

def KB(mu, n, om, l, eta, cdd=None, conj=False):
    T = (n+1)*om/2
    th = 1j*np.pi*mu/T
    c, r, rs = cproj_sp(n, om, l, eta, th, conj=conj)
    ma, mb = mu - (T-1)/2, mu + (T-1)/2
    a_a, b_a = amps(ma, n, om, l, eta, cdd); a_b, b_b = amps(mb, n, om, l, eta, cdd)
    F = ST(T - 2*mu, n, om)
    return (a_a*b_b if not conj else b_a*a_b)*F*c

def blockfit(f, T, mus, dmax=8, tol=1e-6):
    ws = [np.exp(1j*np.pi*m/T) for m in mus]
    for d in range(0, dmax+1):
        if 2*d+2 > len(mus): break
        num, den = ratfit(f[:2*d+2], ws[:2*d+2], d)
        res = max(abs(np.polyval(num[::-1], w)/np.polyval(den[::-1], w) - q) for w, q in zip(ws, f))
        if res < tol:
            g = lambda w: np.round(T*np.log(w+0j)/(1j*np.pi), 4)
            return d, res, [g(z) for z in np.roots(num[::-1])], [g(p) for p in np.roots(den[::-1])], (num[-1]/den[-1] if d else num[0])
    return None

if __name__ == '__main__':
    import sys
    om = 0.37
    for (n, l) in [(1, 1), (2, 1), (2, 0), (3, 1), (3, 2), (3, 0)]:
        T = (n+1)*om/2
        mus = [0.04*k - 0.35 + 0.05j*((k*k) % 7 - 3) for k in range(1, 22)]
        for eta in (0.3+0.4j,):
            f1 = [KB(m, n, om, l, eta) for m in mus]
            f2 = [KB(m, n, om, l, eta, conj=True) for m in mus]
            pr = [x*y for x, y in zip(f1, f2)]
            print(n, l, 'T', T, 'om', om, 'eps', np.round(params(n, om, l, eta)[1], 4), 'epsb', np.round(params(n, om, l, eta)[2], 4))
            for name, f in (('B->Bbar', f1), ('Bbar->B', f2), ('product', pr)):
                print('   ', name, blockfit(f, T, mus))
