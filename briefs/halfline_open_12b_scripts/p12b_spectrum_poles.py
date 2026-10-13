"""Part (b): spectrum rule and pole bookkeeping.
Generates all boundary states |n1,...,nk> allowed by  pi/2 >= nu_{n1} > w_{n2} > nu_{n3} > ... > 0  (BPTT sec. 3.8; prp-bsg-spectrum), their energies
M sum cos(.), the sector (parity of the sum of labels at vt = 0, BPT sec. 4.3.2), and, for the closed-form soliton factor
Q^S = R_0 sigma(x_S) Q_0 prod a_ni, the exact order of every pole in the physical strip 0 < u < pi/2 (p12b_poles.py, checked against numerical
orders in p12b_poles_check.py).  Each strip pole is classed
   N : position of a NEW state S+(m)  (nu_m for even depth, w_m for odd depth, below the last label's position)
   L : at a position nu_j or w_j with j >= 0 (the families discussed by BPTT: sigma(eta), sigma(etabar), numerators of a_n)
   U : anything else, in particular nu_{-j}, j >= 1 (the mirror tower of the denominators of a_n(etabar) and the b-type poles of sigma(etabar))
usage:  python3 p12b_spectrum_poles.py lam eta_over_etaN maxdepth
"""
import sys, math
from p12b_poles import total_events, label_position

def states(eta, lam, maxdepth):
    nu_ = lambda j: eta/lam - (2*j+1)*math.pi/(2*lam)
    etab = math.pi*(lam+1) - eta
    w_ = lambda j: etab/lam - (2*j+1)*math.pi/(2*lam)
    out = [((), None)]
    def rec(S, bound, depth):
        if depth == maxdepth: return
        k = len(S); f = nu_ if k % 2 == 0 else w_
        for m in range(0, 400):
            p = f(m)
            if p <= 1e-9: break
            if p <= math.pi/2 + 1e-12 and (k == 0 or p < bound - 1e-12):
                out.append((S + (m,), p)); rec(S + (m,), p, depth + 1)
    rec((), math.pi/2 + 1e-12, 0)
    return out, nu_, w_

def classify(S, eta, lam, ev, nu_, w_):
    k = len(S); f = nu_ if k % 2 == 0 else w_
    last = None
    if S:
        g = w_ if k % 2 == 0 else nu_      # type of the last label
        last = g(S[-1])
    rows = []
    for u, o in sorted(ev.items(), reverse=True):
        lab = label_position(u, eta, lam)
        new = any(abs(f(m) - u) < 1e-7 for m in range(0, 400)) and (last is None or u < last - 1e-12) and u > 0
        if new: cls = 'N'
        elif any(abs(nu_(j) - u) < 1e-7 or abs(w_(j) - u) < 1e-7 for j in range(0, 400)): cls = 'L'
        else: cls = 'U'
        rows.append((lab, o, cls))
    return rows

if __name__ == '__main__':
    lam = float(sys.argv[1]); frac = float(sys.argv[2]); depth = int(sys.argv[3])
    eta0 = math.pi*(lam+1)/2; eta = frac*eta0
    sts, nu_, w_ = states(eta, lam, depth)
    print(f'lambda={lam}, eta={eta/math.pi:.4f} pi = {frac} eta_N, eta_N = {eta0/math.pi:.4f} pi;  strip edge pi/2')
    print('nu_j (j=-3..6):', [round(nu_(j)/math.pi, 4) for j in range(-3, 7)], ' (units of pi)')
    print('w_j  (j=-3..9):', [round(w_(j)/math.pi, 4) for j in range(-3, 10)])
    for S, p in sts:
        k = len(S)
        E = 0.0
        for i, n in enumerate(S, start=1):
            E += math.cos(nu_(n) if i % 2 == 1 else w_(n))
        par = sum(S) % 2
        ev = total_events(S, eta, lam)
        rows = classify(S, eta, lam, ev, nu_, w_)
        cnt = {c: sum(1 for r in rows if r[2] == c) for c in 'NLU'}
        print(f'|{",".join(map(str, S))}>  E/M={E:.4f} sector={"odd" if par else "even"}  poles: ' +
              ', '.join(f'{lab}:{int(o) if abs(o-round(o))<1e-9 else o}{c}' for lab, o, c in rows))
