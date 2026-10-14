"""B -> 0 limits of the continued SP amplitude L_1 and the continued CDD contributions in book blocks."""
import mpmath as mp
from t5_particles import make_L, blk, CDD, I, PI
mp.mp.dps = 30
ths = [mp.mpc('0.37', '0.11'), mp.mpc('-0.83', '0.29')]
for n in (1, 2, 3, 4):
    h = n + 1
    H = mp.mpf(h)
    for B in (mp.mpf('0.317'), mp.mpf('1.43')):
        base = make_L(n, B, 0, mp.mpc('0.7', '0.2'), 'none')[1]
        sA = make_L(n, B, 0, mp.mpc('0.7', '0.2'), 'sigma_A')[1]
        sB = make_L(n, B, 0, mp.mpc('0.7', '0.2'), 'sigma_B')[1]
        eA = max(abs(sA(t)/base(t) - (blk(B/2, t, h)/blk(H + B/2, t, h))**2) for t in ths)
        eB = max(abs(sB(t)/base(t) - blk(B/2, t, h)**2*blk(B/2 + 1, t, h)*blk(B/2 - 1, t, h)
                     / (blk(H + B/2, t, h)**2*blk(H - 1 + B/2, t, h)*blk(H + 1 + B/2, t, h))) for t in ths)
        print(f'n={n} B={B}: sigma_A contribution = (B/2)^2/(h+B/2)^2 [{mp.nstr(eA, 2)}]; sigma_B contribution = '
              f'(B/2)^2(B/2+1)(B/2-1)/((h+B/2)^2(h-1+B/2)(h+1+B/2)) [{mp.nstr(eB, 2)}]')
    B = mp.mpf('1e-12')
    for cdd in CDD:
        for l, phi in ((0, mp.mpc('0.7', '0.2')), (1, mp.pi/2)):
            L1 = make_L(n, B, l, phi, cdd)[1]
            vals = [L1(t) for t in ths]
            s = [mp.sinh(t) for t in ths]
            g = lambda x, a: (x + I*mp.sin(PI*a/h))/(x - I*mp.sin(PI*a/h))
            E0 = [(x - I*mp.sin(phi))/(x + I*mp.sin(phi + 2*PI*l/h - PI)) for x in s]
            pred = [g(x, H/2)*g(x, H/2 - 1)*e for x, e in zip(s, E0)]
            print(f'   n={n} cdd={cdd:15s} l={l}: L_1(B->0) at two theta: {mp.nstr(vals[0], 8)}, {mp.nstr(vals[1], 8)};'
                  f' ratio to g_(h/2) g_(h/2-1) E: {mp.nstr(vals[0]/pred[0], 8)}, {mp.nstr(vals[1]/pred[1], 8)}')
