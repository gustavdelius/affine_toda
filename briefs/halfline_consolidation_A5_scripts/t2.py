"""Task 1 sweep: every regular chi-interval of the continued pair, a_2..a_5, vs DG (5.5) at B=0."""
import sys, io, contextlib, itertools, time
sys.path.insert(0, '/home/gustav/Git/affine_toda/code')
import numpy as np
import mpmath as mp
from a5lib import *
with contextlib.redirect_stdout(io.StringIO()):
    import boundary_vacua as BV
src = open('/home/gustav/Git/affine_toda/code/particle_reflection.py').read()
PR = {}; exec(src.split('\n# ---------------------------------------------------------------- 1. blocks')[0], PR)
mp.mp.dps = 30

def Kdg(c, b, n, th): return (PR['K_uniform'](c, n)*PR['M_sol'](c, b, n))(th, 0)

def dg_class(C):
    h = len(C); n = h-1
    Cdg = [C[n]] + [C[i-1] for i in range(1, h)]
    ks = [k for k in range(1, h) if Cdg[k] == -1]
    a = len(ks); b = sum((-1)**(a-i)*k for i, k in enumerate(ks, 1))
    return min(b, h-b)

def bc_pattern(bg):
    h = bg.th.h; t = bg.taus(mp.mpf(0)); eh = [mp.sqrt(v) for v in bg.psi(0)]
    dphi = [t[(k-1) % h][1]/t[(k-1) % h][0] - t[k][1]/t[k][0] for k in range(h)]
    return [C for C in itertools.product([1, -1], repeat=h)
            if max(abs(dphi[k] + sum(C[j]*eh[j]*((k == j) - (k == (j+1) % h)) for j in range(h))) for k in range(h)) < mp.mpf(10)**(-25)]

worst = {'off': 0, 'diag': 0, 'unit': 0, 'res': 0, 'pde': 0, 'bv': 0}
t0 = time.time()
for n in (2, 3, 4, 5):
    T = Theory(n); h = n+1
    for a in range(1, (h+1)//2 + (1 if h % 2 == 0 else 0)):
        s = mp.sin(mp.pi*a/h)
        sing = sorted(set(round(float((mp.pi - 2*mp.pi*j*a/h) % (2*mp.pi)), 12) for j in range(h)))
        if 2*a == h:
            cases = [('d', mp.mpf(d)) for d in ('0.3', '0.7', '0.95', '-0.5')]
        else:
            gaps = [(sing[k], sing[k+1] if k+1 < len(sing) else sing[0]+2*np.pi) for k in range(len(sing))]
            cases = [('chi', mp.mpf(lo) + mp.mpf(fr)*(mp.mpf(hi)-mp.mpf(lo))) for lo, hi in gaps for fr in ('0.17', '0.5', '0.81')]
        for kind, par in cases:
            if kind == 'chi':
                d = 1/(1+s); sols = [(a, d*mp.exp(I*par)), (h-a, d*mp.exp(-I*par))]
                chi = float(par)
            else:
                sols = [(a, par)]; chi = float(mp.acos(par))
            bg = Background(T, sols)
            # same tau as boundary_vacua.pair_taus (chi or d = cos chi for a = h/2)
            ref = BV.pair_taus(n, a, chi, np.array([-0.7, -0.2, 0.0]))
            worst['bv'] = max(worst['bv'], max(abs(complex(bg.tau[j](mp.mpf(x))) - ref[j][0][i]) for j in range(h) for i, x in enumerate((-0.7, -0.2, 0.0))))
            pats = bc_pattern(bg)
            assert len(pats) == 1, pats
            C = pats[0]
            # boundary_vacua: pattern of the soliton charge
            if kind == 'chi':
                lam = BV.as_weight(BV.soliton_charge(n, a, chi), h)
                assert BV.pattern(np.array(lam, float), BV.book_roots(n)) == C and sum(lam) == a
            b = dg_class(C)
            assert b == min(a, h-a), (C, b, a)
            line = []
            for lab, Om in (('th1', T.mass(1)*mp.cosh(mp.mpf('0.6'))), ('th2', T.mass(2)*mp.cosh(mp.mpf('0.6'))), ('all', mp.mpf('2.3'))):
                opn, K, Kf, res = kmatrix(bg, C, Om)
                off = max([abs(K[(x, y)]) for x in opn for y in opn if x != y] + [0])
                dd = max(abs(K[(c, c)] - Kdg(c, b, n, mp.acosh(Om/T.mass(c)))) for c in opn)
                u = unitarity_defect(opn, Kf)
                c0 = opn[0]
                pr = max(pde_residual(bg, bg.wave(c, sg), mp.mpf(x), Om) for c, sg in ((c0, I*mp.sqrt(Om**2-T.mass(c0)**2)), (c0, -I*mp.sqrt(Om**2-T.mass(c0)**2))) for x in ('-0.37', '-1.9'))
                for k_, v in (('off', off), ('diag', dd), ('unit', u), ('res', res), ('pde', pr)):
                    worst[k_] = max(worst[k_], v)
                line.append(f'{lab}: open {opn} off {mp.nstr(off,2)} |K-DG| {mp.nstr(dd,2)}')
            print(f'a_{n} a={a} {kind}={mp.nstr(par,6)} C={"".join("+" if x>0 else "-" for x in C)} b={b} | ' + ' ; '.join(line), flush=True)
print({k: mp.nstr(v, 3) for k, v in worst.items()}, 'time', round(time.time()-t0, 1))
