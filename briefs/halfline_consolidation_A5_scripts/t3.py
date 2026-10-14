"""Task 2: limits delta -> 0 of the linearized K (target boundary Hessian) about regular backgrounds approaching a
boundary-singular endpoint."""
import sys, itertools
sys.path.insert(0, '/home/gustav/Git/affine_toda/textbook_code')
import mpmath as mp
from a5lib import *
src = open('/home/gustav/Git/affine_toda/textbook_code/particle_reflection.py').read()
PR = {}; exec(src.split('\n# ---------------------------------------------------------------- 1. blocks')[0], PR)
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
def Kdg(c, b, n, th): return (PR['K_uniform'](c, n)*PR['M_sol'](c, b, n))(th, 0)
blk = PR['blk']
S = lambda C: ''.join('+' if x > 0 else '-' for x in C)

def bc_patterns(bg, tol):
    h = bg.th.h; t = bg.taus(mp.mpf(0)); eh = [mp.sqrt(v) for v in bg.psi(0)]
    dphi = [t[(k-1) % h][1]/t[(k-1) % h][0] - t[k][1]/t[k][0] for k in range(h)]
    return [C for C in itertools.product([1, -1], repeat=h)
            if max(abs(dphi[k] + sum(C[j]*eh[j]*((k == j) - (k == (j+1) % h)) for j in range(h))) for k in range(h)) < tol]

def show(bg, C, Om, n, b=None):
    T = bg.th
    opn, K, Kf, res = kmatrix(bg, C, Om)
    off = max([abs(K[(x, y)]) for x in opn for y in opn if x != y] + [0])
    u = unitarity_defect(opn, Kf)
    out = f'off {mp.nstr(off, 6)} K11 {mp.nstr(K[(1,1)], 8)} K{n}1 {mp.nstr(K[(n,1)], 8)} unit {mp.nstr(u, 2)} res {mp.nstr(res, 2)}'
    if b is not None:
        out += f' |K-DG(b={b})| {mp.nstr(max(abs(K[(c,c)] - Kdg(c, b, n, mp.acosh(Om/T.mass(c)))) for c in opn), 3)}'
    return out, K

# ---- a_5, species-2 pair at chi_end = 5 pi/3 (tau_2(0) = tau_5(0) = 0)
n = 5; T = Theory(n); h = 6; a = 2; s = mp.sin(mp.pi*a/h); d = 1/(1+s)
chi_end = 5*mp.pi/3
th = mp.mpf('0.6')
for Om, lab in ((T.mass(1)*mp.cosh(th), 'Omega = m_1 cosh 0.6'), (mp.mpf('2.3'), 'Omega = 2.3 (all open)')):
    print('==== a_5 species-2 pair,', lab, '; DG b=2 K_1 =', mp.nstr(Kdg(1, 2, 5, mp.acosh(Om/T.mass(1))), 10))
    for side in (+1, -1):
        for dl in ('1e-1', '1e-2', '1e-3', '1e-4', '1e-6', '1e-8'):
            chi = chi_end + side*mp.mpf(dl)
            bg = Background(T, [(a, d*mp.exp(I*chi)), (h-a, d*mp.exp(-I*chi))])
            reg = bc_patterns(bg, mp.mpf(10)**(-14))
            assert len(reg) == 1
            t0 = [bg.tau[j](0).real for j in range(h)]
            for C in [reg[0], (1, 1, -1, -1, -1, -1), (1, 1, -1, -1, 1, -1), (1, 1, -1, 1, -1, -1), (1, 1, -1, 1, 1, -1)]:
                o, K = show(bg, C, Om, n, b=2 if C in (reg[0], (1, 1, -1, -1, -1, -1)) else None)
                print(f'side {side:+d} delta {dl} tau2(0) {mp.nstr(t0[2], 2)} {S(C)}{"(reg)" if C == reg[0] else "     "} {o}', flush=True)

# ---- self-conjugate (middle) soliton endpoint d -> 1: a_3 and a_5
for n in (3, 5):
    T = Theory(n); h = n+1
    Om = T.mass(1)*mp.cosh(th)
    targets = [(-1,)*h] + ([(1, -1, 1, -1), (1, -1, -1, -1)] if n == 3 else [(1, -1, 1, -1, -1, -1), (1, -1, -1, -1, -1, -1), (1, -1, 1, -1, 1, -1)])
    print(f'==== a_{n} middle soliton tau_j = 1 + (-1)^j d e^(2x), d -> 1, Omega = m_1 cosh 0.6; DG b={h//2} K_1 =',
          mp.nstr(Kdg(1, h//2, n, th), 10))
    for e in ('1e-1', '1e-2', '1e-3', '1e-4', '1e-6', '1e-8'):
        bg = Background(T, [(h//2, 1 - mp.mpf(e))])
        for C in targets:
            o, K = show(bg, C, Om, n, b=h//2)
            print(f'1-d {e} {S(C)} {o}', flush=True)

# ---- a_2, prod C = -1: pair endpoint chi = pi (tau_0(0) = 0)
n = 2; T = Theory(n); h = 3; a = 1; s = mp.sin(mp.pi/3); d = 1/(1+s)
Om = T.mass(1)*mp.cosh(th)
print('==== a_2 pair at chi = pi (tau_0(0)=0), Omega = m_1 cosh 0.6; DG b=1 K_1 =', mp.nstr(Kdg(1, 1, 2, th), 10))
for side in (+1, -1):
    for dl in ('1e-1', '1e-3', '1e-6'):
        chi = mp.pi + side*mp.mpf(dl)
        bg = Background(T, [(1, d*mp.exp(I*chi)), (2, d*mp.exp(-I*chi))])
        reg = bc_patterns(bg, mp.mpf(10)**(-14))
        for C in [reg[0], (-1, 1, 1), (-1, -1, -1)]:
            o, K = show(bg, C, Om, n)
            print(f'side {side:+d} delta {dl} {S(C)}{"(reg)" if C == reg[0] else "     "} {o}', flush=True)
