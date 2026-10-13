"""Referee (2026-10-12): same exact tau-method as ref0601_a5_tauK.py, for the a_5 middle-soliton family
tau_j = 1 + (-1)^j d e^{2x} (regular for d<1, exact pattern C=-1 everywhere), approaching d->1 where tau_odd(0)=0.
Target patterns: C=-1 (control), +-+--- (missing, b=3), +----- and +-+-+- (prod C=-1); BP Table 2 gives the C=-1 values
K_1=-1/((1)(5)), K_2=-1/((2)(4)), K_3=-1/(3)^2 for all four (row 1 of the r=5 block, with A=-C)."""
import mpmath as mp, sys
mp.mp.dps = 50
h = 6; PI = mp.pi; w = mp.exp(2j*PI/h); m = [None] + [2*mp.sin(PI*a/h) for a in range(1, h)]
def Aint(c, sigma, b): return (sigma/m[c] - mp.cos(PI*(c-b)/h))/(sigma/m[c] - mp.cos(PI*(c+b)/h))
def tau_bg(j, x, d): return 1 + (-1)**j*d*mp.exp(2*x)
def dtau(j, x, d, c, sg): return w**(j*c)*mp.exp(sg*x)*(1 + Aint(c, sg, 3)*(-1)**j*d*mp.exp(2*x))
def psi_bg(j, x, d): return mp.log(tau_bg((j-1) % h, x, d)*tau_bg((j+1) % h, x, d)/tau_bg(j, x, d)**2)
def dpsi(j, x, d, c, sg):
    f = lambda jj: dtau(jj % h, x, d, c, sg)/tau_bg(jj % h, x, d)
    return f(j-1) + f(j+1) - 2*f(j)
def check_wave(d, c, sg, Om, x=mp.mpf('-0.3')):
    r = 0
    for j in range(h):
        lhs = mp.diff(lambda xx: dpsi(j, xx, d, c, sg), x, 2) + Om**2*dpsi(j, x, d, c, sg)
        e = lambda jj: mp.exp(psi_bg(jj % h, x, d))*dpsi(jj % h, x, d, c, sg)
        r = max(r, abs(lhs - (2*e(j) - e(j-1) - e(j+1))))
    return r
def bc_rows(d, C, c, sg):
    x0 = mp.mpf(0); eh = [mp.exp(psi_bg(j, x0, d)/2) for j in range(h)]
    dp = [dpsi(j, x0, d, c, sg) for j in range(h)]; dpp = [mp.diff(lambda xx: dpsi(j, xx, d, c, sg), x0) for j in range(h)]
    ad = lambda k, j: 2 if k == j else (-1 if (k-j) % h in (1, h-1) else 0)
    return [dpp[k] + sum(C[j]*ad(k, j)*eh[j]*dp[j]/2 for j in range(h)) for k in range(h)]
def Kmatrix(d, C, Om):
    op = [c for c in range(1, h) if Om > m[c]]; cl = [c for c in range(1, h) if Om <= m[c]]
    p = {c: mp.sqrt(Om**2-m[c]**2) for c in op}; kap = {c: mp.sqrt(m[c]**2-Om**2) for c in cl}
    unk = [(c, -1j*p[c]) for c in op] + [(c, kap[c]) for c in cl]
    Mx = mp.matrix(h, len(unk))
    for i, (c, sg) in enumerate(unk):
        col = bc_rows(d, C, c, sg)
        for k in range(h): Mx[k, i] = col[k]
    nf = lambda c: -m[c]**2/(1 - w**c); K = mp.matrix(len(op), len(op)); res = 0
    for ia, a in enumerate(op):
        sol, r = mp.qr_solve(Mx, mp.matrix([-x for x in bc_rows(d, C, a, 1j*p[a])])); res = max(res, r)
        for ib, b in enumerate(op): K[ib, ia] = sol[ib]*nf(b)/nf(a)
    return K, op, res
def blk(x, th): return mp.sinh(th/2 + 1j*PI*x/(2*h))/mp.sinh(th/2 - 1j*PI*x/(2*h))
def Kflat_m(a, th): return -1/(blk(a, th)*blk(h-a, th))
s = lambda C: ''.join('+' if c > 0 else '-' for c in C)
if __name__ == '__main__':
    Om0 = mp.cosh(mp.mpf('0.6'))
    print('wave self-check:', [mp.nstr(check_wave(mp.mpf('0.7'), c, sg, Om0), 3) for c, sg in [(1, 1j*mp.sqrt(Om0**2-1)), (3, mp.sqrt(4-Om0**2))]])
    for Om in (Om0, mp.mpf('2.3')):
        op = [c for c in range(1, h) if Om > m[c]]; ths = {a: mp.acosh(Om/m[a]) for a in op}
        print('\n== Omega=%s open %s; C=-1 values %s' % (mp.nstr(Om, 6), op, {a: mp.nstr(Kflat_m(a, ths[a]), 8) for a in op}))
        for C in [(-1,)*6, (1, -1, 1, -1, -1, -1), (1, -1, -1, -1, -1, -1), (1, -1, 1, -1, 1, -1)]:
            for eta in ('0.3', '1e-2', '1e-4', '1e-6'):
                d = 1 - mp.mpf(eta); K, o, r = Kmatrix(d, C, Om); n = len(o)
                off = max([abs(K[i, j]) for i in range(n) for j in range(n) if i != j] + [0])
                dg = max(abs(K[i, i] - Kflat_m(o[i], ths[o[i]])) for i in range(n))
                print('%s 1-d=%s max|offdiag|=%s max|Kdiag - C=-1 value|=%s lsq-resid %s' % (s(C), eta, mp.nstr(off, 3), mp.nstr(dg, 3), mp.nstr(r, 2)), flush=True)
