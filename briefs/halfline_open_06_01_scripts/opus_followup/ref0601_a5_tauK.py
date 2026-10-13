"""Referee check (2026-10-12), 06_01 item 1: linearized reflection matrix about the a_5 species-2 stationary pair
near and at the singular point chi_end = 5pi/3 (tau_2(0)=tau_5(0)=0), computed EXACTLY by the tau-function method
(BP sec. 5): each linear wave of channel c (exponent e^{sigma x - i Omega t}) is a "soliton" of species c dressed by the
static background solitons (species 2,4) with the Hirota interaction
   A(c,sigma; b) = [sigma/m_c - cos(pi(c-b)/h)] / [sigma/m_c - cos(pi(c+b)/h)]   (static-static limit = worker's A_ab),
   dtau_j = omega^{jc} e^{sigma x} [1 + A_c2 c2 w^{2j} E + A_c4 c4 w^{4j} E + A_c2 A_c4 A_24 c2 c4 E^2],  E = e^{sqrt3 x}.
dpsi_j = dtau_{j-1}/tau_{j-1} + dtau_{j+1}/tau_{j+1} - 2 dtau_j/tau_j. The BC (book): psi_k'(0) = -sum_j C_j (alpha_k.alpha_j) e^{psi_j/2},
linearized: dpsi_k'(0) + sum_j C_j (alpha_k.alpha_j) e^{psi_j(0)/2} dpsi_j(0)/2 = 0.  Everything at x=0 in mpmath (dps 50):
no ODE, so no ill-conditioning from the wall layer.  K in the phi-channel basis (channel c = Fourier vector omega^{kc}).
Self-checks: linearized field equation residual of each dressed wave; regular pattern -> BP/DG diagonal values."""
import mpmath as mp, itertools, sys
mp.mp.dps = 50
h = 6; PI = mp.pi; w = mp.exp(2j*PI/h)
m = [None] + [2*mp.sin(PI*a/h) for a in range(1, h)]
s3 = mp.sqrt(3); d = 4 - 2*s3
A24 = mp.sin(PI*(2-4)/(2*h))**2/mp.sin(PI*(2+4)/(2*h))**2

def bg(chi):
    c2 = d*mp.exp(1j*chi); c4 = d*mp.exp(-1j*chi)
    return c2, c4

def tau_bg(j, x, chi):
    c2, c4 = bg(chi); E = mp.exp(s3*x)
    return 1 + c2*w**(2*j)*E + c4*w**(4*j)*E + A24*c2*c4*E**2

def Aint(c, sigma, b):
    return (sigma/m[c] - mp.cos(PI*(c-b)/h))/(sigma/m[c] - mp.cos(PI*(c+b)/h))

def dtau(j, x, chi, c, sigma):
    c2, c4 = bg(chi); E = mp.exp(s3*x)
    a2 = Aint(c, sigma, 2); a4 = Aint(c, sigma, 4)
    return w**(j*c)*mp.exp(sigma*x)*(1 + a2*c2*w**(2*j)*E + a4*c4*w**(4*j)*E + a2*a4*A24*c2*c4*E**2)

def psi_bg(j, x, chi):
    return mp.log(tau_bg((j-1) % h, x, chi)*tau_bg((j+1) % h, x, chi)/tau_bg(j, x, chi)**2)

def dpsi(j, x, chi, c, sigma):
    f = lambda jj: dtau(jj % h, x, chi, c, sigma)/tau_bg(jj % h, x, chi)
    return f(j-1) + f(j+1) - 2*f(j)

def check_wave(chi, c, sigma, Om, x=mp.mpf('-0.3')):
    """linearized Toda: dpsi_j'' + Om^2 dpsi_j = 2 e^{psi_j} dpsi_j - e^{psi_{j-1}} dpsi_{j-1} - e^{psi_{j+1}} dpsi_{j+1}"""
    res = 0
    for j in range(h):
        lhs = mp.diff(lambda xx: dpsi(j, xx, chi, c, sigma), x, 2) + Om**2*dpsi(j, x, chi, c, sigma)
        e = lambda jj: mp.exp(psi_bg(jj % h, x, chi))*dpsi(jj % h, x, chi, c, sigma)
        rhs = 2*e(j) - e(j-1) - e(j+1)
        res = max(res, abs(lhs - rhs)/max(1, abs(lhs)))
    return res

def bc_rows(chi, C, c, sigma):
    """linearized BC residual vector (k=0..5) for the wave (c,sigma) at x=0"""
    x0 = mp.mpf(0)
    eh = [mp.exp(psi_bg(j, x0, chi)/2) for j in range(h)]
    dp = [dpsi(j, x0, chi, c, sigma) for j in range(h)]
    dpp = [mp.diff(lambda xx: dpsi(j, xx, chi, c, sigma), x0) for j in range(h)]
    def ad(k, j): return 2 if k == j else (-1 if (k-j) % h in (1, h-1) else 0)
    return [dpp[k] + sum(C[j]*ad(k, j)*eh[j]*dp[j]/2 for j in range(h)) for k in range(h)]

def Kmatrix(chi, C, Om):
    open_ = [c for c in range(1, h) if Om > m[c]]; closed = [c for c in range(1, h) if Om <= m[c]]
    p = {c: mp.sqrt(Om**2 - m[c]**2) for c in open_}; kap = {c: mp.sqrt(m[c]**2 - Om**2) for c in closed}
    unk = [(c, -1j*p[c]) for c in open_] + [(c, kap[c]) for c in closed]          # out-waves and decaying waves
    Mcols = [bc_rows(chi, C, c, sg) for c, sg in unk]
    Mx = mp.matrix(h, len(unk))
    for i, col in enumerate(Mcols):
        for k in range(h): Mx[k, i] = col[k]
    K = mp.matrix(len(open_), len(open_))
    # phi-channel normalisation: phi-amplitude = -m_c^2 sqrt(h)/(1-w^c) * tau-wave amplitude
    nf = lambda c: -m[c]**2/(1 - w**c)
    for ia, a in enumerate(open_):
        r = bc_rows(chi, C, a, 1j*p[a])
        b = mp.matrix([-x for x in r])
        sol, resid = mp.qr_solve(Mx, b)
        for ib, bb in enumerate(open_):
            K[ib, ia] = sol[ib]*nf(bb)/nf(a)
    return K, open_, resid

def blk(x, th): return mp.sinh(th/2 + 1j*PI*x/(2*h))/mp.sinh(th/2 - 1j*PI*x/(2*h))
def K_bp_diag(a, th):
    if a in (1, 5): return blk(2, th)*blk(4, th)/(blk(1, th)*blk(5, th))
    if a in (2, 4): return -blk(1, th)*blk(3, th)**2*blk(5, th)/(blk(2, th)*blk(4, th))
    return -blk(2, th)**2*blk(4, th)**2/blk(3, th)**2

def s(C): return ''.join('+' if c > 0 else '-' for c in C)

if __name__ == '__main__':
    chi_end = 5*PI/3
    reg = {+1: (1, -1, -1, 1, -1, -1), -1: (-1, 1, -1, -1, 1, -1)}
    targets = [(1, 1, -1, -1, -1, -1), (1, 1, -1, -1, 1, -1), (1, 1, -1, 1, -1, -1), (1, 1, -1, 1, 1, -1)]
    Om0 = mp.cosh(mp.mpf('0.6'))
    print('wave self-check (linearized Toda residual) at chi=chi_end+0.1, Omega=cosh 0.6:',
          [mp.nstr(check_wave(chi_end+mp.mpf('0.1'), c, sg, Om0), 3) for c, sg in
           [(1, 1j*mp.sqrt(Om0**2-1)), (1, -1j*mp.sqrt(Om0**2-1)), (2, mp.sqrt(3-Om0**2)), (3, mp.sqrt(4-Om0**2))]])
    deltas = [mp.mpf(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [mp.mpf(x) for x in ('0.1', '1e-2', '1e-3', '1e-5', '1e-8')]
    for Om in (Om0, mp.mpf('2.3')):
        open_ = [c for c in range(1, h) if Om > m[c]]
        ths = {a: mp.acosh(Om/m[a]) for a in open_}
        print('\n== Omega=%s open %s  BP diag: %s' % (mp.nstr(Om, 6), open_, {a: mp.nstr(K_bp_diag(a, ths[a]), 8) for a in open_}))
        for side in (+1, -1):
            for C in [reg[side]] + targets:
                for dl in deltas:
                    chi = chi_end + side*dl
                    K, op, resid = Kmatrix(chi, C, Om)
                    nO = len(op)
                    off = max([abs(K[i, j]) for i in range(nO) for j in range(nO) if i != j] + [0])
                    dg = max(abs(K[i, i] - K_bp_diag(op[i], ths[op[i]])) for i in range(nO))
                    Kd = mp.matrix(K); uni = mp.norm(Kd.H*Kd - mp.eye(nO), p=mp.inf)
                    print('side %+d %s%s delta=%s tau2(0)=%s | max|offdiag|=%s max|Kdiag-BP|=%s unitarity %s lsq-resid %s' % (
                        side, s(C), '(reg)' if C == reg[side] else '     ', mp.nstr(dl, 2), mp.nstr(abs(tau_bg(2, 0, chi)), 2),
                        mp.nstr(off, 3), mp.nstr(dg, 3), mp.nstr(uni, 2), mp.nstr(resid, 2)), flush=True)
                    if C == reg[side] and dl < mp.mpf('0.05') and side == -1: pass
