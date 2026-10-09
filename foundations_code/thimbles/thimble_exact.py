"""Exact thimble integrals for n=2 in asymptotic linear coordinates, integrating the deviation
delta = q - sigma with expm1 so that tiny initial deviations keep full relative precision.
Chart: q(c) = sigma + delta_T(c), delta_0 = V e^{-D T} c, flowed upward for time T (speed cap gcap)."""
import numpy as np

def dgrad(ch, sig, D):
    """grad S(sigma + D) - grad S(sigma), D (B, L, r), accurate for tiny D."""
    P = np.concatenate([np.zeros_like(D[:, :1]), D, np.zeros_like(D[:, :1])], axis=1)
    lap = 2*D - P[:, :-2] - P[:, 2:]
    Es = np.exp(1j*(sig @ ch.alpha.T))  # (L,h)
    Em = np.expm1(1j*np.einsum('jr,bxr->bxj', ch.alpha, D))
    gW = -1j*np.einsum('bxj,jr->bxr', Es[None]*Em*ch.kac, ch.alpha)
    return lap + ch.kappa*gW

def field(ch, sig, Dflat, sgn, gcap):
    B = Dflat.shape[0]
    G = ch.phase*dgrad(ch, sig, Dflat.reshape(B, ch.L, ch.r)).reshape(B, -1)
    F = sgn*np.conj(G)
    nrm = np.linalg.norm(F, axis=1)
    fac = np.where(nrm > gcap, gcap/np.maximum(nrm, 1e-300), 1.0)
    return F*fac[:, None]

def rk4(ch, sig, D, dt, sgn, gcap):
    k1 = field(ch, sig, D, sgn, gcap); k2 = field(ch, sig, D+0.5*dt*k1, sgn, gcap)
    k3 = field(ch, sig, D+0.5*dt*k2, sgn, gcap); k4 = field(ch, sig, D+dt*k3, sgn, gcap)
    return D + dt/6*(k1+2*k2+2*k3+k4)

def thimble_integral(ch, sig, hbar, V, mu, nrho=32, nphi=64, dt=0.005, gcap=50.0, start_eps=1e-8, rho_cut=40.0):
    assert ch.n == 2
    s0 = sig.reshape(-1)
    Ss = ch.phase*ch.action(sig)
    N = (2*np.pi*hbar)**(-ch.r*(ch.L+1)/2)
    rho_max = np.sqrt(2*rho_cut*hbar/mu.min())
    T = np.log(rho_max/start_eps)/mu.min()
    x, w = np.polynomial.legendre.leggauss(nrho)
    rho = 0.5*rho_max*(x+1); wr = 0.5*rho_max*w
    ph = np.linspace(0, 2*np.pi, nphi, endpoint=False)
    R, P = np.meshgrid(rho, ph, indexing='ij')
    C = np.stack([R*np.cos(P), R*np.sin(P)], -1).reshape(-1, 2)
    hc = 1e-5
    Cs = np.concatenate([C, C+[hc, 0], C-[hc, 0], C+[0, hc], C-[0, hc]])
    D = (Cs*np.exp(-mu*T)[None, :]) @ V.T
    nst = int(np.ceil(T/dt)); h = T/nst
    M = C.shape[0]
    active = np.ones(M, bool)
    with np.errstate(all='ignore'):
        for k in range(nst):
            idx = np.where(active)[0]
            if len(idx) == 0: break
            sel = np.concatenate([idx + j*M for j in range(5)])
            D[sel] = rk4(ch, sig, D[sel], h, +1, gcap)
            if k % 10 == 0:
                qa = s0[None, :] + D[idx]
                dSa = (ch.phase*ch.action(qa.reshape(len(idx), ch.L, ch.r)) - Ss).real/hbar
                active[idx[(dSa > rho_cut + 20) | ~np.isfinite(dSa)]] = False
    d1 = (D[M:2*M]-D[2*M:3*M])/(2*hc); d2 = (D[3*M:4*M]-D[4*M:])/(2*hc)
    Jc = d1[:, 0]*d2[:, 1]-d1[:, 1]*d2[:, 0]
    q0 = s0[None, :] + D[:M]
    dS = ch.phase*ch.action(q0.reshape(M, ch.L, ch.r)) - Ss
    f = np.exp(-dS/hbar)*Jc
    f[~active] = 0
    f[~np.isfinite(f)] = 0
    f = f.reshape(nrho, nphi)
    I = np.sum(f*(wr*rho)[:, None])*(2*np.pi/nphi)
    return N*np.exp(-Ss/hbar)*I, np.nanmax(np.abs(dS.imag[active]))
