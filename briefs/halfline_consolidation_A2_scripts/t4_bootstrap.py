"""Task 4: SP soliton bootstrap V^(x n) -> Vbar with book code (rcheck, an_vector, intertwiners).

In-state  V(th'_1) (x) ... (x) V(th'_n),  mu'_k = mu - T + k omega  (the only ordering containing Vbar(th)).
Boundary on the right; particles reflect from the right, each incoming one passes the already reflected ones
(R-check_VV at th_k + th_j) and reflects with K_V.  Out-state V(-th'_n) (x) ... (x) V(-th'_1): same rapidity pattern,
so the same embedding matrix iota (checked).  Matrix part:  M iota = P(mu) iota K_Vbar(th).
Physical:  bt(mu) = P(mu) prod_k a(mu'_k) prod_{j<k} F11(mu'_j + mu'_k) * [CDD]   ->  Q = bt / (P prod a prod F11).
Usage: python3 -u t4_bootstrap.py n [omega]"""
import sys
import numpy as np
import mpmath as mp
from sp_lib import an_vector, q_phys, tensor_rep, nodes_of, K_V, K_Vbar, etabar, Rc
from boundary_common import intertwiners
from sp_mp import SP, F11, blk


def embedding(n, om, th):
    N, T, q = n + 1, (n + 1)*om/2, q_phys(om)
    mu = -1j*T*th/np.pi
    mus = [mu - T + k*om for k in range(1, n + 1)]
    W = tensor_rep([an_vector(n, q, np.exp(om*1j*np.pi*m/T)) for m in mus])
    Ub = an_vector(n, q, np.exp(om*th), conj=True)
    null, rel = intertwiners([(Ub[(kd, j)], W[(kd, j)]) for j in nodes_of(Ub) for kd in 'efk'], tol=1e-9)
    assert len(null) == 1, rel[:3]
    return null[0], mus


def apply_two(v, M, p, n, N):
    m = v.shape[1]
    t = v.reshape((N**p, N*N, N**(n - p - 2), m))
    return np.einsum('ab,ibjm->iajm', M, t).reshape(N**n, m)


def apply_one(v, K, p, n, N):
    m = v.shape[1]
    t = v.reshape((N**p, N, N**(n - p - 1), m))
    return np.einsum('ab,ibjm->iajm', K, t).reshape(N**n, m)


def boot_matrix(n, om, l, eta, th):
    N = n + 1
    T = N*om/2
    io_in, mus = embedding(n, om, th)
    io_out, _ = embedding(n, om, -th)
    i0 = np.unravel_index(np.argmax(np.abs(io_in)), io_in.shape)
    io_in, io_out = io_in/io_in[i0], io_out/io_out[i0]
    same = np.abs(io_in - io_out).max()
    ths = [1j*np.pi*m/T for m in mus]
    v = io_in.copy()
    pos = list(range(n))                      # particle index at each position; reflected ones marked by sign
    refl = [False]*n
    for k in range(n - 1, -1, -1):
        p = k
        while p < n - 1:
            j = pos[p + 1]
            v = apply_two(v, Rc(n, om, 'V', 'V', np.exp(om*(ths[k] + ths[j]))), p, n, N)
            pos[p], pos[p + 1] = pos[p + 1], pos[p]
            p += 1
        v = apply_one(v, K_V(n, om, l, eta, ths[k]), n - 1, n, N)
        refl[k] = True
    target = io_out @ K_Vbar(n, om, l, etabar(n, om, l, eta), th)
    P = np.vdot(target.ravel(), v.ravel())/np.vdot(target.ravel(), target.ravel())
    res = np.linalg.norm(v - P*target)/np.linalg.norm(v)
    return P, res, same, mus


if __name__ == '__main__':
    n = int(sys.argv[1])
    om = float(sys.argv[2]) if len(sys.argv) > 2 else 0.37
    mp.mp.dps = 30
    N, T = n + 1, (n + 1)*om/2
    for l in range(N + 1):
        for eta in (0.41 + 0.23j, -0.3 + 1.1j):
            S = SP(n, mp.mpf(om), l, mp.mpc(eta.real, eta.imag))
            for mu in (0.11 + 0.19j, -0.23 + 0.07j, 0.31 - 0.12j):
                th = 1j*np.pi*mu/T
                P, res, same, mus = boot_matrix(n, om, l, eta, th)
                mmu = mp.mpc(mu.real, mu.imag)
                mmus = [mmu - S.T + k*S.om for k in range(1, n + 1)]
                rhs = mp.mpc(P.real, P.imag)
                for k in range(n):
                    rhs *= S.a(mmus[k])
                for j in range(n):
                    for k in range(j + 1, n):
                        rhs *= F11(mmus[j] + mmus[k], n, S.om)
                Q = S.bt(mmu)/rhs
                # exponential family m = +-1 for comparison
                Qm = [S.bt(mmu, m)/(rhs/mp.fprod([S.a(x) for x in mmus])*mp.fprod([S.a(x, m) for x in mmus]))
                      for m in (1, -1)]
                print(f'n={n} om={om} l={l} eta={eta} mu={mu}: matrix res {res:.1e} (iota_out=iota_in {same:.1e}) '
                      f'P={complex(P):.6f}  Q={mp.nstr(Q, 12)}  Q(m=+1)={mp.nstr(Qm[0], 8)}  Q(m=-1)={mp.nstr(Qm[1], 8)}',
                      flush=True)
