"""Task 3: B_1 -> Bbar_1 amplitude of the SP matrices.
(a) double precision, book code: zero-weight block of Rc_VVbar against the closed form obtained from crossing and the
    V V entries r_d; singlet projection with rcheck and the singlets of an_vector (as in an_nondiagonal_scalar check 2).
(b) 30 digits: the same projection with mpmath matrices, times a(mu_a) F(2mu) bt(mu_b).
Usage: python3 -u t3_breather.py"""
import numpy as np
import mpmath as mp
from sp_lib import an_vector, q_phys, cop, nodes_of, K_V, K_Vbar, etabar, Rc, kron
from sp_mp import SP, blk

mp.mp.dps = 30


# ------------------------------------------------------------ closed form of Rc_VV diagonal entries and Rc_VVbar(0-weight)
def rvv_diag(n, q, z, a, b):
    """Rc_VV(z)[aN+b, aN+b]: 1 for a = b, r_d = (q - 1/q) z^(N-d)/(q Z - 1/q) for d = (b - a) mod N, Z = z^N."""
    N = n + 1
    if a == b:
        return 1
    d = (b - a) % N
    return (q - 1/q)*z**(N - d)/(q*z**N - 1/q)


def t_of(q, Z):
    return (Z - 1)/(q*Z - 1/q)


def rvw0(n, q, om, th, a, b, lib=np):
    """Zero-weight block of Rc_VVbar(th): [aN+a, bN+b] = Rc_VV(i pi - th)[aN+b, aN+b] / t(i pi - th)."""
    zc = lib.exp(om*(1j*lib.pi - th))
    return rvv_diag(n, q, zc, a, b)/t_of(q, zc**(n + 1))


def singlet_np(A, B):
    d = A['dim']*B['dim']
    M = np.vstack([cop(A, B, kd, j) - (np.eye(d) if kd == 'k' else 0) for j in nodes_of(A) for kd in 'efk'])
    _, sv, Vh = np.linalg.svd(M)
    return Vh[-1].conj(), sv[-1]/sv[0]


def proj_np(n, om, l, eta, mu, bar=False):
    """Matrix coefficient c (and residual): B_1 (V(ta) x Vbar(tb)) -> Bbar_1 (Vbar(-tb) x V(-ta)), or the conjugate.
    Normalization: in-singlet and out-singlet have coefficient 1 on (hw of left factor) (x) (lw of right factor)."""
    N, T, q = n + 1, (n + 1)*om/2, q_phys(om)
    th = 1j*np.pi*mu/T
    u = np.pi*(1 - 1/T)
    ta, tb = th - 0.5j*u, th + 0.5j*u
    eb = etabar(n, om, l, eta)
    KV = lambda t: K_V(n, om, l, eta, t)
    KW = lambda t: K_Vbar(n, om, l, eb, t)
    A_in, B_in = an_vector(n, q, np.exp(om*ta), conj=bar), an_vector(n, q, np.exp(om*tb), conj=not bar)
    A_out, B_out = an_vector(n, q, np.exp(-om*tb), conj=not bar), an_vector(n, q, np.exp(-om*ta), conj=bar)
    w_in, r1 = singlet_np(A_in, B_in)
    w_out, r2 = singlet_np(A_out, B_out)
    i_in = A_in['hw']*N + (N - 1 - B_in['hw'])          # hw (x) lw  (lw index = n - hw index)
    i_out = A_out['hw']*N + (N - 1 - B_out['hw'])
    w_in, w_out = w_in/w_in[i_in], w_out/w_out[i_out]
    Kb, Ka = (KW(tb), KV(ta)) if not bar else (KV(tb), KW(ta))
    R = Rc(n, om, 'V' if not bar else 'W', 'W' if not bar else 'V', np.exp(om*(ta + tb)))
    v = kron(np.eye(N), Ka) @ (R @ (kron(np.eye(N), Kb) @ w_in))
    c = np.vdot(w_out, v)/np.vdot(w_out, w_out)
    return c, np.linalg.norm(v - c*w_out)/np.linalg.norm(v), max(r1, r2), w_in, w_out


# ------------------------------------------------------------ mpmath versions
def an_vector_mp(n, q, y, conj=False):
    N = n + 1
    g = {}
    for j in range(N):
        a, b = (j - 1) % N, j % N
        if j == 0:
            a, b = n, 0
        E = mp.zeros(N, N); E[a, b] = 1
        F = E.T
        H = [0]*N; H[a] += 1; H[b] -= 1
        K = mp.diag([q**(-h if conj else h) for h in H])
        g[('e', j)] = y*(F if conj else E)
        g[('f', j)] = (E if conj else F)/y
        g[('k', j)] = K
    g['dim'] = N
    g['hw'] = n if conj else 0
    return g


def kron_mp(A, B):
    ra, ca, rb, cb = A.rows, A.cols, B.rows, B.cols
    C = mp.zeros(ra*rb, ca*cb)
    for i in range(ra):
        for j in range(ca):
            if A[i, j] != 0:
                for k in range(rb):
                    for m in range(cb):
                        C[i*rb + k, j*cb + m] = A[i, j]*B[k, m]
    return C


def singlet_mp(A, B, i_norm):
    """Zero-weight singlet sum_i x_i a_i (x) b_i' in A (x) B, normalized x[i_norm] = 1 (30 digits)."""
    N = A['dim']
    I = mp.eye(N)
    rows = []
    for j in range(N):
        e = kron_mp(A[('e', j)], I) + kron_mp(A[('k', j)], B[('e', j)])
        f = kron_mp(A[('f', j)], B[('k', j)]**-1) + kron_mp(I, B[('f', j)])
        rows += [e, f]
    # zero-weight subspace: basis vectors with the same index in both factors (V x Vbar and Vbar x V)
    cols = [i*N + i for i in range(N)]
    M = mp.matrix(sum(r.rows for r in rows), N)
    r0 = 0
    for r in rows:
        for i in range(r.rows):
            for c_, col in enumerate(cols):
                M[r0 + i, c_] = r[i, col]
        r0 += r.rows
    # fix x[i_norm] = 1 and solve least squares
    others = [c_ for c_ in range(N) if c_ != i_norm]
    Ared = mp.matrix(M.rows, N - 1)
    rhs = mp.matrix(M.rows, 1)
    for i in range(M.rows):
        rhs[i] = -M[i, i_norm]
        for k, c_ in enumerate(others):
            Ared[i, k] = M[i, c_]
    sol = mp.lu_solve(Ared.H*Ared, Ared.H*rhs)
    x = [mp.mpf(0)]*N
    x[i_norm] = mp.mpf(1)
    for k, c_ in enumerate(others):
        x[c_] = sol[k]
    res = mp.norm(M*mp.matrix(x))
    return x, res


def proj_mp(S, mu, bar=False):
    n, om, N, T, q = S.n, S.om, S.N, S.T, S.q
    th = 1j*mp.pi*mu/T
    u = mp.pi*(1 - 1/T)
    ta, tb = th - 0.5j*u, th + 0.5j*u
    kV = lambda t, j: (mp.exp(2*j*om*t)*mp.sinh(S.eta + T*t) if j < S.l
                       else mp.exp((2*j - N)*om*t)*mp.sinh(S.eta - T*t))
    kW = lambda t, j: (mp.exp(-2*j*om*t)*mp.sinh(S.etab - T*t) if j < S.l
                       else mp.exp(-(2*j - N)*om*t)*mp.sinh(S.etab + T*t))
    A_in, B_in = an_vector_mp(n, q, mp.exp(om*ta), bar), an_vector_mp(n, q, mp.exp(om*tb), not bar)
    A_out, B_out = an_vector_mp(n, q, mp.exp(-om*tb), not bar), an_vector_mp(n, q, mp.exp(-om*ta), bar)
    x, r1 = singlet_mp(A_in, B_in, A_in['hw'])
    y, r2 = singlet_mp(A_out, B_out, A_out['hw'])
    kb = [(kW if not bar else kV)(tb, j) for j in range(N)]
    ka = [(kV if not bar else kW)(ta, j) for j in range(N)]
    vals = []
    for a in range(N):
        if not bar:       # Rc_VVbar(ta + tb) zero-weight block
            s = sum(rvw0(n, q, om, ta + tb, a, b, mp)*kb[b]*x[b] for b in range(N))
        else:             # Rc_VbarV: by form 1 of the crossing, [aN+a, bN+b] = Rc_VV(i pi - th)[bN+a, bN+a]/t
            zc = mp.exp(om*(1j*mp.pi - (ta + tb)))
            s = sum(rvv_diag(n, q, zc, b, a)/t_of(q, zc**N)*kb[b]*x[b] for b in range(N))
        vals.append(ka[a]*s/y[a])
    spread = max(abs(v/vals[0] - 1) for v in vals)
    return vals[0], spread, max(r1, r2)


if __name__ == '__main__':
    # (a) closed form of the zero-weight blocks against rcheck
    dmax = 0
    for n in (1, 2, 3, 4):
        N = n + 1
        for om in (0.37, 1.29):
            q = q_phys(om)
            for th in (0.31 + 0.12j, -0.45 + 0.27j):
                RVW = Rc(n, om, 'V', 'W', np.exp(om*th))
                RWV = Rc(n, om, 'W', 'V', np.exp(om*th))
                zc = np.exp(om*(1j*np.pi - th))
                for a in range(N):
                    for b in range(N):
                        dmax = max(dmax, abs(RVW[a*N + a, b*N + b] - rvw0(n, q, om, th, a, b)),
                                   abs(RWV[a*N + a, b*N + b] - rvv_diag(n, q, zc, b, a)/t_of(q, zc**N)))
    print(f'(a) zero-weight blocks of Rc_VVbar, Rc_VbarV: closed form vs rcheck, n = 1..4: {dmax:.1e}', flush=True)
    # (b) projection: double precision (book singlets and rcheck) vs 30-digit mpmath
    for n in (1, 2, 3):
        for l in range(n + 2):
            for eta in (0.41 + 0.23j, -0.3 + 1.1j):
                S = SP(n, mp.mpf('0.37'), l, mp.mpc(eta.real, eta.imag))
                for mu in (0.11 + 0.19j, -0.23 + 0.07j):
                    for bar in (False, True):
                        c, res, sres, _, _ = proj_np(n, 0.37, l, eta, mu, bar)
                        cm, spread, r = proj_mp(S, mp.mpc(mu.real, mu.imag), bar)
                        print(f'n={n} l={l} eta={eta} mu={mu} {"Bbar->B" if bar else "B->Bbar"}: np res {res:.1e} '
                              f'c_np/c_mp - 1 = {abs(c/complex(cm) - 1):.1e}  mp spread {mp.nstr(spread, 2)} '
                              f'singlet res {mp.nstr(r, 2)}', flush=True)
