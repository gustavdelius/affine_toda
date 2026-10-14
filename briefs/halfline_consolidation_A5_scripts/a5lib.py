"""A5 library: exact linearized reflection matrix about static Hirota backgrounds of real-coupling a_n^(1) on x <= 0.

Conventions (book ch. 31): m = 1, alpha_j = e_j - e_{j+1 mod h} (j = 0..n), phi = -(1/b) sum_j alpha_j ln tau_j,
e^{b alpha_j.phi} = tau_{j-1} tau_{j+1}/tau_j^2.  Linearized (rescaled b = 1):
   eps'' = M^2(x) eps - Omega^2 eps,   M^2 = sum_j alpha_j alpha_j^T e^{psi_j},
   BC  eps'(0) + H eps(0) = 0,         H = (1/2) sum_j C_j alpha_j alpha_j^T e^{psi_j(0)/2}.
Linear waves: add to the N static background solitons (species a_i, coefficient c_i, rapidity 0) one more Hirota
soliton of species c with E = e^{sigma x - i Omega t}, sigma = m_c cosh(theta_k), to first order in its amplitude:
   dtau_j = w^{jc} e^{sigma x} sum_{mu} prod_{i<k} A_ik^{mu_i mu_k} prod_i (A(c, a_i; sigma/m_c) c_i w^{j a_i} e^{m_{a_i} x})^{mu_i},
   A(a, b; ch) = (ch - cos(pi(a-b)/h))/(ch - cos(pi(a+b)/h))   (Hirota interaction; verified below by the PDE residual),
   eps = -sum_j alpha_j dtau_j/tau_j.
Every expression is a finite sum of exponentials, differentiated exactly; nothing is integrated.
"""
import itertools
import mpmath as mp

I = mp.mpc(0, 1)


class Theory:
    def __init__(self, n):
        self.n, self.h = n, n + 1
        self.w = mp.exp(2 * I * mp.pi / self.h)

    def mass(self, c):
        return 2 * mp.sin(mp.pi * c / self.h)

    def A(self, a, b, ch):
        h = self.h
        return (ch - mp.cos(mp.pi * (a - b) / h)) / (ch - mp.cos(mp.pi * (a + b) / h))


class ExpSum:
    """sum_r coef_r e^{rate_r x}"""
    def __init__(self, terms):
        self.t = terms

    def __call__(self, x, k=0):
        return mp.fsum(c * r**k * mp.exp(r * x) for c, r in self.t)


class Background:
    """static Hirota background: list of (species, coefficient)."""
    def __init__(self, th, sols):
        self.th, self.sols = th, sols
        T = th; h = T.h
        N = len(sols)
        self.Ast = {(i, k): T.A(sols[i][0], sols[k][0], mp.mpf(1)) for i in range(N) for k in range(N) if i < k}
        self.tau = []
        for j in range(h):
            terms = []
            for mu in itertools.product([0, 1], repeat=N):
                coef, rate = mp.mpc(1), mp.mpf(0)
                for i in range(N):
                    if mu[i]:
                        a, ci = sols[i]
                        coef *= ci * T.w**(j * a); rate += T.mass(a)
                for (i, k), Aik in self.Ast.items():
                    if mu[i] and mu[k]:
                        coef *= Aik
                terms.append((coef, rate))
            self.tau.append(ExpSum(terms))

    def wave(self, c, sigma):
        """dtau_j of the linear wave of species c with x-exponent sigma."""
        T = self.th; h = T.h; sols = self.sols; N = len(sols)
        ch = sigma / T.mass(c)
        Aw = [T.A(c, a, ch) for a, _ in sols]
        out = []
        for j in range(h):
            terms = []
            for mu in itertools.product([0, 1], repeat=N):
                coef, rate = T.w**(j * c), sigma
                for i in range(N):
                    if mu[i]:
                        a, ci = sols[i]
                        coef *= Aw[i] * ci * T.w**(j * a); rate += T.mass(a)
                for (i, k), Aik in self.Ast.items():
                    if mu[i] and mu[k]:
                        coef *= Aik
                terms.append((coef, rate))
            out.append(ExpSum(terms))
        return out

    def taus(self, x):
        return [(t(x), t(x, 1), t(x, 2)) for t in self.tau]

    def psi(self, x):
        """e^{psi_j} = tau_{j-1} tau_{j+1}/tau_j^2"""
        h = self.th.h; tv = [t(x) for t in self.tau]
        return [tv[(j - 1) % h] * tv[(j + 1) % h] / tv[j]**2 for j in range(h)]


def eps_of(bg, dtau, x):
    """eps, eps', eps'' (vectors of length h) of the linear wave dtau at x."""
    h = bg.th.h
    q = []
    for j in range(h):
        t0, t1, t2 = (bg.tau[j](x, k) for k in range(3))
        d0, d1, d2 = (dtau[j](x, k) for k in range(3))
        q0 = d0 / t0
        q1 = d1 / t0 - d0 * t1 / t0**2
        q2 = d2 / t0 - 2 * d1 * t1 / t0**2 - d0 * t2 / t0**2 + 2 * d0 * t1**2 / t0**3
        q.append((q0, q1, q2))
    # eps = -sum_j alpha_j q_j, alpha_j = e_j - e_{j+1}: eps_k = q_{k-1} - q_k
    return [[q[(k - 1) % h][r] - q[k][r] for k in range(h)] for r in range(3)]


def matvec_alpha(h, coeffs, v):
    """sum_j coeffs_j alpha_j (alpha_j . v)"""
    out = [mp.mpc(0)] * h
    for j in range(h):
        s = v[j] - v[(j + 1) % h]
        out[j] += coeffs[j] * s
        out[(j + 1) % h] -= coeffs[j] * s
    return out


def pde_residual(bg, dtau, x, Om):
    h = bg.th.h
    e0, e1, e2 = eps_of(bg, dtau, x)
    Mv = matvec_alpha(h, bg.psi(x), e0)
    r = [e2[k] - Mv[k] + Om**2 * e0[k] for k in range(h)]
    return max(abs(v) for v in r) / max(abs(v) for v in e2 + e0)


def bc_vec(bg, dtau, C):
    h = bg.th.h
    e0, e1, _ = eps_of(bg, dtau, mp.mpf(0))
    eh = [mp.sqrt(v) for v in bg.psi(mp.mpf(0))]
    Hv = matvec_alpha(h, [C[j] * eh[j] / 2 for j in range(h)], e0)
    return [e1[k] + Hv[k] for k in range(h)]


def kmatrix(bg, C, Om):
    """K[c'][c] in the Fourier polarization basis rho_c = (w^{jc})_j, for open channels (m_c < Om);
    returns (open, K, Kflux, residual). Kflux = sqrt(p_c'/p_c) K (unitary if the problem is self-adjoint)."""
    T = bg.th; h = T.h; n = T.n
    opn = [c for c in range(1, h) if T.mass(c) < Om]
    cls = [c for c in range(1, h) if T.mass(c) >= Om]
    p = {c: mp.sqrt(Om**2 - T.mass(c)**2) for c in opn}
    kap = {c: mp.sqrt(T.mass(c)**2 - Om**2) for c in cls}
    unk = [(c, -I * p[c]) for c in opn] + [(c, kap[c]) for c in cls]
    cols = [bc_vec(bg, bg.wave(c, s), C) for c, s in unk]
    M = mp.matrix(n, n)
    for i, col in enumerate(cols):
        for k in range(n):
            M[k, i] = col[k]
    # asymptotic polarization of a wave with unit Hirota amplitude: eps -> -(1 - w^{-c}) rho_c e^{sigma x}
    nf = lambda c: -(1 - T.w**(-c))
    K = {}; res = 0
    for c in opn:
        r = bc_vec(bg, bg.wave(c, I * p[c]), C)
        rhs = mp.matrix([-r[k] for k in range(n)])
        sol = mp.lu_solve(M, rhs)
        # residual on the full h rows
        full = [r[k] + mp.fsum(cols[i][k] * sol[i] for i in range(n)) for k in range(h)]
        res = max(res, max(abs(v) for v in full))
        for i, c2 in enumerate(opn):
            K[(c2, c)] = sol[i] * nf(c2) / nf(c)
    Kf = {(c2, c): K[(c2, c)] * mp.sqrt(p[c2] / p[c]) for (c2, c) in K}
    return opn, K, Kf, res


def unitarity_defect(opn, Kf):
    m = len(opn)
    U = mp.matrix(m, m)
    for i, a in enumerate(opn):
        for j, b in enumerate(opn):
            U[i, j] = Kf[(a, b)]
    return mp.mnorm(U.H * U - mp.eye(m), 1)
