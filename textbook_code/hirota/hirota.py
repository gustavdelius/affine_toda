"""Static Hirota solitons and exact fluctuation solutions for simply-laced imaginary ATFT.

Bilinear form (derived from the field equation of Section 2.2 with
phi = (i/beta) sum_j alpha_j ln tau_j, |alpha_j|^2 = 2, m = 1):

    (1/2)(D_x^2 - D_t^2) tau_j . tau_j = n_j ( tau_j^2 - prod_{l ~ j} tau_l ).

Static soliton: tau_j = sum_p t_{jp} E^p, E = exp(mu x + xi).
Fluctuation: tau_j + eps exp(i k x - i w t) S_j(E), S_j = sum_q s_{jq} E^q, w^2 = k^2 + m_b^2 = lambda.
Linearized equation (static background):
    T F'' - 2 T' F' + T'' F + w^2 T F = n_j (2 T_j F_j - sum_{l~j} F_l prod_{l'~j, l' != l} T_l').
Coefficient of exp(ikx) E^K:  sum_{p+q=K} t_p s_q [ (ik + (q-p) mu)^2 + w^2 ] = ...
"""
import numpy as np
import mpmath as mp


class Poly:
    pass


def polymul_coef(polys, K, skip=None):
    """coefficient of E^K in product of coefficient lists; polys: list of lists"""
    # straightforward convolution
    acc = [1.0 + 0j] + [0j] * K
    for P in polys:
        new = [0j] * (K + 1)
        for i in range(K + 1):
            if acc[i] == 0:
                continue
            for q in range(0, K + 1 - i):
                if q < len(P):
                    new[i + q] += acc[i] * P[q]
        acc = new
    return acc[K]


def conv(P, Q, K):
    s = 0
    for p in range(max(0, K - len(Q) + 1), min(K, len(P) - 1) + 1):
        s += P[p] * Q[K - p]
    return s


def prod_trunc(polys, K, ctx):
    """full product truncated at degree K"""
    acc = [ctx.one] + [ctx.zero] * K
    for P in polys:
        new = [ctx.zero] * (K + 1)
        for i in range(K + 1):
            if acc[i] == 0:
                continue
            for q in range(0, min(len(P), K + 1 - i)):
                new[i + q] += acc[i] * P[q]
        acc = new
    return acc


class NumCtx:
    """numpy complex128 context"""
    one = 1.0 + 0j
    zero = 0j

    @staticmethod
    def solve(A, b):
        return np.linalg.solve(np.array(A, complex), np.array(b, complex))

    @staticmethod
    def num(x):
        return complex(x)


class MpCtx:
    def __init__(self, dps=40):
        mp.mp.dps = dps
        self.one = mp.mpc(1)
        self.zero = mp.mpc(0)

    def solve(self, A, b):
        return list(mp.lu_solve(mp.matrix(A), mp.matrix(b)))

    def num(self, x):
        return mp.mpc(x)


def soliton(alg, delta, mu, Kmax=None, ctx=None, tol=1e-10, NC=None):
    """Order-by-order solution of the static Hirota equations with first-order data delta.
    Returns list of coefficient lists t[j][p] (trimmed), and the max |c_K| beyond termination."""
    ctx = ctx or NumCtx
    r1 = alg.r + 1
    NCm = alg.NC if NC is None else NC
    n = alg.n
    nb = alg.nbrs
    if Kmax is None:
        Kmax = int(4 * max(n) * 6) + 10
    t = [[ctx.one, ctx.num(delta[j])] for j in range(r1)]
    mu2 = ctx.num(mu) ** 2 if ctx is not NumCtx else mu * mu
    norms = [1.0, max(abs(complex(d)) for d in delta)]
    zero_run = 0
    for K in range(2, Kmax + 1):
        rhs = []
        for j in range(r1):
            tj = t[j] + [ctx.zero]  # c_K = 0 placeholder
            # LHS nonlinear: (mu^2/2) sum_{0<p,q<K} c_p c_q (p-q)^2
            lhs = ctx.zero
            sq = ctx.zero
            for p in range(1, K):
                q = K - p
                cp = tj[p] if p < len(tj) else ctx.zero
                cq = tj[q] if q < len(tj) else ctx.zero
                lhs += cp * cq * (p - q) ** 2
                sq += cp * cq
            lhs = lhs * mu2 / 2
            prodK = prod_trunc([t[l] + [ctx.zero] for l in nb[j]], K, ctx)[K]
            rhs.append(n[j] * (sq - prodK) - lhs)
        A = [[(mu2 * K * K if i == j else 0) - NCm[i][j] for j in range(r1)] for i in range(r1)]
        cK = ctx.solve(A, rhs)
        for j in range(r1):
            t[j].append(cK[j])
        nrm = max(abs(complex(c)) for c in cK)
        norms.append(nrm)
        if nrm < tol * max(norms):
            zero_run += 1
        else:
            zero_run = 0
        if zero_run > 12 and K > 6:
            break
    # trim
    big = max(norms)
    out = []
    for j in range(r1):
        cs = t[j]
        d = len(cs) - 1
        while d > 0 and abs(complex(cs[d])) < tol * big:
            d -= 1
        out.append(cs[:d + 1])
    tail = max([abs(complex(c)) for j in range(r1) for c in t[j][len(out[j]):]] + [0])
    return out, tail / big


def check_soliton(alg, t, mu, xs=None, xi=0.3 + 0.7j):
    """residual of Hirota equations at sample points (relative)"""
    if xs is None:
        xs = np.linspace(-6, 6, 13) / mu
    res = 0
    for x in xs:
        E = np.exp(mu * x + xi)
        tau = [sum(complex(c) * E ** p for p, c in enumerate(tj)) for tj in t]
        d1 = [sum(complex(c) * p * mu * E ** p for p, c in enumerate(tj)) for tj in t]
        d2 = [sum(complex(c) * (p * mu) ** 2 * E ** p for p, c in enumerate(tj)) for tj in t]
        for j in range(alg.r + 1):
            lhs = tau[j] * d2[j] - d1[j] ** 2
            pr = np.prod([tau[l] for l in alg.nbrs[j]])
            rhs = alg.n[j] * (tau[j] ** 2 - pr)
            scale = abs(tau[j]) ** 2 + abs(pr) + abs(d1[j]) ** 2
            res = max(res, abs(lhs - rhs) / scale)
    return res


def fluctuation(alg, t, mu, deltab, mb2, z, ctx=None, extra=6, NC=None):
    """Exact fluctuation solution about static soliton t (E = e^{mu x}) in channel with first-order data
    deltab (eigvec of NC, eigenvalue mb2), at z = i k / m_b.  Returns s[j][q] and tail norm."""
    ctx = ctx or NumCtx
    r1 = alg.r + 1
    NCm = alg.NC if NC is None else NC
    n = alg.n
    nb = alg.nbrs
    mb = np.sqrt(mb2) if ctx is NumCtx else mp.sqrt(ctx.num(mb2))
    ik = z * mb
    D = [len(tj) - 1 for tj in t]
    Kmax = max(D) + extra
    s = [[ctx.num(deltab[j])] for j in range(r1)]
    w2 = mb2 - z * z * mb2 if ctx is NumCtx else ctx.num(mb2) * (1 - ctx.num(z) ** 2)  # w^2 = k^2 + m_b^2 = m_b^2 (1 - z^2)
    tails = 0
    for K in range(1, Kmax + 1):
        rhs = []
        for j in range(r1):
            tj = t[j]
            sj = s[j] + [ctx.zero]
            acc = ctx.zero
            sm = ctx.zero
            for p in range(1, min(K, len(tj) - 1) + 1):
                q = K - p
                if q >= len(sj):
                    continue
                coef = (ik + (q - p) * mu) ** 2 + w2
                acc += tj[p] * sj[q] * coef
                sm += tj[p] * sj[q]
            # neighbour term, excluding the s_K part (which is linear, in matrix)
            nsum = ctx.zero
            for l in nb[j]:
                others = [t[l2] for l2 in nb[j] if l2 != l]
                # handle repeated neighbours (none for single bonds)
                Pl = prod_trunc(others, K, ctx) if others else [ctx.one] + [ctx.zero] * K
                sl = s[l]
                for q in range(0, min(K - 1, len(sl) - 1) + 1):
                    nsum += sl[q] * Pl[K - q]
            # equation: acc + c_K s_jK = n_j(2 sm + 2 s_jK - nsum - sum_l s_lK)
            rhs.append(n[j] * (2 * sm - nsum) - acc)
        cK = (ik + K * mu) ** 2 + w2
        A = [[(cK if i == j else 0) - NCm[i][j] for j in range(r1)] for i in range(r1)]
        sK = ctx.solve(A, rhs)
        for j in range(r1):
            s[j].append(sK[j])
    for j in range(r1):
        for q in range(D[j] + 1, len(s[j])):
            tails = max(tails, abs(complex(s[j][q])))
    big = max(abs(complex(c)) for sj in s for c in sj)
    return s, tails / big


def plus_infinity_vector(alg, t, s):
    """f_j(+inf) = s_{j,D_j}/t_{j,D_j}; returns delta-phi direction sum_j alpha_j f_j(+inf)"""
    f = np.array([complex(s[j][len(t[j]) - 1]) / complex(t[j][-1]) for j in range(alg.r + 1)])
    return f, alg.alpha.T @ f


# ---------------------------------------------------------------------------------------------
# Robust versions: compute up to a given degree, then verify the polynomial identities exactly
# (all orders) in high precision.

def polymul(P, Q):
    out = [0] * (len(P) + len(Q) - 1)
    for i, p in enumerate(P):
        if p == 0:
            continue
        for j, q in enumerate(Q):
            out[i + j] += p * q
    return out


def polyadd(P, Q):
    n = max(len(P), len(Q))
    return [(P[i] if i < len(P) else 0) + (Q[i] if i < len(Q) else 0) for i in range(n)]


def soliton_deg(alg, delta, mu, degs=None, ctx=None, NC=None):
    """order-by-order up to max(degs); coefficients beyond degs[j] dropped; returns t and the
    max relative size of dropped coefficients and of the residual of the full identities."""
    ctx = ctx or MpCtx(40)
    r1 = alg.r + 1
    if degs is None:
        degs = [int(round(x)) for x in alg.n]
    Kmax = max(degs)
    NCm = alg.NC if NC is None else NC
    NCm = [[ctx.num(NCm[i][j]) for j in range(r1)] for i in range(r1)]
    n = [ctx.num(x) for x in alg.n]
    nb = alg.nbrs
    mu = ctx.num(mu)
    mu2 = mu * mu
    t = [[ctx.one, ctx.num(delta[j])] for j in range(r1)]
    resonances = []
    for K in range(2, Kmax + 1):
        rhs = []
        for j in range(r1):
            tj = t[j] + [ctx.zero]
            lhs = ctx.zero
            sq = ctx.zero
            for p in range(1, K):
                q = K - p
                lhs += tj[p] * tj[q] * (p - q) ** 2
                sq += tj[p] * tj[q]
            lhs = lhs * mu2 / 2
            prodK = prod_trunc([t[l] + [ctx.zero] for l in nb[j]], K, ctx)[K]
            rhs.append(n[j] * (sq - prodK) - lhs)
        cK = resonant_solve(alg, mu2 * K * K, rhs, ctx, info=resonances)
        for j in range(r1):
            t[j].append(cK[j])
    big = max(abs(c) for tj in t for c in tj)
    dropped = 0
    out = []
    for j in range(r1):
        dropped = max([dropped] + [abs(c) / big for c in t[j][degs[j] + 1:]])
        out.append(t[j][:degs[j] + 1])
    return out, dropped, soliton_residual(alg, out, mu, ctx)


def soliton_residual(alg, t, mu, ctx):
    """max |coefficient| of  tau tau'' - tau'^2 - n_j (tau^2 - prod tau_l), relative"""
    mu = ctx.num(mu)
    res = 0
    big = 0
    for j in range(alg.r + 1):
        tj = t[j]
        d1 = [c * p * mu for p, c in enumerate(tj)]
        d2 = [c * (p * mu) ** 2 for p, c in enumerate(tj)]
        lhs = polyadd(polymul(tj, d2), [-c for c in polymul(d1, d1)])
        pr = [ctx.one]
        for l in alg.nbrs[j]:
            pr = polymul(pr, t[l])
        rhs = [ctx.num(alg.n[j]) * c for c in polyadd(polymul(tj, tj), [-c for c in pr])]
        diff = polyadd(lhs, [-c for c in rhs])
        res = max([res] + [abs(c) for c in diff])
        big = max([big] + [abs(c) for c in lhs] + [abs(c) for c in rhs])
    return res / big


def fluct_deg(alg, t, mu, deltab, mb2, z, ctx=None, NC=None):
    """fluctuation polynomial S_j with deg S_j <= deg T_j, by recursion; returns s, dropped, residual"""
    ctx = ctx or MpCtx(40)
    r1 = alg.r + 1
    NCm = alg.NC if NC is None else NC
    NCm = [[ctx.num(NCm[i][j]) for j in range(r1)] for i in range(r1)]
    n = [ctx.num(x) for x in alg.n]
    nb = alg.nbrs
    D = [len(tj) - 1 for tj in t]
    Kmax = max(D)
    mu = ctx.num(mu)
    mb2 = ctx.num(mb2)
    z = ctx.num(z)
    mb = mp.sqrt(mb2)
    ik = z * mb
    w2 = mb2 * (1 - z * z)
    s = [[ctx.num(deltab[j])] for j in range(r1)]
    for K in range(1, Kmax + 1):
        rhs = []
        for j in range(r1):
            tj = t[j]
            sj = s[j]
            acc = ctx.zero
            sm = ctx.zero
            for p in range(1, min(K, len(tj) - 1) + 1):
                q = K - p
                if q >= len(sj):
                    continue
                acc += tj[p] * sj[q] * ((ik + (q - p) * mu) ** 2 + w2)
                sm += tj[p] * sj[q]
            nsum = ctx.zero
            for l in nb[j]:
                others = [t[l2] for l2 in nb[j] if l2 != l]
                Pl = prod_trunc(others, K, ctx) if others else [ctx.one] + [ctx.zero] * K
                for q in range(0, min(K - 1, len(s[l]) - 1) + 1):
                    nsum += s[l][q] * Pl[K - q]
            rhs.append(n[j] * (2 * sm - nsum) - acc)
        cK = (ik + K * mu) ** 2 + w2
        A = [[(cK if i == j else 0) - NCm[i][j] for j in range(r1)] for i in range(r1)]
        sK = ctx.solve(A, rhs)
        for j in range(r1):
            s[j].append(sK[j])
    big = max(abs(c) for sj in s for c in sj)
    dropped = 0
    out = []
    for j in range(r1):
        dropped = max([dropped] + [abs(c) / big for c in s[j][D[j] + 1:]])
        out.append(s[j][:D[j] + 1])
    return out, dropped, fluct_residual(alg, t, out, mu, mb2, z, ctx)


def fluct_residual(alg, t, s, mu, mb2, z, ctx):
    """residual of the linearized bilinear equation as polynomial identity in E (factor e^{ikx} removed)"""
    mb = mp.sqrt(mb2)
    ik = z * mb
    w2 = mb2 * (1 - z * z)
    res = 0
    big = 0
    for j in range(alg.r + 1):
        tj, sj = t[j], s[j]
        L = [ctx.zero] * (len(tj) + len(sj) - 1)
        for p, tp in enumerate(tj):
            for q, sq in enumerate(sj):
                L[p + q] += tp * sq * ((ik + (q - p) * mu) ** 2 + w2)
        R = [2 * c for c in polymul(tj, sj)]
        for l in alg.nbrs[j]:
            pr = list(s[l])
            for l2 in alg.nbrs[j]:
                if l2 != l:
                    pr = polymul(pr, t[l2])
            R = polyadd(R, [-c for c in pr])
        R = [ctx.num(alg.n[j]) * c for c in R]
        diff = polyadd(L, [-c for c in R])
        res = max([res] + [abs(c) for c in diff])
        big = max([big] + [abs(c) for c in L] + [abs(c) for c in R])
    return res / big


def resonant_solve(alg, shift, rhs, ctx, info=None, tol=1e-12):
    """solve (shift - NC) c = rhs.  If shift is (numerically) an eigenvalue m_c^2 of NC (a resonance
    K m_a = m_c), return the solution with no component along the resonant eigenvector (N^{-1}-orthogonal
    complement) and record the size of the resonant component of rhs (must vanish for consistency)."""
    r1 = alg.r + 1
    sq = [mp.sqrt(mp.mpf(int(round(x)))) for x in alg.n]
    Cs = mp.matrix([[sq[i] * int(round(alg.C[i, j])) * sq[j] for j in range(r1)] for i in range(r1)])
    ev, V = mp.eigsy(Cs)
    y = [rhs[i] / sq[i] for i in range(r1)]
    out = [mp.mpc(0)] * r1
    for k in range(r1):
        proj = sum(V[i, k] * y[i] for i in range(r1))
        den = shift - ev[k]
        if abs(den) < tol * (1 + abs(shift)):
            if info is not None:
                info.append((complex(shift), float(abs(proj))))
            continue
        for i in range(r1):
            out[i] += V[i, k] * proj / den
    return [out[i] * sq[i] for i in range(r1)]
