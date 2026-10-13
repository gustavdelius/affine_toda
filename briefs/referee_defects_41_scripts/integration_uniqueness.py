"""Integration of the ch. 41 referee report, item 1: uniqueness of the a_r Weyl-module intertwiners.

Written from the chapter's definitions (eq-dt-weyl, an_vector conventions), not from the referee's scripts.

(A) Exact certificate (rational arithmetic, fractions.Fraction): the space Hom_{U_q(b)}(P (x) W, W (x) P') restricted to
    one charge sector Q is cut out by the loop conditions L_B T_Q = T_Q L_A (commutator loops of two nodes and the loop
    through all nodes once).  Every intertwiner satisfies them and is determined by T_Q (the sector maps of Delta(f_j),
    Delta(e_j) are invertible), so nullity 1 of the loop system bounds dim Hom by 1.  At a rational point; rank is lower
    semicontinuous, so the bound holds for generic parameters.
(B) The all-rank argument: Omega = Delta(f_j)Delta(f_{j+1}) - q^-1 Delta(f_{j+1})Delta(f_j) has no pure-defect part;
    on V (x) W its image in each sector is the coordinate line of v_{j+2}; on Vbar (x) W its kernel is the coordinate
    hyperplane of the vbar_i, i != j+2.  Same for b^+ with e (Omega+ = Delta(e_j)Delta(e_{j+1}) - q Delta(e_{j+1})Delta(e_j)).
    Checked numerically at random complex parameters for h = 3..6.
(C) Determinants of the closed forms on a sector (sympy, symbolic q, X, offset): independent of the sector.
(D) Root of unity: at q^5 = 1, S|N> = |N + 5 alpha_1> commutes with the module.
(E) Reversed orientation: the two-entry Vbar intertwiner of the module with q^{-(N+s)_{j+1}}, and its coefficient qX in the
    basis q^{-|n|^2/2}|n>.
"""
from fractions import Fraction as Fr
import itertools
import numpy as np
import sympy as sp
import time

T0 = time.time()


# ------------------------------------------------------------------ generic exact linear algebra
def rank_exact(rows, ncols):
    """Rank of a list of rows (lists of Fractions) by Gaussian elimination."""
    M = [list(r) for r in rows if any(x != 0 for x in r)]
    rank = 0; col = 0; nrows = len(M)
    while rank < nrows and col < ncols:
        piv = next((i for i in range(rank, nrows) if M[i][col] != 0), None)
        if piv is None:
            col += 1; continue
        M[rank], M[piv] = M[piv], M[rank]
        p = M[rank][col]
        M[rank] = [x/p for x in M[rank]]
        for i in range(nrows):
            if i != rank and M[i][col] != 0:
                c = M[i][col]
                M[i] = [a - c*b for a, b in zip(M[i], M[rank])]
        rank += 1; col += 1
    return rank


def mat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum((A[i][k]*B[k][j] for k in range(m) if A[i][k] != 0 and B[k][j] != 0), Fr(0)) for j in range(p)]
            for i in range(n)]


def mat_inv(A):
    n = len(A)
    M = [list(A[i]) + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        p = M[c][c]; M[c] = [x/p for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f*b for a, b in zip(M[i], M[c])]
    return [row[n:] for row in M]


# ------------------------------------------------------------------ representations (exact, q rational, offsets t_i = q^{s_i})
def qpow(q, n):
    return q**n if n >= 0 else 1/(q**(-n))


def alpha(h, j):
    """alpha_j = e_j - e_{j+1} (1-based, j = 1..r), alpha_0 = e_h - e_1; as integer vectors (0-based components)."""
    a = [0]*h
    if j == 0:
        a[h - 1] += 1; a[0] -= 1
    else:
        a[j - 1] += 1; a[j] -= 1
    return tuple(a)


def particle(h, q, y, conj=False):
    """V(y) (conj False) or Vbar(y): basis 0..h-1 (v_1..v_h), weights +-e_i.
    V: f_j v_j = v_{j+1}/y, e_j v_{j+1} = y v_j, k_j = q^{alpha_j . e_i}.  Vbar: f_j vb_{j+1} = vb_j/y, e_j vb_j = y vb_{j+1},
    k_j = q^{-alpha_j . e_i}.  Node j moves v_j <-> v_{j+1} (1-based; node 0: v_h <-> v_1)."""
    rep = {'dim': h, 'wt': [tuple((1 if c == i else 0)*(-1 if conj else 1) for c in range(h)) for i in range(h)]}
    for j in range(h):
        a, b = (j - 1) % h, j % h            # 0-based v_j, v_{j+1}  (j = 0: v_h, v_1)
        if j == 0:
            a, b = h - 1, 0
        F = [[Fr(0)]*h for _ in range(h)]; E = [[Fr(0)]*h for _ in range(h)]
        if not conj:
            F[b][a] = 1/y; E[a][b] = y
        else:
            F[a][b] = 1/y; E[b][a] = y
        al = alpha(h, j)
        K = [[Fr(0)]*h for _ in range(h)]
        for i in range(h):
            K[i][i] = qpow(q, sum(al[c]*rep['wt'][i][c] for c in range(h)))
        rep[('f', j)] = F; rep[('e', j)] = E; rep[('k', j)] = K
    return rep


def tensor(A, B):
    """Coproduct Delta(e) = e x 1 + k x e, Delta(f) = f x k^-1 + 1 x f, Delta(k) = k x k on A (x) B; basis (a, b) -> a*dB + b."""
    dA, dB = A['dim'], B['dim']
    def kron(X, Y):
        return [[X[i1][j1]*Y[i2][j2] for j1 in range(dA) for j2 in range(dB)] for i1 in range(dA) for i2 in range(dB)]
    IA = [[Fr(int(i == j)) for j in range(dA)] for i in range(dA)]
    IB = [[Fr(int(i == j)) for j in range(dB)] for i in range(dB)]
    rep = {'dim': dA*dB, 'wt': [tuple(x + y for x, y in zip(wa, wb)) for wa in A['wt'] for wb in B['wt']]}
    h = len(A['wt'][0])
    for j in range(h):
        kinvB = mat_inv(B[('k', j)])
        e = kron(A[('e', j)], IB); e2 = kron(A[('k', j)], B[('e', j)])
        f = kron(A[('f', j)], kinvB); f2 = kron(IA, B[('f', j)])
        rep[('e', j)] = [[a + b for a, b in zip(r1, r2)] for r1, r2 in zip(e, e2)]
        rep[('f', j)] = [[a + b for a, b in zip(r1, r2)] for r1, r2 in zip(f, f2)]
        rep[('k', j)] = kron(A[('k', j)], B[('k', j)])
    return rep


class Weyl:
    """eq-dt-weyl on labels n = N + s, N in Z^h with sum 0; q^{n_i} = t_i q^{N_i} with t_i = q^{s_i} (prod t_i = 1).
    b^- (sign -1): f_j|N> = m q^{orient n_{j+1}} |N - alpha_j>;  b^+ (sign +1): e_j|N> = m q^{orient n_{j+1}} |N + alpha_j>;
    k_j|N> = q^{alpha_j . n}.  (n_{j+1} is the component on the second end of alpha_j.)"""
    def __init__(self, h, q, m, t, sign, orient=1):
        self.h, self.q, self.m, self.t, self.sign, self.orient = h, q, m, t, sign, orient

    def qn(self, N, c):
        return self.t[c]*qpow(self.q, N[c])

    def k(self, j, N):
        al = alpha(self.h, j); out = Fr(1)
        for c in range(self.h):
            if al[c]:
                out *= self.qn(N, c)**al[c]
        return out

    def x(self, j, N):
        """(coefficient, new label) of the Borel generator x_j = f_j (b^-) or e_j (b^+)."""
        al = alpha(self.h, j)
        second = [c for c in range(self.h) if al[c] == -1][0]
        coef = self.m*self.qn(N, second)**self.orient
        return coef, tuple(n + self.sign*a for n, a in zip(N, al))


def sector_basis(P, Q, front):
    """Basis of the sector of total charge Q (integer vector, sum may be nonzero; labels N = Q - wt): list of (p, N)."""
    return [(p, tuple(qc - w for qc, w in zip(Q, P['wt'][p]))) for p in range(P['dim'])]


def sector_map(P, W, j, Q, particle_first):
    """Matrix of Delta(x_j) (x = f for b^-, e for b^+) from sector Q to sector Q +- alpha_j, on P (x) W (particle_first) or
    W (x) P, in the sector bases (index p)."""
    h = W.h; sgn = W.sign; al = alpha(h, j)
    Qn = tuple(qc + sgn*a for qc, a in zip(Q, al))
    src = sector_basis(P, Q, particle_first); tgt = sector_basis(P, Qn, particle_first)
    look = {st: i for i, st in enumerate(tgt)}
    d = P['dim']; M = [[Fr(0)]*d for _ in range(d)]
    X = 'f' if sgn == -1 else 'e'
    kinvP = mat_inv(P[('k', j)])
    for i, (p, N) in enumerate(src):
        cW, N2 = W.x(j, N)
        kW = W.k(j, N)
        if X == 'f':
            # Delta(f) = f x k^-1 + 1 x f
            if particle_first:     # P (x) W
                for p2 in range(d):
                    if P[('f', j)][p2][p] != 0:
                        M[look[(p2, N)]][i] += P[('f', j)][p2][p]/kW
                M[look[(p, N2)]][i] += cW
            else:                  # W (x) P
                M[look[(p, N2)]][i] += cW*kinvP[p][p]
                for p2 in range(d):
                    if P[('f', j)][p2][p] != 0:
                        M[look[(p2, N)]][i] += P[('f', j)][p2][p]
        else:
            # Delta(e) = e x 1 + k x e
            if particle_first:
                for p2 in range(d):
                    if P[('e', j)][p2][p] != 0:
                        M[look[(p2, N)]][i] += P[('e', j)][p2][p]
                M[look[(p, N2)]][i] += P[('k', j)][p][p]*cW
            else:
                M[look[(p, N2)]][i] += cW
                for p2 in range(d):
                    if P[('e', j)][p2][p] != 0:
                        M[look[(p2, N)]][i] += kW*P[('e', j)][p2][p]
    return M, Qn


def path(P, W, Q, nodes, particle_first):
    """Product of the sector maps along the nodes (first node first)."""
    d = P['dim']; M = [[Fr(int(i == j)) for j in range(d)] for i in range(d)]; Qc = Q
    for j in nodes:
        A, Qc = sector_map(P, W, j, Qc, particle_first)
        M = mat_mul(A, M)
    return M, Qc


def loops(P, W, Q, particle_first):
    h = W.h; out = []
    for j, k in itertools.combinations(range(h), 2):
        Pjk, Q1 = path(P, W, Q, [j, k], particle_first)
        Pkj, Q2 = path(P, W, Q, [k, j], particle_first)
        assert Q1 == Q2
        out.append(mat_mul(mat_inv(Pkj), Pjk))
    C, Qc = path(P, W, Q, list(range(h)), particle_first)
    assert Qc == Q
    out.append(C)
    return out


def loop_nullity(Psrc, Ptgt, W, Q):
    """Nullity of { T : L_B T = T L_A for all loops } with A on Psrc (x) W, B on W (x) Ptgt (sector bases)."""
    LA = loops(Psrc, W, Q, True); LB = loops(Ptgt, W, Q, False)
    d = Psrc['dim']
    rows = []
    for A, B in zip(LA, LB):
        # (B T - T A)_{ab} = sum_c B_ac T_cb - T_ac A_cb ; unknown T_cb at index c*d + b
        for a in range(d):
            for b in range(d):
                r = [Fr(0)]*(d*d)
                for c in range(d):
                    r[c*d + b] += B[a][c]
                    r[a*d + c] -= A[c][b]
                rows.append(r)
    rk = rank_exact(rows, d*d)
    return d*d - rk, rk, d*d


q = Fr(2); m = Fr(3, 7); y = Fr(5, 4); y2 = Fr(11, 3)
offsets = {3: [(Fr(1), Fr(1), Fr(1)), (Fr(3, 2), Fr(5, 7), Fr(14, 15))],
           4: [(Fr(1),)*4, (Fr(3, 2), Fr(5, 7), Fr(2, 3), Fr(7, 5))]}
print('(A) exact loop certificates at q = 2, m = mu/y_d or mu y_d = 3/7, y = 5/4 (offsets t_i = q^{s_i}):')
for h in (3, 4):
    Q = tuple([1] + [0]*(h - 1))
    for t in offsets[h]:
        for sign in (-1, 1):
            W = Weyl(h, q, m, t, sign)
            for conj in (False, True):
                P = particle(h, q, y, conj)
                nul, rk, n = loop_nullity(P, P, W, Q)
                print(f'  h={h} t={"generic" if t[0] != 1 else "s=0"} {"b+" if sign == 1 else "b-"} {"Vbar" if conj else "V   "}: '
                      f'unknowns {n}, rank {rk}, nullity {nul}')
# hypothesis (ii) for a_2: V(y1) (x) V(y2) (x) W -> W (x) V(y2) (x) V(y1)
h = 3
for t in offsets[3]:
    W = Weyl(h, q, m, t, -1)
    P12 = tensor(particle(h, q, y), particle(h, q, y2)); P21 = tensor(particle(h, q, y2), particle(h, q, y))
    t1 = time.time()
    nul, rk, n = loop_nullity(P12, P21, W, (1, 1, -2))
    print(f'  hypothesis (ii), a_2, b-, {"generic" if t[0] != 1 else "s=0"}: V(y1)V(y2)W -> W V(y2)V(y1): unknowns {n}, rank {rk}, '
          f'nullity {nul}  ({time.time() - t1:.1f} s)')


# ------------------------------------------------------------------ (A') existence at the rational point: the closed forms satisfy the
# sector equations T_{Q'} A_j(Q) = B_j(Q) T_Q exactly (b^-, V and Vbar; b^+, Vbar), so the bound is attained.
def closed_form(h, q, X, t, Q, kind):
    """eq-dt-tv (kind 'V'), eq-dt-tvbar ('Vbar'), eq-dt-tplus ('Vbar+', X = X_+) on the sector Q: matrix (target p) x (source p)."""
    def qn(N, c):
        return t[c]*qpow(q, N[c])
    T = [[Fr(0)]*h for _ in range(h)]
    for i in range(h):
        if kind == 'V':
            N = tuple(Q[c] - (1 if c == i else 0) for c in range(h))
            T[i][i] = 1/qn(N, i); T[(i - 1) % h][i] = X
        elif kind == 'Vbar':
            N = tuple(Q[c] + (1 if c == i else 0) for c in range(h))
            pw = Fr(1)
            for j in range(h):
                b = (i + j) % h; pw *= qn(N, b)
                T[b][i] = X**j*pw
        else:
            N = tuple(Q[c] + (1 if c == i else 0) for c in range(h))
            T[i][i] = 1/qn(N, i); b = (i - 1) % h
            T[b][i] = X*qn(N, i)/qn(N, b)/q
    return T


ok = True
for h in (3, 4, 5):
    t = (Fr(3, 2), Fr(5, 7)) + (Fr(1),)*(h - 3) + (Fr(14, 15),)
    for kind, sign, conj in (('V', -1, False), ('Vbar', -1, True), ('Vbar+', 1, True)):
        W = Weyl(h, q, m, t, sign); P = particle(h, q, y, conj)
        X = -(q - 1/q)*m*y if sign == -1 else (q - 1/q)*m/y
        for Q in (tuple([1] + [0]*(h - 1)), tuple([2, -1] + [0]*(h - 3) + [0])):
            for j in range(h):
                A, Qn = sector_map(P, W, j, Q, True); B, _ = sector_map(P, W, j, Q, False)
                lhs = mat_mul(closed_form(h, q, X, t, Qn, kind), A); rhs = mat_mul(B, closed_form(h, q, X, t, Q, kind))
                ok &= lhs == rhs
print(f"(A') closed forms eq-dt-tv, eq-dt-tvbar (b-) and eq-dt-tplus (b+) satisfy the sector equations exactly, h = 3, 4, 5, "
      f"generic offset: {ok}   [X = -(q - 1/q) m y for b-, X_+ = (q - 1/q) m/y for b+]")


# ------------------------------------------------------------------ (B) structure of Omega (numerical, random complex parameters)
def num_rep(rep):
    out = {}
    for k_, v in rep.items():
        out[k_] = np.array(v, complex) if isinstance(v, list) and k_ not in ('wt',) else v
    return out


def omega_structure(h, sign, conj, rng):
    qc = np.exp(1j*rng.uniform(0.3, 1.2))*rng.uniform(0.8, 1.2)
    mc = complex(rng.normal(), rng.normal()); yc = complex(rng.normal(), rng.normal())
    tc = np.exp(rng.normal(size=h) + 1j*rng.normal(size=h)); tc = tc/np.prod(tc)**(1/h)
    # complex version of the sector maps: reuse the exact code with python complex numbers
    W = Weyl(h, qc, mc, tuple(tc), sign)
    P = particle(h, qc, yc, conj)
    res = []
    Q = tuple(rng.integers(-2, 3, size=h - 1).tolist()); Q = Q + (-sum(Q),)
    Q = tuple(c + (1 if i == 0 else 0) for i, c in enumerate(Q))
    for j in range(h):
        k = (j + 1) % h
        first = np.array(path(P, W, Q, [k, j], True)[0], complex)     # x_j x_{j+1}: node j+1 acts first
        second = np.array(path(P, W, Q, [j, k], True)[0], complex)    # x_{j+1} x_j
        Om = first - (second/qc if sign == -1 else qc*second)
        sv = np.linalg.svd(Om, compute_uv=False)
        rank1 = sv[1]/sv[0] < 1e-12
        nzr = [i for i in range(h) if np.abs(Om[i]).max() > 1e-12*np.abs(Om).max()]       # image support
        nzc = [i for i in range(h) if np.abs(Om[:, i]).max() > 1e-12*np.abs(Om).max()]    # columns acting (kernel = complement)
        res.append((j, rank1, nzr, nzc))
    return res


rng = np.random.default_rng(7)
print('(B) Omega_j on sectors (rank 1; image rows / acting columns, 0-based, j = node):')
for h in (3, 4, 5, 6):
    for sign in (-1, 1):
        for conj in (False, True):
            rr = omega_structure(h, sign, conj, rng)
            # predicted: b-: V image row j+2 (0-based index (j+1)%h for v_{j+2}) ; Vbar acting column only vbar_{j+2}
            okr = True
            for j, r1, nzr, nzc in rr:
                tgt = (j + 1) % h        # 0-based index of v_{j+2} (v_1 is index 0, node j moves v_j -> v_{j+1})
                if sign == -1:
                    okr &= r1 and ((nzr == [tgt]) if not conj else (nzc == [tgt]))
                else:
                    okr &= r1 and ((nzc == [tgt]) if not conj else (nzr == [tgt]))
            print(f'  h={h} {"b+" if sign == 1 else "b-"} {"Vbar" if conj else "V   "}: rank 1 and coordinate image/kernel as '
                  f'predicted: {okr}   e.g. j=0: rows {rr[0][2]}, cols {rr[0][3]}')

# ------------------------------------------------------------------ (C) determinants on a sector, symbolic
qs, Xs = sp.symbols('q X', nonzero=True)
print('(C) determinants of the closed forms on the sector Q (symbolic q, X, offsets t_i, prod t_i = 1):')
for h in (3, 4):
    ts = sp.symbols(f't1:{h + 1}', nonzero=True)
    subs_last = {ts[-1]: 1/sp.Mul(*ts[:-1])}
    for kind in ('V', 'Vbar', 'Vbar+'):
        dets = set()
        for Q in (tuple([1] + [0]*(h - 1)), tuple([2, -1] + [0]*(h - 2)), tuple([0, 1, 1] + [-1]*(h - 3) + ([] if h == 3 else []))[:h]):
            Q = tuple(Q[:h]); Q = Q[:-1] + (1 - sum(Q[:-1]),)
            T = sp.zeros(h, h)
            for i in range(h):
                if kind == 'V':
                    N = [Q[c] - (1 if c == i else 0) for c in range(h)]
                    T[i, i] = 1/(ts[i]*qs**N[i]); T[(i - 1) % h, i] = Xs
                elif kind == 'Vbar':
                    N = [Q[c] + (1 if c == i else 0) for c in range(h)]
                    pw = 1
                    for j in range(h):
                        b = (i + j) % h; pw *= ts[b]*qs**N[b]; T[b, i] = Xs**j*pw
                else:
                    N = [Q[c] + (1 if c == i else 0) for c in range(h)]
                    T[i, i] = 1/(ts[i]*qs**N[i]); b = (i - 1) % h
                    T[b, i] = Xs*ts[i]*qs**N[i]/(ts[b]*qs**N[b])/qs
            dets.add(sp.factor(sp.simplify(T.det().subs(subs_last))))
        print(f'  h={h} {kind:5s}: {dets}')

# ------------------------------------------------------------------ (D) root of unity
h = 3; qr = np.exp(2j*np.pi/5); mr = 0.7 - 0.2j
Wr = Weyl(h, qr, mr, (1.0, 1.0, 1.0), -1)
res = 0
for N in [(1, -1, 0), (0, 2, -2), (3, 0, -3)]:
    for j in range(h):
        c1, N1 = Wr.x(j, N)
        NS = tuple(n + 5*a for n, a in zip(N, alpha(h, 1)))
        c2, N2 = Wr.x(j, NS)
        res = max(res, abs(c1 - c2), abs(Wr.k(j, N) - Wr.k(j, NS)))
        assert N2 == tuple(n + 5*a for n, a in zip(N1, alpha(h, 1)))
print(f'(D) q = e^(2 pi i/5) (omega = 0.6 gives q = -e^(-0.6 i pi) = {np.exp(1j*np.pi*0.4):.4f}): S|N> = |N + 5 alpha_1> commutes with '
      f'k_j and f_j [{res:.1e}]; (S x 1) T is a second intertwiner')

# ------------------------------------------------------------------ (E) reversed orientation
print('(E) orientation -1 (f_j|n> = m q^{-n_{j+1}}|n - alpha_j>): Vbar intertwiner and its CZ2007 gauge')
for h in (3, 4, 5):
    t = (Fr(3, 2), Fr(5, 7)) + (Fr(1),)*(h - 3) + (Fr(14, 15),)
    W = Weyl(h, q, m, t, -1, orient=-1); P = particle(h, q, y, True)
    X = -(q - 1/q)*m*y
    def Tm(Q):
        T = [[Fr(0)]*h for _ in range(h)]
        for i in range(h):
            N = tuple(Q[c] + (1 if c == i else 0) for c in range(h))
            qn = lambda c: t[c]*qpow(q, N[c])
            T[i][i] = qn(i); T[(i + 1) % h][i] = X*qn(i)/qn((i + 1) % h)
        return T
    ok = True
    for Q in (tuple([1] + [0]*(h - 1)), tuple([2, -1] + [0]*(h - 2))):
        for j in range(h):
            A, Qn = sector_map(P, W, j, Q, True); B, _ = sector_map(P, W, j, Q, False)
            ok &= mat_mul(Tm(Qn), A) == mat_mul(B, Tm(Q))
    print(f'  h={h}: M^-_Vbar(vb_i|n>) = q^(n_i)|n>vb_i + X q^(n_i - n_(i+1))|n - alpha_i>vb_(i+1) intertwines exactly: {ok}')
# gauge |n>' = q^{-|n|^2/2}|n>: converting coefficient X q^{n_i-n_{i+1}} q^{(|n-alpha_i|^2-|n|^2)/2} = X q^{n_i-n_{i+1}} q^{-(n_i-n_{i+1})+1}
nsym = sp.symbols('n1:5'); qq, XX = sp.symbols('q X')
for h in (3, 4):
    for i in range(h):
        n = sp.Matrix(nsym[:h]); al = sp.zeros(h, 1); al[i] = 1; al[(i + 1) % h] = -1
        expo = sp.expand(((n - al).dot(n - al) - n.dot(n))/2 + n[i] - n[(i + 1) % h])
        assert expo == 1, expo
print('  in the basis q^(-|n|^2/2)|n> the converting coefficient is q X (exponent identity checked for h = 3, 4, all i); '
      'the diagonal q^(n_i) = Q_CZ^(-alpha.l_i) is unchanged: CZ2007 (5.10) with xhat = q X.  To keep xhat = X: mu -> mu/q.')
print(f'total {time.time() - T0:.1f} s')
