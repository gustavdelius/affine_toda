"""Defects meet boundaries (Part VIII, ch. 42: sec-doubled-theory, sec-diagonal-qsp, sec-jump-coideal,
sec-defect-boundary-fusion, sec-tensor-k).

Book ids supported: prp-reflected-copy (a), prp-diagonal-coideal (b), cnj-twisted-diagonal and prp-diagonal-borel (c),
prp-defect-fusion-k (d), prp-fused-bdg (e), sec-tensor-k (f).  Exercise exr-dbd-mixed-twist uses (c).

Conventions (conventions brief part8_defects_conventions.md; boundary_common.py, defect_common.py):
  Part VI coproduct Delta(e) = e (x) 1 + k (x) e, Delta(f) = f (x) k^-1 + 1 (x) f; q = q_phys = -exp(-i pi lambda);
  y = exp(lambda theta); T(theta): V (x) W -> W (x) V, Tbar(theta) = T(-theta)^-1 : W (x) V(-theta) -> V(-theta) (x) W.
  Doubled theory: M(theta) = V(theta) + V'(theta); the reflected copy V'(u) carries pi_{V(-u)} = pi_{V(1/y)} and the
  algebra acts on it through Delta^op, so the doubled bulk algebra is A = U (x) U^cop with
      Delta_A(g (x) h) = sum (g_(1) (x) h_(2)) (x) (g_(2) (x) h_(1)).
  g (x) h acts on M as diag(pi(g) eps(h), eps(g) pi'(h)).  Doubled R-check: Rc on V V, P Rc P on V' V', flip on mixed
  pairs (eq-doubled-k, sec-doubling-trick).  Doubled K (eq-doubled-k) on (V + V') (x) W, ordered (V (x) W) + (V' (x) W):
  V (x) W -> V' (x) W is P T(theta), V' (x) W -> V (x) W is Tbar(theta) P.  The boundary module W sits on the right, so a
  symmetry is a LEFT coideal B of A, acting on M (x) W by (sigma (x) chi) Delta_A.
  Windows: W = span{|n>, |n| <= N}; equations are imposed on labels |n| <= N - 2 only.

Checks:
 a. Reflected copy (prp-reflected-copy): the doubled R-check intertwines Delta_A of U (x) U^cop for sine-Gordon and the
    a_2 vector solitons, and fails with U (x) U; V(1/y) o omega_q = Vbar(y) exactly (omega_q(e) = f, omega_q(k) = k^-1),
    and P Rc P = Rc_{Vbar Vbar}: the reflected soliton is the conjugate multiplet.  The equivalence "doubled reflection
    equation <=> triangle relations" (prp-doubling item 1) is checked in defect_axioms.py (check d) and not repeated.
 b. Transparent defect and invertible lines (prp-diagonal-coideal): in the space X spanned by g (x) h, g at node j and h at
    node j' in {k^-1, 1, k, e, f}, total degree <= 1 (21 elements), the commutant of the doubled K on one-particle states
    and on the two-particle doubled reflection is larger than a coideal; its largest left coideal (iterate
    S -> {x in S : Delta_A(x) in A (x) S}) is exactly (1 (x) psi) Delta(U) cap X, for W = C and T = 1 (psi = id;
    sine-Gordon with spin 1/2 and spin 1 probes, a_2 vector), T = exp(i alpha mu) (psi = torus), T = C (psi = node swap tau).
 c. Jump defect (cnj-twisted-diagonal, prp-diagonal-borel), sine-Gordon types I_- (Weyl module of U_q(b^-)) and I_+:
    the diagonal Borel Delta(U_q(b^-+)) is a symmetry with W the Weyl module (s = 0 and 0.3); for the generators of every
    twisted diagonal (1 (x) psi)Delta(U), psi in {torus, tau o torus} (all Hopf automorphisms, Kolb 1207.6036 Thm. 3.2),
    a least-squares search over theta-independent weight-graded operators on W (three rapidities, windows N = 5, 7, 9)
    finds: type I_-: e/id and f/tau fail, f/id and e/tau hold; type I_+: the mirror image.  Exactly (sympy,
    noncommutative symbols for arbitrary operators on W): with L = P T = D + y O (I_-) or D + O/y (I_+), one power of y
    in a single generator equation has a block that contains no unknown and does not vanish (prp-diagonal-borel (3)).
    The second Borel survives twisted by tau, with an explicit second module structure rho~ on W (relations and the
    intertwining of all four generators checked numerically on windows and exactly for symbolic n).  The two Borel
    symmetries do not combine: rho(k_1) and rho~(k_1) q^4-commute although k_1 (x) k_1 and k_1 (x) k_0 commute, and any
    coideal containing both contains e_1 (x) 1 (type I_-) or 1 (x) f_0 (type I_+), which are not symmetries (numerically,
    and exactly by a power of y).
    Limits (prp-diagonal-borel (5)): the off-diagonal part of L is (y/y_d) x const for I_- and (y_d/y) x const for I_+,
    so the diagonal limit is theta_d -> +inf for I_- and theta_d -> -inf for I_+ (single-generator residuals decay like
    e^{-lambda |theta_d|} in these directions and tend to 1 in the opposite ones).  In neither limit is there a twisted
    diagonal on all of W: the e_1 and f_1 equations force chi(Delta k_1) and chi(Delta k_1^-1) to be the same operator
    t^-1 q^(-+(2n+s)), not inverse ones, and in the off-diagonal limit chi(k_1 (x) k_0) and its would-be inverse both raise
    n by 2 (exact).  Each charge sector C|n> of the diagonal limit is the twist line t_1 = 1/t_0 = q^(-2n-s).
    Control: a particle used as a defect (W = V(theta_d), T = Rc) has the full diagonal Delta(U).
 d. Fusion onto a Dirichlet boundary (prp-defect-fusion-k): K_f = Tbar (1 (x) K_D) T on V (x) W solves the reflection
    equation on windows N = 5, 7, 9 (perturbed control fails); its charge-sector blocks have label-independent gauge
    invariants (Corrigan-Zambon 1202.6016 (2.18)) and equal the Ghoshal-Zamolodchikov matrix eq-bsg-k-matrix with
    (eta, vartheta) = (eta_D, lambda theta_d); with the scalar factors (eq-sgd-scalar; R_0, sigma of eq-bsg-k-full) the
    identity holds entry by entry.  Sinh-Gordon: Tbar K_D T = K_shG(E, F = 2 i theta_d/pi) (eq-shg-reflection) re-checked.
    Corrigan-Zambon 1202.6016 in their conventions (rows = incoming soliton, Q^2 = 1/q): direct fusion (2.9) of their
    general transmission matrix (2.14) with the Dirichlet matrix (2.16) agrees with their (2.17); for the type-II
    parameters it agrees with (2.19) only after two corrections, h = -(|b_-|^2 + |b_+|^2) Q^-2 and the opposite overall
    sign of s_- (each alone fails); the charge dependence of the type-II result survives.
 e. Classical (prp-fused-bdg): K_phi(lam) = sigma(L_d(1/lam)^-1) K_psi(lam) L_d(lam)^-1 (sigma = transposition) satisfies the
    dynamical compatibility d_t K = K a_t(lam) + sigma(a_t(1/lam)) K on the fused system (defect + psi-boundary with
    field-independent K_psi: GZ for a_1, Durham C = +-1 and Neumann for a_2, a_3); the undressed or wrongly dressed K fails.
    The sine-Gordon defect potential is the Baseilhac-Delius-George H_osc (eq-bsg-oscillator) at
    (p^, q^) = (beta psi/2 - i theta_d, beta psi/2 + i theta_d) (sympy); the BDG flow does not preserve p^ - q^.
 f. Tensor K (sec-tensor-k): for the sine-Gordon coideal with generic eps-hat and W = V(theta_d) (a particle frozen as a
    defect), Tbar (1 (x) K) T intertwines the coideal on V(theta) (x) W, the space of such intertwiners is one-dimensional
    (so it is Appel-Vlaar's trigonometric tensor K-matrix), and it solves the reflection equation.
Runtime: about 15 s (20 checks).
"""
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.linalg import expm, null_space as sp_null
from boundary_common import report, q_phys, sl2_spin, an_vector, rcheck, prop_residual, cop, coideal_gen, \
    coideal_gen_tensor, intertwiners
from defect_common import fock_module_b, borel_l_operator, flip_VW, sg_tmatrix, an_data, dmatrix, \
    defect_potential, a_t, an_lax_rep

kron = np.kron
LAM = 1.37
q = q_phys(LAM)
rng = np.random.default_rng(42)
THETAS = [0.4 + 0.1j, -0.7 + 0.05j, 0.9 - 0.08j]
PAIRS = [(0.64 + 0.05j, 0.21 - 0.08j), (0.93 + 0.02j, -0.37 + 0.06j)]
sg = lambda y: sl2_spin(1, q, y)


def flip(d1, d2):
    """P: C^d1 (x) C^d2 -> C^d2 (x) C^d1."""
    P = np.zeros((d1*d2, d1*d2))
    for a in range(d1):
        for b in range(d2):
            P[b*d1 + a, a*d2 + b] = 1
    return P


def bdiag(A, B):
    Z = np.zeros((A.shape[0] + B.shape[0],)*2, complex)
    Z[:A.shape[0], :A.shape[0]] = A
    Z[A.shape[0]:, A.shape[0]:] = B
    return Z


# ================================================================ the doubled algebra U (x) U^cop on M = V + V'
class Doubled:
    def __init__(self, rep_fn):
        self.f, self.cache = rep_fn, {}

    def rep(self, y):
        k = complex(y)
        if k not in self.cache:
            self.cache[k] = self.f(y)
        return self.cache[k]

    def sig(self, y, which, letter):
        """first (which=1) or second (which=2) tensor factor of A, acting on M(y)."""
        V, Vp = self.rep(y), self.rep(1/y)
        d = V['dim']
        eps = np.eye(d) if letter[0] == 'k' else np.zeros((d, d))
        return bdiag(V[letter], eps) if which == 1 else bdiag(eps, Vp[letter])

    def tens(self, y1, y2, which, letter, cop_kind='std'):
        """(sigma (x) sigma) of Delta (which=1) or Delta^op (which=2; cop_kind='wrong' uses Delta) on M(y1) (x) M(y2)."""
        s1 = lambda L: self.sig(y1, which, L)
        s2 = lambda L: self.sig(y2, which, L)
        kind, j = letter
        K1, K2 = s1(('k', j)), s2(('k', j))
        I1, I2 = np.eye(len(K1)), np.eye(len(K2))
        op = which == 2 and cop_kind == 'std'
        if kind == 'k':
            return kron(K1, K2)
        if kind == 'e':
            return kron(I1, s2(letter)) + kron(s1(letter), K2) if op else kron(s1(letter), I2) + kron(K1, s2(letter))
        return kron(np.linalg.inv(K1), s2(letter)) + kron(s1(letter), I2) if op else \
            kron(s1(letter), np.linalg.inv(K2)) + kron(I1, s2(letter))

    def R(self, y1, y2):
        """doubled R-check M(y1) (x) M(y2) -> M(y2) (x) M(y1)."""
        d = self.rep(y1)['dim']
        D = 2*d
        blocks = {0: rcheck(self.rep(y1), self.rep(y2))[0],
                  1: flip(d, d) @ rcheck(self.rep(1/y2), self.rep(1/y1))[0] @ flip(d, d)}
        out = np.zeros((D*D, D*D), complex)
        for i in range(D):
            for j in range(D):
                ci, cj, li, lj = i//d, j//d, i % d, j % d
                if ci == cj:
                    for k in range(d):
                        for l in range(d):
                            out[(ci*d + k)*D + ci*d + l, i*D + j] = blocks[ci][k*d + l, li*d + lj]
                else:
                    out[j*D + i, i*D + j] = 1
        return out


# ================================================================ a. the reflected copy
def bulk_intertwining(DD, nodes, y1, y2, cop_kind='std'):
    Rd = DD.R(y1, y2)
    worst = 0
    for which in (1, 2):
        for j in nodes:
            for kind in 'efk':
                X12 = DD.tens(y1, y2, which, (kind, j), cop_kind)
                X21 = DD.tens(y2, y1, which, (kind, j), cop_kind)
                worst = max(worst, np.abs(Rd @ X12 - X21 @ Rd).max()/np.abs(X12).max())
    return worst


y1, y2 = np.exp(LAM*(0.71 + 0.13j)), np.exp(LAM*(-0.29 + 0.05j))
a2 = lambda y: an_vector(2, q, y)
res_ok = max(bulk_intertwining(Doubled(sg), (0, 1), y1, y2), bulk_intertwining(Doubled(a2), (0, 1, 2), y1, y2))
res_bad = min(bulk_intertwining(Doubled(sg), (0, 1), y1, y2, 'wrong'),
              bulk_intertwining(Doubled(a2), (0, 1, 2), y1, y2, 'wrong'))
report(f'doubled R-check (Rc on V V, P Rc P on V\'V\', flip on mixed pairs) intertwines Delta_A of U (x) U^cop for '
       f'sine-Gordon and a_2 [{res_ok:.0e}]; with U (x) U instead it fails [{res_bad:.2f}]', res_ok < 1e-12 and res_bad > 0.1)


def omega_pull(rep, nodes):
    return {('e', j): rep[('f', j)] for j in nodes} | {('f', j): rep[('e', j)] for j in nodes} | \
        {('k', j): np.linalg.inv(rep[('k', j)]) for j in nodes}


yy = np.exp(LAM*0.4)
W_a2 = omega_pull(a2(1/yy), (0, 1, 2))
Vb = an_vector(2, q, yy, conj=True)
d_om = max(np.abs(W_a2[k] - Vb[k]).max() for k in W_a2)
Cx = np.array([[0, 1], [1, 0]])
W_sg = omega_pull(sg(1/yy), (0, 1))
d_om_sg = max(np.abs(W_sg[k] - Cx @ sg(yy)[k] @ Cx).max() for k in W_sg)
d_R = max(prop_residual(flip(3, 3) @ rcheck(a2(1/y2), a2(1/y1))[0] @ flip(3, 3),
                        rcheck(an_vector(2, q, y1, conj=True), an_vector(2, q, y2, conj=True))[0])[0],
          prop_residual(flip(2, 2) @ rcheck(sg(1/y2), sg(1/y1))[0] @ flip(2, 2),
                        kron(Cx, Cx) @ rcheck(sg(y1), sg(y2))[0] @ kron(Cx, Cx))[0])
report(f'reflected copy = conjugate multiplet: V(1/y) o omega_q = Vbar(y) exactly (a_2 [{d_om:.0e}]; sine-Gordon: '
       f'C V(y) C [{d_om_sg:.0e}]), and P Rc(1/y2, 1/y1) P = Rc on Vbar(y1) (x) Vbar(y2) up to a scalar [{d_R:.0e}]',
       max(d_om, d_om_sg, d_R) < 1e-12)


# ================================================================ b. transparent defect, invertible lines: maximal coideal
def basis_X(j1, j2):
    """g letters at node j1 (first factor), h letters at node j2 (second factor), total degree <= 1: 21 monomials.
    A monomial (kind, j, a) is x_j k_j^a with x in {1, e, f}."""
    G = [('1', j1, a) for a in (-1, 0, 1)]
    H = [('1', j2, a) for a in (-1, 0, 1)]
    return ([(g, h) for g in G for h in H] + [(('e', j1, 0), h) for h in H] + [(('f', j1, 0), h) for h in H]
            + [(g, ('e', j2, 0)) for g in G] + [(g, ('f', j2, 0)) for g in G])


def mono_act(DD, y, which, mono, two=None):
    kind, j, a = mono
    if two is None:
        K = DD.sig(y, which, ('k', j))
        X = DD.sig(y, which, (kind, j)) if kind != '1' else None
    else:
        K = DD.tens(y, two, which, ('k', j))
        X = DD.tens(y, two, which, (kind, j)) if kind != '1' else None
    Ka = np.linalg.matrix_power(K if a >= 0 else np.linalg.inv(K), abs(a))
    return Ka if X is None else X @ Ka


def commutant(DD, Kfun, B):
    """null space of: K x_in = x_out K (one particle, THETAS) and KK X_in = X_out KK (two particles, PAIRS)."""
    D = 2*DD.rep(1.0)['dim']
    KKs = []
    for t1, t2 in PAIRS:
        Y1, Y2 = np.exp(LAM*t1), np.exp(LAM*t2)
        KKs.append(kron(np.eye(D), Kfun(Y2)) @ DD.R(Y2, 1/Y1) @ kron(np.eye(D), Kfun(Y1)) @ DD.R(Y1, Y2))
    cols1, cols2 = [], []
    for g, h in B:
        e1, e2 = [], []
        for th in THETAS:
            y = np.exp(LAM*th)
            si = mono_act(DD, y, 1, g) @ mono_act(DD, y, 2, h)
            so = mono_act(DD, 1/y, 1, g) @ mono_act(DD, 1/y, 2, h)
            e1.append((Kfun(y) @ si - so @ Kfun(y)).ravel())
        for (t1, t2), KK in zip(PAIRS, KKs):
            Y1, Y2 = np.exp(LAM*t1), np.exp(LAM*t2)
            Ti = mono_act(DD, Y1, 1, g, Y2) @ mono_act(DD, Y1, 2, h, Y2)
            To = mono_act(DD, 1/Y1, 1, g, 1/Y2) @ mono_act(DD, 1/Y1, 2, h, 1/Y2)
            e2.append((KK @ Ti - To @ KK).ravel())
        cols1.append(np.concatenate(e1))
        cols2.append(np.concatenate(e2))
    E1, E2 = np.array(cols1).T, np.array(cols2).T
    return nullsp(E1), nullsp(np.vstack([E1, E2]))


def nullsp(M, rtol=1e-9):
    _, s, Vh = np.linalg.svd(M)
    return Vh[int((s > rtol*s.max()).sum()):].conj().T


def cop_mono(m):
    kind, j, a = m
    if kind == '1':
        return [(1, ('1', j, a), ('1', j, a))]
    if kind == 'e':
        return [(1, ('e', j, a), ('1', j, a)), (1, ('1', j, a + 1), ('e', j, a))]
    return [(1, ('f', j, a), ('1', j, a - 1)), (1, ('1', j, a), ('f', j, a))]


def max_left_coideal(S, B, tol=1e-9):
    """largest subspace of span(S) (columns = coefficient vectors over B) with Delta_A(S) in A (x) S."""
    idx = {b: i for i, b in enumerate(B)}
    firsts = {}
    for i, (g, h) in enumerate(B):
        for c1, g1, g2 in cop_mono(g):
            for c2, h1, h2 in cop_mono(h):
                firsts.setdefault((g1, h2), np.zeros((len(B), len(B)), complex))[idx[(g2, h1)], i] += c1*c2
    while S.shape[1]:
        Q, _ = np.linalg.qr(S)
        Pc = np.eye(len(B)) - Q @ Q.conj().T
        M = np.vstack([Pc @ Y @ S for Y in firsts.values()])
        _, s, Vh = np.linalg.svd(M)
        rank = int((s > tol*max(1.0, s.max())).sum())
        if rank == 0:
            break
        S = S @ Vh[rank:].conj().T
    return S


def twisted_delta(B, j1, j2, t):
    """(1 (x) psi) Delta(U) cap X, psi(e_j1) = t e_j2, psi(f_j1) = f_j2/t, psi(k_j1) = k_j2."""
    idx = {b: i for i, b in enumerate(B)}

    def v(terms):
        c = np.zeros(len(B), complex)
        for coef, b in terms:
            c[idx[b]] += coef
        return c
    return np.array([v([(1, (('1', j1, a), ('1', j2, a)))]) for a in (-1, 0, 1)]
                    + [v([(1, (('e', j1, 0), ('1', j2, 0))), (t, (('1', j1, 1), ('e', j2, 0)))]),
                       v([(1, (('f', j1, 0), ('1', j2, -1))), (1/t, (('1', j1, 0), ('f', j2, 0)))])]).T


def same_span(A, B):
    if A.shape[1] != B.shape[1]:
        return np.inf
    Q, _ = np.linalg.qr(A)
    return np.linalg.norm(B - Q @ (Q.conj().T @ B))/np.linalg.norm(B)


def Kswap(T, Tb):
    d = T.shape[0]
    K = np.zeros((2*d, 2*d), complex)
    K[d:, :d], K[:d, d:] = T, Tb
    return K


al = 0.37
Ttw = np.diag([np.exp(1j*al), np.exp(-1j*al)])
cases = [('transparent, sine-Gordon spin 1/2', sg, lambda y: Kswap(np.eye(2), np.eye(2)), [(0, 0, 1), (1, 1, 1)]),
         ('transparent, sine-Gordon spin 1 probe', lambda y: sl2_spin(2, q, y),
          lambda y: Kswap(np.eye(3), np.eye(3)), [(1, 1, 1)]),
         ('transparent, a_2 vector', a2, lambda y: Kswap(np.eye(3), np.eye(3)), [(0, 0, 1), (1, 1, 1), (2, 2, 1)]),
         ('twist line exp(i alpha mu)', sg, lambda y: Kswap(Ttw, np.linalg.inv(Ttw)),
          [(0, 0, np.exp(-2j*al)), (1, 1, np.exp(2j*al))]),
         ('charge conjugation line C', sg, lambda y: Kswap(Cx, Cx), [(0, 1, 1), (1, 0, 1)])]
summary = []
ok_b = True
for name, fn, Kfun, nodes in cases:
    DD = Doubled(fn)
    for j1, j2, t in nodes:
        B = basis_X(j1, j2)
        N1, N12 = commutant(DD, Kfun, B)
        S1, S12 = max_left_coideal(N1, B), max_left_coideal(N12, B)
        Ex = twisted_delta(B, j1, j2, t)
        dev = max(same_span(S1, Ex), same_span(S12, Ex))
        ok_b &= dev < 1e-8
        summary.append((name, N1.shape[1], N12.shape[1], S12.shape[1], dev))
dims = sorted({(s[1], s[2], s[3]) for s in summary})
report('largest left coideal inside the commutant of the doubled K (X: 21 monomials of degree <= 1 per node pair): '
       'transparent defect -> Delta(U) cap X (sine-Gordon spin 1/2 and spin 1 probes, a_2); twist line -> '
       '(1 (x) t)Delta(U) cap X, t_1 = 1/t_0 = e^{2 i alpha}; charge conjugation -> (1 (x) tau)Delta(U) cap X '
       f'[max deviation {max(s[4] for s in summary):.0e}; commutant dims (one particle, with two particles, coideal) {dims}]',
       ok_b)


# ================================================================ c. the jump defect
def T_I(theta, theta_d, N, sign, s=0.0):
    """type I_- (sign=-1, b^-) / I_+ (sign=+1, b^+) transmission matrix part on the Weyl window, and its module."""
    W = fock_module_b(sign, q, np.exp(LAM*theta_d), N, s=s, kind='weyl')
    return flip_VW(2, W['dim']) @ borel_l_operator(sign, np.exp(LAM*theta), W), W


def doubled_K_W(Tf, theta, dW):
    n = 2*dW
    K = np.zeros((2*n, 2*n), complex)
    K[n:, :n] = flip_VW(2, dW).T @ Tf(theta)
    K[:n, n:] = np.linalg.inv(Tf(-theta)) @ flip_VW(2, dW)
    return K


def interior(W, margin=2):
    m = np.tile(np.abs(W['n']) <= W['n'].max() - margin, 2)
    return np.concatenate([m, m])


def borel_residual(Tf, W, sign, theta):
    dW = W['dim']
    m = interior(W)
    blk = np.ix_(m, m)
    K = doubled_K_W(Tf, theta, dW)
    worst = 0
    for kd in (('e', 'k') if sign == 1 else ('f', 'k')):
        for j in (0, 1):
            def X(th):
                y = np.exp(LAM*th)
                return bdiag(cop(sg(y), W, kd, j), flip_VW(2, dW).T @ cop(W, sg(1/y), kd, j) @ flip_VW(2, dW))
            worst = max(worst, np.abs((K @ X(theta) - X(-theta) @ K)[blk]).max()/np.abs((K @ X(theta))[blk]).max())
    return worst


worst_b = 0
for sign in (-1, 1):
    for N in (6, 8):
        for s in (0.0, 0.3):
            W = T_I(0.1, 0.15, N, sign, s)[1]
            worst_b = max(worst_b, borel_residual(lambda th: T_I(th, 0.15, N, sign, s)[0], W, sign, 0.4 + 0.1j))
report(f'prp-diagonal-borel: the doubled K of the type-I_-+ defect commutes with the diagonal Borel Delta(U_q(b^-+)) '
       f'acting on M (x) W through Delta_A, W the Weyl module (windows N = 6, 8; s = 0, 0.3) [{worst_b:.0e}]',
       worst_b < 1e-12)

DE, DF = {1: 1, 0: -1}, {1: -1, 0: 1}          # shift of the defect label n under e_j, f_j


def gmat(kind, j, p=1):
    return lambda V: np.linalg.matrix_power(V[(kind, j)] if p >= 0 else np.linalg.inv(V[(kind, j)]), abs(p))


ONE = lambda V: np.eye(2)


def generator_terms(kind, j, jt):
    """(sigma (x) chi) Delta_A of x = (1 (x) psi)Delta(x_j), psi(x_j) ~ x_jt (torus factors absorbed in kap):
    e: on V (x) W: pi(e_j)(x)1 + pi(k_j)(x)Y,     on V' (x) W: 1(x)Y + pi'(e_jt)(x)kap,
    f: on V (x) W: pi(f_j)(x)kap + 1(x)Y,        on V' (x) W: pi'(k_jt^-1)(x)Y + pi'(f_jt)(x)1.
    Entries: (matrix on V, unknown name or None for the identity on W, shift of n)."""
    if kind == 'e':
        return ([(gmat('e', j), None, 0), (gmat('k', j), 'Y', DE[j])],
                [(ONE, 'Y', DE[j]), (gmat('e', jt), 'kap', DE[j] - DE[jt])])
    return ([(gmat('f', j), 'kap', DF[jt] - DF[j]), (ONE, 'Y', DF[jt])],
            [(gmat('k', jt, -1), 'Y', DF[jt]), (gmat('f', jt), None, 0)])


def ls_symmetry(Tf, W, termsV, termsP, thetas=THETAS):
    """min over theta-independent shift operators U_t on W of |K X(theta) - X(-theta) K| / |inhomogeneous part|."""
    dW = W['dim']
    n = 2*dW
    m = interior(W)
    blk = np.ix_(m, m)
    labels = list(W['n'])
    names = []
    for _, nm, d in termsV + termsP:
        if nm is not None and (nm, d) not in names:
            names.append((nm, d))

    def X(yy, which):
        Z = np.zeros((2*n, 2*n), complex)
        for cp, terms, yv in ((0, termsV, yy), (1, termsP, 1/yy)):
            V = sg(yv)
            for fn, nm, d in terms:
                if which is None and nm is None:
                    Z[cp*n:(cp + 1)*n, cp*n:(cp + 1)*n] += kron(fn(V), np.eye(dW))
                elif which is not None and nm == which[0]:
                    Z[cp*n:(cp + 1)*n, cp*n:(cp + 1)*n] += kron(fn(V), which[1])
        return Z
    A_rows, b_rows = [], []
    for th in thetas:
        K = doubled_K_W(Tf, th, dW)
        y = np.exp(LAM*th)
        b_rows.append(-(K @ X(y, None) - X(1/y, None) @ K)[blk].ravel())
        cols = []
        for nm, d in names:
            for c, lab in enumerate(labels):
                if lab + d in labels:
                    E = np.zeros((dW, dW))
                    E[labels.index(lab + d), c] = 1
                    cols.append((K @ X(y, (nm, E)) - X(1/y, (nm, E)) @ K)[blk].ravel())
        A_rows.append(np.array(cols).T)
    A, bvec = np.vstack(A_rows), np.concatenate(b_rows)
    sol = np.linalg.lstsq(A, bvec, rcond=1e-12)[0]
    return np.linalg.norm(A @ sol - bvec)/np.linalg.norm(bvec)


table = {}
for sign in (-1, 1):
    for N in (5, 7, 9):
        W = T_I(0.1, 0.15, N, sign)[1]
        Tf = lambda th: T_I(th, 0.15, N, sign)[0]
        for kind in 'ef':
            for tw in ('id', 'tau'):
                r = max(ls_symmetry(Tf, W, *generator_terms(kind, j, j if tw == 'id' else 1 - j)) for j in (0, 1))
                table.setdefault((sign, kind, tw), []).append(r)
hold = {(-1, 'f', 'id'), (-1, 'e', 'tau'), (1, 'e', 'id'), (1, 'f', 'tau')}
ok_hold = max(max(table[k]) for k in hold) < 1e-12
fail_min = min(min(v) for k, v in table.items() if k not in hold)
lab = lambda k: f"{'I-' if k[0] == -1 else 'I+'} {k[1]}/{k[2]}"
report('cnj-twisted-diagonal, least-squares test: for no twisted diagonal (1 (x) psi)Delta(U) (psi = torus or '
       'tau o torus) do both its e- and f-generators commute with the doubled jump-defect K for a theta-independent '
       'action on W of the weight allowed by charge conservation; minimal residuals over windows N = 5, 7, 9 (they '
       'depend on the sampled rapidities): '
       + ', '.join(f'{lab(k)} {min(v):.2f}-{max(v):.2f}' for k, v in table.items() if k not in hold),
       fail_min > 0.5)
report('what holds instead (same test): type I_- keeps f/id (diagonal b^-) and e/tau (tau-twisted b^+); type I_+ keeps '
       f'e/id and f/tau [{max(max(table[k]) for k in hold):.0e}, all windows]', ok_hold)

# exact obstruction (prp-diagonal-borel (3), (4), (5)): arbitrary operators on W as noncommuting symbols.
# Basis v_+ (x) W, v_- (x) W; L = P T = Dg + y Og (I_-) or Dg + Og/y (I_+), Dg = diag(A, D), Og = [[0, B], [C, 0]],
# with A, B, C, D theta-independent and invertible.  pi_y(e_1) = y E+, pi_y(f_1) = E-/y, pi_y(e_0) = y E-,
# pi_y(f_0) = E+/y, pi(k_1) = diag(q, 1/q) = pi(k_0)^-1; the reflected copy at -theta carries pi_y as well.
nA, nB, nC, nD, nY, nK, nZ = sp.symbols('A B C D Y kappa Z', commutative=False)
qs, ys, ts = sp.symbols('q y t', nonzero=True)
nE = {'+': sp.Matrix([[0, 1], [0, 0]]), '-': sp.Matrix([[0, 0], [1, 0]])}
I2s = sp.eye(2)
kmat = lambda op, p=1: sp.Matrix([[qs**p*op, 0], [0, qs**-p*op]])      # pi(k_1)^p (x) op
Dg, Og = sp.Matrix([[nA, 0], [0, nD]]), sp.Matrix([[0, nB], [nC, 0]])


def ycoef(M, p):
    """coefficient of y^p of a matrix polynomial in y, y^-1"""
    return M.applyfunc(lambda e: sp.expand(e*ys**4).coeff(ys, p + 4))


def gen_eq(L, kind, tw):
    """L X_V - X_V' L for (1 (x) psi)Delta(x_1), psi = t id or t tau (eq-dbd-coproduct; generator_terms)"""
    if kind == 'e':          # V: pi(e_1)(x)1 + pi(k_1)(x)Y ;  V': 1(x)Y + t pi(e_psi)(x)kappa
        XV = ys*nE['+'] + kmat(nY)
        XP = I2s*nY + ts*ys*(nE['+'] if tw == 'id' else nE['-'])*nK
    else:                    # V: pi(f_1)(x)kappa + 1(x)Y ;  V': pi(k_psi^-1)(x)Y + t^-1 pi(f_psi)(x)1
        XV = nE['-']/ys*nK + I2s*nY
        XP = kmat(nY, -1 if tw == 'id' else 1) + (nE['-'] if tw == 'id' else nE['+'])/(ts*ys)
    return (L*XV - XP*L).applyfunc(sp.expand)


Lm_s, Lp_s = Dg + ys*Og, Dg + Og/ys
cases_ex = [('I_- e/id', gen_eq(Lm_s, 'e', 'id'), 2, (1, 1), nC), ('I_- f/tau', gen_eq(Lm_s, 'f', 'tau'), -1, (0, 1), -nD/ts),
            ('I_+ f/id', gen_eq(Lp_s, 'f', 'id'), -2, (1, 1), -nB/ts), ('I_+ e/tau', gen_eq(Lp_s, 'e', 'tau'), 1, (0, 1), nA),
            ('I_- e_1 (x) 1', (Lm_s*(ys*nE['+'] + kmat(nZ)) - I2s*nZ*Lm_s), 2, (1, 1), nC),
            ('I_+ 1 (x) f_0', (Lp_s*I2s*nZ - (kmat(nZ) + nE['+']/ys)*Lp_s), -2, (0, 0), -nC),
            # in the limits the other twist fails already: diagonal limit L = Dg (tau), off-diagonal limit L = Og (id)
            ('diag. limit e/tau', gen_eq(Dg, 'e', 'tau'), 1, (0, 1), nA), ('off-diag. limit e/id', gen_eq(Og, 'e', 'id'), 1, (1, 1), nC)]
ex_ok, ex_txt = True, []
for name, M, p, (i, j), expect in cases_ex:
    blk = ycoef(M.applyfunc(sp.expand), p)[i, j]
    ex_ok &= sp.simplify(blk - expect) == 0
    ex_txt.append(f'{name}: y^{p} block {"+-"[i]}{"+-"[j]} = {blk}')
# diagonal limit L = Dg: the e_1 (y^1, block +-) and f_1 (y^-1, block -+) equations of the untwisted diagonal are
# A - t kappa D = 0 and D kappa' - A/t = 0: they force kappa = A D^-1/t = chi(Delta k_1), kappa' = D^-1 A/t = chi(Delta k_1^-1)
be = ycoef(gen_eq(Dg, 'e', 'id'), 1)[0, 1]
bf = ycoef(gen_eq(Dg, 'f', 'id'), -1)[1, 0]
lim_forms = sp.simplify(be - (nA - ts*nK*nD)) == 0 and sp.simplify(bf - (nD*nK - nA/ts)) == 0
# off-diagonal limit L = Og, tau-twisted: C - t kappa B = 0 and B kappa' - C/t = 0, so kappa = C B^-1/t and
# kappa' = B^-1 C/t; C raises n by 1 and B lowers it by 1, so both raise n by 2 and cannot be inverse to each other
be_t = ycoef(gen_eq(Og, 'e', 'tau'), 1)[1, 1]
bf_t = ycoef(gen_eq(Og, 'f', 'tau'), -1)[0, 0]
lim_forms &= sp.simplify(be_t - (nC - ts*nK*nB)) == 0 and sp.simplify(bf_t - (nB*nK - nC/ts)) == 0
ns_, ss_ = sp.symbols('n s')
prod_m = sp.simplify((qs**-ns_/qs**(ns_ + ss_))**2/ts**2)      # I_-: A = q^-n, D = q^(n+s) commute
prod_p = sp.simplify((qs**ns_/qs**(-ns_ - ss_))**2/ts**2)      # I_+: A = q^n, D = q^(-n-s)
report('exact obstruction for arbitrary theta-independent operators Y, kappa, Z on W (noncommuting symbols): '
       + '; '.join(ex_txt) + '.  Diagonal limit: the e_1 and f_1 equations of (1 (x) t)Delta(U) force chi(Delta k_1) = '
       f'A D^-1/t and chi(Delta k_1^-1) = D^-1 A/t, the same operator; their product is {prod_m} (I_-), {prod_p} (I_+), '
       'not 1; off-diagonal limit: chi(k_1 (x) k_0) = C B^-1/t and chi(k_1^-1 (x) k_0^-1) = B^-1 C/t both raise n by 2',
       ex_ok and lim_forms and sp.simplify(prod_m.subs(ns_, 1) - prod_m.subs(ns_, 0)) != 0
       and sp.simplify(prod_p.subs(ns_, 1) - prod_p.subs(ns_, 0)) != 0)


# explicit second module rho~ for type I_- (sg_tmatrix normalization A_n = q^-n, D_n = q^(n+s), ...)
def rho_tilde(yd, N, s):
    ns = np.arange(-N, N + 1)
    d = len(ns)
    a, b = 1j*yd*q**(-s)/(q**2 - 1), 1j*yd*q**s/(q**2 - 1)
    g = {k: np.zeros((d, d), complex) for k in (('e', 1), ('e', 0), ('k', 1), ('k', 0))}
    for i, n in enumerate(ns):
        if i + 1 < d:
            g[('e', 1)][i + 1, i] = a*q**(-2*n)
        if i - 1 >= 0:
            g[('e', 0)][i - 1, i] = b*q**(2*n)
        if i + 2 < d:
            g[('k', 1)][i + 2, i] = q**(-2*n - 2 - s)
        if i - 2 >= 0:
            g[('k', 0)][i - 2, i] = q**(2*n - 2 + s)
    g.update(dim=d, n=ns)
    return g


N, s_, thd = 8, 0.3, 0.15
rt = rho_tilde(np.exp(LAM*thd), N, s_)
ns = rt['n']
inner = np.abs(ns) <= N - 3
blk = np.ix_(inner, inner)
e1, e0, k1, k0 = rt[('e', 1)], rt[('e', 0)], rt[('k', 1)], rt[('k', 0)]
rel = max(np.abs((k1 @ k0 - np.eye(len(ns)))[blk]).max(),
          np.abs((k1 @ e1 - q**2*e1 @ k1)[blk]).max(), np.abs((k1 @ e0 - q**-2*e0 @ k1)[blk]).max(),
          np.abs((k0 @ e0 - q**2*e0 @ k0)[blk]).max(), np.abs((e1 @ e0 - e0 @ e1)[blk]).max())
cen = np.abs(np.diag(e1 @ e0)[inner] - np.diag(e1 @ e0)[inner][0]).max()
worst_int = 0
for th in THETAS:
    T, V, W = sg_tmatrix(th, LAM, q, thd, N, s=s_, scalar=False)
    m_in = np.tile(np.abs(ns) <= N - 3, 1)
    rows = np.repeat(np.abs(ns) <= N - 3, 2)              # W (x) V output
    cols = np.tile(np.abs(ns) <= N - 3, 2)                # V (x) W input
    for j in (0, 1):
        jt = 1 - j
        for kd in 'ek':
            lhs = T @ (kron(V[(kd, j)], np.eye(len(ns))) + kron(V[('k', j)], rt[('e', j)]) if kd == 'e'
                       else kron(V[('k', j)], rt[('k', j)]))
            rhs = ((kron(rt[('e', j)], np.eye(2)) + kron(rt[('k', j)], V[('e', jt)])) if kd == 'e'
                   else kron(rt[('k', j)], V[('k', jt)])) @ T
            worst_int = max(worst_int, np.abs((lhs - rhs)[np.ix_(rows, cols)]).max()/np.abs(lhs[np.ix_(rows, cols)]).max())
# the same exactly, for symbolic n, q, s, y, y_d: operators on W as {shift m: coefficient(n)}, |n> -> c(n)|n+m>
nn, yds = sp.symbols('n y_d')


def s_mul(P, R):
    out = {}
    for m, c in R.items():
        for m2, c2 in P.items():
            out[m + m2] = out.get(m + m2, 0) + c*c2.subs(nn, nn + m)
    return out


def s_add(P, R, c=1):
    out = dict(P)
    for m, v in R.items():
        out[m] = out.get(m, 0) + c*v
    return out


def b_mul(P, R):
    return [[s_add(s_mul(P[i][0], R[0][j]), s_mul(P[i][1], R[1][j])) for j in range(2)] for i in range(2)]


def b_zero(P):
    return all(sp.simplify(v) == 0 for row in P for d in row for v in d.values())


def on_V(M, op=None):            # (2x2 matrix on V) (x) (operator on W, identity if None)
    op = {0: sp.Integer(1)} if op is None else op
    return [[{k: v*M[i][j] for k, v in op.items()} if M[i][j] != 0 else {} for j in range(2)] for i in range(2)]


aS, bS = sp.I*yds*qs**-ss_/(qs**2 - 1), sp.I*yds*qs**ss_/(qs**2 - 1)
rtS = {('e', 1): {1: aS*qs**(-2*nn)}, ('k', 1): {2: qs**(-2*nn - 2 - ss_)},
       ('e', 0): {-1: bS*qs**(2*nn)}, ('k', 0): {-2: qs**(2*nn - 2 + ss_)}}
LmS = [[{0: qs**-nn}, {-1: -sp.I*(ys/yds)*qs**(nn + ss_)}], [{1: -sp.I*(ys/yds)*qs**-nn}, {0: qs**(nn + ss_)}]]
piS = {('e', 1): [[0, ys], [0, 0]], ('e', 0): [[0, 0], [ys, 0]], ('k', 1): [[qs, 0], [0, 1/qs]], ('k', 0): [[1/qs, 0], [0, qs]]}
ex_int = True
for j in (1, 0):
    jt = 1 - j
    XV = [[s_add(a_, b_) for a_, b_ in zip(r1, r2)] for r1, r2 in zip(on_V(piS[('e', j)]), on_V(piS[('k', j)], rtS[('e', j)]))]
    XP = [[s_add(a_, b_) for a_, b_ in zip(r1, r2)] for r1, r2 in zip(on_V([[1, 0], [0, 1]], rtS[('e', j)]),
                                                                    on_V(piS[('e', jt)], rtS[('k', j)]))]
    lhs, rhs = b_mul(LmS, XV), b_mul(XP, LmS)
    ex_int &= b_zero([[s_add(lhs[i][k], rhs[i][k], -1) for k in range(2)] for i in range(2)])
    lhs, rhs = b_mul(LmS, on_V(piS[('k', j)], rtS[('k', j)])), b_mul(on_V(piS[('k', jt)], rtS[('k', j)]), LmS)
    ex_int &= b_zero([[s_add(lhs[i][k], rhs[i][k], -1) for k in range(2)] for i in range(2)])
r_ = lambda x_, y_: s_mul(rtS[x_], rtS[y_])
ex_rel = b_zero([[s_add(r_(('k', 0), ('k', 1)), {0: sp.Integer(1)}, -1), s_add(r_(('k', 1), ('e', 1)), r_(('e', 1), ('k', 1)), -qs**2)],
                 [s_add(r_(('k', 1), ('e', 0)), r_(('e', 0), ('k', 1)), -qs**-2), s_add(r_(('e', 1), ('e', 0)), r_(('e', 0), ('e', 1)), -1)]])
ex_rel &= sp.simplify(r_(('e', 1), ('e', 0))[0] - aS*bS*qs**2) == 0
report('second module rho~ of U_q(b^+) on W (type I_-): rho~(e_1)|n> = a q^-2n |n+1>, rho~(k_1)|n> = q^(-2n-2-s)|n+2>, '
       'rho~(e_0)|n> = b q^2n |n-1>, rho~(k_0)|n> = q^(2n-2+s)|n-2>, a = i y_d q^-s/(q^2-1), b = i y_d q^s/(q^2-1): '
       f'k_0 k_1 = 1, k e k^-1 = q^(+-2) e, e_1 e_0 = e_0 e_1 central [{max(rel, cen):.0e}]; and T intertwines '
       f'(pi (x) rho~)Delta(x) with (rho~ (x) pi o tau)Delta(x) for x = e_j, k_j (s = 0.3) [{worst_int:.0e}]; '
       f'exactly, for symbolic n, q, s, y, y_d: relations {ex_rel}, the four generator equations {ex_int}',
       max(rel, cen) < 1e-12 and worst_int < 1e-12 and ex_rel and ex_int)

rho_k1 = np.diag(q**(2.0*ns + s_))
qc = np.abs((rho_k1 @ k1 - q**4*k1 @ rho_k1)[blk]).max()
plain = np.abs((rho_k1 @ k1 - k1 @ rho_k1)[blk]).max()
W5 = T_I(0.1, 0.15, 7, -1)[1]
r_e1 = max(ls_symmetry(lambda th: T_I(th, 0.15, 7, -1)[0], W5,
                       [(gmat('e', j), None, 0), (gmat('k', j), 'Z', DE[j])], [(ONE, 'Z', DE[j])]) for j in (0, 1))
# type I_+: Delta(U_q(b^+)) and (1 (x) tau)Delta(U_q(b^-)); conjugating f_1 (x) k_0^-1 + 1 (x) f_0 by k_1 (x) k_1 isolates
# 1 (x) f_0, acting as 1 (x) Z on V (x) W and as pi(k_0^-1) (x) Z + pi(f_0) (x) 1 on V' (x) W
W5p = T_I(0.1, 0.15, 7, 1)[1]
r_f0 = ls_symmetry(lambda th: T_I(th, 0.15, 7, 1)[0], W5p, [(ONE, 'Z', DF[0])], [(gmat('k', 0, -1), 'Z', DF[0]),
                                                                                 (gmat('f', 0), None, 0)])
qn_ = complex(q)
Ek, Kk2 = np.array([[0, 1], [0, 0]]), np.diag([qn_, 1/qn_])
ya_, yb_ = np.exp(LAM*(0.3 + 0.1j)), np.exp(LAM*(-0.5 + 0.2j))
Gk = kron(Kk2, Kk2)
iso = lambda x_: (Gk @ x_ @ np.linalg.inv(Gk) - x_/qn_**2)/(qn_**2 - qn_**-2)
alg_ok = (np.allclose(iso(kron(ya_*Ek, np.eye(2)) + kron(Kk2, yb_*Ek.T)), kron(ya_*Ek, np.eye(2)))
          and np.allclose(iso(kron(Ek.T/ya_, Kk2) + kron(np.eye(2), Ek/yb_)), kron(np.eye(2), Ek/yb_)))
report(f'the two Borel symmetries do not combine: rho(k_1) rho~(k_1) = q^4 rho~(k_1) rho(k_1) '
       f'[{qc:.0e}; commutator {plain:.2f}], while k_1 (x) k_1 and k_1 (x) k_0 commute; e_1 (x) 1 = (G x G^-1 - q^-2 x)/'
       f'(q^2 - q^-2) with G = k_1 (x) k_1, x = e_1 (x) 1 + k_1 (x) e_0, and 1 (x) f_0 likewise from x = f_1 (x) k_0^-1 + '
       f'1 (x) f_0 (in pi (x) pi: {alg_ok}); neither is a symmetry: e_j (x) 1 for type I_- [residual {r_e1:.2f}], '
       f'1 (x) f_0 for type I_+ [{r_f0:.2f}]', qc < 1e-10 and plain > 0.1 and r_e1 > 0.5 and r_f0 > 0.5 and alg_ok)

# limits theta_d -> +-inf (prp-diagonal-borel (5)): residuals of the single failing generator of (3), N = 6
TW = {'e': 'id', 'f': 'tau'}                          # the failing generators of (3): I_- e/id, f/tau; I_+ f/id, e/tau
TW_P = {'e': 'tau', 'f': 'id'}
runs = [(sg_, k_, d_*t_) for sg_, k_, d_ in ((-1, 'e', 1), (-1, 'f', -1), (1, 'f', -1), (1, 'e', 1))
        for t_ in (1.5, 3.0, 4.5)]                     # towards the limit in which the generator's own equation is solvable
runs += [(-1, 'e', -4.5), (-1, 'f', 4.5), (1, 'f', 4.5), (1, 'e', -4.5)]   # opposite directions
lim = {}
for sign, kind, thd_ in runs:
    tw = (TW if sign == -1 else TW_P)[kind]
    W6 = T_I(0.1, thd_, 6, sign)[1]
    Tf = lambda th: T_I(th, thd_, 6, sign)[0]
    lim[(sign, kind, thd_)] = ls_symmetry(Tf, W6, *generator_terms(kind, 1, 1 if tw == 'id' else 0))
seq = {k: [lim[(k[0], k[1], t)] for t in (k[2]*1.5, k[2]*3.0, k[2]*4.5)]
       for k in ((-1, 'e', 1), (-1, 'f', -1), (1, 'f', -1), (1, 'e', 1))}
dec_ok = all(v[0] > v[1] > v[2] and v[2] < 0.02 and v[1]/v[2] > 4 for v in seq.values())
wrong = [lim[(-1, 'e', -4.5)], lim[(-1, 'f', 4.5)], lim[(1, 'f', 4.5)], lim[(1, 'e', -4.5)]]
nm_ = {(-1, 'e', 1): 'I_- e/id, theta_d > 0', (-1, 'f', -1): 'I_- f/tau, theta_d < 0', (1, 'f', -1): 'I_+ f/id, theta_d < 0',
       (1, 'e', 1): 'I_+ e/tau, theta_d > 0'}
report('limits (N = 6, |theta_d| = 1.5, 3, 4.5): single-generator residuals decay like e^{-lambda |theta_d|} towards the '
       'diagonal limit (twist lines) of I_- at +inf and of I_+ at -inf, and towards the off-diagonal limit (conjugation '
       'with a charge shift) in the other direction: ' + '; '.join(f'{nm_[k]}: ' + ', '.join(f'{x:.1e}' for x in v)
                                                                    for k, v in seq.items())
       + f'; in the opposite directions they stay of order 1 [{min(wrong):.2f}-{max(wrong):.2f} at |theta_d| = 4.5]. '
       'They measure the approach of the matrix, not a module: see the exact check above', dec_ok and min(wrong) > 0.5)
report('each charge sector C|n> of the diagonal limit is the twist line of prp-diagonal-coideal: diag(q^-n, q^(n+s)) '
       'conjugates pi(e_1), pi(e_0), pi(k_1) into t pi(e_1), pi(e_0)/t, pi(k_1), t = q^(-2n-s) (I_-), and diag(q^n, '
       'q^(-n-s)) with t = q^(2n+s) (I_+), so the sector twist depends on n',
       all(np.allclose(np.diag([a_, d_]) @ sg(1.3)[('e', 1)] @ np.diag([1/a_, 1/d_]), (a_/d_)*sg(1.3)[('e', 1)])
           and np.allclose(np.diag([a_, d_]) @ sg(1.3)[('e', 0)] @ np.diag([1/a_, 1/d_]), (d_/a_)*sg(1.3)[('e', 0)])
           and abs(a_/d_ - q**(-sgn*(2*n_ + 0.3))) < 1e-12
           for sgn in (1, -1) for n_ in (-2, 0, 3) for a_, d_ in [(q**(-sgn*n_), q**(sgn*(n_ + 0.3)))]))

# control: a particle as a defect has the full diagonal Delta(U) (W = V(theta_d), T = Rc, chi = pi_W)
Wp = sl2_spin(1, q, np.exp(LAM*(0.15 + 0.2j)))
worst_c = 0
for th in THETAS:
    y = np.exp(LAM*th)
    T = rcheck(sg(y), Wp)[0]
    Tm = rcheck(sg(1/y), Wp)[0]
    K = np.zeros((8, 8), complex)
    K[4:, :4] = flip(2, 2).T @ T
    K[:4, 4:] = np.linalg.inv(Tm) @ flip(2, 2)
    for kd in 'efk':
        for j in (0, 1):
            X = lambda yy: bdiag(cop(sg(yy), Wp, kd, j), flip(2, 2).T @ cop(Wp, sg(1/yy), kd, j) @ flip(2, 2))
            worst_c = max(worst_c, np.abs(K @ X(y) - X(1/y) @ K).max()/np.abs(K @ X(y)).max())
report(f'control: a particle frozen as a defect (W = V(theta_d), T = Rc) is invariant under the full diagonal Delta(U), '
       f'W acting by pi_W [{worst_c:.0e}]', worst_c < 1e-12)


# ================================================================ d. fusing the defect onto a Dirichlet boundary
def KD(y, eta):
    return np.diag([y*np.exp(1j*eta) + np.exp(-1j*eta)/y, y*np.exp(-1j*eta) + np.exp(1j*eta)/y])


def K_fused(theta, theta_d, eta, N, scalars=False):
    T, V, W = sg_tmatrix(theta, LAM, q, theta_d, N, s=0.0, scalar=scalars)
    Tm = sg_tmatrix(-theta, LAM, q, theta_d, N, s=0.0, scalar=scalars)[0]
    Kb = KD(np.exp(LAM*theta), eta)
    if scalars:
        u = -1j*theta
        Kb = complex(mp.exp(log_R0(u) + log_sigma(eta, u))/(2*mp.cos(eta)))*Kb
    return np.linalg.inv(Tm) @ kron(np.eye(W['dim']), Kb) @ T, W


def K_gz(y, eta, vt):
    """eq-bsg-k-matrix with eq-bsg-eta-theta."""
    e1, e0 = -np.cos(eta - 1j*vt)/np.sin(np.pi*LAM), -np.cos(eta + 1j*vt)/np.sin(np.pi*LAM)
    kap = q - 1/q
    return np.array([[kap*(e1*y + e0/y), y**2 - y**-2], [y**2 - y**-2, kap*(e0*y + e1/y)]])


def blocks(K, W, margin=2):
    ns = list(W['n'])
    N = max(ns)
    dW = W['dim']
    out = []
    for n in ns:
        if abs(n) <= N - margin and abs(n + 1) <= N - margin:
            ix = [ns.index(n), dW + ns.index(n + 1)]
            out.append(K[np.ix_(ix, ix)])
    return out


inv2 = lambda B: np.array([B[0, 0]/B[1, 1], B[0, 1]*B[1, 0]/(B[0, 0]*B[1, 1])])
eta_D, thd = 0.37, 0.21
dev_n, dev_gz, dev_alt = 0, 0, np.inf
for N in (5, 7, 9):
    for th in THETAS:
        Kf, W = K_fused(th, thd, eta_D, N)
        invs = np.array([inv2(B) for B in blocks(Kf, W)])
        dev_n = max(dev_n, np.abs(invs - invs[0]).max())
        y = np.exp(LAM*th)
        dev_gz = max(dev_gz, np.abs(inv2(K_gz(y, eta_D, LAM*thd)) - invs[0]).max())
        dev_alt = min(dev_alt, np.abs(inv2(K_gz(y, eta_D, -LAM*thd)) - invs[0]).max(),
                      np.abs(inv2(K_gz(y, eta_D, thd)) - invs[0]).max())
report('defect + Dirichlet: the charge-sector blocks of K_f = Tbar (1 (x) K_D) T have label-independent gauge invariants '
       f'(Corrigan-Zambon (2.18)) [{dev_n:.0e}] and equal Ghoshal-Zamolodchikov with (eta, vartheta) = (eta_D, lambda theta_d) '
       f'[{dev_gz:.0e}]; vartheta = -lambda theta_d or theta_d fail [{dev_alt:.2f}] (windows N = 5, 7, 9)',
       dev_n < 1e-12 and dev_gz < 1e-12 and dev_alt > 1e-3)


def re_res(Kf_th, t1, t2, W, margin=2):
    dW = W['dim']
    y1_, y2_ = np.exp(LAM*t1), np.exp(LAM*t2)
    Rc = lambda a, b: rcheck(sg(a), sg(b))[0]
    I2, IW = np.eye(2), np.eye(dW)
    K1, K2 = Kf_th(t1), Kf_th(t2)
    L = kron(Rc(1/y2_, 1/y1_), IW) @ kron(I2, K1) @ kron(Rc(y1_, 1/y2_), IW) @ kron(I2, K2)
    R = kron(I2, K2) @ kron(Rc(y2_, 1/y1_), IW) @ kron(I2, K1) @ kron(Rc(y1_, y2_), IW)
    m = np.abs(np.tile(W['n'], 4)) <= W['n'].max() - margin
    return prop_residual(L[np.ix_(m, m)], R[np.ix_(m, m)])[0]


re_w, re_c = [], []
for N in (5, 7, 9):
    W = K_fused(0.1, thd, eta_D, N)[1]
    re_w.append(re_res(lambda th: K_fused(th, thd, eta_D, N)[0], *PAIRS[0], W))
    pert = np.diag(1 + 0.05*np.cos(np.arange(2*W['dim'])))
    re_c.append(re_res(lambda th: K_fused(th, thd, eta_D, N)[0] @ pert, *PAIRS[0], W))
report('the infinite-dimensional K_f solves the reflection equation on V (x) V (x) W, interior labels, windows N = 5, 7, 9 '
       f'[{", ".join(f"{v:.0e}" for v in re_w)}]; a label-dependent perturbation fails [{min(re_c):.2f}]',
       max(re_w) < 1e-12 and min(re_c) > 1e-3)

# full scalar identity: Tbar K_D T = K_GZ(eta_D, lambda theta_d) including R_0, sigma (BPTT form, eq-bsg-k-full) and Phi
mp.mp.dps = 20
LL = mp.mpf(LAM)


def log_R0(u):
    u = mp.mpmathify(u)

    def t(l, uu):
        v = 2*LL*uu/mp.pi
        return (mp.loggamma(4*l*LL - v) + mp.loggamma(4*LL*(l - 1) + 1 - v)
                - mp.loggamma((4*l - 3)*LL - v) - mp.loggamma((4*l - 1)*LL + 1 - v))
    return mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])


def log_sigma(x, u):
    x, u = mp.mpmathify(x), mp.mpmathify(u)

    def t(l, uu):
        v, a = LL*uu/mp.pi, x/mp.pi
        return (mp.loggamma(0.5 + a + (2*l - 1)*LL - v) + mp.loggamma(0.5 - a + (2*l - 1)*LL - v)
                - mp.loggamma(0.5 - a + (2*l - 2)*LL - v) - mp.loggamma(0.5 + a + 2*l*LL - v))
    return mp.log(mp.cos(x)/mp.cos(x + LL*u)) + mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])


def K_full(theta, eta, vt):
    u = -1j*mp.mpmathify(theta)
    P0 = lambda sg_: mp.cos(LL*u)*mp.cos(eta)*mp.cosh(vt) - sg_*mp.sin(LL*u)*mp.sin(eta)*mp.sinh(vt)
    pref = mp.exp(log_R0(u) + log_sigma(eta, u) + log_sigma(1j*vt, u))/(mp.cos(eta)*mp.cosh(vt))
    return complex(pref)*np.array([complex(P0(1)), complex(P0(-1))])


dev_sc = 0
for th in (0.33, 0.21 + 0.17j):
    Kf, W = K_fused(th, thd, eta_D, 5, scalars=True)
    Kg = K_full(th, eta_D, LAM*thd)
    for B in blocks(Kf, W):
        dev_sc = max(dev_sc, abs(B[0, 0]/Kg[0] - 1), abs(B[1, 1]/Kg[1] - 1))
dev_inf = abs(complex(mp.exp(log_sigma(24j, -0.33j))) - 1)
report('with the scalar factors (Phi of eq-sgd-scalar for T and Tbar; Dirichlet K_D = R_0 sigma(eta,u)/(2 cos eta) '
       'diag(y e^{i eta} + e^{-i eta}/y, ...), the vartheta -> inf limit of eq-bsg-k-full since sigma(i vartheta,u) -> 1 '
       f'[{dev_inf:.0e} at vartheta = 24]) the diagonal entries of Tbar K_D T equal those of eq-bsg-k-full at '
       f'(eta_D, lambda theta_d) [{dev_sc:.0e}; the accuracy of the Richardson-summed Phi at complex rapidity is ~1e-8]',
       dev_sc < 1e-7 and dev_inf < 1e-6)

# sinh-Gordon: eq-shg-reflection, blocks (x) = sinh(th/2 + i pi x/4)/sinh(th/2 - i pi x/4)
Bc = 0.43
blk_ = lambda x, th: np.sinh(th/2 + 1j*np.pi*x/4)/np.sinh(th/2 - 1j*np.pi*x/4)
Kshg = lambda th, E, F: (blk_(1, th)*blk_(2 - Bc/2, th)*blk_(1 + Bc/2, th)
                         / (blk_(1 - E, th)*blk_(1 + E, th)*blk_(1 - F, th)*blk_(1 + F, th)))
KDsh = lambda th, E: blk_(1, th)*blk_(2 - Bc/2, th)*blk_(1 + Bc/2, th)/(blk_(1 - E, th)*blk_(1 + E, th))
Tcl = lambda th: -1j*np.sinh(th/2 - 1j*np.pi/4)/np.sinh(th/2 + 1j*np.pi/4)
dev_sh = max(abs(Tcl(1j*np.pi - th - thd)*KDsh(th, 0.31)*Tcl(th - thd)/Kshg(th, 0.31, 2j*thd/np.pi) - 1)
             for th in (0.3, 1.1 + 0.2j, -0.6))
report(f'sinh-Gordon (Bajnok-Simon): Tbar_shG K_D(E) T_shG = K_shG(E, F = 2 i theta_d/pi) of eq-shg-reflection [{dev_sh:.0e}]',
       dev_sh < 1e-12)

# Corrigan-Zambon 1202.6016 in their conventions: T^{b beta}_{a alpha}, rows = incoming (a, alpha), a + alpha = b + beta;
# (2.9) is the row-incoming product T(theta) R(theta) T(-theta)^-1; x = e^{gamma theta}; Q^alpha = q^{-alpha/2}.
# On the sector {(+, alpha), (-, alpha+2)}: T = [[a(alpha), b(alpha)], [c(alpha+2), d(alpha+2)]] of (2.14), and the fused
# matrix is [[r_+(alpha), s_+(alpha)], [s_-(alpha+2), r_-(alpha+2)]] up to a factor common to all sectors.
Qc, xc, rc_ = np.exp(-0.355j), np.exp(0.41), np.exp(0.23j)
qc_ = Qc**-2
Xc = xc**2 - xc**-2
RDc = np.diag([rc_*xc + 1/(rc_*xc), rc_/xc + xc/rc_])                     # (2.16)


def cz_T(xx, A, p):
    ap, am, bp, bm, cp, cm, dp, dm = p
    A2 = A*Qc**2                                                       # Q^(alpha+2)
    return np.array([[ap*A/xx + am*xx/A, bp*A + bm/A], [cp*A2 + cm/A2, dp*A2*xx + dm/(A2*xx)]])


def cz_217(A, p):
    ap, am, bp, bm, cp, cm, dp, dm = p
    qa = A**-2                                                         # q^alpha
    g, h = am*dp/qc_ + ap*dm*qc_, bm*cp/qc_ + bp*cm*qc_
    rp = Xc*(am*dm*qa*qc_*rc_*xc - ap*dp/(qa*qc_*rc_*xc)) + xc*(g*rc_ - h/rc_) + (g/rc_ - h*rc_)/xc
    rm = Xc*(ap*dp*qc_/qa*xc/rc_ - am*dm*qa/qc_*rc_/xc) + xc*(g/rc_ - h*rc_) + (g*rc_ - h/rc_)/xc
    sp_ = Xc*(ap*bm/rc_ - am*bp*rc_ - am*bm*qa*rc_ + ap*bp/(qa*rc_))
    sm_ = Xc*(dm*cp*rc_ - dp*cm/rc_ + dm*cm*qa*rc_ - dp*cp/(qa*rc_))
    return rp, rm, sp_, sm_


def cz_219(A, b, fix_h=False, fix_sign=False):
    bp, bm = b
    Bp, Bm = np.conj(bp), np.conj(bm)
    g = (abs(bp)**2*abs(bm)**2 + 1)*Qc**-2
    h = -(abs(bm)**2 + abs(bp)**2)*Qc**-2 if fix_h else -(abs(bm)**2*Qc**-2 + abs(bp)**2)*Qc**-2
    rp = Xc*(bp*Bm*A**2*Qc**-2/(rc_*xc) - Bp*bm*A**-2*Qc**-2*rc_*xc) + xc*(g*rc_ - h/rc_) + (g/rc_ - h*rc_)/xc
    rm = Xc*(-bp*Bm*A**2*Qc**-6*xc/rc_ + Bp*bm*A**-2*Qc**2*rc_/xc) + xc*(g/rc_ - h*rc_) + (g*rc_ - h/rc_)/xc
    sp_ = Xc*(bm/rc_ + abs(bp)**2*bm*rc_ + bp*A**2/rc_ + bm**2*Bp*A**-2*rc_)
    sm_ = (-1 if fix_sign else 1)*Xc*(Bm*Qc**-4*rc_ + abs(bp)**2*Bm*Qc**-4/rc_ + Bm**2*bp*A**2*Qc**-8/rc_ + Bp*A**-2*rc_)
    return rp, rm, sp_, sm_


def cz_dev(p, formula):
    """max relative spread of fused/formula over all entries of the sectors alpha = 0, 2, 4"""
    ratios = []
    for A in (1, Qc**2, Qc**4):
        F = cz_T(xc, A, p) @ RDc @ np.linalg.inv(cz_T(1/xc, A, p))
        e0, e2 = formula(A), formula(A*Qc**2)
        M = np.array([[e0[0], e0[2]], [e2[3], e2[1]]])
        ratios += list((F/M).ravel())
    ratios = np.array(ratios)
    return np.abs(ratios - ratios[0]).max()/abs(ratios[0])


rng_cz = np.random.default_rng(7)
cpx = lambda: complex(rng_cz.normal(), rng_cz.normal())
ap_, am_, bp_, bm_, cp_, cm_ = (cpx() for _ in range(6))
gen_p = (ap_, am_, bp_, bm_, cp_, cm_, bp_*cp_/ap_, bm_*cm_/am_)                  # a_+- d_+- = b_+- c_+-
dev217 = cz_dev(gen_p, lambda A: cz_217(A, gen_p))
bII = (0.37 + 0.19j, -0.28 + 0.13j)
Bp_, Bm_ = np.conj(bII[0]), np.conj(bII[1])
pII = (1, -Bp_*bII[1], bII[0], bII[1], -Bm_*qc_**2, -Bp_, -bII[0]*Bm_*qc_**2, 1)  # type II, CZ p. 5
devII = {(fh, fs): cz_dev(pII, lambda A, fh=fh, fs=fs: cz_219(A, bII, fh, fs)) for fh in (False, True) for fs in (False, True)}
dev217II = cz_dev(pII, lambda A: cz_217(A, pII))
FII = [cz_T(xc, A, pII) @ RDc @ np.linalg.inv(cz_T(1/xc, A, pII)) for A in (1, Qc**2, Qc**4)]
ratII = [B[0, 0]/B[1, 1] for B in FII]
report('Corrigan-Zambon 1202.6016: direct fusion (2.9) of their (2.14) (generic, a d = b c) with (2.16) equals (2.17) '
       f'up to a factor common to all sectors [{dev217:.0e}], also at the type-II parameters [{dev217II:.0e}]; (2.19) as '
       f'printed fails [{devII[(False, False)]:.2f}], with only h = -(|b_-|^2 + |b_+|^2) Q^-2 [{devII[(True, False)]:.2f}] or '
       f'only the sign of s_- reversed [{devII[(False, True)]:.2f}] it still fails, with both it holds '
       f'[{devII[(True, True)]:.0e}]; the type-II diagonal ratio depends on the charge (alpha = 0, 2, 4: '
       + ', '.join(f'{v:.3f}' for v in ratII) + ')',
       dev217 < 1e-12 and dev217II < 1e-12 and devII[(True, True)] < 1e-12
       and min(devII[(False, False)], devII[(True, False)], devII[(False, True)]) > 1e-2
       and abs(ratII[0] - ratII[1]) > 1e-2)


# ================================================================ e. classical: the K-matrix dressed by L_d; BDG
beta, m = 0.7, 1.0


def Ld(lam, phi, psi, sigma, r):
    """type I_- defect matrix (eq-defect-lax-matrix) in the defining representation, real coupling."""
    U = an_data(r)['U']
    D, E = dmatrix(r, 1)
    Lt = np.eye(r + 1) + (sigma/lam)*sum(Ej.T for Ej in an_lax_rep(r))
    return expm(-(beta/2)*np.diag(U @ (E @ phi + D @ psi))) @ Lt @ expm((beta/2)*np.diag(U @ (D.T @ phi - E @ psi)))


def K_bdry(lam, C, r):
    """field-independent K solving eq-boundary-k-equation."""
    d = an_data(r)
    U, Es, n = d['U'], an_lax_rep(r), r + 1
    rows = []
    for Cj, a, Ej in zip(C, d['alpha'], Es):
        h = np.diag(U @ a)
        rows.append(np.kron((-Cj*h + lam*Ej - Ej.T/lam).T, np.eye(n)) - np.kron(np.eye(n), Cj*h + lam*Ej - Ej.T/lam))
    return sp_null(np.vstack(rows))[:, 0].reshape(n, n).T


def dressed_compat(r, C, sigma, lam, dress):
    phi, psi, psit = rng.normal(size=r)*0.6, rng.normal(size=r)*0.6, rng.normal(size=r)
    _, gphi, gpsi, _ = defect_potential(phi, psi, sigma, r, beta, m, 'real', 1)
    gB = (m/beta)*sum(Cj*a*np.exp(beta*(a @ psi)/2) for Cj, a in zip(C, an_data(r)['alpha']))
    D, E = dmatrix(r, 1)
    phit = np.linalg.solve(D.T, -gB + E @ psit - gpsi)      # psi_x(0) from the defect = -grad B_psi (fused limit)
    phix = E @ phit + D @ psit - gphi
    Kp = K_bdry(lam, C, r)
    if dress == 'inv':
        Kf = lambda ph, ps: np.linalg.inv(Ld(1/lam, ph, ps, sigma, r)).T @ Kp @ np.linalg.inv(Ld(lam, ph, ps, sigma, r))
    elif dress == 'direct':
        Kf = lambda ph, ps: Ld(1/lam, ph, ps, sigma, r).T @ Kp @ Ld(lam, ph, ps, sigma, r)
    else:
        Kf = lambda ph, ps: Kp
    h = 1e-6
    dK = (Kf(phi + h*phit, psi + h*psit) - Kf(phi - h*phit, psi - h*psit))/(2*h)
    K0 = Kf(phi, psi)
    rhs = K0 @ a_t(lam, phi, phix, r, beta, m) + a_t(1/lam, phi, phix, r, beta, m).T @ K0
    return np.abs(dK - rhs).max()/np.abs(rhs).max()


cl_cases = [(1, [0.7, -1.3]), (1, [0, 0]), (2, [1, -1, 1]), (2, [0, 0, 0]), (3, [1, -1, 1, -1])]
cl_ok = max(dressed_compat(r, C, 0.8, 0.9 + 0.3j, 'inv') for r, C in cl_cases)
cl_bad = min(min(dressed_compat(r, C, 0.8, 0.9 + 0.3j, 'direct'), dressed_compat(r, C, 0.8, 0.9 + 0.3j, 'bare'))
             for r, C in cl_cases if any(C))
report('classical fusion: K_phi(lam) = [L_d(1/lam)^-1]^T K_psi(lam) L_d(lam)^-1 satisfies d_t K = K a_t(lam) + a_t(1/lam)^T K '
       'on the fused defect + boundary system (a_1 with GZ C = (0.7, -1.3) and Neumann; a_2 with C = (1,-1,1) and Neumann; '
       f'a_3 with C = (1,-1,1,-1)) [{cl_ok:.0e}, finite differences]; L_d instead of L_d^-1, or K_psi undressed, fails '
       f'[{cl_bad:.2f}]', cl_ok < 1e-8 and cl_bad > 0.1)

# BDG: sine-Gordon defect potential (type I_-) = H_osc at complex oscillator values
f_, g_, bs, m0, td = sp.symbols('phi psi beta_SG m_0 theta_d')
sig_ = sp.exp(-td)
Bd = -(2*m0/bs**2)*(sig_*sp.cos(bs*(f_ + g_)/2) + sp.cos(bs*(f_ - g_)/2)/sig_)
ph_, qh_ = bs*g_/2 - sp.I*td, bs*g_/2 + sp.I*td
Hosc = -(2*m0/bs**2)*(sp.exp(-sp.I*bs*f_/2)*sp.cos(ph_) + sp.exp(sp.I*bs*f_/2)*sp.cos(qh_))
id_ok = sp.simplify((Bd - Hosc).rewrite(sp.exp)) == 0
C1 = sp.cos(qh_).subs(g_, sp.Symbol('phi_D'))
C1_ok = sp.simplify((C1 - sp.cosh(td - sp.I*bs*sp.Symbol('phi_D')/2)).rewrite(sp.exp)) == 0
# BDG flow (dot p = -(b^2/4) dH/dq, dot q = (b^2/4) dH/dp): d(p - q)/dt on the surface p - q = -2 i theta_d
P_, Q_ = sp.symbols('P Q')
Hgen = -(2*m0/bs**2)*(sp.exp(-sp.I*bs*f_/2)*sp.cos(P_) + sp.exp(sp.I*bs*f_/2)*sp.cos(Q_))
dpq = (-(bs**2/4)*sp.diff(Hgen, Q_) - (bs**2/4)*sp.diff(Hgen, P_)).subs({P_: ph_, Q_: qh_})
val = complex(dpq.subs({f_: 0.3, g_: 0.7, bs: 1.1, m0: 1.0, td: 0.4}))
report('BDG: the sine-Gordon defect potential (type I_-) equals H_osc of eq-bsg-oscillator at (p^, q^) = '
       '(beta_SG psi/2 - i theta_d, beta_SG psi/2 + i theta_d); at psi = phi_D, C_1 = cos q^ = cosh(theta_d - i beta_SG phi_D/2); '
       f'the BDG flow does not preserve p^ - q^ = -2 i theta_d [d(p^-q^)/dt = {abs(val):.2f} at a sample point]',
       id_ok and C1_ok and abs(val) > 0.1)


# ================================================================ f. tensor K-matrix: a particle fused onto a coideal boundary
eps = (0.37 - 0.21j, -0.52 + 0.4j)
cq = 1/q
Wd = sl2_spin(1, q, np.exp(LAM*(0.2 + 0.15j)))


def Kbare(y):
    null, _ = intertwiners([(coideal_gen(sg(y), j, cq, eps[j]), coideal_gen(sg(1/y), j, cq, eps[j])) for j in (0, 1)])
    return null[0]


def Ktensor(th):
    y = np.exp(LAM*th)
    return np.linalg.inv(rcheck(sg(1/y), Wd)[0]) @ kron(np.eye(2), Kbare(y)) @ rcheck(sg(y), Wd)[0]


worst_f, dims_f = 0, []
for th in THETAS:
    y = np.exp(LAM*th)
    Kt = Ktensor(th)
    pairs = [(coideal_gen_tensor(sg(y), Wd, j, cq, eps[j]), coideal_gen_tensor(sg(1/y), Wd, j, cq, eps[j])) for j in (0, 1)]
    worst_f = max(worst_f, max(np.abs(Kt @ A - B @ Kt).max() for A, B in pairs)/np.abs(Kt).max())
    dims_f.append(len(intertwiners(pairs, tol=1e-8)[0]))
t1, t2 = PAIRS[0]
Y1, Y2 = np.exp(LAM*t1), np.exp(LAM*t2)
Rc = lambda a, b: rcheck(sg(a), sg(b))[0]
I2, I2w = np.eye(2), np.eye(2)
Lre = kron(Rc(1/Y2, 1/Y1), I2w) @ kron(I2, Ktensor(t1)) @ kron(Rc(Y1, 1/Y2), I2w) @ kron(I2, Ktensor(t2))
Rre = kron(I2, Ktensor(t2)) @ kron(Rc(Y2, 1/Y1), I2w) @ kron(I2, Ktensor(t1)) @ kron(Rc(Y1, Y2), I2w)
re_t = prop_residual(Lre, Rre)[0]
report('tensor K-matrix: for the sine-Gordon coideal (generic eps-hat) and W = V(theta_d), Tbar (1 (x) K) T with T = Rc '
       f'intertwines Delta(b_j) on V(theta) (x) W [{worst_f:.0e}], the intertwiner space is {set(dims_f)}-dimensional '
       f'(Appel-Vlaar Rem. 6.5.3(b)), and the reflection equation holds [{re_t:.0e}]',
       worst_f < 1e-10 and set(dims_f) == {1} and re_t < 1e-10)
