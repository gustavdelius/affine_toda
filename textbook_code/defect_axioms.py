"""Axioms of defect scattering (Part VIII, ch. 36: sec-defect-kinematics, sec-rt-algebra, sec-triangle-relations,
sec-defect-unitarity-crossing, sec-doubling-trick, sec-defect-bound-states, sec-defect-circle, sec-topological-lines,
sec-defect-lee-yang).

Conventions (book, ch. 36): the defect operator D sits to the right of the particles on its left;
  A_a(theta) D_alpha = R_{a alpha}^{b beta}(theta) A_b(-theta) D_beta + TT_{a alpha}^{beta b}(theta) D_beta A_b(theta),
  D_alpha A_a(-theta) = Rbar(theta) D A(theta) + TTbar_{alpha a}^{b beta}(theta) A_b(-theta) D_beta      (theta > 0),
TT(theta): V (x) W -> W (x) V (left to right), TTbar(theta): W (x) V -> V (x) W (right to left), tensor order = spatial
order, TT(v_a (x) w_alpha) = sum TT_{a alpha}^{beta b} w_beta (x) v_b.  Bulk: A_a(t1) A_b(t2) = S_ab^cd(t12) A_d(t2) A_c(t1),
braid form Scheck(v_a (x) v_b) = sum S_ab^cd v_d (x) v_c, R = P Rcheck on fixed slots (book, sec-trigonometric-r).
Sine-Gordon R-check from boundary_common (Part VI coproduct, q = q_phys, principal gradation y = exp(lambda theta)).

Checks:
 a. Triangle relations (eq-triangle-st, eq-triangle-braid, eq-triangle-rtt, prp-triangle-forms): for a random
    charge-conserving TT on W = spin 1 the residuals of the index (S/TT) form, the braid form
    Rc_23 TT_12(t1) TT_23(t2) = TT_12(t2) TT_23(t1) Rc_12 and the RTT form R_12 L_13(t1) L_23(t2) = L_23(t2) L_13(t1) R_12
    (L = P TT) coincide entry by entry after relabelling; TT = Rc_{V,W}(y/y_d) (a spin-1 'particle' as defect) solves
    all three, a perturbation of it solves none.
 b. RT algebra (eq-defect-zf, prp-rt-nogo): a symbolic ZF engine reduces A_a(t1) A_b(t2) D in the two orders; for a
    diagonal S the differences are R_b(t2) TT_a(t1) [S_ab(t1+t2) - S_ab(t1-t2)] and R_a(t1) TT_b(t2) [1 - S_ab(t1-t2)
    S_ba(t1+t2)]; the both-transmitted coefficients agree identically, the both-reflected ones if S_ab = S_ba; with
    a = b and R, TT != 0 these force S^2 = 1 (Groebner basis).  sin chi = -g/(1 + g^2/4), cos chi = (1 - g^2/4)/(1 + g^2/4).  The Delfino-Mussardo-Simonetti Ising energy line (hep-th/9409076 (3.5)) satisfies
    unitarity (2.8) = eq-defect-unitarity with reflection, crossing (2.9), the RT relations with S = -1, the duality
    (3.12) and T = 0 at g = +-2.
 c. Lee-Yang defect (eq-ly-transmission, prp-ly-transmission): TT_b = [b+1][b-1] (Bajnok-Simon 0712.4292 (22)) =
    f_{2/3}(theta - i pi(3-b)/6), TTbar_b = [5-b][-5-b] = f_{2/3}(theta + i pi(3-b)/6); unitarity, crossing, bootstrap;
    fusion TTbar K_(1) TT = K_b of eq-ly-reflection; pole of TTbar at i pi(b+1)/6, excited defect TT^{|1>} =
    S(theta + iu) TT_b = TT_{b+2}, TTbar^{|1>} = S(theta - iu) TTbar_b; energies m sin(b pi/6), m sin((b+2)pi/6);
    pure phase for b = 3 + i a and b = -3 + i a; a CDD-type modification violating crossing is caught.  The
    sinh-Gordon defect factor (Bajnok-Simon (16)) satisfies unitarity, crossing, TT(theta) TT(theta + i pi) = 1 (the
    defect CDD condition, exr-defect-cdd) and TT^2 = 1/f_(1/2)(theta + B vartheta).
 d. Doubling trick (eq-doubled-k, prp-doubling): doubled K = (0, TTbar P; P TT, 0) and doubled bulk
    (Rc on the first copy, P Rc P on the reflected copy, flip between copies) satisfy the braid reflection equation
    eq-reflection-braid exactly when TT satisfies the triangle relations; the phi phi -> psi' psi' block of the
    reflection-equation residual equals the triangle residual entry by entry; K(theta) K(-theta) = 1 is defect unitarity;
    for Lee-Yang the doubled crossing-unitarity eq-boundary-crossing-unitarity holds and is equivalent to defect crossing.
 e. Twist defect (prp-twist-defect): dim W = 1, TT = exp(i alpha.mu) on a_2 vector solitons satisfies the triangle
    relation; the defect Bethe-Yang transfer matrices with the defect anywhere on the circle commute and have, sector by
    sector, the spectrum of the Bloch-twisted matrices exp(i alpha.mu_1) T_1 of prp-bloch-sectors (N = 3); with the
    defect at the seam T_1^D = exp(i alpha.mu_1) T_1 exactly; applying exp(i alpha.mu_j) after each cyclic T_j (a
    different seam for each j) gives the right spectra but non-commuting matrices.
 f. Modified crossing (sec-topological-lines, Copetti-Cordova-Komatsu 2403.04835 (18), (21), (22)): for the
    phi_{1,3} kinks of M_5 the literature S-matrix obeys standard crossing, the proposal (21) the modified crossing (22),
    and they differ by (d_a d_c/(d_b d_d))^{i theta/2 pi}, a rapidity-dependent normalization of the kink states.
Runtime: a few seconds.
"""
import itertools
import numpy as np
import mpmath as mp
import sympy as sp
from boundary_common import report, q_phys, sl2_spin, an_vector, rcheck

mp.mp.dps = 30
rng = np.random.default_rng(7)
kron = np.kron


def op_on(M, slots, dims):
    """Matrix of the operator M, acting on the tensor product of the factors `slots` (in that order) and mapping it to
    the same factors, on the full tensor product with factor dimensions `dims` (numpy kron ordering)."""
    n = len(dims)
    rest = [i for i in range(n) if i not in slots]
    order = list(slots) + rest
    D = int(np.prod(dims))
    dord = [dims[i] for i in order]
    Pm = np.zeros((D, D))
    for idx in itertools.product(*[range(d) for d in dims]):
        Pm[np.ravel_multi_index([idx[i] for i in order], dord), np.ravel_multi_index(idx, dims)] = 1
    big = kron(M, np.eye(int(np.prod([dims[i] for i in rest])) if rest else 1))
    return Pm.T @ big @ Pm


def flipmat(d1, d2):
    """P: C^d1 (x) C^d2 -> C^d2 (x) C^d1."""
    P = np.zeros((d1*d2, d1*d2))
    for i in range(d1):
        for j in range(d2):
            P[j*d1 + i, i*d2 + j] = 1
    return P


# ================================================================ a. three forms of the triangle relations
LAM = 1.37
q = q_phys(LAM)
V = lambda y: sl2_spin(1, q, y)
W = lambda y: sl2_spin(2, q, y)
Rc = lambda y1, y2: rcheck(V(y1), V(y2))[0]                       # V(y1) (x) V(y2) -> V(y2) (x) V(y1)
yd = np.exp(LAM*(0.21 + 0.13j))                                    # defect rapidity theta_d
TT_sol = lambda y: rcheck(V(y), W(yd))[0]                          # V(y) (x) W -> W (x) V(y)
dV, dW = 2, 3
y1, y2 = np.exp(LAM*(0.64 + 0.05j)), np.exp(LAM*(-0.33 + 0.11j))
wV = np.array([1, -1]); wW = np.array([2, 0, -2])                  # weights (2 x spin component)
mask = np.array([[1.0 if wW[be] + wV[b] == wV[a] + wW[al] else 0.0
                  for a in range(dV) for al in range(dW)] for be in range(dW) for b in range(dV)])
Nrand = lambda: mask*(rng.normal(size=mask.shape) + 1j*rng.normal(size=mask.shape))


def three_forms(TTf):
    R = Rc(y1, y2)
    T1, T2 = TTf(y1), TTf(y2)
    # braid form, V (x) V (x) W -> W (x) V (x) V
    IW, IV = np.eye(dW), np.eye(dV)
    Lb = kron(IW, R) @ kron(T1, IV) @ kron(IV, T2)
    Rb = kron(T2, IV) @ kron(IV, T1) @ kron(R, IW)
    res_b = (Lb - Rb).reshape(dW, dV, dV, dV, dV, dW)              # [g, e, f, a, b, al]
    # index form: S_ab^cd = Rc[(d,c),(a,b)];  TT_{a al}^{be b} = T[(be,b),(a,al)]
    S = R.reshape(dV, dV, dV, dV).transpose(2, 3, 1, 0)           # S[a,b,c,d]
    t1 = T1.reshape(dW, dV, dV, dW).transpose(2, 3, 0, 1)          # t[a,al,be,b]
    t2 = T2.reshape(dW, dV, dV, dW).transpose(2, 3, 0, 1)
    lhs = np.einsum('abcd,cABf,dBGe->abAGef', S, t1, t2)           # S(t12) TT(t1) TT(t2)   (path B)
    rhs = np.einsum('bABd,aBGc,cdfe->abAGef', t2, t1, S)           # TT(t2) TT(t1) S(t12)   (path A)
    res_i = (rhs - lhs).transpose(3, 4, 5, 0, 1, 2)                # [g, e, f, a, b, al]
    # RTT form on V1 (x) V2 (x) W, R = P Rc, L = P TT : V (x) W -> V (x) W
    Rfix = flipmat(dV, dV) @ R
    L1, L2 = flipmat(dW, dV) @ T1, flipmat(dW, dV) @ T2
    dims = [dV, dV, dW]
    R12 = op_on(Rfix, [0, 1], dims)
    Lr = R12 @ op_on(L1, [0, 2], dims) @ op_on(L2, [1, 2], dims)  # R12 L13(t1) L23(t2)   (path A)
    Rr = op_on(L2, [1, 2], dims) @ op_on(L1, [0, 2], dims) @ R12   # L23(t2) L13(t1) R12   (path B)
    res_r = (Lr - Rr).reshape(dV, dV, dW, dV, dV, dW)              # [f, e, g, a, b, al]
    res_r = res_r.transpose(2, 1, 0, 3, 4, 5)                      # [g, e, f, a, b, al]
    nrm = np.abs(res_b).max()
    return nrm, np.abs(res_b - res_i).max(), np.abs(res_b - res_r).max(), Lb


nrm0, di0, dr0, Lb0 = three_forms(TT_sol)
report(f'triangle relation, TT = Rc_(V,W)(y/y_d) with W = spin 1: braid, index (S/TT) and RTT forms all hold '
       f'[{nrm0/np.abs(Lb0).max():.0e}]', nrm0/np.abs(Lb0).max() < 1e-12 and di0 < 1e-12 and dr0 < 1e-12)
Nfix = Nrand()
nrm1, di1, dr1, _ = three_forms(lambda y: TT_sol(y) + 0.05*Nfix)
report(f'perturbed TT (+ 0.05 random charge-conserving): all three forms fail [{nrm1:.2f}], with residuals equal entry '
       f'by entry [{max(di1, dr1):.0e}]', nrm1 > 1e-3 and max(di1, dr1) < 1e-12)
Ra, Rb_ = Nrand(), Nrand()
nrm2, di2, dr2, _ = three_forms(lambda y: Ra + y*Rb_)
report(f'random charge-conserving TT(y) = A + y B: residuals of the three forms coincide entry by entry '
       f'[{max(di2, dr2):.0e}] (size {nrm2:.1f})', nrm2 > 1e-3 and max(di2, dr2) < 1e-12)

# ================================================================ b. RT algebra and the no-go theorem
t1s, t2s = sp.symbols('theta1 theta2', positive=True)
Sf, Rf, Tf, Rbf, Tbf = (sp.Function(n) for n in ('S', 'R', 'TT', 'Rbar', 'TTbar'))
num = {t1s: 2, t2s: 1}                                             # theta1 > theta2 > 0


def val(x):
    return float(x.subs(num))


def step(word, coeff, where):
    """Apply one exchange rule at the leftmost ('L') or rightmost ('R') applicable position.
    Tokens: ('A', species, rapidity) or ('D',). Returns list of (word, coeff) or None if in normal form."""
    pos = range(len(word) - 1) if where == 'L' else reversed(range(len(word) - 1))
    for i in pos:
        u, v = word[i], word[i + 1]
        if u[0] == 'A' and v[0] == 'A' and val(u[2]) > val(v[2]):
            return [(word[:i] + (v, u) + word[i + 2:], coeff*Sf(u[1], v[1], u[2] - v[2]))]
        if u[0] == 'A' and v[0] == 'D' and val(u[2]) > 0:
            a, x = u[1], u[2]
            return [(word[:i] + (('A', a, -x), v) + word[i + 2:], coeff*Rf(a, x)),
                    (word[:i] + (v, u) + word[i + 2:], coeff*Tf(a, x))]
        if u[0] == 'D' and v[0] == 'A' and val(v[2]) < 0:
            a, x = v[1], v[2]
            return [(word[:i] + (u, ('A', a, -x)) + word[i + 2:], coeff*Rbf(a, -x)),
                    (word[:i] + (v, u) + word[i + 2:], coeff*Tbf(a, -x))]
    return None


def normal_form(word, where):
    todo, out = [(word, sp.Integer(1))], {}
    while todo:
        w, c = todo.pop()
        nxt = step(w, c, where)
        if nxt is None:
            out[w] = out.get(w, 0) + c
        else:
            todo += nxt
    return out


a_s, b_s = sp.symbols('a b')
start = (('A', a_s, t1s), ('A', b_s, t2s), ('D',))
NA, NB = normal_form(start, 'R'), normal_form(start, 'L')
diffs = {w: sp.factor(NA.get(w, 0) - NB.get(w, 0)) for w in set(NA) | set(NB)}
key = lambda *toks: tuple(toks)
w_rr = key(('A', a_s, -t1s), ('A', b_s, -t2s), ('D',))
w_tt = key(('D',), ('A', b_s, t2s), ('A', a_s, t1s))
w_bt = key(('A', b_s, -t2s), ('D',), ('A', a_s, t1s))
w_at = key(('A', a_s, -t1s), ('D',), ('A', b_s, t2s))
S_ = lambda x, y, z: Sf(x, y, z)
e1 = Rf(b_s, t2s)*Tf(a_s, t1s)*(S_(a_s, b_s, t1s + t2s) - S_(a_s, b_s, t1s - t2s))
e2 = Rf(a_s, t1s)*Tf(b_s, t2s)*(1 - S_(a_s, b_s, t1s - t2s)*S_(b_s, a_s, t1s + t2s))
e0 = Rf(a_s, t1s)*Rf(b_s, t2s)*(S_(a_s, b_s, t1s + t2s)*S_(b_s, a_s, t1s - t2s)
                                 - S_(a_s, b_s, t1s - t2s)*S_(b_s, a_s, t1s + t2s))
same = lambda x, y: sp.simplify(x - y) == 0 or sp.simplify(x + y) == 0
ok_struct = (len(diffs) == 4 and diffs[w_tt] == 0 and same(diffs[w_rr], e0) and same(diffs[w_bt], e1)
             and same(diffs[w_at], e2) and sp.simplify(diffs[w_rr].subs(b_s, a_s)) == 0)
report('ZF associativity for A_a(t1) A_b(t2) D (diagonal S): both-transmitted coefficients agree identically; '
       'both-reflected ones differ by R_a R_b [S_ab(t1+t2) S_ba(t12) - S_ab(t12) S_ba(t1+t2)], zero if S_ab = S_ba '
       '(e.g. a = b); the mixed ones by R_b TT_a [S_ab(t1+t2) - S_ab(t12)] and R_a TT_b [1 - S_ab(t12) S_ba(t1+t2)]',
       ok_struct)
su, sv = sp.symbols('s_u s_v')
G = sp.groebner([sv - su, su*sv - 1], sv, su, order='lex')
report(f'with a = b and R, TT not identically zero: S(u) = S(v), S(u) S(v) = 1 for all 0 < u < v give the Groebner '
       f'basis {list(G.exprs)}, so S^2 = 1', su**2 - 1 in [sp.expand(g) for g in G.exprs])
# parity-free: only the left rules were used
used = set(f.func.__name__ for N_ in (NA, NB) for c in N_.values() for f in c.atoms(sp.Function))
report(f'the derivation used only the left exchange relation: functions appearing {sorted(used)} (no Rbar, TTbar, '
       f'so no parity assumption)', used == {'S', 'R', 'TT'})

# DMS Ising energy line
beta, g = sp.symbols('beta g')
sinx, cosx = -g/(1 + g**2/4), (1 - g**2/4)/(1 + g**2/4)          # sin chi (DMS after (3.4)), continuous cos chi
T_is = cosx*sp.sinh(beta)/(sp.sinh(beta) - sp.I*sinx)
R_is = sp.I*sinx*sp.cosh(beta)/(sp.sinh(beta) - sp.I*sinx)
pts_g = [(0.37, 0.8), (-1.3, 0.3), (0.9, -2.7)]
num_ = lambda e, b_, g_: complex(e.subs({beta: b_, g: g_}).evalf(20))
u1 = max(abs(num_(R_is*R_is.subs(beta, -beta) + T_is*T_is.subs(beta, -beta) - 1, b_, g_)) for b_, g_ in pts_g)
u2 = max(abs(num_(R_is*T_is.subs(beta, -beta) + T_is*R_is.subs(beta, -beta), b_, g_)) for b_, g_ in pts_g)
c1 = max(abs(num_(T_is - T_is.subs(beta, sp.I*sp.pi - beta), b_, g_)) for b_, g_ in pts_g)
c2 = max(abs(num_(R_is - (-1)*R_is.subs(beta, sp.I*sp.pi - beta), b_, g_)) for b_, g_ in pts_g)
d1 = max(abs(num_(T_is.subs(g, 4/g) + T_is, b_, g_)) + abs(num_(R_is.subs(g, 4/g) - R_is, b_, g_)) for b_, g_ in pts_g)
z2 = max(abs(num_(T_is, b_, 2)) + abs(num_(T_is, b_, -2)) for b_, _ in pts_g)
report(f'DMS Ising energy line: R R(-) + T T(-) = 1 [{u1:.0e}], R T(-) + T R(-) = 0 [{u2:.0e}] (DMS (2.8), parity: '
       f'Rbar = R, TTbar = T); crossing T(b) = T(i pi - b), R(b) = S(2b) R(i pi - b) with S = -1 [{max(c1, c2):.0e}] '
       f'(DMS (2.9)); RT relations hold since S = -1 is constant; duality T(4/g) = -T, R(4/g) = R [{d1:.0e}]; '
       f'T = 0 at g = +-2 [{z2:.0e}]', max(u1, u2, c1, c2, d1, z2) < 1e-15)

# ================================================================ c. Lee-Yang defect
fx = lambda x: (lambda th: (mp.sinh(th) + 1j*mp.sin(mp.pi*x))/(mp.sinh(th) - 1j*mp.sin(mp.pi*x)))
S_ly = fx(mp.mpf(2)/3)
blk = lambda x: (lambda th: mp.sinh(th/2 + 1j*mp.pi*x/12)/mp.sinh(th/2 - 1j*mp.pi*x/12))   # (x), h = 6 (ch. 28)
brk = lambda x: (lambda th: 1j*mp.sinh(th/2 + 1j*mp.pi*x/12)/mp.sinh(th/2 + 1j*mp.pi*x/12 - 1j*mp.pi/2))  # BS [x]
TT_ly = lambda b: (lambda th: S_ly(th - 1j*mp.pi*(3 - b)/6))
TTb_ly = lambda b: (lambda th: S_ly(th + 1j*mp.pi*(3 - b)/6))
K1 = lambda th: blk(1)(th)*blk(3)(th)/blk(4)(th)                                          # eq-ly-reflection
Kb = lambda b: (lambda th: K1(th)*blk(1 + b)(th)*blk(5 - b)(th)/(blk(1 - b)(th)*blk(5 + b)(th)))
pts = [mp.mpc('0.37', '0.11'), mp.mpc('-0.8', '0.3'), mp.mpc('1.3', '-0.2')]
bs = [mp.mpf('0.4'), mp.mpf('-1.3'), mp.mpf('1.7'), mp.mpc('0.6', '0.8')]
e_form = max(abs(TT_ly(b)(t) - brk(b + 1)(t)*brk(b - 1)(t)) + abs(TTb_ly(b)(t) - brk(5 - b)(t)*brk(-5 - b)(t))
             for b in bs for t in pts)
report(f'Bajnok-Simon [b+1][b-1] = f_(2/3)(theta - i pi(3-b)/6) and [5-b][-5-b] = f_(2/3)(theta + i pi(3-b)/6) '
       f'[{float(e_form):.0e}]', e_form < 1e-25)
e_ax = 0
for b in bs:
    T, Tb = TT_ly(b), TTb_ly(b)
    for t in pts:
        e_ax = max(e_ax, abs(T(t)*Tb(-t) - 1), abs(T(t) - Tb(1j*mp.pi - t)),
                   abs(T(t + 1j*mp.pi/3)*T(t - 1j*mp.pi/3) - T(t)), abs(Tb(t + 1j*mp.pi/3)*Tb(t - 1j*mp.pi/3) - Tb(t)),
                   abs(Tb(t)*K1(t)*T(t) - Kb(b)(t)))
report(f'Lee-Yang TT_b: unitarity TT(theta) TTbar(-theta) = 1, crossing TT(theta) = TTbar(i pi - theta), bootstrap at '
       f'u = 2 pi/3 for TT and TTbar, fusion TTbar K_(1) TT = K_b (eq-ly-reflection) [{float(e_ax):.0e}] '
       f'(b = 0.4, -1.3, 1.7, 0.6+0.8i)', e_ax < 1e-25)
# a crossing-violating modification is caught (unitarity kept)
gb = blk(mp.mpf('1.3'))
T_bad, Tb_bad = (lambda t: TT_ly(mp.mpf('0.4'))(t)*gb(t)), (lambda t: TTb_ly(mp.mpf('0.4'))(t)*gb(t))
eu = max(abs(T_bad(t)*Tb_bad(-t) - 1) for t in pts); ec = max(abs(T_bad(t) - Tb_bad(1j*mp.pi - t)) for t in pts)
report(f'control: TT (x) with a block (1.3) keeps unitarity [{float(eu):.0e}] and violates crossing [{float(ec):.2f}]',
       eu < 1e-25 and ec > 1e-3)
# defect bound state and excited defect
ok = True
for b in (mp.mpf('-0.5'), mp.mpf('0.6'), mp.mpf('1.4')):
    u = (b + 1)*mp.pi/6
    ok &= abs(1/TTb_ly(b)(1j*u + mp.mpf('1e-20'))) < 1e-15
    ok &= max(abs(S_ly(t + 1j*u)*TT_ly(b)(t) - TT_ly(b + 2)(t)) for t in pts) < 1e-25
    ok &= max(abs(S_ly(t - 1j*u)*TTb_ly(b)(t) - TTb_ly(b + 2)(t)) for t in pts) < 1e-25
    ok &= abs(mp.sin(b*mp.pi/6) + mp.cos(u) - mp.sin((b + 2)*mp.pi/6)) < 1e-25
    ok &= abs(mp.cos(mp.pi*(3 - b)/6) - mp.sin(b*mp.pi/6)) < 1e-25
report('TTbar_b has a pole at i pi(b+1)/6; excited defect TT = S(theta + iu) TT_b = TT_{b+2}, TTbar = S(theta - iu) '
       'TTbar_b = TTbar_{b+2}; m sin(b pi/6) + m cos u = m sin((b+2)pi/6); m cos(pi(3-b)/6) = m sin(b pi/6) '
       '(b = -0.5, 0.6, 1.4)', ok)
ph = max(abs(abs(TT_ly(bb)(t)) - 1) + abs(abs(TTb_ly(bb)(t)) - 1)
         for bb in (mp.mpc(3, '0.7'), mp.mpc(-3, '-1.9')) for t in (mp.mpf('0.3'), mp.mpf('-1.7'), mp.mpf('2.2')))
report(f'TT_b and TTbar_b are pure phases for real theta at b = 3 + 0.7i and b = -3 - 1.9i [{float(ph):.0e}]; at b = 3 '
       f'TT = S(theta) (a standing particle)', ph < 1e-25 and max(abs(TT_ly(3)(t) - S_ly(t)) for t in pts) < 1e-25)

# sinh-Gordon defect (Bajnok-Simon (16)): a defect CDD factor, the square root of a bulk one
cc = mp.mpf('0.37')                                                # B vartheta/2
Tm = lambda t: -1j*mp.sinh(t/2 - 1j*mp.pi/4 + cc)/mp.sinh(t/2 + 1j*mp.pi/4 + cc)
Tp = lambda t: 1j*mp.sinh(t/2 - 1j*mp.pi/4 - cc)/mp.sinh(t/2 + 1j*mp.pi/4 - cc)
e_sh = max(max(abs(Tp(-t)*Tm(t) - 1), abs(Tm(t) - Tp(1j*mp.pi - t)), abs(Tm(t)*Tm(t + 1j*mp.pi) - 1),
               abs(Tm(t)**2*fx(mp.mpf(1)/2)(t + 2*cc) - 1)) for t in pts)
report(f'sinh-Gordon defect TT = -i sinh(theta/2 - i pi/4 + c)/sinh(theta/2 + i pi/4 + c) (Bajnok-Simon (16)): '
       f'unitarity, crossing, the defect CDD condition TT(theta) TT(theta + i pi) = 1, and TT(theta)^2 = '
       f'1/f_(1/2)(theta + 2c) [{float(e_sh):.0e}]', e_sh < 1e-25)

# ================================================================ d. doubling trick
def doubled_S(R):
    """Doubled bulk on (V+V') (x) (V+V'): R on the first copy, P R P on the reflected copy, flip between copies."""
    d = 2*dV
    SD = np.zeros((d*d, d*d), complex)
    Rt = R.reshape(dV, dV, dV, dV)                                 # [k, l, i, j]: (k,l) <- (i,j)
    for i in range(d):
        for j in range(d):
            ci, cj = i // dV, j // dV; li, lj = i % dV, j % dV
            if ci == cj == 0:
                for k in range(dV):
                    for l in range(dV):
                        SD[k*d + l, i*d + j] = Rt[k, l, li, lj]
            elif ci == cj == 1:
                for k in range(dV):
                    for l in range(dV):
                        SD[(dV + k)*d + dV + l, i*d + j] = Rt[l, k, lj, li]     # P R P
            else:
                SD[j*d + i, i*d + j] = 1
    return SD


def doubled_K(T, Tb):
    """K on (V+V') (x) W: phi -> psi' with P TT, psi' -> phi with TTbar P (eq-doubled-k)."""
    d = 2*dV
    K = np.zeros((d*dW, d*dW), complex)
    PT = flipmat(dW, dV) @ T                                       # V (x) W -> V (x) W
    TbP = Tb @ flipmat(dV, dW)
    for a in range(dV):
        for al in range(dW):
            for b in range(dV):
                for be in range(dW):
                    K[(dV + b)*dW + be, a*dW + al] = PT[b*dW + be, a*dW + al]
                    K[b*dW + be, (dV + a)*dW + al] = TbP[b*dW + be, a*dW + al]
    return K


def re_doubled(TTf):
    TbF = lambda y: np.linalg.inv(TTf(1/y))                        # TTbar(theta) = TT(-theta)^-1
    d = 2*dV
    SD_m, SD_p = doubled_S(Rc(y1, y2)), doubled_S(Rc(y1*y2, 1.0))
    K1_, K2_ = doubled_K(TTf(y1), TbF(y1)), doubled_K(TTf(y2), TbF(y2))
    I = np.eye(dW); Id = np.eye(d)
    L = kron(SD_m, I) @ kron(Id, K1_) @ kron(SD_p, I) @ kron(Id, K2_)
    Rr = kron(Id, K2_) @ kron(SD_p, I) @ kron(Id, K1_) @ kron(SD_m, I)
    res = (L - Rr).reshape(d, d, dW, d, d, dW)
    # phi phi -> psi' psi' block; output psi'(a) psi'(b) w  <->  w psi(b) psi(a) of the triangle relation
    blk_ = res[dV:, dV:, :, :dV, :dV, :]                          # [a', b', g, a, b, al]
    tri = (kron(np.eye(dW), Rc(y1, y2)) @ kron(TTf(y1), np.eye(dV)) @ kron(np.eye(dV), TTf(y2))
           - kron(TTf(y2), np.eye(dV)) @ kron(np.eye(dV), TTf(y1)) @ kron(Rc(y1, y2), np.eye(dW)))
    tri = tri.reshape(dW, dV, dV, dV, dV, dW)                      # [g, e, f, a, b, al]
    blk_tri = blk_.transpose(2, 1, 0, 3, 4, 5)                     # psi'(a') psi'(b') -> w psi(b') psi(a')
    KK = max(np.abs(doubled_K(TTf(yy), TbF(yy)) @ doubled_K(TTf(1/yy), TbF(1/yy)) - np.eye(d*dW)).max()
             for yy in (y1, y2))
    return np.abs(res).max()/np.abs(L).max(), np.abs(blk_tri - tri).max(), np.abs(tri).max(), KK


r0, m0, t0, kk0 = re_doubled(TT_sol)
report(f'doubling trick: doubled K of a triangle solution satisfies the braid reflection equation [{r0:.0e}] and '
       f'K(theta) K(-theta) = 1 [{kk0:.0e}]', r0 < 1e-12 and kk0 < 1e-12)
r1, m1, t1_, kk1 = re_doubled(lambda y: TT_sol(y) + 0.05*Nfix)
report(f'perturbed TT: doubled reflection equation fails [{r1:.2f}]; its phi phi -> psi\' psi\' block equals the '
       f'triangle residual entry by entry [{m1:.0e}, size {t1_:.2f}]; unitarity still holds [{kk1:.0e}]',
       r1 > 1e-3 and m1 < 1e-12 and t1_ > 1e-3 and kk1 < 1e-10)


def doubled_cu(T, Tb):
    """Doubled crossing-unitarity k(theta) = Scheck_D(2 theta) k(-theta), k = sum K_{b abar}(i pi/2 - theta) v_a v_b;
    species phi (0), psi' (1), both self-conjugate; K_{psi' phi} = TT, K_{phi psi'} = TTbar."""
    def Kmat(th):
        return mp.matrix([[0, Tb(th)], [T(th), 0]])               # K[b, a] = amplitude a -> b

    def kvec(th):
        K = Kmat(1j*mp.pi/2 - th)
        return [K[b, a] for a in range(2) for b in range(2)]      # index a*2 + b

    def SD(th):                                                    # braid S on (phi, psi') pairs
        M = mp.zeros(4, 4)
        for a in range(2):
            for b in range(2):
                M[b*2 + a, a*2 + b] = S_ly(th) if a == b else 1
        return M
    err = 0
    for t in pts:
        kv = mp.matrix(kvec(-t))
        rhs = SD(2*t)*kv
        err = max(err, max(abs(kvec(t)[i] - rhs[i]) for i in range(4)))
    return err


ecu = doubled_cu(TT_ly(mp.mpf('0.4')), TTb_ly(mp.mpf('0.4')))
ecb = doubled_cu(T_bad, Tb_bad)
report(f'Lee-Yang: the doubled boundary satisfies crossing-unitarity eq-boundary-crossing-unitarity [{float(ecu):.0e}]; '
       f'the crossing-violating control fails it [{float(ecb):.2f}]', ecu < 1e-25 and ecb > 1e-3)

# ================================================================ e. twist defect = Bloch sector
OM = 0.71
qa = q_phys(OM)
Va = lambda y: an_vector(2, qa, y)
Ra = lambda y1_, y2_: rcheck(Va(y1_), Va(y2_))[0]
mu = np.eye(3) - 1/3                                               # weights eps_i of the vector representation
alpha = np.array([0.83, -0.27, 0.0])
Om = np.diag(np.exp(1j*mu @ alpha))
d3 = 3
ys = [np.exp(OM*t) for t in (0.52, -0.31, 1.07)]
Nn = 3
dims = [d3]*Nn
Rfix = lambda i, j: flipmat(d3, d3) @ Ra(ys[i], ys[j])            # S_ij on fixed slots (particle i in slot i)
S_ = {(i, j): op_on(Rfix(i, j), [i, j], dims) for i in range(Nn) for j in range(Nn) if i != j}
Omj = {j: op_on(Om, [j], dims) for j in range(Nn)}
tri_tw = np.abs(Ra(ys[0], ys[1]) @ kron(Om, Om) - kron(Om, Om) @ Ra(ys[0], ys[1])).max()
report(f'twist defect TT = exp(i alpha.mu) (dim W = 1) satisfies the triangle relation for the a_2 vector R-matrix '
       f'[{tri_tw:.0e}]; TTbar = exp(-i alpha.mu) gives unitarity and crossing (mu_abar = -mu_a)', tri_tw < 1e-12)


def T_def(j, k):
    """Defect Bethe-Yang transfer matrix of particle j, defect between positions k and k+1 (k = Nn-1: wrap point).
    Particle j crosses j+1, ..., Nn-1, then 0, ..., j-1 (cyclically), and the defect when it passes position k."""
    M = np.eye(d3**Nn, dtype=complex)
    seq = [(j + s) % Nn for s in range(1, Nn)]
    pos = j
    for other in seq:
        if pos == k:
            M = Omj[j] @ M
        M = S_[(j, other)] @ M
        pos = other
    if pos == k:
        M = Omj[j] @ M
    return M


charges = np.array([sum(mu[i] for i in idx) for idx in itertools.product(range(d3), repeat=Nn)])
sectors = {}
for n_, c in enumerate(np.round(charges, 6)):
    sectors.setdefault(tuple(c), []).append(n_)
ok_spec, ok_comm = True, True
for k in range(Nn):
    T0 = T_def(0, k)
    Tbl = Omj[0] @ T_def(0, -1)                                   # Bloch: exp(i alpha.mu_1) T_1 (no defect)
    for idx in sectors.values():
        ev1 = np.sort_complex(np.linalg.eigvals(T0[np.ix_(idx, idx)]))
        ev2 = np.sort_complex(np.linalg.eigvals(Tbl[np.ix_(idx, idx)]))
        ok_spec &= np.abs(ev1 - ev2).max() < 1e-9
    for i in range(Nn):
        for j in range(i + 1, Nn):
            Ti, Tj = T_def(i, k), T_def(j, k)
            ok_comm &= np.abs(Ti @ Tj - Tj @ Ti).max() < 1e-10
leak = max(np.abs(T_def(0, k)[np.ix_(a, b)]).max() for k in range(Nn)
           for a in sectors.values() for b in sectors.values() if a is not b)
report(f'twist defect on a circle, N = 3 a_2 solitons: for every defect position the transfer matrices T_j^D commute, '
       f'preserve the charge sectors [{leak:.0e}], and T_1^D has in each sector the spectrum of the Bloch-twisted '
       f'exp(i alpha.mu_1) T_1 (prp-bloch-sectors)', ok_spec and ok_comm and leak < 1e-12)

Tcyc = {j: Omj[j] @ T_def(j, -5) for j in range(Nn)}               # twist applied last on each cyclic T_j
ok_sp = all(np.abs(np.sort_complex(np.linalg.eigvals(Tcyc[j][np.ix_(idx, idx)]))
                   - np.sort_complex(np.linalg.eigvals(T_def(j, Nn - 1)[np.ix_(idx, idx)]))).max() < 1e-9
            for j in range(Nn) for idx in sectors.values())
ncomm = max(np.abs(Tcyc[i] @ Tcyc[j] - Tcyc[j] @ Tcyc[i]).max() for i in range(Nn) for j in range(i + 1, Nn))
report(f'the twist must sit at a common seam: exp(i alpha.mu_j) applied after each cyclic T_j gives the right spectra '
       f'sector by sector but matrices that do not commute [{ncomm:.2f}]; T_1^D with the defect at the seam is exactly '
       f'exp(i alpha.mu_1) T_1 [{np.abs(T_def(0, Nn - 1) - Tcyc[0]).max():.0e}]',
       ok_sp and ncomm > 1e-3 and np.abs(T_def(0, Nn - 1) - Tcyc[0]).max() < 1e-14)

# ================================================================ f. modified crossing of Copetti-Cordova-Komatsu
n_ = 5
vac = [mp.mpf(i)/2 for i in range(n_ - 1)]                          # a = 0, 1/2, ..., n/2 - 1
dq = lambda a: mp.sin((2*a + 1)*mp.pi/n_)/mp.sin(mp.pi/n_)
adj = lambda x, y: abs(x - y) == mp.mpf(1)/2


def S18(a, b, c, d, th):
    if not (adj(a, b) and adj(b, c) and adj(c, d) and adj(d, a)):
        return None
    r = dq(a)*dq(c)/(dq(b)*dq(d))
    return r**(1j*th/(2*mp.pi))*(mp.sqrt(r)*mp.sinh(th/n_)*(b == d) + mp.sinh((1j*mp.pi - th)/n_)*(a == c))


def S21(a, b, c, d, th):
    if not (adj(a, b) and adj(b, c) and adj(c, d) and adj(d, a)):
        return None
    r = dq(a)*dq(c)/(dq(b)*dq(d))
    return mp.sqrt(r)*mp.sinh(th/n_)*(b == d) + mp.sinh((1j*mp.pi - th)/n_)*(a == c)


e18 = e21 = eg = 0
e21std = 0
for a, b, c, d in itertools.product(vac, repeat=4):
    for t in pts:
        s18 = S18(a, b, c, d, t)
        if s18 is None:
            continue
        r = dq(a)*dq(c)/(dq(b)*dq(d))
        e18 = max(e18, abs(s18 - S18(b, c, d, a, 1j*mp.pi - t)))
        e21 = max(e21, abs(S21(a, b, c, d, t) - mp.sqrt(r)*S21(b, c, d, a, 1j*mp.pi - t)))
        e21std = max(e21std, abs(S21(a, b, c, d, t) - S21(b, c, d, a, 1j*mp.pi - t)))
        eg = max(eg, abs(s18 - r**(1j*t/(2*mp.pi))*S21(a, b, c, d, t)))
report(f'M_5 + phi_13 kinks (CCK (18), (21), d_a = sin((2a+1)pi/5)/sin(pi/5)): (18) obeys standard crossing '
       f'[{float(e18):.0e}], (21) the modified crossing (22) [{float(e21):.0e}] but not the standard one '
       f'[{float(e21std):.2f}]; (18) = (d_a d_c/(d_b d_d))^(i theta/2 pi) (21) [{float(eg):.0e}]',
       max(e18, e21, eg) < 1e-25 and e21std > 1e-3)
