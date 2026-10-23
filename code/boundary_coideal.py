"""Coideal subalgebras and soliton reflection matrices (Part VII, chapters 29 and 32).

Part VI conventions throughout (boundary_common.py): Delta(e)=e(x)1+k(x)e, Delta(f)=f(x)k^-1+1(x)f,
physical q = -exp(-i pi omega), principal gradation with y = x^(1/h).  The coideal generators are
b_j = e_j + q^-1 f_j k_j + eps_j k_j (parity-symmetric normalization), and a K-matrix is an
intertwiner K b_j = b_j K.

Checks:
 1. coproduct identity Delta(b_j) = (b_j - eps_j k_j)(x)1 + k_j(x)b_j on tensor products; the
    parity-symmetric form q^-1 f_j k_j = k_j^(1/2) f_j k_j^(1/2); d_n^(1) vector rep relations.
 2. sine-Gordon (worked example): K: V(y) -> V(1/y) from the linear equations; a solution exists
    only if c_0 c_1 = q^-2; closed form with c_j = q^-1.
 3. a_n^(1) vector solitons, n = 2..7: K: V(y) -> Vbar(1/y) exists iff all eps_j = 0 (then K = 1)
    or all eps_j^2 = 1/((1-q)(1-q^-1)) = 1/(4 cos^2(pi omega/2)); not for |eps_j| = 1 (DM (4.21)).
    Closed form (corrected DM (4.21)) and the substitution that maps DM (4.21) onto it; gauge
    equivalence of sign patterns.  For n = 2, 3, 4: twisted reflection equation for all four (mu, nu); failure
    of the printed K at -q (where the antisoliton is Vbar(-y) and eps = 0 has a solution only for odd n).
 4. Identification of the solutions: the non-diagonal s=+1 solution is Gandenberger's K^+ for a_2
    (hep-th/9806003 (3.7), his Neumann candidate), s=-1 is his K^-; eps = 0 gives DG (3.7) (K = 1).
    Boundary crossing-unitarity (GZ (3.35)) in book form k(theta) = f R_{Vbar Vbar}(2theta) k(-theta),
    k = sum_ab K_ba(i pi/2 - theta) v_a (x) v_b, holds for the matrix part (n = 2, 3).
 5. The 2x2 systems behind eq-an-eps-condition (symbolic, a_3); DM's own corrected (4.15)-(4.16) are solved by
    their printed (4.21) only after q_(4.21) -> -(-q)^(1/2) (n = 2..7).
 6. Hypothesis (ii) of thm-k-reflection: dim End_B(V_mu (x) V_nu) = 1 (a_1, a_2, a_3).
 7. prp-boundary-module against BPTT: one-dimensional B-submodules of V(y) at y = exp(-i lambda nu_n) carry the
    beta_j of the excited boundary states.
 8. B = Lbar (K (x) 1) L (eq-b-matrix) at first order in 1/y: leading terms of the universal R-matrix against the
    exact R-check; for the a_n vector representation (n <= 6) the coefficients are exactly the b_j of
    eq-coideal-generators, with eps = 0 (K = 1) or eps = (s eps_*, eps_*, ..., eps_*) (eq-an-k-matrix).
 9. Sine-Gordon eps_j near the Neumann point from eq-bsg-eta-theta and cnj-bsg-uv-ir (prp-dm-first-order).
Runtime: about a minute.
"""
import numpy as np
import sympy as sp
from boundary_common import (report, q_phys, an_vector, sl2_spin, dn_vector, check_relations, cartan_a,
                             cartan_d, coideal_gen, coideal_gen_tensor, cop, k_intertwiner, rcheck,
                             re_twisted, re_preserving, dual, kron, intertwiners, nodes_of)

rng = np.random.default_rng(7)
OMEGA = 0.37
q = q_phys(OMEGA)
EPS_STAR = 1/np.sqrt((1 - q)*(1 - 1/q))          # = 1/(2 cos(pi omega/2))

# ---------------------------------------------------------------- 1. coproduct identity
ok = True
for A, B in [(sl2_spin(1, q, 0.7+0.2j), sl2_spin(2, q, 1.3-0.4j)), (an_vector(2, q, 0.8j+0.3), an_vector(2, q, 1.1, True))]:
    for j in range(len([k for k in A if isinstance(k, tuple) and k[0] == 'e'])):
        c, e = 1/q, rng.normal() + 1j*rng.normal()
        lhs = cop(A, B, 'e', j) + c*cop(A, B, 'f', j) @ cop(A, B, 'k', j) + e*cop(A, B, 'k', j)
        bA, bB = coideal_gen(A, j, c, e), coideal_gen(B, j, c, e)
        rhs = kron(bA - e*A[('k', j)], np.eye(B['dim'])) + kron(A[('k', j)], bB)
        ok &= np.allclose(lhs, rhs) and np.allclose(lhs, coideal_gen_tensor(A, B, j, c, e))
report('Delta(b_j) = (b_j - eps_j k_j) (x) 1 + k_j (x) b_j  (left coideal)', ok)
ok = True
for V in (sl2_spin(1, q, 0.7), sl2_spin(2, q, 0.7), an_vector(3, q, 0.7), an_vector(3, q, 0.7, True)):
    for j in (0, 1):
        k = V[('k', j)]; kh = np.diag(np.sqrt(np.diag(k)))
        # k_j^(1/2) with the branch q^(m/2) fixed by the weight: compare up to that branch
        lhs = V[('f', j)] @ k/q
        rhs = kh @ V[('f', j)] @ kh
        ok &= np.allclose(lhs, rhs) or np.allclose(lhs, -rhs)
report('parity-symmetric normalization: q^-1 f_j k_j = k_j^(1/2) f_j k_j^(1/2)', ok)
report('d_n^(1) vector representation (n=4,5) satisfies the q-Serre relations',
       max(check_relations(dn_vector(n, q, 0.6+0.5j), cartan_d(n), q) for n in (4, 5)) < 1e-12)
report('a_n vector and conjugate vector, spin-j: Drinfeld-Jimbo relations',
       max([check_relations(an_vector(n, q, 0.6+0.5j, cj), cartan_a(n), q) for n in (2, 3, 4) for cj in (0, 1)] +
           [check_relations(sl2_spin(j2, q, 0.6+0.5j), cartan_a(1), q) for j2 in (1, 2, 3)]) < 1e-12)
report('Vbar is the antisoliton: V(z)* = Vbar(-z/q) with C = 1 (principal gradation)',
       all(max(np.abs(dual(an_vector(n, q, 0.8+0.1j))[k] - an_vector(n, q, -(0.8+0.1j)/q, True)[k]).max()
               for k in dual(an_vector(n, q, 1.0)) if isinstance(k, tuple)) < 1e-14 for n in (2, 3, 4)))

# DM section 4 factor 2: (4.26) Q_i Qbar_i - q^{-alpha_i.alpha_i} Qbar_i Q_i = (q^{2T_i}-1)/(q_i^2-1), |alpha|^2=1,
# q_i = q^(1/2), on the DM representation (4.13) Q_i = x e_{i+1,i}, Qbar_i = x^-1 e_{i,i+1}: fails with
# T_i = -e_ii + e_{i+1,i+1} as printed, holds with T_i halved
qd = 0.3 + 0.9j; Nd = 4; xd = 1.3
def dm_rel(halve):
    res = 0
    for i in range(Nd):
        a, b = i, (i + 1) % Nd
        Q = np.zeros((Nd, Nd), complex); Q[b, a] = xd
        Qb = np.zeros((Nd, Nd), complex); Qb[a, b] = 1/xd
        T = np.zeros(Nd); T[a] -= 1; T[b] += 1
        T = T/2 if halve else T
        res = max(res, np.abs(Q @ Qb - qd**-1*Qb @ Q - np.diag((qd**(2*T) - 1)/(qd - 1))).max())
    return res
report(f'DM (4.26) on DM (4.13): fails as printed [{dm_rel(False):.2f}], holds with T_i -> T_i/2 [{dm_rel(True):.0e}]',
       dm_rel(False) > 1e-3 and dm_rel(True) < 1e-14)

# ---------------------------------------------------------------- 2. sine-Gordon from the linear equations
qs, ys, e0, e1, c0, c1 = sp.symbols('q y epsilon_0 epsilon_1 c_0 c_1', nonzero=True)
Es = sp.Matrix([[0, 1], [0, 0]]); Fs = Es.T; Ks = sp.diag(qs, 1/qs)


def bsg(j, yy, c, e):
    return yy*Es + c*Fs*Ks/yy + e*Ks if j == 1 else yy*Fs + c*Es*Ks.inv()/yy + e*Ks.inv()


Km = sp.Matrix(2, 2, sp.symbols('k0:4'))
eqs = []
for j, c, e in ((0, c0, e0), (1, c1, e1)):
    eqs += list(Km*bsg(j, ys, c, e) - bsg(j, 1/ys, c, e)*Km)
A_, _ = sp.linear_eq_to_matrix(eqs, list(Km))
A0 = sp.Matrix(A_).subs({e0: 0, e1: 0})
rank_good = A0.subs(c1, 1/(qs**2*c0)).rank(simplify=True)
rank_bad = A0.subs(c1, 2/(qs**2*c0)).rank(simplify=True)
report('sine-Gordon, eps=0: a nonzero K exists iff c_0 c_1 = q^-2 (it is then purely off-diagonal)',
       rank_good < 4 and rank_bad == 4)
eqs = []
for j, e in ((0, e0), (1, e1)):
    eqs += list(Km*bsg(j, ys, 1/qs, e) - bsg(j, 1/ys, 1/qs, e)*Km)
sol = sp.solve(eqs, list(Km), dict=True)[0]
Ksol = Km.subs(sol)
Ksol = sp.simplify(Ksol/Ksol[0, 1])
Kclosed = sp.Matrix([[(qs - 1/qs)*(e1*ys + e0/ys), ys**2 - ys**-2],
                     [ys**2 - ys**-2, (qs - 1/qs)*(e0*ys + e1/ys)]])/(ys**2 - ys**-2)
report('sine-Gordon K(theta) = [[(q-q^-1)(e1 y+e0/y), y^2-y^-2],[y^2-y^-2, (q-q^-1)(e0 y+e1/y)]] (unique)',
       sp.simplify(Ksol - Kclosed) == sp.zeros(2) and len(sp.solve(eqs, list(Km), dict=True)) == 1)

# ---------------------------------------------------------------- 3. a_n vector solitons
def k_an(n, y, eps, qq=q, ybar_sign=1, conj_in=False):
    """K: V(y) -> Vbar(ybar_sign/y) (or Vbar -> V if conj_in)."""
    N = n + 1
    Vin = an_vector(n, qq, y, conj_in); Vout = an_vector(n, qq, ybar_sign/y, not conj_in)
    null, rel = k_intertwiner(Vin, Vout, [1/qq]*N, eps)
    return null, rel


def k_closed(n, y, s, qq=q):
    """Closed form for eps_1 = ... = eps_n = eps_*, eps_0 = s eps_* (gauge-fixed; s = sign of prod eps_j/eps_*).
    X = (-q)^(1/2) y = exp(omega(theta - i pi/2)), w = (-q)^(1/2) = exp(-i pi omega/2), N = n+1:
    K_ij = X^(j-i-N/2) (j > i),  K_ij = s X^(j-i+N/2) (j < i),
    K_ii = (w X^(-N/2) - s w^-1 X^(N/2)) / (eps_* (q^-1 - q))."""
    N = n + 1
    w = np.sqrt(-qq + 0j)
    es = 1/np.sqrt((1 - qq)*(1 - 1/qq))
    X = w*y
    K = np.zeros((N, N), complex)
    for i in range(N):
        for j in range(N):
            if i == j:
                K[i, j] = (w*X**(-N/2) - s/w*X**(N/2))/(es*(1/qq - qq))
            elif j > i:
                K[i, j] = X**(j - i - N/2)
            else:
                K[i, j] = s*X**(j - i + N/2)
    return K


def k_dm(n, xdm, qdm, epsdm):
    """Delius-MacKay (4.21) as printed, K^i_j with indices 1..n+1 (eps_hat = prod eps_i)."""
    N = n + 1
    eh = np.prod(epsdm)
    K = np.zeros((N, N), complex)
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            if i == j:
                K[i-1, j-1] = (1/qdm*(-qdm*xdm)**(N/2) - eh*qdm*(-qdm*xdm)**(-N/2))/(1/qdm - qdm)
            elif j > i:
                K[i-1, j-1] = np.prod(epsdm[i:j])*(-qdm*xdm)**(i - j + N/2)          # eps_i ... eps_{j-1}
            else:
                K[i-1, j-1] = np.prod(epsdm[j:i])*eh*(-qdm*xdm)**(i - j - N/2)      # K^i_j, i > j
    return K


def prop(A, B):
    c = np.vdot(B.ravel(), A.ravel())/np.vdot(B.ravel(), B.ravel())
    return np.linalg.norm(A - c*B)/np.linalg.norm(A)


y0 = 0.9 + 0.4j
for n in (2, 3, 4, 5, 6, 7):
    N = n + 1
    cases = {'generic eps': (rng.normal(size=N) + 1j*rng.normal(size=N), 0),
             'eps = 0': (np.zeros(N), 1),
             'all eps_j = +eps_*': (EPS_STAR*np.ones(N), 1),
             'eps_j = +-eps_* (random signs)': (EPS_STAR*rng.choice([-1, 1], size=N), 1),
             'all |eps_j| = 1 (DM normalization), eps_j = 1': (np.ones(N), 0),
             'all |eps_j| = 1, random phases': (np.exp(2j*np.pi*rng.random(N)), 0),
             'one eps_j = 0, others eps_*': (np.r_[0, EPS_STAR*np.ones(n)], 0)}
    ok = True
    for name, (e, want) in cases.items():
        null, rel = k_an(n, y0, e)
        ok &= len(null) == want
    report(f'a_{n}: K: V(y)->Vbar(1/y) exists (uniquely) iff eps = 0 or all eps_j^2 = 1/((1-q)(1-q^-1)); '
           'not for |eps_j| = 1', ok)
    null, _ = k_an(n, y0, np.zeros(N))
    report(f'a_{n}: eps = 0 gives K = 1 (not an arbitrary diagonal matrix)',
           np.allclose(null[0]/null[0][0, 0], np.eye(N)))
    ok = True
    for s in (1, -1):
        for yy in (y0, 1.7 - 0.3j):
            e = EPS_STAR*np.ones(N); e[0] *= s
            null, _ = k_an(n, yy, e)
            ok &= prop(null[0], k_closed(n, yy, s)) < 1e-11
    report(f'a_{n}: closed form K(theta) (corrected DM (4.21)), s = +1 and -1', ok)
    # DM (4.21) as printed vs the book K: equal (s=+1) or diagonal-gauge equivalent (s=-1) after basis
    # reversal, q_DM -> -(-q)^(1/2) = -exp(-i pi omega/2), x_DM -> y (so -q_DM x_DM = X), eps_i -> eps_i/eps_*
    J = np.eye(N)[::-1]; qdm = -np.sqrt(-q + 0j)
    Kp = J @ k_dm(n, y0, qdm, np.ones(N)) @ J
    Km_ = J @ k_dm(n, y0, qdm, np.r_[np.ones(n), -1]) @ J
    ratio = k_closed(n, y0, -1)/Km_
    sv = np.linalg.svd(ratio, compute_uv=False)
    report(f'a_{n}: DM (4.21) = book K after basis reversal, q_DM -> -(-q)^(1/2), x_DM -> y, eps_i -> eps_i/eps_* '
           '(s=+1 equal, s=-1 up to a diagonal gauge)', prop(Kp, k_closed(n, y0, 1)) < 1e-12 and sv[1]/sv[0] < 1e-12)
    # gauge: different sign patterns with the same product give the same diagonal and |entries|
    e1_ = EPS_STAR*np.ones(N); e1_[1] *= -1; e1_[2] *= -1
    K1 = k_an(n, y0, e1_)[0][0]; K0 = k_closed(n, y0, 1)
    K1 = K1/K1[0, 0]*K0[0, 0]
    report(f'a_{n}: sign patterns with equal product s = prod sign(eps_j) are gauge equivalent',
           np.allclose(np.diag(K1), np.diag(K0)) and np.allclose(np.abs(K1), np.abs(K0)))

# twisted reflection equation, all four (mu, nu), physical q
for n in (2, 3, 4):
    N = n + 1
    cache = {}

    def Rf(a, b, z, qq=q, sb=1):
        key = (a, b, complex(z), qq, sb)
        if key not in cache:
            A = an_vector(n, qq, z if a == 'V' else sb*z, a == 'W')
            B = an_vector(n, qq, 1.0 if b == 'V' else sb*1.0, b == 'W')
            cache[key] = rcheck(A, B)[0]
        return cache[key]
    res_all = []
    for label, e in (('eps=0', np.zeros(N)), ('s=+1', EPS_STAR*np.ones(N)), ('s=-1', EPS_STAR*np.r_[-1, np.ones(n)])):
        def Kf(a, yy, e=e):
            return k_an(n, yy, e, conj_in=(a == 'W'))[0][0]
        pairs = [('V', 'V'), ('V', 'W'), ('W', 'V'), ('W', 'W')] if n < 4 else [('V', 'V'), ('V', 'W')]
        for mu, nu in pairs:
            r, c = re_twisted(Rf, Kf, mu, nu, 1.3 + 0.2j, 0.7 - 0.4j)
            res_all.append((label, mu, nu, r, c))
    worst = max(r for *_, r, c in res_all)
    report(f'a_{n}: twisted reflection equation (all eps solutions, (mu,nu) in {{V,Vbar}}^2) [{worst:.1e}]', worst < 1e-9)
    cmax = max(abs(c - 1) for (lab, mu, nu, r, c) in res_all if mu == nu)
    report(f'a_{n}: ... with R-checks normalized on highest weights the constant is 1 for mu = nu [{cmax:.0e}]',
           cmax < 1e-8)
    # -q: the printed (physical) K with the R-matrices of the other sign, antisoliton Vbar(-y) there
    qm = -q
    def Kprinted(a, yy):
        K = k_closed(n, yy, 1)
        return K if a == 'V' else np.linalg.inv(k_closed(n, 1/yy, 1))
    Rm = lambda a, b, z: Rf(a, b, z, qm, -1)
    r_m, _ = re_twisted(Rm, Kprinted, 'V', 'V', 1.3 + 0.2j, 0.7 - 0.4j)
    r_p, _ = re_twisted(Rf, Kprinted, 'V', 'V', 1.3 + 0.2j, 0.7 - 0.4j)
    report(f'a_{n}: printed K (s=+1) solves the twisted RE at q [{r_p:.0e}] and fails at -q [{r_m:.2f}]',
           r_p < 1e-9 and r_m > 1e-3)
    # at -q (antisoliton Vbar(-y)) the coideal solutions are different: eps=0 exists only for odd n
    null, _ = k_an(n, y0, np.zeros(N), qq=qm, ybar_sign=-1)
    report(f'a_{n}: at -q (antisoliton Vbar(-y)) an eps=0 solution exists iff n is odd', (len(null) == 1) == (n % 2 == 1))

# ---------------------------------------------------------------- 4. identification with Gandenberger and DG
# Gandenberger hep-th/9806003 (3.7) with g = h = 1: off-diagonal e^{+-i pi(mu/3 - lambda/4)}, diagonal kappa_+-(mu),
# mu = -i h lambda theta/(2 pi), lambda = omega, h = 3; kappa_+ = sin(pi(mu - lambda/4))/sin(pi lambda/2),
# kappa_- = i cos(pi(mu - lambda/4))/sin(pi lambda/2).
ok = True
for th in (0.31, -0.7 + 0.2j):
    mu = -1j*3*OMEGA*th/(2*np.pi); lam = OMEGA
    ep, em = np.exp(1j*np.pi*(mu/3 - lam/4)), np.exp(-1j*np.pi*(mu/3 - lam/4))
    for s, kap in ((1, np.sin(np.pi*(mu - lam/4))/np.sin(np.pi*lam/2)),
                   (-1, 1j*np.cos(np.pi*(mu - lam/4))/np.sin(np.pi*lam/2))):
        G = np.array([[kap, ep, s*em], [s*em, kap, ep], [ep, s*em, kap]])
        Kb = k_an(2, np.exp(OMEGA*th), s*EPS_STAR*np.ones(3))[0][0]
        ok &= min(prop(Kb, G), prop(Kb.T, G)) < 1e-12
report("a_2: s = +1 solution = Gandenberger's K^+ (his Neumann candidate), s = -1 = K^-, at the physical q", ok)

# boundary crossing-unitarity, matrix part: k(theta) = f(theta) R_{Vbar Vbar}(y(2theta)) k(-theta),
# k = sum_ab K_{b a}(i pi/2 - theta) v_a (x) v_b in Vbar (x) Vbar (C = 1 in the principal gradation)
ok, gs = True, []
for n in (2, 3):
    N = n + 1
    for s_ in (1, -1):
        for th in (0.31, 0.83):
            ym, yp = np.exp(OMEGA*(1j*np.pi/2 - th)), np.exp(OMEGA*(1j*np.pi/2 + th))
            km = k_closed(n, ym, s_).T.reshape(N*N); kp = k_closed(n, yp, s_).T.reshape(N*N)
            R2 = rcheck(an_vector(n, q, np.exp(2*OMEGA*th), True), an_vector(n, q, 1.0, True))[0]
            v = R2 @ kp
            ok &= np.linalg.norm(np.outer(km, v) - np.outer(v, km))/(np.linalg.norm(km)*np.linalg.norm(v)) < 1e-11
report('a_n (n=2,3): boundary crossing-unitarity k(theta) = f(theta) R_{Vbar Vbar}(2theta) k(-theta) holds for '
       'the matrix part of both non-diagonal K (eps = 0: trivially)', ok)

# ---------------------------------------------------------------- 5. why eps_j^2 = eps_*^2: the 2x2 systems
# For node j (alpha_j = eps_a - eps_b) and c not in {a, b}, the (a, c) and (b, c) entries of
# K pi_y(b_j) = pibar_{1/y}(b_j) K involve only K_ac, K_bc, and the (c, a), (c, b) entries only K_ca, K_cb; in both
# cases the 2x2 determinant is +-(eps_j^2 (1-q)(1-q^-1) - 1).  Symbolic, a_3, all nodes.
qs_, ys_, es_ = sp.symbols('q y epsilon', nonzero=True)
n = 3; N = n + 1
Ksym = sp.Matrix(N, N, sp.symbols(f'K0:{N*N}'))          # row-major: K[r, c] = symbol r*N + c
ok = True
for j in range(N):
    a, b = (n, 0) if j == 0 else (j - 1, j)
    Ej = sp.zeros(N); Ej[a, b] = 1; Fj = Ej.T
    kj = sp.diag(*[qs_ if i == a else (1/qs_ if i == b else 1) for i in range(N)]); kbj = kj.inv()
    bV = ys_*Ej + Fj*kj/(qs_*ys_) + es_*kj                 # pi_y(b_j) on V(y)
    bW = Fj/ys_ + ys_*Ej*kbj/qs_ + es_*kbj                 # pibar_{1/y}(b_j): e_j = F_j/y, f_j = y E_j, k_j -> k_j^-1
    Eq = Ksym*bV - bW*Ksym
    for c in range(N):
        if c in (a, b):
            continue
        for pair, idx in (((Eq[a, c], Eq[b, c]), (a*N + c, b*N + c)), ((Eq[c, a], Eq[c, b]), (c*N + a, c*N + b))):
            A2, _ = sp.linear_eq_to_matrix(list(pair), list(Ksym))
            ok &= {i for i in range(N*N) if any(A2[:, i])} == set(idx)
            det = A2[:, list(idx)].det()
            ok &= sp.simplify(det**2 - (es_**2*(1 - qs_)*(1 - 1/qs_) - 1)**2) == 0
report('a_3: for each node j the off-block entries of rows/columns a_j, b_j obey 2x2 systems with determinant '
       '+-(eps_j^2 (1-q)(1-q^-1) - 1) (origin of eq-an-eps-condition)', ok)


# DM's own equations (4.15)-(4.16), with the constant term of q^{T_i} restored, in DM's basis and at the q of their
# representation: the printed (4.21) solves them only with q_(4.21) = -(-q)^(1/2), not with q_(4.21) = q (n = 2..7)
def dm_system(n, qd, x, eps):
    N = n + 1; I = np.eye(N); rows = []
    def E(a, b):
        M = np.zeros((N, N), complex); M[a % N, b % N] = 1; return M
    for i in range(N):
        qT = I + (1/qd - 1)*E(i, i) + (qd - 1)*E(i + 1, i + 1); qTb = I + (qd - 1)*E(i, i) + (1/qd - 1)*E(i + 1, i + 1)
        M = x*E(i + 1, i) + E(i, i + 1)/x + eps[i]*qT; Mp = E(i, i + 1)/x + x*E(i + 1, i) + eps[i]*qTb
        rows.append(np.kron(M.T, I) - np.kron(I, Mp))
    return np.vstack(rows)


xdm = 0.8*np.exp(0.3j)
res_good, res_bad = 0, np.inf
for n in range(2, 8):
    N = n + 1
    for signs in [(1,)*N, (-1,) + (1,)*n, tuple(rng.choice([1, -1], N))]:
        A = dm_system(n, q, xdm, EPS_STAR*np.array(signs))
        for qdm, good in ((-np.sqrt(-q + 0j), True), (q, False)):
            # dm_system labels the node acting on basis vectors (i, i+1) by i, k_dm by i+1 (DM's labels)
            v = k_dm(n, xdm, qdm, np.roll(np.array(signs, float), 1)).reshape(-1, order='F')
            r = np.linalg.norm(A @ v)/np.linalg.norm(v)
            res_good, res_bad = (max(res_good, r), res_bad) if good else (res_good, min(res_bad, r))
report(f'DM (4.21) solves DM\'s corrected (4.15)-(4.16) at their q only after q_(4.21) -> -(-q)^(1/2), n = 2..7 '
       f'[{res_good:.0e}; with q_(4.21) = q: {res_bad:.2f}]', res_good < 1e-12 and res_bad > 1e-2)

# ---------------------------------------------------------------- 6. hypothesis (ii) of thm-k-reflection
# Schur needs dim End_B(V_mu(y1) (x) V_nu(y2)) = 1, B acting by the coproduct; computed at generic rapidities.
y1c, y2c = np.exp(0.41 + 0.2j), np.exp(-0.23 + 0.1j)


def comm_dim(A, B, eps):
    gens = [coideal_gen_tensor(A, B, j, 1/q, eps[j]) for j in range(len(eps))]
    return len(intertwiners([(G, G) for G in gens])[0])


dims = []
for n in (2, 3):
    N = n + 1
    for eps in (np.zeros(N), EPS_STAR*np.ones(N), EPS_STAR*np.r_[-1, np.ones(n)], rng.normal(size=N)):
        for cA in (False, True):
            for cB in (False, True):
                dims.append(comm_dim(an_vector(n, q, y1c, cA), an_vector(n, q, y2c, cB), eps))
for j2 in (1, 2):
    for eps in ([0, 0], [0.3, -0.7]):
        dims.append(comm_dim(sl2_spin(1, q, y1c), sl2_spin(j2, q, y2c), eps))
report(f'dim End_B(V_mu (x) V_nu) = 1: a_2, a_3 with (mu, nu) in {{V, Vbar}}^2 at eps = 0, s = +-1 and generic eps '
       f'(where no K exists); a_1 spin 1/2 (x) 1/2 and 1/2 (x) 1 [{len(dims)} cases]', set(dims) == {1})

# ---------------------------------------------------------------- 7. one-dimensional boundary modules vs BPTT
# At y = exp(-i lambda nu_n), nu_n = eta/lambda - (2n+1) pi/(2 lambda), pi_y(b_0) and pi_y(b_1) of eq-sg-coideal-rep
# (eps_j from eq-bsg-eta-theta) have a common eigenvector: V(y) (x) C_eps has a one-dimensional B-submodule.  Its
# eigenvalues beta_j are eq-bsg-eta-theta with eta -> etabar = pi(lambda+1) - eta and s <-> sbar (theta -> -theta),
# i.e. the matrix part of BPTT's excited states eq-bsg-excited-k, for every n.  At y = exp(+i lambda nu_n) the
# one-dimensional submodule carries the complementary eigenvalues eps_j (q + q^-1) - beta_j.  lambda = 4.3, eta = 8.1:
# nu_0, nu_1, nu_2 all in the physical strip.
lam_, eta_, vth = 4.3, 8.1, 0.37
qx = q_phys(lam_); sx = np.sin(np.pi*lam_)
epsx = {1: -np.cos(eta_ - 1j*vth)/sx, 0: -np.cos(eta_ + 1j*vth)/sx}
etab = np.pi*(lam_ + 1) - eta_
pred = {1: -np.cos(etab + 1j*vth)/sx, 0: -np.cos(etab - 1j*vth)/sx}


def common_eig(y):
    """(residual, beta_1, beta_0) for the eigenvector of pi_y(b_1) that is closest to an eigenvector of pi_y(b_0)."""
    V = sl2_spin(1, qx, y)
    b1, b0 = coideal_gen(V, 1, 1/qx, epsx[1]), coideal_gen(V, 0, 1/qx, epsx[0])
    w, U = np.linalg.eig(b1)
    out = []
    for k in range(2):
        v = U[:, k]; beta0 = np.vdot(v, b0 @ v)/np.vdot(v, v)
        out.append((np.linalg.norm(b0 @ v - beta0*v)/np.linalg.norm(v), w[k], beta0))
    return min(out, key=lambda t: t[0])


ok, nus = True, []
for nn in range(3):
    nu = eta_/lam_ - (2*nn + 1)*np.pi/(2*lam_); nus.append(nu)
    r, be1, be0 = common_eig(np.exp(-1j*lam_*nu))
    ok &= r < 1e-12 and abs(be1 - pred[1]) < 1e-12 and abs(be0 - pred[0]) < 1e-12
    r, be1, be0 = common_eig(np.exp(1j*lam_*nu))
    trq = qx + 1/qx
    ok &= r < 1e-12 and abs(be1 - (epsx[1]*trq - pred[1])) < 1e-12 and abs(be0 - (epsx[0]*trq - pred[0])) < 1e-12
ok &= all(0 < nu < np.pi/2 for nu in nus) and common_eig(np.exp(0.3 + 0.2j))[0] > 1e-3
report(f'sine-Gordon: at y = exp(-i lambda nu_n), n = 0, 1, 2, V(y) has a 1-dim B-submodule with beta_(1,0) = '
       f'-cos(etabar +- i theta)/sin(pi lambda) (BPTT excited states), the same for all n; complementary values at '
       f'exp(+i lambda nu_n); none at generic y', ok)

# ---------------------------------------------------------------- 8. B = Lbar (K (x) 1) L to first order in 1/y
# (a) Leading terms of the universal R-matrix for the Part VI coproduct, in the principal gradation:
#     R = q^{sum G h (x) h} (1 + sum_j (q_j - q_j^-1) f_j (x) e_j + ...), a power series in y2/y1.  Compared with the exact
#     R-check: (pi_A(y1) (x) pi_B(y2))(R) = q^{(lam_A, lam_B)} P R-check (normalized to 1 on highest weights; the scalar
#     factor is 1 + O((y2/y1)^2)), at first order by central differences.  a_1 spin 1/2 (x) 1 and 1 (x) 1/2 separate
#     the Cartan placement; a_2 with all four (mu, nu) in {V, Vbar}^2.
hR = 2e-4
cR = q - 1/q


def perm(dA, dB):
    """P: A (x) B -> B (x) A."""
    P = np.zeros((dA*dB, dA*dB))
    for i in range(dA):
        for j in range(dB):
            P[j*dA + i, i*dB + j] = 1
    return P


def r_expansion(mkA, mkB, cart):
    """(M0, M1): (pi_A(y1) (x) pi_B(y2))(R) = M0 + t M1 + O(t^2), t = y2/y1."""
    vals = {}
    for sg in (1, -1):
        A, B = mkA(sg/np.sqrt(hR)), mkB(np.sqrt(hR))
        Rc = rcheck(A, B)[0]
        i0 = A['hw']*B['dim'] + B['hw']
        vals[sg] = cart[i0, i0]*perm(B['dim'], A['dim']) @ Rc      # P: B (x) A -> A (x) B
    return (vals[1] + vals[-1])/2, (vals[1] - vals[-1])/(2*hR)


def rel(a, b):
    return np.linalg.norm(a - b)/np.linalg.norm(b)


ok = True
cases = [(lambda y, jA=jA: sl2_spin(jA, q, y), lambda y, jB=jB: sl2_spin(jB, q, y),
          np.diag([q**(a_*b_/2) for a_ in range(jA, -jA - 1, -2) for b_ in range(jB, -jB - 1, -2)]), True)
         for jA, jB in ((1, 2), (2, 1))]
W2 = np.eye(3) - 1/3                                       # projected weights eps_i of the a_2 vector representation
for cA in (False, True):
    for cB in (False, True):
        WA, WB = (-W2 if cA else W2), (-W2 if cB else W2)
        cases.append((lambda y, cA=cA: an_vector(2, q, y, cA), lambda y, cB=cB: an_vector(2, q, y, cB),
                      np.diag([q**(WA[i] @ WB[j]) for i in range(3) for j in range(3)]), False))
for mkA, mkB, C, is_sl2 in cases:
    M0, M1 = r_expansion(mkA, mkB, C)
    A, B = mkA(1.0), mkB(1.0)
    nodes = nodes_of(A)
    FE = sum(kron(A[('f', j)], B[('e', j)]) for j in nodes)
    FEk = sum(kron(A[('f', j)] @ A[('k', j)], np.linalg.inv(B[('k', j)]) @ B[('e', j)]) for j in nodes)
    ok &= rel(M0, C) < 1e-6 and rel(M1, cR*C @ FE) < 1e-6 and rel(M1, cR*FEk @ C) < 1e-6
    ok &= rel(M1, -cR*C @ FE) > 1e-2
    if is_sl2:                                             # placement variants differ only on a_1 spin 1/2 (x) 1
        ok &= rel(M1, cR*FE @ C) > 1e-2 and rel(M1, cR*C @ FEk) > 1e-2
report('universal R for the Part VI coproduct = q^{G h(x)h}(1 + sum_j (q-q^-1) f_j (x) e_j + ...) at first order in y2/y1, '
       'against the exact R-check (a_1 spin 1/2 (x) 1, a_2 all {V,Vbar} pairs); sign and Cartan-placement variants fail', ok)


# (b) B-matrix of eq-b-matrix for the a_n vector representation at first order in u = 1/y, with the second tensor factor
#     kept abstract.  Elements coeff * X * K_lam with X in {1, e_j, f_j} and K_lam group-like (K_lam e_j K_lam^-1 =
#     q^{(lam, alpha_j)} e_j, K_{alpha_j} = k_j) are stored as {(X, lam): coeff}; lam in the eps basis.
#     L = C_V (1 + c u sum_j F_j (x) e_j),  Lbar = C_Vbar (1 + c u sum_j F_j (x) f_j),  C_V = diag K_{eps_i}, C_Vbar = diag K_{-eps_i},
#     K = kappa (1 + u K1 + O(u^2)) from eq-an-k-matrix (Laurent coefficients by FFT) or K = 1.
def b_first_order(n, s, qq):
    N = n + 1; c = qq - 1/qq; es = 1/np.sqrt((1 - qq)*(1 - 1/qq))
    alpha, Fm = {}, {}
    for j in range(N):
        a, b = (n, 0) if j == 0 else (j - 1, j)
        al = [0]*N; al[a] += 1; al[b] -= 1; alpha[j] = tuple(al)
        Fm[j] = np.zeros((N, N)); Fm[j][b, a] = 1
    eps_v = [tuple(int(k == i) for k in range(N)) for i in range(N)]
    pair = lambda u_, v_: sum(x_*y_ for x_, y_ in zip(u_, v_))
    wt = lambda X: alpha[X[1]] if X[0] == 'e' else tuple(-x_ for x_ in alpha[X[1]])
    z0 = (0,)*N

    def add(d, key, v):
        d[key] = d.get(key, 0) + v

    def mul(x_, y_):
        out = {}
        for (l1, m1), c1 in x_.items():
            for (l2, m2), c2 in y_.items():
                m = tuple(u1 + u2 for u1, u2 in zip(m1, m2))
                add(out, (l1, m) if l2 is None else (l2, m), c1*c2*(1 if l2 is None else qq**pair(m1, wt(l2))))
        return out

    def mat(f):
        return [[f(r_, s_) for s_ in range(N)] for r_ in range(N)]

    def _sum(ds):
        out = {}
        for d in ds:
            for k, v in d.items():
                add(out, k, v)
        return out

    def matmul(A, B):
        return mat(lambda r_, s_: _sum([mul(A[r_][t], B[t][s_]) for t in range(N)]))
    neg = lambda lam: tuple(-x_ for x_ in lam)
    C_V = mat(lambda r_, s_: {(None, eps_v[r_]): 1} if r_ == s_ else {})
    C_W = mat(lambda r_, s_: {(None, neg(eps_v[r_])): 1} if r_ == s_ else {})
    L1 = mat(lambda r_, s_: _sum([{(('e', j), eps_v[r_]): c*Fm[j][r_, s_]*qq**pair(eps_v[r_], alpha[j])}
                                  for j in range(N) if Fm[j][r_, s_]]))
    Lb1 = mat(lambda r_, s_: _sum([{(('f', j), neg(eps_v[r_])): c*Fm[j][r_, s_]*qq**pair(eps_v[r_], alpha[j])}
                                   for j in range(N) if Fm[j][r_, s_]]))
    if s is None:
        K1 = np.zeros((N, N)); epshat = [0]*N
    else:                                                  # K / X^(N/2) is a polynomial in u of degree <= N
        w = np.sqrt(-qq + 0j); Mf = 2*N + 2; us = 0.5*np.exp(2j*np.pi*np.arange(Mf)/Mf)
        P = [k_closed(n, 1/u_, s, qq)/(w/u_)**(N/2) for u_ in us]
        coef = lambda m: sum(P[k]*us[k]**(-m) for k in range(Mf))/Mf
        K0c, K1c = coef(0), coef(1)
        assert np.abs(K0c - K0c[0, 0]*np.eye(N)).max() < 1e-12*abs(K0c[0, 0])
        K1 = K1c/K0c[0, 0]; epshat = [s*es] + [es]*n
    Ks1 = mat(lambda r_, s_: {(None, z0): K1[r_, s_]} if abs(K1[r_, s_]) > 0 else {})
    B0 = matmul(C_W, C_V)
    B1 = [[_sum([x_, y_, z_]) for x_, y_, z_ in zip(*rows)] for rows in
          zip(matmul(Lb1, C_V), matmul(matmul(C_W, Ks1), C_V), matmul(C_W, L1))]
    want = mat(lambda r_, s_: _sum([{(('e', j), z0): c*Fm[j][r_, s_], (('f', j), alpha[j]): c*Fm[j][r_, s_]/qq,
                                     (None, alpha[j]): c*Fm[j][r_, s_]*epshat[j]} for j in range(N) if Fm[j][r_, s_]]))
    clean = lambda d: {k: v for k, v in d.items() if abs(v) > 1e-12}
    ok0 = all(clean(B0[r_][s_]) == ({(None, z0): 1} if r_ == s_ else {}) for r_ in range(N) for s_ in range(N))
    ok1 = all(clean(_sum([B1[r_][s_], {k: -v for k, v in want[r_][s_].items()}])) == {}
              for r_ in range(N) for s_ in range(N))
    return ok0 and ok1


ok = all(b_first_order(n, s, qq) for n in range(1, 7) for s in (None, 1, -1) for qq in (q, 0.7 + 0.3j))
report('B = Lbar (K(x)1) L, a_n vector rep (n = 1..6), K = 1 and eq-an-k-matrix (s = +-1), physical and generic q: '
       'B/kappa = 1 + (q-q^-1) y^-1 sum_j pi_V(f_j)|_{y=1} (x) (e_j + q^-1 f_j k_j + eps_j k_j) + O(y^-2), '
       'eps = 0 resp. (s eps_*, eps_*, ..., eps_*)', ok)

# ---------------------------------------------------------------- 9. sine-Gordon: eps near the Neumann point
# eq-bsg-eta-theta with cnj-bsg-uv-ir: cos((eta -+ i vartheta)/(lambda+1)) = mu_+-, mu_+- = (mu~/mu~_crit) e^{+-i beta phi_0/2}.
# So eps_1 = -cos((lambda+1) arccos mu_+)/sin(pi lambda), eps_0 the same with mu_-, and near mu~ = 0
# eps = eps_* - (lambda+1)/(2 sin(pi lambda/2)) mu_+- + O(mu^2); eps = 0 at phi_0 = 0, mu~/mu~_crit = cos(pi/(2(lambda+1))).
lamS = 1.37
sS = np.sin(np.pi*lamS)
eps_sg = lambda mu: -np.cos((lamS + 1)*np.arccos(mu + 0j))/sS
ok = True
for m_, ph in ((0.31, 0.4), (0.8, -1.1), (1.7, 0.25)):       # m_ = mu~/mu~_crit, ph = beta phi_0/2
    Am, Ap = (lamS + 1)*np.arccos(m_*np.exp(1j*ph) + 0j), (lamS + 1)*np.arccos(m_*np.exp(-1j*ph) + 0j)
    eta_s, vth_s = (Am + Ap)/2, (Ap - Am)/2j                   # Am = eta - i vartheta, Ap = eta + i vartheta
    ok &= abs(eta_s.imag) < 1e-12 and abs(vth_s.imag) < 1e-12
    e_, t_ = eta_s.real/(lamS + 1), vth_s.real/(lamS + 1)
    ok &= (abs(np.cos(e_)*np.cosh(t_) - m_*np.cos(ph)) < 1e-12 and abs(np.sin(e_)*np.sinh(t_) - m_*np.sin(ph)) < 1e-12)
    ok &= (abs(-np.cos(eta_s.real - 1j*vth_s.real)/sS - eps_sg(m_*np.exp(1j*ph))) < 1e-12 and
           abs(-np.cos(eta_s.real + 1j*vth_s.real)/sS - eps_sg(m_*np.exp(-1j*ph))) < 1e-12)
hS, zS = 1e-5, np.exp(0.7j)
dS = (eps_sg(hS*zS) - eps_sg(-hS*zS))/(2*hS*zS)
kapS = -(lamS + 1)/(2*np.sin(np.pi*lamS/2))
ok &= abs(eps_sg(0) - 1/(2*np.cos(np.pi*lamS/2))) < 1e-12 and abs(dS - kapS) < 1e-8
ok &= abs(eps_sg(np.cos(np.pi/(2*(lamS + 1))))) < 1e-12
report(f'sine-Gordon: eq-bsg-eta-theta + cnj-bsg-uv-ir give eps_(1,0) = -cos((lambda+1) arccos mu_+-)/sin(pi lambda), '
       f'= eps_* - (lambda+1)/(2 sin(pi lambda/2)) mu_+- + O(mu^2) near Neumann; eps = 0 at mu_+- = cos(pi/(2(lambda+1)))', ok)
