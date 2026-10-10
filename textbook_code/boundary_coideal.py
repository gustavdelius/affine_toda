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
 3. a_n^(1) vector solitons, n = 2, 3, 4: K: V(y) -> Vbar(1/y) exists iff all eps_j = 0 (then K = 1)
    or all eps_j^2 = 1/((1-q)(1-q^-1)) = 1/(4 cos^2(pi omega/2)); not for |eps_j| = 1 (DM (4.21)).
    Closed form (corrected DM (4.21)) and the substitution that maps DM (4.21) onto it; gauge
    equivalence of sign patterns; twisted reflection equation for all four (mu, nu); failure of the
    printed K at -q (where the antisoliton is Vbar(-y) and eps = 0 has a solution only for odd n).
 4. Identification of the solutions: the non-diagonal s=+1 solution is Gandenberger's K^+ for a_2
    (hep-th/9806003 (3.7), his Neumann candidate), s=-1 is his K^-; eps = 0 gives DG (3.7) (K = 1).
    Boundary crossing-unitarity (GZ (3.35)) in book form k(theta) = f R_{Vbar Vbar}(2theta) k(-theta),
    k = sum_ab K_ba(i pi/2 - theta) v_a (x) v_b, holds for the matrix part (n = 2, 3).
Runtime: about a minute.
"""
import numpy as np
import sympy as sp
from boundary_common import (report, q_phys, an_vector, sl2_spin, dn_vector, check_relations, cartan_a,
                             cartan_d, coideal_gen, coideal_gen_tensor, cop, k_intertwiner, rcheck,
                             re_twisted, re_preserving, dual, kron)

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
for n in (2, 3, 4):
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
