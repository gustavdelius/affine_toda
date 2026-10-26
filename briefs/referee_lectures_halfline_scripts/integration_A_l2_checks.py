#!/usr/bin/env python3
"""Integration agent A: independent checks of the findings on half-line Lecture 2.

Generators built from elementary matrices; no book or referee helper is imported.
"""
import itertools
import numpy as np
import sympy as sp


def rep(name, ok, extra=''):
    print(('PASS ' if ok else 'FAIL ') + name + (f'  [{extra}]' if extra else ''))


def e(n, i, j):
    m = np.zeros((n, n), complex); m[i, j] = 1; return m


# ---------------------------------------------------------------- generators
def an_gens(N):
    """a_{N-1}^(1) defining rep: node j (1..N-1) links v_j, v_{j+1}; node 0 links v_N -> v_1 (E_0 = e_{N,1})."""
    g = []
    E0 = e(N, N - 1, 0); g.append((E0 @ E0.T - E0.T @ E0, E0, E0.T))
    for j in range(N - 1):
        Ej = e(N, j, j + 1); g.append((Ej @ Ej.T - Ej.T @ Ej, Ej, Ej.T))
    return g


def node_system(gens, Cs, lam):
    d = gens[0][0].shape[0]; Id = np.eye(d)
    rows = []
    for (hj, Ej, Fj), Cj in zip(gens, Cs):
        Xm = -Cj*hj + lam*Ej - Fj/lam
        Xp = Cj*hj + lam*Ej - Fj/lam
        rows.append(np.kron(Id, Xm.T) - np.kron(Xp, Id))      # row-major vec(K): K Xm - Xp K
    return np.vstack(rows)


def smin(gens, Cs, lam):
    M = node_system(gens, Cs, lam)
    s = np.linalg.svd(M, compute_uv=False)
    return s[-1]/s[0]


lam = 1.27 + 0.31j

# ---------------------------------------------------------------- L2-3 a_1 matrix, invertibility
L, C0, C1 = sp.symbols('lambda C0 C1')
Ka1 = sp.Matrix([[L**2 - L**-2, -2*(C0*L - C1/L)], [-2*(C1*L - C0/L), L**2 - L**-2]])
h1 = sp.diag(1, -1); E1 = sp.Matrix([[0, 1], [0, 0]]); F1 = E1.T
nodes = [(-h1, F1, E1, C0), (h1, E1, F1, C1)]
ok = all(sp.simplify(Ka1*(-C*hh + L*EE - FF/L) - (C*hh + L*EE - FF/L)*Ka1) == sp.zeros(2) for hh, EE, FF, C in nodes)
rep('L2-3: eq-lec2-k-a1 solves both node equations for symbolic C0, C1', ok)
det = sp.expand(Ka1.det()*L**4)
rep('L2-3: L^4 det K is a polynomial in L with leading term L^8 for all C0, C1 (generically invertible)',
    sp.Poly(det, L).LC() == 1, f'{sp.factor(det)}')
rep('L2-3: at C0 = C1 = 0, lambda = 1 the displayed matrix is zero', Ka1.subs({C0: 0, C1: 0, L: 1}) == sp.zeros(2))

# ---------------------------------------------------------------- L2-4 mixed branches for a_2
g3 = an_gens(3)        # order: node 0, node 1, node 2
vals = [0.0, 1.0, -1.0, 0.43 - 0.21j, 1.7]
allowed, forbidden = [], []
for Cs in itertools.product(vals, repeat=3):
    s = smin(g3, Cs, lam)
    good = all(c == 0 for c in Cs) or all(abs(c**2 - 1) < 1e-12 for c in Cs)
    (allowed if good else forbidden).append(s)
rep('L2-4: a_2 node equations have a nonzero K only for C = 0 or all C_j^2 = 1 (125 patterns incl. mixed)',
    max(allowed) < 1e-12 and min(forbidden) > 1e-3, f'max allowed {max(allowed):.1e}, min forbidden {min(forbidden):.3f}')

# symbolic propagation: node 1 with C1^2 != 1 kills k13,k23,k31,k32; then node 2's spectator equations kill k12,k21
k = sp.Matrix(3, 3, sp.symbols('k11 k12 k13 k21 k22 k23 k31 k32 k33'))
Cs_sym = sp.symbols('c0 c1 c2')
Lsym = sp.symbols('lambda', nonzero=True)
gsym = []
for (hj, Ej, Fj) in g3:
    gsym.append((sp.Matrix(hj.real.astype(int)), sp.Matrix(Ej.real.astype(int)), sp.Matrix(Fj.real.astype(int))))
eqs = {}
for (hj, Ej, Fj), Cj, nm in zip(gsym, Cs_sym, ('0', '1', '2')):
    eqs[nm] = list(k*(-Cj*hj + Lsym*Ej - Fj/Lsym) - (Cj*hj + Lsym*Ej - Fj/Lsym)*k)
sub = {k[0, 2]: 0, k[1, 2]: 0, k[2, 0]: 0, k[2, 1]: 0}
# node-1 spectator system: column 3 rows 1,2 and row 3 columns 1,2
spec1 = [eqs['1'][0*3 + 2], eqs['1'][1*3 + 2], eqs['1'][2*3 + 0], eqs['1'][2*3 + 1]]
Mspec = sp.Matrix([[sp.diff(q_, v_) for v_ in (k[0, 2], k[1, 2], k[2, 0], k[2, 1])] for q_ in spec1])
rep('L2-4: node-1 spectator equations have determinant (1 - c1^2)^2 (only k13,k23,k31,k32 enter)',
    sp.factor(Mspec.det()) == sp.factor((1 - Cs_sym[1]**2)**2) and
    all(sp.diff(q_, v_) == 0 for q_ in spec1 for v_ in k if v_ not in (k[0, 2], k[1, 2], k[2, 0], k[2, 1])))
node2 = [sp.simplify(q_.subs(sub)) for q_ in eqs['2']]
sol = sp.solve(node2, [k[0, 1], k[1, 0]], dict=True)
rep('L2-4: with k13=k23=k31=k32=0, node 2 forces k12 = k21 = 0 for every c2',
    sol == [] or all(s_.get(k[0, 1], None) == 0 and s_.get(k[1, 0], None) == 0 for s_ in sol)
    or sp.solve([q_ for q_ in node2 if q_ != 0], [k[0, 1], k[1, 0]]) == {k[0, 1]: 0, k[1, 0]: 0},
    f'{[q_ for q_ in node2 if q_ != 0]}')
diag_sub = dict(sub); diag_sub.update({k[0, 1]: 0, k[1, 0]: 0})
rest = [sp.simplify(q_.subs(diag_sub)) for nm in '012' for q_ in eqs[nm]]
rest = [q_ for q_ in rest if q_ != 0]
print('     L2-4: remaining equations for diagonal K:', rest)

# ---------------------------------------------------------------- L2-5 figure slices
s_or = [smin(g3, (c, 1.0, -1.0), 1.25 + 0.2j) for c in (0.0, 1.0, -1.0)]
s_gr = [smin(g3, (c, c, c), 1.25 + 0.2j) for c in (0.0, 1.0, -1.0)]
rep('L2-5: orange slice (C,+1,-1) has a solution at C = +-1 but not at C = 0', s_or[0] > 1e-3 and max(s_or[1:]) < 1e-12,
    f'{s_or[0]:.3f}')
rep('L2-5: green slice (C,C,C) has solutions at C = 0, +-1', max(s_gr) < 1e-12)

# ---------------------------------------------------------------- L2-6 Dirichlet: Cartan condition
for N in (2, 3, 4):
    Hs = [np.diag([1.0 if a == i else (-1.0 if a == i + 1 else 0.0) for a in range(N)]) for i in range(N - 1)]
    Id = np.eye(N)
    Mc = np.vstack([np.kron(Id, H.T) + np.kron(H, Id) for H in Hs])     # {K, H} = 0
    nul = Mc.shape[1] - np.linalg.matrix_rank(Mc)
    Mcomm = np.vstack([np.kron(Id, H.T) - np.kron(H, Id) for H in Hs])  # [K, H] = 0
    nulc = Mcomm.shape[1] - np.linalg.matrix_rank(Mcomm)
    print(f'     L2-6: sl_{N}: dim{{K: {{K,H}}=0 for all H}} = {nul}; dim{{K: [K,H]=0}} = {nulc}')
    if N == 2:
        rep('L2-6: sl_2: anticommutator leaves the antidiagonal K (dimension 2)', nul == 2)
    else:
        rep(f'L2-6: sl_{N}: anticommutator {{K,H}} = 0 for all Cartan H forces K = 0', nul == 0)
# derivation of the Cartan term: K a_t(lam) + a_t(1/lam)^T K, Cartan part (beta/2) dphi.H on both sides
H = np.diag([0.3, -1.1, 0.8]); K = np.random.default_rng(3).normal(size=(3, 3))
rep('L2-6: transpose keeps H, so the Cartan part of K a_t + a_t(1/lam)^T K is K H + H K',
    np.allclose(K @ H + H.T @ K, K @ H + H @ K))

# ---------------------------------------------------------------- L2-8 d_4 vector representation
def d4_vector():
    """so(8) vector rep, basis e_1..e_4, e_-1..e_-4 (indices 0..3, 4..7).  E_{e_i - e_j} = e_{i,j} - e_{-j,-i},
    E_{e_i + e_j} = e_{i,-j} - e_{j,-i}; E_{-e_i - e_j} = transpose.  Simple roots alpha_0 = -e1-e2,
    alpha_1 = e1-e2, alpha_2 = e2-e3, alpha_3 = e3-e4, alpha_4 = e3+e4; Kac labels (1,1,2,1,1)."""
    n = 8
    p = lambda i: i - 1; m = lambda i: 3 + i
    def Emin(i, j):   # e_i - e_j
        return e(n, p(i), p(j)) - e(n, m(j), m(i))
    def Eplus(i, j):  # e_i + e_j
        return e(n, p(i), m(j)) - e(n, p(j), m(i))
    Es = [Eplus(1, 2).T, Emin(1, 2), Emin(2, 3), Emin(3, 4), Eplus(3, 4)]
    gens = []
    for Ej in Es:
        Fj = Ej.T
        hj = Ej @ Fj - Fj @ Ej
        gens.append((hj, Ej, Fj))
    return gens


g8 = d4_vector()
okrel = all(np.allclose(hj @ Ej - Ej @ hj, 2*Ej) for hj, Ej, Fj in g8)
rep('L2-8: d_4 vector Chevalley generators: [h_j, E_j] = 2 E_j', okrel)
for j, (hj, Ej, Fj) in enumerate(g8):
    act = sorted(set(np.nonzero(Ej)[0]) | set(np.nonzero(Ej)[1]))
    print(f'     L2-8: node {j}: E_j acts on basis vectors {act} (two doublets), spectators {sorted(set(range(8)) - set(act))}')
nj = [1, 1, 2, 1, 1]
lam8 = 1.317 + 0.07j
norm_s = [smin(g8, list(sg), lam8) for sg in itertools.product((1, -1), repeat=5)]
alt_s = [smin(g8, [c/np.sqrt(n_) for c, n_ in zip(sg, nj)], lam8) for sg in itertools.product((1, -1), repeat=5)]
alt2_s = [smin(g8, [c/np.sqrt(n_) for c, n_ in zip((sg[0], sg[1], sg[2]*np.sqrt(2), sg[3], sg[4]), nj)], lam8)
          for sg in itertools.product((1, -1), repeat=5)]
rep('L2-8: normalized potential: K exists for all 32 sign patterns C_j = +-1 (d_4 vector)', max(norm_s) < 1e-10,
    f'{max(norm_s):.1e}')
rep('L2-8: unnormalized potential (C_j/sqrt n_j in the node equation): no K for any |C_j| = 1', min(alt_s) > 1e-3,
    f'{min(alt_s):.3f}')
rep('L2-8: unnormalized potential with C_2^2 = 2 (others +-1): K exists', max(alt2_s) < 1e-10, f'{max(alt2_s):.1e}')
alt_node2 = [smin(g8, [1, 1, c/np.sqrt(2), 1, 1], lam8) for c in (1.0, 1.2)]
rep('L2-8: equation stays node-local: only node 2 coefficient changes (C_2/sqrt2); kernel iff (C_2/sqrt2)^2 = 1',
    min(alt_node2) > 1e-3 and smin(g8, [1, 1, 1.0, 1, 1], lam8) < 1e-10)

# ---------------------------------------------------------------- L2-10 coweight flips and inversion (a_2)
def Ka2(lm, c0, c1, c2):
    s = c0*c1*c2; d = s*lm**2 - 1/lm
    return np.array([[d, 2*c1, -2*c1*c2*lm], [-2*c0*c2*lm, d, 2*c2], [2*c0, -2*c0*c1*lm, d]], complex)


def solves(K, Cs, lm):
    return max(np.abs(K @ (-C*hj + lm*Ej - Fj/lm) - (C*hj + lm*Ej - Fj/lm) @ K).max() for (hj, Ej, Fj), C in zip(g3, Cs))


pats = list(itertools.product((1, -1), repeat=3))
rep('L2-10: eq-lec2-k-a2 solves the node equations for all 8 patterns', max(solves(Ka2(lam, *p), p, lam) for p in pats) < 1e-12)
pairs = {0: (2, 0), 1: (0, 1), 2: (1, 2)}     # node -> (a, b) linked basis indices
orbit_ok = True
for dvec in itertools.product((1, -1), repeat=3):
    D = np.diag(dvec)
    t = [dvec[pairs[j][0]]*dvec[pairs[j][1]] for j in range(3)]
    for p in pats:
        newp = tuple(t[j]*p[j] for j in range(3))
        Kp = D @ Ka2(lam, *p) @ D
        orbit_ok &= solves(Kp, newp, lam) < 1e-12          # same lambda
rep('L2-10: D K D^{-1} solves the node equations with C_j -> t_j C_j, t_j = d_a d_b, at the same lambda', orbit_ok)
tsets = {tuple(dvec[pairs[j][0]]*dvec[pairs[j][1]] for j in range(3)) for dvec in itertools.product((1, -1), repeat=3)}
rep('L2-10: the flips t have product +1 (preserve s) and form 4 elements', all(np.prod(tt) == 1 for tt in tsets) and len(tsets) == 4)
orbits = []
for p in pats:
    o = {tuple(tt[j]*p[j] for j in range(3)) for tt in tsets}
    if o not in orbits:
        orbits.append(o)
rep('L2-10: coweight flips alone give two classes of four (s = +1 and s = -1)',
    len(orbits) == 2 and all(len(o) == 4 for o in orbits))
inv_ok = all(solves(np.linalg.inv(Ka2(lam, *p)), tuple(-c for c in p), lam) < 1e-12 for p in pats)
rep('L2-10: K^{-1} solves with C_j -> -C_j (inversion flips s), so all eight form one orbit', inv_ok)
Kinv = np.linalg.inv(Ka2(lam, 1, 1, 1)); Kmmm = Ka2(lam, -1, -1, -1)
cfit = np.vdot(Kmmm.ravel(), Kinv.ravel())/np.vdot(Kmmm.ravel(), Kmmm.ravel())
rep('L2-10: K_{(+++)}^{-1} is proportional to K_{(---)}', np.abs(Kinv - cfit*Kmmm).max() < 1e-12)
splus = [p for p in pats if np.prod(p) == 1]
rep('L2-10: the s = +1 class is (+++) plus three patterns with two minus signs', sorted(splus) == sorted([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]))

# ---------------------------------------------------------------- L2-2 K_- = 1 at the vacuum
rng = np.random.default_rng(11)
cj = [1.0, 1.0, 1.0]
at = lambda lm: 0.5*sum(c*(lm*Ej - Fj/lm) for c, (hj, Ej, Fj) in zip(cj, g3))
rep('L2-2: at the vacuum a_t(lam)^T + a_t(1/lam) = 0, so K_- = 1 satisfies the x_- condition (sigma = T)',
    np.abs(at(lam).T + at(1/lam)).max() < 1e-14)
# imaginary-coupling vacua: e^{i beta alpha_j . phi/2} = +-1 per node; still cancels
sg = [1, -1, -1]
at2 = lambda lm: 0.5*sum(s_*(lm*Ej - Fj/lm) for s_, (hj, Ej, Fj) in zip(sg, g3))
rep('L2-2: also at a vacuum with node signs (imaginary coupling)', np.abs(at2(lam).T + at2(1/lam)).max() < 1e-14)

# ---------------------------------------------------------------- L2-7 coweight extension of the soliton-preserving wall
beta, mm = 0.9, 1.3
roots = [np.array([-1, 0, 1.]), np.array([1, -1, 0.]), np.array([0, 1, -1.])]   # alpha_0, alpha_1, alpha_2 in eps basis
Hvec = [np.diag([1., 0, 0]), np.diag([0, 1., 0]), np.diag([0, 0, 1.])]          # eps-basis Cartan (traceless part irrelevant)
nu = np.array([2/3, -1/3, -1/3])           # fundamental coweight omega_1 (alpha_j . nu in Z)
ok_all = True
for trial in range(3):
    ImPhi = rng.normal(size=3); ImPhi -= ImPhi.mean()
    lattice = np.array([1, 0, -1.])*rng.integers(-2, 3)           # alpha_j . lattice in Z: coroot lattice element
    RePhi = (np.pi/beta)*nu + (2*np.pi/beta)*lattice
    phi = RePhi + 1j*ImPhi
    lr = 1.37 + 0.1*trial                                          # real lambda
    w = [np.exp(1j*beta*np.dot(a, phi)/2) for a in roots]
    steps = lambda lm, ww: 0.5*mm*sum(wj*(lm*Ej - Fj/lm) for wj, (hj, Ej, Fj) in zip(ww, g3))
    Kp = np.diag(np.exp(1j*np.pi*nu))                              # e^{i pi nu.H}
    resid = Kp @ steps(lr, w) + steps(1/lr, w).conj().T @ Kp      # Neumann on Im phi removes the Cartan part
    ok_all &= np.abs(resid).max() < 1e-12
    # without K_+: fails
    resid1 = steps(lr, w) + steps(1/lr, w).conj().T
    ok_all &= np.abs(resid1).max() > 1e-3
rep('L2-7: K_+ = e^{i pi nu.H} cancels the step terms when Re phi(0) in (pi/beta) nu + (2 pi/beta) coweights; K_+ = 1 does not',
    ok_all)
