"""Independent checks by integration agent C for Lectures 5 and 6 of the half-line course.

Written from the lecture and book definitions; imports no referee script and no lecture script.
Each check prints one PASS/FAIL line.
"""
import itertools
import numpy as np
import mpmath as mp
from scipy.integrate import quad

RESULTS = []


def report(name, ok, info=""):
    RESULTS.append(ok)
    print(("PASS " if ok else "FAIL ") + name + (f"  [{info}]" if info else ""))


# ---------------------------------------------------------------- a_2 representations
def elem(N, a, b):
    M = np.zeros((N, N), dtype=complex)
    M[a, b] = 1.0
    return M


def an_rep(q, y, N, conj=False):
    """V(y) (conj=False) or Vbar(y) (conj=True) of U_q(a_{N-1}^(1)) in the principal gradation.
    Node j (1..N-1) has alpha_j = eps_j - eps_{j+1}; node 0 has alpha_0 = eps_N - eps_1 (0-based indices)."""
    rep = {"dim": N}
    for j in range(N):
        a, b = ((N - 1, 0) if j == 0 else (j - 1, j))
        E, F = elem(N, a, b), elem(N, b, a)
        H = E @ F - F @ E
        k = np.diag(q ** np.diag(H))
        if not conj:
            rep["e", j], rep["f", j], rep["k", j] = y * E, F / y, k
        else:
            rep["e", j], rep["f", j], rep["k", j] = y * F, E / y, np.linalg.inv(k)
    rep["nodes"] = N
    return rep


def tensor(A, B):
    IA, IB = np.eye(A["dim"]), np.eye(B["dim"])
    T = {"dim": A["dim"] * B["dim"], "nodes": A["nodes"]}
    for j in range(A["nodes"]):
        T["e", j] = np.kron(A["e", j], IB) + np.kron(A["k", j], B["e", j])
        T["f", j] = np.kron(A["f", j], np.linalg.inv(B["k", j])) + np.kron(IA, B["f", j])
        T["k", j] = np.kron(A["k", j], B["k", j])
    return T


def bgen(rep, j, eps, q):
    return rep["e", j] + (1 / q) * rep["f", j] @ rep["k", j] + eps * rep["k", j]


def intertwiner_nullity(q, y, eps, N, tol=1e-9):
    Vin, Vout = an_rep(q, y, N), an_rep(q, 1 / y, N, conj=True)
    I = np.eye(N)
    rows = [np.kron(I, bgen(Vin, j, eps[j], q).T) - np.kron(bgen(Vout, j, eps[j], q), I) for j in range(N)]
    s = np.linalg.svd(np.vstack(rows), compute_uv=False)
    return s[-1], int(np.sum(s < tol))


omega = np.sqrt(2.0)
q = -np.exp(-1j * np.pi * omega)
w = np.exp(-0.5j * np.pi * omega)
eps_star = 1 / (2 * np.cos(np.pi * omega / 2))
y0 = np.exp(omega * (0.37 + 0.21j))

# V0: eps_* identities used in Lecture 5
report("eps_*^2 = 1/((1-q)(1-1/q)), 1-q = w/eps_*, eps_*(1/q-q) = w-1/w",
       abs(eps_star ** 2 * (1 - q) * (1 - 1 / q) - 1) < 1e-14 and abs(1 - q - w / eps_star) < 1e-14
       and abs(eps_star * (1 / q - q) - (w - 1 / w)) < 1e-14)

# V1 (L5-2): a_2 dichotomy over 64 patterns, incl. mixed ones
vals = [0.0, eps_star, -eps_star, (np.e / 4) * eps_star]
allowed, bad = [], []
for pat in itertools.product(range(4), repeat=3):
    eps = [vals[i] for i in pat]
    smin, nul = intertwiner_nullity(q, y0, eps, 3)
    expected = all(i == 0 for i in pat) or all(i in (1, 2) for i in pat)
    (allowed if nul > 0 else bad).append(pat)
    if (nul > 0) != expected:
        print("  mismatch", pat, smin)
ok = len(allowed) == 9 and all(all(i == 0 for i in p) or all(i in (1, 2) for i in p) for p in allowed)
report("a_2: K exists for exactly 9 of 64 node patterns (all zero, or all +-eps_*)", ok, f"{len(allowed)} allowed")

# V1b: the propagation argument: if node 1 has eps^2 != eps_*^2, the column-3/row-3 entries vanish and
# then node-2 and node-0 spectator equations force the other off-diagonals to vanish whatever their eps.
# (checked implicitly above; here the hand argument: each 2x2 spectator system has y^{+-1} off the diagonal)
M2 = lambda e, y: np.array([[e * (1 / q - 1), y], [1 / y, e * (q - 1)]])
okp = all(abs(M2(e, y0)[1, 0]) > 0 and abs(M2(e, y0)[0, 1]) > 0 for e in vals)
report("spectator 2x2 systems have nonzero off-diagonal y, 1/y (one vanishing unknown forces the other)", okp)

# V2 (L5-4): d_ij formula and a_2 spectral argument
W = tensor(tensor(an_rep(q, 0.7 + 0.4j, 3), an_rep(q, 1.3 - 0.2j, 3, conj=True)), an_rep(q, 0.5 + 0.9j, 3))
err = 0.0
rng = np.random.default_rng(7)
for (i, j) in [(1, 2), (2, 1), (0, 1), (2, 0)]:
    ei, ej = rng.normal(size=2) + 1j * rng.normal(size=2)
    bi, bj = bgen(W, i, ei, q), bgen(W, j, ej, q)
    F = bi @ bi @ bj - (q + 1 / q) * bi @ bj @ bi + bj @ bi @ bi
    d = ej * ((1 - q) * (1 - 1 / q) * ei ** 2 - 1)
    err = max(err, np.linalg.norm(F - bj - d * W["k", i] @ W["k", i] @ W["k", j]))
report("F_ij(b_i,b_j) - b_j = d_ij k_i^2 k_j with d_ij = eps_j((1-q)(1-1/q)eps_i^2 - 1) (c=1/q)", err < 1e-12, f"{err:.1e}")

V, Vb = an_rep(q, 1.0, 3), an_rep(q, 1.0, 3, conj=True)
okspec, sylv = True, []
for (i, j) in [(1, 2), (2, 1), (2, 0), (0, 2), (0, 1), (1, 0)]:
    kv = V["k", i] @ V["k", i] @ V["k", j]
    kvb = Vb["k", i] @ Vb["k", i] @ Vb["k", j]
    ev, evb = np.diag(kv), np.diag(kvb)
    gap = min(abs(a - b) for a in ev for b in evb)
    S = np.kron(np.eye(3), kv.T) - np.kron(kvb, np.eye(3))  # K kv - kvb K
    sylv.append(np.linalg.svd(S, compute_uv=False)[-1])
    okspec &= gap > 0.1
report("a_2: spectra of k_i^2 k_j on V and Vbar are disjoint, so K kappa_V = kappa_Vbar K forces K=0",
       okspec and min(sylv) > 0.1, f"min sigma={min(sylv):.2f}")
spec12 = np.round(np.diag(V["k", 1] @ V["k", 1] @ V["k", 2]) / np.array([q ** 2, 1 / q, 1 / q]), 12)
report("a_2: spectrum of k_1^2 k_2 on V is {q^2, 1/q, 1/q}, on Vbar its inverse", np.allclose(spec12, 1))
# a_3: the same operator is singular (eigenvalue 1 on the fourth basis vector), so the simple argument is a_2-specific
V4, V4b = an_rep(q, 1.0, 4), an_rep(q, 1.0, 4, conj=True)
kv4 = V4["k", 1] @ V4["k", 1] @ V4["k", 2]
kv4b = V4b["k", 1] @ V4b["k", 1] @ V4b["k", 2]
S4 = np.kron(np.eye(4), kv4.T) - np.kron(kv4b, np.eye(4))
report("a_3: the same spectral operator is singular (simple argument does not extend)", np.linalg.svd(S4, compute_uv=False)[-1] < 1e-12)

# V3 (L5-5): q-Dolan-Grady relation for a_1^(1) coideal, arbitrary eps
W1 = tensor(tensor(an_rep(q, 0.8 + 0.3j, 2), an_rep(q, 1.1 - 0.5j, 2)), an_rep(q, 0.6 + 0.2j, 2))
three = q ** 2 + 1 + q ** -2
rho = (q + 1 / q) ** 2
errdg = 0.0
for e0, e1 in [(0.31 - 0.72j, 1.17 + 0.4j), (eps_star, 0.0), (0.0, 0.0)]:
    b = {0: bgen(W1, 0, e0, q), 1: bgen(W1, 1, e1, q)}
    for i, j in [(0, 1), (1, 0)]:
        bi, bj = b[i], b[j]
        lhs = bi @ bi @ bi @ bj - three * bi @ bi @ bj @ bi + three * bi @ bj @ bi @ bi - bj @ bi @ bi @ bi
        errdg = max(errdg, np.linalg.norm(lhs - rho * (bi @ bj - bj @ bi)) / np.linalg.norm(lhs))
report("a_1^(1): q-Dolan-Grady b_i^3b_j-[3]b_i^2b_jb_i+[3]b_ib_jb_i^2-b_jb_i^3 = (q+1/q)^2[b_i,b_j], any eps", errdg < 1e-12, f"{errdg:.1e}")

# V4 (L5-7): character twist (id x chi) Delta(b_j) = b_j + beta_j k_j, on a representation
A = an_rep(q, 0.9 + 0.1j, 3)
# Delta(b) = (e + q^-1 f k) x 1 + k x b ; apply chi to the second factor: b -> beta
errchi = 0.0
for j in range(3):
    beta = eps_star * (-1) ** j
    twisted = (A["e", j] + (1 / q) * A["f", j] @ A["k", j]) + A["k", j] * beta
    errchi = max(errchi, np.linalg.norm(twisted - bgen(A, j, beta, q)))
# and the image satisfies F_ij(b_i,b_j) = b_j at eps = +-eps_* (q-Onsager relation)
errons = 0.0
for (i, j) in [(1, 2), (2, 0), (0, 1)]:
    bi, bj = bgen(W, i, eps_star, q), bgen(W, j, -eps_star, q)
    errons = max(errons, np.linalg.norm(bi @ bi @ bj - (q + 1 / q) * bi @ bj @ bi + bj @ bi @ bi - bj))
report("character twist maps b_j(0) to b_j(beta); at beta=+-eps_* F_ij(b_i,b_j)=b_j", errchi < 1e-14 and errons < 1e-12)


# V5 (L5-10, L5-11): Toeplitz matrix, Exercise 5.1 entries and Exercise 5.3 eigenvalues
def toeplitz(y, s, N=3):
    X = w * y
    K = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            if i == j:
                K[i, j] = (w * X ** (-N / 2) - s / w * X ** (N / 2)) / (eps_star * (1 / q - q))
            elif j > i:
                K[i, j] = X ** (j - i - N / 2)
            else:
                K[i, j] = s * X ** (j - i + N / 2)
    return K


th = 0.43 + 0.17j
y = np.exp(omega * th)
X = np.exp(omega * (th - 0.5j * np.pi))
ok5 = abs(w * y - X) < 1e-12
for s in (1, -1):
    K = toeplitz(y, s)
    eps = [s * eps_star, eps_star, eps_star]
    Vin, Vout = an_rep(q, y, 3), an_rep(q, 1 / y, 3, conj=True)
    res = max(np.linalg.norm(K @ bgen(Vin, j, eps[j], q) - bgen(Vout, j, eps[j], q) @ K) for j in range(3))
    ok5 &= res < 1e-12
    # Exercise 5.1 part 1 entries (0-based indices)
    ok5 &= abs(K[2, 0] - s * X ** -0.5) < 1e-12 and abs(K[1, 0] - s * X ** 0.5) < 1e-12 and abs(K[2, 1] - s * X ** 0.5) < 1e-12
    # node 0 spectator (c=2 -> index 1): column pair (K32,K12), row pair (K23,K21): relations
    ok5 &= abs(K[2, 1] - s * X * K[0, 1]) < 1e-12 and abs(K[1, 0] - s * X * K[1, 2]) < 1e-12
    # node 1 spectator row 3: K31 = K32 / X
    ok5 &= abs(K[2, 0] - K[2, 1] / X) < 1e-12
    # (1,1) entry from node 1
    K11 = (s * y * X ** 0.5 - X ** -0.5 / y) / (eps_star * (q - 1 / q))
    ok5 &= abs(K[0, 0] - K11) < 1e-12
report("Exercise 5.1: K31=sX^-1/2, K21=K32=sX^1/2 (node 0 gives K32,K21; K31 needs node 1 or 2), K11 formula", ok5)

ok53 = True
for s in (1, -1):
    Z = np.array([[0, 1, 0], [0, 0, 1], [s, 0, 0]], dtype=complex)
    ok53 &= np.allclose(Z.T @ Z, np.eye(3)) and np.allclose(np.linalg.matrix_power(Z, 3), s * np.eye(3))
    ok53 &= (s == 1) or not np.allclose(Z.T, s * np.linalg.inv(Z))  # the printed s Z^-1 is wrong for s=-1
    K = toeplitz(y, s)
    kap = K[0, 0]
    ok53 &= np.allclose(K, kap * np.eye(3) + X ** -0.5 * Z + X ** 0.5 * Z @ Z)
    for zeta in np.roots([1, 0, 0, -s]):
        lam = kap + X ** -0.5 * zeta + X ** 0.5 * zeta ** 2
        formula = (1 - s * X ** 3) * (X * zeta + q) / (w * X ** 1.5 * (w - 1 / w) * (X * zeta - 1))
        ok53 &= abs(lam - formula) < 1e-11 * abs(lam)
    Km = toeplitz(np.exp(-omega * th), s)
    T = 1.5 * omega
    U = (s * np.cosh(2 * T * th) - np.cos(np.pi * T)) / (2 * np.sin(np.pi * omega / 2) ** 2)
    ok53 &= np.allclose(Km.T @ K, U * np.eye(3))
report("Exercise 5.3: Z^T Z = 1 (not Z^T = sZ^-1), eigenvalue kappa(zeta), unitarity U_s", ok53)

# V6 (L5-9): conjectured UV-IR curve
lam = 2 * np.sqrt(3.0)
f = lambda mu: -np.cos((lam + 1) * np.arccos(mu)) / np.sin(np.pi * lam)
es = 1 / (2 * np.cos(np.pi * lam / 2))
h_ = 1e-6
slope = (f(h_) - f(-h_)) / (2 * h_)
mu0 = np.cos(np.pi / (2 * (lam + 1)))
report("UV-IR: f(0)=eps_*, f'(0)=-(lam+1)/(2 sin(pi lam/2)), zero at mu0=cos(pi/(2(lam+1))) not at mu=1",
       abs(f(0) - es) < 1e-12 and abs(slope + (lam + 1) / (2 * np.sin(np.pi * lam / 2))) < 1e-6
       and abs(f(mu0)) < 1e-12 and abs(f(1.0) + 1 / np.sin(np.pi * lam)) < 1e-12,
       f"mu0={mu0:.4f}")


# V7 (L6-4): crossed boundary vector of K=1 is sum_i vbar_i x vbar_i, fixed by normalized R_{Vbar Vbar}
def rcheck(A, B):
    d = A["dim"] * B["dim"]
    AB, BA = tensor(A, B), tensor(B, A)
    rows = []
    for j in range(A["nodes"]):
        for g in ("e", "f", "k"):
            rows.append(np.kron(np.eye(d), AB[g, j].T) - np.kron(BA[g, j], np.eye(d)))
    _, s, vh = np.linalg.svd(np.vstack(rows))
    R = vh[-1].conj().reshape(d, d)
    return R / R[0, 0], s[-1]


R, smin = rcheck(an_rep(q, np.exp(2 * omega * 0.31j + 0.2), 3, conj=True), an_rep(q, 1.0, 3, conj=True))
kvec = np.eye(3).reshape(9)
report("R_{Vbar Vbar} (normalized on |11>) fixes sum_i vbar_i x vbar_i (K=1 crossing vector)",
       np.linalg.norm(R @ kvec - kvec) < 1e-12 and smin < 1e-12)


# ---------------------------------------------------------------- blocks
def blk(x, theta, h):
    return mp.sinh(theta / 2 + 1j * mp.pi * x / (2 * h)) / mp.sinh(theta / 2 - 1j * mp.pi * x / (2 * h))


def blkI(a, mu, h, om):
    return mp.sin(mp.pi * (mu + a) / (h * om)) / mp.sin(mp.pi * (mu - a) / (h * om))


mp.mp.dps = 40
# V8 (L6-6): continuation of the three blocks of K_{B1} and the period 2h
h = 3
n = h - 1
B = mp.mpf(1) / mp.sqrt(5)
th = mp.mpc(0.63, 0.11)
om = -2 / B
T = h * om / 2
mu = -1j * T * th / mp.pi
KB1 = blkI(-n * om / 2, mu, h, om) * blkI(mp.mpf(1) / 2 - h * om / 2, mu, h, om) * blkI(-om / 2 - mp.mpf(1) / 2, mu, h, om)
K1 = blk(h - B / 2, th, h) / (blk(h - 1, th, h) * blk(1 - B / 2, th, h))
okc = abs(KB1 - K1) < mp.mpf(10) ** -30
okc &= abs(blk(-h - B / 2, th, h) - blk(h - B / 2, th, h)) < mp.mpf(10) ** -30
okc &= abs(blkI(mp.mpf(1) / 2 - h * om / 2, mu, h, om) - blk(-h - B / 2, th, h)) < mp.mpf(10) ** -30
report("continuation: (a)_I -> (2a/om) at om=-2/B; middle label -h-B/2 = h-B/2 by period 2h; K_B1 -> K_1", okc)


# V9 (L6-7): CDD factor P_x
def Px(mu, x, h, om):
    T = h * om / 2
    return blkI(T / 2 - x, mu, h, om) * blkI(T / 2 + x, mu, h, om)


okcdd = True
boot = {}
for hh in (3, 4, 5, 6):
    nn = hh - 1
    omg = mp.sqrt(2) + mp.mpf(1) / 7
    T = hh * omg / 2
    for xn, x in (("om/2", omg / 2), ("(om+1)/2", (omg + 1) / 2)):
        m1 = mp.mpc(0.23, -0.17)
        okcdd &= abs(Px(m1, x, hh, omg) * Px(-m1, x, hh, omg) - 1) < 1e-25
        okcdd &= abs(Px(T / 2 - m1, x, hh, omg) - Px(T / 2 + m1, x, hh, omg)) < 1e-25
        prod = mp.mpf(1)
        for k in range(1, nn + 1):
            prod *= Px(m1 + T - k * omg, x, hh, omg)
        boot[nn, xn] = abs(Px(m1, x, hh, omg) - prod) < 1e-20
report("P_{om/2} and P_{(om+1)/2}: unitary and crossing-symmetric (n=2..5)", okcdd)
# DG's soliton-bootstrap condition sigma(mu) = prod_k sigma(mu + T - k om): holds for odd n, fails for even n
report("soliton-bootstrap condition for P_x holds for n=3,5 and fails for n=2,4 (both x) [book checks only a_3]",
       all(boot[n_, x_] == (n_ % 2 == 1) for (n_, x_) in boot), str(boot))


def D(chi, th, h, B):
    return (blk(h + chi + B / 2, th, h) * blk(chi - B / 2, th, h) * blk(h - chi + B / 2, th, h)
            * blk(-chi - B / 2, th, h))


okD = True
for hh in (3, 4):
    for x_of in (lambda o: o / 2, lambda o: (o + 1) / 2):
        Bv = mp.mpf(1) / mp.sqrt(7)
        om = -2 / Bv
        T = hh * om / 2
        thv = mp.mpc(0.71, 0.13)
        muv = -1j * T * thv / mp.pi
        s_ = T / 2 - mp.mpf(1) / 2
        x = x_of(om)
        lhs = Px(muv - s_, x, hh, om) * Px(muv + s_, x, hh, om)
        okD &= abs(lhs - D(2 * x / om, thv, hh, Bv)) < 1e-25
report("CDD factor P_x enters the continued breather amplitude as D(chi), chi=2x/om (h=3,4)", okD)

thv = mp.mpc(0.71, 0)
ratios = []
for Bv in (mp.mpf("1e-2"), mp.mpf("1e-3"), mp.mpf("1e-4")):
    r = (D(1 - Bv / 2, thv, 3, Bv) / D(1, thv, 3, Bv) - 1) / Bv ** 2
    ratios.append(r)
okO = abs(D(1, thv, 3, 0) - 1) < 1e-30 and abs(D(mp.mpf("0.37"), thv, 3, 0) - 1) < 1e-30
okO &= abs(ratios[2] - ratios[1]) < 1e-2 * abs(ratios[2]) and abs(ratios[2]) > 0.1
report("D(chi)=1 at B=0; P_{(om+1)/2} instead of P_{om/2} changes K at O(B^2) (nonzero)", okO,
       ", ".join(mp.nstr(abs(r), 8) for r in ratios))


def Kunif(a, th, h, B):
    p = mp.mpf(1)
    for c in range(1, a + 1):
        p *= blk(c - 1, th, h) * blk(h + 1 - c - B / 2, th, h) / (blk(h - c, th, h) * blk(c - B / 2, th, h))
    return p


def KN(a, th, h, B):
    p = mp.mpf(1)
    for c in range(1, a + 1):
        p *= blk(c - 1, th, h) * blk(h - c + B / 2, th, h) / (blk(h - c, th, h) * blk(c - 1 + B / 2, th, h))
    return p


th = mp.mpc(0.61, -0.19)
Bv = mp.mpf(3) / mp.sqrt(17)
okdual = all(abs(KN(a, th, hh, Bv) - Kunif(a, th, hh, 2 - Bv)) < 1e-30 for hh in (3, 4, 5) for a in range(1, hh))
# K_1^N from continued three blocks
hh = 4
okdual &= abs(blk(-(hh - 1), th, hh) * blk(hh - 1 + Bv / 2, th, hh) * blk(-Bv / 2, th, hh) - KN(1, th, hh, Bv)) < 1e-30
# alternative CDD choice breaks exact duality
Kalt = KN(1, th, 3, Bv) * D(1 - Bv / 2, th, 3, Bv) / D(1, th, 3, Bv)
okdual &= abs(Kalt - Kunif(1, th, 3, 2 - Bv)) > 1e-3
report("K_a^N(B) = K_a(2-B) exactly; (-n)(n+B/2)(-B/2) = K_1^N; the P_{(om+1)/2} choice is not dual", okdual)

# K^(-) relation
okm = all(abs(-blk(a, th, hh) * blk(hh - a, th, hh) * Kunif(a, th, hh, Bv)
              - Kunif(a, th, hh, Bv) / (-1 / (blk(a, th, hh) * blk(hh - a, th, hh)))) < 1e-30
          for hh in (3, 4) for a in range(1, hh))
report("-(a)(h-a) K_a = K_a / K_a^cl with K_a^cl = -1/((a)(h-a))", okm)


# V10 (L6-12): classical moment of inertia of the a_2 solitonic vacuum
def inertia(chi):
    s = np.sqrt(3) / 2  # sin(pi/3); m_1 = 1, beta = 1
    alpha = np.array([[-1, 0, 1], [1, -1, 0], [0, 1, -1]], dtype=float)  # alpha_0, alpha_1, alpha_2

    def integrand(x):
        u = np.exp(x)
        r = []
        for j in range(3):
            cj = np.cos(chi + 2 * np.pi * j / 3)
            tau = 1 + 2 * cj / (1 + s) * u + (1 - s) / (1 + s) * u ** 2
            dtau = -2 * np.sin(chi + 2 * np.pi * j / 3) / (1 + s) * u
            r.append(dtau / tau)
        v = -(alpha.T @ np.array(r))  # d phi / d chi (beta = 1)
        return float(v @ v)

    pts = [-40, -10, -3, -1, -0.3, -0.05, 0]
    return sum(quad(integrand, pts[i], pts[i + 1], limit=200, epsabs=1e-13, epsrel=1e-12)[0] for i in range(len(pts) - 1))


Is = [inertia(c) for c in (0.0, 0.31, -0.77, 0.9, 1.7, 2.6)]
report("a_2 solitonic vacuum: I = int |d phi/d chi|^2 dx = 6/(beta^2 m_1), chi-independent", max(abs(v - 6) for v in Is) < 1e-8,
       ", ".join(f"{v:.10f}" for v in Is))
# rotor-like limit: m_1 sin^2(n delta)/sin(delta) -> m_1 delta n^2 = beta^2 m_1 n^2/12 with delta = pi B/6, B ~ beta^2/(2 pi)
bet2 = 1e-4
Bsm = 2 * bet2 / (bet2 + 4 * np.pi)
dl = np.pi * Bsm / 6
okrot = all(abs(np.sin(nq * dl) ** 2 / np.sin(dl) / (nq ** 2 / (2 * 6 / bet2)) - 1) < 1e-3 for nq in (1, 2, 3))
report("diagonal tower m_1 sin^2(n delta)/sin delta -> n^2/(2I) as beta -> 0", okrot)

# V11 (L6-13): Coleman-Thun residue identity at B = 9/20 with high precision
mp.mp.dps = 50
h = 3
Bct = mp.mpf(9) / 20


def S11(th, B):
    return blk(0, th, 3) * blk(2, th, 3) / (blk(B, th, 3) * blk(2 - B, th, 3))


def Kex1(n, m, th, B):
    K = Kunif(1, th, 3, B)
    f = lambda x: blk(x, th, 3)
    return (K * f(mp.mpf(1) / 2 - (2 * n + 1) * B / 2) * f(mp.mpf(5) / 2 + (2 * m + 1) * B / 2)
            * f(mp.mpf(1) / 2 - (2 * n - 1) * B / 2) * f(mp.mpf(5) / 2 + (2 * m - 1) * B / 2)
            * f(mp.mpf(3) / 2 + (2 * n - 1) * B / 2) * f(mp.mpf(3) / 2 - (2 * m - 1) * B / 2)
            * f(-mp.mpf(5) / 2 + (2 * n + 1) * B / 2) * f(-mp.mpf(1) / 2 - (2 * m + 1) * B / 2))


def Kex2(n, m, th, B):
    return Kex1(m, n, th, B)


def res(fun, th0):
    vals = [e * fun(th0 + e) for e in (mp.mpf(10) ** -20, mp.mpf(10) ** -25)]
    return vals[1], abs(vals[1] - vals[0])


thp = 1j * mp.pi / 3 * (mp.mpf(1) / 2 + 3 * Bct / 2)
thc = 1j * mp.pi / 3 * (mp.mpf(3) / 2 - 3 * Bct / 2)
thR = 1j * mp.pi / 3 * (3 * Bct / 2 - mp.mpf(1) / 2)
rK, e1 = res(lambda t: Kex1(-2, 2, t, Bct), thp)
rC, e2 = res(lambda t: Kex1(-1, 1, t, Bct), thc)
rS, e3 = res(lambda t: S11(t, Bct), 2j * mp.pi / 3)
dK2 = mp.diff(lambda t: Kex2(-1, 1, t, Bct), thR)
zero_val = abs(Kex2(-1, 1, thR, Bct))
g2 = 2 * rC / 1j
G2 = rS / 1j
rhs = g2 * G2 / 2 * dK2
okct = abs(rK - rhs) < mp.mpf(10) ** -15 * abs(rK) and zero_val < mp.mpf(10) ** -30 and abs(dK2) > 0.1
okct &= abs(mp.re(-1j * rK) - mp.mpf("-1.1874897238607920702")) < 1e-15 and abs(mp.im(-1j * rK)) < 1e-20
report("CT: Res K_1^{b_-2,2} = (1/2)g^2 Gamma^2 K_2^{b_-1,1}'(theta_R), simple zero, -iRes = -1.18749 (B=9/20)", okct,
       f"-iRes={mp.nstr(mp.re(-1j*rK), 12)}, rel diff={mp.nstr(abs(rK-rhs)/abs(rK), 3)}")
# order of pole: eps*K bounded and eps^2*K -> 0 ; also counting 6 internal lines, 2 loops for the described graph
Vtx, Edg = 5, 6  # vertices E,F,R,S,A ; internal lines EF, FR, RS, SA (bulk) and ER, RA (boundary)
report("CT graph: 6 internal lines (4 bulk + 2 boundary), 2 loops; naive order 6-2*2=2, one zero -> 1 (5-4-1=0 printed is wrong)",
       Edg - Vtx + 1 == 2 and Edg - 2 * (Edg - Vtx + 1) - 1 == 1 and 5 - 2 * 2 - 1 == 0)

# V12 (L6-15): a_2 fusion angle and Exercise 6.1 identities
mp.mp.dps = 40
u = mp.acos(-mp.mpf(1) / 2)
okang = abs(u - 2 * mp.pi / 3) < 1e-30 and abs(mp.sqrt(2 * (1 + mp.cos(mp.pi / 3))) - mp.sqrt(3)) < 1e-30
th = mp.mpc(0.71, 0.05)
Bv = mp.sqrt(2) / 3
okang &= all(abs(blk(x, 2 * th, 3) + blk(x / 2, th, 3) * blk(x / 2 + 3, th, 3)) < 1e-30 for x in (Bv, 2 - Bv, mp.mpf(2), mp.mpf("0.3")))
okang &= all(abs(blk(x, th + 1j * mp.pi / 3, 3) * blk(x, th - 1j * mp.pi / 3, 3) - blk(x + 1, th, 3) * blk(x - 1, th, 3)) < 1e-30
             for x in (Bv, mp.mpf(2)))
K2 = Kunif(1, th + 1j * mp.pi / 3, 3, Bv) * S11(2 * th, Bv) * Kunif(1, th - 1j * mp.pi / 3, 3, Bv)
okang &= abs(K2 - Kunif(1, th, 3, Bv)) < 1e-30
report("a_2 fusion angle u_11^2 = 2pi/3 (pi/3 gives mass sqrt3 m); doubling and shift identities; K_2 = K_1", okang)

# V13 (L6-16): one loop, false second expression, Kim's [1/2][3/2]
th = mp.mpf("0.71")
K1N = lambda B: blk(2 + B / 2, th, 3) / (blk(2, th, 3) * blk(B / 2, th, 3))
alt = lambda B: -blk(2 + B / 2, th, 3) / (blk(2, th, 3) * blk(3 - B / 2, th, 3))
sq = lambda x, B: blk(x - mp.mpf(1) / 2, th, 3) * blk(x + mp.mpf(1) / 2, th, 3) / (
    blk(x - mp.mpf(1) / 2 + B / 2, th, 3) * blk(x + mp.mpf(1) / 2 - B / 2, th, 3))
Kim = lambda B: sq(mp.mpf(1) / 2, B) * sq(mp.mpf(3) / 2, B)
Bt = mp.mpf("0.43")
ok13 = abs(K1N(Bt) - alt(Bt)) > 0.1
ok13 &= abs(Kim(Bt) - Kim(2 - Bt)) < 1e-30
d1N = mp.diff(K1N, 0)
d1K = mp.diff(Kim, 0)
d2N = mp.diff(K1N, 0, 2)
d2K = mp.diff(Kim, 0, 2)
# beta^2 coefficient: B = beta^2/(2pi) + O(beta^4)
coef = d1N / (2 * mp.pi)
target = -mp.mpf(1) / 4 * mp.cos(th / 2j) / mp.sin(3 * th / 2j)
ok13 &= abs(coef - target) < 1e-25 and abs(d1N - d1K) < 1e-25 and abs(d2N - d2K) > 0.1
# first-order block formula
x0 = mp.mpf("0.8")
dblk = mp.diff(lambda B: blk(x0 + B / 2, th, 3) / blk(x0, th, 3), 0)
ok13 &= abs(dblk - 1j * mp.pi / 6 * mp.sinh(th) / (mp.cosh(th) - mp.cos(mp.pi * x0 / 3))) < 1e-25
report("Ex 6.2: O(beta^2) coefficient of K_1^N; second printed form false; Kim [1/2][3/2] self-dual, same O(B), different O(B^2)", ok13,
       f"|second-derivative diff|={mp.nstr(abs(d2N-d2K), 8)}")

# V14 (L6-17): energies, degeneracy, ranges
Fk = lambda k, d: np.sin(k * d) * np.cos(np.pi / 6 - k * d) / np.sin(d)
d = np.pi * 0.37 / 6
ok14 = all(abs(sum(np.cos(np.pi / 6 - (2 * k + 1) * d) for k in range(nk)) - Fk(nk, d)) < 1e-13 for nk in range(1, 8))
ok14 &= abs(Fk(-1, d) + Fk(1, d) - np.sin(d)) < 1e-14
d = np.pi / 12
ok14 &= abs(Fk(-2, d) + Fk(2, d) - Fk(1, d)) < 1e-14 and abs(Fk(1, d) - np.cos(np.pi / 12)) < 1e-14
# ranges: creation pole of b_{n+1,m} at p = 1/2-(2n+1)B/2 > 0 iff n+1 < 1/(2B)+1/2;
# creation pole of b_{-n-1,n+1} at p = 3/2-(2n+1)B/2 > 0 iff n+1 < 3/(2B)+1/2
for Bv in (0.23, 0.45, 0.61):
    for nk in range(0, 6):
        ok14 &= ((0.5 - (2 * nk + 1) * Bv / 2 > 0) == (nk + 1 < 1 / (2 * Bv) + 0.5))
        ok14 &= ((1.5 - (2 * nk + 1) * Bv / 2 > 0) == (nk + 1 < 3 / (2 * Bv) + 0.5))
report("Ex 6.3: telescoping F(n), F(-1)+F(1)=sin delta, degeneracy at B=1/2; positive-label range bounds from strip poles", ok14)

print(f"{sum(RESULTS)}/{len(RESULTS)} PASS")
