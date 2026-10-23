"""Referee check of prp-soliton-preserving-half (ch. 30) and the three-soliton
images of ch. 31, written from scratch (no code imports).

Conventions (book, ch. 10 and ch. 30): a_n^(1), h=n+1, alpha_j = e_j - e_{j+1}
(alpha_0 = e_h - e_1), |alpha|^2 = 2, m = beta = 1.
  phi = i sum_j alpha_j ln tau_j,
  tau_j = sum_S prod_{k in S} omega^{j a_k} E_k prod_{k<l} A_{a_k a_l}(theta_k - theta_l),
  E = exp(m_a (x cosh th - t sinh th) + xi),  m_a = 2 sin(pi a/h).
Imaginary-coupling Lax pair:
  a_t(lam) = (i/2) d_x phi.H + (1/2) sum_j c_j w_j (lam E_j - E_{-j}/lam),
  w_j = exp(i alpha_j.phi/2).
Checks
 1. 'if' direction: K = exp(i pi nu.H) satisfies K a_t(lam) + a_t(1/lam)^dag K = 0
    for random d_x Re phi, random Im phi(0), d_x Im phi = 0 and
    Re phi(0) = pi nu + 2 pi mu, in the defining rep of a_2, a_3, in Lambda^2
    of a_3 and in the vector rep of d_4.
 2. 'only if' (a_n defining rep): the solution space for K is 1-dim on the
    half lattice, 0-dim off it.
 3. Three-soliton images: tau_j(x)/conj(tau_j(-x)) = c_j F(x,t) with |c_j| = 1
    (so phi(x)+conj(phi(-x)) is constant), and the class of nu.
 4. Which classes of Lambda_W / 2 Lambda_W are reached by the topological
    charges of static single solitons (a_2, a_3, a_4, a_5).
"""
import itertools
import numpy as np

rng = np.random.default_rng(5)


def roots_an(n):
    h = n + 1
    al = []
    for j in range(h):
        v = np.zeros(h)
        if j == 0:
            v[h - 1] += 1; v[0] -= 1
        else:
            v[j - 1] += 1; v[j] -= 1
        al.append(v)
    return np.array(al)        # rows alpha_0..alpha_n in R^h (sum-zero plane)


def defining_rep(n):
    h = n + 1
    Es = []
    for j in range(h):
        M = np.zeros((h, h), dtype=complex)
        if j == 0:
            M[h - 1, 0] = 1      # E_{alpha_0} = e_{h,1}: raises e_1 -> weight e_h? we use H=diag
        else:
            M[j - 1, j] = 1
        Es.append(M)
    # H_a = diag(e_a) in the R^h embedding: phi.H = diag(phi)
    return Es


def Hdot(v):
    return np.diag(v).astype(complex)


def check_compat(K, Es, Hfun, al, phi0, dphi, lam):
    """K a_t(lam) + a_t(1/lam)^dag K at imaginary coupling"""
    def at(l):
        A = 0.5j * Hfun(dphi)
        for j, E in enumerate(Es):
            w = np.exp(0.5j * al[j] @ phi0)
            A = A + 0.5 * w * (l * E - E.conj().T / l)
        return A
    return np.linalg.norm(K @ at(lam) + at(1 / lam).conj().T @ K)


def fundamental_coweights(n):
    h = n + 1
    out = []
    for k in range(1, h):
        v = np.zeros(h); v[:k] = 1; v -= k / h
        out.append(v)
    return np.array(out)


def check1():
    worst = 0
    for n in (2, 3, 4):
        al = roots_an(n)
        Es = defining_rep(n)
        cw = fundamental_coweights(n)
        for coeffs in itertools.product([0, 1], repeat=n):
            nu = np.array(coeffs) @ cw
            mu = rng.integers(-2, 3, size=n) @ cw
            re0 = np.pi * nu + 2 * np.pi * mu
            im0 = rng.normal(size=n + 1); im0 -= im0.mean()
            dre = rng.normal(size=n + 1); dre -= dre.mean()
            phi0 = re0 + 1j * im0
            K = np.diag(np.exp(1j * np.pi * nu))
            lam = 0.83
            worst = max(worst, check_compat(K, Es, Hdot, al, phi0, dre + 0j, lam))
    return worst


def lambda2_rep(n):
    """Lambda^2 of the defining rep of a_n, with H, E as matrices"""
    h = n + 1
    pairs = [(i, j) for i in range(h) for j in range(i + 1, h)]
    idx = {p: k for k, p in enumerate(pairs)}
    N = len(pairs)

    def lift(M):
        R = np.zeros((N, N), dtype=complex)
        for (i, j), k in idx.items():
            # M acts on e_i ^ e_j
            for a in range(h):
                if M[a, i] != 0:
                    p = (a, j)
                    if a != j:
                        s = 1 if a < j else -1
                        R[idx[tuple(sorted(p))], k] += s * M[a, i]
                if M[a, j] != 0:
                    p = (i, a)
                    if a != i:
                        s = 1 if i < a else -1
                        R[idx[tuple(sorted(p))], k] += s * M[a, j]
        return R
    Es = [lift(E) for E in defining_rep(n)]
    Hfun = lambda v: lift(np.diag(v).astype(complex))
    return Es, Hfun


def check1_lambda2():
    n = 3
    al = roots_an(n)
    Es, Hfun = lambda2_rep(n)
    cw = fundamental_coweights(n)
    worst = 0
    for coeffs in itertools.product([0, 1], repeat=n):
        nu = np.array(coeffs) @ cw
        re0 = np.pi * nu
        im0 = rng.normal(size=n + 1); im0 -= im0.mean()
        dre = rng.normal(size=n + 1); dre -= dre.mean()
        K = Hfun(np.zeros(n + 1))
        K = __import__('scipy.linalg', fromlist=['expm']).expm(1j * np.pi * Hfun(nu))
        worst = max(worst, check_compat(K, Es, Hfun, al, re0 + 1j * im0, dre + 0j, 1.21))
    return worst


def check1_d4():
    # d_4 vector rep in basis (e1..e4,-e4..-e1); H = diag(phi, -phi reversed)
    N = 8

    def E(i, j):
        m = np.zeros((N, N), dtype=complex); m[i, j] = 1; return m
    Es1 = [E(0, 1) - E(6, 7), E(1, 2) - E(5, 6), E(2, 3) - E(4, 5), E(2, 4) - E(3, 5)]
    # E_{alpha_0}, alpha_0 = -e1-e2: maps e1 -> -e2 component etc. lowest root generator
    E0 = E(7, 1) - E(6, 0)
    Es = [E0] + Es1
    al = np.array([[-1, -1, 0, 0], [1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [0, 0, 1, 1]], float)
    Hfun = lambda v: np.diag(np.concatenate([v, -v[::-1]])).astype(complex)
    # check E's are root vectors
    for a, Em in zip(al, Es):
        v = rng.normal(size=4)
        assert np.linalg.norm(Hfun(v) @ Em - Em @ Hfun(v) - (a @ v) * Em) < 1e-12
    # coweights of d4 (simply laced: weights). Fundamental weights:
    w = np.array([[1, 0, 0, 0], [1, 1, 0, 0], [.5, .5, .5, -.5], [.5, .5, .5, .5]])
    worst = 0
    for coeffs in itertools.product([0, 1], repeat=4):
        nu = np.array(coeffs) @ w
        im0 = rng.normal(size=4); dre = rng.normal(size=4)
        K = np.diag(np.exp(1j * np.pi * np.concatenate([nu, -nu[::-1]])))
        worst = max(worst, check_compat(K, Es, Hfun, al, np.pi * nu + 1j * im0, dre + 0j, 0.77))
    return worst


def check2():
    """solution space of K (any h x h matrix), defining rep, for given boundary data:
       Dirichlet on Re phi (d_x Re phi free), Neumann on Im phi."""
    res = []
    for n in (2, 3):
        al = roots_an(n); Es = defining_rep(n); h = n + 1
        cw = fundamental_coweights(n)
        for case in ('half', 'off'):
            nu = np.array([1] + [0] * (n - 1)) @ cw
            re0 = np.pi * nu if case == 'half' else np.pi * nu + 0.3 * rng.normal(size=h)
            re0 = re0 - re0.mean()
            im0 = rng.normal(size=h); im0 -= im0.mean()
            cols = []
            for k in range(h * h):
                K = np.zeros(h * h, complex); K[k] = 1; K = K.reshape(h, h)
                rows = []
                for trial in range(h):
                    dre = np.zeros(h); dre[trial] = 1; dre -= dre.mean()
                    # linear in K: evaluate the residual matrix (not norm)
                    def at(l, dphi):
                        A = 0.5j * Hdot(dphi)
                        for j, E in enumerate(Es):
                            w = np.exp(0.5j * al[j] @ (re0 + 1j * im0))
                            A = A + 0.5 * w * (l * E - E.conj().T / l)
                        return A
                    lam = 0.9
                    rows.append((K @ at(lam, dre) + at(1 / lam, dre).conj().T @ K).reshape(-1))
                cols.append(np.concatenate(rows))
            A = np.array(cols).T
            s = np.linalg.svd(A, compute_uv=False)
            null = int(np.sum(s < 1e-10 * s[0])) + max(0, h * h - len(s))
            res.append((n, case, null))
    return res


# ---------------- solitons ----------------
def Acoef(a, b, th, h):
    return (np.cosh(th) - np.cos(np.pi * (a - b) / h)) / (np.cosh(th) - np.cos(np.pi * (a + b) / h))


def taus(sols, x, t, h):
    """sols: list of (a, theta, xi); returns array tau_j, j=0..h-1"""
    om = np.exp(2j * np.pi / h)
    N = len(sols)
    E = [np.exp(2 * np.sin(np.pi * a / h) * (x * np.cosh(th) - t * np.sinh(th)) + xi) for a, th, xi in sols]
    out = np.zeros(h, complex)
    for r in range(N + 1):
        for S in itertools.combinations(range(N), r):
            pref = 1.0 + 0j
            for k in S:
                pref *= E[k]
            for k, l in itertools.combinations(S, 2):
                pref *= Acoef(sols[k][0], sols[l][0], sols[k][1] - sols[l][1], h)
            for j in range(h):
                out[j] += pref * om ** (j * sum(sols[k][0] for k in S))
    return out


def image_check(n, a, b, th, xi1, im3, xs, ts):
    h = n + 1
    A_ab = Acoef(a, b, th, h)
    xi3 = -np.log(A_ab) + 1j * im3
    xi2 = -np.conj(xi1) - np.log(Acoef(a, a, 2 * th, h)) - np.log(A_ab)
    sols = [(a, th, xi1), (a, -th, xi2), (b, 0.0, xi3)]
    cs = []
    for x in xs:
        for t in ts:
            tp = taus(sols, x, t, h)
            tm = np.conj(taus(sols, -x, t, h))
            r = tp / tm
            cs.append(r / r[0])
    cs = np.array(cs)
    spread = np.max(np.abs(cs - cs[0]))
    c = cs[0]
    al = roots_an(n)
    # Re phi(0,t) from the continuous branch of ln tau_j along x in (-L, 0]
    nus = []
    for t in ts:
        xg = np.linspace(-14.0, 0.0, 28001)
        T = np.array([taus(sols, x, t, h) for x in xg])
        ph = np.unwrap(np.angle(T), axis=0)[-1]
        nus.append(-(ph @ al) / np.pi)          # beta Re phi(0)/pi
    nus = np.array(nus)
    nu = nus[0]
    tdep = np.max(np.abs(nus - nu))
    pair = al[1:] @ nu
    integral = np.max(np.abs(pair - np.round(pair)))
    even = np.all(np.round(pair).astype(int) % 2 == 0)
    # without the correction -ln A_ab in xi2 and xi3 the image fails:
    sols_bad = [(a, th, xi1), (a, -th, -np.conj(xi1) - np.log(Acoef(a, a, 2 * th, h))), (b, 0.0, 1j * im3)]
    cb = []
    for x in xs[:4]:
        for t in ts[:3]:
            r = taus(sols_bad, x, t, h) / np.conj(taus(sols_bad, -x, t, h))
            cb.append(r / r[0])
    cb = np.array(cb)
    return spread, np.max(np.abs(np.abs(c) - 1)), max(integral, tdep), even, np.round(pair).astype(int), np.max(np.abs(cb - cb[0]))


def weight_class(nu, n):
    al = roots_an(n)
    p = np.round(al[1:] @ nu).astype(int) % 2
    return tuple(p)


def static_charges(n, b, ny=721):
    """classes of the topological charges of a static species-b soliton, Re xi=0, as Im xi varies"""
    h = n + 1
    om = np.exp(2j * np.pi / h)
    al = roots_an(n)
    classes = set()
    charges = set()
    for y in np.linspace(-np.pi, np.pi, ny)[1:-1] + 1e-3:
        theta = np.angle(om ** (b * np.arange(h)) * np.exp(1j * y))   # Arg in (-pi,pi)
        if np.min(np.abs(np.abs(theta) - np.pi)) < 1e-6:
            continue
        mu = -(theta @ al) / (2 * np.pi)
        charges.add(tuple(np.round(mu, 6)))
        classes.add(weight_class(mu, n))
    return classes, charges


def all_classes(n):
    h = n + 1
    out = {}
    for b in range(1, h):
        for S in itertools.combinations(range(h), b):
            v = np.zeros(h); v[list(S)] = 1; v -= b / h
            c = weight_class(v, n)
            out.setdefault(c, set()).add(b)
    return out


if __name__ == '__main__':
    ok = True
    w1 = check1(); w2 = check1_lambda2(); w3 = check1_d4()
    print(f"1. compat. with K=exp(i pi nu.H): a_2..a_4 defining {w1:.1e}, a_3 Lambda^2 {w2:.1e}, d_4 vector {w3:.1e}")
    ok &= max(w1, w2, w3) < 1e-12
    r = check2()
    print("2. nullity of K (n, case, nullity):", r)
    ok &= all((nl == 1) if c == 'half' else (nl == 0) for _, c, nl in r)
    xs = np.linspace(-2.0, 2.0, 9); ts = np.linspace(-1.5, 1.5, 4)
    print("3. three-soliton images (n,a,b): spread of c_j, ||c_j|-1|, non-integrality or t-dependence of nu, nu in 2Lambda?, (alpha_k.nu), spread without corrections")
    for (n, a, b, th, xi1, im3) in [(2, 1, 1, 0.7, 0.3 + 0.4j, 0.5), (2, 1, 2, 0.5, -0.2 + 1.1j, -0.9),
                                    (3, 1, 1, 0.6, 0.1 + 0.2j, 0.3), (3, 2, 1, 0.8, 0.4 - 0.7j, 1.2),
                                    (3, 1, 2, 0.6, 0.2 - 2.5j, -1.5), (4, 2, 1, 0.9, 0.2 + 0.0j, -1.0)]:
        # phases of the last two chosen regular with a2_regular_scan.py (min |tau_j| / scale > 0.18)
        sp, mod, integ, even, pair, bad = image_check(n, a, b, th, xi1, im3, xs, ts)
        print(f"   a_{n} a={a} b={b}: {sp:.1e} {mod:.1e} {integ:.1e} {even} {pair} {bad:.2f}")
        ok &= sp < 1e-9 and mod < 1e-9 and integ < 1e-9 and (not even) and bad > 1e-3
    print("4. classes of Lambda_W/2Lambda_W reached by static single solitons:")
    for n in (2, 3, 4, 5):
        allc = all_classes(n)
        reached = set()
        for b in range(1, n + 1):
            cl, ch = static_charges(n, b)
            reached |= cl
        nontriv = [c for c in allc if any(c)]
        missing = [c for c in nontriv if c not in reached]
        print(f"   a_{n}: {len(nontriv)} nontrivial classes, {len(missing)} not reached by a static single soliton;"
              f" missing (alpha_1..alpha_n parities): {missing}")
        if n == 3:
            ok &= len(missing) == 2
    print("PASS" if ok else "FAIL", "half-coweight soliton-preserving boundaries (see 4 for the classes without a single half-soliton)")
