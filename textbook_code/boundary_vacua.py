"""Vacua of real-coupling a_n^(1) affine Toda theory on the half-line x <= 0, and linearized reflection (Part VII, ch. 31).

Book conventions (sec-boundary-potential): real coupling b > 0,
  V = (m^2/b^2) sum_j (e^{b alpha_j.phi} - 1),  B = (2m/b^2) sum_j C_j (e^{b alpha_j.phi/2} - 1)   (a_n: n_j = 1, alpha^2 = 2),
  d_x phi(0) = -(m/b) sum_j C_j alpha_j e^{b alpha_j.phi(0)/2}.
Continued Hirota solutions: phi = -(1/b) sum_j alpha_j ln tau_j with the tau_j of eq-n-soliton (imaginary-coupling beta = -i b).
Checks (ids of ch. 31: sec-real-coupling-vacua, prp-solitonic-vacua, prp-dg-weights, sec-boundary-missing-charges,
sec-classical-reflection-factors, eq-linear-reflection):
1. Delius-Gandenberger hep-th/9904002 Theorem 4.1 for n = 1..7: lambda -> C_j = (-1)^{alpha_j.lambda} maps the weights of the
   fundamental representations 2-to-1 (lambda, -lambda) onto the 2^n - 1 sign patterns with prod_{j=0}^n C_j = 1, C != (+...+);
   formula (4.8)-(4.9) gives the positive preimage and its representation; the number of C_j = -1 is twice the number of
   cyclic blocks of the weight.  Delius hep-th/9807189 (1.15)-(1.16) is the same statement with A_i = -C_i.
2. Which solitonic patterns have a single-soliton vacuum (charges realized by single solitons, sec-topological-charges):
   all for a_1, a_2; a_3 misses 2 of 7 (the alternating patterns), a_4 misses 5 of 15, a_5 21 of 31.
3. The continued static soliton-mirror pair for a_2 (and a_3, a_4): real, regular on x <= 0 for chi in open intervals,
   solves the static field equation and the boundary condition with C_j = (-1)^{alpha_j.lambda} (lambda = charge of the soliton),
   energy E = (2m/b^2)(2k - h m_a/m) (k = number of C_j = -1) independent of the modulus chi; the parity partner
   (d = 1/(1-s)) is singular on x < 0.
4. a_{2k-1}, C = -1: the static self-conjugate soliton tau_j = 1 + (-1)^j e^{2mx+xi} satisfies the boundary condition for every
   xi (sympy); for real xi < 0 it is real and regular on x <= 0 with energy 0, the energy of phi = 0 (a modulus).
5. Linearized reflection about phi = 0 (uniform C): K_a = (sinh th - i C s_a)/(sinh th + i C s_a) = -((a)(h-a))^{-C};
   boundary bound states for C = -1 at omega = m_a cos(pi a/h) (Corrigan-Dorey-Rietdijk hep-th/9407148 (4.7)-(4.9), (5.3));
   B -> 0 limit of Delius-Gandenberger (3.16) is the C = +1 value.
"""
import itertools
import numpy as np
import sympy as sp
from scipy.integrate import quad


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


m, b = 1.0, 0.8

# ---------------------------------------------------------------- 1. Delius-Gandenberger Theorem 4.1


def dg_roots(n):
    """DG labels: alpha_i = e_i - e_{i+1} (i = 1..n), alpha_0 = e_{n+1} - e_1, in R^{n+1}; returns list indexed 0..n."""
    h = n + 1; e = np.eye(h)
    return [e[h - 1] - e[0]] + [e[i - 1] - e[i] for i in range(1, h)]


def pattern(lam, roots):
    return tuple(int(round((-1) ** round(float(a @ lam)))) for a in roots)


ok_map = ok_formula = ok_blocks = True
for n in range(1, 8):
    h = n + 1; al = dg_roots(n)
    weights = [np.array(v, float) for i in range(1, h) for v in itertools.product([0, 1], repeat=h) if sum(v) == i]
    pre = {}
    for v in weights:
        C = pattern(v, al)
        ok_map &= np.prod(C) == 1 and C != (1,) * h
        pre.setdefault(C, []).append(v)
        blocks = sum(1 for j in range(h) if v[j] == 1 and v[(j + 1) % h] == 0)
        ok_blocks &= C.count(-1) == 2 * blocks
    target = {C for C in itertools.product([1, -1], repeat=h) if np.prod(C) == 1 and C != (1,) * h}
    ok_map &= set(pre) == target and len(target) == 2 ** n - 1
    ok_map &= all(len(vs) == 2 and np.allclose(vs[0] + vs[1], 1) for vs in pre.values())   # lambda and -lambda = complement
    lam_fund = [None] + [np.concatenate([np.ones(i), np.zeros(h - i)]) for i in range(1, h)]   # lambda_i = eps_1+...+eps_i
    for C in target:
        P = [k for k in range(1, h) if C[k] == -1]
        a = len(P)
        lam = sum((-1) ** (a - j) * lam_fund[k] for j, k in enumerate(P, 1))
        irr = sum((-1) ** (a - j) * k for j, k in enumerate(P, 1))
        ok_formula &= (pattern(lam, al) == C and set(np.unique(lam)) <= {0.0, 1.0} and lam.sum() == irr
                       and lam[-1] == 0 and any(np.allclose(lam, v) for v in pre[C]))
report('DG Theorem 4.1, n = 1..7: lambda -> C_j = (-1)^{alpha_j.lambda} maps the fundamental weights 2-to-1 (lambda, -lambda) onto the '
       '2^n - 1 patterns with prod_{j=0}^n C_j = 1, C != (+...+)', ok_map)
report('DG (4.8)-(4.9): lambda = sum_j (-1)^{a-j} lambda_{k_j} is the positive preimage and lies in Lambda_i, i = sum_j (-1)^{a-j} k_j', ok_formula)
report('number of C_j = -1 equals twice the number of cyclic blocks of 1s of lambda in the eps basis', ok_blocks)

# Delius 98a (1.15)-(1.16): A_i = -e^{i pi alpha_i.lambda}, prod A_i = (-1)^{n+1}; with A = -C this is prod C = 1
ok = True
for n in range(1, 6):
    al = dg_roots(n)
    for v in itertools.product([0, 1], repeat=n + 1):
        if 0 < sum(v) < n + 1:
            A = [-np.exp(1j * np.pi * (a @ np.array(v, float))) for a in al]
            ok &= abs(np.prod(A) - (-1) ** (n + 1)) < 1e-12 and abs(np.prod([-x for x in A]) - 1) < 1e-12
report('Delius 98a (1.16) prod_i A_i^(lambda) = (-1)^{n+1} is DG (4.2) prod_i C_i = 1 under A_i = -C_i', ok)

# ---------------------------------------------------------------- 2. which patterns have a single-soliton vacuum


def book_roots(n):
    """book/boundary_images labels: alpha_j = e_j - e_{j+1 mod h}, j = 0..n (cyclic)."""
    h = n + 1; e = np.eye(h)
    return [e[j] - e[(j + 1) % h] for j in range(h)]


def soliton_charge(n, a, imxi):
    """(beta/2pi) phi(+infty) of the single soliton (eq-an-soliton), sec-topological-charges, as a vector in R^h."""
    h = n + 1; al = book_roots(n)
    th = np.array([2 * np.pi * j * a / h + imxi for j in range(h)])
    red = (th + np.pi) % (2 * np.pi) - np.pi
    return -sum(al[j] * red[j] for j in range(h)) / (2 * np.pi)


def as_weight(w, h):
    """match a sum-zero vector to a 0/1 vector v with w = v - (|v|/h) e; returns tuple v or None."""
    for i in range(1, h):
        v = w + i / h
        if np.allclose(v, np.round(v), atol=1e-9) and set(np.round(v)) <= {0.0, 1.0} and round(v.sum()) == i:
            return tuple(int(x) for x in np.round(v))
    return None


missing = {}
ok = True
for n in range(1, 8):
    h = n + 1; al = book_roots(n)
    realized = set()
    for a in range(1, h):
        ws = set()
        for imxi in np.linspace(0, 2 * np.pi, 7 * h, endpoint=False) + 0.0123:
            v = as_weight(soliton_charge(n, a, imxi), h)
            ok &= v is not None and sum(v) == a
            ws.add(v)
        ok &= len(ws) == h // np.gcd(a, h)
        realized |= ws
    pats = {pattern(np.array(v, float), al) for v in realized}
    allpats = {C for C in itertools.product([1, -1], repeat=h) if np.prod(C) == 1 and C != (1,) * h}
    missing[n] = sorted(allpats - pats)
counts = {n: (2 ** n - 1 - len(missing[n]), 2 ** n - 1) for n in missing}
print('   solitonic patterns with a single-soliton vacuum (have, total): ' + ', '.join(f'a_{n}: {c[0]}/{c[1]}' for n, c in counts.items()))
print(f'   a_3 patterns (C_0, C_1, C_2, C_3) without one: {missing[3]}')
report('single-soliton vacua: species a realizes h/gcd(a,h) weights; a_1, a_2 cover all solitonic patterns, a_3 misses exactly the '
       'two alternating patterns (+-+-), (-+-+)', ok and missing[1] == [] and missing[2] == []
       and missing[3] == [(-1, 1, -1, 1), (1, -1, 1, -1)])

# ---------------------------------------------------------------- 3. the continued static pair (Bowcock; Delius 98a sec. 4)


def pair_taus(n, a, chi, x, d=None):
    """tau_j = 1 + 2 d cos(chi + 2 pi j a/h) e^{m_a x} + A d^2 e^{2 m_a x}, d = 1/(1+s), A = cos^2(pi a/h); and x-derivatives."""
    h = n + 1; s = np.sin(np.pi * a / h); ma = 2 * m * s; A = 1 - s**2
    d = 1 / (1 + s) if d is None else d
    E = np.exp(ma * x)
    out = []
    for j in range(h):
        c = 2 * d * np.cos(chi + 2 * np.pi * j * a / h)
        out.append((1 + c * E + A * d**2 * E**2, ma * c * E + 2 * ma * A * d**2 * E**2, ma**2 * c * E + 4 * ma**2 * A * d**2 * E**2))
    return out


def fields(n, taus):
    al = book_roots(n); h = n + 1
    phi = -sum(al[j][:, None] * np.log(taus[j][0])[None, :] for j in range(h)) / b
    dphi = -sum(al[j][:, None] * (taus[j][1] / taus[j][0])[None, :] for j in range(h)) / b
    ddphi = -sum(al[j][:, None] * (taus[j][2] / taus[j][0] - (taus[j][1] / taus[j][0]) ** 2)[None, :] for j in range(h)) / b
    return phi, dphi, ddphi


def energy(n, phi_of_x, C):
    """E = int_{-inf}^0 (phi'^2/2 + V) dx + B(phi(0)) for a static field."""
    al = book_roots(n)
    def dens(x):
        phi, dphi, _ = phi_of_x(x)
        return 0.5 * dphi @ dphi + (m**2 / b**2) * sum(np.exp(b * (a_ @ phi)) - 1 for a_ in al)
    bulk = quad(dens, -np.inf, 0, epsabs=1e-12, epsrel=1e-12, limit=400)[0]
    phi0 = phi_of_x(0.0)[0]
    return bulk + (2 * m / b**2) * sum(c * (np.exp(b * (a_ @ phi0) / 2) - 1) for c, a_ in zip(C, al))


ok_real = ok_eom = ok_bc = ok_E = ok_mod = True
lines = []
for n, a in [(2, 1), (3, 1), (4, 1), (4, 2)]:
    h = n + 1; al = book_roots(n); s = np.sin(np.pi * a / h); ma = 2 * m * s
    xs = np.linspace(-15, 0, 3001)
    # chi in an open interval between the zeros cos(chi + 2 pi j a/h) = -1 of tau_j(0); pick points inside
    sing = sorted(((np.pi - 2 * np.pi * j * a / h) % (2 * np.pi)) for j in range(h))
    gaps = [(sing[k], sing[k + 1] if k + 1 < len(sing) else sing[0] + 2 * np.pi) for k in range(len(sing))]
    Es = {}
    for lo, hi in gaps:
        if hi - lo < 1e-9:
            continue
        for frac in (0.3, 0.5, 0.8):
            chi = lo + frac * (hi - lo)
            T = pair_taus(n, a, chi, xs)
            mins = min(np.min(t[0]) for t in T)
            ok_real &= mins > 1e-6
            phi, dphi, ddphi = fields(n, T)
            rhs = (m**2 / b) * sum(al[j][:, None] * np.exp(b * (al[j] @ phi)) for j in range(h))
            ok_eom &= np.max(np.abs(ddphi - rhs)) < 1e-8
            # which sign pattern does it satisfy?
            p0, dp0 = phi[:, -1], dphi[:, -1]
            good = [C for C in itertools.product([1, -1], repeat=h)
                    if np.max(np.abs(dp0 + (m / b) * sum(c * a_ * np.exp(b * (a_ @ p0) / 2) for c, a_ in zip(C, al)))) < 1e-9]
            lam = as_weight(soliton_charge(n, a, chi), h)          # Im xi of the soliton = chi
            pred = pattern(np.array(lam, float), al)
            ok_bc &= good == [pred] and np.prod(pred) == 1
            k = pred.count(-1)
            Eth = (2 * m / b**2) * (2 * k - h * ma / m)
            En = energy(n, lambda x: tuple(v[:, 0] for v in fields(n, pair_taus(n, a, chi, np.array([x])))), pred)
            ok_E &= abs(En - Eth) < 1e-7
            Es.setdefault(pred, []).append(En)
        ok_mod &= all(np.ptp(v) < 1e-7 for v in Es.values())
    lines.append(f'a_{n} a={a}: patterns {sorted(Es)} E b^2/m = ' + ', '.join(f'{np.mean(v) * b**2 / m:+.6f}' for v in Es.values()))
    # parity partner d = 1/(1-s): singular on x < 0 for every (generic) chi; at isolated chi the zero is a double zero
    for chi in np.linspace(0, 2 * np.pi, 37, endpoint=False) + 0.05:
        T = pair_taus(n, a, chi, np.linspace(-15, -1e-3, 6001), d=1 / (1 - s))
        ok_real &= min(np.min(t[0]) for t in T) < 0
for l in lines:
    print('   ' + l)
report('continued static pair (d = 1/(1+s)): real and regular on x <= 0 for chi inside each interval; the partner d = 1/(1-s) is singular on x<0', ok_real)
report('continued static pair solves the static real-coupling field equation phi\'\' = (m^2/b) sum alpha_j e^{b alpha_j.phi}  (a_2, a_3, a_4)', ok_eom)
report('it satisfies exactly one sign pattern, C_j = (-1)^{alpha_j.lambda} with lambda the charge of the soliton, prod C_j = 1', ok_bc)
report('energy (book normalization of B) = (2m/b^2)(2k - h m_a/m), k = #{C_j = -1}; independent of the modulus chi', ok_E and ok_mod)

# ---------------------------------------------------------------- 4. self-conjugate static soliton, C = -1, h even
x_, xi_ = sp.symbols('x xi')
ms_, bs_ = sp.symbols('m b', positive=True)
ok = True
for n in (1, 3, 5):
    h = n + 1; al = [sp.Matrix([int(v) for v in r]) for r in book_roots(n)]
    E_ = sp.exp(2 * ms_ * x_ + xi_)
    tau = [1 + (-1) ** j * E_ for j in range(h)]
    phi = -sum((al[j] * sp.log(tau[j]) for j in range(h)), sp.zeros(h, 1)) / bs_
    dphi = phi.diff(x_)
    # e^{b alpha_k.phi/2} = (tau_{k-1} tau_{k+1})^{1/2}/tau_k = tau_{k+1}/tau_k here
    ex = [tau[(k + 1) % h] / tau[k] for k in range(h)]
    bc = dphi + (ms_ / bs_) * sum((-1 * al[k] * ex[k] for k in range(h)), sp.zeros(h, 1))
    ok &= all(sp.simplify(c) == 0 for c in bc)
    ok &= all(sp.simplify(sp.exp(bs_ * (al[k].T * phi)[0] / 2).subs({x_: 0, xi_: -1}) - ex[k].subs({x_: 0, xi_: -1})) == 0 for k in range(h))
report('a_1, a_3, a_5: the static self-conjugate soliton tau_j = 1 + (-1)^j e^{2mx+xi} satisfies the C = -1 boundary condition for every xi', ok)
ok = True; Evals = []
for n in (1, 3):
    h = n + 1
    for xi in (-0.3, -1.0, -3.0):
        def fx(x, xi=xi, n=n):
            h = n + 1; E = np.exp(2 * m * np.atleast_1d(x) + xi)
            T = [(1 + (-1) ** j * E, (-1) ** j * 2 * m * E, (-1) ** j * 4 * m**2 * E) for j in range(h)]
            return tuple(v[:, 0] for v in fields(n, T))
        if n == 1:
            continue
        En = energy(n, fx, (-1,) * h)
        Evals.append(En)
        ok &= np.exp(xi) < 1 and abs(En) < 1e-8
report(f'a_3, C = -1: for real xi < 0 the static soliton is real and regular on x <= 0 with energy 0 = E(phi=0), independent of xi  '
       f'[E = {max(abs(np.array(Evals))):.1e}]', ok)

# ---------------------------------------------------------------- 5. linearized reflection about phi = 0


def block(x, th, h):
    """book eq-blocks: (x) = sinh(th/2 + i pi x/2h)/sinh(th/2 - i pi x/2h)."""
    return np.sinh(th / 2 + 1j * np.pi * x / (2 * h)) / np.sinh(th / 2 - 1j * np.pi * x / (2 * h))


ok_K = ok_bs = ok_dg = True
for n in (2, 3, 4, 5):
    h = n + 1; al = book_roots(n)
    M2 = m**2 * sum(np.outer(a_, a_) for a_ in al)
    w, V = np.linalg.eigh(M2)
    for C in (1, -1):
        # linearized BC: eps' = -(m/2) C sum alpha alpha^T eps = -(C/2m) M2 eps
        for th in (0.37, 1.3, 0.2 + 0.4j):
            for k in range(1, h):                                  # eigenvalue m_a^2 with a determined from w
                ma2 = w[k]; ma = np.sqrt(ma2); a = int(round(h / np.pi * np.arcsin(ma / (2 * m))))
                p = ma * np.sinh(th)
                # eps = v (e^{ipx} + K e^{-ipx}): ip(1 - K) = -(C ma^2/2m)(1 + K)
                K = (1j * p + C * ma2 / (2 * m)) / (1j * p - C * ma2 / (2 * m))
                s_a = np.sin(np.pi * a / h)
                Kf = (np.sinh(th) - 1j * C * s_a) / (np.sinh(th) + 1j * C * s_a)
                Kb = -(block(a, th, h) * block(h - a, th, h)) ** (-C)
                ok_K &= abs(K - Kf) < 1e-10 and abs(K - Kb) < 1e-10
        # bound states for C = -1: pole at th = i pi a/h (a < h/2), omega = m_a cos(pi a/h)
    for a in range(1, (h + 1) // 2):
        th = 1j * np.pi * a / h; ma = 2 * m * np.sin(np.pi * a / h)
        den = np.sinh(th) + 1j * (-1) * np.sin(np.pi * a / h)
        ok_bs &= abs(den) < 1e-12 and abs(ma * np.cosh(th) - ma * np.cos(np.pi * a / h)) < 1e-12
    # Delius-Gandenberger (3.16): K_a = prod_{c=1}^a (c-1)(c-n-1)(-c+B/2)(-c-n-B/2); B -> 0 gives the C = +1 value
    for a in range(1, h):
        for th in (0.41, 1.7):
            Bsm = 1e-7
            Kdg = np.prod([block(c - 1, th, h) * block(c - n - 1, th, h) * block(-c + Bsm / 2, th, h) * block(-c - n - Bsm / 2, th, h)
                           for c in range(1, a + 1)])
            ok_dg &= abs(Kdg - (-1 / (block(a, th, h) * block(h - a, th, h)))) < 1e-5
report('linearized reflection about phi = 0: K_a = (sinh th - i C s_a)/(sinh th + i C s_a) = -((a)(h-a))^{-C} (CDR (4.7); a_2..a_5)', ok_K)
report('C = -1: poles at th = i pi a/h, boundary bound states omega_a = m_a cos(pi a/h) = m sin(2 pi a/h) (CDR (4.9))', ok_bs)
th = 0.77
report('CDR (5.3): a_2, C = -1: (ip - 3m/2)/(ip + 3m/2) = -(1)(2)',
       abs((1j * np.sqrt(3) * np.sinh(th) - 1.5) / (1j * np.sqrt(3) * np.sinh(th) + 1.5) - (-(block(1, th, 3) * block(2, th, 3)))) < 1e-12)
report('Delius-Gandenberger (3.16) at B -> 0 equals the C = +1 linearized reflection factor -1/((a)(h-a)) (a_2..a_5)', ok_dg)
