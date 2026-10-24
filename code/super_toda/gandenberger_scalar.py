"""Gandenberger's relations between the scalar factors of soliton S-matrices and
minimal supersymmetric S-matrices (Gandenberger 1998, Sec. 5; sec-known-smatrices-notes).

Notation of Gandenberger: mu = theta/(2 pi i), h the Coxeter number, lambda the
coupling function of his chapter 3, and lambda = 1/h the special coupling.

Checks:
 1. N = 2 (his (5.6), (3.28), (5.10)).  For a_N^(1), h = N + 1, his scalar factor
    (3.28) with a = b = 1 at lambda = 1/h is
        F_11 = prod_j G(mu+j+1/h)G(mu+j-1/h)G(-mu+j)^2 / [G(-mu+j+1/h)G(-mu+j-1/h)G(mu+j)^2],
    and the minimal N = 2 factor (5.6) differs from it only in G(-mu+j)^2/G(mu+j)^2 ->
    G(-mu+j-1)G(-mu+j+1)/[G(mu+j-1)G(mu+j+1)].  The ratio of the partial products up to
    j = M is -(M - mu)/(M + mu), exactly, so G^(N=2)_11 = -F_11|_(lambda=1/h): the identity
    (5.10) holds up to an overall sign -1.  Checked for N = 2, 4 at complex mu.
 2. Fusion.  His (3.12) and (5.7) build F_ab and G^(N=2)_ab from ab copies of the (1,1)
    factor with the same shifts (lambda/2)(a + b - 2j - 2k + 2), lambda = 1/h.  The partial
    products then differ by prod_jk [-(M - mu_jk)/(M + mu_jk)], so G_ab = (-1)^(ab) F_ab
    at lambda = 1/h, not -F_ab.  Checked for N = 4 and (a, b) = (1, 2), (2, 1), (2, 2), (1, 3).
 3. N = 1 (his (4.4), (5.14), (5.17)).  For a_2n^(2), h = 2n + 1, the infinite product of
    (4.4) at lambda = 1/h equals the infinite product in (5.14) factor by factor, so (5.17)
    holds with the stated trigonometric prefactor.  Checked for n = 1, 2.
"""
import sys
import mpmath as mp

mp.mp.dps = 30
FAILS = []
lg = mp.loggamma


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        FAILS.append(name)


def F11_an(mu, N, lam, M):
    """log of the partial product (j <= M) of Gandenberger (3.28), a = b = 1, for a_N^(1)."""
    h = N + 1
    s = 0
    for j in range(1, M + 1):
        x = j * h * lam
        s += (lg(mu + x - lam / 2 * (2 * N) + 1) + lg(mu + x - lam)
              - lg(-mu + x - lam / 2 * (2 * N) + 1) - lg(-mu + x - lam)
              + lg(-mu + x - lam / 2 * (2 * N + 2) + 1) + lg(-mu + x)
              - lg(mu + x - lam / 2 * (2 * N + 2) + 1) - lg(mu + x))
    return s


def G11_N2(mu, h, M):
    """log of the partial product of Gandenberger (5.6)."""
    s = 0
    for j in range(1, M + 1):
        s += (lg(mu + j + mp.mpf(1) / h) + lg(mu + j - mp.mpf(1) / h) + lg(-mu + j - 1) + lg(-mu + j + 1)
              - lg(-mu + j + mp.mpf(1) / h) - lg(-mu + j - mp.mpf(1) / h) - lg(mu + j - 1) - lg(mu + j + 1))
    return s


ok = True
for N in (2, 4):
    h = N + 1
    for mu in (mp.mpc('0.13', '0.07'), mp.mpc('0.31', '-0.2')):
        for M in (10, 60):
            r = mp.exp(G11_N2(mu, h, M) - F11_an(mu, N, mp.mpf(1) / h, M))
            ok &= abs(r + (M - mu) / (M + mu)) < 1e-20
report('N = 2: the partial products of (5.6) and of (3.28) at lambda = 1/h differ by the factor '
       '-(M - mu)/(M + mu), so G^(N=2)_11 = -F_11 (Gandenberger (5.10) holds up to the sign -1)', ok)

ok = True
N, M = 4, 40
h = N + 1
for a, b in ((1, 2), (2, 1), (2, 2), (1, 3)):
    for mu in (mp.mpc('0.13', '0.07'), mp.mpc('0.31', '-0.2')):
        shifts = [mu + mp.mpf(a + b - 2 * j - 2 * k + 2) / (2 * h) for j in range(1, a + 1) for k in range(1, b + 1)]
        logG = sum(G11_N2(s, h, M) for s in shifts)                       # (5.7)
        logF = sum(F11_an(s, N, mp.mpf(1) / h, M) for s in shifts)       # (3.12) at lambda = 1/h
        finite_M = mp.fprod((M - s) / (M + s) for s in shifts)
        ok &= abs(mp.exp(logG - logF) - (-1) ** (a * b) * finite_M) < 1e-20
report('N = 2, fused: by (3.12) and (5.7), G^(N=2)_ab = (-1)^(ab) F_ab at lambda = 1/h '
       '((a, b) = (1, 2), (2, 1), (2, 2), (1, 3); N = 4)', ok)


def F11_a2n2(mu, n, lam, M):
    """log of the partial product of Gandenberger (4.4) for a_2n^(2)."""
    h = 2 * n + 1
    s = 0
    for j in range(1, M + 1):
        x = j * h * lam
        s += (lg(mu + x - lam) + lg(mu + x - 2 * n * lam + 1) - lg(-mu + x - lam) - lg(-mu + x - 2 * n * lam + 1)
              + lg(mu + x - (n + mp.mpf(1) / 2) * lam) + lg(mu + x - (n + mp.mpf(1) / 2) * lam + 1)
              - lg(-mu + x - (n + mp.mpf(1) / 2) * lam) - lg(-mu + x - (n + mp.mpf(1) / 2) * lam + 1)
              + lg(-mu + x - (2 * n + 1) * lam) + lg(-mu + x + 1) - lg(mu + x - (2 * n + 1) * lam) - lg(mu + x + 1)
              + lg(-mu + x - (n + mp.mpf(3) / 2) * lam) + lg(-mu + x - (n - mp.mpf(1) / 2) * lam + 1)
              - lg(mu + x - (n + mp.mpf(3) / 2) * lam) - lg(mu + x - (n - mp.mpf(1) / 2) * lam + 1))
    return s


def G11_N1_product(mu, h, M):
    """log of the infinite-product part of Gandenberger (5.14)."""
    s = 0
    a, half = mp.mpf(1) / h, mp.mpf(1) / 2
    for j in range(1, M + 1):
        s += (lg(mu + j + a) + lg(mu + j - a) - lg(-mu + j + a) - lg(-mu + j - a)
              + lg(-mu + j - 1) + lg(-mu + j + 1) - lg(mu + j - 1) - lg(mu + j + 1)
              + lg(mu + j - half) + lg(mu + j + half) - lg(-mu + j - half) - lg(-mu + j + half)
              + lg(-mu + j - a - half) + lg(-mu + j + a + half) - lg(mu + j - a - half) - lg(mu + j + a + half))
    return s


ok = True
for n in (1, 2):
    h = 2 * n + 1
    for mu in (mp.mpc('0.13', '0.07'), mp.mpc('0.31', '-0.2')):
        for M in (10, 60):
            ok &= abs(F11_a2n2(mu, n, mp.mpf(1) / h, M) - G11_N1_product(mu, h, M)) < 1e-20
report('N = 1: at lambda = 1/h the product (4.4) for a_2n^(2) equals the product in (5.14) factor by '
       'factor, so (5.17) holds as printed (n = 1, 2)', ok)

sys.exit(1 if FAILS else 0)
