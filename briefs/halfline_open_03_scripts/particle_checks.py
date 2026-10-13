"""Particle-level consequences (continued breather amplitudes), using the blocks of particle_reflection.py."""
from fractions import Fraction as Fr
import mpmath as mp
from pr_core import Prod, S, K_uniform, K_neumann, axioms, close, blk, I, PI, report, THS
# continuation (a)_I -> book block: a = c om/2 - d  ->  (c + d B)
def cont(terms, h, sign=1):
    """terms: list of ((c, d), e) meaning DG block (c om/2 - d)_I^e."""
    return Prod(h, {(c, d): e for (c, d), e in terms}, sign)
ok1 = ok2 = True
for n in range(1, 8):
    h = n + 1
    # K_B1^(+) = (-n om/2)(n om/2 - 1/2)(1/2): c, d with a = c om/2 - d
    Kp = cont([((-n, 0), 1), ((n, Fr(1, 2)), 1), ((0, Fr(-1, 2)), 1)], h)
    ok1 &= close(Kp, K_neumann(1, n))
    # K_B1^(-) = DG(3.11) sigma_1 = (-n om/2)(1/2 - T)(-om/2 - 1/2) * (-(om/2)(T - om/2)),  T = h om/2
    Km = cont([((-n, 0), 1), ((-h, Fr(-1, 2)), 1), ((-1, Fr(1, 2)), 1), ((1, 0), 1), ((h - 1, 0), 1)], h, -1)
    Km_closed = Prod(h, {(1, 0): 1, (h - 1, 0): 1}, -1)*K_uniform(1, n)
    ok2 &= close(Km, Km_closed)
report('continued K_B1^(+) = K^N_1 (eq-neumann-amplitude, a = 1), n = 1..7', ok1)
report('continued K_B1^(-) = -(1)(h-1) K^unif_1 = -(1)(h - B/2)(-1 + B/2), n = 1..7', ok2)
# K^(-)_a = -(a)(h-a) K^unif_a: axioms
def Kmin(a, n):
    return Prod(n + 1, sign=-1).mul((a, 0)).mul((n + 1 - a, 0))*K_uniform(a, n)   # (a)(h-a), a = h-a allowed
res = [axioms(lambda a, n=n: Kmin(a, n), n, 'minus') for n in range(1, 8)]
report('K^(-)_a = -(a)(h-a) K^unif_a: unitarity, crossing (both forms), K_a = K_{h-a}, bulk bootstrap, n = 1..7',
       all(all(r) for r in res))
report('K^(-)_a -> 1 as B -> 0, n = 1..7',
       all(abs(Kmin(a, n)(th, 0) - 1) < 1e-20 for n in range(1, 8) for a in range(1, n + 1) for th in THS))
report('n = 2: K^(-)_1 = Gandenberger (6.4) -(1)(-1+B/2)(3-B/2)',
       close(Kmin(1, 2), Prod(3, {(1, 0): 1, (-1, Fr(1, 2)): 1, (3, Fr(-1, 2)): 1}, -1)))
# the alternative for odd n without g: breather factor divided by g(mu_a)g(mu_b), continued = k_c, c = h/2 - 1
for n in (3, 5, 7):
    h = n + 1; c = Fr(h, 2) - 1
    kc = Prod(h, {(Fr(h, 2) + c, Fr(1, 2)): 1, (Fr(h, 2) - c, Fr(-1, 2)): 1, (-Fr(h, 2) + c, Fr(-1, 2)): 1,
                  (-Fr(h, 2) - c, Fr(1, 2)): 1})
    kg = Prod(h, {(-1, Fr(-1, 2)): 1, (h - 1, Fr(1, 2)): 1, (1, Fr(-1, 2)): 1, (h + 1, Fr(1, 2)): 1})
    eps = mp.mpf('1e-12'); th = mp.mpc('0.4', '0.1')
    dB = (kg(th, eps) - 1)/eps
    report(f'n = {n}: g(mu_a)g(mu_b) continues to DG k_c with c = h/2 - 1 (3.25); its O(B) term is '
           f'{mp.nstr(dB, 6)} (nonzero, so Kim\'s one loop separates the two choices for n = 3)',
           close(kg, kc) and abs(dB) > 1e-3)
res = [axioms(lambda a, n=n: K_neumann(a, n), n, 'N') for n in range(1, 8)]
report('K^N_a (eq-neumann-amplitude) generated from K^N_1: unitarity, crossing (both forms), K_a = K_{h-a}, '
       'bootstrap at every fusing, n = 1..7', all(all(r) for r in res))
report('K^N_a -> 1 as B -> 0, n = 1..7',
       all(abs(K_neumann(a, n)(th, 0) - 1) < 1e-20 for n in range(1, 8) for a in range(1, n + 1) for th in THS))
