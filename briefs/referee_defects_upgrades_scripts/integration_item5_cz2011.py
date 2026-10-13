"""Item 5: CZ2011 (1012.4186) (3.12) against its own (3.8), (3.10), with the spectral variable of (3.2).
Written from the paper's printed formulas (no book or referee code).
 a = pi/(3 xi); q = e^{i 4 pi^2/beta^2} (3.6) = i e^{i pi a} since xi = (2/3) pi beta^2/(8 pi - beta^2);
 z = i(theta - nu)/xi, so theta -> theta - i pi is z -> z + 3a and theta -> theta -+ i pi/3 is z -> z +- a;
 epsilon = tau~ = 1, epsilon~ = tau = e^{-pi nu/xi} (after (3.10)); rho = rho^ e^{-pi(theta-nu)/xi}/(2 pi), e^{-pi(theta-nu)/xi} = e^{i pi z}.
 (3.8):  rho(z) rho(z + 3a) (1 + tau^2 q^-4 x)(1 - tau^2 q^-2 x) = 1;   (3.10): rho(z) = (1 + tau^2 x) rho(z - a) rho(z + a).
 Reading A (corrected (3.2)): x = e^{2 pi theta/xi}, so tau^2 x = e^{2 pi (theta - nu)/xi} = e^{-2 pi i z}.
 Reading B (printed (3.2)):   x = q^{2 pi theta/xi}, so tau^2 x = e^{-2 pi nu/xi} q^{2 pi theta/xi}, theta = nu - i xi z.
Checks:
 1. exact finite identity: printed partial product P_K (outside factors + both products to K) / P_K(0)
    = prod_{k<K} F_k(z)/F_k(0) * Gamma(A_{K,1}-z) Gamma(A_{K,2}-z) / (Gamma(A_{K,1}) Gamma(A_{K,2})), A_{k,l} = 1/2 + (6k+l) a;
    and the unnormalized P_K(0) itself does not converge (z-independent divergence), so 'drift like K^{-2z}' needs the normalization;
 2. with H(z) = prod_k F_k(z)/F_k(0) (equal k) and rho = C e^{i pi z} H, C = 1/(2 H(a)^2), reading A satisfies (3.10) and (3.8);
 3. reading B fails both, at the same constant (for real nu, theta);
 4. pole positions: Gamma(1/2 - z + a), Gamma(1/2 - z + 2a) have poles at theta = nu - i(xi/2 + pi/3), nu - i(xi/2 + 2 pi/3);
    at the printed theta = nu - i(xi + pi/3), nu - i(xi + 2pi/3) their arguments are -1/2.
"""
import mpmath as mp
mp.mp.dps = 30

def setup(xi):
    a = mp.pi/(3*xi)
    beta2 = 8*mp.pi*(3*xi/2)/(mp.pi + 3*xi/2)          # invert xi = (2/3) pi b^2/(8 pi - b^2)
    q = mp.exp(1j*4*mp.pi**2/beta2)
    return a, beta2, q

def logF(z, k, a):
    A = lambda l: mp.mpf(1)/2 + (6*k + l)*a
    s = 0
    for l, sg in ((1, 1), (2, 1), (4, -1), (5, -1)):
        s += sg*(mp.loggamma(A(l) + z) + mp.loggamma(A(l) - z))
    return s

def logH(z, a, K=3000):
    # sum_k [logF_k(z) - logF_k(0)], terms O(k^-2): Richardson on partial sums at K, 2K, 4K
    terms = [logF(z, k, a) - logF(0, k, a) for k in range(4*K)]
    c = []
    acc = 0
    for t in terms:
        acc += t; c.append(acc)
    S1, S2, S4 = c[K-1], c[2*K-1], c[4*K-1]
    R1, R2 = 2*S2 - S1, 2*S4 - S2
    return (4*R2 - R1)/3

ok = True
# 1. finite identity (printed grouping) and divergence of the unnormalized product
xi = mp.mpf('0.93'); a, beta2, q = setup(xi)
def logP_printed(z, K):
    s = mp.loggamma(mp.mpf(1)/2 - z + a) + mp.loggamma(mp.mpf(1)/2 - z + 2*a)
    for k in range(1, K + 1):
        s += (mp.loggamma(0.5 + z + 6*k*a - 4*a) + mp.loggamma(0.5 + z + 6*k*a - 5*a)
              - mp.loggamma(0.5 + z + 6*k*a - 2*a) - mp.loggamma(0.5 + z + 6*k*a - a))
        s += (mp.loggamma(0.5 - z + 6*k*a + 2*a) + mp.loggamma(0.5 - z + 6*k*a + a)
              - mp.loggamma(0.5 - z + 6*k*a - 2*a) - mp.loggamma(0.5 - z + 6*k*a - a))
    return s
z0 = mp.mpc('0.21', '0.37')
for K in (5, 40):
    lhs = logP_printed(z0, K) - logP_printed(0, K)
    rhs = sum(logF(z0, k, a) - logF(0, k, a) for k in range(K))
    AK = lambda l: mp.mpf(1)/2 + (6*K + l)*a
    rhs += (mp.loggamma(AK(1) - z0) + mp.loggamma(AK(2) - z0) - mp.loggamma(AK(1)) - mp.loggamma(AK(2)))
    d = abs(mp.exp(lhs - rhs) - 1)
    print(f'1. K = {K}: |P_K(z)/P_K(0) / (equal-k partial product x endpoint) - 1| = {mp.nstr(d, 3)}')
    ok &= d < 1e-20
print('   log P_K(0) (unnormalized, z-independent) at K = 50, 100, 200:',
      [mp.nstr(mp.re(logP_printed(0, K)), 6) for K in (50, 100, 200)])
print('   endpoint log-modulus vs -2 Re(z) log(6aK) at K = 1000:',
      mp.nstr(mp.re(mp.loggamma(0.5 + 6*1000*a + a - z0) + mp.loggamma(0.5 + 6*1000*a + 2*a - z0)
                    - mp.loggamma(0.5 + 6*1000*a + a) - mp.loggamma(0.5 + 6*1000*a + 2*a)), 6),
      mp.nstr(-2*mp.re(z0)*mp.log(6*a*1000), 6))

# 2./3. functional equations
for xi in (mp.mpf('0.93'), mp.mpf('1.73')):
    a, beta2, q = setup(xi)
    ok &= abs(q - 1j*mp.exp(1j*mp.pi*a)) < 1e-25
    Ha = mp.exp(logH(a, a))
    C = 1/(2*Ha**2)
    rho = lambda z: C*mp.exp(1j*mp.pi*z)*mp.exp(logH(z, a))
    nu = mp.mpf('0.4')
    for theta in (mp.mpf('-0.3'), mp.mpf('0.55')):
        z = 1j*(theta - nu)/xi
        t2xA = mp.exp(-2j*mp.pi*z)
        t2xB = mp.exp(-2*mp.pi*nu/xi)*q**(2*mp.pi*theta/xi)
        for name, t2x in (('A (x = e^{2 pi theta/xi})', t2xA), ('B (printed x = q^{2 pi theta/xi})', t2xB)):
            e310 = rho(z)/((1 + t2x)*rho(z - a)*rho(z + a)) - 1
            e38 = rho(z)*rho(z + 3*a)*(1 + t2x*q**-4)*(1 - t2x*q**-2) - 1
            print(f'2/3. xi = {mp.nstr(xi,3)}, theta = {mp.nstr(theta,3)}, reading {name}: |(3.10) - 1| = {mp.nstr(abs(e310),3)}, |(3.8) - 1| = {mp.nstr(abs(e38),3)}')
            if name.startswith('A'):
                ok &= abs(e310) < 1e-7 and abs(e38) < 1e-7
            else:
                ok &= abs(e310) > 1e-3 and abs(e38) > 1e-3
    print(f'   tau^2 x = e^(-2 pi i z) holds identically in reading A; C = {mp.nstr(C, 10)}')

# 4. pole positions
xi = mp.mpf('0.93'); a = mp.pi/(3*xi); nu = mp.mpf('0.4')
for shift in (mp.pi/3, 2*mp.pi/3):
    th_true = nu - 1j*(xi/2 + shift)
    th_print = nu - 1j*(xi + shift)
    l = 1 if shift == mp.pi/3 else 2
    arg_true = mp.mpf(1)/2 - 1j*(th_true - nu)/xi + l*a
    arg_print = mp.mpf(1)/2 - 1j*(th_print - nu)/xi + l*a
    print(f'4. Gamma(1/2 - z + {l}a): argument at nu - i(xi/2 + {l}pi/3) = {mp.nstr(arg_true, 5)}, at printed nu - i(xi + {l}pi/3) = {mp.nstr(arg_print, 5)}')
    ok &= abs(arg_true) < 1e-25 and abs(arg_print + mp.mpf(1)/2) < 1e-25
print('PASS' if ok else 'FAIL', ': item 5 - (3.12) solves (3.8), (3.10) only with x = e^{2 pi theta/xi}; printed (3.2) x = q^{...} fails')
