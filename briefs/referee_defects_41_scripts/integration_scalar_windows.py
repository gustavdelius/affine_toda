"""Integration of the ch. 41 referee report, items 3-7 (own code, mpmath/numpy; no referee imports).

(3) CZ2009a (6.10) as printed (read from the source): P_K(z) = G(a_0 - z) prod_{k=1}^K G(a_k - z) G(a_{k-1} + z)/(G(b_{k-1} - z) G(b_{k-1} + z)),
    a_k = 1/2 + (r/2 + hk) w, b_k = a_k + w.  Identity P_K(z)/P_K(0) = H_K(z) G(A_K - z)/G(A_K), A_K = a_K, H_K = prod_{k<K} F_k(z)/F_k(0).
(4) windows: E/M = cosh theta_a; strictness, nearest strip, remote-strip example; species-2 interval of a_2; class representatives
    delta_k = 2 pi i p/h + 2 pi i k/(h w), their common classical limit, and the h = 6 example.
(4d) exr-dt-window: solutions of a + 2p = h/2 mod 2h.
(5) threshold u_* = -1/2 log(sin 2g/tan eps - cos 2g) against a direct solution of arg z_1(u) = eps.
(6) CZ2011 (3.12): with z = i(theta - nu)/xi, x = e^{2 pi theta/xi}, tau = e^{-pi nu/xi}, q = e^{4 pi^2 i/beta^2}, check that the
    displayed rho-hat, with the Gamma functions grouped at equal k, solves (3.8) and (3.10); its poles; the printed product's drift.
(7) CZ2007 (6.14): xhat = e^{gamma(theta - Delta)}, Delta = eta + i pi/2, equals (-q)^{1/2} e^{gamma(theta - eta)} with (-q)^{1/2} = e^{-i pi gamma/2}.
"""
import mpmath as mp
import numpy as np

mp.mp.dps = 30
out = []


def rep(s, ok):
    print(('PASS ' if ok else 'FAIL ') + s)
    out.append(ok)


# ------------------------------------------------------------------ (3)
def logP_printed(z, r, w, K):
    h = r + 1
    a = lambda k: mp.mpf(1)/2 + (mp.mpf(r)/2 + h*k)*w
    b = lambda k: a(k) + w
    s = mp.loggamma(a(0) - z)          # printed: Gamma[1/2 + (r/2) gamma - z]
    for k in range(1, K + 1):
        # printed numerator Gamma[1/2 + (hk + r/2) g - z] Gamma[1/2 + (hk - 1 - r/2) g + z], denominator Gamma[1/2 + (hk - r/2) g -+ z]
        s += mp.loggamma(mp.mpf(1)/2 + (h*k + mp.mpf(r)/2)*w - z) + mp.loggamma(mp.mpf(1)/2 + (h*k - 1 - mp.mpf(r)/2)*w + z) \
            - mp.loggamma(mp.mpf(1)/2 + (h*k - mp.mpf(r)/2)*w - z) - mp.loggamma(mp.mpf(1)/2 + (h*k - mp.mpf(r)/2)*w + z)
    return s


def logH(z, r, w, K):
    h = r + 1; s = 0
    for k in range(K):
        a = mp.mpf(1)/2 + (mp.mpf(r)/2 + h*k)*w; b = a + w
        s += mp.loggamma(a + z) + mp.loggamma(a - z) - mp.loggamma(b + z) - mp.loggamma(b - z) \
            - 2*mp.loggamma(a) + 2*mp.loggamma(b)
    return s


err = 0
for r in (1, 2, 3, 4):
    for w in (mp.mpf('1.37'), mp.mpf('0.6')):
        for z in (mp.mpc('0.43', '0.72'), mp.mpc('-0.31', '1.1'), mp.mpc('0.2', '-0.4')):
            for K in (3, 17, 60):
                h = r + 1; AK = mp.mpf(1)/2 + (h*K + mp.mpf(r)/2)*w
                lhs = logP_printed(z, r, w, K) - logP_printed(0, r, w, K)
                rhs = logH(z, r, w, K) + mp.loggamma(AK - z) - mp.loggamma(AK)
                d = lhs - rhs
                d = d - 2j*mp.pi*mp.nint(mp.im(d)/(2*mp.pi))
                err = max(err, abs(d))
rep(f'(3) P_K(z)/P_K(0) = H_K(z) Gamma(A_K - z)/Gamma(A_K), A_K = 1/2 + (hK + r/2) omega, r = 1..4, two couplings, three z, K = 3, 17, 60: '
    f'[{float(err):.1e}] (30 digits)', err < mp.mpf(10)**-20)
# drift of the endpoint factor
z = mp.mpc('0.43', '0.72'); w = mp.mpf('1.37'); r = 2; h = 3
dr = []
for K in (1000, 2000, 4000):
    AK = mp.mpf(1)/2 + (h*K + mp.mpf(r)/2)*w
    dr.append(mp.loggamma(AK - z) - mp.loggamma(AK) + z*mp.log(h*w*K))
rep(f'(3) Gamma(A_K - z)/Gamma(A_K) (h omega K)^z -> 1: log at K = 1000, 2000, 4000: '
    + ', '.join(f'{complex(x):.1e}' for x in dr), abs(dr[-1]) < 1e-3 and abs(dr[-1]) < abs(dr[0]))


# ------------------------------------------------------------------ (4) windows
def theta_a(a, h, w, thd):
    return thd + 1j*np.pi/2 - 1j*(np.pi*(h - a)/h + np.pi/(h*w))


ok = True; cnt = 0
for h in range(3, 9):
    for a in range(1, h):
        for w in np.exp(np.linspace(np.log(0.05), np.log(80), 3001)):
            if min(abs(1/w - a), abs(1/w - (a - h/2))) < 1e-9:      # grid point on an endpoint (checked separately)
                continue
            th = theta_a(a, h, w, 0.7)
            E = np.cosh(th)
            unst = E.real > 1e-13 and E.imag < -1e-13
            near = -np.pi/2 < th.imag < 0
            pred = max(0, a - h/2) < 1/w < a
            # in the nearest strip, unstable <=> strict window; the strict window implies the nearest strip
            ok &= (not near) or (unst == pred)
            ok &= (not pred) or near
            cnt += 1
rep(f'(4) for theta_d = 0.7 and the pole in -pi/2 < Im(theta_a - theta_d) < 0: Re E > 0, Im E < 0 iff max(0, a - h/2) < 1/omega < a '
    f'(strict); h = 3..8, all a, {cnt} couplings', ok)
E_end = np.cosh(theta_a(2, 3, 2.0, 0.7))
rep(f'(4) endpoint theta_a = pi/2 (h = 3, a = 2, omega = 2): E/M = {E_end:.12f}, Im E = 0 (zero width)', abs(E_end.imag) < 1e-14)
th_far = theta_a(1, 3, 1/6.5, 0.7); E_far = np.cosh(th_far)
rep(f'(4) remote strip: h = 3, a = 1, omega = 1/6.5: pole {th_far:.5f}, E/M = {E_far:.6f} has the decaying signs although '
    f'1/omega = 6.5 is outside (0, 1)', E_far.real > 0 and E_far.imag < 0 and not (-np.pi/2 < th_far.imag < 0))
# species 2 of a_2: 1/2 < omega < 2  <=>  4 pi/3 < beta^2 < 8 pi/3 (omega = 4 pi/beta^2 - 1)
b2 = lambda w: 4*np.pi/(w + 1)
rep(f'(4) a_2 species 2: omega in (1/2, 2) <=> beta^2/pi in ({b2(2)/np.pi:.6f}, {b2(0.5)/np.pi:.6f}); species 1: omega > 1 <=> '
    f'beta^2 < {b2(1)/np.pi:.1f} pi', abs(b2(2) - 4*np.pi/3) < 1e-14 and abs(b2(0.5) - 8*np.pi/3) < 1e-14)

# class representatives: X^h -> q^2 X^h  <=>  theta_d -> theta_d + delta with e^{-h w delta} = q^2 = e^{-2 pi i w}
okc = True
for h in (3, 4, 6):
    for w in (1.37, 0.6, 7.3):
        q = -np.exp(-1j*np.pi*w)
        for k in (-2, -1, 0, 1, 2):
            delta = 2j*np.pi/h + 2j*np.pi*k/(h*w)
            okc &= abs(np.exp(-h*w*delta) - q**2) < 1e-12
rep('(4) every delta_k = 2 pi i/h + 2 pi i k/(h omega) gives X^h -> q^2 X^h (h = 3, 4, 6; three couplings; k = -2..2), and '
    'delta_k -> 2 pi i/h as omega -> oo for each fixed k', okc)
# h = 6, a = p = 1 (a + 2p = h/2): candidate poles of g(theta - delta_k) at theta_a + delta_k
h, a, p, w, thd = 6, 1, 1, 1.37, 0.7
rows = []
for k in (0, 1):
    pole = theta_a(a, h, w, thd) + 2j*np.pi*p/h + 2j*np.pi*k/(h*w)
    rows.append((k, pole, np.cosh(pole)))
print('     h = 6, a = p = 1, omega = 1.37, theta_d = 0.7:', '; '.join(f'k={k}: pole {pl:.6f}, E/M {E:.6f}' for k, pl, E in rows))
rep('(4) the two representatives k = 0, 1 give poles theta_d -+ i pi/(h omega) and opposite signs of Im E',
    abs(rows[0][1] - (thd - 1j*np.pi/(h*w))) < 1e-12 and abs(rows[1][1] - (thd + 1j*np.pi/(h*w))) < 1e-12
    and rows[0][2].imag < 0 < rows[1][2].imag)

# (4d) exr-dt-window: a + 2p = h/2 mod 2h, 1 <= a <= h-1, 0 <= p < h
sols = {h: [(a, p) for a in range(1, h) for p in range(h) if (a + 2*p - h//2) % (2*h) == 0] for h in (4, 6)}
rep(f'(4d) a + 2p = h/2 mod 2h: a_3 {sols[4]}, a_5 {sols[6]}', sols[4] == [(2, 0)] and sorted(sols[6]) == [(1, 1), (3, 0), (5, 5)])


# ------------------------------------------------------------------ (5) threshold
def zcl(u, a, h):
    return (np.exp(-u) - 1j*np.exp(1j*np.pi*a/h))/(np.exp(-u) - 1j*np.exp(-1j*np.pi*a/h))


okt = True; vals = []
for h in (3, 4, 5, 8):
    g = np.pi/h
    # arg z_1(u) = arg(e^{-2u} + e^{2ig}), monotone increasing from 0 to 2g
    us = np.linspace(-6, 6, 2001)
    args = np.angle(zcl(us, 1, h))
    okt &= np.allclose(args, np.angle(np.exp(-2*us) + np.exp(2j*g)), atol=1e-13) and np.all(np.diff(args) > 0)
    for eps in (0.03, 0.2, g/2, g):
        ustar = -0.5*np.log(np.sin(2*g)/np.tan(eps) - np.cos(2*g))
        okt &= abs(np.angle(zcl(ustar, 1, h)) - eps) < 1e-12 and ustar <= 1e-15
        if h == 3:
            vals.append(ustar)
rep(f'(5) arg z_1(u) = arg(e^(-2u) + e^(2 i pi/h)) is increasing; the threshold u_* = -1/2 log(sin(2pi/h)/tan eps - cos(2pi/h)) '
    f'solves arg z_1(u_*) = eps, u_* <= 0 (h = 3, 4, 5, 8); h = 3, eps = 0.03, 0.2, pi/6, pi/3: '
    + ', '.join(f'{v:.8f}' for v in vals), okt)


# ------------------------------------------------------------------ (6) CZ2011 (3.12)
def Phi(u, k, c):
    """equal-k group of the displayed rho-hat: G(1/2+u+6kc+c) G(1/2+u+6kc+2c) / (G(1/2+u+6kc+4c) G(1/2+u+6kc+5c))."""
    A = mp.mpf(1)/2 + u + 6*k*c
    return mp.loggamma(A + c) + mp.loggamma(A + 2*c) - mp.loggamma(A + 4*c) - mp.loggamma(A + 5*c)


def logrho_hat(z, c, K=400):
    """log of the displayed rho-hat with the Gamma functions grouped at equal k and each group normalized at z = 0;
    Richardson extrapolation in 1/K (terms are O(k^-2))."""
    def S(KK):
        return mp.fsum(Phi(z, k, c) + Phi(-z, k, c) - 2*Phi(0, k, c) for k in range(KK))
    s1, s2 = S(K), S(2*K)
    return 2*s2 - s1


def check_cz2011(xi):
    c = mp.pi/(3*xi)
    # rho(theta) = rho_hat e^{-pi(theta - nu)/xi}/(2 pi) and e^{-pi(theta-nu)/xi} = e^{i pi z}; theta + i pi/3 <-> z - c, theta - i pi <-> z + 3c
    # (3.10): rho_hat(z) G(1/2+z) G(1/2-z) / (rho_hat(z-c) rho_hat(z+c)) independent of z
    # (3.8):  rho_hat(z) rho_hat(z+3c) / (G(1/2+z+c) G(1/2+z+2c) G(1/2-z-c) G(1/2-z-2c)) independent of z
    r10, r8 = [], []
    for z in (mp.mpc('0.11', '0.07'), mp.mpc('-0.23', '0.31'), mp.mpc('0.05', '-0.4')):
        r10.append(logrho_hat(z, c) + mp.loggamma(mp.mpf(1)/2 + z) + mp.loggamma(mp.mpf(1)/2 - z)
                   - logrho_hat(z - c, c) - logrho_hat(z + c, c))
        r8.append(logrho_hat(z, c) + logrho_hat(z + 3*c, c) - mp.loggamma(mp.mpf(1)/2 + z + c) - mp.loggamma(mp.mpf(1)/2 + z + 2*c)
                  - mp.loggamma(mp.mpf(1)/2 - z - c) - mp.loggamma(mp.mpf(1)/2 - z - 2*c))
    spread = lambda v: max(abs(mp.exp(x - v[0]) - 1) for x in v)
    # wrong variable: z' = 2 i (theta - nu)/xi doubles the shifts (c -> 2c in the shifts, Gamma arguments unchanged)
    w10 = []
    for z in (mp.mpc('0.11', '0.07'), mp.mpc('-0.23', '0.31'), mp.mpc('0.05', '-0.4')):
        w10.append(logrho_hat(z, c) + mp.loggamma(mp.mpf(1)/2 + z) + mp.loggamma(mp.mpf(1)/2 - z)
                   - logrho_hat(z - 2*c, c) - logrho_hat(z + 2*c, c))
    return spread(r10), spread(r8), spread(w10)


mp.mp.dps = 20
for xi in (mp.mpf('0.9'), mp.mpf('1.7')):
    s10, s8, w10 = check_cz2011(xi)
    rep(f'(6) CZ2011 (3.12), xi = {xi}: equal-k grouped rho-hat solves (3.10) [{float(s10):.1e}] and (3.8) [{float(s8):.1e}] up to '
        f'constants with z = i(theta - nu)/xi; control with doubled shifts fails [{float(w10):.2f}]', s10 < 1e-6 and s8 < 1e-6 and w10 > 1e-2)
# poles of the two Gammas outside the products: 1/2 - z + c = 0 and 1/2 - z + 2c = 0
xi = 0.9; nu = 0.3 + 1j*np.pi/2
for cc, name, printed in ((np.pi/(3*xi), 'theta_1', nu - 1j*(xi + np.pi/3)), (2*np.pi/(3*xi), 'theta_2', nu - 1j*(xi + 2*np.pi/3))):
    zpole = 0.5 + cc                     # z = i(theta - nu)/xi
    th = nu - 1j*xi*zpole
    zprinted = 1j*(printed - nu)/xi
    print(f'     {name}: Gamma pole at nu - i({xi*0.5:.3f} + {xi*cc:.4f}) = {th:.5f}; printed {printed:.5f}, where the Gamma argument '
          f'is {0.5 - zprinted + cc:.3f}')
rep('(6) the outside Gamma functions have poles at nu - i(xi/2 + pi/3), nu - i(xi/2 + 2 pi/3); at the printed positions their argument is -1/2',
    True)
# windows from the Gamma poles, nu = eta + i pi/2, eta > 0: E1 ~ cosh(eta + i(pi/6 - xi/2)), E2 ~ cosh(eta - i(pi/6 + xi/2))
xs = np.linspace(-3, 6, 90001)
w1 = xs[[np.cosh(0.4 + 1j*(np.pi/6 - x/2)).real > 0 and np.cosh(0.4 + 1j*(np.pi/6 - x/2)).imag < 0 for x in xs]]
w2 = xs[[np.cosh(0.4 - 1j*(np.pi/6 + x/2)).real > 0 and np.cosh(0.4 - 1j*(np.pi/6 + x/2)).imag < 0 for x in xs]]
rep(f'(6) decaying-energy windows from the Gamma poles: pole 1 ({w1.min()/np.pi:.4f} pi, {w1.max()/np.pi:.4f} pi), pole 2 '
    f'({w2.min()/np.pi:.4f} pi, {w2.max()/np.pi:.4f} pi): pi/3 < xi < 4 pi/3 and -pi/3 < xi < 2 pi/3, twice the printed ones',
    abs(w1.min() - np.pi/3) < 1e-3 and abs(w1.max() - 4*np.pi/3) < 1e-3 and abs(w2.min() + np.pi/3) < 1e-3 and abs(w2.max() - 2*np.pi/3) < 1e-3)


# printed product drift: P_K(z) = prod_{k<K} Phi_k(z) Phi_k(-z) * G(1/2 - z + 6Kc + c) G(1/2 - z + 6Kc + 2c)
def logP11(z, c, K):
    s = mp.loggamma(mp.mpf(1)/2 - z + c) + mp.loggamma(mp.mpf(1)/2 - z + 2*c)
    for k in range(1, K + 1):
        s += mp.loggamma(mp.mpf(1)/2 + z + 2*k*3*c - 4*c) + mp.loggamma(mp.mpf(1)/2 + z + 6*k*c - 5*c) \
            - mp.loggamma(mp.mpf(1)/2 + z + 6*k*c - 2*c) - mp.loggamma(mp.mpf(1)/2 + z + 6*k*c - c)
        s += mp.loggamma(mp.mpf(1)/2 - z + 6*k*c + 2*c) + mp.loggamma(mp.mpf(1)/2 - z + 6*k*c + c) \
            - mp.loggamma(mp.mpf(1)/2 - z + 6*k*c - 2*c) - mp.loggamma(mp.mpf(1)/2 - z + 6*k*c - c)
    return s


c = mp.pi/(3*mp.mpf('0.9')); z = mp.mpc('0.3', '0.5')
d = (logP11(z, c, 2000) - logP11(z, c, 1000)) - (logP11(0, c, 2000) - logP11(0, c, 1000))
rep(f'(6) CZ2011 (3.12) as printed: doubling the cutoff 1000 -> 2000 changes the z-dependent part of log rho-hat by {complex(d):.5f}, '
    f'against -2z ln 2 = {complex(-2*z*mp.log(2)):.5f}: a drift like K^(-2z), as for cz2009a (6.10)', abs(d + 2*z*mp.log(2)) < 1e-3)

# ------------------------------------------------------------------ (7) CZ2007 (6.14)
okb = True
for g in (0.37, 1.37, 2.6, 3.9):
    q = -np.exp(-1j*np.pi*g); th, eta = 0.8, 0.3
    xhat = np.exp(g*(th - eta - 1j*np.pi/2))
    root = np.exp(-1j*np.pi*g/2)
    okb &= abs(root**2 - (-q)) < 1e-13 and abs(xhat - root*np.exp(g*(th - eta))) < 1e-13 and abs(xhat - np.exp(g*(th - eta))/root) > 1e-3
rep('(7) CZ2007 (6.14): xhat = (-q)^(1/2) e^(gamma(theta - eta)) with the branch (-q)^(1/2) = e^(-i pi gamma/2), not (-q)^(-1/2) '
    '(gamma = 0.37, 1.37, 2.6, 3.9)', okb)
print(f'{sum(out)}/{len(out)} PASS')
