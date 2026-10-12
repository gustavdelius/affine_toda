"""Solitons and particles through an a_r^(1) jump defect, book normalization (Part VIII, ch. 39).

Conventions: defect_common.py (phi on x < 0, psi on x > 0, sigma = e^{-theta_d}; type I_- is BCZ's defect, type I_+ its
time reverse with B_d^{I+}(sigma) = B_d^{I-}(1/sigma)).  Solitons at imaginary coupling (b = i beta), eq-an-soliton and
eq-n-soliton: phi = (i/beta) sum_j alpha_j ln tau_j, tau_j = 1 + omega^{ja} E, E = exp(m_a(x cosh th - t sinh th) + xi),
m_a = 2 m sin(pi a/h), soliton mass M_a = 2 h m_a/beta^2; Mhat = 2 h m/beta^2.  The basic delay is eq-dc-ar-delay,
    z_a(u) = (e^{-u} - i e^{i pi a/h})/(e^{-u} - i e^{-i pi a/h})     (cz2009a (4.10); defect_classical_an.py check 8).
Jump class p in Z_h: (beta/2pi)(psi(0) - phi(0)) = lambda_p mod the root lattice (lambda_p a fundamental weight).
Linear waves at real coupling: e^{-i(w t - p x)}, w = m_k cosh th, p = m_k sinh th, on the eigenvector rho_k of the
rotation pi (pi rho_k = omega^k rho_k), particle k of mass m_k; T = transmitted/incident amplitude.

Checks (book ids they support):
 1. Delays, a_2 and a_3, all species, both types, jump classes p: the soliton on x < 0 and the soliton with E -> z E on
    x > 0 satisfy the defect conditions at 7 times (real and complex rapidities, complex xi) with
       type I_-: z = z_a(theta - theta_d - 2 pi i p/h),   type I_+: z = z_a(theta_d - 2 pi i p/h - theta);
    independent least-squares fits of z (type I_+, class 1); time reversal maps I_-(sigma) solutions to I_+(1/sigma)
    solutions; B_d(phi, psi + 2pi lambda_p/beta; sigma) = B_d(phi, psi; sigma e^{-+2 pi i p/h}) + const (root shifts:
    invariant); dictionary: BCZ hep-th/0401020 (6.2) at theta -> -theta, cz2007 (2.20) at a -> h - a.
                                                                              [eq-an-defect-delay, prp-an-defect-delay]
 2. Stored charge and selection rule: static coweight jumps solve the conditions, class-p energy and momentum
    Mhat[(cosh, sinh) theta_d - (cosh, sinh)(theta_d +- 2 pi i p/h)] (= cz2007 (2.17)-(2.18), cz2009a (4.8) with their
    constants; mass shell Mhat; real iff 2p = 0 mod h); a transmitted soliton changes the jump by lambda_in - lambda_out
    (continuation in t at x = 0 vs. continuation in x), jumps lie in the weight lattice; arg z_a is monotone with range
    (0, 2 pi a/h) [a < h/2] or (-2 pi abar/h, 0) [a > h/2]; charge boundaries in Im xi are 2 pi gcd(a,h)/h apart (cz2009a;
    cz2007's 2 gamma_a only if a | h); the outgoing charge is k steps along the cycle, 0 <= k <= a/gcd(a,h), maximum
    attained (a_2, a_3, a_4); theta - theta_d -> +inf: lambda_out = pi lambda_in.              [prp-stored-charge]
 3. Absorption and emission: at theta_abs = theta_d -+ i pi/2 +- i pi (a + 2p)/h (types I_-+), (soliton, const) and
    (const, soliton) solve the conditions (theta_d of both signs; a_2, a_3, a_5); Re cosh theta_abs has the sign of
    sin(pi(a + 2p)/h): if positive the class changes p -> p + a and the defect stores (Delta B_d, Delta P_d) =
    M_a (cosh, sinh) theta_abs; if negative the profile is inverted (species abar at theta_abs - i pi, classes p + a -> p);
    real theta_abs iff a + 2p = h/2 mod 2h (h <= 8);
    E + B_d and P + P_d constant by quadrature through real-rapidity absorptions (a_3 a = 2 p = 0; a_5 a = 1 p = 1);
    envelope velocity tanh theta_d; theta_d -> -inf limit of type I_- is the rotation psi = pi phi.     [prp-ds-absorption]
 4. Several solitons: two-soliton solutions (a_2, a_3; both types; class 1) pass with E_k -> z(theta_k) E_k and A_ab
    unchanged; two defects in a row (I_- then I_+) multiply the delays; classical bootstrap
    z_a(u + i pi b/h) z_b(u - i pi a/h) = z_{a+b}(u) (= cz2007 (2.25)), h = 3..8.                    [sec-defect-multi]
 5. Reality: z_abar(conj theta; conj theta_d) = conj z_a; a Theta-invariant a_2 breather passes a real-sigma defect and
    stays Theta-invariant; the Theta-image of a class-p solution solves the conditions in class -p (h - p); class-p
    and class-(h-p) energies are conjugate; h even, p = h/2 is sigma -> -sigma; at real rapidity |z| > 1 iff
    cos(2 pi p/h) > 0.                                                                               [sec-defect-reality]
 6. Linear waves (real coupling), a_1..a_4, both types, every particle k: R = 0, no mixing, |T| = 1,
    T^{I-}_k(theta) = (e^{-u} + e^{i pi k/h})/(e^{-u} + e^{-i pi k/h}) = z_k(u + i pi/2), u = theta - theta_d,
    T^{I+}_k(theta) = z_k(theta_d - theta - i pi/2); from the right Tbar_k(theta) = T_kbar(i pi - theta), and
    T_k(theta) Tbar_k(-theta) = 1; a_1 gives eq-shg-linear-transmission; limits 1 and omega^k; T_k(theta_d) = e^{i pi k/h},
    pole at theta_d - i(pi - pi k/h) and zero at theta_d + i(pi - pi k/h) (I_-), -e^{i pi k/h} and pole/zero at
    theta_d -+ i pi k/h (I_+); cz2007 in their conventions
    (D -> D^T, rho_k = -sum alpha_l omega^{lk}, m = 1): the corrected (7.11) (e^{-theta} in the denominator) is their
    result and equals (7.5) for a_2, the printed (7.11) is not; about a class-p defect at imaginary coupling the factor is
    T_k at theta_d + 2 pi i p/h, unimodular only for 2p = 0 mod h.    [eq-an-linear-transmission, prp-an-linear-transmission]
 7. Frozen soliton (pre-registered: does a pair of defects act on every species as a soliton does?):
    X_ab(i sinh u) = A_ab(u + i pi/2) and T_b(u) = z_b(u + i pi/2) (the same continuation); two defects of the same type
    reproduce kappa A_cb^{+-1} (any parameters) iff b = h/2 (zero/pole matching, h = 2..8, and a grid search);
    the mixed pair I_- at theta_c + i(pi/2 + pi c/h), I_+ at theta_c - i(pi/2 + pi c/h) gives omega^b A_cb(theta - theta_c)
    and omega^b X_cb for every b and c, carries M_c (cosh, sinh) theta_c with BCZ's constants, and is realized at real
    sigma by the classes 4p = h + 2c mod 4h when h + 2c = 0 mod 4.                          [prp-defect-frozen-soliton]
Runtime: about 10 seconds.
"""
import itertools
import time
from math import gcd
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import simpson
from defect_common import an_data, dmatrix, defect_potential, defect_dx

BETA, M = 0.8, 1.1
T_START = time.time()
rng = np.random.default_rng(39)


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


# ======================================================================== helpers
def zdel(a, h, u):
    """eq-dc-ar-delay."""
    x = np.exp(-u); g = np.pi*a/h
    return (x - 1j*np.exp(1j*g))/(x - 1j*np.exp(-1j*g))


def zdelay(a, h, theta, thd, typ=1, p=0):
    """eq-an-defect-delay: type I_- (typ = 1) or I_+ (typ = -1), jump class p."""
    if typ == 1:
        return zdel(a, h, theta - thd - 2j*np.pi*p/h)
    return zdel(a, h, thd - 2j*np.pi*p/h - theta)


def mass(a, h, m=M):
    return 2*m*np.sin(np.pi*a/h)


def Aint(a, b, h, th):                     # eq-interaction-coefficient
    return (np.cosh(th) - np.cos(np.pi*(a - b)/h))/(np.cosh(th) - np.cos(np.pi*(a + b)/h))


def soliton_terms(h, sols):
    """tau_j = sum_S prod_{k in S} omega^{j a_k} E_k prod_{k<l in S} A, as a list of terms (c_j, kx, wt):
    tau_j = sum c_j exp(kx x - wt t).  sols: list of (a, theta, xi), xi complex (delays go into xi)."""
    om = np.exp(2j*np.pi/h); js = np.arange(h); terms = []
    for k in range(len(sols) + 1):
        for S in itertools.combinations(range(len(sols)), k):
            c = np.ones(h, complex); kx = 0j; wt = 0j
            for i in S:
                a, th, xi = sols[i]
                c = c*om**(js*a)*np.exp(xi); kx = kx + mass(a, h)*np.cosh(th); wt = wt + mass(a, h)*np.sinh(th)
            for i1, i2 in itertools.combinations(S, 2):
                c = c*Aint(sols[i1][0], sols[i2][0], h, sols[i1][1] - sols[i2][1])
            terms.append((c, kx, wt))
    return terms


def const_terms(h, a):
    """tau_j = omega^{ja}: the right asymptote of a species-a soliton (a static coweight configuration)."""
    return [(np.exp(2j*np.pi*np.arange(h)*a/h), 0j, 0j)]


def field(r, terms, shift=None):
    """(phi, phi_x, phi_t)(x, t) of the Hirota form, plus (2 pi/beta) shift; x may be an array (then r x n arrays)."""
    al = an_data(r)['alpha']
    sh = np.zeros(r) if shift is None else (2*np.pi/BETA)*np.asarray(shift)

    def f(x, t):
        x = np.asarray(x, dtype=complex)
        ex = [np.exp(k*x - w*t) for _, k, w in terms]
        tau = sum(np.multiply.outer(c, e) for (c, _, _), e in zip(terms, ex))
        tx = sum(np.multiply.outer(c*k, e) for (c, k, _), e in zip(terms, ex))
        tt = sum(np.multiply.outer(-c*w, e) for (c, _, w), e in zip(terms, ex))
        ph = (1j/BETA)*np.tensordot(al, np.log(tau), axes=(0, 0))
        phx = (1j/BETA)*np.tensordot(al, tx/tau, axes=(0, 0))
        pht = (1j/BETA)*np.tensordot(al, tt/tau, axes=(0, 0))
        return ph + (sh if ph.ndim == 1 else sh[:, None]), phx, pht
    return f


def resid(r, F, G, thd, typ=1, ts=np.linspace(-2.5, 2.5, 7), x=0.0):
    """max |residual| of the defect conditions (imaginary coupling) at the points (x, t)."""
    out = 0.0
    for t in ts:
        ph, phx, pht = F(x, t); ps, psx, pst = G(x, t)
        _, gp, gq, _ = defect_potential(ph, ps, np.exp(-thd), r, BETA, M, 'imag', typ)
        px, sx = defect_dx(pht, pst, gp, gq, r, typ)
        out = max(out, np.abs(phx - px).max(), np.abs(psx - sx).max())
    return out


def lam(r, p):
    """fundamental weight lambda_p (p mod h; lambda_0 = 0) in the coordinates of an_data(r)."""
    return an_data(r)['lam'][p % (r + 1)]


def pair(r, a, theta, xi, z, p=0):
    """phi = soliton (a, theta, xi) on x < 0; psi = soliton with E -> z E plus the jump (2 pi/beta) lambda_p."""
    h = r + 1
    return field(r, soliton_terms(h, [(a, theta, xi)])), field(r, soliton_terms(h, [(a, theta, xi + np.log(z))]), lam(r, p))


# ======================================================================== 1. delays
worst = 0.0; fitdev = 0.0; worst_bad = np.inf
for r in (2, 3):
    h = r + 1
    for a in range(1, h):
        for typ in (1, -1):
            for p in range(h):
                for theta, thd in ((0.9, 0.3), (-0.4, 0.3), (1.2 + 0.3j, -0.5), (0.2, -0.8)):
                    xi = 0.2 + 0.6j
                    z = zdelay(a, h, theta, thd, typ, p)
                    F, G = pair(r, a, theta, xi, z, p)
                    worst = max(worst, resid(r, F, G, thd, typ))
                    # control: the class-0 delay in class p != 0 must fail
                    if p % h:
                        F, G = pair(r, a, theta, xi, zdelay(a, h, theta, thd, typ, 0), p)
                        worst_bad = min(worst_bad, resid(r, F, G, thd, typ))
# independent determination of z by least squares (type I_+, class 1, and type I_-, class 2 of a_3)
for r, a, typ, p in ((2, 1, -1, 1), (3, 2, -1, 0), (3, 1, 1, 2)):
    h = r + 1; theta, thd, xi = 0.7, -0.2, 0.3 + 0.5j

    def fun(v, r=r, a=a, typ=typ, p=p, theta=theta, thd=thd, xi=xi):
        F, G = pair(r, a, theta, xi, v[0] + 1j*v[1], p)
        out = []
        for t in np.linspace(-2, 2, 5):
            ph, phx, pht = F(0.0, t); ps, psx, pst = G(0.0, t)
            _, gp, gq, _ = defect_potential(ph, ps, np.exp(-thd), r, BETA, M, 'imag', typ)
            px, sx = defect_dx(pht, pst, gp, gq, r, typ)
            out += list(phx - px) + list(psx - sx)
        out = np.array(out)
        return np.r_[out.real, out.imag]
    z0 = zdelay(a, h, theta, thd, typ, p)
    sol = least_squares(fun, [z0.real + 0.3, z0.imag - 0.2], xtol=1e-15, ftol=1e-15, gtol=1e-15)
    fitdev = max(fitdev, abs(sol.x[0] + 1j*sol.x[1] - z0))
# time reversal: an I_-(sigma) solution, reversed in t, solves the I_+(1/sigma) conditions
tr = 0.0
for r, a in ((2, 1), (3, 2)):
    h = r + 1; theta, thd, xi = 0.8, 0.25, 0.1 + 0.4j
    F, G = pair(r, a, theta, xi, zdelay(a, h, theta, thd, 1))

    def rev(H):
        return lambda x, t: (lambda v: (v[0], v[1], -v[2]))(H(x, -t))
    tr = max(tr, resid(r, rev(F), rev(G), -thd, -1))
# coweight shift of B_d: class p acts as theta_d -> theta_d +- 2 pi i p/h; root shifts do nothing
cw = 0.0
for r in (2, 3, 4):
    h = r + 1; al = an_data(r)['alpha']
    for typ in (1, -1):
        for p in range(1, h):
            vals = []
            for _ in range(3):
                z1 = rng.normal(size=r) + 0.3j*rng.normal(size=r); z2 = rng.normal(size=r) + 0.3j*rng.normal(size=r)
                sig = 0.7*np.exp(0.2j)
                B1, g1, q1, _ = defect_potential(z1, z2 + (2*np.pi/BETA)*lam(r, p), sig, r, BETA, M, 'imag', typ)
                B2, g2, q2, _ = defect_potential(z1, z2, sig*np.exp(-2j*np.pi*p*typ/h), r, BETA, M, 'imag', typ)
                B3, g3, q3, _ = defect_potential(z1, z2 + (2*np.pi/BETA)*al[1], sig, r, BETA, M, 'imag', typ)
                B0, g0, q0, _ = defect_potential(z1, z2, sig, r, BETA, M, 'imag', typ)
                vals.append(B1 - B2)
                cw = max(cw, np.abs(g1 - g2).max(), np.abs(q1 - q2).max(), abs(B3 - B0), np.abs(g3 - g0).max())
            cw = max(cw, np.ptp(np.real(vals)), np.ptp(np.imag(vals)))
# dictionary for type I_-, class 0
dic = 0.0
for h in (3, 4, 5):
    for a in range(1, h):
        for theta, thd in ((0.4, -0.3), (1.1, 0.6)):
            g = np.pi*a/h; sig = np.exp(-thd)
            bcz = (1j*np.exp(-theta) + sig*np.exp(1j*g))/(1j*np.exp(-theta) + sig*np.exp(-1j*g))        # (6.2), theta -> -theta
            gb = np.pi*(h - a)/h
            cz7 = (np.exp(-(theta - thd)) + 1j*np.exp(-1j*gb))/(np.exp(-(theta - thd)) + 1j*np.exp(1j*gb))  # (2.20), a -> h - a
            dic = max(dic, abs(bcz - zdel(a, h, theta - thd)), abs(cz7 - zdel(a, h, theta - thd)))
report(f'delays (a_2, a_3, all species, both types, every class p): psi = soliton with E -> z E, '
       f'z = z_a(theta - theta_d - 2 pi i p/h) (I_-), z_a(theta_d - 2 pi i p/h - theta) (I_+), at 7 times [res {worst:.0e}; '
       f'class-0 delay in class p != 0: >= {worst_bad:.2f}]; least-squares z agrees to {fitdev:.0e}; time reversal '
       f'I_-(sigma) -> I_+(1/sigma) [{tr:.0e}]; B_d(phi, psi + 2 pi lambda_p/beta; sigma) = B_d(phi, psi; sigma e^(-+2 pi i p/h)) '
       f'+ const, root shifts invariant [{cw:.0e}]; = BCZ (6.2) at theta -> -theta = cz2007 (2.20) at a -> h - a [{dic:.0e}]',
       worst < 1e-11 and worst_bad > 1e-2 and fitdev < 1e-8 and tr < 1e-11 and cw < 1e-10 and dic < 1e-13)


# ======================================================================== 2. stored charge, selection rule
def charge(r, a, imxi):
    """topological charge (beta/2pi) phi(+inf) of the species-a soliton (sec-topological-charges)."""
    h = r + 1; al = an_data(r)['alpha']
    th = 2*np.pi*np.arange(h)*a/h + imxi
    th = (th + np.pi) % (2*np.pi) - np.pi
    return -(al.T @ th)/(2*np.pi)


ok_static = True; ok_cz = True; ok_shell = True; ok_real = True
for r in (2, 3, 4, 5):
    h = r + 1; Mh = 2*h*M/BETA**2
    for typ in (1, -1):
        thd = 0.45
        for p in range(h):
            ps = (2*np.pi/BETA)*lam(r, p)
            B, gp, gq, Pd = defect_potential(np.zeros(r), ps, np.exp(-thd), r, BETA, M, 'imag', typ)
            ok_static &= np.abs(gp).max() < 1e-12 and np.abs(gq).max() < 1e-12
            s = typ*2j*np.pi*p/h
            ok_static &= np.isclose(B, Mh*(np.cosh(thd) - np.cosh(thd + s))) and np.isclose(Pd, Mh*(np.sinh(thd) - np.sinh(thd + s)))
            # BCZ / CZ normalization (no -1's): E + const = -Mhat (cosh, sinh)(theta_d + s), on the mass shell Mhat
            E0, P0 = B - Mh*np.cosh(thd), Pd - Mh*np.sinh(thd)
            ok_shell &= np.isclose(E0**2 - P0**2, Mh**2) and np.isclose(E0, -Mh*np.cosh(thd + s))
            if typ == 1:     # cz2009a (4.8): (u, v) = 2pi(lambda_a, lambda_b)/beta, P_CZ = -P, here a = 0, b = p
                ok_cz &= np.isclose(E0, -Mh*np.cosh(thd - 2*(0 - p)*np.pi*1j/h)) and np.isclose(-P0, Mh*np.sinh(thd - 2*(0 - p)*np.pi*1j/h))
            ok_real &= (abs(np.imag(B)) < 1e-9) == ((2*p) % h == 0)
# charges of transmitted solitons: continuation in x (charge) and in t at x = 0 (jump)
ok_charge = True; jump_dev = 0.0; ok_lat = True


def unwrap_log(vals):
    return np.log(np.abs(vals)) + 1j*np.unwrap(np.angle(vals), axis=-1)


for r, a, theta, imxi in ((2, 1, 1.3, 0.4), (2, 2, 0.9, -1.1), (3, 1, 1.6, 2.0), (3, 2, 0.8, 0.5), (4, 2, 2.2, 0.3)):
    h = r + 1; al = an_data(r)['alpha']; thd = 0.2; xi = 0.1 + 1j*imxi
    z = zdelay(a, h, theta, thd)
    om = np.exp(2j*np.pi*np.arange(h)*a/h)
    E = lambda x, t, zz=1.0: zz*np.exp(mass(a, h)*(x*np.cosh(theta) - t*np.sinh(theta)) + xi)
    xs = np.linspace(-60, 60, 40001)/mass(a, h)
    Ui, Uo = np.imag(unwrap_log(1 + np.outer(om, E(xs, 0.0)))), np.imag(unwrap_log(1 + np.outer(om, E(xs, 0.0, z))))
    lin = -(al.T @ (Ui[:, -1] - Ui[:, 0]))/(2*np.pi)
    lout = -(al.T @ (Uo[:, -1] - Uo[:, 0]))/(2*np.pi)
    ok_charge &= np.allclose(lin, charge(r, a, imxi)) and np.allclose(lout, charge(r, a, imxi + np.angle(z)))
    T = 60/(mass(a, h)*abs(np.sinh(theta)))
    ts = np.linspace(-T, T, 40001)
    Lphi = unwrap_log(1 + np.outer(om, E(0.0, ts))); Lpsi = unwrap_log(1 + np.outer(om, E(0.0, ts, z)))
    dphi = (1j/BETA)*(al.T @ (Lphi[:, -1] - Lphi[:, 0])); dpsi = (1j/BETA)*(al.T @ (Lpsi[:, -1] - Lpsi[:, 0]))
    jump_dev = max(jump_dev, np.abs((BETA/(2*np.pi))*(dpsi - dphi) - (lin - lout)).max())
    for t in (-T, T):                     # the jump itself lies in the weight lattice
        F, G = pair(r, a, theta, xi, z)
        jp = (BETA/(2*np.pi))*(G(0.0, t)[0] - F(0.0, t)[0])
        c = al[1:] @ jp
        ok_lat &= np.allclose(c, np.round(c.real), atol=1e-8)
# selection rule
ok_arg = True; steps = {}; ok_pi = True; spacing_ok = True
for r in (2, 3, 4):
    h = r + 1; pim = an_data(r)['pi']
    for a in range(1, h):
        g = gcd(a, h); s = 2*np.pi*g/h; ncyc = h//g
        # boundaries between charge sectors in Im xi
        grid = np.linspace(0, 2*np.pi, 2001, endpoint=False)
        ch = np.array([charge(r, a, v) for v in grid])
        jumps = grid[1:][np.abs(np.diff(ch, axis=0)).max(axis=1) > 1e-9]
        spacing_ok &= np.allclose(np.diff(jumps), s, atol=2*np.pi/2000 + 1e-12) and len(jumps) in (ncyc - 1, ncyc)
        ths = np.linspace(-6, 6, 241)
        args = np.angle(zdel(a, h, ths))
        if 2*a < h:
            ok_arg &= np.all(np.diff(args) > 0) and args.min() > 0 and args.max() < 2*np.pi*a/h
        elif 2*a > h:
            ok_arg &= np.all(np.diff(args) < 0) and args.max() < 0 and args.min() > -2*np.pi*(h - a)/h
        ks = set()
        for th in ths[::8]:
            for imxi in np.linspace(0.013, 2*np.pi, 37):
                lout = charge(r, a, imxi + np.angle(zdel(a, h, th)))
                k = [j for j in range(ncyc) if np.allclose(charge(r, a, imxi + j*s), lout)]
                kk = k[0] if 2*a <= h else (ncyc - k[0]) % ncyc
                ks.add(kk)
        steps[(h, a)] = (max(ks), (a if 2*a <= h else h - a)//g)
        for imxi in (0.3, 1.7):
            ok_pi &= np.allclose(charge(r, a, imxi + 2*np.pi*a/h), pim @ charge(r, a, imxi))
ok_steps = all(v[0] == v[1] for v in steps.values())
report(f'stored charge: static coweight jumps solve the conditions; class p stores Mhat[(cosh, sinh)theta_d - (cosh, sinh)'
       f'(theta_d +- 2 pi i p/h)], with BCZ constants -Mhat (cosh, sinh)(theta_d +- 2 pi i p/h) on the mass shell Mhat '
       f'(= cz2009a (4.8)), real iff 2p = 0 mod h (a_2..a_5, both types); transmitted solitons: charges by x-continuation '
       f'= phase formula, jump change (t-continuation at x = 0) = lambda_in - lambda_out [{jump_dev:.0e}], jumps in the weight '
       f'lattice; arg z_a monotone in (0, 2 pi a/h) or (-2 pi abar/h, 0); sector boundaries 2 pi gcd(a,h)/h apart; max steps '
       f'along the cycle = a/gcd(a,h) {sorted(steps.items())}; theta >> theta_d: lambda_out = pi lambda_in',
       ok_static and ok_cz and ok_shell and ok_real and ok_charge and jump_dev < 1e-6 and ok_lat and ok_arg and ok_steps
       and ok_pi and spacing_ok)
report('cz2007 (2.22): the sector boundaries are 2 gamma_a apart only when a divides h (e.g. a_4, a = 2: spacing 2 pi/5, '
       'not 4 pi/5; a species-2 soliton can move two steps)', steps[(5, 2)] == (2, 2) and steps[(3, 1)] == (1, 1))


# ======================================================================== 3. absorption and emission
def theta_abs(a, h, p, thd, typ=1):
    return thd - typ*1j*np.pi/2 + typ*1j*np.pi*(a + 2*p)/h


ok_conf = True; worst_ae = 0.0; bal = 0.0
for r, a, p in ((2, 1, 0), (2, 2, 1), (3, 1, 0), (3, 2, 0), (3, 3, 2), (5, 1, 1), (5, 2, 3)):
    h = r + 1
    for typ in (1, -1):
        for thd in (0.6, -0.6):
            ta = theta_abs(a, h, p, thd, typ)
            u = (ta - thd - 2j*np.pi*p/h) if typ == 1 else (thd - 2j*np.pi*p/h - ta)
            ok_conf &= abs(np.exp(-u) - 1j*np.exp(-1j*np.pi*a/h)) < 1e-12           # pole of z
            xi = 0.2 + 0.5j
            F = field(r, soliton_terms(h, [(a, ta, xi)])); G = field(r, const_terms(h, a), lam(r, p))
            worst_ae = max(worst_ae, resid(r, F, G, thd, typ))                     # (soliton, const)
            F2 = field(r, [(np.ones(h, complex), 0j, 0j)]); G2 = field(r, soliton_terms(h, [(a, ta, xi)]), lam(r, p))
            worst_ae = max(worst_ae, resid(r, F2, G2, thd, typ))                   # (const, soliton)
            # stored energy and momentum: B_d(+T) - B_d(-T) = M_a (cosh, sinh) theta_abs (for theta_d > 0 the
            # soliton arrives from the left; the envelope moves with velocity tanh theta_d)
            if thd > 0:
                T = 40/(mass(a, h)*np.sinh(thd)*abs(np.cos(np.imag(ta))) + 1e-300)
                Bm, _, _, Pm = defect_potential(F(0.0, -T)[0], G(0.0, -T)[0], np.exp(-thd), r, BETA, M, 'imag', typ)
                Bp, _, _, Pp = defect_potential(F(0.0, T)[0], G(0.0, T)[0], np.exp(-thd), r, BETA, M, 'imag', typ)
                # orientation: Re cosh theta_abs has the sign of sin(pi(a + 2p)/h).  If it is positive the object is a
                # species-a soliton (energy M_a cosh theta_abs) and the class goes p -> p + a; if negative the profile
                # is inverted, the object is the species-abar soliton at theta_abs - i pi (energy -M_a cosh theta_abs),
                # and the class goes p + a -> p.
                Ma = 2*h*mass(a, h)/BETA**2; sg = np.sign(np.sin(np.pi*(a + 2*p)/h))
                ok_conf &= np.sign(np.real(np.cosh(ta))) == sg
                bal = max(bal, abs((Bp - Bm) - sg*Ma*np.cosh(ta))/Ma, abs((Pp - Pm) - sg*Ma*np.sinh(ta))/Ma)
                Mh = 2*h*M/BETA**2
                cls = lambda q: Mh*(np.cosh(thd) - np.cosh(thd + typ*2j*np.pi*q/h))
                pin, pout = (p, p + a) if sg > 0 else (p + a, p)
                bal = max(bal, abs(Bm - cls(pin))/Ma, abs(Bp - cls(pout))/Ma)
            # the zero of the class-(p+a) delay sits at the same rapidity
            ok_conf &= abs(zdelay(a, h, ta, thd, typ, p + a)) < 1e-12
# envelope velocity
ok_vel = True
for (a, h, p, thd) in ((1, 3, 0, 0.7), (2, 5, 1, -0.4)):
    ta = theta_abs(a, h, p, thd)
    E = lambda x, t: np.exp(mass(a, h)*(x*np.cosh(ta) - t*np.sinh(ta)))
    # the curve |E| = 1 moves with velocity tanh theta_d
    ok_vel &= abs(abs(E(np.tanh(thd)*3.0, 3.0)) - 1) < 1e-12
# real-rapidity absorption rule
ok_rule = True
for h in range(2, 9):
    for a in range(1, h):
        for p in range(h):
            im = np.imag(theta_abs(a, h, p, 0.3)) % (2*np.pi)
            real = min(im, 2*np.pi - im) < 1e-12
            ok_rule &= real == ((a + 2*p - h/2) % (2*h) == 0)
ok_rule &= (theta_abs(1, 6, 1, 0.3) == 0.3) and np.isclose(theta_abs(2, 4, 0, 0.3), 0.3)


# quadrature: E + B_d and P + P_d through real-rapidity absorptions
def bulk_totals(r, F, t, L=40.0, n=48001):
    xs = np.linspace(-L, 0.0, n)
    ph, phx, pht = F(xs, t)
    V = sum((M**2/(1j*BETA)**2)*(np.exp(1j*BETA*(al_ @ ph)) - 1) for al_ in an_data(r)['alpha'])
    e = 0.5*(pht*pht).sum(0) + 0.5*(phx*phx).sum(0) + V
    pdens = -(pht*phx).sum(0)
    return simpson(e, x=xs), simpson(pdens, x=xs)


cons = 0.0
for r, a, p in ((3, 2, 0), (5, 1, 1)):
    h = r + 1; thd = 0.5; ta = theta_abs(a, h, p, thd).real
    F = field(r, soliton_terms(h, [(a, ta, 0.3j)])); G = field(r, const_terms(h, a), lam(r, p))
    tot = []
    for t in np.linspace(-30, 15, 10):
        Eb, Pb = bulk_totals(r, F, t)
        B, _, _, Pd = defect_potential(F(0.0, t)[0], G(0.0, t)[0], np.exp(-thd), r, BETA, M, 'imag', 1)
        tot.append((Eb + B, Pb + Pd))
    tot = np.array(tot); Ma = 2*h*mass(a, h)/BETA**2
    cons = max(cons, np.abs(tot - tot[0]).max()/Ma)
# the theta_d -> -infinity limit of type I_-: z_a -> omega^a, i.e. tau'_j = tau_{j+1}, psi = pi phi
ok_twist = True
for r in (2, 3):
    h = r + 1; pim = an_data(r)['pi']
    for a in range(1, h):
        ok_twist &= abs(zdel(a, h, 25.0) - np.exp(2j*np.pi*a/h)) < 1e-10 and abs(zdel(a, h, -25.0) - 1) < 1e-10
        F = field(r, soliton_terms(h, [(a, 0.4, 0.2 + 0.3j)]))
        G = field(r, soliton_terms(h, [(a, 0.4, 0.2 + 0.3j + 2j*np.pi*a/h)]))
        for x, t in ((0.3, 0.2), (-1.1, 0.7)):
            ok_twist &= np.allclose(G(x, t)[0], pim @ F(x, t)[0])
report(f'absorption/emission at theta_abs = theta_d -+ i pi/2 +- i pi(a + 2p)/h (types I_-+, theta_d = +-0.6; a_2, a_3, a_5): '
       f'pole of z; (soliton, const) and (const, soliton) solve the conditions [res {worst_ae:.0e}]; Re cosh theta_abs has the '
       f'sign of sin(pi(a + 2p)/h); if positive class p -> p + a with (Delta B_d, Delta P_d) = M_a (cosh, sinh) theta_abs, if '
       f'negative the abar soliton at theta_abs + i pi, classes p + a -> p [rel {bal:.0e}]; same rapidity = zero of the '
       f'class-(p+a) delay; '
       f'envelope velocity tanh theta_d; theta_abs real iff a + 2p = h/2 mod 2h (h <= 8); E + B_d, P + P_d constant by '
       f'quadrature through the real absorptions a_3 (a = 2, p = 0), a_5 (a = 1, p = 1) [rel {cons:.0e}]; '
       f'theta - theta_d -> +inf: z_a -> omega^a, psi = pi phi (the rotation line)',
       ok_conf and worst_ae < 1e-11 and bal < 1e-9 and ok_vel and ok_rule and cons < 1e-6 and ok_twist)

# ======================================================================== 4. several solitons
w2 = 0.0
for r, (a1, a2), p in ((2, (1, 2), 0), (2, (1, 1), 1), (3, (1, 2), 1), (3, (2, 3), 0)):
    h = r + 1
    for typ in (1, -1):
        thd = 0.3
        s1, s2 = (a1, 1.1, 0.2 + 0.7j), (a2, -0.5, -0.4 + 0.3j)
        z1, z2 = zdelay(a1, h, s1[1], thd, typ, p), zdelay(a2, h, s2[1], thd, typ, p)
        F = field(r, soliton_terms(h, [s1, s2]))
        G = field(r, soliton_terms(h, [(a1, s1[1], s1[2] + np.log(z1)), (a2, s2[1], s2[2] + np.log(z2))]), lam(r, p))
        w2 = max(w2, resid(r, F, G, thd, typ, ts=np.linspace(-4, 4, 9)))
# two defects in a row: I_- at x = 0 (theta_d1), I_+ at x = L (theta_d2): the delays multiply
wdd = 0.0
for r, a in ((2, 1), (3, 2)):
    h = r + 1; th, xi, t1, t2, L = 0.8, 0.1 + 0.4j, 0.2, -0.6, 1.7
    z1 = zdelay(a, h, th, t1, 1); z2 = zdelay(a, h, th, t2, -1)
    F = field(r, soliton_terms(h, [(a, th, xi)]))
    C = field(r, soliton_terms(h, [(a, th, xi + np.log(z1))]))
    G = field(r, soliton_terms(h, [(a, th, xi + np.log(z1*z2))]))
    wdd = max(wdd, resid(r, F, C, t1, 1), resid(r, C, G, t2, -1, x=L))
# classical bootstrap
boot = 0.0
for h in range(3, 9):
    for a in range(1, h):
        for b in range(1, h - a):
            for u in (0.4 + 0.1j, -1.3, 2.0 - 0.5j):
                boot = max(boot, abs(zdel(a, h, u + 1j*np.pi*b/h)*zdel(b, h, u - 1j*np.pi*a/h) - zdel(a + b, h, u)))
report(f'several solitons: two-soliton solutions (a_2: (1,2), (1,1); a_3: (1,2), (2,3); both types; classes 0, 1) pass with '
       f'E_k -> z(theta_k) E_k and A_ab unchanged [res {w2:.0e}]; two defects (I_- then I_+) multiply the delays [{wdd:.0e}]; '
       f'z_a(u + i pi b/h) z_b(u - i pi a/h) = z_(a+b)(u) (cz2007 (2.25)), h = 3..8 [{boot:.0e}]',
       w2 < 1e-10 and wdd < 1e-11 and boot < 1e-12)

# ======================================================================== 5. reality
ok_th = True
for h in (3, 4, 5, 6):
    for a in range(1, h):
        for p in range(h):
            for theta, thd in ((0.4 + 0.2j, -0.3), (1.1, 0.5 + 0.3j)):
                for typ in (1, -1):
                    ok_th &= np.isclose(zdelay(h - a, h, np.conj(theta), np.conj(thd), typ, (h - p) % h),
                                        np.conj(zdelay(a, h, theta, thd, typ, p)))
# Theta-invariant a_2 breather through a real-sigma defect
r, h = 2, 3; th1, xi1, thd = 0.3 + 0.8j, 0.2 + 0.4j, 0.5
zb1, zb2 = zdelay(1, h, th1, thd), zdelay(2, h, np.conj(th1), thd)
ok_th &= np.isclose(zb2, np.conj(zb1))
F = field(r, soliton_terms(h, [(1, th1, xi1), (2, np.conj(th1), np.conj(xi1))]))
G = field(r, soliton_terms(h, [(1, th1, xi1 + np.log(zb1)), (2, np.conj(th1), np.conj(xi1) + np.log(zb2))]))
res_b = resid(r, F, G, thd, 1, ts=np.linspace(-3, 3, 7) + 0.137)
# Theta-image of a class-p solution solves the conditions in class -p (real sigma); class-(h-p) energy = conjugate
res_t = 0.0; ok_conj = True
for r, a, p in ((2, 1, 1), (3, 1, 1), (3, 2, 3)):
    h = r + 1; theta, thd, xi = 0.9, 0.25, 0.3 + 0.2j
    F, G = pair(r, a, theta, xi, zdelay(a, h, theta, thd, 1, p), p)
    Th = lambda H: (lambda x, t: tuple(-np.conj(v) for v in H(x, t)))
    res_t = max(res_t, resid(r, Th(F), Th(G), thd, 1))
    B1 = defect_potential(np.zeros(r), (2*np.pi/BETA)*lam(r, p), np.exp(-thd), r, BETA, M, 'imag', 1)[0]
    B2 = defect_potential(np.zeros(r), (2*np.pi/BETA)*lam(r, h - p), np.exp(-thd), r, BETA, M, 'imag', 1)[0]
    ok_conj &= np.isclose(B2, np.conj(B1))
# real rapidities: |z| > 1 (a delay) iff cos(2 pi p/h) > 0; in classes with cos < 0 solitons are advanced
ok_adv = True
for h in (3, 4, 5, 6, 7):
    for a in range(1, h):
        for p in range(h):
            if abs(np.cos(2*np.pi*p/h)) < 1e-12:
                continue
            for th in (-2.0, -0.3, 0.4, 1.9):
                ok_adv &= (abs(zdelay(a, h, th, 0.1, 1, p)) > 1) == (np.cos(2*np.pi*p/h) > 0)
                ok_adv &= (abs(zdelay(a, h, th, 0.1, -1, p)) > 1) == (np.cos(2*np.pi*p/h) > 0)
# h even, p = h/2: sigma -> -sigma
ok_odd = all(np.isclose(np.exp(-(0.3 + 2j*np.pi*(h//2)/h)), -np.exp(-0.3)) for h in (2, 4, 6))
report(f'reality: z_abar(conj theta; conj theta_d) = conj z_a in class h - p; a Theta-invariant a_2 breather passes a '
       f'real-sigma defect with conjugate delays [res {res_b:.0e}]; Theta-image of a class-p solution solves the conditions '
       f'(class -p) [{res_t:.0e}]; class h - p stores the conjugate energy; h even, p = h/2 is sigma -> -sigma; real '
       f'rapidities: delayed (|z| > 1) iff cos(2 pi p/h) > 0, advanced otherwise',
       ok_th and res_b < 1e-10 and res_t < 1e-11 and ok_conj and ok_odd and ok_adv)


# ======================================================================== 6. linear waves
def hessian(r, sigma, typ, m=M, Dpot=None):
    """Quadratic part of B_d about (0, 0) (both couplings): blocks (H_pp, H_pq, H_qq); Dpot is the D in f_j, g_j."""
    D0 = dmatrix(r, 1)[0] if Dpot is None else Dpot
    s = sigma if typ == 1 else 1/sigma
    Hpp = Hpq = Hqq = 0
    for al_ in an_data(r)['alpha']:
        u1, u2 = D0 @ al_, D0.T @ al_
        Hpp = Hpp + (m/4)*(s*np.outer(u1, u1) + np.outer(u2, u2)/s)
        Hpq = Hpq + (m/4)*(s*np.outer(u1, u2) - np.outer(u2, u2)/s)
        Hqq = Hqq + (m/4)*(s*np.outer(u2, u2) + np.outer(u2, u2)/s)
    return Hpp, Hpq, Hqq


def modes(r, rho=None):
    d = an_data(r); h = r + 1
    if rho is None:
        ev, V = np.linalg.eig(d['pi'])
        ks = np.round(np.angle(ev)*h/(2*np.pi)).astype(int) % h
        rho = {k: V[:, i] for i, k in enumerate(ks)}
    return rho


def linear_T(r, k, theta, thd, typ=1, side='L', m=M, D=None, E=None, H=None, rho=None):
    """Transmission/reflection of a plane wave of particle k (mode rho_k); side 'L': rapidity theta from the left,
    'R': rapidity -theta from the right.  Returns (R coefficients, T coefficients) in the mode basis."""
    h = r + 1
    if D is None:
        D, E = dmatrix(r, typ)
    Hpp, Hpq, Hqq = hessian(r, np.exp(-thd), typ, m) if H is None else H
    rho = modes(r, rho)
    ks = sorted(rho)
    w = mass(k, h, m)*np.cosh(theta); pk = mass(k, h, m)*np.sinh(theta)
    pj = {j: np.sqrt(complex(w**2 - mass(j, h, m)**2)) for j in ks}
    pj = {j: (v if v.imag > 0 or (abs(v.imag) < 1e-14 and v.real*np.real(pk) >= 0) else -v) for j, v in pj.items()}
    pj[k] = pk
    n = len(ks); Rm = np.array([rho[j] for j in ks]).T
    # unknowns: c (reflected, on the incident side), d (transmitted)
    rows = []; rhs = []
    sgn = 1 if side == 'L' else -1
    inc = rho[k]
    P = np.diag([pj[j] for j in ks])
    if side == 'L':      # phi = inc e^{i pk x} + sum c_j rho_j e^{-i p_j x};  psi = sum d_j rho_j e^{i p_j x}
        phi0 = lambda c, d: inc + Rm @ c;  phix = lambda c, d: 1j*pk*inc - 1j*Rm @ (P @ c)
        psi0 = lambda c, d: Rm @ d;        psix = lambda c, d: 1j*Rm @ (P @ d)
    else:                # psi = inc e^{-i pk x} + sum c_j rho_j e^{i p_j x};  phi = sum d_j rho_j e^{-i p_j x}
        psi0 = lambda c, d: inc + Rm @ c;  psix = lambda c, d: -1j*pk*inc + 1j*Rm @ (P @ c)
        phi0 = lambda c, d: Rm @ d;        phix = lambda c, d: -1j*Rm @ (P @ d)

    def eqs(c, d):
        ph, ps = phi0(c, d), psi0(c, d)
        pht, pst = -1j*w*ph, -1j*w*ps
        e1 = phix(c, d) - (E @ pht + D @ pst - (Hpp @ ph + Hpq @ ps))
        e2 = psix(c, d) - (D.T @ pht - E @ pst + (Hpq.T @ ph + Hqq @ ps))
        return np.r_[e1, e2]
    z0 = np.zeros(n, complex)
    b = -eqs(z0, z0)
    A = np.array([eqs(*np.split(np.eye(2*n)[i], 2)) - eqs(z0, z0) for i in range(2*n)]).T
    sol = np.linalg.solve(A, b)
    c, d = sol[:n], sol[n:]
    return {j: c[i] for i, j in enumerate(ks)}, {j: d[i] for i, j in enumerate(ks)}


def T_closed(k, h, theta, thd, typ=1):
    return zdel(k, h, theta - thd + 1j*np.pi/2) if typ == 1 else zdel(k, h, thd - theta - 1j*np.pi/2)


# Hessian against finite differences of defect_potential
hd = 0.0
for r in (2, 3):
    for typ in (1, -1):
        Hpp, Hpq, Hqq = hessian(r, 0.7, typ)
        e = 1e-5
        for i in range(r):
            dv = np.zeros(r); dv[i] = e
            _, gp1, gq1, _ = defect_potential(dv, np.zeros(r), 0.7, r, BETA, M, 'real', typ)
            _, gp2, gq2, _ = defect_potential(-dv, np.zeros(r), 0.7, r, BETA, M, 'real', typ)
            hd = max(hd, np.abs((gp1 - gp2)/(2*e) - Hpp[:, i]).max(), np.abs((gq1 - gq2)/(2*e) - Hpq.T[:, i]).max())
            _, gp1, gq1, _ = defect_potential(np.zeros(r), dv, 0.7, r, BETA, M, 'real', typ)
            _, gp2, gq2, _ = defect_potential(np.zeros(r), -dv, 0.7, r, BETA, M, 'real', typ)
            hd = max(hd, np.abs((gp1 - gp2)/(2*e) - Hpq[:, i]).max(), np.abs((gq1 - gq2)/(2*e) - Hqq[:, i]).max())
lin = 0.0; refl = 0.0; mix = 0.0; unim = 0.0; cross = 0.0
for r in (1, 2, 3, 4):
    h = r + 1
    for typ in (1, -1):
        for theta, thd in ((0.4, -0.3), (1.3, 0.5), (0.7, 0.7)):
            for k in range(1, h):
                Rc, Tc = linear_T(r, k, theta, thd, typ)
                refl = max(refl, max(abs(v) for v in Rc.values()))
                mix = max([mix] + [abs(v) for j, v in Tc.items() if j != k])
                lin = max(lin, abs(Tc[k] - T_closed(k, h, theta, thd, typ)))
                unim = max(unim, abs(abs(Tc[k]) - 1))
                Rb, Tb = linear_T(r, k, theta, thd, typ, side='R')
                refl = max(refl, max(abs(v) for v in Rb.values()))
                cross = max(cross, abs(Tb[k] - T_closed(h - k, h, 1j*np.pi - theta, thd, typ)),
                            abs(T_closed(k, h, theta, thd, typ)*T_closed(h - k, h, 1j*np.pi + theta, thd, typ) - 1))
# value at theta = theta_d, pole at theta_d - i(pi - pi k/h), zero at theta_d + i(pi - pi k/h) (I_-)
pz = 0.0
for h in (2, 3, 4, 5):
    for k in range(1, h):
        g = np.pi*k/h; thd = 0.3
        tp = thd - 1j*(np.pi - g)              # the denominator e^{-(theta - theta_d)} + e^{-i pi k/h} vanishes there
        pz = max(pz, abs(T_closed(k, h, thd, thd, -1) + np.exp(1j*g)), abs(T_closed(k, h, thd + 1j*g, thd, -1)),
                 1e-6*abs(1/T_closed(k, h, thd - 1j*g + 1e-9, thd, -1)))   # I_+: -e^{i g}, zero at +i g, pole at -i g
        pz = max(pz, abs(T_closed(k, h, thd, thd) - np.exp(1j*g)), abs(T_closed(k, h, thd + 1j*(np.pi - g), thd)),
                 abs(np.exp(-(tp - thd)) + np.exp(-1j*g)), abs(np.exp(-(tp - thd)) + np.exp(1j*g)) < 1e-3)
# a_1: eq-shg-linear-transmission; limits
u = 0.37
a1 = abs(T_closed(1, 2, u + 0.2, 0.2) - (-1j*np.sinh(u/2 - 1j*np.pi/4)/np.sinh(u/2 + 1j*np.pi/4)))
lim = max(abs(T_closed(k, h, 0.3, 30.0) - 1) + abs(T_closed(k, h, 0.3, -30.0) - np.exp(2j*np.pi*k/h))
          for h in (3, 4, 5) for k in range(1, h))
# cz2007 conventions: D -> D^T, E -> 1 - D^T, potential with D^T, rho_k = -sum alpha_l omega^{lk}, m = 1
czv = 0.0; czprint = np.inf; cz75 = 0.0
for r in (2, 3):
    h = r + 1; D0 = dmatrix(r, 1)[0]; al = an_data(r)['alpha']
    Dc = D0.T; Ec = np.eye(r) - Dc
    rho_cz = {k: -sum(al[l]*np.exp(2j*np.pi*l*k/h) for l in range(h)) for k in range(1, h)}
    for theta, eta in ((0.6, 0.1), (1.4, -0.5)):
        sig = np.exp(-eta)
        H = hessian(r, sig, 1, m=1.0, Dpot=Dc)
        for k in range(1, h):
            _, Tc = linear_T(r, k, theta, eta, 1, m=1.0, D=Dc, E=Ec, H=H, rho=rho_cz)
            mk = 2*np.sin(np.pi*k/h); om = np.exp(2j*np.pi/h)
            corr = (1j*np.exp(-theta)*mk - sig*(1 - om**(-k)))/(1j*np.exp(-theta)*mk + sig*(1 - om**k))
            prnt = (1j*np.exp(-theta)*mk - sig*(1 - om**(-k)))/(1j*np.exp(theta)*mk + sig*(1 - om**k))
            czv = max(czv, abs(Tc[k] - corr), abs(Tc[k] - T_closed(h - k, h, theta, eta)))
            czprint = min(czprint, abs(Tc[k] - prnt))
            if h == 3:
                v = (theta - eta)/2
                ref = (np.exp(-1j*np.pi/3)*np.sinh(v - 1j*np.pi/6)/np.sinh(v + 1j*np.pi/6) if k == 1 else
                       np.exp(1j*np.pi/3)*np.sinh(v + 2j*np.pi/3)/np.sinh(v + 1j*np.pi/3))
                cz75 = max(cz75, abs(Tc[k] - ref))
# about a class-p defect at imaginary coupling (numerical Hessian at (0, 2 pi lambda_p/beta))
chg = 0.0; ok_mod = True
for r, p in ((2, 1), (3, 1), (3, 2)):
    h = r + 1; thd = 0.35; sig = np.exp(-thd); e = 1e-5; q0 = (2*np.pi/BETA)*lam(r, p)
    Hpp = np.zeros((r, r), complex); Hpq = np.zeros((r, r), complex); Hqq = np.zeros((r, r), complex)
    for i in range(r):
        dv = np.zeros(r); dv[i] = e
        _, gp1, gq1, _ = defect_potential(dv, q0, sig, r, BETA, M, 'imag', 1)
        _, gp2, gq2, _ = defect_potential(-dv, q0, sig, r, BETA, M, 'imag', 1)
        Hpp[:, i] = (gp1 - gp2)/(2*e)
        _, gp1, gq1, _ = defect_potential(np.zeros(r), q0 + dv, sig, r, BETA, M, 'imag', 1)
        _, gp2, gq2, _ = defect_potential(np.zeros(r), q0 - dv, sig, r, BETA, M, 'imag', 1)
        Hpq[:, i] = (gp1 - gp2)/(2*e); Hqq[:, i] = (gq1 - gq2)/(2*e)
    for k in range(1, h):
        _, Tc = linear_T(r, k, 0.8, thd, 1, H=(Hpp, Hpq, Hqq))
        chg = max(chg, abs(Tc[k] - T_closed(k, h, 0.8, thd + 2j*np.pi*p/h)))
        ok_mod &= (abs(abs(Tc[k]) - 1) < 1e-6) == ((2*p) % h == 0)
report(f'linear waves (real coupling; a_1..a_4, both types, all particles): R = 0, no mixing, |T| = 1 '
       f'[{max(refl, mix, unim):.0e}]; T^(I-)_k = (e^-u + e^(i pi k/h))/(e^-u + e^(-i pi k/h)) = z_k(u + i pi/2), '
       f'T^(I+)_k(theta) = z_k(theta_d - theta - i pi/2) [{lin:.0e}]; from the right Tbar_k(theta) = T_kbar(i pi - theta), '
       f'T_k(theta) T_kbar(theta + i pi) = 1 [{cross:.0e}]; a_1 = eq-shg-linear-transmission [{a1:.0e}]; limits 1 '
       f'(theta_d -> +inf) and omega^k (theta_d -> -inf) [{lim:.0e}]; T_k(theta_d) = e^(i pi k/h), pole at theta_d - '
       f'i(pi - pi k/h), zero at theta_d + i(pi - pi k/h) (I_-), -e^(i pi k/h), poles/zeros at theta_d -+ i pi k/h (I_+) '
       f'[{pz:.0e}]; Hessian vs finite differences [{hd:.0e}]',
       max(refl, mix, unim) < 1e-10 and lin < 1e-10 and cross < 1e-10 and a1 < 1e-14 and lim < 1e-10 and hd < 1e-7
       and pz < 1e-12)
report(f'cz2007 in their conventions (D^T, rho_k, m = 1): their linearized system gives the corrected (7.11) with e^-theta '
       f'in the denominator [{czv:.0e}], = book T_(h-k), = (7.5) for a_2 [{cz75:.0e}]; the printed (7.11) misses by '
       f'>= {czprint:.2f}; about a class-p defect (imaginary coupling) T_k at theta_d + 2 pi i p/h [{chg:.0e}], unimodular '
       f'only if 2p = 0 mod h', czv < 1e-10 and cz75 < 1e-10 and czprint > 1e-2 and chg < 1e-6 and ok_mod)

# ======================================================================== 7. frozen soliton (pre-registered)
# (i) the continuation u -> u + i pi/2 relates soliton factors and linear factors, in the bulk and at the defect
cont = 0.0
for h in (3, 4, 5, 6):
    for a in range(1, h):
        for b in range(1, h):
            for u in (0.3, -1.2, 0.8 + 0.2j):
                X = (1j*np.sinh(u) - np.cos(np.pi*(a - b)/h))/(1j*np.sinh(u) - np.cos(np.pi*(a + b)/h))   # eq-an-transmission
                cont = max(cont, abs(X - Aint(a, b, h, u + 1j*np.pi/2)))
        for u in (0.3, -1.2):
            cont = max(cont, abs(T_closed(a, h, u + 0.4, 0.4) - zdel(a, h, u + 1j*np.pi/2)))


# (ii) same-type pairs: zero/pole matching.  An I_- factor z_b(u + i alpha) has (zero, pole) at angles
# (pi/2 + g_b + alpha, pi/2 - g_b + alpha) on |e^-u| = 1, difference +2 g_b; an I_+ factor difference -2 g_b.
# A_cb has zeros at +-(g_c - g_b) and poles at +-(g_c + g_b) in e^-u (all on the unit circle, so alpha is real).
def match(h, b, c, signs):
    gb, gc = np.pi*b/h, np.pi*c/h
    Zt, Pt = [gc - gb, -(gc - gb)], [gc + gb, -(gc + gb)]
    found = False
    for Zs, Ps in ((Zt, Pt), (Pt, Zt)):                  # A or 1/A
        for perm in ((0, 1), (1, 0)):
            diffs = [Zs[i] - Ps[perm[i]] for i in range(2)]
            want = [s*2*gb for s in signs]
            for wperm in ((0, 1), (1, 0)):
                if all(abs(np.exp(1j*(diffs[i] - want[wperm[i]])) - 1) < 1e-12 for i in range(2)):
                    found = True
    return found


ok_same = True; ok_mixed = True
for h in range(2, 9):
    for b in range(1, h):
        for c in range(1, h):
            ok_same &= match(h, b, c, (1, 1)) == (2*b == h) and match(h, b, c, (-1, -1)) == (2*b == h)
            ok_mixed &= match(h, b, c, (1, -1))
# grid search confirmation for a_2, b = c = 1, same type: no constant ratio
h, b, c = 3, 1, 1
us = np.linspace(-3, 3, 60)                 # avoids the zero of A_11 at u = 0
best = np.inf
for a1_ in np.linspace(0, 2*np.pi, 90, endpoint=False):
    for a2_ in np.linspace(0, 2*np.pi, 90, endpoint=False):
        rat = zdel(b, h, us + 1j*a1_)*zdel(b, h, us + 1j*a2_)/Aint(c, b, h, us)
        dev = np.abs(rat/rat[0] - 1).max()
        if np.isfinite(dev):
            best = min(best, dev)
# (iii) the mixed pair, as a function and on fields, and its energy-momentum with BCZ's constants
mixed = 0.0; en = 0.0
for h in range(2, 9):
    for b in range(1, h):
        for c in range(1, h):
            thc = 0.37; vm = thc + 1j*(np.pi/2 + np.pi*c/h); vp = np.conj(vm)
            for th in (-1.2, 0.3, 0.9 + 0.2j, 2.1):
                mixed = max(mixed, abs(zdelay(b, h, th, vm, 1)*zdelay(b, h, th, vp, -1)
                                       - np.exp(2j*np.pi*b/h)*Aint(c, b, h, th - thc)))
                X = lambda z: (z - np.cos(np.pi*(c - b)/h))/(z - np.cos(np.pi*(c + b)/h))
                mixed = max(mixed, abs(T_closed(b, h, th, vm, 1)*T_closed(b, h, th, vp, -1)
                                       - np.exp(2j*np.pi*b/h)*X(1j*np.sinh(th - thc))))
        Mh = 2*h*M/BETA**2; Mc = 2*h*mass(c, h)/BETA**2; thc = 0.37; vm = thc + 1j*(np.pi/2 + np.pi*c/h)
        # BCZ normalization: each defect carries -Mhat (cosh, sinh)(effective parameter)
        en = max(en, abs(-Mh*(np.cosh(vm) + np.cosh(np.conj(vm))) - Mc*np.cosh(thc))/Mc,
                 abs(-Mh*(np.sinh(vm) + np.sinh(np.conj(vm))) - Mc*np.sinh(thc))/Mc)
# realization at real sigma through classes, and the field check (a_3, c = 2, p = 2; a_5, c = 1, p = 2)
real_cls = 0.0; fld = 0.0
for r, c in ((1, 1), (3, 2), (5, 1), (5, 3)):
    h = r + 1
    assert (h + 2*c) % 4 == 0
    p = ((h + 2*c)//4) % h
    thc = 0.37; Mh = 2*h*M/BETA**2; Mc = 2*h*mass(c, h)/BETA**2
    for b in range(1, h):
        for th in (0.8, -0.4):
            real_cls = max(real_cls, abs(zdelay(b, h, th, thc, 1, p)*zdelay(b, h, th, thc, -1, p)
                                         - np.exp(2j*np.pi*b/h)*Aint(c, b, h, th - thc)))
    q0 = (2*np.pi/BETA)*lam(r, p)
    Bm, _, _, Pm = defect_potential(np.zeros(r), q0, np.exp(-thc), r, BETA, M, 'imag', 1)
    Bp, _, _, Pp = defect_potential(np.zeros(r), q0, np.exp(-thc), r, BETA, M, 'imag', -1)
    real_cls = max(real_cls, abs(Bm + Bp - 2*Mh*np.cosh(thc) - Mc*np.cosh(thc))/Mc,
                   abs(Pm + Pp - 2*Mh*np.sinh(thc) - Mc*np.sinh(thc))/Mc)
    # fields: soliton b through I_- (class p) at x = 0 and I_+ (class p) at x = L: total factor omega^b A_cb
    if r > 1:
        b = 1; th, xi, L = 0.9, 0.2 + 0.3j, 1.3
        z1 = zdelay(b, h, th, thc, 1, p); z2 = zdelay(b, h, th, thc, -1, p)
        F = field(r, soliton_terms(h, [(b, th, xi)]))
        C = field(r, soliton_terms(h, [(b, th, xi + np.log(z1))]), lam(r, p))
        G = field(r, soliton_terms(h, [(b, th, xi + np.log(z1*z2))]), lam(r, 2*p))
        fld = max(fld, resid(r, F, C, thc, 1), resid(r, C, G, thc, -1, x=L),
                  abs(z1*z2 - np.exp(2j*np.pi*b/h)*Aint(c, b, h, th - thc)))
report(f'frozen soliton: X_ab(i sinh u) = A_ab(u + i pi/2) and T_b(u) = z_b(u + i pi/2) [{cont:.0e}]; two defects of '
       f'the same type give kappa A_cb^(+-1) iff b = h/2 (zero/pole matching, h = 2..8; a_2 grid search: best deviation '
       f'{best:.2f}); the mixed pair I_-(theta_c + i(pi/2 + pi c/h)) I_+(theta_c - i(pi/2 + pi c/h)) gives omega^b A_cb and '
       f'omega^b X_cb for all b, c [{mixed:.0e}] and carries M_c (cosh, sinh) theta_c with BCZ constants [{en:.0e}]; at real '
       f'sigma via classes 4p = h + 2c mod 4h (a_1 c = 1, a_3 c = 2, a_5 c = 1, 3) [{real_cls:.0e}], on fields [{fld:.0e}]',
       cont < 1e-12 and ok_same and ok_mixed and 0.05 < best < np.inf and mixed < 1e-12 and en < 1e-12 and real_cls < 1e-10
       and fld < 1e-10)
print(f'   runtime {time.time() - T_START:.0f} s')
