"""Sine-Gordon on the half-line: soliton reflection matrix in book conventions (Part VII, ch. 29/32).

Part VI conventions (boundary_common.py): coproduct Delta(e)=e(x)1+k(x)e, physical q = -exp(-i pi lambda),
principal gradation y = exp(lambda theta); spin-1/2 basis (soliton, antisoliton).  Coideal generators
b_j = e_j + q^-1 f_j k_j + eps_j k_j.  GZ = Ghoshal-Zamolodchikov hep-th/9306002, BPT = Bajnok-Palla-
Takacs-Toth hep-th/0106070 (their (3.2) parametrization of GZ).

Checks:
 1. K: V(y) -> V(1/y) from the coideal: unique; solves the reflection equation (braid form) with the
    book's R-check; the reflection equation is invariant under q -> -q (S_T -> -S_T) for the same K,
    and K(eps; -q) = K(-eps; q): the sign of q is invisible in the reflection equation.
 2. Dictionary to GZ (5.12): eps_1 = eps e^{i xi}, eps_0 = eps e^{-i xi}, k = 1/(eps sin(pi lambda));
    to BPT (3.2): eps_1 = -cos(eta - i vartheta)/sin(pi lambda), eps_0 = -cos(eta + i vartheta)/sin(pi lambda);
    consistency with GZ (5.25).  Neumann (GZ xi=0, k=1/sin(pi lambda/2); BPT eta=pi(lambda+1)/2) is
    eps_0 = eps_1 = 1/(2cos(pi lambda/2)) = (2-q-q^-1)^(-1/2), the a_n special value; eps = 0 is eta = pi/2.
 3. Boundary crossing-unitarity in book form  k(theta) = S(2 theta) k(-theta),
    k = sum_ab K_{b abar}(i pi/2 - theta) v_a (x) v_b, S = S_0 R-check (eq-sg-s0, eq-sg-smatrix):
    the matrix part reduces it to a scalar equation at the physical q with C = sigma_x; at -q only with
    C' = [[0,1],[-1,0]] (the bulk kappa = -1 conjugation of sec-crossing-sign).
 4. Scalar factor: GZ's a(u) = sin(lambda(pi-u)) rho(u) equals the book S_0; GZ/BPT R_0 and sigma
    satisfy their functional equations; the full BPT K satisfies boundary unitarity and boundary
    crossing-unitarity with the book S_0 (numerically, two parameter points).
 5. Breather reflection factors: the B_1, B_2 amplitudes obtained from the soliton K by the boundary
    bootstrap equal Ghoshal (3.5)-(3.10); B_1 factor unitary and crossing-unitary with S_11 (eq-sg-s11).
 6. Continuation lambda -> -2/B: Ghoshal B_1 factor = Corrigan-Delius (2.5) with E = B eta/pi, F = i B vt/pi;
    Dirichlet and Neumann values.
 7. Excited boundary states (Bajnok-Palla-Takacs-Toth): bootstrap identity for Q on |0> (eta -> etabar),
    poles nu_0, nu_1 of sigma(eta,u), K on |0> (s <-> sbar, eta -> etabar) unitary and crossing-unitary.
 8. Limits of the UV-IR relation (BPTT sec. 4; in Al. B. Zamolodchikov's normalization BPT hep-th/0108157 (3.2), (3.3)):
    Neumann, Dirichlet, the eps-hat = 0 point, and the classical value of M_crit (the C = 1 boundary).
 9. Poles (sec-bsg-k-matrix): R_0 has simple poles at u = n pi/(2 lambda), none at pi/2 - n pi/(2 lambda); at these
    points Q is regular and P has the pole, except for vt = 0 and odd n; Neumann tower: P and Q have simple poles at
    u = pi/2 - pi/lambda and u = pi/2.
10. GZ (5.22)/(5.24) normalization (notes): with sigma normalized as GZ (5.24) prints it, the prefactor 1/cos(xi) of
    (5.22) violates (5.20) by k^2/cos^2 xi; the prefactor 1/k (equivalently, the BPTT normalization) satisfies it.
11. prp-bsg-breather-k: the residue of R-check(theta_1 - theta_2) at theta_12 = i u_n is rank 1, with image the breather
    singlet in V(theta - i u_n/2) (x) V(theta + i u_n/2).
12. BPTT (3.5): a_n = a(u+nu_n)a(u-nu_n)/(a(u+nu_0)a(u-nu_0)) = prod_{l=1}^n {2(eta/pi - l)} for n = 1, 2, 3; the bootstrap
    a(u - nu_n) Q(eta) b(u + nu_n) = Q(etabar) a_n holds for n = 0, 1, 2 (eq-bsg-excited-k, prp-bsg-excited).
13. BPTT (3.9): B_1, B_2 on |0>, |1> by the full matrix bootstrap equal K_{B_m}(etabar) b_mk,
    b_mk = a_k(u + pi/2 - m pi/2lam) a_k(u - pi/2 + m pi/2lam); and (exr-bsg-sign-q) Pi K(eps) Pi ~ K(-eps).
Runtime: about 50 seconds.
"""
import numpy as np
import mpmath as mp
from boundary_common import report, q_phys, sl2_spin, rcheck, k_intertwiner, re_preserving

mp.mp.dps = 20
LAM = 1.37                         # lambda = 8 pi/beta_SG^2 - 1 (attractive regime, two breathers)
q = q_phys(LAM)
C = np.array([[0, 1], [1, 0]])     # charge conjugation s <-> sbar


def Rsg(z, qq=q):
    return rcheck(sl2_spin(1, qq, z), sl2_spin(1, qq, 1.0))[0]


def K_coideal(y, e0, e1, qq=q):
    null, _ = k_intertwiner(sl2_spin(1, qq, y), sl2_spin(1, qq, 1/y), [1/qq, 1/qq], [e0, e1])
    assert len(null) == 1
    return null[0]


def K_closed(y, e0, e1, qq=q):
    d = y**2 - y**-2
    return np.array([[(qq - 1/qq)*(e1*y + e0/y), d], [d, (qq - 1/qq)*(e0*y + e1/y)]])/d


def prop(A, B):
    c = np.vdot(B.ravel(), A.ravel())/np.vdot(B.ravel(), B.ravel())
    return np.linalg.norm(A - c*B)/np.linalg.norm(A)


# ---------------------------------------------------------------- 1. K from the coideal, reflection equation
e0, e1 = 0.43 - 0.2j, -0.71 + 0.35j
report('K(theta) from the coideal equals the closed form (unique)',
       max(prop(K_coideal(y, e0, e1), K_closed(y, e0, e1)) for y in (0.8+0.3j, 1.9)) < 1e-12)
r, c = re_preserving(Rsg, lambda y: K_coideal(y, e0, e1), 1.4 + 0.3j, 0.6 - 0.2j)
report(f'reflection equation (braid form) at the physical q [{r:.0e}]', r < 1e-10)
r2, _ = re_preserving(lambda z: Rsg(z, -q), lambda y: K_coideal(y, e0, e1), 1.4 + 0.3j, 0.6 - 0.2j)
report(f'the same K also solves it with the R-check at -q [{r2:.0e}]: the RE does not see the sign of q', r2 < 1e-10)
report('K(eps_0, eps_1; -q) = K(-eps_0, -eps_1; q)',
       prop(K_coideal(0.8+0.3j, e0, e1, -q), K_coideal(0.8+0.3j, -e0, -e1, q)) < 1e-12)
R_q, R_mq = Rsg(1.7+0.2j), Rsg(1.7+0.2j, -q)
Pi = np.diag([1, -1])
report('R-check(-q) = (Pi (x) 1) R-check(q) (Pi (x) 1), Pi = diag(1,-1): S_T -> -S_T',
       np.allclose(R_mq, np.kron(Pi, np.eye(2)) @ R_q @ np.kron(Pi, np.eye(2))))

# ---------------------------------------------------------------- 2. dictionary to GZ and BPT
ok = True
for xi, kk in ((0.4, 0.7), (1.1 + 0.2j, 2.3)):
    eps = 1/(kk*np.sin(np.pi*LAM))
    for th in (0.37, -0.81 + 0.1j):
        u = -1j*th; y = np.exp(LAM*th)
        GZ = np.array([[np.cos(xi + LAM*u), kk/2*np.sin(2*LAM*u)], [kk/2*np.sin(2*LAM*u), np.cos(xi - LAM*u)]])
        ok &= prop(K_closed(y, eps*np.exp(-1j*xi), eps*np.exp(1j*xi)), GZ) < 1e-12
report('GZ (5.12): K = R(u)[[cos(xi+lambda u), k/2 sin 2lambda u],[.., cos(xi-lambda u)]] with '
       'eps_1 = eps e^{i xi}, eps_0 = eps e^{-i xi}, k = 1/(eps sin pi lambda)', ok)
ok = True
for eta, vt in ((0.6, 0.3), (2.0, 1.1)):
    E1 = -np.cos(eta - 1j*vt)/np.sin(np.pi*LAM); E0 = -np.cos(eta + 1j*vt)/np.sin(np.pi*LAM)
    for th in (0.37, -0.81 + 0.1j):
        u = -1j*th; y = np.exp(LAM*th)
        Pp = np.cos(LAM*u)*np.cos(eta)*np.cosh(vt) - np.sin(LAM*u)*np.sin(eta)*np.sinh(vt)
        Pm = np.cos(LAM*u)*np.cos(eta)*np.cosh(vt) + np.sin(LAM*u)*np.sin(eta)*np.sinh(vt)
        Q0 = -np.sin(LAM*u)*np.cos(LAM*u)
        ok &= prop(K_closed(y, E0, E1), np.array([[Pp, Q0], [Q0, Pm]])) < 1e-12
    # GZ (5.25): cos eta cosh vt = -cos(xi)/k, cos^2 eta + cosh^2 vt = 1 + 1/k^2
    eps = np.sqrt(E0*E1); xi = -1j*np.log(E1/eps); kk = 1/(eps*np.sin(np.pi*LAM))
    ok &= abs(np.cos(eta)*np.cosh(vt) + np.cos(xi)/kk) < 1e-12 and abs(np.cos(eta)**2 + np.cosh(vt)**2 - 1 - 1/kk**2) < 1e-12
report('BPT (3.2): eps_1 = -cos(eta - i vt)/sin(pi lambda), eps_0 = -cos(eta + i vt)/sin(pi lambda); '
       'consistent with GZ (5.25)', ok)
eps_N = 1/(2*np.cos(np.pi*LAM/2))
report('Neumann (GZ xi=0, k=1/sin(pi lambda/2); BPT eta=pi(lambda+1)/2, vt=0): eps_0 = eps_1 = 1/(2cos(pi lambda/2)) '
       '= (2-q-1/q)^(-1/2)',
       abs(1/((1/np.sin(np.pi*LAM/2))*np.sin(np.pi*LAM)) - eps_N) < 1e-14 and
       abs(-np.cos(np.pi*(LAM+1)/2)/np.sin(np.pi*LAM) - eps_N) < 1e-14 and abs(eps_N**2*(2 - q - 1/q) - 1) < 1e-12)
eps_bpt = lambda eta, vt: (-np.cos(eta + 1j*vt)/np.sin(np.pi*LAM), -np.cos(eta - 1j*vt)/np.sin(np.pi*LAM))
report('dictionary: eps_0 = eps_1 = 0 <-> BPT (eta, vt) = (pi/2, 0); BPT (0, 0) <-> eps_0 = eps_1 = -1/sin(pi lambda)',
       np.allclose(eps_bpt(np.pi/2, 0), 0) and np.allclose(eps_bpt(0, 0), -1/np.sin(np.pi*LAM)))


# ---------------------------------------------------------------- 3./4. scalar factor and crossing-unitarity
def lg(z):
    return mp.loggamma(z)


L = mp.mpf(LAM)


def log_rho(u):
    """GZ (5.7)."""
    u = mp.mpmathify(u)
    def logF(l, uu):
        v = L*uu/mp.pi
        return lg(2*l*L - v) + lg(1 + 2*l*L - v) - lg((2*l+1)*L - v) - lg(1 + (2*l-1)*L - v)
    s = mp.nsum(lambda l: logF(l, u) + logF(l, mp.pi - u) - logF(l, 0) - logF(l, mp.pi), [1, mp.inf])
    return mp.log(-1/mp.pi) + lg(L) + lg(1 - L*u/mp.pi) + lg(1 - L + L*u/mp.pi) + s


def a_gz(theta):
    u = -1j*mp.mpmathify(theta)
    return mp.sin(L*(mp.pi - u))*mp.exp(log_rho(u))


def S0_book(theta):
    """eq-sg-s0, xi_SG = pi/lambda, real theta."""
    xs = mp.pi/L
    f = lambda t: mp.sinh((mp.pi - xs)*t/2)*mp.sin(theta*t)/(t*mp.sinh(xs*t/2)*mp.cosh(mp.pi*t/2))
    return -mp.exp(-1j*mp.quad(f, [0, 1, 5, 20, mp.inf]))


def log_R0(u):
    """BPT (3.2) R_0 (= GZ (5.21))."""
    u = mp.mpmathify(u)
    def t(l, uu):
        v = 2*L*uu/mp.pi
        return lg(4*l*L - v) + lg(4*L*(l-1) + 1 - v) - lg((4*l-3)*L - v) - lg((4*l-1)*L + 1 - v)
    return mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])


def log_sigma(x, u):
    """BPT form of GZ's sigma: sigma(x,u) sigma(x,-u) = cos^2 x/(cos(x+lambda u)cos(x-lambda u))."""
    x, u = mp.mpmathify(x), mp.mpmathify(u)
    def t(l, uu):
        v = L*uu/mp.pi; a = x/mp.pi
        return (lg(0.5 + a + (2*l-1)*L - v) + lg(0.5 - a + (2*l-1)*L - v)
                - lg(0.5 - a + (2*l-2)*L - v) - lg(0.5 + a + 2*l*L - v))
    return mp.log(mp.cos(x)/mp.cos(x + L*u)) + mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])


# GZ (5.23) as printed: the normalizer Pi^2(x, pi/2) Pi^2(-x, pi/2) diverges, log Pi(x, pi/2) ~ (lambda/4) log(l_max)
def logPi_partial(x, u, lmax):
    return sum(lg(0.5 + (2*l+0.5)*L + x/mp.pi - L*u/mp.pi) + lg(0.5 + (2*l+1.5)*L + x/mp.pi)
               - lg(0.5 + (2*l+1.5)*L + x/mp.pi - L*u/mp.pi) - lg(0.5 + (2*l+0.5)*L + x/mp.pi) for l in range(lmax))
d1 = mp.re(logPi_partial(0.7, mp.pi/2, 400) - logPi_partial(0.7, mp.pi/2, 100))
report(f'GZ (5.23) normalizer diverges: log Pi(x,pi/2) grows by (lambda/4) log 4 per factor 4 in l_max '
       f'[{float(d1):.4f} vs {float(L/4*mp.log(4)):.4f}]; BPT form used instead', abs(d1 - L/4*mp.log(4)) < 1e-3)
th_list = [0.23, 0.61]
dev = max(abs(a_gz(th) - S0_book(th)) for th in th_list)
report(f"GZ a(theta) = sin(lambda(pi-u)) rho(u) equals the book S_0 (eq-sg-s0) [{float(dev):.0e}]", dev < 1e-10)
x0 = mp.mpf('0.7')
dev1 = max(abs(mp.exp(log_sigma(x0, u) + log_sigma(x0, -u))*mp.cos(x0 + L*u)*mp.cos(x0 - L*u)/mp.cos(x0)**2 - 1)
           for u in (0.3j, 0.2 + 0.1j))
dev2 = max(abs(mp.exp(log_sigma(x0, mp.pi/2 - u) - log_sigma(x0, mp.pi/2 + u)) - 1) for u in (0.3j, 0.2 + 0.1j))
report(f'sigma(x,u)sigma(x,-u) = cos^2x/(cos(x+lambda u)cos(x-lambda u)) and sigma(x,pi/2-u) = sigma(x,pi/2+u) '
       f'(GZ (5.24), BPT normalization) [{float(max(dev1, dev2)):.0e}]', max(dev1, dev2) < 1e-10)


def K_full(theta, eta, vt):
    """BPT (3.2): K = R_0 sigma(eta,u) sigma(i vt,u)/(cos eta cosh vt) [[P+, Q0],[Q0, P-]], u = -i theta."""
    u = -1j*mp.mpmathify(theta)
    P0 = lambda s: mp.cos(L*u)*mp.cos(eta)*mp.cosh(vt) - s*mp.sin(L*u)*mp.sin(eta)*mp.sinh(vt)
    Q0 = -mp.sin(L*u)*mp.cos(L*u)
    pref = mp.exp(log_R0(u) + log_sigma(eta, u) + log_sigma(1j*vt, u))/(mp.cos(eta)*mp.cosh(vt))
    M = np.array([[complex(P0(1)), complex(Q0)], [complex(Q0), complex(P0(-1))]])
    return complex(pref)*M


def kvec(Kmat):
    """k[(a,b)] = K_{b abar}: K^{ab} of GZ (3.32)-(3.33) with the plain conjugation abar."""
    Kc = Kmat @ C            # (K C)_{b a} = K_{b abar}
    return Kc.T.reshape(4)   # index a*2 + b


# matrix part of crossing-unitarity: at physical q it reduces to a scalar with the plain conjugation
# C = sigma_x (the C for which the bulk S is crossing symmetric, kappa = +1); at -q it fails with sigma_x
# and holds only with C' = [[0,1],[-1,0]], the conjugation of the kappa = -1 bulk amplitude (sec-crossing-sign)
def cu_res(qq, Cm):
    res = []
    for eta, vt in ((0.6, 0.3), (2.0, 1.1)):
        E1 = -np.cos(eta - 1j*vt)/np.sin(np.pi*LAM); E0 = -np.cos(eta + 1j*vt)/np.sin(np.pi*LAM)
        for th in (0.29, 0.77):
            km = (K_closed(np.exp(LAM*(1j*np.pi/2 - th)), E0, E1) @ Cm).T.reshape(4)
            kp = (K_closed(np.exp(LAM*(1j*np.pi/2 + th)), E0, E1) @ Cm).T.reshape(4)
            v = Rsg(np.exp(2*LAM*th), qq) @ kp
            res.append(np.linalg.norm(np.outer(km, v) - np.outer(v, km))/(np.linalg.norm(km)*np.linalg.norm(v)))
    return max(res)


Cp = np.array([[0, 1], [-1, 0]])
r_q, r_mq, r_mq2 = cu_res(q, C), cu_res(-q, C), cu_res(-q, Cp)
report(f'crossing-unitarity k(theta) = f(theta) R(2theta) k(-theta): matrix part consistent at q with C = sigma_x '
       f'[{r_q:.0e}], inconsistent at -q [{r_mq:.2f}], restored at -q only by C = [[0,1],[-1,0]] [{r_mq2:.0e}]',
       r_q < 1e-12 and r_mq > 1e-3 and r_mq2 < 1e-12)

# full scalar factor: boundary unitarity and crossing-unitarity with the book S_0
devU, devC = 0.0, 0.0
for eta, vt in ((mp.mpf('0.6'), mp.mpf('0.3')), (mp.pi*(L+1)/2, mp.mpf(0))):     # generic; Neumann
    for th in (0.29,):
        devU = max(devU, np.abs(K_full(th, eta, vt) @ K_full(-th, eta, vt) - np.eye(2)).max())
        km = kvec(K_full(1j*np.pi/2 - th, eta, vt)); kp = kvec(K_full(1j*np.pi/2 + th, eta, vt))
        S2 = complex(S0_book(2*th))*Rsg(np.exp(2*LAM*th))
        devC = max(devC, np.abs(km - S2 @ kp).max()/np.abs(km).max())
report(f'boundary unitarity K(theta)K(-theta) = 1 with R_0 sigma(eta) sigma(i vt)/(cos eta cosh vt) [{devU:.0e}]',
       devU < 1e-9)
report(f'boundary crossing-unitarity k(theta) = S_0(2theta) R(2theta) k(-theta) with the book S_0 [{devC:.0e}] '
       '(generic point and Neumann)', devC < 1e-9)

# ---------------------------------------------------------------- 5. breather reflection factors by the bootstrap (ch. 29)
# The breather B_n sits in V(theta - i u_n/2) (x) V(theta + i u_n/2) (u_n = pi - n pi/lambda) as the U_q-invariant
# vector w (in Part VI conventions the singlet exists in this ordering only).  Reflecting the two constituents,
# (1 (x) K(theta_a)) S(theta_a + theta_b) (1 (x) K(theta_b)) w = R_B(theta) w, with the full scalar factors.
from boundary_common import cop


def inv_vec(ya, yb):
    A, B = sl2_spin(1, q, ya), sl2_spin(1, q, yb)
    M = np.vstack([cop(A, B, kd, j) - (np.eye(4) if kd == 'k' else 0) for j in (0, 1) for kd in 'efk'])
    _, s, Vh = np.linalg.svd(M)
    w = Vh[-1].conj()
    return s[-1]/s[0], w/w[1]                       # normalize the s(x)sbar component to 1


def ghoshal_breather(n, u, eta, vt, lam=None):
    """Ghoshal hep-th/9310188 (3.5)-(3.10) for n = 1, 2."""
    lam = L if lam is None else lam
    pre = ((-1)**(n+1)*mp.cos(u/2 + n*mp.pi/(4*lam))*mp.cos(u/2 - mp.pi/4 - n*mp.pi/(4*lam))*mp.sin(u/2 + mp.pi/4)
           / (mp.cos(u/2 - n*mp.pi/(4*lam))*mp.cos(u/2 + mp.pi/4 + n*mp.pi/(4*lam))*mp.sin(u/2 - mp.pi/4)))
    for l in range(1, n):
        pre *= (mp.sin(u + l*mp.pi/(2*lam))*mp.cos(u/2 - mp.pi/4 - l*mp.pi/(4*lam))**2
                / (mp.sin(u - l*mp.pi/(2*lam))*mp.cos(u/2 + mp.pi/4 + l*mp.pi/(4*lam))**2))

    def Sx(x):
        if n == 1:
            return (mp.cos(x/lam) - mp.sin(u))/(mp.cos(x/lam) + mp.sin(u))
        c1, c2 = mp.cos(x/lam - mp.pi/(2*lam)), mp.cos(x/lam + mp.pi/(2*lam))
        return (mp.sin(u) - c1)*(mp.sin(u) - c2)/((mp.sin(u) + c1)*(mp.sin(u) + c2))
    return pre*Sx(eta)*Sx(1j*vt)


eta, vt = mp.mpf('0.6'), mp.mpf('0.3')
devs, ok_order = [], True
for n in (1, 2):
    un = np.pi - n*np.pi/LAM
    for th in (0.3, -0.55):
        ta, tb = th - 1j*un/2, th + 1j*un/2
        r_in, w = inv_vec(np.exp(LAM*ta), np.exp(LAM*tb))
        ok_order &= r_in < 1e-12 and inv_vec(np.exp(LAM*tb), np.exp(LAM*ta))[0] > 1e-3
        S2 = complex(a_gz(ta + tb))*Rsg(np.exp(LAM*(ta + tb)))
        v = np.kron(np.eye(2), K_full(ta, eta, vt)) @ S2 @ np.kron(np.eye(2), K_full(tb, eta, vt)) @ w
        RB = v[1]
        devs.append(np.linalg.norm(v - RB*w)/np.linalg.norm(v))
        devs.append(abs(RB/complex(ghoshal_breather(n, -1j*mp.mpf(th), eta, vt)) - 1))
report(f'breather bootstrap from the BPT soliton K: B_1, B_2 reflect diagonally with exactly Ghoshal (3.5)-(3.10) '
       f'[{max(devs):.0e}]; singlet in V(theta - i u_n/2) (x) V(theta + i u_n/2) only', max(devs) < 1e-10 and ok_order)
S11 = lambda th: (mp.sinh(th) + 1j*mp.sin(mp.pi/L))/(mp.sinh(th) - 1j*mp.sin(mp.pi/L))     # eq-sg-s11, xi = pi/lambda
R1B = lambda th: ghoshal_breather(1, -1j*th, eta, vt)
pts = [mp.mpc('0.31', '0.07'), mp.mpc('-0.6', '0.2')]
dev = max(max(abs(R1B(t)*R1B(-t) - 1), abs(R1B(1j*mp.pi/2 - t) - S11(2*t)*R1B(1j*mp.pi/2 + t))) for t in pts)
report(f'B_1 reflection factor: unitarity and crossing-unitarity with S_11 = f_(xi/pi) [{float(dev):.0e}]', dev < 1e-20)

# ---------------------------------------------------------------- 6. continuation to sinh-Gordon (Corrigan-Delius (2.5))
Bs = mp.mpf('0.71')
lam_shg = -2/Bs                                         # lambda -> -2/B under beta_SG^2 -> -beta_SG^2
cdb = lambda x, th: mp.sinh(th/2 + 1j*mp.pi*x/4)/mp.sinh(th/2 - 1j*mp.pi*x/4)        # CD (2.3) = eq-blocks, h = 2
dev = 0
for eta_, vt_ in ((mp.mpf('0.9'), mp.mpf('0.4')), (mp.mpf('2.3'), mp.mpf('0'))):
    E, F = Bs*eta_/mp.pi, 1j*Bs*vt_/mp.pi
    for th in (mp.mpc('0.4', '0.1'), mp.mpc('-1.1', '0.3')):
        cd = (cdb(1, th)*cdb(2 - Bs/2, th)*cdb(1 + Bs/2, th)
              / (cdb(1 - E, th)*cdb(1 + E, th)*cdb(1 - F, th)*cdb(1 + F, th)))
        dev = max(dev, abs(ghoshal_breather(1, -1j*th, eta_, vt_, lam_shg)/cd - 1))
report(f'Ghoshal B_1 factor at lambda = -2/B equals Corrigan-Delius (2.5) with E = B eta/pi, F = i B vt/pi [{float(dev):.0e}]',
       dev < 1e-20)
th = mp.mpc('0.4', '0.1')
KD_shg = cdb(2 - Bs/2, th)*cdb(1 + Bs/2, th)/cdb(1, th)
lim = ghoshal_breather(1, -1j*th, mp.mpf(0), mp.mpf(60), lam_shg)
EN = Bs*(mp.pi*(lam_shg + 1)/2)/mp.pi
report(f'Dirichlet phi_0 = 0 (eta = 0, vt -> infinity) gives K_D = (2-B/2)(1+B/2)/(1) [{float(abs(lim/KD_shg - 1)):.0e}]; '
       f'Neumann eta = pi(lambda+1)/2 gives E = B/2 - 1, i.e. |E| = 1 - B/2', abs(lim/KD_shg - 1) < 1e-8 and
       abs(EN - (Bs/2 - 1)) < 1e-25)

# ---------------------------------------------------------------- 7. excited boundary states (BPTT hep-th/0106070)
L = mp.mpf('3.3')                                       # lambda large enough for two boundary states nu_0, nu_1


eta = mp.mpf('5.5'); etab = mp.pi*(L + 1) - eta
nu = lambda n: eta/L - (2*n + 1)*mp.pi/(2*L)
a_u = lambda u: a_gz(1j*u)
b_u = lambda u: mp.sin(L*u)/mp.sin(L*(mp.pi - u))*a_gz(1j*u)
dev = 0
for u in (mp.mpc('0.3', '0.2'), mp.mpc('-0.5', '0.1')):
    lhs = a_u(u - nu(0))*b_u(u + nu(0))
    rhs = mp.exp(log_sigma(etab, u) - log_sigma(eta, u))*mp.cos(eta)/mp.cos(etab)
    dev = max(dev, abs(lhs/rhs - 1))
ok_ph = 0 < nu(1) < nu(0) < mp.pi/2 and nu(2) < 0 and eta <= mp.pi*(L + 1)/2
report(f'BPTT bootstrap on the first excited state: a(u - nu_0) Q(eta, vt, u) b(u + nu_0) = Q(etabar, vt, u), '
       f'etabar = pi(lambda+1) - eta [{float(dev):.0e}]; nu_n = eta/lambda - (2n+1)pi/(2lambda)', dev < 1e-15 and ok_ph)
# pole of the ground-state P^+ at u = nu_n (from sigma(eta, u)); simple pole: (u - nu) sigma finite and nonzero
for n in (0, 1):
    d = mp.mpf('1e-8')
    r1 = abs(mp.exp(log_sigma(eta, nu(n) + d))*d); r2 = abs(mp.exp(log_sigma(eta, nu(n) + 2*d))*2*d)
    ok_ph &= abs(r1/r2 - 1) < 1e-5 and r1 > 1e-12
report('sigma(eta, u) has simple poles at u = nu_0, nu_1 in the physical strip (boundary bound states |0>, |1>)', ok_ph)
L = mp.mpf(LAM)
Kx = np.array([[0, 1], [1, 0]])
eta0, vt0 = mp.mpf('2.9'), mp.mpf('0.3')
etab0 = mp.pi*(L + 1) - eta0
K0 = lambda th: Kx @ K_full(th, etab0, vt0) @ Kx       # K on |0>: s <-> sbar and eta -> etabar (BPTT (3.6))
th = 0.29
devU = np.abs(K0(th) @ K0(-th) - np.eye(2)).max()
km = kvec(K0(1j*np.pi/2 - th)); kp = kvec(K0(1j*np.pi/2 + th))
devC = np.abs(km - complex(S0_book(2*th))*Rsg(np.exp(2*LAM*th)) @ kp).max()/np.abs(km).max()
report(f'K on the first excited boundary |0> satisfies unitarity [{devU:.0e}] and crossing-unitarity [{devC:.0e}]',
       devU < 1e-9 and devC < 1e-9)

# ---------------------------------------------------------------- 8. limits of Al. B. Zamolodchikov's UV-IR relation (BPT hep-th/0108157 (3.2))
# cos(eta/(lambda+1)) cosh(vt/(lambda+1)) = r cos(b phi0/2),  sin(eta/(lambda+1)) sinh(vt/(lambda+1)) = r sin(b phi0/2),
# r = M_0/M_crit, beta_SG^2/8pi = 1/(lambda+1), M_crit = sqrt(2 mu/sin(beta_SG^2/8)).
def uvir_lhs(eta, vt, lam=LAM):
    return mp.cos(eta/(lam+1))*mp.cosh(vt/(lam+1)), mp.sin(eta/(lam+1))*mp.sinh(vt/(lam+1))


ok = True
c, s_ = uvir_lhs(mp.pi*(LAM + 1)/2, mp.mpf(0))
ok &= abs(c) < 1e-20 and abs(s_) < 1e-20                                               # Neumann <-> M_0 = 0
bphi = mp.mpf('0.4'); V = mp.mpf(40)
c, s_ = uvir_lhs((LAM + 1)*bphi, (LAM + 1)*V)
ok &= abs(s_/c - mp.tan(bphi)) < 1e-20 and c > 1e15                                     # Dirichlet: r -> oo, eta = 4pi phi_0/beta
c, s_ = uvir_lhs(mp.pi/2, mp.mpf(0))
ok &= abs(c - mp.cos(mp.pi/(2*(LAM + 1)))) < 1e-20 and abs(s_) < 1e-20                  # eps-hat = 0 point, phi_0 = 0
bs2 = mp.mpf('1e-4'); m0 = mp.mpf(1); mu = m0**2/bs2                                     # classical limit: mu = m_0^2/beta^2
Mcrit = mp.sqrt(2*mu/mp.sin(bs2/8))
ok &= abs(Mcrit/(4*m0/bs2) - 1) < 1e-8                                                  # = M_0 of the C = 1 boundary
report('Zamolodchikov UV-IR relation: M_0 -> 0 gives Neumann (eta, vt) = (pi(lambda+1)/2, 0); M_0 -> infinity gives '
       'the Dirichlet eta = 4 pi phi_0/beta_SG; (pi/2, 0) <-> M_0/M_crit = cos(pi/(2(lambda+1))); classically '
       'M_crit -> 4 m_0/beta_SG^2, the C_0 = C_1 = 1 boundary', ok)

# ---------------------------------------------------------------- 9. poles of R_0, P and Q (sec-bsg-k-matrix, "Poles")
# R_0 has simple poles at u = n pi/(2 lambda) (soliton + image -> B_n at zero momentum; GZ after (5.21), BPTT sec. 3.4),
# not at pi/2 - n pi/(2 lambda) (that is R_0^(n) of the breathers).  sin(lambda u) in Q_0 removes them from Q; for
# vt = 0, cos(lambda u) removes the odd ones from P.  Neumann tower: P and Q have poles at u = pi/2 - n pi/lambda.
L = mp.mpf('3.3')


def order(vals):
    """pole order p from |f(u0 + d)|/|f(u0 + 2d)| = 2^p; vals = (f(u0 + d), f(u0 + 2d))."""
    return float(mp.log(abs(vals[0])/abs(vals[1]))/mp.log(2))


def PQ(u, eta, vt):
    pref = mp.exp(log_R0(u) + log_sigma(eta, u) + log_sigma(1j*vt, u))
    return (pref*(mp.cos(L*u)*mp.cos(eta)*mp.cosh(vt) - mp.sin(L*u)*mp.sin(eta)*mp.sinh(vt)),
            -pref*mp.sin(L*u)*mp.cos(L*u))


def orders_PQ(u0, eta, vt, d=mp.mpf('1e-7')):
    v1, v2 = PQ(u0 + d, eta, vt), PQ(u0 + 2*d, eta, vt)
    return order((v1[0], v2[0])), order((v1[1], v2[1]))


d7 = mp.mpf('1e-7')
R0f = lambda u: mp.exp(log_R0(u))
o_new = [order((R0f(n*mp.pi/(2*L) + d7), R0f(n*mp.pi/(2*L) + 2*d7))) for n in (1, 2, 3)]
o_old = [order((R0f(mp.pi/2 - n*mp.pi/(2*L) + d7), R0f(mp.pi/2 - n*mp.pi/(2*L) + 2*d7))) for n in (1, 2, 3)]
report(f'R_0 has simple poles at u = n pi/(2 lambda), n = 1, 2, 3 < lambda = 3.3 [orders {[round(o, 3) for o in o_new]}], '
       f'and none at pi/2 - n pi/(2 lambda) [{[round(o, 3) for o in o_old]}]',
       all(abs(o - 1) < 1e-3 for o in o_new) and all(abs(o) < 1e-3 for o in o_old))
gen = [orders_PQ(n*mp.pi/(2*L), mp.mpf('2.1'), mp.mpf('0.3')) for n in (1, 2, 3)]
vt0 = [orders_PQ(n*mp.pi/(2*L), mp.mpf('2.1'), mp.mpf(0)) for n in (1, 2, 3)]
ok_gen = all(abs(p - 1) < 1e-3 and abs(qq_) < 1e-3 for p, qq_ in gen)
ok_vt0 = all(abs(vt0[i][0] - (1 if n % 2 == 0 else 0)) < 1e-3 and abs(vt0[i][1]) < 1e-3 for i, n in enumerate((1, 2, 3)))
report('at u = n pi/(2 lambda): Q is regular, P has the simple pole for generic (eta, vt) = (2.1, 0.3); for vt = 0 '
       'P keeps it only for even n (n = 1, 2, 3)', ok_gen and ok_vt0)
etaN = mp.pi*(L + 1)/2
oN = [orders_PQ(u0, etaN, mp.mpf(0)) for u0 in (mp.pi/2 - mp.pi/L, mp.pi/2)]
report(f'Neumann point: P and Q have simple poles at u = pi/2 - pi/lambda (energy M sin(pi/lambda)) and at u = pi/2 '
       f'[orders {[(round(p, 3), round(qq_, 3)) for p, qq_ in oN]}]', all(abs(p - 1) < 1e-3 and abs(qq_ - 1) < 1e-3 for p, qq_ in oN))

# ---------------------------------------------------------------- 10. GZ (5.22)/(5.24) normalization (notes)
# With sigma normalized as GZ (5.24) prints it, sigma(x,u) sigma(x,-u) = 1/(cos(x+lambda u) cos(x-lambda u)), GZ's
# prefactor 1/cos(xi) in (5.22) violates the unitarity condition (5.20) by k^2/cos^2 xi; the prefactor 1/k works,
# i.e. (5.22) is right with the BPTT normalization cos^2 x/(cos(x+lambda u) cos(x-lambda u)).
L = mp.mpf(LAM)
eta, vt = mp.mpf('0.6'), mp.mpf('0.3')
E1 = -mp.cos(eta - 1j*vt)/mp.sin(mp.pi*L); E0 = -mp.cos(eta + 1j*vt)/mp.sin(mp.pi*L)
epsB = mp.sqrt(E0*E1); xiB = -1j*mp.log(E1/epsB); kB = 1/(epsB*mp.sin(mp.pi*L))
u = mp.mpc('0.3', '0.1')
sGZ = lambda x, uu: mp.exp(log_sigma(x, uu))/mp.cos(x)
prodGZ = sGZ(eta, u)*sGZ(eta, -u)*sGZ(1j*vt, u)*sGZ(1j*vt, -u)
s2 = mp.sin(L*u)**2
D520 = mp.cos(xiB)**2 - (1 + kB**2)*s2 + kB**2*s2**2
Mgz = lambda uu: mp.matrix([[mp.cos(xiB + L*uu), kB/2*mp.sin(2*L*uu)], [kB/2*mp.sin(2*L*uu), mp.cos(xiB - L*uu)]])
MM = Mgz(u)*Mgz(-u)
devM = max(abs(MM[0, 0] - D520), abs(MM[1, 1] - D520), abs(MM[0, 1]), abs(MM[1, 0]))
r_printed = prodGZ*D520/mp.cos(xiB)**2
r_k = prodGZ*D520/kB**2
report(f'GZ (5.22) with sigma normalized by (5.24): R(u)R(-u) M(u)M(-u) = k^2/cos^2 xi = {mp.nstr(mp.re(kB**2/mp.cos(xiB)**2), 8)} '
       f'[{float(abs(r_printed - kB**2/mp.cos(xiB)**2)):.0e}] instead of 1 (M(u)M(-u) = (5.20) denominator [{float(devM):.0e}]); '
       f'with prefactor 1/k it is 1 [{float(abs(r_k - 1)):.0e}]',
       abs(r_printed - kB**2/mp.cos(xiB)**2) < 1e-15 and abs(kB**2/mp.cos(xiB)**2 - 1) > 0.1 and abs(r_k - 1) < 1e-15
       and devM < 1e-15)

# ---------------------------------------------------------------- 11. breather singlet = image of the R-check residue
# prp-bsg-breather-k: R-check(theta_1 - theta_2): V(theta_1) (x) V(theta_2) -> V(theta_2) (x) V(theta_1) has a rank-1
# residue at the physical pole theta_12 = i u_n, whose image is the singlet in V(theta - i u_n/2) (x) V(theta + i u_n/2).
LAM3 = 3.3
q3 = q_phys(LAM3)


def inv_vec_q(ya, yb, qq):
    A, B = sl2_spin(1, qq, ya), sl2_spin(1, qq, yb)
    Mx = np.vstack([cop(A, B, kd, j) - (np.eye(4) if kd == 'k' else 0) for j in (0, 1) for kd in 'efk'])
    _, sv, Vh = np.linalg.svd(Mx)
    w = Vh[-1].conj()
    return sv[-1]/sv[0], w/w[1]


ok_res, ovs = True, []
for n in (1, 2):
    un = np.pi - n*np.pi/LAM3
    Rn = lambda z: Rsg(np.exp(LAM3*z), q3)
    dd = 1e-6
    res = (Rn(1j*un + dd)*dd + Rn(1j*un + 2*dd)*2*dd)/2
    U, sv, _ = np.linalg.svd(res)
    r_out, w_out = inv_vec_q(np.exp(LAM3*(0.37 - 1j*un/2)), np.exp(LAM3*(0.37 + 1j*un/2)), q3)
    ovs.append(abs(np.vdot(w_out, U[:, 0]))/np.linalg.norm(w_out))
    ok_res &= sv[1]/sv[0] < 1e-4 and r_out < 1e-10
report(f'the residue of R-check(theta_12) at theta_12 = i u_n is rank 1 and its image is the breather singlet in '
       f'V(theta - i u_n/2) (x) V(theta + i u_n/2) (n = 1, 2, lambda = 3.3) [overlap - 1: {max(abs(o - 1) for o in ovs):.0e}]',
       ok_res and max(abs(o - 1) for o in ovs) < 1e-9)

# ---------------------------------------------------------------- 12. BPTT a_n (eq. (3.5)) and the bootstrap for n = 1, 2
# a_n = a(u + nu_n) a(u - nu_n)/(a(u + nu_0) a(u - nu_0)) = prod_{l=1}^n {2(eta/pi - l)}, BPTT blocks
# (x) = sin(u/2 + pi x/2)/sin(u/2 - pi x/2) (eq-blocks with h = 1) and {y} = (y+1/2lam)(y-1/2lam)/((y+1/2lam - 1)(y-1/2lam + 1));
# a(u - nu_n) Q(eta, vt, u) b(u + nu_n) = Q(etabar, vt, u) a_n(eta, u), a = S_0, b = S_T.
L = mp.mpf('3.3')
eta, vt = mp.mpf('5.5'), mp.mpf('0.3')
etab = mp.pi*(L + 1) - eta
nu = lambda n: eta/L - (2*n + 1)*mp.pi/(2*L)
blk1 = lambda x, uu: mp.sin(uu/2 + x*mp.pi/2)/mp.sin(uu/2 - x*mp.pi/2)
curly = lambda y, uu: (blk1((y + 1)/(2*L), uu)*blk1((y - 1)/(2*L), uu)
                       / (blk1((y + 1)/(2*L) - 1, uu)*blk1((y - 1)/(2*L) + 1, uu)))
a_n = lambda n, uu: mp.fprod([curly(2*(eta/mp.pi - l), uu) for l in range(1, n + 1)]) if n else mp.mpf(1)
dev_ratio, dev_boot = 0, 0
for u in (mp.mpc('0.3', '0.2'), mp.mpc('-0.5', '0.1')):
    av = {(n, s): a_gz(1j*(u + s*nu(n))) for n in (0, 1, 2, 3) for s in (1, -1)}
    dev_ratio = max([dev_ratio] + [abs(av[(n, 1)]*av[(n, -1)]/(av[(0, 1)]*av[(0, -1)])/a_n(n, u) - 1) for n in (1, 2, 3)])
    Qrat = mp.exp(log_sigma(eta, u) - log_sigma(etab, u))*mp.cos(etab)/mp.cos(eta)      # Q(eta)/Q(etabar)
    bv = lambda n: mp.sin(L*(u + nu(n)))/mp.sin(L*(mp.pi - u - nu(n)))*av[(n, 1)]
    dev_boot = max([dev_boot] + [abs(av[(n, -1)]*bv(n)*Qrat/a_n(n, u) - 1) for n in (0, 1, 2)])
report(f'BPTT a_n: a(u+nu_n)a(u-nu_n)/(a(u+nu_0)a(u-nu_0)) = prod_(l=1)^n {{2(eta/pi - l)}} for n = 1, 2, 3 '
       f'[{float(dev_ratio):.0e}]; bootstrap a(u - nu_n) Q(eta) b(u + nu_n) = Q(etabar) a_n for n = 0, 1, 2 '
       f'[{float(dev_boot):.0e}]', dev_ratio < 1e-15 and dev_boot < 1e-15)

# ---------------------------------------------------------------- 13. breathers on the excited boundaries |0>, |1>
# Soliton K on |k> = sigma_x K(etabar) sigma_x a_k (eq-bsg-excited-k); the full matrix bootstrap for B_m gives
# K_{B_m}(etabar) b_mk with b_mk = a_k(u + pi/2 - m pi/2lam) a_k(u - pi/2 + m pi/2lam) (BPTT (3.9)), m = 1, 2, k = 0, 1.
X = np.array([[0, 1], [1, 0]])
devs = []
for n in (1, 2):
    un = np.pi - n*np.pi/LAM3
    th = 0.3
    ta, tb = th - 1j*un/2, th + 1j*un/2
    _, w = inv_vec_q(np.exp(LAM3*ta), np.exp(LAM3*tb), q3)
    S2 = complex(a_gz(ta + tb))*Rsg(np.exp(LAM3*(ta + tb)), q3)
    Ka, Kb = X @ K_full(ta, etab, vt) @ X, X @ K_full(tb, etab, vt) @ X
    uth, ua, ub = -1j*mp.mpf(th), -1j*mp.mpmathify(ta), -1j*mp.mpmathify(tb)
    for k in (0, 1):
        v = np.kron(np.eye(2), Ka*complex(a_n(k, ua))) @ S2 @ np.kron(np.eye(2), Kb*complex(a_n(k, ub))) @ w
        RB = v[1]
        bmk = a_n(k, uth + mp.pi/2 - n*mp.pi/(2*L))*a_n(k, uth - mp.pi/2 + n*mp.pi/(2*L))
        pred = complex(ghoshal_breather(n, uth, etab, vt)*bmk)
        devs.append(max(np.linalg.norm(v - RB*w)/np.linalg.norm(v), abs(RB/pred - 1)))
report(f'B_1, B_2 on the excited boundaries |0>, |1>: the full matrix bootstrap gives K_(B_m)(etabar) b_mk, '
       f'b_mk = a_k(u + pi/2 - m pi/2lam) a_k(u - pi/2 + m pi/2lam) (BPTT (3.9)) [{max(devs):.0e}]', max(devs) < 1e-10)

# exr-bsg-sign-q: Pi K(eps) Pi is proportional to K(-eps)
yy = 0.8 + 0.3j
report('Pi K(eps_0, eps_1) Pi is proportional to K(-eps_0, -eps_1), Pi = diag(1,-1) (exr-bsg-sign-q)',
       prop(Pi @ K_closed(yy, 0.43 - 0.2j, -0.71 + 0.35j) @ Pi, K_closed(yy, -0.43 + 0.2j, 0.71 - 0.35j)) < 1e-12)
L = mp.mpf(LAM)
