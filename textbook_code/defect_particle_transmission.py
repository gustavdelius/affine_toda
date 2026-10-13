"""Chapter 43 (sec-particle-transmission): breather and particle transmission amplitudes of the a_r^(1) defects.

Book ids supported:
  prp-sg-breather-transmission (re-check, a), prp-lightest-breather-universal (b), prp-an-particle-transmission and
  eq-dpt-particle (c, d), cnj-transmission-continuation and sec-transmission-limits (d, e), prp-dpt-vacuum (g),
  prp-dpt-doubled (f: the ch. 42 expectation for a_r, sec-jump-coideal), tbl-dpt-spectrum (a_r particle rows, c);
  exercises exr-dpt-telescope, exr-dpt-cz, exr-dpt-vacuum, exr-dpt-bs, exr-dpt-limits, exr-dpt-cdd, exr-dpt-doubled.

Conventions (briefs/part8_defects_conventions.md; ch. 41): Part VI coproduct Delta(e) = e (x) 1 + k (x) e,
Delta(f) = f (x) k^-1 + 1 (x) f; q = q_phys = -e^{-i pi omega}, omega = 4 pi/beta^2 - 1 (the coupling, Corrigan and
Zambon's gamma); y = e^{omega theta}; V = vector solitons (species 1, weights e_i), Vbar = species r (weights -e_i),
boundary_common.an_vector; T: V (x) W -> W (x) V, tensor order = spatial order; W = the a_r Weyl module eq-dt-weyl of
U_q(b^-) (type I_-), defect_common.an_weyl_module, mu = -e^{-i pi omega/2}/(q - q^-1), X = e^{omega(theta-theta_d-i pi/2)};
matrix parts eq-dt-tv, eq-dt-tvbar (an_T_vector), scalar factors g (eq-dt-g, an_scalar_g) and gbar = prod g
(eq-dt-bootstrap).  The root of unity of ch. 39 is written e^{2 pi i/h}; gamma_k = pi k/h.  Classical linear factors
(eq-an-linear-transmission): type I_-  T^cl_k(D) = (e^{-D} + e^{i gamma_k})/(e^{-D} + e^{-i gamma_k}), D = theta - theta_d;
type I_+  T^cl_{+,k}(D) = (e^{D} - e^{i gamma_k})/(e^{D} - e^{-i gamma_k}).

Checks (PASS lines):
 (a) sine-Gordon: the lowest breather, the singlet of V(theta - i u/2) (x) V(theta + i u/2), u = pi - pi/lambda,
     crosses eq-sg-transmission-matrix (with eq-sgd-scalar) with the factor T^cl of eq-shg-linear-transmission, which is
     T^cl_1 at h = 2; lambda = 1.37, 2.6, three defect labels.
 (b) prp-lightest-breather-universal.
     (b1) from the closed forms with their scalar factors (windows of an_weyl_module): the U_q(a_r^(1)) singlet of
          (species r)(theta - i ub) (x) (species 1)(theta + i ub), ub = pi/2 - pi/(h omega), crosses the defect as
          T^cl_1 times the identity on W, and (species 1) (x) (species r) as T^cl_r; with the opposite imaginary shifts
          there is no singlet; a_2, a_3, a_4.  cz2007 (7.5) is (T^cl_2, T^cl_1) in book labels.
     (b2) every species: B^(a)_1 = (species h-a)(theta - i ub) (x) (species a)(theta + i ub), each built from vector
          solitons at theta_c + i pi(2m - n + 1)/h (Lambda^n V, eq-dt-bootstrap), so that h vector solitons cross with
          scalar g each: the singlet is unique in the weight-zero subspace of V^{(x) h}; exact propagation through
          eq-dt-tv (no window); factor = T^cl_a for h = 3, 4, 5, all a, omega = 0.6, 1.37, 3.1, real and complex theta;
          independent of the defect label N and of the offset s (classes p = 0, 1, 2).
 (c) prp-an-particle-transmission: T^cl_a (both types) satisfy defect unitarity and crossing,
     T_a(theta) T_abar(theta + i pi) = 1, |T_a| = 1 for real theta and theta_d, Hermitian analyticity
     Tbar_a(theta) = conj T_a(-theta) with Tbar_a(theta) = T_abar(i pi - theta), and the bulk bootstrap
     T_c(theta) = T_a(theta + i ub_a) T_b(theta - i ub_b) at every fusing of a_r^(1), r <= 7, with the fusing angles
     located as poles of eq-an-smatrix; defect CDD factors: S_ab(theta - theta_0) (a standing particle) and a second
     defect satisfy the same equations; the factors have no pole in the physical strip for real theta_d.
 (d) classical limit: the linearized defect conditions (defect_common.defect_potential, dmatrix) about phi = psi = 0
     give R = 0 and T^cl_k (I_-) and T^cl_{+,k} (I_+) on the eigenvectors of pi (D = 1 + i cot gamma_k), r = 2, 3, 4;
     limits theta_d -> +-infinity: 1 and e^{2 pi i k/h} (type I_-), the reverse for I_+.
 (e) sinh-Gordon (Bajnok-Simon 0712.4292): their classical factor -(m sinh mu - i omega)/(m cosh mu - i k) and their
     (16) are T^cl with mu = -theta_d^cl, B_BS vartheta = -theta_d; their one-loop result (Sec. 4.2)
     T(theta, B vartheta(mu)) = T(theta, mu) + (mu b^2/8 pi)/(1 - i sinh(theta + mu)) + O(b^4) is the O(b^2) term of
     T^cl(theta - (1 - B/2) theta_d^cl) with B/2 = b^2/(8 pi + b^2) (cnj-defect-parameter).
 (f) the ch. 42 expectation for a_r, type I_-: for theta-independent unknown operators Y, kappa on W,
        e-type  T [pi(e_j) (x) 1 + pi(k_j) (x) Y] = [Y (x) 1 + kappa (x) pi(e_psi(j))] T,
        f-type  T [pi(f_j) (x) kappa + 1 (x) Y] = [Y (x) pi(k_psi(j))^-1 + 1 (x) pi(f_psi(j))] T
     (the block P T of eq-doubled-k with Delta_A of eq-dbd-coproduct), solved jointly for V and Vbar at several
     rapidities by least squares: e-type with psi(e_j) = e_{j-1} (the rotation pi) and f-type with psi = id hold;
     e-type with psi = id, pi^-1 and f-type with psi = pi, pi^-1 fail (a_2 all nodes, a_3 nodes 0, 1).
 (g) real-coupling vacuum: B_d^real >= 0 for real fields and sigma > 0 (both types), with equality only at
     phi = psi = 0; for sigma < 0 it is unbounded below (r = 2, 3).
Runtime: about one minute (62 s alone, 120 s under load; check f takes most of it).
"""
import time
import itertools
import numpy as np
import scipy.sparse as sps
from boundary_common import report, q_phys, an_vector, cop, prop_residual
from defect_common import (an_weyl_module, an_T_vector, an_scalar_g, an_scalar_gbar, an_scalar_constant,
                           sg_tmatrix, invariant_vector, an_data, dmatrix, defect_potential)

T0 = time.time()
TD = 0.25                       # defect rapidity theta_d used throughout


def Tcl(D, k, h, typ=1):
    """eq-an-linear-transmission: type I_- (typ = 1) and type I_+ (typ = -1), D = theta - theta_d."""
    g = np.pi*k/h
    if typ == 1:
        return (np.exp(-D) + np.exp(1j*g))/(np.exp(-D) + np.exp(-1j*g))
    return (np.exp(D) - np.exp(1j*g))/(np.exp(D) - np.exp(-1j*g))


def mu_of(q, om):
    return -np.exp(-1j*np.pi*om/2)/(q - 1/q)


# ---------------------------------------------------------------- (a) sine-Gordon lowest breather
def sg_breather(th, lam, n0):
    qq = q_phys(lam); u = np.pi - np.pi/lam
    t1, t2 = th - 1j*u/2, th + 1j*u/2
    rel, w = invariant_vector(np.exp(lam*t1), np.exp(lam*t2), qq)
    T1, V, W = sg_tmatrix(t1, lam, qq, TD, 6)
    T2 = sg_tmatrix(t2, lam, qq, TD, 6)[0]
    L = np.kron(T1, np.eye(2)) @ np.kron(np.eye(2), T2)
    e = np.zeros(W['dim']); e[list(W['n']).index(n0)] = 1
    res, c = prop_residual(L @ np.kron(w, e), np.kron(e, w))
    return c, max(rel, res)


devA, fitA = 0, 0
for lam in (1.37, 2.6):
    for th in (0.2, -0.5):
        for n0 in (-1, 0, 1):
            c, fit = sg_breather(th, lam, n0)
            D = th - TD
            shg = -1j*np.sinh(D/2 - 1j*np.pi/4)/np.sinh(D/2 + 1j*np.pi/4)
            devA = max(devA, abs(c - shg), abs(shg - Tcl(D, 1, 2))); fitA = max(fitA, fit)
report(f'(a) sine-Gordon: the singlet of V(theta - iu/2) (x) V(theta + iu/2) crosses eq-sg-transmission-matrix as a '
       f'scalar, independent of the defect label, equal to T^cl (eq-shg-linear-transmission) = T^cl_1 at h = 2, for '
       f'lambda = 1.37, 2.6 [{devA:.0e}; proportionality {fitA:.0e}]', devA < 1e-8 and fitA < 1e-10)


# ---------------------------------------------------------------- (b1) a_r: lightest breathers from the closed forms
def singlet(r, q, y1, c1, y2, c2):
    """Common null vector of Delta(e_j), Delta(f_j), Delta(k_j) - 1 on (species of c1)(y1) (x) (species of c2)(y2)."""
    A = an_vector(r, q, y1, conj=c1); B = an_vector(r, q, y2, conj=c2)
    h = r + 1
    M = np.vstack([cop(A, B, kd, j) - (np.eye(h*h) if kd == 'k' else 0) for j in range(h) for kd in 'efk'])
    _, s, Vh = np.linalg.svd(M, full_matrices=False)
    return s[-1]/s[0], Vh[-1].conj()


def breather_closed(r, om, th, cL, cR, sign=1, R2=14, s=None):
    """Two-particle bound state (species cL at th - i sign ub on the left, cR at th + i sign ub on the right) crossing
    the type-I_- defect with T_V = g M_V, T_Vbar = gbar M_Vbar.  Returns (singlet rel. sv, proportionality, factor)."""
    h = r + 1
    q = q_phys(om)
    W = an_weyl_module(r, -1, q, np.exp(om*TD), R2, s=s, mu=mu_of(q, om))
    dW = W['dim']
    ub = np.pi/2 - np.pi/(h*om)
    thL, thR = th - sign*1j*ub, th + sign*1j*ub
    rel, inv = singlet(r, q, np.exp(om*thL), cL, np.exp(om*thR), cR)
    if rel > 1e-10:
        return rel, None, None
    c = an_scalar_constant(r, om, q)

    def Tfull(t, conj):
        M = an_T_vector(W, np.exp(om*(t - TD - 1j*np.pi/2)), conj)
        return (an_scalar_gbar(t - TD, r, om, 400, c) if conj else an_scalar_g(t - TD, r, om, 400, c))*M
    TL, TR = Tfull(thL, cL), Tfull(thR, cR)
    k = np.where(W['deep'])[0][0]
    v = np.kron(inv, np.eye(dW)[k]).reshape(h, h*dW)
    v = (TR @ v.T).T.reshape(h*dW, h)            # [(a1, N), b2]
    v = TL @ v                                   # [(N, b1), b2]
    res, ratio = prop_residual(v.reshape(-1), np.kron(np.eye(dW)[k], inv))
    return rel, res, ratio


out_b1 = []; okb1 = True
for r, R2 in ((2, 14), (3, 14), (4, 8)):
    h = r + 1
    for th in (0.4, -0.3 + 0.2j):
        D = th - TD
        for (cL, cR), kexp in (((True, False), 1), ((False, True), r)):
            rel, res, val = breather_closed(r, 1.37, th, cL, cR, R2=R2, s=np.linspace(-.2, .2, h))
            relx, _, _ = breather_closed(r, 1.37, th, cL, cR, sign=-1, R2=R2)
            dev = abs(val - Tcl(D, kexp, h))
            out_b1.append((r, dev, res, relx))
            okb1 &= dev < 1e-7 and res < 1e-12 and relx > 1e-3
# cz2007 (7.5), their eta = theta_d, book labels: T_1^b1 = T^cl_2, T_2^b1 = T^cl_1 (h = 3)
Ds = np.array([0.3, -1.1, 0.4 + 0.3j])
cz75_1 = np.exp(-1j*np.pi/3)*np.sinh(Ds/2 - 1j*np.pi/6)/np.sinh(Ds/2 + 1j*np.pi/6)
cz75_2 = np.exp(1j*np.pi/3)*np.sinh(Ds/2 + 2j*np.pi/3)/np.sinh(Ds/2 + 1j*np.pi/3)
devcz = max(abs(cz75_1 - Tcl(Ds, 2, 3)).max(), abs(cz75_2 - Tcl(Ds, 1, 3)).max())
report(f'(b1) prp-lightest-breather-universal from eq-dt-tv, eq-dt-tvbar with g, gbar: the singlet of Vbar(theta - i ub) '
       f'(x) V(theta + i ub) crosses as T^cl_1 times the identity on W, that of V (x) Vbar as T^cl_r, a_2, a_3, a_4, '
       f'omega = 1.37, generic offset [max dev {max(o[1] for o in out_b1):.0e}, proportionality '
       f'{max(o[2] for o in out_b1):.0e}]; with the shifts reversed there is no singlet [smallest rel. sv '
       f'{min(o[3] for o in out_b1):.2f}]; cz2007 (7.5) = (T^cl_2, T^cl_1) in book labels [{devcz:.0e}]',
       okb1 and devcz < 1e-14)
print(f'# (a), (b1) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (b2) every species, exact propagation
def _node(h, j):
    return ((j - 1) % h, j % h) if j >= 1 else (h - 1, 0)


def _kdiag(h, q, j, i):
    a, b = _node(h, j)
    return q**((i == a) - (i == b))


def singlet_weight0(h, q, ys):
    """U_q(a_r^(1))-invariant vectors of V(y_1) (x) ... (x) V(y_h) in the weight-zero subspace (the h! permutations of
    v_1, ..., v_h), with Delta^(h)(e_j) = sum_p k_j^{(x) p-1} (x) e_j (x) 1 ..., Delta^(h)(f_j) = sum_p 1 ... (x) f_j (x)
    k_j^-1 ....  Returns (relative singular values, last right singular vector, permutations)."""
    perms = list(itertools.permutations(range(h)))
    col = {p: i for i, p in enumerate(perms)}
    rows, ent = {}, []
    for j in range(h):
        a, b = _node(h, j)
        for p in perms:
            for pos in range(h):
                if p[pos] == b:                                   # e_j v_b = y v_a
                    coef = np.prod([_kdiag(h, q, j, p[i]) for i in range(pos)])*ys[pos]
                    key = ('e', j, p[:pos] + (a,) + p[pos + 1:])
                elif p[pos] == a:                                 # f_j v_a = v_b / y
                    coef = np.prod([1/_kdiag(h, q, j, p[i]) for i in range(pos + 1, h)])/ys[pos]
                    key = ('f', j, p[:pos] + (b,) + p[pos + 1:])
                else:
                    continue
                rows.setdefault(key, len(rows))
                ent.append((rows[key], col[p], coef))
    M = np.zeros((len(rows), len(perms)), complex)
    for i_, j_, v in ent:
        M[i_, j_] += v
    try:
        _, s, Vh = np.linalg.svd(M, full_matrices=False)
    except np.linalg.LinAlgError:                      # LAPACK gesdd can fail to converge; gesvd is slower but robust
        from scipy.linalg import svd
        _, s, Vh = svd(M, full_matrices=False, lapack_driver='gesvd')
    return s/s[0], Vh[-1].conj(), perms


def breather_species(h, om, th, a, s=None, N0=None):
    """B^(a)_1 from h vector solitons crossing the type-I_- defect one after the other, the rightmost first, with
    M_V of eq-dt-tv propagated exactly on the labels N (no window) and the scalar g for each."""
    q = q_phys(om); r = h - 1
    s = np.zeros(h) if s is None else np.asarray(s, float)
    N0 = tuple([0]*h) if N0 is None else tuple(N0)
    ub = np.pi/2 - np.pi/(h*om)
    raps = [th - 1j*ub + 1j*np.pi*(2*m - (h - a) + 1)/h for m in range(h - a)] + \
           [th + 1j*ub + 1j*np.pi*(2*m - a + 1)/h for m in range(a)]
    rel, inv, perms = singlet_weight0(h, q, [np.exp(om*t) for t in raps])
    state = {(N0, p): c for p, c in zip(perms, inv) if abs(c) > 1e-14}
    I = np.eye(h, dtype=int)
    for pos in range(h - 1, -1, -1):
        X = np.exp(om*(raps[pos] - TD - 1j*np.pi/2))
        new = {}
        for (N, p), amp in state.items():
            i = p[pos]
            new[(N, p)] = new.get((N, p), 0) + amp*q**(-(N[i] + s[i]))
            b = (i - 1) % h
            key = (tuple(np.array(N) + I[i] - I[b]), p[:pos] + (b,) + p[pos + 1:])
            new[key] = new.get(key, 0) + amp*X
        state = new
    c = an_scalar_constant(r, om, q)
    G = np.prod([an_scalar_g(t - TD, r, om, 400, c) for t in raps])
    target = {(N0, p): cc for p, cc in zip(perms, inv)}
    keys = list(set(state) | set(target))
    res, ratio = prop_residual(np.array([state.get(k, 0) for k in keys]), np.array([target.get(k, 0) for k in keys]))
    return rel, res, ratio*G


devB, resB, nullB = 0, 0, True
for h in (3, 4, 5):
    for om in (0.6, 1.37, 3.1):
        for th in (0.4, -0.3 + 0.2j):
            for a in range(1, h):
                rel, res, val = breather_species(h, om, th, a, s=np.linspace(-0.3, 0.3, h))
                devB = max(devB, abs(val - Tcl(th - TD, a, h))); resB = max(resB, res)
                nullB &= rel[-1] < 1e-12 and rel[-2] > 1e-3
devS = 0
for s in (np.zeros(4), np.array([.75, -.25, -.25, -.25]), np.array([.5, .5, -.5, -.5])):
    for N0 in ((0, 0, 0, 0), (1, -1, 2, -2)):
        for a in (1, 2):
            _, res, val = breather_species(4, 1.37, 0.4, a, s=s, N0=N0)
            devS = max(devS, abs(val - Tcl(0.4 - TD, a, 4)), res)
report(f'(b2) prp-lightest-breather-universal, every species: B^(a)_1 from h vector solitons (singlet unique in the '
       f'weight-zero subspace) crosses the type-I_- defect as T^cl_a(theta - theta_d) times the identity, h = 3, 4, 5, '
       f'all a, omega = 0.6, 1.37, 3.1, real and complex theta [dev {devB:.0e}, proportionality {resB:.0e}]; the same for '
       f'the classes p = 0, 1, 2 (offset s = lambda_p) and other defect labels (a_3) [{devS:.0e}]',
       nullB and devB < 1e-7 and resB < 1e-12 and devS < 1e-7)
print(f'# (b2) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (c) the particle factors
def blk(x, h, th):
    return np.sinh(th/2 + 1j*np.pi*x/(2*h))/np.sinh(th/2 - 1j*np.pi*x/(2*h))


def brace(x, h, B, th):
    return blk(x - 1, h, th)*blk(x + 1, h, th)/(blk(x - 1 + B, h, th)*blk(x + 1 - B, h, th))


def S_an(a, b, h, B, th):
    """eq-an-smatrix."""
    out = 1
    for p in range(abs(a - b) + 1, a + b, 2):
        out = out*brace(p, h, B, th)
    return out


def fusings(h):
    """(a, b, c, ub_a, ub_b, u) for every fusing a b -> c of a_r^(1): c = a + b (a + b < h), c = a + b - h (a + b > h)."""
    out = []
    for a in range(1, h):
        for b in range(1, h):
            if a + b < h:
                out.append((a, b, a + b, np.pi*b/h, np.pi*a/h, np.pi*(a + b)/h))
            elif a + b > h:
                out.append((a, b, a + b - h, np.pi*(h - b)/h, np.pi*(h - a)/h, np.pi*(2*h - a - b)/h))
    return out


rng = np.random.default_rng(43)
ths = [0.31 + 0.2j, -0.77 + 1.1j, 1.3 - 0.4j]
reals = [0.3, -1.2, 2.1]
okc = True; devU = devH = devF = devP = devCDD = 0; poles_ok = True
for h in range(2, 9):
    for typ in (1, -1):
        for k in range(1, h):
            kb = h - k
            for th in ths:
                devU = max(devU, abs(Tcl(th - TD, k, h, typ)*Tcl(th + 1j*np.pi - TD, kb, h, typ) - 1))
            for th in reals:
                devU = max(devU, abs(abs(Tcl(th - TD, k, h, typ)) - 1))
                # Hermitian analyticity: Tbar_k(theta) = T_kbar(i pi - theta) (as functions of theta + theta_d)
                tbar = Tcl(1j*np.pi - th - TD, kb, h, typ)
                devH = max(devH, abs(tbar - np.conj(Tcl(-th - TD, k, h, typ))))
            # no pole in the physical strip 0 < Im theta < pi/2 for real theta_d: pole of T^cl at theta_d -+ i(pi - g)
            g = np.pi*k/h
            pole = TD - 1j*(np.pi - g) if typ == 1 else TD - 1j*g
            poles_ok &= not (0 < pole.imag < np.pi/2) and abs(1/Tcl(pole - TD + 1e-9, k, h, typ)) < 1e-6
    for (a, b, c, ua, ub_, u) in fusings(h):
        # the fusing angle is a pole of eq-an-smatrix, and the masses close
        eps = 1e-6
        Sp = S_an(a, b, h, 0.37, 1j*u + eps)
        mass = 2*np.sin(np.pi*a/h)*np.cos(ua) + 2*np.sin(np.pi*b/h)*np.cos(ub_) - 2*np.sin(np.pi*c/h)
        poles_ok &= abs(Sp) > 1e3 and abs(mass) < 1e-12
        for typ in (1, -1):
            for th in ths:
                lhs = Tcl(th + 1j*ua - TD, a, h, typ)*Tcl(th - 1j*ub_ - TD, b, h, typ)
                devF = max(devF, abs(lhs - Tcl(th - TD, c, h, typ)))
                # defect CDD factors: a standing particle e (S_ae(theta - theta0)) and a second defect at theta_d'
                for e in (1, h - 1):
                    t0 = 0.4 - 0.2j
                    cdd = lambda t, x: S_an(x, e, h, 0.37, t - t0)
                    devCDD = max(devCDD, abs(cdd(th + 1j*ua, a)*cdd(th - 1j*ub_, b) - cdd(th, c)),
                                 abs(cdd(th, a)*cdd(th + 1j*np.pi, h - a) - 1))
report(f'(c) prp-an-particle-transmission: T^cl_a of both types, h = 2..8: T_a(theta) T_abar(theta + i pi) = 1 and '
       f'|T_a| = 1 for real theta [{devU:.0e}], Hermitian analyticity Tbar_a(theta) = conj T_a(-theta) [{devH:.0e}], the '
       f'bulk bootstrap at every fusing a b -> a+b, a+b-h [{devF:.0e}] (fusing angles are poles of eq-an-smatrix, '
       f'masses close); defect CDD factors S_ae(theta - theta_0) satisfy the same equations [{devCDD:.0e}]; no pole in '
       f'the physical strip for real theta_d', devU < 1e-12 and devH < 1e-12 and devF < 1e-12 and devCDD < 1e-10
       and poles_ok)
print(f'# (c) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (d) classical limit: linearized defect conditions
def linear_RT(r, typ, k, th, thd, m=1.0):
    """Plane wave of mode k (eigenvector of pi with eigenvalue e^{2 pi i k/h}) incident from the left on the
    linearized defect conditions eq-defect-conditions; returns (R, T vectors, rho, D-eigenvalue)."""
    d = an_data(r); h = r + 1
    Dm, E = dmatrix(r, typ)
    sig = np.exp(-thd)
    Hs = {}
    for which in ('phi', 'psi'):
        cols_p, cols_s = [], []
        for i in range(r):
            ei = np.zeros(r); ei[i] = 1e-20
            z = np.zeros(r)
            f, s_ = (1j*ei, z) if which == 'phi' else (z, 1j*ei)
            _, gphi, gpsi, _ = defect_potential(f, s_, sig, r, 1.0, m, 'real', typ)
            cols_p.append(np.imag(gphi)/1e-20); cols_s.append(np.imag(gpsi)/1e-20)
        Hs[which] = (np.array(cols_p).T, np.array(cols_s).T)   # d(grad_phi)/d(which), d(grad_psi)/d(which)
    Hpp, Hsp = Hs['phi']; Hps, Hss = Hs['psi']
    ev, U = np.linalg.eig(d['pi'])
    i = np.argmin(abs(ev - np.exp(2j*np.pi*k/h)))
    rho = U[:, i]
    mk = 2*m*np.sin(np.pi*k/h)
    nu, p = mk*np.cosh(th), mk*np.sinh(th)
    I = np.eye(r)
    A = np.block([[-1j*p*I + 1j*nu*E + Hpp, 1j*nu*Dm + Hps],
                  [1j*nu*Dm.T - Hsp, 1j*p*I - 1j*nu*E - Hss]])
    rhs = -np.concatenate([(1j*p*I + 1j*nu*E + Hpp) @ rho, (1j*nu*Dm.T - Hsp) @ rho])
    sol = np.linalg.solve(A, rhs)
    Dev = (rho.conj() @ dmatrix(r, 1)[0] @ rho)/(rho.conj() @ rho)
    return sol[:r], sol[r:], rho, Dev


devD = 0; okD = True
for r in (2, 3, 4):
    h = r + 1
    for typ in (1, -1):
        for k in range(1, h):
            for th in (0.3, -0.8, 1.4):
                R, Tv, rho, Dev = linear_RT(r, typ, k, th, TD)
                okD &= abs(Dev - (1 + 1j/np.tan(np.pi*k/h))) < 1e-12
                devD = max(devD, np.abs(R).max(), np.abs(Tv - Tcl(th - TD, k, h, typ)*rho).max())
lim = max(abs(Tcl(0.3 - 40, k, h) - 1) + abs(Tcl(0.3 + 40, k, h) - np.exp(2j*np.pi*k/h))
          + abs(Tcl(0.3 - 40, k, h, -1) - np.exp(2j*np.pi*k/h)) + abs(Tcl(0.3 + 40, k, h, -1) - 1)
          for h in (3, 4, 5) for k in range(1, h))
report(f'(d) classical limit: the linearized defect conditions (defect_common, D = 1 + i cot gamma_k on the mode of '
       f'pi-eigenvalue e^{{2 pi i k/h}}) give R = 0 and T = T^cl_k (type I_-), T^cl_+,k (type I_+), r = 2, 3, 4 '
       f'[{devD:.0e}], which by (b) is the quantum B^(k)_1 factor at every coupling; theta_d -> +inf, -inf give '
       f'1, e^{{2 pi i k/h}} (I_-) and the reverse (I_+) [{lim:.0e}]', okD and devD < 1e-10 and lim < 1e-12)


# ---------------------------------------------------------------- (e) sinh-Gordon: Bajnok-Simon
def T_bs(th, mu):
    return -1j*np.sinh(th/2 - 1j*np.pi/4 + mu/2)/np.sinh(th/2 + 1j*np.pi/4 + mu/2)


devE1 = devE2 = 0
for th in (0.3, -0.9, 0.5 + 0.2j):
    for thcl in (0.7, -0.4):
        mu = -thcl
        cl = -(np.sinh(mu) - 1j*np.cosh(th))/(np.cosh(mu) - 1j*np.sinh(th))
        devE1 = max(devE1, abs(cl - T_bs(th, mu)), abs(T_bs(th, mu) - Tcl(th - thcl, 1, 2)))
        b2 = 1e-4
        Bhalf = b2/(8*np.pi + b2)
        exact = Tcl(th - (1 - Bhalf)*thcl, 1, 2)
        oneloop = T_bs(th, mu) + (mu*b2/(8*np.pi))/(1 - 1j*np.sinh(th + mu))
        devE2 = max(devE2, abs(exact - oneloop)/b2)
report(f'(e) sinh-Gordon: Bajnok-Simon\'s classical factor and their (16) equal T^cl with mu = -theta_d^cl '
       f'[{devE1:.0e}]; their one-loop result (Sec. 4.2) is the O(b^2) term of T^cl(theta - (1 - B/2) theta_d^cl), '
       f'B/2 = b^2/(8 pi + b^2) [remainder/b^2 = {devE2:.0e} at b^2 = 1e-4]', devE1 < 1e-13 and devE2 < 1e-3)
print(f'# (d), (e) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (f) the ch. 42 expectation for a_r
def borel_twist_residual(r, kind, psi, R2, thetas, s, nodes):
    """Least-squares residual of the e-type or f-type equation (docstring, check f), jointly for V and Vbar and all
    rapidities, over theta-independent operators Y and kappa of the required weights on W."""
    h = r + 1
    om = 1.37; q = q_phys(om)
    W = an_weyl_module(r, -1, q, np.exp(om*TD), R2, s=s, mu=mu_of(q, om))
    dW = W['dim']; labs = W['labels']; look = W['look']
    I = np.eye(h)
    alpha = [I[(j - 1) % h] - I[j % h] for j in range(h)]
    inter = W['interior']
    rsel = np.where(np.repeat(inter, h))[0]; csel = np.where(np.tile(inter, h))[0]
    Ih = sps.identity(h, format='csr'); IW = sps.identity(dW, format='csr')
    worst = 0
    for j in nodes:
        jp = {'id': j, 'pi': (j - 1) % h, 'pi^-1': (j + 1) % h}[psi]
        wY = alpha[j] if kind == 'e' else -alpha[jp]
        wK = alpha[j] - alpha[jp]
        As, bs = [], []
        for conj in (False, True):
            for th in thetas:
                V = an_vector(r, q, np.exp(om*th), conj=conj)
                T = sps.csr_matrix(an_T_vector(W, np.exp(om*(th - TD - 1j*np.pi/2)), conj))
                Vc = {key: sps.csr_matrix(V[key]) for key in V if isinstance(key, tuple)}
                if kind == 'e':
                    const = T @ sps.kron(Vc[('e', j)], IW)
                    termY = lambda E: T @ sps.kron(Vc[('k', j)], E) - sps.kron(E, Ih) @ T
                    termK = lambda E: -sps.kron(E, Vc[('e', jp)]) @ T
                else:
                    const = -sps.kron(IW, Vc[('f', jp)]) @ T
                    kinv = sps.csr_matrix(np.linalg.inv(V[('k', jp)]))
                    termY = lambda E: T @ sps.kron(Ih, E) - sps.kron(E, kinv) @ T
                    termK = lambda E: T @ sps.kron(Vc[('f', j)], E)
                cm = const.tocsr()[rsel][:, csel].toarray().ravel()
                cols = []
                for w, fn in ((wY, termY), (wK, termK)):
                    for kk in range(dW):
                        t = tuple(map(int, np.round(labs[kk] + w)))
                        if t not in look:
                            continue
                        E = sps.csr_matrix(([1.0], ([look[t]], [kk])), shape=(dW, dW))
                        cols.append(fn(E).tocsr()[rsel][:, csel].toarray().ravel())
                A = np.array(cols).T
                keep = (np.abs(A).sum(1) + np.abs(cm)) > 0
                As.append(A[keep]); bs.append(-cm[keep])
        A = np.vstack(As); b = np.concatenate(bs)
        sol = np.linalg.lstsq(A, b, rcond=None)[0]
        worst = max(worst, np.linalg.norm(A @ sol - b)/np.linalg.norm(b))
    return worst


resF = {}
for r, R2, nodes, thetas in ((2, 18, range(3), (0.3, -0.4, 0.9)), (3, 14, (0, 1), (0.3, -0.4))):
    s = np.linspace(-0.2, 0.2, r + 1)
    for kind, psi in (('e', 'pi'), ('f', 'id'), ('e', 'id'), ('e', 'pi^-1')) + ((('f', 'pi'), ('f', 'pi^-1')) if r == 2 else ()):
        resF[(r, kind, psi)] = borel_twist_residual(r, kind, psi, R2, thetas, s, nodes)
hold = [v for (r, kd, ps), v in resF.items() if (kd, ps) in (('e', 'pi'), ('f', 'id'))]
fail = [v for (r, kd, ps), v in resF.items() if (kd, ps) not in (('e', 'pi'), ('f', 'id'))]
report(f'(f) ch. 42 expectation for a_r, type I_-: the doubled boundary keeps Delta(U_q(b^-)) (f-type, psi = id) and '
       f'the pi-twisted (1 (x) pi) Delta(U_q(b^+)) (e-type, psi(e_j) = e_(j-1)) with theta-independent Y, kappa on W '
       f'[{max(hold):.0e}]; untwisted and pi^-1-twisted e, and pi^+-1-twisted f, fail [min {min(fail):.2f}] (a_2 all '
       f'nodes, 3 rapidities; a_3 nodes 0, 1, 2 rapidities; V and Vbar jointly)', max(hold) < 1e-12 and min(fail) > 0.3)
print(f'# (f) done after {time.time() - T0:.0f} s')


# ---------------------------------------------------------------- (g) real-coupling vacuum
okG = True; minpos = np.inf; minneg = 0
for r in (2, 3):
    for typ in (1, -1):
        for sig in (0.4, 2.3):
            for _ in range(400):
                phi, psi = rng.normal(size=r)*rng.uniform(0.01, 3), rng.normal(size=r)*rng.uniform(0.01, 3)
                Bd = defect_potential(phi, psi, sig, r, 1.0, 1.0, 'real', typ)[0].real
                minpos = min(minpos, Bd/(np.dot(phi, phi) + np.dot(psi, psi)))
            okG &= abs(defect_potential(np.zeros(r), np.zeros(r), sig, r, 1.0, 1.0, 'real', typ)[0]) < 1e-15
            # sigma < 0: along a direction with one f_j growing the potential decreases without bound
            vals = []
            for _ in range(50):
                v = rng.normal(size=2*r)
                vals.append(min(defect_potential(t*v[:r], t*v[r:], -sig, r, 1.0, 1.0, 'real', typ)[0].real
                                for t in (5.0, 10.0)))
            minneg = min(minneg, min(vals))
report(f'(g) real-coupling vacuum: B_d^real(phi, psi) >= 0 for real fields and sigma > 0, both types, r = 2, 3, with '
       f'B_d/|fields|^2 >= {minpos:.2e} on random samples and B_d(0, 0) = 0; for sigma < 0 it is unbounded below '
       f'(min sample {minneg:.1e})', okG and minpos > 0 and minneg < -1e3)
print(f'# total runtime {time.time() - T0:.0f} s')
