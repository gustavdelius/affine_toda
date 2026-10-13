"""Chapter 43 (sec-particle-transmission): breather and particle transmission amplitudes of the a_r^(1) defects.

Book ids supported:
  prp-sg-breather-transmission (re-check, a), prp-lightest-breather-universal (b), prp-an-particle-transmission and
  eq-dpt-particle (c, d), cnj-transmission-continuation and sec-transmission-limits (d, e), prp-dpt-vacuum (g),
  prp-dpt-doubled and eq-dpt-second-module (f, h: the ch. 42 expectation for a_r, sec-jump-coideal), eq-dpt-soliton-cdd
  (i), tbl-dpt-spectrum and sec-defect-spectrum (c, j); exercises exr-dpt-telescope, exr-dpt-cz, exr-dpt-vacuum,
  exr-dpt-bs, exr-dpt-limits, exr-dpt-imaginary, exr-dpt-cdd, exr-dpt-doubled.

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
          independent of the defect label N and of the offset s (classes p = 0, 1, 2).  The lowest breather is a physical
          bound state only for h omega > 2; h = 3, omega = 0.6 tests the algebraic identity only.
     (b3) optional, off by default (environment variable DPT_A5=1, two to three minutes): (b2) for a_5 (h = 6), all five
          species, omega = 1.37, one complex rapidity.
 (c) prp-an-particle-transmission: T^cl_a (both types) satisfy defect unitarity and crossing,
     T_a(theta) T_abar(theta + i pi) = 1, |T_a| = 1 for real theta and theta_d, Hermitian analyticity
     Tbar_a(theta) = conj T_a(-theta) with Tbar_a(theta) = T_abar(i pi - theta), and the bulk bootstrap
     T_c(theta) = T_a(theta + i ub_a) T_b(theta - i ub_b) at every fusing of a_r^(1): exactly with sympy, for symbolic
     u = e^{i gamma_a}, v = e^{i gamma_b} (both fusings, both types; every h), and numerically for r <= 7, with the
     fusing angles located as poles of eq-an-smatrix; defect CDD factors: S_ae(theta - theta_0) (a particle of species e
     at rapidity theta_0, standing for theta_0 = 0) satisfies the same equations; the factors have no pole in the
     physical strip for real theta_d.
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
     phi = psi = 0; for sigma < 0 it is unbounded below (r = 2, 3), along random directions and along phi = psi = t v,
     where it equals (m sigma/beta^2)(sum_j e^{t beta alpha_j.v} - h).
 (h) the second module eq-dpt-second-module (proof of prp-dpt-doubled), on exactly propagated labels (no window):
     rho~(e_j)|N> = e^{omega(theta_d + i pi/2)}/(q - q^-1) q^{-(N+s)_j}|N + alpha_j>, rho~(k_j)|N> = |N + alpha_j - alpha_{j-1}>
     solves the e-type equation (psi = pi) and the Cartan equation T(pi(k_j) (x) kappa) = (kappa (x) pi(k_{j-1})) T for V
     and Vbar, h = 3, 4, 5, every node, three rapidities; adding c_j q^{(N+s)_j} to rho~(e_j) still solves them but breaks
     k_j e_j k_j^-1 = q^2 e_j; the U_q(b^+) relations (Cartan conjugation, q-commutation, q-Serre, prod k_j = 1);
     the sparse inverse of eq-dt-tvbar used in the all-rank proof of the Vbar half, H M_Vbar = 1 - X^h/q with
     H|i,N> = q^{-(N+s)_i}|i,N> - (X/q)|i+1, N - e_i + e_{i+1}>; for h = 2, in the basis tau^-n q^{n(n-1)/2 + ns/2}|n>, tau = -i e^{i pi omega/2}, q^{s/2} M_V is
     eq-sg-transmission-matrix and rho~ is eq-dbd-second-module; the y-power obstructions of the other twists.
 (i) the soliton CDD condition eq-dpt-soliton-cdd: with g -> g f and gbar from eq-dt-bootstrap, the two crossing
     relations of ch. 41 are prod_{k<h} f(D - i pi + 2 pi i k/h) = 1 (at D and at a shifted D); the simpler condition
     f(D) f(D + i pi) = 1 is not necessary for h > 2 (T^cl_1 satisfies the string condition for h = 2..8 but not the
     simpler one), it is sufficient for even h (the string pairs k with k + h/2, arguments i pi apart:
     (sinh D + 0.3i)/(sinh D - 0.3i), which satisfies the simpler condition, passes the string condition at h = 4, 6),
     and it is not sufficient for odd h (the same factor fails at h = 3, 5); the multiplier of B^(2)_1 of a_3 for
     f = T^cl_1 differs at omega = sqrt2 and sqrt5.
 (j) tbl-dpt-spectrum: on every charge sector the residue of T_V at the species-1 pole has rank h - 1 and that of
     T_Vbar at the species-r pole (X^h = q) rank 1, h = 3..6, random offsets; the pole of gbar there is simple and not
     cancelled; the Lee-Yang bound-state poles are in Tbar_b (eq-ly-transmission), not in T_b; for theta_d = i v the
     particle factor is a phase for real theta iff v in pi Z, and at v = pi (inside the window for a < h/2) it is a
     phase with its pole at i gamma_a.
Runtime: about one minute (check f takes most of it).
"""
import os
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
if os.environ.get('DPT_A5'):
    devA5, resA5, nullA5 = 0, 0, True
    for a in range(1, 6):
        rel, res, val = breather_species(6, 1.37, -0.3 + 0.2j, a, s=np.linspace(-0.3, 0.3, 6))
        devA5 = max(devA5, abs(val - Tcl(-0.3 + 0.2j - TD, a, 6))); resA5 = max(resA5, res)
        nullA5 &= rel[-1] < 1e-12 and rel[-2] > 1e-3
    report(f'(b3) optional (DPT_A5=1): a_5, all five species, omega = 1.37, theta = -0.3 + 0.2i: B^(a)_1 crosses as '
           f'T^cl_a [dev {devA5:.0e}, proportionality {resA5:.0e}], singlet unique', nullA5 and devA5 < 1e-7 and resA5 < 1e-12)
    print(f'# (b3) done after {time.time() - T0:.0f} s')


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
# exact identities (proof of prp-an-particle-transmission): x = e^{-(theta - theta_d)} (I_-), y = e^{theta - theta_d} (I_+),
# u = e^{i gamma_a}, v = e^{i gamma_b}; theta -> theta + i c multiplies x by e^{-ic} and y by e^{ic}
import sympy as sp
xs, ys, ws = sp.symbols('x y w', positive=True)
us, vs = sp.symbols('u v', nonzero=True)
Tm = lambda X, U: (X + U)/(X + 1/U)
Tp = lambda Y, U: (Y - U)/(Y - 1/U)
exact = [Tm(xs/vs, us)*Tm(xs*us, vs) - Tm(xs, us*vs),          # I_-, a + b < h, e^{i gamma_c} = u v
         Tm(-xs*vs, us)*Tm(-xs/us, vs) - Tm(xs, -us*vs),       # I_-, a + b > h, e^{i gamma_c} = -u v
         Tp(ys*vs, us)*Tp(ys/us, vs) - Tp(ys, us*vs),          # I_+, a + b < h
         Tp(-ys/vs, us)*Tp(-ys*us, vs) - Tp(ys, -us*vs),       # I_+, a + b > h
         Tm(xs, us)*Tm(-xs, -1/us) - 1, Tp(ys, us)*Tp(-ys, -1/us) - 1,     # unitarity/crossing, e^{i gamma_abar} = -1/u
         Tm(-ws, -1/us) - Tm(ws, 1/us), Tp(-1/ws, -1/us) - Tp(1/ws, 1/us)]  # Hermitian analyticity, w = e^{theta+theta_d}
exact_ok = all(sp.simplify(e) == 0 for e in exact)
report(f'(c) prp-an-particle-transmission: the bootstrap at both fusings, unitarity/crossing and Hermitian analyticity '
       f'are exact rational identities in u = e^(i gamma_a), v = e^(i gamma_b), both types (sympy: {exact_ok}); '
       f'numerically, T^cl_a of both types, h = 2..8: T_a(theta) T_abar(theta + i pi) = 1 and '
       f'|T_a| = 1 for real theta [{devU:.0e}], Hermitian analyticity Tbar_a(theta) = conj T_a(-theta) [{devH:.0e}], the '
       f'bulk bootstrap at every fusing a b -> a+b, a+b-h [{devF:.0e}] (fusing angles are poles of eq-an-smatrix, '
       f'masses close); defect CDD factors S_ae(theta - theta_0) satisfy the same equations [{devCDD:.0e}]; no pole in '
       f'the physical strip for real theta_d', exact_ok and devU < 1e-12 and devH < 1e-12 and devF < 1e-12
       and devCDD < 1e-10 and poles_ok)
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
            # the explicit direction of the proof: phi = psi = t v gives (m sigma'/beta^2)(sum_j e^{t beta alpha_j.v} - h)
            al = an_data(r)['alpha']; v = rng.normal(size=r); v = v/max(al @ v)     # max_j alpha_j.v = 1
            for t in (3.0, 12.0):
                Bt = defect_potential(t*v, t*v, -sig, r, 1.0, 1.0, 'real', typ)[0].real
                sp_ = -sig if typ == 1 else -1/sig
                okG &= abs(Bt - sp_*(np.exp(t*(al @ v)).sum() - (r + 1))) < 1e-9*abs(Bt)
            okG &= Bt < -1e4
report(f'(g) real-coupling vacuum: B_d^real(phi, psi) >= 0 for real fields and sigma > 0, both types, r = 2, 3, with '
       f'B_d/|fields|^2 >= {minpos:.2e} on random samples and B_d(0, 0) = 0; for sigma < 0 it is unbounded below '
       f'(min sample {minneg:.1e}; along phi = psi = t v it is (m sigma/beta^2)(sum_j e^(t beta alpha_j.v) - h))',
       okG and minpos > 0 and minneg < -1e3)


# ---------------------------------------------------------------- (h) the second module eq-dpt-second-module
# Exact propagation on dicts of labels (0-based: alpha_j = e_{j-1} - e_j, first end a = j-1, second end b = j).
OMH = np.sqrt(2); QH = q_phys(OMH); TH = np.exp(OMH*(TD + 1j*np.pi/2))


def _u(h, i):
    v = np.zeros(h, dtype=int); v[i % h] = 1; return v


def _alpha(h, j):
    return _u(h, j - 1) - _u(h, j)


def _acc(d, k, v):
    d[k] = d.get(k, 0) + v


def _e_rep(h, conj, y, j, i):
    """pi(e_j) on v_i: V: e_j v_j = y v_{j-1}; Vbar: e_j vb_{j-1} = y vb_j (an_vector)."""
    a, b = (j - 1) % h, j % h
    if not conj:
        return [(a, y)] if i == b else []
    return [(b, y)] if i == a else []


def _k_rep(h, conj, j, i):
    a, b = (j - 1) % h, j % h
    H = (i == a) - (i == b)
    return QH**(-H if conj else H)


def _T_dict(h, s, X, conj, state):
    """eq-dt-tv / eq-dt-tvbar on a dict {(i, N): amp} of V (x) W; returns {(N, i): amp} on W (x) V."""
    out = {}
    for (i, N), c in state.items():
        Nv = np.array(N); Ns = Nv + s
        if not conj:
            _acc(out, (N, i), c*QH**(-Ns[i]))
            _acc(out, (tuple(Nv + _u(h, i) - _u(h, i - 1)), (i - 1) % h), c*X)
        else:
            kk = np.zeros(h)
            for jj in range(h):
                b = (i + jj) % h
                kk[b] += 1
                _acc(out, (tuple(Nv - _u(h, i) + _u(h, b)), b), c*X**jj*QH**(Ns @ kk))
    return out


def _Y(h, s, j, c=0.0):
    """rho~(e_j) of eq-dpt-second-module, plus c q^{(N+s)_j} (the free constant of the y^1 equations)."""
    a = (j - 1) % h
    return lambda N: [(tuple(np.array(N) + _alpha(h, j)), TH/(QH - 1/QH)*QH**(-(N[a] + s[a])) + c*QH**(N[a] + s[a]))]


def _K(h, j, sign=1):
    w = _alpha(h, j) - _alpha(h, j - 1)
    return lambda N: [(tuple(np.array(N) + sign*w), 1.0)]


def _rel(d1, d2):
    keys = set(d1) | set(d2)
    num = max([abs(d1.get(k, 0) - d2.get(k, 0)) for k in keys] + [0])
    return num/max([abs(v) for v in list(d1.values()) + list(d2.values())] + [1e-300])


def e_type_residual(h, s, th, conj, j, jpsi, c, labels):
    """max relative residual of T[pi(e_j) (x) 1 + pi(k_j) (x) Y] = [Y (x) 1 + kappa (x) pi(e_psi(j))] T and of the
    Cartan equation T(pi(k_j) (x) kappa) = (kappa (x) pi(k_psi(j))) T on v_i (x) |N>."""
    y = np.exp(OMH*th); X = np.exp(OMH*(th - TD - 1j*np.pi/2))
    Y, K = _Y(h, s, j, c), _K(h, j)
    worst = 0
    for N in labels:
        for i in range(h):
            L = {}
            for (i2, v) in _e_rep(h, conj, y, j, i):
                _acc(L, (i2, N), v)
            for (N2, v) in Y(N):
                _acc(L, (i, N2), _k_rep(h, conj, j, i)*v)
            L = _T_dict(h, s, X, conj, L)
            R = {}
            for (N2, i2), c0 in _T_dict(h, s, X, conj, {(i, N): 1.0}).items():
                for (N3, v) in Y(N2):
                    _acc(R, (N3, i2), c0*v)
                for (i3, v) in _e_rep(h, conj, y, jpsi, i2):
                    _acc(R, (K(N2)[0][0], i3), c0*v)
            worst = max(worst, _rel(L, R))
            Lk = _T_dict(h, s, X, conj, {(i, K(N)[0][0]): _k_rep(h, conj, j, i)})
            Rk = {}
            for (N2, i2), c0 in _T_dict(h, s, X, conj, {(i, N): 1.0}).items():
                _acc(Rk, (K(N2)[0][0], i2), c0*_k_rep(h, conj, jpsi, i2))
            worst = max(worst, _rel(Lk, Rk))
    return worst


def _compose(*ops):
    def f(N):
        cur = [(N, 1.0)]
        for op in reversed(ops):
            cur = [(M2, c*v) for (M, c) in cur for (M2, v) in op(M)]
        return cur
    return f


def _opdiff(A, B, labels, scale=None):
    w = 0
    for N in labels:
        a, b = {}, {}
        for k, v in A(N):
            _acc(a, k, v)
        for k, v in B(N):
            _acc(b, k, v)
        if scale is None:
            w = max(w, _rel(a, b))
        else:
            w = max(w, max([abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b)] + [0])/scale(N))
    return w


def _H_dict(h, s, X, state):
    """The sparse inverse of eq-dt-tvbar up to its scalar (proof of prp-dpt-doubled): H = A^-1 (1 - X U),
    H|i, N> = q^{-(N+s)_i}|i, N> - (X/q)|i+1, N - e_i + e_{i+1}>, on {(N, i): amp} of W (x) Vbar; returns {(i, N): amp}."""
    out = {}
    for (N, i), c in state.items():
        Nv = np.array(N)
        _acc(out, (i, N), c*QH**(-(Nv + s)[i]))
        _acc(out, ((i + 1) % h, tuple(Nv - _u(h, i) + _u(h, i + 1))), -c*X/QH)
    return out


rngh = np.random.default_rng(431)
thsH = (-0.32, 0.41, 0.82 + 0.11j)
resV = resVb = resC = cartc = rels = devH = 0; obstr_ok = True
for h in (3, 4, 5):
    s = rngh.normal(size=h)*0.3; s -= s.mean()
    labs = [tuple(np.zeros(h, dtype=int)), tuple(_u(h, 0) - _u(h, 2)), tuple(2*_u(h, 1) - 2*_u(h, h - 1))]
    for th in thsH:                                    # H M_Vbar = 1 - X^h/q on every input
        X = np.exp(OMH*(th - TD - 1j*np.pi/2))
        for N in labs:
            for i in range(h):
                devH = max(devH, _rel(_H_dict(h, s, X, _T_dict(h, s, X, True, {(i, N): 1.0})), {(i, N): 1 - X**h/QH}))
    for j in range(h):
        jp = (j - 1) % h
        for th in thsH:
            resV = max(resV, e_type_residual(h, s, th, False, j, jp, 0.0, labs))
            resVb = max(resVb, e_type_residual(h, s, th, True, j, jp, 0.0, labs))
            resC = max(resC, e_type_residual(h, s, th, False, j, jp, 0.3, labs),
                       e_type_residual(h, s, th, True, j, jp, 0.3, labs))
        # the free constant breaks k_j e_j k_j^-1 = q^2 e_j
        cartc = max(cartc, _opdiff(_compose(_K(h, j), _Y(h, s, j, 0.3), _K(h, j, -1)),
                                   lambda N, j=j: [(M, QH**2*v) for (M, v) in _Y(h, s, j, 0.3)(N)], labs))
        # U_q(b^+) relations of rho~
        for i in range(h):
            aij = 2 if i == j else (-1 if (i - j) % h in (1, h - 1) else 0)
            rels = max(rels, _opdiff(_compose(_K(h, i), _Y(h, s, j), _K(h, i, -1)),
                                     lambda N, j=j, aij=aij: [(M, QH**aij*v) for (M, v) in _Y(h, s, j)(N)], labs))
            if (i - j) % h == 1:                       # adjacent i = j + 1: e_j e_{j+1} = q^-1 e_{j+1} e_j, q-Serre
                Yj, Yi = _Y(h, s, j), _Y(h, s, i)
                rels = max(rels, _opdiff(_compose(Yj, Yi), lambda N: [(M, v/QH) for (M, v) in _compose(Yi, Yj)(N)], labs))
                for (A, B) in ((Yj, Yi), (Yi, Yj)):
                    serre = lambda N, A=A, B=B: (_compose(A, A, B)(N) + [(M, -(QH + 1/QH)*v) for (M, v) in _compose(A, B, A)(N)]
                                                 + _compose(B, A, A)(N))
                    rels = max(rels, _opdiff(serre, lambda N: [], labs, scale=lambda N, A=A, B=B: abs(_compose(A, A, B)(N)[0][1])))
            elif i != j and (j - i) % h != 1:
                rels = max(rels, _opdiff(_compose(_Y(h, s, i), _Y(h, s, j)), _compose(_Y(h, s, j), _Y(h, s, i)), labs))
    rels = max(rels, _opdiff(_compose(*[_K(h, j) for j in range(h)]), lambda N: [(N, 1.0)], labs))
    # y-power obstructions of the other twists, on V (0-based: first end a = j-1, second end b = j of alpha_j).
    # e-type, order y^2, input v_b: the left side is pi(e_j) followed by the X-part of eq-dt-tv, v_b -> v_a -> v_{a-1},
    # non-zero; the right side is the X-part followed by kappa (x) pi(e_psi(j)), v_b -> v_{b-1} = v_a -> e_psi(j) v_a,
    # which vanishes for every kappa unless psi(j) = j - 1.  f-type, order y^-1, input v_{a'}, a' = first end of psi(j):
    # the right side (1 (x) pi(f_psi(j))) times the diagonal part is non-zero; the left side contains f_j v_{a'} = 0
    # unless psi(j) = j.  So exactly psi = pi (e-type) and psi = id (f-type) escape the obstruction.
    for j in range(h):
        b, a = j % h, (j - 1) % h
        for jp in range(h):
            lhs = bool(_e_rep(h, False, 1.0, j, b)); rhs = bool(_e_rep(h, False, 1.0, jp, (b - 1) % h))
            obstr_ok &= (lhs and not rhs) == (jp != (j - 1) % h)
            obstr_ok &= (((jp - 1) % h) != a) == (jp != j)
# h = 2 against ch. 37 and ch. 42
sS = 0.37; s2 = np.array([sS/2, -sS/2]); tau = -1j*np.exp(1j*np.pi*OMH/2); ydH = np.exp(OMH*TD)
dB = lambda n: tau**(-n)*QH**(n*(n - 1)/2 + n*sS/2)
devT2 = devR2 = 0; shifts_ok = True
for th in thsH:
    y = np.exp(OMH*th); X = np.exp(OMH*(th - TD - 1j*np.pi/2))
    for n in range(-3, 4):
        N = (n, -n)
        o1 = _T_dict(2, s2, X, False, {(0, N): 1.0}); o2 = _T_dict(2, s2, X, False, {(1, N): 1.0})
        A_ = QH**(sS/2)*o1[(N, 0)]; C_ = QH**(sS/2)*o1[((n + 1, -n - 1), 1)]*dB(n)/dB(n + 1)
        D_ = QH**(sS/2)*o2[(N, 1)]; B_ = QH**(sS/2)*o2[((n - 1, -n + 1), 0)]*dB(n)/dB(n - 1)
        sg37 = (QH**-n, -1j*(y/ydH)*QH**(n + sS), -1j*(y/ydH)*QH**-n, QH**(n + sS))   # A, B, C, D, eq-sg-transmission-matrix
        devT2 = max(devT2, max(abs(u_ - v_)/abs(v_) for u_, v_ in zip((A_, B_, C_, D_), sg37)))
for n in range(-3, 4):
    N = (n, -n)
    e1 = _Y(2, s2, 1)(N)[0]; k1 = _K(2, 1)(N)[0]; e0 = _Y(2, s2, 0)(N)[0]; k0 = _K(2, 0)(N)[0]
    aa = 1j*ydH*QH**(-sS)/(QH**2 - 1); bb = 1j*ydH*QH**sS/(QH**2 - 1)
    got = [o[1]*dB(n)/dB(o[0][0]) for o in (e1, k1, e0, k0)]
    want = (aa*QH**(-2*n), QH**(-2*n - 2 - sS), bb*QH**(2*n), QH**(2*n - 2 + sS))     # eq-dbd-second-module
    shifts_ok &= tuple(o[0][0] - n for o in (e1, k1, e0, k0)) == (1, 2, -1, -2)
    devR2 = max(devR2, max(abs(u_ - v_)/abs(v_) for u_, v_ in zip(got, want)))
report(f'(h) eq-dpt-second-module solves the e-type (psi = pi) and Cartan equations for V [{resV:.0e}] and Vbar '
       f'[{resVb:.0e}], h = 3, 4, 5, every node, 3 rapidities, exact labels; the sparse inverse of eq-dt-tvbar used in '
       f'the all-rank proof, H M_Vbar = 1 - X^h/q [{devH:.0e}]; with c_j q^(N+s)_j added it still solves them '
       f'[{resC:.0e}] but breaks k_j e_j k_j^-1 = q^2 e_j [{cartc:.2f}]; U_q(b^+) relations [{rels:.0e}]; the other twists '
       f'have a block free of unknowns [{obstr_ok}]; h = 2: q^(s/2) M_V = eq-sg-transmission-matrix [{devT2:.0e}] and '
       f'rho~ = eq-dbd-second-module [{devR2:.0e}] in the basis tau^-n q^(n(n-1)/2 + ns/2)|n>',
       resV < 1e-13 and resVb < 1e-13 and devH < 1e-13 and resC < 1e-13 and cartc > 0.1 and rels < 1e-13 and obstr_ok
       and devT2 < 1e-13 and devR2 < 1e-13 and shifts_ok)


# ---------------------------------------------------------------- (i) the soliton CDD condition eq-dpt-soliton-cdd
def _string(f, D, h):
    return np.prod([f(D - 1j*np.pi + 2j*np.pi*k/h) for k in range(h)], axis=0)


def _Fbar(f, D, h):
    return np.prod([f(D + 1j*np.pi*(2*m - h + 2)/h) for m in range(h - 1)], axis=0)


DsI = np.array([0.41, -0.7 + 0.2j, 1.3 - 0.1j])
devI = 0; okI = True
for h in range(2, 9):
    for f in (lambda D, h=h: Tcl(D, 1, h), lambda D: (np.sinh(D) + 0.3j)/(np.sinh(D) - 0.3j)):
        # crossing (ch. 41) with g -> g f, gbar -> gbar Fbar: (II) Fbar(D) f(D - i pi) = 1, (I) f(D) Fbar(D - i pi) = 1
        devI = max(devI, abs(_Fbar(f, DsI, h)*f(DsI - 1j*np.pi) - _string(f, DsI, h)).max(),
                   abs(f(DsI)*_Fbar(f, DsI - 1j*np.pi, h) - _string(f, DsI + 1j*np.pi - 2j*np.pi/h, h)).max())
    okI &= abs(_string(lambda D: Tcl(D, 1, h), DsI, h) - 1).max() < 1e-12
    if h > 2:
        okI &= abs(Tcl(DsI, 1, h)*Tcl(DsI + 1j*np.pi, 1, h) - 1).max() > 0.1
sinhf = lambda D: (np.sinh(D) + 0.3j)/(np.sinh(D) - 0.3j)
fail3 = abs(_string(sinhf, DsI, 3) - 1).max()
fail5 = abs(_string(sinhf, DsI, 5) - 1).max()
pass46 = max(abs(_string(sinhf, DsI, 4) - 1).max(), abs(_string(sinhf, DsI, 6) - 1).max())
okI &= abs(sinhf(DsI)*sinhf(DsI + 1j*np.pi) - 1).max() < 1e-14 and fail3 > 0.1 and fail5 > 0.1 and pass46 < 1e-12


def _constituents(h, a, om, D):
    ub = np.pi/2 - np.pi/(h*om)
    return ([D - 1j*ub + 1j*np.pi*(2*m - h + a + 1)/h for m in range(h - a)]
            + [D + 1j*ub + 1j*np.pi*(2*m - a + 1)/h for m in range(a)])


mult = [np.prod([Tcl(c, 1, 4) for c in _constituents(4, 2, om, 0.41)]) for om in (np.sqrt(2), np.sqrt(5))]
okI &= abs(mult[0] - (-0.2200 - 0.9755j)) < 1e-4 and abs(mult[1] - (0.2004 - 0.9797j)) < 1e-4
report(f'(i) eq-dpt-soliton-cdd: with g -> g f, both crossing relations of ch. 41 are the string condition '
       f'prod_k f(D - i pi + 2 pi i k/h) = 1 [{devI:.0e}]; T^cl_1 satisfies it for h = 2..8 but not f(D) f(D + i pi) = 1 '
       f'for h > 2 (the simpler condition is not necessary); (sinh D + 0.3i)/(sinh D - 0.3i) satisfies the latter, passes '
       f'the string condition at h = 4, 6 [{pass46:.0e}] (sufficient for even h) and fails it at h = 3, 5 by {fail3:.2f}, '
       f'{fail5:.2f} (not sufficient for odd h); f = T^cl_1 multiplies the factor of B^(2)_1 of a_3 at D = 0.41 by {mult[0]:.4f} at omega = sqrt2 and '
       f'{mult[1]:.4f} at sqrt5', okI and devI < 1e-12)


# ---------------------------------------------------------------- (j) tbl-dpt-spectrum: ranks, Lee-Yang, v = pi
def _sector_blocks(h, om, X, s):
    """Sector blocks of eq-dt-tv (basis v_i (x) |Q - e_i>, sum Q = 1) and eq-dt-tvbar (vb_i (x) |Qb + e_i>, sum Qb = -1)."""
    q = q_phys(om); I = np.eye(h)
    Q = I[1]; Qb = -I[0]
    MV = np.zeros((h, h), complex); MVb = np.zeros((h, h), complex)
    for i in range(h):
        N = Q - I[i]
        MV[i, i] += q**(-(N + s)[i]); MV[(i - 1) % h, i] += X
        N = Qb + I[i]; kk = np.zeros(h)
        for jj in range(h):
            b = (i + jj) % h; kk[b] += 1
            MVb[b, i] += X**jj*q**((N + s) @ kk)
    return MV, MVb


rngj = np.random.default_rng(4313)
okJ = True; worstJ = 0
for h in (3, 4, 5, 6):
    r = h - 1
    for om in (np.sqrt(2), 1.37, 3.1):
        s = rngj.normal(size=h)*0.3; s -= s.mean()
        X1 = np.exp(-1j*np.pi*om*r/h - 1j*np.pi/h)          # species-1 pole theta_1 (eq-dt-pole)
        Xr = np.exp(-1j*np.pi*om/h - 1j*np.pi/h)            # species-r pole, X^h = q
        sv1 = np.linalg.svd(_sector_blocks(h, om, X1, s)[0], compute_uv=False)
        svr = np.linalg.svd(_sector_blocks(h, om, Xr, s)[1], compute_uv=False)
        okJ &= sv1[-1]/sv1[0] < 1e-12 and sv1[-2]/sv1[0] > 0.1 and svr[1]/svr[0] < 1e-12
        worstJ = max(worstJ, sv1[-1]/sv1[0], svr[1]/svr[0])
resid = []
for h in (3, 4, 5):
    r = h - 1; om = 1.37; c = an_scalar_constant(r, om, q_phys(om))
    Dr = 1j*np.pi/2 - 1j*(np.pi/h + np.pi/(h*om))
    v1, v2 = [eps*an_scalar_gbar(Dr + eps, r, om, 400, c) for eps in (1e-5, 5e-6)]
    okJ &= abs(v1) > 1e-2 and abs(v1 - v2) < 1e-4*abs(v1)
    resid.append(v1)
fLY = lambda th: (np.sinh(th) + 1j*np.sin(2*np.pi/3))/(np.sinh(th) - 1j*np.sin(2*np.pi/3))
for b in (-0.5, 0.5, 1.5):
    for y0 in ((b + 1)/6,) + (((b - 1)/6,) if b > 1 else ()):
        th = 1j*np.pi*y0 + 1e-7
        okJ &= abs(1e-7*fLY(th + 1j*np.pi*(3 - b)/6)) > 1 and abs(1e-7*fLY(th - 1j*np.pi*(3 - b)/6)) < 1e-5
devPi = 0; okPi = True
for h in (3, 4, 5, 6):
    for a in range(1, h):
        g = np.pi*a/h
        devPi = max(devPi, max(abs(abs(Tcl(t - 1j*np.pi, a, h)) - 1) for t in np.linspace(-3, 3, 31)))
        okPi &= max(abs(abs(Tcl(t - 1j*(np.pi - 0.3), a, h)) - 1) for t in np.linspace(-3, 3, 31)) > 0.05
        okPi &= ((np.pi - g < np.pi < 1.5*np.pi - g) == (2*a < h)) and abs(1/Tcl(1j*g - 1j*np.pi + 1e-9, a, h)) < 1e-8
report(f'(j) tbl-dpt-spectrum: on a charge sector the residue of T_V at the species-1 pole has rank h - 1 and that of '
       f'T_Vbar at the species-r pole (X^h = q) rank 1, h = 3..6, three couplings, random offsets [{worstJ:.0e}]; gbar has '
       f'a simple uncancelled pole there [eps gbar = {abs(resid[0]):.3f}, {abs(resid[1]):.3f}, {abs(resid[2]):.3f}]; the '
       f'Lee-Yang bound-state poles are in Tbar_b, not T_b; for theta_d = i v, T_a is a phase for real theta at v = pi '
       f'[{devPi:.0e}] but not at v = pi - 0.3, and v = pi lies in the window iff a < h/2, pole at i gamma_a',
       okJ and okPi and devPi < 1e-12)

print(f'# total runtime {time.time() - T0:.0f} s')
