"""Shared helpers for the defect scripts (Part VIII, ch. 36-43); prints nothing.

CLASSICAL CONVENTIONS (book normalization: eq-imaginary-potential, eq-lax-pair, longest roots alpha^2 = 2;
checked in defect_classical_sg.py and defect_classical_an.py)

  fields       phi on x < 0, psi on x > 0 (both in R^r); the defect sits at x = 0.
  action       S = int dt { int_{-inf}^0 L(phi) dx + int_0^inf L(psi) dx
                            + [ 1/2 phi.E phi_t + phi.D psi_t + 1/2 psi.E psi_t - B_d(phi, psi) ]_{x=0} }
               type I_-:  D + D^T = 2, E = 1 - D antisymmetric (a_1: D = 1, E = 0).
               type I_+: D -> -D, E -> -E (the time reverse, equivalently the parity image, of type I_-).
  conditions   phi_x = E phi_t + D psi_t - grad_phi B_d,     psi_x = D^T phi_t - E psi_t + grad_psi B_d   (x = 0)
               (variation of S; matches BCZ hep-th/0305022 (1.5)-(1.6) and hep-th/0401020 (3.2) with F = E).
  a_r data     h = r + 1;  alpha_j = e_j - e_{j+1} in R^h (j = 0..r, indices mod h), so alpha_0 = -sum_{j>=1} alpha_j;
               lambda_i (i = 1..r) fundamental weights, alpha_i.lambda_j = delta_ij;  lambda_0 = lambda_{r+1} = 0;
               coordinates: an orthonormal basis U of the sum-zero subspace obtained by Gram-Schmidt on alpha_1..alpha_r
               (so for a_1: alpha_1 = sqrt2, alpha_0 = -sqrt2 as in sec-sg-dictionary).
               D = 2 sum_{i=1}^r (lambda_i - lambda_{i+1}) lambda_i^T,   D alpha_j = 2(lambda_j - lambda_{j+1}),
               D alpha_0 = -2 lambda_1,  D = 2(1 - pi)^{-1} with pi alpha_j = alpha_{j-1} (pi e_k = e_{k-1} on R^h),
               D^T alpha_j = -D alpha_{j-1},  alpha_i.D alpha_j = 2 delta_ij - 2 delta_{i,j+1}.
  potential    (real coupling; B_d^imag = B_d^real(beta -> i beta) as for eq-imaginary-potential, eq-boundary-potential)
               f_j = exp(b alpha_j.(D^T phi + D psi)/2),  g_j = exp(b alpha_j.D (phi - psi)/2),   b = beta or i beta,
               B_d^{I-}(sigma)  = (m/b^2) sum_{j=0}^r [ sigma (f_j - 1) + sigma^{-1} (g_j - 1) ],
               P_d^{I-}(sigma)  = -(m/b^2) sum_j [ sigma (f_j - 1) - sigma^{-1} (g_j - 1) ],
               B_d^{I+}(sigma) = B_d^{I-}(1/sigma),  P_d^{I+}(sigma) = -P_d^{I-}(1/sigma)   (D in f_j, g_j is always the type-I_- D).
               E + B_d and P + P_d are conserved, with P = -int phi_t.phi_x dx (momentum M sinh theta for eq-an-soliton).
               sigma = exp(-theta_d) for both types: a sine-Gordon defect carrying a 2 pi/beta_SG jump stores
               energy M cosh theta_d and momentum M sinh theta_d.
               sine-Gordon (a_1, imaginary coupling, beta_SG = sqrt2 beta, m0 = 2m):
               B_d = -(2 m0/beta_SG^2) [ sigma (cos(beta_SG(phi+psi)/2) - 1) + sigma^{-1} (cos(beta_SG(phi-psi)/2) - 1) ].
  integrability   1/2 |grad_phi B_d|^2 - 1/2 |grad_psi B_d|^2 = V(phi) - V(psi)  (sigma-independent).
  Lax          a_t = (A_+ + A_-)/2 of eq-boundary-lax-at; defect matrix L_d with  d_t L_d = L_d a_t(psi) - a_t(phi) L_d,
               L_d = exp(-b H.(E phi + D psi)/2) Ltilde exp(b H.(D^T phi - E psi)/2), and in the defining representation
               type I_-:  Ltilde = 1 + sigma lambda^{-1} sum_j E_{-alpha_j}   (BCZ hep-th/0401020 (3.39), d = sigma, c = 1),
               type I_+:  Ltilde^{-1} ~ 1 - sigma^{-1} lambda sum_j E_{alpha_j}  (Ltilde = sum_{k<h} (lambda/sigma)^k P^k)
               (checked in defect_classical_an.py; for a soliton of rapidity theta the Lax parameter is lambda = -+e^{-theta}).
  charges      both types conserve modified spin +-1 AND spin +-2 charges (a_2, a_3 checked; contrary to BCZ 2004 sec. 5).
  quantum      type I_- <-> U_q(b^-) Weyl module (BCZ hep-th/0506169 (3.51), Konik-LeClair); type I_+ <-> U_q(b^+).

QUANTUM CONVENTIONS (Part VI, boundary_common.py: Delta(e) = e(x)1 + k(x)e, Delta(f) = f(x)k^-1 + 1(x)f,
q = q_phys = -exp(-i pi lambda), y = exp(lambda theta) for sine-Gordon; checked in defect_sg_transmission.py)

  defect space W = span{|n>, n in Z}, topological charge 2n + s (k_1 = q^(2n+s), k_0 = k_1^-1).
  Weyl modules of the Borel subalgebras, one parameter y_d = exp(lambda theta_d) (plus the offset s):
     U_q(b^+) = <e_0, e_1, k^{+-1}>:  e_1 |n> = mu y_d |n+1>,  e_0 |n> = mu y_d |n-1>,
     U_q(b^-) = <f_0, f_1, k^{+-1}>:  f_1 |n> = mu/y_d |n-1>,  f_0 |n> = mu/y_d |n+1>,
     mu = 1/(2 sin(pi lambda)) = i/(q - q^-1) at q = q_phys (only mu^2 is basis independent; this value makes
     theta_d real and equal to BCZ's eta).  e_0 e_1 = e_1 e_0 (resp. f's) is central, so the q-Serre relations hold.
  transmission T(theta - theta_d): V(y) (x) W -> W (x) V(y) (soliton passes from left to right), intertwining the
     Part VI coproduct of the Borel subalgebra; triangle relation on V(y1) (x) V(y2) (x) W -> W (x) V(y2) (x) V(y1):
        (1 (x) Rc(y1/y2)) (T(y1) (x) 1) (1 (x) T(y2)) = (T(y2) (x) 1) (1 (x) T(y1)) (Rc(y1/y2) (x) 1),
     Rc = boundary_common.rcheck.  Matrix elements: T (v_+ (x) |n>) = A_n |n> (x) v_+ + C_n |n+1> (x) v_-,
     T (v_- (x) |n>) = B_n |n-1> (x) v_+ + D_n |n> (x) v_-.

Windows: W is truncated to n in [-N, N]; every module dict carries 'n' (the labels) so that equations can be imposed
only where the truncation is invisible (interior rows/columns).
"""
import numpy as np
from boundary_common import cop, prop_residual

kron = np.kron


# ======================================================================== classical: root data and D
def _t(typ):
    """'I-' (BCZ's defect, default) -> +1, 'I+' (its time reverse) -> -1; integers pass through."""
    return {'I-': 1, 'I+': -1}.get(typ, typ)


def an_data(r):
    """Root data of a_r^(1) in R^h (h = r+1) and in orthonormal coordinates on the sum-zero subspace.
    Returns dict: h, alpha_h (h x h, row j = alpha_j in R^h, j = 0..r), U (h x r), alpha (r+1 x r, rows alpha_j),
    lam (r+2 x r, rows lambda_0..lambda_{r+1}, lambda_0 = lambda_{r+1} = 0), pi (r x r)."""
    h = r + 1
    I = np.eye(h)
    alpha_h = np.array([I[j] - I[(j + 1) % h] for j in range(h)])
    Q, R = np.linalg.qr(alpha_h[1:].T)          # Gram-Schmidt on alpha_1..alpha_r
    Q = Q*np.sign(np.diag(R))                     # alpha_1 has positive first coordinate
    U = Q                                         # h x r
    alpha = alpha_h @ U                           # (r+1) x r
    simple = alpha[1:]                            # alpha_1..alpha_r
    lam_s = np.linalg.inv(simple).T               # rows lambda_i: alpha_i.lambda_j = delta_ij (alpha^2 = 2)
    lam = np.vstack([np.zeros(r), lam_s, np.zeros(r)])
    P = np.array([I[(k - 1) % h] for k in range(h)]).T        # P e_k = e_{k-1}
    pi = U.T @ P @ U
    return dict(h=h, alpha_h=alpha_h, U=U, alpha=alpha, lam=lam, pi=pi)


def dmatrix(r, typ=1):
    """(D, E) of the defect action in the orthonormal coordinates of an_data(r); typ=+1 type I_-, typ=-1 type I_+.
    a_1: D = 1, E = 0 (type I_-)."""
    typ = _t(typ)
    d = an_data(r)
    lam = d['lam']
    D = 2*sum(np.outer(lam[i] - lam[i + 1], lam[i]) for i in range(1, r + 1))
    E = np.eye(r) - D
    return (D, E) if typ == 1 else (-D, -E)


# ======================================================================== classical: potential and conditions
def defect_potential(phi, psi, sigma, r, beta=1.0, m=1.0, coupling='real', typ=1):
    """B_d, grad_phi B_d, grad_psi B_d, P_d for a_r^(1) (r = 1: sinh-/sine-Gordon), fields in the coordinates of
    an_data(r).  The exponents always use the type-I_- D; type I_+ has B_d^{I+}(sigma) = B_d^{I-}(1/sigma),
    P_d^{I+}(sigma) = -P_d^{I-}(1/sigma)."""
    typ = _t(typ)
    d = an_data(r)
    D, _ = dmatrix(r, 1)
    b = beta if coupling == 'real' else 1j*beta
    s = sigma if typ == 1 else 1/sigma
    phi = np.asarray(phi); psi = np.asarray(psi)
    B = 0; gphi = 0; gpsi = 0; Pd = 0
    for a in d['alpha']:
        Da, DTa = D @ a, D.T @ a
        f = np.exp(b*(Da @ phi + DTa @ psi)/2)
        g = np.exp(b*(DTa @ (phi - psi))/2)
        B = B + (m/b**2)*(s*(f - 1) + (g - 1)/s)
        Pd = Pd - (m/b**2)*(s*(f - 1) - (g - 1)/s)
        gphi = gphi + (m/(2*b))*(s*f*Da + g*DTa/s)
        gpsi = gpsi + (m/(2*b))*(s*f*DTa - g*DTa/s)
    if typ != 1:
        Pd = -Pd
    return B, gphi, gpsi, Pd


def bulk_potential(phi, r, beta=1.0, m=1.0, coupling='real'):
    """V and grad V of eq-imaginary-potential for a_r^(1) (n_j = 1)."""
    d = an_data(r)
    b = beta if coupling == 'real' else 1j*beta
    V = 0; gV = 0
    for a in d['alpha']:
        e = np.exp(b*(a @ phi))
        V = V + (m**2/b**2)*(e - 1)
        gV = gV + (m**2/b)*a*e
    return V, gV


def defect_dx(phit, psit, gphi, gpsi, r, typ=1):
    """phi_x, psi_x at x = 0 from the defect conditions."""
    D, E = dmatrix(r, _t(typ))
    return E @ phit + D @ psit - gphi, D.T @ phit - E @ psit + gpsi


def defect_condition_residual(phi, psi, phix, psix, phit, psit, sigma, r, beta=1.0, m=1.0, coupling='real', typ=1):
    """max |residual| of the two defect conditions at x = 0."""
    _, gphi, gpsi, _ = defect_potential(phi, psi, sigma, r, beta, m, coupling, typ)
    px, sx = defect_dx(np.asarray(phit), np.asarray(psit), gphi, gpsi, r, typ)
    return max(np.abs(np.asarray(phix) - px).max(), np.abs(np.asarray(psix) - sx).max())


# ======================================================================== classical: Lax pair (defining rep)
def an_lax_rep(r):
    """Defining rep of a_r^(1): E_j = e_{j,j+1} (j mod h) for alpha_j = e_j - e_{j+1}; H_k = e_kk."""
    h = r + 1
    Es = []
    for j in range(h):
        M = np.zeros((h, h)); M[j, (j + 1) % h] = 1.0
        Es.append(M)
    return Es


def a_t(lam, phi, phix, r, beta=1.0, m=1.0, coupling='real'):
    """a_t = (A_+ + A_-)/2 (eq-boundary-lax-at) in the defining rep; phi, phix in an_data(r) coordinates."""
    d = an_data(r)
    b = beta if coupling == 'real' else 1j*beta
    U = d['U']
    out = (b/2)*np.diag(U @ phix).astype(complex)
    for a, E in zip(d['alpha'], an_lax_rep(r)):
        out = out + (m/2)*np.exp(b*(a @ phi)/2)*(lam*E - E.T/lam)
    return out


# ======================================================================== quantum: Borel Weyl modules
def weyl_module_b(sign, q, yd, N, s=0.0, mu=1.0, kappa=1.0):
    """Shift (Weyl) module of U_q(b^+) (sign=+1: e_0, e_1, k_0, k_1) or U_q(b^-) (sign=-1: f_0, f_1, k_0, k_1)
    of U_q(sl2^), on W = span{|n>, n = -N..N}, in the boundary_common dict format.
      k_1 |n> = q^(2n+s) |n>,  k_0 = kappa k_1^-1  (kappa = k_0 k_1 is central; kappa = 1 is level zero);
      b^+: e_1 |n> = mu yd |n+1>,  e_0 |n> = mu yd |n-1>;      b^-: f_1 |n> = mu/yd |n-1>,  f_0 |n> = mu/yd |n+1>.
    Extra keys: 'n' (labels), 'interior' (|n| <= N-1), 'sign'."""
    ns = np.arange(-N, N + 1)
    d = len(ns)
    up = np.zeros((d, d), complex); dn = np.zeros((d, d), complex)
    for i in range(d - 1):
        up[i + 1, i] = 1           # |n> -> |n+1>
        dn[i, i + 1] = 1           # |n> -> |n-1>
    K1 = np.diag(q**(2.0*ns + s)).astype(complex)
    g = {('k', 1): K1, ('k', 0): kappa*np.linalg.inv(K1)}
    if sign == 1:
        g[('e', 1)] = mu*yd*up; g[('e', 0)] = mu*yd*dn
    else:
        g[('f', 1)] = mu*dn/yd; g[('f', 0)] = mu*up/yd
    g.update(dim=d, hw=None, n=ns, interior=np.abs(ns) <= N - 1, sign=sign)
    return g


def labels_VW(V, W):
    """defect label n of each basis vector of V (x) W and of W (x) V (V index slow / fast)."""
    nVW = np.repeat(W['n'][None, :], V['dim'], axis=0).ravel()
    nWV = np.repeat(W['n'][:, None], V['dim'], axis=1).ravel()
    return nVW, nWV


def borel_kinds(sign):
    return [('e', 0), ('e', 1), ('k', 0), ('k', 1)] if sign == 1 else [('f', 0), ('f', 1), ('k', 0), ('k', 1)]


def masked_intertwiners(pairs, row_mask, col_mask, tol=1e-7):
    """All X with (X A - B X)[r, c] = 0 for row_mask[r] and col_mask[c], for every (A, B) in pairs.
    Null space from the eigenvectors of the Gram matrix M^H M (relative singular values below tol).
    Returns (list of X, relative singular values)."""
    A0, B0 = pairs[0]
    dU, dW = A0.shape[0], B0.shape[0]
    IU, IW = np.eye(dU), np.eye(dW)
    sel = (np.outer(col_mask, row_mask)).ravel()            # column-major vec: index c*dW + r
    G = 0
    for A, B in pairs:
        Mg = (kron(A.T, IW) - kron(IU, B))[sel]
        G = G + Mg.conj().T @ Mg
    ev, Uv = np.linalg.eigh(G)
    rel = np.sqrt(np.clip(ev, 0, None)/ev[-1])
    vecs = [Uv[:, i] for i in range(len(ev)) if rel[i] < tol]
    return [v.reshape(dU, dW).T for v in vecs], rel


def transmission_window(V, W, tol=1e-7):
    """Borel intertwiners T: V (x) W -> W (x) V on the window, equations imposed on interior rows and columns.
    Returns (T, rank): T is the unique (up to scalar) solution restricted to the interior block (other entries 0),
    rank = dimension of the solution space after restriction to the interior block (absolute threshold 1e-5 on the
    singular values of the unit null vectors' interior parts); (None, 0) if there is none."""
    kinds = borel_kinds(W['sign'])
    pairs = [(cop(V, W, kd, j), cop(W, V, kd, j)) for kd, j in kinds]
    nVW, nWV = labels_VW(V, W)
    N = W['n'].max()
    rows, cols = np.abs(nWV) <= N - 1, np.abs(nVW) <= N - 1
    null, _ = masked_intertwiners(pairs, rows, cols, tol)
    blk = np.ix_(rows, cols)
    if not null:
        return None, 0
    stack = np.array([X[blk].ravel() for X in null])           # null vectors have unit norm: junk supported on the
    sv = np.linalg.svd(stack, compute_uv=False)                  # window edges has interior projection ~ 0
    rank = int((sv > 1e-5).sum())
    if rank == 0:
        return None, 0
    _, _, Vh = np.linalg.svd(stack)
    T = np.zeros_like(null[0])
    T[blk] = Vh[0].reshape(rows.sum(), cols.sum())
    return T, rank


def tmatrix_elements(T, V, W):
    """Read off A_n, B_n, C_n, D_n (dicts n -> value) from T: V (x) W -> W (x) V (V basis v_+ = 0, v_- = 1)."""
    dV, dW = V['dim'], W['dim']
    iVW = lambda a, k: a*dW + k
    iWV = lambda k, a: k*dV + a
    A, B, C, D = {}, {}, {}, {}
    for k, n in enumerate(W['n']):
        A[n] = T[iWV(k, 0), iVW(0, k)]
        D[n] = T[iWV(k, 1), iVW(1, k)]
        if k + 1 < dW:
            C[n] = T[iWV(k + 1, 1), iVW(0, k)]
        if k - 1 >= 0:
            B[n] = T[iWV(k - 1, 0), iVW(1, k)]
    return A, B, C, D


def tmatrix_from_elements(Af, Bf, Cf, Df, V, W):
    """Build T: V (x) W -> W (x) V from functions n -> A_n, B_n, C_n, D_n (window truncation)."""
    dV, dW = V['dim'], W['dim']
    T = np.zeros((dW*dV, dV*dW), complex)
    for k, n in enumerate(W['n']):
        T[k*dV + 0, 0*dW + k] = Af(n)
        T[k*dV + 1, 1*dW + k] = Df(n)
        if k + 1 < dW:
            T[(k + 1)*dV + 1, 0*dW + k] = Cf(n)
        if k - 1 >= 0:
            T[(k - 1)*dV + 0, 1*dW + k] = Bf(n)
    return T


def triangle_residual_braid(Rcheck_fn, T_fn, y1, y2, W, dV=2, margin=2):
    """Triangle relation on V(y1) (x) V(y2) (x) W -> W (x) V(y2) (x) V(y1):
         (1 (x) Rc(y1/y2)) (T(y1) (x) 1) (1 (x) T(y2))  ~  (T(y2) (x) 1) (1 (x) T(y1)) (Rc(y1/y2) (x) 1),
    compared (up to a scalar, prop_residual) on inputs and outputs with |n| <= N - margin.
    Rcheck_fn(z): V (x) V -> V (x) V at ratio z;  T_fn(y): V(y) (x) W -> W (x) V(y) on the window."""
    dW = W['dim']; N = W['n'].max()
    R = Rcheck_fn(y1/y2)
    T1, T2 = T_fn(y1), T_fn(y2)
    I2, IW = np.eye(dV), np.eye(dW)
    L = kron(IW, R) @ kron(T1, I2) @ kron(I2, T2)
    Rr = kron(T2, I2) @ kron(I2, T1) @ kron(R, IW)
    nin = np.tile(W['n'], dV*dV)                               # V (x) V (x) W: W index fastest
    nout = np.repeat(W['n'], dV*dV)                            # W (x) V (x) V: W index slowest
    blk = np.ix_(np.abs(nout) <= N - margin, np.abs(nin) <= N - margin)
    return prop_residual(L[blk], Rr[blk])


# ======================================================================== sine-Gordon transmission: scalar factors
# (appended for ch. 37; numpy/scipy versions of Gamma products, summed with three-level Richardson extrapolation)
def richardson_sum(terms, K):
    """sum_{k>=0} t_k from the partial sums at K, 2K, 4K (terms has length >= 4K, t_k = O(1/k^2)):
    removes the 1/K and 1/K^2 tails."""
    c = np.cumsum(terms)
    S1, S2, S4 = c[K - 1], c[2*K - 1], c[4*K - 1]
    R1, R2 = 2*S2 - S1, 2*S4 - S2
    return (4*R2 - R1)/3


def sg_transmission_scalar(D, lam, K=2000):
    """Scalar factor Phi(theta - theta_d) of the sine-Gordon soliton transmission matrix (eq-sg-transmission-matrix),
    s = 0: the Konik-LeClair minimal solution in the form of BCZ hep-th/0506169 (3.43)-(3.44) with gamma = lambda,
    theta~ = D:  Phi = e^{i pi (1+lam)/4} r(D)/((1 + i e^{lam D}) rbar(D)),
    r = prod_{k>=0} G(k lam + 1/4 - i w/2) G((k+1) lam + 3/4 - i w/2)/[G((k+1/2) lam + 1/4 - i w/2) G((k+1/2) lam + 3/4 - i w/2)],
    w = lam D/pi, rbar the same with -i w -> +i w (paired as in (3.44), which converges; (3.42)/(3.47) do not).
    |Phi|^2 (1 + e^{2 lam D}) = 1 for real D and Phi(D) Phi(D + i pi)(1 + e^{2 lam D}) = 1."""
    from scipy.special import loggamma as LG
    k = np.arange(4*K)
    w = lam*D/np.pi

    def logr(a):
        return richardson_sum(LG(k*lam + .25 + a) + LG((k + 1)*lam + .75 + a)
                              - LG((k + .5)*lam + .25 + a) - LG((k + .5)*lam + .75 + a), K)
    return np.exp(1j*np.pi*(1 + lam)/4 + logr(-1j*w/2) - logr(1j*w/2))/(1 + 1j*np.exp(lam*D))


def sg_S0(theta, lam, K=2000):
    """Sine-Gordon S_0 (eq-sg-s0) for complex theta, as Ghoshal-Zamolodchikov's a(u) = sin(lam(pi - u)) rho(u),
    u = -i theta (hep-th/9306002 (5.7); equal to the book S_0, checked in boundary_sine_gordon.py)."""
    from scipy.special import loggamma as LG
    u = -1j*theta
    l = np.arange(1, 4*K + 1)

    def logF(uu):
        v = lam*uu/np.pi
        return LG(2*l*lam - v) + LG(1 + 2*l*lam - v) - LG((2*l + 1)*lam - v) - LG(1 + (2*l - 1)*lam - v)
    s = richardson_sum(logF(u) + logF(np.pi - u) - logF(0) - logF(np.pi), K)
    logrho = np.log(-1/np.pi + 0j) + LG(lam) + LG(1 - lam*u/np.pi) + LG(1 - lam + lam*u/np.pi) + s
    return np.sin(lam*(np.pi - u))*np.exp(logrho)


def sg_tmatrix(theta, lam, q, theta_d, N, s=0.0, scalar=True):
    """The U_q(b^-) transmission matrix (type I_-, eq-sg-transmission-matrix) on the window n = -N..N:
    A_n = q^-n, D_n = q^(n+s), C_n = -i Y q^-n, B_n = -i Y q^(n+s), Y = e^{lam(theta - theta_d)}, times Phi (s = 0 form)
    if scalar.  Returns (T, V, W)."""
    from boundary_common import sl2_spin
    y, yd = np.exp(lam*theta), np.exp(lam*theta_d)
    V = sl2_spin(1, q, y)
    W = weyl_module_b(-1, q, yd, N, s=s)
    T = tmatrix_from_elements(lambda n: q**-n, lambda n: -1j*(y/yd)*q**(n + s), lambda n: -1j*(y/yd)*q**-n,
                              lambda n: q**(n + s), V, W)
    return (sg_transmission_scalar(theta - theta_d, lam) if scalar else 1.0)*T, V, W


def invariant_vector(ya, yb, q):
    """The vector of V(ya) (x) V(yb) (spin 1/2) annihilated by Delta(e_j), Delta(f_j), Delta(k_j) - 1 (the breather
    singlet of prp-bsg-breather-k); returns (relative smallest singular value, vector)."""
    from boundary_common import sl2_spin
    A, B = sl2_spin(1, q, ya), sl2_spin(1, q, yb)
    M = np.vstack([cop(A, B, kd, j) - (np.eye(4) if kd == 'k' else 0) for j in (0, 1) for kd in 'efk'])
    _, sv, Vh = np.linalg.svd(M)
    return sv[-1]/sv[0], Vh[-1].conj()


# ======================================================================== ch. 40: q-oscillators, Fock modules, L-operators
# (appended for ch. 40; checked in defect_universal_r.py and defect_transfer.py)
def q_oscillator_module(kind, q, N, s=0.0):
    """Representations of the q-oscillator algebra eq-q-oscillator (BLZ hep-th/9604044 (2.27) normalization, book q):
         q a+ a- - q^-1 a- a+ = 1/(q - q^-1),   q^h a+- q^-h = q^(+-2) a+-,    q^h |n> = q^(2n+s) |n>,
      'weyl': n = -N..N,  a+ |n> = |n+1>/(q-q^-1),  a- |n> = |n-1>/(q-q^-1)            (Z = a+ a- - (q-q^-1)^-2 = 0),
      'fock': n = 0..N,   a+ |n> = |n+1>/(q-q^-1),  a- |n> = (1 - q^(2n)) |n-1>/(q-q^-1)  (a- |0> = 0, Z |n> = -q^(2n)/(q-q^-1)^2).
    Returns dict ap, am, h (array, the eigenvalues of h), n (labels), interior (labels unaffected by the truncation:
    |n| <= N-1 for 'weyl', n <= N-1 for 'fock'), kind."""
    kap = q - 1/q
    ns = np.arange(-N, N + 1) if kind == 'weyl' else np.arange(0, N + 1)
    d = len(ns)
    ap = np.zeros((d, d), complex); am = np.zeros((d, d), complex)
    for i in range(d - 1):
        ap[i + 1, i] = 1/kap
        am[i, i + 1] = (1/kap) if kind == 'weyl' else (1 - q**(2*ns[i + 1]))/kap
    interior = (np.abs(ns) <= N - 1) if kind == 'weyl' else (ns <= N - 1)
    return dict(ap=ap, am=am, h=2.0*ns + s, n=ns, interior=interior, kind=kind)


def fock_module_b(sign, q, yW, N, s=0.0, kind='fock'):
    """U_q(b^+) (sign=+1) or U_q(b^-) (sign=-1) module obtained from q_oscillator_module through the homomorphisms
         b^+:  e_1 = i yW a+,  e_0 = i yW a-,      b^-:  f_1 = (i/yW) a-,  f_0 = (i/yW) a+,      k_1 = q^h = k_0^-1
    (eq-weyl-module, eq-fock-module).  kind='fock': the Fock module (BLZ's Q-operators); kind='weyl' reproduces
    weyl_module_b(sign, q, yW, N, s) exactly (mu = i/(q - q^-1)).  Dict in the boundary_common format plus the keys
    n, interior, sign, h, ap, am, yW, q."""
    O = q_oscillator_module(kind, q, N, s)
    K1 = np.diag(q**O['h']).astype(complex)
    g = {('k', 1): K1, ('k', 0): np.linalg.inv(K1)}
    if sign == 1:
        g[('e', 1)] = 1j*yW*O['ap']; g[('e', 0)] = 1j*yW*O['am']
    else:
        g[('f', 1)] = 1j*O['am']/yW; g[('f', 0)] = 1j*O['ap']/yW
    g.update(dim=len(O['n']), hw=None, n=O['n'], interior=O['interior'], sign=sign, h=O['h'], ap=O['ap'], am=O['am'],
             yW=yW, q=q)
    return g


def borel_l_operator(sign, y, W):
    """Closed-form L-operator on V(y) (x) W (V index slow), W from fock_module_b (Weyl or Fock), up to a scalar
    (eq-l-operator):  sign=+1: (pi_y (x) rho)(R), the U_q(b^+) intertwiner;  sign=-1: (pi_y (x) rho)(R_21^-1), U_q(b^-).
    With kappa = q - q^-1 and Z = a+ a- - kappa^-2 (Z = 0 on Weyl modules):
      b^+:  [[q^(h/2), i kappa x q^(h/2) a-], [i kappa x q^(-h/2) a+, q^(-h/2) (1 - kappa^2 x^2 Z)]],         x = yW/y,
      b^-:  [[(1 - kappa^2 x^2 q^2 Z) q^(-h/2), -i kappa x a- q^(h/2)], [-i kappa x a+ q^(-h/2), q^(h/2)]],     x = y/yW.
    On Weyl modules these are q^(+-s/2) times the matrices of weyl-module intertwiners in sg_tmatrix / the conventions
    brief; T = flip . L : V (x) W -> W (x) V is the transmission matrix."""
    q = W['q']; kap = q - 1/q
    ap, am = W['ap'], W['am']
    I = np.eye(W['dim'])
    Z = ap @ am - I/kap**2
    qh, qhi = np.diag(q**(W['h']/2)), np.diag(q**(-W['h']/2))
    if sign == 1:
        x = W['yW']/y
        blocks = [[qh, 1j*kap*x*qh @ am], [1j*kap*x*qhi @ ap, qhi @ (I - kap**2*x**2*Z)]]
    else:
        x = y/W['yW']
        blocks = [[(I - kap**2*x**2*q**2*Z) @ qhi, -1j*kap*x*am @ qh], [-1j*kap*x*ap @ qhi, qh]]
    E = lambda a, b: np.outer(np.eye(2)[a], np.eye(2)[b])
    return sum(kron(E(a, b), blocks[a][b]) for a in range(2) for b in range(2))


def flip_VW(dV, dW):
    """The flip V (x) W -> W (x) V (V index slow on the input, W index slow on the output)."""
    P = np.zeros((dW*dV, dV*dW))
    for a in range(dV):
        for k in range(dW):
            P[k*dV + a, a*dW + k] = 1
    return P
