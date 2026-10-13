"""Dirichlet-type boundary conditions of a_n^(1) affine Toda theory on the half-line (Part VII, ch. 30,
sec-soliton-preserving, prp-soliton-preserving, eq-theta-p-reflection, eq-boundary-compatibility; ch. 31,
sec-images-preserving; ch. 35, sec-boundary-open-problems, items 7 and 11; appendix D, item Delius 1998b (1.4)).
Written in the consolidation of half-line open problem 11 (classical claims 4 and 5), 2026-10-13.

Book conventions: a_n^(1) in the defining representation with cyclic labels (boundary_classical.an_rep:
alpha_j = e_j - e_{j+1}, E_j = e_{j,j+1}, E_{-alpha_j} = E_j^T, H_k = e_kk, fields in the sum-zero subspace of R^h);
a_t of eq-boundary-lax-at (boundary_classical.a_t; imaginary coupling by beta -> i beta); a_x of the same equation;
compatibility eq-boundary-compatibility K a_t(lam) + sigma(a_t(1/lam)) K = 0 at x = 0; double-row transfer matrix
eq-double-row with K_- = 1 at a vacuum; solitons phi = (i/beta) sum_j alpha_j ln tau_j (boundary_images.Sol).
A coweight class nu in Lambda/2Lambda is represented by eps_S = sum_{k in S} e_k - |S|/h, S a subset of {1,...,n};
D_nu = e^{i pi nu.H} = diag(e^{i pi nu_k}).
Checks:
 1. sigma = Hermitian conjugation, imaginary coupling, Re phi(0) = (pi/beta) nu + (2 pi/beta) x coweight, Im phi(0)
    free, d_x Re phi(0) free, d_x Im phi(0) = 0: K_+ = D_nu solves the compatibility condition for all 2^n classes;
    K_+ = 1 only for nu in 2 Lambda (prp-soliton-preserving). n = 1, 2, 3.
 2. Classification for these boundary data: the constant K (any h x h matrix) is unique up to a scalar and equals D_nu
    whenever Re phi(0) lies in (pi/beta) Lambda; no K exists when Re phi(0) is moved off (pi/beta) Lambda, or without
    the Neumann condition on Im phi. n = 1, 2, 3.
 3. D_nu (all classes) solves the classical reflection equation eq-classical-reflection with the book's r-matrix at
    real spectral parameters (r real, so r^sigma = r^t); so does any constant diagonal K_+ (the sign pattern is selected by
    the compatibility condition alone); a random non-diagonal constant K does not. N = 2, 3, 4.
 4. Static self-mirror solitons (species b, theta = 0, xi = i y): phi(x) + conj phi(-x) = 2 Re phi(0) for all x, with
    beta Re phi(0)/pi a coweight NOT in 2 Lambda, and d_x Im phi(0) = 0; D_nu is compatible at x = 0, K = 1 is not.
    Re of the half-line energy is half the full-line energy, which is the soliton mass 2 h m_b/beta^2. a_1 (y = pi/2:
    the real sine-Gordon half-kink, phi(0) = pi/beta_SG, the maximum of the potential), a_2 (b = 1, 2), a_3 (b = 1, 2).
 5. Integrability on exact solutions: the three-soliton image solutions (a, th, xi1), (a, -th, xi2), (b, 0, xi3) with
    Re xi3 = -ln A_ab(th), xi2 = -conj(xi1) - ln A_aa(2th) - ln A_ab(th) satisfy the boundary condition of check 4 at
    all times, and tau(lam) = tr(K_+ T(lam) T(1/lam)^dagger) on [-16, 0] is conserved with K_+ = D_nu and not with
    K_+ = 1. Control: the book's Theta P pair (eq-thetap-mirror; trivial class) conserves tau with K_+ = 1.
    a_1 (including the real kink reflected at the sine-Gordon Dirichlet boundary at the maximum), a_2, a_3 (regular phases).
 6. Real coupling, Dirichlet phi(0,t) = phi_0 on all components, field-independent K, defining representation: an
    invertible K exists for a_1 (sinh-Gordon) and, for a_2, a_3, a_4, for no linear sigma (transposition, antipode, and
    both composed with every Z_h rotation and with the diagram reflection), at generic phi_0 and at phi_0 = 0.
    On real data antilinear sigma (dagger, -conj) acts like the linear ones: no K either.
 7. Other representations: adjoint (a_2, a_3), Lambda^2 (a_3), V + Vbar (a_2) and the vector of d_4^(1): no K;
    a_1 adjoint: K exists.
 8. Imaginary coupling, complex fields, Dirichlet phi(0) = phi_0 (complex) on all components: no K for transposition,
    antipode, dagger or -conj (a_2, a_3); for a_1, K exists for the linear ones only.
 9. The obstruction is the cycle condition: at lam = e^{i pi/h} (lam^(2h) = 1) a K exists (antipode, defining rep;
    both types in the self-conjugate adjoint(a_2) and Lambda^2(a_3)); at lam = e^{i pi/(2h)} it does not.
10. Independent of Sklyanin's formalism: the spin-(+-2) current of a_n^(1), n = 2, 3, 4, is unique, T = A(p,p,p) + B(q,p),
    Theta = -B(p, grad V) (p = d_+ phi, q = d_+^2 phi; ansatz of Corrigan, Dorey, Rietdijk, Sasaki hep-th/9404108 sec. 2),
    and d_- T = d_+ Theta holds on an exact two-soliton solution. At a Dirichlet wall (u = d_x phi(0)) the
    energy-like flux 2A(u,u,u) - 4B(u, grad V) has no u-dot term and is nonzero, and the momentum-like flux
    4B(u-dot, u) is not a time derivative (B antisymmetric, nonzero): no combination of the spin +-2 charges is
    conserved. a_1 has no spin-2 current.
11. Delius's own conditions (5.1)-(5.2) of hep-th/9809140 (his Lax pair (3.1), field phi_D = i beta phi), with his
    image condition (1.4) for an arbitrary coweight lambda (Im phi_D(0) = -pi lambda, d_x Re phi_D(0) = 0): solved by
    K_+ = e^{i pi lambda.H} for every coweight, by his K_+ = 1 only for lambda in 2 Lambda, which is when (1.4) implies (1.5).
    n = 1, 2, 3.
Runtime: about two minutes.
"""
import importlib.util
import itertools
import os

import numpy as np
import scipy.linalg as sla

from boundary_common import report


# ---------------------------------------------------------------- book definitions (the book scripts run their checks on import)
def _book_source(name):
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), name + '.py')
    if os.path.exists(here):
        return open(here).read()
    return open(importlib.util.find_spec(name).origin).read()


def _block(src, start, stop):
    i = 0 if start is None else src.index(start)
    return src[i:src.index(stop, i + 1)]


def load_classical():
    src = _book_source('boundary_classical')
    code = _block(src, None, '\n# check the labels')
    code += _block(src, '# ---------------------------------------------------------------- 2. K-matrix compatibility', '\nsigns = lambda r')
    code += _block(src, '# ---------------------------------------------------------------- 8. Sklyanin', '\npairs = [(0.7, 1.3)')
    ns = {'__name__': 'boundary_classical_defs'}
    exec(compile(code, 'boundary_classical.py', 'exec'), ns)
    return ns


def load_images():
    src = _book_source('boundary_images')
    ns = {'__name__': 'boundary_images_defs'}
    exec(compile(_block(src, None, '\n# ---------------------------------------------------------------- 1. boundary conditions'),
                 'boundary_images.py', 'exec'), ns)
    return ns


BC = load_classical()                 # beta = 0.83, m = 1.27
BS = load_classical()                 # a second copy, set to the soliton conventions of boundary_images
IM = load_images()                    # beta = 0.9, m = 1.0
BS['beta'], BS['m'] = IM['beta'], IM['m']
algebra, a_t, eij = BC['algebra'], BC['a_t'], BC['eij']
Sol, interaction, roots_an = IM['Sol'], IM['interaction'], IM['roots_an']
rng = np.random.default_rng(2026)
trap = getattr(np, 'trapezoid', None) or np.trapz


def eps_S(S, h):
    v = -len(S) / h * np.ones(h); v[list(S)] += 1
    return v


def classes(n):
    return [S for k in range(n + 1) for S in itertools.combinations(range(1, n + 1), k)]


def solve_K(As, Bs):
    """nullity of {K : K A + B K = 0 for all (A, B)}, the last right-singular vector, and whether a generic null K is invertible."""
    N = As[0].shape[0]
    rows = []
    for A, B in zip(As, Bs):
        cols = []
        for p in range(N):
            for q in range(N):
                K = np.zeros((N, N), complex); K[p, q] = 1
                cols.append((K @ A + B @ K).ravel())
        rows.append(np.array(cols).T)
    U, sv, Vh = np.linalg.svd(np.vstack(rows))
    nul = int((sv < 1e-10 * sv[0]).sum())
    if nul == 0:
        return 0, None, False
    Kr = sum(rng.normal() * Vh[-1 - k].conj() for k in range(nul)).reshape(N, N)
    s = np.linalg.svd(Kr, compute_uv=False)
    return nul, Vh[-1].conj().reshape(N, N), bool(s[-1] > 1e-7 * s[0])


# ---------------------------------------------------------------- 1, 2. sigma = dagger: Dirichlet on Re phi, Neumann on Im phi
def sp_data(n, v, lam, dIm=0.0, ns=BC):
    roots, Es, Hs, nj, sub = ns['algebra'](f'a{n}(1)')
    w, u, s = (sub @ rng.normal(size=sub.shape[1]) for _ in range(3))
    phi, dx = v + 1j * w, u + 1j * dIm * s
    A = ns['a_t'](lam, phi, dx, roots, Es, Hs, nj, coupling='imag')
    B = ns['a_t'](1 / lam, phi, dx, roots, Es, Hs, nj, coupling='imag').conj().T
    return A, B


def resid(K, A, B):
    return np.linalg.norm(K @ A + B @ K) / (np.linalg.norm(K) * (np.linalg.norm(A) + np.linalg.norm(B)))


beta = BC['beta']
ok1 = True; worstD = 0; best1 = np.inf; info2 = []; ok2 = True
for n in (1, 2, 3):
    h = n + 1
    for S in classes(n):
        nu = eps_S(S, h)
        v = (np.pi / beta) * nu + (2 * np.pi / beta) * eps_S(tuple(rng.choice(h, size=rng.integers(0, h), replace=False)), h)
        D = np.diag(np.exp(1j * np.pi * nu))
        rD = max(resid(D, *sp_data(n, v, lam)) for lam in (0.6, 1.7))
        r1 = max(resid(np.eye(h), *sp_data(n, v, lam)) for lam in (0.6, 1.7))
        worstD = max(worstD, rD)
        if S:
            best1 = min(best1, r1)
        else:
            ok1 &= r1 < 1e-13
        # 2. classification
        nul, K, inv = solve_K(*zip(*[sp_data(n, v, 0.6) for _ in range(6)]))
        ok2 &= nul == 1 and inv and np.linalg.norm(K / K[0, 0] - D / D[0, 0]) < 1e-12
        pert = rng.normal(size=h); pert -= pert.mean()
        nul_p = solve_K(*zip(*[sp_data(n, v + 0.05 * pert, 0.6) for _ in range(6)]))[0]
        nul_d = solve_K(*zip(*[sp_data(n, v, 0.6, dIm=1.0) for _ in range(6)]))[0]
        ok2 &= nul_p == 0 and nul_d == 0
report(f'sigma = dagger, Re phi(0) = (pi/beta) nu (mod (2pi/beta) Lambda), Neumann on Im phi: K_+ = e^{{i pi nu.H}} is compatible for all '
       f'2^n classes, K_+ = 1 only for nu in 2 Lambda (n = 1, 2, 3)  [residuals: D {worstD:.1e}; K = 1 off 2 Lambda >= {best1:.2f}]',
       ok1 and worstD < 1e-13 and best1 > 0.1)
report('classification (sigma = dagger, these boundary data, n = 1, 2, 3): the constant K is unique up to scale and equals e^{i pi nu.H} '
       'for Re phi(0) in (pi/beta) Lambda; no K for Re phi(0) off (pi/beta) Lambda or without Neumann on Im phi', ok2)

# ---------------------------------------------------------------- 3. classical reflection equation
RE_plus, r_toda = BC['RE_plus'], BC['r_toda']
worst3 = 0; bad3 = np.inf
for N in (2, 3, 4):
    rf = lambda z, N=N: r_toda(N, z)
    for S in itertools.chain.from_iterable(itertools.combinations(range(N), k) for k in range(N)):
        D = np.diag(np.exp(1j * np.pi * eps_S(S, N)))
        for l1, l2 in [(0.7, 1.3), (1.9, 0.45), (2.2, -0.6)]:
            worst3 = max(worst3, np.linalg.norm(RE_plus(lambda l: D, rf, l1, l2, N)) / np.linalg.norm(np.kron(D, D)))
    Dr = np.diag(rng.normal(size=N) + 1j * rng.normal(size=N))           # any constant diagonal K_+ solves it as well
    worst3 = max(worst3, np.linalg.norm(RE_plus(lambda l: Dr, rf, 0.7, 1.3, N)) / np.linalg.norm(np.kron(Dr, Dr)))
    Kr = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    bad3 = min(bad3, np.linalg.norm(RE_plus(lambda l: Kr, rf, 0.7, 1.3, N)) / np.linalg.norm(np.kron(Kr, Kr)))
report(f'e^{{i pi nu.H}} (all classes), and indeed any constant diagonal K_+, solves eq-classical-reflection with the book r-matrix at real '
       f'spectral parameters, N = 2, 3, 4; a random non-diagonal constant K does not  [{worst3:.1e}; random {bad3:.2f}]', worst3 < 1e-13 and bad3 > 0.1)


# ---------------------------------------------------------------- solitons: fields, Lax pair, transport
bS, mS = IM['beta'], IM['m']


def taus_t(sol, x, t):
    """d tau_j/dt, same expansion as Sol.taus."""
    x = np.atleast_1d(np.asarray(x, float))
    k = sol.ma[:, None] * np.cosh(sol.th)[:, None]; w = sol.ma[:, None] * np.sinh(sol.th)[:, None]
    expo = k * x[None, :] - w * t + sol.xi[:, None]
    dT = np.zeros((sol.h, x.size), complex)
    for S, c in zip(sol.subsets, sol.coef):
        if S:
            e = c * np.exp(expo[list(S)].sum(0)); dw = -w[list(S)].sum(0)
            for j in range(sol.h):
                dT[j] += sol.om ** (j * sum(sol.sp[i] for i in S)) * dw * e
    return dT


def fields(sol, x, t):
    ph, dph, T = sol.phi(x, t)
    dT = taus_t(sol, x, t)
    al = roots_an(sol.n)
    return ph, dph, (1j / bS) * sum(al[j][:, None] * (dT[j] / T[j])[None, :] for j in range(sol.h)), T


def a_x(lam, phi, dtphi, roots, Es, Hs, nj):
    """x-component of eq-boundary-lax-at at imaginary coupling (beta -> i beta)."""
    b = 1j * bS
    out = (b / 2) * sum(dtphi[k] * Hs[k] for k in range(len(Hs)))
    for a, E, n in zip(roots, Es, nj):
        out = out + (mS / 2) * np.sqrt(n * (a @ a) / 2) * np.exp(b * (a @ phi) / 2) * (lam * E + E.T / lam)
    return out


def transport(sol, t, lam, L=16.0, N=1600):
    """T(lam) from x = -L to 0, dT/dx = -a_x T (book convention), fourth-order Magnus with two Gauss points per step."""
    roots, Es, Hs, nj, sub = algebra(f'a{sol.n}(1)')
    hs = L / N; xs = -L + hs * np.arange(N); g = hs / (2 * np.sqrt(3))
    pts = np.empty(2 * N); pts[0::2] = xs + hs / 2 - g; pts[1::2] = xs + hs / 2 + g
    ph, _, dtph, _ = fields(sol, pts, t)
    T = np.eye(sol.h, dtype=complex)
    for k in range(N):
        A1 = -a_x(lam, ph[:, 2 * k], dtph[:, 2 * k], roots, Es, Hs, nj)
        A2 = -a_x(lam, ph[:, 2 * k + 1], dtph[:, 2 * k + 1], roots, Es, Hs, nj)
        T = sla.expm((hs / 2) * (A1 + A2) + (np.sqrt(3) * hs**2 / 12) * (A2 @ A1 - A1 @ A2)) @ T
    return T


def compat_at0(sol, t, lam, K):
    roots, Es, Hs, nj, sub = BS['algebra'](f'a{sol.n}(1)')
    ph, dph, _, _ = fields(sol, np.linspace(-25, 0, 5001), t)
    A = BS['a_t'](lam, ph[:, -1], dph[:, -1], roots, Es, Hs, nj, coupling='imag')
    B = BS['a_t'](1 / lam, ph[:, -1], dph[:, -1], roots, Es, Hs, nj, coupling='imag').conj().T
    return resid(K, A, B)


def nu_of(v, n):
    al = roots_an(n)
    c = np.array([bS * (al[k] @ v) / np.pi for k in range(n + 1)])          # alpha_k . nu
    return bS * v / np.pi, np.max(np.abs(c - np.round(c))), bool(np.all(np.round(c).astype(int) % 2 == 0))


def V_imag(phi, n):
    return -(mS**2 / bS**2) * sum(np.exp(1j * bS * np.einsum('i,i...->...', a, phi)) - 1 for a in roots_an(n))


# ---------------------------------------------------------------- 4. static self-mirror solitons (half-solitons at the wall)
xs_full = np.linspace(-25, 25, 20001); i0 = 10000; dxs = xs_full[1] - xs_full[0]
ok4 = True; info4 = []
for n, b, y in [(1, 1, np.pi / 2), (1, 1, 0.7), (2, 1, 0.9), (2, 2, -1.3), (3, 1, 0.4), (3, 2, 0.8)]:
    s = Sol(n, [b], [0.0], [1j * y])
    ph, dph, T = s.phi(xs_full, 0.0)
    c = ph + np.conj(ph[:, ::-1]); v = ph[:, i0].real
    nu, off, triv = nu_of(v, n)
    D = np.diag(np.exp(1j * np.pi * nu))
    img = max(np.max(np.abs(c - c[:, [i0]])), np.max(np.abs(c[:, i0] - 2 * v)))
    cD, c1 = compat_at0(s, 0.0, 0.7, D), compat_at0(s, 0.0, 0.7, np.eye(n + 1))
    dens = 0.5 * np.einsum('ix,ix->x', dph, dph) + V_imag(ph, n)              # static energy density (eq-boundary-energy, B = 0)
    Efull = trap(dens, dx=dxs); Ehalf = trap(dens[:i0 + 1], dx=dxs)
    ok4 &= img < 1e-12 and np.max(np.abs(dph[:, i0].imag)) < 1e-12 and off < 1e-12 and not triv and cD < 1e-13 and c1 > 0.1
    ok4 &= abs(Efull.imag) < 1e-8 * abs(Efull) and abs(Ehalf.real - Efull.real / 2) < 1e-8 * abs(Efull)
    ok4 &= abs(Efull.real - 2 * (n + 1) * 2 * mS * np.sin(np.pi * b / (n + 1)) / bS**2) < 1e-6 * abs(Efull)   # M_b = 2 h m_b/beta^2
    if n == 1 and abs(y - np.pi / 2) < 1e-12:                                 # the real sine-Gordon half-kink
        ok4 &= np.max(np.abs(ph.imag)) < 1e-12 and abs(Ehalf.imag) < 1e-8
    info4.append(f'a{n} b={b} y={y:.2f}: nu={np.round(nu, 3)}, E={Efull.real:.3f}, Im E(x<0)={Ehalf.imag:.3f}')
report('static self-mirror solitons: phi(x) + conj phi(-x) = 2 Re phi(0) with beta Re phi(0)/pi a coweight outside 2 Lambda, '
       'd_x Im phi(0) = 0; e^{i pi nu.H} compatible, K = 1 not; Re E(x<0) = E/2, E real (half-soliton at the wall; '
       f'a1 = sine-Gordon half-kink)  [{"; ".join(info4)}]', ok4)


# ---------------------------------------------------------------- 5. tau conservation along exact three-soliton image solutions
def three(n, a, b, th, xi1, y3):
    h = n + 1
    Aab = interaction(a, b, th, 0.0, h).real; Aaa = interaction(a, a, th, -th, h).real
    return Sol(n, [a, a, b], [th, -th, 0.0], [xi1, -np.conj(xi1) - np.log(Aaa) - np.log(Aab), -np.log(Aab) + 1j * y3])


def tau_variation(sol, Ks, lams=(0.6, 1.5), ts=(-6.0, -2.0, 2.0, 6.0)):
    """max over lam of max_t |tau(t) - tau(t_0)|/|tau(t_0)|, for each K_+ in Ks (transports computed once)."""
    out = np.zeros(len(Ks))
    for lam in lams:
        TT = [transport(sol, t, lam) @ transport(sol, t, 1 / lam).conj().T for t in ts]
        for i, K in enumerate(Ks):
            vals = np.array([np.trace(K @ M) for M in TT])
            out[i] = max(out[i], np.max(np.abs(vals - vals[0])) / abs(vals[0]))
    return out


ok5 = True; cons = 0; noncons = np.inf; worst_img = 0
xs_img = np.linspace(-25, 25, 10001); j0 = 5000
for (n, a, b, th, xi1, y3) in [(1, 1, 1, 0.5, 1.0 + 0.5j * np.pi, 0.5 * np.pi), (1, 1, 1, 0.5, 1.0 + 0.7j, 1.1), (2, 1, 1, 0.5, 1.0 + 2.094j, 2.094), (2, 2, 1, 0.6, 1.0 + 0.0j, 4.189),
                               (3, 1, 2, 0.5, 1.0 + 5.76j, 1.571), (3, 2, 1, 0.45, 1.0 + 1.571j, 5.76)]:
    s = three(n, a, b, th, xi1, y3)
    vals = []
    for t in np.linspace(-8, 8, 9):
        ph, dph, T = s.phi(xs_img, t)
        c = ph + np.conj(ph[:, ::-1])
        worst_img = max(worst_img, np.max(np.abs(c - c[:, [j0]])), np.max(np.abs(dph[:, j0].imag)))
        vals.append(ph[:, j0].real)
    nu, off, triv = nu_of(vals[0], n)
    ok5 &= np.ptp(np.array(vals), axis=0).max() < 1e-12 and off < 1e-12 and not triv
    D = np.diag(np.exp(1j * np.pi * nu))
    vD, v1 = tau_variation(s, [D, np.eye(n + 1)])
    cons = max(cons, vD); noncons = min(noncons, v1)
ctrl = 0
for (n, a, xi1) in [(2, 1, 0.4 + 0.25j), (3, 2, 0.1 + 0.6j)]:
    th = 0.6; Aaa = interaction(a, a, th, -th, n + 1).real
    ctrl = max(ctrl, tau_variation(Sol(n, [a, a], [th, -th], [xi1, -np.conj(xi1) - np.log(Aaa)]), [np.eye(n + 1)])[0])
report(f'exact three-soliton image solutions (soliton + mirror + half-soliton; a1, a2, a3) obey the new condition at all times, and tau(lam) is '
       f'conserved with K_+ = e^{{i pi nu.H}}, not with K_+ = 1; control: Theta P pair, K_+ = 1  '
       f'[image {worst_img:.1e}; rel. variation D {cons:.1e}, K=1 >= {noncons:.1e}; control {ctrl:.1e}]',
       ok5 and worst_img < 1e-12 and cons < 1e-6 and noncons > 1e-2 and ctrl < 1e-6)


# ---------------------------------------------------------------- 6-9. Dirichlet on all components
SIG = {'T': lambda X: X.T, 'S': lambda X: -X, 'dag': lambda X: X.conj().T, 'nconj': lambda X: -X.conj()}


def dirichlet_K(rep, phi0, sig, lam, coupling='real', P=None, cplx=False, nsamp=6):
    """nullity / invertibility of {K : K a_t(lam) + sigma(a_t(1/lam)) K = 0} at phi(0) = phi0, d_x phi(0) free."""
    roots, Es, Hs, nj, sub = rep
    s0 = SIG[sig]
    s = s0 if P is None else (lambda X: P @ s0(X) @ np.linalg.inv(P))
    As, Bs = [], []
    for _ in range(nsamp):
        dx = sub @ rng.normal(size=sub.shape[1])
        if cplx:
            dx = dx + 1j * (sub @ rng.normal(size=sub.shape[1]))
        As.append(a_t(lam, phi0, dx, roots, Es, Hs, nj, coupling)); Bs.append(s(a_t(1 / lam, phi0, dx, roots, Es, Hs, nj, coupling)))
    nul, _, inv = solve_K(As, Bs)
    return nul if inv else -nul


def cyclic(h, k):
    P = np.zeros((h, h))
    for i in range(h):
        P[(i + k) % h, i] = 1
    return P


def reflection(h):
    J = np.zeros((h, h))
    for i in range(h):
        J[i, h - 1 - i] = (-1) ** i
    return J


ok6 = True; tab6 = []
for n in (1, 2, 3, 4):
    rep = algebra(f'a{n}(1)'); sub = rep[4]; h = n + 1
    res = []
    for phi0 in (sub @ rng.normal(size=sub.shape[1]), np.zeros(h)):
        for sig in ('T', 'S'):
            res.append(dirichlet_K(rep, phi0, sig, 0.83))
            res += [dirichlet_K(rep, phi0, sig, 0.83, P=cyclic(h, k)) for k in range(1, h)]
            res.append(dirichlet_K(rep, phi0, sig, 0.83, P=reflection(h)))
        res += [dirichlet_K(rep, phi0, sig, 0.83) for sig in ('dag', 'nconj')]
    tab6.append(f'a{n}: {sorted(set(res))}')
    ok6 &= all(r > 0 for r in res) if n == 1 else all(r == 0 for r in res)
report(f'real coupling, Dirichlet on all components (generic phi_0 and phi_0 = 0): an invertible K exists for a1 and for no linear sigma '
       f'(transposition, antipode, composed with all Z_h rotations and the reflection; dagger, -conj on real data) for a2, a3, a4  '
       f'[nullities {"; ".join(tab6)}]', ok6)


def adjoint(n):
    """adjoint of sl(h) in a Hilbert-Schmidt orthonormal basis, so that rho(X^T) = rho(X)^T."""
    roots, Es, Hs, nj, sub = algebra(f'a{n}(1)'); h = n + 1
    basis = [eij(h, i, k) for i in range(h) for k in range(h) if i != k]
    for k in range(1, h):
        d = np.zeros(h); d[:k] = 1; d[k] = -k; basis.append(np.diag(d / np.sqrt(k * (k + 1))))
    Bm = np.array([x.ravel() for x in basis])
    rho = lambda X: Bm @ np.array([(X @ y - y @ X).ravel() for y in basis]).T
    return roots, [rho(E) for E in Es], [rho(H) for H in Hs], nj, sub


def wedge2(n):
    roots, Es, Hs, nj, sub = algebra(f'a{n}(1)'); h = n + 1
    P = []
    for i, k in itertools.combinations(range(h), 2):
        v = np.zeros((h, h)); v[i, k] = 1 / np.sqrt(2); v[k, i] = -1 / np.sqrt(2); P.append(v.ravel())
    P = np.array(P)
    rho = lambda X: P @ (np.kron(X, np.eye(h)) + np.kron(np.eye(h), X)) @ P.T
    return roots, [rho(E) for E in Es], [rho(H) for H in Hs], nj, sub


def v_vbar(n):
    roots, Es, Hs, nj, sub = algebra(f'a{n}(1)'); Z = np.zeros((n + 1, n + 1))
    rho = lambda X: np.block([[X, Z], [Z, -X.T]])
    return roots, [rho(E) for E in Es], [rho(H) for H in Hs], nj, sub


def rep_ok(rep):
    roots, Es, Hs, nj, sub = rep
    ok = all(np.allclose((E @ E.T - E.T @ E) @ E - E @ (E @ E.T - E.T @ E), 2 * E) for E in Es)
    for w in sub.T:
        pH = sum(w[k] * Hs[k] for k in range(len(Hs)))
        ok &= all(np.allclose(pH @ E - E @ pH, (a @ w) * E) for a, E in zip(roots, Es))
    return ok


reps = {'a2 adjoint': adjoint(2), 'a3 adjoint': adjoint(3), 'a3 Lambda^2': wedge2(3), 'a2 V+Vbar': v_vbar(2), 'd4 vector': algebra('d4(1)'),
        'a1 adjoint': adjoint(1)}
ok7 = True; tab7 = []
for name, rep in reps.items():
    phi0 = rep[4] @ rng.normal(size=rep[4].shape[1])
    r = (dirichlet_K(rep, phi0, 'T', 0.83), dirichlet_K(rep, phi0, 'S', 0.83))
    tab7.append(f'{name} {r}')
    ok7 &= rep_ok(rep) and (all(x > 0 for x in r) if name == 'a1 adjoint' else all(x == 0 for x in r))
report(f'other representations, Dirichlet on all components: no K (transposition, antipode) in adjoint(a2, a3), Lambda^2(a3), V+Vbar(a2), '
       f'vector of d4^(1); K exists for a1 adjoint  [{"; ".join(tab7)}]', ok7)

ok8 = True; tab8 = []
for n in (1, 2, 3):
    rep = algebra(f'a{n}(1)'); sub = rep[4]
    phi0 = sub @ (rng.normal(size=sub.shape[1]) + 1j * rng.normal(size=sub.shape[1]))
    r = {s: dirichlet_K(rep, phi0, s, 0.83, coupling='imag', cplx=True) for s in ('T', 'S', 'dag', 'nconj')}
    tab8.append(f'a{n} {r}')
    ok8 &= (r['T'] > 0 and r['S'] > 0 and r['dag'] == 0 and r['nconj'] == 0) if n == 1 else all(x == 0 for x in r.values())
report(f'imaginary coupling, complex Dirichlet phi(0) = phi_0 on all components: no K for any sigma (a2, a3); a1: linear sigma only  '
       f'[{"; ".join(tab8)}]', ok8)

ok9 = True; tab9 = []
for name, rep, h, sigs in [('a2 defining', algebra('a2(1)'), 3, ('S',)), ('a3 defining', algebra('a3(1)'), 4, ('S',)),
                           ('a2 adjoint', reps['a2 adjoint'], 3, ('T', 'S')), ('a3 Lambda^2', reps['a3 Lambda^2'], 4, ('T', 'S'))]:
    phi0 = rep[4] @ rng.normal(size=rep[4].shape[1])
    at_root = [dirichlet_K(rep, phi0, s, np.exp(1j * np.pi / h)) for s in sigs]
    generic = [dirichlet_K(rep, phi0, s, np.exp(1j * np.pi / (2 * h))) for s in sigs]
    tab9.append(f'{name} {sigs}: {at_root} vs {generic}')
    ok9 &= all(x > 0 for x in at_root) and all(x == 0 for x in generic)
report(f'the obstruction is the cycle condition: K exists at lam = e^(i pi/h) (lam^(2h) = 1), not at lam = e^(i pi/2h)  [{"; ".join(tab9)}]', ok9)


# ---------------------------------------------------------------- 10. spin +-2 charges at a Dirichlet wall (no Sklyanin)
def spin2_current(n, bt):
    roots, Es, Hs, nj, sub = algebra(f'a{n}(1)')
    r = [sub.T @ a for a in roots]
    s3 = list(itertools.combinations_with_replacement(range(n), 3)); a2 = list(itertools.combinations(range(n), 2))
    def unpack(vec):
        A = np.zeros((n, n, n), complex); B = np.zeros((n, n), complex)
        for v, idx in zip(vec[:len(s3)], s3):
            for P in set(itertools.permutations(idx)):
                A[P] = v
        for v, (i, j) in zip(vec[len(s3):], a2):
            B[i, j], B[j, i] = v, -v
        return A, B
    cols = []
    for k in range(len(s3) + len(a2)):
        A, B = unpack(np.eye(len(s3) + len(a2))[k])
        cols.append(np.concatenate([((3 / bt) * np.einsum('xyc,c->xy', A, rj) - np.outer(B @ rj, rj) - np.outer(rj, B @ rj))[np.triu_indices(n)]
                                    for rj in r]))
    U, sv, Vh = np.linalg.svd(np.array(cols).T)
    nul = int((sv < 1e-10 * sv[0]).sum()) + max(0, len(s3) + len(a2) - len(sv))
    A, B = unpack(Vh[-1].conj())
    return nul, A, B, r, sub


def gradV(phi, r, bt, mm):
    return (mm**2 / bt) * sum(rj * np.exp(bt * (rj @ phi)) for rj in r)


ok10 = True; tab10 = []
for n in (1, 2, 3, 4):
    nul, A, B, r, sub = spin2_current(n, 0.83)
    if n == 1:
        ok10 &= nul == 0; tab10.append('a1: none'); continue
    Fp = max(abs(2 * np.einsum('abc,a,b,c->', A, u, u, u) - 4 * u @ B @ gradV(phi0, r, 0.83, 1.27))
             for u, phi0 in ((rng.normal(size=n), rng.normal(size=n)) for _ in range(10)))
    ok10 &= nul == 1 and Fp > 1.0 and np.linalg.norm(B + B.T) < 1e-12 and np.linalg.norm(B) > 0.1 and np.linalg.norm(A) > 0.1
    tab10.append(f'a{n}: dim 1, max |energy flux| {Fp:.1f}, |A| {np.linalg.norm(A):.2f}, |B| {np.linalg.norm(B):.2f}')
# the current is conserved on an exact solution: book solitons at imaginary coupling solve the real-coupling equations with beta -> i beta
nul, A, B, r, sub = spin2_current(2, 1j * bS)
sol2 = Sol(2, [1, 2], [0.5, -0.3], [0.2 + 0.4j, -0.6 + 0.9j]); hh = 1e-3; kk = 2e-3


def phi_at(x, t):
    ph, _, _ = sol2.phi(np.concatenate([np.linspace(-30, x - 0.01, 3000), [x]]), t)
    return sub.T @ ph[:, -1]


def TTh(x, t):
    f = lambda a, b: phi_at(x + a, t + b)
    p = sum(c * f(k * hh, k * hh) for k, c in [(-2, 1 / 12), (-1, -8 / 12), (1, 8 / 12), (2, -1 / 12)]) / hh
    q = sum(c * f(k * hh, k * hh) for k, c in [(-2, -1 / 12), (-1, 16 / 12), (0, -30 / 12), (1, 16 / 12), (2, -1 / 12)]) / hh**2
    return np.einsum('abc,a,b,c->', A, p, p, p) + q @ B @ p, -p @ B @ gradV(f(0, 0), r, 1j * bS, mS)


st = [(-2, 1 / 12), (-1, -8 / 12), (1, 8 / 12), (2, -1 / 12)]
cl = 0; sc = 0
for x0 in (-1.5, 0.0, 1.5):
    dmT = sum(c * TTh(x0 - j * kk, 0.2 + j * kk)[0] for j, c in st) / kk
    dpTh = sum(c * TTh(x0 + j * kk, 0.2 + j * kk)[1] for j, c in st) / kk
    cl = max(cl, abs(dmT - dpTh)); sc = max(sc, abs(dmT))
ok10 &= cl < 1e-6 * sc
report(f'spin +-2 charges (no Sklyanin): unique current for a2, a3, a4, conserved on an exact solution [rel. {cl / sc:.1e}]; at a Dirichlet wall the '
       f'energy-like flux has no u-dot term and is nonzero, the momentum-like flux 4B(u-dot, u) is not a time derivative (B antisymmetric)  '
       f'[{"; ".join(tab10)}]', ok10)


# ---------------------------------------------------------------- 11. Delius's own equations (hep-th/9809140 (5.1), (5.2)) with his (1.4)
# Delius's field phi_D = i beta phi; his a_t (3.1): (1/2) d_x phi_D.H + sum_i sqrt(eta_i alpha_i^2/8) (lam E_i - E_-i/lam) e^{alpha_i.phi_D/2};
# [A, B]_+- = AB +- BA.  (1.4) phi_D^*(-x) = phi_D(x) + 2 pi i lambda gives Im phi_D(0) = -pi lambda; (1.5) requires Im phi_D(0)/2pi in Lambda.
ok11 = True; w51 = 0; best11 = np.inf
for n in (1, 2, 3):
    roots, Es, Hs, nj, sub = algebra(f'a{n}(1)'); h = n + 1
    for S in classes(n):
        lam_cw = eps_S(S, h) + 2 * eps_S(tuple(rng.choice(h, size=rng.integers(0, h), replace=False)), h)
        for K in (np.diag(np.exp(1j * np.pi * lam_cw)), np.eye(h)):
            worst = 0
            for lam in (0.6, 1.7):
                phiD = sub @ rng.normal(size=sub.shape[1]) - 1j * np.pi * lam_cw          # Re free, Im = -pi lambda
                dRe, dIm = np.zeros(h), sub @ rng.normal(size=sub.shape[1])            # Neumann on Re phi_D
                X = [0.5 * (lam * E - E.T / lam) for E in Es]
                ee = [np.exp(a @ phiD / 2) for a in roots]
                H1 = 0.5 * sum(dRe[k] * Hs[k] for k in range(h)); H2 = 0.5 * sum(dIm[k] * Hs[k] for k in range(h))
                Y1 = sum(x * e.real for x, e in zip(X, ee)); Y2 = sum(x * e.imag for x, e in zip(X, ee))
                e51 = (K @ H1 + H1 @ K) + (K @ Y1 - Y1 @ K)
                e52 = (K @ H2 - H2 @ K) + (K @ Y2 + Y2 @ K)
                worst = max(worst, np.linalg.norm(e51) + np.linalg.norm(e52))
            if K is not None and np.allclose(K, np.eye(h)):
                if S:
                    best11 = min(best11, worst)
                else:
                    ok11 &= worst < 1e-12
            else:
                w51 = max(w51, worst)
report(f'Delius hep-th/9809140 (5.1)-(5.2), with his image condition (1.4) for every coweight lambda: solved by K_+ = e^{{i pi lambda.H}}; '
       f'K_+ = 1 (his proof) only for lambda in 2 Lambda, the case in which (1.4) implies his (1.5)  [{w51:.1e}; K = 1 off 2 Lambda >= {best11:.2f}]',
       ok11 and w51 < 1e-12 and best11 > 0.1)
