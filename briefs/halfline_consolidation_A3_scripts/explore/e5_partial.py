"""(a) Partial Dirichlet from involutive automorphisms of the affine a_n diagram (real coupling): image phi(x) = M phi(-x)
       => P_- phi(0) = 0 (Dirichlet at 0 on M-odd part), P_+ d_x phi(0) = 0 (Neumann on M-even part). Sklyanin K?
   (b) transposition-type sigma in self-conjugate reps at lam^(2h) = 1 (mechanism of the obstruction)
   (c) book's classical reflection equation (RE_plus, r_toda, real spectral parameters) for constant diagonal K_+ = e^{i pi nu.H}
"""
import itertools
import numpy as np
from booklib import load_classical
import types
_e3 = open('e3_claim5.py').read()
E3 = types.SimpleNamespace()
_ns = {'__name__': 'e3defs'}
exec(compile(_e3[:_e3.index('lam_gen = 0.83')], 'e3_claim5.py[defs]', 'exec'), _ns)
E3.__dict__.update(_ns)

C = load_classical()
algebra, a_t = C['algebra'], C['a_t']
rng = np.random.default_rng(17)
SIG = E3.SIG


def nullity_partial(rep, M, sig, lam, swap=False, c_odd=0.0, nsamp=8):
    roots, Es, Hs, nj, sub = rep
    h = len(roots[0]); N = Es[0].shape[0]
    Pp = (np.eye(h) + M) / 2; Pm = (np.eye(h) - M) / 2
    if swap:
        Pp, Pm = Pm, Pp
    odd_dir = sub @ rng.normal(size=sub.shape[1]); odd_dir = Pm @ odd_dir
    rows = []
    for _ in range(nsamp):
        phi0 = Pp @ (sub @ rng.normal(size=sub.shape[1])) + c_odd * odd_dir    # Dirichlet: odd part fixed (0 or c)
        dx = Pm @ (sub @ rng.normal(size=sub.shape[1]))                           # Neumann: even part of d_x phi = 0
        A = a_t(lam, phi0, dx, roots, Es, Hs, nj, 'real')
        B = SIG[sig](a_t(1 / lam, phi0, dx, roots, Es, Hs, nj, 'real'))
        cols = []
        for p in range(N):
            for q in range(N):
                K = np.zeros((N, N), complex); K[p, q] = 1
                cols.append((K @ A + B @ K).ravel())
        rows.append(np.array(cols).T)
    Mt = np.vstack(rows)
    U, sv, Vh = np.linalg.svd(Mt)
    nul = int((sv < 1e-10 * sv[0]).sum())
    if nul == 0:
        return 0, None
    K = Vh[-1].conj().reshape(N, N)
    inv = np.linalg.svd(K, compute_uv=False); ok = inv[-1] > 1e-7 * inv[0]
    return (nul if ok else -nul), K


def involutions(h):
    out = []
    if h % 2 == 0:
        M = np.zeros((h, h))
        for i in range(h):
            M[(i + h // 2) % h, i] = 1
        out.append(('rotation by h/2', M))
    for c in range(h):
        M = np.zeros((h, h))
        for i in range(h):
            M[(c + 1 - i) % h, i] = -1
        out.append((f'reflection j -> {c}-j', M))
    return out


print('=== (a) partial Dirichlet from involutive diagram automorphisms (real coupling, defining rep, lam = 0.83)')
for n in (1, 2, 3):
    rep = algebra(f'a{n}(1)'); roots, Es, Hs, nj, sub = rep; h = n + 1
    for name, M in involutions(h):
        assert all(any(np.allclose(M @ a, b) for b in roots) for a in roots) and np.allclose(M @ M, np.eye(h))
        Pm = (np.eye(h) - M) / 2
        dD = int(round(np.trace(sub.T @ Pm @ sub)))
        res = {s: nullity_partial(rep, M, s, 0.83)[0] for s in ('T', 'S')}
        sw = {s: nullity_partial(rep, M, s, 0.83, swap=True)[0] for s in ('T', 'S')}
        cnz = {s: nullity_partial(rep, M, s, 0.83, c_odd=0.4)[0] for s in ('T', 'S')} if dD > 0 else '-'
        fixed = [j for j in range(h) if np.allclose(M @ roots[j], roots[j])]
        print(f'a_{n} {name:22s} (fixed nodes {fixed}; Dirichlet on {dD} of {n} components): K for T/S {res}; '
              f'swapped D/N {sw}; odd value c != 0: {cnz}')

print('=== (b) transposition-type sigma, self-conjugate reps, lam^(2h) = 1 vs generic')
for name, rep, h in [('a2 adjoint', E3.adjoint(2), 3), ('a3 Lambda^2', E3.wedge2(3), 4), ('a2 V+Vbar', E3.v_plus_vbar(2), 3)]:
    roots, Es, Hs, nj, sub = rep
    phi0 = sub @ rng.normal(size=sub.shape[1])
    out = {f'lam={np.round(l, 3)}': (E3.nullity(rep, phi0, 'T', l)[0], E3.nullity(rep, phi0, 'S', l)[0])
           for l in (np.exp(1j * np.pi / h), np.exp(1j * np.pi / (2 * h)), 0.83)}
    print(f'{name}: (T, S) nullities {out}')

print('=== (c) classical reflection equation (book RE_plus with r_toda at real lam) for K_+ = e^{i pi nu.H}')
RE_plus, r_toda = C['RE_plus'], C['r_toda']
for N in (2, 3, 4):
    rf = lambda z, N=N: r_toda(N, z)
    worst = 0
    for S in itertools.chain.from_iterable(itertools.combinations(range(N), k) for k in range(N)):
        nu = -len(S) / N * np.ones(N); nu[list(S)] += 1
        D = np.diag(np.exp(1j * np.pi * nu))
        for l1, l2 in [(0.7, 1.3), (1.9, 0.45), (2.2, -0.6)]:
            worst = max(worst, np.linalg.norm(RE_plus(lambda l: D, rf, l1, l2, N)) / np.linalg.norm(np.kron(D, D)))
    Kr = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    bad = np.linalg.norm(RE_plus(lambda l: Kr, rf, 0.7, 1.3, N)) / np.linalg.norm(np.kron(Kr, Kr))
    print(f'N={N}: all {2**(N-1)} classes: max relative residual {worst:.1e}; random constant K: {bad:.1e}')
