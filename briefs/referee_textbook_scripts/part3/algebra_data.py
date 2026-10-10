"""Data for appendix A (sec-app-algebra-data): classical masses of the untwisted affine Toda theories.

Simply-laced: eigenvalues of M^2 = m^2 sum_j n_j alpha_j alpha_j^T (m=1).
Non-simply-laced untwisted: fold the simply-laced parent by a finite-diagram automorphism that
fixes alpha_0 (c_n <- a_{2n-1}, b_n <- d_{n+1}, g_2 <- d_4, f_4 <- e_6) and restrict M^2 to the
invariant subspace (sec-affine-folding). Prints the masses and checks the closed forms of the appendix.
"""
import itertools
import numpy as np

exec(open('toda_classical.py').read().split('# ---------------------------------------------------------------- root data')[1]
     .split('ok_pf = True')[0])


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def masses(typ, n):
    C, al, roots, h, A, nj, M2 = data(typ, n)
    return h, nj, np.sort(np.sqrt(np.abs(np.linalg.eigvalsh(M2))))


def folded(typ, n, perm):
    """perm: permutation of finite nodes 0..n-1 (a Dynkin automorphism fixing alpha_0)."""
    C, al, roots, h, A, nj, M2 = data(typ, n)
    P = np.zeros((n, n))
    # linear map sigma with sigma(alpha_i) = alpha_perm(i); invariant subspace of phi: alpha_i.phi = alpha_perm(i).phi
    Ms = np.array([al[i] - al[perm[i]] for i in range(n)])
    # null space of Ms = invariant subspace
    u, s, vt = np.linalg.svd(Ms)
    rank = int((s > 1e-9).sum())
    Vinv = vt[rank:].T                       # columns: orthonormal basis of invariant subspace
    Mr = Vinv.T @ M2 @ Vinv
    return np.sort(np.sqrt(np.abs(np.linalg.eigvalsh(Mr))))


ok = True
for n in range(1, 9):
    h, nj, ms = masses('a', n) if n > 1 else (2, None, np.array([2.0]))
    if n > 1:
        ok &= np.allclose(ms, np.sort(2*np.sin(np.pi*np.arange(1, n+1)/(n+1))))
report('a_n (n<=8): m_a = 2 m sin(pi a/h)', ok)
ok = True
for n in range(4, 9):
    h, nj, ms = masses('d', n)
    exp = np.sort(np.concatenate([2*np.sqrt(2)*np.sin(np.pi*np.arange(1, n-1)/h), [np.sqrt(2)]*2]))
    ok &= np.allclose(ms, exp)
report('d_n (4<=n<=8): m_a = 2 sqrt2 m sin(pi a/h), a<=n-2, two spinors sqrt2 m', ok)
for n in (6, 7, 8):
    h, nj, ms = masses('e', n)
    print(f'   e_{n}: h={h}, labels sum={nj.sum()}, masses/m = {np.round(ms, 6)}')

# untwisted non-simply-laced by folding
ok = True
for n in range(2, 7):    # c_n <- a_{2n-1}: reflection i -> 2n-2-i of the chain (0-based)
    perm = [2*n - 2 - i for i in range(2*n - 1)]
    ms = folded('a', 2*n - 1, perm)
    ok &= np.allclose(ms, np.sort(2*np.sin(np.pi*np.arange(1, n+1)/(2*n))))
report('c_n (n<=6) by folding a_{2n-1}: m_a = 2 m sin(pi a/2n), a=1..n', ok)
ok = True
for n in range(3, 8):    # b_n <- d_{n+1}: swap the two spinor nodes (last two)
    N = n + 1; perm = list(range(N)); perm[N-2], perm[N-1] = perm[N-1], perm[N-2]
    ms = folded('d', N, perm)
    exp = np.sort(np.concatenate([2*np.sqrt(2)*np.sin(np.pi*np.arange(1, n)/(2*n)), [np.sqrt(2)]]))
    ok &= np.allclose(ms, exp)
report('b_n (3<=n<=7) by folding d_{n+1}: m_a = 2 sqrt2 m sin(pi a/2n), a<n, and sqrt2 m', ok)
# g2 <- d4 triality: nodes (0-based) chain 0-1, 1-2? our d4 cartan: chain 0-1-2? check: d_n chain 0..n-2, node n-1 attached to n-3
# for d4: chain 0-1-2, node 3 attached to node 1: centre 1, ends 0,2,3
ms = folded('d', 4, [2, 1, 3, 0])
report(f'g_2 by folding d_4 (triality): masses {np.round(ms, 6)} = sqrt2 m, sqrt6 m', np.allclose(ms, [np.sqrt(2), np.sqrt(6)]))
# f4 <- e6: e6 chain 0-1-2-3-4, node 5 attached to node 2; reflection 0<->4, 1<->3
ms = folded('e', 6, [4, 3, 2, 1, 0, 5])
print(f'   f_4 by folding e_6: masses/m = {np.round(ms, 6)}; ratios to lightest {np.round(ms/ms[0], 6)}')
r = ms/ms[0]
report('f_4: mass ratios 1, 2cos(pi/4), 2cos(pi/12), 1+sqrt3 (the classical pattern of sec-sm-f4)',
       np.allclose(r, [1, 2*np.cos(np.pi/4), 2*np.cos(np.pi/12), 1 + np.sqrt(3)]))
