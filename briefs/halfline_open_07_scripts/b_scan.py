"""Problem 7(b): unimodularity of the boundary Bethe-Yang transfer-matrix eigenvalues, a_2^(1), by Z_3 sector,
and the periodic (bulk) comparison.  Output: max_j max ||Lambda|-1| over random real rapidities."""
import numpy as np, sys
from b_transfer import Tdr, Zshift, kron_all, embed2
from p07_common import Rc
rng = np.random.default_rng(5)
n = 2; d = 3
zetas = np.exp(2j*np.pi*np.arange(3)/3)

def joint(Ts, Z):
    G = sum(rng.normal()*T for T in Ts) + 0.37*Z
    w, P = np.linalg.eig(G); Pi = np.linalg.inv(P)
    lam = np.array([np.diag(Pi @ T @ P) for T in Ts])        # (Nparticles, dim)
    zz = np.diag(Pi @ Z @ P)
    return lam, zz

def bulk_T(types, th, om):
    """periodic Bethe-Yang transfer matrix of particle 0 (product of R-checks, then cyclic shift back)."""
    Np = len(types); y = np.exp(om*np.asarray(th)); M = np.eye(d**Np, dtype=complex)
    order = list(range(Np))
    for p in range(Np - 1):
        a, b = order[p], order[p + 1]
        M = embed2(Rc(n, types[a], types[b], y[a]/y[b], om), p, Np, d) @ M
        order[p], order[p + 1] = b, a
    # particle 0 now at the right end; cyclic shift moves it back to position 0
    S = np.zeros((d**Np, d**Np))
    for idx in np.ndindex(*([d]*Np)):
        new = (idx[-1],) + idx[:-1]
        S[np.ravel_multi_index(new, [d]*Np), np.ravel_multi_index(idx, [d]*Np)] = 1
    return S @ M

def scan(cfg, Np, om, nsamp=12, types=None):
    types = types or ['V']*Np
    Z = kron_all([Zshift(n, -1 if cfg == 'MM' else 1)]*Np)
    worst = {}
    for _ in range(nsamp):
        th = np.sort(rng.uniform(0.05, 2.5, Np))[::-1]
        Ts = [Tdr(j, types, th, cfg[0], cfg[1], n, om) for j in range(Np)]
        lam, zz = joint(Ts, Z)
        dev = np.abs(np.abs(lam) - 1).max(axis=0)
        for z, dv in zip(zz, dev):
            key = int(np.round(np.angle(z)/(np.pi/3))) % 6 if cfg != 'NM' else 0
            worst[key] = max(worst.get(key, 0), dv)
    return worst

def lab(k):
    return {0: 'z=1', 2: 'z=w', 4: 'z=w^2', 1: 'z=e^{i pi/3}', 3: 'z=-1', 5: 'z=e^{-i pi/3}'}[k]

oms = [0.37, 0.83, 1.43, 2.6]
for cfg in ('UU', 'UN', 'NN', 'MM', 'NM'):
    for Np in (1, 2, 3):
        row = []
        for om in oms:
            w = scan(cfg, Np, om, 8 if Np == 3 else 16)
            row.append(' '.join(f'{lab(k) if cfg != "NM" else "all"}:{v:.2g}' for k, v in sorted(w.items())))
        print(f'{cfg} N={Np}: ' + ' | '.join(f'om={om}: {r}' for om, r in zip(oms, row)))
        sys.stdout.flush()
# mixed soliton/antisoliton content, uniform walls
for types in (['V', 'W'], ['V', 'V', 'W']):
    row = []
    for om in oms:
        w = scan('UU', len(types), om, 10, types)
        row.append(' '.join(f'{lab(k)}:{v:.2g}' for k, v in sorted(w.items())))
    print(f'UU types={"".join(types)}: ' + ' | '.join(f'om={om}: {r}' for om, r in zip(oms, row)))
# bulk periodic comparison: T_1 eigenvalue moduli
for types in (['V', 'V'], ['V', 'W'], ['V', 'V', 'V'], ['V', 'V', 'W']):
    row = []
    for om in oms:
        dv = 0
        for _ in range(10):
            th = rng.uniform(-2, 2, len(types))
            dv = max(dv, np.abs(np.abs(np.linalg.eigvals(bulk_T(types, th, om))) - 1).max())
        row.append(f'{dv:.2g}')
    print(f'bulk periodic types={"".join(types)}: max||s|-1| ' + ' | '.join(f'om={om}: {r}' for om, r in zip(oms, row)))
