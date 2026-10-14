"""Fast part of b2_dn_classical.py: K-matrix comparison and pole positions only (no bound-state scan)."""
import numpy as np, cmath
import b2_dn_classical as b2
sm = b2.load_smatrix()
for n in (4, 5, 6):
    al, mk = b2.dn_roots(n)
    h, bl, pfm, _ = sm.setup('d', n)
    masses = np.sqrt(np.linalg.eigvalsh((al.T*mk) @ al)); mmin = masses.min()
    node_mass = {a: pfm[a]*mmin for a in range(n)}
    psi, p0 = b2.vacuum(n)
    vres, bres = b2.check_vacuum(n, psi, p0)
    def prodS(m, th):
        nodes = [a for a in range(n) if abs(node_mass[a]-m) < 1e-8]
        vals = [np.prod([sm.S(bl[(a, b)], h, 1.0, th) for b in range(n)]) for a in nodes]
        assert max(abs(v-vals[0]) for v in vals) < 1e-10
        return vals[0]
    thr = sorted(set(np.round(masses, 10)))
    Oms = [0.5*(thr[0]+thr[1]), 0.5*(thr[-2]+thr[-1]), thr[-1]*1.3, thr[-1]*2.1]
    wd = wo = 0
    for Om in Oms:
        K, opn, ms, th, U = b2.kmatrix(n, psi, p0, Om)
        for i, a in enumerate(opn):
            wd = max(wd, abs(K[i, i]-prodS(ms[a], th[a])))
            for j in range(len(opn)):
                if j != i: wo = max(wo, abs(K[i, j]))
    print(f"d_{n}: vacuum residuals {vres:.1e} {bres:.1e}; Omega {np.round(Oms,3)}: |K_aa - prod_b S_ab(th,1)| {wd:.1e}, off-diag {wo:.1e}", flush=True)
    poles = {}
    for a in range(n):
        e = {}
        for c in range(n):
            for x, k in bl[(a, c)].items():
                for y, m in ((x-1, k), (x+1, k), (x, -2*k)):
                    e[y % (2*h)] = e.get(y % (2*h), 0) + m
        for z in range(1, h):
            order = e.get(z % (2*h), 0) - e.get((-z) % (2*h), 0)
            if order > 0 and 2*z < h:
                poles.setdefault(round(2*node_mass[a]*np.sin(np.pi*z/h), 6), []).append((round(node_mass[a], 4), z, order))
    ms = sorted(set(np.round(masses, 6)))
    print(f"   poles: 2|p| -> (mass, z, order): {poles}; masses {ms}", flush=True)
