"""Checks of chapter 5 (05-root-systems.qmd): tbl-simple-lie, Coxeter element (Steinberg), PF vector, Casimir."""
import numpy as np, itertools
from lie import *
rng = np.random.default_rng(1)

table = {  # (d, h, hv, exponents) as printed in tbl-simple-lie
 'A': lambda n: (n*(n+2), n+1, n+1, list(range(1, n+1))),
 'B': lambda n: (n*(2*n+1), 2*n, 2*n-1, list(range(1, 2*n, 2))),
 'C': lambda n: (n*(2*n+1), 2*n, n+1, list(range(1, 2*n, 2))),
 'D': lambda n: (n*(2*n-1), 2*n-2, 2*n-2, sorted(list(range(1, 2*n-2, 2)) + [n-1])),
 'E': lambda n: {6: (78, 12, 12, [1,4,5,7,8,11]), 7: (133, 18, 18, [1,5,7,9,11,13,17]),
                 8: (248, 30, 30, [1,7,11,13,17,19,23,29])}[n],
 'F': lambda n: (52, 12, 9, [1,5,7,11]),
 'G': lambda n: (14, 6, 4, [1,5]),
}
centre = {'A': lambda n: n+1, 'B': lambda n: 2, 'C': lambda n: 2, 'D': lambda n: 4,
          'E': lambda n: {6: 3, 7: 2, 8: 1}[n], 'F': lambda n: 1, 'G': lambda n: 1}
cases = [('A', n) for n in range(2, 7)] + [('B', n) for n in range(2, 6)] + [('C', n) for n in range(3, 6)] + \
        [('D', n) for n in range(4, 8)] + [('E', 6), ('E', 7), ('E', 8), ('F', 4), ('G', 2)]
allok = True
for typ, n in cases:
    R = RootSystem(typ, n)
    d, h, hv, ex = table[typ](n)
    nroots = len(R.roots)
    ok = (nroots + n == d) and R.h == h and abs(R.hv - hv) < 1e-9 and nroots == n*h
    ok &= sum(ex) == len(R.pos)
    # Coxeter elements in random orders
    perms = list(itertools.permutations(range(n))) if n <= 5 else [list(rng.permutation(n)) for _ in range(30)]
    cp = None
    for p in perms:
        w = R.coxeter(p)
        ok &= order_of(w) == h
        poly = np.round(np.poly(w), 8)
        if cp is None: cp = poly
        ok &= np.allclose(poly, cp)
    ok &= np.allclose(exponents_of(R.coxeter(), h), ex)
    # bicoloured Coxeter element and Steinberg plane
    W_, B_ = R.bicolour()
    sW = np.eye(n); sB = np.eye(n)
    for i in W_: sW = sW @ refl(R.al[i])
    for i in B_: sB = sB @ refl(R.al[i])
    w = sB @ sW
    ev, V = np.linalg.eig(w)
    k = np.argmin(abs(ev - np.exp(2j*np.pi/h)))
    z = V[:, k]; P = np.linalg.qr(np.column_stack([z.real, z.imag]))[0]
    Wp = P.T @ w @ P
    rot = np.allclose(w @ P, P @ Wp) and np.allclose(Wp @ Wp.T, np.eye(2)) and \
          abs(abs(np.arctan2(Wp[1, 0], Wp[0, 0])) - 2*np.pi/h) < 1e-9
    # orbits of w on roots
    Rt = np.array(R.roots)
    perm = [int(np.argmin(np.linalg.norm(Rt - w @ r, axis=1))) for r in Rt]
    assert all(np.allclose(Rt[perm[i]], w @ Rt[i]) for i in range(len(Rt)))
    seen = set(); orbits = []
    for i in range(len(Rt)):
        if i in seen: continue
        o = [i]; j = perm[i]
        while j != i: o.append(j); j = perm[j]
        seen |= set(o); orbits.append([Rt[k] for k in o])
    sizes = sorted(len(o) for o in orbits)
    polygons = True; radii = []
    for o in orbits:
        pr = np.array([P.T @ v for v in o]); rad = np.linalg.norm(pr, axis=1)
        ang = np.diff(np.unwrap(np.arctan2(pr[:, 1], pr[:, 0])))
        polygons &= np.ptp(rad) < 1e-9 and rad[0] > 1e-6 and np.allclose(abs(ang), 2*np.pi/h)
        radii.append(rad[0])
    ok_st = rot and sizes == [h]*n and polygons
    # Casimir of adjoint: theta.(theta+2 rho) = 2 hv ; rho = sum lambda_i = half sum positive roots
    rho = R.lam.sum(0)
    ok_cas = np.allclose(rho, 0.5*sum(R.pos)) and abs(R.theta @ (R.theta + 2*rho) - 2*R.hv) < 1e-9
    ok_det = abs(np.linalg.det(R.A) - centre[typ](n)) < 1e-9
    line = f'{typ}{n}: |Phi|={nroots} d={d} h={R.h} hv={R.hv:g} marks={list(R.marks)} comarks={[float(x) for x in R.comarks]} ' \
           f'table={ok} Steinberg(plane rotation 2pi/h, {len(orbits)} orbits of size h, regular h-gons)={ok_st} ' \
           f'Casimir C(theta)=2hv & rho={ok_cas} detA=|centre|={ok_det}'
    if typ in 'ADE':
        ce = np.sort(np.linalg.eigvalsh(R.A))
        ok_ev = np.allclose(ce, np.sort(4*np.sin(np.pi*np.array(ex)/(2*h))**2))
        pf = np.abs(np.linalg.eigh(R.A)[1][:, 0]); 
        line += f' Cartan-eigs={ok_ev} PF/min={np.round(np.sort(pf)/pf.min(), 6)} radii/min={np.round(np.sort(radii)/min(radii), 6)}'
        ok_ev &= np.allclose(np.sort(pf)/pf.min(), np.sort(radii)/min(radii))
        allok &= ok_ev
    allok &= ok and ok_st and ok_cas and ok_det
    print(line)
print('ALL', 'PASS' if allok else 'FAIL')
