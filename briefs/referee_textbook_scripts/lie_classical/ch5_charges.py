"""sec-perron-frobenius: 'charge eigenvalues chi_a^(s) are eigenvectors of the Cartan matrix'.
Test: fusings (a,b,c) from Dorey's rule (all incoming); triangle directions from the PF masses alone;
for each spin s=1..h-1 solve the conservation laws sum_{x in abc} q_x e^{i s phi_x}=0 for real q,
and compare with the eigenspace of the Cartan matrix with eigenvalue 4 sin^2(pi s/2h)."""
import numpy as np, itertools
from lie import *

def fusings(R):
    n, h = R.n, R.h
    W_, B_ = R.bicolour()
    sW = np.eye(n); sB = np.eye(n)
    for i in W_: sW = sW @ refl(R.al[i])
    for i in B_: sB = sB @ refl(R.al[i])
    w = sB @ sW; wi = np.linalg.inv(w)
    gam = [(np.eye(n) - wi) @ R.lam[a] for a in range(n)]
    orb = [[np.linalg.matrix_power(w, p) @ gam[a] for p in range(h)] for a in range(n)]
    out = []
    for a, b, c in itertools.combinations_with_replacement(range(n), 3):
        if any(np.allclose(x + y + z, 0) for x in orb[a] for y in orb[b] for z in orb[c]):
            out.append((a, b, c))
    return out

def directions(ma, mb, mc):
    # interior angles opposite sides a,b,c
    A = np.arccos((mb**2 + mc**2 - ma**2)/(2*mb*mc)); B = np.arccos((ma**2 + mc**2 - mb**2)/(2*ma*mc))
    Cc = np.pi - A - B
    return 0.0, np.pi - Cc, -(np.pi - B)

for typ, n in [('A', 4), ('D', 4), ('D', 5), ('D', 6), ('E', 6), ('E', 7), ('E', 8)]:
    R = RootSystem(typ, n); h = R.h
    pf = np.abs(np.linalg.eigh(R.A)[1][:, 0]); m = pf/pf.min()
    F = fusings(R)
    ex = [int(round(x)) for x in exponents_of(R.coxeter(), h)]
    res = []
    allok = True
    for s in range(1, h):
        rows = []
        for (a, b, c) in F:
            ph = directions(m[a], m[b], m[c])
            assert abs(m[a]*np.exp(1j*ph[0]) + m[b]*np.exp(1j*ph[1]) + m[c]*np.exp(1j*ph[2])) < 1e-9
            row = np.zeros(n, complex)
            for x, p in zip((a, b, c), ph): row[x] += np.exp(1j*s*p)
            rows += [row.real, row.imag]
        M = np.array(rows)
        u, sv, vt = np.linalg.svd(M); null = vt[(sv > 1e-9).sum():]          # real solutions q
        lam = 4*np.sin(np.pi*s/(2*h))**2
        E = np.linalg.svd(R.A - lam*np.eye(n))
        eig = E[2][(E[1] > 1e-9).sum():]                                     # Cartan eigenspace
        mult = ex.count(s)
        # does the Cartan eigenspace lie in the conserved space, and do they coincide?
        inside = len(eig) == 0 or np.allclose(M @ eig.T, 0, atol=1e-9)
        same = len(null) == len(eig) and inside
        allok &= (len(eig) == mult) and same
        res.append(f's={s}:dimCons={len(null)},mult={mult},{"=" if same else ("sub" if inside else "NO")}')
    print(f'{typ}{n} (h={h}, {len(F)} fusings): ' + ' '.join(res) + f'  => {"PASS" if allok else "FAIL"}')
