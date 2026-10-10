"""Ch.9 checks: e8 masses; Dorey rule + area rule in charge-eigenvector basis for degenerate cases;
Coxeter-plane projection statements.  Uses foundations_code/hirota/lie.py (copied)."""
import itertools, numpy as np
from lie import Algebra

def report(n, ok): print(('PASS ' if ok else 'FAIL ') + n)

# ---- e8 masses
A = Algebra('e', 8)
m = np.array(sorted(A.mass.values())); r = m / m[0]
c = np.cos; p = np.pi
m2 = 2*c(p/5)
book = [1, 2*c(p/5), 2*c(p/30)]
known = [1, m2, 2*c(p/30), 2*m2*c(7*p/30), 2*m2*c(2*p/15), 2*m2*c(p/30), 4*m2*c(p/5)*c(7*p/30), 4*m2*c(p/5)*c(2*p/15)]
print('e8 mass ratios      ', np.round(r, 6))
print('closed forms (Zam.) ', np.round(known, 6))
report('e8: book table ratios 1 : 2cos(pi/5) : 2cos(pi/30) are the three lightest', np.allclose(r[:3], book, atol=1e-10))
report('e8: all eight ratios match the standard closed forms', np.allclose(r, known, atol=1e-10))
print('   sum m^2 =', round(float((m**2).sum()), 10), ' 2h =', 2*A.h)

def bicolour_coxeter(Al):
    r = Al.r; al = Al.alpha[1:]; C = Al.Cfin
    col = {0: 1}; st = [0]
    while st:
        i = st.pop()
        for j in range(r):
            if j != i and C[i, j] < 0 and j not in col:
                col[j] = -col[i]; st.append(j)
    refl = lambda v: np.eye(r) - np.outer(v, v)
    sW = np.eye(r); sB = np.eye(r)
    for i in range(r):
        if col[i] == 1: sW = refl(al[i]) @ sW
        else: sB = refl(al[i]) @ sB
    w = sB @ sW
    lam = np.linalg.inv(C) @ al            # fundamental weights (rows), node a -> row a-1
    gam = [(np.eye(r) - np.linalg.inv(w)) @ lam[a] for a in range(r)]
    orbits = {a+1: [np.linalg.matrix_power(w, k) @ gam[a] for k in range(Al.h)] for a in range(r)}
    return w, orbits

def roots_of(Al):
    al = list(Al.alpha[1:]); R = {tuple(np.round(a, 9)) for a in al}; fr = list(al)
    while fr:
        new = []
        for v in fr:
            for a in al:
                w_ = v - (v @ a)*a; t = tuple(np.round(w_, 9))
                if t not in R: R.add(t); new.append(w_)
        fr = new
    return [np.array(t) for t in R]

def check(typ, n, swap_conj=False):
    Al = Algebra(typ, n); h = Al.h; nodes = sorted(Al.u)
    u = dict(Al.u)
    if swap_conj:     # relabel: species a <- eigenvector of conj(a) (global charge conjugation)
        u = {a: Al.u[Al.conj[a]] for a in nodes}
    w, orb = bicolour_coxeter(Al)
    roots = roots_of(Al)
    allorb = {tuple(np.round(v, 6)) for a in nodes for v in orb[a]}
    part = len(allorb) == len(roots) == n*h
    cub = lambda a, b, c: sum(Al.n[j]*(Al.alpha[j] @ u[a])*(Al.alpha[j] @ u[b])*(Al.alpha[j] @ u[c]) for j in range(n+1))
    ok_rule = ok_area = True; consts = []; nz = 0; bad = []
    for a, b, c in itertools.combinations_with_replacement(nodes, 3):
        val = cub(a, b, c)
        ma, mb, mc = Al.mass[a], Al.mass[b], Al.mass[c]
        s = (ma+mb+mc)/2; ar2 = s*(s-ma)*(s-mb)*(s-mc)
        rule = any(np.allclose(x+y+z, 0, atol=1e-8) for x in orb[a] for y in orb[b] for z in orb[c])
        if (abs(val) > 1e-8) != rule: ok_rule = False; bad.append((a, b, c, abs(val), rule))
        if abs(val) > 1e-8:
            nz += 1
            if ar2 <= 1e-12: ok_area = False
            else: consts.append(abs(val)/np.sqrt(ar2))
    ok_area &= bool(consts) and np.ptp(consts) < 1e-8 and abs(consts[0] - 4/np.sqrt(h)) < 1e-8
    # Coxeter plane: real 2-plane of eigenvalue exp(2 pi i/h)
    ev, V = np.linalg.eig(w)
    k = np.argmin(abs(ev - np.exp(2j*np.pi/h)))
    P = np.array([V[:, k].real, V[:, k].imag]); P = np.linalg.qr(P.T)[0].T   # orthonormal basis rows
    proj = {a: [P @ v for v in orb[a]] for a in nodes}
    lens = {a: np.linalg.norm(proj[a][0]) for a in nodes}
    const_len = all(np.allclose([np.linalg.norm(x) for x in proj[a]], lens[a]) for a in nodes)
    ratio = np.array([lens[a]/Al.mass[a] for a in nodes]); ok_len = const_len and np.ptp(ratio) < 1e-8*ratio.max()
    # angles between all projected roots: multiples of pi/h ?
    angs = []
    for a in nodes:
        for b in nodes:
            x, y = proj[a][0], proj[b][0]
            ang = np.arctan2(x[0]*y[1]-x[1]*y[0], x @ y)
            angs.append(ang/(np.pi/h))
    ok_ang = np.allclose(angs, np.round(angs), atol=1e-8)
    tag = f'{typ}{n}' + (' (conj-swapped labels)' if swap_conj else '')
    print(f'   {tag}: h={h}, nonzero C_abc: {nz}, |C|/area in [{min(consts):.8f},{max(consts):.8f}], 4/sqrt(h)={4/np.sqrt(h):.8f}; mismatches {bad[:4]}')
    return part, ok_rule, ok_area, ok_len, ok_ang

for typ, n in [('a', 4), ('a', 5), ('d', 4), ('d', 5), ('d', 6), ('e', 6), ('e', 7), ('e', 8)]:
    part, ok_rule, ok_area, ok_len, ok_ang = check(typ, n)
    report(f'{typ}{n}: orbits partition roots; Dorey rule iff; |C|=(4/sqrt h)*area; |proj root| prop. to m_a; angles in (pi/h)Z',
           part and ok_rule and ok_area and ok_len and ok_ang)
    if not (part and ok_rule and ok_area and ok_len and ok_ang): print('     ', part, ok_rule, ok_area, ok_len, ok_ang)
for typ, n in [('a', 4), ('d', 5), ('e', 6)]:
    part, ok_rule, ok_area, ok_len, ok_ang = check(typ, n, swap_conj=True)
    report(f'{typ}{n} with globally conjugated labels: Dorey rule still iff', ok_rule)
