"""eq-tpg-rule beyond sl2^: vector reps of U_q(a_2^(1)), U_q(c_2^(1)), U_q(b_2^(1)), U_q(c_3^(1)), U_q(b_3^(1)).

Builds e_i,f_i,k_i (i=0..n) explicitly in the book's conventions, checks all
defining relations (incl. q-Serre with q_i = q^{alpha_i^2/2}), solves Jimbo's
equations numerically in the homogeneous gradation at random q, x, y, and
compares the eigenvalues of Rc(z), normalized to 1 on v1(x)v1, with
eq-tpg-rule: rho_mu/rho_nu = <(C(nu)-C(mu))/4>, <a> = (1 - z q^{2a})/(z - q^{2a}).
"""
import numpy as np
from qaff import check_relations, solve_jimbo, qnum

rng = np.random.default_rng(1)


def Eij(N, i, j):
    M = np.zeros((N, N), dtype=complex); M[i-1, j-1] = 1; return M


def br(a, z, q):
    return (1 - z*q**(2*a))/(z - q**(2*a))


def report(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)


def match(ev, want, tol=1e-7):
    """ev: eigenvalues; want: list of (value, multiplicity). True if multisets agree."""
    ev = list(ev); ok = True
    for val, m in want:
        hits = [i for i, v in enumerate(ev) if abs(v - val) < tol*max(1, abs(val))]
        if len(hits) != m:
            ok = False
        for i in sorted(hits, reverse=True):
            ev.pop(i)
    return ok and not ev


def run(name, rep, A, qi, q, want_fn, N):
    err, bad = check_relations(*rep, A, qi, verbose=True)
    report(f'{name}: all U_q relations incl. q-Serre (max residual {err:.1e})', not bad)
    for trial in range(3):
        x = complex(rng.normal(), rng.normal()); y = complex(rng.normal(), rng.normal()); z = x/y
        R, sv = solve_jimbo(rep, x, y)
        R = R/R[0, 0]
        ev = np.linalg.eigvals(R)
        want = want_fn(z, q)
        ok = match(ev, want)
        report(f'{name}: solution space 1-dim (last sing. values {abs(sv[0]):.1e},{abs(sv[1]):.1e}) and eigenvalues = TPG prediction, trial {trial}',
               ok and abs(sv[0]) > 1e-4 and abs(sv[1]) < 1e-6)
        if not ok and trial == 0:
            print('   eigenvalues:', np.round(np.sort_complex(ev), 6))
            print('   predicted  :', [(np.round(v, 6), m) for v, m in want])
    return R


# random generic q = t^2 (t = q^{1/2})
t = complex(0.83, 0.41); q = t*t

# ---------------- a_2^(1) vector -----------------
N = 3
e = [Eij(N, 3, 1), Eij(N, 1, 2), Eij(N, 2, 3)]
f = [m.T.copy() for m in e]
h = [np.diag([-1, 0, 1]), np.diag([1, -1, 0]), np.diag([0, 1, -1])]
k = [np.diag(q**np.diag(hh).astype(complex)) for hh in h]
A = [[2, -1, -1], [-1, 2, -1], [-1, -1, 2]]
run('a_2^(1) vector', (e, f, k), A, [q, q, q], q,
    lambda z, q: [(1, 6), (br(1, z, q), 3)], N)

# ---------------- c_2^(1) vector (C^4) -----------------
# weights e1,e2,-e2,-e1 with eps_i.eps_j = delta/2; alpha_1=eps1-eps2 short, alpha_2=2eps2 long, alpha_0=-2eps1 long
N = 4
e = [Eij(N, 4, 1), Eij(N, 1, 2) + Eij(N, 3, 4), Eij(N, 2, 3)]
f = [m.T.copy() for m in e]
h = [np.diag([-1, 0, 0, 1]), np.diag([1, -1, 1, -1]), np.diag([0, 1, -1, 0])]
d = [1, 0.5, 1]                       # alpha_i^2/2
qi = [t**(2*di) for di in d]          # q_i = q^{alpha_i^2/2}
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(3)]
A = [[2, -1, 0], [-2, 2, -2], [0, -1, 2]]
# Casimirs (long roots length^2 2): C(2l1)=2n+2=6, C(l2)=2n=4, C(0)=0  (n=2)
run('c_2^(1) vector', (e, f, k), A, qi, q,
    lambda z, q: [(1, 10), (br(0.5, z, q), 5), (br(1.5, z, q), 1)], N)

# ---------------- b_2^(1) vector (C^5) -----------------
# weights e1,e2,0,-e2,-e1 orthonormal; alpha_1=e1-e2 long, alpha_2=e2 short, alpha_0=-e1-e2
N = 5
q2 = t                                  # q_2 = q^{1/2}
e = [Eij(N, 5, 2) + Eij(N, 4, 1), Eij(N, 1, 2) + Eij(N, 4, 5), qnum(2, q2)*Eij(N, 2, 3) + Eij(N, 3, 4)]
f = [Eij(N, 2, 5) + Eij(N, 1, 4), Eij(N, 2, 1) + Eij(N, 5, 4), Eij(N, 3, 2) + qnum(2, q2)*Eij(N, 4, 3)]
h = [np.diag([-1, -1, 0, 1, 1]), np.diag([1, -1, 0, 1, -1]), np.diag([0, 2, 0, -2, 0])]
qi = [q, q, q2]
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(3)]
A = [[2, 0, -1], [0, 2, -1], [-2, -2, 2]]
# C(2l1)=4n+2=10, C(l2)=4n-2=6, C(0)=0 (n=2): <1>, then <(2n-1)/2>=<3/2>
run('b_2^(1) vector', (e, f, k), A, qi, q,
    lambda z, q: [(1, 14), (br(1, z, q), 10), (br(1, z, q)*br(1.5, z, q), 1)], N)

# ---------------- c_3^(1) vector (C^6) -----------------
# basis 1..6 = eps1,eps2,eps3,-eps3,-eps2,-eps1
N = 6
e = [Eij(N, 6, 1), Eij(N, 1, 2) + Eij(N, 5, 6), Eij(N, 2, 3) + Eij(N, 4, 5), Eij(N, 3, 4)]
f = [m.T.copy() for m in e]
h = [np.diag([-1, 0, 0, 0, 0, 1]), np.diag([1, -1, 0, 0, 1, -1]), np.diag([0, 1, -1, 1, -1, 0]), np.diag([0, 0, 1, -1, 0, 0])]
qi = [q, t, t, q]
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(4)]
A = [[2, -1, 0, 0], [-2, 2, -1, 0], [0, -1, 2, -2], [0, 0, -1, 2]]
# C(2l1)=2n+2=8, C(l2)=2n=6, C(0)=0: rho_l2=<1/2>, rho_0=<(n+1)/2>=<2>
run('c_3^(1) vector', (e, f, k), A, qi, q,
    lambda z, q: [(1, 21), (br(0.5, z, q), 14), (br(2, z, q), 1)], N)

# ---------------- b_3^(1) vector (C^7) -----------------
# basis 1..7 = e1,e2,e3,0,-e3,-e2,-e1 ; alpha_1=e1-e2, alpha_2=e2-e3 long, alpha_3=e3 short, alpha_0=-e1-e2
N = 7
e = [Eij(N, 7, 2) + Eij(N, 6, 1), Eij(N, 1, 2) + Eij(N, 6, 7), Eij(N, 2, 3) + Eij(N, 5, 6), qnum(2, t)*Eij(N, 3, 4) + Eij(N, 4, 5)]
f = [Eij(N, 2, 7) + Eij(N, 1, 6), Eij(N, 2, 1) + Eij(N, 7, 6), Eij(N, 3, 2) + Eij(N, 6, 5), Eij(N, 4, 3) + qnum(2, t)*Eij(N, 5, 4)]
h = [np.diag([-1, -1, 0, 0, 0, 1, 1]), np.diag([1, -1, 0, 0, 0, 1, -1]), np.diag([0, 1, -1, 0, 1, -1, 0]), np.diag([0, 0, 2, 0, -2, 0, 0])]
qi = [q, q, q, t]
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(4)]
A = [[2, 0, -1, 0], [0, 2, -1, 0], [-1, -1, 2, -1], [0, 0, -2, 2]]
# C(2l1)=4n+2=14, C(l2)=4n-2=10: <1>, <(2n-1)/2>=<5/2>
run('b_3^(1) vector', (e, f, k), A, qi, q,
    lambda z, q: [(1, 27), (br(1, z, q), 21), (br(1, z, q)*br(2.5, z, q), 1)], N)

# ---------------- negative controls: wrong normalisations must fail -----------------
from qaff import solve_jimbo as sj
N = 4
e = [Eij(N, 4, 1), Eij(N, 1, 2) + Eij(N, 3, 4), Eij(N, 2, 3)]
f = [m.T.copy() for m in e]
h = [np.diag([-1, 0, 0, 1]), np.diag([1, -1, 1, -1]), np.diag([0, 1, -1, 0])]
qi = [q, t, q]
k = [np.diag(qi[i]**np.diag(h[i]).astype(complex)) for i in range(3)]
x, y = 0.3+0.9j, -1.2+0.4j; z = x/y
R, _ = sj((e, f, k), x, y); ev = np.linalg.eigvals(R/R[0, 0])
for label, want in (('C with long roots length^2 4 (<1>,<3>)', [(1, 10), (br(1, z, q), 5), (br(3, z, q), 1)]),
                    ('rule with z -> 1/z', [(1, 10), (br(0.5, 1/z, q), 5), (br(1.5, 1/z, q), 1)]),
                    ('book conventions (<1/2>,<3/2>)', [(1, 10), (br(0.5, z, q), 5), (br(1.5, z, q), 1)])):
    print(f'   control c_2^(1): eigenvalues match {label}: {match(ev, want)}')
