"""Universal R-matrix of U_q(sl2), eq-universal-r-sl2 (chapter 7).

Checks R Delta(a) = Delta^op(a) R for the coproduct eq-coproduct on spin 1/2
and spin 1, for R = q^{H(x)H/2} sum_n q^{n(n-1)/2}(q-1/q)^n/[n]! E^n (x) F^n,
and that the reversed ordering fails on spin 1. Also prints the spin-1/2 matrix.
"""
import sympy as sp

q = sp.symbols('q', positive=True)
qn = lambda n: (q**n - q**-n)/(q - 1/q)
kp = sp.kronecker_product


def rep(twoj):
    d = twoj + 1
    E, F = sp.zeros(d), sp.zeros(d)
    for k in range(d-1):
        F[k+1, k] = 1
    for k in range(1, d):
        E[k-1, k] = qn(k)*qn(twoj+1-k)  # eq-spin-j
    Hd = [twoj - 2*k for k in range(d)]
    return E, F, sp.diag(*[q**h for h in Hd]), Hd


def check(twoj):
    E, F, K, Hd = rep(twoj); d = twoj + 1; I = sp.eye(d)
    ok_rep = sp.simplify(E*F - F*E - (K - K.inv())/(q - 1/q)) == sp.zeros(d)
    qHH = sp.diag(*[q**sp.Rational(a*b, 2) for a in Hd for b in Hd])
    fact = lambda n: sp.prod([qn(j) for j in range(1, n+1)])
    Ssum = sp.zeros(d*d)
    for n in range(d):
        Ssum += (q**sp.Rational(n*(n-1), 2)*(q - 1/q)**n/fact(n))*kp(E**n, F**n)
    Pm = sp.zeros(d*d)
    for i in range(d):
        for j in range(d):
            Pm[d*j+i, d*i+j] = 1
    D = [kp(E, K) + kp(I, E), kp(F, I) + kp(K.inv(), F), kp(K, K)]
    good = lambda R: all(sp.simplify(R*M - Pm*M*Pm*R) == sp.zeros(d*d) for M in D)
    return ok_rep, good(qHH*Ssum), good(Ssum*qHH), qHH*Ssum


for twoj in (1, 2):
    ok_rep, ok, rev, R = check(twoj)
    print(('PASS ' if ok_rep else 'FAIL ') + f'spin {twoj}/2 is a representation')
    print(('PASS ' if ok else 'FAIL ') + f'universal R intertwines on spin {twoj}/2 (x) spin {twoj}/2')
    if twoj == 2:
        print(('PASS ' if not rev else 'FAIL ') + 'reversed ordering fails on spin 1')
    if twoj == 1:
        print('spin-1/2 matrix / q^(1/2):'); sp.pprint(sp.simplify(R/sp.sqrt(q)))
