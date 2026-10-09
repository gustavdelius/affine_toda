import numpy as np
from a22 import rep, W
from rsolve import intertwiner
P = np.zeros((9, 9))
for i in range(3):
    for j in range(3): P[j*3+i, i*3+j] = 1
def cubic(x, q):
    # TW (5.5): (x-q^4)(x+q^6) l^3 + q^6(2+q^2)(l^2+x^2 l) + x(q^2-1)(1-3q^4+q^8)(l^2+l) - q^2(1+2q^2)(x^2 l^2 + l) + (1-q^4 x)(1+q^6 x) = 0
    c3 = (x - q**4)*(x + q**6)
    c2 = q**6*(2 + q**2) + x*(q**2 - 1)*(1 - 3*q**4 + q**8) - q**2*(1 + 2*q**2)*x**2
    c1 = q**6*(2 + q**2)*x**2 + x*(q**2 - 1)*(1 - 3*q**4 + q**8) - q**2*(1 + 2*q**2)
    c0 = (1 - q**4*x)*(1 + q**6*x)
    return np.roots([c3, c2, c1, c0])
for om in (0.13, 0.61, 0.87):
    q = np.exp(-1j*np.pi*om)
    x = 2.0
    R = intertwiner(rep(x, q), rep(1.0, q), (0, 1), W, W)[0]; R /= R[0, 0]
    e = np.linalg.eigvals(P @ R)
    print("omega", om, "mine |e|", np.round(np.sort(np.abs(e)), 4))
    for nm, Q in (("q", q), ("1/q", 1/q)):
        for X in (x, 1/x):
            r = cubic(X, Q)
            ok = all(np.min(np.abs(e - z)) < 1e-6 for z in r)
            print("   TW cubic with", nm, "x" if X == x else "1/x", ":", np.round(r, 4), "|r|", np.round(np.abs(r), 4), "in my spectrum" if ok else "")
    # TW claimed doubly degenerate eigenvalues
    for nm, Q in (("q", q), ("1/q", 1/q)):
        for X in (x, 1/x):
            c = [(1 - Q**2*np.sqrt(X))/(Q**2 - np.sqrt(X)), (1 + Q**2*np.sqrt(X))/(Q**2 + np.sqrt(X))]
            print("   TW pair", nm, "x" if X == x else "1/x", np.round(c, 4), [np.min(np.abs(e - z)) < 1e-6 for z in c])
