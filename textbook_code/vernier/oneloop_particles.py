"""One-loop particle mass ratios of real-coupling Toda theory, in the foundations normalization
   (L = 1/2 (d phi)^2 - (m^2/beta^2) sum_j n_j e^{beta alpha_j.phi}, longest roots |alpha|^2 = 2).
   Extract dH/d(beta^2) of the floating Coxeter number from each mass ratio.
   Imaginary coupling: beta^2 -> -beta^2, so H = h - (dH/dbeta^2) beta^2 + O(beta^4)."""
import numpy as np
from scipy.integrate import quad
def toda(roots, kac):
    A = np.array(roots, float); n = np.array(kac, float)
    assert np.allclose(n @ A, 0)
    M2 = sum(ni*np.outer(a, a) for ni, a in zip(n, A))
    m2, V = np.linalg.eigh(M2); o = np.argsort(m2); m2, V = m2[o], V[:, o]
    Ae = A @ V
    C3 = np.einsum('i,ia,ib,ic->abc', n, Ae, Ae, Ae)
    return m2, C3
def J(pa2, mb, mc):
    return quad(lambda x: 1/(x*mb**2 + (1-x)*mc**2 - x*(1-x)*pa2), 0, 1, limit=200)[0]/(4*np.pi)
def dm2(m2, C3):
    m = np.sqrt(m2); r = len(m2)
    return np.array([-0.5*sum(C3[a,b,c]**2*J(m2[a], m[b], m[c]) for b in range(r) for c in range(r) if abs(C3[a,b,c]) > 1e-10) for a in range(r)])
def dH_from(m2, d, ratio_of_H, h, pairs):
    """pairs: list of (i, j, f) meaning m_i/m_j = f(H). Returns dH/dbeta^2 from each."""
    out = []
    for i, j, f in pairs:
        R = np.sqrt(m2[i]/m2[j]); assert abs(R - f(h)) < 1e-9, (R, f(h))
        dR = 0.5*R*(d[i]/m2[i] - d[j]/m2[j])            # dR/dbeta^2
        e = 1e-6; slope = (f(h+e) - f(h-e))/(2*e)
        out.append(dR/slope)
    return out
s2 = 1/np.sqrt(2)
def a2n1_roots(n):
    e = np.eye(n)
    roots = [-s2*(e[0]+e[1])] + [s2*(e[i]-e[i+1]) for i in range(n-1)] + [np.sqrt(2)*e[n-1]]
    kac = [1, 1] + [2]*(n-2) + [1]
    return roots, kac
def cn_roots(n):
    e = np.eye(n)
    roots = [-np.sqrt(2)*e[0]] + [s2*(e[i]-e[i+1]) for i in range(n-1)] + [np.sqrt(2)*e[n-1]]
    kac = [1] + [2]*(n-1) + [1]
    return roots, kac
print("Calibration, c_n^(1) (companion: H = 2n + B, B = 2b^2/(8pi+b^2) at real coupling => dH/db^2 = 1/(4pi) = %.6f):" % (1/(4*np.pi)))
for n in [2, 3, 4]:
    m2, C3 = toda(*cn_roots(n)); d = dm2(m2, C3)
    # classical c_n masses m_a ~ sin(a pi/2n), a = 1..n ; sorted ascending = a order
    pairs = [(a-1, 0, (lambda a: lambda H: np.sin(a*np.pi/H)/np.sin(np.pi/H))(a)) for a in range(2, n+1)]
    print(f"  n={n}: dH/dbeta^2 =", np.round(dH_from(m2, d, None, 2*n, pairs), 6))
print("\na_{2n-1}^(2), foundations normalization (long root |alpha_n|^2 = 2):")
for n in [2, 3, 4, 5]:
    m2, C3 = toda(*a2n1_roots(n)); d = dm2(m2, C3)
    h = 2*n - 1; m = np.sqrt(m2)
    # identify particle n: the one not of the form 2 m_n sin(a pi/h)
    ok = None
    for k in range(n):
        targ = sorted([2*m[k]*np.sin(a*np.pi/h) for a in range(1, n)] + [m[k]])
        if np.allclose(targ, m, atol=1e-9): ok = k
    idx = {}
    for a in range(1, n):
        idx[a] = int(np.argmin(abs(m - 2*m[ok]*np.sin(a*np.pi/h))))
    pairs = [(idx[a], ok, (lambda a: lambda H: 2*np.sin(a*np.pi/H))(a)) for a in range(1, n)]
    r = dH_from(m2, d, None, h, pairs)
    print(f"  n={n}: particle n = eigen-index {ok}; dH/dbeta^2 from m_a/m_n, a=1..n-1: {np.round(r, 7)};  x 8pi = {np.round(np.array(r)*8*np.pi, 6)}")
