"""Problem 11 core: soliton-preserving diagonal K in book conventions (principal gradation)."""
import numpy as np
from boundary_common import q_phys, an_vector
from nk_common import fast_rcheck

def KV(n, om, th, l, eta):
    """Diagonal SP matrix V(th)->V(-th): block S={0..l-1}.  K = diag(e^{a_j om th}) diag(sinh(eta +- T th))."""
    N = n + 1; T = N*om/2
    d = []
    for j in range(N):
        if j < l:
            d.append(np.exp(2*j*om*th)*np.sinh(eta + T*th))
        else:
            d.append(np.exp((2*j - N)*om*th)*np.sinh(eta - T*th))
    return np.diag(d)

def R(A, B):
    return fast_rcheck(A, B)

def Vy(n, q, th, om, conj=False):
    return an_vector(n, q, np.exp(om*th), conj)

def one(K, d):
    return np.kron(np.eye(d), K)

def prop(A, B):
    c = np.vdot(B.ravel(), A.ravel())/np.vdot(B.ravel(), B.ravel())
    return np.linalg.norm(A - c*B)/np.linalg.norm(A), c

def KVbar(n, om, th, l, eta):
    """Antisoliton SP matrix Vbar(th)->Vbar(-th): diag(e^{-a_j om th}) diag(sinh(etab -+ T th)), etab = eta + i pi om (l - N/2)."""
    N = n + 1; T = N*om/2; etab = eta + 1j*np.pi*om*(l - N/2)
    d = []
    for j in range(N):
        if j < l:
            d.append(np.exp(-2*j*om*th)*np.sinh(etab - T*th))
        else:
            d.append(np.exp(-(2*j - N)*om*th)*np.sinh(etab + T*th))
    return np.diag(d)

def Phis(n, om, l, eta, th):
    """Crossing functions: kB(th) = Phi1 R_{V Vbar}(2th) kA(-th),  kA(th) = Phi2 R_{Vbar V}(2th) kB(-th)."""
    from boundary_common import q_phys
    q = q_phys(om); N = n + 1; a = 0.5j*np.pi
    V = lambda t: Vy(n, q, t, om); W = lambda t: Vy(n, q, t, om, True)
    E = np.eye(N)
    kA = lambda t: sum(KVbar(n, om, a - t, l, eta)[i, i]*np.kron(E[i], E[i]) for i in range(N))
    kB = lambda t: sum(KV(n, om, a - t, l, eta)[i, i]*np.kron(E[i], E[i]) for i in range(N))
    r1, c1 = prop(kB(th), R(V(th), W(-th)) @ kA(-th))
    r2, c2 = prop(kA(th), R(W(th), V(-th)) @ kB(-th))
    return c1, c2, max(r1, r2)
