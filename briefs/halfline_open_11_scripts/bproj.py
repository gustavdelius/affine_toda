"""B1 singlet projection for the SP diagonal K: B1 (in V(ta) x Vbar(tb)) -> Bbar1 (in Vbar(-tb) x V(-ta)), and conjugate."""
import numpy as np
from boundary_common import q_phys, an_vector, cop, nodes_of, kron
from nk_common import fast_rcheck
from sp_core import KV, KVbar

def singlet(A, B):
    d = A['dim']*B['dim']
    M = np.vstack([cop(A, B, kd, j) - (np.eye(d) if kd == 'k' else 0) for j in nodes_of(A) for kd in 'efk'])
    _, sv, Vh = np.linalg.svd(M); w = Vh[-1].conj()
    return sv[-1]/sv[0], w/w[0]           # normalise the coefficient of basis vector 0 (x) 0 to 1

def cproj_sp(n, om, l, eta, th, conj=False):
    q = q_phys(om); T = (n+1)*om/2; u = np.pi*(1 - 1/T); N = n+1
    lam = np.exp(-1j*np.pi*om*(1 - l))
    ta, tb = th - 0.5j*u, th + 0.5j*u
    Va = an_vector(n, q, np.exp(om*ta), conj); Vb = an_vector(n, q, np.exp(om*tb), not conj)
    r_in, w = singlet(Va, Vb)
    Vbo = an_vector(n, q, np.exp(-om*tb), not conj); Vao = an_vector(n, q, np.exp(-om*ta), conj)
    r_out, wo = singlet(Vbo, Vao)
    Kv = lambda t: KV(n, om, t, l, eta); Kw = lambda t: lam*KVbar(n, om, t, l, eta)
    Kb = Kw(tb) if not conj else Kv(tb)
    Ka = Kv(ta) if not conj else Kw(ta)
    v = kron(np.eye(N), Kb) @ w
    R = fast_rcheck(Va, Vbo)
    v = kron(np.eye(N), Ka) @ (R @ v)
    c = np.vdot(wo, v)/np.vdot(wo, wo)
    return c, np.linalg.norm(v - c*wo)/np.linalg.norm(v), max(r_in, r_out)

if __name__ == '__main__':
    for n in (1, 2, 3):
        for l in range(0, n+1):
            for th in (0.31+0.07j, -0.2+0.15j):
                c, r, rs = cproj_sp(n, 0.37, l, 0.3+0.4j, th)
                c2, r2, rs2 = cproj_sp(n, 0.37, l, 0.3+0.4j, th, conj=True)
                print(n, l, th, 'c', np.round(c, 6), 'res %.1e sing %.1e' % (r, rs), '| conj c', np.round(c2, 6), 'res %.1e' % r2)
