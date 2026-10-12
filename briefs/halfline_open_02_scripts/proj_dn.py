"""Matrix part of the lowest-breather reflection for d_n vector solitons: (1 (x) K(th_a)) Rc(th_a + th_b) (1 (x) K(th_b))
applied to the B_1 singlet of V(th_a) (x) V(th_b), th_a,b = th -+ i u/2, u = pi (1 - 1/T); projected on the outgoing
singlet of V(-th_b) (x) V(-th_a)."""
from dn_common import *


def singlet(A, B):
    d = A['dim']*B['dim']
    M = np.vstack([cop(A, B, kd, j) - (np.eye(d) if kd == 'k' else 0) for j in nodes_of(A) for kd in 'efk'])
    _, sv, Vh = np.linalg.svd(M)
    w = Vh[-1].conj()
    return sv[-1]/sv[0], w/w[int(np.argmax(np.abs(w) > 1e-8))]


def cproj(n, om, th, Kf):
    q = q_phys(om); T = (n - 1)*om; u = np.pi*(1 - 1/T); N = 2*n
    ta, tb = th - 0.5j*u, th + 0.5j*u
    Va = dn_vector(n, q, np.exp(om*ta)); Vb = dn_vector(n, q, np.exp(om*tb))
    r_in, w = singlet(Va, Vb)
    r_out, wo = singlet(dn_vector(n, q, np.exp(-om*tb)), dn_vector(n, q, np.exp(-om*ta)))
    v1 = kron(np.eye(N), Kf(tb)) @ w
    R = fast_rcheck(Va, dn_vector(n, q, np.exp(-om*tb)))
    v = kron(np.eye(N), Kf(ta)) @ (R @ v1)
    c = np.vdot(wo, v)/np.vdot(wo, wo)
    return c, np.linalg.norm(v - c*wo)/np.linalg.norm(v), max(r_in, r_out)


if __name__ == '__main__':
    for n in (4, 5):
        for om in (0.37, 1.31):
            q = q_phys(om)
            J = Kdn(n, q, 0.9 + 0.2j, np.zeros(n + 1)); J = J/J[2*n - 1, 0]
            for th in (0.31 + 0.07j, -0.5 + 0.3j):
                c, r, rs = cproj(n, om, th, lambda t: J)
                print(n, om, th, 'c =', np.round(c, 10), f'res {r:.1e} singlet res {rs:.1e}')
