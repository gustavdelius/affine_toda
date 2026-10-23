"""Referee check of prp-dn-nondiagonal-scalar (i), (ii), (iv) and the sign rule, at the matrix level.

Uses the book's library boundary_common (representations, coproduct, dense intertwiner solver)
as allowed by the brief; the assembly (normalizations, crossing vector, breather singlet) is my own.

Conventions: q = -exp(-i pi omega), y = exp(omega theta), T = (n-1) omega, mu = -i T theta / pi.
Coideal b_j = e_j + q^-1 f_j k_j + eps_j k_j, eps_0 = s eps_*, eps_{j>=1} = eps_*,
eps_* = 1/(2 cos(pi omega/2));  K(y): V(y) -> V(1/y).
Normalization: K_11 = -sin(pi T/2)/sin(pi mu) (varsigma=+1), -i sin(pi T/2)/cos(pi mu) (varsigma=-1),
varsigma = (-1)^n s.
Checks
 (a) the corner entry of K vanishes at y^h = varsigma, not at y^h = -varsigma;
 (i) K(1/y) K(y) = U_varsigma 1;
 (ii) crossing vector: k(theta) = Psi Rcheck_{VbarVbar}(2 theta) k(-theta), with Psi = Psi_varsigma,
      control: eps-hat = 0 gives Psi = 1;
 (iv) the singlet of V(theta_a) (x) V(theta_b) at the breather pole is mapped by
      (1 x K(theta_a)) Rcheck (1 x K(theta_b)) to c^d_varsigma times the singlet of V(-theta_b) (x) V(-theta_a);
      control eps-hat = 0 (K = J): coefficient 1.
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
import boundary_common as bc   # noqa: E402

pi = np.pi


def dn_bar(n, q, y):
    V1 = bc.dn_vector(n, q, 1.0)
    g = {}
    for j in bc.nodes_of(V1):
        E = V1[('e', j)]; F = V1[('f', j)]; K = V1[('k', j)]
        g[('e', j)] = y * F; g[('f', j)] = E / y; g[('k', j)] = np.linalg.inv(K)
    rep = dict(g); rep['hw'] = 2 * n - 1; rep['dim'] = 2 * n
    return rep


def Kmat(n, om, th, s, eps0=False, raw=False):
    q = bc.q_phys(om)
    y = np.exp(om * th)
    es = 1 / (2 * np.cos(pi * om / 2))
    eps = [0.0] * (n + 1) if eps0 else [s * es] + [es] * n
    c = [1 / q] * (n + 1)
    null, rel = bc.k_intertwiner(bc.dn_vector(n, q, y), bc.dn_vector(n, q, 1 / y), c, eps)
    assert len(null) == 1, (len(null), rel[:3])
    K = null[0]
    if raw:
        return K
    T = (n - 1) * om
    mu = -1j * T * th / pi
    if eps0:
        return K / K[2 * n - 1, 0]          # J normalized to +1 on v_1 -> v_-1
    vs = (-1) ** n * s
    K11 = -np.sin(pi * T / 2) / np.sin(pi * mu) if vs > 0 else -1j * np.sin(pi * T / 2) / np.cos(pi * mu)
    return K / K[0, 0] * K11


def U_s(mu, T, om, vs):
    c2 = np.cos(2 * pi * mu)
    return vs * (c2 - vs * np.cos(pi * T)) * (c2 - vs * np.cos(pi * om)) / (2 * np.sin(pi * om / 2) ** 2 * (c2 - vs))


def Psi_s(mu, om, vs):
    if vs > 0:
        return np.sin(pi * (om / 2 - mu)) / np.sin(pi * (om / 2 + mu))
    return np.cos(pi * (mu - om / 2)) / np.cos(pi * (mu + om / 2))


def c_d(mu, T, om, vs):
    s2 = np.sin(pi * om / 2) ** 2
    if vs > 0:
        c = -np.cos(pi * mu) * np.cos(pi * (mu + om / 2)) * np.cos(pi * (mu - T)) / (s2 * np.cos(pi * (mu - om / 2)))
        e = np.cos(pi * (mu - T / 2 + om / 2)) * np.cos(pi * (mu - T / 2 - om / 2)) / np.cos(pi * (mu - T / 2)) ** 2
    else:
        c = np.sin(pi * mu) * np.sin(pi * (mu + om / 2)) * np.sin(pi * (mu - T)) / (s2 * np.sin(pi * (mu - om / 2)))
        e = np.sin(pi * (mu - T / 2 + om / 2)) * np.sin(pi * (mu - T / 2 - om / 2)) / np.sin(pi * (mu - T / 2)) ** 2
    return c * e


def conj_map(n, q):
    V, Vb = bc.dn_vector(n, q, 1.3), dn_bar(n, q, 1.3)
    pairs = [(V[(k, j)], Vb[(k, j)]) for j in bc.nodes_of(V) for k in 'efk']
    null, rel = bc.intertwiners(pairs)
    assert len(null) == 1
    C = null[0]
    # check y-independence
    V2, Vb2 = bc.dn_vector(n, q, 0.6), dn_bar(n, q, 0.6)
    assert max(np.abs(C @ V2[(k, j)] - Vb2[(k, j)] @ C).max() for j in bc.nodes_of(V) for k in 'efk') < 1e-10
    return C


def weights(rep):
    """integer weight vectors (H_1..H_n eigenvalues) of the basis vectors of a dn_vector/dn_bar rep"""
    nodes = [j for j in bc.nodes_of(rep) if j != 0]
    cols = []
    for j in nodes:
        e, f = rep[('e', j)], rep[('f', j)]
        cols.append(np.round(np.real(np.diag(e @ f - f @ e))).astype(int))
    return [tuple(v) for v in np.array(cols).T]


def rcheck_w(A, B):
    """R-check A (x) B -> B (x) A for the Part VI coproduct, unknowns restricted to
    weight-preserving entries; normalized to 1 on hw_A (x) hw_B -> hw_B (x) hw_A."""
    import scipy.sparse as sps
    wA, wB = weights(A), weights(B)
    dA, dB = A['dim'], B['dim']
    win = [tuple(np.add(wA[i], wB[j])) for i in range(dA) for j in range(dB)]
    wout = [tuple(np.add(wB[k], wA[l])) for k in range(dB) for l in range(dA)]
    byw = {}
    for o, w in enumerate(wout):
        byw.setdefault(w, []).append(o)
    unk = [(o, i) for i, w in enumerate(win) for o in byw.get(w, [])]
    col = {u: c for c, u in enumerate(unk)}
    rows, cols, vals = [], [], []
    r0 = 0
    D = dA * dB
    for j in bc.nodes_of(A):
        for kind in 'ef':
            X = bc.cop(A, B, kind, j); Y = bc.cop(B, A, kind, j)
            # equation (R X - Y R)[o, i] = 0 for all o, i
            for (o, i), c in col.items():
                # R[o,i] contributes to (R X)[o, i'] with X[i, i'] and to (Y R)[o', i] with Y[o', o]
                for ip in np.nonzero(X[i])[0]:
                    rows.append(r0 + o * D + ip); cols.append(c); vals.append(X[i, ip])
                for op in np.nonzero(Y[:, o])[0]:
                    rows.append(r0 + op * D + i); cols.append(c); vals.append(-Y[op, o])
            r0 += D * D
    M = sps.csr_matrix((vals, (rows, cols)), shape=(r0, len(unk)))
    G = (M.conj().T @ M).toarray()
    ev, U = np.linalg.eigh(G)
    assert ev[0] / ev[-1] < 1e-14 and ev[1] / ev[-1] > 1e-8, (ev[:3] / ev[-1])
    R = np.zeros((D, D), complex)
    for (o, i), c in col.items():
        R[o, i] = U[c, 0]
    i_in = A['hw'] * dB + B['hw']; i_out = B['hw'] * dA + A['hw']
    return R / R[i_out, i_in], 1


def singlet(A, B):
    """null vector of Delta(e_j), Delta(f_j), Delta(k_j)-1 on A (x) B"""
    M = []
    for j in bc.nodes_of(A):
        M.append(bc.cop(A, B, 'e', j)); M.append(bc.cop(A, B, 'f', j))
        M.append(bc.cop(A, B, 'k', j) - np.eye(A['dim'] * B['dim']))
    M = np.vstack(M)
    _, s, Vh = np.linalg.svd(M)
    assert s[-1] / s[0] < 1e-10 and s[-2] / s[0] > 1e-6, (s[-2:] / s[0])
    return Vh[-1].conj()


def run():
    ok = True
    for n in (4, 5):
        h = 2 * n - 2
        for om in (0.37, 0.71):
            q = bc.q_phys(om)
            T = (n - 1) * om
            C = conj_map(n, q)
            for s in (+1, -1):
                vs = (-1) ** n * s
                # (a) corner entries at y^h = +-1
                out = {}
                for sign, label in ((1, '+1'), (-1, '-1')):
                    # y^h = +-1 exactly is degenerate (x = 1/x): approach it, the corner entry is O(dtheta)
                    th = (2j * pi if sign > 0 else 1j * pi) / (h * om) + 1e-6
                    K = Kmat(n, om, th, s, raw=True)
                    out[label] = max(abs(K[0, 2 * n - 1]), abs(K[2 * n - 1, 0])) / abs(K[0, 0])
                zero_at = '+1' if out['+1'] < 1e-4 < out['-1'] else ('-1' if out['-1'] < 1e-4 < out['+1'] else 'none')
                ok &= (zero_at == ('+1' if vs > 0 else '-1'))
                # (i) unitarity
                th = 0.31 + 0.17j
                mu = -1j * T * th / pi
                K, Km = Kmat(n, om, th, s), Kmat(n, om, -th, s)
                U = Km @ K
                resU = np.abs(U - U_s(mu, T, om, vs) * np.eye(2 * n)).max() / abs(U_s(mu, T, om, vs))
                # (ii) crossing vector
                def kvec(t):
                    Kc = C @ Kmat(n, om, 1j * pi / 2 - t, s)      # V -> Vbar
                    return Kc.T.reshape(-1)                       # component (a,b) = Kc[b,a]
                R, nul = rcheck_w(dn_bar(n, q, np.exp(om * th)), dn_bar(n, q, np.exp(-om * th)))
                kp, km = kvec(th), R @ kvec(-th)
                i0 = np.argmax(abs(km))
                psi = kp[i0] / km[i0]
                resP = max(np.abs(kp - psi * km).max() / np.abs(kp).max(), abs(psi / Psi_s(mu, om, vs) - 1))
                # (iv) breather singlet
                tha = th - 1j * pi / 2 + 1j * pi / (2 * T)
                thb = th + 1j * pi / 2 - 1j * pi / (2 * T)
                ya, yb = np.exp(om * tha), np.exp(om * thb)
                sv = singlet(bc.dn_vector(n, q, ya), bc.dn_vector(n, q, yb))
                sv = sv / sv[2 * n - 1]                           # normalize BEFORE mapping (v_1 (x) v_-1 -> 1)
                Rab, _ = rcheck_w(bc.dn_vector(n, q, ya), bc.dn_vector(n, q, 1 / yb))
                w = bc.one_k(Kmat(n, om, tha, s), 2 * n) @ Rab @ bc.one_k(Kmat(n, om, thb, s), 2 * n) @ sv
                sv2 = singlet(bc.dn_vector(n, q, 1 / yb), bc.dn_vector(n, q, 1 / ya))
                i1 = 0 * 2 * n + (2 * n - 1)                      # v_1 (x) v_-1
                sv2 = sv2 / sv2[i1]
                cc = w[i1]
                resS = max(np.abs(w - cc * sv2).max() / np.abs(w).max(), abs(cc / c_d(mu, T, om, vs) - 1))
                print(f"n={n} om={om} s={s:+d} (varsigma={vs:+d}): corner zero at y^h={zero_at} "
                      f"[{out['+1']:.1e}, {out['-1']:.1e}]; U {resU:.1e}; Psi {resP:.1e}; c^d {resS:.1e}")
                ok &= max(resU, resP, resS) < 1e-8
            # controls at eps-hat = 0
            th = 0.31 + 0.17j
            def kvec0(t):
                return (C @ Kmat(n, om, 1j * pi / 2 - t, 1, eps0=True)).T.reshape(-1)
            R, _ = rcheck_w(dn_bar(n, q, np.exp(om * th)), dn_bar(n, q, np.exp(-om * th)))
            kp, km = kvec0(th), R @ kvec0(-th)
            i0 = np.argmax(abs(km)); psi0 = kp[i0] / km[i0]
            tha = th - 1j * pi / 2 + 1j * pi / (2 * T); thb = th + 1j * pi / 2 - 1j * pi / (2 * T)
            ya, yb = np.exp(om * tha), np.exp(om * thb)
            sv = singlet(bc.dn_vector(n, q, ya), bc.dn_vector(n, q, yb))
            Rab, _ = rcheck_w(bc.dn_vector(n, q, ya), bc.dn_vector(n, q, 1 / yb))
            J = Kmat(n, om, 0.2, 1, eps0=True)
            w = bc.one_k(J, 2 * n) @ Rab @ bc.one_k(J, 2 * n) @ sv
            sv2 = singlet(bc.dn_vector(n, q, 1 / yb), bc.dn_vector(n, q, 1 / ya))
            i1 = 2 * n - 1
            c0 = (w / sv[i1])[i1] / (sv2 / sv2[i1])[i1]
            print(f"   control eps-hat=0: Psi = {psi0:.12f} (expect 1), singlet coefficient {c0:.12f} (expect 1)")
            ok &= abs(psi0 - 1) < 1e-9 and abs(c0 - 1) < 1e-9
    print("PASS" if ok else "FAIL", "d_n matrix data: sign rule, U, Psi, singlet coefficient (n=4,5)")


if __name__ == '__main__':
    run()
