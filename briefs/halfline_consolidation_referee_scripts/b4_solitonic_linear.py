"""Referee check of prp-solitonic-linear (ch. 31) and of the singular-endpoint bullets of
sec-boundary-missing-charges, by direct ODE integration (a different method from the book's
script, which uses Bowcock's exact Hirota dressing).

a_n^(1) real coupling, m = beta = 1, phi in R^h (sum zero), alpha_j = e_j - e_{j+1} (cyclic),
phi = -sum_j alpha_j ln tau_j, e^{alpha_k.phi} = tau_{k-1} tau_{k+1} / tau_k^2.
Linearization eps'' = (M^2(x) - Omega^2) eps, M^2 = sum_k alpha_k alpha_k^T e^{alpha_k.phi},
boundary eps'(0) = -R eps(0), R = 1/2 sum_k C_k alpha_k alpha_k^T e^{alpha_k.phi(0)/2}.
Particle c <-> v_c = (omega^{ic})_i / sqrt(h); K_{c'c} = outgoing v_{c'} per unit incoming v_c.
"""
import itertools
import sys
import cmath
import numpy as np
from scipy.integrate import solve_ivp

pi = np.pi


def pair_background(n, a, chi, d=None):
    """static species-a pair; d = None: the continued vacuum (d = 1/(1+s));
    returns tau(x) (h-vector) and tau'(x)"""
    h = n + 1; ma = 2 * np.sin(pi * a / h); s = np.sin(pi * a / h)
    if d is None:
        d = 1 / (1 + s)
    c2 = d * d * (1 - s * s)
    js = np.arange(h)
    cs = np.cos(chi + 2 * pi * js * a / h)

    def tau(x):
        E = np.exp(ma * x)
        return 1 + 2 * d * cs * E + c2 * E * E, ma * (2 * d * cs * E + 2 * c2 * E * E)
    return tau


def selfconj_background(n, d):
    h = n + 1; js = np.arange(h)

    def tau(x):
        E = np.exp(2 * x)
        return 1 + (-1.0) ** js * d * E, 2 * (-1.0) ** js * d * E
    return tau


def roots(n):
    h = n + 1
    al = np.zeros((h, h))
    for k in range(h):
        al[k, k] += 1; al[k, (k + 1) % h] -= 1
    return al            # row k: alpha_k = e_k - e_{k+1}, node k (k = 0 is alpha_0 = e_0 - e_1 in this cyclic labelling)


def expo(t):
    h = len(t)
    return np.array([t[(k - 1) % h] * t[(k + 1) % h] / t[k] ** 2 for k in range(h)])


def pattern(n, tau):
    t, tp = tau(0.0)
    h = n + 1
    r = tp / t
    e = np.sqrt(expo(t))
    best = None
    for C in itertools.product([1, -1], repeat=h):
        v = r - np.array(C) * e
        spread = np.ptp(v)
        if best is None or spread < best[0]:
            best = (spread, C)
    return best[1], best[0]


def kmatrix(n, tau, C, Om, L=16.0):
    h = n + 1
    al = roots(n)
    om = np.exp(2j * pi / h)
    V = np.array([[om ** (i * c) for i in range(h)] for c in range(1, h)]).T / np.sqrt(h)   # columns v_c, c=1..n
    masses = np.array([2 * np.sin(pi * c / h) for c in range(1, h)])
    t0, _ = tau(0.0)
    R = 0.5 * (al.T * (np.array(C) * np.sqrt(expo(t0)))) @ al

    def rhs(x, y):
        t, _ = tau(x)
        M2 = (al.T * expo(t)) @ al
        return np.concatenate([y[h:], (M2 - Om ** 2 * np.eye(h)) @ y[:h]])
    opn = [i for i in range(n) if Om > masses[i]]
    cls = [i for i in range(n) if Om <= masses[i]]
    finals = {}
    for i in opn:
        k = np.sqrt(Om ** 2 - masses[i] ** 2)
        for sg, lab in ((1, 'in'), (-1, 'out')):
            amp = np.exp(1j * sg * k * (-L))
            y0 = np.concatenate([V[:, i] * amp, 1j * sg * k * V[:, i] * amp])
            sol = solve_ivp(rhs, (-L, 0.0), y0, method='DOP853', rtol=1e-12, atol=1e-14)
            yf = sol.y[:, -1]
            finals[(lab, i)] = V.conj().T @ (yf[h:] + R @ yf[:h])
    for i in cls:
        kap = np.sqrt(masses[i] ** 2 - Om ** 2)
        amp = np.exp(-kap * L)
        y0 = np.concatenate([V[:, i] * amp, kap * V[:, i] * amp]).astype(complex)
        sol = solve_ivp(rhs, (-L, 0.0), y0, method='DOP853', rtol=1e-12, atol=1e-14)
        yf = sol.y[:, -1]
        finals[('dec', i)] = V.conj().T @ (yf[h:] + R @ yf[:h])
    unknown = [('out', i) for i in opn] + [('dec', i) for i in cls]
    Mx = np.array([finals[u] for u in unknown]).T
    K = np.zeros((len(opn), len(opn)), complex)
    for col, i in enumerate(opn):
        K[:, col] = np.linalg.solve(Mx, -finals[('in', i)])[:len(opn)]
    th = {i: np.arccosh(Om / masses[i]) for i in opn}
    return K, opn, th


def blk(x, th, h):
    return cmath.sinh(th / 2 + 1j * pi * x / (2 * h)) / cmath.sinh(th / 2 - 1j * pi * x / (2 * h))


def DG(c, b, th, h, B=0.0):
    K = 1
    for k in range(1, c + 1):
        K *= blk(k - 1, th, h) * blk(h + 1 - k - B / 2, th, h) / (blk(h - k, th, h) * blk(k - B / 2, th, h))
    M = 1
    for k in range(1, c + 1):
        M *= (blk(c + b + h / 2 - 2 * k + B / 2, th, h) * blk(b - c - h / 2 + 2 * k - B / 2, th, h)
              / (blk(c + b - h / 2 - 2 * k + B / 2, th, h) * blk(b - c + h / 2 + 2 * k - B / 2, th, h)))
    return K * M


def regular(tau, xs=np.linspace(-12, 0, 2401)):
    return min(np.min(tau(x)[0]) for x in xs) > 0


def main():
    ok = True
    if len(sys.argv) > 1 and sys.argv[1] == 'p3':
        return part3()
    print("1. linearized K about continued vacua vs DG K^(b)_c at B = 0 (b = a)")
    for n in (2, 3, 4, 5):
        h = n + 1
        for a in range(1, (h + 1) // 2):
            if 2 * a == h:
                continue
            ends = sorted(((pi - 2 * pi * j * a / h) % (2 * pi)) for j in range(h))
            ends = sorted(set(np.round(ends, 12)))
            mids = [(ends[i] + ends[(i + 1) % len(ends)] + (2 * pi if i == len(ends) - 1 else 0)) / 2 for i in range(len(ends))]
            for chi in mids[:3]:
                tau = pair_background(n, a, chi)
                if not regular(tau):
                    print(f"   a_{n} a={a} chi={chi:.3f}: not regular (skipped)")
                    continue
                C, spread = pattern(n, tau)
                worst_d = worst_o = 0
                for Om in (2 * np.sin(pi / h) * np.cosh(0.83), 1.9, 2.7):
                    K, opn, th = kmatrix(n, tau, C, Om)
                    for p, i in enumerate(opn):
                        worst_d = max(worst_d, abs(K[p, p] - DG(i + 1, a, th[i], h)))
                        for q, j in enumerate(opn):
                            if q != p:
                                worst_o = max(worst_o, abs(K[p, q]))
                print(f"   a_{n} a={a} chi={chi:.3f} pattern {''.join('+' if c > 0 else '-' for c in C)} (bc spread {spread:.0e}):"
                      f" |K_cc - DG| {worst_d:.1e}, off-diagonal {worst_o:.1e}")
                ok &= worst_d < 1e-8 and worst_o < 1e-8 and spread < 1e-10
    print("   self-conjugate family tau_j = 1 + (-1)^j d e^{2x}, C = -1, b = h/2")
    for n in (3, 5):
        h = n + 1
        for d in (0.3, 0.8):
            tau = selfconj_background(n, d)
            C, spread = pattern(n, tau)
            worst_d = worst_o = 0
            for Om in (2 * np.sin(pi / h) * np.cosh(0.83), 1.9, 2.7):
                K, opn, th = kmatrix(n, tau, C, Om)
                for p, i in enumerate(opn):
                    worst_d = max(worst_d, abs(K[p, p] - DG(i + 1, h // 2, th[i], h)))
                    for q, j in enumerate(opn):
                        if q != p:
                            worst_o = max(worst_o, abs(K[p, q]))
            print(f"   a_{n} d={d} pattern {''.join('+' if c > 0 else '-' for c in C)}: |K - DG| {worst_d:.1e}, off-diagonal {worst_o:.1e}")
            ok &= worst_d < 1e-8 and worst_o < 1e-8
    print("2. a_3 middle-soliton endpoint d -> 1 with the alternating pattern's boundary Hessian")
    n, h = 3, 4
    for e in (1e-2, 1e-3, 1e-4):
        tau = selfconj_background(n, 1 - e)
        C = (1, -1, 1, -1)
        Om = 2 * np.sin(pi / h) * np.cosh(0.6)
        K, opn, th = kmatrix(n, tau, C, Om)
        target = -blk(1, th[0], h) * blk(h - 1, th[0], h)
        print(f"   1-d = {e:.0e}: |K_11 - (-(1)(3))| = {abs(K[0, 0] - target):.2e}, max off-diagonal {np.max(np.abs(K - np.diag(np.diag(K)))):.2e}")
    for C in ((-1, 1, -1, 1),):
        tau = selfconj_background(n, 1 - 1e-4)
        K, opn, th = kmatrix(n, tau, C, 2 * np.sin(pi / h) * np.cosh(0.6))
        print(f"   other alternating pattern (-+-+), 1-d = 1e-4: |K_11 + (1)(3)| = {abs(K[0, 0] + blk(1, th[0], h) * blk(3, th[0], h)):.2e}")
    part3()
    print("PASS" if ok else "FAIL", "linearized reflection about the regular continued vacua = DG (5.7) at B = 0 (a_2..a_5)")


def part3():
    print("3. a_5 pattern (+,+,-,-,-,-) at the species-2 pair endpoint chi = 5pi/3: path dependence")
    n, h, a = 5, 6, 2
    C = (1, 1, -1, -1, -1, -1)
    Om = 2 * np.sin(pi / h) * np.cosh(0.6)
    s = np.sin(pi * a / h)
    for dl in (1e-2, 1e-3, 1e-4):
        out = []
        for sgn in (+1, -1):
            tau = pair_background(n, a, 5 * pi / 3 + sgn * dl)
            assert regular(tau)
            K, opn, th = kmatrix(n, tau, C, Om)
            i1, i5 = opn.index(0), opn.index(4)          # particles c = 1 and c = 5 = h-1
            out.append((K[i1, i1], abs(K[i5, i1]), abs(K[i1, i5])))
        tau = pair_background(n, a, 5 * pi / 3, d=(1 - dl) / (1 + s))
        Kd, opnd, _ = kmatrix(n, tau, C, Om)
        i1, i5 = opnd.index(0), opnd.index(4)
        print(f"   delta = {dl:.0e}: chi path K_11 = {out[0][0]:.6f} / {out[1][0]:.6f}, |K_15| = {out[0][1]:.6f} / {out[1][1]:.6f};"
              f"  d path K_11 = {Kd[i1, i1]:.6f}, |K_15| = {abs(Kd[i5, i1]):.6f}")
    print(f"   DG b=2 diagonal value (2)(4)/((1)(5)) = {blk(2, 0.6, h) * blk(4, 0.6, h) / (blk(1, 0.6, h) * blk(5, 0.6, h)):.6f}")


if __name__ == '__main__':
    main()
