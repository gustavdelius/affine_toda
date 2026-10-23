"""Soliton-preserving reflection matrices of a_n^(1) in book conventions (Part VII, ch. 33, sec-preserving-conjugating;
half-line open problem 11).

Book conventions: physical q = -exp(-i pi omega), principal gradation y = exp(omega theta), V and Vbar = an_vector,
R-check = rcheck of boundary_common (Part VI coproduct). N = n+1, T = N omega/2 = h omega/2.
Checks:
 1. K_V(theta) = diag(e^{a_j omega theta}) diag(sinh(eta +- T theta)), with + and a_j = 2j on the block S = {0..l-1}
    and - and a_j = 2j - N off S, solves the braid reflection equation (eq-reflection-braid) for n = 2, 3 and every
    l = 0..N; flipping the exponents off S breaks it. In the homogeneous gradation, K_hom = diag(y^{-2j}) K_V, this
    is de Vega-Gonzalez-Ruiz's K_a = e^u sinh(xi+u) (a <= l), e^{-u} sinh(xi-u) (a > l), u = T theta, up to a scalar.
 2. The antisoliton matrix fixed by the mixed reflection equation
      R_{Vbar V}(y1/y2) K_V(y1)_2 R_{V Vbar}(y1 y2) K_Vbar(y2)_2 = K_Vbar(y2)_2 R_{Vbar V}(y1 y2) K_V(y1)_2 R_{V Vbar}(y1/y2)
    is K_Vbar = diag(e^{-a_j omega theta}) diag(sinh(etabar -+ T theta)) with etabar = eta + i pi omega (l - N/2).
 3. Unitarity is scalar: K_V(theta) K_V(-theta) = sinh(eta + T theta) sinh(eta - T theta) 1.
 4. No member of the family is Z_N-covariant for n >= 2: the ratios k_j/k_{j+1} (cyclically) are not all equal.
Runtime: about 2 minutes.
"""
import numpy as np
from boundary_common import an_vector, rcheck, q_phys, one_k, prop_residual, report

OM = 0.37
ETA = 0.41 + 0.23j
Y1, Y2 = np.exp(OM*(0.31 + 0.12j)), np.exp(OM*(-0.57 + 0.05j))


def K_V(n, l, eta, theta, flip=False):
    N, T = n + 1, (n + 1)*OM/2
    s = -1 if flip else 1
    return np.diag([np.exp(2*j*OM*theta)*np.sinh(eta + T*theta) if j < l
                    else np.exp(s*(2*j - N)*OM*theta)*np.sinh(eta - T*theta) for j in range(N)])


def K_Vbar(n, l, etab, theta):
    N, T = n + 1, (n + 1)*OM/2
    return np.diag([np.exp(-2*j*OM*theta)*np.sinh(etab - T*theta) if j < l
                    else np.exp(-(2*j - N)*OM*theta)*np.sinh(etab + T*theta) for j in range(N)])


def make_R(n):
    q, cache = q_phys(OM), {}
    def R(a, b, z):
        key = (a, b, complex(np.round(z, 12)))
        if key not in cache:
            cache[key] = rcheck(an_vector(n, q, z, conj=(a == 'W')), an_vector(n, q, 1.0, conj=(b == 'W')))[0]
        return cache[key]
    return R


th1, th2 = np.log(Y1)/OM, np.log(Y2)/OM
res1, ctrl1, res2, ctrl2, dvgr = [], [], [], [], []
for n in (2, 3):
    N, R = n + 1, make_R(n)
    for l in range(N + 1):
        for flip, out in ((False, res1), (True, ctrl1)):
            K1, K2 = K_V(n, l, ETA, th1, flip), K_V(n, l, ETA, th2, flip)
            L = R('V', 'V', Y1/Y2) @ one_k(K1, N) @ R('V', 'V', Y1*Y2) @ one_k(K2, N)
            Rr = one_k(K2, N) @ R('V', 'V', Y1*Y2) @ one_k(K1, N) @ R('V', 'V', Y1/Y2)
            if not (flip and l == N):           # for l = N nothing is off S, so the control is void
                out.append(prop_residual(L, Rr)[0])
        # gradation map to de Vega-Gonzalez-Ruiz
        T = N*OM/2
        Khom = np.diag(Y1**(-2*np.arange(N))) @ K_V(n, l, ETA, th1)
        dv = np.array([np.exp(T*th1)*np.sinh(ETA + T*th1) if j < l else np.exp(-T*th1)*np.sinh(ETA - T*th1)
                       for j in range(N)])
        dvgr.append(prop_residual(np.diag(Khom), dv)[0])
        # antisoliton partner
        for etab, out in ((ETA + 1j*np.pi*OM*(l - N/2), res2), (ETA, ctrl2)):
            KV1, KW2 = K_V(n, l, ETA, th1), K_Vbar(n, l, etab, th2)
            L = R('W', 'V', Y1/Y2) @ one_k(KV1, N) @ R('V', 'W', Y1*Y2) @ one_k(KW2, N)
            Rr = one_k(KW2, N) @ R('W', 'V', Y1*Y2) @ one_k(KV1, N) @ R('V', 'W', Y1/Y2)
            if out is res2 or l not in (0, N) and 2*l != N:    # control only where etabar != eta matters
                out.append(prop_residual(L, Rr)[0])
report(f'soliton-preserving K_V solves the braid reflection equation, n = 2, 3, l = 0..N [{max(res1):.0e}]; '
       f'flipped exponents fail [min {min(ctrl1):.1e}]; homogeneous gradation gives de Vega-Gonzalez-Ruiz '
       f'[{max(dvgr):.0e}]', max(res1) < 1e-13 and min(ctrl1) > 1e-3 and max(dvgr) < 1e-13)
report(f'mixed reflection equation: K_Vbar with etabar = eta + i pi omega (l - N/2) [{max(res2):.0e}]; '
       f'etabar = eta fails [min {min(ctrl2):.1e}]', max(res2) < 1e-13 and min(ctrl2) > 1e-3)

dev3, minspread = 0, np.inf
for n in (2, 3, 4):
    N, T = n + 1, (n + 1)*OM/2
    for l in range(N + 1):
        P = K_V(n, l, ETA, th1) @ K_V(n, l, ETA, -th1)
        dev3 = max(dev3, np.linalg.norm(P - np.sinh(ETA + T*th1)*np.sinh(ETA - T*th1)*np.eye(N)))
        k = np.diag(K_V(n, l, ETA, th1))
        ratios = k/np.roll(k, -1)
        minspread = min(minspread, np.max(np.abs(ratios - ratios[0])))
report(f'unitarity K_V(theta) K_V(-theta) = sinh(eta + T theta) sinh(eta - T theta) [{dev3:.0e}]; no member is '
       f'Z_N-covariant for n = 2, 3, 4 (spread of k_j/k_(j+1): min {minspread:.1e})', dev3 < 1e-13 and minspread > 1e-3)
