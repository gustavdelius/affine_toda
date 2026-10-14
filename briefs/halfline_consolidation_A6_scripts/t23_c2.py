"""Task 3: c_2 by folding a_3. (i) classical: linearized reflection at phi = 0 for C = +-1 in c_2 itself (roots: long
alpha_0, alpha_2 of length^2 2, short alpha_1 of length^2 1, Kac labels (1,2,1); book V and B with sqrt(2n_j/alpha_j^2));
(ii) quantum: the c_2 S-matrix of the heaviest particle (book eq-cn-snn, H = 4 + B) is not a_3's S_22 = {1}{3}."""
import numpy as np, cmath
# folded a_3 roots in the invariant plane: a_3 roots e_j - e_{j+1} in R^4 (sum-zero), invariant under j -> 4 - j
E = np.eye(4)
a3 = [E[3] - E[0], E[0] - E[1], E[1] - E[2], E[2] - E[3]]           # alpha_0..alpha_3
# invariant subspace of the folding sigma: alpha_1 <-> alpha_3, alpha_0, alpha_2 fixed
rows = [a3[1] - a3[3], np.ones(4)]
_, s, Vt = np.linalg.svd(np.array(rows)); P = Vt[2:].T            # 2-dim subspace orthogonal to both
al = [P.T @ a3[0], P.T @ (a3[1] + a3[3])/2, P.T @ a3[2]]             # c_2 roots (long, short, long)
nj = np.array([1, 2, 1.0])
print('root lengths^2', [round(a @ a, 12) for a in al], ' sum n_j alpha_j =', np.round(sum(n*a for n, a in zip(nj, al)), 12))
M2 = sum(n*np.outer(a, a) for n, a in zip(nj, al))
Rb = sum(np.sqrt(2*n/(a @ a))*np.outer(a, a) for n, a in zip(nj, al))/2
ev, U = np.linalg.eigh(M2); ma = np.sqrt(ev)
print('masses', ma, ' BC matrix = M2/2:', np.allclose(Rb, M2/2))
blk = lambda x, th, h: cmath.sinh(th/2 + 1j*cmath.pi*x/(2*h))/cmath.sinh(th/2 - 1j*cmath.pi*x/(2*h))
worst = 0
for C in (1, -1):
    for k in range(2):
        for th in (0.37, 1.1 + 0.2j):
            p = ma[k]*cmath.sinh(th)
            # eps = rho(e^{ipx} + K e^{-ipx}), eps'(0) = -C Rb eps(0): ip(1 - K) = -C (m^2/2)(1 + K)
            K = (1j*p + C*ma[k]**2/2)/(1j*p - C*ma[k]**2/2)
            a = 1 if abs(ma[k] - np.sqrt(2)) < 1e-9 else 2
            pred = -(blk(a, th, 4)*blk(4 - a, th, 4))**(-C)
            worst = max(worst, abs(K - pred))
print('c_2 linearized reflection at phi = 0 = folded a_3 targets -[(1)(3)]^{-C}, -[(2)^2]^{-C}:', worst)
# quantum: S_22 of c_2 vs a_3
br = lambda x, h, B, th: blk(x - 1, th, h)*blk(x + 1, th, h)/(blk(x - 1 + B, th, h)*blk(x + 1 - B, th, h))
for B in (0.4, 1.0):
    H = 4 + B
    th = 0.3 + 0.2j
    Sc2 = br(1, H, B, th)*br(H - 1, H, B, th)*br(3, H, B, th)*br(H - 3, H, B, th)
    Sa3 = br(1, 4, B, th)*br(3, 4, B, th)
    print(f'B = {B}: S_22(c_2)/S_22(a_3) = {Sc2/Sa3:.6f}; mass ratio m_1/m_2 = {np.sin(np.pi/H)/np.sin(2*np.pi/H):.6f} '
          f'(a_3: {1/np.sqrt(2):.6f})')
