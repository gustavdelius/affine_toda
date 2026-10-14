"""Normalisation-free test:  delta H (one-loop floating Coxeter shift) / B_tree (tree-level S-matrix strength), from the Toda Lagrangian."""
import numpy as np
from scipy.integrate import quad
def toda_data(roots, kac):
    A = np.array(roots, float); n = np.array(kac, float)
    assert np.allclose(n @ A, 0)
    M2 = sum(ni*np.outer(a, a) for ni, a in zip(n, A))
    m2, V = np.linalg.eigh(M2); order = np.argsort(m2); m2, V = m2[order], V[:, order]
    Ae = A @ V                                                 # roots in mass eigenbasis
    C3 = np.einsum('i,ia,ib,ic->abc', n, Ae, Ae, Ae)           # cubic couplings (m = beta = 1)
    C4 = np.einsum('i,ia,ib,ic,id->abcd', n, Ae, Ae, Ae, Ae)   # quartic couplings
    return m2, C3, C4
def J(ma2, mb, mc):
    return quad(lambda x: 1/(x*mb**2 + (1 - x)*mc**2 - x*(1 - x)*ma2), 0, 1)[0]/(4*np.pi)
def one_loop_dm2(m2, C3):
    m = np.sqrt(m2); r = len(m2)
    return np.array([-0.5*sum(C3[a, b, c]**2*J(m2[a], m[b], m[c]) for b in range(r) for c in range(r) if abs(C3[a, b, c]) > 1e-9) for a in range(r)])
def tree_over_4m2sinh(m2, C3, C4, a, th):
    s, u, t = 2*m2[a]*(1 + np.cosh(th)), 2*m2[a]*(1 - np.cosh(th)), 0.0
    M = -C4[a, a, a, a] - sum(C3[a, a, c]**2*(1/(s - m2[c]) + 1/(u - m2[c]) + 1/(t - m2[c])) for c in range(len(m2)))
    return M/(4*m2[a]*np.sinh(th))
def blk(x, H, th): return np.sinh(th/2 + 1j*np.pi*x/(2*H))/np.sinh(th/2 - 1j*np.pi*x/(2*H))
def brace(x, B, H, th, nu=0.0):
    return blk(x - nu*B - 1, H, th)*blk(x + nu*B + 1, H, th)/(blk(x + nu*B + B - 1, H, th)*blk(x - nu*B - B + 1, H, th))
def analyse(name, roots, kac, h, S11, ratio_formula, expected):
    m2, C3, C4 = toda_data(roots, kac)
    dm2 = one_loop_dm2(m2, C3)
    # floating Coxeter shift from each mass ratio (relative to the lightest particle 0)
    dH = []
    for a in range(1, len(m2)):
        dR = (dm2[a]*m2[0] - m2[a]*dm2[0])/m2[0]**2
        eps = 1e-6; slope = (ratio_formula(a, h + eps) - ratio_formula(a, h - eps))/(2*eps)
        if abs(slope) > 1e-8: dH.append(dR/slope)
    # tree strength: log S_11 ~ B g(theta) ; compare with i M/(4 m^2 sinh)
    eps = 1e-7; ths = [0.4, 0.9, 1.5, 2.3]
    g = [np.log(S11(h, eps, th)).imag/eps for th in ths]
    tr = [tree_over_4m2sinh(m2, C3, C4, 0, th) for th in ths]
    kap = np.array(tr)/np.array(g)
    print(f"{name}: classical mass ratios {np.round(np.sqrt(m2/m2[0]), 4)}")
    print(f"   one-loop delta H from each mass ratio: {np.round(dH, 6)}")
    print(f"   tree amplitude / g(theta) at theta = {ths}: {np.round(kap, 6)}  (constant => same tree shape)")
    print(f"   delta H / B_tree = {np.mean(dH)/np.mean(kap):.6f}   (expected {expected})")
# c_3^(1): long roots length^2 2
s2 = 1/np.sqrt(2); e = np.eye(3)
c3 = [s2*(e[0]-e[1]), s2*(e[1]-e[2]), s2*2*e[2], -s2*2*e[0]]
dgz11 = lambda H, B, th: brace(1, B, H, th)*brace(H - 1, B, H, th)
analyse("c_3^(1)", c3, [2, 2, 1, 1], 6, lambda h, B, th: dgz11(h + B, B, th), lambda a, H: np.sin((a+1)*np.pi/H)**2/np.sin(np.pi/H)**2, "1 (H = h + B)")
# a_5^(2): finite c_3 roots, alpha_0 = -theta_short = -(e1+e2); long roots length^2 2
a5 = [s2*(e[0]-e[1]), s2*(e[1]-e[2]), s2*2*e[2], -s2*(e[0]+e[1])]
def a5ratio(a, H):            # particles 1,2 float as 2 sin(a pi/H), particle 3 fixed (m_n^2 = 2 relative units) -> use the lightest as reference
    ms = [2*np.sin(np.pi/H), 2*np.sin(2*np.pi/H), 1.0]; ms = sorted(ms); return ms[a]**2/ms[0]**2
analyse("a_5^(2)", a5, [1, 2, 1, 1], 5, lambda h, B, th: dgz11(h + B/2, B, th), a5ratio, "1/2 (H = h + B/2)")
# f_4^(1): Bourbaki roots, long length^2 2
E4 = np.eye(4)
f4 = [E4[1]-E4[2], E4[2]-E4[3], E4[3], 0.5*(E4[0]-E4[1]-E4[2]-E4[3]), -(E4[0]+E4[1])]
def cdsS11(H, B, th):
    return brace(1, B, H, th)*brace(H - 1, B, H, th)*brace(H/3 + 1, B, H, th)*brace(H - H/3 - 1, B, H, th)
def cdsS11_param(h, B, th): return cdsS11(12 + 3*B, B, th)
def f4ratio(a, H):
    Hp = 6*H/(H - 6)
    ms = [np.sin(np.pi/H)*np.sin(2*np.pi/Hp), np.sin(3*np.pi/H)*np.sin(np.pi/Hp), np.sin(2*np.pi/H)*np.sin(2*np.pi/Hp), np.sin(3*np.pi/H)*np.sin(2*np.pi/Hp)]
    ms = sorted(ms); return ms[a]**2/ms[0]**2
analyse("f_4^(1)", f4, [2, 3, 4, 2, 1], 12, cdsS11_param, f4ratio, "3 if CDS (H = 12 + 3B), 9/2 if the soliton construction")
