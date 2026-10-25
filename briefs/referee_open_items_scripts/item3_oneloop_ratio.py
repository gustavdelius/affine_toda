"""Investigation and verification of Item 3: a_5^{(2)} in code/rmatrix/smatrix/oneloop_ratio.py."""
import sys
import numpy as np

# Run code/rmatrix/smatrix/oneloop_ratio.py and inspect its exact output and internals:
sys.path.insert(0, 'code/rmatrix/smatrix')
import oneloop_ratio as olr

print("\n--- Diagnosing a_5^(2) in oneloop_ratio.py ---")
s2 = 1/np.sqrt(2)
e = np.eye(3)
a5 = [s2*(e[0]-e[1]), s2*(e[1]-e[2]), s2*2*e[2], -s2*(e[0]+e[1])]
kac = [1, 2, 1, 1]
h = 5
m2, C3, C4 = olr.toda_data(a5, kac)
dm2 = olr.one_loop_dm2(m2, C3)
print("m2 (sorted):", m2)
print("m (sorted):", np.sqrt(m2))
print("Note: index 0 has m = 1.0 (which is particle n=3 of a_5^(2)!), while index 1 has m = 2*sin(pi/5) = 1.17557 (particle 1)!")

# Check dH from mass ratios:
dH = []
for a in range(1, len(m2)):
    dR = (dm2[a]*m2[0] - m2[a]*dm2[0])/m2[0]**2
    eps = 1e-6
    slope = (olr.a5ratio(a, h + eps) - olr.a5ratio(a, h - eps))/(2*eps)
    dH.append(dR/slope)
print("dH (one-loop shift of H per beta^2):", dH, "-> equals 1/(8*pi) =", 1/(8*np.pi))

# Why did tree_over_4m2sinh / g give the wrong B_tree?
# Because in oneloop_ratio.py:
# 1) tree_over_4m2sinh(m2, C3, C4, 0, th) evaluates the tree amplitude of index a=0 (which is particle n=3, mass 1.0),
# 2) while S11 was passed as dgz11(h + B/2, B, th) = brace(1, B, H, th)*brace(H-1, B, H, th) (which is the c_n^(1) S_11 amplitude with nu=0, NOT the a_5^(2) S_{33} or S_{11} amplitude!).
# Wait: what ARE the exact DGZ S-matrix amplitudes S_{33} (for particle n=3, index 0) and S_{11} (for particle 1, index 1) of a_5^(2)?
# In Ch. 25 (eq-a2n1-smatrix, 25-a2n-solitons.qmd line 50) and Ch. 14 (14-nonsimply-laced.qmd):
# For particle n=3 (index 0, the lightest particle of a_5^(2)):
#   S_{33}(theta) = prod_{p=1}^{n-1} {2p - 1}_{1/2} = {1}_{1/2} {3}_{1/2}
#   where {x}_{1/2} in oneloop_ratio.py's brace(x, B, H, th, nu) has nu = -0.5!
# Let's check nu = -0.5 vs nu = +0.5 and compare with tree_over_4m2sinh(m2, C3, C4, 0, th)!
ths = [0.4, 0.9, 1.5, 2.3]
eps = 1e-7
tr0 = [olr.tree_over_4m2sinh(m2, C3, C4, 0, th) for th in ths]
tr1 = [olr.tree_over_4m2sinh(m2, C3, C4, 1, th) for th in ths]

for nu in (-0.5, 0.0, 0.5):
    S33_fn = lambda h_val, B, th: olr.brace(1, B, h_val + B/2, th, nu=nu) * olr.brace(3, B, h_val + B/2, th, nu=nu)
    g33 = [np.log(S33_fn(h, eps, th)).imag / eps for th in ths]
    kap33 = np.array(tr0) / np.array(g33)
    print(f"Particle n=3 (index 0) with S_33 = {{1}}_{{nu}}{{3}}_{{nu}} (nu={nu:+.1f}):")
    print(f"  tr0 / g33 = {np.round(kap33, 6)}  (1/(4*pi) = {1/(4*np.pi):.6f})")
    print(f"  delta H / B_tree = {np.mean(dH)/np.mean(kap33):.6f}")

# And what about Particle 1 (index 1, m_1 = 2*sin(pi/5))?
# For particle 1 of a_{2n-1}^{(2)}, S_{11}(theta) = {1}_0 {2}_{1/2} {H-2}_{1/2} {H-1}_0 -> wait: let's check {1}_0 {H-1}_0!
for nu_pair in ((0.0, -0.5), (0.0, 0.5), (0.0, 0.0)):
    nu1, nu2 = nu_pair
    S11_a5 = lambda h_val, B, th: (olr.brace(1, B, h_val + B/2, th, nu=nu1) *
                                   olr.brace(h_val + B/2 - 1, B, h_val + B/2, th, nu=nu2))
    g11 = [np.log(S11_a5(h, eps, th)).imag / eps for th in ths]
    kap11 = np.array(tr1) / np.array(g11)
    print(f"Particle 1 (index 1) with brace(1, nu={nu1})*brace(H-1, nu={nu2}):")
    print(f"  tr1 / g11 = {np.round(kap11, 6)}")
    print(f"  delta H / B_tree = {np.mean(dH)/np.mean(kap11):.6f}")
    if nu_pair == (0.0, 0.0):
        ratio_fixed = np.mean(dH) / np.mean(kap11)
        spread_fixed = np.ptp(kap11)
        print("PASS a5_2_particle1_index_fix" if abs(ratio_fixed - 0.5) < 1e-7 and spread_fixed < 1e-8 else "FAIL a5_2_particle1_index_fix")


