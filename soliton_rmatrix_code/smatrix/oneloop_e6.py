import numpy as np
exec(open('oneloop_ratio.py').read().split("# c_3^(1): long roots length^2 2")[0])
E4 = np.eye(4)
f4 = [E4[1]-E4[2], E4[2]-E4[3], E4[3], 0.5*(E4[0]-E4[1]-E4[2]-E4[3]), -(E4[0]+E4[1])]
e62 = [2*np.array(a)/np.dot(a, a) for a in f4]              # e_6^(2) Toda roots = coroots of f_4^(1)
kac_e6 = [1, 2, 3, 2, 1]                                     # order: alpha_1..alpha_4, alpha_0 coroots
# find the Kac labels n with sum n_i alpha_i^vee = 0
A = np.array(e62); ns = np.linalg.svd(A.T)[2][-1]; ns = ns/ns[np.argmin(np.abs(ns))]; ns = ns/np.min(np.abs(ns))
print("e_6^(2) Kac labels (on alpha_1..alpha_4, alpha_0):", np.round(ns, 6))
def cdsS11(H, B, th):
    return brace(1, B, H, th)*brace(H - 1, B, H, th)*brace(H/3 + 1, B, H, th)*brace(H - H/3 - 1, B, H, th)
# near the e_6^(2) end: B = 2 - Bt, H = 18 - 3 Bt ; log S ~ Bt * gt(theta)
S11_e6 = lambda h, Bt, th: cdsS11(18 - 3*Bt, 2 - Bt, th)
def f4ratio(a, H):
    Hp = 6*H/(H - 6)
    ms = [np.sin(np.pi/H)*np.sin(2*np.pi/Hp), np.sin(3*np.pi/H)*np.sin(np.pi/Hp), np.sin(2*np.pi/H)*np.sin(2*np.pi/Hp), np.sin(3*np.pi/H)*np.sin(2*np.pi/Hp)]
    ms = sorted(ms); return ms[a]**2/ms[0]**2
analyse("e_6^(2)", e62, list(np.round(ns).astype(int)), 18, S11_e6, f4ratio, "-3 if CDS (H = 18 - 3 Bt)")
