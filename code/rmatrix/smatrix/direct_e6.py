import numpy as np
exec(open('direct_lagrangian.py').read().split("# ---- f_4^(1)")[0])
def e6_ours(w):                                    # soliton-derived e_6^(2) breather amplitude = CDS S_11 at H = 18 + 6/w, T = 9w + 3
    T = 9*w + 3; H = 18 + 6/w; Bv = (H - 12)/3
    def brk(x): return [(x - 1, 1), (x + 1, 1), (x - 1 + Bv, -1), (x + 1 - Bv, -1)]
    ys = []
    for x in (1, H/3 + 1, H - 1, H - H/3 - 1):
        for z, s in brk(x): ys.append(z*T/H if s > 0 else -z*T/H)
    return ys, T, 1
def e6_ratio(a):
    def r(w):
        H = 18 + 6/w; Hp = 6*H/(H - 6)
        ms = sorted([np.sin(np.pi/H)*np.sin(2*np.pi/Hp), np.sin(3*np.pi/H)*np.sin(np.pi/Hp), np.sin(2*np.pi/H)*np.sin(2*np.pi/Hp), np.sin(3*np.pi/H)*np.sin(2*np.pi/Hp)])
        return ms[a]**2/ms[0]**2
    return r
E4 = np.eye(4)
f4 = [E4[1]-E4[2], E4[2]-E4[3], E4[3], 0.5*(E4[0]-E4[1]-E4[2]-E4[3]), -(E4[0]+E4[1])]
e62 = [2*np.array(a)/np.dot(a, a) for a in f4]          # e_6^(2) roots = coroots of f_4^(1) (long length^2 2, short 4)
ths = [0.25, 0.5, 0.8, 1.1, 1.5, 2.0, 2.6, 3.3]
run("e_6^(2)", e62, [2, 3, 2, 1, 1], e6_ours, {1: e6_ratio(1), 2: e6_ratio(2), 3: e6_ratio(3)}, ths)
print("gradation assumed for e_6^(2): w = 2 pi/beta^2 - 1, i.e. 1/w ~ beta^2/(2 pi)")
