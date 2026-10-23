import numpy as np
from fractions import Fraction as Fr
exec(open('oneloop_ratio.py').read().split("def blk(x, H, th)")[0])
E4 = np.eye(4)
f4 = [E4[1]-E4[2], E4[2]-E4[3], E4[3], 0.5*(E4[0]-E4[1]-E4[2]-E4[3]), -(E4[0]+E4[1])]
m2, C3, C4 = toda_data(f4, [2, 3, 4, 2, 1]); dm2 = one_loop_dm2(m2, C3)
dR = {a: (dm2[a]*m2[0] - m2[a]*dm2[0])/m2[0]**2 for a in (1, 2)}         # per beta_r^2 (Lagrangian)
ths = [0.25, 0.5, 0.8, 1.1, 1.5, 2.0, 2.6, 3.3]
Atree = np.array([tree_over_4m2sinh(m2, C3, C4, 0, th) for th in ths])
def val(y, w): return float(y[0])*w + float(y[1])
def bootstrap_BB(SB, sgn, Tf, s):
    """S[B,B](theta) = S[B,27](theta + i s) S[B,27](theta - i s): symmetric shift of blocks"""
    return [(y[0] + s[0], y[1] + s[1]) for y in SB] + [(y[0] - s[0], y[1] - s[1]) for y in SB], sgn*sgn
def weak(blocks, sgn, Tf, th):
    out = []
    for w in (2000.0, 4000.0):
        T = val(Tf, w)
        v = sgn*np.prod([np.sinh(th/2 + 1j*np.pi*val(y, w)/(2*T))/np.sinh(th/2 - 1j*np.pi*val(y, w)/(2*T)) for y in blocks])
        out.append(w*np.log(v).imag)
    return 2*out[1] - out[0]
def test(name, SB, sgn, Tf, a378, a79):
    s = ((Tf[0])/2, (Tf[1] - 1)/2)                      # tB/2 with tB = T - 1
    BB, sg = bootstrap_BB(SB, sgn, Tf, s)
    f = np.array([weak(BB, sg, Tf, th) for th in ths])
    kap = Atree/f; k = np.mean(kap)
    print(f"{name}: T = {Tf[0]}w + {Tf[1]}")
    print(f"   tree shape: Lagrangian/amplitude ratio at 8 rapidities = {np.round(kap, 7)}  (spread {np.ptp(kap):.1e})")
    print(f"   tree coupling map: 1/w = beta^2 / {1/abs(k)/np.pi:.4f} pi      (gradation implies beta^2 / 4 pi)")
    for a, (num, lab) in {2: (a378, 'm_3'), 1: (a79, 'm_2')}.items():
        r = lambda w: (2*np.cos(np.pi*val(num, w)/(2*val(Tf, w))))**2
        slope = (r(1e7) - r(5e6))/(1e-7 - 2e-7)
        pred_tree = slope*k; pred_grad = slope*(-1/(4*np.pi))
        print(f"   one-loop ({lab}/m_1)^2:  amplitude/Lagrangian = {pred_tree/dR[a]:.6f} (tree map),  {pred_grad/dR[a]:.6f} (gradation map)")
F = Fr
# (a)  T = 6w + 9/2, zeros w+1, 3w+5/2, 4w+3  (S[B,27] identified with J = 4000)
SBa = [(0, F(3, 4)), (1, F(5, 4)), (3, F(7, 4)), (3, F(11, 4)), (5, F(13, 4)), (6, F(15, 4)), (-2, F(-3, 4)), (-4, F(-15, 4))]
test("(a)", SBa, 1, (6, F(9, 2)), (1, 1), (3, F(5, 2)))
# (b)  T = 6w - 3/2, zeros w, 3w-1/2, 4w-1: convert the stored S[B,B] blocks to exact forms
BBn, sgn_b, Tb, Ab = np.load("f4BB_0,-0.5,-1,-1.5.npy", allow_pickle=True)
def exact(y, w=2.37):
    for a in range(-24, 25):
        b = y - a*w
        if abs(b*4 - round(b*4)) < 1e-6 and abs(b) < 30: return (a, F(round(b*4), 4))
    raise ValueError(y)
BBb = [exact(y) for y in BBn]
# test() takes S[B,27]; for (b) feed S[B,B] directly by bypassing the bootstrap step
def test_BB(name, BB, sg, Tf, a378, a79):
    f = np.array([weak(BB, sg, Tf, th) for th in ths]); kap = Atree/f; k = np.mean(kap)
    print(f"{name}: T = {Tf[0]}w + {Tf[1]}")
    print(f"   tree shape: Lagrangian/amplitude ratio at 8 rapidities = {np.round(kap, 7)}  (spread {np.ptp(kap):.1e})")
    print(f"   tree coupling map: 1/w = beta^2 / {1/abs(k)/np.pi:.4f} pi      (gradation implies beta^2 / 4 pi)")
    for a, (num, lab) in {2: (a378, 'm_3'), 1: (a79, 'm_2')}.items():
        r = lambda w: (2*np.cos(np.pi*val(num, w)/(2*val(Tf, w))))**2
        slope = (r(1e7) - r(5e6))/(1e-7 - 2e-7)
        print(f"   one-loop ({lab}/m_1)^2:  amplitude/Lagrangian = {slope*k/dR[a]:.6f} (tree map),  {slope*(-1/(4*np.pi))/dR[a]:.6f} (gradation map)")
test_BB("(b)", BBb, int(sgn_b), (6, F(-3, 2)), (1, 0), (3, F(-1, 2)))
# reference: the original choice T = 6w + 3/2
BBo, sgo, To, Ao = np.load("f4BB_1,1.5,1,1.5.npy", allow_pickle=True)
test_BB("(original)", [exact(y) for y in BBo], int(sgo), (6, F(3, 2)), (1, 1), (3, F(3, 2)))
