import numpy as np
exec(open('f4gamma.py').read().split("Ta = (6, F(9, 2)); s =")[0])
def W_of(A, extra, T, th):
    FA = Fmaker(A, T); FE = Fmaker(A + extra, T); bet = np.pi*(T - 1)/T
    z = [FE(th + 1j*bet)/FA(th + 1j*bet), FE(th)/FA(th), FE(th - 1j*bet)/FA(th - 1j*bet)]
    return z[0]*z[1]**2*z[2], z[1]
for label, Tf, Aex, ours, cds in (("(a)", (6, F(9, 2)), [(1, 1), (3, F(5, 2)), (4, 3), (6, F(9, 2))], None, None),):
    T = val(Tf, w0); A = [val(a, w0) for a in Aex]
    print(f"{label}: T = {T:.3f}, zeros {np.round(A, 3)}")
    for e in (w0 + 0.5, 3*w0 + 2.0, 2*w0 + 2.0, 0.5):
        th = 0.31 + 0.17j
        W, Z = W_of(A, [e, T - e], T, th)
        print(f"   single pair e = {e:.3f}:  Z_E(theta) = {Z:.4f}   W_E(theta) = {W:.6f}")
    for sgn_, name in ((-1, "a_j - 1/2"), (+1, "a_j + 1/2")):
        extra = []
        for a in A[:-1]:
            e = a + sgn_*0.5; extra += [e, T - e]
        e = A[-1] + sgn_*0.5; extra += [e, T - e]
        Ws = [W_of(A, extra, T, th)[0] for th in ths]
        print(f"   full doubling (copies at {name}, with crossing partners): W_E at 4 rapidities = {np.round(Ws, 5)}")
