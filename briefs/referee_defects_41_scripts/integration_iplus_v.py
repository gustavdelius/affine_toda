"""Type I_+, V: an invertible intertwiner for every r from the Chevalley involution omega (e <-> f, k -> k^-1), which is a
coalgebra anti-automorphism.  Claim: on the sector Q of V(y) (x) W^+ (W^+: e_j|n> = M q^{n_{j+1}}|n + alpha_j>, offset s),
T_Q = (M'_{-Q})^{-1}, where M' is the two-entry intertwiner of the reversed-orientation b^- module (f_j|n'> = M q^{-n'_{j+1}}|n' - alpha_j>,
offset -s) for Vbar(1/y):  M'(vb_i|n'>) = q^{n'_i}|n'>vb_i + X' q^{n'_i - n'_{i+1}}|n' - alpha_i>vb_{i+1}, X' = -(q - 1/q) M/y.
Checked exactly (Fractions) against the I_+ sector equations of V (x) W^+, h = 3..6."""
from fractions import Fraction as Fr
from integration_uniqueness import Weyl, particle, sector_map, mat_mul, mat_inv, qpow
q = Fr(2); M = Fr(3, 7); y = Fr(5, 4)
for h in (3, 4, 5, 6):
    t = (Fr(3, 2), Fr(5, 7)) + (Fr(1),)*(h - 3) + (Fr(14, 15),)
    tinv = tuple(1/x for x in t)
    W = Weyl(h, q, M, t, +1); P = particle(h, q, y, False)
    Xp = -(q - 1/q)*M/y
    def TQ(Q):
        Qm = tuple(-c for c in Q)
        Mp = [[Fr(0)]*h for _ in range(h)]
        for i in range(h):
            n = tuple(Qm[c] + (1 if c == i else 0) for c in range(h))     # Vbar weight -e_i: n' = -Q + e_i
            qn = lambda c: tinv[c]*qpow(q, n[c])
            Mp[i][i] = qn(i); Mp[(i + 1) % h][i] = Xp*qn(i)/qn((i + 1) % h)
        return mat_inv(Mp)
    ok = True; full = True
    for Q in (tuple([1] + [0]*(h - 1)), tuple([2, -1] + [0]*(h - 2)), tuple([0, 0, 1] + [0]*(h - 3))):
        for j in range(h):
            A, Qn = sector_map(P, W, j, Q, True); B, _ = sector_map(P, W, j, Q, False)
            ok &= mat_mul(TQ(Qn), A) == mat_mul(B, TQ(Q))
        full &= all(x != 0 for row in TQ(Q) for x in row)
    print(f'h={h}: T_Q = (M\'_-Q)^-1 satisfies the type-I_+ equations for V exactly: {ok}; all h^2 entries nonzero: {full}')
