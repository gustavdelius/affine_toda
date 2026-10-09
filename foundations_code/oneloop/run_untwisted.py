import numpy as np
from fold import folded, describe
pi = np.pi
# b_n^(1) from d_{n+1}^(1): sigma = spinor swap (r-1 r)
for r in [4, 5, 6, 7, 8]:
    perm = list(range(r+1)); perm[r-1], perm[r] = r, r-1
    F = folded('d', r, perm, f"b_{r-1}^(1)")
    describe(F)
# g_2^(1) from d_4^(1): triality (1 3 4)
perm = [0, 3, 2, 4, 1]
F = folded('d', 4, perm, "g_2^(1)"); describe(F)
# f_4^(1) from e_6^(1): (1 6)(3 5)
perm = [0, 6, 2, 5, 4, 3, 1]
F = folded('e', 6, perm, "f_4^(1)"); describe(F)
