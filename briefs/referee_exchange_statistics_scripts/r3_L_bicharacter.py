# Extract L with Rcheck(-q) = L Rcheck(q) L (L diagonal +-1), and test whether L is, up to a function of the total weight
# and a symmetric one-soliton gauge d(a)d(b), a bicharacter (-1)^{b_a^T M b_b} in the bit labels b (1 = minus sign).
import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sectors', 'bloch'))
os.environ['QSIGN'] = '1'
from spin_common import *
def Lmatrix(n, alg, om):
    D, Wn, setup = make(n, alg)
    q = np.exp(-1j*np.pi*om); Rp, Rm = setup(q), setup(-q)
    P = flip(D, D)
    nzall = np.zeros((D*D, D*D), bool); L = {}
    for x in (1.7+0.4j, 0.3-1.1j, 2.9):
        Cp, Cm = P@Rp(x), P@Rm(x)        # Rcheck = P R
        nz = np.abs(Cp) > 1e-10; assert np.array_equal(nz, np.abs(Cm) > 1e-10)
        r = np.where(nz, Cm/np.where(nz, Cp, 1), 0)
        assert np.allclose(np.abs(r[nz]), 1, atol=1e-9) and np.allclose(r[nz].imag, 0, atol=1e-9)
        nzall |= nz
    sgn = np.sign(r.real).astype(int)
    for seed in range(D*D):
        if seed in L: continue
        L[seed] = 1; st = [seed]
        while st:
            i = st.pop()
            for j in np.nonzero(nzall[i] | nzall[:, i])[0]:
                s = sgn[i, j] if nzall[i, j] else sgn[j, i]
                if j not in L: L[j] = s*L[i]; st.append(j)
                elif L[i]*L[j] != s: raise Exception('not diagonal similarity')
    bits = np.array([[1 if w < 0 else 0 for w in Wn[i]] for i in range(D)])
    return D, Wn, bits, np.array([L[i] for i in range(D*D)]).reshape(D, D)
for n, alg in ((2, 'c'), (3, 'c'), (2, 'a'), (3, 'a')):
    D, Wn, bits, L = Lmatrix(n, alg, 1.61)
    tot = {}
    sols = []
    for Mf in itertools.product((0, 1), repeat=n*n):
        M = np.array(Mf).reshape(n, n)
        for v in itertools.product((0, 1), repeat=n):
            v = np.array(v); ok = True; tot = {}
            for a in range(D):
                for b in range(D):
                    val = L[a, b]*(-1)**(bits[a]@M@bits[b] + v@bits[a] + v@bits[b])
                    key = tuple(bits[a] + bits[b])
                    if tot.setdefault(key, val) != val: ok = False; break
                if not ok: break
            if ok: sols.append((M, v))
    print(f"{alg}_{n}: L symmetric-pairs? {np.array_equal(L, L.T)};  bicharacter fits found: {len(sols)}")
    for M, v in sols[:4]:
        A = (M - M.T) % 2
        print('   M =', M.tolist(), ' v =', v.tolist(), '  M+M^T mod 2 =', ((M + M.T) % 2).tolist())
