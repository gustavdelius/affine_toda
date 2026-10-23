"""delta H (in units of beta^2, foundations normalization, longest root length^2 = 2) extracted from each
one-loop soliton mass ratio, for every folded family, with the counting rule and with the single count."""
import numpy as np
from fold import folded
from analysis import soliton_data, dspin, drev, ds1, drho, arefl, D4_Z3, E6_Z2, E6_Z3, E7_Z2
pi = np.pi

def fam(kind, n=None):
    if kind == 'sin':            # c_n / d_{n+1}^(2) family: m_a = sin(a pi/H)
        return lambda H: np.array([np.sin(a*pi/H) for a in range(1, n+1)])
    if kind == 'bn':             # b_n / a_{2n-1}^(2) family: 2 sin(a pi/H) (a<n), 1
        return lambda H: np.array([2*np.sin(a*pi/H) for a in range(1, n)] + [1.0])
    if kind == 'g2':
        return lambda H: np.array([1.0, 2*np.cos(pi/H)])
    if kind == 'f4':             # CDS (f_4^(1), e_6^(2)): 1, 2cos(pi/H+pi/6), 2cos(pi/H), 1+2cos(2pi/H)
        return lambda H: np.array([1.0, 2*np.cos(pi/H + pi/6), 2*np.cos(pi/H), 1 + 2*np.cos(2*pi/H)])

def dH(F, massfun, H0, f=1.0):
    S = soliton_data(F)
    Ms = np.array([s['Mcl'] for s in S]); o = np.argsort(Ms); S = [S[i] for i in o]; Ms = Ms[o]
    m0 = massfun(H0); om = np.argsort(m0); m0s = m0[om]
    assert np.allclose(Ms/Ms[0], m0s/m0s[0], atol=1e-6), (Ms/Ms[0], m0s/m0s[0])
    eps = 1e-6
    lr = lambda H: np.log(np.sort(massfun(H))/np.sort(massfun(H))[0]) if False else np.log(massfun(H)[om]/massfun(H)[om][0])
    slope = (lr(H0 + eps) - lr(H0 - eps))/(2*eps)
    out = {}
    for key in ('r', 'rMW'):
        d = np.array([S[i][key] - S[0][key] for i in range(len(S))])
        out[key] = [d[i]/slope[i] for i in range(1, len(S)) if abs(slope[i]) > 1e-9]
    return out

rows = []
def row(name, F, massfun, H0, ref, refval):
    o = dH(F, massfun, H0)
    rows.append((name, H0, o['r'], o['rMW'], ref, refval))
    print(f"{name:10s} H0={H0:5g}: counting rule dH/beta^2 = {np.round(o['r'], 6)};  single count = {np.round(o['rMW'], 6)};  exact: {ref} = {refval:.6f}")

if __name__ == "__main__":
    for n in [3, 4, 5]:
        row(f"c_{n}^(1)", folded('a', 2*n-1, arefl(2*n-1), "", verbose=False), fam('bn', n), 2*n, "-1/(4pi) [DG; companion]", -1/(4*pi))
    for n in [3, 4, 5, 6, 7]:
        row(f"b_{n}^(1)", folded('d', n+1, dspin(n+1), "", verbose=False), fam('sin', n), 2*n, "+1/(4pi) [GMW]", 1/(4*pi))
    row("g_2^(1)", folded('d', 4, D4_Z3, "", verbose=False), fam('g2'), 6, "+1/(2pi) [Takacs: 6+2/w]", 1/(2*pi))
    row("f_4^(1)", folded('e', 6, E6_Z2, "", verbose=False), fam('f4'), 12, "breather CDS: -3/(4pi)", -3/(4*pi))
    for n in [2, 3, 4]:
        row(f"a_{2*n-1}^(2)", folded('d', 2*n, drev(2*n), "", verbose=False), fam('bn', n), 2*n-1, "-1/(8pi) [companion, beta_c^2=beta^2/2]", -1/(8*pi))
    for n in [3, 4, 5, 6]:
        row(f"d_{n+1}^(2)", folded('d', n+2, ds1(n+2), "", verbose=False), fam('sin', n), 2*n+2, "+1/(2pi) [GM, beta_GM^2=beta^2/2]", 1/(2*pi))
    row("d_4^(3)", folded('e', 6, E6_Z3, "", verbose=False), fam('g2'), 12, "+3/(2pi) [Takacs 12+6/w, beta_c^2=beta^2/3]", 3/(2*pi))
    row("e_6^(2)", folded('e', 7, E7_Z2, "", verbose=False), fam('f4'), 18, "+3/(2pi) [companion 18+6/w, beta_c^2=beta^2/2]", 3/(2*pi))
