"""Referee check (2026-10-12): states with <=3 labels allowed by the BPTT rule nu_{n1}>w_{n2}>nu_{n3}, all in (0,pi/2),
at BPT's point lambda=7, eta=0.9 eta_N, vartheta=0; C-parity sector = (-1)^(sum of labels) (BPT sec. 4.3.2).
Compared with BPT (hep-th/0108157) eqs. (4.7), (4.8) (uv_ir3.tex l.1113-1120) and the Fig. 4.3 captions."""
import math, itertools
lam = 7.0; etaN = math.pi*(lam+1)/2; eta = 0.9*etaN
nu = lambda n: eta/lam - (2*n+1)*math.pi/(2*lam)
w = lambda n: math.pi - eta/lam - (2*n-1)*math.pi/(2*lam)
print('nu_n/pi:', [round(nu(n)/math.pi, 4) for n in range(-1, 5)], ' w_n/pi:', [round(w(n)/math.pi, 4) for n in range(0, 5)])
def ok(S):
    p = [nu(n) if i % 2 == 0 else w(n) for i, n in enumerate(S)]
    return all(0 < x < math.pi/2 for x in p) and all(p[i] > p[i+1] for i in range(len(p)-1))
allowed = [S for k in (1, 2, 3) for S in itertools.product(range(0, 6), repeat=k) if ok(S)]
even = [S for S in allowed if sum(S) % 2 == 0]; odd = [S for S in allowed if sum(S) % 2 == 1]
E = lambda S: sum(math.cos(nu(n)) if i % 2 == 0 else math.cos(w(n)) for i, n in enumerate(S))
bpt47 = [(0,), (2,), (0, 2), (1, 3), (0, 1, 1), (0, 1, 3), (1, 2, 3), (2, 3, 3)]
bpt48 = [(1,), (3,), (1, 2), (2, 3), (0, 1, 2), (1, 2, 2)]
cap_even = [(0,), (2,), (0, 2), (0, 1, 1)]; cap_odd = [(1,), (3,), (0, 1), (0, 3), (1, 2), (0, 1, 2)]
for name, mine, bpt, cap in (('even', even, bpt47, cap_even), ('odd', odd, bpt48, cap_odd)):
    print(f'{name}: allowed {len(mine)}:', ' '.join(f'{S}[E/M={E(S):.3f}]' for S in sorted(mine, key=E)))
    print(f'   in BPT list but not allowed/wrong sector: {[S for S in bpt if S not in mine]}')
    print(f'   allowed but missing from BPT list: {[S for S in mine if S not in bpt]}')
    print(f'   caption (TCSA-identified) states not in BPT list: {[S for S in cap if S not in bpt]};  caption states not allowed: {[S for S in cap if S not in mine]}')
print('|1,1> allowed?', ok((1, 1)), ' (w_1>nu_1:', w(1) > nu(1), ')')
