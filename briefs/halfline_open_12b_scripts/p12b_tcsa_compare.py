"""Part (b): comparison with the numbers/lists that BPT (hep-th/0108157) actually quote.
(1) BPT Table 5.2: first Dirichlet excited state |0> in both boundaries: E1-E0 = 2 M cos((eta - pi/2)/lambda), eta = (lambda+1) pi x, x = beta phi_0/(2 pi).
(2) BPT sec. 4.3.2 (lambda = 7, eta = 0.9 eta_N, vartheta = 0): sequence nu_0 > w_1 > nu_1 > w_2 > nu_2 > w_3 > nu_3 and the low-lying states in the even/odd
    (C-parity = parity of the sum of labels) sectors, eqs. (4.7), (4.8), against the spectrum rule of prp-bsg-spectrum."""
import math
from p12b_spectrum_poles import states
def report(name, ok): print(('PASS ' if ok else 'FAIL ') + name)
tab = [(31, 0.495, 0.032428), (17, 0.485, 0.099750), (7, 0.48, 0.14349), (7, 0.45, 0.35711), (41/8, 0.44, 0.44675), (41/8, 0.36, 1.0035), (17/8, 0.4, 0.89148)]
ok = True
for lam, x, pred in tab:
    eta = (lam+1)*math.pi*x
    E = 2*math.cos((eta - math.pi/2)/lam)
    print(f'lambda={lam:7.4f} x={x}: 2cos(...) = {E:.5f}, BPT predicted {pred}')
    ok &= abs(E - pred) < 6e-5*max(1, pred)
report('BPT Table 5.2 "predicted" column reproduced by E_1-E_0 = 2 M cos((eta-pi/2)/lambda)', ok)

lam = 7.0; eta = 0.9*math.pi*(lam+1)/2
sts, nu_, w_ = states(eta, lam, 3)
order = sorted([(nu_(j), f'nu_{j}') for j in range(0, 4)] + [(w_(j), f'w_{j}') for j in range(1, 4)], reverse=True)
print('sequence in the strip:', ' > '.join(n for _, n in order))
seq_ok = [n for _, n in order] == ['nu_0', 'w_1', 'nu_1', 'w_2', 'nu_2', 'w_3', 'nu_3']
report('BPT (4.6): nu_0 > w_1 > nu_1 > w_2 > nu_2 > w_3 > nu_3', seq_ok)
def E(S): return sum(math.cos(nu_(n) if i % 2 == 0 else w_(n)) for i, n in enumerate(S))
even = sorted([(E(S), S) for S, _ in sts if S and sum(S) % 2 == 0]); odd = sorted([(E(S), S) for S, _ in sts if S and sum(S) % 2 == 1])
print('even sector (E/M, labels):', [(round(e, 3), ','.join(map(str, S))) for e, S in even])
print('odd  sector (E/M, labels):', [(round(e, 3), ','.join(map(str, S))) for e, S in odd])
bpt_even = {(0,), (2,), (0, 2), (1, 3), (0, 1, 1), (0, 1, 3), (1, 2, 3), (2, 3, 3)}
bpt_odd = {(1,), (3,), (1, 2), (2, 3), (0, 1, 2), (1, 2, 2)}
ev = {S for _, S in even}; od = {S for _, S in odd}
report('every state listed in BPT (4.7) (even) and (4.8) (odd) is allowed by the rule and in the right C-parity sector', bpt_even <= ev and bpt_odd <= od)
print('allowed (<= 3 labels) but not in the BPT lists:', sorted(ev - bpt_even), sorted(od - bpt_odd))
report('|1,1> is excluded (BPT: no TCSA line, explained by a Coleman-Thun diagram because w_1 > nu_1)', (1, 1) not in ev | od)

# (3) Zamolodchikov's boundary energy, BPT (3.4):  E(eta,vt) = -M/(2cos(pi/2lam)) [cos(eta/lam) + cosh(vt/lam) - cos(pi/2lam)/2 + sin(pi/2lam)/2 - 1/2].
# BPTT's conjecture m_{|0>}(eta,vt) = m_{|>}(etabar,vt) is then EXACTLY E(etabar,vt) - E(eta,vt) = M cos(nu_0) with nu_0 = (eta - pi/2)/lam, for all eta, vt.
def Ebdy(eta, vt, lam):
    return -1/(2*math.cos(math.pi/(2*lam)))*(math.cos(eta/lam) + math.cosh(vt/lam) - 0.5*math.cos(math.pi/(2*lam)) + 0.5*math.sin(math.pi/(2*lam)) - 0.5)
dev = 0
for lam, eta, vt in ((3.3, 5.5, 0.3), (7.0, 2.0, 0.0), (12.7, 3.1, 4.0), (1.7, 0.9, 1.2)):
    etab = math.pi*(lam + 1) - eta
    dev = max(dev, abs(Ebdy(etab, vt, lam) - Ebdy(eta, vt, lam) - math.cos((eta - math.pi/2)/lam)))
report(f'E(etabar,vt) - E(eta,vt) = M cos nu_0 exactly (Zamolodchikov boundary energy = BPTT m_|0>(eta) = m_|>(etabar)) [{dev:.0e}]', dev < 1e-13)
