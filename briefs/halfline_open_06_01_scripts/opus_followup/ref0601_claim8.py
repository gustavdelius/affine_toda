"""Referee check (2026-10-12) of claim 8 (06_01): at each singular endpoint chi_end of the stationary-pair family
(species a of a_n, tau_j = 1 + 2d cos(chi+2 pi j a/h) e^{m_a x} + cos^2(pi a/h) d^2 e^{2 m_a x}, d=1/(1+sin(pi a/h))),
find the exact (regular) patterns on both sides, the vanishing nodes, and the set of patterns sharing the endpoint
(C=-1 on vanishing nodes, any sign on the other nodes, all of which have e^{psi/2} -> 0 when they neighbour a vanishing node)."""
import numpy as np, itertools
def s(C): return ''.join('+' if c > 0 else '-' for c in C)
for n in (2, 3, 4, 5):
    h = n+1; AL = np.array([np.eye(h)[j]-np.eye(h)[(j+1) % h] for j in range(h)])
    for a in range(1, (h+1)//2):
        if 2*a == h: continue
        m = 2*np.sin(np.pi*a/h); d = 1/(1+np.sin(np.pi*a/h)); A = np.cos(np.pi*a/h)**2
        def tau(chi):
            j = np.arange(h); c = np.cos(chi+2*np.pi*j*a/h)
            return 1+2*d*c+A*d*d, 2*d*c*m+2*A*d*d*m
        def pat(chi):
            t, tp = tau(chi); eh = np.sqrt(np.roll(t, 1)*np.roll(t, -1))/t; phid = -(AL.T@(tp/t))
            return min((np.abs(phid+AL.T@(np.array(C)*eh)).max(), C) for C in itertools.product([1, -1], repeat=h))
        ends = sorted({round((np.pi-2*np.pi*j*a/h) % (2*np.pi), 12) for j in range(h)})
        for ce in ends:
            t0 = tau(ce)[0]; van = [j for j in range(h) if abs(t0[j]) < 1e-9]
            (r1, C1), (r2, C2) = pat(ce+0.03), pat(ce-0.03)
            nb = sorted({(v+1) % h for v in van} | {(v-1) % h for v in van})
            free = [j for j in range(h) if j not in van]
            ok_nb = set(free) <= set(nb)
            opp = lambda C: all(C[(v-1) % h]*C[(v+1) % h] == -1 for v in van)
            shared = [C for C in itertools.product([1, -1], repeat=h) if all(C[v] == -1 for v in van)]
            print(f'a_{n} a={a} chi_end={ce:.4f} vanishing {van} (all other nodes are neighbours: {ok_nb}); '
                  f'regular sides {s(C1)} (res {r1:.0e}), {s(C2)} (res {r2:.0e}); opposite-neighbour on both: {opp(C1)},{opp(C2)}; '
                  f'flipped nodes {[j for j in range(h) if C1[j] != C2[j]]}')
            prodp = [C for C in shared if np.prod(C) == 1]; prodm = [C for C in shared if np.prod(C) == -1]
            print(f'      patterns sharing the endpoint: {len(shared)} ({len(prodp)} with prod C=+1, {len(prodm)} with -1); '
                  f'opposite-neighbour ones: {[s(C) for C in shared if opp(C)]}; prodC=-1 ones all non-opposite: {all(not opp(C) for C in prodm)}')
        break_after = True
