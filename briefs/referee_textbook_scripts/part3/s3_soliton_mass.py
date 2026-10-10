"""eq-soliton-mass M_a = 2 h m_a / beta^2 for d_n, e_n single solitons (static, complex xi), by direct integration of
int (1/2 phi'.phi' + V_imag) dx (no complex conjugation), and the total-derivative formula
E = (2/beta^2) sum_j (2/alpha_j^2) [d_x ln tau_j]_{-inf}^{+inf} = (2/beta^2) mu sum_j (2/alpha_j^2) deg tau_j ."""
import numpy as np, mpmath as mp
from lie import Algebra, delta_mp
from hirota import soliton, MpCtx
mp.mp.dps = 30
def report(n, ok): print(('PASS ' if ok else 'FAIL ') + n)
beta = mp.mpf('0.9')
def energy(Al, a):
    delta, lam, res = delta_mp(Al, a, dps=40)
    mp.mp.dps = 30
    mu = mp.sqrt(lam)
    t, tail = soliton(Al, delta, mu, ctx=MpCtx(30), Kmax=14, tol=mp.mpf(10)**-22)
    degs = [len(c)-1 for c in t]
    r1 = Al.r + 1; alpha = [[mp.mpf(float(x)) for x in row] for row in Al.alpha]
    # exact-ish roots: refine alpha via Cholesky in mp
    Cf = mp.matrix([[int(round(Al.Cfin[i, j])) for j in range(Al.r)] for i in range(Al.r)])
    L = mp.cholesky(Cf)
    fin = [[L[i, k] for k in range(Al.r)] for i in range(Al.r)]
    marks = [int(round(x)) for x in Al.n[1:]]
    th = [sum(marks[i]*fin[i][k] for i in range(Al.r)) for k in range(Al.r)]
    alpha = [[-x for x in th]] + fin
    best = None
    for im in np.linspace(0.05, 3.1, 40):           # choose Im xi with no zero of tau_j on the real line
        xi = mp.mpc(0, im)
        xs = np.linspace(-12, 12, 241)/float(mu)
        mn = min(abs(sum(c*mp.e**(p*(mu*x+xi)) for p, c in enumerate(cs))) for x in xs for cs in t)
        if best is None or mn > best[0]: best = (mn, xi)
    xi = best[1]
    nbr = Al.nbrs; n = [int(round(x)) for x in Al.n]
    def dens(x):
        E = mp.e**(mu*x + xi)
        T = [sum(c*E**p for p, c in enumerate(cs)) for cs in t]
        dT = [sum(c*p*mu*E**p for p, c in enumerate(cs)) for cs in t]
        dP = [1j/beta*sum(alpha[j][k]*dT[j]/T[j] for j in range(r1)) for k in range(Al.r)]
        V = -(1/beta**2)*sum(n[j]*(mp.fprod([T[l] for l in nbr[j]])/T[j]**2 - 1) for j in range(r1))
        return sum(d*d for d in dP)/2 + V
    Lx = 50/mu
    E = mp.quad(dens, mp.linspace(-Lx, Lx, 21))
    pred = 2*Al.h*mu/beta**2
    formula = (2/beta**2)*mu*sum(degs)            # (2/alpha_j^2)=1 simply laced
    return E, pred, formula, degs, n, float(best[0]), xi
for typ, r in [('d', 4), ('d', 5), ('e', 6)]:
    Al = Algebra(typ, r)
    for a in sorted(Al.mass):
        E, pred, formula, degs, n, mn, xi = energy(Al, a)
        ok = abs(E - pred) < 1e-15*abs(pred) and degs == n
        print(f'   {typ}{r} species {a}: m_a={float(Al.mass[a]):.6f}, Im xi={float(mp.im(xi)):.3f}, min|tau|={mn:.2e}, deg tau_j={degs} (n_j={n}), E={mp.nstr(E, 16)}, 2h m_a/beta^2={mp.nstr(pred, 16)}')
        report(f'{typ}{r} species {a}: E = 2 h m_a/beta^2 and deg tau_j = n_j', ok)
