"""Task 3 at 30 digits: K_{1->n}(mu) = a(mu_a) F(2mu) bt(mu_b) c(mu) (minimal solution; the sign s cancels) against
the worker's fit -K0 E_eta, with E_eta = (sinh th - i sin(pi A/T))/(sinh th + i sin(pi Abar/T)),
A = eps + (T-1)/2, Abar = A + omega l - T, K0 = (T/2)(T/2-1/2)(d)(1/2-d)(1/2+d)(T-1/2-d), d = (T-omega)/2.
Normalization of the singlets: coefficient 1 on (hw of left factor) (x) (lw of right factor) [mine] and, for
comparison, on index 0 (x) index 0 [worker].  Usage: python3 -u t3_amp.py n"""
import sys
import mpmath as mp
from sp_mp import SP, blk
from t3_breather import proj_mp, an_vector_mp, singlet_mp

mp.mp.dps = 30
n = int(sys.argv[1])
om = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.37')


def K0(mu, T, om):
    d = (T - om)/2
    h = mp.mpf(1)/2
    return (blk(T/2, mu, T)*blk(T/2 - h, mu, T)*blk(d, mu, T)*blk(h - d, mu, T)*blk(h + d, mu, T)
            * blk(T - h - d, mu, T))


def Eeta(mu, S, bar=False):
    T = S.T
    th = 1j*mp.pi*mu/T
    A = S.eps + (T - 1)/2
    Ab = A + S.om*S.l - T
    if bar:
        A, Ab = Ab, A
    return (mp.sinh(th) - 1j*mp.sin(mp.pi*A/T))/(mp.sinh(th) + 1j*mp.sin(mp.pi*Ab/T))


def worker_conv(S):
    """y_0 of the singlet in Vbar(1) (x) V(e^{i omega u}) normalized on Vbar_hw (x) V_lw = vbar_n (x) v_n.
    This is the out-singlet of B -> Bbar and the in-singlet of Bbar -> B (theta-independent).  With the worker's
    normalization (index 0 (x) 0 everywhere): K_w = K * y_0 for B -> Bbar and Kbar_w = Kbar / y_0."""
    A = an_vector_mp(S.n, S.q, mp.mpf(1), True)
    B = an_vector_mp(S.n, S.q, mp.exp(S.om*1j*mp.pi*(1 - 1/S.T)), False)
    y, _ = singlet_mp(A, B, A['hw'])
    return y[0]


for l in range(n + 2):
    for eta in ('0.41+0.23j', '-0.3+1.1j'):
        S = SP(n, om, l, mp.mpmathify(eta))
        T = S.T
        h = mp.mpf(1)/2
        conv = worker_conv(S)
        rows = []
        for mu in (mp.mpc('0.11', '0.19'), mp.mpc('-0.23', '0.07'), mp.mpc('0.05', '-0.3')):
            ma, mb = mu - (T - 1)/2, mu + (T - 1)/2
            c, sp1, _ = proj_mp(S, mu)
            cb, sp2, _ = proj_mp(S, mu, bar=True)
            K = S.a(ma)*S.F(2*mu)*S.bt(mb)*c
            Kb = S.bt(ma)*S.F(2*mu)*S.a(mb)*cb
            # unitarity K_{1->n}(mu) K_{n->1}(-mu) = 1
            c2, _, _ = proj_mp(S, -mu, bar=True)
            Kb_m = S.bt(-mu - (T - 1)/2)*S.F(-2*mu)*S.a(-mu + (T - 1)/2)*c2
            r = K/(-K0(mu, T, om)*Eeta(mu, S))
            rb = Kb/(-K0(mu, T, om)*Eeta(mu, S, bar=True))
            rows.append((r, rb, K*Kb_m, r*conv, rb/conv))
        print(f'n={n} om={om} l={l} eta={eta}:')
        for r, rb, U, rw, rbw in rows:
            print(f'   K/(-K0 E) = {mp.nstr(r, 15)}   Kbar/(-K0 Ebar) = {mp.nstr(rb, 15)}   unitarity-1 = '
                  f'{mp.nstr(U - 1, 2)}   [worker normalization: {mp.nstr(rw, 10)}, {mp.nstr(rbw, 10)}]', flush=True)
