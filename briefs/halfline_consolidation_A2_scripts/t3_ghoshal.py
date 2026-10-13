"""n = 1: the SP breather amplitude (30 digits) times sigma_A(mu_a) sigma_A(mu_b) against Ghoshal's K_B1
(eq-bsg-b1-reflection) at Dirichlet (vartheta -> infinity, eta_B = xi = pi eps - pi/2), lambda = omega.
Also: for n = 1 the in- and out-singlets (normalized hw (x) lw) are the same vector under J: Vbar ~ V."""
import mpmath as mp
from sp_mp import SP, blk
from t3_breather import proj_mp, an_vector_mp, singlet_mp

mp.mp.dps = 30
# singlets: in V(a) (x) Vbar(b) [x], out Vbar(-b) (x) V(-a) [y]; J swaps v0 <-> v1
for om in (mp.mpf('0.37'), mp.mpf('0.83')):
    q = -mp.exp(-1j*mp.pi*om)
    r = mp.exp(om*1j*mp.pi*(1 - 1/om))         # y_right / y_left for both (pattern of the B_1 constituents)
    x, _ = singlet_mp(an_vector_mp(1, q, 1, False), an_vector_mp(1, q, r, True), 0)
    y, _ = singlet_mp(an_vector_mp(1, q, 1, True), an_vector_mp(1, q, r, False), 1)
    # x: coefficients of v_i (x) vbar_i -> (1 x J): v_i (x) v_{1-i};  y: vbar_i (x) v_i -> (J x 1): v_{1-i} (x) v_i
    print(f'om={om}: in-singlet as V(x)V: v0v1 {mp.nstr(x[0], 10)}, v1v0 {mp.nstr(x[1], 10)};  '
          f'out-singlet: v0v1 {mp.nstr(y[1], 10)}, v1v0 {mp.nstr(y[0], 10)}')


def ghoshal_dirichlet(u, lam, xi):
    pre = (mp.cos(u/2 + mp.pi/(4*lam))*mp.cos(u/2 - mp.pi/4 - mp.pi/(4*lam))*mp.sin(u/2 + mp.pi/4)
           / (mp.cos(u/2 - mp.pi/(4*lam))*mp.cos(u/2 + mp.pi/4 + mp.pi/(4*lam))*mp.sin(u/2 - mp.pi/4)))
    return pre*(mp.cos(xi/lam) - mp.sin(u))/(mp.cos(xi/lam) + mp.sin(u))


for om in ('0.37', '0.83'):
    for eta in ('0.41+0.23j', '-0.3+1.1j', '0.05+1.4j'):
        S = SP(1, mp.mpf(om), 1, mp.mpmathify(eta))
        T = S.T
        xi = mp.pi*S.eps - mp.pi/2
        out = []
        for mu in (mp.mpc('0.11', '0.19'), mp.mpc('-0.23', '0.07'), mp.mpc('0.05', '-0.3')):
            ma, mb = mu - (T - 1)/2, mu + (T - 1)/2
            c, _, _ = proj_mp(S, mu)
            K = S.a(ma)*S.F(2*mu)*S.bt(mb)*c
            sA = lambda m: blk(-T/2, m, T)**2
            u = mp.pi*mu/S.om
            out.append(K*sA(ma)*sA(mb)/ghoshal_dirichlet(u, S.om, xi))
        print(f'om={om} eta={eta}: K_SP sigma_A sigma_A / K_B1(Ghoshal, Dirichlet) = '
              + ', '.join(mp.nstr(v, 12) for v in out))
