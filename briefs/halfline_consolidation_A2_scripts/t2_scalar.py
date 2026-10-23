"""Task 2: the SP scalar factors at 30 digits.  Run: python3 -u t2_scalar.py > t2_scalar.log"""
import mpmath as mp
from sp_mp import SP, F11, mp as _mp

mp.mp.dps = 30

# F11 against the book's an_nondiagonal_scalar.scalar_data (exec of its definitions only)
src = open('/home/gustav/Git/affine_toda/code/an_nondiagonal_scalar.py').read()
defs = src[src.index('mp.mp.dps = 30'):src.index('d3 = 0')]
ns = {'mp': mp}
exec(defs, ns)
dF = 0
for n, om, m in ((2, '0.53', '0.11+0.19j'), (3, '0.61', '-0.27+0.13j'), (1, '0.37', '0.2+0.1j')):
    T, Ap, F11b = ns['scalar_data'](n, mp.mpf(om))
    dF = max(dF, abs(F11(mp.mpmathify(m), n, mp.mpf(om))/F11b(mp.mpmathify(m)) - 1))
print('F11 (mine) vs book scalar_data F11:', mp.nstr(dF, 3))

cases = [(1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (3, 1), (3, 2), (3, 3), (3, 4)]
for om in ('0.37', '0.83'):
    for n, l in cases:
        for eta in ('0.41+0.23j', '-0.3+1.1j'):
            S = SP(n, om, l, mp.mpmathify(eta))
            T = S.T
            worst = {}
            for mu in (mp.mpc('0.11', '0.19'), mp.mpc('-0.23', '0.07')):
                worst['rho unit'] = max(worst.get('rho unit', 0), abs(S.rho(mu)*S.rho(-mu) - 1))
                worst['rho cross'] = max(worst.get('rho cross', 0),
                                         abs(S.rho(T/2 - mu)/S.rho(T/2 + mu)*S.s(mu)/(-S.F(2*mu)) - 1))
                worst['a unit'] = max(worst.get('a unit', 0), abs(S.a(mu)*S.a(-mu)*S.U(mu) - 1))
                worst['b unit'] = max(worst.get('b unit', 0), abs(S.bt(mu)*S.bt(-mu)*S.Ub(mu) - 1))
                worst['cross1'] = max(worst.get('cross1', 0),
                                      abs(S.a(T/2 - mu)*S.Phi1(mu)/(S.F(2*mu)*S.bt(T/2 + mu)) - 1))
                worst['cross2'] = max(worst.get('cross2', 0),
                                      abs(S.bt(T/2 - mu)*S.Phi2(mu)/(S.F(2*mu)*S.a(T/2 + mu)) - 1))
                worst['beta0=alpha0(T-.)'] = max(worst.get('beta0=alpha0(T-.)', 0),
                                                 abs(S.beta0(mu)/S.alpha0(T - mu) - 1))
            print(f'om={om} n={n} l={l} eta={eta}: ' + '  '.join(f'{k} {mp.nstr(v, 2)}' for k, v in worst.items()),
                  flush=True)
