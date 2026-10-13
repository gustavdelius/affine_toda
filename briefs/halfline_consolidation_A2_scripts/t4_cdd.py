"""Which CDD factors are compatible with the SP soliton bootstrap?  sigma_a on the soliton, sigma_b on the antisoliton.
Allowed SP CDD pairs: sigma_a 2T-periodic with sigma_a(mu)sigma_a(-mu) = 1, sigma_b(mu) = sigma_a(T - mu).
Bootstrap (the minimal solution satisfies it up to the sign (-1)^(n+1)):  sigma_b(mu) = prod_k sigma_a(mu - T + k om)."""
import mpmath as mp
from sp_mp import blk
mp.mp.dps = 30
mus = [mp.mpc('0.11', '0.19'), mp.mpc('-0.23', '0.07'), mp.mpc('0.31', '-0.12')]
for om in (mp.mpf('0.37'), mp.mpf('0.83')):
    for n in (1, 2, 3, 4, 5):
        N = n + 1; T = N*om/2
        cands = {
            'sigma_A = (-T/2)^2': lambda m: blk(-T/2, m, T)**2,
            'sigma_B': lambda m: -1/(blk(T/2, m, T)**2*blk(T/2 - om/2, m, T)*blk(T/2 + om/2, m, T)),
            'sigma_c=1 (DG)': lambda m: -blk(om/2, m, T)*blk(T - om/2, m, T),
            '(om)_I on a only': lambda m: blk(om, m, T),
            '(om/2)_I on a only': lambda m: blk(om/2, m, T),
            '(T/2)_I on a only': lambda m: blk(T/2, m, T),
        }
        out = []
        for name, sa in cands.items():
            sb = lambda m: sa(T - m)
            # crossing symmetric common factor? (sigma_b == sigma_a)
            common = max(abs(sb(m)/sa(m) - 1) for m in mus)
            boot = max(abs(sb(m)/mp.fprod([sa(m - T + k*om) for k in range(1, n + 1)]) - 1) for m in mus)
            out.append(f'{name}: common {mp.nstr(common, 2)} boot {mp.nstr(boot, 2)}')
        print(f'om={om} n={n}: ' + ' | '.join(out))
