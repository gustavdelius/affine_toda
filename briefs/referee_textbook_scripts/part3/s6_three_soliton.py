"""eq-n-soliton for N = 3 (a_3^(1), a_4^(1)); reuses toda_classical.py functions."""
import itertools, numpy as np, mpmath as mp
src = open('toda_classical.py').read()
exec('def report' + src.split('def report')[1].split('# ----')[0])
exec('mp.mp.dps = 30\n' + src.split('mp.mp.dps = 30')[1].split('r1 = max(')[0])
r3 = max(field_residual(n, sp_, [mp.mpf('0.5'), mp.mpf('-0.3'), mp.mpf('0.1')], [mp.mpc('0.2', '0.9'), mp.mpc('-0.1', '0.4'), mp.mpc('0.3', '-1.1')], pts=2)
         for n, sp_ in [(3, [1, 2, 3]), (4, [1, 3, 2])])
report(f'three-soliton tau (eq-n-soliton, N=3) solves the field equation  [res {float(r3):.1e}]', r3 < 1e-15)
