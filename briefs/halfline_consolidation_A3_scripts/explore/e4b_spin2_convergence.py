"""Convergence of the spin-2 conservation law d_- T = d_+ Theta on an exact two-soliton solution (4th-order outer stencil).
Log: ../logs/e4b_spin2_convergence.log (run inline on 2026-10-13 with identical code)."""
import io, contextlib
import numpy as np
src = open('e4_spin2.py').read()
cut = src.index('# conservation law d_- T = d_+ Theta checked')
ns = {'__name__': 'e4b'}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src[:cut], 'e4_spin2.py[head]', 'exec'), ns)
phi_at, deriv_jet, TTh, x0s, t0 = ns['phi_at'], ns['deriv_jet'], ns['TTh'], ns['x0s'], ns['t0']
for hh in (1e-3, 5e-4):
    ns['hh'] = hh
    for k in (8e-3, 4e-3, 2e-3):
        worst = 0; scale = 0
        for x0 in x0s[::2]:
            TT = lambda a, b: TTh(phi_at(x0 + a, t0 + b), deriv_jet(x0 + a, t0 + b), 'p', 'q')
            st = [(-2, 1 / 12), (-1, -8 / 12), (1, 8 / 12), (2, -1 / 12)]
            dmT = sum(c * TT(-j * k, j * k)[0] for j, c in st) / k
            dpTh = sum(c * TT(j * k, j * k)[1] for j, c in st) / k
            worst = max(worst, abs(dmT - dpTh)); scale = max(scale, abs(dmT))
        print(f'inner hh={hh:.0e}, outer k={k:.0e} (4th order): |d_- T - d_+ Theta| = {worst:.2e}  (|d_- T| ~ {scale:.1e})')
