"""Test closed forms of the d_n singlet coefficients c_sigma(mu) (K-hat normalization)."""
import sys
import numpy as np
from t10_fitc import cval, n

S = lambda z: np.sin(np.pi*z)


def cform(sig, mu, om):
    T = (n - 1)*om
    if sig == 1:
        m = mu - 0.5
    else:
        m = mu
    return (S(m - T/2 + om/2)*S(m - T/2 - om/2)*S(m - T)*S(m + om/2)*S(m)
            / (S(m - T/2)**2*S(m - om/2)))


rng = np.random.default_rng(5)
for om in (0.53, 0.71, 1.37):
    for name, s in (('s+', 1), ('s-', -1)):
        sig = s*(-1)**n
        ratios = []
        for m in [complex(a, b) for a, b in rng.uniform(-0.6, 0.6, size=(4, 2))]:
            ratios.append(cval(name, om, m, sig)/cform(sig, m, om))
        ratios = np.array(ratios)
        print(f'om={om} {name} sig={sig}: ratio c/cform = {ratios[0]:.10f}, spread {np.ptp(np.abs(ratios - ratios[0])):.1e}; '
              f'1/sin^2(pi om/2) = {1/S(om/2)**2:.10f}; ratio*sin^2 = {ratios[0]*S(om/2)**2:.10f}', flush=True)
