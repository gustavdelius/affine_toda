# Independent check of the mass formula: evaluate (8.1) via the heat trace of Prop. 8.5 part 2
# and the Mellin representation, doing the k-integral first at fixed t (the defining order).
import numpy as np
from scipy.integrate import quad
from massformula import dM_channels
def heat_mass(channels):
    def F(t):
        tot = 0.0
        for mb, roots in channels:
            c = 2*mb*sum(-o*w for w, o in roots)
            tot += sum(o*np.exp(-t*mb*mb*(1-w*w)) for w, o in roots if w.real > 0)
            g = lambda k: sum(o/(np.pi*1j)/(k + 1j*mb*w) for w, o in roots)
            I = quad(lambda k: (np.exp(-t*(mb*mb+k*k))*g(k)).real, 0, np.inf, limit=400)[0]
            tot += I + c/(2*np.sqrt(np.pi))*np.sqrt(t)*np.exp(-t*mb*mb)
        return tot
    return -(1/(4*np.sqrt(np.pi)))*quad(lambda t: t**-1.5*F(t), 0, np.inf, limit=400)[0]
for ch in ([(1.0, [(1.0+0j, 1), (-1.0+0j, -1)])],
           [(1.3, [(0.37+0j, 1), (-0.81+0j, 1), (-0.37+0j, -1), (0.81+0j, -1)])],
           [(0.7, [(0.5+0j, 2), (-0.5+0j, -2)]), (1.9, [(0.2+0j, 1), (-0.9+0j, 1), (-0.2+0j, -1), (0.9+0j, -1)])]):
    print("heat-trace:", heat_mass(ch), "  closed form:", dM_channels(ch).real)
