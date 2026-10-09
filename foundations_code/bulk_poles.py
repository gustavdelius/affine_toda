"""Bulk vacuum-energy density of simply-laced ATFT in terms of the bare coupling, continued to imaginary coupling.
Fateev's mass-mu relation and the bulk energy f = mbar^2 sin(pi/h)/(8 sin(pi B/h) sin(pi(1-B)/h))
(as quoted in Ahn-Fateev-Kim-Rim-Yang, hep-th/0002213, eqs. (mmu),(fg), with H=h), with b^2 -> -kappa,
kappa = beta^2/4pi in the foundations normalization.  At fixed bare coupling z,
   f  propto  (pi z gamma(1-kappa))^{1/(1-kappa)} * G(kappa),
   G = Gamma(1-1/(h(1-k))) Gamma(1+k/(h(1-k))) / ( B Gamma(1/(h(1-k))) Gamma(1-k/(h(1-k))) ),  B=-k/(1-k).
Poles of G at kappa = 1 - 1/(k h): the neutral collapse thresholds beta^2 = 4 pi (1 - 1/(k h))."""
import numpy as np
from scipy.special import gamma as G
def Gfun(k,h):
    a=1/(h*(1-k)); c=k/(h*(1-k)); B=-k/(1-k)
    return G(1-a)*G(1+c)/(B*G(a)*G(1-c))
for h in [2,3,4,6,12]:
    poles=[1-1/(kk*h) for kk in (1,2,3)]
    kap=np.linspace(0.02,0.97,200000); v=np.abs(Gfun(kap,h))
    # local maxima of |G| that blow up
    idx=[i for i in range(1,len(v)-1) if v[i]>v[i-1] and v[i]>v[i+1] and v[i]>1e3]
    print(f"h={h}: predicted poles {np.round(poles,4)}; numerical blow-ups at {np.round(kap[idx][:3],4)}")
