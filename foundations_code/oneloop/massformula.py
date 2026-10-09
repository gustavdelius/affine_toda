"""One-loop soliton mass from transmission factors (counting rule of Section 8).

Convention (foundations paper, Section 8.5): channel-b fluctuation e^{kappa x} e_b at -inf,
X_b e^{kappa x} e_b at +inf, z = kappa/m_b = i k/m_b; bound states at zeros with Re z > 0,
omega = m_b sqrt(1-z^2).

Derived formula (see REPORT.md):
    Delta M = sum_b m_b * sum_{roots w of X_b} ord_w * phi(w)
    phi(w)  = 1/2 sqrt(1-w^2) [Re w>0]  + (1-w^2) L(w)/(2 pi i)  - w/(2 pi)
    L(w)    = int_0^inf dk / ((k + i w) sqrt(k^2+1))
For channels grouped with their conjugates (X_bbar(z) = 1/X_b(-z), Lemma 8.4) and real roots,
this reduces to
    Delta M = 1/2 sum_b m_b sum_{roots} ord * psi(arcsin w),   psi(t) = (t cos t - sin t)/pi.
"""
import numpy as np
from scipy.integrate import quad

def psi(w):
    w = np.asarray(w, dtype=complex)
    t = np.arcsin(w)
    return (t*np.cos(t) - np.sin(t))/np.pi

def L(w):
    f = lambda u: 1.0/(np.sinh(u) + 1j*w)
    re = quad(lambda u: f(u).real, 0, np.inf, limit=400)[0]
    im = quad(lambda u: f(u).imag, 0, np.inf, limit=400)[0]
    return re + 1j*im

def phi_general(w):
    """single-root contribution (per unit m_b), general complex root, not paired."""
    bound = 0.5*np.sqrt(1 - w*w + 0j) if np.real(w) > 0 else 0.0
    return bound + (1 - w*w)*L(w)/(2j*np.pi) - w/(2*np.pi)

def dM_channels(channels, paired=True):
    """channels: list of (m_b, [(root, order), ...]) ; order>0 zero, <0 pole."""
    tot = 0.0
    for mb, roots in channels:
        if paired:
            tot += 0.5*mb*sum(o*psi(w) for w, o in roots)
        else:
            tot += mb*sum(o*phi_general(w) for w, o in roots)
    return complex(tot)

def dM_direct(channels, kmax=np.inf):
    """Direct evaluation: bound states + continuum phase-shift integral + counterterm,
    Delta M = sum_b [ 1/2 sum ord omega + 1/2 int_0^inf E (g - c/(2 pi E^2)) dk + c/(4 pi) ],
    g = (1/(pi i)) d/dk log X_b(k),  c_b = 2 m_b (sum poles - sum zeros)."""
    tot = 0.0
    for mb, roots in channels:
        c = 2*mb*sum(-o*w for w, o in roots)
        bound = sum(o*0.5*mb*np.sqrt(1 - w*w + 0j) for w, o in roots if np.real(w) > 0)
        def integrand(k):
            g = sum(o/(np.pi*1j)/(k + 1j*mb*w) for w, o in roots)
            E = np.sqrt(k*k + mb*mb)
            return 0.5*E*(g - c/(2*np.pi*E*E))
        re = quad(lambda k: integrand(k).real, 0, kmax, limit=800)[0]
        im = quad(lambda k: integrand(k).imag, 0, kmax, limit=800)[0]
        tot += bound + re + 1j*im + c/(4*np.pi)
    return complex(tot)

if __name__ == "__main__":
    # sine-Gordon kink, A = -d^2 + 1 - 2 sech^2: X = (z-1)/(z+1), m = 1  ->  -1/pi
    sg = [(1.0, [(1.0, 1), (-1.0, -1)])]
    print("SG paired  :", dM_channels(sg), " expected", -1/np.pi)
    print("SG general :", dM_channels(sg, paired=False))
    print("SG direct  :", dM_direct(sg))
    # heat-trace route for SG (independent check of Mellin + counterterm):
    from scipy.special import erf
    F = lambda t: erf(np.sqrt(t)) - 2*np.sqrt(t/np.pi)*np.exp(-t)
    val = -(1/(4*np.sqrt(np.pi)))*quad(lambda t: t**-1.5*F(t), 0, np.inf, limit=400)[0]
    print("SG heat-trace Mellin:", val)
    # random self-conjugate channel test: X(z) = (z-a)(z+b)/((z+a)(z-b)), a,b in (0,1)
    a, b = 0.37, 0.81
    ch = [(1.3, [(a, 1), (-b, 1), (-a, -1), (b, -1)])]
    print("test paired :", dM_channels(ch), " general:", dM_channels(ch, paired=False), " direct:", dM_direct(ch))
