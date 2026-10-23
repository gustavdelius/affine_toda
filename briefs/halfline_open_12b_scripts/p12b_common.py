"""Shared helpers for the part (b) scripts (higher boundary sine-Gordon states).
Book conventions of code/boundary_sine_gordon.py: u = -i theta, BPT/BPTT (3.1)-(3.8).
All functions take lambda (lam) explicitly."""
import mpmath as mp
mp.mp.dps = 20
lg = mp.loggamma
PI = mp.pi

def a_S0(u, lam):
    """BPTT (3.1) a(u) = S_{++}, written as in the book (GZ rho form): a = sin(lam(pi-u)) rho(u)."""
    L = mp.mpf(lam); u = mp.mpmathify(u)
    def logF(l, uu):
        v = L*uu/PI
        return lg(2*l*L - v) + lg(1 + 2*l*L - v) - lg((2*l+1)*L - v) - lg(1 + (2*l-1)*L - v)
    s = mp.nsum(lambda l: logF(l, u) + logF(l, PI - u) - logF(l, 0) - logF(l, PI), [1, mp.inf])
    lr = mp.log(-1/PI) + lg(L) + lg(1 - L*u/PI) + lg(1 - L + L*u/PI) + s
    return mp.sin(L*(PI - u))*mp.exp(lr)

def b_S0(u, lam):
    L = mp.mpf(lam)
    return mp.sin(L*u)/mp.sin(L*(PI - u))*a_S0(u, lam)

def log_R0(u, lam):
    L = mp.mpf(lam); u = mp.mpmathify(u)
    def t(l, uu):
        v = 2*L*uu/PI
        return lg(4*l*L - v) + lg(4*L*(l-1) + 1 - v) - lg((4*l-3)*L - v) - lg((4*l-1)*L + 1 - v)
    return mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])

def log_sigma(x, u, lam):
    L = mp.mpf(lam); x, u = mp.mpmathify(x), mp.mpmathify(u)
    def t(l, uu):
        v = L*uu/PI; a = x/PI
        return (lg(0.5 + a + (2*l-1)*L - v) + lg(0.5 - a + (2*l-1)*L - v)
                - lg(0.5 - a + (2*l-2)*L - v) - lg(0.5 + a + 2*l*L - v))
    return mp.log(mp.cos(x)/mp.cos(x + L*u)) + mp.nsum(lambda l: t(l, u) - t(l, -u), [1, mp.inf])

def blk(x, u):                     # (x) = sin(u/2 + pi x/2)/sin(u/2 - pi x/2)
    return mp.sin(u/2 + x*PI/2)/mp.sin(u/2 - x*PI/2)

def curly(y, u, lam):
    L = mp.mpf(lam)
    return blk((y+1)/(2*L), u)*blk((y-1)/(2*L), u)/(blk((y+1)/(2*L) - 1, u)*blk((y-1)/(2*L) + 1, u))

def a_n(n, eta, u, lam):
    """BPTT (3.5) product form."""
    r = mp.mpf(1)
    for l in range(1, n+1):
        r *= curly(2*(mp.mpf(eta)/PI - l), u, lam)
    return r

def nu(n, eta, lam):
    return mp.mpf(eta)/lam - (2*n+1)*PI/(2*lam)

def etabar(eta, lam):
    return PI*(lam+1) - eta

def Qfun(eta, u, lam):
    """Q(eta, vt, u) without the vt-dependent factors (they cancel in every bootstrap ratio): 
    R0 sigma(eta,u)/cos(eta) * Q0 ; Q0 = -sin(lam u) cos(lam u).  Positive-vt factors sigma(i vt,u)/cosh vt are eta-independent."""
    L = mp.mpf(lam)
    return mp.exp(log_R0(u, lam) + log_sigma(eta, u, lam))/mp.cos(eta)*(-mp.sin(L*u)*mp.cos(L*u))
