"""Shared helpers for problem 7 (a), (b).  Imports an unchanged copy of code/boundary_common.py.
Conventions: Part VI coproduct, q = -exp(-i pi omega), principal gradation, y = exp(omega theta) for a_n vectors.
V = an_vector(conj=False), W = Vbar = an_vector(conj=True), both in the basis v_1..v_N (code indices 0..n)."""
import numpy as np
from boundary_common import an_vector, rcheck, q_phys, k_intertwiner, nodes_of

def Cmat(N):
    return np.eye(N)[::-1].astype(complex)

def nu(j, N):            # diagram automorphism fixing node 0: alpha_j -> alpha_{N-j}
    return (N - j) % N

_rc = {}
def Rc(n, a, b, z, omega):
    """R-check a(z) (x) b(1) -> b (x) a, a,b in {'V','W'}, normalized to 1 on hw (x) hw."""
    key = (n, a, b, complex(z), omega)
    if key not in _rc:
        q = q_phys(omega)
        A = an_vector(n, q, z, a == 'W'); B = an_vector(n, q, 1.0, b == 'W')
        _rc[key] = rcheck(A, B)[0]
    return _rc[key]

def eps_star(omega):     # analytic branch (report of problem 3): 1/(2 cos(pi omega/2))
    return 1/(2*np.cos(np.pi*omega/2))

def k_book(n, theta, s, omega):
    """eq-an-k-matrix on the analytic branch X^p = exp(p omega (theta - i pi/2)); V(y) -> Vbar(1/y)."""
    N = n + 1; q = q_phys(omega); w = np.exp(-1j*np.pi*omega/2)
    Xp = lambda p: np.exp(p*omega*(theta - 1j*np.pi/2))
    K = np.zeros((N, N), complex)
    for i in range(N):
        for j in range(N):
            if i == j:
                K[i, j] = (w*Xp(-N/2) - s/w*Xp(N/2))/(eps_star(omega)*(1/q - q))
            elif j > i:
                K[i, j] = Xp(j - i - N/2)
            else:
                K[i, j] = s*Xp(j - i + N/2)
    return K

def cK(n, theta, s, omega):
    """Krein constant K^dag C K = cK C (claimed closed form)."""
    T = (n + 1)*omega/2
    return (np.cosh(2*T*theta) - s*np.cos(np.pi*T))/(2*np.sin(np.pi*omega/2)**2)

def A_min(n, theta, s, omega):
    """Analytic prefactor of the problem-3 scalar factor: |A_min|^2 = 1/cK on the real axis."""
    T = (n + 1)*omega/2; mu = -1j*T*theta/np.pi
    if s == 1:
        return np.sin(np.pi*omega/2)/np.sin(np.pi*(mu - T/2))
    return -1j*np.sin(np.pi*omega/2)/np.cos(np.pi*(mu - T/2))
