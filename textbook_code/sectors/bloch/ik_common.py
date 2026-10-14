"""Izergin-Korepin kink amplitude (conventions of foundations_code/sectors), twisted and graded transfer matrices."""
import sys, os, numpy as np
SECTORS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, SECTORS)
from lib import *
src = open(os.path.join(SECTORS, 'a22three.py')).read().split('q = np.exp(-1j*0.77)')[0]
exec(src)                     # defines setup(q): x -> particle-labelled R = P Rcheck (9x9), and P9
M = np.array([1, 0, -1])      # topological charge m of basis states 0,1,2
G = np.abs(M) % 2             # grading: charged kinks odd
Qg = np.diag([(-1.0)**(G[i]*G[j]) for i in range(3) for j in range(3)])
Pg3 = Qg @ flip(3, 3)
def q_of_xi(xi_pi):
    xi = xi_pi*np.pi; return np.conj(1j*np.exp(1j*np.pi**2/(3*xi)))
def psi_of_xi(xi_pi):
    p = (2/(3*xi_pi)) % 2
    return p - 2 if p > 1 else p      # psi/pi in (-1,1]
def transfer_graded(Rg, ls, D=3):
    Nn = len(ls) + 1; T = np.eye(D**Nn, dtype=complex)
    for j, l in enumerate(ls, start=1):
        v = T
        for k in range(j, 1, -1): v = apply_pair(Pg3, v, k-1, k, D, Nn)
        v = apply_pair(Rg(np.exp(l)), v, 0, 1, D, Nn)
        for k in range(2, j+1): v = apply_pair(Pg3, v, k-1, k, D, Nn)
        T = v
    return T
def charges(Nn):
    import itertools
    return np.array([sum(M[list(c)]) for c in itertools.product(range(3), repeat=Nn)])
def m_first(Nn):
    import itertools
    return np.array([M[c[0]] for c in itertools.product(range(3), repeat=Nn)])
