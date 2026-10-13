import numpy as np
from gz3 import X, T, lam
from rhonp import blk
for mu in (0.13+0.21j,-0.3+0.1j,0.05-0.4j):
    print(X(mu)/blk(T/2,mu,T)**2)
