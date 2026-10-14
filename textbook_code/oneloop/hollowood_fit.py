import numpy as np
from massformula import dM_channels
from calib_an_cn import Xa, mass
pi=np.pi
hs = np.arange(2, 30)
v = np.array([dM_channels([(mass(h,b), Xa(h,1,b)) for b in range(1,h)]).real/mass(h,1) for h in hs])
basis = lambda h: [h, 1.0, 1/np.tan(pi/h), 1/np.sin(pi/h), h/np.tan(pi/h), 1/h, np.tan(pi/(2*h)), 1/np.tan(pi/(2*h))]
A = np.array([basis(h) for h in hs])
c, res, *_ = np.linalg.lstsq(A, v*pi, rcond=None)
print("coeffs (x 1/pi):", np.round(c, 8), "resid", np.abs(A@c - v*pi).max())
