"""Exact real-cycle integral Z = N int_{R^{rL}} exp(-e^{-i eps} S/hbar) d^{rL}u by grid quadrature (r=2)."""
import numpy as np
from scipy.signal import fftconvolve

def Z_grid(ch, hbar, h=0.02, pad=None):
    assert ch.r == 2
    L = ch.L
    uL = ch.uL
    if pad is None:
        pad = 9*np.sqrt(hbar*(L+1)/4) + 1.0
    lo = np.minimum(0, uL) - pad; hi = np.maximum(0, uL) + pad
    xs = np.arange(lo[0], hi[0]+h, h); ys = np.arange(lo[1], hi[1]+h, h)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    U = np.stack([X, Y], -1)
    E = np.exp(1j*np.einsum('jr,xyr->xyj', ch.alpha, U))
    W = -np.sum((E-1)*ch.kac, axis=-1)
    f = np.exp(-ch.phase*ch.kappa*W/hbar)
    N1 = (2*np.pi*hbar)**(-1.0)  # per link normalization for r=2
    def kern(D2):
        return N1*np.exp(-ch.phase*D2/(2*hbar))
    v = kern(X**2+Y**2) * f
    # kernel on grid for convolution
    kx = np.arange(-int(pad/h)-1, int(pad/h)+2)*h
    KX, KY = np.meshgrid(kx, kx, indexing='ij')
    K = kern(KX**2+KY**2)
    for x in range(2, L+1):
        v = fftconvolve(v, K, mode='same')*h*h * f
    D2 = (X-uL[0])**2+(Y-uL[1])**2
    return np.sum(v*kern(D2))*h*h
