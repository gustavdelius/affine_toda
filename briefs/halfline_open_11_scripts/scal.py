"""Scalar factors for the soliton-preserving diagonal K (problem 11).
mu = -i T th/pi.  a: soliton (with KV), b: antisoliton (with lam*KVbar, lam = (-q)^(1-l)).
Equations:  a(mu)a(-mu) = 1/U, U = sinh(eta+i pi mu)sinh(eta-i pi mu);  b(mu)b(-mu) = 1/Ub, Ub = lam^2 sinh(etab+..)sinh(etab-..)
            a(T/2-mu) = Y(mu) b(T/2+mu),  b(T/2-mu) = Y(mu) a(T/2+mu),  Y = -S^T(T-2mu)/s(mu), s = sin pi(2mu-T)/sin pi(2mu+T)."""
import numpy as np
from scipy.special import loggamma as lg

def rich(f, K):
    return 2*f(2*K) - f(K)

def logE(x, c, P, K=2000):
    """E(x+P)/E(x) = const*sin pi(x+c) (regularized Gamma product)."""
    def part(KK):
        k = np.arange(1, KK+1)
        return np.sum(lg(x + c + k*P) - lg(1 - x - c + k*P) - (2*x + 2*c - 1)*np.log(k*P))
    return lg(x + c) + rich(part, K)

def F11(mu, h, lam, K=2000):
    def part(KK):
        j = np.arange(1, KK + 1); a = h*j*lam
        return np.sum(lg(mu + a - (h-1)*lam + 1) + lg(mu + a - lam) + lg(-mu + a - h*lam + 1) + lg(-mu + a)
                      - lg(-mu + a - (h-1)*lam + 1) - lg(-mu + a - lam) - lg(mu + a - h*lam + 1) - lg(mu + a))
    return -np.exp(rich(part, K))

def ST(mu, n, om):
    h = n + 1
    return np.sin(np.pi*mu)/np.sin(np.pi*(om - mu))*F11(mu, h, om)

def Y(mu, n, om):
    T = (n+1)*om/2
    s = np.sin(np.pi*(2*mu - T))/np.sin(np.pi*(2*mu + T))
    return -ST(T - 2*mu, n, om)/s

def logRt(nu, n, om):
    T = (n+1)*om/2; P = 4*T; x = 2*nu
    return logE(x, 2*T+om, P) + logE(x, 2*T-om, P) - logE(x, 0, P) - logE(x, 4*T, P)

def params(n, om, l, eta):
    N = n+1; T = N*om/2
    eps = eta/(1j*np.pi); epsb = eps + om*l - T
    lam = np.exp(-1j*np.pi*om*(1 - l))
    return T, eps, epsb, lam

def logg_a(nu, n, om, l, eta):
    T, eps, epsb, lam = params(n, om, l, eta)
    return (logE(nu, eps, 2*T) + logE(nu, 1-eps, 2*T) - logE(nu, T+epsb, 2*T) - logE(nu, T+1-epsb, 2*T)
            - (nu/T)*np.log(lam))

def logg_b(nu, n, om, l, eta):
    T, eps, epsb, lam = params(n, om, l, eta)
    return (logE(nu, epsb, 2*T) + logE(nu, 1-epsb, 2*T) - logE(nu, T+eps, 2*T) - logE(nu, T+1-eps, 2*T)
            + (nu/T)*np.log(lam))

def U(mu, eta):
    return np.sinh(eta + 1j*np.pi*mu)*np.sinh(eta - 1j*np.pi*mu)

def logL(x, c, P):
    """odd-unitary block: L(x)L(-x)=1, L(x+P)/L(x) = sin pi(x+c) sin pi(c-P-x)."""
    return logE(x, c, P) - logE(-x, c, P)

def logR(mu, n, om):
    T = (n+1)*om/2
    return logL(2*mu, 2*T+om, 4*T) - logL(2*mu, 0, 4*T)

def logga(mu, n, om, l, eta, ratio):
    T, eps, epsb, lam = params(n, om, l, eta)
    return (mu/T)*np.log(ratio) + logL(mu, eps+2*T, 2*T) - logL(mu, epsb+T, 2*T)

def loggb(mu, n, om, l, eta, ratio):
    T, eps, epsb, lam = params(n, om, l, eta)
    return -(mu/T)*np.log(ratio) + logL(mu, epsb+2*T, 2*T) - logL(mu, eps+T, 2*T)

def blk(x, mu, T):
    """DG block (x)_I = sin(pi(mu+x)/2T)/sin(pi(mu-x)/2T)."""
    return np.sin(np.pi*(mu+x)/(2*T))/np.sin(np.pi*(mu-x)/(2*T))

def gamma_shift(n, om, l, eta):
    """k such that gamma -> gamma e^{k mu}, gammabar -> gammabar e^{-k mu} makes Qg = 1 (minimal branch)."""
    T, eps, epsb, lam = params(n, om, l, eta)
    ca = 1j; cb = 1j/lam; ratio = (cb/ca)**2
    mu = 0.11+0.05j
    Z = (cb/ca)*np.sin(np.pi*(eps+T/2-mu))/np.sin(np.pi*(epsb+T/2+mu))
    qg = np.exp(logga(T/2-mu, n, om, l, eta, ratio) - loggb(T/2+mu, n, om, l, eta, ratio))/Z
    return -np.log(qg)/T, ratio

def amps(mu, n, om, l, eta, cdd=None):
    """Return (a, b): soliton and antisoliton scalar factors (multiplying KV and lam*KVbar)."""
    T, eps, epsb, lam = params(n, om, l, eta)
    k, ratio = gamma_shift(n, om, l, eta)
    base = logR(mu, n, om) + np.log(blk(om/2, mu, T))
    a = 1j*np.exp(base + logga(mu, n, om, l, eta, ratio) + k*mu)/np.sin(np.pi*(eps+mu))
    b = (1j/lam)*np.exp(base + loggb(mu, n, om, l, eta, ratio) - k*mu)/np.sin(np.pi*(epsb+mu))
    if cdd is not None:
        c = cdd(mu); a *= c; b *= c
    return a, b
