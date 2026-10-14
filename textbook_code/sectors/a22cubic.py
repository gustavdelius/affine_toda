import numpy as np
def cubic(x, q):
    c3 = (x - q**4)*(x + q**6)
    c2 = q**6*(2 + q**2) + x*(q**2 - 1)*(1 - 3*q**4 + q**8) - q**2*(1 + 2*q**2)*x**2
    c1 = q**6*(2 + q**2)*x**2 + x*(q**2 - 1)*(1 - 3*q**4 + q**8) - q**2*(1 + 2*q**2)
    c0 = (1 - q**4*x)*(1 + q**6*x)
    return np.roots([c3, c2, c1, c0])
def allev(x, q):
    s = np.sqrt(x + 0j)
    return np.concatenate([cubic(x, q), [(1 - q**2*s)/(q**2 - s), (1 + q**2*s)/(q**2 + s), 1]])
lxs = np.concatenate([np.linspace(1e-3, 6, 1200), np.linspace(6, 60, 200)])
print("xi/pi   q=i e^{i pi^2/3xi} ; maxdev x>0 [logx at max] ; maxdev x<0")
for xi_pi in [0.1, 0.2, 0.25, 0.3, 1/3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.62, 0.64, 0.66, 2/3, 0.68, 0.7, 0.75, 0.8, 0.9, 1.0, 1.2, 1.5, 2, 3, 5, 10]:
    xi = xi_pi*np.pi; q = 1j*np.exp(1j*np.pi**2/(3*xi))
    dp = [np.abs(np.abs(allev(np.exp(l), q)) - 1).max() for l in lxs]
    dm = [np.abs(np.abs(allev(-np.exp(l), q)) - 1).max() for l in lxs]
    i = int(np.argmax(dp))
    print(f"{xi_pi:6.3f}  {dp[i]:.2e} [{lxs[i]:.2f}]   {max(dm):.2e}")
