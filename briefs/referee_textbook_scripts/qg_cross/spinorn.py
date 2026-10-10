import numpy as np, itertools, os, time
from scipy.special import loggamma
from rsolve import intertwiner
N = int(os.environ.get('NSPIN', '4')); omega = float(os.environ.get('OMEGA', '2.37')); J = int(os.environ.get('JTRUNC', '2000'))
q = np.exp(-1j*np.pi*omega); T = N*omega + (N - 1)/2; mu = 2*T; D = 2**N
sp, sm, I2 = np.array([[0,1],[0,0]],complex), np.array([[0,0],[1,0]],complex), np.eye(2)
def site(op, i):
    mats = [I2]*N; mats[i] = op; out = mats[0]
    for m in mats[1:]: out = np.kron(out, m)
    return out
Wn = np.array([[0.5 if b == 0 else -0.5 for b in bits] for bits in itertools.product([0, 1], repeat=N)])
def rep(x):
    e, f, k, qi = {}, {}, {}, {}
    for i in range(1, N):
        e[i] = site(sp, i-1) @ site(sm, i); f[i] = e[i].T.copy(); k[i] = np.diag(q**(2*(Wn[:, i-1] - Wn[:, i]))).astype(complex); qi[i] = q**2
    e[N] = site(sp, N-1); f[N] = e[N].T.copy(); k[N] = np.diag(q**(2*Wn[:, N-1])).astype(complex); qi[N] = q
    e[0] = x*site(sm, 0); f[0] = site(sp, 0)/x; k[0] = np.diag(q**(-2*Wn[:, 0])).astype(complex); qi[0] = q
    return e, f, k, qi
def br(x, l, s): return (x - s*q**l)/(1 - s*x*q**l)
def rho(x):
    out = [1.0]
    for j in range(1, N+1): out.append(out[-1]*br(x, 2*j, (-1)**(j+1)))
    return out
cache = f"Pk_n{N}_w{omega}.npy"
if os.path.exists(cache): Pk = list(np.load(cache))
else:
    x0 = 1.37+0.41j; Rc = intertwiner(rep(x0), rep(1.0), range(N+1), Wn, Wn)[0]; Rc = Rc/Rc[0, 0]; lam = rho(x0); Pk = []
    for kk in range(N+1):
        M = np.eye(D*D, dtype=complex)
        for jj in range(N+1):
            if jj != kk: M = M @ (Rc - lam[jj]*np.eye(D*D))/(lam[kk] - lam[jj])
        Pk.append(M)
    np.save(cache, np.array(Pk))
Rmat = lambda x: sum(r*Pm for r, Pm in zip(rho(x), Pk))
_off = os.environ.get('AOFF')
A = [j*omega + (float(_off.split(',')[j-1]) if _off else (j - 1)/2) for j in range(1, N+1)]
def c_of(t): return np.prod([np.sin(np.pi*(t - a)) for a in A], axis=0)
def lf1(z):
    z = np.asarray(z, dtype=complex); return sum(loggamma(z - a) + loggamma(1 + z + a) for a in A) - len(A)*np.log(np.pi)
js = np.arange(1, J+1); Aj = 2*T*(js-1)
def logf_rel(t):
    d = lambda z, z0: lf1(z) - lf1(z0)
    r = 0.371                                            # generic reference shift (avoids a_j = T coincidences, e.g. n = 2)
    return np.sum(d(t + Aj, Aj + r) - d(t + Aj + T, Aj + T + r) + d(-t + Aj + T, Aj + T + r) - d(-t + Aj + 2*T, Aj + 2*T + r))
_t1 = 0.123
_K = c_of(_t1)*c_of(-_t1)*np.exp(logf_rel(_t1) + logf_rel(-_t1))
_norm = np.sqrt(_K + 0j)
def _Fraw(t): return c_of(t)*np.exp(logf_rel(t))/_norm
_sgn = np.sign(_Fraw(1e-6).real)
def F(theta):
    t = mu*theta/(2j*np.pi); return _sgn*_Fraw(t)
def S(theta): return F(theta)*Rmat(np.exp(mu*theta))
P = np.zeros((D*D, D*D))
for i in range(D):
    for j in range(D): P[j*D+i, i*D+j] = 1
if __name__ == "__main__":
    print(f"n = {N}, omega = {omega}, T = {T}, projector ranks {[int(round(np.trace(Pm).real)) for Pm in Pk]}")
    th = 0.21 + 0.13j
    print(f"unitarity S(th)S(-th) = 1: {np.abs(S(th) @ S(-th) - np.eye(D*D)).max():.1e}")
    rng = np.random.default_rng(0); X = rng.normal(size=(D*D, 60)) + 1j*rng.normal(size=(D*D, 60))
    def probe(t0):
        n_ = [np.linalg.norm(S(1j*np.pi*(t0 + d)/T) @ X) for d in (1e-5, 1e-7)]
        M = 1e-7*S(1j*np.pi*(t0 + 1e-7)/T); sv = np.linalg.svd(M, compute_uv=False)
        return np.log10(n_[1]/n_[0])/2, int(np.sum(sv > 1e-6*sv[0]))
    H = None
    for j in range(1, N):
        od, rk = probe(A[j-1]); print(f"  t = a_{j} = {A[j-1]:7.4f}  {N}+{N} -> soliton {N-j:<2d}  pole order {od:.2f}, residue rank {rk}")
    for t0 in (T - 1, T - 2):
        od, rk = probe(t0); print(f"  t = {t0:7.4f}  breather            pole order {od:.2f}, residue rank {rk}")
    Hh = 2*T/(omega + 0.5)
    M = {N - j: 2*np.cos(np.pi*A[j-1]/(2*T)) for j in range(1, N)}
    print("soliton masses from fusion angles M_a/M_n:", {a: round(v, 6) for a, v in sorted(M.items())})
    print("          2 sin(a pi/H), H = 2T/(w+1/2):   ", {a: round(2*np.sin(a*np.pi/Hh), 6) for a in sorted(M)})
