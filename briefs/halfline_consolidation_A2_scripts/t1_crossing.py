"""Task 1: bulk crossing V V <-> Vbar V (normalization of the soliton-antisoliton S-matrix) and the SP crossing
function, with the book's rcheck.  Run: python3 -u t1_crossing.py > t1_crossing.log"""
import numpy as np
from sp_lib import K_V, K_Vbar, etabar, Rc, q_phys

# (a) my K_V, K_Vbar agree with the book's soliton_preserving.py definitions (exec up to the checks)
src = open('/home/gustav/Git/affine_toda/textbook_code/soliton_preserving.py').read().split('\nth1, th2')[0]
ns = {}
exec(src, ns)
d = 0
for n in (1, 2, 3):
    for l in range(n + 2):
        for th in (0.3 + 0.1j, -0.7 + 0.2j):
            d = max(d, np.abs(K_V(n, 0.37, l, ns['ETA'], th) - ns['K_V'](n, l, ns['ETA'], th)).max(),
                    np.abs(K_Vbar(n, 0.37, l, etabar(n, 0.37, l, ns['ETA']), th)
                           - ns['K_Vbar'](n, l, etabar(n, 0.37, l, ns['ETA']), th)).max())
print(f'(a) K_V, K_Vbar agree with soliton_preserving.py: {d:.1e}')

# (b) bulk crossing.  Braid matrix elements S[(out_left,out_right),(in_left,in_right)].
#  form 1: S_VV(th)[(d,c),(a,b)] = S_VbarV(i pi - th)[(c,bbar),(dbar,a)]
#  form 2: S_VV(th)[(d,c),(a,b)] = S_VVbar(i pi - th)[(abar,d),(b,cbar)]
#  with S_VV = F11 Rc_VV, S_VbarV = F_VbarV Rc_VbarV; hypothesis F_VbarV(x) = S^T(T - x), S^T = t F11, so
#  Rc_VV(mu) = t(mu) Rc_VbarV(T - mu) entrywise (identity pairing v_i <-> vbar_i).
print('(b) bulk crossing, max relative deviation of Rc_VV(mu) - t(mu) Rc_(VbarV / VVbar)(T-mu) (crossed entries)')
for n in (1, 2, 3):
    N = n + 1
    for om in (0.37, 1.29):
        T = N*om/2
        for th in (0.31 + 0.12j, -0.45 + 0.27j):
            mu = -1j*T*th/np.pi
            t = np.sin(np.pi*mu)/np.sin(np.pi*(om - mu))
            RVV = Rc(n, om, 'V', 'V', np.exp(om*th))
            RWV = Rc(n, om, 'W', 'V', np.exp(om*(1j*np.pi - th)))
            RVW = Rc(n, om, 'V', 'W', np.exp(om*(1j*np.pi - th)))
            e1 = e2 = 0
            scale = np.abs(RVV).max()
            for a in range(N):
                for b in range(N):
                    for c in range(N):
                        for dd in range(N):
                            lhs = RVV[dd*N + c, a*N + b]
                            e1 = max(e1, abs(lhs - t*RWV[c*N + b, dd*N + a])/scale)
                            e2 = max(e2, abs(lhs - t*RVW[a*N + dd, b*N + c])/scale)
            # transmission entries (i != j) of Rc_VbarV and Rc_VVbar
            tr = max(abs(RWV[j*N + i, i*N + j] - 1) for i in range(N) for j in range(N) if i != j)
            tr2 = max(abs(RVW[j*N + i, i*N + j] - 1) for i in range(N) for j in range(N) if i != j)
            print(f'  n={n} om={om} th={th}: form1 {e1:.1e}  form2 {e2:.1e}  transmission entries -1: {tr:.1e} {tr2:.1e}')

# (c) SP boundary crossing: kB(th) = Phi1 Rc_VVbar(2th) kA(-th), kA(th) = Phi2 Rc_VbarV(2th) kB(-th)
print('(c) SP crossing functions')
for n in (1, 2, 3):
    N = n + 1
    E = np.eye(N)
    for om in (0.37, 1.29):
        T, q = N*om/2, q_phys(om)
        for l in range(N + 1):
            for eta in (0.41 + 0.23j, -0.3 + 1.1j):
                eb = etabar(n, om, l, eta)
                kA = lambda t: sum(K_Vbar(n, om, l, eb, 0.5j*np.pi - t)[i, i]*np.kron(E[i], E[i]) for i in range(N))
                kB = lambda t: sum(K_V(n, om, l, eta, 0.5j*np.pi - t)[i, i]*np.kron(E[i], E[i]) for i in range(N))
                worst = 0
                for th in (0.31 + 0.12j, -0.45 + 0.27j):
                    mu = -1j*T*th/np.pi
                    v1 = Rc(n, om, 'V', 'W', np.exp(2*om*th)) @ kA(-th)
                    v2 = Rc(n, om, 'W', 'V', np.exp(2*om*th)) @ kB(-th)
                    c1 = np.vdot(v1, kB(th))/np.vdot(v1, v1)
                    c2 = np.vdot(v2, kA(th))/np.vdot(v2, v2)
                    r1 = np.linalg.norm(kB(th) - c1*v1)/np.linalg.norm(kB(th))
                    r2 = np.linalg.norm(kA(th) - c2*v2)/np.linalg.norm(kA(th))
                    s = np.sin(np.pi*(2*mu - T))/np.sin(np.pi*(2*mu + T))
                    lam = (-q)**(1 - l)
                    f1, f2 = -lam*s, -s/lam
                    worst = max(worst, r1, r2, abs(c1/f1 - 1), abs(c2/f2 - 1))
                print(f'  n={n} om={om} l={l} eta={eta}: proportionality and closed form '
                      f'Phi1 = -(-q)^(1-l) s, Phi2 = -(-q)^(l-1) s: {worst:.1e}')
