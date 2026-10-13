import sys
sys.path.insert(0, '/home/gustav/Git/affine_toda/textbook_code')
exec(open('base.py').read())

# ---- 8. singlet projection, s = -1
dP = 0
for n, om, mu in ((2, 0.53, 0.11 + 0.19j), (3, 0.61, 0.17 + 0.23j), (4, 0.47, -0.21 + 0.13j), (2, 1.37, 0.05 - 0.12j)):
    h, N = n + 1, n + 1
    T, q = h*om/2, q_phys(om)
    th = 1j*np.pi*mu/T
    u = np.pi*(1 - 2/(h*om))
    ta, tb = th - 0.5j*u, th + 0.5j*u
    Va, Vb = an_vector(n, q, np.exp(om*ta)), an_vector(n, q, np.exp(om*tb), True)
    w_in = singlet(Va, Vb)
    Wa = an_vector(n, q, np.exp(-om*tb))
    w_out = singlet(Wa, an_vector(n, q, np.exp(-om*ta), True))
    R = rcheck(Va, Wa)[0]
    v = kron(np.eye(N), kmat(n, om, ta, -1)) @ (R @ (kron(np.eye(N), kmat(n, om, tb, -1).T) @ w_in))
    c = np.vdot(w_out, v)/np.vdot(w_out, w_out)
    cm = (np.sin(np.pi*mu)*np.sin(np.pi*(mu + om/2))*np.sin(np.pi*(mu - T))
          / (np.sin(np.pi*om/2)**2*np.sin(np.pi*(mu - om/2))))
    print(n, om, 'resid', np.linalg.norm(v - c*w_out)/np.linalg.norm(v), 'c/c_-', c/cm)
    dP = max(dP, np.linalg.norm(v - c*w_out)/np.linalg.norm(v), abs(c/cm - 1))
print('check 8', dP)


# ---- 9. A^(-)
def Am_fun(n, om):
    h = n + 1
    T, L = h*om/2, h*om
    half = mp.mpf(1)/2

    def Am(mu):
        a1 = mp.exp(logprod(2*L, G(2*mu, -1.5*L + 1) + G(2*mu, -0.5*L) + G(2*mu, -1.5*n*om - 0.5*om + 1, -1)
                            + G(2*mu, -0.5*n*om - 1.5*om, -1)))
        a0 = mp.sin(mp.pi*om/2)/mp.cos(mp.pi*(mu - T/2))*mp.exp(
            logprod(L, G(mu, om/2 - 1.5*T + half) + G(mu, half - om/2 - T/2) + G(mu, om/2 - T/2 + half, -1)
                    + G(mu, half - om/2 - 1.5*T, -1)))
        return -1j*blk(-T/2, mu, T)*a0*a1
    return Am


d9 = 0
for n, om, mu in ((1, mp.mpf('0.71'), mp.mpc('0.13', '0.08')), (2, mp.mpf('0.53'), mp.mpc('0.11', '0.19')),
                  (3, mp.mpf('0.61'), mp.mpc('0.17', '0.23')), (4, mp.mpf('0.47'), mp.mpc('-0.21', '0.13')),
                  (2, mp.mpf('1.37'), mp.mpc('0.05', '-0.12'))):
    T, _, F11 = scalar_data(n, om)
    Am = Am_fun(n, om)
    U = -mp.cos(mp.pi*(T/2 + mu))*mp.cos(mp.pi*(T/2 - mu))/mp.sin(mp.pi*om/2)**2
    Psi = mp.cos(mp.pi*(mu - om/2))/mp.cos(mp.pi*(mu + om/2))
    cm = (mp.sin(mp.pi*mu)*mp.sin(mp.pi*(mu + om/2))*mp.sin(mp.pi*(mu - T))
          / (mp.sin(mp.pi*om/2)**2*mp.sin(mp.pi*(mu - om/2))))
    half = mp.mpf(1)/2
    lhs = Am(mu - T/2 + half)*F11(2*mu)*Am(mu + T/2 - half)*cm
    rhs = -blk(om/2, mu, T)*blk(half - T, mu, T)*blk(-om/2 - half, mu, T)
    e = [abs(Am(mu)*Am(-mu)*U - 1), abs(Am(T/2 - mu)/Am(T/2 + mu)*Psi/F11(2*mu) - 1), abs(lhs/rhs - 1)]
    print(n, om, [mp.nstr(x, 3) for x in e])
    d9 = max(d9, *e)
print('check 9', d9)
