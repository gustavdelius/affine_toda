"""Independent referee calculations; imports no book code. Run with python3.

Matrix indices are (out_particle, out_defect, in_particle, in_defect).
All finite-window comparisons exclude the artificial endpoints.
"""
import itertools
import numpy as np
from scipy.special import loggamma, zeta
from scipy.integrate import quad
import sympy as sp

PI = np.pi


def log_r(z, lam, count):
    k = np.arange(count, dtype=float) * lam
    return np.sum(loggamma(k + .25 - z)
                  + loggamma(k + lam + .75 - z)
                  - loggamma(k + lam / 2 + .25 - z)
                  - loggamma(k + lam / 2 + .75 - z))


def phi(theta, lam, count=96):
    z = 1j * lam * theta / (2 * PI)
    # Euler--Maclaurin tail from the asymptotic log-Gamma expansion. The 1/k
    # coefficient cancels in the paired ratio; retain terms k^-2 ... k^-8.
    paired = log_r(z, lam, count) - log_r(-z, lam, count)
    for m, coefficient in tail_coefficients:
        paired += coefficient(z,lam) * zeta(m,count) / lam**m
    return np.exp(1j * PI * (1 + lam) / 4 + paired) / (1 + 1j * np.exp(lam * theta))


zz,ll=sp.symbols('z l')
shifts=[(sp.Rational(1,4),1),(ll+sp.Rational(3,4),1),
        (ll/2+sp.Rational(1,4),-1),(ll/2+sp.Rational(3,4),-1)]
tail_coefficients=[]
for mm in range(2,9):
    coefficient=sp.expand(sum(sign*(sp.bernoulli(mm+1,a-zz)-sp.bernoulli(mm+1,a+zz))
                              for a,sign in shifts)*(-1)**(mm+1)/(mm*(mm+1)))
    tail_coefficients.append((mm,sp.lambdify((zz,ll),coefficient,'numpy')))


def bulk_scalar(theta,lam):
    xi=PI/lam
    def integrand(t):
        if t<1e-8:
            return theta*(PI-xi)/xi
        ratio=2*np.sign(PI-xi)*np.exp(-min(PI,xi)*t)*(-np.expm1(-abs(PI-xi)*t))
        ratio/=(-np.expm1(-xi*t))*(1+np.exp(-PI*t))
        return ratio*np.sin(theta*t)/t
    return -np.exp(-1j*quad(integrand,0,np.inf,epsabs=2e-12)[0])


def matrix(theta, lam, q, radius=5, scalar=True):
    ns = list(range(-radius, radius + 1))
    w = len(ns)
    ans = np.zeros((2, w, 2, w), complex)
    y = np.exp(lam * theta)
    for i, n in enumerate(ns):
        ans[0, i, 0, i] = q ** (-n)
        ans[1, i, 1, i] = q ** n
        if i + 1 < w:
            ans[1, i + 1, 0, i] = -1j * y * q ** (-n)
        if i > 0:
            ans[0, i - 1, 1, i] = -1j * y * q ** n
    return ans.reshape(2 * w, 2 * w) * (phi(theta, lam) if scalar else 1)


def conjugate_particle_transpose(mat, c):
    w = mat.shape[0] // 2
    a = mat.reshape(2, w, 2, w).transpose(2, 1, 0, 3).reshape(2*w, 2*w)
    return np.kron(c, np.eye(w)) @ a @ np.kron(np.linalg.inv(c), np.eye(w))


def six_vertex(theta, lam, q):
    y = np.exp(lam * theta)
    den = q * y - 1 / (q * y)
    b = (y - 1 / y) / den
    c = (q - 1 / q) / den
    return np.array([[1,0,0,0], [0,c,b,0], [0,b,c,0], [0,0,0,1]], complex)


def item_1_2():
    print('\nITEMS 1, 2: non-diagonal doubled boundary and q sign')
    cx = np.array([[0,1], [1,0]], complex)
    ca = np.array([[0,1], [-1,0]], complex)
    radius, w = 5, 11
    interior = [a*w+i for a in range(2) for i in range(2, w-2)]
    ix = np.ix_(interior, interior)
    flip = np.eye(4)[[0,2,1,3]]
    for lam, td, theta in [(1.37,.31,.42), (2.6,-.23,.73), (.73,.17,.29)]:
        q = -np.exp(-1j*PI*lam)
        u = 1j*PI/2-theta
        v = 1j*PI/2+theta
        def doubled_k(arg, sign=1):
            t = matrix(arg-td,lam,sign*q,radius)
            tb = np.linalg.inv(matrix(-arg-td,lam,sign*q,radius))
            zero = np.zeros_like(t)
            return np.block([[zero,tb],[t,zero]]).reshape(4,w,4,w)
        # Include the physical bulk scalar even though k has support only in
        # mixed-copy blocks, so same-copy scattering drops out identically.
        r = bulk_scalar(2*theta,lam)*six_vertex(2*theta,lam,q)
        s = np.zeros((4,4,4,4), complex)
        for a,b in itertools.product(range(4), repeat=2):
            if a//2 != b//2:
                s[b,a,a,b] = 1
            else:
                block = r if a//2 == 0 else flip @ r @ flip
                for c,d in itertools.product(range(2), repeat=2):
                    s[2*(a//2)+c,2*(a//2)+d,a,b] = block[2*c+d,2*(a%2)+b%2]
        def pair(k):
            return np.stack([k[:, :, a^1, :] for a in range(4)], axis=0)
        lhs = pair(doubled_k(u))
        rhs = np.einsum('abij,ijmn->abmn',s,pair(doubled_k(v)))
        residual = np.max(np.abs((lhs-rhs)[:,:,2:-2,2:-2]))
        bad = pair(doubled_k(v)).copy()
        bad[:,2:] *= 1.07
        negative = np.max(np.abs((lhs-np.einsum('abij,ijmn->abmn',s,bad))[:,:,2:-2,2:-2]))
        print(f'lambda={lam}: doubled crossing={residual:.3g}, negative control={negative:.3g}')
        for sign in (1,-1):
            t = matrix(theta-td,lam,sign*q)
            inv = np.linalg.inv(matrix(theta-1j*PI-td,lam,sign*q))
            res_c = np.max(np.abs((t-conjugate_particle_transpose(inv,cx))[ix]))
            res_ca = np.max(np.abs((t-conjugate_particle_transpose(inv,ca))[ix]))
            unit = np.max(np.abs((t.conj().T @ t-np.eye(2*w))[ix]))
            print(f'  q sign {sign:+}: Hermitian unitarity={unit:.3g}, Cx={res_c:.3g}, Ca={res_ca:.3g}')
        # A joint, non-factorized Klein redefinition maps T(-q) to T(q):
        # U_in(+n)=U_in(-n)=(-1)^n, U_out=I. This changes the ZF exchange
        # conventions, not a common single-particle/single-defect basis.
        parity = np.diag([(-1.)**n for n in range(-radius,radius+1)])
        tq = matrix(theta-td,lam,q)
        tm = matrix(theta-td,lam,-q)
        print('  exact one-sided defect parity identity:',
              f'{np.max(np.abs(tm-tq@np.kron(np.eye(2),parity))):.3g}')
    kp=np.array([[0,1],[1,0]],complex)
    km=np.diag([1,-1]).astype(complex)
    j=kp@km
    residual=0
    for n in range(-5,6):
        dn=np.linalg.matrix_power(j,n)
        for ka,kb,beta in [(kp,kp,n),(km,km,n),(kp,km,n+1),(km,kp,n-1)]:
            residual=max(residual,np.max(abs(ka@dn-(-1.)**n*np.linalg.matrix_power(j,beta)@kb)))
    print('  explicit anticommuting Klein generators, all four transmissions:',residual)


def item_4_5():
    print('\nITEMS 4, 5: renormalization and scalar')
    g = sp.symbols('g', positive=True)  # g=b^2/(8*pi)
    print('BS theta_d/theta_cl =', sp.simplify((g/(1+g))/g))
    print('BS scale M/Mclass =', sp.series(sp.sqrt(sp.pi*g/sp.sin(sp.pi*g)),g,0,4))
    print('PP scale =', sp.series(sp.sqrt(1+g),g,0,4))
    for lam in (.73,1.37,2.6):
        ztheta = .41+.21j
        crossing = abs(phi(ztheta,lam)*phi(ztheta+1j*PI,lam)*(1+np.exp(2*lam*ztheta))-1)
        uni = abs(abs(phi(.41,lam))**2*(1+np.exp(2*lam*.41))-1)
        print(f'lambda={lam}: scalar crossing={crossing:.3g}, unitarity={uni:.3g}')
    lam, theta = 1.37,.61
    def kernel(t):
        if t < 1e-7:
            return lam*theta*(lam+1)/(4*PI)
        h = np.exp(-t/4)*(-np.expm1(-(lam+1)*t/2))/(-np.expm1(-t)*(1+np.exp(-lam*t/2)))
        return 2*np.sin(lam*theta*t/(2*PI))*h/t
    phase = quad(kernel,0,np.inf,epsabs=2e-12)[0]
    integral_phi = np.exp(1j*PI*(1+lam)/4+1j*phase)/(1+1j*np.exp(lam*theta))
    print('Independent integral vs paired Gamma product:',f'{abs(integral_phi-phi(theta,lam)):.3g}')
    z = 1j*lam*theta/(2*PI)
    for n in (100,1000,10000):
        growth = np.exp(log_r(z,lam,2*n)-log_r(z,lam,n))
        k = np.arange(1,n+1)*lam
        def compact(count):
            k = np.arange(1,count+1)*lam
            return np.sum(loggamma(.5+2*k-2*z)+loggamma(.5+2*k-lam+2*z)
                          -loggamma(.5+2*k+2*z)-loggamma(.5+2*k-lam-2*z))
        drift = compact(2*n)-compact(n)
        print(f'N={n}: r(2N)/r(N)={growth:.9g}; compact phase increment={drift.imag:.9g}')
    print('Predicted limits:',2**((lam+1)/4),-lam*theta/PI*np.log(2))
    for eps in (1e-3,1e-5,1e-7):
        p = -1j*PI/(2*lam)
        print(f'eps={eps:g}: lower pole residue={eps*phi(p+eps,lam):.9g}; upper cancelled value={phi(-p+eps,lam):.9g}')
    for lam in (.3,.73,1.37,2.6):
        e = np.cosh(.7-1j*PI/(2*lam))
        print(f'resonance energy lambda={lam}: {e:.9g}')


def item_6():
    print('\nITEM 6: actual CCK formulas, ordinary and modified crossing')
    n=5
    ds={a:np.sin((a+1)*PI/n)/np.sin(PI/n) for a in range(n-1)}
    def face(a,b,c,d,t,ordinary=False):
        ratio=ds[a]*ds[c]/(ds[b]*ds[d])
        # Common crossing-symmetric Z is removed, so this tests tensor crossing.
        val=np.sqrt(ratio)*np.sinh(t/n)*(b==d)+np.sinh((1j*PI-t)/n)*(a==c)
        return val*ratio**(1j*t/(2*PI)) if ordinary else val
    allowed=[f for f in itertools.product(range(n-1),repeat=4)
             if all(abs(f[i]-f[(i+1)%4])==1 for i in range(4))]
    t=.73+.12j
    old_res=new_res=ordinary_failure=0.
    for a,b,c,d in allowed:
        ratio=ds[a]*ds[c]/(ds[b]*ds[d])
        old_res=max(old_res,abs(face(a,b,c,d,t,True)-face(b,c,d,a,1j*PI-t,True)))
        new_res=max(new_res,abs(face(a,b,c,d,t)-np.sqrt(ratio)*face(b,c,d,a,1j*PI-t)))
        ordinary_failure=max(ordinary_failure,abs(face(a,b,c,d,t)-face(b,c,d,a,1j*PI-t)))
    print(f'{len(allowed)} allowed faces: ordinary={old_res:.3g}, modified={new_res:.3g}, ungauged ordinary failure={ordinary_failure:.3g}')
    # Pair phase obtained directly from four ONE-kink rephasings.
    a,b,c,d=0,1,2,1
    t1,t2=.8,-.3
    f=lambda x,y,t:(ds[x]/ds[y])**(1j*t/(2*PI))
    state_factor=f(a,b,t1)*f(b,c,t2)/(f(a,d,t2)*f(d,c,t1))
    expected=(ds[a]*ds[c]/(ds[b]*ds[d]))**(1j*(t1-t2)/(2*PI))
    print('four one-kink phases vs amplitude factor:',abs(state_factor-expected))
    # Three-kink braid YBE on every A_4 path, both gauges; common scalars cancel.
    paths=[p for p in itertools.product(range(n-1),repeat=4)
           if all(abs(p[i]-p[i+1])==1 for i in range(3))]
    def op(position,t,ordinary):
        out=np.zeros((len(paths),len(paths)),complex)
        for j,p in enumerate(paths):
            for h in range(n-1):
                pp=list(p);pp[position]=h;pp=tuple(pp)
                if pp in paths:
                    # CCK opposite vacua b,d are path endpoints.
                    out[paths.index(pp),j]=face(p[position],p[position+1],h,p[position-1],t,ordinary)
        return out
    for ordinary in (False,True):
        x,y=.42,.71
        lhs=op(1,x,ordinary)@op(2,x+y,ordinary)@op(1,y,ordinary)
        rhs=op(2,y,ordinary)@op(1,x+y,ordinary)@op(2,x,ordinary)
        unit=op(1,x,ordinary)@op(1,-x,ordinary)
        scalar=np.sinh((x-1j*PI)/n)*np.sinh((-x-1j*PI)/n)
        print(f'ordinary={ordinary}: YBE={np.max(abs(lhs-rhs)):.3g}, algebraic unitarity={np.max(abs(unit-scalar*np.eye(len(paths)))):.3g}')


def ly(theta):
    return (np.sinh(theta)+1j*np.sin(2*PI/3))/(np.sinh(theta)-1j*np.sin(2*PI/3))


def ly_t(theta,b,bar=False):
    return ly(theta+(1 if bar else -1)*1j*PI*(3-b)/6)


def item_7_8():
    print('\nITEMS 7, 8: half-soliton and Lee-Yang')
    tt=lambda t,d:-1j*np.sinh((t-d)/2-1j*PI/4)/np.sinh((t-d)/2+1j*PI/4)
    for lam in (1.37,2.6,5.1):
        t,d=.83,.27
        shift=1j*PI*(1+1/(2*lam))
        prod=tt(t,d+shift)*tt(t,d-shift)
        bulk=(np.sinh(t-d)+1j*np.cos(PI/(2*lam)))/(np.sinh(t-d)-1j*np.cos(PI/(2*lam)))
        print(f'lambda={lam}: two-defect breather identity={abs(prod-bulk):.3g}')
    for b in (-3,3,9,-3+.7j,.5,1.8):
        t=.43+.17j
        ha=abs(ly_t(t,b,True)-ly_t(-t.conjugate(),b).conjugate())
        mirror=abs(ly_t(1j*PI/2-.43,b,True))
        print(f'b={b}: Hermitian analyticity={ha:.3g}; mirror |Tbar|={mirror:.9g}')
    for b in (-.4,.6,1.6):
        t=.43
        u=PI*(b+1)/6
        v=PI*(b-1)/6
        left1=ly_t(t,b)*ly(t+1j*u)
        right1=ly_t(t,b,True)*ly(t-1j*u)
        ident=max(abs(left1-ly_t(t,b+2)),abs(right1-ly_t(t,b+2,True)))
        print(f'b={b}: first excitation shift identity={ident:.3g}')
        if b>1:
            left2=ly_t(t,b)*ly(t+1j*v)
            via1=left1*ly(t-1j*PI*(3-b)/6)
            energy=abs(np.cos(u)+np.cos(PI*(3-b)/6)-np.cos(v))
            print(f'  second excitation two routes={abs(left2-via1):.3g}; energy triangle={energy:.3g}')
            for name,bar,height in [('2+ right simple',True,b-1),('2+ left simple',False,3-b),('2+ right double',True,b+1)]:
                def excited(z):
                    return ly_t(z,b,bar)*ly(z+(-1 if bar else 1)*1j*v)
                p=1j*PI*height/6
                f1=abs(excited(p+1e-5));f2=abs(excited(p+5e-6))
                print(f'  {name}: pole order={np.log2(f2/f1):.7f}')
    # Free-energy expansion: positive mirror factor at b=-3 leads negative E.
    b=-3
    for length in (6,10):
        def integrand(t):
            return np.cosh(t)*ly_t(1j*PI/2-t,b,True)*np.exp(-length*np.cosh(t))/(2*PI)
        val=-quad(lambda t:integrand(t).real,-12,12,epsabs=1e-12)[0]
        print(f'b=-3, m=1, L={length}: Luscher correction={val:.12g}')


if __name__=='__main__':
    item_1_2()
    item_4_5()
    item_6()
    item_7_8()
