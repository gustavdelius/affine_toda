"""Independent checks of Lecture 3 conventions, momentum, Lax matrix and solitons.

Uses ambient roots e_j-e_{j+1}; no book or referee helper is imported.
Run from the repository root. One PASS/FAIL per group, nonzero on failure.
"""
import numpy as np
from scipy.integrate import quad

checks = []
def check(name, error, tol=2e-9):
    ok = bool(error < tol)
    checks.append(ok)
    print(f'{"PASS" if ok else "FAIL"} {name}: max residual {error:.3g}')

def geometry(h):
    eye = np.eye(h)
    pi = np.roll(eye, -1, axis=0)  # pi e_j = e_{j-1}
    # Rows must be e_j - e_{j+1}.
    roots = np.array([eye[j]-eye[(j+1)%h] for j in range(h)])
    jmat = np.ones((h,h))/h
    q = eye-jmat
    d = 2*np.linalg.inv(eye-pi+jmat)-2*jmat
    return roots, pi, q, d, q-d

def fg(phi, psi, roots, d, b):
    f = np.exp(b/2*roots@(d.T@phi+d@psi))
    g = np.exp(b/2*roots@d@(phi-psi))
    return f,g

def gradients(f,g,roots,d,b,m,sigma):
    bp=m/(2*b)*(sigma*d@roots.T@f+1/sigma*d.T@roots.T@g)
    bs=m/(2*b)*(sigma*d.T@roots.T@f-1/sigma*d.T@roots.T@g)
    return bp,bs

def at(phi, phix, roots, b, m, lam):
    h=len(phi)
    mat=-b/2*np.diag(phix.astype(complex))
    for j in range(h):
        ep=np.exp(b/2*roots[j]@phi)
        mat[j,(j+1)%h] -= m/2*lam*ep
        mat[(j+1)%h,j] += m/2/lam*ep
    return mat

rng=np.random.default_rng(20261025)
errs={k:0. for k in ['necklace','Fourier','momentum','Lax','Hirota','sewing','delay forms','mass','poles','cubic','frozen pair']}
for h in range(2,7):
    roots,pi,q,d,e=geometry(h)
    expected=2*np.eye(h)-2*np.roll(np.eye(h),1,axis=0)
    errs['necklace']=max(errs['necklace'],np.max(abs(roots@d@roots.T-expected)),np.max(abs(d+d.T-2*q)))
    for k in range(1,h):
        v=sum(np.exp(2j*np.pi*k*j/h)*roots[j] for j in range(h))
        dk=1+1j/np.tan(np.pi*k/h)
        errs['Fourier']=max(errs['Fourier'],np.max(abs(d@v-dk*v)))
    for b in [0.8,0.8j]:
        phi=q@rng.normal(size=h)*.2
        psi=q@rng.normal(size=h)*.2
        vt=q@rng.normal(size=h)
        wt=q@rng.normal(size=h)
        m,sigma,lam=1.2,.83,1.27+.3j
        f,g=fg(phi,psi,roots,d,b)
        bp,bs=gradients(f,g,roots,d,b,m,sigma)
        phix=e@vt+d@wt-bp
        psix=d.T@vt-e@wt+bs
        pot=lambda x:m*m/b**2*np.sum(np.exp(b*roots@x)-1)
        pdp=e@bp-d@bs
        pds=-d.T@bp-e@bs
        residual=(-vt@vt-phix@phix+wt@wt+psix@psix)/2+pot(phi)-pot(psi)+pdp@vt+pds@wt
        errs['momentum']=max(errs['momentum'],abs(residual),np.max(abs(e@e+d@d.T-q)))
        left=np.exp(-b/2*(e@phi+d@psi))
        right=np.exp(b/2*(d.T@phi-e@psi))
        st=pi.T
        ld=left[:,None]*(np.eye(h)+sigma/lam*st)*right[None,:]
        ld_dot=(-b/2*(e@vt+d@wt))[:,None]*ld+ld*(b/2*(d.T@vt-e@wt))[None,:]
        res=ld_dot-at(phi,phix,roots,b,m,lam)@ld+ld@at(psi,psix,roots,b,m,lam)
        errs['Lax']=max(errs['Lax'],np.max(abs(res)))
    for a in range(1,h):
        gamma=np.pi*a/h
        theta,theta_d,beta,m=.93,.37,.81,1.13
        sigma=np.exp(-theta_d)
        u=theta-theta_d
        z=(np.exp(-u)-1j*np.exp(1j*gamma))/(np.exp(-u)-1j*np.exp(-1j*gamma))
        z2=np.exp(1j*gamma)*np.sinh((u+1j*np.pi/2+1j*gamma)/2)/np.sinh((u+1j*np.pi/2-1j*gamma)/2)
        mod=(np.exp(-2*u)+2*np.exp(-u)*np.sin(gamma)+1)/(np.exp(-2*u)-2*np.exp(-u)*np.sin(gamma)+1)
        errs['delay forms']=max(errs['delay forms'],abs(z-z2),abs(abs(z)**2-mod))
        for c in range(1,h):
            gc=np.pi*c/h
            zfun=lambda v:(np.exp(-v)-1j*np.exp(1j*gamma))/(np.exp(-v)-1j*np.exp(-1j*gamma))
            paired=zfun(u-1j*(np.pi/2+gc))*zfun(-u-1j*(np.pi/2+gc))
            target=np.exp(2j*gamma)*(np.cosh(u)-np.cos(gc-gamma))/(np.cosh(u)-np.cos(gc+gamma))
            errs['frozen pair']=max(errs['frozen pair'],abs(paired-target))
        ma=2*m*np.sin(gamma)
        phase=np.exp(2j*np.pi*a*np.arange(h)/h)
        for ev in [.13+.07j,.8+.21j,2.1+.31j]:
            tau=1+phase*ev
            ttau=1+phase*z*ev
            phi=1j/beta*roots.T@np.log(tau)
            psi=1j/beta*roots.T@np.log(ttau)
            phix=1j/beta*roots.T@(phase*ev/tau)*ma*np.cosh(theta)
            psix=1j/beta*roots.T@(phase*z*ev/ttau)*ma*np.cosh(theta)
            phit=-np.tanh(theta)*phix
            psit=-np.tanh(theta)*psix
            f,g=fg(phi,psi,roots,d,1j*beta)
            ft=np.roll(tau,-1)*np.roll(ttau,1)/(tau*ttau)
            gt=np.roll(tau,1)*ttau/(tau*np.roll(ttau,1))
            errs['Hirota']=max(errs['Hirota'],np.max(abs(f-ft)),np.max(abs(g-gt)),np.max(abs(np.exp(1j*beta*roots@phi)-np.roll(tau,1)*np.roll(tau,-1)/tau**2)))
            # Check the bulk equation too, independently of the defect.
            wave=-1j/beta*ma**2*roots.T@(phase*ev/tau**2)
            bulk_gradient=m*m/(1j*beta)*roots.T@np.exp(1j*beta*roots@phi)
            errs['Hirota']=max(errs['Hirota'],np.max(abs(wave+bulk_gradient)))
            bp,bs=gradients(f,g,roots,d,1j*beta,m,sigma)
            errs['sewing']=max(errs['sewing'],np.max(abs(phix-e@phit-d@psit+bp)),np.max(abs(psix-d.T@phit+e@psit-bs)))
        # Integrate the static energy density independently, including both kinetic and potential.
        def energy(s):
            ev=np.exp(s+.317j)
            tau=1+phase*ev
            phix=1j/beta*roots.T@(phase*ev/tau)*ma
            exponent=np.roll(tau,1)*np.roll(tau,-1)/tau**2
            return ((phix@phix)/2-m*m/beta**2*np.sum(exponent-1))/ma
        en=quad(lambda s:energy(s).real,-35,35,epsabs=1e-9,limit=300)[0]+1j*quad(lambda s:energy(s).imag,-35,35,epsabs=1e-9,limit=300)[0]
        errs['mass']=max(errs['mass'],abs(en-2*h*ma/beta**2))
        pole=-1j*(np.pi/2-gamma)
        zero=-1j*(np.pi/2+gamma)
        errs['poles']=max(errs['poles'],abs(np.exp(-pole)-1j*np.exp(-1j*gamma)),abs(np.exp(-zero)-1j*np.exp(1j*gamma)))
    for k1 in range(1,h):
        for k2 in range(1,h):
            k3=(-k1-k2)%h
            if k3:
                product=np.prod([1+1j/np.tan(np.pi*k/h) for k in (k1,k2,k3)])
                errs['cubic']=max(errs['cubic'],abs(product.real))
for name,error in errs.items():
    check(name,error)
raise SystemExit(0 if all(checks) else 1)
