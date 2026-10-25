"""Independent differentiation of Hirota subset polynomials and Robin scattering.

No imports/copies of book helpers. Fields are differentiated with mpmath.diff
from the finite tau polynomial; the fluctuation is the coefficient of one
additional infinitesimal soliton. All computations use 90 decimal digits.
"""
import itertools
import mpmath as m
m.mp.dps=90

class Field:
    def __init__(self,h,solitons):
        self.h=h; self.solitons=solitons
        self.roots=m.eye(h)
        for j in range(h): self.roots[j,(j+1)%h]=-1
    def mass(self,a): return 2*m.sin(m.pi*a/self.h)
    def interaction(self,a,b,c):
        return (c-m.cos(m.pi*(a-b)/self.h))/(c-m.cos(m.pi*(a+b)/self.h))
    def tau(self,x,wave=None):
        out=[]
        for j in range(self.h):
            total=m.mpc(0)
            for size in range(len(self.solitons)+1):
                for subset in itertools.combinations(range(len(self.solitons)),size):
                    term=m.mpc(1)
                    for k in subset:
                        a,z=self.solitons[k]
                        term*=z*m.exp(2j*m.pi*j*a/self.h+self.mass(a)*x)
                        if wave:
                            c,sigma=wave
                            term*=self.interaction(a,c,sigma/self.mass(c))
                    for k,l in itertools.combinations(subset,2):
                        term*=self.interaction(self.solitons[k][0],self.solitons[l][0],1)
                    total+=term
            if wave: total*=m.exp(2j*m.pi*j*wave[0]/self.h+wave[1]*x)
            out.append(total)
        return m.matrix(out)
    def phi(self,x): return -self.roots.T*m.matrix([m.log(v) for v in self.tau(x)])
    def wave(self,x,c,sigma):
        t=self.tau(x); dt=self.tau(x,(c,sigma))
        return -self.roots.T*m.matrix([dt[j]/t[j] for j in range(self.h)])
    def hessian(self,x,pattern=None):
        t=self.tau(x)
        e=[t[(j-1)%self.h]*t[(j+1)%self.h]/t[j]**2 for j in range(self.h)]
        if pattern: e=[pattern[j]*m.sqrt(e[j])/2 for j in range(self.h)]
        return self.roots.T*m.diag(e)*self.roots
    def robin(self,c,sigma,P):
        f=lambda x:self.wave(x,c,sigma)
        return m.diff(f,0)+self.hessian(0,P)*f(0)
    def reflect(self,P):
        energy=self.mass(1)*m.cosh(m.mpf('.6'))
        rates=[-1j*m.sqrt(energy**2-self.mass(c)**2) if energy>self.mass(c)
               else m.sqrt(self.mass(c)**2-energy**2) for c in range(1,self.h)]
        cols=[self.robin(c,rates[c-1],P) for c in range(1,self.h)]
        incoming=self.robin(1,-rates[0],P)
        mat=m.matrix([[col[k] for col in cols] for k in range(self.h-1)])
        z=m.lu_solve(mat,-incoming[:self.h-1,:])
        assert max(abs(sum(cols[k][j]*z[k] for k in range(self.h-1))+incoming[j])
                   for j in range(self.h))<m.mpf('1e-40')
        # Far-left coefficient in component 0 is exp(-2pi i c/h)-1.
        K=[z[c-1]*(m.exp(-2j*m.pi*c/self.h)-1)/(m.exp(-2j*m.pi/self.h)-1)
           for c in range(1,self.h)]
        reverse_in=self.robin(self.h-1,-rates[-1],P)
        reverse=m.lu_solve(mat,-reverse_in[:self.h-1,:])
        reverse_off=reverse[0]*(m.exp(-2j*m.pi/self.h)-1)/(m.exp(-2j*m.pi*(self.h-1)/self.h)-1)
        assert abs(abs(reverse_off)-abs(K[-1]))<m.mpf('1e-40')
        return K[0],abs(K[-1])
    def validate(self):
        x=m.mpf('-.7'); Om=self.mass(1)*m.cosh(m.mpf('.6'))
        sig=1j*m.sqrt(Om**2-self.mass(1)**2)
        f=lambda xx:self.wave(xx,1,sig)
        residual=m.diff(f,x,2)-(self.hessian(x)-Om**2*m.eye(self.h))*f(x)
        assert max(abs(v) for v in residual)<m.mpf('1e-65')

def pair(chi,d=None):
    if d is None: d=1/(1+m.sin(m.pi/3))
    return Field(6,[(2,d*m.exp(1j*chi)),(4,d*m.exp(-1j*chi))])
P=[1,1,-1,-1,-1,-1]; chi=5*m.pi/3
for eps in ['1e-4','1e-6','1e-8']:
    for sign in [-1,1]:
        f=pair(chi+sign*m.mpf(eps)); f.validate(); a,b=f.reflect(P)
        print('a5 regular',sign,eps,m.nstr(a,16),m.nstr(b,16),flush=True)
        assert abs(b-m.mpf('.7191932055'))<m.mpf('1e-6')
print('PASS O2 regular-pair limit 0.719193',flush=True)
for eps in ['1e-6','1e-8','1e-10']:
    f=pair(chi,(1-m.mpf(eps))/(1+m.sin(m.pi/3))); f.validate(); a,b=f.reflect(P)
    print('a5 d path',eps,m.nstr(a,16),m.nstr(b,16),flush=True)
assert abs(b-m.mpf('.9374372024'))<m.mpf('1e-7')
print('PASS O2 d-path limit 0.937437',flush=True)

def exceptional(d1,d2):
    z=d1*m.exp(1j*9*m.pi/5); w=d2*m.exp(1j*18*m.pi/5)
    return Field(5,[(1,z),(4,m.conj(z)),(2,w),(3,m.conj(w))])
def equations(d1,d2):
    f=exceptional(d1,d2); t=f.tau(0); dt=m.diff(f.tau,0)
    return m.re(t[2]),m.re(dt[0]/t[0]-dt[3]/t[3])
d1,d2=m.findroot(equations,(m.mpf('.258'),m.mpf('.564')))
print('endpoint d1,d2',m.nstr(d1,25),m.nstr(d2,25),flush=True)
P=[1,1,-1,1,-1]
for eps in ['1e-5','1e-7','1e-9']:
    for path in [1,2]:
        ee=m.mpf(eps)
        f=exceptional(d1*(1-ee) if path==1 else d1,d2*(1-ee) if path==2 else d2)
        f.validate(); a,b=f.reflect(P)
        print('a4 inward',path,eps,m.nstr(a,16),m.nstr(b,16),flush=True)
assert abs(b-m.mpf('.14713'))<m.mpf('1e-5')
print('PASS O2 exceptional two-inward-path limit 0.14713',flush=True)
