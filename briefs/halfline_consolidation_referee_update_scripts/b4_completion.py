"""Independent B4: differentiate the general Hirota subset sum.

No imports from book or earlier referee scripts. Input definitions are
eq-n-soliton, eq-interaction-coefficient, and the root-space field/boundary
equations. Probe waves have spatial exponent k with k^2=m_c^2-Omega^2.
"""
from itertools import combinations, product
from fractions import Fraction
from math import comb
import sys
import mpmath as mp


class Background:
    def __init__(self, h, solitons):
        self.h, self.solitons = h, solitons
        self.z = mp.exp(2j*mp.pi/h)
        self.m = [2*mp.sin(mp.pi*c/h) for c in range(h)]
        self.terms = []
        for bits in product((0, 1), repeat=len(solitons)):
            ids = [j for j, bit in enumerate(bits) if bit]
            rate = sum(self.m[solitons[j][0]] for j in ids)
            charge = sum(solitons[j][0] for j in ids)
            coeff = mp.fprod(solitons[j][1] for j in ids)
            for a, b in combinations(ids, 2):
                coeff *= self.interaction(solitons[a][0], solitons[b][0], 1)
            self.terms.append((ids, rate, [coeff*self.z**(j*charge) for j in range(h)]))
        self.V = mp.matrix([[self.z**(j*c)/mp.sqrt(h) for c in range(1,h)] for j in range(h)])

    def interaction(self, a, b, cosh):
        return (cosh-mp.cos(mp.pi*(a-b)/self.h))/(cosh-mp.cos(mp.pi*(a+b)/self.h))

    def jets(self, x, order=2, probe=None):
        out = [[mp.mpc(0) for j in range(self.h)] for k in range(order+1)]
        for ids, rate, coeffs in self.terms:
            factor = mp.exp(rate*x)
            if probe is not None:
                c, k = probe
                factor *= mp.exp(k*x)*mp.fprod(self.interaction(self.solitons[j][0], c, k/self.m[c]) for j in ids)
                rate += k
            for j, coeff in enumerate(coeffs):
                value = coeff*factor
                if probe is not None:
                    value *= self.z**(j*probe[0])
                for r in range(order+1):
                    out[r][j] += rate**r*value
        return out

    def weights(self, x):
        t = self.jets(x,0)[0]
        assert max(abs(mp.im(v)) for v in t) < mp.mpf('1e-35')
        assert min(mp.re(v) for v in t) > 0
        return [mp.sqrt(mp.re(t[j-1]*t[(j+1)%self.h]))/mp.re(t[j]) for j in range(self.h)]

    def root_hessian(self, coefficients):
        H = mp.zeros(self.h)
        for j, v in enumerate(coefficients):
            k = (j+1)%self.h
            H[j,j] += v; H[k,k] += v
            H[j,k] -= v; H[k,j] -= v
        return H

    def wave(self, c, k, x=0, order=2):
        t, u = self.jets(x,order), self.jets(x,order,probe=(c,k))
        quotient = []
        for r in range(order+1):
            quotient.append([(u[r][j]-sum(comb(r,l)*t[l][j]*quotient[r-l][j]
                                          for l in range(1,r+1)))/t[0][j]
                             for j in range(self.h)])
        norm = mp.sqrt(self.h)*(1-self.z**(-c))
        return [mp.matrix([(v[j]-v[j-1])/norm for j in range(self.h)]) for v in quotient]

    def pattern(self):
        t = self.jets(0,1); w = self.weights(0)
        candidates = []
        for C in product((1,-1),repeat=self.h):
            residual = [t[1][j]/t[0][j]-C[j]*w[j] for j in range(self.h)]
            candidates.append((max(abs(v-residual[0]) for v in residual), C))
        return min(candidates)

    def boundary_defect(self,C):
        t,tp=self.jets(0,1); w=self.weights(0)
        vals=[tp[j]/t[j]-C[j]*w[j] for j in range(self.h)]
        return max(abs(v-vals[0]) for v in vals)

    def reflect(self, C, energy):
        R = self.root_hessian([c*w/2 for c,w in zip(C,self.weights(0))])
        opened = [c for c in range(1,self.h) if energy > self.m[c]]
        incoming, unknown = [], []
        for c in opened:
            k = 1j*mp.sqrt(energy**2-self.m[c]**2)
            e, ep = self.wave(c,k,order=1)
            incoming.append(self.V.H*(ep+R*e))
            e, ep = self.wave(c,-k,order=1)
            unknown.append(self.V.H*(ep+R*e))
        for c in range(1,self.h):
            if c not in opened:
                k = mp.sqrt(self.m[c]**2-energy**2)
                e, ep = self.wave(c,k,order=1)
                unknown.append(self.V.H*(ep+R*e))
        M = mp.matrix([[v[j] for v in unknown] for j in range(self.h-1)])
        A = mp.matrix([[v[j] for v in incoming] for j in range(self.h-1)])
        columns = [mp.lu_solve(M,-A[:,j]) for j in range(len(opened))]
        X = mp.matrix([[v[j] for v in columns] for j in range(self.h-1)])
        residual = mp.norm(M*X+A)/(1+mp.norm(A))
        return X[:len(opened),:], opened, residual

    def field_residual(self, x):
        t = self.jets(x,2)
        second = [t[2][j]/t[0][j]-(t[1][j]/t[0][j])**2 for j in range(self.h)]
        w = self.weights(x)
        # phi'' = -sum alpha_j (log tau_j)'', grad V=sum alpha_j w_j^2.
        vals = [-second[j]-w[j]**2 for j in range(self.h)]
        return max(abs(v-vals[0]) for v in vals)


def pair(h,a,chi,d=None):
    if d is None:
        d = 1/(1+mp.sin(mp.pi*a/h))
    return Background(h,[(a,d*mp.exp(1j*chi)),(h-a,d*mp.exp(-1j*chi))])


def block(x,theta,h):
    return mp.sinh(theta/2+1j*mp.pi*x/(2*h))/mp.sinh(theta/2-1j*mp.pi*x/(2*h))


def dg(a,b,theta,h):
    v = lambda x: block(x,theta,h)
    result = mp.fprod(v(k-1)*v(h+1-k)/(v(h-k)*v(k)) for k in range(1,a+1))
    return result*mp.fprod(v(a+b+mp.mpf(h)/2-2*k)*v(b-a-mp.mpf(h)/2+2*k)/
                          (v(a+b-mp.mpf(h)/2-2*k)*v(b-a+mp.mpf(h)/2+2*k))
                          for k in range(1,a+1))


def regular():
    worst = dict(field=mp.mpf(0),wave=mp.mpf(0),boundary=mp.mpf(0),reflection=mp.mpf(0),solve=mp.mpf(0))
    count = 0
    for h in range(3,7):
        backgrounds = []
        for a in range(1,(h+1)//2):
            cuts = sorted({(1-Fraction(2*j*a,h))%2 for j in range(h)})
            # Use two interior points in every distinct phase interval.
            for i,left in enumerate(cuts):
                right = cuts[i+1] if i+1<len(cuts) else cuts[0]+2
                left, right = (mp.mpf(v.numerator)/v.denominator for v in (left,right))
                for fraction in (mp.mpf('.37'),mp.mpf('.61')):
                    backgrounds.append((pair(h,a,mp.pi*(left+fraction*(right-left))),a))
        if h%2==0:
            backgrounds += [(Background(h,[(h//2,mp.mpf(d))]),h//2) for d in ('.23','.73')]
        for bg,b in backgrounds:
            boundary,C = bg.pattern(); worst['boundary']=max(worst['boundary'],boundary)
            worst['field']=max(worst['field'],bg.field_residual(mp.mpf('-.37')))
            for energy in (bg.m[1]*mp.cosh(mp.mpf('.83')),mp.mpf('1.9'),mp.mpf('2.7')):
                K,opened,res = bg.reflect(C,energy); worst['solve']=max(worst['solve'],res)
                for i,c in enumerate(opened):
                    theta = mp.acosh(energy/bg.m[c])
                    worst['reflection']=max(worst['reflection'],abs(K[i,i]-dg(c,b,theta,h)),
                                            max([abs(K[j,i]) for j in range(len(opened)) if j!=i]+[mp.mpf(0)]))
                # Check the analytic wave in the differential equation, independently
                # of the wall solve and expected reflection formula.
                c=1; k=1j*mp.sqrt(energy**2-bg.m[c]**2); x=mp.mpf('-.37')
                e,ep,epp=bg.wave(c,k,x)
                H=bg.root_hessian([w*w for w in bg.weights(x)])
                worst['wave']=max(worst['wave'],mp.norm(epp-(H-energy**2*mp.eye(h))*e))
            count+=1
        print(f'a{h-1} completed regular backgrounds',flush=True)
    print('regular total',count,{k:mp.nstr(v,5) for k,v in worst.items()},flush=True)
    assert max(worst.values()) < mp.mpf('1e-40')
    print('PASS regular a2-a5 field, fluctuation, boundary and DG reflection identities',flush=True)


def pairs():
    h=6; a=2; C=(1,1,-1,-1,-1,-1); C2=(1,1,-1,-1,1,-1)
    energy=2*mp.sin(mp.pi/h)*mp.cosh(mp.mpf('.6'))
    records=[]
    for exponent in (3,5,7,9):
        delta=mp.mpf(10)**(-exponent); row=[]; second=[]; defects=[]
        for sign in (1,-1):
            bg=pair(h,a,5*mp.pi/3+sign*delta)
            K,opened,res=bg.reflect(C,energy)
            row.append(abs(K[0,1]))
            defects.append(bg.boundary_defect(C))
            K2,_,_=bg.reflect(C2,energy); second.append(abs(K2[0,1]))
        bg=pair(h,a,5*mp.pi/3,(1-delta)/(1+mp.sin(mp.pi*a/h)))
        K,_,res=bg.reflect(C,energy); row.append(abs(K[0,1]))
        defects.append(bg.boundary_defect(C))
        records.append(row)
        print('pair delta='+mp.nstr(delta,3),'chi+, chi-, d:',[mp.nstr(v,14) for v in row],
              'other pattern:',[mp.nstr(v,12) for v in second],flush=True)
        print('  zeroth-order boundary defects:',[mp.nstr(v,5) for v in defects],flush=True)
    assert abs(records[-1][0]-records[-2][0]) < mp.mpf('1e-10')
    assert abs(records[-1][0]-records[-1][1]) < mp.mpf('1e-10')
    assert abs(records[-1][2]-records[-2][2]) < mp.mpf('1e-6')
    assert abs(records[-1][0]-records[-1][2]) > mp.mpf('.2')
    assert max(defects)<mp.mpf('.001')
    print('PASS stationary-pair limits differ by path; asymmetric target differs by side',flush=True)


def middle():
    for h,C in ((4,(1,-1,1,-1)),(6,(1,-1,1,-1,-1,-1))):
        values=[]
        for exponent in (3,5,7,9):
            delta=mp.mpf(10)**(-exponent); row=[]
            for admixture in (False,True):
                solitons=[(h//2,1-delta)]
                if admixture:
                    solitons += [(1,delta**2*mp.exp(mp.mpf('.31')*1j)),
                                 (h-1,delta**2*mp.exp(-mp.mpf('.31')*1j))]
                bg=Background(h,solitons)
                K,opened,res=bg.reflect(C,bg.m[1]*mp.cosh(mp.mpf('.6')))
                target=dg(1,h//2,mp.mpf('.6'),h)
                row.append(max(abs(K[0,0]-target),abs(K[0,len(opened)-1])))
            values.append(row)
            print(f'middle a{h-1} delta={mp.nstr(delta,3)} errors:',[mp.nstr(v,9) for v in row],flush=True)
        assert max(values[-1]) < mp.mpf('2e-8')
        assert max(values[-1]) < max(values[-2])/50
    print('PASS middle-soliton limits including species-1 admixtures',flush=True)


def exceptional():
    h=5; phase=9*mp.pi/5; C=(1,1,-1,1,-1)
    def build(d1,d2):
        return Background(h,[(1,d1*mp.exp(1j*phase)),(4,d1*mp.exp(-1j*phase)),
                             (2,d2*mp.exp(2j*phase)),(3,d2*mp.exp(-2j*phase))])
    def equations(d1,d2):
        t,tp=build(d1,d2).jets(0,1)
        return mp.re(t[2]),mp.re(tp[0]/t[0]-tp[3]/t[3])
    # At the singular wall, tau2=tau4=0. On the other nodes the
    # exponential boundary term vanishes, leaving equal log derivatives.
    d1,d2=mp.findroot(equations,(mp.mpf('.258'),mp.mpf('.564')))
    bg=build(d1,d2); t,tp=bg.jets(0,1)
    constraint=max(abs(t[2]),abs(t[4]),abs(tp[0]/t[0]-tp[3]/t[3]))
    assert constraint<mp.mpf('1e-60')
    # The regular-node boundary terms decay as sqrt(-x), not as x.
    x=-mp.mpf('1e-30'); tx,tpx=bg.jets(x,1); w=bg.weights(x)
    boundary=[tpx[j]/tx[j]-C[j]*w[j] for j in range(h)]
    finite=max(abs(v-boundary[0]) for v in boundary)
    assert finite<mp.mpf('1e-14')
    print('PASS exceptional background from wall equations: d1='+mp.nstr(d1,30),
          'd2='+mp.nstr(d2,30),'finite boundary residual='+mp.nstr(finite,5),flush=True)
    final=[]
    for exponent in (4,6,8,10):
        delta=mp.mpf(10)**(-exponent); row=[]
        for direction in (0,1):
            b=build(d1*(1-delta if direction==0 else 1),
                    d2*(1-delta if direction==1 else 1))
            assert b.field_residual(mp.mpf('-.37'))<mp.mpf('1e-60')
            K,opened,res=b.reflect(C,b.m[1]*mp.cosh(mp.mpf('.6')))
            assert res<mp.mpf('1e-55')
            row.append(abs(K[0,len(opened)-1]))
        print(f'exceptional delta={mp.nstr(delta,3)} |K14|:',[mp.nstr(v,14) for v in row],flush=True)
        final=row
    assert abs(final[0]-final[1])<mp.mpf('1e-8')
    assert abs(final[0]-mp.mpf('.1471365'))<mp.mpf('1e-6')
    print('PASS exceptional a4 non-diagonal limit along two specified paths',flush=True)


def precision():
    rows=[]
    for digits in (80,120):
        with mp.workdps(digits):
            delta=mp.mpf('1e-9'); h=6; a=2
            energy=2*mp.sin(mp.pi/h)*mp.cosh(mp.mpf('.6'))
            result=[]
            for bg in (pair(h,a,5*mp.pi/3+delta),
                       pair(h,a,5*mp.pi/3,(1-delta)/(1+mp.sin(mp.pi*a/h)))):
                K,opened,res=bg.reflect((1,1,-1,-1,-1,-1),energy)
                result.extend([K[0,0],K[0,1]])
                # Test all open and closed fluctuation species in the wall layer.
                x=-mp.mpf('1e-5')
                H=bg.root_hessian([w*w for w in bg.weights(x)])
                for c in range(1,h):
                    k=mp.sqrt(bg.m[c]**2-energy**2)
                    e,ep,epp=bg.wave(c,k,x)
                    relative=mp.norm(epp-(H-energy**2*mp.eye(h))*e)/(1+mp.norm(epp))
                    assert relative<mp.mpf('1e-35')
            rows.append(result)
    error=max(abs(a-b) for a,b in zip(*rows))
    print('precision comparison:',mp.nstr(error,10),flush=True)
    assert error<mp.mpf('1e-20')
    print('PASS pair wall-layer fluctuation equations and 80/120-digit refinement; '
          'max amplitude change='+mp.nstr(error,6),flush=True)


if __name__=='__main__':
    mp.mp.dps=80
    {'regular':regular,'pairs':pairs,'middle':middle,'exceptional':exceptional,
     'precision':precision}[sys.argv[1]]()
