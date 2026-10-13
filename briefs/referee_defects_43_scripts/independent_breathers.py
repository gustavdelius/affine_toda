"""Independent ch43 audit. No imports of textbook or previous referee code.

Chevalley generators are applied to tensor basis words directly; transmission
is propagated on the full set of reached lattice labels. The Gamma product
uses Bernoulli/Hurwitz-zeta tails, not the book's Richardson extrapolation.
"""
import itertools
from collections import defaultdict
from functools import lru_cache
import numpy as np
import mpmath as mp

mp.mp.dps = 38

def scalar_raw(D, h, w, cutoff=32):
    w, D = mp.mpf(w), mp.mpc(D)
    z = 1j*h*w*(D-1j*mp.pi/2)/(2*mp.pi)
    a0 = (1+(h-1)*w)/2
    shifts = [(z,1),(-z,1),(w+z,-1),(w-z,-1),(0,-2),(w,2)]
    val = -w*(D-1j*mp.pi/2)/2
    for k in range(cutoff):
        val += sum(sign*mp.loggamma(a0+k*h*w+c) for c,sign in shifts)
    # Generalized Stirling expansion of log Gamma(A+c). n=1 cancels.
    for n in range(2,19):
        coeff = (-1)**(n+1)*sum(sign*mp.bernpoly(n+1,c) for c,sign in shifts)/(n*(n+1))
        val += coeff*mp.zeta(n,cutoff+a0/(h*w))/(h*w)**n
    return mp.exp(val)

@lru_cache(None)
def normalization(h,w,cutoff=32):
    D = mp.mpf('0.193')
    q = -mp.exp(-1j*mp.pi*w)
    X = mp.exp(w*(D-1j*mp.pi/2))
    v = scalar_raw(D-1j*mp.pi,h,w,cutoff)*(1-q*X**h)
    for m in range(h-1):
        v *= scalar_raw(D+1j*mp.pi*(2*m-h+2)/h,h,w,cutoff)
    assert abs(mp.im(v)) < mp.mpf('1e-22')*abs(v), v
    return mp.re(v)**(-mp.mpf(1)/h)

def scalar(D,h,w,cutoff=32):
    return scalar_raw(D,h,w,cutoff)*normalization(h,w,cutoff)

def roots(h,j):
    return (j-1)%h,j

def invariant(h,w,rapidities,conjugates):
    """Full common kernel, restricted only by the exact zero-weight condition."""
    q=-np.exp(-1j*np.pi*w)
    signs=np.array([-1 if c else 1 for c in conjugates])
    length=len(rapidities)
    if length==h and not any(conjugates):
        words=list(itertools.permutations(range(h)))
    else:
        words=[]
        for word in itertools.product(range(h),repeat=length):
            weight=sum((signs[p]*np.eye(h,dtype=int)[i] for p,i in enumerate(word)),np.zeros(h,dtype=int))
            if np.all(weight==weight[0]): words.append(word)
    rows={}
    for col,word in enumerate(words):
        for j in range(h):
            a,b=roots(h,j)
            cartan=[signs[p]*((i==a)-(i==b)) for p,i in enumerate(word)]
            for kind in ('e','f'):
                for p,i in enumerate(word):
                    start,end=(b,a) if kind=='e' else (a,b)
                    if conjugates[p]: start,end=end,start
                    if i!=start: continue
                    target=list(word);target[p]=end
                    if kind=='e': c=np.exp(w*rapidities[p])*q**sum(cartan[:p])
                    else: c=np.exp(-w*rapidities[p])*q**(-sum(cartan[p+1:]))
                    key=(kind,j,tuple(target))
                    if key not in rows: rows[key]=np.zeros(len(words),complex)
                    rows[key][col]+=c
    A=np.stack(list(rows.values()))
    _,sv,Vh=np.linalg.svd(A,full_matrices=False)
    return words,Vh[-1].conj(),sv/sv[0]

def transmit(i,N,s,theta,h,w,conjugate=False,td=.27):
    q=-np.exp(-1j*np.pi*w)
    X=np.exp(w*(theta-td-1j*np.pi/2))
    if not conjugate:
        yield i,N,q**(-N[i]-s[i])
        j=(i-1)%h
        out=list(N);out[i]+=1;out[j]-=1
        yield j,tuple(out),X
    else:
        for step in range(h):
            j=(i+step)%h
            out=list(N);out[i]-=1;out[j]+=1
            exponent=sum(N[(i+k)%h]+s[(i+k)%h] for k in range(step+1))
            yield j,tuple(out),X**step*q**exponent

def propagate(words,v,rapidities,conjugates,N,s,h,w):
    state={(word,tuple(N)):c for word,c in zip(words,v)}
    for p in reversed(range(len(rapidities))):
        new=defaultdict(complex)
        for (word,label),c in state.items():
            for j,label2,amp in transmit(word[p],label,s,rapidities[p],h,w,conjugates[p]):
                target=list(word);target[p]=j
                new[(tuple(target),label2)]+=c*amp
        state=new
    expected={(word,tuple(N)):c for word,c in zip(words,v)}
    ratio=sum(np.conj(c)*state.get(key,0) for key,c in expected.items())
    residual=np.sqrt(sum(abs(c-ratio*expected.get(key,0))**2 for key,c in state.items()))/max(1,abs(ratio))
    return ratio,residual

def classical(D,a,h):
    return (np.exp(-D)+np.exp(1j*np.pi*a/h))/(np.exp(-D)+np.exp(-1j*np.pi*a/h))

def run():
    worstscalar=0
    for h in (3,4,5):
        for w in (np.sqrt(2),np.sqrt(5),.6):
            q=-mp.exp(-1j*mp.pi*w)
            D=mp.mpc('.31','.17')
            gs=scalar(D-1j*mp.pi,h,w)
            for k in range(h-1):gs*=scalar(D+1j*mp.pi*(2*k-h+2)/h,h,w)
            cross=abs(gs*(1-q*mp.exp(h*w*(D-1j*mp.pi/2)))-1)
            conv=abs(scalar(D,h,w,32)/scalar(D,h,w,64)-1)
            worstscalar=max(worstscalar,float(cross),float(conv))
            print(f'scalar h={h} omega={w:.12g} crossing={float(cross):.3e} cutoff32/64={float(conv):.3e}',flush=True)
            dev=prop=svmin=0;gap=1
            for a in range(1,h):
                theta=.413+.127j
                u=np.pi/2-np.pi/(h*w)
                ts=[theta-1j*u+1j*np.pi*(2*k-(h-a)+1)/h for k in range(h-a)]
                ts += [theta+1j*u+1j*np.pi*(2*k-a+1)/h for k in range(a)]
                words,v,sv=invariant(h,w,ts,[False]*h)
                svmin=max(svmin,sv[-1]);gap=min(gap,sv[-2])
                fullscalar=complex(mp.fprod(scalar(t-.27,h,w) for t in ts))
                offsets=[np.zeros(h),np.linspace(-.24,.24,h)]
                if h==4:
                    offsets += [np.array([1]*p+[0]*(h-p))-p/h for p in (1,2)]
                for offset in offsets:
                    for N in ([0]*h,[1,-1]+[0]*(h-2)):
                        ratio,residual=propagate(words,v,ts,[False]*h,N,offset,h,w)
                        prop=max(prop,residual)
                        dev=max(dev,abs(ratio*fullscalar-classical(theta-.27,a,h)))
            print(f'all species h={h} omega={w:.12g}: kernel={svmin:.3e} gap={gap:.3e} leakage={prop:.3e} factor_error={dev:.3e}',flush=True)
            assert max(svmin,prop,dev)<1e-9
        # Independent two-module check of Vbar/V and V/Vbar, including reversed shifts.
        w=np.sqrt(2); theta=-.213+.083j;u=np.pi/2-np.pi/(h*w)
        for cs,a in (([True,False],1),([False,True],h-1)):
            ts=[theta-1j*u,theta+1j*u]
            words,v,sv=invariant(h,w,ts,cs)
            reverse=invariant(h,w,ts[::-1],cs)[2][-1]
            ratio,residual=propagate(words,v,ts,cs,[1,-1]+[0]*(h-2),np.linspace(-.17,.17,h),h,w)
            sf=1
            for t,c in zip(ts,cs):
                sf*=complex(mp.fprod(scalar(t-.27+1j*mp.pi*(2*k-h+2)/h,h,w) for k in range(h-1))) if c else complex(scalar(t-.27,h,w))
            error=abs(ratio*sf-classical(theta-.27,a,h))
            print(f'two modules h={h} a={a}: kernel={sv[-1]:.3e} reversed={reverse:.3e} leakage={residual:.3e} factor_error={error:.3e}',flush=True)
            assert max(sv[-1],residual,error)<1e-9 and reverse>.01
    print(f'max scalar validation error={worstscalar:.3e}')

if __name__=='__main__':run()
