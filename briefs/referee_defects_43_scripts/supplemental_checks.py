"""Pole residues, Lee-Yang direction, analytic bootstrap and citation coverage."""
from pathlib import Path
import re,json,hashlib
import numpy as np
import sympy as sp
from independent_breathers import scalar,transmit
from independent_particles import T

x,u,v=sp.symbols('x u v',nonzero=True)
first=((x/v+u)/(x/v+1/u))*((x*u+v)/(x*u+1/v))-(x+u*v)/(x+1/(u*v))
second=((-x*v+u)/(-x*v+1/u))*((-x/u+v)/(-x/u+1/v))-(x-u*v)/(x-1/(u*v))
print('exact type-I-minus bootstrap residuals:',sp.factor(first),sp.factor(second))
# Same calculation for type I-plus, whose x=e^D has the opposite shift rule.
firstp=((x*v-u)/(x*v-1/u))*((x/u-v)/(x/u-1/v))-(x-u*v)/(x-1/(u*v))
secondp=((-x/v-u)/(-x/v-1/u))*((-x*u-v)/(-x*u-1/v))-(x+u*v)/(x+1/(u*v))
print('exact type-I-plus bootstrap residuals:',sp.factor(firstp),sp.factor(secondp))

for h in (3,4,5):
    w=np.sqrt(5)
    Dp=1j*(np.pi*(h-1)/h-np.pi/2-np.pi/(h*w))
    vals=[]
    for eps in (1e-5,5e-6):
        p=complex(eps)
        for k in range(h-1):p*=complex(scalar(Dp+eps+1j*np.pi*(2*k-h+2)/h,h,w))
        vals.append(p)
    print(f'conjugate scalar simple pole h={h}: epsilon*gbar={vals}, relative change={abs(vals[0]/vals[1]-1):.3e}',flush=True)

# LY Tbar, not T, supports the two quoted states in the fundamental b interval.
def ly(z):return (np.sinh(z)+1j*np.sin(2*np.pi/3))/(np.sinh(z)-1j*np.sin(2*np.pi/3))
b=.5;shift=1j*np.pi*(3-b)/6;pole=1j*np.pi*(b+1)/6
for eps in (1e-5,5e-6):
    print(f'LY b=.5 epsilon={eps}: eps*T={eps*ly(pole+eps-shift)}, eps*Tbar={eps*ly(pole+eps+shift)}')

# Source/anchor coverage independently parsed from the chapter and bibliography.
chapter=Path('part8-defects/43-particle-transmission.qmd').read_text()
entries={m.group(1):m.group(2) for m in re.finditer(r'(?ms)^@\w+\{([^,]+),(.*?)(?=^@|\Z)',Path('references.bib').read_text())}
refs=sorted(set(re.findall(r'@([A-Za-z0-9_]+)',chapter)))
keys=[x for x in refs if x in entries]
missing=[x for x in refs if '_' in x and x not in entries]
print('citation keys:',len(keys),'missing bibliography keys:',missing)
print('theorem labels:',[(i+1,line) for i,line in enumerate(chapter.splitlines()) if '[Theorem]' in line])
print('exercise count:',len(re.findall(r'\{#exr-',chapter)))
manifest=[]
for key in keys:
    match=re.search(r'eprint\s*=\s*\{([^}]+)',entries[key])
    arxiv=match.group(1) if match else {'gandenberger1995':'hep-th/9501136'}[key]
    source=Path('referee_scratch/defects_43/sources')/(arxiv.replace('/','_')+'.txt')
    assert source.exists(),source
    cited=[i+1 for i,l in enumerate(chapter.splitlines()) if '@'+key in l]
    manifest.append({'key':key,'url':'https://arxiv.org/abs/'+arxiv,'source':str(source),'cited_lines':cited,'sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
Path('referee_scratch/defects_43/citation_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('all cited primary texts present; manifest written')
