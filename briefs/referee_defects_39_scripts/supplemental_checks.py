"""Independent targeted counterexamples and exercise/bibliography audit."""
import hashlib
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]  # repository root
import numpy as np
from independent_checks import geometry, fields, tau_data, delay, charge, residual

for theta,td in [(0.8,-1.4),(-0.8,-3.0)]:
    h,a=5,2
    z=delay(h,a,theta,td)
    values=[]
    for time in [-30.,30.]:
        L=fields(h,tau_data(h,[a],np.array([theta]),np.array([.55j]),0,time))[0]
        R=fields(h,tau_data(h,[a],np.array([theta]),np.array([.55j]),0,time,np.array([z])))[0]
        values.append((R-L)/(2*np.pi))
    difference=charge(h,a,.55)-charge(h,a,.55+np.angle(z))
    jump=values[1]-values[0]
    print(f'3 theta={theta:+g}: jump={np.round(jump.real,6)}, lambda_left-lambda_right={np.round(difference,6)}')
    print(f'  residual with sign(theta)={max(abs(jump-np.sign(theta)*difference)):.3e}')

# At cos(chi)=0 the envelope has no x,t dependence, although the sewing
# equations continue to hold. This is a failure of the physical theorem,
# not a failure of the meromorphic pole identity.
h,a,p,td=4,2,1,.6
th=td+1j*np.pi/2
for x,t in [(0.,-1.),(0.,2.),(.7,.3)]:
    envelope=np.exp((2*np.sin(np.pi*a/h)*(x*np.cosh(th)-t*np.sinh(th))).real)
    print(f'2 degenerate envelope x={x:g},t={t:g}: |E|={envelope:.16f}')

# The pole at h=3,a=1,p=0 has an unavoidable real spacetime tau-zero.
# Choose xi=0, j=0, phase=-pi. At x=t*tanh(td),
# Im log E = -m_a*sin(chi)*t/cosh(td).
h,a,td=3,1,.6
chi=-np.pi/6;ma=2*np.sin(np.pi*a/h)
time=np.pi*np.cosh(td)/(ma*np.sin(chi))
x=time*np.tanh(td)
data=tau_data(h,[a],np.array([td+1j*chi]),np.array([0j]),x,time)
print(f'2 complex-pole singularity h=3 a=1: x={x:.9f}, t={time:.9f}, |tau_0|={abs(data[0][0]):.3e}')

# Negative control: do not accidentally compare the formula only with itself.
h,a,theta,td,p,kind=3,1,.8,.2,1,1
L=fields(h,tau_data(h,[a],np.array([theta]),np.array([.23+.37j]),0,.12))
correct=delay(h,a,theta,td,p,kind)
wrong=delay(h,a,theta,td,0,kind)
for label,z in [('correct charged delay',correct),('wrong class-zero delay',wrong)]:
    R=fields(h,tau_data(h,[a],np.array([theta]),np.array([.23+.37j]),0,.12,np.array([z])),geometry(h)[1][p])
    print(f'1 negative control {label}: sewing residual={residual(h,L,R,td,kind):.9g}')

chapter=ROOT/'part8-defects/39-defect-solitons.qmd'
text=chapter.read_text()
book_files=[ROOT/'index.qmd']+list(ROOT.glob('part*/*.qmd'))
id_pairs=[(id_,str(file)) for file in book_files for id_ in re.findall(r'\{#([^\s}]+)',file.read_text())]
all_ids={id_ for id_,_ in id_pairs}
bib=(ROOT/'references.bib').read_text()
keys=set(re.findall(r'@\w+\{([^,]+),',bib))
refs=set(re.findall(r'@([A-Za-z][\w:.-]*)',text))
missing=sorted(refs-all_ids-keys)
print(f'9 metadata: {len(refs)} distinct chapter targets/citations; unresolved={missing}')
print('9 chapter citation keys:',', '.join(sorted(refs&keys)))
print('9 chapter sha256:',hashlib.sha256(chapter.read_bytes()).hexdigest())
print('9 exercises:',re.findall(r'\{#(exr-[^}]+)\}',text))
