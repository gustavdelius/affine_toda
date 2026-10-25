"""Reproducible source/provenance audit. Usage: python scope_audit.py SNAPSHOT_ROOT."""
from pathlib import Path
import re,sys,hashlib
root=Path(sys.argv[1]).resolve()
qmd={str(p.relative_to(root)):p.read_text() for p in root.rglob('*.qmd') if '_site' not in p.parts}
ids={}
for path,txt in qmd.items():
    for n,line in enumerate(txt.splitlines(),1):
        for key in re.findall(r'\{#([^ }]+)',line): ids[key]=(path,n)
for key in ['prp-soliton-preserving-half','prp-sp-scalar','prp-dn-uniform-vacuum',
            'prp-dn-uniform-classical','prp-dn-nondiagonal-scalar','prp-twist-defect',
            'prp-sg-spin-three','prp-dc-gauge','prp-defect-theta','prp-ds-absorption',
            'prp-defect-frozen-soliton','prp-fused-bdg','prp-an-particle-transmission','prp-dpt-vacuum']:
    p,n=ids[key]; print('\nCANDIDATE',key,p,n)
    print('\n'.join(qmd[p].splitlines()[n-1:n+7]))
for key in ['prp-cn-smatrix','prp-a2n-smatrix','prp-e6-smatrix','prp-f4-r','prp-f4-breathers',
            'prp-dn-vector-r','prp-sector-table','eq-boundary-anomaly','eq-boundary-topological']:
    print('OLD REPORT TARGET',key,ids.get(key,'ABSENT'))
for path,needle in [('part5-semiclassics/20-one-loop-masses.qmd','The shifts of all mass ratios vanish'),
                    ('preface.qmd','Every result in the book carries a status label')]:
    print('OLD REPORT QUOTE',path,repr(needle),'PRESENT' if needle in qmd[path] else 'ABSENT')
total=0
for path,txt in sorted(qmd.items()):
    if re.match(r'lecture-notes/defects/\d\d-',path):
        count=len(re.findall(r'^::: \{#exr-',txt,re.M)); total+=count
        print('EXERCISES',path,count)
print('TOTAL CLAUDE EXERCISES',total)
assert total==84
print('PASS exercise census: 84 in nine existing lectures')
lua=(root/'status-tags.lua').read_text()
assert 'Proposition = true' in lua and 'Refuted = true' not in lua
print('PASS Proposition is supported; Refuted is not a coloured status')
for fn in ['out_ha.txt','out_qsign/irf_ha_q1.txt','out_qsign/irf_ha_q-1.txt',
           'out_qsign/irf_tm2_q1.txt','out_qsign/irf_tm2_q-1.txt']:
    p=root/'code/irf'/fn
    for n,line in enumerate(p.read_text().splitlines(),1):
        match=re.search(r'W3\((\d+),(\d+)\).*breathers=(yes|no)',line)
        if match:
            a,b=int(match[1]),int(match[2]); expected='yes' if 3*(b-a)>2*a else 'no'
            if match[3]!=expected: print('STALE',str(p.relative_to(root)),n,match[0],'=>',expected)
print('PASS stale output flags independently classified by integer inequality')
for p in ['part7-half-line/33-boundary-quantum-groups.qmd','part9-foundations/48-continuum-window.qmd',
          'part5-semiclassics/20-one-loop-masses.qmd','appendices/E-original-results.qmd','preface.qmd']:
    print('SHA256',p,hashlib.sha256((root/p).read_bytes()).hexdigest())
