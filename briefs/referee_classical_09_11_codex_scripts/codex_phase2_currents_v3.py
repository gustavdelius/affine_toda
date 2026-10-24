"""Construct A2 affine-Toda currents by solving the conservation identity, no book code."""
import sympy as s
import json
from pathlib import Path
jets=[[s.Symbol(f'{field}{k}') for k in range(1,7)] for field in ['u','v']]
vars=[jets[a][k-1] for k in range(1,6) for a in range(2)]
weights=[k for k in range(1,6) for a in range(2)]
def D(expr):
    return s.expand(sum(s.diff(expr,jets[a][k])*jets[a][k+1] for a in range(2) for k in range(5)))
def monomials(weight):
    out=[]
    def rec(i,left,term):
        if i==len(vars):
            if left==0: out.append(term)
            return
        for power in range(left//weights[i]+1):rec(i+1,left-power*weights[i],term*vars[i]**power)
    rec(0,weight,s.Integer(1))
    return out
b=[(1,-1),(1,2),(-2,-1)]; f=[(-1,1),(0,-1),(1,0)]
result={}
for spin in [2,4]:
    tm=monomials(spin+1); xm=monomials(spin-1)
    tc=s.symbols(f't0:{len(tm)}'); xc=s.symbols(f'x0:{3*len(xm)}')
    T=sum(c*m for c,m in zip(tc,tm)); X=[sum(xc[j*len(xm)+i]*m for i,m in enumerate(xm)) for j in range(3)]
    equations=[]
    for j in range(3):
        bp=b[j][0]*jets[0][0]+b[j][1]*jets[1][0]
        bell=[s.Integer(1)]
        for k in range(1,spin+1):bell.append(s.expand(D(bell[-1])+bp*bell[-1]))
        dplusT=sum(s.diff(T,jets[a][k])*f[j][a]*bell[k] for a in range(2) for k in range(spin+1))
        residual=s.Poly(s.expand(dplusT+D(X[j])+bp*X[j]),*vars)
        equations.extend(residual.coeffs())
    print('built equations',spin,len(tm),len(xm),flush=True)
    mat,_=s.linear_eq_to_matrix(equations,tc+xc)
    print('matrix',mat.shape,flush=True)
    from sympy.polys.matrices import DomainMatrix
    ns=DomainMatrix.from_Matrix(mat).to_field().nullspace().to_Matrix()
    null=[ns[i,:].T for i in range(ns.rows)]
    print('nullspace done',spin,len(null),flush=True)
    # Find a conservation law not a total D derivative (which has zero integrated T contribution
    # after its full current is adjusted). Test T modulo D(polynomials of weight spin).
    exp_index={s.Poly(m,*vars).monoms()[0]:i for i,m in enumerate(tm)}
    def coeffvec(expr):
        v=s.zeros(len(tm),1)
        for ex,c in s.Poly(s.expand(expr),*vars).terms():v[exp_index[ex]]=c
        return v
    improvements=s.Matrix.hstack(*[coeffvec(D(m)) for m in monomials(spin)])
    rank=improvements.rank()
    print('improvements done',spin,rank,flush=True)
    primitive=next(v for v in null if improvements.row_join(v[:len(tm),:]).rank()>rank)
    primitive=primitive/s.gcd_list(list(primitive))
    lawT=s.expand(T.subs(dict(zip(tc+xc,primitive))))
    lawX=[s.expand(e.subs(dict(zip(tc+xc,primitive)))) for e in X]
    # Strip leading total derivatives to choose a compact representative using rref pivot columns of improvements.
    # Raw representative is reproducible and enough; print it exactly.
    print(f'PASS spin {spin}: conservation nullspace={len(null)}, improvement rank={rank}',flush=True)
    print('T =',lawT,flush=True)
    for j,e in enumerate(lawX):print(f'X{j} =',e,flush=True)
    result[str(spin)]={'T':str(lawT),'X':[str(e) for e in lawX]}
Path('referee_scratch/classical_09_11/codex/codex_phase2_currents_v3.json').write_text(json.dumps(result,indent=2)+'\n')
