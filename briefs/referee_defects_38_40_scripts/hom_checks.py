"""Exact finite obstruction to non-scalar endomorphisms of an infinite module.
No artificial top/bottom boundary of the bilateral Weyl module is introduced.
"""
import itertools
import sympy as s

q=s.Rational(2); ys=[s.Rational(3,2),s.Rational(5,4)];yd=s.Rational(7,6)
states=list(itertools.product((1,-1),repeat=2)); mu=s.I/(q-1/q)

def gen(l,borel,node):
    M=s.zeros(4); tj=1 if node==1 else -1
    # positive e_j raises charge by 2*tj; negative f_j lowers by 2*tj
    direction=tj*borel
    for col,ab in enumerate(states):
        n=l-sum(ab)//2
        for i in range(3):
            out=list(ab)
            if i<2:
                if ab[i]!=-direction:continue
                out[i]=direction
                if borel==1:amp=ys[i]*q**(tj*sum(ab[:i]))
                else:amp=q**(-tj*(sum(ab[i+1:])+2*n))/ys[i]
            else:
                amp=mu*yd*q**(tj*sum(ab)) if borel==1 else mu/yd
            M[states.index(tuple(out)),col]+=amp
    return M

if __name__=='__main__':
    for borel in (1,-1):
        raisenode=1 if borel==1 else 0; lowernode=1-raisenode
        A=lambda l:gen(l,borel,raisenode)
        B=lambda l:gen(l,borel,lowernode)
        transport=s.eye(4);constraints=[]
        for l in range(3):
            C=s.simplify(transport.inv()*B(l+1)*A(l)*transport)
            constraints.append(s.kronecker_product(C.T,s.eye(4))-s.kronecker_product(s.eye(4),C))
            transport=A(l)*transport
        rank=s.Matrix.vstack(*constraints).rank()
        print('Borel',borel,'loop-commutant rank',rank,'dimension',16-rank,'det raising A0',s.factor(A(0).det()),flush=True)
        assert rank==15
    # Invertibility of raising matrices is structural: particle flips are
    # strictly triangular in the particle weight, and the W-shift term is
    # nonzero diagonal. Every X_l is propagated uniquely from X_0.
