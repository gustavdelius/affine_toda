"""Bowcock-Perkins hep-th/9909174 Table 1 (vacua) and Table 2 (classical K) translated to book conventions C_j=-A_j.
Energies: E_book = -2(n+1) C_BP - 2 sum_j C_j (m/beta^2).  Each BP row stands for its cyclic orbit."""
import numpy as np, itertools
s3=3**.5; s5=5**.5
# (A string, solitons, C_BP, singular, K diag?)  rows in the order of BP tables
ROWS={2:[('+++','none',1,False,'d'),('-++','1 1b',s3-1,False,'d'),('--+','1 1b',s3-1,True,'nd'),('---','none',-1,False,'d')],
3:[('++++','none',1,False,'d'),('-+++','2',1,True,'d'),('-+-+','2',1,True,'d'),('--++','1 1b',2**.5-1,False,'d'),('---+','1 1b',2**.5-1,True,'nd'),('----','none',-1,False,'d')],
4:[('+++++','none',1,False,'d'),('-++++','2 2b',np.sqrt(5+s5)/2**.5-1,False,'d'),('-+-++','2 2b',np.sqrt(5+s5)/2**.5-1,True,'nd'),
   ('--+++','1 1b 2 2b',np.sqrt(5+2*s5)-s5,True,'nd'),('--+-+','1 1b 2 2b',np.sqrt(5+2*s5)-s5,True,'nd'),
   ('---++','1 1b',np.sqrt(5-s5)/2**.5-1,False,'d'),('----+','1 1b',np.sqrt(5-s5)/2**.5-1,True,'nd'),('-----','none',-1,False,'d')],
5:[('++++++','none',1,False,'d'),('-+++++','3',1,True,'d'),('-+-+++','3',1,True,'d'),('-+-+-+','3',1,True,'d'),
   ('-++-++','2 2b',s3-1,False,'d'),('--++++','2 2b',s3-1,True,'d'),('--++-+','2 2b',s3-1,True,'nd'),('--+-++','2 2b',s3-1,True,'nd'),('--+--+','2 2b',s3-1,True,'nd'),
   ('---+++','1 1b 2 2b 3',s3+2-np.sqrt(2*(2**(1/3)+2**(2/3)+2)),True,'nd'),('---+-+','1 1b 2 2b 3',s3+2-np.sqrt(2*(2**(1/3)+2**(2/3)+2)),True,'nd'),
   ('----++','1 1b',0,False,'d'),('-----+','1 1b',0,True,'nd'),('------','none',-1,False,'d')]}
def book_C(A): return tuple(-1 if c=='+' else 1 for c in A)
def class_b(C):
    """DG class b of a solitonic pattern: b=sum (-1)^{a-i} k_i ; returns min(b,h-b)"""
    h=len(C); ks=[k for k in range(1,h) if C[k]==-1]; a=len(ks)
    b=sum((-1)**(a-i)*k for i,k in enumerate(ks,1)); return b
def bp_E(n,Cbp,C): return -2*(n+1)*Cbp-2*sum(C)
if __name__=='__main__':
    for n,rows in ROWS.items():
        print('a_%d'%n)
        for A,sol,Cbp,sing,kd in rows:
            C=book_C(A); orb=len({tuple(np.roll(C,k)) for k in range(n+1)})
            print('  A=%s -> C=%s orbit %d prod C=%+d  sol=%-10s C_BP=%.5f E_book=%.5f  singular@0=%s K:%s %s'%(A,''.join('+' if c>0 else '-' for c in C),orb,np.prod(C),sol,Cbp,bp_E(n,Cbp,C),sing,kd, ('b=%d'%class_b(C) if np.prod(C)==1 and C!=(1,)*(n+1) else '')))
        tot=sum(len({tuple(np.roll(book_C(A),k)) for k in range(n+1)}) for A,*_ in rows); print('  total patterns covered',tot,'of',2**(n+1))
