"""Check Lemma 6.13: for a rooted tree with scales s_j <= s_parent, root <= t,
   int prod_j s_j * max(s_j,b)^{-e_j} ds  <=  t^{2h} max(t,b)^{-C} prod_j 2/min(E_j,2h_j)
   whenever all cumulative E_j = 2h_j - C_j > 0.  Exact recursion on a grid (piecewise power laws)."""
import numpy as np, itertools, random
from scipy.integrate import quad

def I(node, t, b, tree, e):
    # integral over subtree of node with s_node <= t ; tree[node] = list of children (internal only)
    f = lambda s: s*max(s,b)**(-e[node])*np.prod([I(c, s, b, tree, e) for c in tree[node]]) if tree[node] else s*max(s,b)**(-e[node])
    pts=[b] if b<t else []
    return quad(f, 0, t, points=pts, limit=200)[0]

def stats(node, tree, e):
    h=1; C=e[node]; out=[]
    for c in tree[node]:
        hc,Cc,oc=stats(c,tree,e); h+=hc; C+=Cc; out+=oc
    out.append((h,C)); return h,C,out

random.seed(1); worst=0; n=0
trees=[{0:[1],1:[]},{0:[1,2],1:[],2:[]},{0:[1],1:[2],2:[]},{0:[1,2],1:[3],2:[],3:[]}]
for tree in trees:
    for trial in range(60):
        e={j:random.uniform(-3,3) for j in tree}
        h,C,out=stats(0,tree,e)
        if any(2*hh-CC<=0.05 for hh,CC in out): continue
        for b,t in [(0.01,1.0),(0.3,1.0),(1.0,0.5),(0.05,0.07)]:
            lhs=I(0,t,b,tree,e)
            rhs=t**(2*h)*max(t,b)**(-C)*np.prod([2/min(2*hh-CC,2*hh) for hh,CC in out])
            worst=max(worst,lhs/rhs); n+=1
print(n,"cases; max lhs/rhs =",worst)
