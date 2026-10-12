"""Scan over (lambda, eta/eta_N): for which parameters does any valid state (depth <= 3) have a strip pole of class U (not at nu_j, w_j, j>=0)?
Compare with the analytic condition  nu_{-1} < pi/2  <=>  eta < eta_N - pi  (and |0> exists: eta > pi/2).  Also the first state in which U appears."""
import math
from p12b_spectrum_poles import states, classify
from p12b_poles import total_events
rows=[]
print('lambda  eta/eta_N  eta/pi  nstates  first state with U-poles (labels) : U positions   | predicted nu_{-1}<pi/2')
for lam in (1.53, 2.71, 3.77, 5.31, 7.43, 12.67):
    eta_N = math.pi*(lam+1)/2
    for frac in (0.31, 0.53, 0.67, 0.83, 0.911, 0.973):
        eta = frac*eta_N
        sts, nu_, w_ = states(eta, lam, 3)
        firstU=None; nU=0
        for S,p in sts:
            ev = total_events(S, eta, lam)
            r = classify(S, eta, lam, ev, nu_, w_)
            u=[x for x in r if x[2]=='U']
            if u:
                nU+=1
                if firstU is None: firstU=(S,[(x[0],x[1]) for x in u])
        pred = nu_(-1) < math.pi/2 and eta>math.pi/2
        print(f'{lam:5.1f}  {frac:5.2f}  {eta/math.pi:6.3f}  {len(sts):3d}   U in {nU} states; first {firstU}   | {pred}')
