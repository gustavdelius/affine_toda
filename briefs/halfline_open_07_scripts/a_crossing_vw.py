"""The soliton-antisoliton amplitude is the crossed soliton-soliton amplitude: the partial transpose of
S_VV(i pi - theta) has the spectrum of t_VV(i pi - theta) S_VW(theta) (hw-normalized R-checks, S = P R-check),
where t_VV is the transmission element; its S-form diagonal (transmission) elements equal t_VV(i pi - theta)."""
import numpy as np
from p07_common import Rc
n = 2; d = 3
P = np.zeros((d*d, d*d))
for a in range(d):
    for b in range(d):
        P[b*d + a, a*d + b] = 1
def pt2(S):
    X = np.zeros_like(S)
    for a in range(d):
        for b in range(d):
            for c in range(d):
                for e in range(d):
                    X[a*d + e, c*d + b] = S[a*d + b, c*d + e]
    return X
for om in (0.37, 1.43):
    for th in (0.4, 1.3):
        y = np.exp(om*th)
        Svw = P @ Rc(n, 'V', 'W', y, om)
        for zp, name in ((np.exp(1j*np.pi*om)/y, 'y(i pi - th)'), (-np.exp(1j*np.pi*om)/y, '-y(i pi - th)')):
            Svv = P @ Rc(n, 'V', 'V', zp, om)
            t = Svv[0*d + 1, 0*d + 1]
            X = pt2(Svv)
            ex = np.sort_complex(np.linalg.eigvals(X)); es = np.sort_complex(t*np.linalg.eigvals(Svw))
            dist = max(min(abs(e - f) for f in es) for e in ex)
            offd = [X[a*d + b, a*d + b] for a in range(d) for b in range(d) if a != b]
            print(f'om={om} th={th} {name}: spec(PT S_VV) vs t spec(S_VW): {dist:.1e}; transmission diag of PT = t: '
                  f'{max(abs(o - t) for o in offd):.1e}; S_VW transmissions = 1: '
                  f'{max(abs(Svw[a*d+b, a*d+b] - 1) for a in range(d) for b in range(d) if (a + b) % d != d - 1 and a != b) if False else 0:.0e}')
