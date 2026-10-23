"""A_2^(1) RSOS (IRF) face weights for kinks in the 3 and 3bar of U_q(sl_3), symmetric gauge.

Conventions
  heights: shifted weights h = a + rho = (m1, m2) in Dynkin labels, m1, m2 >= 1, m1 + m2 <= p - 1
           (level k = p - 3 alcove).  Generic (unrestricted) use: pass p=None.
  colours: weights of the 3:  e[0]=(1,0), e[1]=(-1,1), e[2]=(0,-1);  3bar: -e[mu].
  a_{mu nu}(h) = (h, e_mu - e_nu):  a_01 = m1, a_12 = m2, a_02 = m1 + m2.
  [x] = sin(gam x)/sin(gam),  gam = pi*lam,  lam = (p' - p)/p = 4 pi/beta^2 - 1  (Takacs-Watts lambda).
  q-dimension d(h) = [m1][m2][m1+m2]/[2].
  face W_AB(a,b,c,d|u):  a left, b bottom (in), c right, d top (out);
     in kinks K_ab (type A, rapidity th1), K_bc (type B, th2) -> out K_ad (type B, th2), K_dc (type A, th1).
  additive spectral parameter u = 3 i theta/(2 pi): theta = i pi  <->  u = -3/2 (crossing).
"""
import numpy as np, itertools

E = [(1, 0), (-1, 1), (0, -1)]
EB = [(-1, 0), (1, -1), (0, 1)]


def add(h, w):
    return (h[0] + w[0], h[1] + w[1])


def sub(h, w):
    return (h[0] - w[0], h[1] - w[1])


class A2RSOS:
    def __init__(self, p, pp=None, lam=None):
        """p: restriction (heights m1+m2<=p-1); pp: p'; or give lam directly (then gam = pi lam)."""
        self.p = p
        if lam is None:
            lam = (pp - p) / p
        self.lam = lam
        self.gam = np.pi * lam
        if p is not None:
            self.H = [(m1, m2) for m1 in range(1, p) for m2 in range(1, p) if m1 + m2 <= p - 1]
        else:
            self.H = None
        self.Hset = set(self.H) if self.H else None

    # q-number with complex argument
    def br(self, x):
        return np.sin(self.gam * x) / np.sin(self.gam)

    def allowed(self, h):
        if h[0] < 1 or h[1] < 1:
            return False
        if self.p is None:
            return True
        return h[0] + h[1] <= self.p - 1

    def qdim(self, h):
        m1, m2 = h
        return (self.br(m1) * self.br(m2) * self.br(m1 + m2) / self.br(2)).real

    @staticmethod
    def amn(h, mu, nu):
        # (h, e_mu - e_nu) with Dynkin inner product: (h, alpha_i) = m_i
        lab = {0: 0, 1: 0, 2: 0}
        # e0-e1 = alpha1, e1-e2 = alpha2, e0-e2 = alpha1+alpha2
        m1, m2 = h
        tab = {(0, 1): m1, (1, 2): m2, (0, 2): m1 + m2}
        if (mu, nu) in tab:
            return tab[(mu, nu)]
        return -tab[(nu, mu)]

    @staticmethod
    def colour(diff, bar=False):
        L = EB if bar else E
        return L.index(diff) if diff in L else None

    # ---- the basic 3 x 3 face weight (JMO, symmetric gauge), unnormalized ----
    def W33(self, a, b, c, d, u, sqrt_branch=None):
        mu = self.colour(sub(b, a)); nu = self.colour(sub(c, b))
        mp = self.colour(sub(d, a)); nup = self.colour(sub(c, d))
        if None in (mu, nu, mp, nup):
            return 0.0
        br = self.br
        if mu == nu:
            return br(1 + u) / br(1) if d == b else 0.0
        A = self.amn(a, mu, nu)
        if d == b:
            return br(A - u) / br(A)
        # d = a + e_nu
        s = (br(A + 1) * br(A - 1)).real / br(A).real ** 2
        r = np.sqrt(complex(s))
        if sqrt_branch is not None:
            r = sqrt_branch(r, a, b, c, d)
        return r * br(u) / br(1)


def steps(t):
    return EB if t == 1 else E


def paths(model, types, start=None, periodic=False):
    """all height paths h0 -> h1 -> ... -> hN with h_{i}-h_{i-1} in wt(types[i-1])."""
    out = []
    starts = model.H if start is None else [start]
    for h0 in starts:
        stack = [(h0,)]
        while stack:
            pth = stack.pop()
            i = len(pth) - 1
            if i == len(types):
                if (not periodic) or pth[-1] == pth[0]:
                    out.append(pth)
                continue
            for w in steps(types[i]):
                h = add(pth[-1], w)
                if model.allowed(h):
                    stack.append(pth + (h,))
    return sorted(out)


def local_op(model, Wfun, types, i, u, src, dst_index):
    """X_i(u): acts on height h_i (between kinks i and i+1, 1-based kinks i,i+1 -> positions i-1,i in types).
    Wfun(A,B,a,b,c,d,u).  returns matrix dst x src."""
    A, B = types[i - 1], types[i]
    M = np.zeros((len(dst_index), len(src)), complex)
    for col, pth in enumerate(src):
        a, b, c = pth[i - 1], pth[i], pth[i + 1]
        for w in steps(B):
            d = add(a, w)
            if not model.allowed(d):
                continue
            if sub(c, d) not in steps(A):
                continue
            new = pth[:i] + (d,) + pth[i + 1:]
            val = Wfun(A, B, a, b, c, d, u)
            if val != 0:
                M[dst_index[new], col] += val
    return M


def swapped(types, i):
    t = list(types); t[i - 1], t[i] = t[i], t[i - 1]; return tuple(t)


def ybe_error(model, Wfun, types, u12, u23, periodic=False):
    u13 = u12 + u23
    P0 = paths(model, types)
    def idx(ts):
        ps = paths(model, ts); return ps, {p: k for k, p in enumerate(ps)}
    t1 = swapped(types, 1); P1, I1 = idx(t1)
    t12 = swapped(t1, 2); P12, I12 = idx(t12)
    t121 = swapped(t12, 1); P121, I121 = idx(t121)
    L = local_op(model, Wfun, t12, 1, u23, P12, I121) @ local_op(model, Wfun, t1, 2, u13, P1, I12) @ local_op(model, Wfun, types, 1, u12, P0, I1)
    t2 = swapped(types, 2); Q2, J2 = idx(t2)
    t21 = swapped(t2, 1); Q21, J21 = idx(t21)
    t212 = swapped(t21, 2); Q212, J212 = idx(t212)
    assert t212 == t121
    Rm = local_op(model, Wfun, t21, 2, u12, Q21, J212) @ local_op(model, Wfun, t2, 1, u13, Q2, J21) @ local_op(model, Wfun, types, 2, u23, P0, J2)
    # align bases (same type sequence, same path list)
    assert Q212 == P121
    return np.abs(L - Rm).max(), np.abs(L).max()


# ---------------- symmetric-gauge weights for all four sectors ----------------
def make_weights(m):
    """returns Wfun(A,B,a,b,c,d,u); A,B in {0 (=3), 1 (=3bar)}; per-height square roots s(h)=sqrt(d_h)."""
    s = {h: np.sqrt(complex(m.qdim(h))) for h in m.H}
    br = lambda r, a, b, c, d: s[b] * s[d] / (s[a] * s[c])
    W33 = lambda a, b, c, d, u: m.W33(a, b, c, d, u, sqrt_branch=br)
    kap = lambda a, b, c, d: s[b] * s[d] / (s[a] * s[c])
    def W(A, B, a, b, c, d, u):
        if (A, B) == (0, 0): return W33(a, b, c, d, u)
        if (A, B) == (1, 0): return kap(a, b, c, d) * W33(b, c, d, a, -1.5 - u)
        if (A, B) == (0, 1): return kap(a, b, c, d) * W33(d, a, b, c, -1.5 - u)
        return W33(c, d, a, b, u)
    return W


def uof(theta):
    return 1.5j * theta / np.pi


def block(m, W, A, B, a, c, u):
    """2-body IRF S-matrix block (fixed left a, right c): rows = out middle d, cols = in middle b."""
    ins = [add(a, w) for w in steps(A) if m.allowed(add(a, w)) and sub(c, add(a, w)) in steps(B)]
    outs = [add(a, w) for w in steps(B) if m.allowed(add(a, w)) and sub(c, add(a, w)) in steps(A)]
    M = np.array([[W(A, B, a, b, c, d, u) for b in ins] for d in outs], complex)
    return ins, outs, M


def periodic_paths(m, types):
    return [p for p in paths(m, types) if p[-1] == p[0]]


def transfer(m, W, types, thetas, scal=None):
    """Bethe-Yang transfer matrix of kink 1 on periodic paths of the given types.
    scal(A,B,theta): scalar prefactor of S_AB.  Returns (paths, T)."""
    N = len(types)
    cur_types = tuple(types)
    P0 = periodic_paths(m, types)
    cur = P0
    T = np.eye(len(P0), dtype=complex)
    for j in range(1, N):
        # kink 1 is at position j (1-based), scatters with kink at j+1
        nt = swapped(cur_types, j)
        nxt = sorted({q for q in paths(m, nt) if q[-1] == q[0]})
        Inx = {q: k for k, q in enumerate(nxt)}
        u = uof(thetas[0] - thetas[j])
        X = local_op(m, W, cur_types, j, u, cur, Inx)
        if scal is not None:
            X = X * scal(cur_types[j - 1], cur_types[j], thetas[0] - thetas[j])
        T = X @ T
        cur, cur_types = nxt, nt
    # cyclic relabel: (g0, g1',...,g_{N-1}', g0) -> (g_{N-1}', g0, g1', ..., g_{N-1}')
    I0 = {q: k for k, q in enumerate(P0)}
    Cm = np.zeros((len(P0), len(cur)))
    for k, q in enumerate(cur):
        new = (q[N - 1],) + q[:N - 1] + (q[N - 1],)
        Cm[I0[new], k] = 1
    return P0, Cm @ T


# ---------------- faster periodic transfer matrices ----------------
class Fast:
    def __init__(self, m, W=None):
        self.m = m; self.W = W if W is not None else make_weights(m)
        self.nb = {t: {h: [add(h, w) for w in steps(t) if m.allowed(add(h, w))] for h in m.H} for t in (0, 1)}
        self.cache = {}

    def per(self, types):
        types = tuple(types)
        if types not in self.cache:
            out = []
            for h0 in self.m.H:
                stack = [(h0,)]
                while stack:
                    pth = stack.pop(); i = len(pth) - 1
                    if i == len(types):
                        if pth[-1] == pth[0]: out.append(pth)
                        continue
                    for h in self.nb[types[i]][pth[-1]]:
                        stack.append(pth + (h,))
            out = sorted(out)
            self.cache[types] = (out, {p: k for k, p in enumerate(out)})
        return self.cache[types]

    def scal(self, A, B, u):
        br = self.m.br
        if A == B:
            return br(1) / br(1 + u)          # F_33 = Phi/rho with |Phi| = 1 (phase dropped)
        return 1.0 / np.sqrt((br(1.5 + u) * br(1.5 - u) / br(1) ** 2).real)  # |F_33bar| = sigma^{-1/2}

    def X(self, types, i, u, normalize=True):
        A, B = types[i - 1], types[i]
        src, _ = self.per(types)
        nt = swapped(types, i); dst, Id = self.per(nt)
        M = np.zeros((len(dst), len(src)), complex)
        for col, pth in enumerate(src):
            a, b, c = pth[i - 1], pth[i], pth[i + 1]
            for d in self.nb[B][a]:
                if sub(c, d) not in steps(A): continue
                v = self.W(A, B, a, b, c, d, u)
                if v != 0:
                    M[Id[pth[:i] + (d,) + pth[i + 1:]], col] += v
        if normalize:
            M = M * self.scal(A, B, u)
        return M, nt

    def transfer(self, types, thetas):
        types = tuple(types); N = len(types)
        P0, I0 = self.per(types)
        T = np.eye(len(P0), dtype=complex); cur = types
        for j in range(1, N):
            Xm, cur = self.X(cur, j, uof(thetas[0] - thetas[j]))
            T = Xm @ T
        Pc, _ = self.per(cur)
        Cm = np.zeros((len(P0), len(Pc)))
        for k, q in enumerate(Pc):
            Cm[I0[(q[N - 1],) + q[:N - 1] + (q[N - 1],)], k] = 1
        return Cm @ T
