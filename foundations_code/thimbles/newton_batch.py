"""Batched damped Newton search for critical points of the chain action."""
import numpy as np

def hess_batch(ch, U):
    B, L, r = U.shape
    E = np.exp(1j*np.einsum('jr,bxr->bxj', ch.alpha, U)) * ch.kac
    H = np.zeros((B, L, r, L, r), dtype=complex)
    blk = np.einsum('bxj,ja,jc->bxac', E, ch.alpha, ch.alpha)*ch.kappa
    I = np.eye(r)
    for x in range(L):
        H[:, x, :, x, :] = 2*I + blk[:, x]
        if x+1 < L:
            H[:, x, :, x+1, :] = -I
            H[:, x+1, :, x, :] = -I
    return H.reshape(B, L*r, L*r)

def newton_batch(ch, U0, maxit=80, tol=1e-12, maxstep=1.5):
    U = U0.copy()
    B = U.shape[0]
    done = np.zeros(B, bool); bad = np.zeros(B, bool)
    with np.errstate(all='ignore'):
        for it in range(maxit):
            act = np.where(~done & ~bad)[0]
            if len(act) == 0: break
            G = ch.grad(U[act]).reshape(len(act), -1)
            nG = np.linalg.norm(G, axis=1)
            conv = nG < tol
            done[act[conv]] = True
            bad[act[~np.isfinite(nG)]] = True
            act2 = act[~conv & np.isfinite(nG)]
            if len(act2) == 0: continue
            H = hess_batch(ch, U[act2])
            G2 = G[~conv & np.isfinite(nG)]
            try:
                d = np.linalg.solve(H, -G2[..., None])[..., 0]
            except np.linalg.LinAlgError:
                d = np.zeros_like(G2)
                for i in range(len(act2)):
                    try: d[i] = np.linalg.solve(H[i], -G2[i])
                    except np.linalg.LinAlgError: bad[act2[i]] = True
            nd = np.linalg.norm(d, axis=1)
            fac = np.minimum(1.0, maxstep/np.maximum(nd, 1e-300))
            U[act2] += (fac[:, None]*d).reshape(len(act2), ch.L, ch.r)
            bad[act2[np.max(np.abs(U[act2].imag), axis=(1, 2)) > 25]] = True
    G = ch.grad(U).reshape(B, -1)
    ok = np.isfinite(np.linalg.norm(G, axis=1)) & (np.linalg.norm(G, axis=1) < 1e-10)
    return U, ok

def find_cps(ch, nstart=20000, seed=0, re_spread=np.pi, im_spread=1.5, batch=2000, tol_dup=1e-6, extra=None):
    from model import lin_interp
    rng = np.random.default_rng(seed)
    base = lin_interp(ch)
    found = []; counts = []
    chunks = []
    for s in range(0, nstart, batch):
        m = min(batch, nstart - s)
        chunks.append(base[None] + rng.uniform(-re_spread, re_spread, (m,)+base.shape) + 1j*rng.normal(0, im_spread, (m,)+base.shape))
    if extra is not None:
        chunks.append(np.asarray(extra))
    for U0 in chunks:
        U, ok = newton_batch(ch, U0)
        for u in U[ok]:
            hit = False
            for i, v in enumerate(found):
                if np.max(np.abs(v - u)) < tol_dup:
                    counts[i] += 1; hit = True; break
            if not hit:
                found.append(u); counts.append(1)
    return found, counts
