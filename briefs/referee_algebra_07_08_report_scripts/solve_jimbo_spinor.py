import numpy as np
import scipy.linalg
import itertools

from spinor_jimbo import make_spinor_rep, Qs

def solve_jimbo(n, z_val=1.4+0.3j, Q=Qs, conv='7'):
    rep, cartan, qi = make_spinor_rep(n, Q)
    dim = 2**n
    dim2 = dim * dim
    
    e_x, f_x, k_x = rep(z_val)
    e_y, f_y, k_y = rep(1.0)
    
    def coprod(e1, f1, k1, e2, f2, k2):
        k1_inv = np.linalg.inv(k1)
        k2_inv = np.linalg.inv(k2)
        if conv == '7':
            de = np.kron(e1, k2) + np.kron(np.eye(dim), e2)
            df = np.kron(f1, np.eye(dim)) + np.kron(k1_inv, f2)
            dk = np.kron(k1, k2)
        else:
            de = np.kron(e1, np.eye(dim)) + np.kron(k1, e2)
            df = np.kron(f1, k2_inv) + np.kron(np.eye(dim), f2)
            dk = np.kron(k1, k2)
        return de, df, dk

    basis = list(itertools.product([1, -1], repeat=n))
    weights = []
    for a in range(dim):
        for b in range(dim):
            w = tuple(basis[a][i] + basis[b][i] for i in range(n))
            weights.append(w)
            
    weight_classes = {}
    for idx, w in enumerate(weights):
        weight_classes.setdefault(w, []).append(idx)
        
    num_unknowns = sum(len(indices)**2 for indices in weight_classes.values())
    print(f"\n================ n={n} ================")
    print(f"dim={dim}, dim^2={dim2}. Weight-preserving unknowns = {num_unknowns} (out of {dim2**2})")
    
    allowed_pairs = []
    for indices in weight_classes.values():
        for r in indices:
            for c in indices:
                allowed_pairs.append((r, c))
    pair_map = {p: i for i, p in enumerate(allowed_pairs)}
    
    eq_rows = []
    for i in range(n+1):
        de_xy, df_xy, _ = coprod(e_x[i], f_x[i], k_x[i], e_y[i], f_y[i], k_y[i])
        de_yx, df_yx, _ = coprod(e_y[i], f_y[i], k_y[i], e_x[i], f_x[i], k_x[i])
        
        for D_xy, D_yx in ((de_xy, de_yx), (df_xy, df_yx)):
            for r in range(dim2):
                for c in range(dim2):
                    row_dict = {}
                    col_c = np.nonzero(D_xy[:, c])[0]
                    for m in col_c:
                        if (r, m) in pair_map:
                            var = pair_map[(r, m)]
                            row_dict[var] = row_dict.get(var, 0.0) + D_xy[m, c]
                    row_r = np.nonzero(D_yx[r, :])[0]
                    for m in row_r:
                        if (m, c) in pair_map:
                            var = pair_map[(m, c)]
                            row_dict[var] = row_dict.get(var, 0.0) - D_yx[r, m]
                            
                    if row_dict:
                        eq_rows.append(row_dict)
                        
    M = np.zeros((len(eq_rows), num_unknowns), complex)
    for r_idx, row_dict in enumerate(eq_rows):
        for col_idx, val in row_dict.items():
            M[r_idx, col_idx] = val
            
    u, s, vh = np.linalg.svd(M)
    null_mask = s < 1e-9 * s[0]
    num_null = np.sum(null_mask) + (M.shape[1] - len(s))
    print(f"SVD: smallest singular value = {s[-1]:.2e}, second smallest = {s[-2]:.2e}")
    print(f"Dimension of solution space = {num_null}")
    
    sol = vh[-1, :]
    Rcheck = np.zeros((dim2, dim2), complex)
    for (r, c), idx in pair_map.items():
        Rcheck[r, c] = sol[idx]
        
    top_idx = 0
    Rcheck = Rcheck / Rcheck[top_idx, top_idx]
    
    evals = np.linalg.eigvals(Rcheck)
    
    clusters = []
    for ev in evals:
        matched = False
        for cl in clusters:
            if abs(ev - cl['val']) < 1e-5:
                cl['mult'] += 1
                matched = True
                break
        if not matched:
            clusters.append({'val': ev, 'mult': 1})
            
    clusters.sort(key=lambda cl: cl['mult'], reverse=True)
    print(f"Distinct eigenvalues of Rcheck(z) at z={z_val}:")
    for cl in clusters:
        print(f"  val = {cl['val']:.6f}, multiplicity = {cl['mult']}")
        
    print("\nComparing with theoretical formula: rho_{\\Lambda^{n-j}} = prod_{i=1}^j <i/2>_{(-1)^i}")
    for q_candidate, q_label in [(Q, "Kac's Q (short root length^2=2)"), (Q**0.5, "Q^(1/2) (book's long root length^2=2)")]:
        print(f"\n--- Testing with q = {q_label} ---")
        cur = 1.0 + 0j
        theo_vals = [cur]
        for j in range(1, n+1):
            sgn = (-1)**j
            a = j / 2.0
            factor = (1.0 + sgn * z_val * q_candidate**(2*a)) / (z_val + sgn * q_candidate**(2*a))
            cur *= factor
            theo_vals.append(cur)
        for j, val in enumerate(theo_vals):
            matching_cl = [cl for cl in clusters if abs(val - cl['val']) < 1e-5]
            mult = matching_cl[0]['mult'] if matching_cl else 0
            err = min(abs(val - cl['val']) for cl in clusters)
            print(f"  j={j} (Lambda^{n-j}): theo = {val:.6f}, matched mult = {mult}, residual = {err:.2e}")

solve_jimbo(2)
solve_jimbo(3)
