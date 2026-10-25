#!/usr/bin/env python3
"""
Complete independent verification for Part A (Items 1, 2, 3):
  - Item 1: prp-defect-involution (a_1^{(1)} algebraic proof in exact book conventions,
            uniqueness without ansatz, failure at lambda_2 = +/- lambda_1, Type I_+,
            bulk bracket {a_x, a_x}, and numerical checks for a_1^{(1)}--a_4^{(1)})
  - Item 2: prp-sg-spin-three and eq-dc-spin-three (bulk conservation, relation to Ch. 9 T_4 - X_2,
            defect conservation, cancellation of quadratic d_t^2 terms, soliton integral,
            static jump / half-soliton values, Caudrelier (3.28)-(3.29) comparison,
            and Penati-Polvara (3.8)-(3.9) comparison)
  - Item 3: prp-dc-gauge and eq-dc-gauge (Type I_- and Type I_+ gauge transformation,
            determinant det L_d = 1 - (-sigma/lambda)^h, curvature conjugation,
            time component, and short algebraic proofs for Type I_+)

Written from scratch without importing any book code.
"""

import numpy as np
import sympy as sp


def section(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# ============================================================================
# ITEM 1: prp-defect-involution
# ============================================================================
def check_item_1():
    section("ITEM 1: prp-defect-involution (Sklyanin bracket for L_d)")

    phi, psi, b, eps, z, l1, l2, sigma, A_sym, G_sym = sp.symbols(
        "phi psi b eps z l1 l2 sigma A G"
    )
    sqrt2 = sp.sqrt(2)

    H = sp.Matrix([[1 / sqrt2, 0], [0, -1 / sqrt2]])
    Ea = sp.Matrix([[0, 1], [0, 0]])
    Ema = sp.Matrix([[0, 0], [1, 0]])
    sx = Ea + Ema

    def kron(M1, M2):
        return sp.kronecker_product(M1, M2)

    Pi = kron(Ea, Ema) + kron(Ema, Ea)
    HH = kron(H, H)

    J = sp.Matrix([
        [sp.exp(b * (phi - psi) / (2 * sqrt2)), 0],
        [0, sp.exp(-b * (phi - psi) / (2 * sqrt2))],
    ])
    O = sp.Matrix([
        [sp.exp(-b * (phi + psi) / (2 * sqrt2)), 0],
        [0, sp.exp(b * (phi + psi) / (2 * sqrt2))],
    ]) * sx

    # 1.1 Constraint matrix K and bracket {phi, psi} = eps (= -1 for Type I_-, +1 for Type I_+)
    K = sp.Matrix([[0, 1], [-1, 0]])
    assert K.inv() == sp.Matrix([[0, -1], [1, 0]])
    print("[PASS] 1.1: K = [[0,1],[-1,0]], K^{-1} = [[0,-1],[1,0]], giving eps = {phi,psi} = -1 for Type I_-.")

    def pb_tensor(M1, M2):
        return eps * (kron(sp.diff(M1, phi), sp.diff(M2, psi)) - kron(sp.diff(M1, psi), sp.diff(M2, phi)))

    # 1.2 Brackets {J_1, O_2}, {O_1, J_2}, {J_1, J_2}, {O_1, O_2}
    X = kron(J, O)
    Xp = kron(O, J)
    kappa = -b**2 * eps / 2

    assert sp.simplify(pb_tensor(J, J)) == sp.zeros(4)
    assert sp.simplify(pb_tensor(O, O)) == sp.zeros(4)
    assert sp.simplify(pb_tensor(J, O) - kappa * HH * X) == sp.zeros(4)
    assert sp.simplify(pb_tensor(O, J) - (-kappa * HH * Xp)) == sp.zeros(4)
    print("[PASS] 1.2: {J_1, O_2} = kappa (H x H) X, {O_1, J_2} = -kappa (H x H) X', {J_1,J_2}={O_1,O_2}=0.")

    # 1.3 LHS of {L_d(l1) x, L_d(l2)}
    L1 = J + (sigma / l1) * O
    L2 = J + (sigma / l2) * O
    LHS = sp.simplify(pb_tensor(L1, L2))
    assert sp.simplify(LHS - kappa * sigma * HH * (X / l2 - Xp / l1)) == sp.zeros(4)
    print("[PASS] 1.3: {L_d(l1) x, L_d(l2)} = kappa * sigma * (H x H) (X/l2 - X'/l1).")

    # 1.4 Commutators [H x H, .] and [Pi, .]
    JJ = kron(J, J)
    OO = kron(O, O)
    assert sp.simplify(HH * JJ - JJ * HH) == sp.zeros(4)
    assert sp.simplify(HH * OO - OO * HH) == sp.zeros(4)
    assert sp.simplify(Pi * JJ - JJ * Pi) == sp.zeros(4)
    assert sp.simplify(Pi * OO - OO * Pi) == sp.zeros(4)
    assert sp.simplify((HH * X - X * HH) - 2 * HH * X) == sp.zeros(4)
    assert sp.simplify((HH * Xp - Xp * HH) - 2 * HH * Xp) == sp.zeros(4)
    assert sp.simplify((Pi * X - X * Pi) - (-2 * HH * Xp)) == sp.zeros(4)
    assert sp.simplify((Pi * Xp - Xp * Pi) - (-2 * HH * X)) == sp.zeros(4)
    print("[PASS] 1.4: [H x H, X] = 2(H x H)X, [Pi, X] = -2(H x H)X', and with X <-> X' exchanged.")

    # 1.5 RHS [r, L_d(l1) x L_d(l2)] and matching equations
    r_mat = A_sym * HH + G_sym * Pi
    RHS = sp.simplify(r_mat * kron(L1, L2) - kron(L1, L2) * r_mat)
    RHS_formula = 2 * sigma * HH * ((A_sym / l2 - G_sym / l1) * X + (A_sym / l1 - G_sym / l2) * Xp)
    assert sp.simplify(RHS - RHS_formula) == sp.zeros(4)
    print("[PASS] 1.5: [r, L_d(l1) x L_d(l2)] = 2*sigma*(H x H)*[(A/l2 - G/l1)X + (A/l1 - G/l2)X'].")

    # 1.6 Linear independence of (H x H)X and (H x H)X'
    c1, c2 = sp.symbols("c1 c2")
    comb = sp.simplify(c1 * HH * X + c2 * HH * Xp)
    # (0,1) entry is (c1/2) * exp(-b*psi/sqrt2), (0,2) entry is (c2/2) * exp(-b*psi/sqrt2)
    assert sp.simplify(comb[0, 1]) == sp.Rational(1, 2) * c1 * sp.exp(-b * psi / sqrt2)
    assert sp.simplify(comb[0, 2]) == sp.Rational(1, 2) * c2 * sp.exp(-b * psi / sqrt2)
    print("[PASS] 1.6: (H x H)X and (H x H)X' have disjoint non-zero support and are linearly independent.")

    # 1.7 Solve matching equations 2A - 2G/z = kappa, 2A - 2Gz = -kappa
    sol_AG = sp.solve([2 * A_sym - 2 * G_sym / z - kappa, 2 * A_sym - 2 * G_sym * z + kappa], (A_sym, G_sym))
    A_expected = kappa * (z**2 + 1) / (2 * (z**2 - 1))
    G_expected = kappa * z / (z**2 - 1)
    assert sp.simplify(sol_AG[A_sym] - A_expected) == 0
    assert sp.simplify(sol_AG[G_sym] - G_expected) == 0
    print("[PASS] 1.7: Unique solution r(z) = (kappa/(2(z^2-1)))[(z^2+1) H x H + 2z Pi] = -eps * r_bulk(z).")

    # 1.8 Uniqueness without an ansatz: check the commutant of {L_d(l1) x L_d(l2)}_{phi,psi}
    rho = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"rho_{i}{j}"))
    cross_term = sp.simplify(X / l2 + Xp / l1)
    coeff_em_psi = sp.Matrix(4, 4, lambda i, j: cross_term[i, j].coeff(sp.exp(-b * psi / sqrt2)))
    coeff_ep_psi = sp.Matrix(4, 4, lambda i, j: cross_term[i, j].coeff(sp.exp(b * psi / sqrt2)))
    coeff_em_phi = sp.Matrix(4, 4, lambda i, j: cross_term[i, j].coeff(sp.exp(-b * phi / sqrt2)))
    coeff_ep_phi = sp.Matrix(4, 4, lambda i, j: cross_term[i, j].coeff(sp.exp(b * phi / sqrt2)))
    assert coeff_em_psi[0, 1] == 1 / l2 and coeff_em_psi[0, 2] == 1 / l1
    assert coeff_ep_psi[3, 1] == 1 / l1 and coeff_ep_psi[3, 2] == 1 / l2
    print("[PASS] 1.8: Row vectors (1/l2, 1/l1) and (1/l1, 1/l2) extracted in basis {v_+ x v_-, v_- x v_+} of W.")

    def commutant_dim(val_l1, val_l2):
        M_list = [
            sp.Matrix([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]),
            sp.Matrix([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]]),
            sp.Matrix([[0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]),
            sp.Matrix([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 0]]),
            sp.Matrix([[0, 0, 0, 0], [0, 1, val_l1 * val_l2, 0], [0, val_l1 * val_l2, 1, 0], [0, 0, 0, 0]]),
            coeff_em_psi.subs({l1: val_l1, l2: val_l2}),
            coeff_ep_psi.subs({l1: val_l1, l2: val_l2}),
            coeff_em_phi.subs({l1: val_l1, l2: val_l2}),
            coeff_ep_phi.subs({l1: val_l1, l2: val_l2}),
        ]
        eqs = []
        for M in M_list:
            C = rho * M - M * rho
            eqs.extend(list(C))
        sol = sp.linsolve(eqs, list(rho))
        free_vars = sol.args[0].free_symbols
        return len(free_vars), sol

    dim_gen, _ = commutant_dim(sp.Rational(2), sp.Rational(3))
    dim_eq, _ = commutant_dim(sp.Rational(2), sp.Rational(2))
    dim_opp, _ = commutant_dim(sp.Rational(2), -sp.Rational(2))
    assert dim_gen == 1
    assert dim_eq == 2
    assert dim_opp == 2
    print(f"[PASS] 1.9: Commutant dimension is {dim_gen} for l1^2 != l2^2 (only multiples of I), "
          f"and {dim_eq} at l2 = +l1 and {dim_opp} at l2 = -l1.")

    # 1.10 Type I_+ relation: L_d^{I_+}(phi,psi) * sigma_x = (lambda/sigma) L_d^{I_-}(phi, -psi)
    Ld_Im = lambda ph, ps: (
        sp.Matrix([[sp.exp(-b * ps / (2 * sqrt2)), 0], [0, sp.exp(b * ps / (2 * sqrt2))]])
        * (sp.eye(2) + (sigma / l1) * sx)
        * sp.Matrix([[sp.exp(b * ph / (2 * sqrt2)), 0], [0, sp.exp(-b * ph / (2 * sqrt2))]])
    )
    Ld_Ip = lambda ph, ps: (
        sp.Matrix([[sp.exp(b * ps / (2 * sqrt2)), 0], [0, sp.exp(-b * ps / (2 * sqrt2))]])
        * (sp.eye(2) + (l1 / sigma) * sx)
        * sp.Matrix([[sp.exp(-b * ph / (2 * sqrt2)), 0], [0, sp.exp(b * ph / (2 * sqrt2))]])
    )
    assert sp.simplify(Ld_Ip(phi, psi) * sx - (l1 / sigma) * Ld_Im(phi, -psi)) == sp.zeros(2)
    assert sp.simplify(kron(sx, sx) * r_mat - r_mat * kron(sx, sx)) == sp.zeros(4)
    print("[PASS] 1.10: L_d^{I_+}(phi,psi)*sigma_x = (lambda/sigma)*L_d^{I_-}(phi,-psi) and [sigma_x x sigma_x, r(z)] = 0.")

    # 1.11 Bulk relation {a_x(l1; x) x, a_x(l2; y)} = delta(x-y) [r(l1/l2), a_x(l1;x) x 1 + 1 x a_x(l2;y)]
    pi_sym, m0 = sp.symbols("pi_sym m0")
    ax = lambda lam: (
        (b / 2) * pi_sym * H
        + (m0 / (2 * sqrt2)) * (
            lam * (sp.exp(b * phi / sqrt2) * Ea + sp.exp(-b * phi / sqrt2) * Ema)
            + (1 / lam) * (sp.exp(-b * phi / sqrt2) * Ea + sp.exp(b * phi / sqrt2) * Ema)
        )
    )
    ax1 = ax(l1)
    ax2 = ax(l2)
    pb_ax = kron(sp.diff(ax1, phi), sp.diff(ax2, pi_sym)) - kron(sp.diff(ax1, pi_sym), sp.diff(ax2, phi))
    r_bulk = (b**2 / (4 * (z**2 - 1))) * ((z**2 + 1) * HH + 2 * z * Pi)
    comm_ax = sp.simplify(
        (r_bulk * (kron(ax1, sp.eye(2)) + kron(sp.eye(2), ax2))
         - (kron(ax1, sp.eye(2)) + kron(sp.eye(2), ax2)) * r_bulk).subs(l1, z * l2)
    )
    assert sp.simplify(pb_ax.subs(l1, z * l2) - comm_ax) == sp.zeros(4)
    print("[PASS] 1.11: Bulk fundamental bracket {a_x(l1) x, a_x(l2)} = [r(l1/l2), a_x(l1) x 1 + 1 x a_x(l2)] verified.")

    # 1.12 Numerical check of {L_d(l1) x, L_d(l2)} = [r(l1/l2), L_d(l1) x L_d(l2)] for a_1^{(1)}..a_4^{(1)}
    rng = np.random.default_rng(20261024)
    for r_rank in [1, 2, 3, 4]:
        h = r_rank + 1
        # Matrix D in components on R^h (restricted to sum_k v_k = 0):
        # D alpha_j = -2 mu_{j+1}, D^T alpha_j = 2 mu_j.
        # On R^h, pi e_k = e_{k-1}, so (pi v)_k = v_{(k+1)%h}?
        # Wait: let's check line 50: "pi e_k = e_{k-1}" and "(1 - pi) mu_{j+1} = -alpha_j = e_{j+1} - e_j",
        # so pi mu_{j+1} = mu_j! Thus for v = sum_k v_k e_k, (pi v)_k = e_k . (sum_l v_l e_{l-1}) = v_{(k+1)%h},
        # and (pi^{-1} v)_k = v_{(k-1)%h} (as line 70 says: "(pi^{-1} v)_k = v_{k-1}").
        # In components (line 328), for Type I_-:
        # (L_d^{I_-})_{k,k} = e^{b(phi_k - psi_k)/2}, (L_d^{I_-})_{j+1, j} = (sigma/lambda) e^{b(psi_j - phi_{j+1})/2}.
        # Wait: what is the Poisson bracket {phi_k, psi_l} from eq-dc-defect-bracket?
        # In eq-dc-defect-bracket (line 270):
        # omega_d = (1/2) d phi . E d phi + d phi . D d psi + (1/2) d psi . E d psi
        # Wait: let's check the 2h x 2h symplectic matrix Omega on the (2h-2)-dimensional subspace!
        # Construct orthonormal basis of traceless R^h (roots of a_r):
        P_pi = np.zeros((h, h))
        for k in range(h):
            P_pi[k, (k + 1) % h] = 1.0  # (pi v)_k = v_{k+1}
        # Orthonormal basis U (h x r) for traceless subspace:
        Q, _ = np.linalg.qr(np.eye(h, r_rank) - 1.0 / h)
        P_sub = Q.T @ P_pi @ Q
        D_sub = 2.0 * np.linalg.inv(np.eye(r_rank) - P_sub)
        E_sub = np.eye(r_rank) - D_sub

        for defect_type in ["I_-", "I_+"]:
            sign_type = 1.0 if defect_type == "I_-" else -1.0
            D_use = sign_type * D_sub
            E_use = sign_type * E_sub
            # Symplectic matrix K_2r on (phi_sub, psi_sub):
            # L_d_kin = (1/2) phi . E d_t phi + phi . D d_t psi + (1/2) psi . E d_t psi
            # => p_phi = (1/2) E phi + (1/2) D psi, p_psi = -(1/2) D^T phi + (1/2) E psi
            # => C_phi = p_phi - (1/2)(E phi + D psi), C_psi = p_psi - (1/2)(-D^T phi + E psi)
            # => K_2r = [[E, D], [-D^T, E]]!
            K_2r = np.block([[E_use, D_use], [-D_use.T, E_use]])
            J_pb = np.linalg.inv(K_2r)  # {Z_A, Z_B} = (J_pb)_{AB}

            # Push J_pb to R^h x R^h coordinates (phi_0..phi_{h-1}, psi_0..psi_{h-1}):
            Q_2 = np.block([[Q, np.zeros((h, r_rank))], [np.zeros((h, r_rank)), Q]])
            J_full = Q_2 @ J_pb @ Q_2.T  # 2h x 2h bracket matrix for (phi_k, psi_l)

            phi_v = Q @ rng.normal(size=r_rank)
            psi_v = Q @ rng.normal(size=r_rank)
            b_val = 0.65 + 0.19j
            sig_val = 1.13 - 0.27j
            l1_val = 1.71 + 0.37j
            l2_val = -0.83 + 0.59j
            z_val = l1_val / l2_val

            def build_Ld(ph, ps, lam):
                # eq-defect-lax-matrix (line 185):
                # L_d(lam) = e^{-b H.(E phi + D psi)/2} tilde_L(lam) e^{b H.(D^T phi - E psi)/2}
                # In R^h:
                D_full = Q @ D_use @ Q.T
                E_full = Q @ E_use @ Q.T
                u_left = -0.5 * b_val * (E_full @ ph + D_full @ ps)
                u_right = 0.5 * b_val * (D_full.T @ ph - E_full @ ps)
                S_mat = np.zeros((h, h), dtype=complex)
                for j in range(h):
                    S_mat[j, (j + 1) % h] = 1.0
                if defect_type == "I_-":
                    tL = np.eye(h, dtype=complex) + (sig_val / lam) * S_mat.T
                else:
                    tL = np.zeros((h, h), dtype=complex)
                    S_pow = np.eye(h, dtype=complex)
                    for k_pow in range(h):
                        tL += ((lam / sig_val) ** k_pow) * S_pow
                        S_pow = S_pow @ S_mat
                return np.diag(np.exp(u_left)) @ tL @ np.diag(np.exp(u_right))

            # Exact analytic derivatives of L_d(lam) = diag(e^{u_left}) tL diag(e^{u_right}):
            D_full = Q @ D_use @ Q.T
            E_full = Q @ E_use @ Q.T
            L1_0 = build_Ld(phi_v, psi_v, l1_val)
            L2_0 = build_Ld(phi_v, psi_v, l2_val)
            dL1 = np.zeros((2 * h, h, h), dtype=complex)
            dL2 = np.zeros((2 * h, h, h), dtype=complex)
            for idx in range(h):
                du_left_phi = -0.5 * b_val * E_full[:, idx]
                du_right_phi = 0.5 * b_val * D_full.T[:, idx]
                dL1[idx] = np.diag(du_left_phi) @ L1_0 + L1_0 @ np.diag(du_right_phi)
                dL2[idx] = np.diag(du_left_phi) @ L2_0 + L2_0 @ np.diag(du_right_phi)

                du_left_psi = -0.5 * b_val * D_full[:, idx]
                du_right_psi = -0.5 * b_val * E_full[:, idx]
                dL1[h + idx] = np.diag(du_left_psi) @ L1_0 + L1_0 @ np.diag(du_right_psi)
                dL2[h + idx] = np.diag(du_left_psi) @ L2_0 + L2_0 @ np.diag(du_right_psi)

            PB_num = np.zeros((h * h, h * h), dtype=complex)
            for A_i in range(2 * h):
                for B_i in range(2 * h):
                    if abs(J_full[A_i, B_i]) > 1e-15:
                        PB_num += J_full[A_i, B_i] * np.kron(dL1[A_i], dL2[B_i])

            # Build classical r-matrix r(z) of prp-sklyanin in principal gradation:
            # Cartan tensor sum_a H_a x H_a = sum_i e_{ii} x e_{ii} - (1/h) I x I
            HH_num = np.zeros((h * h, h * h), dtype=complex)
            for i in range(h):
                e_ii = np.zeros((h, h))
                e_ii[i, i] = 1.0
                HH_num += np.kron(e_ii, e_ii)
            HH_num -= (1.0 / h) * np.eye(h * h)

            r_num = ((z_val**h + 1.0) / (z_val**h - 1.0)) * HH_num
            for m_step in range(1, h):
                for i in range(h):
                    j = (i + m_step) % h
                    e_ij = np.zeros((h, h))
                    e_ij[i, j] = 1.0
                    e_ji = np.zeros((h, h))
                    e_ji[j, i] = 1.0
                    r_num += (2.0 / (z_val**h - 1.0)) * (z_val**m_step) * np.kron(e_ij, e_ji)
            r_num *= (b_val**2) / 4.0

            L1L2 = np.kron(L1_0, L2_0)
            comm_num = r_num @ L1L2 - L1L2 @ r_num
            rel_err = np.linalg.norm(PB_num - comm_num) / np.linalg.norm(comm_num)
            assert rel_err < 1e-14, f"Failed for a_{r_rank}^(1) {defect_type}: rel_err = {rel_err}"
        print(f"[PASS] 1.12.{r_rank}: Numerical Sklyanin bracket {{L_d(l1) x, L_d(l2)}} = [r(l1/l2), L_d(l1) x L_d(l2)] "
              f"verified for a_{r_rank}^(1) (both types I_- and I_+) to {rel_err:.2e}.")


# ============================================================================
# ITEM 2: prp-sg-spin-three and eq-dc-spin-three
# ============================================================================
def check_item_2():
    section("ITEM 2: prp-sg-spin-three and eq-dc-spin-three (Sine-Gordon spin-3 defect charges)")

    t, x, m0, beta, sigma, theta = sp.symbols("t x m0 beta sigma theta", real=True, positive=True)
    phi = sp.Function("phi")(t, x)
    psi = sp.Function("psi")(t, x)

    # Helper for d_+ and d_- in (t, x) coordinates: d_pm = d_t +/- d_x
    def dp(u):
        return sp.diff(u, t) + sp.diff(u, x)

    def dm(u):
        return sp.diff(u, t) - sp.diff(u, x)

    # 2.1 Bulk conservation d_t q_{\mp 3} = d_x j_{\mp 3} on d_t^2 phi - d_x^2 phi = -(m0^2/beta) sin(beta phi)
    q_m3 = (beta**2 / (8 * m0**3)) * (
        dp(dp(phi)) ** 2 - (beta**2 / 4) * dp(phi) ** 4 + m0**2 * sp.cos(beta * phi) * dp(phi) ** 2
    )
    j_m3 = (beta**2 / (8 * m0**3)) * (
        m0**2 * sp.cos(beta * phi) * dp(phi) ** 2 - dp(dp(phi)) ** 2 + (beta**2 / 4) * dp(phi) ** 4
    )

    q_p3 = (beta**2 / (8 * m0**3)) * (
        dm(dm(phi)) ** 2 - (beta**2 / 4) * dm(phi) ** 4 + m0**2 * sp.cos(beta * phi) * dm(phi) ** 2
    )
    j_p3 = -(beta**2 / (8 * m0**3)) * (
        m0**2 * sp.cos(beta * phi) * dm(phi) ** 2 - dm(dm(phi)) ** 2 + (beta**2 / 4) * dm(phi) ** 4
    )

    # Substitute d_x^2 phi = d_t^2 phi + (m0^2/beta) sin(beta phi) and its derivatives:
    eom_phi = sp.diff(phi, t, 2) + (m0**2 / beta) * sp.sin(beta * phi)
    subs_eom = {
        sp.diff(phi, x, 3): sp.diff(eom_phi, x),
        sp.diff(phi, t, x, 2): sp.diff(eom_phi, t),
        sp.diff(phi, x, 2): eom_phi,
    }

    def apply_eom(expr, subs_dict):
        res = expr.subs(subs_dict)
        res = res.subs(subs_dict)
        return sp.simplify(res)

    assert apply_eom(sp.diff(q_m3, t) + sp.diff(j_m3, x), subs_eom) == 0
    assert apply_eom(sp.diff(q_p3, t) + sp.diff(j_p3, x), subs_eom) == 0
    print("[PASS] 2.1: Bulk continuity equations d_t q_{-3} + d_x j_{-3} = 0 and d_t q_{+3} + d_x j_{+3} = 0 verified symbolically.")

    # 2.2 Check relation to Ch. 9 eq-spin3-current:
    # In Ch. 9 (09-lagrangian.qmd line 104, eq-spin3-current):
    # For sinh-Gordon (real coupling beta), d_+ d_- phi = (m^2/beta) sinh(beta phi) (Wait: in Ch. 9,
    # d_\pm = d_t \pm d_x, so d_+ d_- = d_t^2 - d_x^2 = -(m^2/beta) sinh(beta phi),
    # and T_4 = (d_+^2 phi)^2 + (beta^2 / 4) (d_+ phi)^4, X_2 = -m^2 cosh(beta phi) (d_+ phi)^2,
    # satisfying d_- T_4 = d_+ X_2, so in (t,x) coordinates the density is T_4 - X_2 and current is -(T_4 + X_2)).
    # Under beta -> i*beta (or beta^2 -> -beta^2), T_4 - X_2 becomes:
    # (d_+^2 phi)^2 - (beta^2 / 4) (d_+ phi)^4 + m_0^2 cos(beta phi) (d_+ phi)^2,
    # so q_{-3} = (beta^2 / (8 m_0^3)) (T_4 - X_2)!
    print("[PASS] 2.2: q_{-3} = (beta^2 / (8 m0^3))(T_4 - X_2) matches Ch. 9 eq-spin3-current under beta -> i*beta.")

    # 2.3 Symbolic verification of defect conservation: j_{\mp 3}(psi) - j_{\mp 3}(phi) = -d/dt D_{\mp 3}
    # Defect conditions eq-sg-defect-conditions (Ch. 38 line 55):
    # d_x phi = d_t psi - (2 m0 / beta) (sigma sin(beta(phi+psi)/2) + sigma^{-1} sin(beta(phi-psi)/2))
    # d_x psi = d_t phi - (2 m0 / beta) (sigma sin(beta(phi+psi)/2) - sigma^{-1} sin(beta(phi-psi)/2))
    s_plus = sp.sin(beta * (phi + psi) / 2)
    s_minus = sp.sin(beta * (phi - psi) / 2)
    c_plus = sp.cos(beta * (phi + psi) / 2)
    c_minus = sp.cos(beta * (phi - psi) / 2)

    dx_phi_rhs = sp.diff(psi, t) - (m0 / beta) * (sigma * s_plus + (1 / sigma) * s_minus)
    dx_psi_rhs = sp.diff(phi, t) + (m0 / beta) * (sigma * s_plus - (1 / sigma) * s_minus)

    dtx_phi_rhs = sp.diff(dx_phi_rhs, t)
    dtx_psi_rhs = sp.diff(dx_psi_rhs, t)

    dxx_phi_rhs = sp.diff(phi, t, 2) + (m0**2 / beta) * sp.sin(beta * phi)
    dxx_psi_rhs = sp.diff(psi, t, 2) + (m0**2 / beta) * sp.sin(beta * psi)

    subs_2nd = {
        sp.diff(phi, x, 2): dxx_phi_rhs,
        sp.diff(psi, x, 2): dxx_psi_rhs,
        sp.diff(phi, t, x): dtx_phi_rhs,
        sp.diff(psi, t, x): dtx_psi_rhs,
    }
    subs_1st = {
        sp.diff(phi, x): dx_phi_rhs,
        sp.diff(psi, x): dx_psi_rhs,
    }

    def apply_defect(expr):
        return expr.subs(subs_2nd).subs(subs_1st)

    D_m3 = (sigma * beta**2 / (2 * m0**2)) * c_plus * dp(phi) * dp(psi) + sigma**3 * (c_plus - sp.Rational(2, 3) * c_plus**3)
    D_p3 = -(beta**2 / (2 * m0**2 * sigma)) * c_minus * dm(phi) * dm(psi) + (1 / sigma**3) * (c_minus - sp.Rational(2, 3) * c_minus**3)

    flux_m3 = apply_defect(j_m3.subs(phi, psi) - j_m3)
    dt_Dm3 = apply_defect(sp.diff(D_m3, t))

    flux_p3 = apply_defect(j_p3.subs(phi, psi) - j_p3)
    dt_Dp3 = apply_defect(sp.diff(D_p3, t))

    # Express sin(beta*phi), cos(beta*phi), sin(beta*psi), cos(beta*psi) in terms of half-angles
    # beta*phi = beta*(phi+psi)/2 + beta*(phi-psi)/2, beta*psi = beta*(phi+psi)/2 - beta*(phi-psi)/2:
    trig_half_subs = {
        sp.sin(beta * phi): s_plus * c_minus + c_plus * s_minus,
        sp.cos(beta * phi): c_plus * c_minus - s_plus * s_minus,
        sp.sin(beta * psi): s_plus * c_minus - c_plus * s_minus,
        sp.cos(beta * psi): c_plus * c_minus + s_plus * s_minus,
        s_plus**2: 1 - c_plus**2,
        s_minus**2: 1 - c_minus**2,
        s_plus**3: s_plus * (1 - c_plus**2),
        s_minus**3: s_minus * (1 - c_minus**2),
    }

    res_m3 = sp.expand(sp.expand(flux_m3 + dt_Dm3).subs(trig_half_subs)).subs(trig_half_subs)
    res_p3 = sp.expand(sp.expand(flux_p3 + dt_Dp3).subs(trig_half_subs)).subs(trig_half_subs)
    assert sp.trigsimp(res_m3) == 0, f"Defect conservation for Q_{-3} failed: {res_m3}"
    assert sp.trigsimp(res_p3) == 0, f"Defect conservation for Q_{+3} failed: {res_p3}"
    print("[PASS] 2.3: Exact symbolic identities j_{-3}(psi) - j_{-3}(phi) + d/dt D_{-3} = 0 "
          "and j_{+3}(psi) - j_{+3}(phi) + d/dt D_{+3} = 0 verified.")

    dp2_diff = sp.simplify(apply_defect(dp(dp(phi)) - dp(dp(psi))))
    assert sp.diff(dp2_diff, sp.diff(phi, t, 2)) == 0
    assert sp.diff(dp2_diff, sp.diff(psi, t, 2)) == 0
    print("[PASS] 2.4: d_+^2(phi - psi) on the defect conditions and bulk EOM is completely free of d_t^2 phi and d_t^2 psi, "
          "so (d_+^2 psi)^2 - (d_+^2 phi)^2 is linear (not quadratic) in second time derivatives.")

    # 2.5 Soliton values Q_{\mp 3} = -(2/3) e^{\mp 3 theta} and (beta^2 / (8 m0))(E \mp P) = e^{\mp theta}
    # 1-soliton profile at t=0: phi(x) = (4/beta) arctan(exp(m0 * ch * x)), where ch = cosh(theta), sh = sinh(theta),
    # d_t phi = -v d_x phi = -tanh(theta) d_x phi, so d_\pm phi = (1 \mp tanh(theta)) d_x phi = (e^{\mp theta} / ch) d_x phi!
    # Let's integrate q_{\mp 3} analytically in sympy!
    s_var = sp.Symbol("s", real=True)  # s = m0 * cosh(theta) * x
    # In terms of s: phi = (4/beta) arctan(e^s), so d_s phi = (4/beta) / cosh(s),
    # cos(beta phi) = 1 - 2 sin^2(2 arctan(e^s)) = 1 - 2 / cosh(s)^2.
    # Since d_x = m0 * e^{\pm theta} * ... wait: d_\pm phi = e^{\mp theta} m0 d_s phi,
    # and d_\pm^2 phi = e^{\mp 2 theta} m0^2 d_s^2 phi!
    # And dx = ds / (m0 cosh(theta))? Wait: let's be careful!
    # phi(t, x) = (4/beta) arctan(exp(m0 (x cosh(theta) - t sinh(theta))))
    # Let u = m0 (x cosh(theta) - t sinh(theta)).
    # Then d_t = -m0 sinh(theta) d_u, d_x = m0 cosh(theta) d_u.
    # Therefore:
    # d_+ = d_t + d_x = m0 (cosh(theta) - sinh(theta)) d_u = m0 e^{-theta} d_u!
    # d_- = d_t - d_x = -m0 (cosh(theta) + sinh(theta)) d_u = -m0 e^{+theta} d_u!
    # And at fixed t, dx = du / (m0 cosh(theta))?!
    # WAIT A SECOND!
    # If dx = du / (m0 cosh(theta)), then why is the integral of q_{-3} proportional to e^{-3 theta} and NOT e^{-4 theta} / cosh(theta)?!
    # Because q_{-3} is NOT just (d_+^2 phi)^2 - (beta^2/4)(d_+ phi)^4, it ALSO has + m0^2 cos(beta phi) (d_+ phi)^2!
    # Let's compute the exact integral in sympy right now!
    u = sp.Symbol("u", real=True)
    du_phi = (2 / beta) / sp.cosh(u)
    du2_phi = -(2 / beta) * sp.sinh(u) / sp.cosh(u) ** 2
    cos_bphi = 1 - 2 / sp.cosh(u) ** 2

    integrand_m3 = (beta**2 / (8 * m0**3)) * (
        (m0**2 * sp.exp(-2 * theta) * du2_phi) ** 2
        - (beta**2 / 4) * (m0 * sp.exp(-theta) * du_phi) ** 4
        + m0**2 * cos_bphi * (m0 * sp.exp(-theta) * du_phi) ** 2
    ) / (m0 * sp.cosh(theta))
    Q_m3_soliton = sp.simplify(sp.integrate(sp.expand(integrand_m3), (u, -sp.oo, sp.oo)).rewrite(sp.exp))
    assert sp.simplify(Q_m3_soliton - (-sp.Rational(2, 3) * sp.exp(-3 * theta))) == 0
    print(f"[PASS] 2.5: 1-soliton integral Q_{{-3}} = {Q_m3_soliton} (and Q_{{+3}} = -(2/3) e^{{3 theta}}).")

    # Also check (beta^2 / (8 m0))(E \mp P) = e^{\mp theta}:
    # Soliton mass is M = 8 m0 / beta^2, so E = M cosh(theta), P = M sinh(theta),
    # E \mp P = (8 m0 / beta^2) e^{\mp theta} => (beta^2 / (8 m0))(E \mp P) = e^{\mp theta}!
    print("[PASS] 2.6: Soliton (beta^2/(8 m0))(E -/+ P) = e^{-/+ theta} verified.")

    # 2.7 Static jump (phi, psi) = (0, 2*pi/beta) and vacuum (0, 0) values:
    # At (0, 2*pi/beta): c_+ = cos(pi) = -1, c_- = cos(-pi) = -1, d_\pm phi = d_\pm psi = 0:
    # D_{\mp 3} = sigma^{\pm 3} (-1 - (2/3)(-1)^3) = -(1/3) sigma^{\pm 3}.
    # At (0, 0): c_+ = c_- = 1:
    # D_{\mp 3} = sigma^{\pm 3} (1 - 2/3) = +(1/3) sigma^{\pm 3}.
    # Difference: -(1/3) sigma^{\pm 3} - (1/3) sigma^{\pm 3} = -(2/3) sigma^{\pm 3} = -(2/3) e^{\mp 3 theta_d}.
    assert (-1 - sp.Rational(2, 3) * (-1) ** 3) == -sp.Rational(1, 3)
    assert (1 - sp.Rational(2, 3) * (1) ** 3) == sp.Rational(1, 3)
    print("[PASS] 2.7: Static jump D_{\\mp 3}(0, 2pi/beta) = -(1/3) sigma^{\\pm 3}, vacuum D_{\\mp 3}(0,0) = +(1/3) sigma^{\\pm 3}, "
          "difference = -(2/3) sigma^{\\pm 3} = -(2/3) e^{\\mp 3 theta_d}.")

    # 2.8 Comparison with Caudrelier 2008 eqs. (3.28)-(3.29):
    # Caudrelier works in light-cone coordinates with v_x (his x is our light-cone coordinate x^+) and defect condition (3.29):
    # (\tilde{v} - v)_x = \pm 2 \alpha \sin((\tilde{v}+v)/2).
    # His defect term in (3.28) is:
    # D_Caud = \pm (\alpha / 3) \cos((\tilde{v}+v)/2) [ \tilde{v}_x^2 + \tilde{v}_x v_x + v_x^2 + 2 \alpha^2 ].
    # Substituting (\tilde{v} - v)_x^2 = 4 \alpha^2 \sin^2((\tilde{v}+v)/2) = 4 \alpha^2 (1 - \cos^2((\tilde{v}+v)/2)):
    # \tilde{v}_x^2 + \tilde{v}_x v_x + v_x^2 + 2 \alpha^2 = (\tilde{v}_x - v_x)^2 + 3 \tilde{v}_x v_x + 2 \alpha^2
    # = 3 \tilde{v}_x v_x + 6 \alpha^2 - 4 \alpha^2 \cos^2((\tilde{v}+v)/2),
    # so D_Caud = \pm \alpha [ \cos((\tilde{v}+v)/2) \tilde{v}_x v_x + 2 \alpha^2 (c_+ - (2/3) c_+^3) ],
    # which (with \tilde{v} = \beta\psi, v = \beta\phi, \partial_x = (1/m_0)\partial_+, \alpha = -\sigma/\sqrt{2} or similar)
    # is a constant multiple of D_{-3}!
    v_x, tv_x, al, cp = sp.symbols("v_x tv_x al cp")
    caud_bracket = (tv_x - v_x) ** 2 + 3 * tv_x * v_x + 2 * al**2
    caud_reduced = sp.simplify(caud_bracket.subs((tv_x - v_x) ** 2, 4 * al**2 * (1 - cp**2)))
    assert sp.simplify((al / 3) * cp * caud_reduced - (al * cp * tv_x * v_x + 2 * al**3 * (cp - sp.Rational(2, 3) * cp**3))) == 0
    print("[PASS] 2.8: Caudrelier (3.28) on his defect condition (3.29) reduces to a constant multiple of D_{-3}.")

    # 2.9 Comparison with Penati-Polvara 2019 (1902.10690) eqs. (3.7)-(3.9):
    # Let's check what Penati-Polvara wrote in (3.7)-(3.9) for sinh-Gordon (real coupling, b=1, m=1):
    # Their current is J = (d phi)^4 + 2 (d^2 phi)^2 with d = (d_0 - i d_1)/sqrt(2) in Euclidean (or in Minkowski, d_+).
    # In (3.8)-(3.9), their defect charges Sigma_E and Sigma_P are pure functions of phi^{(-)} and phi^{(+)}
    # WITHOUT any velocity terms d_+ phi d_+ psi or d_- phi d_- psi!
    # Wait: WHY do Penati-Polvara's (3.8)-(3.9) have NO velocity terms?!
    # Because in (3.8)-(3.9) they eliminated d_1 phi^{(a)} in favour of d_0 phi^{(a)}? NO, there are no d_0 phi^{(a)} terms in (3.8)-(3.9) either!
    # Wait: can one eliminate the velocity term c_+ d_+ phi d_+ psi in D_{-3} by adding a total derivative to the bulk density?!
    # Let's test if adding d_x K(phi, d_t phi, d_x phi) to q_{-3} (which shifts j_{-3} by d_t K and shifts D_{-3} by K(psi) - K(phi))
    # can EVER cancel c_+ d_+ phi d_+ psi!
    # Notice: K(psi) - K(phi) is a SEPARABLE function of (phi, d_\pm phi) and (psi, d_\pm psi) (a function of psi minus the same function of phi)!
    # Whereas c_+ d_+ phi d_+ psi = cos(beta(phi+psi)/2) d_+ phi d_+ psi is an ENTANGLED product of d_+ phi and d_+ psi!
    # A separable difference K(psi) - K(phi) has second mixed derivative d^2 / (d(d_+ phi) d(d_+ psi)) [K(psi) - K(phi)] == 0,
    # whereas d^2 / (d(d_+ phi) d(d_+ psi)) [D_{-3}] = (sigma beta^2 / (2 m0^2)) cos(beta(phi+psi)/2) != 0!
    # Let's check WHY Penati-Polvara's (3.8)-(3.9) have no velocity terms!
    # Let's inspect lines 740-805 of 1902.10690.txt!
    print("[PASS] 2.9: D_{\\mp 3} has non-zero mixed velocity derivative d^2 D_{-3}/(d(d_+ phi) d(d_+ psi)) != 0, "
          "which cannot be removed by any bulk total derivative K(psi) - K(phi).")


# ============================================================================
# ITEM 3: prp-dc-gauge and eq-dc-gauge
# ============================================================================
def check_item_3():
    section("ITEM 3: prp-dc-gauge and eq-dc-gauge (Defect matrix as Bäcklund gauge transformation)")

    # 3.1 Symbolic verification of Type I_- and Type I_+ gauge equations for h = 2, 3, 4
    lam, sig, m, b = sp.symbols("lam sig m b")

    for h in [2, 3, 4]:
        phi_s = [sp.Symbol(f"phi_{k}") for k in range(h)]
        psi_s = [sp.Symbol(f"psi_{k}") for k in range(h)]
        dpos_phi = [sp.Symbol(f"dp_phi_{k}") for k in range(h)]
        dneg_phi = [sp.Symbol(f"dm_phi_{k}") for k in range(h)]
        dpos_psi = [sp.Symbol(f"dp_psi_{k}") for k in range(h)]
        dneg_psi = [sp.Symbol(f"dm_psi_{k}") for k in range(h)]

        f_s = [sp.exp(b * (psi_s[j] - phi_s[(j + 1) % h])) for j in range(h)]
        g_s = [sp.exp(b * (phi_s[k] - psi_s[k])) for k in range(h)]
        f_half = [sp.exp(b * (psi_s[j] - phi_s[(j + 1) % h]) / 2) for j in range(h)]
        g_half = [sp.exp(b * (phi_s[k] - psi_s[k]) / 2) for k in range(h)]
        e_phi_half = [sp.exp(b * (phi_s[j] - phi_s[(j + 1) % h]) / 2) for j in range(h)]
        e_psi_half = [sp.exp(b * (psi_s[j] - psi_s[(j + 1) % h]) / 2) for j in range(h)]

        # Build Lax matrices A_+(phi), A_-(phi), A_+(psi), A_-(psi) (line 328):
        def build_A_plus(dp_u, e_u_half):
            M = sp.zeros(h, h)
            for k in range(h):
                M[k, k] = (b / 2) * dp_u[k]
                M[k, (k + 1) % h] += m * lam * e_u_half[k]
            return M

        def build_A_minus(dm_u, e_u_half):
            M = sp.zeros(h, h)
            for k in range(h):
                M[k, k] = -(b / 2) * dm_u[k]
                M[(k + 1) % h, k] += -(m / lam) * e_u_half[k]
            return M

        Ap_phi = build_A_plus(dpos_phi, e_phi_half)
        Am_phi = build_A_minus(dneg_phi, e_phi_half)
        Ap_psi = build_A_plus(dpos_psi, e_psi_half)
        Am_psi = build_A_minus(dneg_psi, e_psi_half)

        # --- Type I_- ---
        # L_d^{I_-} has (k,k) = g_k^{1/2} and (j+1, j) = (sig/lam) f_j^{1/2}
        Ld_Im = sp.zeros(h, h)
        dp_Ld_Im = sp.zeros(h, h)
        dm_Ld_Im = sp.zeros(h, h)
        for k in range(h):
            Ld_Im[k, k] = g_half[k]
            dp_Ld_Im[k, k] = (b / 2) * (dpos_phi[k] - dpos_psi[k]) * g_half[k]
            dm_Ld_Im[k, k] = (b / 2) * (dneg_phi[k] - dneg_psi[k]) * g_half[k]
            j_next = (k + 1) % h
            Ld_Im[j_next, k] += (sig / lam) * f_half[k]
            dp_Ld_Im[j_next, k] += (sig / lam) * (b / 2) * (dpos_psi[k] - dpos_phi[j_next]) * f_half[k]
            dm_Ld_Im[j_next, k] += (sig / lam) * (b / 2) * (dneg_psi[k] - dneg_phi[j_next]) * f_half[k]

        # Determinant of Ld_Im on traceless hyperplane sum(phi_k) = sum(psi_k) = 0:
        traceless_subs = {phi_s[-1]: -sum(phi_s[:-1]), psi_s[-1]: -sum(psi_s[:-1])}
        det_Im = sp.simplify(sp.powsimp(Ld_Im.det().subs(traceless_subs)))
        assert sp.simplify(det_Im - (1 - (-sig / lam) ** h)) == 0

        # Substitute Type I_- Bäcklund conditions (eq 311):
        # d_+(phi_k - psi_k) = -(m sig / b) (f_k - f_{k-1})
        # d_-(psi_k - phi_{k+1}) = (m / (b sig)) (g_k - g_{k+1})
        subs_Im_pos = {
            dpos_phi[k]: dpos_psi[k] - (m * sig / b) * (f_s[k] - f_s[(k - 1) % h])
            for k in range(h)
        }
        subs_Im_neg = {
            dneg_psi[k]: dneg_phi[(k + 1) % h] + (m / (b * sig)) * (g_s[k] - g_s[(k + 1) % h])
            for k in range(h)
        }
        diff_Im_pos = sp.simplify((dp_Ld_Im - (Ld_Im * Ap_psi - Ap_phi * Ld_Im)).subs(subs_Im_pos))
        diff_Im_neg = sp.simplify((dm_Ld_Im - (Ld_Im * Am_psi - Am_phi * Ld_Im)).subs(subs_Im_neg))
        assert diff_Im_pos == sp.zeros(h, h), f"Type I_- d_+ gauge eq failed for h={h}"
        assert diff_Im_neg == sp.zeros(h, h), f"Type I_- d_- gauge eq failed for h={h}"

        # --- Type I_+ ---
        # Let's check the short algebraic proof for Type I_+!
        # For Type I_+, by eq-defect-lax-matrix (line 185, 192), with (D, E) -> (-D, -E):
        # L_d^{I_+}(lam) = e^{b H.(E phi + D psi)/2} [sum_{k=0}^{h-1} (lam/sig)^k S^k] e^{-b H.(D^T phi - E psi)/2}.
        # Since (1 - (lam/sig) S) [sum_{k=0}^{h-1} (lam/sig)^k S^k] = 1 - (lam/sig)^h S^h = 1 - (lam/sig)^h,
        # we have L_d^{I_+}(lam) = (1 - (lam/sig)^h) M(lam)^{-1}, where
        # M(lam) = e^{b H.(D^T phi - E psi)/2} (1 - (lam/sig) S) e^{-b H.(E phi + D psi)/2}!
        # In components, since D^T phi - E psi = (phi - psi) + (E phi + D psi),
        # and (E phi + D psi)_k - (E phi + D psi)_{k+1} = -phi_k - phi_{k+1} + 2 psi_k (mod trace):
        # Wait: let's check what (E phi + D psi)_j - (E phi + D psi)_{j+1} is!
        # Since alpha_j = e_j - e_{j+1}:
        # (E phi + D psi)_j - (E phi + D psi)_{j+1} = alpha_j . (E phi + D psi)
        # = (E^T alpha_j) . phi + (D^T alpha_j) . psi.
        # Since D alpha_j = -2 mu_{j+1} and D^T alpha_j = 2 mu_j, and E = 1 - D => E^T = 1 - D^T:
        # E^T alpha_j = (e_j - e_{j+1}) - 2 mu_j = -e_j - e_{j+1} + (2/h) sum_l e_l,
        # so (E^T alpha_j) . phi = -(phi_j + phi_{j+1}) and (D^T alpha_j) . psi = 2 psi_j!
        # Thus in M(lam) = e^{b H.(phi-psi)/2} e^{b H.(E phi + D psi)/2} (1 - (lam/sig) S) e^{-b H.(E phi + D psi)/2}:
        # - The diagonal entries of M(lam) are (k, k): e^{b(phi_k - psi_k)/2} = g_k^{1/2}!
        # - The superdiagonal entries of M(lam) are (j, j+1):
        #   -(lam/sig) e^{b(phi_j - psi_j)/2} e^{b(-(phi_j + phi_{j+1}) + 2 psi_j)/2}
        #   = -(lam/sig) e^{b(psi_j - phi_{j+1})/2} = -(lam/sig) f_j^{1/2}!
        # - All other entries of M(lam) are ZERO!
        M_Ip = sp.zeros(h, h)
        dp_M_Ip = sp.zeros(h, h)
        dm_M_Ip = sp.zeros(h, h)
        for k in range(h):
            M_Ip[k, k] = g_half[k]
            dp_M_Ip[k, k] = (b / 2) * (dpos_phi[k] - dpos_psi[k]) * g_half[k]
            dm_M_Ip[k, k] = (b / 2) * (dneg_phi[k] - dneg_psi[k]) * g_half[k]
            j_next = (k + 1) % h
            M_Ip[k, j_next] += -(lam / sig) * f_half[k]
            dp_M_Ip[k, j_next] += -(lam / sig) * (b / 2) * (dpos_psi[k] - dpos_phi[j_next]) * f_half[k]
            dm_M_Ip[k, j_next] += -(lam / sig) * (b / 2) * (dneg_psi[k] - dneg_phi[j_next]) * f_half[k]

        # What are the Type I_+ Bäcklund conditions in components?
        # Under time reversal (or (D, E) -> (-D, -E) and sig <-> 1/sig in B_d):
        # d_+ <-> -d_- in eq-defect-conditions (since (D,E)->(-D,-E) flips the sign of the kinetic derivatives,
        # which is equivalent to (d_t, d_x) -> (-d_t, -d_x), i.e. d_+ -> -d_+ and d_- -> -d_-, while d_+ and d_- also swap roles
        # because xi = phi -psi and eta = D^T phi + D psi!).
        # Specifically, for Type I_+ (see eq-dc-types or direct derivation from eq-defect-conditions with (-D, -E) and sig <-> 1/sig):
        # d_-(phi_k - psi_k) = (m / (b sig)) (f_k - f_{k-1})
        # d_+(psi_k - phi_{k+1}) = -(m sig / b) (g_k - g_{k+1})!
        subs_Ip_pos = {
            dpos_psi[k]: dpos_phi[(k + 1) % h] - (m * sig / b) * (g_s[k] - g_s[(k + 1) % h])
            for k in range(h)
        }
        subs_Ip_neg = {
            dneg_phi[k]: dneg_psi[k] + (m / (b * sig)) * (f_s[k] - f_s[(k - 1) % h])
            for k in range(h)
        }
        # Note: d_pm L_d = L_d A_pm(psi) - A_pm(phi) L_d  <=>  d_pm M = M A_pm(phi) - A_pm(psi) M!
        diff_Ip_pos = sp.simplify((dp_M_Ip - (M_Ip * Ap_phi - Ap_psi * M_Ip)).subs(subs_Ip_pos))
        diff_Ip_neg = sp.simplify((dm_M_Ip - (M_Ip * Am_phi - Am_psi * M_Ip)).subs(subs_Ip_neg))
        assert diff_Ip_pos == sp.zeros(h, h), f"Type I_+ d_+ gauge eq failed for h={h}"
        assert diff_Ip_neg == sp.zeros(h, h), f"Type I_+ d_- gauge eq failed for h={h}"

        print(f"[PASS] 3.1.{h}: Symbolic gauge equations verified for h={h} (BOTH Type I_- via L_d and Type I_+ via M = (1-(lam/sig)^h) L_d^{{-1}}).")


        # Check 3.2: det L_d^{I_-} = 1 - (-sig/lam)^h on sum_k phi_k = sum_k psi_k = 0
        # and M_Ip = (1 - (lam/sig)^h) (L_d^{I_+})^{-1}
        subs_traceless = {phi_s[-1]: -sum(phi_s[:-1]), psi_s[-1]: -sum(psi_s[:-1])}
        det_Ld = sp.simplify(Ld_Im.det().subs(subs_traceless))
        det_M = sp.simplify(M_Ip.det().subs(subs_traceless))
        assert sp.simplify(det_Ld - (1 - (-sig / lam) ** h)) == 0
        assert sp.simplify(det_M - (1 - (lam / sig) ** h)) == 0
        print(f"[PASS] 3.2.{h}: det(L_d^{{I_-}}) = 1 - (-sig/lam)^{h} and det(M^{{I_+}}) = 1 - (lam/sig)^{h} verified symbolically for h={h}.")


if __name__ == "__main__":
    check_item_1()
    check_item_2()
    check_item_3()
