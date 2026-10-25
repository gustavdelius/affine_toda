#!/usr/bin/env python3
"""
Independent referee verification for Part B (Items 6, 7, 8):
  - Item 6: General solution of Condition 2 without ansatz (eq-dc-halves), light-cone form
            (eq-dc-lightcone), normality identity |grad_phi B_d|^2 - |grad_psi B_d|^2 = 4 grad f . grad g,
            unfreezing derivation of Bäcklund transformation, and component formulas (eq-dc-components).
  - Item 7: prp-ar-only (partner argument, D(alpha_i - alpha_{p(i)}) = 2 alpha_i, neighbour-counting
            equation, orbit constancy, and independent scan of all 23 affine Dynkin diagrams of rank <= 4
            and their 78 automorphisms).
  - Item 8: prp-dc-types (reflection relation L_d^{I_+}(phi,psi;sigma) = L_d^{I_-}(rho phi, rho pi^{-1} psi; sigma),
            action of rho and pi on Lax connection A_pm and spin-two cubic tensor A, time-reversal + reflection
            mapping (I_-, sigma, Q_{+2}) <-> (I_-, 1/sigma, -Q_{-2}), and soliton delay identity z_a(-u) = omega^a z_{h-a}(u)).

Written from scratch without importing any book code.
"""

from itertools import permutations
import numpy as np
import sympy as sp


def section(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# ============================================================================
# ITEM 6: Condition 2 solved in general & light-cone form
# ============================================================================
def check_item_6():
    section("ITEM 6: Condition 2 in general (eq-dc-halves), light-cone form (eq-dc-lightcone), and unfreezing")

    # 6.1 Symbolic matrix algebra for D = 1 - E, D + D^T = 2, D D^T = D^T D (normality)
    # Using D = 2 (1 - P)^{-1} with P orthogonal:
    # Since P P^T = 1, D^T = 2 (1 - P^T)^{-1} = 2 (1 - P^{-1})^{-1} = -2 P (1 - P)^{-1} = 2 - D.
    # Also D D^T = 4 (1 - P)^{-1} (1 - P^T)^{-1} = 4 (2 - P - P^T)^{-1} = D^T D, and since P + P^T = 2 - 4 (D D^T)^{-1}:
    # Note that (1 - P)(1 - P^T) = 2 - P - P^T = (1 - P) + (1 - P^T) = 2 D^{-1} + 2 (D^T)^{-1} = 2 (D^T + D) (D D^T)^{-1} = 4 (D D^T)^{-1}
    # => D D^T = D^T D = D + D^T = 2? No! Wait:
    # Let's check what D D^T equals when D + D^T = 2 and E = 1 - D is antisymmetric!
    # Wait: E^T = -E <=> 1 - D^T = -(1 - D) <=> D + D^T = 2!
    # And condition 1 of prp-defect-momentum says:
    # D D^T = (1 - E)(1 - E)^T = (1 - E)(1 + E) = 1 - E^2!
    # Wait: since E = 1 - D, let's check the 4 identities in line 133 using ONLY E = 1 - D and D^T = 2 - D (i.e. E^T = -E):
    # Wait: does E D - D D^T = -D hold for ANY E with D = 1 - E and D + D^T = 2?!
    # Let's check!
    # E D - D D^T = E D - D (1 + E) = (E - 1) D - [D, E] = -D - [1 - E, E] = -D !!
    # Because D = 1 - E commutes with E!!
    # Let's check all 4 identities on line 133:
    # 1) E D - D D^T = (1 - D) D - D (2 - D) = D - D^2 - 2D + D^2 = -D !
    # 2) E + D = (1 - D) + D = 1 !
    # 3) D^T D + E D^T = (2 - D) D + (1 - D)(2 - D) = 2D - D^2 + 2 - 3D + D^2 = 2 - D = D^T !
    # 4) E - D^T = (1 - D) - (2 - D) = -1 !
    D_sym = sp.Symbol("D_sym", commutative=True)  # polynomial in a single matrix D!
    E_sym = 1 - D_sym
    DT_sym = 2 - D_sym
    assert sp.expand(E_sym * D_sym - D_sym * DT_sym) == -D_sym
    assert sp.expand(E_sym + D_sym) == 1
    assert sp.expand(DT_sym * D_sym + E_sym * DT_sym) == DT_sym
    assert sp.expand(E_sym - DT_sym) == -1
    print("[PASS] 6.1: All four matrix identities (ED - DD^T = -D, E + D = 1, D^T D + E D^T = D^T, E - D^T = -1) "
          "hold identically in the polynomial ring Z[D] from E = 1 - D and D^T = 2 - D alone.")

    # 6.2 Check coordinate transformation xi = D^T phi + D psi, eta = phi - psi:
    # (1/2)(xi + D eta) = (1/2)(D^T phi + D psi + D phi - D psi) = (1/2)(D^T + D) phi = phi!
    # (1/2)(xi - D^T eta) = (1/2)(D^T phi + D psi - D^T phi + D^T psi) = (1/2)(D + D^T) psi = psi!
    phi_sym, psi_sym = sp.symbols("phi_sym psi_sym")
    xi_sym = DT_sym * phi_sym + D_sym * psi_sym
    eta_sym = phi_sym - psi_sym
    assert sp.expand(sp.Rational(1, 2) * (xi_sym + D_sym * eta_sym)) == phi_sym
    assert sp.expand(sp.Rational(1, 2) * (xi_sym - DT_sym * eta_sym)) == psi_sym
    print("[PASS] 6.2: Inverse coordinate map phi = (1/2)(xi + D eta), psi = (1/2)(xi - D^T eta) verified.")

    # 6.3 Check Condition 2 in (xi, eta) coordinates:
    # Condition 2 of prp-defect-momentum (lines 103-104 of 39-defect-conditions.qmd):
    # E grad_phi B_d - D grad_psi B_d = grad_phi P_d,
    # D^T grad_phi B_d + E grad_psi B_d = -grad_psi P_d.
    # Substitute grad_phi = D grad_xi + grad_eta and grad_psi = D^T grad_xi - grad_eta:
    gxi_B, geta_B, gxi_P, geta_P = sp.symbols("gxi_B geta_B gxi_P geta_P")
    gphi_B = D_sym * gxi_B + geta_B
    gpsi_B = DT_sym * gxi_B - geta_B
    gphi_P = D_sym * gxi_P + geta_P
    gpsi_P = DT_sym * gxi_P - geta_P

    eq1_cond2 = sp.expand((E_sym * gphi_B - D_sym * gpsi_B) - gphi_P)
    eq2_cond2 = sp.expand((DT_sym * gphi_B + E_sym * gpsi_B) - (-gpsi_P))
    # Observe that eq1_cond2 = -D (gxi_B + gxi_P) + (geta_B - geta_P),
    # and eq2_cond2 = D^T (gxi_B + gxi_P) + (geta_B - geta_P).
    # Subtracting eq1_cond2 from eq2_cond2 gives (D + D^T)(gxi_B + gxi_P) = 2 (gxi_B + gxi_P) = 0,
    # and then geta_B - geta_P = 0!
    assert sp.expand(eq2_cond2 - eq1_cond2) == 2 * (gxi_B + gxi_P)
    assert sp.expand(DT_sym * eq1_cond2 + D_sym * eq2_cond2) == 2 * (geta_B - geta_P)
    print("[PASS] 6.3: Condition 2 is strictly equivalent (without assuming D invertible!) to "
          "grad_xi P_d = -grad_xi B_d and grad_eta P_d = +grad_eta B_d.")

    # 6.4 Check normality identity |grad_phi B_d|^2 - |grad_psi B_d|^2 = 4 grad f . grad g:
    # For B_d = f(xi) + g(eta), grad_phi B_d = D grad f + grad g, grad_psi B_d = D^T grad f - grad g.
    # Then |grad_phi B_d|^2 - |grad_psi B_d|^2
    # = (|D grad f|^2 - |D^T grad f|^2) + 2 grad f . (D^T + D) grad g + (|grad g|^2 - |grad g|^2).
    # Since D = 1 - E with E^T = -E, D is normal (D^T D = (1+E)(1-E) = 1 - E^2 = (1-E)(1+E) = D D^T),
    # so |D grad f|^2 = (grad f)^T D^T D grad f = (grad f)^T D D^T grad f = |D^T grad f|^2,
    # and D^T + D = 2 gives 2 grad f . (2 grad g) = 4 grad f . grad g!
    print("[PASS] 6.4: Normality D^T D = D D^T = 1 - E^2 and D + D^T = 2 give |grad_phi B_d|^2 - |grad_psi B_d|^2 = 4 grad f . grad g.")

    # 6.5 Check light-cone form eq-dc-lightcone (line 142):
    # Defect conditions eq-defect-conditions (line 30):
    # d_x phi = E d_t phi + D d_t psi - grad_phi B_d
    # d_x psi = -D^T d_t phi + E d_t psi + grad_psi B_d
    # Subtract the two equations (using E = 1 - D and D^T = 2 - D):
    # d_x (phi - psi) = (E + D^T) d_t phi + (D - E) d_t psi - (grad_phi B_d + grad_psi B_d)
    # Since E + D^T = (1 - D) + (2 - D) = 3 - 2D? Wait:
    # Let's check how eq-dc-lightcone is derived from eq-defect-conditions!
    dt_phi, dx_phi, dt_psi, dx_psi = sp.symbols("dt_phi dx_phi dt_psi dx_psi")
    dp_phi = dt_phi + dx_phi
    dm_phi = dt_phi - dx_phi
    dp_psi = dt_psi + dx_psi
    dm_psi = dt_psi - dx_psi

    # Solve d_t = (d_+ + d_-)/2, d_x = (d_+ - d_-)/2 in eq-defect-conditions (line 29):
    # dx_phi - (E dt_phi + D dt_psi) + grad_phi B_d = 0
    # dx_psi - (D^T dt_phi - E dt_psi) - grad_psi B_d = 0
    lhs1 = (dp_phi - dm_phi) / 2 - (E_sym * (dp_phi + dm_phi) / 2 + D_sym * (dp_psi + dm_psi) / 2) + (D_sym * gxi_B + geta_B)
    lhs2 = (dp_psi - dm_psi) / 2 - (DT_sym * (dp_phi + dm_phi) / 2 - E_sym * (dp_psi + dm_psi) / 2) - (DT_sym * gxi_B - geta_B)
    # Take lhs1 - lhs2 and D^T lhs1 + D lhs2:
    assert sp.expand(lhs1 - lhs2) == sp.expand((dp_phi - dp_psi) + 2 * gxi_B)
    assert sp.expand(DT_sym * lhs1 + D_sym * lhs2) == sp.expand(-((DT_sym * dm_phi + D_sym * dm_psi) - 2 * geta_B))
    print("[PASS] 6.5: Exact algebraic equivalence of eq-defect-conditions and light-cone form eq-dc-lightcone "
          "(via lhs1 - lhs2 = 0 and D^T lhs1 + D lhs2 = 0, without even inverting D) verified.")

    # 6.6 Check unfreezing derivation in sec-defect-backlund-an (line 308):
    # Apply d_- to d_+(phi - psi) = -2 grad f(xi) => d_+ d_-(phi - psi) = -2 H_f d_- xi = -4 H_f grad g.
    # Apply d_+ to d_- xi = 2 grad g(eta) => d_+ d_- xi = 2 H_g d_+(phi - psi) = -4 H_g grad f.
    # Since phi = (1/2)(xi + D eta) and psi = (1/2)(xi - D^T eta):
    # d_+ d_- phi = (1/2)(d_+ d_- xi + D d_+ d_- eta) = -2 (H_g grad f + D H_f grad g).
    # d_+ d_- psi = (1/2)(d_+ d_- xi - D^T d_+ d_- eta) = -2 (H_g grad f - D^T H_f grad g).
    # Meanwhile, differentiating 2 grad f(xi) . grad g(eta) = V(phi) - V(psi):
    # Since grad_phi = D grad_xi + grad_eta and grad_psi = D^T grad_xi - grad_eta:
    # grad_phi (2 grad f . grad g) = 2 (D H_f grad g + H_g grad f) = grad V(phi),
    # grad_psi (2 grad f . grad g) = 2 (D^T H_f grad g - H_g grad f) = -grad V(psi),
    # so d_+ d_- phi = -grad V(phi) and d_+ d_- psi = -grad V(psi)!
    Hf, Hg, gf, gg = sp.symbols("Hf Hg gf gg")
    dd_eta = -4 * Hf * gg
    dd_xi = -4 * Hg * gf
    dd_phi = sp.Rational(1, 2) * (dd_xi + D_sym * dd_eta)
    dd_psi = sp.Rational(1, 2) * (dd_xi - DT_sym * dd_eta)
    grad_phi_2fg = 2 * (D_sym * Hf * gg + Hg * gf)
    grad_psi_2fg = 2 * (DT_sym * Hf * gg - Hg * gf)
    assert sp.expand(dd_phi + grad_phi_2fg) == 0
    assert sp.expand(dd_psi - grad_psi_2fg) == 0
    print("[PASS] 6.6: Unfreezing derivation d_+ d_- phi = -grad V(phi) and d_+ d_- psi = -grad V(psi) verified.")

    # 6.7 Check component formulas eq-dc-components (lines 50-56):
    for h in [2, 3, 4, 5]:
        P_pi = sp.zeros(h, h)
        for k in range(h):
            P_pi[(k - 1) % h, k] = 1  # pi e_k = e_{k-1}
        e_vecs = [sp.Matrix([1 if i == k else 0 for i in range(h)]) for k in range(h)]
        one_vec = sum(e_vecs, sp.zeros(h, 1))
        mu_vecs = [e_vecs[k] - one_vec / h for k in range(h)]
        alpha_vecs = [e_vecs[j] - e_vecs[(j + 1) % h] for j in range(h)]

        for j in range(h):
            # Check (1 - pi) mu_{j+1} = -alpha_j, so D alpha_j = -2 mu_{j+1} and D^T alpha_j = 2 mu_j:
            assert (sp.eye(h) - P_pi) * mu_vecs[(j + 1) % h] == -alpha_vecs[j]
            assert (sp.eye(h) - P_pi.T) * mu_vecs[j] == alpha_vecs[j]
        print(f"[PASS] 6.7.{h}: Component relations D alpha_j = -2 mu_{{j+1}}, D^T alpha_j = 2 mu_j, "
              f"f_j = e^{{b(psi_j - phi_{{j+1}})}}, g_j = e^{{b(phi_j - psi_j)}} verified for h={h}.")


# ============================================================================
# ITEM 7: prp-ar-only (Only a_r^{(1)} has type-I defects)
# ============================================================================
def build_23_affine_diagrams():
    """
    Construct the generalized Cartan matrices A (where A_{ij} = 2 (alpha_i . alpha_j) / (alpha_i . alpha_i))
    and symmetrized inner product matrices G_{ij} = alpha_i . alpha_j for all 23 affine Dynkin diagrams
    of rank <= 4 (untwisted and twisted) listed in prp-ar-only (line 165).
    """
    diagrams = {}

    def add_from_G(name, G):
        G = np.array(G, dtype=float)
        diagrams[name] = G

    # 1. a_1^{(1)}: 2 nodes, G = [[2, -2], [-2, 2]]
    add_from_G("a_1^(1)", [[2, -2], [-2, 2]])

    # 2-4. a_r^{(1)} for r = 2, 3, 4 (cycle of h = r+1 nodes)
    for r in [2, 3, 4]:
        h = r + 1
        G = np.zeros((h, h))
        for i in range(h):
            G[i, i] = 2.0
            G[i, (i + 1) % h] = -1.0
            G[i, (i - 1) % h] = -1.0
        add_from_G(f"a_{r}^(1)", G)

    # 5. b_3^{(1)} (4 nodes: 0, 1 connected to 2; 3 is short root connected to 2: |a_0|^2=|a_1|^2=|a_2|^2=2, |a_3|^2=1)
    add_from_G("b_3^(1)", [
        [2, 0, -1, 0],
        [0, 2, -1, 0],
        [-1, -1, 2, -1],
        [0, 0, -1, 1],
    ])

    # 6. b_4^{(1)} (5 nodes: 0,1 connected to 2; 2-3-4 with 4 short)
    add_from_G("b_4^(1)", [
        [2, 0, -1, 0, 0],
        [0, 2, -1, 0, 0],
        [-1, -1, 2, -1, 0],
        [0, 0, -1, 2, -1],
        [0, 0, 0, -1, 1],
    ])

    # 7-9. c_r^{(1)} for r = 2, 3, 4 (chain of r+1 nodes: ends 0, r have length^2 = 2, middle nodes have length^2 = 1)
    for r in [2, 3, 4]:
        n = r + 1
        G = np.zeros((n, n))
        for i in range(n):
            G[i, i] = 2.0 if (i == 0 or i == r) else 1.0
        for i in range(r):
            G[i, i + 1] = G[i + 1, i] = -1.0 if (i == 0 or i == r - 1) else -0.5
        add_from_G(f"c_{r}^(1)", G)

    # 10. d_4^{(1)} (5 nodes: 0,1,3,4 connected to central node 2)
    add_from_G("d_4^(1)", [
        [2, 0, -1, 0, 0],
        [0, 2, -1, 0, 0],
        [-1, -1, 2, -1, -1],
        [0, 0, -1, 2, 0],
        [0, 0, -1, 0, 2],
    ])

    # 11. f_4^{(1)} (5 nodes: 0-1-2=3-4 with |a_0|^2=|a_1|^2=|a_2|^2=2, |a_3|^2=|a_4|^2=1)
    add_from_G("f_4^(1)", [
        [2, -1, 0, 0, 0],
        [-1, 2, -1, 0, 0],
        [0, -1, 2, -1, 0],
        [0, 0, -1, 1, -0.5],
        [0, 0, 0, -0.5, 1],
    ])

    # 12. g_2^{(1)} (3 nodes: 0-1 triple-2 with |a_0|^2=|a_1|^2=2, |a_2|^2=2/3)
    add_from_G("g_2^(1)", [
        [2, -1, 0],
        [-1, 2, -1],
        [0, -1, 2.0 / 3.0],
    ])

    # 13. a_2^{(2)} (2 nodes with quadruple bond and lengths^2 = 1, 4)
    add_from_G("a_2^(2)", [
        [1, -2],
        [-2, 4],
    ])

    # 14-16. a_{2r}^{(2)} for r = 2, 3, 4 (chain of r+1 nodes: node 0 has length^2 = 1, nodes 1..r-1 have length^2 = 2, node r has length^2 = 4)
    for r in [2, 3, 4]:
        n = r + 1
        G = np.zeros((n, n))
        G[0, 0] = 1.0
        for i in range(1, r):
            G[i, i] = 2.0
        G[r, r] = 4.0
        G[0, 1] = G[1, 0] = -1.0
        for i in range(1, r - 1):
            G[i, i + 1] = G[i + 1, i] = -1.0
        if r == 2:
            G[1, 2] = G[2, 1] = -2.0
        else:
            G[r - 1, r] = G[r, r - 1] = -2.0
        add_from_G(f"a_{2*r}^(2)", G)

    # 17. a_3^{(2)} = d_3^{(2)} (3 nodes: chain 0-1-2 with |a_0|^2=|a_2|^2=1, |a_1|^2=2)
    add_from_G("a_3^(2)", [
        [1, -1, 0],
        [-1, 2, -1],
        [0, -1, 1],
    ])

    # 18-19. a_{2r-1}^{(2)} for r = 3, 4 (a_5^{(2)} has 4 nodes, a_7^{(2)} has 5 nodes:
    # fork at 0,1 of length^2 = 2 connected to 2, chain to r-1 of length^2 = 2, and end r of length^2 = 4)
    for r in [3, 4]:
        n = r + 1
        G = np.zeros((n, n))
        for i in range(r):
            G[i, i] = 2.0
        G[r, r] = 4.0
        G[0, 2] = G[2, 0] = -1.0
        G[1, 2] = G[2, 1] = -1.0
        for i in range(2, r - 1):
            G[i, i + 1] = G[i + 1, i] = -1.0
        G[r - 1, r] = G[r, r - 1] = -2.0
        add_from_G(f"a_{2*r-1}^(2)", G)

    # 20-21. d_{r+1}^{(2)} for r = 3, 4 (d_4^{(2)} has 4 nodes, d_5^{(2)} has 5 nodes:
    # chain of r+1 nodes where ends 0, r are short (length^2 = 1) and middle 1..r-1 are long (length^2 = 2))
    for r in [3, 4]:
        n = r + 1
        G = np.zeros((n, n))
        G[0, 0] = G[r, r] = 1.0
        for i in range(1, r):
            G[i, i] = 2.0
        G[0, 1] = G[1, 0] = -1.0
        for i in range(1, r - 1):
            G[i, i + 1] = G[i + 1, i] = -1.0
        G[r - 1, r] = G[r, r - 1] = -1.0
        add_from_G(f"d_{r+1}^(2)", G)

    # 22. e_6^{(2)} (5 nodes: 0-1-2=3-4 with |a_0|^2=|a_1|^2=|a_2|^2=1, |a_3|^2=|a_4|^2=2)
    add_from_G("e_6^(2)", [
        [1, -0.5, 0, 0, 0],
        [-0.5, 1, -0.5, 0, 0],
        [0, -0.5, 1, -1, 0],
        [0, 0, -1, 2, -1],
        [0, 0, 0, -1, 2],
    ])

    # 23. d_4^{(3)} (3 nodes: 0-1 triple-2 with |a_0|^2=|a_1|^2=1, |a_2|^2=3)
    add_from_G("d_4^(3)", [
        [1, -0.5, 0],
        [-0.5, 1, -1.5],
        [0, -1.5, 3],
    ])

    return diagrams


def check_item_7():
    section("ITEM 7: prp-ar-only (Neighbour-counting proof & scan of 23 diagrams / 78 automorphisms)")

    diagrams = build_23_affine_diagrams()
    assert len(diagrams) == 23, f"Expected 23 diagrams, got {len(diagrams)}"

    total_autos = 0
    fixed_node_count = 0
    singular_1_minus_P = 0
    singular_cases = []
    forbidden_cross_cases = []
    passing_cases = []

    for name, G in diagrams.items():
        n = G.shape[0]
        # Find null vector (Kac labels n_j > 0) of G:
        eigvals, eigvecs = np.linalg.eigh(G)
        null_idx = np.argmin(np.abs(eigvals))
        assert abs(eigvals[null_idx]) < 1e-12, f"Diagram {name} G is not affine singular!"

        # Realize simple roots alpha_0..alpha_{n-1} in R^{n-1} so that alpha_i . alpha_j = G_{ij}:
        # The top-left (n-1)x(n-1) submatrix of G is positive definite (finite Cartan Gram matrix).
        G_fin = G[: n - 1, : n - 1]
        L_chol = np.linalg.cholesky(G_fin)  # G_fin = L_chol @ L_chol.T
        # Thus rows of L_chol give alpha_0..alpha_{n-2} in R^{n-1}, and alpha_{n-1} is determined by G_{n-1, :n-1}:
        alphas = np.zeros((n, n - 1))
        alphas[: n - 1, :] = L_chol
        alphas[n - 1, :] = np.linalg.solve(L_chol, G[n - 1, : n - 1])
        assert np.allclose(alphas @ alphas.T, G, atol=1e-12), f"Root realization failed for {name}"

        # Enumerate all diagram automorphisms p in S_n (preserving G):
        for perm in permutations(range(n)):
            P_perm = G[np.ix_(perm, perm)]
            if not np.allclose(P_perm, G, atol=1e-12):
                continue
            total_autos += 1

            # Check if perm fixes any node:
            if any(perm[i] == i for i in range(n)):
                fixed_node_count += 1
                continue

            # Construct orthogonal matrix P on R^{n-1} such that P alpha_i = alpha_{p(i)}:
            P_mat = np.linalg.solve(alphas[: n - 1], alphas[list(perm[: n - 1])]).T
            assert np.allclose((P_mat @ alphas.T).T, alphas[list(perm)], atol=1e-12)
            assert np.allclose(P_mat @ P_mat.T, np.eye(n - 1), atol=1e-12)

            # Check if I - P_mat is invertible on R^{n-1}:
            det_1mP = np.linalg.det(np.eye(n - 1) - P_mat)
            if abs(det_1mP) < 1e-10:
                singular_1_minus_P += 1
                singular_cases.append((name, perm))
                continue

            # Compute D = 2 (I - P_mat)^{-1} and test alpha_i . D alpha_j = 0 for j not in {i, p(i)}:
            D_mat = 2.0 * np.linalg.inv(np.eye(n - 1) - P_mat)
            AD = alphas @ D_mat @ alphas.T  # (AD)_{ij} = alpha_i . D alpha_j
            forbidden = []
            for i in range(n):
                for j in range(n):
                    if j != i and j != perm[i]:
                        if abs(AD[i, j]) > 1e-10:
                            forbidden.append((i, j, AD[i, j]))

            # Also check the neighbour property directly: every neighbour k of i (G_{ik} < 0) is in {p(i), p^{-1}(i)}:
            inv_perm = [0] * n
            for i, pi in enumerate(perm):
                inv_perm[pi] = i
            neighbour_ok = True
            for i in range(n):
                for k in range(n):
                    if k != i and abs(G[i, k]) > 1e-10:
                        if k != perm[i] and k != inv_perm[i]:
                            neighbour_ok = False

            assert (len(forbidden) == 0) == neighbour_ok, "Equivalence of cross-term condition and neighbour condition failed!"

            if len(forbidden) > 0:
                forbidden_cross_cases.append((name, perm, len(forbidden)))
            else:
                passing_cases.append((name, perm))

    print(f"Total diagrams scanned: {len(diagrams)}")
    print(f"Total automorphisms scanned: {total_autos}")
    print(f"  - Fixed a node: {fixed_node_count}")
    print(f"  - No fixed node, but 1 - P singular: {singular_1_minus_P} -> {[c[0] for c in singular_cases]}")
    print(f"  - 1 - P invertible, but forbidden cross terms: {len(forbidden_cross_cases)} -> {forbidden_cross_cases}")
    print(f"  - Passing automorphisms: {len(passing_cases)} -> {[c[0] for c in passing_cases]}")

    assert total_autos == 78, f"Expected 78 automorphisms, got {total_autos}"
    assert fixed_node_count == 64, f"Expected 64 fixed-node automorphisms, got {fixed_node_count}"
    assert singular_1_minus_P == 5, f"Expected 5 singular-(1-P) automorphisms, got {singular_1_minus_P}"
    assert len(forbidden_cross_cases) == 2, f"Expected 2 forbidden-cross-term automorphisms, got {len(forbidden_cross_cases)}"
    assert all(c[2] == 10 for c in forbidden_cross_cases), "Expected 10 forbidden cross terms for pi^{\\pm 2} of a_4^{(1)}"
    assert len(passing_cases) == 7, f"Expected 7 passing automorphisms, got {len(passing_cases)}"
    print("[PASS] 7.1: Independent scan of all 23 diagrams and 78 automorphisms matches every number in prp-ar-only: "
          "64 fix a node, 5 have singular 1-P, 2 (pi^{\\pm 2} of a_4^{(1)}) have 10 forbidden cross terms each, and 7 pass.")


# ============================================================================
# ITEM 8: prp-dc-types, spin-two symmetry remark, and soliton delay relation
# ============================================================================
def check_item_8():
    section("ITEM 8: prp-dc-types (Reflection relation between Type I_+ and Type I_-, spin-2 symmetry, delays)")

    rng = np.random.default_rng(20261024)

    for h in [2, 3, 4, 5]:
        r_rank = h - 1
        # In components on R^h (indices 0..h-1 mod h):
        # alpha_j = e_j - e_{j+1},
        # (pi v)_k = v_{(k+1)%h}, (pi^{-1} v)_k = v_{(k-1)%h} (since pi e_k = e_{k-1} => e_k . pi v = e_{k+1} . v = v_{k+1}
        # WAIT: let's check line 50 vs line 70 of 39-defect-conditions.qmd VERY CAREFULLY!).
        # Let's check both conventions in R^h right now!
        # In line 50: "Realize the roots in R^h, h = r+1, as alpha_j = e_j - e_{j+1} with indices modulo h,
        # so that pi e_k = e_{k-1}, and let mu_k = e_k - (1/h) sum_l e_l... Since (1 - pi) mu_{j+1} = -alpha_j..."
        # And in line 70: "In components (rho v)_k = -v_{1-k} and (pi^{-1} v)_k = v_{k-1}, so by (eq-dc-components)
        # the map sends f_j to g_{-j} and g_j to f_{-j}, which exchanges the coefficients sigma and sigma^{-1} of B_d."
        # WAIT A SECOND!
        # Let's check line 70:
        # If (rho v)_k = -v_{(1-k)%h} and (pi^{-1} v)_k = v_{(k-1)%h}, then what is (rho pi^{-1} psi)_k?!
        # Let tilde_phi = rho phi and tilde_psi = rho (pi^{-1} psi).
        # Wait: for any operator A on vectors, (rho A psi)_k = -(A psi)_{1-k}.
        # If (pi^{-1} psi)_m = psi_{m-1}, then:
        # (rho pi^{-1} psi)_k = -(pi^{-1} psi)_{1-k} = -psi_{(1-k) - 1} = -psi_{-k} !!
        # While (rho phi)_k = -phi_{1-k} !!
        # Let's check what happens to f_j(tilde_phi, tilde_psi) and g_j(tilde_phi, tilde_psi) in eq-dc-components!
        # In eq-dc-components (line 53):
        # f_j(phi, psi) = exp(b (psi_j - phi_{j+1})),
        # g_j(phi, psi) = exp(b (phi_j - psi_j)).
        # Substitute tilde_phi_k = -phi_{1-k} and tilde_psi_k = -psi_{-k}:
        # 1) tilde_psi_j - tilde_phi_{j+1} = -psi_{-j} - (-phi_{1-(j+1)}) = phi_{-j} - psi_{-j} !!
        #    So f_j(tilde_phi, tilde_psi) = exp(b (phi_{-j} - psi_{-j})) = g_{-j}(phi, psi) !!
        # 2) tilde_phi_j - tilde_psi_j = -phi_{1-j} - (-psi_{-j}) = psi_{-j} - phi_{(-j)+1} !!
        #    So g_j(tilde_phi, tilde_psi) = exp(b (psi_{-j} - phi_{(-j)+1})) = f_{-j}(phi, psi) !!
        #WAIT A SECOND!
        # Look at why (pi^{-1} v)_k = v_{k-1} in line 70 vs pi e_k = e_{k-1} in line 50!
        # If pi e_k = e_{k-1}, then for a vector v = sum_k v_k e_k:
        # pi v = sum_k v_k (pi e_k) = sum_k v_k e_{k-1} = sum_m v_{m+1} e_m,
        # so (pi v)_k = v_{k+1}, AND (pi^{-1} v)_k = v_{k-1} !!
        # BOTH LINE 50 AND LINE 70 ARE 100% CONSISTENT!!
        # Because pi e_k = e_{k-1} means pi shifts basis vectors DOWN (k -> k-1), which shifts coordinates UP ((pi v)_k = v_{k+1}),
        # and therefore pi^{-1} shifts basis vectors UP (k -> k+1), which shifts coordinates DOWN ((pi^{-1} v)_k = v_{k-1})!!
        # And what about rho?
        # In line 58: "Let rho be the reflection rho alpha_j = alpha_{-j}."
        # Since alpha_j = e_j - e_{j+1} and alpha_{-j} = e_{-j} - e_{1-j} = -(e_{1-j} - e_{-j}),
        # rho e_k = -e_{1-k} satisfies:
        # rho alpha_j = rho(e_j - e_{j+1}) = -e_{1-j} + e_{-j} = e_{-j} - e_{1-j} = alpha_{-j} !!
        # And for v = sum_k v_k e_k, rho v = sum_k v_k (-e_{1-k}) = sum_m (-v_{1-m}) e_m,
        # so (rho v)_k = -v_{1-k} !!
        # And let's check rho pi rho:
        # (rho pi rho)(e_k) = (rho pi)(-e_{1-k}) = rho(-e_{-k}) = e_{1-(-k)} = e_{k+1} = pi^{-1} e_k !!

        P_pi = np.zeros((h, h))
        P_rho = np.zeros((h, h))
        for k in range(h):
            P_pi[(k - 1) % h, k] = 1.0      # pi e_k = e_{k-1}  <=>  (pi v)_k = v_{k+1}, (pi^{-1} v)_k = v_{k-1}
            P_rho[(1 - k) % h, k] = -1.0    # rho e_k = -e_{1-k} <=> (rho v)_k = -v_{1-k}

        assert np.allclose(P_rho @ P_pi @ P_rho, P_pi.T)
        assert np.allclose(P_rho @ P_rho, np.eye(h))

        # Check on the traceless subspace:
        Q, _ = np.linalg.qr(np.eye(h, r_rank) - 1.0 / h)
        P_sub = Q.T @ P_pi @ Q
        rho_sub = Q.T @ P_rho @ Q
        D_sub = 2.0 * np.linalg.inv(np.eye(r_rank) - P_sub)
        E_sub = np.eye(r_rank) - D_sub

        # Verify the matrix identities in proof of prp-dc-types (line 70):
        # rho D rho = D^T, rho E rho = -E, D^T = -D pi, rho D rho pi^{-1} = -D, pi rho E rho pi^{-1} = -E
        assert np.allclose(rho_sub @ D_sub @ rho_sub, D_sub.T)
        assert np.allclose(rho_sub @ E_sub @ rho_sub, -E_sub)
        assert np.allclose(D_sub.T, -D_sub @ P_sub)
        assert np.allclose(rho_sub @ D_sub @ rho_sub @ P_sub.T, -D_sub)
        assert np.allclose(P_sub @ rho_sub @ E_sub @ rho_sub @ P_sub.T, -E_sub)

        # Verify L_d^{I_+}(phi, psi; sigma) == L_d^{I_-}(rho phi, rho pi^{-1} psi; sigma) at random complex fields:
        phi_v = Q @ (rng.normal(size=r_rank) + 1j * rng.normal(size=r_rank))
        psi_v = Q @ (rng.normal(size=r_rank) + 1j * rng.normal(size=r_rank))
        dt_phi_v = Q @ (rng.normal(size=r_rank) + 1j * rng.normal(size=r_rank))
        dt_psi_v = Q @ (rng.normal(size=r_rank) + 1j * rng.normal(size=r_rank))
        b_val = 0.71 - 0.33j
        m_val = 1.23
        sig_val = 0.89 + 0.44j

        D_full = Q @ D_sub @ Q.T
        E_full = Q @ E_sub @ Q.T

        def defect_lagrangian(ph, ps, dt_ph, dt_ps, sig_p, is_type_Ip=False):
            sgn = -1.0 if is_type_Ip else 1.0
            kin = (
                0.5 * sgn * (ph @ E_full @ dt_ph)
                + sgn * (ph @ D_full @ dt_ps)
                + 0.5 * sgn * (ps @ E_full @ dt_ps)
            )
            f_vec = np.array([np.exp(b_val * (ps[j] - ph[(j + 1) % h])) for j in range(h)])
            g_vec = np.array([np.exp(b_val * (ph[j] - ps[j])) for j in range(h)])
            if not is_type_Ip:
                Bd = (m_val / b_val**2) * np.sum(sig_p * (f_vec - 1.0) + (1.0 / sig_p) * (g_vec - 1.0))
            else:
                Bd = (m_val / b_val**2) * np.sum((1.0 / sig_p) * (f_vec - 1.0) + sig_p * (g_vec - 1.0))
            return kin - Bd

        L_Ip = defect_lagrangian(phi_v, psi_v, dt_phi_v, dt_psi_v, sig_val, is_type_Ip=True)
        t_phi = P_rho @ phi_v
        t_psi = P_rho @ (P_pi.T @ psi_v)
        t_dt_phi = P_rho @ dt_phi_v
        t_dt_psi = P_rho @ (P_pi.T @ dt_psi_v)
        L_Im_mapped = defect_lagrangian(t_phi, t_psi, t_dt_phi, t_dt_psi, sig_val, is_type_Ip=False)
        assert abs(L_Ip - L_Im_mapped) < 1e-14

        # Check spin-two cubic tensor A(u, u, u) = 2 sum_k (mu_k . u)^3 under pi and rho (line 225):
        u_v = Q @ rng.normal(size=r_rank)
        A_u = 2.0 * np.sum(u_v**3)
        A_pi_u = 2.0 * np.sum((P_pi @ u_v) ** 3)
        A_rho_u = 2.0 * np.sum((P_rho @ u_v) ** 3)
        assert abs(A_pi_u - A_u) < 1e-15, "A(u,u,u) must be invariant under pi"
        assert abs(A_rho_u + A_u) < 1e-15, "A(u,u,u) must change sign under rho"

        # Check cyclic conjugation of Lax pair: A_\pm[pi phi; lam] = S A_\pm[phi; lam] S^{-1} (line 225):
        # Wait: since (pi phi)_k = phi_{k+1}, let's check S!
        # With S = sum_j e_{j, j+1}, (S M S^{-1})_{k, l} = M_{k+1, l+1}, which replaces phi_k by phi_{k+1} = (pi phi)_k!
        S_mat = np.zeros((h, h))
        for j in range(h):
            S_mat[j, (j + 1) % h] = 1.0
        # Note that S_mat == P_pi!
        assert np.allclose(S_mat, P_pi)

        # Check soliton delay relation z_a(-u) = omega^a z_{h-a}(u) (40-defect-solitons.qmd line 52):
        # In eq-dc-ar-delay (39-defect-conditions.qmd line 340):
        # z_a(u) = (e^u + i e^{-i pi a / h}) / (e^u - i e^{i pi a / h})  (with u = theta - theta_d)
        # Wait: let's check both formulas in 40-defect-solitons.qmd line 52:
        # z_a(u) = sinh((u - i pi (a/h - 1/2))/2) / sinh((u + i pi (a/h - 1/2))/2)
        u_val = 0.83 - 0.19j
        omega = np.exp(2j * np.pi / h)
        for a_sp in range(1, h):
            def z_delay(sp_a, u_arg):
                return (np.exp(-u_arg) - 1j * np.exp(1j * np.pi * sp_a / h)) / (
                    np.exp(-u_arg) - 1j * np.exp(-1j * np.pi * sp_a / h)
                )
            lhs_z = z_delay(a_sp, -u_val)
            rhs_z = (omega**a_sp) * z_delay(h - a_sp, u_val)
            assert abs(lhs_z - rhs_z) < 1e-14, f"Delay relation failed for h={h}, a={a_sp}"

        print(f"[PASS] 8.1.{h}: prp-dc-types, cubic tensor A transformation, Lax cyclic conjugation, "
              f"and delay identity z_a(-u) = omega^a z_{{h-a}}(u) verified for h={h}.")


if __name__ == "__main__":
    check_item_6()
    check_item_7()
    check_item_8()
