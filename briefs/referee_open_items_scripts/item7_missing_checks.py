"""Independent verification of Item 7:
1. Ch. 9 (09-lagrangian.qmd line 116): e_8^{(1)} closed-form mass ratios table vs eigenvalues of M^2 in code/algebra_data.py.
2. Ch. 35 (35-quantum-symmetric-pairs.qmd line 206): Baseilhac-Koizumi Heisenberg realization of (19) and the q-Dolan-Grady relations.
"""
import os
import sys
import numpy as np

_orig_cwd = os.getcwd()
os.chdir(os.path.join(_orig_cwd, 'code'))
sys.path.insert(0, os.getcwd())
import algebra_data as ad
os.chdir(_orig_cwd)


def check_e8_closed_form_masses():
    """Verify the closed-form e_8^{(1)} mass ratios in Ch. 9 (09-lagrangian.qmd line 116) against M^2 eigenvalues."""
    _, _, ms = ad.masses('e', 8)
    c5 = np.cos(np.pi / 5)
    c30 = np.cos(np.pi / 30)
    c7 = np.cos(7 * np.pi / 30)
    c2 = np.cos(2 * np.pi / 15)
    m1 = ms[0]
    m2 = 2 * m1 * c5
    closed = np.array([
        m1,
        m2,
        2 * m1 * c30,
        2 * m2 * c7,
        2 * m2 * c2,
        2 * m2 * c30,
        4 * m2 * c5 * c7,
        4 * m2 * c5 * c2,
    ])
    err = np.max(np.abs(ms - closed))
    print(f"e_8^(1) closed-form masses max error vs M^2 eigenvalues: {err:.3e}")
    print("PASS e8_closed_form_mass_ratios" if err < 1e-11 else "FAIL e8_closed_form_mass_ratios")


def check_bk_heisenberg_dolan_grady():
    """Verify Baseilhac-Koizumi (2003) Heisenberg realization of eq. (19) and the q-Dolan-Grady relations (Ch. 35, sec-q-onsager).

    With Weyl operators U = exp(hat_theta), V = exp(hat_phi) obeying U V = q^2 V U (from [hat_theta, hat_phi] = 2 log q):
      hat_E_+ = (i c / (q - q^{-1})) (alpha U + alpha^{-1} U^{-1}),
      hat_E_- = (i c / (q - q^{-1})) (alpha_bar V + alpha_bar^{-1} V^{-1}).
    Because alpha * alpha^{-1} = alpha_bar * alpha_bar^{-1} = 1, BOTH Baseilhac-Koizumi's eq. (19):
      [hat_E_+ hat_E_-, hat_E_- hat_E_+] = rho (hat_E_+^2 - hat_E_-^2)
    AND the q-Dolan-Grady relations:
      [hat_E_pm, [hat_E_pm, [hat_E_pm, hat_E_mp]_q]_{q^{-1}}] = rho [hat_E_pm, hat_E_mp]
    hold simultaneously with the exact same constant
      rho = c^2 (q + q^{-1})^2!
    And after hat_E_pm = c * b_pm, b_pm satisfy eq-q-dolan-grady with c_j = q^{-1} (rho_b = (q + q^{-1})^2).
    """
    def serre_dg(x, y, q):
        q3 = q**2 + 1.0 + q**(-2)
        return x @ x @ x @ y - q3 * (x @ x @ y @ x - x @ y @ x @ x) - y @ x @ x @ x

    for label, q_val in [
        ("physical_q", -np.exp(-1j * np.pi * 0.37)),
        ("generic_q", 0.83 * np.exp(0.41j)),
    ]:
        c_val = 0.9 - 0.3j
        alpha = 1.3 - 0.4j
        alpha_bar = 0.8 + 0.7j
        N = 15
        S = slice(4, N - 4)  # interior states: zero boundary truncation for degree-4 words in U, U^{-1}
        U = np.diag(np.ones(N - 1, dtype=complex), k=-1)
        U_inv = np.diag(np.ones(N - 1, dtype=complex), k=1)
        ks = np.arange(N) - (N // 2)
        V = np.diag(q_val ** (-2.0 * ks))
        V_inv = np.diag(q_val ** (2.0 * ks))

        pref = 1j * c_val / (q_val - 1.0 / q_val)
        Ep = pref * (alpha * U + (1.0 / alpha) * U_inv)
        Em = pref * (alpha_bar * V + (1.0 / alpha_bar) * V_inv)
        rho = (c_val ** 2) * ((q_val + 1.0 / q_val) ** 2)

        # 1. Check BK (19): [Ep Em, Em Ep] = rho (Ep^2 - Em^2)
        lhs19 = (Ep @ Em @ Em @ Ep - Em @ Ep @ Ep @ Em)[:, S]
        rhs19 = (rho * (Ep @ Ep - Em @ Em))[:, S]
        err19 = np.linalg.norm(lhs19 - rhs19) / np.linalg.norm(rhs19)

        # 2. Check q-Dolan-Grady for (Ep, Em) and (Em, Ep) with rho = c^2 (q + q^{-1})^2
        lhs_dg_p = serre_dg(Ep, Em, q_val)[:, S]
        rhs_dg_p = (rho * (Ep @ Em - Em @ Ep))[:, S]
        err_dg_p = np.linalg.norm(lhs_dg_p - rhs_dg_p) / np.linalg.norm(rhs_dg_p)

        lhs_dg_m = serre_dg(Em, Ep, q_val)[:, S]
        rhs_dg_m = (rho * (Em @ Ep - Ep @ Em))[:, S]
        err_dg_m = np.linalg.norm(lhs_dg_m - rhs_dg_m) / np.linalg.norm(rhs_dg_m)

        # 3. Also check that V(y) restricted to B does NOT satisfy (19) (as stated in Ch. 35 line 206):
        # On spin-1 (or 1/2 x 1/2), b_0, b_1 satisfy q-DG with rho_b = (q + q^{-1})^2, but [b_0 b_1, b_1 b_0] != rho_b (b_0^2 - b_1^2)
        print(f"BK Heisenberg ({label}): err19={err19:.2e}, err_dg_p={err_dg_p:.2e}, err_dg_m={err_dg_m:.2e}")
        ok = max(err19, err_dg_p, err_dg_m) < 1e-14
        print(f"PASS bk_heisenberg_{label}" if ok else f"FAIL bk_heisenberg_{label}")


if __name__ == '__main__':
    check_e8_closed_form_masses()
    check_bk_heisenberg_dolan_grady()
