# Exchange statistics and Bloch sectors

Checks behind Propositions `prp-exchange-twist` and `prp-bloch-sectors` and the Bloch-sector rows of the table in Section `sec-sector-status` (`part7-foundations/32-scattering.qmd`). The scripts import `lib.py`, `a22three.py`, `a22.py`, `rsolve.py` and `dtw.py` from the parent directory `sectors/`. Run them from this directory; outputs are in `out/`.

Conventions: $q^2=-e^{i\psi}$, $\psi=2\pi^2/3\xi$ reduced to $(-\pi,\pi]$, $q_{\rm code}=\bar q_{\rm TW}$, as in `../a22three.py`. Couplings are offset by $0.0013\pi$ in $\xi$ to avoid roots of unity. Spinor amplitudes use $q=e^{-i\pi\omega}$, as in `../spinor_scan.py`. The twist in the Bloch sector $\alpha$ is $e^{i\alpha\cdot\mu_1}$ on soliton 1; its sign is immaterial. "Unbroken" means $\max\lvert\lvert s\rvert-1\rvert<10^{-6}$ over the rapidities sampled.

| Script (arguments) | Output | What it checks |
|---|---|---|
| `ik_common.py`, `spin_common.py` | | shared setup: Izergin–Korepin and spinor R-matrices, graded and twisted transfer matrices |
| `t1_graded_is_twist.py` | `out_t1_graded_is_twist.txt` | $T_1^{\rm gr}=(-1)^{g_1(G-1)}T_1$ for $a_2^{(2)}$ kinks, $N\le4$ |
| `t2_twist_scan.py 2 40` | `out_t2_N2.txt` | two kinks: $\max\lvert\lvert s\rvert-1\rvert$ and Krein pairing by charge $Q$ and twist $\alpha$ |
| `t3_threshold_N2.py` | `out_t3_threshold_N2.txt` | two kinks, $Q=0$: critical twist $\alpha_c(\psi)$ by bisection; compare $\min(2\lvert\psi\rvert-\pi,\pi-\lvert\psi\rvert)$ |
| `t4_twist_N.py N ntr alphas xis` | `out_t4_N3.txt`, `out_t4_N4.txt`, `out_t4_N5.txt` | $N=3,4,5$ kinks by $Q$ and $\alpha$; run as `t4_twist_N.py 3 80 0,0.125,0.25,0.333,0.5,0.667,0.75,1 0.25,0.3,0.45,0.6,0.9,1.3,1.6,2.0,2.5,3.0`, `4 40 0,0.5,1 0.45,0.6,0.9,1.3,1.6,2.5,3.0,4.0`, `5 25 0,1 3.7,3.03,2.22` |
| `t5_spin_graded_twist.py` | `out_t5_spin_graded_twist.txt` | the identity of `t1` for the $c_2^{(1)}$, $c_3^{(1)}$, $a_5^{(2)}$ spinors; breaking at $\alpha=0$ by total weight, two solitons |
| `t6_spin_torus.py n c steps` | `out_t6_c2.txt` (`2 c 24`), `out_t6_c3.txt` (`3 c 12`) | spinors, two solitons, $Q=0$ and $Q=e_1$: scan of the Bloch torus |
| `t7_krein_ops.py` | `out_t7_krein_ops.txt` | $R_{12}(x)^\dagger=R_{12}(1/x)$ and $\mathcal CR_{12}\mathcal C=R_{21}$ for the Izergin–Korepin and spinor R-matrices |
| `t8_transpose.py` | `out_t8_transpose.txt` | $R_{12}^{\mathsf T}=R_{21}$ (holds for the spinors, fails for Izergin–Korepin in this gradation) |
| `t9_ik_gauge.py` | `out_t9_ik_gauge.txt` | no diagonal regauging $x^{am}$ restores $R^{\mathsf T}=R_{21}$ for Izergin–Korepin |
| `t10_crossC.py omega` | `out_t10_crossC_w2.37.txt`, `out_t10_crossC_w0.002.txt` | $c_n^{(1)}$ spinor charge conjugation: $C^{-1}C^{\mathsf T}$ is diagonal; as $q\to1$, $C^{\mathsf T}=\epsilon_nC$ with $\epsilon_n=-1,-1,+1,+1$ |
| `t11_psi_pi.py` | `out_t11_psi_pi.txt` | at $\psi=\pi$ ($\xi=2\pi/3$, $q^2=1$) the kink amplitude is diagonal and unbroken |
| `t12_excited_sectors.py` | `out_t12_excited_sectors.txt` | soliton–excited soliton ($\check R(-z)$) by $Q$ at $\alpha=0,\pi$ |
| `t13_Qpm1_all_alpha.py` | `out_t13_Qpm1_all_alpha.txt` | two kinks, $Q=\pm1$: unimodular for every $\alpha\in[-\pi,\pi]$ |
