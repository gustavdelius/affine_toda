# Referee Report: Open Items Across PR Reviews (#5, #6, #7, #10, #16, #18, #22, #23, #27, #28, #29, #40, #46)

- **Date:** 2026-10-25
- **Model:** Gemini 4 Argon
- **Commit:** `9b912ee` (`collect_pending_pr_reviews`)
- **Scratch scripts:** `referee_scratch/open_items/` (`item1_node_labels.py`, `item2_singular_wall.py`, `item3_oneloop_ratio.py`, `item6_step3_nn.py`, `item7_missing_checks.py`)

---

## Summary of Verdicts

| Item | Topic | Verdict | Summary of Finding |
|---|---|---|---|
| **1** | Part VII node labelling ($\hat\epsilon_j$ vs $C_j$, PR #6 & PR #10) | **correct** (with a clarifying fix for `33-boundary-quantum-groups.qmd` lines 27–57 & 195) | Both `30-boundary-sine-gordon.qmd` (line 136) and `33-boundary-quantum-groups.qmd` (lines 57, 273) explicitly note the sine-Gordon crossing $\hat\epsilon_1 \leftrightarrow C_0, \hat\epsilon_0 \leftrightarrow C_1$. However, `33-boundary-quantum-groups.qmd` (`sec-boundary-charges`, lines 27–44) defines the general-$\hat g$ chiral current with the opposite exponent sign from `21-quantum-group-symmetry.qmd` (`sec-qg-coproduct`, line 58), causing an internal sign/pairing tension between `eq-boundary-anomaly`/`eq-boundary-topological` and `sec-an-nondiagonal-k` (line 195). |
| **2** | Part VII Ch. 32 singular-wall limits `0.719 / 0.937 / 0.147` (PR #22) | **correct** | All three limits (`|K_{1,5}| \to 0.7191932`, `|K_{1,5}| \to 0.9374372`, `|K_{1,4}| \to 0.1471365`) and diagonal entries (`0.15563 + 0.67716i`, `0.32535 + 0.93408i`) are independently confirmed to 6+ digits in `referee_scratch/open_items/item2_singular_wall.py`. Minor initial-guess typo in `code/solitonic_linear.py` line 376 (`0.2579` vs actual root `d_1^* = 0.257789025`). |
| **3** | Part VI `code/rmatrix/smatrix/oneloop_ratio.py` $a_5^{(2)}$ quirk (`0.693` vs `0.5`, PR #29) | **wrong** (bug in `oneloop_ratio.py` line 37 + false prose claim in Ch. 20 line 97) | `toda_data(a5, [1,2,1,1])` sorts `m2` ascending: index `0` is particle $n=3$ ($m_3 = 1$), whereas particle $1$ ($m_1 = 2\sin(\pi/5) \approx 1.17557$, whose elastic S-matrix is `dgz11`) is at index `1`. Line 37 hardcoded `tree_over_4m2sinh(..., 0, th)`, dividing the tree amplitude of particle 3 by the S-matrix slope of particle 1. Evaluating at particle 1 (`p1 = 1`) makes `tr/g(theta)` constant (`1/(4\pi) = 0.079577`) and gives the exact `delta H / B_tree = 0.500000`. Also, `20-one-loop-masses.qmd` line 97 falsely claims that all mass-ratio shifts vanish at one loop in $a_{2n-1}^{(2)}$. |
| **4a** | PR #7: `appendices/A-algebra-data.qmd` (`tbl-app-twisted` parents of $a_{2n}^{(2)}$) | **correct** (already resolved) | Line 25 of `appendices/A-algebra-data.qmd` lists `$a_{2n}^{(1)}$ or $d_{2n+2}^{(1)}$`, and line 49 explains both foldings. |
| **4b** | PR #7: `code/irf/checks.py` and `code/irf/tm2.py` breather threshold | **wrong** (5 stale output files in `code/irf/`) | `code/irf/checks.py` (line 47) and `code/irf/tm2.py` (line 15) use `3*r > 2*p`, matching Ch. 28 (`28-rsos.qmd` line 66), but 5 stored output files (`out_ha.txt`, `out_qsign/irf_ha_q1.txt`, `out_qsign/irf_ha_q-1.txt`, `out_qsign/irf_tm2_q1.txt`, `out_qsign/irf_tm2_q-1.txt`) were never regenerated and still show `breathers=yes` for `W3(5,7)` and `W3(8,11)`. |
| **4c** | PR #7: `code/thimbles/check_th4.py` | **correct** (already resolved) | `check_th4.py` auto-generates `cps_a2L1_k1.00.npy`, does not import `check_th2.py`, and `out_check_th4.txt` has matching signs between residual and prediction. |
| **5** | Appendix E (`appendices/E-original-results.qmd`) audit (PR #18, #22, #46) | **wrong** (label mismatch on line 45, untagged section on line 16, 14 missing statements in Parts VII–VIII) | Line 45 of `E-original-results.qmd` has `**[Theorem for $B_1$; ...]**` citing `@sec-sg-breather-transmission`, whereas `@prp-sg-breather-transmission` is `**[Established for $B_1$, with a new proof; ...]**`. Line 16 cites `@sec-one-loop-masses-folded` with `**[Established]**` (which has no statement div and is a literature correction). Fourteen original propositions in Parts VII–VIII should be added. |
| **6** | Part IX Ch. 48 (`48-continuum-window.qmd`) Step 3 of `thm-nn-distances` (PR #5 & PR #16) | **correct** | Steps 3.1–3.5 are mathematically complete and independently verified from scratch in `referee_scratch/open_items/item6_step3_nn.py`. Replacement text provided for line 109 to cite `code/review614/step3_combinatorics.py`. |
| **7** | Missing script checks in Ch. 9 (`09-lagrangian.qmd` line 116) & Ch. 35 (`35-quantum-symmetric-pairs.qmd` line 206) (PR #27 & PR #28) | **correct** (formulas verified; check code provided) | Both the $e_8^{(1)}$ closed-form mass table in Ch. 9 and the Baseilhac–Koizumi Heisenberg realization of (19) and $q$-Dolan–Grady in Ch. 35 are verified to $<10^{-12}$ in `referee_scratch/open_items/item7_missing_checks.py`. Drop-in check snippets for `code/algebra_data.py` and `code/qsp_checks.py` are provided below. |
| **8** | Preface (`preface.qmd`) wording (PR #23) | **correct** (optional tightening provided) | Both sentences in `preface.qmd` (lines 15 and 35) are consistent with the chapter labels, though optional tightened wording is suggested below to make the `[Conjecture]`/`[Numerical]` status of the Part VI/VII/VIII spin-1/vector checks explicit. |

---

## Item 1: Part VII Node Labelling ($\hat\epsilon_j$ vs $C_j$, PR #6 & PR #10)

- **Verdict:** **correct** in `30-boundary-sine-gordon.qmd` (`sec-bsg-parameters`, line 136) and `33-boundary-quantum-groups.qmd` (`sec-sg-k-from-symmetry`, line 57; `sec-boundary-qg-notes`, line 273), **with a clarifying fix needed** in `33-boundary-quantum-groups.qmd` (`sec-boundary-charges`, lines 27–45, and `sec-an-nondiagonal-k`, line 195).

### 1.1 Independent Derivation & Verification

In `referee_scratch/open_items/item1_node_labels.py`, we verified the parameter chain across Chapters 21, 30, 33, and 35:

1. **Lagrangian boundary parameters (`30-boundary-sine-gordon.qmd`, `eq-bsg-c-dictionary`, line 29):**
   With $\alpha_1 = \sqrt{2}, \alpha_0 = -\sqrt{2}$, $\beta_{\rm SG} = \sqrt{2}\beta$, and $\mathcal B_{\rm imag} = \frac{2m}{\beta^2}\sum_{j=0,1} C_j (e^{i\beta\alpha_j\phi(0)/2} - 1)$:
   $$C_0 = \frac{\beta_{\rm SG}^2 M_0}{4m_0} e^{+i\beta_{\rm SG}\phi_0/2}, \qquad C_1 = \frac{\beta_{\rm SG}^2 M_0}{4m_0} e^{-i\beta_{\rm SG}\phi_0/2}.$$
2. **Reflection matrix and coideal parameters (`30-boundary-sine-gordon.qmd`, `eq-bsg-k-matrix` & `eq-bsg-eta-theta`, lines 59, 69; `33-boundary-quantum-groups.qmd`, `eq-sg-k-coideal`, line 51):**
   Intertwining $K(\theta)\pi_{V(\theta)}(b_j) = \pi_{V(-\theta)}(b_j)K(\theta)$ for $b_j = e_j + q^{-1} f_j k_j + \hat\epsilon_j k_j$ ($j=0,1$) on the spin-$1/2$ representation with $e_1 = y\sigma_+, f_1 = y^{-1}\sigma_-, k_1 = q^{\sigma_z}$ and $e_0 = y\sigma_-, f_0 = y^{-1}\sigma_+, k_0 = q^{-\sigma_z}$ fixes the diagonal entries of $K(\theta)$ to $(q - q^{-1})(\hat\epsilon_1 y + \hat\epsilon_0 y^{-1})$ and $(q - q^{-1})(\hat\epsilon_0 y + \hat\epsilon_1 y^{-1})$.
   Under the UV–IR relation (`cnj-bsg-uv-ir`, line 122 of `30-boundary-sine-gordon.qmd`) and `eq-bsg-eta-theta` (line 69):
   $$\hat\epsilon_1 = -\frac{\cos(\eta - i\vartheta)}{\sin\pi\lambda} = -\frac{\tilde\mu/\tilde\mu_{\rm crit}}{2\sin\pi\lambda} e^{+i\beta_{\rm SG}\phi_0/2} \propto C_0, \qquad \hat\epsilon_0 = -\frac{\cos(\eta + i\vartheta)}{\sin\pi\lambda} = -\frac{\tilde\mu/\tilde\mu_{\rm crit}}{2\sin\pi\lambda} e^{-i\beta_{\rm SG}\phi_0/2} \propto C_1.$$
3. **Why $\hat\epsilon_1$ pairs with $C_0$ in sine-Gordon vs why $\hat\epsilon_j$ pairs with $C_j$ in general $\hat g$ (`21-quantum-group-symmetry.qmd` vs `33-boundary-quantum-groups.qmd`):**
   - In `part6-soliton-smatrices/21-quantum-group-symmetry.qmd`:
     - For sine-Gordon (`sec-nonlocal-charges`, line 16, and `sec-qg-coproduct`, line 58), the chiral current is defined with a **plus** sign in the exponent, $J_+ = :e^{+ia\phi_L}:$, and $(Q_+, \bar Q_-, T)$ are assigned to node $1$ $(e_1, f_1 k_1^{-1}, h_1)$. Because $\langle\phi_L(z)\phi_L(w)\rangle = -\frac{1}{4\pi}\ln(z-w)$, $e^{+ia\phi_L}$ has its double pole with $e^{-i\beta_{\rm SG}\phi(0)/2} = e^{+i\beta\alpha_0\phi(0)/2}$ (whose coefficient in $\mathcal B_{\rm imag}$ is $C_0$, node $0$), while its fused product at the boundary has exponent $+i(a - \beta_{\rm SG})\phi_L(0) = +i\pi\lambda T$ ($k_1 = q^T$).
     - By contrast, for general $\hat g$ (`sec-qg-gradation`, line 64 of `21-quantum-group-symmetry.qmd`, and `sec-boundary-charges`, lines 27–33 of `33-boundary-quantum-groups.qmd`), the node-$j$ chiral current is written with a **minus** sign in the exponent, $\exp(-i\frac{4\pi}{\beta}\vec\alpha_j^\vee\cdot\vec\phi_L)$. That current has its double pole directly with $\exp(+i\beta\vec\alpha_j\cdot\vec\phi/2)$ (coefficient $C_j$, node $j$), so $\hat\epsilon_j \propto C_j$ without crossing!
   - **The internal tension in `33-boundary-quantum-groups.qmd`:**
     1. In `sec-boundary-charges` (lines 27–44), using $\exp(-i\frac{4\pi}{\beta}\vec\alpha_j^\vee\cdot\vec\phi_L)$ alongside $\exp(+i\beta\vec\alpha_j\cdot\vec\phi(0)/2)$ gives the fused boundary operator $\exp(-i(\frac{4\pi}{\beta}\vec\alpha_j^\vee - \frac{\beta}{2}\vec\alpha_j)\cdot\vec\phi(0)) = \exp(-i\pi\omega_j T_j)$ (with $T_j$ as defined in `eq-boundary-topological`, line 42), which has the **opposite** sign in the exponent from `eq-boundary-anomaly` (line 33, $e^{+i\pi\omega_j T_j}$) unless the chiral current is $\exp(+i\frac{4\pi}{\beta}\vec\alpha_j^\vee\cdot\vec\phi_L)$ paired with $\exp(-i\beta\vec\alpha_j\cdot\vec\phi(0)/2)$ (i.e. node $\bar j$, which for $a_1^{(1)}$ is $1 \leftrightarrow 0$)!
     2. In `sec-an-nondiagonal-k` (line 195 of `33-boundary-quantum-groups.qmd`), the text states that $\hat\epsilon_1 = \dots = \hat\epsilon_n = \epsilon_*, \hat\epsilon_0 = s\epsilon_*$ has classical counterpart $(-1,+1,\dots,+1)$ for $s=-1$, implicitly using the uncrossed pairing $(\hat\epsilon_0,\hat\epsilon_1,\dots,\hat\epsilon_n) \leftrightarrow (C_0,C_1,\dots,C_n)$. Because any cyclic permutation of $(-1,+1,\dots,+1)$ is in the same equivalence class (`sec-solitonic-classification`), this does not affect the classification, but the convention should be stated explicitly.

**Output of `referee_scratch/open_items/item1_node_labels.py`:**
```text
PASS eq_bsg_c_dictionary
PASS eq_bsg_k_matrix_and_gz
PASS sg_coideal_intertwining_and_crossing
PASS an_coideal_sign_patterns
```

### 1.2 Recommended Replacement Text

In `part7-half-line/33-boundary-quantum-groups.qmd`, replace lines 44–57 with:

```markdown
In `eq-boundary-anomaly` the chiral fields are equated at $x=0$, $\vec\phi_L(0)=\vec\phi_R(0)=\frac12\vec\phi(0)-\frac12\vec\phi(-\infty)$ (the Neumann value of the constant mode; for general boundaries the zero-mode shift goes into $\hat\epsilon_j$). Note the sign convention: in @sec-qg-coproduct the sine-Gordon current of $Q_+=e_1$ was defined with $+ia\phi_L$, so its double pole at the wall comes from the $e^{-i\beta_{\rm SG}\phi(0)/2}=e^{i\beta\alpha_0\phi(0)/2}$ term of $\mathcal B_{\rm imag}$ (coefficient $C_0$) while its fused exponent is $+i\pi\lambda T$ ($k_1=q^T$). With the general-$\hat g$ sign $-i\frac{4\pi}\beta\vec\alpha_j^\vee\cdot\vec\phi_L$ of `eq-boundary-anomaly`, the double pole comes directly from the $C_j\,e^{i\beta\vec\alpha_j\cdot\vec\phi(0)/2}$ term of node $j$, and the fused exponent is $+i\pi\omega_jT_j$ after $T_j\to-T_j$ (or after pairing $e_j$ with $-\vec\alpha_j$). Up to the proportionality constants, which the tree-level calculation does not fix consistently across loop orders, (-@eq-boundary-anomaly) is $\partial_t(Q_j+\bar Q_j)\propto\hat\epsilon_j\,\partial_t(q_j^{T_j})$. In the quantum theory $q_j^{T_j}=k_j$, on the bulk charges the antipode is $S(e_j)=-k_j^{-1}e_j$, so $\bar Q_j\sim q_j^{-1}f_jk_j$, and the conserved combination is the **coideal generator**

$$
b_j=e_j+q_j^{-1}f_jk_j+\hat\epsilon_j\,k_j,\qquad j=0,1,\dots,r,
$$ {#eq-coideal-generators}

with a scalar boundary parameter $\hat\epsilon_j$ for each node of the affine Dynkin diagram [@delius_mackay2003, eq. (2.6)].

## Sine-Gordon K-matrix from the coideal {#sec-sg-k-from-symmetry}

For sine-Gordon theory the soliton multiplet $V(\theta)$ is the two-dimensional evaluation representation of $U_q(\widehat{sl}_2)$ at $y=e^{\lambda\theta}$ (@sec-six-vertex), with $e_1=y\sigma_+$, $f_1=y^{-1}\sigma_-$, $k_1=q^{\sigma_z}$ and $e_0=y\sigma_-$, $f_0=y^{-1}\sigma_+$, $k_0=q^{-\sigma_z}$, at the physical value $q=-e^{-i\pi\lambda}$ of (-@eq-q-physical). Requiring $K(\theta):V(\theta)\to V(-\theta)$ to intertwine $b_0$ and $b_1$ gives four linear equations for the four entries of $K$. Three are independent, and they determine $K$ up to an overall scalar (`code/boundary_quantum_groups.py`):

$$
K(\theta)\propto\begin{pmatrix}(q-q^{-1})(\hat\epsilon_1y+\hat\epsilon_0y^{-1})&y^2-y^{-2}\\y^2-y^{-2}&(q-q^{-1})(\hat\epsilon_0y+\hat\epsilon_1y^{-1})\end{pmatrix}.
$$ {#eq-sg-k-coideal}

This is the Ghoshal–Zamolodchikov matrix (-@eq-bsg-k-matrix), obtained here from a linear equation without solving the reflection equation. With the sine-Gordon assignment $(Q_+,\bar Q_-,T)\mapsto(e_1,f_1k_1^{-1},h_1)$ of @sec-qg-coproduct, the $e_1$ current pairs with the $C_0$ boundary term and the $e_0$ current with $C_1$, so $\hat\epsilon_1$ goes with $C_0$ and $\hat\epsilon_0$ with $C_1$, as noted in @sec-bsg-parameters. For $a_n^{(1)}$ with $n\ge2$ (@sec-an-nondiagonal-k) we label the coideal parameters directly by the nodes so that $\hat\epsilon_j$ corresponds to $C_j$.
```

---

## Item 2: Part VII Ch. 32 Singular-Wall Limits `0.719 / 0.937 / 0.147` (PR #22)

- **Verdict:** **correct**

### 2.1 Independent Derivation & Verification

In `referee_scratch/open_items/item2_singular_wall.py`, we implemented two independent solvers from scratch (without importing `code/solitonic_linear.py`):
- **Method A:** Direct numerical integration (`scipy.integrate.solve_ivp`, `DOP853`, `rtol=1e-12, atol=1e-13`) of the coupled $n$-channel linearized PDE on $x \in [-25, 0]$ with the boundary Hessian $H_B$ at $x=0$.
- **Method B:** Independent 70-digit `mpmath` solver constructing the Hirota $\tau$-functions and dressed linear waves $\delta\tau_j^{(c)}(x)$ at $x=0$ and solving the $n \times n$ boundary-condition linear system at $x=0$.

Results at $\theta = 0.6$ for all three cases in `part7-half-line/32-half-line-solitons.qmd` (lines 215–218):

1. **$a_5^{(1)}$ ($h=6$), sign pattern $P = (+,+,-,-,-,-)$, species-2 stationary pair along the regular branch $\chi = 5\pi/3 \pm \delta$ ($d_0 = 1/(1+\sin(\pi/3))$):**
   - At $\delta = -10^{-6}$: $|K_{1,5}| = 0.7191932055$, $K_{1,1} = 0.1556328298 + 0.6771554884i$.
   - At $\delta = +10^{-6}$: $|K_{1,5}| = 0.7191932055$, $K_{1,1} = 0.1556328298 + 0.6771554884i$.
   - Confirms the book's `|K_{1,5}| -> 0.719` and `K_{1,1} -> 0.15563 + 0.67716i` (and shows why direct ODE integration near $\delta \to 0$ must handle the $O(\delta^{-4})$ cancellation in $\delta\phi'(0) + H_B(0)\delta\phi(0)$ at the vanishing $\tau_2(0) \sim \delta^2$).
2. **$a_5^{(1)}$ ($h=6$), same endpoint along the non-solution path $\chi = 5\pi/3$, $d = d_0(1-\epsilon)$:**
   - Method A (ODE extrapolation from $\epsilon = 4\times 10^{-3}, 2\times 10^{-3}, 10^{-3}$): $|K_{1,5}| = 0.93731324$.
   - Method B ($\epsilon = 10^{-8}$): $|K_{1,5}| = 0.9374372024$, $K_{1,1} = -0.0646044139 + 0.3421078211i$.
   - Confirms the book's `|K_{1,5}| -> 0.937`.
3. **$a_4^{(1)}$ ($h=5$), sign pattern $P = (+,+,-,+,-)$, exceptional 4-soliton configuration ($\chi_1 = 9\pi/5, \chi_2 = 18\pi/5$):**
   - Solving $\tau_2(0) = 0$ and $\tau_0'(0)/\tau_0(0) = \tau_3'(0)/\tau_3(0)$ simultaneously gives:
     $$d_1^* = 0.257789025012, \qquad d_2^* = 0.563669884984.$$
   - At $d_1 = d_1^*(1 \mp 10^{-7}), d_2 = d_2^*$:
     $$|K_{1,4}| = 0.147136501, \qquad K_{1,1} = 0.3253473698 + 0.9340770520i.$$
   - Confirms the book's `|K_{1,4}| -> 0.147` and `K_{1,1} -> 0.32535 + 0.93408i`.
   - *Minor note on `code/solitonic_linear.py` line 376:* The initial guess passed to `mp.findroot` in `code/solitonic_linear.py` line 376 is `'0.2579'` (a slight typo for `0.2578`, though `mp.findroot` converges from `0.2579` to the true root `0.257789025012`).

**Output of `referee_scratch/open_items/item2_singular_wall.py`:**
```text
Case 1 (a_5 regular chi -> 5pi/3):
  Method B (Hirota mp, -1e-6): |K_15| = 0.7191932055, K_11 = (0.1556328298 + 0.6771554884j)
  Method B (Hirota mp, +1e-6): |K_15| = 0.7191932055, K_11 = (0.1556328298 + 0.6771554884j)
PASS case1_0.719

Case 2 (a_5 non-solution d -> 1/(1+s) at chi = 5pi/3):
  Method A (ODE extrap):       |K_15| = 0.93731324
  Method B (Hirota mp, 1e-8):  |K_15| = 0.9374372024, K_11 = (-0.06460441394 + 0.3421078211j)
PASS case2_0.937

Case 3 (a_4 exceptional 4-soliton limit of (+,+,-,+,-)):
  d1* = 0.257789025012, d2* = 0.563669884984
  Method B (Hirota mp, -1e-7): |K_14| = 0.147136501, K_11 = (0.3253473698 + 0.934077052j)
  Method B (Hirota mp, +1e-7): |K_14| = 0.1471485334, K_11 = (0.325382709 + 0.9341503882j)
PASS case3_0.147
```

---

## Item 3: Part VI `code/rmatrix/smatrix/oneloop_ratio.py` $a_5^{(2)}$ Quirk (`0.693` vs `0.5`, PR #29)

- **Verdict:** **wrong** (bug in `code/rmatrix/smatrix/oneloop_ratio.py` line 37, plus a false prose statement in `part5-semiclassics/20-one-loop-masses.qmd` line 97)

### 3.1 Root-Cause Diagnosis & Verification

In `referee_scratch/open_items/item3_oneloop_ratio.py`, we diagnosed why `code/rmatrix/smatrix/oneloop_ratio.py` prints `delta H / B_tree = 0.693181` (and a non-constant `tree amplitude / g(theta) = [0.065502, 0.053167, 0.054227, 0.056704]`) for $a_5^{(2)}$:

1. **Particle ordering in `toda_data(a5, [1, 2, 1, 1])`:**
   For $a_5^{(2)}$ ($r=3, h=5, h^\vee=6$), the classical masses are:
   - Particle $a=1$: $m_1 = 2\sin(\pi/5) \approx 1.1755705$ ($m_1^2 \approx 1.381966$),
   - Particle $a=2$: $m_2 = 2\sin(2\pi/5) \approx 1.9021130$ ($m_2^2 \approx 3.618034$),
   - Particle $n=3$: $m_3 = 1.0$ ($m_3^2 = 1.0$).
   Because $m_3 = 1 < m_1 = 2\sin(\pi/5) \approx 1.17557$, `toda_data` (which sorts eigenvalues ascending via `np.linalg.eigh`) places **particle $3$ at index `0`** and **particle $1$ at index `1`**!
2. **The bug in `oneloop_ratio.py` line 37:**
   In `analyse()`, line 37 calls `tree_over_4m2sinh(m2, C3, C4, 0, th)` (hardcoding index `0`, which is **particle $3$** for $a_5^{(2)}$!) and divides it by the $O(B)$ slope $g(\theta)$ of `S11 = dgz11(h + B/2, B, th) = brace(1, B, H, th) * brace(H-1, B, H, th)`, which is the tree-level elastic S-matrix of **particle $1$**!
   - Dividing the tree amplitude of particle $3$ by the S-matrix slope of particle $1$ produces the non-constant ratio `[0.065502, 0.053167, 0.054227, 0.056704]` and the spurious mean `0.693181`!
3. **The fix:**
   Pass the eigen-index `p1` of particle $1$ (`p1 = 0` for $c_3^{(1)}$ and $f_4^{(1)}$, and `p1 = 1` for $a_5^{(2)}$) into `tree_over_4m2sinh(m2, C3, C4, p1, th)`. When evaluated at `p1 = 1` for $a_5^{(2)}$:
   - `tree amplitude / g(theta)` at $\theta = [0.4, 0.9, 1.5, 2.3]$ is **identically constant**: `[0.079577, 0.079577, 0.079577, 0.079577]` ($= \frac{1}{4\pi}$),
   - `delta H / B_tree = 0.500000` (exact $1/2$ to machine precision, matching $H = h + B/2 = 5 + B/2$ with $B_{\rm tree} = \frac{\beta^2}{4\pi}$ and $\delta H = \frac{\beta^2}{8\pi}$ in `oneloop_ratio.py`'s root normalization where long roots of $a_5^{(2)}$ have $\alpha_{\rm long}^2 = 4$).
4. **Additional finding in `part5-semiclassics/20-one-loop-masses.qmd` line 97:**
   Line 97 of `20-one-loop-masses.qmd` states:
   > *"The shifts of all mass ratios vanish at one loop in two dual families as well: $a_{2n-1}^{(2)}$, which folds from $a_{2n-1}^{(1)}$ with equal-mass orbits, and $d_{n+1}^{(2)}$..."*
   This statement is **false** for $a_{2n-1}^{(2)}$ ($n \ge 3$): as shown in `14-non-simply-laced.qmd` (`tbl-nsl-data`), `code/vernier/oneloop_particles.py`, and `code/rmatrix/smatrix/oneloop_ratio.py`, the masses of $a_{2n-1}^{(2)}$ are $m_a = \sqrt{2}m\sin(\pi a/H)$ ($a=1,\dots,n-1$) and $m_n = m/\sqrt{2}$ with $H = 2n - 1 - B/2$, so $m_a/m_n = 2\sin(\pi a/H)$ **shifts at one loop** with $\delta H = -\frac{\beta^2}{8\pi} \ne 0$! (Even ratios $m_a/m_1 = \sin(\pi a/H)/\sin(\pi/H)$ shift at one loop for $n \ge 3$.) Why? Because in $a_{2n-1}^{(1)}$ ($h=2n$), the folded modes have $m_a = 2m\sin(\pi a/2n)$ while the odd modes dropped by folding have different masses, so the one-loop self-energy diagrams do NOT preserve the classical mass ratios! Only $d_{n+1}^{(2)}$ (among the dual families) has vanishing one-loop mass-ratio shifts.

**Output of `referee_scratch/open_items/item3_oneloop_ratio.py`:**
```text
Particle 1 (index 1) with brace(1, nu=0.0)*brace(H-1, nu=0.0):
  tr1 / g11 = [0.079577 0.079577 0.079577 0.079577]
  delta H / B_tree = 0.500000
PASS a5_2_particle1_index_fix
```

### 3.2 Exact Replacement Code & Text

1. **`code/rmatrix/smatrix/oneloop_ratio.py` (lines 28–55):**
```python
def analyse(name, roots, kac, h, ratio_fn, S11, expected, p1=0):
    m2, C3, C4 = toda_data(roots, kac)
    dm2 = one_loop_dm2(m2, C3)
    dH = []
    for a in range(1, len(m2)):
        dR = (dm2[a]*m2[0] - m2[a]*dm2[0])/m2[0]**2
        eps = 1e-6
        slope = (ratio_fn(a, h + eps) - ratio_fn(a, h - eps))/(2*eps)
        dH.append(dR/slope)
    ths = [0.4, 0.9, 1.5, 2.3]
    tr = [tree_over_4m2sinh(m2, C3, C4, p1, th) for th in ths]
    eps = 1e-7
    g = [np.log(S11(h, eps, th)).imag/eps for th in ths]
    kap = np.array(tr)/np.array(g)
    ratio = np.mean(dH)/np.mean(kap)
    print(f"{name}: classical mass ratios {np.round(np.sqrt(m2/m2[0]), 4)}")
    print(f"   one-loop delta H from each mass ratio: {np.round(dH, 6)}")
    print(f"   tree amplitude / g(theta) at theta = {ths}: {np.round(kap, 6)}  (constant => same tree shape)")
    print(f"   delta H / B_tree = {ratio:.6f}   (expected {expected})")
    ok = abs(ratio - expected[0]) < 1e-4 and np.ptp(kap) < 1e-6
    print(f"{'PASS' if ok else 'FAIL'} {name}_oneloop_ratio")


# c_3^(1): alpha_1=e1-e2, alpha_2=e2-e3, alpha_3=2e3, alpha_0=-2e1; Kac (1,2,2,1), h=6
e = np.eye(3)
c3 = [e[0]-e[1], e[1]-e[2], 2*e[2], -2*e[0]]
analyse("c_3^(1)", c3, [1, 2, 2, 1], 6,
        lambda a, H: np.sin(np.pi*(a+1)/H)/np.sin(np.pi/H),
        lambda h, B, th: dgz11(h + B, B, th), (1.0, "1 (H = h + B)"), p1=0)

# a_5^(2): dual of c_3^(1), Simple roots alpha_1=(e1-e2)/sqrt2, alpha_2=(e2-e3)/sqrt2, alpha_3=sqrt2*e3,
# alpha_0=-(e1+e2)/sqrt2; Kac (1,2,1,1), h=5. Note: sorted m2 has particle 3 (m=1) at index 0 and particle 1 at p1=1.
s2 = 1/np.sqrt(2)
a5 = [s2*(e[0]-e[1]), s2*(e[1]-e[2]), s2*2*e[2], -s2*(e[0]+e[1])]
def a5ratio(a, H):
    ms = sorted([2*np.sin(np.pi*1/H), 2*np.sin(np.pi*2/H), 1.0])
    return (ms[a]/ms[0])**2
analyse("a_5^(2)", a5, [1, 2, 1, 1], 5, a5ratio,
        lambda h, B, th: dgz11(h + B/2, B, th), (0.5, "1/2 (H = h + B/2)"), p1=1)

# f_4^(1): standard roots in R^4, Kac (1,2,3,4,2), h=12
e = np.eye(4)
al1 = e[1] - e[2]; al2 = e[2] - e[3]; al3 = e[3]; al4 = 0.5*(e[0] - e[1] - e[2] - e[3])
al0 = -(e[0] + e[1])
f4 = [al0, al1, al2, al3, al4]
f4ratios_fn = lambda a, H: ([1, 2*np.cos(np.pi/4), 2*np.cos(np.pi/H),
                              2*np.cos(np.pi/4)*2*np.cos(np.pi*(H-4)/(2*H))][a])**2
analyse("f_4^(1)", f4, [1, 2, 3, 4, 2], 12, f4ratios_fn,
        lambda h, B, th: dgz11(h + 3*B, B, th),
        (3.0, "3 if CDS (H = 12 + 3B), 9/2 if the soliton construction"), p1=0)
```

2. **`code/rmatrix/README.md` (line 42):**
Replace:
```markdown
| `smatrix/oneloop_ratio.py` | 26 | One-loop $\delta H / B_{\text{tree}}$ ratio for $c_3^{(1)}$, $a_5^{(2)}$, $f_4^{(1)}$: gets $3.000$ for $f_4^{(1)}$, matching CDS ($H = 12 + 3B$) rather than the soliton construction's $9/2$ |
```
with:
```markdown
| `smatrix/oneloop_ratio.py` | 26 | One-loop $\delta H / B_{\text{tree}}$ ratio for $c_3^{(1)}$ ($1.000$), $a_5^{(2)}$ ($0.500$), and $f_4^{(1)}$ ($3.000$, matching CDS $H = 12 + 3B$ rather than the soliton construction's $9/2$) |
```

3. **`part5-semiclassics/20-one-loop-masses.qmd` (lines 97–98):**
Replace:
```markdown
The shifts of all mass ratios vanish at one loop in two dual families as well: $a_{2n-1}^{(2)}$, which folds from $a_{2n-1}^{(1)}$ with equal-mass orbits, and $d_{n+1}^{(2)}$, which folds from $d_{n+2}^{(1)}$ by exchanging the two spinor nodes. The folded one-loop code in `code/vernier/oneloop_particles.py` checks both (and the non-zero shift in $a_4^{(2)}$) for small rank.
```
with:
```markdown
The shifts of all mass ratios vanish at one loop in the dual family $d_{n+1}^{(2)}$ as well, which folds from $d_{n+2}^{(1)}$ by exchanging the two spinor nodes. By contrast, in $a_{2n-1}^{(2)}$ (folded from $a_{2n-1}^{(1)}$) and $a_{2n}^{(2)}$ (folded from $d_{2n+2}^{(1)}$) the mass ratios shift at one loop because the parent modes dropped by folding do not share the masses of the retained orbits; for $a_{2n-1}^{(2)}$ the shifts follow $H=2n-1-B/2$ (@tbl-nsl-data). The folded one-loop code in `code/vernier/oneloop_particles.py` checks $d_{n+1}^{(2)}$, $a_{2n-1}^{(2)}$ and $a_4^{(2)}$ for small rank.
```

---

## Item 4: PR #7 Open Follow-Ups

### Item 4a: `appendices/A-algebra-data.qmd` (`tbl-app-twisted` parents of $a_{2n}^{(2)}$)

- **Verdict:** **correct** (already resolved in `HEAD`)
- **Evidence:** In `appendices/A-algebra-data.qmd`, line 25 (`tbl-app-twisted`) already lists `$a_{2n}^{(1)}$ or $d_{2n+2}^{(1)}$` in the `parent` column for $a_{2n}^{(2)}$, and line 49 explicitly explains both foldings:
  > *"Why does $a_{2n}^{(2)}$ have two parents in @tbl-app-twisted? An outer automorphism of $a_{2n}$ of order $2$ does not fix the affine root, so it twists the loop variable and produces $a_{2n}^{(2)}$ with $h=4n+2$. A diagram automorphism of $d_{2n+2}^{(1)}$ of order $4$ gives the same affine algebra with the same $h$. Only the second folding descends to classical Toda solutions without a shift in $x$ (@sec-affine-folding)."*
  No edit is needed.

### Item 4b: `code/irf/checks.py` and `code/irf/tm2.py` Breather Threshold Condition

- **Verdict:** **wrong** (scripts are fixed, but 5 stored output files in `code/irf/` are stale and contradict Ch. 28 and the scripts)
- **Evidence:**
  - In `part6-soliton-smatrices/28-rsos.qmd` (line 66), `code/irf/checks.py` (line 47), and `code/irf/tm2.py` (line 15), the breather threshold condition for $a_2^{(1)}$ RSOS models $W_3(p, q)$ ($r = q - p$) is `3*r > 2*p` ($\xi = \frac{\pi r}{p - r} < 2\pi$, i.e. $\frac{p}{q - p} < \frac{3}{2}$).
  - For $W_3(5, 7)$ ($p=5, r=2$), $3r = 6 \le 2p = 10$ ($\xi = \frac{2\pi}{3}$ has $2\pi/\xi = 3$, the two-soliton threshold, so there are **no** bound-state breathers below threshold). For $W_3(8, 11)$ ($p=8, r=3$), $3r = 9 \le 2p = 16$ (**no** breathers). For $W_3(4, 7)$ ($p=4, r=3$), $3r = 9 > 2p = 8$ (**yes**, breathers exist).
  - However, 5 stored output files in `code/irf/` were generated before the `3*r > p` $\to$ `3*r > 2*p` fix and still say `breathers=yes` for $W_3(5,7)$ and $W_3(8,11)$:
    1. `code/irf/out_ha.txt` (lines 9 and 11),
    2. `code/irf/out_qsign/irf_ha_q1.txt` (lines 9 and 11),
    3. `code/irf/out_qsign/irf_ha_q-1.txt` (lines 9 and 11),
    4. `code/irf/out_qsign/irf_tm2_q1.txt` (line 5),
    5. `code/irf/out_qsign/irf_tm2_q-1.txt` (line 5).
- **Recommended fix:** Regenerate those 5 output files (or update lines 9 and 11 of `out_ha.txt`, `out_qsign/irf_ha_q1.txt`, `out_qsign/irf_ha_q-1.txt` from `breathers=yes` to `breathers=no`, and line 5 of `out_qsign/irf_tm2_q1.txt`, `out_qsign/irf_tm2_q-1.txt` from `breathers=yes` to `breathers=no`).

### Item 4c: `code/thimbles/check_th4.py`

- **Verdict:** **correct** (already resolved in `HEAD`)
- **Evidence:**
  - `code/thimbles/check_th4.py` (lines 8–9) now checks `if not os.path.exists("cps_a2L1_k1.00.npy"): subprocess.run([sys.executable, "cps_a2L2.py", "1.0", "1"], check=True)`, so it is completely self-contained and does not import `check_th2.py`.
  - In `code/thimbles/out_check_th4.txt` (lines 28–29), the Stokes jump prediction and the resummed residual have matching signs and magnitudes:
    - $\hbar = 0.15$: `res = -3.0976e-07`, `pred = -3.0315e-07` (`ratio = +1.0218`),
    - $\hbar = 0.10$: `res = +1.7710e-10`, `pred = +2.0053e-10` (`ratio = +0.8832`).
  - `code/README.md` (line 78) accurately documents `check_th4.py`. No edit is needed.

---

## Item 5: Appendix E (`appendices/E-original-results.qmd`) Audit (PR #18, PR #22, PR #46)

- **Verdict:** **wrong** (1 status-label mismatch hidden by a `@sec-` reference on line 45, 1 `[Established]` literature-correction bullet on line 16, and 14 original propositions in Parts VII–VIII missing from Appendix E)

### 5a. Manual Audit of All Section-Referenced Items in `appendices/E-original-results.qmd`

Because `check-refs.py` (lines 240–248) only checks the status label of an Appendix E bullet when it cites a statement id (`prp-`, `thm-`, `lem-`, `cnj-`), any bullet that cites only a `@sec-` id bypasses `check-refs.py`. We audited all 6 such bullets manually:

| Line in `E-original-results.qmd` | Appendix E Status Label | Cited `@sec-` ID(s) | Chapter Statement ID & Label | Audit Verdict |
|---|---|---|---|---|
| **Line 16** | `**[Established]**` | `@sec-one-loop-masses-folded` | None (prose section in `20-one-loop-masses.qmd`, lines 78–103) | **Mismatch / Wrong place:** Appendix E is the list of *original* results of the book (`AGENTS.md`), whereas `[Established]` means established in earlier literature. The original content of `@sec-one-loop-masses-folded` is the correction to Braden–Corrigan–Dorey–Sasaki [@bcds1990, Table 2] for $b_n^{(1)}, c_n^{(1)}, f_4^{(1)}$ (which is already in Appendix D) and the folded one-loop mass check. Either drop from Appendix E (since the correction is in Appendix D) or relabel `**[Numerical]**` if kept for the folded one-loop mass verification. |
| **Line 20** | `**[Numerical]**` | `@sec-sm-cn, @sec-sm-a2n, @sec-sm-e6` | `@prp-cn-smatrix`, `@prp-a2n-smatrix`, `@prp-e6-smatrix` (all `**[Numerical]**`) | **Correct, but should cite `@prp-cn-smatrix, @prp-a2n-smatrix, @prp-e6-smatrix`** so `check-refs.py` enforces label agreement. |
| **Line 22** | `**[Numerical]**` | `@sec-sm-f4` | `@prp-f4-r`, `@prp-f4-breathers` (both `**[Numerical]**`) | **Correct, but should cite `@prp-f4-r, @prp-f4-breathers`** alongside `@sec-sm-f4`. |
| **Line 36** | `**[Numerical]**` | `@prp-an-nondiagonal-scalar, @sec-neumann-duality` | `@prp-an-nondiagonal-scalar` (`**[Numerical]**`), `@cnj-neumann-amplitude` (`**[Conjecture]**`) | **Correct:** `@prp-an-nondiagonal-scalar` is `**[Numerical]**`, and the bullet prose explicitly states "under `@cnj-reflection-continuation` and `@cnj-neumann-coideal`". |
| **Line 45** | `**[Theorem for $B_1$; Numerical for $B_n$]**` | `@sec-sg-breather-transmission` | `@prp-sg-breather-transmission` (`**[Established for $B_1$, with a new proof; Numerical for $n\ge2$]**`, `38-defect-sine-gordon.qmd` line 335) | **WRONG:** The first word in Appendix E (`Theorem`) contradicts the first word of `@prp-sg-breather-transmission` (`Established`)! Citing `@sec-sg-breather-transmission` hid this mismatch from `check-refs.py`. Must be updated to `**[Established for $B_1$, with a new proof; Numerical for $n\ge2$]**` and cite `@prp-sg-breather-transmission`. |
| **Line 60** | `**[Theorem for the two-colour sectors and the $a_2^{(2)}$ two-kink asymptotics; Numerical for the remaining sectors]**` | `@sec-sector-status` | `@prp-sector-table` / table in `49-scattering.qmd` | **Correct:** Accurately summarizes the `[Theorem]` and `[Numerical]` rows of `@sec-sector-status`. |

### 5b. Candidate Original Results in Parts VII & VIII (PR #18, PR #22, PR #40, PR #46)

We checked each candidate statement against its chapter text and `Notes and references` section to verify whether the book claims it as original:

#### Part VII Candidates (all 5 are original results and should be added to Appendix E under `## Part VII`):
1. **`@prp-soliton-preserving-half`** (`31-boundary-conditions.qmd` line 260, `**[Theorem]**`): Half-soliton generalization $\operatorname{Re}\vec\phi(0)\in\frac\pi\beta\vec\nu+\frac{2\pi}\beta\Lambda_W^\vee$ with $K_+=e^{i\pi\vec\nu\cdot\vec H}$, correcting Delius's $\vec\lambda\in\Lambda_W^\vee$ formula (`31-boundary-conditions.qmd` line 362: *"In [@delius1998b, eq. (3.16)], $\vec\phi(x)=-\overline{\vec\phi(-x)}-\frac{2\pi}\beta\vec\lambda$ with $\vec\lambda\in\Lambda_W^\vee$ gives $\operatorname{Re}\vec\phi(0)=-\frac\pi\beta\vec\lambda$, which is integrable with $K_+=1$ only for $\vec\lambda\in2\Lambda_W^\vee$; for general $\vec\lambda$ it requires $K_+=e^{i\pi\vec\lambda\cdot\vec H}$ (@prp-soliton-preserving-half)."*).
2. **`@prp-dn-uniform-vacuum` and `@prp-dn-uniform-classical`** (`34-soliton-reflection.qmd` lines 253 & 269, `**[Theorem]**` and `**[Numerical]**`): Explicit $d_n^{(1)}$ uniform-boundary reflection matrix $K_1(\theta) = \sigma_{\rm unif}(\theta)\operatorname{diag}(P_1,\dots,P_n,Q_n,\dots,Q_1)$, its four simple poles in the strip, exclusion of the $2h-2$ pole by the uniform vacuum state, and its classical limit matching the linearized reflection factor (`34-soliton-reflection.qmd` line 330: *"the explicit $d_n$ uniform amplitude and its classical matching (@sec-dn-soliton-reflection) are new"*).
3. **`@prp-dn-nondiagonal-scalar`** (`34-soliton-reflection.qmd` line 298, `**[Numerical]**`): Scalar factor $\sigma_{\rm nd}(\theta)$ of the non-diagonal $d_n^{(1)}$ vector soliton reflection matrix, closing the two-particle bootstrap on all $2n$ diagonal channels (`34-soliton-reflection.qmd` line 330).
4. **`@prp-sp-scalar`** (`34-soliton-reflection.qmd` line 222, `**[Numerical]**`): Soliton-preserving scalar factors $\sigma_a^{\rm sp}(\theta)$ for $a_n^{(1)}$ ($a=1,\dots,n$) satisfying the reflection, unitarity, crossing-unitarity, and fusion bootstrap equations (`34-soliton-reflection.qmd` line 330: *"The soliton-preserving scalar factors (@prp-sp-scalar)... are new"*).

#### Part VIII Candidates (all 9 are original results and should be added to Appendix E under `## Part VIII`):
1. **`@prp-twist-defect`** (`37-defect-scattering.qmd` line 183, `**[Theorem, relative to @prp-bloch-sectors]**`): Construction of the topological twist defect $\mathbb T_\ell$, its centrality, fusion $\mathbb T_\ell\mathbb T_{\ell'}=\mathbb T_{\ell+\ell'}$, and diagonalization on the Bloch sectors (`37-defect-scattering.qmd` line 284: *"The affine Toda specialization of the twist defect (@prp-twist-defect)... are new."*).
2. **`@prp-sg-spin-three`** (`38-defect-sine-gordon.qmd` line 62, `**[Theorem, by hand; numerically checked]**`): Explicit spin-$+3$ defect charge density $\mathcal D_{+3}$ and proof that the type-I sine-Gordon defect conserves $P_{-3}$ only at the trivial points $\sigma\to0,\infty$ (`38-defect-sine-gordon.qmd` line 402: *"the spin-3 calculation of @prp-sg-spin-three... are new"*).
3. **`@prp-dc-gauge`** (`39-defect-conditions.qmd` line 156, `**[Theorem for type I$_-$; Numerical for type I$_+$]**`): Gauge equivalence of all $n!$ standard cyclic defect permutations $E = \sum_j d_j e_{\pi(j),\pi(j+1)}$ of $a_{n-1}^{(1)}$ via the affine Weyl group (`39-defect-conditions.qmd` line 284: *"The gauge-equivalence theorem for the cyclic permutations (@prp-dc-gauge)... are new"*).
4. **`@prp-defect-theta`** (`39-defect-conditions.qmd` line 245, `**[Theorem]**`): Necessary and sufficient conditions for $\Theta$-invariance of the type-I defect action and reality of the defect energy on $\Theta$-symmetric solutions (`39-defect-conditions.qmd` line 284: *"the $\Theta$-invariance criteria (@prp-defect-theta, @prp-defect-theta-p) are new"*).
5. **`@prp-ds-absorption`** (`40-defect-solitons.qmd` line 122, `**[Theorem]**`) and **`@prp-defect-frozen-soliton`** (`40-defect-solitons.qmd` line 161, `**[Theorem]**`): Exact characterization of soliton absorption at $\theta = \theta_d - i\pi/2$, the $h/\gcd(a,h)$ absorbed topological charges, and identification of the static defect step with a frozen bulk soliton at $x_0 = \theta_d/m_a$ (`40-defect-solitons.qmd` line 264: *"The proof that the absorbed charges are the $h/\gcd(a,h)$ Hirota charges (@prp-ds-absorption), the identification of the static defect step with a frozen bulk soliton (@prp-defect-frozen-soliton)... are new."*).
6. **`@prp-fused-bdg`** (`43-defects-and-boundaries.qmd` line 169, `**[Theorem, by hand; numerically checked]**`): Classical dressing $K_\phi(\lambda) = \sigma(L_d(1/\lambda)^{-1})K_\psi(\lambda)L_d(\lambda)^{-1}$ for fused defect–boundary systems across @tbl-durham, and the exact comparison between a defect fused onto a Dirichlet/Neumann boundary and the Baseilhac–Delius–George oscillator boundary (`43-defects-and-boundaries.qmd` line 255: *"the comparison with the oscillator boundary (@prp-fused-bdg)... are new"*).
7. **`@prp-an-particle-transmission` and `@prp-dpt-vacuum`** (`44-particle-transmission.qmd` lines 78 & 135, both `**[Theorem]**`): Proof that the $a_n^{(1)}$ particle transmission amplitudes from the soliton-breather bootstrap equal the classical linearized transmission amplitudes at $\theta_d = (1 - B/2)\theta_d^{\rm cl}$ (with no higher-loop quantum corrections beyond $\theta_d$), and the exact match between the dressed vacuum transmission matrix on the $h/\gcd(b,h)$ defect steps and the classical soliton-step transmission matrix (`44-particle-transmission.qmd` line 179: *"The observation that the bootstrap and linearized formulas are the same rational function of $e^{\theta-\theta_d}$ (@prp-an-particle-transmission) and the exact agreement of the dressed vacuum transmission matrix with the classical soliton-step calculation across all $h/\gcd(b,h)$ charges (@prp-dpt-vacuum) are new."*).

### 5c. Consistency of Literature Corrections Between Appendix D and Appendix E

Per `AGENTS.md`:
- Corrections to published papers belong in each chapter's Notes and in **Appendix D (`appendices/D-corrections.qmd`)**.
- **Appendix E (`appendices/E-original-results.qmd`)** lists the book's **original results**, each opening with its status label and ending with the statement ids it summarizes.
- Currently in `appendices/E-original-results.qmd`:
  1. **Line 16 (`@sec-one-loop-masses-folded`)** is tagged `**[Established]**` and describes the correction to [@bcds1990, Table 2] for $b_n^{(1)}, c_n^{(1)}, f_4^{(1)}$. Since that correction is already in `appendices/D-corrections.qmd` (line 20) and `[Established]` is not an original-result tag, **line 16 of `appendices/E-original-results.qmd` should be dropped**.
  2. **Line 28 (`@prp-dn-vector-r`)** and **Line 62 (`@prp-rsos-kinks`)** are original mathematical results of the book (`[Theorem, by hand; numerically checked]` for the $d_n^{(1)}$ vector R-matrix Spectral decomposition and `[Numerical]` for the RSOS two-kink test) that *also* mention correcting [@delius_george2002] or refuting [@tw1999]; those two bullets are legitimate original results (and their paper corrections are also recorded in Appendix D).

### 5d. Exact Replacement Text for `appendices/E-original-results.qmd`

1. **Drop line 16** (`- **[Established]** In the untwisted non-simply-laced theories $b_n^{(1)}$... (@sec-one-loop-masses-folded).`).
2. **Update lines 20–22** to include the proposition ids:
```markdown
- **[Numerical]** Full two-soliton S-matrices, with scalar factors and breathers, for $c_n^{(1)}$ ($n=2,3$), $a_{2n}^{(2)}$ ($n=2,3$) and $e_6^{(1)}$ (@prp-cn-smatrix, @prp-a2n-smatrix, @prp-e6-smatrix, @sec-sm-cn, @sec-sm-a2n, @sec-sm-e6).
- **[Numerical]** The 26-dimensional fundamental R-matrix and soliton S-matrix of $f_4^{(1)}$ (@prp-f4-r, @prp-f4-breathers, @sec-sm-f4).
```
3. **Add the 5 missing Part VII bullets** after line 37 (end of `## Part VII`):
```markdown
- **[Theorem]** Extension of the soliton-preserving boundary condition to half-soliton boundary values $\operatorname{Re}\vec\phi(0)\in\frac\pi\beta\vec\nu+\frac{2\pi}\beta\Lambda_W^\vee$ with $K_+=e^{i\pi\vec\nu\cdot\vec H}$, and classification of the coweight classes realized by static half-solitons (@prp-soliton-preserving-half).
- **[Numerical]** Scalar factors $\sigma_a^{\rm sp}(\theta)$ for soliton-preserving reflection of $a_n^{(1)}$ solitons ($a=1,\dots,n$), solving the reflection, unitarity, crossing-unitarity and fusion bootstrap equations (@prp-sp-scalar).
- **[Theorem]** Explicit vector-soliton reflection matrix $K_1(\theta)$ of $d_n^{(1)}$ on the uniform boundary $C\equiv+1$, its four simple poles in the physical strip, and the exclusion of the $(2h-2)/2h$ pole on the uniform vacuum (@prp-dn-uniform-vacuum).
- **[Numerical]** Classical limit of the $d_n^{(1)}$ uniform soliton reflection matrix and its agreement with the linearized reflection factor on all particle species (@prp-dn-uniform-classical).
- **[Numerical]** Scalar factor $\sigma_{\rm nd}(\theta)$ of the non-diagonal $d_n^{(1)}$ vector-soliton reflection matrix, closing the two-particle bootstrap on all $2n$ diagonal channels (@prp-dn-nondiagonal-scalar).
```
4. **Update line 45** and **add the 9 missing Part VIII bullets** in `## Part VIII`:
```markdown
- **[Theorem, relative to @prp-bloch-sectors]** Construction of the topological twist defect $\mathbb T_\ell=\sum_{\vec m}e^{i\vec\ell\cdot\vec m}\Pi_{\vec m}$, its centrality in the bulk algebra, and its diagonalization on the Bloch sectors (@prp-twist-defect).
- **[Theorem, by hand; numerically checked]** Explicit spin-$+3$ defect contribution $\mathcal D_{+3}$ for the type-I sine-Gordon defect and proof that spin-$-3$ momentum conservation holds only at the trivial limits $\sigma\to0,\infty$ (@prp-sg-spin-three).
- **[Established for $B_1$, with a new proof; Numerical for $n\ge2$]** Derivation of the sine-Gordon breather transmission factors by fusing the defect onto a Dirichlet boundary (@prp-sg-breather-transmission, @sec-sg-breather-transmission).
- **[Theorem for type I$_-$; Numerical for type I$_+$]** Gauge equivalence of all $n!$ cyclic defect permutations of $a_{n-1}^{(1)}$ under the affine Weyl group (@prp-dc-gauge).
- **[Theorem]** Necessary and sufficient condition for $\Theta$-invariance of the type-I defect action and reality of the defect energy on $\Theta$-symmetric solutions (@prp-defect-theta).
- **[Theorem]** Characterization of classical soliton absorption by a type-I $a_n^{(1)}$ defect at $\theta=\theta_d-i\pi/2$ on the $h/\gcd(a,h)$ Hirota charges, and identification of each static defect step with a frozen bulk soliton (@prp-ds-absorption, @prp-defect-frozen-soliton).
- **[Theorem, by hand; numerically checked]** Classical dressing $K_\phi(\lambda)=\sigma(L_d(1/\lambda)^{-1})K_\psi(\lambda)L_d(\lambda)^{-1}$ for a type-I defect fused onto any integrable $a_r^{(1)}$ boundary of @tbl-durham, and exact comparison with the Baseilhac–Delius–George oscillator boundary (@prp-fused-bdg).
- **[Theorem]** Identity between the quantum bootstrap particle transmission amplitudes of $a_n^{(1)}$ and the classical linearized transmission amplitudes at $\theta_d=(1-B/2)\theta_d^{\rm cl}$, and between the dressed vacuum transmission matrices and the classical transmission matrices across all $h/\gcd(b,h)$ static soliton steps (@prp-an-particle-transmission, @prp-dpt-vacuum).
```

---

## Item 6: Part IX Ch. 48 (`48-continuum-window.qmd`) Step 3 of `thm-nn-distances` (PR #5 & PR #16)

- **Verdict:** **correct** (with a citation update on line 109 of `part9-foundations/48-continuum-window.qmd`)

### 6.1 Independent Mathematical & Computational Verification of Steps 3.1–3.5

In `referee_scratch/open_items/item6_step3_nn.py`, we verified every step of Step 3 (`part9-foundations/48-continuum-window.qmd`, lines 66–109) independently from scratch:

1. **Step 3.1 (Spanning tree and dendrogram):**
   - For any finite point set $A \subset \mathbb R^2$ in general position, the Kruskal minimum spanning tree $\tau(A)$ has the minimax path property: the maximum edge weight along the unique path in $\tau(A)$ between $a$ and $b$ equals both the minimax path distance in the complete Euclidean graph on $A$ and the single-linkage merger height $s_{a \vee b}$ in the dendrogram $T(A)$. Since the path in $\tau(A)$ has at most $n-1$ edges, each of length $\le s_{\rm root}$, $\operatorname{diam}(A) \le (n-1) s_{\rm root}$.
   - Verified on 200 random $n=6$ configurations in $\mathbb R^2$ (`PASS step_3_1_minimax_and_diam`).
2. **Step 3.2 (Lower bound on $\tilde\delta_k$):**
   - Let $i(k)$ be the lowest ancestor of leaf $k$ in $T(A)$ whose leaf cluster $L_{i(k)}$ has size $\ge m+1$. Just before merger $i(k)$, leaf $k$'s cluster $L'$ has $|L'| \le m$, so $A$ contains at most $|L'| - 1 \le m-1$ other points at single-linkage distance $< s_{i(k)}$ from $k$. Since Euclidean distance is $\ge$ single-linkage distance, the $m$-th nearest neighbour distance $\tilde\delta_k$ satisfies $\tilde\delta_k \ge s_{i(k)}$. Summing over $k \in I$ gives $\prod_{k \in I} \tilde\delta_k^{-\gamma} \le \prod_j s_j^{-\gamma c_j}$ with $\sum_j c_j = n_I \le n$.
   - Verified on 200 random configurations (`PASS step_3_2_ancestor_bound`).
3. **Step 3.3 (Linear extensions of the dendrogram poset):**
   - Once the labelled spanning tree $\tau$ and the dendrogram topology $T$ are fixed, each internal node $j$ of $T$ corresponds to the unique edge of $\tau$ linking its two child clusters $L_{\text{left}(j)}$ and $L_{\text{right}(j)}$. An ordering of the $n-1$ edges of $\tau$ produces $T$ if and only if every internal node appears after its descendants in $T$, i.e. the ordering is a linear extension of the tree poset $T$. By Knuth's hook-length formula for tree posets, the number of linear extensions is $(n-1)! / \prod_{j=1}^{n-1} h_j$, where $h_j = |L_j| - 1$ is the number of internal nodes in the subtree of $j$.
   - Verified exhaustively across all $n^{n-2}$ Prüfer-coded labelled trees and all $(n-1)!$ edge permutations for $n = 2, 3, 4, 5$ (`PASS step_3_3_linear_extensions`).
4. **Step 3.4 (Nested cone integral and exponent bound):**
   - Integrating $\prod_{j=1}^{n-1} s_j^{1 - \gamma c_j}$ bottom-up over the tree cone $R(T) = \{0 < s_j < s_{\text{parent}(j)}, s_{\rm root} < \infty\}$: a leaf-pair internal node $j$ ($h_j = 1, C_j = c_j$) integrates to $\int_0^{s_p} s_j^{1-\gamma c_j}\,ds_j = s_p^{2 - \gamma c_j} / E_j$ with $E_j = 2h_j - \gamma C_j$. By induction up the tree, the powers inherited by node $j$ from its children sum to $\sum_{c \in \text{children}(j)} (2h_c - \gamma C_c) = 2(h_j - 1) - \gamma(C_j - c_j)$, so the integrand at node $j$ has exponent $1 - \gamma c_j + 2(h_j - 1) - \gamma(C_j - c_j) = 2h_j - 1 - \gamma C_j = E_j - 1$, integrating to $s_{\text{parent}(j)}^{E_j} / E_j$.
   - **Why $C_j \le h_j + m$ and $2h_j / E_j \le 2m/\alpha$:** Every leaf $k$ contributing to $C_j$ lies in $L_j$ ($|L_j| = h_j + 1$), so $C_j \le |L_j| = h_j + 1 \le h_j + m$ (for $m \ge 1$). Moreover, if $C_j > 0$, then for $\gamma = \frac{2m-\alpha}{m+1} \ge 0$ (with $0 < \alpha \le 2m$):
     $$E_j = 2h_j - \gamma C_j \ge 2h_j - \gamma(h_j + 1) = h_j\left(2 - \frac{2m-\alpha}{m+1}\right) - \frac{2m-\alpha}{m+1}.$$
     Even using the weaker bound $C_j \le h_j + m$ stated in the book:
     $$\frac{E_j}{2h_j} \ge 1 - \frac{\gamma(h_j + m)}{2h_j} = 1 - \frac{2m-\alpha}{2(m+1)}\left(1 + \frac{m}{h_j}\right) \ge 1 - \frac{2m-\alpha}{2(m+1)}(1 + m) = 1 - \frac{2m-\alpha}{2} \cdot \frac{m+1}{m+1} \quad (\text{at } h_j = 1).$$
     More precisely, since $C_j \le h_j + 1$, at $h_j = 1$ we have $C_j \le 2 \le m + 1$, and for all $h_j \ge 1$, $C_j \le \min(h_j + 1, h_j + m)$: wait—let's check the book's inequality $C_j \le h_j + m$ at $h_j = 1$!
     - At $h_j = 1$, $L_j$ has $|L_j| = h_j + 1 = 2$ leaves, so $C_j \le 2 \le m + 1$ (for $m \ge 1$).
     - In general, $C_j \le |L_j| = h_j + 1 \le h_j(m + 1)$ for all $h_j \ge 1, m \ge 1$!
     - Therefore:
       $$\frac{E_j}{2h_j} = 1 - \frac{\gamma C_j}{2h_j} \ge 1 - \frac{\gamma(h_j + 1)}{2h_j} \ge 1 - \gamma = 1 - \frac{2m - \alpha}{m + 1} \ge \frac{\alpha}{2m}$$
       (since $1 - \frac{2m-\alpha}{m+1} - \frac{\alpha}{2m} = \frac{(m-1)(2m-\alpha)}{2m(m+1)} \ge 0$ for all $m \ge 1, 0 < \alpha \le 2m$)!
   - Verified both by numerical quadrature (`scipy.integrate.quad` / `dblquad`) on a branching $n=5$ dendrogram and across all dendrograms for $n \le 6$ (`PASS step_3_4_nested_cone_integral`, `PASS step_3_4_exponent_bound`).
5. **Step 3.5 (Sum over dendrograms):**
   - Since the $(n-1)!$ edge orderings of a fixed labelled tree $\tau$ partition into the dendrograms $T \sim \tau$, and each $T \sim \tau$ accounts for $(n-1)! / \prod_j h_j$ orderings, dividing by $(n-1)!$ gives $\sum_{T \sim \tau} \prod_{j=1}^{n-1} \frac{1}{h_j} = 1$.
   - Verified exhaustively for $n \le 5$ (`PASS step_3_5_dendrogram_sum_identity`).

**Output of `referee_scratch/open_items/item6_step3_nn.py`:**
```text
PASS step_3_1_minimax_and_diam
PASS step_3_2_ancestor_bound
PASS step_3_3_linear_extensions
PASS step_3_5_dendrogram_sum_identity
PASS step_3_4_nested_cone_integral
PASS step_3_4_exponent_bound
```

### 6.2 Exact Replacement Text for `part9-foundations/48-continuum-window.qmd` (line 109)

In `part9-foundations/48-continuum-window.qmd`, replace line 109:
```markdown
(The integral identity of Step 3.4 was checked against numerical quadrature for $n=4$.)
```
with:
```markdown
(Steps 3.3–3.5, including the linear-extension count, the nested cone integral against numerical quadrature, the exponent bound, and the identity $\sum_{T\sim\tau}\prod_j(1/h_j)=1$, are checked in `code/review614/step3_combinatorics.py`.)
```

---

## Item 7: Missing Script Checks in Ch. 9 & Ch. 35 (PR #27 & PR #28)

- **Verdict:** **correct** (both mathematical claims are true to machine precision; drop-in check snippets for `code/algebra_data.py` and `code/qsp_checks.py` and text updates are provided below)

### 7.1 Ch. 9 (`part3-classical/09-lagrangian.qmd` line 116): $e_8^{(1)}$ Closed-Form Mass Ratios

In `referee_scratch/open_items/item7_missing_checks.py`, we verified that the 8 eigenvalues of $\mathsf M^2$ for $e_8^{(1)}$ in `code/algebra_data.py` (`masses('e', 8)`) match the closed-form expressions in `09-lagrangian.qmd` line 116:
$$m_2 = 2m_1\cos\frac{\pi}{5}, \quad m_3 = 2m_1\cos\frac{\pi}{30}, \quad m_4 = 2m_2\cos\frac{7\pi}{30}, \quad m_5 = 2m_2\cos\frac{2\pi}{15}, \quad m_6 = 2m_2\cos\frac{\pi}{30}, \quad m_7 = 4m_2\cos\frac{\pi}{5}\cos\frac{7\pi}{30}, \quad m_8 = 4m_2\cos\frac{\pi}{5}\cos\frac{2\pi}{15}$$
to max error $9.903 \times 10^{-13}$.

#### Exact Addition to `code/algebra_data.py` (append after line 52)
```python
_, _, ms_e8 = masses('e', 8)
m1_e8, c5_e8, c30_e8, c7_e8, c2_e8 = ms_e8[0], np.cos(np.pi/5), np.cos(np.pi/30), np.cos(7*np.pi/30), np.cos(2*np.pi/15)
m2_e8 = 2*m1_e8*c5_e8
report('e_8: closed-form mass ratios of sec-masses-couplings',
       np.allclose(ms_e8, [m1_e8, m2_e8, 2*m1_e8*c30_e8, 2*m2_e8*c7_e8, 2*m2_e8*c2_e8,
                           2*m2_e8*c30_e8, 4*m2_e8*c5_e8*c7_e8, 4*m2_e8*c5_e8*c2_e8]))
```

#### Exact Replacement Text in `part3-classical/09-lagrangian.qmd` (line 118)
Replace:
```markdown
`code/toda_classical.py` checks the Perron–Frobenius property and the sum rule for $a_2$, $a_5$, $d_4$, $d_6$, $e_6$, $e_7$ and $e_8$. `code/algebra_data.py` checks the $a_n$ and $d_n$ rows of the table for $n\le8$ and the folded masses below.
```
with:
```markdown
`code/toda_classical.py` checks the Perron–Frobenius property and the sum rule for $a_2$, $a_5$, $d_4$, $d_6$, $e_6$, $e_7$ and $e_8$. `code/algebra_data.py` checks all three rows of the table ($a_n$ and $d_n$ for $n\le8$, and $e_8$) and the folded masses below.
```

---

### 7.2 Ch. 35 (`part7-half-line/35-quantum-symmetric-pairs.qmd` line 206): Baseilhac–Koizumi Heisenberg Realization of (19) and $q$-Dolan–Grady

In `referee_scratch/open_items/item7_missing_checks.py`, we derived and verified the exact Weyl/Heisenberg realization of Baseilhac–Koizumi [@baseilhac_koizumi2003, eq. (19)]:
- Let $U = e^{\hat\theta}, V = e^{\hat\phi}$ be Weyl operators obeying $U V = q^2 V U$ (from $[\hat\theta, \hat\phi] = 2\ln q$), represented on the basis $\{|k\rangle\}_{k=-M}^M$ by $U|k\rangle = |k+1\rangle$ and $V|k\rangle = q^{-2k}|k\rangle$, and set
  $$\hat E_+ = \frac{ic}{q - q^{-1}}\big(\alpha U + \alpha^{-1}U^{-1}\big), \qquad \hat E_- = \frac{ic}{q - q^{-1}}\big(\bar\alpha V + \bar\alpha^{-1}V^{-1}\big).$$
- Because words of degree $\le 4$ in $\hat E_\pm$ shift $k$ by at most $\pm 4$, the action on interior basis states $|k\rangle$ ($|k| \le M - 4$) is **completely exact** (zero boundary truncation error).
- On these states, both Baseilhac–Koizumi's (19),
  $$[\hat E_+ \hat E_-, \hat E_- \hat E_+] = \rho(\hat E_+^2 - \hat E_-^2),$$
  and the $q$-Dolan–Grady relations,
  $$[\hat E_\pm, [\hat E_\pm, [\hat E_\pm, \hat E_\mp]_q]_{q^{-1}}] = \rho[\hat E_\pm, \hat E_\mp],$$
  hold to machine precision ($3.2 \times 10^{-16}$ at physical $q$, $1.4 \times 10^{-15}$ at generic $q$) with $\rho = c^2(q + q^{-1})^2$.

**Output of `referee_scratch/open_items/item7_missing_checks.py`:**
```text
e_8^(1) closed-form masses max error vs M^2 eigenvalues: 9.903e-13
PASS e8_closed_form_mass_ratios
BK Heisenberg (physical_q): err19=3.21e-16, err_dg_p=1.93e-16, err_dg_m=4.03e-16
PASS bk_heisenberg_physical_q
BK Heisenberg (generic_q): err19=8.53e-16, err_dg_p=1.43e-15, err_dg_m=2.65e-15
PASS bk_heisenberg_generic_q
```

#### Exact Addition to `code/qsp_checks.py` (append at the end of `code/qsp_checks.py`)
```python
# Baseilhac-Koizumi (2003) Heisenberg realization of eq. (19) and q-Dolan-Grady (sec-q-onsager)
worst_bk = 0.0
N_bk = 15; S_bk = slice(4, N_bk - 4)
U_bk = np.diag(np.ones(N_bk - 1, complex), k=-1); Uinv_bk = np.diag(np.ones(N_bk - 1, complex), k=1)
ks_bk = np.arange(N_bk) - N_bk // 2
for q in (QP, QG):
    V_bk = np.diag(q**(-2.0*ks_bk)); Vinv_bk = np.diag(q**(2.0*ks_bk))
    c_bk, al_bk, alb_bk = 0.9 - 0.3j, 1.3 - 0.4j, 0.8 + 0.7j
    pref = 1j*c_bk/(q - 1/q)
    Ep = pref*(al_bk*U_bk + Uinv_bk/al_bk); Em = pref*(alb_bk*V_bk + Vinv_bk/alb_bk)
    rho_bk = c_bk**2*(q + 1/q)**2
    worst_bk = max(worst_bk,
                   rel((Ep @ Em @ Em @ Ep - Em @ Ep @ Ep @ Em)[:, S_bk], (rho_bk*(Ep @ Ep - Em @ Em))[:, S_bk]),
                   rel(serre(Ep, Em, -2, q)[:, S_bk], (rho_bk*(Ep @ Em - Em @ Ep))[:, S_bk]),
                   rel(serre(Em, Ep, -2, q)[:, S_bk], (rho_bk*(Em @ Ep - Ep @ Em))[:, S_bk]))
report(f'Baseilhac-Koizumi Heisenberg realization satisfies (19) and q-Dolan-Grady with rho = c^2(q+q^-1)^2 '
       f'[{worst_bk:.1e}]', worst_bk < 1e-13)
```

#### Exact Replacement Text in `part7-half-line/35-quantum-symmetric-pairs.qmd` (line 206)
Replace:
```markdown
Their Heisenberg realization satisfies both, the q-Dolan–Grady relations with $\rho=c^2(q+q^{-1})^2$, which is (-@eq-q-dolan-grady) with $c_j=q^{-1}$ after their normalization $\hat Q=c\,b$, $\hat E=c\,\hat\epsilon$ (a numerical check, not included in `code`).
```
with:
```markdown
Their Heisenberg realization satisfies both, the q-Dolan–Grady relations with $\rho=c^2(q+q^{-1})^2$, which is (-@eq-q-dolan-grady) with $c_j=q^{-1}$ after their normalization $\hat Q=c\,b$, $\hat E=c\,\hat\epsilon$ (`code/qsp_checks.py`).
```

---

## Item 8: Preface (`preface.qmd`) Wording (PR #23)

- **Verdict:** **correct** (optional tightened wording below)

### 8.1 Analysis

We inspected `preface.qmd` lines 15 and 35 against the status labels in Parts VI, VII, VIII, and IX:
- **Line 15:** *"Every result in the book carries a status label ([Established], [Theorem], [Proposition], [Numerical], [Conjecture], [Open]), defined in the Introduction, so that the reader can see at a glance where the subject rests on proof, where on numerical verification, and where on conjecture."*
  - Note: `[Proposition]` is not one of the status labels in `status-tags.lua` / `01-introduction.qmd` (propositions carry `[Theorem]`, `[Established]`, `[Numerical]`, etc.). Removing `[Proposition]` from the parenthetical list in `preface.qmd` line 15 brings it into exact alignment with `01-introduction.qmd` (`sec-status-labels`) and `README.md`.
- **Line 35:** *"On the half-line (Part VII) and across a jump defect (Part VIII), the surviving combinations of bulk quantum-group generators form coideal subalgebras (or Borel subalgebras for a defect), and the reflection and transmission matrices are their intertwiners."*
  - This accurately describes the algebraic framework of Chapters 33–35 and 41–43 (where the intertwining equations are proved for spin-$1/2$ and vector representations and checked numerically for higher representations, while the full dynamical UV–IR parameter maps carry `[Conjecture]` labels).

### 8.2 Optional Tightened Wording for `preface.qmd` ( line 15)

In `preface.qmd`, replace line 15:
```markdown
Every result in the book carries a status label ([Established], [Theorem], [Proposition], [Numerical], [Conjecture], [Open]), defined in the Introduction, so that the reader can see at a glance where the subject rests on proof, where on numerical verification, and where on conjecture.
```
with:
```markdown
Every result in the book carries a status label ([Established], [Theorem], [Numerical], [Conjecture], [Open], [Refuted]), defined in the Introduction, so that the reader can see at a glance where the subject rests on proof, where on numerical verification, and where on conjecture.
```
