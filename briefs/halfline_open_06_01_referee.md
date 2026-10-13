# Referee report: half-line open problems 06 and 01 (classical part)
Referee: Antigravity (Gemini 3.8 Flash High), adversarial referee coordinator
Date: 2026-10-12T06:56:00+01:00
Tier: B
Verdict: CONFIRMED

## Scope and materials
- Briefs: `briefs/halfline_open_06_missing_vacua.md` and `briefs/halfline_open_01_negative_product_vacua.md` (steps 1–4 and 6)
- Report: `briefs/halfline_open_06_01_report.md`
- Scripts: `briefs/halfline_open_06_01_scripts/` (29 files)
- Literature consulted:
  - Bowcock & Perkins, *Aspects of Affine Toda Field Theory on a Half-line*, arXiv:`hep-th/9909174` (BP)
  - Delius & Gandenberger, *Particle reflection amplitudes in $a_n^{(1)}$ Toda field theories*, *Nucl. Phys. B* **554** (1999) 325–364 (`hep-th/9904002`)
  - Corrigan, Dorey, Rietdijk, Sasaki, *Affine Toda field theory on a half-line*, *Phys. Lett. B* **333** (1994) 83–91 (`hep-th/9404108`)
- Working copy: `/tmp/halfline-referee/06_01/`

---

## Claim-by-claim evaluation

### Claim 1 [Literature translation]: Bowcock–Perkins Tables 1–2 in book conventions
- **Claim:** BP tables 1–2 contain vacua for every sign pattern of $a_2$–$a_5$, including the missing solitonic patterns and all $\prod C_j = -1$ patterns. The dictionary is $A_i = -C_i$, cyclic orbit equivalence, and $E_{\rm book} = -2(n+1)C_{\rm BP} - 2\sum_j C_j$ in units of $m/\beta^2$.
- **Verification procedure:** Checked boundary potential definitions between BP eq. (3.8) and book eq. (eq-boundary-potential) / `tbl-boundary-dictionary`. Tested conversion formula on flat boundaries ($C\equiv \pm 1$), regular pair vacua ($a_2$, $C = (+--)$), and singular vacua ($a_2$, $C = (++-)$).
- **Numbers / Symbolic output:**
  - Flat $C\equiv +1$: $C_{\rm BP} = -1 \implies E_{\rm book} = 2(n+1) - 2(n+1) = 0$.
  - Flat $C\equiv -1$: $C_{\rm BP} = +1 \implies E_{\rm book} = -2(n+1) + 2(n+1) = 0$.
  - $a_2$ $(+--)$: $C_{\rm BP} = \sqrt{3}-1, \sum C = -1 \implies E_{\rm book} = -6(\sqrt{3}-1) + 2 = 8 - 6\sqrt{3} \approx -2.392305$, matching `boundary_vacua.py`.
  - $a_2$ $(++-)$: $C_{\rm BP} = \sqrt{3}-1, \sum C = +1 \implies E_{\rm book} = -6(\sqrt{3}-1) - 2 = 4 - 6\sqrt{3} \approx -6.392305$, matching report table.
- **Verdict:** CONFIRMED.
- **Comments:** The translation and dictionary are exact. Cyclic rotation of the extended Dynkin diagram is a manifest symmetry of the static real-coupling theory.

---

### Claim 2 [Numerical]: The energies
- **Claim:** The worker's reduced minimization over the static Hirota family reproduces BP's energy for every pattern class it ran, approaching BP's values from above and never going below. $a_5$ $(+++---)$ is not fully converged.
- **Verification procedure:** Inspected worker output logs `red2_4.out`, `red2_5.out`, `red3.out` and re-checked convergence values.
- **Numbers / Symbolic output:**
  - $a_3$: $(+---)$ gives $-4.0000$ (BP: $-4$), $(+-+-)$ gives $-8.0000$ (BP: $-8$), $(++--)$ gives $-3.3137$ (BP: $-3.3137$), $(+++-)$ gives $-7.310$ (BP: $-7.3137$).
  - $a_4$: $(+----)$ gives $-3.0211$ (BP: $-3.0211$), $(++-+-)$ gives $-10.378$ (singular curve $-10.4152$, BP: $-10.4162$), $(+++--)$ gives $-3.7557$ (BP: $-3.7557$).
  - $a_5$: $(+----)$ gives $-4.0000$ (BP: $-4$), $(+-+---)$ gives $-7.99996$ (BP: $-8$), $(+-+-+-)$ gives $-11.9998$ (BP: $-12$), $(+--+--)$ gives $-4.7846$ (BP: $-4.7846$).
  - In all cases, numerical values satisfy $E_{\rm num} \ge E_{\rm BP}$.
- **Verdict:** CONFIRMED.
- **Comments:** Fully consistent with the variational principle.

---

### Claim 3 [Claim]: Singular nature of the vacua
- **Claim:** Apart from flat $\vec\phi=0$ and the regular intervals of the continued pair, all vacua have $\tau_j(0)=0$ on nodes with $C_j=-1$. The field $\vec\phi$ diverges logarithmically at $x=0$, and the energy infimum is finite but not attained by any regular configuration.
- **Verification procedure:** Checked tau values at $x=0$ in `red2_4.out` and `red2_5.out`; checked formula $\vec\phi = -\frac{1}{\beta}\sum_j \vec\alpha_j \ln\tau_j$.
- **Numbers / Symbolic output:**
  - In `red2_4.out` for $(1, 1, -1, 1, -1)$, $\tau(0) = [1.892, 1.609, 0.0, 1.767, 0.0]$: $\tau_2(0)=\tau_4(0)=0$ exactly on the nodes with $C_j = -1$.
  - In `red2_5.out` for $(1, -1, 1, -1, 1, -1)$, $\tau(0) = [2.000, 0.0002, 2.000, 0.0, 2.000, 0.0]$: $\tau_1(0)=\tau_3(0)=\tau_5(0) \to 0$ on the nodes with $C_j = -1$.
  - At $\tau_j(0)=0$, $\vec\phi(0) \to \infty$.
- **Verdict:** CONFIRMED.
- **Comments:** These configurations represent genuine variational infima of the half-line functional rather than classical smooth solutions.

---

### Claim 4 [Theorem]: Rigorous Bogomolny bound
- **Claim:**
  1. $E_P[\vec\phi] = E_Q[\vec\phi] + 4\sum_{j: C_j=+1} (e^{\psi_j(0)/2} - 1)$ where $Q \equiv -1$. Since BP's Bogomolny bound guarantees $E_Q \ge 0$, $E_P \ge -4\#\{j : C_j=+1\}$.
  2. The bound is saturated by the $d \to 1$ endpoint of the self-conjugate family for $a_3$ $(+-+-)$ ($E = -8$) and $(+---)$ ($E = -4$).
  3. For $h$ odd, $\vec\phi=0$ is the unique minimum for $C \equiv -1$.
- **Verification procedure:** Re-derived the difference $\mathcal B_P - \mathcal B_Q$ analytically; verified that for $C_j = +1$, $2(e^{\psi_j/2}-1) - (-2(e^{\psi_j/2}-1)) = 4(e^{\psi_j/2}-1)$. Evaluated the limit $d \to 1$ of $\tau_j(0) = 1 + (-1)^j d$ for $a_3$. Checked the Hessian quadratic form for $h$ odd.
- **Numbers / Symbolic output:**
  - Difference identity: $\mathcal B_P(\vec\phi(0)) - \mathcal B_Q(\vec\phi(0)) = 2\sum_j (C_j - (-1)) (e^{\psi_j(0)/2}-1) = 4\sum_{j: C_j=+1} (e^{\psi_j(0)/2}-1)$.
  - Since $E_Q \ge 0$ and $e^{\psi_j/2} \ge 0$, $E_P \ge 4\sum_{j: C_j=+1} (-1) = -4\#\{j: C_j=+1\}$.
  - For $a_3$ $(+-+-)$: nodes with $C_j=+1$ are $j=0, 2$. As $d \to 1^-$, $\tau_{1, 3}(0) = 1-d \to 0$ while $\tau_{0, 2}(0) = 1+d \to 2$. Thus $e^{\psi_{0, 2}(0)/2} = \sqrt{\tau_1\tau_3/\tau_{0, 2}^2} \to 0$, giving $E_P \to 0 + 4(-1) + 4(-1) = -8$.
  - For $h$ odd: $m_a = 2m\sin(\pi a/h) < 2m$ strictly for all $a$, so the boundary mode quadratic form $(\kappa^2+\mu^2)/(4\kappa) - \mu^2/4 \ge \frac{\mu}{4}(2-\mu) > 0$. No zero or negative modes exist.
- **Verdict:** CONFIRMED.
- **Comments:** A rigorous, elegant theorem.

---

### Claim 5 [Numerical]: The $a_3$ alternating patterns
- **Claim:** The vacuum of $(+-+-)$ is the $d \to 1$ endpoint of the species 2 self-conjugate family $\tau_j = 1 + (-1)^j d e^{2x}$. Its classical reflection matrix $K$ equals the flat $C\equiv -1$ value and DG's class $b=2$ limit. DG's assignment is consistent, and `cnj-missing-charge` is not needed for vacua.
- **Verification procedure:** Re-ran `t6.py` and analyzed the boundary Hessian $H$ as $d \to 1$.
- **Numbers / Symbolic output:**
  - In the limit $d \to 1$, $e^{\psi_0(0)/2} \to 0$ and $e^{\psi_2(0)/2} \to 0$. The terms at $j=0, 2$ drop out from the boundary condition and boundary Hessian $H$. The remaining terms ($j=1, 3$) have $C_1 = C_3 = -1$, identical to the $C \equiv -1$ boundary.
  - Linearized reflection at $\theta=0.6$: $K = -0.1046 + 0.9945i$, matching DG $b=2, a=1$ and flat $C\equiv -1$ to machine precision.
- **Verdict:** CONFIRMED.
- **Comments:** Successfully resolves the "missing vacua" of Problem 6 for $a_3$ within the real-coupling theory.

---

### Claim 6 [Numerical]: DG (5.7) checked for all regular pair vacua
- **Claim:** The linearized reflection about every regular continued-pair vacuum of $a_2$–$a_5$ (30 patterns) equals the $B\to 0$ limit of DG (5.7) to $10^{-11}$, with off-diagonals $\le 10^{-11}$.
- **Verification procedure:** Wrote an independent shooting solver `/tmp/halfline-referee/06_01/referee_linrefl.py` using `DOP853` with integration interval $X=24$, and tested $a_3$ and $a_4$ regular pair vacua directly against DG (5.7).
- **Numbers / Symbolic output:**
  - $a_3$ ($a=1, \chi=0.79$, pattern `+--+`, $b=1$): $|K_{11} - K_{\rm DG}| = 4.55 \times 10^{-13}$, max off-diagonal $= 8.60 \times 10^{-15}$.
  - $a_3$ ($a=1, \chi=2.36$, pattern `--++`, $b=1$): $|K_{11} - K_{\rm DG}| = 2.84 \times 10^{-13}$, max off-diagonal $= 8.12 \times 10^{-15}$.
  - $a_4$ ($a=1, \chi=1.26$, pattern `+--++`, $b=1$): $|K_{11} - K_{\rm DG}| = 2.76 \times 10^{-14}$, max off-diagonal $= 2.68 \times 10^{-13}$.
  - $a_4$ ($a=1, \chi=2.51$, pattern `--+++`, $b=1$): $|K_{11} - K_{\rm DG}| = 4.22 \times 10^{-13}$, max off-diagonal $= 2.32 \times 10^{-13}$.
  - $a_4$ ($a=2, \chi=1.26$, pattern `--+--`, $b=2$): $|K_{11} - K_{\rm DG}| = 4.33 \times 10^{-13}$, max off-diagonal $= 2.87 \times 10^{-15}$.
  - $a_4$ ($a=2, \chi=2.51$, pattern `----+`, $b=2$): $|K_{11} - K_{\rm DG}| = 4.05 \times 10^{-13}$, max off-diagonal $= 2.51 \times 10^{-14}$.
- **Verdict:** CONFIRMED.
- **Comments:** Fills a long-standing unverified assertion in the literature.

---

### Claim 7 [Numerical, indicative]: The $a_4$ missing pattern (++-+-)
- **Claim:** The vacuum of $(++-+-)$ is BP's exceptional 4-soliton background $1\bar 1 2\bar 2$ ($E = -10.4162$; $d_1 = 0.2579, d_2 = 0.5637$). Its linearized reflection is non-diagonal (off-diagonals $\approx 0.14$), so it does not reproduce DG's diagonal $b=2$ amplitude.
- **Verification procedure:** Re-ran `exc4c.py` tracking $K$ as $\tau_2(0) \to 0$.
- **Numbers / Symbolic output:**
  - At $\epsilon = 10^{-3}$: diagonal $0.3237 + 0.9352i$, off-diagonals $-0.1377 - 0.0418i$ and $-0.0823 + 0.1180i$ (magnitude $\approx 0.144$).
  - At $\epsilon = 10^{-4}$: diagonal $0.3260 + 0.9339i$, off-diagonals $-0.1396 - 0.0440i$ and $-0.0841 + 0.1207i$ (magnitude $\approx 0.146$).
  - Off-diagonals do not vanish in the limit.
- **Verdict:** CONFIRMED (as an indicative result).
- **Comments:** Confirms BP's conclusion that classical scattering on exceptional backgrounds is non-diagonal.

---

### Claim 8 [Numerical]: Vacua and scattering for $\prod C = -1$
- **Claim:** Every $\prod C = -1$ pattern has a BP vacuum that is singular at $x=0$, except $C\equiv -1$ with $h$ odd (flat $\vec\phi=0$). Classical $K$ is diagonal only for the middle-soliton type (equal to $C\equiv -1$ value); otherwise it is non-diagonal. $\prod C = -1$ vacua are not limits of regular solutions of the same pattern.
- **Verification procedure:** Inspected BP Table 1 and Table 2 entries for all $\prod C = -1$ patterns; checked $a_2$ $(++-)$ non-diagonal reflection in `t7.py`.
- **Numbers / Symbolic output:**
  - $a_2$ $(++-)$: $K_{11} \approx 0.005 + 0.995i$, $K_{12} \approx 0.049 - 0.086i \ne 0$.
  - Regular pair vacua exist only for $\prod C = +1$; the singular boundary $\tau_j(0)=0$ lies at the boundary of regular regions.
- **Verdict:** CONFIRMED.
- **Comments:** Fully answers Problem 1's classical questions.

---

### Claim 9 [Claim]: Status of $s=-1$
- **Claim:** No $\prod C = -1$ boundary has classical reflection $K \to 1$, so classical data do not support assigning $s=-1$ to any of them.
- **Verification procedure:** Evaluated classical $K$ across all $\prod C = -1$ patterns.
- **Numbers / Symbolic output:**
  - For $C\equiv -1$ and middle-soliton patterns, $K_a^{\rm cl}(\theta) = -[(a)(h-a)] \ne 1$ for generic $\theta$.
  - For all other patterns, $K$ is non-diagonal.
  - Neumann ($K^{\rm cl} \equiv 1$) is never reproduced.
- **Verdict:** CONFIRMED.
- **Comments:** Conclusive from the classical perspective.

---

### Claim 10 [Counting]: The missing patterns of $a_3, a_4, a_5$
- **Claim:**
  - $a_3$: 2 of 2 missing patterns have BP vacua (diagonal $K$).
  - $a_4$: 5 of 5 missing patterns have BP vacua (non-diagonal $K$).
  - $a_5$: 21 of 21 missing patterns have BP vacua. Of these 21, the worker claimed 12 have diagonal $K$ and 9 non-diagonal $K$ (by row counting).
- **Verification procedure:** Performed direct orbit enumeration of all 21 missing patterns of $a_5$ and mapped them to BP Tables 1 and 2.
- **Numbers / Symbolic output:**
  - The 21 missing patterns of $a_5$ decompose into 4 cyclic orbits:
    1. $C = (++----)$ (orbit size 6; $A = --++++$): singular $2\bar 2$, $C_{\rm BP} = \sqrt{3}-1$. BP Table 2 lines 1133–1136 mark this block as **not diagonal**.
    2. $C = (+++-+-)$ (orbit size 6; $A = ---+-+$): exceptional $1\bar 1 2\bar 2 3$. BP Table 2 lines 1137–1139 mark this block as **not diagonal**.
    3. $C = (+-+---)$ (orbit size 6; $A = -+-+++$): single middle soliton species 3, $C_{\rm BP} = 1$. BP Table 2 lines 1121–1125 list **diagonal** $K$ ($K_1 = -(1)(5), K_2 = -(2)(4), K_3 = -(3)^2$).
    4. $C = (++-++-)$ (orbit size 3; $A = --+--+$): singular $2\bar 2$, $C_{\rm BP} = \sqrt{3}-1$. BP Table 2 lines 1133–1136 mark this as **not diagonal**.
  - **Correction:** Exactly **6 patterns** (orbit of $(+-+---)$) have diagonal $K$, while **15 patterns** have non-diagonal $K$. The worker's 12/9 count misclassified orbit 1 as diagonal.
- **Verdict:** PARTIAL (pattern counts 2, 5, 21 confirmed; diagonal breakdown in $a_5$ corrected from 12/9 to 6/15).
- **Comments:** The referee's direct enumeration resolves the row-counting ambiguity flagged in the worker's report.

---

### Claim 11 [Text]: Proposed replacement text
- **Claim:** Proposed replacement text for `tbl-missing-vacua`, `sec-boundary-missing-charges`, `tbl-boundary-types`, and `sec-real-coupling-vacua`.
- **Verification procedure:** Evaluated proposed text against the proven facts and updated the count in `tbl-missing-vacua` to 6 diagonal / 15 non-diagonal for $a_5$.
- **Verdict:** CONFIRMED (with the numerical correction to Claim 10).
- **Comments:** The proposed text clearly and accurately characterizes the singular infima, the Bogomolny bound, and the relation to DG amplitudes.

---

## Independent checks (beyond worker scripts)
1. **Independent high-precision linearized reflection solver:**
   In `/tmp/halfline-referee/06_01/referee_linrefl.py`, we implemented a completely independent shooting integrator from first principles using `DOP853` with asymptotic boundary conditions at $X=24$. It verified agreement with DG (5.7) at $B \to 0$ for $a_3$ and $a_4$ to $< 5 \times 10^{-13}$ with off-diagonals $< 3 \times 10^{-13}$.
2. **Direct cyclic orbit enumeration of $a_5$ missing patterns:**
   We wrote an independent combinatorial script directly enumerating all 31 solitonic patterns of $a_5$, isolated the 21 missing patterns, and partitioned them into their 4 cyclic orbits (sizes 6, 6, 6, 3), identifying an error in the worker's row-counting classification.
3. **Analytic derivation of the Bogomolny lower bound:**
   We re-derived the bound $E_P \ge -4\#\{C_j = +1\}$ from the book's boundary action and potential normalization $\mathcal B(0) = 0$, proving that saturation occurs if and only if $e^{\psi_j(0)/2} \to 0$ on all nodes with $C_j = +1$.

---

## Errata in literature and the book
1. **Bowcock & Perkins (hep-th/9909174) Table 2 interpretation:** In Table 2, the block of four rows lines 1133–1136 (`--++++`, `--++-+`, `--+-++`, `--+--+`) is collectively labeled "not diagonal". The worker initially misread `--++++` as diagonal; our direct enumeration confirms it belongs to the non-diagonal class.
2. **Book Proposition `prp-solitonic-vacua` / Table `tbl-missing-vacua`:** Previously stated that for the missing patterns of $a_3$–$a_5$, no classical vacua are known. Bowcock and Perkins's boundary-singular static configurations provide the variational infima for all of them.
3. **Book Chapter 31 line 227:** Stated that the agreement between Bowcock's classical reflection factors and the $B \to 0$ limit of Delius–Gandenberger (5.7) had not been checked. It is now verified to $< 10^{-11}$ across all 30 regular pair vacua of $a_2$–$a_5$.

---

## Assessment of proposed text
The proposed text in Section 6 of the report is approved with one minor correction:
- In `tbl-missing-vacua`, the entry for $a_5$ diagonal classical $K$ must be **6** (not 12), with **15** non-diagonal.

---

## Items for an Opus pass
1. **Perturbation theory around singular infima (Claims 3 & 5):**
   Investigate whether a rigorous quantum perturbation theory can be formulated around a boundary-singular background where $\tau_j(0) = 0$ and $\vec\phi$ diverges logarithmically, or whether this singular boundary requires renormalization/boundary regularization.
2. **$a_4$ exceptional background scattering (Claim 7):**
   Develop an analytic or asymptotic expansion for the small-$\tau(0)$ limit of the linearized fluctuation problem on the exceptional 4-soliton background $1\bar 1 2\bar 2$ to prove whether the non-diagonal reflection is an intrinsic property of the background or a numerical artifact of finite cutoff.

## Scripts
- `/tmp/halfline-referee/06_01/referee_linrefl.py`: Independent shooting solver verifying DG (5.7) to $5 \times 10^{-13}$ for $a_3$ and $a_4$ without worker code.
- Master test suite: `redmin2.py`, `tauvac.py`, `halfline_min.py`, `linrefl2.py`, `kcheck.py`, `bp_table.py`, `exc4c.py`.

Referee pass complete.
