# Adversarial Referee Report: Part VIII (Foundations of the Non-Unitary Theory)

**Target:** Part VIII (`part8-foundations/36-twelve-gaps.qmd` through `41-unitarity.qmd`), Chapter 42 (`outlook/42-outlook.qmd`), and accompanying codebase (`foundations_code/`, `lightcone_krein_check.py`).  
**Referee:** Mathematical-physics adversarial referee pass.  
**Date:** 2026-10-11  

---

## 1. Executive Summary & Headline Verdict

Part VIII investigates whether imaginary-coupling affine Toda field theories (ATFTs) exist as valid non-unitary quantum field theories. The author formulates twelve historical gaps (H1–H12) and addresses them through lattice bounds, Coulomb-gas cluster expansions, Krein-space transfer matrices, and quantum affine scattering.

The manuscript presents substantial constructive mathematics—notably the multi-scale cluster expansions in Chapter 39 extending the continuum limit into the window $2\pi \le \beta^2 < \beta^2_{\rm UV}$, and the exact Krein-unitarity of the fundamental soliton transfer matrix in Chapter 40. 

However, an adversarial audit reveals several **overclaims, unproved assertions, mathematical gaps, and an empirical error in the author's numerical claims:**

1. **Empirical Falsification in Chapter 40 (§40.2, line 118):** The author asserts that for fundamental solitons of $a_2^{(1)}$ ($n=3$), $\Theta$-breaking in sectors containing all three colours occurs "at every coupling sampled. The same holds for $n=3$ with $N=4$". Running the author's committed script [`foundations_code/scan_n4.py`](file:///home/gustav/Git/affine_toda/foundations_code/scan_n4.py) and extensive checks shows that at $\mu=2.5$, the spectrum for $n=3, N=4$ is **completely unbroken** across all sectors (deviations $\le 5.98\times 10^{-12}$). The author's own code outputted `9.24e-13` at $\mu=2.5$, yet the text declared universal breaking.
2. **Vacuum Non-Degeneracy Gap ([Theorem 37.4 (`thm-vacuum-energy`)](file:///home/gustav/Git/affine_toda/part8-foundations/37-lattice.qmd#L128-L150)):** Pringsheim's theorem proves that the spectral radius $\rho(T)$ is an eigenvalue and $\operatorname{Re} E \ge E_0$. It **does not prove** non-degeneracy, algebraic simplicity, or strict minimality ($\operatorname{Re} E > E_0$ for all $\lambda \ne \rho(T)$). Pringsheim's theorem is mathematically compatible with complex-conjugate eigenvalue pairs $\lambda_\pm = \rho(T)e^{\pm i\theta}$ ($\theta\ne 0$) on the spectral circle. In Perron–Frobenius theory, strict dominance requires kernel positivity, but the transfer matrix $T$ is not positivity-preserving on $L^2$ due to phase oscillations of $e^{-\frac{a^2}{2}V(\vec\phi)}$.
3. **Overclaim on the Sign Problem (H1, H2, [Theorem 37.2](file:///home/gustav/Git/affine_toda/part8-foundations/37-lattice.qmd#L54-L86)):** The claim that the complex action "never produces a sign problem" holds *exclusively* for the grand-canonical Coulomb gas representation of vertex operators with real charges. For polynomials of $\vec\phi$, correlation functions are sign-indefinite. In the transfer-matrix formulation on the light-cone lattice ([`lightcone_krein_check.py`](file:///home/gustav/Git/affine_toda/lightcone_krein_check.py)), the average phase $\operatorname{Tr} V^m / \operatorname{Tr} W^m$ decays to $0.18$ by $m=8$.
4. **Continuum Limit Inhomogeneity (H3, [Theorem 38.2](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L141-L150), [Theorem 39.7](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L299-L308)):** The theory constructed in finite volume is *not* standard ATFT with constant coupling. [Proposition 38.3(2)](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L115-L135) proves that constant-coupling Wick ordering diverges at the boundary for all $\beta^2 \ge 2\pi$. The finite-volume theory exists solely because pointwise Wick ordering forces the coupling to vanish at the Dirichlet walls ($g_j(u) \sim (2d)^{x_j} \to 0$).
5. **Missing Non-Simply-Laced Algebras in $\beta^2_{\rm UV}$ List ([Ch. 38, line 30](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L30)):** The list of algebras with $\beta^2_{\rm UV} = 4\pi$ omits $f_4^{(1)}$ and all twisted affine algebras ($a_{2n-1}^{(2)}, d_{n+1}^{(2)}, e_6^{(2)}, d_4^{(3)}$).
6. **Asymptotic Base Eigenvalue Slip ([Ch. 40, line 185](file:///home/gustav/Git/affine_toda/part8-foundations/40-scattering.qmd#L185)):** The large-$x$ splitting formula for the two-kink $a_2^{(2)}$ amplitude approaches $q^{-4}$, not $q^4$, as $x \to \infty$.
7. **Crossed-Channel Krein Unitarity Unproved (H9, [Ch. 40, line 120](file:///home/gustav/Git/affine_toda/part8-foundations/40-scattering.qmd#L120)):** Exact Krein-unitarity is proved for fundamental solitons of $a_n^{(1)}$. For crossed (soliton–antisoliton) amplitudes, the required twisted Krein metric is unconstructed.
8. **Verification of Step 3 in [Theorem 39.2](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L41-L107) and [Theorem 39.4](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L213-L295):** The referee independently re-derived and certified the EMST dendrogram induction, nested scale integration, and hook-length summation.

---

## 2. Detailed Findings

### Finding 1: Empirical Falsification in Chapter 40 (§40.2, line 118)
- **Manuscript Text:**
  > "For $n=3$ and $N=3$, with random real rapidities, $\max| |s|-1 |$ in the sector containing all three colours ranges from $0.07$ to $0.47$ for $0.3\le\mu\le 3$, so the breaking occurs at every coupling sampled. The same holds for $n=3$ with $N=4$, and for $n=4$ with $N=3$."
- **Referee Recheck:**
  Running [`foundations_code/scan_n4.py`](file:///home/gustav/Git/affine_toda/foundations_code/scan_n4.py) outputs:
  ```
  n=3 N=4 mu=0.5: max||s|-1| by number of distinct colours: {1: 2.22e-16, 2: 2.34e-11, 3: 0.166}
  n=3 N=4 mu=1.5: max||s|-1| by number of distinct colours: {1: 1.11e-16, 2: 1.31e-13, 3: 0.362}
  n=3 N=4 mu=2.5: max||s|-1| by number of distinct colours: {1: 2.22e-16, 2: 1.19e-12, 3: 9.24e-13}
  ```
  Notice that at $\mu=2.5$, the maximum deviation for 3 distinct colours is $9.24 \times 10^{-13}$.
  A targeted check over 200 random rapidity sets at $\mu=2.5$ confirmed that across every topological sector of $n=3, N=4$ (including the 3-colour sectors $(0,0,1,2)$, $(0,1,1,2)$, $(0,1,2,2)$), the maximum deviation from unimodularity is at most $5.98 \times 10^{-12}$.
  The spectrum at $\mu=2.5$ is **unbroken**.
- **Assessment:**
  The blanket claim "The same holds for $n=3$ with $N=4$" is false. The author had an incentive to claim that $\Theta$-breaking occurs for all couplings whenever $\ge 3$ colours are present, but their own code outputted unimodularity at $\mu=2.5$.
- **Required Fix:** Replace line 118 with:
  > "For $n=3$ with $N=4$, breaking occurs across most couplings sampled (e.g. $\mu=0.5, 1.5$), but at $\mu=2.5$ the eigenvalues remain unimodular to $10^{-12}$, indicating an unbroken window at large $\mu$."

---

### Finding 2: Vacuum Energy Non-Degeneracy Gap ([Theorem 37.4](file:///home/gustav/Git/affine_toda/part8-foundations/37-lattice.qmd#L128-L150))
- **Manuscript Text:**
  > "[Theorem] The spectral radius $\rho(T)$ is positive and is an eigenvalue of $T$. Consequently the lattice vacuum energy $E_0=-a^{-1}\ln\rho(T)$ is real, and every energy $E=-a^{-1}\ln\lambda$ with $\lambda\in\sigma(T)\setminus\{0\}$ satisfies $\operatorname{Re}E\ge E_0$."
  In Chapter 36: "the vacuum is the state of lowest real energy."
- **Referee Analysis:**
  The proof uses Pringsheim's theorem on $f(w) = \sum_{N\ge 1}\operatorname{Tr}(T^N) w^N = \sum_\lambda m_\lambda \frac{\lambda w}{1 - \lambda w}$.
  1. Because $\operatorname{Tr} T^N = Z_\Lambda(N) > 0$, Pringsheim's theorem guarantees that the positive real point on the circle of convergence $w = 1/\rho(T)$ is a singularity.
  2. Because $T$ is trace class, the spectrum outside 0 consists of isolated eigenvalues of finite multiplicity, and poles cannot cancel. Hence $\rho(T) \in \sigma(T)$.
  3. This proves that $E_0 = -a^{-1}\ln\rho(T)$ is real, and for any $\lambda \in \sigma(T)$, $|\lambda| \le \rho(T) \implies \operatorname{Re} E = -a^{-1}\ln|\lambda| \ge -a^{-1}\ln\rho(T) = E_0$.
- **The Gap:**
  Pringsheim's theorem does **not** assert that $w = 1/\rho(T)$ is the *only* singularity on the circle $|w| = 1/\rho(T)$:
  - Any power series $\sum a_N w^N$ with $a_N > 0$ can have multiple singularities on its circle of convergence (e.g. $g(w) = 2/(1-w) + 1/(1+w^2) = \sum_{N\ge 0} (2 + \cos(N\pi/2))w^N$ has poles at $w = 1$ and $w = \pm i$).
  - If $T$ had an eigenvalue pair $\lambda_\pm = \rho(T) e^{\pm i\theta}$ ($\theta \ne 0$), their real energy would satisfy $\operatorname{Re}(-a^{-1}\ln\lambda_\pm) = E_0$.
  - In Perron–Frobenius theory, strict dominance ($|\lambda| < \rho(T)$) requires the transfer kernel to be strictly positive. But $T$ is not positivity-preserving on $L^2$ because $V(\vec\phi)$ is complex.
  - Thus, Pringsheim's theorem proves non-strict minimality $\operatorname{Re} E \ge E_0$, but does not prove algebraic simplicity ($m_{\rho(T)} = 1$) or strict minimality ($\operatorname{Re} E > E_0$ for $\lambda \ne \rho(T)$).
- **Required Fix:** Keep **[Theorem]**, but add a clarifying remark following line 150 documenting that non-degeneracy and strict isolation are open for the field-theoretic transfer matrix, though verified in the zero-mode reduction.

---

### Finding 3: Overclaimed Scope of "No Sign Problem" (H1, H2, [Theorem 37.2](file:///home/gustav/Git/affine_toda/part8-foundations/37-lattice.qmd#L54-L86))
- **Manuscript Text:**
  > "The 'complex action' therefore never produces a sign problem." (Abstract §36.3, line 56; §37.3, line 84).
- **Referee Analysis:**
  The equivalence of $Z_\Lambda$ to a positive Coulomb gas is exact: since $C_a \otimes \mathbb{I} > 0$, $\exp(-\frac{1}{2}\langle J, C_a J\rangle) \in (0,1]$ for real charges.
  However:
  1. The elementary field $\vec\phi(x)$ is not a positive observable. Differentiating with respect to the source $\vec w$ at $\vec w=0$ inserts factors of $i$, so multi-point correlation functions of $\vec\phi$ are sign-indefinite.
  2. In the transfer-matrix formulation on the light-cone lattice ([`lightcone_krein_check.py`](file:///home/gustav/Git/affine_toda/lightcone_krein_check.py)), the average phase $\operatorname{Tr} V^m / \operatorname{Tr} W^m$ decays exponentially with time slice number $m$ ($0.78, 0.66, 0.49, 0.18$ for $m=1,2,4,8$).
  3. Real-time evolution exhibits exponential growth in $\Theta$-broken sectors.
- **Required Fix:** Qualify the claim: the absence of a sign problem applies strictly to the dual Coulomb gas ensemble of vertex operators with real charges, not to the elementary field algebra or real-time dynamics.

---

### Finding 4: Inhomogeneous Nature of the Continuum Limit (H3, [Theorem 38.2](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L141-L150), [Theorem 39.7](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L299-L308))
- **Manuscript Text:**
  Table row H3 states: "proved here in finite volume for $\beta^2 < \beta^2_{\rm UV}$".
- **Referee Analysis:**
  [Proposition 38.3(2)](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L115-L135) proves that site-independent Wick ordering (constant bulk coupling) diverges at the Dirichlet boundary as $a^{1 - \beta^2/2\pi}$ for all $\beta^2 \ge 2\pi$.
  Convergence in the window $2\pi \le \beta^2 < \beta^2_{\rm UV}$ in [Theorem 39.7](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L299-L308) holds only because pointwise Wick ordering imposes a position-dependent coupling profile $g_j(u) = z R_\Omega(u)^{x_j} \sim (2d)^{x_j}$ that vanishes at the walls.
- **Required Fix:** In Chapter 36 (abstract and row H3), state clearly that the finite-volume continuum limit is proved for the theory with this Dirichlet boundary profile.

---

### Finding 5: Omission of Non-Simply-Laced Algebras in $\beta^2_{\rm UV}$ List ([Ch. 38, line 30](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L30))
- **Manuscript Text:**
  > "- For $b_n^{(1)}$, $c_n^{(1)}$, $g_2^{(1)}$ and $a_2^{(2)}$ the minimum is $4\pi$, attained by a pair of adjacent roots with $\vec\alpha_i\cdot\vec\alpha_j=-1$."
- **Referee Analysis:**
  $f_4^{(1)}$ is omitted. In $f_4^{(1)}$, adjacent long and short roots have $\vec\alpha_2\cdot\vec\alpha_3 = -1 \implies 4\pi$. The neutral cluster has labels $(1,2,3,4,2)$, giving threshold $8\pi(11)/18 = 44\pi/9 \approx 4.89\pi > 4\pi$.
  The twisted algebras $a_{2n-1}^{(2)}, d_{n+1}^{(2)}, e_6^{(2)}, d_4^{(3)}$ are also omitted and all share the threshold $4\pi$.
- **Required Fix:** Update line 30 to include $f_4^{(1)}$ and the twisted families.

---

### Finding 6: Slip in Large-$x$ Splitting Formula ([Ch. 40, line 185](file:///home/gustav/Git/affine_toda/part8-foundations/40-scattering.qmd#L185))
- **Manuscript Text:**
  > "At large $x$ two roots approach $q^4$ as $\lambda=q^4(1+\delta)$, with $\delta^2=-16\sin^2\psi\cos\psi/x+o(1/x)$, which proves the breaking for $\cos\psi<0$ [Theorem]."
- **Referee Verification:**
  We computed the large-$x$ eigenvalues of the two-kink $a_2^{(2)}$ amplitude numerically and asymptotically at $\xi = 0.9\pi$:
  - At $x = e^{10} \approx 22026.5$: the roots approach $\lambda \to q^{-4}$, not $q^4$.
  - The relative shifts from $q^{-4}$ are $d_1 \approx +0.01637$ and $d_2 \approx -0.01610$, yielding $d_1 d_2 \approx -2.636 \times 10^{-4}$.
  - The predicted splitting $-16\sin^2\psi\cos\psi/x$ equals $+2.637 \times 10^{-4}$.
  - The formula $\delta^2 = -16\sin^2\psi\cos\psi/x$ is accurate to 4 decimal places, but the base eigenvalue approached as $x \to \infty$ is $q^{-4}$.
- **Required Fix:** Replace $q^4$ with $q^{-4}$ in line 185.

---

### Finding 7: Independent Certification of Step 3 in Theorems 39.2 & 39.4
The author noted in Chapter 39 (lines 109 and 311) that Step 3 of [Theorem 39.2](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L41-L107) and [Theorem 39.4](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L213-L295) had only been tested by AI-assisted checks and requested an independent referee reading.
- **Referee Derivation:**
  1. *EMST tree-path property:* For any points $u_k, u_l$, $|u_k - u_l| \ge \max_{e\in\operatorname{path}_\tau(u_k, u_l)} \operatorname{length}(e)$. Since any point outside the cluster of $u_k$ just before merger $i(k)$ is separated by an edge of length at least $s_{i(k)}$, and the cluster has at most $m$ points, the $m$-th nearest neighbour satisfies $\tilde\delta_k^{(m)} \ge s_{i(k)}$.
  2. *Nested Scale Integration:*
     $$\int_{R(T)} \prod_{j=1}^{n-1} s_j^{1 - \gamma c_j} ds = \prod_{j=1}^{n-1} \frac{1}{E_j},\qquad E_j = 2h_j - \gamma C_j.$$
     Verified by induction from leaves to root.
  3. *Hook-Length Summation:* $\sum_{T\sim\tau} \prod_j 1/h_j = 1$ correctly sums the dendrogram orderings to 1 for every spanning tree $\tau$.
  4. *Exponents in Theorem 39.4:* Checked all three cases:
     - Case (i) (active nodes inside unit): $E_j = 2h_j - p\kappa D(S_j) \ge \eta > 0$.
     - Case (ii) (inactive sub-unit nodes): $C_j = 0 \implies E_j = 2h_j > 0$.
     - Case (iii) (superclusters with $\ge m+1$ points): pair exponents cancel, leaving only $m$-th neighbour exponents $C_j = p\sum \gamma_k$, where $E_j \ge 2m - (m+1)p\bar\gamma > 0$.
  5. The multi-external point cluster bound in Step 3/4 correctly uses the macroscopic separation $d_v$ to bound chains of points connecting different external charges.
- **Verdict:** **Keep [Theorem]**. Step 3 is certified as mathematically sound.

---

### Finding 8: Crossed-Channel Krein Unitarity Unproved (H9)
- [Proposition 40.3 (`prp-krein-fundamental`)](file:///home/gustav/Git/affine_toda/part8-foundations/40-scattering.qmd#L91-L115) rigorously proves Krein-unitarity for fundamental solitons of $a_n^{(1)}$ where the Krein metric is $\mathcal C = C^{\otimes N}$.
- For soliton–antisoliton scattering, crossing introduces $q$-dependent diagonal twists, and the corresponding Krein metric is not proven (§40.2, line 120).
- **Required Fix:** Row H9 in Chapter 36 should state: "closed for fundamental solitons of $a_n^{(1)}$; conjectural for crossed channels".

---

## 3. Formal Verdict Table for Part VIII Theorems

| ID | Location | Claimed Status | Referee Verdict | Key Requirement |
|---|---|---|---|---|
| `thm-bounded-weight` | [Ch. 37, L34](file:///home/gustav/Git/affine_toda/part8-foundations/37-lattice.qmd#L34) | Theorem | **Keep [Theorem]** | Verified; $\operatorname{Re} V \ge 0 \implies \|e^{-a^2\sum V}\| \le 1$. |
| `thm-coulomb-gas` | [Ch. 37, L54](file:///home/gustav/Git/affine_toda/part8-foundations/37-lattice.qmd#L54) | Theorem | **Keep [Theorem]** | Positive Boltzmann weights verified. Qualify sign-problem claims. |
| `thm-krein-transfer` | [Ch. 37, L100](file:///home/gustav/Git/affine_toda/part8-foundations/37-lattice.qmd#L100) | Theorem | **Keep [Theorem]** | Verified; $T^\dagger = PTP$, $\Theta = PK$ commutes with $T$. |
| `thm-vacuum-energy` | [Ch. 37, L128](file:///home/gustav/Git/affine_toda/part8-foundations/37-lattice.qmd#L128) | Theorem | **Keep [Theorem]** | Proves $\rho(T) \in \sigma(T)$ and $\operatorname{Re} E \ge E_0$. Add remark on non-degeneracy gap. |
| `prp-continuum-limit` | [Ch. 38, L82](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L82) | Theorem | **Keep [Theorem]** | Established via Theorems 38.2 and 39.7 with the boundary profile. |
| `prp-necessity` | [Ch. 38, L115](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L115) | Theorem | **Keep [Theorem]** | Boundary divergence for constant coupling at $\beta^2 \ge 2\pi$ verified. |
| `thm-below-2pi` | [Ch. 38, L141](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L141) | Theorem | **Keep [Theorem]** | Nearest-neighbour integral bound fully verified. |
| `lem-newton` | [Ch. 38, L152](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L152) | Theorem | **Keep [Theorem]** | Discrete harmonic exit law and screening verified. |
| `lem-multiple-sites` | [Ch. 38, L170](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L170) | Theorem | **Keep [Theorem]** | Verified. |
| `lem-nn-integral` | [Ch. 38, L201](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L201) | Theorem | **Keep [Theorem]** | Geometric angle bound $> 60^\circ \implies \|A\| \le 5$ verified. |
| `prp-uv-divergence` | [Ch. 38, L235](file:///home/gustav/Git/affine_toda/part8-foundations/38-continuum-thresholds.qmd#L235) | Theorem | **Keep [Theorem]** | Uniform cluster collapse scaling $s^{2n-3-\kappa D(S)}$ verified. |
| `thm-nn-distances` | [Ch. 39, L41](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L41) | Theorem | **Keep [Theorem]** | Step 3 EMST dendrogram induction independently certified. |
| `lem-overlap` | [Ch. 39, L117](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L117) | Theorem | **Keep [Theorem]** | Markov property identity verified. |
| `lem-unit-majorant` | [Ch. 39, L142](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L142) | Theorem | **Keep [Theorem]** | Attractive/repulsive signed pairs bounds verified. |
| `lem-tightness` | [Ch. 39, L166](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L166) | Theorem | **Keep [Theorem]** | Combinatorial tightness bound $\hat R_B \ge c_m \delta_k^{(m)}$ verified. |
| `lem-nested-scales` | [Ch. 39, L194](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L194) | Theorem | **Keep [Theorem]** | Tested against 708 and 13,326 random trees; verified. |
| `thm-cluster-bound` | [Ch. 39, L213](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L213) | Theorem | **Keep [Theorem]** | Multi-scale cluster induction verified. |
| `thm-window` | [Ch. 39, L299](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L299) | Theorem | **Keep [Theorem]** | Vitali/Tannery convergence in the window verified. |
| `lem-stability` | [Ch. 39, L313](file:///home/gustav/Git/affine_toda/part8-foundations/39-continuum-window.qmd#L313) | Theorem | **Keep [Theorem]** | Hölder/Stirling bounds verified. |
| `prp-krein-unitarity`| [Ch. 40, L35](file:///home/gustav/Git/affine_toda/part8-foundations/40-scattering.qmd#L35) | Prop, sketch | **Keep [Prop, sketch]**| Asymptotic completeness in Krein space remains open. |
| `prp-krein-regularization`| [Ch. 40, L60](file:///home/gustav/Git/affine_toda/part8-foundations/40-scattering.qmd#L60) | Theorem | **Keep [Theorem]** | Light-cone circuit verified on $L=6,8$ to machine precision. |
| `prp-krein-fundamental`| [Ch. 40, L91](file:///home/gustav/Git/affine_toda/part8-foundations/40-scattering.qmd#L91) | Theorem | **Keep [Theorem]** | Parts 1–3 rigorously verified. |
| `prp-scalar-breathers` | [Ch. 40, L164](file:///home/gustav/Git/affine_toda/part8-foundations/40-scattering.qmd#L164) | Theorem | **Keep [Theorem]** | Moduli 1 for real shifts verified. |
| `prp-exchange-twist` | [Ch. 40, L191](file:///home/gustav/Git/affine_toda/part8-foundations/40-scattering.qmd#L191) | Theorem | **Keep [Theorem]** | Telescoping identity verified. |
| `prp-bloch-sectors` | [Ch. 40, L216](file:///home/gustav/Git/affine_toda/part8-foundations/40-scattering.qmd#L216) | Prop, sketch | **Keep [Prop, sketch]**| Keep sketch label; properly restricted to bosonic forms. |
| `lem-no-go` | [Ch. 41, L24](file:///home/gustav/Git/affine_toda/part8-foundations/41-unitarity.qmd#L24) | Theorem | **Keep [Theorem]** | Rigorous linear algebra. |

---

## 4. Prior Art and Originality

1. **Multi-Component Coulomb Gas & Thresholds:**
   The two-component Coulomb gas for sine-Gordon ($a_1^{(1)}$) is classical (Fröhlich 1976; Benfatto–Gallavotti–Nicolò 1982). Extending this to general Lie algebras and showing that the non-Hermitian action becomes a positive Coulomb gas is genuine new work.
2. **The Continuum Window ($2\pi \le \beta^2 < \beta^2_{\rm UV}$):**
   The cluster expansion technique preserving exact short-range attractions and repulsions within tight units, combined with the EMST dendrogram bounds, is an original constructive QFT contribution.
3. **Krein-Unitarity of Soliton S-Matrices:**
   Takács and Watts (1999) observed complex-conjugate eigenvalue pairs empirically in finite volume. Identifying the exact algebraic symmetry $\mathcal C$ intertwining the transfer matrix, proving that two-colour subsectors are unitarily similar, and formulating the Krein-unitarity framework are original and settle the structural origin of their observations.

---

## 5. Integration (2026-10-13, branch foundations-referee-fixes)

The report predates the defects part. Its Part VIII is now Part IX, and its chapters 36–42 are now 44–50 (`part9-foundations/44-twelve-gaps.qmd` … `49-unitarity.qmd`, `outlook/50-outlook.qmd`); the line numbers it quotes have moved. Every finding was re-derived or recomputed before any change; the verdicts follow. Scripts of this pass: `foundations_code/scan_n4.py` (rewritten), `foundations_code/uv_thresholds.py` (new), `foundations_code/zero_mode.py` (one check added), and `briefs/referee_part8_scripts/n4_window_independent.py`.

**Finding 1 (no breaking for $n=3$, $N=4$ at $\mu=2.5$): confirmed, and sharpened.** The old `scan_n4.py` does print `9.24e-13` for the three-colour sectors at $\mu=2.5$, so "the same holds for $n=3$ with $N=4$" was false. A scan in $\mu$ shows a sharp window: the three-colour sectors of four solitons are broken for $\mu<2\pi/3$ and unbroken for $2\pi/3<\mu<\pi$, i.e. $\operatorname{Re}q<-\frac12$ (deviation $1.2\times10^{-2}$ at $\mu=2.09$, $1.9\times10^{-8}$ at $2\pi/3$, $4\times10^{-11}$ at $2.0954$). An independent implementation (tensordot, all four $T_j$, which commute to $10^{-14}$) gives the same. More solitons show a further structure. At $q=-1$, $\tilde R$ is $+1$ on equal colours and $-1$ on different ones, so $T_1$ is diagonal with eigenvalue $(-1)^{N-m_c}$. In every case computed ($n=3$ with $N\le7$, $n=4$ with $N\le6$), a sector is unbroken near $q=-1$ exactly when at most two colours have even multiplicity and at most two odd: $(2,1,1)$, $(2,2,1)$, $(3,2,1)$, $(4,1,1)$, $(3,2,2)$, $(4,2,1)$, $(2,2,1,1)$ unbroken, while $(1,1,1)$, $(3,1,1)$, $(2,2,2)$, $(3,3,1)$, $(5,1,1)$, $(1,1,1,1)$, $(2,1,1,1)$, $(3,1,1,1)$ are broken. The pattern $(2,2,1,1)$ was predicted by the rule before it was computed. The onset depends on the pattern: between $\mu=2.0$ and $2.1$ for $(2,2,1)$, between $2.45$ and $2.55$ for $(4,1,1)$. Also, "$n=4$ with $N=3$" is not independent evidence: on three colours the $U_q(\widehat{sl}_n)$ R-matrix is the $U_q(\widehat{sl}_3)$ one. The N=3 range was recomputed as $0.08$–$0.49$ (the old $0.07$–$0.47$ came from a different sample). The referee's replacement text was not used. Changes: the *Numerics* after `prp-krein-fundamental` (ch. 48), the sector-status table (one row split in two), the unbroken-sector list of ch. 49 (U-table), outlook problem 5; `scan_n4.py` now checks every claim.

**Finding 2 (vacuum non-degeneracy): confirmed as a limit of the statement, not an error.** `thm-vacuum-energy` claims only $\rho(T)\in\sigma(T)$ and $\operatorname{Re}E\ge E_0$, and its proof is correct; but neither simplicity of $\rho(T)$ nor the absence of other eigenvalues of modulus $\rho(T)$ follows, and the kernel of $T$ is not positivity preserving, so Jentzsch/Perron–Frobenius does not apply. The text of ch. 44 ("has the smallest real part") read as if strict. Positivity of the traces does give something: with $a_N=\rho^{-N}\operatorname{Tr}T^N>0$ and the non-peripheral part $o(1)$, the Bohr mean of $a_Ne^{-iN\vartheta}$ is $m_{\rho e^{i\vartheta}}$ and is bounded by the mean of $|a_N|=a_N$, which is $m_\rho$; so $m_\lambda\le m_{\rho(T)}$ on the spectral circle. Zero-mode reduction: the vacuum is simple and the next real part is higher by $1.98$, $2.29$, $1.35$, $1.92$ ($a_2$, $g=0.5,1$; $c_2$, $g=0.5,1$). Changes: a paragraph after the theorem and `exr-peripheral` (ch. 45), a check in `zero_mode.py`, the abstract and row H5 (ch. 44, "non-degeneracy open").

**Finding 3 (scope of "no sign problem"): partly confirmed.** Point 1 is right that correlation functions of $\vec\phi$ are not positive: with charges present the Gaussian integral shifts $\vec\phi$ by $i\sum_kq_{j_k}C(\cdot,u_k)$, and $\overline{V(\vec\phi)}=V(-\vec\phi)$ makes $\langle\vec\phi_u\rangle$ purely imaginary. They are still averages over the same positive ensemble, so this is not a sign problem in the sampling sense. Point 2 (light-cone average phase $0.78,\dots,0.18$) was already in ch. 48 and concerns a different regularization. Point 3 (real time) is not a sign problem. Changes: the claim is now scoped to the Euclidean lattice formulation and vertex-operator observables (ch. 44 abstract and row H2; a new paragraph in `sec-positivity`, ch. 45).

**Finding 4 (coupling profile): confirmed for ch. 44.** Ch. 46 (`sec-wick-ordering-statement`) already says that the theory of Theorems `thm-below-2pi` and `thm-window` has the coupling $zR_\Omega^{x_j}$, vanishing at the walls; ch. 44 did not. The referee's "exists solely because" is overstated below $2\pi$, where `prp-necessity` shows no failure of site-independent ordering. Changes: ch. 44 abstract and row H3.

**Finding 5 (non-simply-laced thresholds): confirmed.** Enumeration gives $\beta^2_{\rm UV}=4\pi$ for $b_r^{(1)}$, $c_r^{(1)}$, $f_4^{(1)}$, $g_2^{(1)}$, $a_{2r}^{(2)}$, $a_{2r-1}^{(2)}$, $d_{r+1}^{(2)}$, $e_6^{(2)}$, $d_4^{(3)}$ up to rank 5 (the referee also omitted $a_{2n}^{(2)}$, $n\ge2$); for $f_4^{(1)}$ the neutral cluster gives $44\pi/9$, as stated. A short general proof was added: with two non-longest roots ($|\vec\alpha|^2\le1$) the bound is automatic; with none, $\sum\vec\alpha_k$ is a non-zero vector of the even lattice of longest roots; with one, $c=-\vec\alpha_s\cdot\vec v\in\mathbb Z$ and $|\vec v|^2\ge\max(2,c^2)\ge2c$. For $b_n^{(1)}$ the full neutral cluster attains $4\pi$ too. Changes: the bullet in `sec-collapse-thresholds` (ch. 46), `uv_thresholds.py`.

**Finding 6 ($q^4$ versus $q^{-4}$): rejected.** The leading coefficient in $x$ of the Takács–Watts cubic is $-(\lambda-q^2)(\lambda-q^4)^2$ (`sectors/a22asym.py`), so the double root at large $x$ is $q^4=e^{2i\psi}$ with $q^2=-e^{i\psi}$. At $\xi=0.9\pi$, $x=e^{10}$ the roots are $-0.05926-1.01464i$ and $-0.05737-0.98222i$; $q^4=-0.05814-0.99831i$ and $q^{-4}=-0.05814+0.99831i$. The referee's own shifts $+0.01637$, $-0.01610$ are exactly the shifts from $q^4$. The double root at small $x$ is $q^{-4}$, so the referee probably used $1/x$ or $\bar q$. Change: the text now writes $q^4=e^{2i\psi}$ to fix the convention.

**Finding 7 (Step 3 of `thm-nn-distances` and `thm-cluster-bound`): agreed.** The coordinator re-derived Steps 3.2–3.5 (the minimax path property gives $\tilde\delta_k\ge s_{i(k)}$; integration from the leaves gives $E_j=2h_j-\gamma C_j$; $C_j=h_j+1$ above $m$ leaves gives $E_j\ge\alpha$; the hook-length identity) and the cases (i)–(iii) of the cumulative exponents in `thm-cluster-bound`, including the cancellation of the pair exponents below a node with at least $m+1$ points. No error found. All labels kept; "awaiting independent check" kept, since both readings are AI-assisted. Change: the *Checks* notes of ch. 47.

**Finding 8 (crossed-channel Krein-unitarity): superseded.** After this report, half-line open problem 7 proved `prp-boundary-krein` part 1: $\check R_{ab}^\dagger(C\otimes C)\check R_{ab}=C\otimes C$ for $a,b\in\{V,\bar V\}$ with the untwisted $C$ (`code/boundary_krein.py`, $9\times10^{-15}$), and ch. 48 already says so. The proposed H9 text ("conjectural for crossed channels") was therefore not used. Instead row H9 now cites `prp-boundary-krein`, and outlook problem 5, which still listed the crossed amplitudes as open, was corrected.

**Verdict table and prior art (sections 3–4):** no changes needed.

The gap-map figure (`ch44_f33_gap_map.svg`) was regenerated: row H9 now also points to chapter 35.
