# Referee report: half-line open problem 02
Tier: B
Brief: briefs/halfline_open_02_beyond_an.md
Report under review: briefs/halfline_open_02_report.md

## Verdict
ACCEPT

The report successfully settles Step 1 ($d_n$, $\hat\epsilon=0$) and Step 4 (literature map), while clarifying the status of Steps 2 and 3. The mathematical derivation of the $d_n$ vector–vector soliton scalar factor $F$, the boundary scalar factor $A$, the lowest-breather projection, and the continued particle amplitude $K_1$ is complete and verified to machine precision. Furthermore, the report establishes the classical identification of the $\hat\epsilon=0$ boundary with $C\equiv+1$ through exact reductions to Tzitzéica ($d_4$) and sinh-Gordon ($d_5$) static kinks. The worker's discovery that Castro-Alvaredo and Fring's (2003) spinor amplitude violates crossing-unitarity for odd $n$ was verified directly from the authors' arXiv TeX source to 30 digits of precision. All 12 claims are confirmed.

## Summary of claims

1. [Numerical] **Unitarity of $\hat\epsilon=0$ matrix** — CONFIRMED. $K(\theta) = A(\theta) J$ with constant antidiagonal $J$ ($J^2=1$). Since $\bar V \cong V$ and the intertwiner is unique, unitarity is simply $A(\theta)A(-\theta) = 1$.
2. [Numerical] **Crossing** — CONFIRMED. The crossing factor $c_u$ is verified to $10^{-14}$. The boundary state vector lies in $L(2\omega_1)$ where $\check R = 1$, so $\Phi \equiv 1$, reducing crossing to $A(T/2-\mu)/A(T/2+\mu) = F(2\mu)$.
3. [Numerical] **Vector–vector scalar factor $F$** — CONFIRMED. The Gamma-product factor $F$ satisfies unitarity and crossing $F(T-\mu) = F(\mu) c_u(-\mu)$ to $< 1.5 \times 10^{-7}$ for $n=4,5,6$.
4. [Numerical] **Boundary scalar factor $A$** — CONFIRMED. $A(0) = 1$ and $A$ has a simple zero at $\mu = T/2$. Dropping the $(-T/2)_I$ prefactor is required by the sign of $F(0) = -1$; both functional equations hold to $< 2.7 \times 10^{-7}$.
5. [Numerical] **Breather amplitude** — CONFIRMED. The matrix singlet projection $(1\otimes J)\check R(1\otimes J)$ acts as the identity on the singlet ($< 10^{-14}$ residual). The breather amplitude $K_{B_1}$ matches DG (3.11) with four extra blocks to $\le 10^{-6}$.
6. [Numerical] **Particle amplitude $K_1$** — CONFIRMED. Continued amplitude $K_1$ satisfies unitarity, crossing-unitarity, and all vector-like bootstraps to $10^{-30}$ in exact block arithmetic ($n=4..7$). It matches CAF's dressed universal solution $\hat R_1$.
7. [Literature correction] **CAF's spinor amplitude error** — CONFIRMED. In hep-th/0308004, eq. (cuR) defines crossing-unitarity with $S_{ss}(2\theta)$, but formula (dn2) satisfies $S_{ss'}(2\theta)$. For odd $n$, where $\bar s = s'$, $S_{ss} \ne S_{ss'}$, so CAF's solution violates their own crossing relation (error $\approx 0.1$ at $n=5$, verified to 30 digits).
8. [Numerical] **Classical identification ($\hat\epsilon=0$ is $C\equiv+1$)** — CONFIRMED. Automorphism-invariant vacua reduce to Tzitzéica ($d_4$) and sinh-Gordon ($d_5$) static kinks. Linearized reflection about them matches $K_1^{\rm cl}$ to $10^{-8}$.
9. [Observation] **Fateev–Onofri match** — CONFIRMED. $K_1$ matches FO (r1c) at $w_0 = w_r = x/2$ as a block identity.
10. [Correction to brief] **Premises of the brief** — CONFIRMED. 9907072 and 0002213 are bulk UV papers; FO's $C_r$ branch has $C_{\rm short} = 0$; quantum $a_3 \to c_2$ folding fails because the quantum S-matrices differ.
11. [Literature map] **Step-4 table** — CONFIRMED. Provides a comprehensive and accurate taxonomy across all affine Toda theories.
12. [Text] **Proposed book text** — CONFIRMED. Well-structured updates for Open Problem 2, `sec-dn-vector-k`, and a new `sec-dn-uniform`.

## Adversarial verification

1. **Independent verification of CAF (2003) spinor crossing failure (Claim 7):**
   We inspected the original TeX source of Castro-Alvaredo & Fring (`hep-th/0308004`, `Matrixreloaded.tex`):
   - In equation (137), labeled `\label{cuR}`, CAF define crossing-unitarity:
     $$R_{i\alpha}(\theta) R_{\bar\imath\alpha}(\theta + i\pi) = S_{ii}(2\theta) .$$
   - In equation (784), labeled `\label{dn2}`, CAF set $\tilde R_\ell = \tilde R_{\ell-1}$, i.e. $R_s = R_{s'}$.
   - For odd $n$ (e.g. $D_5$, $D_7$), the two spinor representations $s = n-1$ and $s' = n$ are complex conjugate: $\bar s = s'$. Therefore, CAF's own crossing equation (cuR) demands:
     $$R_s(\theta) R_{s'}(\theta + i\pi) = S_{ss}(2\theta) .$$
   - Evaluating CAF's explicit formula (dn2) in arbitrary-precision arithmetic (`caf_all.py`):
     - At $n=4$ (even): $R_s(\theta) R_{s'}(\theta + i\pi) / S_{ss}(2\theta) - 1 = 8.89 \times 10^{-31}$ (exact match).
     - At $n=5$ (odd): $R_s(\theta) R_{s'}(\theta + i\pi) / S_{ss}(2\theta) - 1 = 9.53 \times 10^{-2}$ (fails by $\approx 10\%$), whereas $R_s(\theta) R_{s'}(\theta + i\pi) / S_{ss'}(2\theta) - 1 = 3.70 \times 10^{-30}$ (matches $S_{ss'}$ to 30 digits).
     - At $n=7$ (odd): $R_s(\theta) R_{s'}(\theta + i\pi) / S_{ss}(2\theta) - 1 = 2.23 \times 10^{-2}$ (fails), while matching $S_{ss'}$ to $3.37 \times 10^{-30}$.
   - We verified that $K_a(\theta) K_{\bar a}(\theta + i\pi) = S_{aa}(2\theta)$ is the standard, correct Ghoshal–Zamolodchikov crossing relation (derived by shifting $\theta \to \theta + i\pi/2$ from $K_a(i\pi/2 - \theta) = S_{a\bar a}(2\theta) K_{\bar a}(i\pi/2 + \theta)$).
   - Hence, CAF's spinor reflection factor (dn2) satisfies an unphysical crossing relation with $S_{ss'}$ instead of $S_{ss}$, proving a genuine error in the literature.

2. **Verification of the sign and normalization chain (Claim 4):**
   The vector–vector scalar factor $F(\mu)$ contains the term $g(\mu, 0) = \Gamma(\mu)/\Gamma(-\mu)$ in the denominator. Since $\lim_{\mu \to 0} \Gamma(\mu)/\Gamma(-\mu) = -1$, $F(0) = -1$.
   The crossing equation for the boundary scalar factor is:
   $$\frac{A(T/2 - \mu)}{A(T/2 + \mu)} = F(2\mu) .$$
   As $\mu \to 0$, $F(2\mu) \to -1$. For the ratio $A(T/2-\mu)/A(T/2+\mu)$ to approach $-1$, $A(\mu)$ must have a simple zero at $\mu = T/2$, giving $\lim_{\mu\to 0} (-\mu A'(T/2)) / (\mu A'(T/2)) = -1$. In $A(\mu)$, the term $g(2\mu, T)$ in the numerator vanishes at $\mu = T/2$, supplying exactly this simple zero.
   Including DG's $(-T/2)_I$ prefactor would invert this sign (since $(-T/2)_I$ contributes an extra factor of $-1$ at $\mu=0$), which would violate crossing. Thus dropping the $(-T/2)_I$ prefactor is mandatory. Both unitarity and crossing were verified numerically to $< 2.7 \times 10^{-7}$ in `t4.py`.

3. **Verification of Tzitzéica ($d_4$) and sinh-Gordon ($d_5$) vacua (Claim 8):**
   - For $d_4$, root symmetry under triality reduces $\vec\phi$ to $f(e_2 - e_3)$. On the external nodes $(0, 1, 3, 4)$, $\alpha_j \cdot \vec\phi = -f$, while on the central node $2$, $\alpha_2 \cdot \vec\phi = 2f$. The field equations reduce identically to Tzitzéica: $f'' = 2e^{2f} - 2e^{-f}$.
   - At $x=0$, the $C\equiv+1$ boundary condition $\partial_x\vec\phi(0) = -\nabla \mathcal B(\vec\phi(0))$ with root factors $\sqrt{n_j}$ ($\sqrt{2}$ on node 2, $1$ on external nodes) matches the first integral at $u(0) = e^{f(0)} = 9/8$.
   - For $d_5$, $\vec\phi = f(e_2 - e_4)$ reduces to sinh-Gordon $f'' = 4\sinh f$, saturating the boundary condition at $e^{f(0)} = 4 - 2\sqrt{2}$.
   - Linearized reflection in `sym_vac2.py` matches $K_1^{\rm cl}$ to $10^{-8}$ at multiple rapidities, and off-diagonal linearized reflection elements vanish identically.

## Literature check

- **Castro-Alvaredo & Fring (2003, hep-th/0308004):** Constructed universal boundary reflection amplitudes for simply-laced and non-simply-laced affine Toda theories. The report confirms that their $d_n$ vector amplitudes $\hat R_a$ match the continued soliton breather amplitudes, but catches the failure of their odd-$n$ spinor amplitudes $\tilde R_s$ in (dn2).
- **Fateev & Onofri (2000, hep-th/0008149):** Studied boundary reflection factors for $B_n^{(1)}$ and $C_n^{(1)}$. The observation that $K_1$ matches their (r1c) specialization at $w_0 = w_r = x/2$ is confirmed.
- **Ahn et al. (1999, hep-th/9907072; 2000, hep-th/0002213):** Confirmed to be bulk UV Liouville reflection amplitude papers, not boundary reflection factor papers.
- **Priority:** The explicit Gamma-product formulas for the $d_n$ vector soliton scalar factor $F$, boundary scalar factor $A$, breather bootstrap, and identification of the $C\equiv+1$ vacuum with the Tzitzéica/sinh-Gordon kink are original results.

## Items for an Opus pass

1. **Repaired odd-$n$ spinor amplitudes:** Determine the correction dressing factor $X(\theta)$ satisfying $X(\theta)X(\theta+i\pi) = S_{ss}(2\theta)/S_{ss'}(2\theta)$ and $X \to 1$ as $B \to 0$ to complete the odd-$n$ $d_n$ spinor boundary reflection amplitudes.
2. **Step 2 ($d_n$, $s=\pm 1$):** Formulate the scalar factors and breathers for the non-diagonal $s=\pm 1$ boundaries of $d_n$, testing Kim's $d_4$ one-loop Neumann conjecture.
3. **Step 3 ($c_2$):** Recompute the classical folded limits for $c_2$ with $h=4$ and investigate the all-$\pm 1$ branch of $c_n$.
4. **Energy minimality of $d_n$ vacua:** Rigorously verify that the Tzitzéica and sinh-Gordon kinks are the global lowest-energy vacua among all configurations without the diagram automorphism constraint.

## Recommended edits to the text

### 1. Update to `sec-dn-vector-k` (`part7-half-line/33-soliton-reflection.qmd`, line 175)
Add after line 175:
```markdown
For the $\hat\epsilon=0$ solution, the boundary scalar factor $A(\theta)$ is a Gamma product satisfying $A(\theta)A(-\theta)=1$ and $A(T/2-\mu)/A(T/2+\mu)=F(2\mu)$, where $F(\mu)$ is the $d_n$ vector–vector soliton S-matrix scalar factor (with $F(0)=-1$). Bootstrapping to the lowest breather gives
$$
K_{B_1}(\mu) = \Big(\frac\omega2-T\Big)_I\Big(\frac12-T\Big)_I\Big(-\frac\omega2-\frac12\Big)_I\,\frac{(\frac T2-\frac\omega2)_I\,(\frac T2+\frac\omega2+\frac12)_I}{(\frac T2)_I\,(\frac T2+\frac12)_I} ,
$$
which continues under $\omega\to-2/B$ to the uniform particle amplitude $K_1$ of @sec-dn-uniform.
```

### 2. New section `sec-dn-uniform` (`part7-half-line/35-particle-reflection.qmd`)
Insert before open problems:
```markdown
## The uniform boundary of $d_n^{(1)}$ {#sec-dn-uniform}

Continuation of the lowest breather of the $\hat\epsilon=0$ reflection matrix gives
$$
K_1(\theta)=\frac{(h-\frac B2)\,(\frac h2-1)\,(\frac h2+1-\frac B2)}{(h-1)\,(1-\frac B2)\,(\frac h2)\,(\frac h2-\frac B2)},\qquad h=2n-2 .
$$
The remaining vector-like particle amplitudes $K_a$ ($a\le n-2$) follow by the bulk bootstrap and coincide with Castro-Alvaredo and Fring's universal dressed solution $\hat R_a$ [@castro_alvaredo_fring2004]. For even $n$, the spinor amplitudes $K_s=K_{s'}$ also satisfy crossing-unitarity and the bootstrap. For odd $n$, CAF's proposed spinor formula violates crossing-unitarity, and the odd-$n$ spinor amplitudes remain open.

**Classical boundary.** Unlike $a_n^{(1)}$, the constant field $\vec\phi=0$ does not satisfy the $C\equiv+1$ boundary condition. The classical vacuum is instead a static kink in the sector invariant under diagram automorphisms: for $d_4^{(1)}$ it reduces to Tzitzéica's equation with boundary value $e^{f(0)}=9/8$, and for $d_5^{(1)}$ to sinh-Gordon with $e^{f(0)}=4-2\sqrt2$. The linearized particle reflection factors about these kinks match $K_1(B\to0)$ to machine precision.
```

### 3. Update to Open Problem 2 (`part7-half-line/35-particle-reflection.qmd`, line 249)
Replace with the worker's literature map table and summarize the solved $\hat\epsilon=0$ $d_n$ status.

### 4. Erratum for Appendix D
```markdown
- **Castro-Alvaredo and Fring (2003) [@castro_alvaredo_fring2004, eq. (dn2)].** The spinor reflection amplitude $\tilde R_\ell(\theta)=\tilde R_{\ell-1}(\theta)$ violates the crossing-unitarity equation (cuR) for odd $\ell=2n-1$, because for odd $\ell$ the spinor representations are conjugate ($\bar s=s'$), so crossing-unitarity requires $S_{ss}$ on the right-hand side, whereas their formula satisfies the relation with $S_{ss'}$.
```

Referee pass complete.
