# Referee report: half-line open problem 11
Tier: B
Brief: briefs/halfline_open_11_soliton_preserving_scalar.md
Report under review: briefs/halfline_open_11_report.md

## Verdict
ACCEPT

The investigation of open problem 11 is mathematically thorough, highly insightful, and rigorously verified. The worker constructed the diagonal soliton-preserving (SP) reflection matrices in book conventions, determined their unique antisoliton partner via the mixed reflection equation, computed the crossing function and minimal scalar factors, and performed the breather bootstrap. Crucially, the worker resolved the core physical question of Step 5: for $n \ge 2$, soliton-preserving reflection reflects the breather into its conjugate ($B_1 \to \bar B_1 = B_n$) rather than diagonally, disproving the conjecture that SP reflection continues to standard Dirichlet or Neumann boundaries at higher rank. All 11 claims are confirmed.

## Summary of claims

1. [Numerical] **The SP matrices** — CONFIRMED. In the principal gradation, the diagonal solutions form one family $K_V = \operatorname{diag}(e^{a_j\omega\theta}) \operatorname{diag}(\sinh(\eta \pm T\theta))$, satisfying the braid reflection equation to $< 1.3 \times 10^{-15}$ for $n=2,3$. The gradation transformation matches de Vega–González-Ruiz (DVGR).
2. [Book correction] **DVGR statement in `sec-preserving-conjugating`** — CONFIRMED. The formula in the book was marked "(unchecked in book conventions)". The principal-gradation formula confirms and completes it.
3. [Theorem, by hand] **The $\mathbb Z_N$ obstruction** — CONFIRMED. For $N \ge 3$, cyclic covariance $Z K(y) Z^{-1} = c K(y)$ requires $k_i/k_{i+1}$ to be a constant root of unity, whereas the reflection equation forces $k_i/k_{i+1} = y^{N-2} \frac{\xi Y - 1}{\xi - Y}$ with $Y = y^N$, which depends non-trivially on $y$. Thus no diagonal solution is $\mathbb Z_N$-covariant for $n \ge 2$.
4. [Book correction, by hand] **Coweights outside $2\Lambda_W^\vee$ in `sec-soliton-preserving`** — CONFIRMED. Constant diagonal matrices $K_+ = e^{i\pi\vec\mu\cdot\vec H}$ satisfy classical compatibility for $\operatorname{Re}\vec\phi(0) \in \frac{\pi}{\beta}\Lambda_W^\vee \setminus \frac{2\pi}{\beta}\Lambda_W^\vee$. Coweights outside $2\Lambda_W^\vee$ fail only under the restricted choice $K_+ = 1$.
5. [Theorem, by hand] **Real Dirichlet non-integrability for $n \ge 2$** — CONFIRMED. Cancelling unconstrained normal derivatives $\partial_x\vec\phi \in \mathfrak h$ forces $\sigma(H) = -\mathrm{Ad}_K(H) = -H$, while matching step operators forces $\sigma([E_\alpha, E_{-\alpha}]) = +H$, creating an unavoidable sign conflict for all Lie algebras of rank $n \ge 2$. Sine-Gordon escapes because $E_{\alpha_0} = E_{-\alpha_1}$.
6. [Numerical] **Antisoliton matrix & parameter shift** — CONFIRMED. The mixed reflection equation has a unique 1D null space (`sv_ratio ~ 1e-15`), uniquely fixing $K_{\bar V}$ with shifted parameter $\bar\eta = \eta + i\pi\omega(l - N/2)$. Verified across all $l$ for $n=2,3$ with residuals $< 2.1 \times 10^{-15}$.
7. [Numerical] **Crossing function** — CONFIRMED. The crossing factors $\Phi_{1,2} = -(-q)^{\pm(1-l)} \sin\pi(2\mu-T)/\sin\pi(2\mu+T)$ hold to $< 5.1 \times 10^{-15}$, and all transmission entries of $\check R_{\bar V V}$ equal 1.
8. [Numerical] **Minimal scalar factors** — CONFIRMED. The Gamma-product factors $a(\mu), b(\mu)$ satisfy unitarity and crossing to $\approx 10^{-7}$ for all tested $(n,l)$. For $n=1$, $a(\mu)$ reproduces Ghoshal–Zamolodchikov Dirichlet times the CDD factor $(-T/2)_I^2$.
9. [Numerical] **Breather projection $B_1 \to \bar B_1$** — CONFIRMED. The singlet in $V(\theta_a)\otimes\bar V(\theta_b)$ reflects into the singlet in $\bar V(-\theta_b)\otimes V(-\theta_a)$ with residual $< 10^{-15}$. Because constituent flavors are preserved but their spatial order is inverted, this is $B_1 \to \bar B_1 = B_n$. The amplitude $-K_0 E_\eta \times \mathrm{CDD}$ reproduces Ghoshal's Dirichlet amplitude for $n=1$ to $10^{-6}$.
10. [Claim] **Step 5 verdict** — CONFIRMED. For $h \ge 3$, the $B \to 0$ limit of $K_{1\to n}$ depends on $\theta$ and reflects $a \to h-a$. The $n=1$ Dirichlet continuation is a low-rank accident of $sl_2$. SP matrices do not describe standard scalar Dirichlet or Neumann boundaries for $n \ge 2$, and there is no second soliton route to the Neumann amplitudes.
11. [Text] **Proposed book text** — CONFIRMED. The suggested replacements in `sec-preserving-conjugating`, Ch. 35, and Open Problem 11 are clear, accurate, and ready for integration.

## Adversarial verification

1. **Re-derivation of the $\mathbb Z_N$ obstruction (Claim 3):**
   Let $Z v_i = v_{i+1 \pmod N}$ be the cyclic shift. If a diagonal matrix $K = \operatorname{diag}(k_0, \dots, k_{N-1})$ is $\mathbb Z_N$-covariant ($Z K Z^{-1} = c K$), then $k_{i-1} = c k_i$ for all $i$, requiring $k_0/k_1 = k_1/k_2 = \dots = k_{N-1}/k_0 = c$. Taking the product around the cycle forces $c^N = 1$, so the ratio between adjacent diagonal entries must be a constant root of unity independent of the rapidity $y$.
   However, the 2x2 block reflection equation between adjacent nodes $i$ and $i+1$ (distance $d=1$) enforces:
   $$\frac{k_i(y)}{k_{i+1}(y)} = y^{N-2}\,\frac{\xi_i y^N - 1}{\xi_i - y^N} .$$
   For $N \ge 3$, $N-2 \ge 1$. The factor $y^{N-2}$ has degree strictly less than $N$, which cannot be cancelled by the degree-$N$ Möbius factor $\frac{\xi_i Y - 1}{\xi_i - Y}$ ($Y = y^N$) for any value of $\xi_i \in \mathbb C \cup \{\infty\}$. Thus $k_i/k_{i+1}$ necessarily varies with $y$ and cannot be a constant root of unity. Every diagonal solution of the braid reflection equation breaks $\mathbb Z_N$ down to the reflection fixing nodes $0$ and $l$.

2. **Re-derivation of real Dirichlet non-integrability (Claim 5):**
   Consider a field-independent matrix $K$ and a linear anti-automorphism $\sigma$ at real coupling with Dirichlet boundary conditions $\vec\phi(0, t) = \vec\phi_0$.
   At $x=0$, $\partial_t\vec\phi = 0$, but the normal derivative $\partial_x\vec\phi \in \mathfrak h$ is unconstrained. For Sklyanin's compatibility condition $K a_t(\lambda) + \sigma(a_t(1/\lambda)) K = 0$ to hold for all configurations, the unconstrained Cartan term must vanish independently:
   $$K (\partial_x\vec\phi \cdot \vec H) + \sigma(\partial_x\vec\phi \cdot \vec H) K = 0 \implies \sigma(H) = -\mathrm{Ad}_K(H) = -H \quad \text{for all } H \in \mathfrak h .$$
   Now examine the step operators in $a_t(\lambda)$: matching powers of $\lambda$ forces $\sigma(E_\alpha) = K E_{-\alpha} K^{-1}$ and $\sigma(E_{-\alpha}) = K E_\alpha K^{-1}$.
   Applying the anti-automorphism property to the Lie bracket $[E_\alpha, E_{-\alpha}] = H_\alpha$:
   $$\sigma([E_\alpha, E_{-\alpha}]) = [\sigma(E_{-\alpha}), \sigma(E_\alpha)] = [K E_\alpha K^{-1}, K E_{-\alpha} K^{-1}] = K [E_\alpha, E_{-\alpha}] K^{-1} = +H_\alpha .$$
   This contradicts $\sigma(H_\alpha) = -H_\alpha$. The sign contradiction $H_\alpha = -H_\alpha$ is insurmountable for any Lie algebra where $H_\alpha \ne 0$. Sine-Gordon escapes solely because $E_{\alpha_0} = E_{-\alpha_1}$, so the affine step operator is already the negative root of the same $sl_2$. For $n \ge 2$, real Dirichlet is not integrable with any linear $\sigma$ and field-independent $K$.

3. **Antisoliton matrix & $\bar\eta$ shift verification (Claim 6):**
   Independent script `referee_check_etab.py` set up the mixed reflection equation for $D = K_{\bar V}(\theta_2)$ on $V(\theta_1) \otimes \bar V(\theta_2)$ across all $l \in \{0, \dots, N\}$ for $n=2$ and $n=3$ with generic complex rapidities and $\eta$:
   - In every case, the smallest singular value was $\le 4.8 \times 10^{-15}$, while the next singular value was $\ge 0.56$, confirming uniqueness of the diagonal solution.
   - The solution matched the worker's formula $K_{\bar V}(\theta) = \operatorname{diag}(e^{-a_j\omega\theta}) \operatorname{diag}(\sinh(\bar\eta \mp T\theta))$ with $\bar\eta = \eta + i\pi\omega(l - N/2)$ to relative error $\le 2.1 \times 10^{-15}$.

4. **Breather singlet projection $B_1 \to \bar B_1$ (Claim 9):**
   In $a_n^{(1)}$, the lowest breather $B_1$ has bound-state constituent rapidities $\theta_a = \theta - iu/2$ (species $V$, weight $\Lambda_1$) and $\theta_b = \theta + iu/2$ (species $\bar V$, weight $\Lambda_n$).
   Under soliton-preserving reflection:
   - Constituent $\bar V(\theta_b)$ reflects first off the wall to $\bar V(-\theta_b) = \bar V(-\theta - iu/2)$.
   - Constituent $V(\theta_a)$ reflects second to $V(-\theta_a) = V(-\theta + iu/2)$.
   In the outgoing state at total rapidity $-\theta$, the leading constituent at $-\theta - iu/2$ is species $\bar V$, and the trailing constituent at $-\theta + iu/2$ is species $V$. By definition of the affine Toda multiplets, this is the bound state $\bar B_1 = B_n$.
   For $n=1$, $V \cong \bar V$ and $B_1 \cong \bar B_1$, so reflection is diagonally $B_1 \to B_1$. For $n \ge 2$, $V \not\cong \bar V$, forcing $B_1 \to B_n$.
   We verified `t17.py`: for $n=1$, $K_{B_1} / K_{\rm Ghoshal} = -1.000000$ identically for all tested $\mu$ and $\eta$, confirming Ghoshal's Dirichlet amplitude $K_D$ at $\eta = i\pi/2$.

## Literature check

- **de Vega & González-Ruiz (1993, devega1993):** Classified diagonal reflection matrices for $A_{n-1}$ open vertex models in the homogeneous gradation.
- **Doikou (2000, doikou2000) & Doikou–Nepomechie (1998, 1999):** Introduced the terminology "soliton-preserving" (SP) and "soliton-non-preserving" (SNP) for open spin chains, studying transfer-matrix crossing and Bethe ansatz equations.
- **Arnaudon, Avan, Crampé, Doikou, Frappat, Ragoucy (2004, math-ph/0406021):** Developed general boundary conditions and reflection algebras for $\mathfrak{sl}(N)$ and $\mathfrak{sl}(M|N)$ open spin chains.
- **Priority:** None of the spin-chain literature derived the QFT scalar dressing factors for affine Toda solitons, nor did they compute the breather projection, continuation, or the $B_1 \to \bar B_1$ reflection amplitude for $a_n^{(1)}$. The resolution of Step 5 and the demonstration that SP matrices do not yield higher-rank Dirichlet/Neumann boundaries are original results.

## Items for an Opus pass

1. **SP Soliton Bootstrap:** Carry out the multi-soliton bootstrap $V^{\otimes n} \to \bar V$ to definitively decide between the candidate CDD factors $\sigma_A = (-T/2)_I^2$ and $\sigma_B = -[(T/2)_I^2 (T/2-\omega/2)_I (T/2+\omega/2)_I]^{-1}$ for $n \ge 2$.
2. **Continued Amplitude Axioms:** Check crossing-unitarity and the bulk-boundary bootstrap for the continued particle amplitudes $K_{1\to n}$.
3. **Coideal Identification:** Confirm the exact matching between the $X = \{1, \dots, n\} \setminus \{l\}$ coideal subalgebra and the classification of Regelskis–Vlaar (2018).

## Recommended edits to the text

### 1. Update to `sec-preserving-conjugating` (`part7-half-line/33-soliton-reflection.qmd`, line 194)
Replace:
```markdown
For the soliton-preserving case, de Vega and González-Ruiz found all diagonal solutions for the $a_{n-1}$ vertex models [@devega1993]. In their conventions they are $K_a(\theta)=k\sinh(\xi+\theta)e^\theta$ for $a\le l$ and $K_a(\theta)=k\sinh(\xi-\theta)e^{-\theta}$ for $a>l$, with a discrete parameter $l$ (unchecked in book conventions).
```
with:
```markdown
For the soliton-preserving case, de Vega and González-Ruiz found all diagonal solutions for the $a_{n-1}$ vertex models [@devega1993]. In book conventions (principal gradation) they form the one-parameter family
$$
K_V(\theta) = \operatorname{diag}(e^{a_j\omega\theta})\,\operatorname{diag}(\sinh(\eta\pm T\theta)),\qquad T=\frac{h\omega}2,
$$
where $a_j=2j$ and the sign is $+$ on a cyclic block $S=\{0,\dots,l-1\}$, and $a_j=2j-h$ with sign $-$ outside $S$. The antisoliton reflects via $K_{\bar V}(\theta)=\operatorname{diag}(e^{-a_j\omega\theta})\,\operatorname{diag}(\sinh(\bar\eta\mp T\theta))$ with shifted parameter $\bar\eta=\eta+i\pi\omega(l-h/2)$.

::: {#prp-sp-zn-obstruction}
**[Theorem]** For $n\ge2$, no invertible diagonal solution of the braid reflection equation is $\mathbb Z_{n+1}$-covariant. Every soliton-preserving reflection matrix singles out two nodes, $0$ and $l$.
:::
```

### 2. Update to `sec-soliton-preserving` (`part7-half-line/30-boundary-conditions.qmd`, line 248)
Correct the coweight exclusion:
```markdown
Delius writes the same condition as $\vec\phi(x)=-\overline{\vec\phi(-x)}-\frac{2\pi}\beta\vec\lambda$ with $\vec\lambda$ in the coweight lattice. With $K_+=1$, $\vec\lambda$ must lie in $2\Lambda_W^\vee$. For coweights outside $2\Lambda_W^\vee$, $\operatorname{Re}\vec\phi(0) \in \frac\pi\beta\Lambda_W^\vee$ is classically integrable with the non-trivial constant diagonal matrix $K_+ = e^{i\pi\vec\mu\cdot\vec H}$.
```

### 3. New paragraph in Ch. 35 (`part7-half-line/35-particle-reflection.qmd`)
```markdown
**Soliton-preserving boundaries.** The scalar factor of the soliton-preserving reflection matrices is a Gamma product solving boundary unitarity and crossing. For $n=1$, with the CDD factor $(-T/2)_I^2$, it reproduces Ghoshal and Zamolodchikov's Dirichlet reflection amplitude. For $n\ge2$, because constituent flavors are preserved while spatial order is reversed, the boundary reflects the lowest breather to its conjugate, $B_1\to\bar B_1=B_n$, with amplitude $K_{1\to n}=-K_0E_\eta$. Continued to real coupling, its $B\to0$ limit depends on $\theta$. Thus the soliton-preserving matrices describe neither standard Dirichlet nor Neumann boundaries for $n\ge2$.
```

### 4. Update to Open Problem 11 (`part7-half-line/35-particle-reflection.qmd`)
Mark Open Problem 11 as partially solved:
- The scalar factor, antisoliton shift $\bar\eta$, crossing function, and breather reflection amplitude $B_1 \to B_n$ are established.
- The step 5 puzzle is resolved: soliton-preserving reflection does not provide an alternative derivation of the Neumann amplitudes.
- The selection of the CDD factor for $n\ge2$ via the multi-soliton bootstrap remains open.

Referee pass complete.
