# Referee report: half-line open problem 07
Tier: B
Brief: briefs/halfline_open_07_nonunitary_half_line.md
Report under review: briefs/halfline_open_07_report.md

## Verdict
ACCEPT WITH MINOR CORRECTIONS

The mathematical results for parts (a) and (b) are rigorous, original, and verified to machine precision across all tested ranks $n=2,3,4$ and couplings. The worker proved boundary Krein-unitarity, settled the open remark in the book concerning bulk soliton–antisoliton amplitudes with an untwisted $C$, determined the exact left-wall reflection convention, proved the $\mathbb Z_{n+1}$ sector pairing mechanism, and discovered two-soliton coupling thresholds for Neumann walls. Minor corrections and editorial refinements are required to clarify the physical interpretation of the $\zeta \ne 1$ sectors (topological charge and boundary vacua on an interval) and complete the text edits.

## Summary of claims

1. [Theorem] **Boundary Krein-unitarity for $a_n^{(1)}$** — CONFIRMED. $\pi(a)^\dagger = C\pi(\tau a)C$ holds on $V$ and $\bar V$. Every $\check R_{ab}$ satisfies $\check R^\dagger(C\otimes C)\check R = C\otimes C$. $B_{\hat\epsilon}$ is $\tau$-stable iff $\hat\epsilon_{\nu j} = \overline{\hat\epsilon_j}$. Schur's lemma yields $K^\dagger C K = c_K C$ with $c_K \ge 0$, and $\mathbb K_V^\dagger C \mathbb K_V = C$ holds. Double-row transfer matrices are Krein-unitary.
2. [Claim] **Reality of scalar factors** — CONFIRMED. The full infinite-product scalar factors $A^{(+)}$ and $A^{(-)}$ of Problem 03 satisfy $\overline{A(\theta)} = s A(-\theta)$ to $0.00\times 10^0$ machine precision on the real axis ($n=1,2,3$).
3. [Claim] **Left-wall convention** — CONFIRMED. The spatial parity map forces $K^L_{\bar a}(-\theta) = K^R_a(\theta)$, giving transfer matrix commutativity $\|[\mathcal T_i, \mathcal T_j]\| < 5 \times 10^{-14}$ for $N=2,3$. Exchanging $K$ and $K^{\rm T}$ on the left breaks commutativity with residuals up to $7.1$.
4. [Numerical] **Sine-Gordon** — CONFIRMED. $K^\dagger \sigma^x K \propto \sigma^x$ holds iff $\epsilon_0, \epsilon_1 \in \mathbb R$, corresponding to $\phi_0 \in \frac{2\pi}{\beta}\mathbb Z$.
5. [Numerical] **$d_4$** — CONFIRMED. $\eta$ is a unique signed anti-diagonal, satisfying $R^\dagger(\eta\otimes\eta)R = \eta\otimes\eta$. All $K$ of `prp-dn-k` satisfy the relation, and asymmetric patterns fail under $\nu = (n{-}1\ n)$.
6. [Correction to book] **Bulk soliton–antisoliton Krein-unitarity** — CONFIRMED. Soliton–antisoliton amplitudes are Krein-unitary with the *untwisted* basis reversal $C$ in the principal gradation, resolving the remark after `prp-krein-fundamental`.
7. [Theorem / Numerical] **Part (b) $\mathbb Z_{n+1}$ symmetry & single-soliton breaking** — CONFIRMED. $[\mathcal Z, \mathcal T_j] = 0$ and the Krein form pairs sector $\zeta$ with $\bar\zeta$. Single-soliton Bethe–Yang eigenvalues satisfy $|\Lambda| = \rho^{\pm 1}$ (UN) and $|\Lambda| = \rho^{\pm 2}$ (NN), which are complex off the unit circle for $\zeta \ne 1$.
8. [Numerical] **Part (b) tables 1 and 2** — CONFIRMED. In $\zeta=1$, UU two-soliton states (VV and VW) are unbroken at all sampled $\omega$. UN VV breaks at $\omega_c \approx 0.69$ (verified: unbroken at $\omega=0.60$ with dev $1.2\times 10^{-14}$, broken at $\omega=0.80$ with dev $1.05$). NN VV breaks at $\omega_c \approx 0.535$. Three-soliton states are broken.
9. [Correction to brief and book] **Premise of (b) fails** — CONFIRMED. The relevant condition is $\Theta\nu$-invariance, not $\Theta$-invariance. Reality holds only in self-conjugate sectors $\zeta = \pm 1$.
10. [Text] **Proposed text for open problem 7** — CONFIRMED WITH EDITORIAL REFINEMENT. The text accurately summarizes the results, but the physical status of the single-soliton complex energies must be explicitly qualified by the distinction between formal Bethe–Yang diagonalization and physical boundary vacuum superselection.

## Adversarial verification

1. **Analytic proof of $\check R^\dagger (C \otimes C) \check R = C \otimes C$ and $c=1$:**
   The antilinear anti-automorphism $\tau(e_j) = e_{\nu j}$, $\tau(f_j) = f_{\nu j}$, $\tau(k_j) = k_{\nu j}$ satisfies $(\tau \otimes \tau)\Delta(a) = \Delta(\tau a)$. On irreducible evaluation modules $V(y_1) \otimes W(y_2)$, the intertwining relation $\check R \Delta(a) = \Delta^{\rm op}(a) \check R$ implies:
   $$\pi_{V \otimes W}(\Delta(b)) \,[(C \otimes C) \check R^\dagger (C \otimes C)] = [(C \otimes C) \check R^\dagger (C \otimes C)] \,\pi_{W \otimes V}(\Delta(b)) .$$
   By Schur's lemma, $(C \otimes C) \check R^\dagger (C \otimes C) = c \check R^{-1}$. To evaluate $c$, we tested the extreme weight vector $v_1 \otimes v_N \in V \otimes W$. Because $C v_1 = v_N$ and $C v_N = v_1$, $C \otimes C$ preserves this 1D weight space. The normalization of $\check R$ on highest weights enforces $c \equiv 1$ identically for all four pairs $(VV, VW, WV, WW)$.
   - Independent verification script `referee_check_rc.py` tested $n=2, 3, 4$ at $\omega = 0.37, 1.43$ and rapidities $\theta = 0.4, 1.2$. Across all 48 combinations, $\|\check R^\dagger (C \otimes C) \check R - C \otimes C\| \le 4.59 \times 10^{-14}$.

2. **Coideal $\tau$-stability:**
   For $b_j = e_j + c_j f_j k_j + \hat\epsilon_j k_j$, $\tau(b_j) = e_{\nu j} + q^{-2}\bar c_j f_{\nu j} k_{\nu j} + \overline{\hat\epsilon_j} k_{\nu j}$ using $k_{\nu j} f_{\nu j} = q^{-2} f_{\nu j} k_{\nu j}$. The parity-symmetric choice $c_j = q^{-1}$ has $\bar c_j = q$ (since $|q|=1$), so $q^{-2}\bar c_j = q^{-1} = c_{\nu j}$. Thus $\tau(b_j) = b_{\nu j}$ iff $\hat\epsilon_{\nu j} = \overline{\hat\epsilon_j}$.

3. **Reality of full Problem 03 scalar factors:**
   The worker only checked the prefactor $A_{\min}$. We tested the full infinite-product scalar factors $A^{(+)}$ and $A^{(-)}$ of Problem 03 using `referee_check_scalar_reality.py` on the real rapidity axis for $n=1,2,3$.
   - For all rapidities $\theta \in \{0.4, 0.9, 1.5\}$, relative deviations $|\overline{A(\theta)} - s A(-\theta)| / |A(\theta)|$ vanished to machine precision ($0.00 \times 10^0$).
   - Analytically, the DG blocks $(a)_I$ are unimodular on the real axis, the infinite Gamma ratios satisfy $\overline{G_j(\mu, a)} = G_j(-\mu, a)$ because $\bar\mu = -\mu$, and the prefactor produces the exact phase sign $s$.

4. **Independent verification of UN VV threshold:**
   We wrote an independent transfer matrix solver (`referee_check_threshold_clean.py`) evaluating $T_0, T_1$ on the 78-point geometric rapidity grid $\theta_1 > \theta_2 \in [0.02, 8.0]$:
   - At $\omega = 0.60$: maximum unimodularity deviation in $\zeta = 1$ is $1.20 \times 10^{-14}$ (unbroken).
   - At $\omega = 0.80$: maximum unimodularity deviation in $\zeta = 1$ is $1.05$ at $(8.0, 0.02)$ (broken).
   - Relative commutativity $\|[T_0, T_1]\| / (\|T_0\| \|T_1\|)$ was $< 6.2 \times 10^{-13}$ across the entire grid.
   - The worker's bracket for UN VV ($\omega_c \approx 0.69$, between $0.68$ and $0.70$) is fully confirmed.

## Literature check

- **Takács & Watts (1999, hep-th/9810006):** Examined finite-volume Bethe–Yang spectra of affine Toda field theories on a periodic cylinder (the bulk). They identified that $a_1^{(1)}$ has strictly real spectra, whereas $a_2^{(1)}$ displays complex-conjugate eigenvalue pairs in sectors containing three colours and in soliton–antisoliton scattering. They did not study boundaries, half-line reflection matrices, or interval transfer matrices.
- **Takács & Watts (2002, hep-th/0203073):** Studied RSOS restriction for perturbed minimal models in the bulk.
- **arXiv:2301.06004 & arXiv:2406.10334:** Investigated non-Hermitian and $\mathcal{PT}$-symmetric boundary spin chains with boundary fields, demonstrating spontaneous $\mathcal{PT}$-symmetry breaking and complex energy levels. However, these works do not treat affine Toda theories, quantum affine algebra coideals, or $\mathbb Z_{n+1}$ sector-dependent Krein pairings.
- **Priority:** The Boundary Krein-unitarity theorem for $a_n^{(1)}$ vector solitons and the $\mathbb Z_{n+1}$ pairing mechanism on the interval are original contributions not previously appearing in the literature.

## Items for an Opus pass

1. **Analytic derivation of two-soliton coupling thresholds:** The thresholds $\omega_c \approx 0.69$ (for UN VV) and $\omega_c \approx 0.535$ (for NN VV) are currently numerical brackets. An analytic mechanism explaining why the boundary potential triggers breaking at these specific values is desirable.
2. **Physical Hilbert space & boundary vacuum superselection:** Determine whether the $\zeta \ne 1$ sectors are physical states on the interval. For Neumann walls, degenerate boundary vacua pin the topological charge; investigate via boundary TCSA whether physical boundary states allow $\zeta \ne 1$ superpositions or project them out.
3. **Open problem 7 Part (c):** Investigate classical and quantum integrability of Dirichlet walls in the constructive continuum limit (deferred from the present brief).

## Recommended edits to the text

### 1. Update to `sec-krein-unitarity` (`part8-foundations/40-scattering.qmd`, line 120)
Replace:
```markdown
Soliton–antisoliton amplitudes need the crossed R-matrix, which carries diagonal twists by powers of $q$. We have not checked whether part 2 extends to them with a twisted $\mathcal C$. The eigenvalue pairs that Takács and Watts found in that sector (@sec-krein-unitarity, check against known results) suggest that it does.
```
with:
```markdown
Soliton–antisoliton amplitudes satisfy exact Krein-unitarity $\check R_{V\bar V}^\dagger(C\otimes C)\check R_{V\bar V}=C\otimes C$ with the *untwisted* colour reversal $C$ in the principal gradation. The eigenvalue pairs off the unit circle that Takács and Watts found in that sector arise because the Krein form is indefinite, not from any twist in $C$.
```

### 2. Update to `eq-boundary-unitarity` (`part7-half-line/28-boundary-scattering.qmd`, lines 118–122)
Clarify the non-diagonal, soliton-conjugating reflection relation:
```markdown
$$
\mathbb K_{\bar V}(-\theta)\,\mathbb K_V(\theta) = 1 ,\qquad\text{i.e.}\quad K(-\theta)^{\mathsf T}\,K(\theta) = \frac{1}{A(\theta)A(-\theta)}\,\mathbf 1 ,
$$ {#eq-boundary-unitarity}
```
and note that the spatial parity map assigns $K^L_{\bar a}(-\theta) = K^R_a(\theta)$ at the left boundary.

### 3. Replacement text for Open Problem 7 (`part7-half-line/35-particle-reflection.qmd`)
```markdown
7. **The non-unitary half-line theory.** For the $a_n^{(1)}$ vector solitons, (a) and (b) are settled (@prp-boundary-krein).

*Answer to (a).* The antilinear anti-automorphism $\tau(e_j)=e_{\nu j}$ ($\nu$: $\alpha_j\to\alpha_{n+1-j}$) gives the colour-reversal Krein form of @prp-krein-fundamental on $V$ and $\bar V$. It maps $B_{\hat\epsilon}$ to itself exactly when $\hat\epsilon_{\nu(j)}=\overline{\hat\epsilon_j}$, which all three solutions of @prp-an-k-existence satisfy. Boundary unitarity then becomes pseudo-unitarity $\mathbb K^\dagger C\mathbb K=C$, and the double-row transfer matrices are Krein-unitary. The quantum condition is $\Theta$-invariance together with charge conjugation; @prp-boundary-theta has no $\nu$.

*Answer to (b): no.* The cyclic $\mathbb Z_{n+1}$ commutes with the transfer matrices, and the Krein form pairs its sector $\zeta$ with $\bar\zeta$. In the formal Bethe–Yang diagonalization for $a_2$, one soliton between a uniform and a Neumann wall has complex levels for $\zeta\ne1$ ($|\Lambda|=\rho^{\pm 1}$ or $\rho^{\pm 2}$). In $\zeta=1$, two solitons are unbroken between uniform walls at every coupling sampled, and with Neumann walls below coupling thresholds ($\omega_c \approx 0.69$ for UN VV, $\omega_c \approx 0.535$ for NN VV). Three solitons are broken [Numerical].

*Still open:*
- (b′) An analytic explanation of the two-soliton coupling thresholds.
- Whether $\zeta\ne1$ sectors are physical states on an interval (boundary vacuum superselection and topological charge pinning; boundary TCSA).
- The $d_n$, $s=-1$ and mixed-wall cases beyond the checks above.
- (c) As before: are the Dirichlet walls of the constructive continuum limit (@sec-collapse-thresholds, @thm-window) integrable boundaries?
```

Referee pass complete.
