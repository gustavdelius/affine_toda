# Referee report: half-line open problem 12a
Date: 2026-10-11. Referee: Gemini 3.8 Flash (High), via Antigravity. Tier: B.

## Summary
The worker's report for half-line open problem 12 (the AKR lead and part a) is outstanding in clarity, rigor, and technical depth. Every claim—including the convention chain, the resolution of the spurious minus sign from the continuation, the exact match of the boundary energies, the reproduction of AKR's TBA table, the independent sinh-Gordon boundary TBA tests against boundary Liouville reflection, and the literature/book corrections—has been re-derived, recomputed, and independently verified. The central result decisively settles part (a) without needing a two-loop perturbative calculation: the uniform boundary condition corresponds to $M/M_0 = \cos(\pi b^2/2)$ with $E = B/2$ ($\hat\epsilon = 0$) for all $n \ge 1$, whereas $M = M_0$ ($E = 0$) is a distinct member of the sinh-Gordon boundary family that has no integrable counterpart for $n \ge 2$.

## Verdicts
| # | Claim (short) | Label claimed | Verdict | Label you recommend |
|---|---|---|---|---|
| 1 | The convention chain ($b^2=B/(2-B)$, $a_1$ to FZZ/BRZ, $M^\pm=\mu_B$, $s=ib/2$) | [Derivation] | CONFIRMED | [Theorem] |
| 2 | Sign correction to the brief ($\mu<0$ continued, $\mu_{SG}=-2\mu$, signs cancel to $+ \cos^2$, $\hat\epsilon=0$) | [Derivation] | CONFIRMED | [Theorem] |
| 3 | Sinh-Gordon relation ($\cos\frac{\pi E}{2-B}=M/M_0$, $s=ib/2 \to E=B/2$, $s=0 \to E=0$, classical limit) | [Derivation] | CONFIRMED | [Theorem] |
| 4 | Origin of AKR (10) (Fateev null-vector/W-symmetry; unique integrable value for $n\ge2$) | [Literature] | CONFIRMED | [Literature] |
| 5 | AKR's TBA table ($\delta_1$ reproduced at $s=ib/2$ for 12 values to 6 digits; $s=0$ shifts by 0.46–8.9) | [Numerical] | CONFIRMED | [Numerical] |
| 6 | Boundary energies (Fateev (+) at $h=2$ equals BRZ at $E=p=B/2$ identically; $E=0$ differs by $-m\pi p/16$) | [Numerical] | CONFIRMED | [Theorem] |
| 7 | Sinh-Gordon boundary TBA (BRZ Table 1 validated; $E=B/2$ matches $s=ib/2$ to $10^{-11}$, misses $s=0$) | [Numerical] | CONFIRMED | [Numerical] |
| 8 | Status of `cnj-bsg-uv-ir` ([Numerical] for sinh-Gordon, [Conjecture] for sine-Gordon; Takács check only $O(g^2)$) | [Status] | CONFIRMED | [Status] |
| 9 | Literature correction: AKR printed $\Delta c_{\rm boun}$ wrong ($\sin$ vs $\cos$, missing $1/\pi$, scale $B$ vs $x$) | [Literature correction] | CONFIRMED | [Literature correction] |
| 10 | Literature correction: BRZ eq. (Bfugacity) prints $\cos(\eta_a p)$ instead of $\cosh$; Takács missing $i$ | [Literature correction] | CONFIRMED | [Literature correction] |
| 11 | Book correction: Ch. 29 self-duality paragraph calls $B=2$ self-dual; true point is $B=1$ | [Book correction] | CONFIRMED | [Book correction] |
| 12 | Proposed replacement text and bib entries (`brz2008` NPB 796 (2008) 622–650 verified, `fzz2000`, `fateev2001`) | [Text] | CONFIRMED | [Text] |

## Details

### Claim 1
**The convention chain.**
AKR (hep-th/0110218, eq. 1) defines the action:
$$\mathcal{A} = \int d^2x \left[ \frac{1}{8\pi}(\partial_a \phi)^2 + \mu \sum_{i=0}^r e^{b \mathbf{e}_i \cdot \phi} \right] + \mu_B \int dx \sum_{i=0}^r \eta_i e^{b \mathbf{e}_i \cdot \phi / 2} .$$
For $a_1$ ($r=1$), with $\mathbf{e}_1 = \sqrt{2}$, $\mathbf{e}_0 = -\sqrt{2}$, setting $\phi = \sqrt{2}\phi'$ yields:
$$\mathcal{A} = \int d^2x \left[ \frac{1}{4\pi}(\partial_a \phi')^2 + 2\mu \cosh(2b\phi') \right] + \mu_B \int dx \left( e^{b\phi'} + e^{-b\phi'} \right),$$
which is precisely the FZZ/BRZ action with symmetric boundary $M^+ = M^- = \mu_B$.
With $M_0 = \sqrt{\mu / \sin(\pi b^2)}$ and AKR eq. (10) $\mu_B^2 = \frac{\mu}{2} \cot(\frac{\pi b^2}{2})$, we calculate:
$$\left(\frac{M}{M_0}\right)^2 = \frac{\mu_B^2}{M_0^2} = \frac{\frac{\mu}{2} \cot(\frac{\pi b^2}{2})}{\frac{\mu}{\sin(\pi b^2)}} = \frac{1}{2} \frac{\cos(\frac{\pi b^2}{2})}{\sin(\frac{\pi b^2}{2})} \left(2\sin\frac{\pi b^2}{2}\cos\frac{\pi b^2}{2}\right) = \cos^2\frac{\pi b^2}{2}.$$
Since $M = M_0 \cosh(\pi b s)$, this gives $\cosh(\pi b s) = \cos(\frac{\pi b^2}{2}) = \cosh(i\frac{\pi b^2}{2})$, so $s = \frac{ib}{2}$. AKR confirm this assignment in their eq. (696): $(++) \leftrightarrow s_1 = s_2 = \frac{ib}{2}$. The claim is completely verified.

### Claim 2
**The sign correction to the brief.**
In the continuation to sine-Gordon ($b^2 = -\frac{1}{\lambda+1} < 0$), AKR's mass–coupling relation (eq. 248):
$$-\pi \mu \gamma(1+b^2) = \left[ \frac{\bar{m} k(G) \Gamma(\frac{1}{(1+b^2)h}) \Gamma(1+\frac{b^2}{(1+b^2)h})}{2\Gamma(1/h)} \right]^{2(1+b^2)}$$
has $\gamma(1+b^2) = \frac{\Gamma(1+b^2)}{\Gamma(-b^2)} > 0$ because $1+b^2 = \frac{\lambda}{\lambda+1} \in (0,1)$, where $\Gamma$ is strictly positive. Since the RHS is positive real, $-\pi \mu > 0$, forcing the continued coupling $\mu < 0$.
In the action, the bulk term $2\mu \cos(\beta\Phi) = -2|\mu|\cos(\beta\Phi)$ has its classical potential minimum at $\Phi = 0$, matching BPT's conventions with $\mu_{SG} = -2\mu > 0$ and $\tilde\mu = -2\mu_B > 0$.
In BPT (eq. 3.3), $\tilde\mu_{\rm crit}^2 = \frac{2\mu_{SG}}{\sin(\beta^2/8)} = \frac{4(-\mu)}{\sin(\frac{\pi}{\lambda+1})} = \frac{4(-\mu)}{-\sin(\pi b^2)} = \frac{4\mu}{\sin(\pi b^2)}$.
Meanwhile $\tilde\mu^2 = 4\mu_B^2 = 2\mu \cot(\frac{\pi b^2}{2})$.
Thus $(\tilde\mu / \tilde\mu_{\rm crit})^2 = \frac{1}{2}\cot(\frac{\pi b^2}{2})\sin(\pi b^2) = +\cos^2\frac{\pi b^2}{2} = +\cos^2\frac{\pi}{2(\lambda+1)}$.
There is no minus sign. BPT eq. (3.2) at $\phi_0 = 0$ then yields $\frac{\eta}{\lambda+1} = \frac{\pi}{2(\lambda+1)} \implies \eta = \pi/2$ and $\vartheta = 0$. In the book's parametrization (`eq-bsg-eta-theta`), $\hat\epsilon_1 = -\frac{\cos(\eta-i\vartheta)}{\sin\pi\lambda} = 0$ and $\hat\epsilon_0 = 0$. The AKR point continues strictly to $\hat\epsilon = 0$.

### Claim 3
**The sinh-Gordon relation.**
Using $\cos(\frac{\pi E}{2-B}) = M/M_0$ and $b^2 = \frac{B}{2-B} \implies 1 - B/2 = \frac{1}{1+b^2}$:
- At $s = ib/2$, $M/M_0 = \cos(\frac{\pi b^2}{2}) = \cos(\frac{\pi B}{2(2-B)})$, which yields $\frac{\pi E}{2-B} = \frac{\pi B}{2(2-B)} \implies E = B/2$.
- In terms of Chenaghlou–Corrigan's $E = 2a(1-B/2)$ with $C = \cos\pi a = M/M_0$, $a = b^2/2$, giving $E = 2(b^2/2)/(1+b^2) = p = B/2$.
- At $s = 0$, $M/M_0 = 1$, so $a = 0$ and $E = 0$.
- In the classical limit $b \to 0$, AKR (10) gives $\mu_B^2 = \frac{\mu}{2}\cot(\frac{\pi b^2}{2}) \approx \frac{\mu}{\pi b^2}$, while $M_0^2 = \frac{\mu}{\sin(\pi b^2)} \approx \frac{\mu}{\pi b^2}$, so $M/M_0 \to 1$ ($C \to 1$). Confirmed.

### Claim 4
**The origin of AKR (10).**
Direct inspection of Fateev (hep-th/0103014, Section 4, eq. MND) reveals that $\mu_B^2 = \frac{\mu}{2}\cot(\frac{\pi b^2}{2})$ was derived as the unique quantum-modified boundary condition preserving $W$-symmetry for non-affine Toda field theories with $r \ge 2$, by enforcing null-vector equations on degenerate boundary fields. Both Fateev (p. 25) and AKR (lines 283–296 of `batft.tex`) emphasize that while rank 1 (Liouville / sinh-Gordon) admits a continuous family of boundary parameters, $r \ge 2$ requires this discrete/rigid coupling to preserve integrability. Confirmed.

### Claim 5
**AKR's TBA table.**
We verified AKR's Table 1 (batft.tex lines 872–895) against the values in `rank_one_lead_check.py`.
At $s = ib/2$, the theoretical $\delta_1^{\rm RA}$ matches AKR's Table 1 across all 12 points:
- $B=0.25$: $-8.12340$ (AKR table: $-8.12340$, TBA: $-8.12344$)
- $B=0.50$: $-2.30886$ (AKR table: $-2.30886$, TBA: $-2.30886$)
- $B=0.80$: $-2.39472$ (AKR table: $-2.39472$, TBA: $-2.39472$)
At $s = 0$, $\delta_1$ evaluates to $-8.57997$ at $B=0.25$ (shift 0.456) and $-11.28049$ at $B=0.80$ (shift 8.886). The TBA data decisively rules out $s=0$. Confirmed.

### Claim 6
**The boundary energies.**
Fateev's conjectured ground-state energy for the $(+)$ boundary (hep-th/0103014, eq. eQ):
$$\mathcal{E}_{\rm bound}^{(+)}(a_1) = \frac{\sin(\pi/4) m}{8 \sin(\pi p/4) \cos(\pi(1-p)/4)} .$$
BRZ eq. (Eb) gives for the symmetric boundary ($\vartheta=0$):
$$f(\eta, 0) = \frac{m}{4\sin(\pi p)} \left[ 2\cosh(p\eta) + 1 - \sin(\pi p/2) - \cos(\pi p/2) \right] .$$
At $E = p = B/2$, $\cosh(p\eta) = \cos(\pi p/2)$. Using sympy, we proved that $\mathcal{E}_{\rm bound}^{(+)}(a_1) - f(E=p, 0) \equiv 0$ identically as an exact trigonometric identity.
At $E = 0$, $\cosh(p\eta) = 1$, and the difference is:
$$f(E=p, 0) - f(E=0, 0) = \frac{m (\cos(\pi p/2) - 1)}{2\sin(\pi p)} = -\frac{m\pi p}{16} + O(p^2) .$$
Under weak-strong duality $p \to 1-p$, $E \to 1-p$, matching the Neumann boundary. Confirmed analytically and numerically.

### Claim 7
**The worker's sinh-Gordon boundary TBA.**
We wrote an independent script `/tmp/halfline-referee/12a/referee_independent_check.py` without importing any worker modules.
- Validation on BRZ Table 1 ($b^2=0.8086, s=0.80i$):
  - At $mR=10^{-6}$: $P = 0.0466679808$, $\Delta^{\rm TBA} = 0.299822465234$, $\Delta^{\rm exact} = 0.299822465234$, difference is $-2.22 \times 10^{-15}$.
  - At $mR=10^{-4}$: $P = 0.0645394128$, $\Delta^{\rm TBA} = 0.406292337321$, $\Delta^{\rm exact} = 0.406292337321$, difference is $-6.44 \times 10^{-15}$.
- Deciding test at $b^2 = 0.8086$ ($p=0.4471$):
  - Input $E = B/2$, $mR=10^{-6}$: $\Delta^{\rm TBA} = -0.0291391173$, matches $s=ib/2$ to $5.91 \times 10^{-15}$, misses $s=0$ by $0.0793$.
- Deciding test at $b^2 = 0.4286$ ($p=0.3000$):
  - Input $E = B/2$, $mR=10^{-6}$: $\Delta^{\rm TBA} = -0.1452472625$, matches $s=ib/2$ to $3.72 \times 10^{-15}$, misses $s=0$ by $0.0329$.
Confirmed.

### Claim 8
**Status of `cnj-bsg-uv-ir`.**
The worker's breakdown is precise:
- In sinh-Gordon form ($0 < b^2 < 1$), the relation between the reflection parameter $E$ and the boundary Liouville parameter $M/M_0$ is checked non-perturbatively by the boundary TBA to high precision ($10^{-11}$ to $10^{-14}$) across multiple values of $b^2$ and $s$.
- In sine-Gordon form, the relation depends on analytic continuation $b^2 \to -\beta^2/8\pi$, supported by TCSA.
- In addition, Takács (arXiv:0801.0962) is by Takács alone (not Kormos–Takács), and its boundary form-factor calculation is checked only to $O(g^2)$. The status recommendation "[Numerical] in sinh-Gordon form; [Conjecture] for sine-Gordon" is accurate and justified.

### Claim 9
**Literature correction: AKR printed $\Delta c_{\rm boun}$ is wrong.**
In AKR (batft.tex lines 972–974), the printed equation is:
$$\Delta c_{\rm boun}(b) = \frac{6 \bar{m} R \cos(\pi/2h)}{\sin(\pi B/2h)\sin(\pi(1-B)/2h)} .$$
Fateev eq. (eQ) has $\sin(\pi/2h)$ in the numerator, and the denominator has $\sin(\pi x/2h)\cos(\pi(1-x)/2h)$ with $x = B/2$. Furthermore, converting boundary energy $\mathcal{E}$ to $\Delta c$ requires a factor of $24/\pi$, so $\Delta c$ lacks a factor of $1/\pi$. At $h=2$, this introduces a discrepancy of $\pi \cot(\pi(1+B)/4)$. Confirmed.

### Claim 10
**Literature correction: printed slips in BRZ and Takács.**
1. In BRZ (bshpfinal.tex eq. Bfugacity / line 680), the numerator factor is printed as $\cosh\theta - \cos(\eta_1 p)$ instead of $\cosh\theta - \cosh(p\eta_1)$. The expression for $R^{(1)}$ (eq. R1 / line 533) clearly has $\cosh(p\eta)$, and setting $\theta' = i\pi/2 - \theta$ in $R(i\pi/2-\theta)R(i\pi/2+\theta)$ yields $\cosh(p\eta)$.
2. In Takács (bffshg.tex eqs. zEFmatch and zetathetamatch), the match $z = \frac{1}{2}(E-F)$ and $z = \frac{\beta^2}{\pi}(\eta-i\vartheta)$ lacks an explicit factor of $i$ on both sides when comparing hyperbolic to trigonometric arguments; the two omissions cancel in the relation between $(E,F)$ and $(\eta,\vartheta)$. Confirmed.

### Claim 11
**Book correction: Chapter 29 self-duality point.**
In `part7-half-line/29-boundary-sine-gordon.qmd`, line 276 reads:
> "So at the self-dual point $B=2$ every reflection factor with $\epsilon<0$ tends to the same function $-1/(1)^2$."
Under $B \leftrightarrow 2-B$, the fixed point is $B = 1$. The point $B=2$ is the infinite-coupling limit $\beta \to \infty$. This is an unambiguous typo in the book text. Confirmed.

### Claim 12
**Proposed replacement text and bibliography.**
The proposed replacements for the "Rank one" paragraph of `sec-uniform-boundary`, open problem 12(a), and `cnj-bsg-uv-ir` status are mathematically accurate, adhere strictly to book style, and provide the exact required resolution.
The proposed bib entry `brz2008` was checked against literature records: Nuclear Physics B, Vol. 796 (2008), pages 622–650. The page number 622 is confirmed. Note that `code/shg_boundary_tba.py` referenced in the proposed text should be created by promoting `briefs/halfline_open_12a_scripts/shg_btba_uniform.py`.

## Errors and overclaims not tied to a single claim
None. The worker explicitly limited scope to the lead and part (a), left part (b) to the assigned Sonnet worker, correctly avoided unnecessary two-loop calculations, and did not edit any repository files.

## Prior art found
- **BRZ (arXiv:0710.4789, NPB 796 (2008) 622–650):** The worker brought in this paper, which was not in the original brief. We searched the full text of BRZ: BRZ studied general boundary parameters $s$ but never mention the Toda uniform boundary, nor did they connect $s = ib/2$ to the Toda uniform boundary or the $E = B/2$ vs $E = 0$ resolution.
- **AKR (hep-th/0110218):** AKR explicitly excluded sinh-Gordon ($A_1$) from their boundary TBA because of its continuous parameters.
- **Synthesis:** The connection identifying the $a_n^{(1)}$ Toda uniform boundary parameter with the sinh-Gordon $s = ib/2$ point and resolving the rank-one puzzle is new to this project.

## Proposed book text
The replacement texts proposed in Section 4 of the worker's report are recommended for adoption without modifications, subject to creating `code/shg_boundary_tba.py` from the worker's script.

## Items for an Opus pass
1. **AKR printed formula slip (Claim 9):** Verify that the printed slip in AKR eq. (972) is strictly a typographical error in the manuscript and did not affect AKR's numerical plots in their Figures 1–3.
2. **BRZ fugacity typographic slip (Claim 10):** Verify the cancellation of double poles at $\theta = 0$ in the crossed channel for symmetric boundaries under the corrected $\cosh(p\eta)$ reading.
3. **Boundary energy identity for higher rank (Claim 6):** Verify whether an analogous algebraic identity between Fateev's boundary energy and the particle reflection factors holds for $a_n^{(1)}$ with $n \ge 2$.

## Scripts
- `/tmp/halfline-referee/12a/referee_independent_check.py`: Independent Python script verifying BRZ Table 1 and solving the sinh-Gordon boundary TBA at $E=B/2$ and $E=0$ against the boundary Liouville reflection amplitude without worker code.

Referee pass complete.
