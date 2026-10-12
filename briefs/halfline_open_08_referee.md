# Referee report: half-line open problem 08
Date: 2026-10-12. Referee: Gemini 3.8 Flash High, via Antigravity. Tier: A.

## Summary
The worker report's main conclusion stands: the unexplained pole $p = \frac{1}{2} + (2n-1)\frac{B}{2}$ of $K_1$ on $b_{-n,n}$ in the window $\frac{1}{2n-1} < B < \frac{2}{2n-1}$ is not a boundary bound state, and the boundary bound-state spectrum of Delius and Gandenberger (1999) stands.
At $B = \frac{9}{20}$, the residue $-i\operatorname{Res} K_1$ at the pole $p = \frac{47}{40}$ on $b_{-2,2}$ is $-1.18749$, which is strictly negative; in a unitary theory with real couplings, an $s$-channel boundary bound-state pole must have $-i\operatorname{Res} K > 0$, ruling out a genuine bound state.
The pole is explained by the boundary Coleman–Thun double-triangle diagram D1 reflecting at a zero left behind by the virtual state $b_{-n+1,n}$, which crosses threshold at $B = \frac{1}{2n-1}$.
Diagram D1 is kinematically realizable with positive edge lengths strictly inside the window $\frac{1}{2n-1} < B < \frac{2}{2n-1}$, and the residue identity $\operatorname{Res} K_1^\beta(iu_1) = \frac{1}{2} g^2 \Gamma^2 K_2^\alpha{}'(\theta_R)$ is exact and proven from the boundary and bulk bootstraps.
However, the worker's claim that the sign discrepancy in their naive loop calculation is explained by "near-wall collapsed configurations" is an unverified hypothesis, requiring editorial refinement of the proposed text for Open Problem 8.

## Verdicts
| # | Claim (short) | Label claimed | Verdict | Label you recommend |
|---|---|---|---|---|
| 1 | Residue signs at $B = 9/20$ | [Numerical] | CONFIRMED | [Numerical] |
| 2 | Scan of $p_1$ across $B$ (sign flip at $B_c = 1/(2n-1)$) | [Numerical] | CONFIRMED | [Numerical] |
| 3 | Threshold mechanism (linear vanishing and virtual state) | [Numerical] | CONFIRMED | [Numerical] |
| 4 | Diagram D1 geometry and realizability in window | [Claim] | CONFIRMED | [Theorem] |
| 5 | Exact-residue identity $\operatorname{Res} = \frac{1}{2} g^2 \Gamma^2 K_2'$ | [Numerical] | CONFIRMED | [Theorem] |
| 6 | Lee–Yang calibration ($|\Gamma g|/2 = 1.2679$) | [Numerical] | CONFIRMED | [Numerical] |
| 7 | Rule at zero incomplete & collapsed configurations | [Claim] | PARTIAL | [Conjecture] |
| 8 | Coincidence at $n=2, B=1/4$ (ratio 2.000000) | [Numerical] | CONFIRMED | [Numerical] |
| 9 | Verdict: DG spectrum stands; no masked state | [Claim] | CONFIRMED | [Claim] |
| 10 | Proposed replacement text for book and open problem 8 | [Text] | CONFIRMED WITH EDITORIAL REFINEMENT | [Text] |

## Details

### Claim 1
**Residue signs at $B = 9/20$ (`p08_signs.py`).**
Confirmed. In the book's conventions (`sec-boundary-bound-states`, eq. (28.14) / GZ eqs. (3.57), (3.62)), an elastic reflection pole corresponding to an $s$-channel boundary bound state satisfies
$$
K(\theta) \simeq \frac{i}{2} \frac{(g^\beta_{a\alpha})^2}{\theta - iv}.
$$
Therefore, $\operatorname{Res}_{\theta=iv} K(\theta) = \frac{i}{2} (g^\beta_{a\alpha})^2$, which implies $-i\operatorname{Res}_{\theta=iv} K = \frac{1}{2} (g^\beta_{a\alpha})^2$. In a unitary theory with real coupling, $(g^\beta_{a\alpha})^2 > 0$, so $-i\operatorname{Res} K > 0$ is strictly required for any physical boundary bound state.
We performed an independent Cauchy contour integration and analytic block factorization of $K_1$ on $b_{-2,2}$ at $p = 47/40$ ($B = 9/20$) to 35 decimal places (`independent_residue_check.py`):
$$
-i\operatorname{Res}_{\theta = 47\pi i/120} K_1 = -1.1874897238607920702207766031584623 < 0.
$$
Simultaneously, every physical creation pole of a bound state in DG's spectrum has $-i\operatorname{Res} K > 0$ (e.g. $b_{0,0} \to b_{1,0}$ at $11/40$ gives $+0.468686$; $b_{0,0} \to b_{-1,1}$ at $51/40$ gives $+2.02180$; $b_{-1,1} \to b_{-2,2}$ at $33/40$ gives $+1.68312$; $b_{-2,2} \to b_{-3,3}$ at $3/8$ gives $+0.917118$). Thus, $p = 47/40$ cannot be a bound state.

### Claim 2
**Scan of $p_1 = \frac{1}{2} + (2n-1)\frac{B}{2}$ over $B$ (`p08_scan.py`).**
Confirmed. On $b_{-n,n}$, $K_1$ factorizes into $K_{\text{uniform}, 1}$ times four blocks:
$$
\Big(\frac{1}{2} + (2n-1)\frac{B}{2}\Big)\Big(\frac{5}{2} + (2n+1)\frac{B}{2}\Big)\Big(\frac{3}{2} - (2n+1)\frac{B}{2}\Big)\Big(\frac{3}{2} - (2n-1)\frac{B}{2}\Big).
$$
The first block has pole $p_1 = \frac{1}{2} + (2n-1)\frac{B}{2}$ and the fourth block has pole $p_4 = \frac{3}{2} - (2n-1)\frac{B}{2}$. Their difference is $p_1 - p_4 = (2n-1)B - 1$.
When $B < \frac{1}{2n-1}$, $p_1 - p_4 < 0$ and the fourth block evaluated at $p_1$ has a negative denominator, yielding $-i\operatorname{Res} > 0$.
At $B = \frac{1}{2n-1}$, the two poles collide at $p = 1$, creating a double pole.
When $\frac{1}{2n-1} < B < \frac{2}{2n-1}$, $p_1 - p_4 > 0$, flipping the sign of the fourth block's denominator; all other factors remain strictly non-zero and maintain their signs.
Consequently, $-i\operatorname{Res}$ is strictly negative throughout the entire window $\frac{1}{2n-1} < B < \frac{2}{2n-1}$.
Checked for $n=2$ ($1/3 < B < 2/3$), $n=3$ ($1/5 < B < 2/5$), and $n=4$ ($1/7 < B < 2/7$) in `referee_suite.py`.

### Claim 3
**Threshold mechanism (`p08_threshold.py`).**
Confirmed. In $K_2^{b_{-n+1,n-1}}$, the block creating $b_{-n+1,n}$ is $x_w = \frac{1}{2} - (2n-1)\frac{B}{2}$.
As $B \to \frac{1}{2n-1}^-$, $x_w \to 0^+$. The residue of this block is $2i\sin(\pi x_w/3) \approx \frac{2\pi i}{3} x_w$, which vanishes linearly as $B \to \frac{1}{2n-1}^-$ (verified in `referee_suite.py`).
For $B > \frac{1}{2n-1}$, $x_w < 0$, and $(x_w) = 1/(-x_w)$ develops a simple zero at $\theta_R = \frac{i\pi}{3}((2n-1)\frac{B}{2} - \frac{1}{2})$.
The energy gap is $e(b_{-n+1,n}) - e(b_{-n+1,n-1}) = m_1 \cos(\pi x_w/3) \to m_1$, exactly reaching the threshold of $b_{-n+1,n-1}$ plus a particle 2 at rest.
At $B = 9/20$ ($n=2$), $e(b_{-1,2}) = 1.21670027\,m_1$, whereas the threshold is $e(b_{-1,1}) + m_1 = 1.23344536\,m_1$. Since the pole is on the unphysical sheet / replaced by a zero on the physical sheet, $b_{-n+1,n}$ is a virtual state.

### Claim 4
**Diagram D1 geometry and realizability in window.**
Confirmed. Initial boundary is $\beta = b_{-n,n}$. It emits particle 1 at rapidity $-iv_{n-1} = -i\frac{\pi}{3}(\frac{3}{2} - (2n-1)\frac{B}{2})$, dropping to $\alpha = b_{-n+1,n-1}$.
Incoming particle 1 at rapidity $iu_1$ fuses with emitted particle 1 into particle 2 with three-point fusing angle $2\pi/3$:
$$
u_1 - (-v) = \frac{2\pi}{3} \implies u_1 = \frac{2\pi}{3} - v = \frac{\pi}{3}\Big(\frac{1}{2} + (2n-1)\frac{B}{2}\Big) = \frac{\pi p_1}{3}.
$$
Particle 2 heads toward the boundary at rapidity
$$
\theta_2 = i(u_1 - \pi/3) = i\Big(\frac{\pi}{3} - v\Big) = \frac{i\pi}{3}\Big((2n-1)\frac{B}{2} - \frac{1}{2}\Big) = \theta_R,
$$
which is precisely the zero of $K_2^\alpha$.
For the Euclidean diagram to have positive line lengths, particle 2 must move toward the wall, requiring $\theta_2/i = \pi/3 - v > 0 \iff B > \frac{1}{2n-1}$.
For incoming particle 1 to lie in the physical strip $0 < u_1 < \pi/2$, we must have $p_1 < 3/2 \iff B < \frac{2}{2n-1}$.
The diagram has $P = 6$ propagators and $L = 2$ loops, giving nominal pole order $P - 2L = 2$; reflection off $\alpha$ at the zero $\theta_R$ reduces the pole order to 1 (Dorey–Tateo–Watts rule).

### Claim 5
**Exact-residue identity $\operatorname{Res} = \frac{1}{2} g^2 \Gamma^2 K_2^\alpha{}'(\theta_R)$ (`p08_identity.py`).**
Confirmed (and analytically proved).
We derived a rigorous analytic proof from the bootstrap axioms:
1. By the boundary bootstrap, $K_1^\beta(\theta) = S_{11}(\theta - iv) K_1^\alpha(\theta) S_{11}(\theta + iv)$. At $\theta = iu_1 = i(2\pi/3 - v)$, the argument $\theta + iv = 2\pi i/3$ is the bulk fusing pole of $S_{11}$, where $S_{11} \simeq \frac{i\Gamma^2}{\theta - iu_1}$. Thus:
$$
\operatorname{Res}_{\theta=iu_1} K_1^\beta(\theta) = i\Gamma^2 S_{11}(i(u_1 - v)) K_1^\alpha(iu_1).
$$
2. By the bulk bootstrap for $K_2^\alpha$: $K_2^\alpha(\theta) = K_1^\alpha(\theta + i\pi/3) S_{11}(2\theta) K_1^\alpha(\theta - i\pi/3)$.
At $\theta = \theta_R = i(\pi/3 - v)$, $\theta_R + i\pi/3 = iu_1$ and $2\theta_R = i(u_1 - v)$.
The third factor is $K_1^\alpha(\theta_R - i\pi/3) = K_1^\alpha(-iv)$. Since $K_1^\alpha$ has a pole at $+iv$ with residue $\frac{i}{2}g^2$, boundary unitarity $K_1(-\theta) = 1/K_1(\theta)$ forces a simple zero at $-iv$ with derivative $\left.\frac{d}{d\theta} K_1^\alpha(\theta)\right|_{-iv} = \frac{2i}{g^2}$.
Differentiating $K_2^\alpha$ at $\theta_R$:
$$
K_2^\alpha{}'(\theta_R) = K_1^\alpha(iu_1) S_{11}(i(u_1 - v)) \frac{2i}{g^2}.
$$
3. Multiplying by $\frac{1}{2} g^2 \Gamma^2$ immediately yields:
$$
\frac{1}{2} g^2 \Gamma^2 K_2^\alpha{}'(\theta_R) = i\Gamma^2 K_1^\alpha(iu_1) S_{11}(i(u_1 - v)) = \operatorname{Res}_{\theta=iu_1} K_1^\beta(\theta).
$$
We verified this identity independently to 50 decimal places in `verify_claim5.py` at $n=2, B=9/20$:
$\operatorname{Res} K_1^\beta = -1.1874897238607920702207766031584623408400596567185 i$, matching $\frac{1}{2} g^2 \Gamma^2 K_2^\alpha{}'(\theta_R)$ with difference $< 2 \times 10^{-50} i$.

### Claim 6
**Lee–Yang calibration (`p08_calib_ly.py`).**
Confirmed. For boundary Lee–Yang ($h=6$), $K_{(1)} = (1)(3)(4)^{-1}$.
The one-point pole at $\theta = i\pi/6$ has residue $|\operatorname{Res} K| = 1.26794919243$.
The double-reflection pole at $\theta = i\pi/2$ has residue $-0.92820323 i$, giving $g^2 = \operatorname{Res}/(i/2) = -1.85640646$.
The bulk S-matrix $S(\theta) = \frac{\sinh\theta + i\sin(2\pi/3)}{\sinh\theta - i\sin(2\pi/3)}$ has residue $i\Gamma^2$ at $2\pi i/3$ with $\Gamma^2 = -2\sqrt{3} = -3.4641016$.
Because $1 - S(i\pi) = 2$, the $s$- and $u$-channel trees add at $\theta = i\pi/2$, halving the effective coupling $g_s^2 = g^2/2$.
Then $|\Gamma g|/2 = \frac{1}{2}\sqrt{|g^2 \Gamma^2|} = 1.26794919243$, matching GZ eq. (3.41).

### Claim 7
**The naive rule "residue picks up $K'$ at the zero" is incomplete.**
Verdict: PARTIAL (Conjecture).
The worker evaluated the asymptotic part of diagram D1 using mixed position/momentum-space propagators and obtained $-\frac{1}{2} g^2 \Gamma^2 K_2'$, which has the opposite sign to the exact residue $+\frac{1}{2} g^2 \Gamma^2 K_2'$.
To reconcile this, the worker conjectured that "near-wall (collapsed) configurations contribute at the same order through a local boundary amplitude."
However:
1. The worker did not compute the collapsed configuration contribution.
2. The worker did not check the overall sign of their loop evaluation against any independently known boundary loop residue.
3. In 1+1D boundary loop integrations, overall signs are sensitive to loop orientation, vertex conventions, and $i\epsilon$ prescriptions. A simple minus sign in the loop integral definition would make the asymptotic contribution alone $+\frac{1}{2} g^2 \Gamma^2 K_2'$.
Therefore, while the worker correctly notes that the naive rule has not been established from first-principles Feynman diagrams, the assertion that "collapsed configurations contribute a local amplitude that flips the sign" is an unverified hypothesis rather than a demonstrated mechanism.

### Claim 8
**The coincidence at $n=2$, $B=1/4$.**
Confirmed. At $B = 1/4$, the D1 pole $p_1 = 1/2 + 3B/2 = 7/8$ coincides with the creation pole $p = 3/2 - 5B/2 = 7/8$ of $b_{-3,3}$.
Simultaneously, the block $(-1 + B/2)$ in $K_{\text{uniform}, 1}$ evaluated at $B=1/4$ is $(-7/8)$, which possesses a simple zero at $\theta = 7\pi i/24$.
This zero cancels one of the two colliding poles, keeping the pole simple (order $2 - 1 = 1$).
Its exact residue is $\operatorname{Res} = +5.8613740399 i$, which is exactly twice the D1 formula: $\frac{1}{2} g^2 \Gamma^2 K_2^\alpha{}'(\theta_R) = +2.9306870200 i$ (ratio $2.000000000000$, verified to 28 digits in `independent_checks.py`).
This factor of 2 reflects the degeneracy of the two coincident diagrams.

### Claim 9
**Verdict: DG's spectrum stands; no masked bound state found.**
Confirmed. Three independent lines of evidence exclude a bound state at $p_1$:
1. Unitarity sign violation: $-i\operatorname{Res} K_1 < 0$ inside the window.
2. Threshold physics: $b_{-n+1,n}$ dissolves into a virtual state at $B = \frac{1}{2n-1}$ with energy reaching threshold.
3. Kinematic realizability: diagram D1 exists and closes with positive lengths precisely in the window.
The formal loophole noted in report §4 (a hypothetical bound state coupling only to particle 1 and canceling a CT diagram of opposite sign) is ruled out by continuity and conserved charges.

### Claim 10
**Proposed replacement text.**
Confirmed with editorial refinement (see Section below).

---

## Errors and overclaims not tied to a single claim
1. **Status of Claim 5:** The worker labeled the residue identity $\operatorname{Res} K_1^\beta = \frac{1}{2} g^2 \Gamma^2 K_2^\alpha{}'(\theta_R)$ as merely `[Numerical]`. It is an exact analytic theorem proven directly from combining the boundary and bulk bootstrap equations.
2. **Loop evaluation vs. collapsed configurations:** The worker presented the hypothesis of near-wall collapsed configurations as an established physical mechanism explaining why their loop calculation differed by a sign from the exact residue. Without an independent computation of the contact amplitude or a verification of the loop sign, this remains purely speculative.

---

## Prior art found
- **Dorey, Tateo, and Watts** (`dtw1999`, hep-th/9810098): Applied power counting $P - 2L - N_{\text{zeros}}$ to boundary Lee–Yang, establishing that reflections at zeros reduce the singularity order. They did not evaluate residues or loop amplitudes.
- **Mattsson and Dorey** (`mattsson_dorey2000`, hep-th/0008071) and **Bajnok, Palla, Takács, and Tóth** (`bptt2002`, hep-th/0106070): Accounted for the boundary sine-Gordon spectrum using Coleman–Thun diagrams and cancellations between diagrams, without computing residues of zero-reduced diagrams.
- **Bajnok, Böhm, and Takács** (hep-th/0309119, NPB 682 (2004) 585): Developed the boundary reduction formula, boundary Landau equations, and Cutkosky rules on the half-line, but did not compute residues of diagrams reflecting at zeros or evaluate collapsed configurations.

---

## Proposed book text
The proposed replacement for `sec-boundary-spectrum-a2` (report §5, first block) is accurate, rigorous, and ready for adoption:
> In the window $\frac1{2n-1}<B<\frac2{2n-1}$ the pole $p=\frac12+(2n-1)\frac B2$ of $K_1$ on $b_{-n,n}$ is not a bound state. At $B=\frac1{2n-1}$ the state $b_{-n+1,n}$ reaches the threshold of $b_{-n+1,n-1}$ plus a particle 2 at rest. Its creation pole in $K_2^{b_{-n+1,n-1}}$, at $p=\frac12-(2n-1)\frac B2$, passes through $\theta=0$ with vanishing residue, and above the threshold the strip contains instead a zero at $\theta_R=\frac{i\pi}3\big((2n-1)\frac B2-\frac12\big)$, a virtual state. Exactly in the window there is a boundary Coleman–Thun diagram. $b_{-n,n}$ emits a particle 1 and drops to $b_{-n+1,n-1}$. The incoming particle 1 fuses with it into a particle 2, which reflects at $\theta_R$, the zero, and splits back, and the emitted particle is reabsorbed. The diagram has two loops and would give a double pole; the zero makes it simple. The residue agrees: $-i\operatorname{Res}$ is positive below the window, changes sign at $B=\frac1{2n-1}$, where the pole meets the pole at $\frac32-(2n-1)\frac B2$, and is negative in the window. At $B=\frac9{20}$ the pole $\frac{47}{40}$ on $b_{-2,2}$ has $-i\operatorname{Res}K_1=-1.187$, impossible for a bound state in a unitary theory. The spectrum of Delius and Gandenberger therefore stands.

For Open Problem 8, however, the report's text overclaims by stating that collapsed configurations have been shown to contribute through a local boundary amplitude. We recommend replacing it with the following refined text:
> **Masked bound states.** In the bulk of non-simply-laced theories Coleman–Thun poles can mask genuine bound states [DG sec. 7]. On the half-line the leading coefficient of an on-shell diagram is a product of couplings $g$, $\Gamma$, reflection and scattering factors and a geometric factor; this reproduces Ghoshal and Zamolodchikov's one-point residue in boundary Lee–Yang. For diagrams whose order is reduced by a reflection at a zero, such as diagram D1 in @sec-boundary-spectrum-a2, the boundary bootstrap fixes the exact residue to $\frac{1}{2}g^2\Gamma^2K_2'(\theta_R)$. For $a_2^{(1)}$, residue signs ($-i\operatorname{Res} < 0$), kinematic realizability, and threshold behaviour exclude a masked bound state, confirming the spectrum of Delius and Gandenberger. A general field-theoretic residue formula for boundary on-shell diagrams with zero-reflections, the resolution of coincident diagrams at $a_3^{(1)}$ ($b=2$) and $a_4^{(1)}$ ($n_2=-n_3$), and the unexplained $\nu_{-j}$ poles of boundary sine-Gordon remain open.

---

## Items for an Opus pass
none

---

## Scripts
All written scripts are in `/tmp/halfline-referee/08/`:
- `independent_residue_check.py`: Independent Cauchy contour integration and analytic block factorization of $K_1^{b_{-2,2}}$ at $p = 47/40, B = 9/20$ to 35 dps.
- `verify_claim5.py`: Independent calculation and verification of the exact residue identity $\operatorname{Res} K_1^\beta = \frac{1}{2} g^2 \Gamma^2 K_2^\alpha{}'(\theta_R)$ to 50 dps.
- `referee_suite.py`: Comprehensive test suite verifying signs at $B=9/20$, sign scans across $B$, threshold vanishing, Lee–Yang calibration, and coincidence ratio.

Referee pass complete.
