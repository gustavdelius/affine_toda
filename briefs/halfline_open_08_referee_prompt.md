# Referee prompt: half-line open problem 08
Tier: A
Brief: briefs/halfline_open_08_masked_bound_states.md
Report: briefs/halfline_open_08_report.md
Scripts: briefs/halfline_open_08_scripts/ (self-contained, mpmath only; each runs in seconds): `p08_common.py`, `p08_signs.py`, `p08_scan.py`, `p08_identity.py`, `p08_calib_ly.py`, `p08_threshold.py`.

Scope:
- **Done:** plan step 3 (the $a_2$ unexplained pole) and the leading-order residue formula with a Lee–Yang calibration.
- **Partial:** step 2 (the rule at a zero) and step 4 (one coincidence).
- **Not done:** the sine-Gordon calibration and step 5 beyond $a_2$.

Tier A: no theorem is claimed and no published paper is corrected. But the result would remove an "unexplained pole" from the book (`sec-boundary-spectrum-a2`), so check it carefully.

## Claims to referee

1. [Numerical] **Residue signs at $B=\frac9{20}$** (`p08_signs.py`).
   - $-i\operatorname{Res}K$ is positive at every creation pole of a state in DG's spectrum.
   - The pole $\frac{47}{40}$ of $K_1$ on $b_{-2,2}$ has $-i\operatorname{Res}=-1.18749$.
2. [Numerical] **Scan of the pole $p_1=\frac12+(2n-1)\frac B2$ over $B$** (`p08_scan.py`). The residue is positive below the window and changes sign exactly at $B=\frac1{2n-1}$, where $p_1$ meets the crossed pole $\frac32-(2n-1)\frac B2$ and the pole is double. It is negative inside the window $\frac1{2n-1}<B<\frac2{2n-1}$, for $n=2,3,4$ (table in report §2).
3. [Numerical] **The threshold mechanism** (`p08_threshold.py`).
   - The creation pole of $b_{-n+1,n}$ in $K_2^{b_{-n+1,n-1}}$, at $p=\frac12-(2n-1)\frac B2$, has a residue that vanishes linearly as $B\to\frac1{2n-1}$ from below.
   - Above that value the strip has a zero at $\theta_R=\frac{i\pi}3((2n-1)\frac B2-\frac12)$, so $b_{-n+1,n}$ becomes a virtual state. At $B=\frac9{20}$ its energy is $1.2167\,m_1$, against a threshold of $1.2334\,m_1$.
4. [Claim] **Diagram D1 is realizable exactly in the window.** The particle-2 angle $\frac\pi3-v>0$ iff $B>\frac1{2n-1}$, and $u_1<\frac\pi2$ gives $B<\frac2{2n-1}$. No internal lines cross, and the pole order 2 is lowered to 1 by the zero.
5. [Numerical] **The exact-residue identity.** $\operatorname{Res}K_1^{b_{-n,n}}(iu_1)=\frac12g^2\Gamma^2K_2^{\alpha}{}'(\theta_R)$, with $\alpha=b_{-n+1,n-1}$ and $g$, $\Gamma$ defined by the residues of $K$ and $S_{11}$. It holds to working precision at six $(n,B)$ points (`p08_identity.py`).
6. [Numerical] **The Lee–Yang calibration** (`p08_calib_ly.py`). With the boundary coupling halved at $\theta=\frac{i\pi}2$ (s- and u-channel trees add, since $1-S(i\pi)=2$), the one-point pole residue is $|\Gamma g|/2=1.26794919$, which is Ghoshal–Zamolodchikov (3.41). Only the magnitude is calibrated.
7. [Claim] **The naive rule "residue picks up $K'$ at the zero" is incomplete.**
   - The worker's own evaluation of the asymptotic part of D1 gives $-\frac12g^2\Gamma^2K'$, the opposite sign to the exact residue.
   - The worker attributes the difference to near-wall (collapsed) configurations that contribute at the same order.
   - The overall sign of the loop evaluation is not independently cross-checked.
8. [Numerical] **The coincidence at $n=2$, $B=\frac14$.** The pole stays simple, with residue twice the D1 structure (ratio 2.000000).
9. [Claim] **Verdict.** DG's spectrum stands; no masked bound state was found. One formal loophole remains (report §4).
10. [Text] **The proposed replacement text** for the paragraph "This range leaves one pole unexplained" in `sec-boundary-spectrum-a2` and for open problem 8.

## Coordinator notes (suspected weak points)

- **Claim 1's sign convention is the foundation of the verdict.** Check the claim "$-i\operatorname{Res}<0$ is impossible for a bound state in a unitary theory" in the book's conventions (`sec-boundary-bound-states`, ch. 28): which sign of $-i\operatorname{Res}K$ corresponds to positive $(g^\beta_{a\alpha})^2$ at real coupling? Recompute $-i\operatorname{Res}K_1$ at $\frac{47}{40}$ on $b_{-2,2}$ yourself, from the book's `code/particle_reflection.py` (`K1nm`, `Kex`) or from DG (6.7) directly, without the worker's code.
- **Claim 7 is self-contradictory as it stands.** If the naive evaluation has the wrong sign, the "D1 explains the pole" story rests on the exact identity (claim 5) plus geometry (claim 4), not on a computed diagram. Decide whether the proposed book text claims more than this. In particular, "the residue agrees" should mean only the sign-change and threshold pattern.
- **Claim 5.** Verify the identity independently at one point, e.g. $n=2$, $B=\frac9{20}$, computing $g^2$, $\Gamma^2$ and $K_2'(\theta_R)$ yourself.
- **Prior art.** Check Dorey–Tateo–Watts (`dtw1999`), Mattsson–Dorey, and Bajnok–Böhm–Takács (hep-th/0309119) for an existing statement of the boundary residue rule at a zero, or of the collapsed-configuration contribution.
