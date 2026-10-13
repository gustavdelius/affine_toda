# Referee prompt: half-line open problem 07
Tier: B
Brief: briefs/halfline_open_07_nonunitary_half_line.md (parts (a) and (b) only)
Report: briefs/halfline_open_07_report.md
Scripts: briefs/halfline_open_07_scripts/
- **Setup.** The scripts import `boundary_common.py`, an unchanged copy of `textbook_code/boundary_common.py`; copy it alongside.
- **Part (a):** `a_krein_algebra.py`, `a_crossing_vw.py`, `a_dn.py`.
- **Part (b):** `b_transfer.py`, `b_checks.py`, `b_scan.py`, `b_bisect*.py`, `b_by_n1.py`. The `*.out` files are the worker's outputs.

Scope: parts (a) and (b) for the $a_n^{(1)}$ vector solitons, with side checks for sine-Gordon and $d_4$. Part (c) was deferred. Tier B because a theorem is claimed (claim 1), together with corrections to the brief's premises and to the book (claims 6, 9).

**Context.** The non-diagonal $s=\pm1$ scalar factors used here come from problem 3 (`briefs/halfline_open_03_report.md`, refereed in `briefs/halfline_open_03_referee.md`). The worker used only the modulus of their prefactor.

## Claims to referee

1. [Theorem] **Boundary Krein-unitarity for $a_n^{(1)}$** (report §1, parts 1–6). Take the antilinear anti-automorphism $\tau(e_j)=e_{\nu j}$, $\tau(f_j)=f_{\nu j}$, $\tau(k_j)=k_{\nu j}$, where $\nu$ is $j\mapsto n+1-j$, and $C$ the basis reversal. Then:
   - (1) $\pi(a)^\dagger=C\pi(\tau a)C$ on $V$ and $\bar V$.
   - (2) Every $\check R_{ab}$ satisfies $\check R^\dagger(C\otimes C)\check R=C\otimes C$, including $V\bar V$.
   - (3) $B_{\hat\epsilon}$ is $\tau$-stable iff $\hat\epsilon_{\nu j}=\overline{\hat\epsilon_j}$, so $K^\dagger CK=c_KC$ with $c_K$ real.
   - (4) $c_K=1$ for $\hat\epsilon=0$, and $c_K=(\cosh2T\theta-s\cos\pi T)/(2\sin^2\frac{\pi\omega}2)\ge0$ for $s=\pm1$.
   - (5) $\mathbb K_V^\dagger C\mathbb K_V=C$ iff $\overline{A(\theta)}=sA(-\bar\theta)$.
   - (6) The double-row transfer matrices are Krein-unitary.
2. [Claim] **The scalar factors have the required reality.** DG's $A_d$ and the problem-3 $A^{(\pm)}$ satisfy $\overline{A(\theta)}=sA(-\bar\theta)$. Only the prefactor $A_{\min}$ was checked, not the full infinite products.
3. [Claim] **The left-wall convention is forced.** The left wall is the image of the right wall under $e\leftrightarrow f$, $k\to k^{-1}$, so $K^L_{\bar a}(-\theta)=K^R_a(\theta)$. With it $[\mathcal T_i,\mathcal T_j]=0$; exchanging $K$ and $K^{\rm T}$ on the left breaks commutativity whenever the left wall is non-diagonal.
4. [Numerical] **Sine-Gordon.** $K^\dagger\sigma^xK\propto\sigma^x$ iff $\epsilon_0,\epsilon_1$ are real, i.e. $\phi_0\in\frac{2\pi}\beta\mathbb Z$. The proportionality constant can be negative when $\epsilon_0\epsilon_1<0$ and $|\epsilon|$ is large.
5. [Numerical] **$d_4$.** The Krein form is a signed anti-diagonal. All $K$ of `prp-dn-k` satisfy the relation, and patterns with $\epsilon_{n-1}\ne\epsilon_n$ fail under the swapped $\nu$.
6. [Correction to book] **Bulk.** The soliton–antisoliton amplitudes are Krein-unitary with the *untwisted* $C$, which settles the remark after `prp-krein-fundamental` and the crossed case of outlook problem 5.
7. [Theorem / Numerical] **Part (b) symmetry.** Every $\mathcal T_j$ commutes with $\mathcal Z=Z^{\otimes N}$, and the Krein form pairs sector $\zeta$ with $\bar\zeta$, so only $\zeta=\pm1$ sectors can be unbroken.
   - For $a_2$, one soliton between U and N walls, or between N and N, has complex energies in $\zeta=e^{\pm2\pi i/3}$. This is in closed form: $|\Lambda|=\rho^{\pm1}$ or $\rho^{\pm2}$, and the Bethe–Yang energies at $mL=10$ are tabulated.
8. [Numerical] **Part (b) tables 1 and 2.**
   - In $\zeta=1$, UU two-soliton states (VV and VW) are unbroken at all sampled $\omega$, although VW is broken in the bulk.
   - With N walls there are thresholds: UN VV breaks at $\omega\approx0.69$, NN VV at $\omega\approx0.535$.
   - Three-soliton states are broken.
9. [Correction to brief and book] **The premise of (b) is false as stated.** "Real for $\Theta$-invariant boundaries" fails. The relevant condition is $\Theta\nu$-invariance, and even then reality holds only in self-conjugate sectors. `prp-boundary-theta` has no $\nu$.
10. [Text] **The proposed text for open problem 7** and the related edits (`sec-krein-unitarity`, `eq-boundary-unitarity`), including a new `prp-boundary-krein`.

## Coordinator notes (suspected weak points)

- **Whether the $\zeta\ne1$ sectors are physical is the central interpretive question** (claims 7, 8). "A single soliton between a uniform and a Neumann wall has complex energies" is a striking physical statement. It holds only if the $\mathbb Z_{n+1}$ sectors with $\zeta\ne1$ are states of the theory on the interval. The worker flags this as unsettled: for Neumann walls they should be Fourier modes of degenerate boundary vacua. Consider topological charge on the interval: which soliton charges can interpolate between the vacua allowed at the two walls, and does that select the sectors? Decide whether the proposed text states this caveat strongly enough.
- **Claim 1 (2), $V\bar V$.** The proof gets $c=1$ from "all transmissions equal the crossed $VV$ transmission". Check that crossing step analytically, not only by the numerical spectrum comparison. Check also that the book's `prp-krein-fundamental` uses the same $C$ (basis reversal) in the same gradation, so that "untwisted $C$" in claim 6 matches the book's statement.
- **Claim 1 (3).** Check the computation of $\tau(b_j)$, and the condition $c_{\nu j}=q^{-2}\bar c_j$ on the coideal's $c_j$, against the book's coideal generators (ch. 32, `sec-epsilon-constraint`). The worker's "parity-symmetric choice $c_j=q^{-1}$" must be the book's.
- **Claim 2.** Verify, at least numerically, that the full problem-3 $A^{(+)}$ (with its infinite products) satisfies $\overline{A(\theta)}=A(-\theta)$ at real $\theta$.
- **Claim 3.** Compare the left-wall construction with the book's double-row transfer matrix (`sec-boundary-bethe-yang`, ch. 28). If the book uses the "swapped" convention, that is a book correction; say so explicitly.
- **Claim 8.** Recompute one threshold yourself, with your own transfer-matrix code: for example UN VV at $\omega=0.60$ (unbroken) and $\omega=0.80$ (broken), on the worker's rapidity grid.
- **Prior art.** Check Takács–Watts (hep-th/9810006, hep-th/0203073) and the $\mathcal{PT}$-symmetric boundary spin-chain papers (arXiv:2301.06004, 2406.10334) for an existing boundary pseudo-unitarity statement, or for sector-dependent reality breaking.
