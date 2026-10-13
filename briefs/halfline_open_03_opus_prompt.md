# Opus referee prompt: half-line open problem 03 (tier-B follow-up)

Written 2026-10-12 by the Claude coordinator, after reading the Gemini verdict `briefs/halfline_open_03_referee.md`.

**What the Gemini pass established.** It independently established the following, so do not redo them:
- the proof of $U_n$ (claim 3, by sympy);
- that the bootstrap fails without the CDD factor $(\frac T2\mp\frac\omega2)_I$ at $n=2$ (claim 8, with its own test);
- the origin of the $-i$ in $A^{(-)}$: Gandenberger's explicit $i$ in $\kappa_-$, `paprev.tex` around line 720.

**Gemini's own Opus items** (the physical-strip pole, the odd-$n$ mechanism, an analytic proof of $\Phi$) are research problems and go to the book's open-problem list. They are not referee items.

Read:
- the brief `briefs/halfline_open_03_neumann_soliton_route.md`;
- the report `briefs/halfline_open_03_report.md`;
- the Gemini verdict;
- the scripts in `briefs/halfline_open_03_scripts/`. They import unchanged copies of `textbook_code/*.py`, so copy those alongside.

## Items

1. **Claim 1 against the book.**
   - Find what the book states about boundary unitarity and crossing-unitarity for a non-diagonal $K$ (`sec-boundary-unitarity-crossing` in ch. 28, `sec-an-nondiagonal-k` in ch. 33, and any use in ch. 34–35), and decide whether it is wrong, incomplete, or already consistent with $A(\theta)A(-\theta)K(-\theta)^{\rm T}K(\theta)=1$.
   - Derive analytically, in the book's conventions (physical $q=-e^{-i\pi\omega}$, the book's coproduct, the book's antiparticle representation $\bar V$), that the antisoliton is reflected by $JKJ$ and that $JKJ=K^{\rm T}$ for (`eq-an-k-matrix`). Is $JKJ=K^{\rm T}$ an identity (persymmetry) or does it hold only on the $\hat\epsilon_j=\hat\epsilon_{n+1-j}$ locus?
   - Gemini checked this only numerically and did not look at the book.
2. **Claim 8, Kim's one loop.** Check against `kim1995` (eqs. (33), (39); fetch the paper) that Kim's $a_3^{(1)}$ Neumann one-loop amplitudes agree with $K^N_a$ including the factor $(\frac T2\mp\frac\omega2)_I$, and disagree without it (the $O(B)$ term of DG's $k_{h/2-1}$). This is the only input that fixes the odd-$n$ CDD choice at $n=3$, and Gemini restated it without showing the check.
3. **Claim 9, independent evaluation.** Evaluate the breather identity $A(\mu_a)F_{11}(2\mu)A(\mu_b)c_+(\mu)=(-\frac{n\omega}2)_I(\frac{n\omega}2-\frac12)_I(\frac12)_I$ at one generic point for $n=3$ and one for $n=4$ to at least 12 digits. Do not use the worker's truncated products: use an integral representation of the Gamma-product blocks, or mpmath with convergence acceleration. This raises the evidence from about $10^{-7}$ to a level the book can quote.
4. **Prior art.** Search seriously (arXiv full-text, INSPIRE citations of Gandenberger hep-th/9806003 and hep-th/9911178, of Delius–Gandenberger hep-th/9904002, and of Ahn–Kim–Rim hep-th/0110218) for any published scalar factor of the non-diagonal $a_n^{(1)}$ K-matrix, or any soliton-side derivation of the Neumann amplitudes for $n\ge3$. Candidates include Doikou, Nepomechie, Avan–Doikou, Kim, Ahn, and Bowcock–Umezawa.

For each item give the verdict (CONFIRMED / CORRECTED / REFUTED / UNCHECKED), the evidence, and any corrected book text. Give script paths.
