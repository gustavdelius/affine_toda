# Brief: referee the proofs of the heat-trace formula and the one-loop mass formula

You are an adversarial mathematical-physics referee, expert in the spectral and scattering theory of non-self-adjoint Schrödinger operators. The text to referee is on branch `worktree-heat-trace-theorems` of the repository `/home/gustav/Git/affine_toda`, checked out at `/home/gustav/Git/affine_toda/.claude/worktrees/heat-trace-theorems`. If that branch has been merged, the same text is on `main`. Do not edit any file in the repository. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files.

## Context

Two results that the one-loop soliton masses of Part V rest on were sketches. They were turned into theorems on 2026-10-10 (brief `briefs/rigorous_heat_trace.md`), and every [Theorem] label below is provisional until it survives this pass. Experience in this repository is that every adversarial referee pass finds real errors, including in proofs that passed numerical checks. The numerics below test the formulas, not the logic of the proofs. So read the proofs line by line.

Find any object with `grep -rn '{#<id>}' part*/`. Run `git -C <worktree> diff main` to see every change. The changes are:

| Id | File | What changed |
|---|---|---|
| `eq-one-loop-mellin` | `part5-semiclassics/17-thimbles.qmd` | new: the Mellin form is made the definition of `eq-one-loop-mass` |
| `prp-net-count` parts 2 and 3 | `part5-semiclassics/18-jordan-chains.qmd` | part 2 was a sketch; part 3 (Pöschl–Teller decomposition, $\sigma(\mathcal A)\subset[0,\infty)$) is new; both now [Theorem], with full proofs |
| paragraphs after those proofs | same | "Numerical checks", "Why no special treatment…", "Which solitons satisfy the hypotheses", "Consequence for (`eq-one-loop-mass`)" |
| `lem-born` | `part5-semiclassics/19-one-loop-masses.qmd` | new lemma: Born asymptotics of $X_b$ |
| `prp-one-loop-mass` | same | was a sketch, now [Theorem] with a full proof; three remarks follow |
| other | chapters 17 and 19, `references.bib` | the $\operatorname{Im}\xi$ sentence in `sec-which-saddles-contribute`; the label of `cnj-counting-rule`; the calibration paragraph and the notes of chapter 19; new entries `simon2005`, `rebhan1997` |

Read first, for the setting: `sec-counting-rule-transmission-factors` (the channel conventions, Hypothesis (R), `lem-wronskian`, `prp-net-count` part 1), Hypothesis (H) and `thm-hirota-factorization`, `prp-hirota-form`, and `sec-counting-rule` (chapter 17). Conventions: the channel-$c$ fluctuation is $e^{ikx}e_c$ at $-\infty$ and $X_c(k)e^{ikx}e_c$ at $+\infty$; $\kappa_c=\sqrt{m_c^2-\lambda}$ with $\operatorname{Re}\kappa_c\ge0$, $\kappa=ik$, $z=\kappa/m$; bound states are zeros with $\operatorname{Im}k<0$, that is $\operatorname{Re}z>0$; $e_b^{\mathsf T}e_c=\delta_{c\bar b}$; $R(\lambda)=(\mathcal A-\lambda)^{-1}$.

## What to attack, in priority order

1. **`prp-net-count` part 2, proof.**
   - Step 1: the numerical-range bound; m-sectoriality; the sign and orientation of the Dunford integral with $R=(\mathcal A-\lambda)^{-1}$; moving Kato's sector contour to the half-strip boundary $\Gamma$.
   - Step 2: enlarging $\Gamma$ to enclose zeros and poles of individual $X_b$ that are not eigenvalues. Does this change the integral?
   - Step 3: are the Hilbert–Schmidt bounds uniform along $\operatorname{Im}\lambda=\pm c$ as $\operatorname{Re}\lambda\to\infty$? Check the identity for $RV_1$, the trace-norm convergence, and the use of part 1 on $\Gamma$, which needs $I+K$ invertible there.
   - Steps 5–7: the conformal map and the preimage $\tilde\Gamma$; the residues; the vanishing connecting segments; the absence of poles on the real axis, including at $k=0$; evenness from (X2).
2. **`lem-born`.**
   - Step 2 identifies the Volterra solution with the (H) solution $\psi_b^-$ at real $k$. Is $\psi_b^-$ defined there, on the boundary of the physical sheet? Is the uniqueness argument for open channels correct when masses are degenerate?
   - Is "(H) with $\psi_b^\mp$ rational in $z_b$" actually established for the $c_n^{(1)}$ solitons and for the foldings, or only asserted? The text cites `sec-example-an1-zeros-alone`, `sec-example-cn1-jordan-blocks` and step 1 of the proof of `prp-hirota-form`.
   - Check the sign and normalization of $c_b$ and the claim that the expansion holds in all directions.
3. **`prp-one-loop-mass`, proof.**
   - Step 1: is $X_\mu(\infty)=1$ properly justified?
   - Step 2: is $(\mathcal A-\mathcal A_0)e^{-t\mathcal A_0}$ trace class, and is its trace as stated?
   - Steps 3–5: the reduction to one root; Tonelli and Frullani; the elementary integral $\int_0^1\log(\cos^2\theta+u^2\sin^2\theta)\,du=2\theta\cot\theta-2$; the case $w=1$.
   - Step 6: the passage between physical-sheet roots and all roots.
4. **`prp-net-count` part 3.** Check the product formula and its sign ambiguity, the erfc integral, the Pöschl–Teller identification, and $\sigma(\mathcal A)\subset[0,\infty)$. For the last, are the analytic Fredholm theorem and part 1 enough to conclude that every non-real or negative spectral point is a zero of $d$?
5. **The definition.** Is `eq-one-loop-mellin` the right reading of `eq-one-loop-mass`, and does it agree with the standard heat-kernel and DHN definitions with normal ordering? The one-loop masses are claimed to agree with the lattice evaluation (Appendix C) and with exact S-matrices. The equivalence of the Mellin form with the continuum limit of the lattice/thimble Gaussian (`prp-thimble-gaussian`) is not proved. Is that gap stated honestly where it matters, and does any [Theorem] label or "closes H8" claim depend on it?
6. **Applicability and status labels.** Does the paragraph "Which solitons satisfy the hypotheses" overclaim? For $d_n^{(1)}$, $e_n^{(1)}$ and the foldings, both (D) and the root location (`prp-hirota-form` part 3) are numerical. Check that the theorems are not presented as unconditional there. Some folded families in the results table of `sec-one-loop-masses-folded` are not among the 19 foldings checked in `sec-hirota-form-dn1-en1`: $b_7^{(1)}$, $d_{n+1}^{(2)}$ with $n=5,6$, and $a_7^{(2)}$. For them the folded transmission factors are assumed to be products of parent factors, $X_b^{\rm fold}=\prod_{a'\in O(a)}X_{a'b}$. Is this stated where the theorems are applied? Check the $\operatorname{Im}\xi$ claim in `sec-which-saddles-contribute`: the theorem covers regular sectors only. Check the new label of `cnj-counting-rule`.
7. **Remarks and prior art.** Check the phase-shift remark (the term $\frac1{4\pi}\int\operatorname{tr}(M-M_0)\,dx$) and the unpaired bookkeeping. Is `rebhan1997` an appropriate citation for it? Is the heat trace of reflectionless potentials, or the Pöschl–Teller decomposition, already in the literature? Look at trace formulas for reflectionless potentials and the spectral methods of Graham, Jaffe and Weigel. If so, the novelty sentence in the notes of chapter 19 must credit it.

## Known pitfalls in this repository

- A status label and a cross-reference in one bracket, such as `[Numerical; @sec-x]`, is parsed by Pandoc as a citation. Write `[Numerical] (@sec-x)`.
- Grouping channels into mass levels is essential: a single channel can have $X_b(0)=0$, so its heat-trace integrand has a $1/k$ singularity.
- The closed form holds only for mass levels, or for channels paired with their conjugates. Unpaired, a root contributes $\tfrac12m_b[\psi(\arcsin w)+\tfrac12\sqrt{1-w^2}]$.

## Required checks

Copy scripts to your scratch directory before running them.

- Rerun `foundations_code/oneloop/heat_trace_checks.py` (copy `oneloop/` and `hirota/` side by side), and compare with `heat_trace_checks.log`. Run `--grid 768` for at least one sector and compare with `heat_trace_grid.log`. At $N=384$ the heavier sectors are off by $10^{-4}$ to $10^{-2}$; this was shown to be discretization error.
- Verify independently, by hand or with a CAS: the erfc integral of part 3; Frullani; the $\log$ integral of step 5; and $\int_0^\infty\frac{dk}{\sqrt{1+k^2}(k^2+w^2)}=\frac{\arccos w}{w\sqrt{1-w^2}}$.
- Test the theorems on a case not in the script: the $\phi^4$ kink ($X=\frac{z-1}{z+1}\cdot\frac{z-1/2}{z+1/2}$), whose one-loop mass is $m\big(\frac1{4\sqrt3}-\frac3{2\pi}\big)$. Compute its heat trace directly, by a grid or the Gel'fand–Yaglom method, and compare with part 3.
- Check that the book renders without unresolved references (`?@`). Render into your scratch directory with `quarto render --output-dir <scratch>/site`; the two warnings about the path configuration are harmless.

## Deliverable

In your final message:

- every problem found, ranked: **fatal** (the statement is false or the proof cannot be repaired), **gap** (the statement is probably true but a step is missing or wrong), **presentation** (unclear, misleading or a wrong reference). For each, give the id and proof step, why it is wrong, and a proposed fix in the book's style: short declarative sentences, LaTeX in `$...$`, numbered proof steps;
- a verdict for each [Theorem] label introduced or changed: keep, keep after the stated fix, or downgrade (to what);
- the checks you ran, with their precision;
- prior art found, with full references.

Usage is limited: work efficiently. Prioritize items 1–3 over the rest.
