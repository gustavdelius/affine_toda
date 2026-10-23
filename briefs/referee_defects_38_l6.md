# Brief: referee pass on the three new [Theorem] labels in chapter 38 (PR #40)

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory, in the worktree `.claude/worktrees/gallant-kilby-b8f7ed` on branch `claude/gallant-kilby-b8f7ed`, which is PR #40.

**Rules.**
- Other agents may be editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Do not spawn subagents.
- Write scratch scripts only under `referee_scratch/defects_38_l6/`.
- Write the report to `briefs/referee_defects_38_l6_report.md`.

## Context

On 2026-10-23 a Claude agent added three results to chapter 38 (`part8-defects/38-defect-conditions.qmd`, `sec-defect-conditions`). They came out of Lecture 6 of the companion lecture notes (`lecture-notes/defects/`), but the book's text and scripts were written and checked independently of the lecture. **No referee has seen them.** Your earlier report on this chapter, `briefs/referee_defects_38_40_report.md`, does not cover them. Do not re-referee anything else in the chapter.

- **Commits.** The additions are commit `62d4bb6`; see them with `git show 62d4bb6 -- part8-defects/38-defect-conditions.qmd`. Commit `ccb4c15` merged `main`, which renamed `textbook_code/` to `code/`. Review the text as it stands at `ccb4c15`. This brief was committed on top of it.
- **Under review** (find each with `grep -n '{#<id>}' part8-defects/38-defect-conditions.qmd`):
  1. `prp-defect-involution` (`sec-defect-r-matrix`). The $a_1^{(1)}$ case was upgraded from [Numerical] to **[Theorem]**, with a new hand proof and a new uniqueness claim. The label is now [Theorem for $a_1^{(1)}$; Numerical for $a_2^{(1)}$–$a_4^{(1)}$]. Also review the paragraph after the proposition, which starts "The normalization is the point".
  2. `prp-sg-spin-three` and `eq-dc-spin-three`, new, **[Theorem]**, at the end of `sec-defect-lax`. Also review the paragraph that introduces them, which starts "The spin-$\pm3$ charges of sine-Gordon theory", and the paragraph after the proof.
  3. `prp-dc-gauge` and `eq-dc-gauge`, new, **[Theorem for type I$_-$; Numerical for type I$_+$]**, in `sec-defect-backlund-an`. Also review the paragraph before it and the one after its proof.
  4. The new sentences in the Notes (`sec-defect-conditions-notes`, the "stated here" paragraph) and the changed learning goal.
- **Check scripts** (rows in `code/README.md`):
  - `code/defect_sklyanin_a1.py`: 7 checks, a few seconds;
  - `code/defect_sg_spin3.py`: 7 checks, about 15 s;
  - `code/defect_ar_gauge.py`: 3 checks, about 10 s.

  The older `code/defect_lax_rmatrix.py` (6 checks) is the numerical check for $r\ge2$.
- **Drafter's flags for the referee.**
  - (a) The $a_1$ [Theorem] of `prp-defect-involution` uses the bulk relation $\{T\overset\otimes,T\}=-[r,T\otimes T]$ of `prp-sklyanin` (ch. 30, [Established]). The new script re-derives the bracket of $a_x$ for $a_1$.
  - (b) The proof of `prp-sg-spin-three` is an exact computer-algebra identity: the prose gives the reduction, the script does the algebra. A precedent is `prp-spin-two` part (3), whose $a_2$ case is a [Theorem] for the same reason.
  - (c) The closed form of $\mathcal D_{-3}$ is credited to Caudrelier, in light-cone coordinates. No Appendix E item was added for any of the three results.
  - (d) Ch. 9, after `eq-spin3-current`, calls the $\partial_+$ charge "spin 3" with eigenvalue $e^{3\theta}$. This clashes with the convention below and is being fixed in a separate session. Confirm or refute the clash in one line under "other observations" and do not report it further.

## Notation and conventions (read before computing)

The book's conventions are in `01-introduction.qmd` (`sec-conventions`) and in `briefs/part8_defects_conventions.md`; read its §§1–5 first. Convention slips are the most common false findings. In particular:

- **Fields and Lagrangian.** The field is $\vec\phi$ on $x<0$ and $\vec\psi$ on $x>0$. The Lagrangian is $\frac12(\partial\vec\phi)^2-V$, with no $1/\beta^2$ in front (`sec-lagrangians`), so the bulk bracket is $\{\phi_a(x),\partial_t\phi_b(y)\}=\delta_{ab}\delta(x-y)$.
- **Coupling and light-cone derivatives.** $b=\beta$ at real coupling and $b=i\beta$ at imaginary coupling. $\partial_\pm=\partial_t\pm\partial_x$.
- **Spin.** Spin $s$ goes with $e^{s\theta}$, so $E+P$ has spin $+1$. Densities built from $\partial_+$ have negative spin.
- **Defect bracket** (`eq-dc-defect-bracket`). $\{z_I,z_J\}=(K^{-1})_{IJ}$ with $K=\begin{pmatrix}E&\mathsf D\\-\mathsf D^{\mathsf T}&E\end{pmatrix}$. For $a_1$, $\{\phi(0),\psi(0)\}=-1$ for type I$_-$ and $+1$ for type I$_+$. The proof calls this number $\epsilon$.
- **r-matrix** (`prp-sklyanin`, ch. 30). $r(z)=\frac{\beta^2}{4(z^h-1)}[(z^h+1)\sum_aH_a\otimes H_a+2\sum_{\alpha>0}\frac{\vec\alpha^2}2(z^{\rm ht\,\alpha}E_\alpha\otimes E_{-\alpha}+z^{h-{\rm ht\,\alpha}}E_{-\alpha}\otimes E_\alpha)]$, with $\beta^2\to b^2$. The local bracket is $\{a_x\overset\otimes,a_x\}=+[r,a_{x1}+a_{x2}]\delta$. With $\partial_xT=-a_xT$ this gives $\{T\overset\otimes,T\}=-[r,T\otimes T]$.
- **$a_1$ in two dimensions.** $H=\operatorname{diag}(1,-1)/\sqrt2$, $E_{\vec\alpha}=\begin{pmatrix}0&1\\0&0\end{pmatrix}$, $\sigma_x=E_{\vec\alpha}+E_{-\vec\alpha}$, $\Pi=E_{\vec\alpha}\otimes E_{-\vec\alpha}+E_{-\vec\alpha}\otimes E_{\vec\alpha}$. The type-I$_-$ defect matrix is $L_d=e^{-bH\psi/2}(1+\frac\sigma\lambda\sigma_x)e^{bH\phi/2}$ (`exr-dc-lax`, `eq-defect-lax-matrix`). The proof writes it as $\mathsf J+\frac\sigma\lambda\mathsf O$; $\mathsf D$ is reserved for the kinetic matrix.
- **$a_r$ components.** In the defining representation, $E_{\vec\alpha_j}=e_{j,j+1}$ with indices mod $h=r+1$ and $H\cdot\vec v=\operatorname{diag}(\vec\mu_k\cdot\vec v)$. The components are $\phi_k=\vec\mu_k\cdot\vec\phi$, with $f_j=e^{b(\psi_j-\phi_{j+1})}$ and $g_j=e^{b(\phi_j-\psi_j)}$ (`eq-dc-components`). The component Bäcklund transformation is the display just before `prp-dc-gauge`'s introductory paragraph.
- **Sine-Gordon** (ch. 37, `sec-defect-sine-gordon`). $\beta=\beta_{\rm SG}=\sqrt2\beta_{\rm book}$ and $m_0=2m$. The soliton mass is $M=8m_0/\beta^2$, and the soliton is $\frac4\beta\arctan e^{m_0(x\cosh\theta-t\sinh\theta)+c}$, with momentum $M\sinh\theta$. $\sigma=e^{-\theta_d}$. The conditions are `eq-sg-defect-conditions`, with the potential `eq-sg-defect-potential`.
- **The lecture notes** use a dimensionless field $\Phi=\beta_{\rm SG}\phi$ and a bracket $c=\{\Phi,\Psi\}=\beta_{\rm SG}^2\epsilon$. Do not use them as evidence.

## Read first

- **Sections of ch. 38:** `sec-defect-lagrangian`, `sec-defect-lax`, `sec-defect-r-matrix`, `sec-defect-backlund-an`.
- **Elsewhere:** ch. 30 `prp-sklyanin`; ch. 37 `prp-sg-defect-momentum`, `prp-half-soliton`, `sec-sg-defect-solitons`; ch. 9 `eq-lax-pair`, `eq-spin3-current`.
- **Sources:**
  - Bowcock–Corrigan–Zambon, hep-th/0305022 (`bcz2004`) and hep-th/0401020 (`bcz2004b`);
  - Caudrelier, 0704.2326 (`caudrelier2008`);
  - Habibullin–Kundu, 0709.4611 (`habibullin_kundu2008`);
  - Caudrelier, 1411.5171 (`caudrelier2015`).

  The texts are in `~/referee-copy/papers/*.txt`, except 0709.4611, which is not there; download it from arXiv if you need it for item 4.
- **Finding an object:** `grep -rn '{#<id>}' part*/ appendices/ outlook/ *.qmd`.

## Evidence rule

- For each item, do **one independent computation**: a derivation written out, or a script of your own.
- You may run the book's scripts for comparison, but they are not evidence.
- A verdict without its own evidence counts as "unclear".
- `briefs/referee_defects_38_40_report.md` and `briefs/referee_defects_upgrades_report.md` show the expected depth.

## Items (most important first)

1. **`prp-defect-involution`: the $a_1^{(1)}$ case, now [Theorem].** Check the proof line by line.
   - **Brackets.** Derive the brackets of $\mathsf J=e^{bH(\phi-\psi)/2}$ and $\mathsf O=e^{-bH(\phi+\psi)/2}\sigma_x$ from $\{\phi(0),\psi(0)\}=\epsilon$: are $\{\overset1{\mathsf J},\overset2{\mathsf O}\}=\kappa\,H\otimes H\,\mathsf J\otimes\mathsf O$ and $\{\overset1{\mathsf O},\overset2{\mathsf J}\}=-\kappa\,H\otimes H\,\mathsf O\otimes\mathsf J$, with $\kappa=-b^2\epsilon/2$?
   - **Commutators.** Check $[H\otimes H,X]=2H\otimes H\,X$, $[\Pi,X]=-2H\otimes H\,X'$, the same with $X$ and $X'$ exchanged, and that $H\otimes H$ and $\Pi$ commute with $\mathsf J\otimes\mathsf J$ and $\mathsf O\otimes\mathsf O$.
   - **Matching.** Solve the matching equations and confirm the result $r=\frac\kappa{2(z^2-1)}[(z^2+1)H\otimes H+2z\Pi]=-\epsilon\,r_{\rm book}$. Reduce the general formula of `prp-sklyanin` to $h=2$ yourself. Is the coefficient exactly 1, with the right sign, for type I$_-$? This sign was wrong once before, in ch. 30.
   - **Uniqueness without an ansatz.** The proof claims that two solutions differ by a matrix commuting with $L_d(\lambda_1)\otimes L_d(\lambda_2)$ for all defect values, hence by a multiple of the identity when $\lambda_1^2\ne\lambda_2^2$. Check the separation into exponentials of $\phi$ and $\psi$, and whether all the exponents named in the proof are distinct. Check the eigenvector argument on $W$, and whether the conclusion fails at $\lambda_2=\pm\lambda_1$.
   - **Type I$_+$.** Check $L_d^{\rm I_+}(\phi,\psi)\sigma_x=\frac\lambda\sigma L_d^{\rm I_-}(\phi,-\psi)$, with $L_d^{\rm I_+}$ from `eq-defect-lax-matrix` and $\tilde L_{\rm I_+}$ for $h=2$. Check that the type-I$_+$ bracket becomes $\{\phi,-\psi\}=-1$.
   - **Bulk relation.** Derive $\{a_x(\lambda_1;x)\overset\otimes,a_x(\lambda_2;y)\}=\delta(x-y)[r,a_x(\lambda_1)\otimes1+1\otimes a_x(\lambda_2)]$ for $a_1$ in the book's normalization, with $a_x$ from `eq-boundary-lax-at`.
   - **Scope.** The statement includes involution, $\{\tau_d,\tau_d\}=0$, under the stated hypothesis (independent defect values that Poisson-commute with the bulk). Is the proof complete for $a_1$ given that hypothesis and `prp-sklyanin`? Is the label justified?
2. **`prp-sg-spin-three`, new [Theorem].** Do the computation yourself, from `eq-sg-defect-conditions` and `eq-sg-defect-potential`.
   - **Conservation.** Check that $Q_{\mp3}$, as defined, is conserved for the type-I$_-$ defect. Eliminate $\partial_x$, $\partial_t\partial_x$ and $\partial_x^2$ of the defect values. Show that $j_{\mp3}(\psi)-j_{\mp3}(\phi)=-\frac d{dt}\mathcal D_{\mp3}$ for both signs, with the normalizations and the sign of the current as written.
   - **The proof's claims.** Check the claim that the terms quadratic in second time derivatives cancel, and the reason the proof gives.
   - **Soliton values.** Check $Q_{\mp3}=-\frac23e^{\mp3\theta}$ for a soliton, and $\frac{\beta^2}{8m_0}(E\mp P)=e^{\mp\theta}$.
   - **Jump values.** Check the values of $\mathcal D_{\mp3}$ at $(0,2\pi/\beta)$ and at $(0,0)$, the stored charge, and the half-soliton statement in the paragraph after the proof, against `prp-half-soliton` part 1.
   - **Spin labels.** Does "spin $\mp3$" agree with the convention above?
   - **Label.** Is [Theorem] justified when the decisive identity is checked by computer algebra? Should the proof show more of the computation? A by-hand derivation of the identity, even for one sign, would settle it.
3. **`prp-dc-gauge`, new [Theorem for type I$_-$; Numerical for type I$_+$].**
   - **Components.** Derive the components of $L_d^{\rm I_-}$ from `eq-defect-lax-matrix`, `eq-d-matrix` and `eq-dc-components`. Derive the component form of $\mathcal A_\pm$ from `eq-lax-pair` with $\beta\to b$ and $c_j=1$.
   - **The orders.** Check each order of $\lambda$ in both light-cone equations. Check that the surviving orders are exactly the component Bäcklund transformation of `sec-defect-backlund-an`, signs and factors $m\sigma/b$ and $m/(b\sigma)$ included, and that no order is missing.
   - **Curvature.** Check the conjugation statement, $\det L_d^{\rm I_-}=1-(-\sigma/\lambda)^h$, and "fails only at finitely many $\lambda$" for type I$_+$ as well.
   - **Type I$_+$.** Is there a short proof, for example through `prp-dc-types`, through $L_d^{-1}$, or by the same order-by-order computation? If so, give it; the label could then be [Theorem] for both types.
   - **Time component.** Check the claim that half the sum of the two equations is the evolution equation of `eq-defect-lax-matrix`.
4. **Prior art and novelty.** Search beyond the sources named, including the editor's own papers. For each item below, quote the source line.
   - **Caudrelier 2008.** Is his (3.28), on his (3.29), a constant multiple of $\mathcal D_{-3}$? Is it fair to call the closed form his, given that his defect sits on a line of constant light-cone coordinate? Does he, or anyone, give $\mathcal D_{+3}$, the defect at rest, or the stored charges?
   - **Bowcock–Corrigan–Zambon 2003.** Do their (2.4)–(2.5) use only the combination $Q_3+Q_{-3}$, without writing its defect term?
   - **Bowcock–Corrigan–Zambon 2004b.** Do they impose only the time equation (3.7) on $K$, and get the Bäcklund transformation by cross-differentiation (Sec. 5.1)?
   - **Habibullin–Kundu 2008.** Is "the same algebraic structure, Cartan exponentials on the diagonal and an off-diagonal partner with a constant bracket" a fair description of their eqs. (5.73)–(5.81)? Does the dependence on the lattice momenta drop out in the continuum limit (Sec. 5.2)?
   - **The Sklyanin relation for $a_1$.** Is the relation, with the bracket taken from the defect term of the action, and with uniqueness of $r$, already in the literature? Look at Avan–Doikou 2012 and Caudrelier 2015, among others.
   - **The gauge property for $a_r$.** Is it known, for example as a Darboux or dressing matrix of the periodic 2D Toda lattice?
   - **Appendix E.** Should any of the three results have an item in `appendices/E-original-results.qmd`? Is the Notes paragraph accurate about what is "stated here"?
5. **Citations.** For every citation added in commit `62d4bb6`, check that the source says what the text says, at the cited equation or section. The citations are:
   - `bcz2004` (2.4)–(2.5);
   - `bcz2004b` (3.7) and Sec. 5.1;
   - `caudrelier2008` (1.12)–(1.13), (3.26), (3.28) and (3.29);
   - `habibullin_kundu2008` (5.73)–(5.81), Sec. 5.2 and Sec. 5.

## Report format (`briefs/referee_defects_38_l6_report.md`)

A header with the date, your model and the commit reviewed (`ccb4c15`). Then one section per item, with:
- the verdict (correct / wrong / unclear);
- your own computation (the script inline or its path, and its output);
- for errors, replacement text ready to paste, with file and line.

End with other observations, labelled as unverified.
