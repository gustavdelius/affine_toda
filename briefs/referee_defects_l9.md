# Brief: referee pass on the book edits from Lecture 9 of the defect lecture notes

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory, in the worktree `.claude/worktrees/jovial-antonelli-69b2cf` on branch `claude/defects-l9-book`.

**Rules.**
- Other agents may be editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Do not spawn subagents.
- Write scratch scripts only under `referee_scratch/defects_l9/`.
- Write the report to `briefs/referee_defects_l9_report.md`.

## Context

On 2026-10-25 a Claude agent carried five findings of Lecture 9 of the companion lecture notes (`lecture-notes/defects/09-folding-fusing.qmd`, PR #56) into Part VIII of the book. The book's text and its check script were written and checked independently of the lecture. **No referee has seen them.** Review only what is listed here.

- **Commit.** The book edits are commit `c444a86`; see them with `git show c444a86`. This brief was committed on top of it.
- **Under review** (find each with `grep -rn '{#<id>}' part8-defects/ appendices/`):
  1. `prp-defect-fusion-k` (ch. 43, `part8-defects/43-defects-and-boundaries.qmd`, `sec-defect-boundary-fusion`). Item (3) was upgraded from [Numerical] to **[Theorem]**; the label is now [Theorem (1), (3); Numerical (2); (4) by hand, given the quantum UV–IR matching]. Its new proof contains `eq-dbd-scalar-identity` and the integral representation `eq-dbd-sigma-integral` of the Ghoshal–Zamolodchikov factor $\sigma(x,u)$. Also review the paragraph after the proof ("We checked (3) on windows") and the changed sentence in the Notes (`sec-defect-boundary-notes`, the "stated here" sentence).
  2. `prp-sg-fused-spin-three` and `eq-sgd-fused-spin-three`, new, **[Theorem]** (ch. 38, `part8-defects/38-defect-sine-gordon.qmd`, in the paragraph "Fusing a defect onto a boundary" of `sec-sg-breather-transmission`), with the two sentences before it and the paragraph after its proof.
  3. The new fourth bullet of the status list after `cnj-defect-parameter` (ch. 38), which starts "A boost gives the factor without a boundary". It is presented as a heuristic.
  4. The new paragraph "Linear waves see the condition" at the end of `sec-defect-folding` (ch. 39, `part8-defects/39-defect-conditions.qmd`).
  5. The new sentences in "Reading the result" after `prp-diagonal-borel` (ch. 43), from "and the halves go with the right chiralities" to "as it rescales the action on the particles".
  6. The new item of the list of original results (`appendices/E-original-results.qmd`, Defects, the item citing `prp-defect-fusion-k` and `prp-sg-fused-spin-three`), and the citation change `[@cz2012, Sec. 2]` → `Sec. 3` in ch. 38.
- **Check script** (row in `code/README.md`): `code/defect_fused_boundary.py`, 9 checks, about 2 minutes. The older `code/defect_doubling.py` (20 checks) is where the scalar identity was first found numerically.
- **Drafter's flags for the referee.**
  - Item 1 is the only label upgrade. The proof uses Malmstén's formula with grouping by $l$ to justify exchanging sum and integral, and the Riemann–Lebesgue lemma for $\sigma(i\vartheta,u)\to1$. Check both steps and the stated strip $|\operatorname{Re}x|+|\operatorname{Re}\lambda u|<\frac\pi2$.
  - Item 5 settles a doubt raised while writing the lecture: the lecture's handoff (`briefs/defects_lectures_handoff.md` on branch `claude/defects-lecture-9`, §7.4, "Left open") thought the assignment of the two Borel halves to the chiralities looked reversed. The drafter concluded that the book's reading is right when chirality is read off from the Lorentz spin, and that the doubt came from `sec-nonlocal-charges` (ch. 21), which builds the charges $Q_j$, identified with the $e_j$ of spin $+\lambda$, from $\phi_L(x^+)$. In the book's conventions a left-mover carries $E-P$, of spin $-1$. Ch. 21 and the single-chirality argument of ch. 41 ("Which non-local charges survive", `sec-perturbed-defects`) were **not** changed. Give your view on them (item 5 below).

## Notation

- $\lambda=8\pi/\beta_{\rm SG}^2-1$ is the coupling, $\beta=\beta_{\rm SG}$ in the spin-3 formulas, $\hbar=1$.
- The physical $q=-e^{-i\pi\lambda}$ (`eq-q-physical`); the coproduct is $\Delta e=e\otimes1+k\otimes e$, $\Delta f=f\otimes k^{-1}+1\otimes f$ (Part VI).
- $\mathbb T(\theta):V(\theta)\otimes\mathcal W\to\mathcal W\otimes V(\theta)$; $y=e^{\lambda\theta}$, $y_d=e^{\lambda\theta_d}$, $Y=y/y_d$; $e_1v_-=yv_+$, $f_1v_+=y^{-1}v_-$.
- $\Phi$ is the scalar factor `eq-sgd-scalar`, with integral form `eq-sgd-scalar-integral` (`prp-sgd-scalar-integral`, ch. 38).
- Reflection matrices use $u=-i\theta$ and the parameters $(\eta,\vartheta)$ of `eq-bsg-k-full` (ch. 30); $\sigma(x,u)$ is the Ghoshal–Zamolodchikov factor printed there, not the defect parameter.
- The defect parameter is $\sigma=e^{-\theta_d^{\rm cl}}$ (`eq-sg-defect-potential`); $\theta_d$ is the quantum rapidity.
- $x^\pm=t\pm x$, $\partial_\pm=\partial_t\pm\partial_x$, and $E+P$ has spin $+1$; a density built from $\partial_+$ has negative spin (ch. 39, before `prp-sg-spin-three`).

## Read first

- `briefs/part8_defects_conventions.md`.
- Ch. 30 `sec-bsg-k-matrix` (the scalar factors $R_0$, $\sigma$), `cnj-bsg-uv-ir`; ch. 38 `prp-sgd-scalar-integral` and its proof, `sec-sg-transmission-matrix`, `cnj-defect-parameter` and its status list; ch. 39 `prp-sg-spin-three`, `prp-dc-types`, `exr-dc-d-matrix`, `sec-type-two-atft`; ch. 40 `eq-an-linear-transmission`; ch. 41 `eq-dq-massless-defect`, `sec-perturbed-defects`, `prp-two-borels`; ch. 43 `prp-diagonal-borel`, `eq-dbd-second-module`, `prp-defect-twist-limit`; ch. 21 `sec-nonlocal-charges`, `sec-qg-gradation`.
- Sources, by arXiv number (texts in `~/referee-copy/papers/*.txt` on the editor's machine): hep-th/9703085 (Konik–LeClair, Sec. 5.1 and 5.2), 1307.4536 (Bajnok–Holló–Watts, eqs. (96)–(102)), 0712.4292 (Bajnok–Simon), 1202.6016 (Corrigan–Zambon 2012), 0908.3126 (Corrigan–Zambon 2009b, Sec. 4.2 and 5), 1304.3129 (Robertson, Sec. 3.1), hep-th/9306002 (Ghoshal–Zamolodchikov).
- To find an object: `grep -rn '{#<id>}' part*/ appendices/ outlook/ *.qmd`.

## Evidence

For each item, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_38_40_report.md` shows the expected depth.

## Items (most important first)

1. **`prp-defect-fusion-k` (3), upgraded to [Theorem].**
   - From the Gamma product of $\sigma(x,u)$ in `eq-bsg-k-full`, derive `eq-dbd-sigma-integral` yourself, including the strip of validity and the justification for exchanging the sum over $l$ and the integral. Check it numerically against the product at complex $(x,u)$.
   - Check that $\sigma(i\vartheta,u)\to1$ as $\vartheta\to\infty$, and the limit $K_D$ of `eq-bsg-k-full` used in the proof.
   - Redo the sector computation: $K_{\rm f}=\bar{\mathbb T}(1\otimes K_D)\mathbb T$ on $\{v_+\otimes|n\rangle,v_-\otimes|n+1\rangle\}$, the identities $k_\pm+y_d^{-2}k_\mp=4y_d^{-1}P_0^\pm$ and $\frac i{y_d}(y^{-2}-y^2)=-4y_d^{-1}Q_0$ at $\vartheta=\lambda\theta_d$, and the conjugating matrix $\Lambda D_n$.
   - Prove or refute `eq-dbd-scalar-identity`, and check it numerically from the Gamma products of $\sigma$ and $\Phi$ at several couplings (use generic $\lambda$, not short decimals) and at complex rapidities. Check that the offset $s$ drops out.
   - Is the [Theorem] label justified? Is anything in (3) still only numerical?
2. **`prp-sg-fused-spin-three`, new [Theorem].**
   - With the densities and defect terms of `prp-sg-spin-three` (`eq-dc-spin-three`) in the book's normalization, check exactly that $Q_3^{\rm f}$ of `eq-sgd-fused-spin-three` is conserved on the fused boundary condition $\partial_x\phi=-\partial_\phi\mathcal B_d(\phi,\phi_D)$ and the field equation, for symbolic $\sigma$ and $\phi_D$. Check the coefficient and the sign of the Dirichlet term, and that the charge fails without it.
   - Check the proof's interval argument: the Dirichlet current $j_{-3}+j_{+3}=-\frac{\beta}{m_0}\sin(\beta\phi_D)\,\partial_t\psi_x$, and the value of $\psi_x$ given by the second defect condition.
   - Check the statement that a boundary can keep only $Q_{-3}+Q_{+3}$, and the attribution to Ghoshal and Zamolodchikov.
3. **The boost heuristic after `cnj-defect-parameter`.**
   - Check $h=\beta_{\rm SG}^2/8\pi=1/(\lambda+1)$ for $e^{i\beta_{\rm SG}\phi_L}$ in the normalization of `sec-nonlocal-charges`, and that the two halves of $\mathcal B_d$ are the chiral fields stated, at the conformal point (`eq-dq-massless-defect`).
   - Check the scaling argument: a boost by $\zeta$ multiplies the couplings by $e^{\mp(1-h)\zeta}$, and with $\sigma=e^{-\theta_d^{\rm cl}}$ this gives $\theta_d^{\rm cl}=(1-h)\theta_d$ up to a constant, i.e. $\lambda\theta_d=(\lambda+1)\theta_d^{\rm cl}$, and $\theta_d=(1-\frac B2)\theta_d^{\rm cl}$ after $\lambda\to-2/B$.
   - Quote Konik–LeClair's statement of the scaling dimension of $g$ (Sec. 5.1) and check it equals $1-h$ in their normalization.
   - Quote Bajnok–Holló–Watts (101)–(102) and their surrounding text, check the equation numbers (97)–(102), and check that with $\theta_b=i\pi(3-b)/6$ (`sec-defect-lee-yang`) their couplings are $e^{\pm\frac65(\theta_b-i\pi)}$ times constants, and that the Lee–Yang defect fields have weight $-\frac15$.
   - Is the bullet clearly a heuristic, and are its assumptions stated? Does it overclaim?
4. **The folding paragraph in `sec-defect-folding`.**
   - Check that the mirror image (`prp-dc-types`, `exr-dc-d-matrix`) transmits particle $k$ with $\mathbb T^{\rm cl}_{h-k}$, from the linearized defect conditions or from the action of the reflection on the modes.
   - Check the ratio $\mathbb T^{\rm cl}_k/\mathbb T^{\rm cl}_{h-k}$, the conclusion $\sigma_2=\pm\sigma_1$ for $2k\ne h$, the transparency at $-\sigma_1$, the product $1/f_{k/h}(u)$ with the block `eq-fx-block`, and the folded soliton delay.
   - Check that for $h=3$ the delay equals the Tzitzéica $z_1$ printed in `sec-type-two-atft` and in Corrigan–Zambon 2009b (4.9)–(4.10), and in which sector.
   - Check the attributions: Robertson Sec. 3.1 for the condition and the reality of the delay, and the remark that Corrigan and Zambon rescale $\sigma$ when they fold (2009b, Sec. 5).
5. **Chirality and the two Borel halves (ch. 43, "Reading the result").**
   - Check, in the book's conventions, which chirality carries negative Lorentz spin, which generators ($e_j$ or $f_j$) have spin $-\lambda$ on $V(\theta)$, and which chirality the type-I$_-$ and type-I$_+$ defects pass at $m_0=0$. Is the new text right that the untwisted half of type I$_-$, $U_q(\mathfrak b^-)$, goes with the left-movers it passes?
   - Check $\rho(f_j)\propto y_d^{-1}$ (Weyl module, `sec-sg-transmission-matrix`) and $\tilde\rho(e_j)\propto y_d$ (`eq-dbd-second-module`).
   - Not under review, but report your view: is `sec-nonlocal-charges` (ch. 21), which builds the charges $Q_\pm$, identified with $e_1,e_0$ of spin $+\lambda$, from $J_\pm=e^{\pm ia\phi_L}$ with $\phi_L$ a function of $x^+=t+x$, consistent with the book's spin convention? Is the sentence of ch. 41 "U_q(b^-) … when the perturbation is left-moving, U_q(b^+) … when it is right-moving" consistent with it? If not, what is the smallest correct change, and does it touch the derivation of the sign of $q$ from the braiding of chiral vertex operators? Answer with your own computation; this is the most likely place for a convention slip.
6. **Prior art, novelty and citations.**
   - The new item of Appendix E claims as original: the scalar identity of the fused Dirichlet boundary, with its proof, and the spin-3 charge of the fused boundary. Search for earlier statements, at least in Bajnok–Simon 0712.4292, Corrigan–Zambon 1202.6016 (which uses the spin-3 argument for a type-II boundary), Bajnok–George hep-th/0404199, and the work of Caudrelier, Avan and Doikou on defects and boundaries. Is the claim of novelty justified?
   - Check that Corrigan–Zambon 1202.6016 put Habibullin's remark in Sec. 3, and each citation of items 3 and 4 against its source.

## Report format

- A header with the date, your model and the commit you reviewed.
- One section per item: the verdict (correct / wrong / unclear); your own computation (script inline or path, and its output); and, for errors, replacement text ready to paste, with file and line.
- End with other observations, labelled as unverified.
