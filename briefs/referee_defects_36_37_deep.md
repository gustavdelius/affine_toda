# Brief: deep referee pass on the high-risk claims of chapters 36–37

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (the worktree `.claude/worktrees/integrable-defects-framework-a54a84`, branch `claude/integrable-defects-framework-a54a84`).

**Rules: other agents are editing files in this worktree at the same time.**
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase, restore).
- Do not run `quarto render`.
- Write scratch scripts only under `referee_scratch/defects_36_37_deep/` (create it).
- Write your report to `briefs/referee_defects_36_37_deep_report.md`.

## Why this pass

Chapters 36 (`part8-defects/36-defect-scattering.qmd`) and 37 (`part8-defects/37-defect-sine-gordon.qmd`) were drafted on 2026-10-11 by Claude agents. They were refereed once, lightly, by Gemini 3.8 Flash; that report and its integration are in `briefs/referee_defects_36_37_report.md`. That pass found two real but small errors. It mostly restated the drafting agents' own checks and did not test them independently. This pass is limited to the claims most likely to be wrong.

For each item, **do one concrete independent computation**: a derivation written out, or a script you write yourself. Do not reuse the book's check scripts as evidence; you may run them for comparison. A verdict without its own evidence counts as "unclear".

**Read first.**
- `briefs/part8_defects_conventions.md` (conventions, verified by script; source dictionary).
- `briefs/part8_defects_status.md` (what is verified and how).
- `briefs/part8_defects_sources.md` (digests of the sources, with their equation numbers).
- Source texts: pdftotext files in `~/referee-copy/papers/`, named by arXiv id with `/` → `_`.
- To find an object: `grep -rn '{#<id>}' part*/ index.qmd`.

**Conventions.**
- The Part VI coproduct is $\Delta e=e\otimes1+k\otimes e$, $\Delta f=f\otimes k^{-1}+1\otimes f$; the physical $q=-e^{-i\pi\lambda}$ (`eq-q-physical`); $y=e^{\lambda\theta}$.
- $\mathbb T:V(\theta)\otimes\mathcal W\to\mathcal W\otimes V(\theta)$, with tensor order equal to spatial order.
- $\sigma=e^{-\theta_d}$.

## Items (most important first)

1. **`prp-doubling` item 3 (doubled crossing-unitarity ⇔ defect crossing) for non-diagonal scattering.** It was checked numerically only for Lee–Yang, which is diagonal.
   - Construct the doubled $K$ of `eq-doubled-k` for the sine-Gordon $\mathbb T$ of ch. 37 (`eq-sg-transmission-matrix` with the scalar `eq-sgd-scalar`; a window of $\mathcal W$).
   - Test `eq-boundary-crossing-unitarity` (ch. 28) for it, with the doubled bulk S-matrix ($\check S$ on $V$, $P\check SP$ on the reflected copy, the flip between copies).
   - Do the index orders in ch. 28 and ch. 36 fit together?
2. **The sign of $q$** (ch. 37 `sec-sg-transmission-matrix`; conventions §9). The claim: unitarity holds at both $\pm q$, and crossing with $C=\sigma_x$ holds at $q_{\rm phys}$ but needs $C'=[[0,1],[-1,0]]$ at $-q_{\rm phys}$.
   - Is this a genuine selection, or could a gauge or Klein-factor redefinition at $-q$ restore $C=\sigma_x$?
   - Relate it to `sec-crossing-sign` (ch. 21), where the same question was settled for the soliton S-matrix.
   - Does the conclusion depend on including the scalar factor?
3. **`prp-rt-nogo` [Theorem]** (ch. 36). Write out the proof's two reductions yourself, with care about which rapidity differences appear.
   - Is "$S_{ab}$ constant" really forced by the mixed coefficient alone, given only meromorphy?
   - Is the hypothesis list complete (associativity, no defect states, diagonal S)?
   - Attempt a counterexample in the stated loophole: non-diagonal S with defect states. Either construct one, or argue why the natural candidates fail.
4. **`cnj-defect-parameter`** (ch. 37). The claimed renormalization is $\theta_d=(1-B/2)\theta_d^{\rm cl}$ for sinh-Gordon (from Bajnok–Simon (17)), continued to sine-Gordon as $\lambda\theta_d=(\lambda+1)\theta_d^{\rm cl}$.
   - Re-derive it from Bajnok–Simon.
   - Check it against Penati–Polvara's $B=\sqrt{1+\gamma/2}\,B^{(0)}$ (their (4.12), with $\gamma/2=\beta_{SG}^2/8\pi$).
   - Are the two statements compatible, contradictory, or about different quantities?
   - Is the chapter's description of the tension accurate?
5. **The scalar factor $\Phi$** (`eq-sgd-scalar`). Independently derive the minimal solution of unitarity and crossing for the ch. 37 matrix, and compare with Konik–LeClair and with BCZ 2005 (3.43)–(3.47).
   - Check the absence of poles in the physical strip.
   - Check the resonance at $\theta_d-i\pi/(2\lambda)$.
   - Check the claim that BCZ (3.42)/(3.47) "diverge" while the KL pairing (3.44) converges.
6. **Modified crossing** (ch. 36 `sec-topological-lines`). Take the RSOS kink S-matrices actually printed in ch. 27 (`sec-rsos`) for a minimal model with $\phi_{13}$ (e.g. $\mathcal M_5$ or $\mathcal M_{4,5}$) and compute the crossing behaviour yourself.
   - Do they satisfy ordinary crossing?
   - Do they differ from the Copetti–Córdova–Komatsu form by exactly $(d_ad_c/d_bd_d)^{i\theta/2\pi}$?
   - Is that factor a legitimate rescaling of kink states, i.e. does it preserve unitarity, Yang–Baxter and the bootstrap?
7. **`prp-half-soliton` and the correction to BCZ's two-defect remark** (ch. 37). Derive $z_{\rm odd}^2=A_{11}(\theta-\theta_d)$ and the quantum identity $\mathbb T_{B_1}(\theta;\vartheta+i\pi+i\xi/2)\,\mathbb T_{B_1}(\theta;\vartheta-i\pi-i\xi/2)=S_{sB_1}(\theta-\vartheta)$, and judge whether "the left-to-right factors are inverse" is a fair reading of BCZ 2005 Sec. 2.
8. **Defect Lee–Yang** (ch. 36 `sec-defect-lee-yang`).
   - The fundamental range $b\in[-3,2]$.
   - "Hermitian analyticity iff Re $b=\pm3$".
   - The excited state and Bajnok–Simon's $|2^+\rangle$ with its Coleman–Thun closure.
   - The defect-TBA driving term and the Lüscher formula, with their signs.

## Report format (`briefs/referee_defects_36_37_deep_report.md`)

A header with the date, model and commit (`git log -1 --format=%h`). Then one section per item, with:
- a verdict: correct / wrong / unclear;
- your own computation, with the script inline or its path and output;
- for errors, corrected text ready to paste, with file and line.

End with a short list of anything else you noticed, labelled as unverified.
