# Brief: draft chapter 39, "Solitons and particles through a defect"

*For a Claude drafting agent. Drafting is done by Claude only, for a consistent style; other models referee.*

You are drafting chapter 39 of Part VIII ("Affine Toda theory with defects") of a Quarto graduate textbook on affine Toda field theory, together with its verification script. You are working in the repository worktree `.claude/worktrees/integrable-defects-framework-a54a84`, branch `claude/integrable-defects-framework-a54a84`.

**Rules: other agents are drafting chapters 38 and 40 in this worktree at the same time.**
- You may edit only `part8-defects/39-defect-solitons.qmd` (replace the skeleton) and create `textbook_code/defect_solitons_an.py`.
- You may also create a report file, `briefs/draft_ch39_report.md`.
- Do not edit any other file: not `defect_common.py` (import from it only), not `references.bib`, not other chapters. List needed edits in your report instead.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase, restore).
- Do not run `quarto render`.

## Read first

1. `briefs/part8_defects_conventions.md` holds the conventions, verified by script. Use them exactly.
   - **Geometry and defect:** fields $\vec\phi$ on $x<0$ and $\vec\psi$ on $x>0$; the defect action and conditions (`eq-defect-action`, `eq-defect-conditions`); the $a_r$ matrix $\mathsf D=2(1-\pi)^{-1}$ (`eq-d-matrix`); the defect potential $\mathcal B_d$ in book normalization (`eq-defect-potential`).
   - **Parameters and types:** $\sigma=e^{-\theta_d}$; types I$_-$ and I$_+$ (named by Borel subalgebra; both keep the spin $\pm2$ charges); type II, Corrigan–Zambon's defect with an auxiliary field Λ.
   - **Lax and momentum:** the defect matrix $L_d$; momentum $P=-\int\phi_t\cdot\phi_x$.
   - **Sections 8–9:** the source dictionary and traps (σ = e^{+η} in cz2009b; the BCZ 2003 rapidity sign).
   - **Clash table:** e.g. $\mathbb T$ for transmission; $X_{ab}$ is Part V's fluctuation transmission factor; "folding" means Dynkin folding.
2. `briefs/part8_defects_status.md` lists what chapters 36–37 established and the **obligations for chapter 39** (under chapter 37):
   - generalize the absorption/emission direction rule to $a_r$;
   - generalize the odd-sector correspondence σ → −σ (sine-Gordon: $\theta_d\to\theta_d+i\pi$ classically) to $a_r$.
3. The drafted chapters. Reuse their results; do not repeat their derivations.
   - Ch. 36, `part8-defects/36-defect-scattering.qmd`.
   - Ch. 37, `part8-defects/37-defect-sine-gordon.qmd`: `eq-sg-delay` (the delay is the Bianchi coefficient), `prp-sg-delay`, `prp-half-soliton`, `eq-shg-linear-transmission`, `sec-type-two`.
   - Chapter 38 (`part8-defects/38-defect-conditions.qmd`) is being drafted now. Cite its ids from the registry at the end of the conventions brief (`eq-defect-conditions`, `eq-defect-potential`, `eq-d-matrix`, `eq-defect-lax-matrix`, `prp-defect-momentum`, `sec-defect-backlund-an`, `sec-type-two-atft`, `prp-defect-theta`, `tbl-defect-dictionary`), but do not describe their content in detail.
4. Sources: `briefs/part8_defects_sources.md` has the digests with source equation numbers.
   - The main ones: bcz2004b (hep-th/0401020, §6 soliton delay (6.2)), cz2007 (0705.1066, delay (2.20), classical selection rule (2.22), complex static energies (2.17)–(2.18), linearized/lightest-breather factors (7.5), (7.11) and its misprint), cz2009a (0902.1307), cz2009b (0908.3126, Tzitzéica delays (4.9)–(4.11)), cz2004 (gates).
   - Source texts: `~/referee-copy/papers/*.txt`. Check every formula against them.
5. Book background:
   - `part3-classical/10-solitons.qmd`: `eq-hirota-ansatz`, `eq-an-soliton`, `eq-n-soliton`, `eq-interaction-coefficient`, `sec-topological-charges`, `eq-soliton-mass`, `sec-multi-solitons` (time delays and Wigner's relation).
   - `part3-classical/11-reality.qmd`: Θ, complex solitons, `sec-missing-charges`.
   - `part5-semiclassics/16-dhn.qmd`: `eq-an-transmission`, Part V's $X_{ab}$.
   - Code: `textbook_code/toda_classical.py` (`tau_funcs`, `field_residual`, `an_alpha`, `soliton_mass`), `textbook_code/defect_common.py` (`dmatrix`, `defect_potential`, `defect_condition_residual`, `an_data`).
6. Model chapter: `part7-half-line/31-half-line-solitons.qmd` (Solitons, breathers and vacua on the half-line). Mirror its structure, depth (about 260–320 lines) and tone. Its method of images is the boundary analogue of what you do here.
7. House style, as in chapters 36–37:
   - learning goals; `##` sections with the skeleton ids;
   - `::: {#prp-…}` with a bold status label first, **[Established]**, **[Numerical]**, **[Proposition, sketch]** or **[Conjecture]** (never a cross-reference inside the bracket);
   - `::: {.proof}`; equations cited as `(-@eq-…)`;
   - 7–8 exercises with ids `exr-ds-…`;
   - Notes with the history and "Corrections and clarifications to the sources made in this chapter:".
   - Cite only keys that exist in `references.bib` (grep it); new ids get the prefix `ds-`.

## Content (keep every skeleton id)

- **`sec-an-soliton-defect`.** An $a_r$ soliton of species $a$ on $x<0$ and a shifted soliton on $x>0$, with τ-functions in book form, satisfy `eq-defect-conditions` for all $t$ when the shift is $z_a(\theta-\theta_d)$.
  - `eq-an-defect-delay`: derive it, compare with bcz2004b (6.2) ($\Lambda=(ie^\theta+de^{i\zeta_a})/(ie^\theta+de^{-i\zeta_a})$) and cz2007 (2.20), and relate it to the Bäcklund/Bianchi picture of ch. 37.
  - Treat both types I$_\mp$.
  - `prp-an-defect-delay` [Established; numerically checked for $a_2$, $a_3$, all species].
- **`sec-defect-topological-charge`.** A soliton's topological charge can change at the defect, which stores charge in the coweight lattice $(2\pi/\beta)\Lambda_W^\vee$. Give cz2007's classical selection rule (2.22) in book terms, and say which members of a multiplet emerge. `prp-stored-charge` [Proposition or Established].
- **`sec-defect-absorption`.** Absorption and emission rapidities for $a_r$, and the direction rule (obligation above); the stored energy and momentum (for sine-Gordon a soliton of rapidity $\theta_d$; what is the $a_r$ analogue?).
- **`sec-defect-multi`.** Several solitons: the delays add (test with two solitons), and a classical consistency check of factorization.
- **`sec-defect-reality`.** Complex solitons at imaginary coupling; Θ-images; complex delays; the odd-sector correspondence σ → −σ (obligation above); real σ (cf. `prp-defect-theta`).
- **`sec-classical-transmission-factors`.** Linearize about the vacuum at real coupling.
  - `eq-an-linear-transmission`: $R=0$, $|\mathbb T_a^{\rm cl}|=1$, the closed form, and crossing $a\leftrightarrow\bar a$. Compare with cz2007 (7.5) and with (7.11), whose denominator should read $ie^{-\theta}m_k+\sigma(1-\omega^k)$ (digest); confirm or refute.
  - `prp-an-linear-transmission` [Numerical].
  - `prp-defect-frozen-soliton` [Proposition, sketch, or Conjecture]: the "frozen soliton" comparison of the defect factors with the soliton interaction coefficient and with Part V's $X_{ab}$ (`eq-an-transmission`). Report the outcome either way.

## Script `textbook_code/defect_solitons_an.py`

Print PASS/FAIL lines in the style of the other `textbook_code/defect_*.py` scripts. The docstring lists the book ids each check supports; the runtime should be under 3 minutes. It checks:

- (a) the delay $z_a$ for $a_2$ and $a_3$, all species, both types: the defect conditions hold at several times (residual about $10^{-12}$), with the closed form compared against bcz2004b (6.2) and cz2007 (2.20);
- (b) the stored charge lies in the coweight lattice, with the change in charge equal to the difference of soliton weights;
- (c) the absorption and emission rapidities;
- (d) two-soliton additivity;
- (e) the linear $\mathbb T_a^{\rm cl}$: $R=0$, $|\mathbb T|=1$, the closed form, $a\leftrightarrow\bar a$ crossing, and the cz2007 (7.11) misprint;
- (f) the frozen-soliton comparison, reported either way.

Run it; every check should PASS. Report honestly anything that fails and what it means.

## Report (`briefs/draft_ch39_report.md`)

1. The chapter summary: sections, results with status labels, exercises, new ids.
2. The script output and runtime.
3. The verification record: by script (which check), by hand, unchecked.
4. Corrections to the sources.
5. Requests: edits needed elsewhere, bib keys, updates to the conventions brief, and obligations for chapters 41–43.
