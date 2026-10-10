# Brief: referee chapters 30–31 (classically integrable boundary conditions; solitons on the half-line)

You are a mathematical-physics referee. The repository is `/home/gustav/Git/affine_toda` (branch `claude/affine-toda-boundary-dae69e`, PR #6), a Quarto book on affine Toda field theory. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files.

## Context

Part VII, "Affine Toda theory on the half-line" (`part7-half-line/`, chapters 28–35), was drafted on 2026-10-10 by agents working from the boundary literature. This brief covers:

- **Chapter 30** (`30-boundary-conditions.qmd`, `sec-boundary-conditions`). It contains:
  - the book's boundary convention, $\mathcal B_{\rm real}=\frac{2m}{\beta^2}\sum_jC_j\sqrt{2n_j/\alpha_j^2}(e^{\beta\alpha_j\cdot\phi/2}-1)$, with $\mathcal B_{\rm imag}$ obtained by $\beta\to i\beta$;
  - Sklyanin's formalism with an anti-automorphism σ, following Delius hep-th/9809140;
  - the node-by-node boundary K-equation;
  - the classification `tbl-durham` and its dictionary `tbl-boundary-dictionary`;
  - the soliton-preserving boundary condition;
  - dynamical boundaries;
  - Θ on the half-line.
- **Chapter 31** (`31-half-line-solitons.qmd`, `sec-half-line-solitons`). It contains:
  - the method of images, with mirror shifts, time delays and virtual collision points;
  - solitons returning as antisolitons;
  - ΘP mirrors;
  - boundary breathers;
  - $c_n$ by folding;
  - real-coupling vacua from continued stationary solitons;
  - the link to the missing-charge problem;
  - linearized reflection factors.

Find any object with `grep -rn '{#<id>}' part*/`. Book conventions are in `sec-conventions` (ch. 1), ch. 9 (`eq-imaginary-potential`, `eq-lax-pair`), ch. 10 (`eq-hirota-ansatz`, `eq-an-soliton`, `eq-n-soliton`) and ch. 11 (Θ). Both chapters translate Delius's φ_D = iβφ and his ε = −C, and the relations ξ_book = −ξ_D and C_j = −A_j (Delius 98b) = C_j (Delius–Gandenberger).

The sources are mostly arXiv papers; fetch them with `curl -sL https://arxiv.org/pdf/<id>`:
- Bowcock–Corrigan–Dorey–Rietdijk hep-th/9501098 (BCDR);
- Corrigan–Dorey–Rietdijk–Sasaki hep-th/9404108; Corrigan–Dorey–Rietdijk hep-th/9407148;
- Bowcock–Corrigan–Rietdijk hep-th/9510071; Bowcock hep-th/9609233;
- Delius hep-th/9807189 and hep-th/9809140;
- Delius–Gandenberger hep-th/9904002 (§4–5);
- Baseilhac–Delius–George nlin/0201007; Baseilhac–Koizumi hep-th/0208005.

Experience in this repository is that every adversarial referee pass finds real errors. The drafting itself corrected several published statements, so look hard at the corrections too.

## Already verified (do not redo unless you find a reason)

**By script** (`textbook_code/`, see its README):

- `boundary_classical.py`:
  - **`tbl-durham`.** An invertible field-independent $K$ solving `eq-boundary-k-equation` exists exactly for the listed $C_j$ in the checked rows: $a_1,a_2,a_3,d_4,c_2,c_3,b_3,g_2,a_2^{(2)},a_4^{(2)},a_6^{(2)},d_3^{(2)}$.
  - **Controls.** These fail: without the $\sqrt{n_j}$ factor, or with $|A_i|=\sqrt{2n_i}$.
  - **Coupling.** Imaginary coupling gives the same classification.
  - **φ = 0.** It solves the $|C|=1$ condition for $a_n$ with equal $C_j$, but for no $d_4$ choice.
  - **`prp-soliton-preserving`.** σ = † with $K=1$ is equivalent to Dirichlet on Re φ (in $(2\pi/\beta)\Lambda_W^\vee$) together with Neumann on Im φ at imaginary coupling. At real coupling the two are swapped.
  - **`prp-boundary-theta`.** Checked for $a_2$.
  - **Dictionaries.** Corrigan–Delius, GZ, Saleur–Skorik–Warner, Delius 98a/98b, Delius–Gandenberger and BDG.
  - **BDG equations.** These are reproduced, and energy conservation holds.
- `boundary_images.py`:
  - **Mirror pairs.** For $a_2$ and $a_3$, mirror pairs satisfy the boundary condition at $C=0,\pm1$ for 13 times, with residual about $10^{-14}$. A same-species mirror fails (`prp-antisoliton-return`).
  - **Asymptotics.** The fitted $\Delta t$ and $d$ equal `eq-boundary-time-delay` and Delius (1.4), (1.7), (1.8).
  - **Coweight shifts.**
  - **ΘP mirror** (`eq-thetap-mirror`).
  - **Boundary breathers.** The regularity criterion needs $u^2>\cot^2$ and $|\sinh\rho_-|$, which corrects Delius (3.36) and (3.37). $C=+1$ breathers vanish at $x=0$; $C=-1$ has no regular breathers.
  - **Two pairs** (`eq-multi-pair-shift`).
  - **$c_2$ by folding.**
- `boundary_vacua.py`:
  - **`prp-dg-weights`.** DG Theorem 4.1 holds for $n\le7$.
  - **`tbl-missing-vacua`.** The $a_3$ alternating patterns have no single-soliton vacuum.
  - **`prp-solitonic-vacua`.** The continued static pair is real, regular and solves the field equation, satisfying exactly one sign pattern with $\prod C_j=1$.
  - **`eq-vacuum-energy`.** $E=(2m/\beta^2)(2k-hm_a/m)$, independent of the modulus.
  - **`prp-self-conjugate-vacuum`.** The static species-$h/2$ soliton satisfies $C\equiv-1$ for every $\xi$, with energy 0. Sympy checks $a_1,a_3,a_5$; quadrature checks $a_3$.
  - **`eq-linear-reflection`.** $K_a^{\rm cl}=-[(a)(h-a)]^{-C}$, matching CDR (4.7), (4.9), (5.3). The $B\to0$ limit of DG (3.16) equals the $C=+1$ value.

**By hand** (worth a second look):
- `prp-sklyanin` part (a), the conservation of τ, is proved in the text.
- The proofs of `prp-soliton-preserving` and `prp-boundary-theta`.
- The general proof of `prp-self-conjugate-vacuum`.
- Θ maps solutions of the oscillator system to solutions.
- The "why" of soliton→antisoliton (parity inverts $E$, so $\omega^{ja}\to\omega^{-ja}$).
- Delius 98a (4.4) is the energy of the singular partner. This comparison with Bowcock's $E_\pm$ was done by hand.

## To be checked

Report on each item: correct / wrong (with the correction) / unclear (with what would settle it). Items are ordered roughly by risk.

1. **Corrections to the sources.** Each is stated in the text; confirm or refute it independently.
   - BCDR App. D for $g_2^{(1)}$: is "$C$ on the short root arbitrary" really false?
   - BCDR App. D for $c_2^{(1)}=b_2^{(1)}$: is there an extra branch, and does $a_3^{(2)}\cong d_3^{(2)}$ satisfy the union of the two rules?
   - BCDR (1.4) and (3.18): the reading $|A_i|=2\sqrt{n_i}$. Check it against the journal version, not the arXiv text extraction.
   - Delius 98a (1.5): it is integrable only for $a_n$ and $c_n$.
   - Delius 98a (3.36) $\tan^2\to\cot^2$ and (3.37) $\cosh\to\sinh$.
   - Delius 98a, p. 7: the energy argument for attractive versus repulsive boundaries.
   - Delius 98a (4.4): the vacuum energy belongs to the singular partner.
   - Delius 98a §1.3: "no vacuum unless $\prod A=(-1)^{n+1}$" fails for $C\equiv-1$ with $h$ odd.
   - Delius 98b (1.4): the coweight must lie in $2\Lambda^\vee$.
2. **Status labels and new claims.**
   - `prp-self-conjugate-vacuum` [Theorem] and its interpretation as a vacuum family linked to CDR's zero mode $\omega_{h/2}=0$. This interpretation is the drafting agent's own. Is it right?
   - `prp-soliton-preserving` [Theorem]: is the proof complete, including the claim that $K=1$ solves the classical reflection equation for σ = †? That part is quoted from 98b.
   - `prp-durham`: is [Numerical] the right label, given the mix of checked and quoted rows?
3. **Unchecked rows of `tbl-durham`.** $e_{6,7,8}^{(1)}$, $f_4^{(1)}$, $e_6^{(2)}$, $d_4^{(3)}$, and $a_{2n-1}^{(2)}$ for $n\ge3$. Extend `boundary_classical.py`'s K-solve to them, or explain why the App. D conditions are sufficient.
4. **Sklyanin's formalism** (`prp-sklyanin` (b), `eq-classical-reflection`). The involution property and the classical reflection equations are quoted from Delius 98b (2.7), (2.8). Check them, including the σ-twisted r-matrices. Also check that the Durham $K$ solves the classical reflection equation.
5. **Real-coupling vacua.**
   - The lowest-energy assumption, which neither Bowcock nor Delius proves.
   - Bowcock's multi-soliton vacua (his §3.2).
   - The classical reflection factors about the solitonic vacua: Bowcock (4.21) and DG (5.7).
   - Delius's $c_n$ real-vacuum condition $C_0C_n=1$.
   - Is `tbl-missing-vacua` right, and is the link to the missing-charge problem (`sec-missing-charges`) stated correctly?
6. **Boundary breathers.** `prp-boundary-breathers` was checked only for $a_2$, $a=1$. Does the criterion hold for other species and ranks?
7. **The heuristic** "σ selects the charges, so the species is conjugated or preserved". It is labelled heuristic in the text. Can it be made precise?
8. **The dictionary `tbl-boundary-dictionary`.** Check every row against its source:
   - Saleur–Skorik–Warner $(\zeta,\eta)$; their time delay (2.15) is unchecked;
   - the BDG oscillator variables;
   - the background-field conditions of Bowcock–Corrigan–Rietdijk (unchecked);
   - Baseilhac–Koizumi.
9. **Exercises.** For every `exr-` in the two chapters: is it well posed and solvable with the text, and is any hint wrong?
10. **Bibliography.** Check authors, titles, journal data and DOIs in `references.bib`, and that each citation supports its sentence. The keys added for these chapters are `cdrs1994`, `cdr1995`, `bcdr1995`, `bcr1996`, `bowcock1998`, `delius1998a`, `delius1998b`, `delius_gandenberger1999`, `bdg2002`, `baseilhac_koizumi2003`, `ssw1995`, `gz1994`, `cherednik1984` and `sklyanin1988`.
11. **Prior art.** Are any corrections presented as found here already in the literature? Examples are the $g_2$ entry, the regularity criteria and the energy (4.4). Search the citations of BCDR and Delius 98a (INSPIRE).

## Deliverable

One entry per item: verdict, justification, and for errors the corrected text ready to paste, with file and line. Include any scripts you wrote, inline or as a list of what they check with their output. Do not mark anything verified on the strength of a citation alone.
