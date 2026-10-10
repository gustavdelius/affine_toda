# Brief: referee chapters 32–33 (quantum-group symmetry on the half-line; soliton reflection matrices)

You are a mathematical-physics referee. The repository is `/home/gustav/Git/affine_toda` (branch `claude/affine-toda-boundary-dae69e`, PR #6), a Quarto book on affine Toda field theory. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files.

## Context

Part VII, "Affine Toda theory on the half-line" (`part7-half-line/`, chapters 28–35), was drafted on 2026-10-10 by agents working from the boundary literature. This brief covers two chapters.

**Chapter 32** (`32-boundary-quantum-groups.qmd`, `sec-boundary-qg`):
- Boundary conformal perturbation theory, to the extent needed.
- The conserved charges, and the coideal generators $b_j=e_j+q_j^{-1}f_jk_j+\hat\epsilon_jk_j$ (`eq-coideal-generators`) with the left-coideal coproduct $\Delta(b_j)=(b_j-\hat\epsilon_jk_j)\otimes1+k_j\otimes b_j$ (`eq-coideal-coproduct`).
- $K$ as an intertwiner, and `thm-k-reflection`: intertwiners solve the reflection equation.
- Boundary bound states as coideal modules.
- The sine-Gordon worked example.
- $B=\bar LKL$ and the reflection-equation algebra.
- Quantum versus classical boundary parameters.

**Chapter 33** (`33-soliton-reflection.qmd`, `sec-soliton-reflection`):
- XXZ reflection matrices of any spin (Delius–Nepomechie), and boundary fusion.
- The diagonal $a_n$ $K$ and its scalar factor (Delius–Gandenberger).
- The non-diagonal $a_n$ $K$ (corrected Delius–MacKay (4.21); Gandenberger's $K^\pm$).
- The $d_n$ vector $\check R$ from the tensor-product graph, and the Delius–George $K$.
- Soliton-preserving versus soliton-conjugating reflection (`tbl-k-types`).
- A semiclassical check.

**Conventions.** These chapters must follow Part VI conventions: the coproduct $\Delta(e)=e\otimes1+k\otimes e$, $\Delta(f)=f\otimes k^{-1}+1\otimes f$ (`sec-qg-coproduct`), the physical $q=-e^{-i\pi\omega}=e^{-4\pi^2i/\beta^2}$ (`eq-q-physical`), and the principal gradation. The sign of $q$ is not a convention, which `sec-crossing-sign` explains. Find any object with `grep -rn '{#<id>}' part*/`.

**Sources** (fetch with `curl -sL https://arxiv.org/pdf/<id>`):
- Delius–MacKay hep-th/0112023 (DM)
- Delius–Nepomechie hep-th/0204076 (DN)
- Delius–George math/0208043
- Delius math/0212317; Delius–George hep-th/0212300 (expositions)
- Gandenberger hep-th/9806003 (G98) and hep-th/9911178
- Delius–Gandenberger hep-th/9904002 §3 (DG)
- de Vega–González-Ruiz hep-th/9211114, hep-th/9306089
- Kulish–Sklyanin hep-th/9209054
- Mezincescu–Nepomechie (fusion, J. Phys. A 25 (1992) 2533)
- Doikou hep-th/0006197
- Delius–MacKay–Short hep-th/0109115
- Kim hep-th/9506031

Experience in this repository is that every adversarial referee pass finds real errors. These chapters overturn several published statements, including some by the book's author. Look hardest at those.

## Already verified (do not redo unless you find a reason)

**By script** (`textbook_code/`, see its README).

`boundary_coideal.py` (37 checks):
- The coproduct identity, and the parity-symmetric normalization $q^{-1}fk=k^{1/2}fk^{1/2}$.
- The $d_n$ vector representation satisfies q-Serre.
- $\bar V$ is the antisoliton at the same rapidity only at the physical $q$: $V(z)^*=\bar V(-z/q)$ with $C=1$.
- DM (4.26) holds on (4.13) only after $T_i\to T_i/2$.
- Sine-Gordon: a nonzero $K$ exists iff $c_0c_1=q^{-2}$; $K$ is unique (`prp-sg-k-unique`).
- $a_n$ for $n=2,3,4$ (`prp-an-k-existence`):
  - $K:V(y)\to\bar V(1/y)$ exists iff $\hat\epsilon=0$ (then $K=1$) or all $\hat\epsilon_j^2=1/((1-q)(1-q^{-1}))$; the case $|\hat\epsilon_j|=1$ fails;
  - the closed form `eq-an-k-matrix`;
  - comparison with DM (4.21) under the dictionary, and gauge equivalence of the sign patterns;
  - the twisted reflection equation `eq-reflection-twisted` for all $(\mu,\nu)$, with the constant equal to 1;
  - the printed $K$ fails at $-q$;
  - $s=\pm1$ equal Gandenberger's $K^\pm$;
  - crossing-unitarity of the matrix part.

`boundary_sine_gordon.py`: the coideal $K$ equals GZ/BPT; the Neumann point is $\hat\epsilon=(2-q-q^{-1})^{-1/2}$.

`reflection_matrices.py` (38 checks):
- DN:
  - the closed form solves their linear equations for $j\le2$ at real and complex $\eta$;
  - at physical imaginary $\eta$ with principal branches it fails for $j\ge1$, as their footnote anticipates;
  - their linear equations equal the book intertwiner under `eq-xxz-dictionary`;
  - the spin-$j$ $K$ solves the reflection equation with $R^{(j,j)}$, $R^{(1/2,j)}$ and $R^{(1,2)}$.
- Boundary fusion (`prp-boundary-fusion`).
- $d_4$, $d_5$:
  - the TPG eigenvalues of $\check R$ (`prp-dn-r`, with $q\to1/q$ for the Part VI coproduct);
  - the Delius–George existence conditions, plus the $\hat\epsilon=0$ solution their (22) omits, and $\hat\epsilon_{DG}^2=\epsilon_*^2$;
  - the book coideal $K$ exists iff $\hat\epsilon=0$ (anti-diagonal) or all $\hat\epsilon_j^2=\epsilon_*^2$;
  - the reflection equation holds;
  - torus conjugates also solve it.
- DG (3.9): the scalar factor satisfies unitarity and (3.8) (`prp-an-diagonal-scalar`, $n=2,3,4$).
- DG (3.24): the CDD factors need the overall minus sign for even $n$.
- Semiclassical limit of $A_1$ against the $(+\dots+)$ time delay (`prp-semiclassical-uniform`).

**By hand** (worth a second look):
- The proof of `thm-k-reflection` under its stated hypotheses (existence, B-irreducibility of $V_\mu\otimes V_\nu$, normalization).
- `prp-boundary-module`.
- `prp-b-coideal`, re-proved in book conventions via `eq-quasitriangular`. With the Part VI coproduct, the universal $R$ is the flip of ch. 7's.
- The resolution of DM's "arbitrary diagonal $K$ at $\hat\epsilon=0$" as torus conjugates intertwining a rescaled coideal.

## To be checked

Report on each item: correct / wrong (with the correction) / unclear (with what would settle it). Items are ordered roughly by risk.

1. **Neumann versus uniform boundary** (`cnj-neumann-coideal`, `prp-dm-first-order`, `sec-epsilon-constraint`). The chapters claim:
   - Neumann is $\hat\epsilon_j=\epsilon_*=(2-q-q^{-1})^{-1/2}$ with $s=+1$, not $\hat\epsilon=0$;
   - the diagonal $K$ ($\hat\epsilon=0$) belongs to the uniform $C\equiv+1$ boundary;
   - DM's first-order boundary perturbation fixes only the shift of $\hat\epsilon$ away from the Neumann value.

   The evidence offered:
   - sine-Gordon (GZ's free boundary, and BPT's Neumann point);
   - G98's $O(\beta^2)$ comparison with Kim (Neumann) and Perkins–Bowcock (+++);
   - the semiclassical limit for $(+\dots+)$;
   - the classical limit of DG (3.16), which is in ch. 31.

   Is the argument sound? Is the proposed mechanism (DM's normalization $c\propto$ bulk coupling, so an $\epsilon_*q^T$ term arises at first order in the bulk coupling) right? Can it be computed? Is [Conjecture] the right label, or is the identification established in later literature?
2. **The DM corrections** listed in ch. 32's notes; confirm or refute each:
   - their $q$ is inverted;
   - the factor 2 in §4;
   - (4.27) prints $\Delta(Q)$ twice;
   - (4.19)–(4.20) should read $\hat\epsilon(q-1)$, with two equations per row;
   - the condition is $\hat\epsilon^2=\epsilon_*^2$, so the residual freedom is a sign, not a phase;
   - "arbitrary diagonal" contradicts (4.18);
   - (4.21) is correct after the dictionary.
3. **The other source corrections** in ch. 33's notes:
   - DN (3.8): the square root covers the first product only; the branch caveat;
   - DG (3.9): the sign of the prefactor block is not fixed by (3.8);
   - DG (3.24): the minus sign;
   - Delius–George (22) omits the $\hat\epsilon=0$ solution; their $\eta^2=q^{-h}$ is a rescaling artefact; their continuous families are torus conjugations.
4. **`thm-k-reflection`** [Theorem]. Is the proof complete? Pay attention to:
   - the irreducibility of tensor products under the coideal for generic rapidities, which is assumed; is it known, and for which representations?
   - the normalization argument;
   - the twisted case with four different $\check R$'s.
5. **`eq-coideal-generators`.** Does the choice $c_j=q_j^{-1}$ (the image of $\bar Q_j$) follow from the charges of DM §3–4 once translated to Part VI conventions? Is the identification of $\hat\epsilon_j$ with the boundary coupling of ch. 30 consistent with `eq-boundary-potential`?
6. **The general-$n$ claims** of `prp-an-k-existence` (a sketch beyond $n=4$) and of the $d_n$ results (checked for $n=4,5$). Also: which classical $d_n$ boundaries correspond to $\hat\epsilon=0$, $s=+1$ and $s=-1$? This is open.
7. **$B=\bar LKL$** (`prp-b-coideal`, `eq-b-matrix`). The $O(x)$ expansion that recovers the generators (DM (4.35)ff.) is quoted, not redone. Redo it in book conventions, or say what fails.
8. **Unchecked entry by entry:** the Kulish–Reshetikhin $R^{(1/2,j)}$ of DN (2.11), Jimbo's/Delius–George's (16)–(17) $d_n$ $\check R$, Gandenberger's (3.3)–(3.4) in his own conventions, and de Vega–González-Ruiz's diagonal soliton-preserving $K$.
9. **`prp-boundary-module`.** Are the $\beta_j$ of the Bajnok–Palla–Takács–Tóth excited boundary states consistent with it?
10. **`tbl-k-types` and the soliton-preserving case.** Is it right which reflection equation (twisted or untwisted) goes with which classical boundary condition of ch. 30? Is Doikou's terminology quoted correctly?
11. **The semiclassical check** `prp-semiclassical-uniform`. Is the comparison set up correctly (which time delay, which limit $\omega\to\infty$)? Is the stated $1/\omega$ convergence the expected rate?
12. **Exercises.** For every `exr-` in the two chapters: is it well posed and solvable with the text, and is any hint wrong?
13. **Bibliography.** Check authors, titles, journal data and DOIs in `references.bib`, and that each citation supports its sentence:
    - `delius_mackay2003`, `delius_nepomechie2002`, `delius_george2002`, `delius_george2003`, `delius2002`, `dms2001`
    - `gandenberger1999a`, `gandenberger1999b`
    - `devega1993`, `devega1994`
    - `kulish_sklyanin1992`, `mezincescu_nepomechie1992`, `mezincescu_nepomechie1998`
    - `doikou2000`, `kim1995` (published 1996)
14. **Prior art and later work.** Has the $a_n$ existence condition, the Neumann identification or the omitted $d_n$ solution already appeared in later work? Search the citations of DM and Delius–George (INSPIRE), including the quantum-symmetric-pair literature (Kolb, Balagović–Kolb, Regelskis–Vlaar), which may already contain the corrected parameter conditions.

## Deliverable

One entry per item: the verdict, the justification, and for errors the corrected text ready to paste, with file and line. Include any scripts you wrote, inline or as a list of what they check with their output. Do not mark anything verified on the strength of a citation alone.
