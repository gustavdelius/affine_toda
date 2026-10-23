# Brief: particle reflection amplitudes beyond $a_n^{(1)}$ — $d_n^{(1)}$, $c_n^{(1)}$ and the rest (half-line open problem 2)

You are a mathematical physicist with expertise in affine Toda theory, quantum affine algebras and the boundary bootstrap. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it.

## The problem

Open problem 2 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`):

> **Beyond $a_n^{(1)}$.** The particle amplitudes are known only for $a_n^{(1)}$ (and sinh-Gordon). For $d_n^{(1)}$, soliton reflection matrices of the vector solitons are now available (`sec-dn-vector-k`), but their scalar factors and the breather bootstrap are not. For $c_n^{(1)}$ continued vacua exist for some patterns with $C_0C_n=1$ (`sec-images-cn`), and for those the method of this chapter should apply. The exceptional and twisted cases are untouched. Generic solutions of the boundary bootstrap for all affine Toda theories were proposed by Castro-Alvaredo and Fring; which boundary conditions they describe is not known.

This is a programme. The brief asks for its first two steps and a map of the rest.

## What the book already has

- $d_n$: the vector R-matrix (`eq-dn-vector-r`, `prp-dn-r`), and three reflection matrices, $\hat\epsilon=0$ (anti-diagonal) and $s=\pm1$ (`prp-dn-k`). It is not known which classical $d_n$ boundaries they describe. $\vec\phi=0$ satisfies none of the $|C_j|=1$ conditions for $d_n$ (`sec-dn-vector-k`, last paragraph).
- $c_n$: folding of $a_{2n-1}$ images (`sec-images-cn`). Only $C\equiv-1$ (the self-conjugate species-$n$ soliton) has a continued vacuum for $c_2$; Delius's claim "whenever $C_0C_n=1$" fails for the alternating patterns.
- The classical boundary table `tbl-durham`, including the continuous parameters of $b_n$, $c_n$, $a_{2n}^{(2)}$, $d_{n+1}^{(2)}$.
- Scripts: `code/reflection_matrices.py` (the $d_n$ R and K), `boundary_classical.py`, `boundary_images.py`, `boundary_vacua.py`, `particle_reflection.py`.

## Literature

- Castro-Alvaredo–Fring, hep-th/0308004 (`castro_alvaredo_fring2004`): "universal" generic boundary bootstrap solutions for all ATFTs, simply- and non-simply-laced, in several representations.
- Kim, hep-th/9506031 (`kim1995`): $d_4^{(1)}$ Neumann boundary at one loop, a "non-minimal" exact proposal, and the cancellation of non-meromorphic terms for multi-mass spectra. This is a one-loop check for a $d_n$ Neumann amplitude.
- Riva, hep-th/0102163 (NPB 604 (2001) 511) and thesis hep-th/0106268 (not in `references.bib`): boundary bound-state structure of the $E_n$ theories from two known families of reflection amplitudes.
- Fateev–Onofri, hep-th/0207152 (IJMPA 18 (2003) 879; not in `references.bib`): **boundary $B_r$, $B_r^\vee$, $C_r$, $C_r^\vee$, $BC_r$ affine Toda theories with continuous boundary parameters $(w_0,w_r)$**, from boundary one-point functions. They propose minimal reflection amplitudes (Fourier-integral form, their (66)–(71)), boundary ground-state energies, and classical vacua as $\tau$-functions (their (81), (93)). They also state that dual pairs $G$, $G^\vee$ of classically different boundary theories are one quantum theory. Only positive boundary signs are treated. This is directly relevant to the $c_n$ row of `tbl-durham` with its free long-end parameters.
- Fateev, hep-th/0103014 (not in `references.bib`): boundary one-point functions in affine Toda theories and conjectured boundary ground-state energies.
- Ahn–Kim–Rim (`akr2002`) and their earlier non-simply-laced and $ADE$ reflection-amplitude papers (hep-th/9907072, hep-th/0002213): the UV quantum-mechanics method; the $D$, $E$ cases are said to be blocked by CDD ambiguities.
- Delius–George (`delius_george2002`): the $d_n$ coideal and $K$.

## Suggested plan of attack

1. **$d_n$, $\hat\epsilon=0$.** Construct the scalar factor of the anti-diagonal vector reflection matrix by the method of `sec-an-diagonal-k`: unitarity and crossing-unitarity with the $d_n$ vector–vector scalar factor (find or derive it from the Part VI $d_n$ soliton S-matrix in book conventions). Bootstrap the breather in the singlet of $V\otimes V$ (the pole $\langle n-1\rangle$ of `eq-dn-vector-r`) and continue to real coupling. Compare with Castro-Alvaredo–Fring's $d_n$ solutions and, for $d_4$, with Kim's one-loop Neumann result. Use the $B\to0$ limit to decide which classical boundary it is: the linearized reflection factors about each candidate vacuum, which first have to be found, since $\vec\phi=0$ is not one for $|C_j|=1$.
2. **$d_n$, $s=\pm1$.** Repeat if time permits. By analogy with $a_n$, $s=+1$ may be Neumann. Kim's $d_4$ Neumann one-loop amplitude is the test.
3. **$c_n$ via folding.** For the $c_2$ patterns with a continued vacuum ($C\equiv-1$), fold the $a_3$ solitonic amplitudes (`eq-solitonic-factor`, $b=2$, the self-conjugate class) to $c_2$ particle amplitudes. Check them against all axioms with the $c_2$ bulk S-matrix of Part IV. Compare with Fateev–Onofri's $C_r$, $C_r^\vee$ amplitudes after translating their $(w_0,w_r)$ to the book's $C_0,C_n$ (their duality statement implies that the $c_n$ and $d_{n+1}^{(2)}$ boundaries are one quantum family).
4. **Map of the rest.** For $b_n$, $g_2$, $f_4$, $e_{6,7,8}$ and the twisted algebras, record what exists: Castro-Alvaredo–Fring solutions, Riva's $E_n$ analysis, Fateev–Onofri's non-simply-laced families, Kim's $d_4$. Record the missing ingredient in each case (soliton $K$, scalar factor, vacuum, or identification). Produce a table to replace the prose of problem 2.

## Known traps

- Twisted-algebra coupling normalization: the literature uses longest root squared length $2k$, so $\beta_{\rm lit}^2=\beta^2/k$.
- The $d_n$ R-matrix as printed by Delius–George has four misprints (`sec-dn-vector-k`); use the book's script.
- Delius–George's shifted centre $\eta$ is a normalization effect; use book generators.
- The $c_n$ boundary normalization: $\sqrt{2n_j/\vec\alpha_j^2}=2$ on short roots (`sec-images-cn`).
- Fateev–Onofri's $b$ and $x=b^2/(1+b^2)$ are not the book's $B$ for non-simply-laced algebras; match the bulk S-matrices first.

## Deliverable

- The $d_n$ $\hat\epsilon=0$ scalar factor, breather amplitudes, continued particle amplitudes, all checks, and the proposed identification of the boundary.
- The $c_2$ amplitudes and the comparison with Fateev–Onofri.
- The table of step 4.
- Proposed text for open problem 2 and, where results are firm, for new subsections in ch. 33 and 35.
- Scripts.

## Feasibility (assessment when the brief was written, 2026-10-11)

**Medium** for step 1: the scalar factor and breather bootstrap are routine, but the identification of the classical boundary is new work. **Medium** for step 3: folding is mechanical, and Fateev–Onofri gives an independent comparison. **Low** for a complete answer across all algebras.
