# Brief: referee pass on chapters 5 and 6 (root systems, affine Kac–Moody algebras)

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (worktree `.claude/worktrees/book-preface-section-reorganize-7f6a81`, branch `claude/referee-briefs-part-ii-bac6f0`).

**Rules.**
- Other agents may be editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Do not spawn subagents.
- Write scratch scripts only under `referee_scratch/algebra_05_06/`.
- Write the report to `briefs/referee_algebra_05_06_report.md`.

## Context

Part II, "Algebraic background" (`part2-algebra/`), is the book's account of the Lie theory, affine algebras, quantum groups and R-matrices that the later parts use. This brief covers its first two chapters:

- **Chapter 5**, `part2-algebra/05-root-systems.qmd` (`sec-root-systems`): roots, Cartan matrix, Dynkin diagrams and `tbl-simple-lie`, the Weyl group and Coxeter element, weights and lattices, the Perron–Frobenius eigenvector and its link to masses and conserved charges, `fig-root-systems`.
- **Chapter 6**, `part2-algebra/06-affine-algebras.qmd` (`sec-affine-algebras`): loop algebras, affine Dynkin diagrams and `tbl-untwisted`, twisted algebras, folding and duality (`sec-affine-folding`, `tbl-foldings`, `fig-folding`), gradations and the principal Heisenberg subalgebra, vertex-operator representations, and the summary `sec-lie-algebraic-data`.

A companion brief, `briefs/referee_algebra_07_08.md`, covers chapters 7 and 8. Stay within chapters 5 and 6, but report any inconsistency with other chapters that you meet.

**History.** The chapters were drafted on 2026-10-09 (commit 4776db5). One referee pass was run on 2026-10-10 by **Claude agents only**, as part of a book-wide pass: brief `briefs/referee_textbook.md`, items 9, 10 and 12; report `briefs/referee_textbook_report.md`, sections "Item 9", "Item 10", "1c" and "Item 12". Its replacement text was integrated in commit 2e306b6. A reading pass then rewrote several passages (commit 790bedf; brief and resolution in `briefs/part2_reading_fixes.md`), and figures were added (commit ffab947). The chapter text is unchanged since main at 9d85ab7.

This is therefore the first **model-diverse** pass on Part II, and much of the current text has had **no referee pass at all**:
- the integrated replacement text of 2026-10-10, which the Claude referee wrote and nobody else checked: the duality bullets, the folding normalization sentence, the `tbl-foldings` caption, the period $kh$ of the exponents, the affine-diagram symmetry sentence, the Olive–Turok–Underwood formula;
- the text added by the reading pass: the paragraph "The parent of $a_{2n}^{(2)}$" and the new row of `tbl-foldings`, the rewritten duality bullet with $\beta^\vee=4\pi\sqrt k/\beta$, the spin-$(m_k+h)$ sentences in `sec-perron-frobenius`, the Killing-form normalization, the $\mathbb Z_{n+1}$ charge sentence;
- the figure captions of `fig-root-systems` and `fig-folding`.

Weight your effort toward these. Neither chapter has labelled statements, items in Appendix E (original results) or Appendix D (corrections), and the text cites no check script.

**Check scripts.** None is cited in the chapters. These scripts check some of their claims; each prints one PASS or FAIL line per check, and each runs in about a second:
- `textbook_code/figures/ch5_f10_root_systems.py` (13 checks): $G_2$ data, $E_8$ Coxeter orbits and radii.
- `textbook_code/figures/ch6_f11_folding.py` (8 checks): the $c_2^{(1)}$ and $g_2^{(1)}$ folds.
- `textbook_code/algebra_data.py` (6 checks): classical masses of the untwisted theories, by folding.

The figure scripts overwrite the committed SVGs in `part2-algebra/figures/` whatever `OUT_DIR` says. If you run them, report it, and do not restore the files yourself; the editor will.

## Notation

- **Root normalization (book).** The longest roots of the finite algebra have squared length 2 (`sec-conventions` in `01-introduction.qmd`, `sec-lie-algebraic-data`). For a twisted algebra $X^{(k)}$, much of the literature instead gives the longest root squared length $2k$; then $\beta_{\rm lit}^2=\beta^2/k$. Chapter 7 uses Kac's normalization, $(\alpha_0,\alpha_0)=2$, for $U_q(\hat g)$ (the other brief).
- **Coroots** $\vec\alpha^\vee=2\vec\alpha/\vec\alpha^2$. **Kac labels** $n_j$, with $\sum_jn_j\vec\alpha_j=0$ and $n_0=1$ (`eq-kac-labels`), $h=\sum_jn_j$. **Dual labels** $n_j^\vee=n_j\vec\alpha_j^2/2$ for untwisted algebras, $n_j^\vee=k\,n_j\vec\alpha_j^2/2$ for twisted ones (`tbl-foldings` caption), $h^\vee=\sum_jn_j^\vee$.
- **Potentials** (`eq-imaginary-potential`, ch. 9): $V_{\rm real}=\frac{m^2}{\beta^2}\sum_jn_j(e^{\beta\vec\alpha_j\cdot\vec\phi}-1)$ and $V_{\rm imag}=-\frac{m^2}{\beta^2}\sum_jn_j(e^{i\beta\vec\alpha_j\cdot\vec\phi}-1)$.
- **Sine-Gordon** ($a_1^{(1)}$): $\beta_{\rm SG}^2=2\beta^2$.
- **Two dualities, easily confused.** $\hat g\to\hat g^\vee$ reverses the arrows of the affine diagram (for example $c_n^{(1)}\leftrightarrow d_{n+1}^{(2)}$). "Lie duality" $X^{(1)}\to(X^\vee)^{(1)}$ dualizes the finite root system (for example $c_n^{(1)}\leftrightarrow b_n^{(1)}$). Chapter 6 uses both; check that each sentence uses the right one.
- **Twisted node labelling.** Several conventions exist for which node of a twisted diagram is node 0, especially for $a_{2n}^{(2)}$. Kac (Table Aff 2) puts the label $a_0=2$ on node 0 of $A_{2l}^{(2)}$, whereas the book's folding gives node 0 label 1 and the longest root. Say which convention each of your computations uses.

## Read first

- `01-introduction.qmd`, section `sec-conventions`.
- Chapter 9, `eq-imaginary-potential` and `sec-theories`; ch. 10, `sec-imaginary-coupling-vacua`; ch. 14, `sec-nsl-duality` and `tbl-floating-h`; ch. 20, `sec-qg-coproduct`; ch. 21, `sec-amplitude-and-gradation`. Only as far as items 1, 2 and 4 need them.
- The 2026-10-10 report sections named above, to see what was checked then. Do not take its verdicts as evidence.
- To find an object: `grep -rn '{#<id>}' part*/ appendices/ outlook/ *.qmd`.

**Sources.** None of the sources of these chapters is in `~/referee-copy/papers/` (or `papers/`). Fetch from arXiv what you need:
- Khastgir–Sasaki, hep-th/9512158 (`khastgir1996b`);
- MacKay–Watts, hep-th/9411169 (`mackay_watts1995`);
- Delius–Grisaru–Zanon, hep-th/9201067 (`dgz1992`);
- Corrigan–Dorey–Sasaki, hep-th/9304065 (`corrigan1993`);
- Kneipp–Olive, hep-th/9404030 (`kneipp1996`);
- Olive–Turok–Underwood (`olive1993a`, `olive1993b`; find by title, 1992–93).

The others are pre-arXiv papers and books: Olive–Turok 1983, Braden–Corrigan–Dorey–Sasaki 1990, Freeman 1991, Dorey 1991, Klassen–Melzer 1990, Steinberg 1959, Kostant 1959, Kac's *Infinite-Dimensional Lie Algebras*. If you cannot read a source, say so and mark that citation "unclear"; do not reconstruct what it probably says.

## How to work

For each item, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_38_40_report.md` shows the expected depth.

In numerical work, avoid values of $q$ or of the coupling that are roots of unity. Short decimals such as 0.3 or 1.37 give a rational $\omega$, and then a root of unity. Use generic values such as $\sqrt2$ or $e/2$.

## Items (most important first)

1. **Duality of the dual pairs, and its normalization** (`sec-affine-folding`, bullet "The real-coupling S-matrices obey a strong–weak duality", and `sec-lie-algebraic-data`).
   - Claim: the theory of $\hat g$ with simple roots $\vec\alpha_j$ at coupling $\beta$ has the same S-matrix as the theory of $\hat g^\vee$ with simple roots $\vec\alpha_j^\vee$ at coupling $4\pi/\beta$. This form is said not to depend on the normalization. In the book's normalization it reads $\beta^\vee=4\pi\sqrt k/\beta$, symmetric between the two members of a pair.
   - Ch. 14 (`sec-nsl-duality`) writes $2-B=B(\beta^\vee_*)$, with $\beta^\vee_*=\beta^\vee$ for $d_{n+1}^{(2)}$, $d_4^{(3)}$ and $e_6^{(2)}$, but $\beta^{\vee2}_*=\beta^{\vee2}/2$ for $a_{2n-1}^{(2)}$. Is that consistent with chapter 6's single rule, or does one of them contain a normalization slip?
   - Computation: take $B(\beta)$ and $H(\beta)$ from `tbl-floating-h`. Check $2-B(\beta)=B^\vee(\beta^\vee)$ and the exchange of $H$ for one pair of each kind: $(c_2^{(1)},d_3^{(2)})$, $(b_2^{(1)},a_3^{(2)})$ and $(g_2^{(1)},d_4^{(3)})$. Derive $\vec\alpha_j^\vee$ and their squared lengths yourself for each pair.
   - Check also the first duality bullet: Lie duality, with the example "for $c_n^{(1)}$ it gives $b_n^{(1)}$, not $d_{n+1}^{(2)}$". Check the second bullet, that the non-local charges generate $U_q(\hat g^\vee)$. Check the closing sentence of `sec-affine-folding` and the summary list of dual pairs in `sec-lie-algebraic-data`.
   - Find the source of each duality statement in DGZ 1992 or CDS 1993 and quote it. The reading pass recorded that it "did not confirm [the coroot form] against DGZ or CDS text".

2. **The parent of $a_{2n}^{(2)}$** (paragraph after `tbl-foldings`, and its row). This text has had no referee pass.
   - Fold $a_{2n}^{(1)}$ by $j\mapsto-j$ and $d_{2n+2}^{(1)}$ by the stated order-4 automorphism, for $n=1,2,3$. In each case, check:
     - that the map is a diagram automorphism;
     - the orbits;
     - the folded squared lengths ($2,1,\dots,1,\frac12$ in both cases);
     - the labels and their common factor;
     - $h=2n+1$;
     - the replacement $m^2\to2m^2$ for $d_{2n+2}^{(1)}$, against no change for $a_{2n}^{(1)}$.
   - Is the folded Cartan matrix that of $a_{2n}^{(2)}$, or of its transpose?
   - Check the claim that both foldings give "the $a_{2n}^{(2)}$ theory with the same classical masses". Compute the folded masses from the parent mass matrices of each, with the stated factor $m^2\to2m^2$.
   - Check the general rule stated before the example: "For the twisted foldings ($\alpha_0$ in an orbit of size $k$) these labels have the common factor $k$ … $m^2$ replaced by $km^2$", and that $a_{2n}^{(2)}$ is the only exception. Do this for every twisted row of `tbl-foldings`.
   - Quote Khastgir–Sasaki on both foldings.
   - Part V uses both foldings (`sec-hirota-form-dn1-en1` in ch. 18, `sec-one-loop-masses-folded` in ch. 19). Check that they label the nodes and normalize $m$ as chapter 6 does, or report the clash.

3. **`tbl-untwisted`, `tbl-foldings` and the normalization of the twisted labels.**
   - Recompute independently every Kac-label row, $h$ and $h^\vee$ of `tbl-untwisted`.
   - For every row of `tbl-foldings`: search the automorphisms of the parent diagram, fold by the rule of `sec-affine-folding`, and compare the folded Cartan matrix with the target and its transpose. Then compute $h$ and $h^\vee$ by the caption's rule $n_j^\vee=k\,n_j\vec\alpha_j^2/2$, and check that this rule gives the coprime integer null vector of the transposed Cartan matrix, as the caption says.
   - Check "duality exchanges the Coxeter number and the dual Coxeter number" for every pair.
   - Check the ranks in `sec-twisted`, and the list of diagram automorphisms ("order 2 for $A_n$, $D_n$, $E_6$"; is $A_1$ an exception?).
   - Find a source for the twisted labels. The text says "we follow [@mackay_watts1995]" (in `sec-twisted` and `sec-lie-algebraic-data`). The 2026-10-10 pass left open whether MacKay–Watts are the right citation, or whether they follow Braden–Corrigan–Dorey–Sasaki or Kac. Quote what MacKay–Watts say.

4. **Gradations, exponents and vertex operators** (`sec-gradations`, `sec-vertex-operators`).
   - Claim: the principal degrees of the principal Heisenberg subalgebra are the exponents of $\hat g$, with period $h$ for $g^{(1)}$ and $kh$ for $g^{(k)}$. Check the examples "$1,5$ mod 6 for $a_2^{(2)}$, $1,5,7,11$ mod 12 for $d_4^{(3)}$". Compute the centralizer of $E_1$ in a matrix realization for $a_2^{(2)}$ and one more twisted algebra.
   - The sentence on $2T=4\pi h/\beta^2-h^\vee$: for which algebras, and in which normalization of $\beta$, does chapter 21 (`sec-amplitude-and-gradation`) state it? Is chapter 6's unqualified statement consistent with that for twisted $\hat g$?
   - Claim: $L(\Lambda_j)$ has level $n_j^\vee$.
   - Claim: level-one representations with $n_j=1$ have vertex-operator realizations, with the attributions "Frenkel–Kac … homogeneous" and "Lepowsky–Wilson … principal". The 2026-10-10 pass noted that the general principal construction is due to Kac–Kazhdan–Lepowsky–Wilson.
   - Compare the $\tau_j$ formula and $W_k$ with Olive–Turok–Underwood, equations (2.27)–(2.29) of the vertex-operator paper or wherever they are. Check the scope: untwisted only, with Kneipp–Olive for twisted algebras.

5. **Symmetries of the affine diagram** (`sec-affine-dynkin`, last paragraph).
   - Claim: the symmetry group is the semidirect product of the finite-diagram symmetries with a subgroup isomorphic to the centre of the simply connected group, which permutes the nodes with $n_j=1$ simply transitively. Check it for every untwisted algebra by enumerating automorphisms; include $d_n$ for $n$ even and odd, and $a_1^{(1)}$.
   - Claim: the $a_n^{(1)}$ particles carry a $\mathbb Z_{n+1}$ charge, and a species-$a$ soliton has topological charge among the weights of $\Lambda^a\mathbb C^{n+1}$, whose class modulo the root lattice is $a\bmod(n+1)$.

6. **Chapter 5, Perron–Frobenius and conserved charges** (`sec-perron-frobenius`).
   - Claim: for every spin $s$ of a conserved charge, the $\chi_a^{(s)}$ of `eq-charge-eigenvalue` (ch. 2) form an eigenvector of the Cartan matrix with eigenvalue $4\sin^2\frac{\pi s}{2h}$.
   - Claim, added by the reading pass: the spin $m_k+h$ gives the eigenvector of the exponent $m_k$ with its sign reversed on one colour of the bicolouring, so that $\sin\frac{a(s+h)\pi}h=(-1)^a\sin\frac{as\pi}h$ for $a_n^{(1)}$. Check this from the three-point conservation laws at Dorey's fusings, not from the Cartan matrix alone. Do it for $a_4$, $d_5$ and $e_6$, for one spin above $h$ each.
   - Check the repeated exponent of $D_{2n}$, and the citations `klassen1990` and `dorey1991` for these statements.
   - Check `eq-cartan-eigenvalues`, the $A_n$ example, the $E_8$ golden ratio, and the Ising sentence, which says that the scaling-Ising S-matrix is the "coupling-independent (minimal) part" of the $e_8^{(1)}$ S-matrix.
   - Check the attribution: Freeman for the mass theorem only, "following" the observation of Braden, Corrigan, Dorey and Sasaki.

7. **Chapter 5, other statements.**
   - The normalization in `sec-roots`: the Killing form rescaled so that long roots have squared length 2, and $[E^{\vec\alpha},E^{-\vec\alpha}]=2\vec\alpha\cdot\vec H/\vec\alpha^2$, against `eq-chevalley-serre` with $h_i=\vec\alpha_i^\vee\cdot\vec H$. Are they consistent, and is the $su(2)$ remark right?
   - `tbl-simple-lie` and `eq-coxeter-numbers`.
   - The Coxeter-element facts: conjugacy, order $h$, eigenvalues, $rh$ roots in $r$ orbits of $h$.
   - Steinberg's statement as quoted, including "$h$-gons of different radii" for non-simply-laced $g$.
   - The coweight lattice and the periodicity of $V_{\rm imag}$.
   - `eq-casimir`; $\vec\rho=\sum_i\vec\lambda_i$.
   - The captions of `fig-root-systems`. Is the shaded wedge of $G_2$ the fundamental chamber for the drawn simple roots? Are the $E_8$ radii proportional to the Perron–Frobenius vector?

8. **Prior art, attributions and citations.**
   - Part II claims no original result, and Appendix E (`appendices/E-original-results.qmd`) has no Part II item. Check that no sentence of chapters 5–6 presents a known result as the book's own, or attributes the book's own work to the literature.
   - Check that every citing sentence says what its source says. Examples: Olive and Turok "also showed that every non-simply-laced theory arises this way"; "the role of duality in the masses and S-matrices was clarified by Delius, Grisaru and Zanon and Corrigan, Dorey and Sasaki"; Kostant's relation of the exponents to the principal three-dimensional subalgebra; `goddard1986` and `kac1990` for the topics named; `khastgir1996b`, `kneipp1996`, `olive1993a`, `olive1993b`.
   - Recall the trap: `dgz1992` is Delius–Grisaru–Zanon, `dgz1996` is Delius–Gould–Zhang.

9. **Exercises and scope.**
   - Work every exercise of both chapters at least in outline: `exr-a2-roots`, `exr-g2`, `exr-d4-masses`, `exr-cartan-eigenvalues`, `exr-roots-count`, `exr-jacobi-central`, `exr-affine-cartan`, `exr-fold-d4`, `exr-fold-mass`, `exr-dual-coxeter`. Is each one true, and solvable with what the chapter provides?
   - Check `exr-fold-mass` for $c_2^{(1)}\subset a_3^{(1)}$, and whether its general claim also holds for the $a_{2n}^{(2)}$ folds of item 2.
   - List any claim in these chapters that goes beyond standard textbook material and therefore needs a status label or a check script.

## Report format (`briefs/referee_algebra_05_06_report.md`)

A header with the date, your model, and the commit you reviewed (`git rev-parse HEAD`). Then one section per item, with:
- the verdict (correct / wrong / unclear), per claim where an item has several;
- your own computation (the script inline or its path, and its output);
- for errors, replacement text ready to paste, with file and line.

Then a list of the claims you could not check, and why. End with other observations, labelled as unverified.
