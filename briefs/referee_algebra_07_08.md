# Brief: referee pass on chapters 7 and 8 (quantum groups, R-matrices)

You are an adversarial mathematical-physics referee. You are working in the repository of a Quarto book on affine Toda field theory (worktree `.claude/worktrees/book-preface-section-reorganize-7f6a81`, branch `claude/referee-briefs-part-ii-bac6f0`).

**Rules.**
- Other agents may be editing files in this worktree at the same time.
- Do not edit any existing file.
- Do not run git commands that change the index or the working tree (no commit, checkout, stash, reset, rebase or restore).
- Do not run `quarto render`.
- Do not spawn subagents.
- Write scratch scripts only under `referee_scratch/algebra_07_08/`.
- Write the report to `briefs/referee_algebra_07_08_report.md`.

## Context

Part II, "Algebraic background" (`part2-algebra/`), is the book's account of the Lie theory, affine algebras, quantum groups and R-matrices that the later parts use. Part VI builds the soliton S-matrices on chapters 7 and 8. This brief covers those two chapters:

- **Chapter 7**, `part2-algebra/07-quantum-groups.qmd` (`sec-quantum-groups`):
  - $U_q(sl_2)$, its Hopf structure and universal R-matrix;
  - $U_q(\hat g)$, its central element, normalizations and evaluation representations (`sec-uqghat`);
  - R-matrices with spectral parameter;
  - the six-vertex model and the sine-Gordon dictionary, with the sign of $q$ (`sec-six-vertex`);
  - crossing from the antipode (`sec-qg-crossing`);
  - Hermitian forms at $\lvert q\rvert=1$ (`sec-qg-hermitian`);
  - roots of unity (`sec-roots-of-unity`).
- **Chapter 8**, `part2-algebra/08-r-matrices.qmd` (`sec-r-matrices`):
  - the history in the introduction;
  - Jimbo's equations and the spectral decomposition;
  - the tensor-product graph, untwisted and twisted, with `fig-tpg`, and where it fails (`sec-tpg`);
  - the direct method (`sec-direct-method`);
  - fusion (`sec-fusion`);
  - vertex and IRF forms, and RSOS restriction (`sec-irf`).

A companion brief, `briefs/referee_algebra_05_06.md`, covers chapters 5 and 6. Stay within chapters 7 and 8, but report any inconsistency with other chapters that you meet.

**History.** The chapters were drafted on 2026-10-09 (commit 4776db5). One referee pass was run on 2026-10-10 by **Claude agents only**, as part of a book-wide pass: brief `briefs/referee_textbook.md`, items 2, 3, 4, 11 and 12; report `briefs/referee_textbook_report.md`, sections "Item 2", "Item 3", "Item 4", "Item 11", "1c" and "Item 12". It found several errors, among them the sign of $q$, the crossing-point paragraph, the twisted tensor-product graph and the central element. Its replacement text was integrated in commit 2e306b6. A reading pass then rewrote several passages (commit 790bedf; brief and resolution in `briefs/part2_reading_fixes.md`), and `fig-tpg` was added (commit ffab947). The chapter text is unchanged since main at 9d85ab7.

This is therefore the first **model-diverse** pass on Part II, and much of the current text has had **no referee pass at all**:
- the integrated replacement text of 2026-10-10, which the Claude referee wrote and nobody else checked: the $q$-sign paragraph of `sec-six-vertex`, the $S^2$ paragraph of `sec-qg-crossing`, the central element, the twisted tensor-product graph and "When it fails", the $\mathbf{27}$ count in `sec-direct-method`, the history in the chapter 8 introduction;
- the text added by the reading pass: the normalization sentences for twisted algebras in `sec-uqghat` and `sec-qg-crossing` (including the exception $a_{2n}^{(2)}$), the tilting-module paragraph of `sec-roots-of-unity`, the Hopf glossary, the gloss of $x_c$, the hints of `exr-spin1` and `exr-faces`;
- the caption of `fig-tpg`.

Weight your effort toward these. Neither chapter has labelled statements, items in Appendix E (original results) or Appendix D (corrections), and the text cites no check script.

**Check scripts.** None is cited in the chapters. These scripts check some of their claims; each prints one PASS or FAIL line per check:
- `textbook_code/sixvertex.py` (8 checks, about 2 s): Jimbo's equations, unitarity, braid Yang–Baxter, eigenvalues, the crossing data, the sine-Gordon dictionary.
- `textbook_code/universal_r.py` (5 checks, about 1 s): `eq-universal-r-sl2` on spin $\frac12$ and 1.
- `textbook_code/spin1_tpg.py` (2 checks, about 6 s): `eq-tpg-rule` and `exr-spin1`.
- `textbook_code/rsos_qdims.py` (2 checks, under 1 s): quantum dimensions at roots of unity.
- `textbook_code/figures/ch8_f12_tpg.py` (6 checks, about 1 s): the data of `fig-tpg`. It overwrites the committed SVG in `part2-algebra/figures/` whatever `OUT_DIR` says. If you run it, report it, and do not restore the file yourself; the editor will.

## Notation

- **Chapter 7 coproduct** (`eq-coproduct`): $\Delta(K)=K\otimes K$, $\Delta(E)=E\otimes K+1\otimes E$, $\Delta(F)=F\otimes1+K^{-1}\otimes F$. Antipode (`eq-antipode`): $S(E)=-EK^{-1}$, $S(F)=-KF$.
- **Part VI coproduct** (`sec-qg-coproduct`, ch. 20): the opposite one, $\Delta(e_i)=e_i\otimes1+k_i\otimes e_i$, $\Delta(f_i)=f_i\otimes k_i^{-1}+1\otimes f_i$. Hence $\check R_{\rm VI}(z;q)=\check R_7(z;q^{-1})$.
- **Physical $q$** (`eq-q-physical`, ch. 20, in Part VI's coproduct): $q=e^{-4\pi^2i/\beta^2}=-e^{-i\pi\omega}$, with $\omega=4\pi/\beta^2-1$. For sine-Gordon, $\beta_{\rm SG}^2=2\beta^2$ and $\lambda=8\pi/\beta_{\rm SG}^2-1=\omega$. In chapter 7's coproduct the same charges give $q^{-1}=-e^{i\pi\lambda}$.
- **The map $q\to-q$** is a Klein-factor twist, not a gauge transformation (`sec-crossing-sign`, ch. 21). "Unchanged under $q\to-q$" usually holds only for moduli or for even soliton number.
- **Spectral parameter.** $z=x/y$ in the homogeneous gradation, where only $e_0$, $f_0$ carry $x$ (`eq-evaluation-rep`). For sine-Gordon, $z=e^{2\lambda\theta}$.
- **Tensor order.** $\check R:V(x)\otimes W(y)\to W(y)\otimes V(x)$ (`eq-r-check`), and $R=P\check R$. The operator orderings of `eq-ybe-braid` and `eq-fusion` are as printed.
- **Root normalization.**
  - Chapters 5–6 and the rest of the book normalize the longest roots of the finite algebra to squared length 2.
  - `sec-uqghat` and `sec-qg-crossing` instead use Kac's normalization, $(\alpha_0,\alpha_0)=2$ (1 for $a_{2n}^{(2)}$), with $q_i=q^{(\alpha_i,\alpha_i)/2}$. For twisted $X^{(k)}$ the longest roots then have squared length $2k$, the literature normalization of `sec-conventions`; there $\beta_{\rm lit}^2=\beta^2/k$.
  - The tensor-product-graph rule `eq-tpg-rule` is stated with long roots of squared length 2.
  - Several conventions exist for which node of a twisted diagram is node 0, especially for $a_{2n}^{(2)}$. Kac puts $a_0=2$ on node 0 of $A_{2l}^{(2)}$; the book's folding in ch. 6 gives node 0 label 1 and the longest root.
  - Say which normalization and labelling each of your computations uses.
- **Bethe–Yang operator.** It is $R$ in some sources and $P\check R$ in others; check which one a formula uses.

## Read first

- `01-introduction.qmd`, section `sec-conventions`.
- Chapter 20, `sec-qg-coproduct` and `eq-q-physical`; chapter 21, `sec-amplitude-and-gradation` and `sec-crossing-sign`. These hold the Part VI side of the dictionary.
- Chapter 6, `sec-twisted`, `sec-affine-folding` and `tbl-foldings`, for the twisted algebras and their labels.
- The 2026-10-10 report sections named above, to see what was checked then. Do not take its verdicts as evidence.
- To find an object: `grep -rn '{#<id>}' part*/ appendices/ outlook/ *.qmd`.

**Sources.** The texts of the arXiv sources are in `~/referee-copy/papers/` (or `papers/`, if you run through `run_referee_agy.sh`), as `<id>.txt` with `/` replaced by `_`:
- Delius–Gould–Zhang, `hep-th_9405030.txt` (`dgz1994`) and `q-alg_9508012.txt` (`dgz1996`);
- Gandenberger–MacKay–Watts, `hep-th_9509007.txt` (`gmw1996`);
- Delius, `hep-th_9503079.txt` (`delius1995`);
- Miramontes, `hep-th_9901145.txt` (`miramontes1999`);
- Hoare–Hollowood–Miramontes, `1303.1447.txt` (`hoare2013`).

The PDFs are next to them; consult a PDF when the extracted text garbles a formula. The other sources are pre-arXiv papers and books: Jimbo 1985 and 1986, Drinfeld 1987, Kulish–Reshetikhin–Sklyanin 1981, MacKay 1991, Zhang–Gould–Bracken 1991, Frenkel–Reshetikhin 1992, Bernard–LeClair 1991, Chari–Pressley, Kassel, Gómez–Ruiz-Altaba–Sierra. If you cannot read a source, say so and mark that citation "unclear"; do not reconstruct what it probably says.

## How to work

For each item, do **one independent computation**: a derivation written out, or a script of your own. You may run the book's scripts for comparison, but they are not evidence. A verdict without its own evidence counts as "unclear". `briefs/referee_defects_38_40_report.md` shows the expected depth.

Numerical hints:
- Use a generic $q$, such as $e^{i\sqrt2}$ or $e/2$. Short decimals such as 0.3 or 1.37 give a rational $\omega$, and then a root of unity, where representations degenerate and R-matrix code silently returns garbage.
- Before trusting a representation, check every defining relation, including the $q$-Serre relations.
- For a tensor product of dimension 64 or more, a dense solve of Jimbo's equations is unusable. Restrict the unknowns to the weight-preserving entries.

## Items (most important first)

1. **Normalizations of $U_q(\hat g)$, the central element, and the $S^2$ shift** (`sec-uqghat`, `sec-qg-crossing`).
   - Claim: $\prod_ik_i^{n_i}$ is central, and it equals $q^{\sum_in_i^\vee h_i}$, "$q$ to the power of the central element". Check this for an untwisted non-simply-laced algebra, a twisted one and $a_{2n}^{(2)}$. Which labels $n_i$ (book or Kac) and which normalization make each equality true?
   - Claim: $S^2$ is conjugation by $q^{2\hat\rho}$, $\hat\rho=\bar\rho+h^\vee\Lambda_0$. Its $\Lambda_0$ part multiplies $e_0$ by $q^{2h^\vee}$, so that $V(x)^{**}\cong V(x\,q^{2h^\vee})$ in the homogeneous gradation. With the Part VI coproduct the shift is $q^{-2h^\vee}$.
     - Derive the shift for a general affine algebra, as a function of $(\alpha_0,\alpha_0)$, the label of node 0 and $h^\vee$.
     - Does the stated shift hold for $a_{2n}^{(2)}$, where the text sets $(\alpha_0,\alpha_0)=1$, and with which node labelled 0?
     - Test numerically: build $\pi(S^2(a))$ and look for the $t$ with $V(x)^{**}\cong V(xt)$. Do it for the 3-dimensional representation of $U_q(a_2^{(2)})$, the spinor of $U_q(d_3^{(2)})$, and the vector of $U_q(c_2^{(1)})$ as an untwisted control. Use both coproducts.
   - Claim: if $V(x)^*\cong V(x\,x_c)$, then $x_c^2=q^{2h^\vee}$, and $S^2$ fixes $x_c$ only up to a sign.
   - Claim: in Part VI's gradation, $x(i\pi)=e^{2\pi iT}$. The power of $q$ is the $S^2$ shift, and $h^\vee$ of the Toda algebra in $2T=4\pi h/\beta^2-h^\vee$ supplies only the sign.
     - Check this against every crossing point that Part VI quotes (`grep -n 'x_c' part6-soliton-smatrices/*.qmd`).
     - Keep track of $\hat g$ against $\hat g^\vee$: the quantum group of the solitons of $\hat g$ is $U_q(\hat g^\vee)$.
   - Check the sentence "For twisted $\hat g=X^{(k)}$ … $q_i=q^k$ on [the longest roots], while $q_i=q$ on the roots of squared length 2", and its twin in `sec-qg-crossing` ("except for $a_{2n}^{(2)}$, whose roots have squared lengths 1, 2 and 4 …"). Are they consistent with Part VI's $U_q(d_{n+1}^{(2)})$ ($q_i=q^2$ on long roots, $q$ on short) and with ch. 20's $q_j=q^{\vec\alpha_j^{\vee2}/2}$?

2. **The sign of $q$ and the dictionary with Part VI** (`sec-six-vertex`, last two paragraphs).
   - Sine-Gordon: with $z=e^{2\lambda\theta}$ and $q=-e^{i\pi\lambda}$, check the following:
     - $b=\sinh\lambda\theta/\sinh\lambda(i\pi-\theta)$ is $S_T/S_0$ of `eq-sg-smatrix` (ch. 3);
     - the off-diagonal entries are $e^{\mp\lambda\theta}\sinh(i\pi\lambda)/\sinh\lambda(i\pi-\theta)$;
     - a diagonal, rapidity-dependent change of basis gives $S_R/S_0$, and it is the change to the principal gradation.
   - Claim: $P\check R(1/z)P$ intertwines the opposite coproduct and equals `eq-six-vertex` with $q\to q^{-1}$, at the same $z$, gradation and normalization. Hence $q=-e^{-i\pi\lambda}=e^{-8\pi^2i/\beta_{\rm SG}^2}$, Bernard and LeClair's value. Quote Bernard–LeClair's value of $q$ and their coproduct.
   - Claim: $q\to-q$ leaves $c,c'$ and $x_c=q^2$ unchanged, but reverses $b$, the combination $(1+b^2-cc')/2b=(q+q^{-1})/2$, and $c_u$. It equals conjugation of $\check R$ by $\mathrm{diag}(1,1,-1,1)$, which is not of the form $A\otimes A$.
     - Prove the "not of the form $A\otimes A$" claim. Is there a one-soliton basis change $A$, even $\theta$-dependent, that maps $\check R(q)$ to $\check R(-q)$?
     - Is the stated combination really basis-independent, under $A\otimes A$ changes of basis?
   - Check that the sentence "it changes the relative exchange sign of soliton and antisoliton, and it is what the crossing test of `sec-crossing-sign` detects" agrees with what `sec-crossing-sign` actually says.

3. **The twisted tensor-product graph and the spinor of $U_q(d_{n+1}^{(2)})$** (`sec-tpg`, paragraph "Twisted algebras", and `fig-tpg`).
   - Claim: under $U_q(b_n)$ the spinor tensor square is $\bigoplus_{k=0}^n\Lambda^k\mathbb C^{2n+1}$, multiplicity-free. The graph is the chain $\Lambda^n-\dots-\Lambda^0$, and $\rho_{\Lambda^{n-j}}=\prod_{i=1}^j\langle i/2\rangle_{(-1)^i}$, with $\langle a\rangle_+=(1+zq^{2a})/(z+q^{2a})$ on equal-parity edges.
   - Construct the spinor representation of $U_q(d_{n+1}^{(2)})$ yourself and check all relations. Solve Jimbo's equations directly for $n=2$ (dimension 16) and $n=3$ (dimension 64; weight-preserving unknowns), and compare the eigenvalues with the formula.
   - **In which $q$ does the formula hold?** Kac's $q$, with $q_i=q$ on the short roots, or the $q$ of the book's normalization, with long roots of squared length 2? Does the text make this clear? If not, give the sentence that would.
   - Check the parities of the components. Check the caption of `fig-tpg`, part (b), for $n=3$: $\mathbf8\otimes\mathbf8=(\Lambda^3\oplus\Lambda^0)_{\rm sym}\oplus(\Lambda^2\oplus\Lambda^1)_{\rm anti}$, with $\langle1\rangle_+$ on the middle edge.
   - Check that $e_0$ transforms in $g_{(1)}$ under $U_q(g_{(0)})$, as stated. Find the rule and the formula in DGZ 1996 and quote them; the citation names eq. (4.31) in ch. 22.

4. **The untwisted rule, and where the graph fails** (`sec-spectral-decomposition`, `sec-tpg`).
   - Check `eq-tpg-rule` (with $\langle a\rangle\langle-a\rangle=1$), the $U_q(\widehat{sl}_{n+1})$ vector example, and part (a) of `fig-tpg`. Do one case with unequal root lengths: the vector of $c_2^{(1)}$, where half-integer arguments occur.
   - Claim, "When it fails": the pairs that $e_0$ actually links form a connected subgraph of the extended graph whenever $V(x)\otimes V(y)$ is irreducible for generic $x/y$, citing DGZ 1994. Quote the source.
   - Claim: the other soliton multiplets of $c_n^{(1)}$ are Kirillov–Reshetikhin modules of $U_q(d_{n+1}^{(2)})$, reducible under $U_q(b_n)$. Examples: $\Lambda^1\oplus\Lambda^0$, and the 29-plet $\Lambda^2\oplus\Lambda^1\oplus\Lambda^0$ for $c_3^{(1)}$. Check the decompositions, and that their tensor squares have multiplicities.

5. **History and prior art** (chapter 8 introduction and notes, `sec-direct-method`).
   - Claim: Delius, Gould and Zhang derived the spinor R-matrices of $U_q(b_n^{(1)})$ and $U_q(d_{n+1}^{(2)})$ for all $n$, and the $\mathbf{26}$ R-matrix of $U_q(f_4^{(1)})$ (`dgz1994`, `dgz1996`). Locate each in the papers and quote it.
   - Claim: the $\mathbf{27}$ of $U_q(e_6^{(2)})$ is reducible under $U_q(f_4)$, but irreducible under the other horizontal subalgebra $U_q(c_4)$. Gandenberger, MacKay and Watts obtained its R-matrix by the tensor-product graph for $U_q(c_4)$ (ch. 22 cites their Sec. 4.4). The 2026-10-10 pass left this "unclear".
     - Check the branching under both subalgebras, and that the $\mathbf{27}\otimes\mathbf{27}$ is multiplicity-free under $U_q(c_4)$.
     - Read GMW and quote what they computed.
   - Check the attributions of the graph: MacKay for rational R-matrices, Zhang–Gould–Bracken for trigonometric ones, DGZ for pairs of representations and for twisted algebras. Check Jimbo 1986: "the R-matrices in the vector representations of all the classical quantum affine algebras, with the affine Toda equations in mind".
   - Appendix E (`appendices/E-original-results.qmd`) has no Part II item, but its Part VI items rest on these R-matrices. Check that nothing in chapters 7–8 or in Appendix E presents a known R-matrix as new.
   - `sec-direct-method`:
     - Check the count of unknowns, $\sum_\mu{\rm mult}_\mu^2$.
     - Check the $\mathbf{27}\otimes\mathbf{27}$ example: $26\otimes26=1\oplus26\oplus52\oplus273\oplus324$, multiplicities 2 and 3 for $\mathbf1$ and $\mathbf{26}$, and 16 unknowns.
     - Claim: imposing the $e_0$ equation on highest-weight vectors, projected onto highest-weight vectors, determines $\check R$. Is that enough in general, or only when it is supplemented by the $U_q(g)$-invariance?

6. **Roots of unity and RSOS restriction** (`sec-roots-of-unity`; `sec-irf`, paragraph "RSOS restriction"). This paragraph was rewritten by the reading pass and has had no referee.
   - For $q^{2p}=1$, check the following claims:
     - the summands of a tensor product of two representations with $2j+1<p$ are irreducible for $j<(p-1)/2$, and tilting modules $T(j)$ for $j\ge(p-1)/2$;
     - $T(j)$ is irreducible with $[p]_q=0$ for $2j+1=p$;
     - for $p<2j+1<2p$, $T(j)$ has dimension $2p$, is glued from spin $j$ and spin $p-1-j$, and has quantum dimension $[2j+1]_q+[2p-2j-1]_q=0$.
   - Check the example $p=4$: the tensor square of spin 1 is $T(2)$ plus spin 0, and $[5]_q=-1$ but $[5]_q+[3]_q=0$. Construct the module explicitly at a root of unity.
   - Which primitive root is meant? Does the sign of $q$ matter here?
   - Claim: summands of quantum dimension zero can be consistently discarded.
   - Claim: the face weights with heights $a\le(p-2)/2$ never couple to larger heights. Check it on the spin-$\frac12$ face weights.
   - Check the restriction to `prp-sg-rsos` (ch. 3) and `sec-rsos` (ch. 27).

7. **The core formulas of chapter 7.** They were checked on 2026-10-10 and by the book's scripts, but only by Claude. Re-derive them quickly and independently:
   - the Hopf axioms for `eq-coproduct` and `eq-antipode`, and $S^2(a)=KaK^{-1}$;
   - `eq-universal-r-sl2` and its printed spin-$\frac12$ matrix;
   - the three consequences listed in `sec-universal-r`;
   - `eq-evaluation-rep`, including the $q$-Serre relations for spin $j\ge1$;
   - `eq-six-vertex`, `eq-r-unitarity`, `eq-ybe-braid` (operator order as printed), and `eq-six-vertex-eigenvalue`;
   - the crossing example: $x_c=q^2$, $C=\left(\begin{smallmatrix}0&1\\-q&0\end{smallmatrix}\right)$, $c_u(z)=q(z-1)/(q^2z-1)$ in `eq-r-crossing`;
   - the $*$-structure, with $(*\otimes*)\circ\Delta=\Delta^{\rm op}\circ*$, and the spin-1 norm $[2]_q$.

   One script covering them all is enough.

8. **Fusion, poles and the IRF form** (`sec-fusion`, `sec-irf`, `sec-qg-hermitian`).
   - Check the operator order of `eq-fusion` and the argument that it is the R-matrix. Is the projector $P_c$ well defined where $V_a\otimes V_b$ is reducible but indecomposable?
   - Claim: the six-vertex $\rho_0$ has a pole at $z=q^2$, with residue proportional to the spin-0 projector, and a zero at $z=q^{-2}$, where $\check R$ projects onto spin 1.
   - Check the general claim that poles and zeros lie at $z=\varepsilon q^l$, with $\varepsilon$ a root of unity.
   - Claim: the vertex and IRF forms "are related by a change of basis, so they have the same eigenvalues and the same Bethe–Yang spectra". Is this accurate as worded, given that the IRF space is a space of highest-weight vectors (multiplicity spaces)?
   - Claim: an S-matrix that is not Hermitian analytic in the vertex basis can be Hermitian analytic in the IRF basis. Quote `hoare2013` for it.
   - Check the statements of `sec-qg-hermitian` and their citation `miramontes1999`.

9. **Citations and exercises.**
   - Check that every citing sentence says what its source says. Examples:
     - `delius1995`: the bootstrap axioms follow from the Hopf structure, and the gradation is fixed by Lorentz spin;
     - `frenkel1992` and `chari1994`: the $V^{**}$ shift;
     - `kulish1981`: fusion;
     - `gomez1996`: vertex–IRF and RSOS;
     - `drinfeld1987` and `jimbo1985`.
   - Recall the trap: `dgz1992` is Delius–Grisaru–Zanon, `dgz1996` is Delius–Gould–Zhang.
   - Work every exercise in outline: `exr-coproduct-homomorphism`, `exr-antipode`, `exr-spin-j`, `exr-six-vertex`, `exr-sg-match`, `exr-hermitian-spin1`, `exr-six-vertex-derive`, `exr-tpg-check`, `exr-tpg-an`, `exr-spin1`, `exr-faces`. Is each one true, and solvable with what the chapter provides?
   - Check the hint of `exr-spin1` in particular: $V_1(x)$ is the image of $\check R(q^{-2})$ acting on $V_{1/2}(xq^{-1})\otimes V_{1/2}(xq)$, and is the $q$-symmetric part of $V_{1/2}(xq)\otimes V_{1/2}(xq^{-1})$.
   - List any claim in these chapters that goes beyond standard textbook material and therefore needs a status label or a check script.

## Report format (`briefs/referee_algebra_07_08_report.md`)

A header with the date, your model, and the commit you reviewed (`git rev-parse HEAD`). Then one section per item, with:
- the verdict (correct / wrong / unclear), per claim where an item has several;
- your own computation (the script inline or its path, and its output);
- for errors, replacement text ready to paste, with file and line.

Then a list of the claims you could not check, and why. End with other observations, labelled as unverified.
