# Brief: fixes and gaps from the Part II reading (Chapters 5–8)

Status: done 2026-10-10 (uncommitted, branch `reading-fixes`; see Resolution at the end). Originally written 2026-10-10 after a line-by-line reading of Part II by a reader with graduate-level mathematical physics and QFT, who has met $su(2)$ and $su(3)$ but has no prior knowledge of Kac–Moody algebras, quantum groups or R-matrices.

Files: `part2-algebra/part.qmd`, `05-root-systems.qmd`, `06-affine-algebras.qmd`, `07-quantum-groups.qmd`, `08-r-matrices.qmd`.

Scope: fix the items below in place. Items in A1–A3 need a check against the literature before any change: state in your report what you checked and what you could not confirm. Do not edit anything else without checking with the author. The Part II results that I checked by hand are listed in section D.

## A. Errors and likely errors

**A1. Folding table, `tbl-foldings` in `06-affine-algebras.qmd` (row $a_{2n}^{(2)}$).** The row gives the parent as $d_{2n+2}^{(1)}$. Two checks suggest this is wrong.
- Node count. $a_2^{(2)}$ ($n=1$) has 2 nodes. $d_4^{(1)}$ has 5 nodes. An automorphism of order 2 has orbits of size at most 2, so folding $d_4^{(1)}$ by one gives at least 3 orbits. (A larger group, such as the $\mathbb Z_2\times\mathbb Z_2$ or $\mathbb Z_4$ actions on the four end nodes, would give 2 orbits, but it is not an order-2 twist and so does not fit the framework of "Twisted affine algebras", `sec-twisted`.)
- Coxeter number. The folded $h=2n+1$ equals the Coxeter number of $a_{2n}^{(1)}$, which is also $2n+1$, so the natural parent is the $a_{2n}^{(1)}$ with a reflection fixing node 0, which has $n+1$ orbits. This is the case with no common factor $k$, which is consistent with the self-duality of $a_{2n}^{(2)}$ in the table. The $d_{2n+2}^{(1)}$ entry has the right $h$ only after dividing by $k=2$, which the other rows do not require.
- Fix: check Olive–Turok (1983) and Kac (1990, Table Aff 2) for the parent of $a_{2n}^{(2)}$. If confirmed, change the parent to $a_{2n}^{(1)}$ and the $k$-column accordingly, and then recheck the table's $h$ and $h^\vee$ columns.

**A2. Spins and Cartan eigenvalues, `sec-perron-frobenius` in `05-root-systems.qmd`.** The text says: "for every spin $s$ of a conserved charge, the eigenvalues $\chi_a^{(s)}$ form an eigenvector of the Cartan matrix with eigenvalue $4\sin^2\frac{\pi s}{2h}$". Spins of conserved charges are $s=m_k+kh$ (`sec-affine-algebras`, `sec-sine-gordon`), and $4\sin^2\frac{\pi s}{2h}$ is not invariant under $s\mapsto s+h$, since $\sin^2$ becomes $\cos^2$. So the claim holds only for the basic period $s=m_k$, or needs a different formula for shifted spins. Fix: state the claim for $s=m_k\in\{1,\dots,h-1\}$, or give the correct eigenvalue for $s=m_k+kh$ (from Bernard–LeClair or Kac). Check which before editing.

**A3. Bicoloured Coxeter element, `sec-weyl-coxeter` in `05-root-systems.qmd`.** "the projections of the roots onto this plane form $r$ regular $h$-gons" holds for simply-laced algebras. For non-simply-laced algebras the orbits have different sizes. Add "for simply-laced $g$".

**A4. Strong–weak duality, `sec-affine-folding` (last bullet of "Duality") in `06-affine-algebras.qmd`.** The sentence says "The real-coupling S-matrices of a dual pair are exchanged by the strong–weak duality $\beta\to4\pi/\beta$" and then "in the normalization of this book it reads $\beta^{\vee2}=16\pi^2k/\beta^2$". These two statements are not reconciled. The first is for a single algebra in the untwisted case, which agrees with `sec-sg-breather-particle` in Part I ($\beta_{\rm SG}\to8\pi/\beta_{\rm SG}$). The $k$-dependence in the second is not derived in the text. Fix: state the untwisted and twisted cases separately, give the normalization of $\beta^\vee$ used for each, and derive the $k$-factor (or cite where it is derived).

**A5. Quantum dimension of summands, `sec-roots-of-unity` in `07-quantum-groups.qmd`.** "the summands with highest weight $j\ge(p-1)/2$ ... have quantum dimension zero." For $2j+1>p$, $[2j+1]_q$ is not zero in general (for example $p=3$, $2j+1=4$, $q^6=1$, $[4]_q=-[2]_q\neq0$). What vanishes is the quantum dimension of the indecomposable tilting modules in that range, not of each simple summand. Rewrite the sentence in terms of the indecomposable tilting modules (or the negligible summands) and check the statement for one small case.

**A6. Spectral parameter in the hint of `exr-spin1` in `08-r-matrices.qmd`.** The hint says $V_1(x)$ is "the $q$-symmetric part of $V_{1/2}(xq)\otimes V_{1/2}(xq^{-1})$, the image of $\check R(q^{-2})$". With $z=x/y$ as in `eq-r-check`, the first factor $V(xq)$ and the second $V(xq^{-1})$ give $z=q^2$, not $q^{-2}$. Also, `sec-six-vertex` shows that $\check R(q^{-2})$ is the spin-1 projector only if the ordering convention is the one used there. Check the ordering and correct one of the two places.

**A7. Level-zero statement, `sec-loop-algebras` in `06-affine-algebras.qmd`.** "the finite-dimensional representations used for solitons and R-matrices have level zero." The solitons use the infinite-dimensional level-one representations $L(\Lambda_j)$ of `sec-vertex-operators`. Only the R-matrices use finite-dimensional, level-zero representations. Rewrite the sentence to say so.

**A8. Twisted normalization, section `sec-qg-crossing` in `07-quantum-groups.qmd` (the sentence "For untwisted algebras this is the normalization of `sec-lie-algebraic-data`; for twisted ones the short roots then have squared length 2").** This conflicts with `sec-twisted`/`tbl-foldings` in `06-affine-algebras.qmd`, where the longest root of $X_n^{(k)}$ has squared length $2k$, and with the conversion $\beta_{\rm lit}^2=\beta^2/k$ in Part I, `01-introduction.qmd`. Check which root has $(\alpha_0,\alpha_0)=2$ for each twisted algebra. Make the sentence say which root is normalized, and check that $q_i=q^{\alpha_i^2/2}$ gives the intended $q_i$.

**A9. Figure of $\mathbf{27}$, `sec-uqghat` in `07-quantum-groups.qmd`.** "the $\mathbf{27}=\mathbf{26}\oplus\mathbf1$ of $U_q(e_6^{(2)})$". The $\mathbf{27}$ is a representation of $U_q(e_6)$, and $26\oplus1$ is its decomposition under the fixed-point subalgebra $f_4$ (the classical part of $e_6^{(2)}$). Rephrase so that the algebra and the decomposition are not conflated.

## B. Gaps for a reader without Part II

**B1. Terms used without definition.** Each needs a one-line gloss or a pointer.
- Chapter 5: "Killing form (suitably normalized)"; "Chevalley–Serre presentation" (defined by its relations, which is fine); "Dorey's fusing rule" (forward pointer, `sec-simply-laced`); "Steinberg's bicolouring" (described, not proved).
- Chapter 6: "principal Heisenberg subalgebra" (not defined; only said to be the centralizer of $E_1$); "Lax pair of `sec-theories`" (forward pointer); "minimal affinization" and "Kirillov–Reshetikhin module" (`sec-quantum-groups`, "Evaluation representations", defined only by name); "Lie duality" (forward pointer, `sec-imaginary-coupling-vacua`).
- Chapter 7: "Hopf algebra", "quasitriangular", "cocommutative" (these are defined in the chapter, but the reader must follow the coproduct and antipode conventions carefully); "$q$-Clebsch–Gordan coefficients" (used in `exr-faces`, not defined); "crossing point $x_c$" (defined by $V(x)^*\cong\bar V(xx_c)$, which is a dense definition).
- Chapter 8: "isotypic component" (used in `sec-direct-method` and `sec-tpg`, not defined); "extended graph" (named, not defined); "face weights" and "paths $a\to b\to d$" (defined in `sec-irf`, but the rules for allowed paths are only implicit); "interaction-round-a-face" (named only).

Fix: add a short glossary box at the start of Chapter 7 for the Hopf-algebra terms, and one sentence each for the other items.

**B2. Forward references in chapter 5.** The Dorey fusing rule and the Lie duality are used as results before Parts III–VI. This is acceptable, but say so in the chapter introduction.

**B3. The "assumed background" statement.** The chapter introduction says "assuming only that the reader has met $su(2)$ and $su(3)$ in physics". Chapter 7 assumes the Hopf-algebra language, which is not in that list. Amend the statement to name the Hopf-algebra background, or to refer the reader to Kassel (already in the references).

**B4. $\mathbb Z_{n+1}$ charge.** `sec-affine-dynkin` in `06-affine-algebras.qmd` says the particles of $a_n^{(1)}$ carry a $\mathbb Z_{n+1}$ charge. The mechanism (the centre of $SU(n+1)$ acting as rotations of the affine diagram) is stated but not connected to the topological charges of `sec-imaginary-coupling-vacua`. Add the connecting sentence or a pointer.

## C. Notation clashes

**C1. $k$** is the level of a Kac–Moody algebra (`sec-loop-algebras`, `sec-affine-bracket`), the order of a twist (`sec-twisted`, `tbl-foldings`), and the $k$ in $X_n^{(k)}$, and the generators $k_i$ of $U_q$ (`sec-uqghat`). Suggestion: keep $k$ for the twist order, call the level $\ell$ (or write "level $k$" with a subscript) in `sec-loop-algebras`, and leave $k_i$.

**C2. $E$** is the root vector $E^{\vec\alpha}$ (chapter 5), the spin-½ generator $E$ of $U_q(sl_2)$ (chapter 7), and the Lax element $E_1$ (chapter 6). Suggestion: rename the Lax element $\mathcal E_1$ or $L_1$.

**C3. $\lambda$** is the loop variable of the loop algebra (`sec-loop-algebras`), the fundamental weights $\vec\lambda_i$ (chapters 5 and 8), and the sine-Gordon parameter $\lambda=8\pi/\beta^2-1$ (chapter 7, `sec-six-vertex`; Part I). In chapter 8 the weight $2\vec\lambda_1$ appears next to the loop parameter. Suggestion: keep $\vec\lambda_i$ for weights and write the sine-Gordon parameter as in Part I, or rename the loop variable.

**C4. $\rho$** is the Weyl vector (chapter 5, `eq-casimir`), the Weyl vector $\hat\rho$ of chapter 7, and the eigenvalues $\rho_\mu(z)$ of chapter 8 and of the six-vertex model (chapter 7). Suggestion: rename the eigenvalues to $r_\mu(z)$ or $\varrho_\mu(z)$.

**C5. $P$** is the flip of tensor factors (chapter 7), the projector $P_\mu$ (chapter 8), and field parity $P$ (Part I, `sec-antilinear`). Suggestion: write the flip as $\mathsf P$ or $\sigma$ in chapter 7, and keep $P_\mu$.

**C6. $d$** is the dimension of $\mathfrak g$ (chapter 5, the table `tbl-simple-lie`) and the derivation $d=\lambda\,d/d\lambda$ (chapter 6, `sec-loop-algebras`). Suggestion: rename the derivation as $\mathsf d$ or $\partial_\lambda$.

**C7. $n$** is the rank of $A_n$, $B_n$, etc., the Kac label $n_j$ and dual label $n_j^\vee$, and the summation index in the universal R-matrix (chapter 7, `eq-universal-r-sl2`). Suggestion: rename the summation index to $m$ or $j$ in `eq-universal-r-sl2`.

**C8. $m$** is an exponent $m_k$ (chapter 5), a loop mode $T^a\otimes\lambda^m$ (chapter 6), and a mass $m$ in the vertex-operator formula (chapter 6). Suggestion: keep $m_k$ for exponents and use $n$ or $p$ for loop modes after resolving C7.

**C9. $F$** is the generator $F$ of $U_q(sl_2)$ (chapter 7), the vertex operators $F^a(z)$ (chapter 6), and the fundamental-representation labels elsewhere. Suggestion: rename the vertex operators $\mathcal F^a(z)$.

**C10. $q$ is also $q_i=q^{\vec\alpha_i^2/2}$** (chapter 7), and $q_i$ is used with the exponent $\vec\alpha_i^2/2$. This is consistent, but the reader must be told when $q_i\neq q$ (short roots of twisted algebras). Add the sentence at first use.

## D. Items checked and found correct (do not re-check unless you change them)

- Chapter 5 table of simple Lie algebras: $d$, $h$, $h^\vee$, exponents for $A_n$–$G_2$; check $d=r(h+1)$ (exr-roots-count).
- Coxeter number as $1+\sum n_i$ and $h^\vee=1+\sum n_i^\vee$ (`eq-coxeter-numbers`); the exponents are symmetric under $m\mapsto h-m$.
- Cartan eigenvalues $4\sin^2(\pi m_k/2h)$ (`eq-cartan-eigenvalues`), and the $A_n$ eigenvector $\sin(a\pi/(n+1))$.
- $a_1^{(1)}$ mass $2m$ and $A_n^{(1)}$ masses $2m\sin(a\pi/h)$ (`sec-perron-frobenius`).
- $E_8$ mass ratio $2\cos(\pi/5)$ and the $D_4$ eigenvector with the central component $\sqrt3$ times the legs (exr-d4-masses).
- Chapter 6 untwisted table: $n_j$ sums give $h$ for every row; $e_8^{(1)}$ labels sum to 30, $e_7^{(1)}$ to 18.
- Chapter 6 twisted list and ranks ($a_{2n}^{(2)}$, $a_{2n-1}^{(2)}$, $d_{n+1}^{(2)}$ rank $n$; $e_6^{(2)}$ rank 4; $d_4^{(3)}$ rank 2).
- Foldings table: all rows except $a_{2n}^{(2)}$ (see A1) have $h$ and $h^\vee$ swapped correctly with their duals; the $a_3^{(1)}\to c_2^{(1)}$ example has labels $1,1,2$ and $\alpha_O^2=1$.
- $d_4^{(1)}\to g_2^{(1)}$ by triality: orbits of sizes $1,1,3$ give labels $1,2,3$ and $\alpha_O^2=2/3$ (short), consistent with `exr-fold-d4`.
- Principal gradation: the exponents $1,5 \bmod 6$ for $a_2^{(2)}$ match the Bullough–Dodd spins $1,5,7,11,\dots$.
- Chapter 7: coproduct is a homomorphism (`exr-coproduct-homomorphism`); antipode axiom and $S^2(a)=KaK^{-1}$; the spin-$\frac12$ R-matrix equals $q^{1/2}\times$ the stated matrix; $\check R(z)\check R(1/z)=1$ and eigenvalues $1$ and $\rho_0(z)$ for the six-vertex matrix; $\rho_0$ pole at $z=q^2$ and zero at $z=q^{-2}$.
- Sine-Gordon substitution: $b=\sinh\lambda\theta/\sinh\lambda(i\pi-\theta)$, $a=e^{\lambda\theta}\sinh(i\pi\lambda)/\sinh\lambda(i\pi-\theta)$ with $z=e^{2\lambda\theta}$, $q=-e^{i\pi\lambda}=e^{8\pi^2 i/\beta_{\rm SG}^2}$.
- $(1+b^2-cc')/(2b)=(q+q^{-1})/2$, so $q\to-q$ reverses it; conjugation by $\mathrm{diag}(1,1,-1,1)$ flips only $b$.
- Crossing: $x_c=q^2$ for $U_q(\widehat{sl}_2)$, and the $T$-formula $2T=4\pi h/\beta^2-h^\vee$ (`sec-amplitude-and-gradation`) gives $x(i\pi)=e^{2\pi i\lambda}=q^2$ for $a_1^{(1)}$, consistent with `sec-qg-crossing`.
- Hermitian form: $E^*=F$, $K^*=K^{-1}$ is compatible with the relations; $\langle Fv_0,Fv_0\rangle=[2]_q$ on spin 1; $(*\otimes*)\Delta=\Delta^{\rm op}*$.
- Chapter 8: the TPG rule reproduces $\rho_0/\rho_1=\langle1\rangle$ for spin $\tfrac12\otimes\tfrac12$ and $\rho_2=1$, $\rho_1=\langle2\rangle$, $\rho_0=\langle2\rangle\langle1\rangle$ for spin $1\otimes1$ (`exr-spin1`), using $C(j)=2j(j+1)$ (`exr-tpg-check`); $\langle a\rangle\langle-a\rangle=1$.
- $C(2\vec\lambda_1)-C(\vec\lambda_2)=4$ for $sl_{n+1}$ (checked for $n=2$: $20/3-8/3$).
- Spinor of $d_{n+1}^{(2)}$: $V\otimes V=\bigoplus_{k=0}^n\Lambda^k\mathbb C^{2n+1}$ has dimension $2^{2n}$, consistent with $\dim V=2^n$.
- $\mathbf{27}\otimes\mathbf{27}$ under $f_4$: multiplicities $\mathbf1\times2$, $\mathbf{26}\times3$, and $\mathbf{52},\mathbf{273},\mathbf{324}$ once each, giving $2^2+3^2+1+1+1=16$.
- Fusion and the bootstrap structure of `eq-fusion`.
- All 70 `@` labels and all 34 citation keys used in Part II resolve, and the citation check was tested against a fake key.

## E. Verification checklist for the agent who makes the fixes

1. Re-run the label check: every `@sec-`, `@eq-`, `@prp-`, `@exr-`, `@tbl-` and `@thm-` in `part2-algebra/` resolves to an anchor elsewhere in the book (exclude `.claude/` and `_site/`).
2. For A1, check the parent of $a_{2n}^{(2)}$ against Olive–Turok (1983) or Kac (1990) before changing the table.
3. For each notation rename (C1–C9), grep the whole book (excluding `.claude/` and `_site/`) for the old symbol in the same meaning, and change only those occurrences. Chapter 6's $\lambda$ and the sine-Gordon $\lambda$ (C3) need care.
4. Report what changed, what was left, and anything you could not verify.

## Resolution (2026-10-10)

Checks are in the scratch scripts `checks.py` (A2, A3, A5, A6) and `labels.py` (label check), run with numpy.

- **A1 (partly changed; the brief's claim was only partly right).** $d_{2n+2}^{(1)}$ is a valid parent, but through an order-4 automorphism that moves the affine node. It reverses the chain and permutes the four end nodes cyclically. The orbits give squared lengths $\frac12,1,\dots,1,2$, labels $4,\dots,4,2$ with common factor 2, and $h=2n+1$. The reflection of $a_{2n}^{(1)}$ is also valid, but its orbit $\{n,n+1\}$ consists of adjacent roots, which gives $\vec\alpha_O^2=\frac12$. Khastgir–Sasaki (hep-th/9512158, eq. (4.7) and Fig. 4, and case "iv) $d_{2n+2}\to a_{2n}^{(2)}$") list both foldings. Part V uses both (ch. 18 the reflection, ch. 19 and App. A the order-4 one). Changes: the table cell now reads "$a_{2n}^{(1)}$ or $d_{2n+2}^{(1)}$ (see text)", and a paragraph after the table works out both foldings. The general folding rule now flags $a_{2n}^{(2)}$ as the exception. Added the bib entry `khastgir1996b` (DOI 10.1143/PTP.95.503). The $h$ and $h^\vee$ columns were rechecked and are correct. Not accessed: Olive–Turok 1983 (paywalled) and Kac Table Aff 2 (book). As a twisted loop algebra, $a_{2n}^{(2)}$ comes from $A_{2n}$, as `sec-twisted` already says.
- **A2 (claim wrong; clarified).** $4\sin^2\frac{\pi s}{2h}$ is correct for every spin. The spin $m_k+h$ gives the eigenvalue of the exponent $h-m_k$, and its eigenvector is that of $m_k$ with the bicolouring sign applied. For $a_n^{(1)}$, $\chi_a^{(s)}\propto\sin(as\pi/h)$ is an eigenvector for all $s$, which was checked numerically. Added two sentences saying this.
- **A3 (claim wrong; clarified).** Every Coxeter-element orbit has exactly $h$ elements for every irreducible $g$, and the projections form regular $h$-gons. In the non-simply-laced case they have two radii. Checked for $B_2$, $B_3$ and $G_2$. Added "of different radii", "quoted without proof" and a one-line gloss of Dorey's rule.
- **A4 (changed).** Rewrote the bullet. $(\hat g,\vec\alpha,\beta)\leftrightarrow(\hat g^\vee,\vec\alpha^\vee,4\pi/\beta)$ holds in any normalization. The simply-laced case is self-dual, which matches $\beta_{\rm SG}\to8\pi/\beta_{\rm SG}$. For the dual pairs the coroots of $\hat g$ give $\hat g^\vee$ with longest roots of squared length $2k$, and rescaling yields $\beta^{\vee2}=16\pi^2k/\beta^2$. The coroot form is the standard statement, but I did not confirm it against DGZ or CDS text.
- **A5 (changed).** Rewritten in terms of tilting modules $T(j)$, whose dimension is $2p$ for $p<2j+1<2p$ and whose quantum dimension is $[2j+1]_q+[2p-2j-1]_q=0$. Added the example $p=4$: $[5]_q=-1$ but $[5]_q+[3]_q=0$, as for $T(2)$ in spin $1\otimes1$. Checked numerically for $p=3,4,5,7$.
- **A6 (claim wrong; clarified).** The image of $\check R(q^{-2})$, acting on $V(xq^{-1})\otimes V(xq)$, lies in $V(xq)\otimes V(xq^{-1})$ and is the $q$-symmetric part. Checked numerically: $\Delta(e_0)$ preserves the $q$-symmetric subspace exactly when $x_1/x_2=q^2$, and $\check R(q^{-2})$ has rank 3 with image equal to that subspace. Added "acting on $V_{1/2}(xq^{-1})\otimes V_{1/2}(xq)$" to the hint.
- **A7 (changed).** Level zero now applies to the R-matrix representations only. The solitons live in the highest-weight $L(\Lambda_j)$, which have nonzero level.
- **A8 (changed; there was no real conflict).** In Kac's normalization, $(\alpha_0,\alpha_0)=2$ means short roots of squared length 2 and longest roots of squared length $2k$, which is the same as Part I's literature normalization. The exception is $a_{2n}^{(2)}$, where $(\alpha_0,\alpha_0)=1$ and the squared lengths are 1, 2 and 4. The sentence now says this, and a matching sentence went into `sec-uqghat` (see C10). Part VI's $U_q(d_{n+1}^{(2)})$ ($q_i=q^2$ on long roots, $q$ on short) is consistent with it.
- **A9 (changed).** Now "a 27-dimensional representation of $U_q(e_6^{(2)})$ that decomposes as $\mathbf{26}\oplus\mathbf1$ under $U_q(f_4)$; classically the $\mathbf{27}$ of $e_6$ restricted to $f_4$". The 27 is a representation of the affine algebra, so the original phrase was not wrong, only compressed.
- **B1 (changed).** Glosses added for: the Killing form, Dorey's fusing rule and Steinberg's argument (ch. 5); the principal Heisenberg subalgebra (with where the name comes from) and the Lax pair (ch. 6); a Hopf-algebra glossary box at the start of ch. 7 and a plain-words gloss of $x_c$ (ch. 7); isotypic component, extended graph, the rule for allowed paths, IRF, and $q$-Clebsch–Gordan coefficients (in `exr-faces`) (ch. 8). Lie duality is already defined inline in `sec-affine-folding`, so nothing was added for it. Minimal affinization and KR modules were left as defined.
- **B2, B3 (changed).** The ch. 5 introduction now names the forward-quoted results and the Hopf-algebra background, with a pointer to Kassel.
- **B4 (changed).** Added a sentence: the $\mathbb Z_{n+1}$ charge of a species-$a$ soliton is the class $a\bmod(n+1)$ of its topological charge, a weight of $\Lambda^a\mathbb C^{n+1}$ (`sec-topological-charges`).
- **C1:** the level was renamed $\ell$ (2 occurrences, ch. 6 only). **C2:** a note was added, not a rename ($E_1$ is also used in Part III). **C3:** a note was added at the loop algebra. The projector subscripts in ch. 8 now carry arrows, $P_{2\vec\lambda_1}$ and $P_{\vec\lambda_2}$. **C4:** a note was added in ch. 8 (the $\rho_\mu$ also appear in Parts VI and VII). **C5:** a note was added in ch. 7 (the flip $P$ is used in Parts VI and VII). **C6:** the derivation was renamed $\mathsf d$ (2 occurrences, ch. 6 and ch. 7). **C7:** the summation index of `eq-universal-r-sl2` was renamed $m$ (2 occurrences). **C8:** no change. The loop modes are local to two displays, and the mass $m$ is book-wide. **C9:** a note was added, not a rename ($F^a$ is also used in ch. 10 and ch. 33). **C10:** added at first use in `sec-uqghat`.
- **Verification.** All `@` labels and citation keys in ch. 5–8 resolve, and the check was tested against fake ids. The scratch `quarto render --to html` finished with no errors, and there is no `?@` anywhere in `_site`.
