# Referee report: chapters 5 and 6

Date: 2026-10-15. Model: GPT-6 (Codex).

Reviewed commit: `5a53d5c7e7eeb9836dedd2be1683e0d7ae987bcf`, in the worktree named in the brief, `.claude/worktrees/book-preface-section-reorganize-7f6a81`. The main checkout was at `2d36c9c2c1848067a87a3c342f2549c43a25fe96`; both chapter files were byte-identical between the two checkouts. File and line references below refer to the reviewed worktree. No book files, existing reports, index entries or SVGs were changed. No author check script or Quarto render was run. No subagents were used.

The main mathematical tables survive this pass. In particular, the apparent conflict between chapters 6 and 14 about the dual coupling is **not an error**. The principal-degree examples and the higher-spin colour reversal also survive independent computations. Corrections are needed to the MacKay–Watts normalization attribution, the unqualified formula for the gradation parameter, the scope of the Kneipp–Olive citation, the assertion about different Coxeter-plane radii, and the assertion that all diagram symmetries act by phases. The comparison of the two parents of $a_{2n}^{(2)}$ needs an explicit identification of their mass parameters.

## Evidence and reproducibility

All independent scripts and their captured output are in `referee_scratch/algebra_05_06/` in this worktree:

- `check_algebra.py` and `check_algebra.out`: finite roots constructed by reflection closure from independently specified Euclidean roots/Gram matrices; affine marks; exhaustive diagram automorphisms; orbit-average folding; primitive right and left null vectors; restriction of parent Hessians; Coxeter spectra; conservation equations at all root-triangle fusings; symbolic dual-coupling checks.
- `check_heisenberg.py` and `check_heisenberg.out`: explicit twisted-loop matrix centralizers in $\mathfrak{sl}_3,\mathfrak{sl}_5,\mathfrak{so}_8$.
- `check_misc.py` and `check_misc.out`: rank-two duality compatibility, particle symmetry action, weight classes, mass Hessians, principal $\mathfrak{sl}_2$ decompositions, the minimal $E_8$ $S_{11}$ block identity, and elementary exercises.

Run with `python3 -u referee_scratch/algebra_05_06/check_algebra.py`, and similarly for the other two scripts. Dependencies are NumPy and SymPy. None imports book code. Assertions check the actual equations; the output also records intermediate data. Finite samples are identified below and are not presented as proofs for arbitrary rank.

Local primary sources were read from `/home/gustav/referee-copy/papers/`. Formula pages of MacKay–Watts, Olive–Turok–Underwood and Khastgir–Sasaki were also inspected as rendered PDFs. The committed chapter SVGs were rasterized into scratch space and inspected without running their generators. The earlier Claude report and reading-pass notes were read for context, not used as evidence.

## 1. Duality and coupling normalization

**Verdicts:** coroot rescaling **correct**; compatibility with chapter 14 **correct**; dual-pair list **correct**; Lie duality **correct for classical masses**, but qualify that scope; dual quantum-group identification **correct at the root-data level**. The exact quantum S-matrix interpretation retains the conjectural status stated in chapter 14. DGZ and CDS should not jointly be cited as having established an all-orders exceptional-family coupling map.

### Independent root and coupling calculation

For $C_2^{(1)}$, choose

\[
\alpha_0=-\sqrt2 e_1,\quad
\alpha_1=(e_1-e_2)/\sqrt2,\quad
\alpha_2=\sqrt2 e_2.
\]

The squared lengths are $(2,1,2)$, the marks are $(1,2,1)$, and the coroot squared lengths are $(2,4,2)$. Dividing the coroots by $\sqrt2$ gives the book-normalized $D_3^{(2)}$ roots, with marks $(1,1,1)$.

For $B_2^{(1)}$, take $(-e_1-e_2,e_1-e_2,e_2)$. Lengths are $(2,2,1)$, marks $(1,1,2)$, and coroot lengths $(2,2,4)$. Dividing the coroots by $\sqrt2$ gives $A_3^{(2)}$, with marks $(1,1,1)$.

For $G_2$, take the short simple root $a=(\sqrt{2/3},0)$, the long one $b=(-\sqrt{3/2},1/\sqrt2)$, and $\alpha_0=-3a-2b$. In the order $(\alpha_0,b,a)$, the lengths are $(2,2,2/3)$ and marks $(1,2,3)$. Coroot lengths are $(2,2,6)$; division by $\sqrt3$ gives $D_4^{(3)}$ with marks $(1,2,1)$.

In each case the coroot Cartan matrix is the transpose, since

\[
\frac{2\alpha_i^\vee\cdot\alpha_j^\vee}{|\alpha_i^\vee|^2}
=\frac{2\alpha_i\cdot\alpha_j}{|\alpha_j|^2}=a_{ji}.
\]

If the coroot coupling is $4\pi/\beta$, replacing the coroots by $\alpha^\vee/\sqrt k$ gives the same exponential arguments at coupling $4\pi\sqrt k/\beta$. Rescaling the original roots by $c$ and $\beta$ by $c^{-1}$ leaves this relation covariant. The overall mass parameter can also be rescaled; it does not enter the dimensionless two-particle S-matrix. The simply-laced and sine-Gordon conversions in chapter 6 are correct. The four-pair formula is not intended to use $k=2$ as the root-length ratio for the self-dual $A_{2n}^{(2)}$: that system has squared-length ratio 4 and gives the book-normalized map $\beta' =8\pi/\beta$.

Write $b=\beta^2$, and $F(u)=2u/(u+4\pi)$. Direct substitution into the book's formulas gives:

| pair | $B_u$ | $B_t$ | common $H$ |
|---|---|---|---|
| $C_2^{(1)},D_3^{(2)}$ | $F(b/2)$ | $F(b^\vee)$ | $4+B_u=6-B_t$ |
| $B_2^{(1)},A_3^{(2)}$ | $F(b)$ | $F(b^\vee/2)$ | $4-B_u/2=3+B_t/2$ |
| $G_2^{(1)},D_4^{(3)}$ | $F(b/3)$ | $F(b^\vee)$ | $6+3B_u=12-3B_t$ |

Here $b^\vee=16\pi^2k/b$. Every row has $B_t=2-B_u$, exactly, not just at one loop. The script simplifies both residuals to zero. At $b=\sqrt2$ it gives respectively

```text
C2/D3t: B=0.1065443154 Bdual=1.8934556846; B residual=0, H residual=0
B2/A3t: B=0.2023110829 Bdual=1.7976889171; B residual=0, H residual=0
G2/D4t: B=0.0723136449 Bdual=1.9276863551; B residual=0, H residual=0
```

The $B_2=C_2$ coincidence is a further useful check, since their two mass formulas use different $H$'s. Symbolically $1/H_B+1/H_C=1/2$, so $2\sin(\pi/H_B)=2\cos(\pi/H_C)$; their mass ratios agree at the same book coupling. The exceptional $F_4/E_6^{(2)}$ pair has the same $B$ algebra as $C_n/D_{n+1}^{(2)}$, with $12+3B_u=18-3B_t$.

### Lie duality and the non-local algebra

An independent $C_3$ check uses the $A_5$ parent. Its masses in units of $m$ are $(1,\sqrt3,2,\sqrt3,1)$. Folding pairs species $1,5$ and $2,4$, leaving species 3 fixed. Additivity of the classical multi-soliton energy therefore gives folded soliton ratios $(2,2\sqrt3,2)$. The $B_3$ Hessian has masses $(\sqrt2,\sqrt6,\sqrt2)$, in the corresponding ordering. These ratios agree. The $D_4^{(2)}$ ratios, proportional to $(\sin(\pi/8),\sin(2\pi/8),\sin(3\pi/8))$, do not. Thus the distinction made in line 120 is essential and correct.

The coroot calculation above also independently identifies the Cartan matrix of the non-local-charge algebra as $A^T$. It checks the claimed dual algebra and the consistency with chapter 20; it does not re-prove all quantum Serre relations or the existence of all charges. Those need the Bernard–LeClair/Felder sources cited there, not DGZ/CDS.

### Source audit and replacement

DGZ, `hep-th/9201067`, §9, equations (9.6)–(9.7), explicitly exchange the $A_{2n-1}^{(2)}$ and $B_n^{(1)}$ amplitudes and say the analogous $C_n^{(1)}/D_{n+1}^{(2)}$ exchange works. The exact phrase is “under the substitution $\beta\to4\pi/\beta$”; this is in their conventions. CDS, `hep-th/9304065`, §1, define duality by the exact replacement “$\alpha_i\to\check\alpha_i=2\alpha_i/\alpha_i^2$”. Their §4 also says: “it is not intended in this article to relate B explicitly to the coupling constant.” These are different statements; the last one prevents using CDS as a derivation of the all-orders $B(\beta)$ conjecture.

MacKay–Watts, `hep-th/9411169`, §6.3.1, give the normalized relation directly: “$\beta^2\beta'^2=16\pi^2kk'$”. Their §2 distinguishes classical Lie duality from affine duality. Thus a much better direct citation for the normalized coupling relation is MacKay–Watts, §6.3.1. No explicit claim about non-local quantum-group generators was located in DGZ/CDS; do not assign that claim to them.

At `part2-algebra/06-affine-algebras.qmd:120`, change the start to:

> The classical soliton masses of an untwisted non-simply-laced theory $g^{(1)}$ are proportional to the classical particle masses of $(g^\vee)^{(1)}$, the untwisted algebra of the dual root system ("Lie duality", @sec-imaginary-coupling-vacua).

At line 122, replace its first sentence by:

> The proposed real-coupling S-matrices obey a strong–weak duality (@sec-strong-weak, @sec-nsl-duality); in longest-root normalization its coupling map is stated in [@mackay_watts1995, sec. 6.3.1].

Retain the subsequent coroot calculation. The list at lines 179–182 is correct as a classification of the root systems used by the Toda potentials.

## 2. The two parents of $a_{2n}^{(2)}$

**Verdicts:** the two specified automorphisms, orbits, lengths, raw labels, common factors and $h$ are **correct**. The sentence “the same classical masses” is **wrong at equal parent $m$** and **correct after identifying the folded mass scales**. The normalization must be stated explicitly. Part V uses the correct node actions but inherits this ambiguity when comparing the two sets of masses.

### Independent folds and Cartan matrices

For the odd cycle, $j\mapsto-j$ preserves its edges and fixes only 0. All two-element orbits consist of orthogonal roots except $\{n,n+1\}$. Thus

\[
\left|\frac{\alpha_n+\alpha_{n+1}}2\right|^2
=\frac{2+2-2}{4}=\frac12.
\]

The labels are $1,2,\ldots,2$; their gcd is 1 and their sum is $2n+1$.

For $D_{2n+2}^{(1)}$, the four-cycle in the brief and the chain reversal preserve every adjacency. Its four end roots are mutually orthogonal, the chain-pair averages have length squared 1, and the middle root has length squared 2. This includes $n=1$: there are no chain pairs, and the central root of $D_4^{(1)}$ is fixed. Raw labels, in order from the four-end orbit to the middle, are $4,4,\ldots,4,2$, with gcd 2. The script checks the actual permutations for $n=1,2,3$, not merely their orbit sizes.

Put the long singleton first in both folds. A Euclidean realization of the common reduced roots is

\[
\alpha_0=-\sqrt2e_1,\quad
\alpha_i=(e_i-e_{i+1})/\sqrt2\ (1\le i<n),\quad
\alpha_n=e_n/\sqrt2.
\]

This gives marks $(1,2,\ldots,2)$, with the longest root designated node 0. For $n=1$ the Cartan matrix is

\[
\begin{pmatrix}2&-1\\-4&2\end{pmatrix};
\]

for $n=2$ it is

\[
\begin{pmatrix}2&-1&0\\-2&2&-1\\0&-2&2\end{pmatrix}.
\]

These are the book's $A_{2n}^{(2)}$ convention. The transposed matrix is also isomorphic by reversing the chain. Consequently “the target or its transpose?” cannot distinguish the unlabelled algebra in this self-dual family. With node 0 and the root lengths fixed it distinguishes conventions. In particular, comparing the unrelabelled short four-end orbit of the $D$ fold with the book's long node 0 would produce a spurious disagreement. The book's $n_0=1$ convention is not Kac's $a_0=2$ convention mentioned in the brief.

### Masses and the missing scale identification

Let $P$ be the orthogonal projection onto the invariant field space. For invariant $u,v$,

\[
u^T\mathsf M^2v
=m_P^2\sum_O |O|n_O(u\cdot\alpha_O)(v\cdot\alpha_O).
\]

If $|O|n_O=c\,n_O^{\rm primitive}$, this is the primitive-label mass matrix at $m_F^2=c\,m_P^2$. The script builds the full parent Hessian and an orthonormal basis of the invariant space, and checks this matrix equality for every folding tested.

The following are squared masses, with the common displayed factor suppressed:

| $n$ | primitive $A_{2n}^{(2)}$, in units of $m_F^2$ | $A_{2n}^{(1)}$ restriction, in units of $m_A^2$ | $D_{2n+2}^{(1)}$ restriction, in units of $m_D^2$ |
|---|---|---|---|
| 1 | 3 | 3 | 6 |
| 2 | 1.3819660113, 3.6180339887 | same | 2.7639320225, 7.2360679775 |
| 3 | 0.7530203963, 2.4450418679, 3.8019377358 | same | 1.5060407926, 4.8900837358, 7.6038754716 |

These are $4\sin^2(\pi a/(2n+1))$, $a=1,\ldots,n$, for the primitive theory. Equality of absolute masses requires

\[
m_F^2=m_A^2=2m_D^2.
\]

The phrase “$m^2$ is replaced by $2m^2$” is itself correct: the primitive child obtained from a fixed $D$ parent has that larger coefficient. It does not make the two restrictions at identical parent $m$ equal.

Replace `part2-algebra/06-affine-algebras.qmd:116` with:

> Both foldings give the same primitive-label $a_{2n}^{(2)}$ potential. If $m_A$ and $m_D$ denote the mass parameters of the two parents and $m_F$ that of the folded theory, then $m_F^2=m_A^2=2m_D^2$. The classical mass ratios agree, and the absolute masses agree with this identification. Part V uses both (@sec-hirota-form-dn1-en1, @sec-one-loop-masses-folded).

### Other twisted rows and Part V

For the other four twisted rows, the script finds respectively the common factors $2,2,2,3$, equal to the twist orders. There is no further exception among the table's rows. Qualify line 89 by “For the twisted foldings listed here, other than $a_{2n}^{(2)}$” if the rule is intended to exclude other noncanonical reductions. The table does not classify every possible automorphism reduction.

Khastgir–Sasaki, `hep-th/9512158`, equation (4.7), explicitly list “$a_{2n}^{(1)}\to a_{2n}^{(2)}$”; Figure 12(iv), printed p. 23, explicitly lists “$d_{2n+2}^{(1)}\to a_{2n}^{(2)}$”. Both source diagrams were checked. Their numbering at the two $D$-fork tips differs harmlessly from the book. These quotations support the two parents; the normalization and spectra above were calculated independently.

Chapter 18, line 350, and `textbook_code/hirota/fold.py` use the same $A$-reflection as chapter 6. Chapter 19's route through `textbook_code/oneloop/run_twisted.py:7` uses the inverse four-cycle to the one printed in chapter 6. An automorphism and its inverse have the same orbits and invariant subspace, so this is correct. The oneloop routines retain parent masses rather than performing an explicit conversion to a shared $m_F$. Their mass-ratio comparison is unaffected. Chapter 18's phrase “the same classical masses” should reference the explicit scale identification above; no one-loop numerical result was re-refereed here.

## 3. Tables, twisted labels and rank

**Verdicts:** every numerical row of `tbl-untwisted` and `tbl-foldings` is **correct**. The twisted dual-label caption is **correct in the stated longest-root normalization**. The attribution “we follow MacKay–Watts” for the labels is **wrong**. The rank statements are **correct as field-space ranks**; the assertion about order-two diagram symmetries needs the $A_1$ exception.

### Untwisted labels

From reflection closure, the highest positive root supplies $n_1,\ldots,n_r$. Adjoining its negative gives a Gram matrix $G_{\rm aff}$. I verified $A n=0$, and independently $A^T n^\vee=0$ with $n_i^\vee=n_iG_{ii}/2$. The output records all marks, dimensions and exponents. The classical families were computed for $A_1$–$A_6$, $B_2$–$B_6$, $C_2$–$C_6$, $D_4$–$D_8$; all exceptional types were computed.

The arbitrary-rank formulas follow directly from the same Euclidean roots: $\theta=e_1-e_{n+1}$ for $A_n$, $e_1+e_2$ for $B_n,D_n$, and $\sqrt2e_1$ for the normalized $C_n$. Their expansions give the marks printed in the book, hence:

| type | $h$ | $h^\vee$ |
|---|---|---|
| $A_n$ | $n+1$ | $n+1$ |
| $B_n$ | $2n$ | $2n-1$ |
| $C_n$ | $2n$ | $n+1$ |
| $D_n$ | $2n-2$ | $2n-2$ |
| $E_6,E_7,E_8$ | $12,18,30$ | $12,18,30$ |
| $F_4$ | 12 | 9 |
| $G_2$ | 6 | 4 |

For $F_4$, in the printed chain order, $n=(1,2,3,4,2)$, lengths $(2,2,2,1,1)$, and $n^\vee=(1,2,3,2,1)$. For $G_2$, $n=(1,2,3)$, lengths $(2,2,2/3)$, and $n^\vee=(1,2,1)$. The $E_6$ three-leg description and the $E_7,E_8$ long-chain descriptions agree with the independently generated marks after relabelling.

### Twisted rows

Automorphisms were enumerated by exact preservation of the parent Gram matrix. Orbit averages were compared with independently constructed targets: coroot Cartan matrices of the untwisted duals, and the explicit $BC_n$ coordinates of item 2. Classical rows were checked at $n=2,3,4$, and the two $A_{2n}^{(2)}$ routes at $n=1,2,3$. Exceptional rows were checked in full. Representative output is:

| fold | orbit sizes | raw labels | primitive labels | squared lengths | primitive dual labels |
|---|---|---|---|---|---|
| $D_6^{(1)}\to A_5^{(2)}$ | 2,2,2,1 | 2,2,4,2 | 1,1,2,1 | 1,1,1,2 | 1,1,2,2 |
| $D_5^{(1)}\to D_4^{(2)}$ | 2,1,1,2 | 2,2,2,2 | 1,1,1,1 | 1,2,2,1 | 1,2,2,1 |
| $E_7^{(1)}\to E_6^{(2)}$ | 2,1,2,2,1 | 2,4,6,4,2 | 1,2,3,2,1 | 1,2,1,1,2 | 1,4,3,2,2 |
| $E_6^{(1)}\to D_4^{(3)}$ | 3,1,3 | 3,3,6 | 1,1,2 | 2/3,2,2/3 | 1,3,2 |
| $A_{2n}^{(1)}\to A_{2n}^{(2)}$ | 1,2,...,2 | 1,2,...,2 | same | 2,1,...,1,1/2 | 2,2,...,2,1 |

The order is the script's orbit order, not a new node convention for the book. In every case $k n_i|\alpha_i|^2/2$ is precisely the primitive positive integer null vector of $A^T$. This is stronger than checking its sum. The other non-simply-laced folds give the target, and are not isomorphic to its transpose; the self-dual $A_{2n}^{(2)}$ caveat is explained above. All $h,h^\vee$ pairs exchange on transposition, exactly as claimed.

For arbitrary rank, the $A_{2n-1}^{(2)}$ primitive labels are two fork labels 1, $n-2$ chain labels 2, and a last label 1. Their sum is $2n-1$; dual labels change the last 1 to 2, giving $2n$. The $D_{n+1}^{(2)}$ labels are all 1, with two end lengths 1 and $n-1$ interior lengths 2: $h=n+1$, $h^\vee=2n$. The $A_{2n}^{(2)}$ vectors in item 2 both sum to $2n+1$. These formulas supply the general-rank checks behind the samples.

### The source-normalization error

MacKay–Watts, §2, printed p. 2, choose the longest root to have squared length 2, but their defining equation is **$\sum_a n_a=kh$**. In §5.1 they explicitly call the resulting labels “no longer minimal”. The book instead divides out the common factor and uses $\sum n_a=h$. Thus the prior report's open attribution question has a definite answer: the book follows their root-length convention, not their twisted-label normalization.

At `part2-algebra/06-affine-algebras.qmd:79`, replace the final sentence with:

> We use coprime positive integer labels. MacKay and Watts [@mackay_watts1995] use the same longest-root normalization, but retain nonminimal labels in their twisted potentials; the resulting change of the mass parameter is explained in @sec-affine-folding.

At line 175, replace the second sentence with:

> The Kac labels $n_j$ are the coprime positive integers satisfying $\sum_j n_j\vec\alpha_j=0$.

At line 73, replace the diagram-symmetry sentence with:

> Nontrivial diagram automorphisms of order 2 exist for $A_n$ with $n\ge2$, $D_n$ with $n\ge4$, and $E_6$; $D_4$ also has automorphisms of order 3 (triality).

There is no order-two automorphism of the one-node finite $A_1$ diagram. Twisted field ranks $n,n,n,4,2$ follow from the number of simple-root orbits of the finite diagram automorphism; the affine Cartan matrices have one more node and corank 1.

There is also a verified wording clash with `01-introduction.qmd:128`: it normalizes the longest roots of the **finite** algebra to squared length 2. In the book's $A_{2n}^{(2)}$ convention, deleting the long affine node leaves lengths $(1,\ldots,1,1/2)$, or just $1/2$ for $n=1$. Its finite subsystem therefore does not have longest-root squared length 2. Chapter 6's summary, which normalizes the longest of all the affine simple roots, correctly describes its folding data. Replace the introduction's first normalization sentence by:

> Finite root systems are normalized so that their longest roots have squared length 2. For the affine Toda potentials, the longest affine simple roots have squared length 2, as in @sec-lie-algebraic-data; for $a_{2n}^{(2)}$ this includes the longest affine root at node 0.

## 4. Gradations, exponents and vertex operators

**Verdicts:** principal degrees and their stated periods/examples **correct**; the unqualified $2T$ formula **wrong for twisted algebras in the book's coupling**; fundamental levels **correct for the canonical central element**; vertex construction **correct in its simply-laced scope**, but the principal attribution is incomplete; OTU time dependence **correct after a species-dependent phase convention for $z$**; the Kneipp–Olive reference is **too broad**, since they exclude $A_{2n}^{(2)}$.

### Explicit matrix centralizers

For $A_2^{(2)}$, use

\[
E=E_{12}+E_{23}+\lambda E_{31},\qquad
\sigma(X)=-JX^TJ^{-1},\quad
J=\begin{pmatrix}0&0&1\\0&-1&0\\1&0&0\end{pmatrix}.
\]

Then $\sigma(E(\lambda))=E(-\lambda)$, $E^3=\lambda I$, and the principal degree of $\lambda$ is 3. The matrix is cyclic over $\mathbb C(\lambda)$; its traceless centralizer is spanned by $E,E^2$ over that field. Since $\sigma(E^r)=(-1)^{r+1}\sigma(E)^r$, the twisted condition on $\lambda^mE^r$ is $m\equiv r+1\pmod2$. Its degree is $3m+r$. The allowed residues are therefore 1 and 5 modulo 6. The companion-matrix cyclic vector also shows that a Laurent-polynomial commuting matrix has Laurent-polynomial coefficients in this basis, so no extra degrees arise from rational denominators.

The identical construction in $\mathfrak{sl}_5$, with $E^5=\lambda I$, gives $A_4^{(2)}$: residues $1,3,7,9\pmod {10}$. These are explicit matrix-realization computations for the two twisted algebras requested by the brief.

I also checked the stated $D_4^{(3)}$ example directly. The script constructs the 28-dimensional $\mathfrak{so}_8$ algebra as $8\times8$ matrices, propagates triality from its Chevalley generators by commutators, and verifies $\sigma^3=1$. With $\deg\lambda=4$, it solves $[E_1,X]=0$ degree by degree, imposing $X_m\in g_{(m\bmod3)}$. Output:

```text
PASS A2^(2) in sl(3): residues=[1, 5] mod 6
PASS A4^(2) in sl(5): residues=[1, 3, 7, 9] mod 10
PASS D4^(3) in twisted so(8) loop:
centralizer dimensions degrees 0..11 = [0,1,0,0,0,1,0,1,0,0,0,1]
```

Multiplication by $\lambda^3$ proves the period 12 in the last calculation. The word “period” need not mean the least possible period of the set of degrees: for some families that set is just the odd integers.

The central-extension statement follows directly from the cocycle: if loop-centralizer elements commute, their affine bracket has only the cocycle term, which pairs opposite degrees. This supports the wording at line 151.

### The gradation formula

From the chapter 20 charge spin $s_j=8\pi/(\beta^2|\alpha_j|^2)-1$, independent summation gives

\[
\sum_j n_j^\vee s_j
=\begin{cases}
4\pi h/\beta^2-h^\vee,&\text{untwisted},\\
4\pi kh/\beta^2-h^\vee,&X^{(k)}\text{ twisted}.
\end{cases}
\]

Chapter 21, line 29, explicitly uses $\beta_{\rm lit}^2=\beta^2/k$ for the second case. Chapter 6, line 153, does not. For $D_4^{(3)}$, the two expressions are $48\pi/\beta^2-6$ and $16\pi/\beta^2-6$, respectively; this is not a harmless additive convention.

Replace `part2-algebra/06-affine-algebras.qmd:153` with:

> For the non-local charges of the imaginary-coupling theory, their Lorentz spins determine the gradation parameter $T$ of @sec-amplitude-and-gradation. For untwisted $\hat g$ this gives $2T=4\pi h/\beta^2-h^\vee$. For twisted $\hat g=X^{(k)}$, the formula uses $\beta_{\rm lit}^2=\beta^2/k$; in the normalization of @sec-lie-algebraic-data it reads $2T=4\pi kh/\beta^2-h^\vee$.

This replacement follows the spectral-parameter convention of chapter 21. A separate issue occurs for the book-labelled $A_{2n}^{(2)}$: its dual label at node 0 is 2, not 1. Chapter 20's sentence normalizing $n_0^\vee=1$ “since $x$ multiplies $e_0$ once” therefore cannot be applied to it without an additional spectral-variable convention. The dual-label calculation in item 3 verifies that obstruction. I have not re-refereed a particular $A_{2n}^{(2)}$ amplitude to choose that convention.

### Levels, attributions and OTU

Define $c=\sum_i n_i^\vee h_i$. Since $A^Tn^\vee=0$, $[c,e_j]=0$; since $\Lambda_j(h_i)=\delta_{ij}$, $\Lambda_j(c)=n_j^\vee$. This is the independent level calculation. For simply-laced untwisted algebras $n_j^\vee=n_j$, so the stated level-one criterion is correct. It assumes integrable dominant fundamental modules, not arbitrary nonintegrable highest-weight modules.

OTU, `hep-th/9305160`, introduction, describe their construction as the “principal vertex operator construction [21], [22]”. Their reference [21] names **Kac, Kazhdan, Lepowsky and Wilson**, and [22] is Kac's book. This directly supports expanding the principal attribution. The authors and title of the 1981 paper are also confirmed by its [Rutgers publication record](https://www.researchwithrutgers.org/en/publications/realization-of-the-basic-representations-of-the-euclidean-lie-alg/). A safe replacement for the final sentence of line 157 is:

> These are the Frenkel–Kac–Segal construction in the homogeneous gradation and the principal construction of Kac, Kazhdan, Lepowsky and Wilson, extending the work of Lepowsky and Wilson [@kac1990].

For the time dependence, OTU equation (2.29a), visually checked on printed p. 9, is

\[
W(i,z,x^\mu)=\exp\!\left(-\mu\{\gamma_i\cdot q(1)zx^+
-\gamma_i\cdot q(1)^*z^{-1}x^-\}\right).
\]

Write $p_i=\gamma_i\cdot q(1)=|p_i|e^{i\varphi_i}$ and define $z_{\rm book}=e^{i\varphi_i}z_{\rm OTU}$. Then both coefficients become $|p_i|$, so this is exactly the chapter's $W$, with $m\mu_i=\mu|p_i|$. OTU (2.30) makes $z_{\rm book}=\mp e^\eta$. The minus sign between the two light-cone terms is correct. Their (2.28) also has a central factor $e^{\mu^2m_jx^+x^-}$; it cancels in the ratio in (2.27), leaving the stated matrix-element tau functions. Thus this omission is legitimate, not a missing spacetime dependence. Nilpotency is explicitly treated in OTU §§5–7 for integrable modules of untwisted affine algebras.

To make the convention and the directly checked source clear, replace the parenthesis in line 171 by:

> (we absorb the species-dependent phase of the Coxeter eigenvector into $z_k$; the relative sign in $W_k$ is that of [@olive1993b, eq. (2.29)])

Kneipp–Olive, `hep-th/9404030`, introduction after (1.2), explicitly say of $A_{2n}^{(2)}$ that it “needs special treatment in both approaches and so will be ignored as we seek general arguments.” The $A_2^{(2)}$ computation above also shows the exceptional unequal root lengths independently. Replace the next scope sentence at line 171 by:

> Olive, Turok and Underwood treat untwisted $\hat g$. Kneipp and Olive [@kneipp1996] extend the construction to the twisted families other than $a_{2n}^{(2)}$, which they exclude from their treatment.

The bibliography/reference check for `olive1993a` itself remains unclear: I could not read that paper. The sign is verified from `olive1993b`, which is enough to correct the current citation.

## 5. Affine-diagram symmetries and charges

**Verdicts:** the semidirect-product and simple-transitivity assertions are **correct**; the sentence that the symmetries act by phases is **wrong for the whole group**; the $A_n$ cyclic charge and weight-class statements are **correct**.

The script enumerates automorphisms of the full directed Cartan matrix. It then constructs a normal abelian subgroup acting regularly on the mark-1 nodes and checks its intersection with the node-0 stabilizer. For even $D_n$, the subgroup is explicitly required to have exponent 2: an arbitrary regular normal subgroup of the same order would not identify the centre correctly.

| family | affine automorphism group | stabilizer of 0 | centre subgroup |
|---|---|---|---|
| $A_1$ | $C_2$ | 1 | $C_2$ |
| $A_n,n\ge2$ | dihedral, order $2(n+1)$ | $C_2$ | $C_{n+1}$ |
| $B_n,C_n$ | $C_2$ | 1 | $C_2$ |
| $D_4$ | $S_4$ | $S_3$ | $C_2\times C_2$ |
| $D_n,n\ge5$, even | order 8 | $C_2$ | $C_2\times C_2$ |
| $D_n,n\ge5$, odd | order 8 | $C_2$ | $C_4$ |
| $E_6$ | $S_3$ | $C_2$ | $C_3$ |
| $E_7$ | $C_2$ | 1 | $C_2$ |
| $E_8,F_4,G_2$ | 1 | 1 | 1 |

These counts were obtained for the ranks listed in item 3. For the classical families at arbitrary rank, the cycle rotations/reflections and the two-fork permutations give the displayed groups directly. The centre subgroup is normal, the stabilizer is a complement, and their orders multiply to the full group order. Thus the statement is $\mathrm{Aut}(\widehat\Delta)\cong Z(G)\rtimes\mathrm{Aut}(\Delta)$, with the normal factor understood in that order. $A_1^{(1)}$ has the two-node exceptional diagram rather than an ordinary cycle of length at least three, but its order-two symmetry still gives the stated centre.

The phase statement fails already for $A_2$. For $N=n+1$, let $v_a(j)=e^{2\pi iaj/N}$. Rotation of the affine nodes sends $v_a$ to a scalar multiple of itself, but reflection $j\mapsto-j$ sends $v_a$ to $v_{N-a}$. For $N\ge3$, those are in general distinct particle species. A nonabelian full diagram group cannot generally be represented only by simultaneously diagonal phases. `check_misc.py` verifies both actions for $N=3,\ldots,8$.

Replace the start of the final part of `part2-algebra/06-affine-algebras.qmd:63` with:

> The centre subgroup can be diagonalized on the particle states, so the particles of $a_n^{(1)}$ carry a $\mathbb Z_{n+1}$ charge. Other diagram symmetries can permute species; for $a_n^{(1)}$ the reflection exchanges $a$ and $n+1-a$.

Retain the soliton sentence. Independently, every weight of $\Lambda^a\mathbb C^N$ has the form

\[
\sum_{j\in I}e_j-\frac aN\sum_{j=1}^Ne_j,\qquad |I|=a.
\]

Differences of such weights lie in the $A_{N-1}$ root lattice, and $\lambda_a-a\lambda_1$ is an integer vector of zero sum. Therefore their common class is $a\pmod N$, exactly as stated. For a species-$a$ Hirota soliton, the phase increments in the tau functions are $2\pi a/N$ modulo $2\pi$. Choosing continuous logarithms yields a zero-sum charge vector whose coordinates have the two values $1-a/N$ and $-a/N$, with $a$ occurrences of the former, up to the sign/orientation of the soliton. Hence it is one of these weights. This does not assert that a single-soliton family attains every weight; the chapter carefully says “among”.

## 6. Perron–Frobenius vectors and higher charges

**Verdicts:** the eigenvalue formula, the above-$h$ colour change and the repeated $D_{2n}$ exponent are **correct in the simply-laced scope of this section**. The $A_n$ and $E_8$ mass statements pass independent checks. The minimal-factor interpretation of the Ising statement is **correct as a block identity**; the precise original-paper attributions are partly **unclear**, as specified below.

### Conservation equations, not just Cartan diagonalization

The independent computation first constructs the bicoloured Coxeter element from Euclidean reflections. With signed simple-root representatives $\gamma_a$, it enumerates every triple of orbit roots $u,v,w$ satisfying $u+v+w=0$. Their Coxeter-plane projections supply the physical triangle angles $\phi_u,\phi_v,\phi_w$. At spin $s$, I solve the real and imaginary parts of

\[
q_a e^{is\phi_u}+q_b e^{is\phi_v}+q_c e^{is\phi_w}=0
\]

over all such triples. These equations are assembled before any Cartan eigenvector is used. Only afterwards is the nullspace compared to the claimed Cartan eigenspace and to the colour-flipped lower-spin space.

| algebra | $h$ | above-$h$ spin | root triangles | nullity | Cartan residual |
|---|---|---|---|---|---|
| $A_4$ | 5 | 6 | 20 | 1 | $2.22\times10^{-15}$ |
| $D_5$ | 8 | 9 | 80 | 1 | $9.42\times10^{-16}$ |
| $E_6$ | 12 | 13 | 240 | 1 | $1.67\times10^{-15}$ |
| $D_4$ | 6 | 9 | 32 | 2 | $7.14\times10^{-16}$ |
| $D_6$ | 10 | 15 | 160 | 2 | $5.40\times10^{-16}$ |

In the first three rows the lower spin is 1; in the last two it is the repeated exponent 3 or 5. The norm of the difference between the upper-spin projector and the colour-conjugated lower-spin projector is below $1.8\times10^{-15}$ in all five cases. Thus the reading pass's above-$h$ correction is confirmed by the requested independent fusing computation.

For completeness, the all-spin mechanism is visible in these equations. Coxeter orbit shifts change a root's projected angle by $2\pi/h$; their $e^{ih\phi}$ factor is unchanged. For the two signed simple-root colours, the corresponding $h$-fold phases differ by a minus sign, up to a common phase. Multiplying $q_a$ by the colour sign therefore maps the complete spin-$s$ conservation system to the spin-$s+h$ system. If $J$ is the colour-sign diagonal matrix, $J(2-G)J=4I-(2-G)$, so the eigenvalue changes from $4\sin^2(\pi s/(2h))$ to $4\sin^2(\pi(s+h)/(2h))$. The two derivations agree. At repeated exponents this is a statement about the charge space; an individual basis is not canonical.

For $A_n$, the elementary sine recurrence gives all $n$ eigenvectors $v_a^{(s)}=\sin(as\pi/(n+1))$, $s=1,\ldots,n$, and immediately gives $v_a^{(s+h)}=(-1)^av_a^{(s)}$. For $D_{2n}$, the exponent list contains $2n-1$ twice: once among the odd exponents and once as the extra $r-1$ exponent. The fusion calculation verifies that the degeneracy is retained in the charge equations for $D_4,D_6$.

### Masses and the Ising comparison

The script checks the entire finite Cartan spectrum against the exponent formula, for all types/ranks in item 3. An independent Hessian calculation also compares classical masses with the positive eigenvector for $A_4,D_5,E_6,E_7,E_8$. For $E_8$ the sorted mass ratios are

```text
1, 1.6180339887, 1.9890437907, 2.4048671724,
2.9562952015, 3.2183404585, 3.8911568233, 4.7833861168
```

The golden-ratio residual is $3.33\times10^{-15}$. For $A_n$, the cyclic Hessian is the discrete Laplacian and has eigenvalues $4m^2\sin^2(\pi a/h)$, proving the printed mass formula, including $A_1$.

The minimal $E_8$ $S_{11}$ expression can be checked without the author's script:

\[
\prod_{x\in\{1,11,19,29\}}(x-1)_{30}(x+1)_{30}
= f_{2/3}(\theta)f_{2/5}(\theta)f_{1/15}(\theta),\quad
f_a=\frac{\sinh\theta+i\sin\pi a}{\sinh\theta-i\sin\pi a}.
\]

Pair the blocks $x$ and $30-x$, using $(x)_{30}(30-x)_{30}=-f_{x/30}$, and $(30)_{30}=-1$. This proves the identity; `check_misc.py` also checks it at three rapidities. It is the standard Ising minimal amplitude. The complete Toda amplitude has additional coupling-dependent factors; its zero-coupling limit is 1, not this minimal amplitude. The chapter's wording preserves that distinction.

CDS §1 explicitly state that masses and charge eigenvalues “are known to be components of the eigenvectors of the adjacency matrix of the simple roots”. Dorey's accessible [sequel, §2](https://arxiv.org/pdf/hep-th/9110058) states that the orbit sets “are disjoint, each has h elements, and thus their union is all of Φ”, and describes the charge-bootstrap discussion in the earlier paper. That is supporting primary evidence, but it is **not a reading of `dorey1991` itself**. The full Klassen–Melzer paper, Freeman 1991, Braden et al. 1990 and Zamolodchikov 1989 were not available in the supplied source set or obtained in readable full text. I therefore leave the exact attribution/priority claims to those papers **unclear**, rather than inherit the earlier referee's certification. No mathematical change to line 166 is recommended. Adding “For these simply-laced theories” at its start would make its scope harder to misread.

A cross-chapter inconsistency encountered in the required reading: chapter 2, line 61, and chapter 9's concluding charge-eigenvector sentence do not retain the simply-laced restriction. DGZ §9 explicitly explain that the quantum non-simply-laced charges no longer obey the undeformed incidence-matrix rule when $h$ is replaced by the floating $H$. Those sentences should have the same scope qualification.

## 7. Other chapter 5 statements and figures

**Verdicts:** normalization, Chevalley relations, tables, weights, Casimir, and the two figure captions are **correct**, subject to the small precision points below. The universal assertion that long- and short-root orbit polygons have different radii is **wrong**. A Coxeter plane also requires rank at least two.

### Killing form, root vectors and Casimir

With the normalized invariant form $B$, invariance gives

\[
B([E_\alpha,E_{-\alpha}],H)
=\alpha(H)B(E_\alpha,E_{-\alpha}).
\]

Choosing $B(E_\alpha,E_{-\alpha})=2/|\alpha|^2$ therefore gives exactly the printed coroot commutator. It is compatible with $h_i=\alpha_i^\vee\cdot H$, since $[h_i,e_j]=a_{ij}e_j$. Root-vector normalization is a choice, not a contradiction with orthonormal Cartan generators.

For $\mathfrak{sl}_2$, $H=\sqrt2J^3$, $E_\pm=J^\pm$ gives root $\sqrt2$, $[E_+,E_-]=\sqrt2H=2J^3$, and $[H,E_\pm]=\pm\sqrt2E_\pm$. Thus the physics $su(2)$ remark is correct up to the implicit rescaling. Literally, inside a **complex** simple algebra the complex span of $E_\pm,h$ is $\mathfrak{sl}_2(\mathbb C)$; its compact real form is $su(2)$. For mathematical precision, replace “span an $su(2)$ subalgebra” at `05-root-systems.qmd:29` by “span an $sl_2(\mathbb C)$ subalgebra, whose compact real form is $su(2)$”.

If the invariant-form Casimir is written in dual root bases, its positive-root contribution on a highest-weight vector is $\sum_{\alpha>0}\alpha\cdot\lambda=2\rho\cdot\lambda$; the Cartan contribution is $\lambda^2$. This derives $C(\lambda)=\lambda\cdot(\lambda+2\rho)$ in precisely the chosen normalization. In $su(2)$ it is $2j(j+1)$, which is consistent with long-root squared length 2. Reflection of the positive-root sum by a simple reflection removes twice that simple root, giving $\alpha_i^\vee\cdot\rho=1$ and hence $\rho=\sum_i\lambda_i$. The script verifies the equality with the half-sum of positive roots and $C(\theta)=2h^\vee$ for every tested type.

### Classification, Coxeter elements and the radius counterexample

`tbl-simple-lie` agrees with all the generated dimensions, exponents and Coxeter numbers. Formula $d=r(h+1)$ holds for every row. The classical Euclidean root lists prove the corresponding dimension formulas for arbitrary rank; exceptional root counts are $72,126,240,48,12$.

For all cases tested, a Coxeter product has order $h$, the stated eigenvalues, and exactly $r$ root orbits of size $h$. The usual all-order conjugacy argument can be seen directly on a tree: commute nonadjacent generators and conjugate a leading simple reflection to the end, which reverses a source to a sink in the induced orientation. All tree orientations are connected by these moves. This establishes conjugacy of the Coxeter products beyond the sampled numerical computations. The exponent symmetry follows by complex conjugation of their eigenvalues.

However, the parenthesis in `part2-algebra/05-root-systems.qmd:119` is false. For $B_3$, take the signed cyclic Coxeter action

\[
e_1\mapsto e_2,\quad e_2\mapsto e_3,\quad e_3\mapsto-e_1.
\]

The omitted one-dimensional eigenspace is spanned by $(1,-1,1)$. Orthogonal projection onto the other two-dimensional plane gives

\[
|P e_1|^2=1-1/3=2/3,
\quad |P(e_1-e_2)|^2=2-4/3=2/3.
\]

Here $e_1$ is short and $e_1-e_2$ is long. Their entire Coxeter orbits have the same radius. Another long-root orbit, represented by $e_1+e_2$, has squared radius 2. This disproves both a strict long/short separation and the reading-pass note's stronger assertion that there are always just two radii corresponding to the two root lengths. The script's unit complex-eigenvector coordinates differ by a common factor $\sqrt2$, and print:

```text
B3 root length2 2: radii {0.57735027, 1.0}
B3 root length2 1: radii {0.57735027}
```

Replace the relevant sentence at line 119 with:

> For rank at least two, the bicoloured Coxeter element acts on its Coxeter plane as a rotation by $2\pi/h$ [@steinberg1959]. Each root orbit projects to a regular $h$-gon. Distinct orbits can have the same radius, even when their roots have different lengths. For $A_1$, the action is a reflection on a line.

The mathematical replacement is established by the computation and the elementary eigenspace construction. The exact historical attribution to Steinberg remains subject to the source-access limitation in item 8.

The principal $\mathfrak{sl}_2$ claim has an additional independent check: group roots by height. If $b_j$ is the number of roots of positive height $j$, the multiplicity of a spin-$j$ summand under the principal $\mathfrak{sl}_2$ is $b_j-b_{j+1}$. For $E_6,E_7,E_8$, this recovers the printed exponents; `check_misc.py` records them. This checks the mathematical relation, but does not replace reading Kostant for the attribution.

### Lattices and figures

For $\ell$ in the coweight lattice, $\alpha_i\cdot\ell\in\mathbb Z$. The extra untwisted root is an integer combination of the finite simples, so every potential term is invariant under $\phi\mapsto\phi+2\pi\ell/\beta$. The same reasoning applies to the book's twisted node convention when the label at the chosen affine node is 1. Thus differences of these lattice vacua give the stated normalized topological charges. This does not assert that every lattice difference is realized by an individual classical soliton.

For the displayed $G_2$ roots, $\alpha_1=(\sqrt{2/3},0)$ and $\alpha_2=(-\sqrt{3/2},1/\sqrt2)$. The chamber inequalities are $x>0$ and $y>\sqrt3x$; its angular range is $60^\circ<\arg(x+iy)<90^\circ$, exactly the shaded wedge in the committed SVG. The highest root is $(0,\sqrt2)=3\alpha_1+2\alpha_2$. All six mirrors and both root lengths agree with the caption.

For $E_8$, the independent bicoloured Coxeter calculation gives eight distinct orbit radii, 30 roots per orbit, and radii proportional to the corresponding positive Cartan eigenvector entries. The numerical ratios are those in item 6. The figure caption is correct. The $C_2$ and $G_2$ folding figure also has the correct orbit labels, bond directions, lengths and sums. There is no need to regenerate either SVG.

## 8. Prior art, attribution and citation audit

**Verdicts:** no claim of an original Part II result was found; the absence of a Part II item in Appendix E is **correct**. Citation correctness is **mixed**, with the verified errors identified above and specific inaccessible sources left **unclear**.

Independent audit: I read every citation-bearing paragraph in the two chapters and searched their text for novelty/status wording. The mathematical content is classification, standard representation theory, or explicitly attributed Toda material; there is no sentence presenting a known formula as a new result of the book. The folding and matrix-element formulas were independently checked in items 1–7. The Appendix E search produced no Part II entry. This is an audit of these chapters, not a priority claim about the entire book.

| citation | result of source comparison |
|---|---|
| `mackay_watts1995` | **Wrong for the label-normalization attribution.** The primary source uses $\sum n=kh$; see item 3. **Correct** source for longest-root normalization, classical Lie duality, and the normalized affine strong–weak coupling relation. |
| `dgz1992` | **Correct** source for the moving mass/S-matrix structure and classical-family strong–weak exchanges, §9. It is Delius–Grisaru–Zanon. It should not be substituted by the Delius–Gould–Zhang R-matrix paper. |
| `corrigan1993` | **Correct** for coroot pairing and exceptional floating S-matrices; **not a derivation of exceptional $B(\beta)$**. Their explicit non-claim is quoted in item 1. |
| `khastgir1996b` | **Correct** for both stated $A_{2n}^{(2)}$ parent diagrams, equation (4.7) and Figure 12(iv). |
| `olive1993b` | **Correct** for untwisted tau functions, time factors, integrable-module nilpotency and principal vertex operators; equations (2.27)–(2.30), §§5–7. The rephasing convention should be stated. |
| `kneipp1996` | **Correct only after excluding $A_{2n}^{(2)}$**. Their introduction explicitly excludes it; see item 4. |
| `olive1993a` | **Unclear:** full text not read. The time-factor signs can instead be cited to the directly checked `olive1993b`. |
| `olive1983` | **Partly checked:** the publisher's abstract states that the reduced equations “are shown to correspond precisely to the generalized Cartan matrices classified in the mathematical literature.” This supports the reduction topic. The precise historical assertion “also showed that every non-simply-laced theory arises this way” remains **unclear** without the full paper. Independently, the table's constructive folds establish the mathematical coverage. [Publisher abstract](https://www.sciencedirect.com/science/article/pii/0550321383902560). |
| `dorey1991`, `klassen1990` | **Unclear as exact citations:** the full cited papers were not read. Dorey's accessible sequel supports the root-orbit and charge-bootstrap framework; it is a different paper. The [Klassen–Melzer abstract](https://www.sciencedirect.com/science/article/pii/055032139090643R) discusses the minimal/Toda distinction but does not certify the specific repeated-exponent sentence. |
| `freeman1991`, `braden1990`, `zamolodchikov1989` | **Unclear as exact attribution/priority checks:** full source passages were not obtained. Independent numerical and algebraic checks support the mass and minimal-factor formulas, not historical priority. The book no longer assigns the higher-charge theorem to Freeman, which is the right separation of claims. |
| `steinberg1959`, `kostant1959` | **Unclear as exact source-line attributions:** original full texts were not obtained; the AMS Steinberg PDF returned HTTP 403. The mathematics was independently checked, including the counterexample to the added radius parenthesis. |
| `kac1990`, `goddard1986` | **Unclear under the brief's source-line standard:** book/review metadata and citations from OTU identify the expected subjects, but I did not obtain the relevant full pages. Do not treat that as a verified quotation from either source. |
| `humphreys1972`, `fuchs1997`, `kassel1995` | **Unclear as source-location checks:** no relevant full pages read. The introductory recommendations do not assert a new result. |

The source-based principal-construction correction is given at item 4. It is supported by OTU's own citation [21], rather than merely repeating the earlier referee's advice. No correction to a published mathematical formula is proposed, and no Appendix D change is warranted by this report.

## 9. Exercises and scope

**Verdict:** all ten exercises are **correct and solvable**, with the folded-mass scale clarification below. Here are independent solutions in outline.

1. **`exr-a2-roots`.** The six roots are $\pm\alpha_1,\pm\alpha_2,\pm(\alpha_1+\alpha_2)$. Fundamental weights are $\lambda_1=(1/\sqrt2,1/\sqrt6)$, $\lambda_2=(0,\sqrt{2/3})$. The Cartan matrix is $\left(\begin{smallmatrix}2&-1\\-1&2\end{smallmatrix}\right)$, and $s_1s_2=\left(\begin{smallmatrix}-1/2&-\sqrt3/2\\sqrt3/2&-1/2\end{smallmatrix}\right)$. Its order is 3 and its eigenvalues correspond to exponents 1,2.
2. **`exr-g2`.** Short-root squared length $2/3$, short-first Cartan matrix $\left(\begin{smallmatrix}2&-3\\-1&2\end{smallmatrix}\right)$, highest root $3\alpha_1+2\alpha_2$, marks $(3,2)$, dual marks $(1,2)$. Thus $h=1+3+2=6$, $h^\vee=1+1+2=4$.
3. **`exr-d4-masses`.** Give the three leaves component $a$ and the centre component $b$. The equations are $(2-\lambda)a=b$, $(2-\lambda)b=3a$. Positivity selects $2-\lambda=\sqrt3$, hence $b=\sqrt3a$. The independent script gives three 1's and $\sqrt3$.
4. **`exr-cartan-eigenvalues`.** For $s=1,\ldots,n$, $v_a=\sin(as\pi/(n+1))$ obeys both zero boundary conditions, and the sine addition formula gives eigenvalue $2-2\cos(s\pi/(n+1))$. The distinct eigenvalues supply a complete basis.
5. **`exr-roots-count`.** Substituting each table row gives $d=r(h+1)$. For the exceptionals: $78=6\cdot13$, $133=7\cdot19$, $248=8\cdot31$, $52=4\cdot13$, $14=2\cdot7$. Subtracting the rank gives $|\Phi|=rh$.
6. **`exr-jacobi-central`.** Set $K(X,Y,Z)=B([X,Y],Z)$. Symmetry and invariance of $B$ make $K$ cyclic. In the Jacobiator of $X\lambda^m,Y\lambda^n,Z\lambda^p$, the central coefficient is $\delta_{m+n+p,0}[(m+n)+(n+p)+(p+m)]K=0$. The noncentral coefficient vanishes by the finite-algebra Jacobi identity. Symmetry is also needed for antisymmetry of the cocycle.
7. **`exr-affine-cartan`.** $A_2^{(1)}$ has diagonal 2 and all off-diagonal entries $-1$, annihilating $(1,1,1)^T$. In the long-short-long order, $C_2^{(1)}$ has
   \[
   A=\begin{pmatrix}2&-1&0\\-2&2&-2\\0&-1&2\end{pmatrix},\qquad A(1,2,1)^T=0.
   \]
   These use the chapter's row-coroot convention, not its transpose.
8. **`exr-fold-d4`.** In standard $D_4$ coordinates take $\alpha_0=-e_1-e_2$, central $\alpha_2=e_2-e_3$, and average the other three end roots: $u=(e_1-e_2+2e_3)/3$. Lengths are $2,2,2/3$, adjacent inner products $-1,-1$, and $\alpha_0+2\alpha_2+3u=0$. Thus the chain Cartan matrix is $\left(\begin{smallmatrix}2&-1&0\\-1&2&-1\\0&-3&2\end{smallmatrix}\right)$, the stated $G_2^{(1)}$.
9. **`exr-fold-mass`.** The bilinear-form derivation in item 2 proves the claim for every orbit fold, including orbits with adjacent roots. In $A_3^{(1)}$, the Hessian eigenvalues are $2,2,4$ times $m^2$; the invariant two-dimensional restriction has $2,4$, giving $C_2^{(1)}$ masses $\sqrt2m,2m$. The mass-subset statement includes multiplicities only for the retained invariant modes. For $A_{2n}^{(2)}$, the same claim holds with raw orbit labels; with primitive labels the folded parameter is $m_F^2=c\,m_P^2$. In particular $c=2$ for its $D$ parent, not 1.
10. **`exr-dual-coxeter`.** For $B_n$, the long-label contribution is $1+1+2(n-2)$ and the short contribution is 1, totaling $2n-1$. For $C_n$, two long ends contribute 2 and $n-1$ short nodes contribute $n-1$, totaling $n+1$. The $F_4,G_2$ sums are 9 and 4 as computed in item 3.

At `part2-algebra/06-affine-algebras.qmd:199`, the following addition would remove the only substantive exercise ambiguity:

> If the orbit labels have a common factor $c$ and are divided by it, use the folded mass parameter $m_F^2=c\,m_P^2$, as in @sec-affine-folding.

There is no new mathematical result here needing an originality label. The exact strong–weak S-matrix statement goes beyond elementary Lie theory and should retain an explicit conjectural/proposed status, consistent with chapter 14. The two-parent normalization calculation, dual-coupling conversion, twisted-period examples and higher-spin colour identity would benefit from cited author-maintained checks under the repository's check-script rule. The existing figure scripts are not cited in these chapters. The editor can have the drafting agent add appropriate paths or independently reproduce these checks in `textbook_code/`; referee scratch paths must not be cited in the book. No [Theorem] upgrade is recommended.

## Claims not fully checked

- Exact source-line attribution/priority for Freeman 1991, Braden–Corrigan–Dorey–Sasaki 1990, Zamolodchikov 1989, Dorey 1991, Klassen–Melzer 1990, Steinberg 1959 and Kostant 1959: readable full source passages were not obtained. The corresponding mathematical statements were checked to the extent detailed above; those checks do not certify attribution.
- The full scope of Olive–Turok 1983: only its primary publisher abstract was accessible. Constructive coverage of the table's families was checked independently.
- Original `olive1993a` sign conventions: the source was unavailable. The displayed time factor was checked directly against `olive1993b` (2.29).
- The exact relevant pages of Kac, Goddard–Olive, Humphreys, Fuchs–Schweigert and Kassel: not read. Standard subject coverage or a bibliographic record is not enough under this brief's quotation standard.
- An all-orders derivation from the Lagrangian of the exceptional $B(\beta)$ conjectures, and a complete proof of quantum non-local Serre relations: not supplied by the checks here. The coroot and coupling algebra is verified.
- A fresh audit of Part V one-loop transmission factors and absolute mass corrections: outside the brief's chapter scope. Its node maps and use of parent masses were inspected; only the classical scale conversion was independently calculated.
- The $A_{2n}^{(2)}$ evaluation-parameter convention needed to reconcile chapter 20's $n_0^\vee=1$ sentence with the book's longest-node-0 labeling: the mismatch is verified, but selecting a corrected amplitude convention requires reviewing the relevant Part VI construction.

## Other observations — unverified

- The word “centre” identifies the abstract finite abelian group correctly. If a later argument needs a canonical identification with a weight-lattice quotient rather than its character group, the pairing should be specified. No error in these chapters is asserted on that basis.
- The generality of the phrase “Physical quantities do not depend on the choice” of gradation was not tested for singular parameter-dependent basis changes or special quantum-group parameters. The ordinary generic evaluation-representation interpretation is the one used in this report.
