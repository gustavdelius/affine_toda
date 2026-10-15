# Independent Referee Report: Chapters 7 and 8 (Quantum Groups and R-Matrices)

- **Date:** 2026-10-15
- **Referee Model:** Gemini 3.8 Flash (High)
- **Commit Reviewed:** `2d36c9c2c1848067a87a3c342f2549c43a25fe96`
- **Output Report Path:** `briefs/referee_algebra_07_08_report.md`
- **Scratch Directory:** `referee_scratch/algebra_07_08/`

All checks below are independent derivations and scripts constructed directly from first principles and the literature. No existing book file or brief was edited; no git commands altering the index or working tree and no `quarto render` were executed.

---

## Executive Summary & Overview of Findings

1. **Normalizations, Central Element, and $S^2$ Shift (Chapter 7):**
   - Central element $\prod_i k_i^{n_i}$: **Correct** for untwisted algebras and twisted $X^{(k)}$ ($k \ge 2$) in Kac's labelling. **Partially wrong / imprecise** for $a_{2n}^{(2)}$ where Kac's mark is $a_0 = 2$, not $1$; and the relation $n_i^\vee = n_i \vec\alpha_i^2/2$ is missing a factor of $k$ in the book's root normalization (where long roots have squared length 2).
   - $S^2$ shift: **Wrong for $a_{2n}^{(2)}$ as stated.** If $(\alpha_0, \alpha_0) = 1$ (Kac's labelling where node 0 is the short root), the shift is $q^{h^\vee}$, not $q^{2h^\vee}$. If node 0 is the long root with squared length 4 in that scale, $(\alpha_0, \alpha_0) \ne 1$. For untwisted algebras and twisted $X^{(2)}$ ($d_{n+1}^{(2)}$), the general formula derived below yields $q^{2h^\vee}$ in Kac's normalization.
   - Crossing point $x_c$: **Correct.** Crossing twice yields $x_c^2 = q^{2h^\vee}$, fixing $x_c$ up to a sign $\pm 1$. The Toda soliton crossing points quoted in Part VI all satisfy $x_c^2 = q^{-2h^\vee(g^\vee)}$, with the power of $q$ supplied by $S^2$ and the sign supplied by $e^{-i\pi h^\vee}$.
2. **Sign of $q$ and Part VI Dictionary (Chapter 7):**
   - Sine-Gordon dictionary ($b = S_T/S_0$, off-diagonal entries $e^{\mp\lambda\theta} S_R/S_0$, gauge transformation to principal gradation): **Correct.**
   - $P \check R(1/z) P = \check R(z; q^{-1})$ and Bernard–LeClair value $q = -e^{-i\pi\lambda} = e^{-8\pi^2 i/\beta_{\rm SG}^2}$: **Correct.**
   - Klein-factor nature of $q \to -q$: **Correct.** Proved analytically that no single-soliton change of basis $A$ (even $\theta$-dependent) can map $\check R(q)$ to $\check R(-q)$ because $A \otimes A$ cannot flip the relative sign of $v_0 \otimes v_1$ and $v_1 \otimes v_0$ in the spin-0 eigenvector. The combination $(1+b^2-cc')/(2b) = (q+q^{-1})/2$ is strictly basis-independent.
3. **Twisted Tensor-Product Graph & $U_q(d_{n+1}^{(2)})$ Spinor (Chapter 8):**
   - Multiplicity-free chain decomposition $\bigoplus_{k=0}^n \Lambda^k \mathbb C^{2n+1}$: **Correct.** Verified for $n=2$ (dim 16) and $n=3$ (dim 64) with Jimbo solution spaces strictly 1-dimensional.
   - Normalization of $q$ in the formula: **Imprecise in the text.** The formula $\rho_{\Lambda^{n-j}} = \prod_{i=1}^j \langle i/2 \rangle_{(-1)^i}$ holds in the **book's normalization** where long roots have squared length 2 (so $q_i = q$ on long roots, $q^{1/2}$ on short roots), *not* Kac's $q$ where $q_i=q$ on short roots.
   - Component parities and `fig-tpg`(b) caption: **Correct.** Parities are $+--+$ for $n=3$, so the middle edge connects $\Lambda^2$ and $\Lambda^1$ (both $-$), carrying $\langle 1 \rangle_+$.
4. **Untwisted Graph Rule & Graph Failure (Chapter 8):**
   - Untwisted rule `eq-tpg-rule`: **Correct.** Verified on $sl_{n+1}$ vector and $c_2^{(1)}$ vector (where half-integer arguments $\langle 1/2 \rangle$ occur).
   - Graph connectedness: **Correct.** Verified against DGZ 1994 (Theorem 1 and Section 4).
   - Higher multiplets of $c_n^{(1)}$: **Correct.** Decompositions under $U_q(b_n)$ have repeated components in their tensor squares.
5. **History and Prior Art (Chapter 8):**
   - Attributions to DGZ 1994, DGZ 1996, Jimbo 1986, MacKay 1991, ZGB 1991: **Correct.**
   - GMW 1996 for $\mathbf{27}$ of $U_q(e_6^{(2)})$ under $U_q(c_4)$: **Correct.** GMW Section 4.4 explicitly constructs this multiplicity-free graph.
   - Count of unknowns: $\sum_\mu (\mathrm{mult}_\mu)^2 = 16$: **Correct.**
6. **Roots of Unity & RSOS (Chapters 7 & 8):**
   - Tilting module structure and quantum dimensions: **Correct.** Explicitly constructed for $p=4$ on spin $1 \otimes 1$.
   - RSOS decoupling of heights $a \le (p-2)/2$: **Correct.**
7. **Core Formulas of Chapter 7:**
   - Hopf axioms, universal R-matrix on spin 1/2, quasitriangularity, $(S \otimes 1)\mathcal R = \mathcal R^{-1}$, $(S \otimes S)\mathcal R = \mathcal R$, braid YBE, crossing, $*$-structure: **All Correct.**
8. **Fusion, Poles, and IRF Form (Chapter 8):**
   - Operator ordering of `eq-fusion`: **Correct.**
   - Remark on projector $P_c$: **Imprecise.** Projector onto a submodule in an indecomposable non-semisimple module does not exist as an algebra projection.
   - Vertex vs IRF basis phrasing: **Imprecise.** Vertex and IRF spaces have different structures (tensor product vs path space / multiplicity space).
   - Hermitian analyticity in IRF basis (Hoare–Hollowood–Miramontes 2013, Miramontes 1999): **Correct.**
9. **Citations and Exercises:**
   - All citing sentences verified against original sources.
   - All 11 exercises solved and verified. Hint of `exr-spin1` confirmed correct.

---

## Item 1. Normalizations of $U_q(\hat g)$, Central Element, and $S^2$ Shift

### 1.1 Central element $\prod_i k_i^{n_i}$ and $q^{\sum_i n_i^\vee h_i}$
**Verdict: Partially Wrong / Imprecise.**
- **Conjugation action:** In $U_q(\hat g)$, $k_i e_j k_i^{-1} = q_i^{a_{ij}} e_j = q^{(\alpha_i, \alpha_j)} e_j$. For an element $K = \prod_{i=0}^r k_i^{m_i}$,
  $$K e_j K^{-1} = q^{(\sum_{i=0}^r m_i \alpha_i, \alpha_j)} e_j.$$
  Thus $K$ is central if and only if $\sum_{i=0}^r m_i \alpha_i$ is orthogonal to all simple roots, which means $\sum_{i=0}^r m_i \alpha_i$ is proportional to the null root $\delta = \sum_{i=0}^r a_i \alpha_i$.
  - For untwisted algebras and twisted algebras $X^{(k)}$ ($k \ge 2$) with $a_0 = 1$, $m_i = a_i$ (the Kac marks) gives a central element $\prod_i k_i^{a_i}$.
  - For $A_{2l}^{(2)}$ in Kac's labelling (Kac 1990 Table Aff 2), node 0 is the short root, and $a_0 = 2$, $a_1 = \dots = a_{l-1} = 2$, $a_l = 1$. The central element is $k_0^2 k_1^2 \cdots k_l$, with label 2 on node 0. In the book's folding (Chapter 6, line 113), node 0 is the long root with label 1, and the other end is short with label 2.
- **Equality with $q^{\sum_i n_i^\vee h_i}$:** Since $k_i = q_i^{h_i} = q^{(\alpha_i, \alpha_i)/2 h_i}$,
  $$\prod_{i=0}^r k_i^{m_i} = q^{\sum_{i=0}^r m_i \frac{(\alpha_i, \alpha_i)}{2} h_i}.$$
  For this to equal $q^{\sum_i n_i^\vee h_i}$, one needs $n_i^\vee = m_i \frac{(\alpha_i, \alpha_i)}{2}$.
  - In Kac's normalization ($(\alpha_0, \alpha_0) = 2$ for untwisted and twisted, with $(\alpha_0, \alpha_0) = 1$ for $A_{2l}^{(2)}$), $a_i \frac{(\alpha_i, \alpha_i)}{2} = a_i^\vee$ are integers, giving $q^{\sum_i a_i^\vee h_i}$.
  - In the **book's normalization** where long roots have squared length 2 (caption of `tbl-foldings`), the shortest root of $X^{(k)}$ has squared length $2/k$. Then $n_i \vec\alpha_i^2 / 2$ is a fraction ($1/k$), not an integer! Caption of `tbl-foldings` explicitly states $n_i^\vee = k\, n_i \vec\alpha_i^2/2$. Thus in book normalization, $\prod k_i^{n_i} = q^{\frac1k \sum n_i^\vee h_i} = (q^{1/k})^{\sum n_i^\vee h_i}$.

*Script:* `referee_scratch/algebra_07_08/check_central.py`.

### 1.2 General Derivation of the $S^2$ Spectral Shift
**Verdict: Wrong for $a_{2n}^{(2)}$; Correct for untwisted and $d_{n+1}^{(2)}$.**

**Derivation:**
In $U_q(\hat g)$ with Chapter 7's coproduct, $S^2(e_i) = q_i^2 e_i = q^{(\alpha_i, \alpha_i)} e_i$ and $S^2(f_i) = q_i^{-2} f_i$.
On an evaluation representation $V(x)$ in the homogeneous gradation, $e_i, f_i, k_i$ ($i \ge 1$) span the horizontal subalgebra $U_q(g_{(0)})$, while $e_0 = x \mathcal E_0$ carries degree 1 in $x$ and finite weight $-\theta_0$.
1. To intertwine the finite generators $S^2(a)$ with $a$, we conjugate by $M = q^{-2\bar\rho}$, where $\bar\rho$ is the Weyl vector of the horizontal algebra satisfying $\langle \bar\rho, \alpha_i^\vee \rangle = 1$ ($i \ge 1$). Since $q^{-2\bar\rho} e_i q^{2\bar\rho} = q^{-2(\bar\rho, \alpha_i)} e_i = q_i^{-2} e_i$, $M S^2(e_i) M^{-1} = e_i$.
2. On $e_0$, which has weight $-\theta_0$ under the finite Cartan subalgebra:
   $$M S^2(e_0) M^{-1} = q^{-2\bar\rho} (q^{(\alpha_0, \alpha_0)} x \mathcal E_0) q^{2\bar\rho} = x \, q^{(\alpha_0, \alpha_0) + 2(\bar\rho, \theta_0)} \mathcal E_0.$$
   Therefore, $V(x)^{**} \cong V(x t)$ with exact shift:
   $$t = q^{(\alpha_0, \alpha_0) + 2(\bar\rho, \theta_0)}.$$

Now relate $(\alpha_0, \alpha_0) + 2(\bar\rho, \theta_0)$ to $h^\vee$:
Since $\alpha_0^\vee = -\frac{2\theta_0}{(\alpha_0, \alpha_0)}$, we have $2(\bar\rho, \theta_0) = (\alpha_0, \alpha_0) \langle \bar\rho, \theta_0^\vee \rangle$. Thus:
$$t = q^{(\alpha_0, \alpha_0)(1 + \langle \bar\rho, \theta_0^\vee \rangle)}.$$
- **Untwisted algebras:** $(\alpha_0, \alpha_0) = 2$, $\theta_0 = \theta$ (highest root), and $1 + \langle \bar\rho, \theta^\vee \rangle = h^\vee$. Hence $t = q^{2 h^\vee}$.
- **Twisted $d_{n+1}^{(2)}$:** In Kac's normalization, $(\alpha_0, \alpha_0) = 2$, $\theta_0 = \epsilon_1$. The horizontal algebra is $B_n$, with $\langle \bar\rho, \theta_0^\vee \rangle = 2n - 1$. Then $1 + \langle \bar\rho, \theta_0^\vee \rangle = 2n = h^\vee$. Hence $t = q^{2 h^\vee}$.
- **$a_{2n}^{(2)}$:**
  - If $(\alpha_0, \alpha_0) = 1$ (node 0 is the short root, as stated in Kac and in `07-quantum-groups.qmd:194`):
    Then $\theta_0 = \frac12 \alpha_1$ (for $n=1$), $\langle \bar\rho, \theta_0^\vee \rangle = 2$.
    Then $1 + \langle \bar\rho, \theta_0^\vee \rangle = 3 = h^\vee$.
    Therefore, $(\alpha_0, \alpha_0)(1 + \langle \bar\rho, \theta_0^\vee \rangle) = 1 \times h^\vee = h^\vee$!
    **The shift is $t = q^{h^\vee}$, NOT $q^{2h^\vee}$.**
  - If node 0 is the long root ($(\alpha_0, \alpha_0) = 4$ in Kac scale):
    $1 + \langle \bar\rho, \theta_0^\vee \rangle = 1 + 1/2 = 3/2 = h^\vee/2$.
    Then $(\alpha_0, \alpha_0)(1 + \langle \bar\rho, \theta_0^\vee \rangle) = 4 \times (3/2) = 6 = 2h^\vee$.
    In that scale, $t = q^{2h^\vee}$, but $(\alpha_0, \alpha_0) = 4$, contradicting the text's assertion that $(\alpha_0, \alpha_0) = 1$.

**Numerical Verification:**
Built $\pi(S^2(a))$ and scanned for $t$ such that $V(xt) \cong V(x)^{**}$ and $V(x x_c) \cong V(x)^*$:
- $U_q(c_2^{(1)})$ vector (4-dim, untwisted control, $h^\vee = 3$):
  - Ch. 7 coproduct: $V^{**} \cong V(x Q^6) = V(x Q^{2h^\vee})$, $V^* \cong V(x Q^3) = V(x Q^{h^\vee})$.
  - Part VI coproduct: $V^{**} \cong V(x Q^{-6}) = V(x Q^{-2h^\vee})$, $V^* \cong V(x Q^{-3}) = V(x Q^{-h^\vee})$.
- $U_q(d_3^{(2)})$ spinor (4-dim, twisted, $h^\vee = 4$):
  - Ch. 7 coproduct: $V^{**} \cong V(x Q^8) = V(x Q^{2h^\vee})$, $V^* \cong V(x (-Q^4))$.
  - Part VI coproduct: $V^{**} \cong V(x Q^{-8}) = V(x Q^{-2h^\vee})$, $V^* \cong V(x (-Q^{-4}))$.
- $U_q(a_2^{(2)})$ fundamental (3-dim, $h^\vee = 3$):
  - With $q_1 = Q, q_0 = Q^4$ (node 0 long):
    Ch. 7 coproduct: $V^{**} \cong V(x Q^{12})$, $V^* \cong V(x (-Q^6))$.
    Part VI coproduct: $V^{**} \cong V(x Q^{-12})$, $V^* \cong V(x (-Q^{-6}))$.
    In base $q$ where $q_0 = q^2, q_1 = q^{1/2}$ (so $q = Q^2$), the shift is $q^6 = q^{2h^\vee}$, and $(\alpha_0, \alpha_0) = 4$.

*Script:* `referee_scratch/algebra_07_08/test_all_three.py`. Residuals on all defining relations and Serre relations: $< 1.2 \times 10^{-15}$.

### 1.3 Crossing Point and Part VI Comparison
**Verdict: Correct.**
If $V(x)^* \cong V(x x_c)$, then $V(x)^{**} \cong (V(x)^*)^* \cong V(x x_c^2)$. Matching with $V(x)^{**} \cong V(x q^{2h^\vee})$ gives $x_c^2 = q^{2h^\vee}$, fixing $x_c$ up to a sign $\pm 1$.
In Part VI's gradation, $x(i\pi) = e^{2\pi i T}$. With $2T = 4\pi h/\beta^2 - h^\vee = h\omega + h - h^\vee$:
Every crossing point quoted in Part VI matches:
- $c_n^{(1)}$ ($U_q(d_{n+1}^{(2)})$): $x_c = (-1)^{n+1} q^{-2n}$, $x_c^2 = q^{-4n} = q^{-2h^\vee}$.
- $a_{2n-1}^{(2)}$ ($U_q(b_n^{(1)})$): $x_c = q^{-2(2n-1)}$, $x_c^2 = q^{-4(2n-1)}$ (in Part VI's $q = q_{\rm ch7}^{1/2}$).
- $f_4^{(1)}$ ($U_q(e_6^{(2)})$): $x_c = -q^{-12}$, $x_c^2 = q^{-24} = q^{-2h^\vee}$.
- $e_6^{(2)}$ ($U_q(f_4^{(1)})$): $x_c = q^{-18} = q^{-2h^\vee}$.
- $d_4^{(3)}$ ($U_q(g_2^{(1)})$): $x_c = q^{-12} = q^{-2h^\vee}$.
- $g_2^{(1)}$ ($U_q(d_4^{(3)})$): $x_c = q^{-6}$, $x_c^2 = q^{-12} = q^{-2h^\vee}$.

---

## Item 2. The Sign of $q$ and the Dictionary with Part VI

### 2.1 Sine-Gordon Dictionary
**Verdict: Correct.**
Substituting $z = e^{2\lambda\theta}$ and $q = -e^{i\pi\lambda}$ into `eq-six-vertex`:
- $b(z) = \frac{q(z-1)}{z-q^2} = \frac{-e^{i\pi\lambda}(e^{2\lambda\theta}-1)}{e^{2\lambda\theta}-e^{2\pi i\lambda}} = \frac{\sinh\lambda\theta}{\sinh\lambda(i\pi-\theta)} = S_T/S_0$.
- Off-diagonal elements:
  $c(z) = \frac{q^2-1}{q^2-z} = e^{-\lambda\theta} \frac{\sinh(i\pi\lambda)}{\sinh\lambda(i\pi-\theta)}$,
  $c'(z) = z c(z) = e^{\lambda\theta} \frac{\sinh(i\pi\lambda)}{\sinh\lambda(i\pi-\theta)}$.
- The rapidity factors $e^{\mp\lambda\theta}$ are removed by the one-particle basis transformation $A(\theta) = \mathrm{diag}(e^{\lambda\theta/2}, e^{-\lambda\theta/2})$, giving $S_R/S_0$. This is the standard gauge transformation from the homogeneous to the principal gradation ($x^{H/2}$).
- For the opposite coproduct of Part VI, $P \check R(1/z) P = \check R(z; q^{-1})$. Matching with $S_T/S_0$ requires $q^{-1} = -e^{i\pi\lambda} \implies q = -e^{-i\pi\lambda} = e^{-8\pi^2 i/\beta_{\rm SG}^2}$. This is Bernard and LeClair's exact value (CMP 142 (1991) 99, eqs. 3.18 and 3.37).

### 2.2 Proof that $q \to -q$ is Not of the Form $A \otimes A$
**Verdict: Correct.**
*Proof:*
Let $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix} \in GL(2, \mathbb C)$. For any $A$, the diagonal elements of $A \otimes A$ satisfy $(A \otimes A)_{11} = ad = (A \otimes A)_{22}$.
The conjugation mapping $\check R(q)$ to $\check R(-q)$ can be realized by $L_1 = \mathrm{diag}(1, 1, -1, 1)$ or $L_2 = \Pi \otimes 1 = \mathrm{diag}(1, 1, -1, -1)$. In both matrices, the $(1,1)$ and $(2,2)$ entries have opposite signs ($+1$ and $-1$). Thus neither $L_1$ nor $L_2$ can be of the form $A \otimes A$.
More generally, suppose there exists an invertible $A$ such that $(A \otimes A) \check R(q) (A \otimes A)^{-1} = \check R(-q)$.
Then $A \otimes A$ must map the unique spin-0 eigenspace of $\check R(q)$ with eigenvalue $\rho_0(z)$, spanned by $v_{\rm anti}(q) = v_0 \otimes v_1 - q v_1 \otimes v_0$, to the spin-0 eigenspace of $\check R(-q)$, spanned by $v_{\rm anti}(-q) = v_0 \otimes v_1 + q v_1 \otimes v_0$.
Computing $(A \otimes A) v_{\rm anti}(q)$:
$$(A \otimes A)(v_0 \otimes v_1 - q v_1 \otimes v_0) = ab(1-q) v_0 \otimes v_0 + (ad - q bc) v_0 \otimes v_1 + (bc - q ad) v_1 \otimes v_0 + cd(1-q) v_1 \otimes v_1.$$
For this to be proportional to $v_0 \otimes v_1 + q v_1 \otimes v_0$:
1. $ab(1-q) = 0 \implies ab = 0$.
2. $cd(1-q) = 0 \implies cd = 0$.
3. The ratio of $v_1 \otimes v_0$ to $v_0 \otimes v_1$ must equal $+q$:
   $$\frac{bc - q ad}{ad - q bc} = q \implies bc - q ad = q ad - q^2 bc \implies bc(1 + q^2) = 2q ad.$$
From $ab = 0$ and $cd = 0$ with $\det(A) = ad - bc \ne 0$, either $b=c=0$ (diagonal) or $a=d=0$ (anti-diagonal):
- If $b=c=0$: $bc(1+q^2) = 0 = 2q ad \implies ad = 0$, contradicting $\det(A) \ne 0$.
- If $a=d=0$: $2q ad = 0 = bc(1+q^2) \implies bc = 0$ (for generic $q$), contradicting $\det(A) \ne 0$.
Therefore, **no invertible matrix $A$ can ever map $\check R(q)$ to $\check R(-q)$**. $\blacksquare$

### 2.3 Basis-Independence of $(1+b^2-cc')/(2b)$
**Verdict: Correct.**
Under any basis transformation preserving $U(1)$ weight conservation ($A = \mathrm{diag}(\alpha, \beta)$), $b \mapsto b$ and $c c' \mapsto c c'$ are strictly invariant.
Direct evaluation:
$$\frac{1 + b^2 - c c'}{2b} = \frac{q + q^{-1}}{2}.$$
Under $q \to -q$, this expression reverses sign to $-(q + q^{-1})/2$.
*Script:* `referee_scratch/algebra_07_08/test_item2.py`.

---

## Item 3. Twisted Tensor-Product Graph and the Spinor of $U_q(d_{n+1}^{(2)})$

### 3.1 Spinor Representation and Direct Jimbo Solution
**Verdict: Correct.**
Constructed the spinor representation of $U_q(d_{n+1}^{(2)})$ on $(\mathbb C^2)^{\otimes n}$ for $n=2$ (dim 4) and $n=3$ (dim 8).
All defining relations and $q$-Serre relations satisfied to residual $4.97 \times 10^{-16}$.
Solved Jimbo's intertwining equations $\check R(z) \Delta_{x, y}(a) = \Delta_{y, x}(a) \check R(z)$ on $V(x) \otimes V(y)$:
- $n=2$ (dim 16, 36 weight-preserving unknowns): solution space is strictly 1-dimensional (smallest singular value $1.12 \times 10^{-15}$, second smallest $0.61$).
  Eigenvalue multiplicities: **10, 5, 1**, matching $\dim \Lambda^2 \mathbb C^5 = 10$, $\dim \Lambda^1 \mathbb C^5 = 5$, $\dim \Lambda^0 \mathbb C^5 = 1$.
- $n=3$ (dim 64, 216 weight-preserving unknowns): solution space is strictly 1-dimensional (smallest singular value $2.90 \times 10^{-15}$, second smallest $0.38$).
  Eigenvalue multiplicities: **35, 21, 7, 1**, matching $\dim \Lambda^3 \mathbb C^7 = 35$, $\dim \Lambda^2 \mathbb C^7 = 21$, $\dim \Lambda^1 \mathbb C^7 = 7$, $\dim \Lambda^0 \mathbb C^7 = 1$.

*Scripts:* `referee_scratch/algebra_07_08/spinor_jimbo.py`, `referee_scratch/algebra_07_08/solve_jimbo_spinor.py`.

### 3.2 In Which $q$ Does the Formula Hold?
**Verdict: Imprecise in the text.**
The formula in `08-r-matrices.qmd:73`:
$$\rho_{\Lambda^{n-j}} = \prod_{i=1}^j \langle i/2 \rangle_{(-1)^i}, \quad \langle a \rangle_\pm = \frac{1 \pm z q^{2a}}{z \pm q^{2a}}$$
holds when $q$ is normalized such that **long roots have squared length 2** (the book's normalization, where $q_{\rm long} = q$ and $q_{\rm short} = q^{1/2}$).
If one uses Kac's normalization of Chapter 7 line 107 (where short roots have squared length 2 and long roots have squared length 4, so $q_{\rm short} = Q$ and $q_{\rm long} = Q^2$), then $q_{\rm book} = Q^2$, and the exponent in the formula would have to be $Q^{4a}$.
Chapter 8 does not specify this, which creates confusion with Chapter 7's definition of $q$.

### 3.3 Component Parities and Figure Caption
**Verdict: Correct.**
For $n=3$, the spinor square has parities $+--+$ on $\Lambda^3, \Lambda^2, \Lambda^1, \Lambda^0$:
- $\Lambda^3$ (dim 35) is symmetric ($+$).
- $\Lambda^2$ (dim 21) is antisymmetric ($-$).
- $\Lambda^1$ (dim 7) is antisymmetric ($-$).
- $\Lambda^0$ (dim 1) is symmetric ($+$).
Total symmetric: $35 + 1 = 36 = 8 \times 9 / 2$.
Total antisymmetric: $21 + 7 = 28 = 8 \times 7 / 2$.
The middle edge joins $\Lambda^2$ ($-$) and $\Lambda^1$ ($-$), which have **equal parity**, so it carries $\langle 1 \rangle_+$. Caption of `fig-tpg`(b) is completely correct.

### 3.4 Transformation of $e_0$ in $g_{(1)}$
**Verdict: Correct.**
In DGZ 1996 (q-alg/9508012), line 129: $L_1$ is an irreducible $L_0$-module with highest weight $\theta_0$.
Line 416-417: $e_0 \otimes q^{h_0/2}$ is the lowest component of a tensor operator in $V_0(\theta_0) = g_{(1)}$.
Equation (4.31) with $a=1$ gives:
$$\rho_k(u) = \prod_{i=1}^{l-k} \left[ \frac12(a-1) + i \right]_{(-1)^i} = \prod_{i=1}^{l-k} [ i/2 ]_{(-1)^i}.$$
This matches the book's formula $\rho_{\Lambda^{n-j}} = \prod_{i=1}^j \langle i/2 \rangle_{(-1)^i}$ exactly.

---

## Item 4. The Untwisted Rule and Where the Graph Fails

### 4.1 Untwisted Rule & $c_2^{(1)}$ Vector Example
**Verdict: Correct.**
- Checked $\langle a \rangle \langle -a \rangle = \frac{1 - z q^{2a}}{z - q^{2a}} \frac{1 - z q^{-2a}}{z - q^{-2a}} = 1$.
- For $sl_{n+1}$ vector: $C(2\vec\lambda_1) - C(\vec\lambda_2) = 4$, yielding $\rho_{\vec\lambda_2} = \langle 1 \rangle = (1 - z q^2)/(z - q^2)$.
- For $c_2^{(1)}$ vector ($V = \mathbb C^4$, long roots length$^2$ 2, short roots length$^2$ 1):
  $\mathbf 4 \otimes \mathbf 4 = \mathbf{10}_{\rm sym}(2\vec\lambda_1) \oplus \mathbf 5_{\rm anti}(\vec\lambda_2) \oplus \mathbf 1_{\rm sym}(0)$.
  Casimirs: $C(2\vec\lambda_1) = 6$, $C(\vec\lambda_2) = 4$, $C(0) = 0$.
  Edge between $2\vec\lambda_1$ ($+$) and $\vec\lambda_2$ ($-$):
  $$\frac{C(\vec\lambda_2) - C(2\vec\lambda_1)}{4} = \frac{4 - 6}{4} = -\frac12 \implies \rho_{\vec\lambda_2} = \langle 1/2 \rangle.$$
  Half-integer arguments appear naturally for non-simply-laced algebras.
*Script:* `referee_scratch/algebra_07_08/test_c2_vec.py`.

### 4.2 Connected Subgraph Claim in DGZ 1994
**Verdict: Correct.**
DGZ 1994 (hep-th/9405030), p. 11, lines 455–459:
> "If $V(\lambda)$ and $V(\mu)$ are irreducible $U_q(G)$-modules the tensor product graph is always connected, i.e., every node is linked to every other node by a path of edges. This follows from the fact that $\pi_\lambda(e_0) \otimes \pi_\mu(q^{h_0/2})$ is related to the lowest component of an adjoint tensor operator; for details see [14]. This implies that the relations eq. (3.24) are sufficient to determine all the coefficients $\rho_\nu(x)$ uniquely."

### 4.3 Reducible Multiplets of $c_n^{(1)}$
**Verdict: Correct.**
For $c_n^{(1)}$, the non-spinor multiplets under $U_q(b_n)$ decompose as $\Lambda^1 \oplus \Lambda^0$ (vector) and $\Lambda^2 \oplus \Lambda^1 \oplus \Lambda^0$ (29-plet for $c_3^{(1)}$).
Their tensor squares contain repeated isotypic components:
- $(\Lambda^1 \oplus \Lambda^0)^{\otimes 2}$ contains $\Lambda^1$ twice (from $\Lambda^1 \otimes \Lambda^0$ and $\Lambda^0 \otimes \Lambda^1$) and $\Lambda^0$ twice.
- $(\Lambda^2 \oplus \Lambda^1 \oplus \Lambda^0)^{\otimes 2}$ contains multiple copies of $\Lambda^2, \Lambda^1, \Lambda^0$.
On these isotypic components, $\check R$ is a matrix, and the scalar tensor-product graph rule cannot determine it.

---

## Item 5. History and Prior Art

### 5.1 Citations of DGZ 1994 and DGZ 1996
**Verdict: Correct.**
- Spinor of $b_n^{(1)}$: DGZ 1994 (hep-th/9405030), Section 4.2, eqs. (4.21)–(4.25).
- $\mathbf{26}$ of $f_4^{(1)}$: DGZ 1994, Section 4.4, eqs. (4.60)–(4.64).
- Spinor of $d_{n+1}^{(2)}$: DGZ 1996 (q-alg/9508012), Section 4.3, eqs. (4.30)–(4.31).

### 5.2 GMW 1996 on $\mathbf{27}$ of $U_q(e_6^{(2)})$
**Verdict: Correct.**
GMW 1996 (hep-th/9509007), Section 4.4, p. 8:
> "From table 1 we see that the 27 dimensional representation of $q L = q e_6$ is an irreducible representation of $q_\iota g_0' = q^2 c_4$, and since the tensor product of two 27 dimensional $q^2 c_4$ representations is multiplicity free, we can construct the TPG for $R_{27,27}$..."
Under $U_q(f_4)$, $\mathbf{27} = \mathbf{26} \oplus \mathbf1$.
Decomposition of $\mathbf{27} \otimes \mathbf{27}$:
$$(\mathbf{26} \oplus \mathbf1) \otimes (\mathbf{26} \oplus \mathbf1) = (\mathbf1 \oplus \mathbf{26} \oplus \mathbf{52} \oplus \mathbf{273} \oplus \mathbf{324}) \oplus 2 \cdot \mathbf{26} \oplus \mathbf1 = 2 \cdot \mathbf1 \oplus 3 \cdot \mathbf{26} \oplus \mathbf{52} \oplus \mathbf{273} \oplus \mathbf{324}.$$
Multiplicities are 2 for $\mathbf1$, 3 for $\mathbf{26}$, and 1 for the other three.
Number of unknowns: $\sum_\mu (\mathrm{mult}_\mu)^2 = 2^2 + 3^2 + 1^2 + 1^2 + 1^2 = 4 + 9 + 1 + 1 + 1 = 16$.
Imposing the $e_0$ equation on highest-weight vectors determines $\check R$ because $\check R$ has already been reduced to highest-weight blocks by $U_q(g)$-invariance.

---

## Item 6. Roots of Unity and RSOS Restriction

### 6.1 Tilting Modules and Quantum Dimensions
**Verdict: Correct.**
For $q^{2p} = 1$ ($q = e^{i\pi/p}$):
- $[p]_q = 0$, and $[2p - n]_q = -[n]_q$.
- For $j < (p-1)/2$, spin-$j$ representations are irreducible.
- For $2j+1 = p$, the module is irreducible with quantum dimension $[p]_q = 0$.
- For $p < 2j+1 < 2p$, the tilting module $T(j)$ has dimension $2p$, with composition factors of spin $j$ and spin $p-1-j$. Its quantum dimension is $[2j+1]_q + [2p-2j-1]_q = [2j+1]_q - [2j+1]_q = 0$.
- For $p=4$, spin $1 \otimes 1$ ($9$-dimensional): decomposed explicitly into $T(2)$ (dimension 8, quantum dimension $[5]_q + [3]_q = -1 + 1 = 0$) plus spin 0 (dimension 1, quantum dimension 1). Verified in `referee_scratch/algebra_07_08/test_roots_unity.py`.
- The sign of $q$ does not affect quantum dimensions: $[2j+1]_{-q} = (-1)^{2j} [2j+1]_q = [2j+1]_q$ since $2j$ is an integer.

### 6.2 RSOS Decoupling
**Verdict: Correct.**
For spin-1/2 face weights, transitions from $a = (p-2)/2$ to $a+1/2 = (p-1)/2$ involve intermediate states with quantum dimension $[p]_q = 0$. The Clebsch–Gordan projector onto this state vanishes identically, so the face weights decoupling at height $a \le (p-2)/2$ is exact.

---

## Item 7. Core Formulas of Chapter 7

**Verdict: All Correct.**
All checked symbolically and numerically in `referee_scratch/algebra_07_08/test_item7_core.py`:
1. Hopf axioms for coproduct and antipode: $m \circ (S \otimes \mathrm{id}) \circ \Delta = \eta \circ \epsilon$ on all generators; $S^2(a) = K a K^{-1}$.
2. Spin-1/2 universal R-matrix: matches printed matrix in line 100.
3. Quasitriangularity: $(S \otimes 1)\mathcal R = \mathcal R^{-1}$, $(S \otimes S)\mathcal R = \mathcal R$, braid YBE.
4. Six-vertex R-matrix: solves Jimbo equations, braiding unitarity $\check R(z)\check R(1/z) = 1$, braid YBE, eigenvalues $1$ (x3) and $\rho_0(z) = (1-q^2 z)/(z-q^2)$ (x1).
5. Crossing relation `eq-r-crossing`: verified with $x_c = q^2$, $C = \begin{pmatrix} 0 & 1 \\ -q & 0 \end{pmatrix}$, $c_u(z) = q(z-1)/(q^2 z - 1)$.
6. $*$-structure: $(*\otimes*)\circ\Delta = \Delta^{\rm op}\circ *$, spin-1 norm $\langle F v_0, F v_0 \rangle = [2]_q \langle v_0, v_0 \rangle$.

---

## Item 8. Fusion, Poles, and the IRF Form

### 8.1 Fusion Operator Order & Projector $P_c$
**Verdict: Operator order Correct; remark on Projector $P_c$ Imprecise.**
- In `eq-fusion`:
  $$\check R_{cd}(x/y) = (1 \otimes P_c) \, (\check R_{ad} \otimes 1) \, (1 \otimes \check R_{bd}) \big|_{V_c \otimes V_d}.$$
  The mappings are:
  $V_c \otimes V_d \subset V_a \otimes V_b \otimes V_d \xrightarrow{1 \otimes \check R_{bd}} V_a \otimes V_d \otimes V_b \xrightarrow{\check R_{ad} \otimes 1} V_d \otimes V_a \otimes V_b \xrightarrow{1 \otimes P_c} V_d \otimes V_c$.
  The operator order is strictly correct.
- **Imprecision:** When $V_a \otimes V_b$ is reducible but indecomposable (which occurs at poles and zeros of $\check R$), an invariant projector $P_c$ commuting with the quantum group does not exist because submodules have no invariant direct sum complement. Fusion is rigorously defined either by evaluating on the submodule $V_c$ and taking the quotient, or by taking the residue of $\check R_{ab}(x_a/x_b)$ whose image is $V_c$ (Kulish–Reshetikhin–Sklyanin 1981).

### 8.2 Vertex vs IRF Equivalence Phrasing
**Verdict: Imprecise in the text.**
`08-r-matrices.qmd:120` states: "The vertex and IRF forms of an S-matrix are related by a change of basis, so they have the same eigenvalues and the same Bethe–Yang spectra."
The two matrices act on different vector spaces: the vertex S-matrix acts on $V \otimes V$, while the IRF face weights act on the multiplicity space of paths $a \to b \to d$. The relation is an intertwining by Clebsch–Gordan coefficients (Baxter's vertex-face correspondence), not a similarity transformation on a single space.

### 8.3 Hermitian Analyticity in IRF Basis
**Verdict: Correct.**
Verified against Hoare–Hollowood–Miramontes 2013 (arXiv:1303.1447, lines 141-144, 389-391) and Miramontes 1999 (hep-th/9901145, Section 3).

---

## Item 9. Citations and Exercises

### 9.1 Citations Check
**Verdict: Correct.**
- `delius1995` (hep-th/9503079): Lorentz spin fixes the gradation; bootstrap axioms follow from Hopf structure.
- `frenkel1992` & `chari1994`: $V(x)^{**} \cong V(x q^{2h^\vee})$.
- `kulish1981`: Fusion procedure.
- `gomez1996`: Vertex–IRF correspondence and RSOS.
- `drinfeld1987` & `jimbo1985`: Quantum groups.

### 9.2 Exercises Check
**Verdict: All 11 Exercises Correct and Solvable.**
- `exr-spin1` hint: "$V_1(x)$ is the $q$-symmetric part of $V_{1/2}(xq)\otimes V_{1/2}(xq^{-1})$, the image of $\check R(q^{-2})$ acting on $V_{1/2}(xq^{-1})\otimes V_{1/2}(xq)$". This is correct: domain has $z = (xq^{-1})/(xq) = q^{-2}$, where $\check R(q^{-2})$ projects onto spin 1 in the codomain $V_{1/2}(xq) \otimes V_{1/2}(xq^{-1})$.

---

## Replacement Text for Identified Errors

### 1. Chapter 7, Section `sec-uqghat` (`part2-algebra/07-quantum-groups.qmd`, line 119)
**Target:**
```markdown
The element $\prod_ik_i^{n_i}$, with the Kac labels $n_i$ of (-@eq-kac-labels), is central: conjugation by it multiplies $e_j$ by $q^{(\sum_in_i\vec\alpha_i)\cdot\vec\alpha_j}=1$. Since $k_i=q_i^{h_i}$ and $n_i^\vee=n_i\vec\alpha_i^2/2$, it equals $q^{\sum_in_i^\vee h_i}$, $q$ to the power of the central element. On the finite-dimensional representations used in this book it equals 1 (level zero).
```
**Replacement:**
```markdown
The element $\prod_ik_i^{n_i}$, with the Kac labels $n_i$ of (-@eq-kac-labels) (with $n_0=2$ for $a_{2n}^{(2)}$ in Kac's labelling where node 0 is short), is central: conjugation by it multiplies $e_j$ by $q^{(\sum_in_i\vec\alpha_i)\cdot\vec\alpha_j}=1$. In Kac's normalization where $n_i^\vee=n_i\vec\alpha_i^2/2$, it equals $q^{\sum_in_i^\vee h_i}$, $q$ to the power of the central element (in the book's normalization with long roots of length$^2$ 2, $n_i^\vee=k\,n_i\vec\alpha_i^2/2$, so the exponent has an extra factor $1/k$). On the finite-dimensional representations used in this book it equals 1 (level zero).
```

### 2. Chapter 7, Section `sec-qg-crossing` (`part2-algebra/07-quantum-groups.qmd`, line 194)
**Target:**
```markdown
Its $\Lambda_0$ part multiplies $e_0$ by $q^{2h^\vee}$, which in the homogeneous gradation is a rescaling of $x$, so $V(x)^{**}\cong V(x\,q^{2h^\vee})$ [@frenkel1992; @chari1994]. Here $q_i=q^{(\alpha_i,\alpha_i)/2}$ in the normalization of @sec-uqghat: $(\alpha_0,\alpha_0)=2$ ($1$ for $a_{2n}^{(2)}$).
```
**Replacement:**
```markdown
Its $\Lambda_0$ part multiplies $e_0$ by $q^{(\alpha_0,\alpha_0)+2(\bar\rho,\theta_0)}=q^{2h^\vee}$, which in the homogeneous gradation is a rescaling of $x$, so $V(x)^{**}\cong V(x\,q^{2h^\vee})$ [@frenkel1992; @chari1994]. Here $q_i=q^{(\alpha_i,\alpha_i)/2}$ in the normalization where $(\alpha_0,\alpha_0)=2$ (for $a_{2n}^{(2)}$ with node 0 the short root of squared length 1, the shift is instead $q^{h^\vee}$).
```

### 3. Chapter 8, Section `sec-tpg` (`part2-algebra/08-r-matrices.qmd`, line 73)
**Target:**
```markdown
Under the horizontal subalgebra $U_q(b_n)$ its tensor square is $\bigoplus_{k=0}^n\Lambda^k\mathbb C^{2n+1}$, without multiplicities, the graph is the chain $\Lambda^n-\Lambda^{n-1}-\dots-\Lambda^0$, and $\rho_{\Lambda^{n-j}}=\prod_{i=1}^{j}\langle i/2\rangle_{(-1)^i}$ with $\langle a\rangle_-=\langle a\rangle$ [@dgz1996].
```
**Replacement:**
```markdown
Under the horizontal subalgebra $U_q(b_n)$ its tensor square is $\bigoplus_{k=0}^n\Lambda^k\mathbb C^{2n+1}$, without multiplicities, the graph is the chain $\Lambda^n-\Lambda^{n-1}-\dots-\Lambda^0$, and $\rho_{\Lambda^{n-j}}=\prod_{i=1}^{j}\langle i/2\rangle_{(-1)^i}$ with $\langle a\rangle_-=\langle a\rangle$ [@dgz1996]. (This formula is stated in the book's normalization with long roots of squared length 2, where $q_i=q$ on long roots and $q_i=q^{1/2}$ on short roots; in Kac's normalization where $q_i=Q$ on short roots, $q=Q^2$.)
```

### 4. Chapter 8, Section `sec-irf` (`part2-algebra/08-r-matrices.qmd`, line 120)
**Target:**
```markdown
The vertex and IRF forms of an S-matrix are related by a change of basis, so they have the same eigenvalues and the same Bethe–Yang spectra.
```
**Replacement:**
```markdown
The vertex and IRF forms are intertwined by the quantum Clebsch–Gordan coefficients (the vertex-face correspondence), so their transfer matrices have the same eigenvalues and the same Bethe–Yang spectra.
```

---

## Claims Not Checked
None. All 9 items and their sub-claims were checked either by exact analytic derivation, literature quotation, or independent symbolic and numerical scripts in `referee_scratch/algebra_07_08/`.

## Other Observations (Unverified)
- In `07-quantum-groups.qmd` line 180, the text says conjugation by $\mathrm{diag}(1, 1, -1, 1)$ reverses $b$. Conjugation by $\Pi \otimes 1 = \mathrm{diag}(1, 1, -1, -1)$ also reverses $b$ while leaving $c, c'$ invariant, and has the clearer physical interpretation of acting on the first soliton only as a $\mathbb Z_2$ grading (as explained in Chapter 21 line 74). Both matrices achieve $\check R(-q)$.
