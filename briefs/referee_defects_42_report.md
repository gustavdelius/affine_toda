# Deep referee report on chapter 42: defects meet boundaries

Date: 2026-10-13 (review begun 2026-10-12). Model: GPT-6 (Codex). Starting commit: 54dc21087e351dcf60f1fce29e2ee93674861187. Reviewed chapter snapshot: 70b57c11a86309bfa62c409c59d707dd25c094fa, saved as referee_scratch/defects_42/ch42_reviewed_snapshot.qmd; SHA-256 dbaaa71b49e71ce56c1bd106ee1bfaa4aec3ba19638f75b6120e9bef10976a13. Concurrent agents advanced HEAD and inserted a figure during this review. Line numbers below refer to that snapshot unless another file is named.

This follows briefs/referee_defects_42.md. No pre-existing file was edited, no index- or worktree-changing Git command was run, and Quarto was not run. All scratch material is under referee_scratch/defects_42/. The book's script was read to audit its hypotheses, but was neither imported nor used as evidence. The independent scripts import only mathematical libraries and each other's new code.

Reproduce the computations with:

~~~bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 referee_scratch/defects_42/independent_checks.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 referee_scratch/defects_42/scalars_and_coideal.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 referee_scratch/defects_42/source_and_limits.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 referee_scratch/defects_42/tensor_and_cz.py
~~~

Complete output is retained in the corresponding .out files. Quantum calculations assume generic \(q\), nonzero finite \(y_d\), and the book's coproduct and slot ordering.

The principal finding is that the finite-rapidity no-go result can be proved without truncation, but the claimed restoration of full diagonal symmetry in the limits is false on the original Weyl space: the separately fitted Cartan operators fail the inverse relation. The fusion matrix, its tested scalar normalization, classical dressing, and crossing renormalization pass. Algebraic scope and several source descriptions need correction. Direct fusion also reveals two errors in the printed source formula CZ (2.19).

## 1. Refutation of cnj-twisted-diagonal and its scope

**Verdict: correct for torus and torus–node-exchange Hopf twists; unclear as an unrestricted conclusion justified solely by the printed numerical test.** The finite-rapidity conclusion survives, with an exact proof and explicit scope.

### What the existing test actually varies

In code/defect_doubling.py, generator_terms and ls_symmetry, both \(Y\) and \(\kappa\) are unknowns. Every allowed matrix entry is fitted independently, using the same entries at three rapidities. Neither operator is fixed to the original Weyl action. However:

* The test permits one prescribed shift of \(n\) for each unknown.
* Generators are fitted separately, without imposing common Cartan, inverse, commutator, or Serre relations.
* Failure rules out the tested equation; success alone does not construct a module.

For the negative conclusion, the charge-grade restriction can be justified. \(L=P\mathbb T\) preserves \(n+\sigma_z/2\). Decomposing arbitrary operators on the algebraic Weyl space into their matrix diagonals decomposes the equation into independent total-charge shifts. Its inhomogeneous term belongs to one shift; the other shifts solve homogeneous equations and cannot cancel it. This missing argument explains why testing the forced shift suffices for nonexistence here.

### Independent doubled matrix and finite-window computation

The new independent_checks.py constructs \(L\) directly from the four printed transmission amplitudes, then constructs

\[
K(y)=\begin{pmatrix}0&L(1/y)^{-1}\\L(y)&0\end{pmatrix}.
\]

It solves \(K(y)X(y)=X(1/y)K(y)\) with both operators unknown, at new rapidities \(-0.35,0.21,0.64\), \(\lambda=1.37,\theta_d=0.15\). The residual norm is divided by the inhomogeneous term.

| Failing equation | \(N=5\) | \(N=7\) | \(N=9\) |
|---|---:|---:|---:|
| I₋, \(e_1\), id | 0.893808 | 0.905240 | 0.910212 |
| I₋, \(f_1\), τ | 0.920097 | 0.922498 | 0.923595 |
| I₊, \(f_1\), id | 0.902643 | 0.913321 | 0.917960 |
| I₊, \(e_1\), τ | 0.857157 | 0.865835 | 0.869740 |

Complementary equations have residuals below \(8\times10^{-12}\); explicit-module checks below give \(4\times10^{-15}\). Different sampling explains the difference from 0.93–0.97; those residuals are not invariants.

Allowing **every matrix entry** of both operators at \(N=3\) gives 0.879842, 0.917374, 0.889577, 0.847138 in the four failing cases. Each equals its matching graded calculation to the printed precision (source_and_limits.py).

The inverse of the truncated \(L\) splits into independent two-dimensional charge blocks, apart from the endpoint singlets. The interior masks exclude those singlets. More decisively, the following obstruction uses no window.

### Independent exact obstruction, without a weight assumption

For I₋ write

\[
L_-(y)=D+yO,\quad
D=\begin{pmatrix}A&0\\0&D_0\end{pmatrix},\quad
O=\begin{pmatrix}0&B\\C&0\end{pmatrix},
\]

where \(A,D_0\) are invertible diagonal operators and \(B,C\) are nonzero bilateral shifts. Put \(E_+=e_{12},E_-=e_{21}\). Scalar transmission factors cancel.

For untwisted \(e_1\), allowing any nonzero torus factor \(t\), the equation is

\[
L_-(y)(yE_+\otimes1+k\otimes Y)
=(1\otimes Y+t yE_+\otimes\kappa)L_-(y).
\]

Its \(y^2\) coefficient requires

\[
\begin{pmatrix}0&0\\0&C\end{pmatrix}
=t\begin{pmatrix}\kappa C&0\\0&0\end{pmatrix},
\]

which is impossible for arbitrary operators \(Y,\kappa\).

For τ-twisted \(f_1\), the \(y^{-1}\) coefficient requires

\[
D(E_-\otimes\kappa)=t(E_+\otimes1)D.
\]

The nonzero inhomogeneous upper-right block has no counterpart on the left.

For I₊ write \(L_+(y)=D_++y^{-1}O_+\). The \(y^{-2}\) coefficient of untwisted \(f_1\) gives incompatible diagonal blocks; the \(y\) coefficient of τ-twisted \(e_1\) gives incompatible off-diagonal blocks. These are the analogous obstructions. One failing node suffices. Thus the no-go holds for arbitrary rapidity-independent actions on the algebraic Weyl space, not just graded finite matrices.

### Are these all relevant twists?

The source collection contains the missing classification citation: **Kolb 2014, Theorem 3.2**, quoting Twietmeyer, identifies the Hopf automorphisms of the standard generic derived Drinfeld–Jimbo algebra \(U_q(\mathfrak g')\) as torus rescalings semidirect diagram automorphisms. For affine \(sl_2\), these are precisely the two classes tested. The chapter should cite this and specify its algebra. Arbitrary algebra automorphisms, different quotients, or completed settings are not automatically covered. In particular \(\omega_q\) reverses the coproduct; it is not another Hopf automorphism of \(U\) with the same coproduct. [Kolb, §3.2](https://arxiv.org/abs/1207.6036).

### Paste-ready correction

At part8-defects/42-defects-and-boundaries.qmd:70 replace the status and statement with:

~~~markdown
**[Conjecture; false for the standard generic derived quantum affine algebra and Hopf twists]** The proposed full twisted-diagonal symmetry fails for both type-I transmission matrices. For torus twists and their compositions with the node exchange, comparison of powers of $y=e^{\lambda\theta}$ in a single generator equation gives an obstruction for arbitrary $\theta$-independent operators on $\mathcal W$. These exhaust the Hopf automorphisms of the standard derived Drinfeld–Jimbo algebra in the setting of [@kolb2014, Thm. 3.2]. The finite-window test below corroborates the obstruction; it is not by itself a classification of other algebra automorphisms or completions.
~~~

At line 75 append the charge-grade argument and coefficient comparison above. Update “refuted as stated by the test” at line 245 and in the status/conventions summaries on integration. Keep the recorded numerical residuals as results of that particular sampling.

## 2. prp-diagonal-borel

**Verdict: wrong as a whole because (5) overstates symmetry restoration; (1), the displayed I₋ instance of (2), and obstruction (4) are correct. The finite-rapidity claim (3) is correct within item 1's scope.**

### (1): diagonal Borel

For \(b\in U_q(\mathfrak b^-)\), set \(\chi(\Delta b)=\rho(b)\). Expanding \(\Delta_A\Delta(b)\) gives \(\Delta(b)\) on the first particle copy and \(\Delta^{\rm op}(b)\) on the reflected copy. The block equation is

\[
L(y)(\pi_y\otimes\rho)\Delta(b)
=(\pi_y\otimes\rho)\Delta^{\rm op}(b)L(y).
\]

The inverse block follows by replacing \(y\) with \(1/y\) and inverting \(L\). Direct substitution of \(\rho(f_1)=\mu y_d^{-1}S^{-1}\), \(\rho(f_0)=\mu y_d^{-1}S\), and \(\rho(k_1)|n\rangle=q^{2n+s}|n\rangle\) verifies this. I₊ uses the positive Borel identically. The proof in the chapter is sufficient.

### (2): the second module

Applying the printed operators to a basis vector gives

\[
\tilde k_0\tilde k_1=1,\quad
\tilde k_1\tilde e_1=q^2\tilde e_1\tilde k_1,\quad
\tilde k_1\tilde e_0=q^{-2}\tilde e_0\tilde k_1,\quad
\tilde e_1\tilde e_0=\tilde e_0\tilde e_1=abq^2.
\]

The cubic Serre polynomial consequently reduces to a common operator times \(1-[3]_q+[3]_q-1=0\). Also

\[
Y_{n+1}-Y_n=a q^{-2n}(q^{-2}-1)=-iy_dq^{-2n-2-s}.
\]

Cartan covariance forces the additive constant in \(Y_n\) to vanish for generic \(q\). Substitution verifies both \(e_j\) and both \(k_j\) intertwining equations, giving a module rather than unrelated fitted operators. Independent output:

~~~text
second-module relations/Serre: 1.4131e-15
explicit diagonal b-minus and twisted b-plus intertwining: 3.3894e-15
~~~

The text already proves the displayed I₋ formula. Keeping the I₊ counterpart Numerical is reasonable until its formula or a fully specified transport is supplied.

### (4): incompatibility of the halves

With \(G=k_1\otimes k_1\), \(x=e_1\otimes1+k_1\otimes e_0\),

\[
e_1\otimes1=\frac{GxG^{-1}-q^{-2}x}{q^2-q^{-2}}.
\]

Any subalgebra containing the two I₋ halves therefore contains this element. Its intertwining equation lacks the \(E_+\otimes\kappa\) output term, leaving the nonzero \(y^2\) coefficient in item 1 uncancelled. This proves the arbitrary-action obstruction.

For I₊, conjugation of \(f_1\otimes k_0^{-1}+1\otimes f_0\) by \(G\) separates its weights \(q^{-2}\) and \(q^2\). This isolates \(1\otimes f_0\), whose equation has an uncancellable \(y^{-2}\) coefficient. The chapter prints only the I₋ proof; this sentence completes the other type.

The relation \(\rho(k_1)\tilde\rho(k_1)=q^4\tilde\rho(k_1)\rho(k_1)\) gives a second obstruction for the specific displayed actions. There is no contradiction with chapter 40's no-extension argument: the second Borel uses a different embedding and different Cartan action, whereas that argument fixes the untwisted Weyl action.

### (5): full symmetry is not restored on the original Weyl space

There are two errors. First, I₊ has off-diagonal entries proportional to \(y_d/y\), so the directions are reversed: it is diagonal at \(\theta_d\to-\infty\), charge-conjugating at \(+\infty\).

Second, even for I₋, separate generator fits fail to give a common full \(U\)-module. In the diagonal limit,

\[
L_\infty=\operatorname{diag}(A,D_0),\quad
A|n\rangle=q^{-n}|n\rangle,\quad D_0|n\rangle=q^{n+s}|n\rangle.
\]

The coefficients of \(y\) in the \(e_1\) equation and \(y^{-1}\) in the \(f_1\) equation respectively force

\[
\chi(\Delta k_1)=AD_0^{-1},\qquad
\chi(\Delta k_1^{-1})=D_0^{-1}A.
\]

These are the **same operator, not inverses**. Their product on \(|n\rangle\) is \(q^{-4n-2s}\), generically not 1. A fixed torus twist \(t_1=t\) changes the product to \(t^{-2}q^{-4n-2s}\), still not 1 for all \(n\). A fixed \(n\) sector instead has its own twist \(t_n=q^{-2n-s}\) and the trivial boundary module. That sectorwise statement is correct.

At the opposite limit \(L_0=\left(\begin{smallmatrix}0&B\\C&0\end{smallmatrix}\right)\), \(B\) lowers \(n\) and \(C\) raises it. For a τ-twisted diagonal the coefficients force

\[
\chi(k_1\otimes k_0)=CB^{-1},\qquad
\chi(k_1^{-1}\otimes k_0^{-1})=B^{-1}C.
\]

Both raise \(n\) by 2, so their product raises it by 4. They cannot be inverses. The shifted conjugation on particle–defect states is not the line with \(\mathcal W=\mathbb C,\mathbb T=C\) in prp-diagonal-coideal.

The individual-generator residuals do decay. At \(|\theta_d|=2,4,6\):

~~~text
I-minus, e/id, theta_d positive: 0.157434, 0.0102998, 0.000665098
I-minus, f/tau, theta_d negative: 0.0910717, 0.00590276, 0.000381150
I-plus, f/id, theta_d negative:  0.156486, 0.0102995, 0.000665098
I-plus, e/tau, theta_d positive: 0.139420, 0.00909241, 0.000587125
~~~

Successive ratios approach \(e^{2\lambda}=15.49\). This measures the matrix limit, not the missing module relations. The incorrect I₊ directions instead approach residuals 0.994 and 1. Source_and_limits.py records the inverse obstruction; scalars_and_coideal.py tests the correct I₊ directions.

### Paste-ready corrections

At part8-defects/42-defects-and-boundaries.qmd:91 replace (5):

~~~markdown
5. For type I$_-$ the transmission matrix tends to a direct sum of diagonal twist matrices as $\theta_d\to+\infty$, and to a charge-conjugating operator with a shift of the defect charge as $\theta_d\to-\infty$. The directions are reversed for type I$_+$. The relative off-diagonal or diagonal corrections are of order $e^{-\lambda|\theta_d|}$ at fixed particle rapidity. Each fixed-charge diagonal limit has the twisted-diagonal symmetry of its own invertible twist line. This does not assert a full untwisted or $\tau$-twisted diagonal action on the original Weyl space: separately fitted generators need not share a Cartan operator satisfying the inverse relation.
~~~

At line 97 replace the interpretation of the limit residuals with: “These are separate-generator residuals and do not test the relations of a full diagonal module.” At lines 13, 99, 103 (figure caption), and 121 (last symmetry-table row), use the same matrix-limit statement. The shifted conjugation limit cannot be identified with a scalar-boundary-module line without an additional reduction of defect degrees of freedom. Add the I₊ obstruction sentence above to the proof at line 95.

## 3. prp-defect-fusion-k

**Verdict: correct for the matrix and tested scalar identification; unclear as a derivation of an exact finite-coupling UV–IR relation from a classical limit. The qualitative type-II comparison is correct, but literal CZ (2.19) has two errors.**

### Independent fusion block

Set \(d=y_d\), \(r=e^{i\eta_D}\), \(H(z)=\left(\begin{smallmatrix}1&-iz\\-iz&1\end{smallmatrix}\right)\), \(D_n=\operatorname{diag}(q^{-n},q^{n+1+s})\). On the sector \(\{v_+|n\rangle,v_-|n+1\rangle\}\), \(L_n(y)=H(y/d)D_n\). For \(k_+=yr+y^{-1}r^{-1}\), \(k_-=yr^{-1}+y^{-1}r\), direct multiplication gives

\[
K_{f,n}^{\rm mat}=
\frac{D_n^{-1}}{1+y^{-2}d^{-2}}
\begin{pmatrix}
k_++d^{-2}k_-&-ir d^{-1}(y^2-y^{-2})\\
-ir^{-1}d^{-1}(y^2-y^{-2})&k_-+d^{-2}k_+
\end{pmatrix}D_n.
\]

The symbolic check in independent_checks.py returns True. A rapidity-independent diagonal similarity removes \(D_n\), removes the \(r\) factors off the diagonal, and reverses their common sign. At \(\vartheta=\log d=\lambda\theta_d\), the numerator then equals \(4/d\) times the BPTT matrix \(P_0^\pm,Q_0\).

~~~text
GZ gauge invariants, maximum error over windows 5,7,9: 4.9760e-16
vartheta=-lambda*theta_d, minimum error: 0.135388
vartheta=theta_d, minimum error: 0.0172450
~~~

Thus the map and the two proposed negative controls pass independently, including complex rapidity. The similarity depends on the charge sector, not on rapidity.

### Scalar normalization

Canceling the common \(R_0(u)\sigma(\eta_D,u)\) factors reduces the full-matrix scalar comparison to

\[
\boxed{\displaystyle
\sigma(i\lambda\theta_d,-i\theta)=
\frac{1+d^{-2}}{1+y^{-2}d^{-2}}
\frac{\Phi(\theta-\theta_d)}{\Phi(-\theta-\theta_d)}.}
\]

Scalars_and_coideal.py independently sums the printed Gamma products with 35-digit arithmetic and common cutoffs \(40,80,160,320,640\), then extrapolates in \(1/N\).

| \(\lambda,\theta,\theta_d\) | Last raw cutoff error | Four-cutoff extrapolation | Five-cutoff extrapolation |
|---|---:|---:|---:|
| \(1.37,0.33,0.21\) | \(8.34\,10^{-5}\) | \(3.87\,10^{-14}\) | \(3.01\,10^{-15}\) |
| \(1.37,0.21+0.17i,0.21\) | \(6.83\,10^{-5}\) | \(2.40\,10^{-14}\) | \(2.39\,10^{-15}\) |
| \(0.73,0.29+0.11i,-0.38\) | \(7.83\,10^{-5}\) | \(1.01\,10^{-13}\) | \(3.38\,10^{-14}\) |

This supports no extra CDD factor **for the specified minimal factors**, at accuracy exceeding the book's check. It checks the entire gauged matrix because the remaining scalar is common to all entries. Finite sampling is not a global meromorphic identity proof or uniqueness theorem against arbitrary CDD modifications; the Numerical label remains appropriate.

### Corrigan–Zambon in their own conventions

Their matrix rows are incoming soliton labels. With \(Q^2=q^{-1}\), \(\nu^{1/2}=v\), the type-I charge block, as a column-action matrix, is

\[
T_I(x)=\begin{pmatrix}Q^\alpha/v&\varepsilon xQ^{\alpha+2}\\
\varepsilon xQ^{-\alpha}&vQ^{-\alpha-2}\end{pmatrix}.
\]

I multiplied \(T_I(1/x)^{-1}\operatorname{diag}(rx+r^{-1}x^{-1},rx^{-1}+r^{-1}x)T_I(x)\). Its diagonal ratio and off-diagonal product agree symbolically with (2.18): output [0, 0]. Translation in the even sector is \(x=y,\ q_{\rm CZ}=q,\ \alpha=2n,\ r=e^{i\eta_D},\ \nu=1,\ \varepsilon=-iq/d\), with a diagonal charge-dependent similarity. In particular

\[
r_+^{\rm CZ}=q\{y(r+d^{-2}r^{-1})+y^{-1}(r^{-1}+d^{-2}r)\}.
\]

CZ do not prove the chapter's scalar identity there; the chapter correctly treats that comparison as additional work.

For type II I used their displayed \(T_{II}\) at
\((q,b_+,b_-,r,x)=(e^{0.71i},0.37+0.19i,-0.28+0.13i,e^{0.23i},e^{0.41})\).
Direct fusion gives diagonal ratios, for \(\alpha=0,2,4,6\),

~~~text
0.992460-0.122569i, 0.943454+0.331503i,
0.909352+0.416027i, 0.998339-0.057617i.
~~~

The genuine charge dependence survives diagonal similarity. However **literal (2.19) fails**. I visually inspected printed PDF pages 5–6 to exclude extraction errors. Substitution of their type-II parameters into (2.17) gives

\[
h=b_-c_+q^{-1}+b_+c_-q
=-(|b_-|^2+|b_+|^2)q
=-(|b_-|^2+|b_+|^2)Q^{-2}.
\]

The printed expression has an extra \(Q^{-2}\) on \(|b_-|^2\). Also every term of \(s_-\) needs the opposite overall sign: its first term, \(d_-c_+r\), is \(-\bar b_-Q^{-4}r\), but is printed positive. Correcting **both** discrepancies makes all entries of direct fusion proportional to the formula, with residuals \(1.76\,10^{-16},2.63\,10^{-16},1.93\,10^{-16}\) at \(\alpha=0,2,4\). Printed-form residuals are 0.260, 0.175, 0.560. Correcting \(h\) alone leaves the off-diagonal-product mismatch. Evidence: source_and_limits.py and tensor_and_cz.py; [CZ source](https://arxiv.org/abs/1202.6016).

### Sinh-Gordon and the UV–IR consequence

Independent substitution of the elementary transmission factor gives

\[
\frac{T(\theta)}{T(-\theta)}
=\frac{\cosh\theta_d+i\sinh\theta}{\cosh\theta_d-i\sinh\theta}
=\frac1{(1-F)(1+F)},\qquad F=2i\theta_d/\pi,
\]

where \((z)=\sinh(\theta/2+i\pi z/4)/\sinh(\theta/2-i\pi z/4)\). Three real/complex samples give errors \(1.57,2.48,1.11\times10^{-16}\). Bajnok–Simon (16), (17), and (19)–(20) support the transmission, quantum parameter matching, and Neumann CDD statements; their \(B_{\rm BS}\vartheta=-\theta_d\).

Combining \(\vartheta=\lambda\theta_d\) with the **conjectural finite-coupling matching** \(\vartheta/(\lambda+1)=\theta_d^{\rm cl}\) gives the chapter-37 relation, also obtained by \(B\mapsto-2/\lambda\). But a classical limit alone fixes only leading asymptotics. Adding \(c\theta_d^{\rm cl}/\lambda^2\) to \((1+1/\lambda)\theta_d^{\rm cl}\) is invisible in that limit. Exact matching additionally requires the normal-ordered quantum boundary-coupling identification. “Consistent exactly when” is too strong as derived.

The general fusion theorem's unitarity also follows independently: for \(F(\theta)=T(-\theta)^{-1}K(\theta)T(\theta)\),
\[
F(-\theta)F(\theta)=T(\theta)^{-1}K(-\theta)K(\theta)T(\theta)=1.
\]
The RTT/mixed identities in item 5 provide the train argument for the reflection equation, corroborated by item 7's independent tensor check.

### Paste-ready corrections

At part8-defects/42-defects-and-boundaries.qmd:138 replace (4), and use the same qualification in its proof at line 142:

~~~markdown
(4) The parameter map in (3) is consistent with the conjectural defect UV–IR relation: if quantum fusion and boundary UV–IR matching identify $\vartheta/(\lambda+1)=\theta_d^{\rm cl}$, then $\lambda\theta_d=(\lambda+1)\theta_d^{\rm cl}$, precisely @cnj-defect-parameter. Classical fusion checks the leading semiclassical limit; the exact finite-coupling relation additionally assumes the normal-ordered quantum parameter matching.
~~~

Change the label for (4) at line 134 to “Conditional on the UV–IR matching”. At line 244 append:

~~~markdown
In the printed type-II formula [@cz2012, eq. (2.19)], direct substitution into their (2.17) instead gives $h=-(|b_-|^2+|b_+|^2)Q^{-2}$ and the negative of the displayed expression for $s_-$. With these two corrections, direct fusion agrees entry by entry up to the common scalar. The irremovable charge dependence remains.
~~~

The boxed scalar identity and independent convergence results can be appended at line 145. No correction to \((\eta,\vartheta)=(\eta_D,\lambda\theta_d)\) is needed.

## 4. prp-reflected-copy and prp-diagonal-coideal

**Verdict: correct for the reflected coproduct and diagonal construction, and supported for the restricted maximal-coideal calculation; wrong if read as unrestricted maximality. A cross-referenced definition of the conjugate module needs a spectral-parameter correction.**

### Independent coproduct calculation

Conjugating
\(\check R_{12}(\pi_1\otimes\pi_2)\Delta(a)=(\pi_2\otimes\pi_1)\Delta(a)\check R_{12}\)
by the source and target flips gives \(P\check RP\) intertwining \(\Delta^{\rm op}\). Reflection reverses the original particle order, so this is the needed action. On mixed copies the algebra factors act independently and scattering is the flip. This establishes \(U\otimes U^{\rm cop}\); replacing it by \(U\otimes U\) fails generically (not at every cocommutative degeneration).

On generators,
\[
(\omega_q\otimes\omega_q)\Delta^{\rm op}(e)
=1\otimes f+f\otimes k^{-1}=\Delta(f).
\]
The other generators work identically. Thus \(F=\omega_q\otimes{\rm id}\) satisfies
\(\Delta_{U\otimes U}F=(F\otimes F)\Delta_A^{\rm op}\).
It reverses the two paired tensor factors and maps left coideals to right coideals.

The fully indexed diagonal formula is
\[
\Delta_A\Delta(a)=\sum(a_{(1)}\otimes a_{(4)})\otimes(a_{(2)}\otimes a_{(3)}),
\]
equivalent to the chapter's shorthand with \(\Delta(a_{(2)})\). Applying \(1\otimes\psi\) proves the twisted identity. On the counit boundary module the two particle actions are \(\pi(a)\) and \(\pi(\psi(a))\), so the implementation by \(\mathbb T\) proves invariance. The chapter's Theorem parts (1)–(2) are proved.

### Conjugate module

For \(\pi_y(e_j)=yE_j,\pi_y(f_j)=y^{-1}F_j\),
\[
(\pi_{1/y}\circ\omega_q)(e_j)=yF_j,\quad
(\pi_{1/y}\circ\omega_q)(f_j)=y^{-1}E_j,\quad
(\pi_{1/y}\circ\omega_q)(k_j)=q^{-H_j}.
\]
These are the actual principal-gradation antisoliton matrices, so chapter 42 is correct. But chapter 33, line 72, literally calls \(\bar V(y)\) the composition \(V(y)\circ\omega_q\), which would instead have \(e_j=y^{-1}F_j\); chapter 34, line 72, uses similar shorthand. The grading shift must be specified.

Paste-ready replacement at part7-half-line/33-soliton-reflection.qmd:72, and the analogous sentence at part7-half-line/34-quantum-symmetric-pairs.qmd:72:

~~~markdown
In the principal gradation the antisoliton at parameter $y$ has $e_j=yF_j$, $f_j=y^{-1}E_j$ and $k_j=q^{-H_j}$; thus $\bar V(y)=V(1/y)\circ\omega_q$. Equivalently, first conjugate the unshifted module and then apply the grading shift.
~~~

### Maximality: exact finite-space check and global counterexample

The function exact_coideal in scalars_and_coideal.py independently constructs the 21 PBW monomials, uses exact rational arithmetic at the generic test value \(q=2\), computes the one-particle commutant, expands the algebraic coproduct, and imposes closure under every second-leg coefficient map:

~~~text
one-particle commutant dimension: 17
after coideal closure:            5
span equals {1,Delta(k),Delta(k^-1),Delta(e),Delta(f)}: True
~~~

A Hopf twist of the second factor transports the result to the torus and node-exchange cases. This supports the restricted calculation. However \(X\) and the surviving five-dimensional space are not subalgebras: they do not contain \(e^2\). The computation concerns a **coideal subspace in \(X\)**, not a classification of all coideal subalgebras.

Unrestricted maximality would be false for the unquotiented affine algebra on the level-zero probes in the chapter. Let \(c=k_0k_1\). The group-like \(c\otimes1\) acts as the identity on both particle copies and every tensor product of them. Adjoining \(c^{\pm1}\otimes1\) to \(\Delta(U)\) gives a larger invariant left coideal subalgebra. Yet \(c\otimes1\notin\Delta(U)\): its two partial counits would require a preimage simultaneously equal to \(c\) and 1. Explicitly quotienting \(c=1\) removes this counterexample, but does not promote the finite-\(X\) result to global maximality. Appel–Vlaar 4.1.3 concerns a universal fixed-point condition, not the commutant in these nonfaithful probes.

At part8-defects/42-defects-and-boundaries.qmd:51 replace the beginning of (4) with:

~~~markdown
(4) Let $X$ be the 21-dimensional span of products $g\otimes h$ with $g,h\in\{k^{-1},1,k,e,f\}$ and at most one occurrence of $e$ or $f$ in total. For the stated probes, the largest left coideal subspace $S\subset X$ whose represented action intertwines the doubled $K$ is $(1\otimes\psi)\Delta(U)\cap X$. This is a calculation in $X$, not a maximality theorem for all coideal subalgebras of the affine algebra.
~~~

Retain the listed five spanning elements. At lines 115–117 change “maximal coideal [Numerical]” to “maximal coideal subspace in \(X\) [Numerical]”.

### Source scope

Bao–Wang explicitly restrict to finite type in §1.2 and identify the twisted quasi-R-matrix in Remark 4.10(2). Their coproduct and \(\omega\) agree with the dictionary. Wang's Example 0.1 and §2 repeat it. Kolb §2.1 and Appel–Vlaar 2022 §3.6 explicitly assume indecomposability, so their stated results do not automatically cover the two-copy affine Cartan matrix. The Nichols meaning of “diagonal type” in [Kolb–Yakimov](https://arxiv.org/abs/1901.00490) is indeed different.

The Notes incorrectly credit Kolb 2014 with universal K-matrices; that paper constructs the coideals. Item 8 supplies the correction.

## 5. RTT inside the doubled reflection-equation algebra

**Verdict: correct as an equivalence of presentations after an off-diagonal quotient and boundary unitarity; wrong or imprecise as a literal embedding or ambient Hopf-subalgebra assertion.**

### Independent algebraic calculation

From \(R_{12}L_1L_2=L_2L_1R_{12}\), adjoining two-sided matrix inverses gives by multiplication:

\[
L_2^{-1}R_{12}L_1=L_1R_{12}L_2^{-1},\qquad
R_{12}L_2^{-1}L_1^{-1}=L_1^{-1}L_2^{-1}R_{12}.
\]

These are the mixed and reverse-order inverse block relations, with the appropriate spectral arguments. Therefore the doubled reflection equation for
\(K(y)=\left(\begin{smallmatrix}0&L(1/y)^{-1}\\L(y)&0\end{smallmatrix}\right)\)
reduces to RTT and its inverse consequences. **Unitarity**, not off-diagonality alone, identifies the two blocks as inverses.

Setting the diagonal generators to zero imposes relations: it takes a quotient of the universal RE algebra. The resulting presentation with unitarity maps to the RTT algebra with matrix inverse adjoined. No embedding of that quotient back into the original algebra is supplied. “Inverse adjoined” can mean a universal two-sided matrix inverse or a specified formal completion; an Ore-localization theorem is not proved.

Also, for \(x=\Delta(e)\),
\[
\Delta_A(x)=(e\otimes1)\otimes1+(k\otimes1)\otimes\Delta(e)
+(k\otimes e)\otimes\Delta(k).
\]
The first legs are not generally in the diagonal Borel. Hence its diagonal image is a coideal, **not a Hopf subalgebra of \(A\)**. It has an abstract transported Hopf structure through its algebra isomorphism with the Borel; that coproduct is not the restriction of \(\Delta_A\).

### Paste-ready replacement

Replace the paragraph at part8-defects/42-defects-and-boundaries.qmd:108:

~~~markdown
**The RTT presentation of the off-diagonal doubled boundary.** Impose vanishing diagonal blocks and boundary unitarity, so that the remaining blocks are $\mathbb L(\theta)=P\mathbb T(\theta)$ and $\mathbb L(-\theta)^{-1}$. The doubled reflection equation is then equivalent to RTT for $\mathbb L$ and the mixed and inverse relations obtained by multiplying RTT by its two-sided inverses. Thus this quotient of the reflection-equation presentation is the RTT algebra with matrix inverses adjoined; no embedding into the unrestricted reflection-equation algebra is asserted. Likewise $\Delta(U_q(\mathfrak b^\mp))$ is a left coideal in $U\otimes U^{\rm cop}$ and is abstractly isomorphic to the Borel algebra. Its transported Hopf coproduct differs from the ambient coproduct $\Delta_A$, under which it is not a Hopf subalgebra.
~~~

Use this quotient formulation in the corresponding chapter-40 table entry on integration.

## 6. Classical dressing and prp-fused-bdg

**Verdict: correct as a conservation theorem. The explicit reservation about involution and quantization is appropriate.**

### Independent finite differences

The function classical in independent_checks.py builds roots and weights in the sum-zero subspace of \(\mathbb R^{r+1}\), constructs \(D,E,L_d\), differentiates the printed defect potential, and independently solves for constant \(K_\psi\) by equating coefficients of the boundary exponentials. It chooses random root-space fields and \(\dot\psi\), and imposes

\[
D^T\dot\phi=-\nabla B_\psi+E\dot\psi-\nabla_\psi B_d,\qquad
\phi_x=E\dot\phi+D\dot\psi-\nabla_\phi B_d.
\]

Central differences differentiate the dressed matrix along that velocity. For \(\beta=0.73,\sigma=0.81,z=0.94+0.28i\), relative residuals are:

| Rank and boundary | \(h=10^{-4}\) | \(h=10^{-5}\) | \(h=10^{-6}\) | Direct \(L\) control | Bare \(K\) control |
|---|---:|---:|---:|---:|---:|
| \(a_1,C=(0.61,-1.17)\) | \(5.98\,10^{-9}\) | \(9.66\,10^{-11}\) | \(9.07\,10^{-10}\) | 5.63 | 1 |
| \(a_1,\) Neumann | \(1.79\,10^{-10}\) | \(1.01\,10^{-10}\) | \(2.64\,10^{-9}\) | 1.77 | 1 |
| \(a_2,C=(1,-1,1)\) | \(1.14\,10^{-8}\) | \(1.20\,10^{-10}\) | \(9.31\,10^{-11}\) | 2.07 | 1 |
| \(a_2,\) Neumann | \(2.30\,10^{-10}\) | \(6.21\,10^{-11}\) | \(5.26\,10^{-10}\) | 1.86 | 1 |

The step dependence is compatible with central-difference error followed by roundoff; neither control passes.

### Independent proof and oscillator comparison

Put \(M(z)=L_d(z)^{-1}\). Then \(\dot M=M a_\phi-a_\psi M\). Differentiating \(K_\phi=M(1/z)^TK_\psi M(z)\) and subtracting \(K_\phi a_\phi(z)+a_\phi(1/z)^TK_\phi\) leaves

\[
-M(1/z)^T\{a_\psi(1/z)^TK_\psi+K_\psi a_\psi(z)\}M(z)=0.
\]

This proves compatibility and trace conservation for every invertible defect Lax matrix with the stated evolution, independently of rank or defect type. The proof is conditional on the already stated Lax and bare-boundary identities, exactly as it should be.

Expanding exponentials independently gives

\[
e^{-t}\cos\frac{\beta(\phi+\psi)}2+e^t\cos\frac{\beta(\phi-\psi)}2
=e^{-i\beta\phi/2}\cos(\beta\psi/2-it)
+e^{i\beta\phi/2}\cos(\beta\psi/2+it).
\]

The symbolic check returns True. Multiplying by \(-2m_0/\beta^2\) establishes the potential identity up to the vacuum constant. Using BDG's \(\dot f=\{H,f\}\),

\[
\frac{d}{dt}(\hat p-\hat q)=
-\frac{m_0}{2}\left[e^{i\beta\phi/2}\sin\hat q+e^{-i\beta\phi/2}\sin\hat p\right].
\]

At \(m_0=\beta=1,\phi=0.31,\psi=0.72,\theta_d=0.43\) it is \(-0.316667\), not zero. Neumann fusion instead gives \(\dot\phi=-\partial_\psi B_d\) and a \(\dot\psi\) term in the boundary condition. Thus equal potentials do not give equal dynamical systems. BDG's original (7), (9), (10) agree. [BDG source](https://arxiv.org/pdf/nlin/0201007).

The printed Theorem proof is sufficient. Line 168 explicitly excludes involution and quantization and should retain that reservation. No change to the inverse/transposition order is needed. Any summary elsewhere calling the higher-rank result “integrable” should instead say “a dynamical boundary with a conserved Lax generating function” until the Poisson reduction and involution are established.

## 7. Tensor K-matrices and non-invertible boundaries

**Verdict: correct for the tensor identities and crossing renormalization; wrong in the reason given for excluding jump modules and in the global description of the \(p=4\) off-diagonal scenario. Tensor-K identification is only up to scalar.**

### Appel–Vlaar source audit and independent tensor calculation

The supplied 2402.16676.txt supports these precise citations:

| Citation | Checked content |
|---|---|
| Theorem 5.6.1 | Tensor factorization, canonical lift to boundary/bulk module categories, identities (5.5)–(5.8). |
| Theorem 6.4.2 | Formal Laurent spectral action on \(M\otimes V(z)\). |
| Theorem 6.5.2 | Rational tensor K from normalized R and trigonometric K for a boundary module restricting an irreducible bulk module; a formal scalar factor remains. |
| Remark 6.5.3(a) | Ito–Terwilliger's restriction result for finite-dimensional irreducible q-Onsager modules. |
| Remark 6.5.3(b) | One-dimensional intertwiner space suffices for normalization. |
| §1.6.2 | Proposed Borel second leg for boundary Q-operator applications. |

The general theorem allows **boundary coideal modules which need not extend to \(U\)**. Line 182's claimed requirement of a full bulk module on the boundary leg contradicts the theorem and line 178. A given Borel module supplies no boundary-coideal action by itself; that is the actual missing input. Infinite dimension or failure to extend to \(U\) is not alone an exclusion theorem.

Tensor_and_cz.py independently solves the bulk coproduct-intertwiner equations for \(\check R\), constructs the fused \(4\times4\) matrix, and solves the boundary-intertwiner system. For generic complex \(\hat\epsilon\), \(\theta_d=0.29\):

~~~text
theta=0.37:         residual 7.84e-16; intertwiner dimension 1
theta=-0.28+0.13i:  residual 9.17e-16; intertwiner dimension 1
tensor reflection equation: 1.19e-15
~~~

This supports the example, but uniqueness identifies it only up to scalar, not with a specified universal normalization.

### Independent boundary-crossing derivation

Let \(A=d_a,C=d_c,B=d_b,D=d_d\), \(t=i\theta/(2\pi)\), and \(F_b(z)=(AC/B^2)^{iz/(2\pi)}\). Then

\[
\frac{F_b(i\pi/2-\theta)}{F_d(i\pi/2+\theta)}
\sqrt{\frac DB}
=\left(\frac{BD}{AC}\right)^{i\theta/\pi}.
\]

The constant factor is \((D^2/B^2)^{-1/4}\sqrt{D/B}=1\). The rapidity factor is \([(AC/B^2)(AC/D^2)]^{-t}\). This is exactly the multiplier of \(\tilde S^{ba}_{cd}(2\theta)\), with the exercise's index order. There is no leftover constant. A numerical test with arbitrary positive dimensions gives \(1.11\times10^{-16}\); reversing the bulk factor gives 0.0294.

For unitarity the product of boundary multipliers, with intermediate label \(c\), is \((d_b/d_d)^{i\theta/(2\pi)}\), independent of \(c\), and equals 1 on \(\delta_{bd}\). In the boundary Yang–Baxter equation the multipliers are ratios of external kink-state normalizations: internal factors telescope along either scattering path. This proves the claimed preservation of both equations.

SY (3.22)–(3.28) give \(C^B_{ab}=1/(d_bF_b^B)\) and the boundary-factor cancellation. Their (2.32) is equivalent to the exercise's bulk convention. Weak symmetry and footnote 13 are accurately quoted.

### Source-description errors

The vacuum labels carry a fusion action, not generally a permutation: \(L_2\otimes L_2=L_1\oplus L_3\). For \(p=4\), SY define

\[
R=U^1_{22}=U^3_{22},\quad P=U^2_{11}=U^2_{33},\quad
V=X^2_{13}=X^2_{31},\quad R=P+V.
\]

Scenario 2 is \(P=0,V=R\ne0\). Only the outer-vacuum-2 block is off-diagonal; outer-vacuum-1,3 reflections stay diagonal and nonzero. Thus “purely off-diagonal reflection” globally is wrong. The \(p\ge5\) statement is correct. Formula (4.29) covers both \(p=4\) scenarios; (4.35) concerns \(p\ge5\).

### Paste-ready corrections

At part8-defects/42-defects-and-boundaries.qmd:180 replace its last sentence:

~~~markdown
This corroborates identification with the trigonometric tensor K-matrix up to scalar; a scalar normalization must still be chosen.
~~~

At line 182 replace the opening sentence:

~~~markdown
For a jump defect the given structure on $\mathcal W$ is a Borel-module structure. Applying Appel and Vlaar's tensor K-matrix theorem would require a compatible boundary-coideal module in their category $\mathcal W_\theta$, which has not been specified here. Their general theorem does not require that boundary module to extend to the full bulk quantum group; the stronger restriction hypothesis belongs to Theorem 6.5.2.
~~~

At line 190 replace “whose \(p-1\) vacua are permuted by” with “whose \(p-1\) vacua carry the fusion action of”. At line 196 replace the last two sentences:

~~~markdown
For $p=4$, the tricritical Ising model, the scenarios are $V=0,\ P=R\ne0$ or $P=0,\ V=R\ne0$ [@shimamori_yamaguchi2025, eqs. (4.26)–(4.27)]. In the second scenario only the outer-vacuum-$2$ block is off-diagonal; the outer-vacuum-$1,3$ amplitudes remain diagonal. Both scenarios use the minimal scalar solution (4.29); equation (4.35) gives the $p\ge5$ solution.
~~~

No correction to eq-dbd-sy-crossing or exr-dbd-boundary-crossing is needed.

## 8. Status labels, exercises, and bibliography

**Verdict: wrong in several scope and attribution details; the core Theorem claims have proofs, with the small I₊ completion in item 2.**

### Proof/status audit

* prp-diagonal-coideal (1)–(2): proved by the coproduct formula and implementation of the twist. Part (4) is a restricted coideal-subspace calculation.
* prp-diagonal-borel (1): proved. Part (4): I₋ is proved; the analogous I₊ proof is supplied above. The displayed I₋ second module is also proved by the text's algebra. Part (5)'s full-module reading is false despite its Numerical label.
* prp-defect-fusion-k (1): the train argument is a short valid proof sketch under the stated independent boundary/defect slots; unitarity is explicit above. Scalar identification is appropriately Numerical. Part (4) must be conditional on quantum UV–IR matching.
* prp-fused-bdg: product differentiation proves compatibility and conservation, not Poisson involution or quantization. The explicit reservation makes the limited Theorem label justified.
* The preregistered conjecture should not be called refuted solely by least squares. Item 1 supplies the missing exact argument under specified hypotheses.

### Independent exercise spot-checks

1. **exr-dbd-mixed-twist, line 217: correct.** Direct substitution gives the stated difference equation and \(\kappa_n\). Item 2 verifies all module relations. Applying both Cartan orders to \(|n\rangle\) gives the \(q^4\) factor.
2. **exr-dbd-e-one, line 221: correct.** Items 1–2 give the exact linear combination isolating \(e_1\otimes1\) and its uncancellable \(y^2\) term, with no truncation.
3. **exr-dbd-invertible, line 213: correct.** The twist is \(t_1=e^{2i\alpha},t_0=t_1^{-1}\); charge conjugation gives \(\tau\). For arbitrary torus factors, a spectral change \(y'=\sqrt{t_0t_1}\,y\) and diagonal gauge ratio \(\sqrt{t_1/t_0}\) implement the action. A fixed-rapidity twist therefore requires \(t_0t_1=1\), with the associated sign/gauge choice.
4. **exr-dbd-classical-dressing, line 229: correct.** For \(a_2\) Neumann fusion the two components of \(\psi(t)\) are boundary variables, and \(D^T\dot\phi-E\dot\psi+\nabla_\psi B_d=0\), \(\phi_x=E\dot\phi+D\dot\psi-\nabla_\phi B_d\). These are exactly the equations used in the new finite differences.
5. **exr-dbd-boundary-crossing, line 233: correct.** Item 7 supplies the full factor cancellation and unitarity/BYBE argument.
6. **exr-dbd-fusion-gz, line 225: matrix part correct; final instruction too strong.** It asks for an exact finite-coupling relation from a classical limit. The additional quantum matching hypothesis is necessary.

Paste-ready replacement for the last sentence of exr-dbd-fusion-gz:

~~~markdown
Check semiclassical compatibility with $C_{0,1}=\cosh(\theta_d^{\rm cl}\pm i\beta_{\rm SG}\phi_D/2)$. Then, assuming normal-ordered quantum UV–IR matching $\vartheta/(\lambda+1)=\theta_d^{\rm cl}$, derive $\lambda\theta_d=(\lambda+1)\theta_d^{\rm cl}$ and compare with @cnj-defect-parameter. Explain why the classical limit alone cannot determine this exact finite-coupling relation.
~~~

These are independent computations for this item, not just an inventory of exercises.

### Citation-by-citation audit

All 19 bibliography keys cited in the chapter exist. The following records content verification, not re-proofs of the cited literature.

| Key | Source inspected and finding |
|---|---|
| bajnok_george2006 | hep-th_0404199.txt, doubling and (19): attribution correct. |
| bajnok_simon2008 | 0712.4292.txt, (5), (7), (16)–(20): fusion, transmission, parameter matching and Neumann CDD claims correct. |
| cz2012 | 1202.6016.txt and PDF printed pages 5–6: (2.9), (2.18), qualitative charge dependence correct; literal (2.19) errors in item 3. |
| bdg2002 | [Original paper](https://arxiv.org/pdf/nlin/0201007), (7), (9), (10), dynamical reflection construction: correct, with its Hamiltonian sign convention. |
| baseilhac_koizumi2003 | [Original paper, §§3–4](https://arxiv.org/html/hep-th/0208005v2): quantum boundary operators, coideal intertwining and reflection matrices support the claim; q-Onsager is later terminology. |
| baseilhac_koizumi2005 | [Original paper](https://arxiv.org/pdf/hep-th/0507053): deformed Onsager symmetry in open XXZ, suitable for the later algebraic identification. |
| bao_wang2018b | 1610.09271.txt, §1.2 and Remark 4.10(2): finite-type scope, diagonal embedding, quasi-R identification correct. |
| wang2023 | 2112.10911.txt, Example 0.1 and §2: correct. |
| bao_wang2018 | 1310.0103.txt: pioneering type-A/coideal construction for type-B Kazhdan–Lusztig theory; not itself the general Kac–Moody quasi-K theorem. Narrow attribution. |
| balagovic_kolb2019 | 1507.06276.txt: general quasi-K and finite-type universal K; no diagonal-type specialization claimed there. |
| kolb2014 | 1207.6036.txt: indecomposability and coideal construction correct; universal-K attribution wrong. Theorem 3.2 supplies the missing Hopf-automorphism citation. |
| appel_vlaar2022 | 2007.09218.txt, §3.6 and universal construction: indecomposable scope and universal-K attribution correct. |
| kolb_yakimov2020 | [Original abstract](https://arxiv.org/abs/1901.00490): different Nichols meaning of diagonal type confirmed. |
| appel_vlaar2024 | 2402.16676.txt: all numbered citations checked in items 4 and 7; category-exclusion and normalization wording need correction. |
| sklyanin1988 | Author-uploaded [original article](https://www.researchgate.net/publication/231073719_Boundary_conditions_for_integrable_quantum_systems), Proposition 2, (20)–(21), printed p. 2379: dressing attribution correct. The article, despite unreliable hosting-page metadata, is the 1988 paper. |
| runkel2010 | 1004.1909.txt, (1.22), (3.38): fusion onto \((1,1)\) gives \((1,s)\) with coupling \(\lambda+\tilde\lambda\); correct. |
| cck2024 | 2403.04835.txt: modified-crossing and Ward-identity attribution correct. |
| cck2025 | 2408.13132.txt: extended scattering/non-invertible-symmetry and modified-crossing attribution correct. |
| shimamori_yamaguchi2025 | 2504.08375.txt: cited equation numbers and footnote 13 checked; crossing factors correct; the \(p=4\) description and “permuted” require correction. |

### Paste-ready bibliographic correction

At part8-defects/42-defects-and-boundaries.qmd:238 replace the clause beginning “the general quasi-K-matrix” with:

~~~markdown
Bao and Wang introduced the quasi-K intertwiner in their pioneering quantum-symmetric-pair construction [-@bao_wang2018]; Balagović and Kolb developed the general quasi-K construction and universal K-matrices in finite type [-@balagovic_kolb2019]. Kolb constructed the Kac–Moody quantum symmetric pair coideals [-@kolb2014], and Appel and Vlaar constructed universal K-matrices of Kac–Moody type [-@appel_vlaar2022].
~~~

At line 242, and the analogous sentence at line 60, replace “so their universal K-matrices do not cover” with “so their coideal and universal-K constructions, respectively, do not directly cover”. This preserves the indecomposability qualification without misattribution.

## Other observations — unverified

* This review does not identify the fused operator with an ordinary boundary as a complete physical theory, including its boundary-state spectrum and all multi-particle identifications. Sectorwise scattering-matrix equivalence is narrower.
* The higher-rank rotation-twisted second module, its representation relations, and the physical interpretation of the I₋ second module as the alternate two-label I₊ reading were not independently audited here.
* No global analytic uniqueness theorem for the scalar, full infinite-dimensional fused-boundary crossing bootstrap, or Poisson reduction/involution theorem for the fused field theory is supplied.
* The figure inserted concurrently was read through its caption; its SVG geometry and plotting code were not independently inspected. Its full-symmetry-limit language needs the same correction as proposition (5).

## Integration (2026-10-13)

The referee's four scripts are in `briefs/referee_defects_42_scripts/` with their outputs; the coordinator reran them and they reproduce the saved outputs byte for byte. They were not used as verification. The integrating session checked every item with its own scripts, written from the chapter's definitions and importing nothing from the referee: `integration_obstructions.py` (sympy: the $y$-power obstructions for arbitrary operators on $\mathcal W$ as noncommuting symbols, both types; the I$_+$ part of (4); $\tilde\rho$ exactly for symbolic $n$, all relations and all four generator equations, and the free constant; the limits, the forced $\chi(\Delta k_1)$, $\chi(\Delta k_1^{-1})$ and the $\tau$-limit operators), `integration_cz2012.py` (sympy, $\bar b_\pm$ as independent symbols: direct fusion (2.9) of the general (2.14) with (2.16) against (2.17), and of the type-II and type-I specializations against (2.19) and (2.18)) and `integration_scalar_identity.py` (the scalar identity and the matrix step $k_++y_d^{-2}k_-=4y_d^{-1}P_0^+$); all are in the same directory with `*_output.txt`. Further checks: $V(1/y)\circ\omega_q=\bar V(y)$ and $V(y)\circ\omega_q=\bar V(1/y)$ with `boundary_common.an_vector` (2e-16), and the $c\otimes1$ counterexample by hand. Source passages read: Kolb 2014 Sec. 2.1 (indecomposable $A$), Sec. 3.1–3.2 and Thm. 3.2 (Twietmeyer); CZ 2012 Sec. 2, (2.4)–(2.19), on the rendered PDF pages 5–7 (printed pp. 4–6); Appel–Vlaar 2024 Sec. 1.4–1.6, Def. 5.2.1, Prop. 5.2.3, Thm. 5.6.1, Thm. 4.1.3, Thm. 6.5.2, Rem. 6.5.3; Shimamori–Yamaguchi (2.1)–(2.6), (4.1), (4.8)–(4.12), (4.24)–(4.35); the abstracts and introductions of Bao–Wang 1310.0103, Balagović–Kolb 1507.06276, Kolb 1207.6036 and Appel–Vlaar 2007.09218 (Sec. 3.6, indecomposable $A$).

- **Item 1. Verified.** The $y$-power obstruction was re-derived by hand and checked with noncommuting symbols: with $\mathbb L=\mathcal D+y\mathcal O$ (I$_-$), the $y^2$ coefficient of the untwisted $e_1$-equation has the block $C$ on $v_-\otimes\mathcal W$ and the $y^{-1}$ coefficient of the $\tau$-twisted $f_1$-equation the block $-D/t$ from $v_-$ to $v_+$, whatever $Y$ and $\kappa$; for I$_+$ ($\mathcal D+y^{-1}\mathcal O$) the blocks $-B/t$ ($y^{-2}$, untwisted $f_1$) and $A$ ($y$, $\tau$-twisted $e_1$). Kolb Thm. 3.2 (quoting Twietmeyer) is as stated, over $\mathbb K(q)$ for $U_q(\mathfrak g')$; the book says it uses it for generic $q$. The charge-grade argument was checked by hand. Changes: `cnj-twisted-diagonal` label **[Conjecture; refuted for Hopf twists, exact argument below]**; the paragraph after it; "The test" (Twietmeyer/Kolb, $\omega_q$ is not a Hopf twist, the charge-shift argument); `prp-diagonal-borel` (3) now states the scope and is proved in the proof; Notes (book-internal bullet); `index.qmd` (only the clause); the ch. 40 table status cell; status file and conventions §12, label-rule line.
- **Item 2. Verified, with one modification.** (a) I$_+$ directions reversed: confirmed (off-diagonal entries $\propto y_d/y$; the I$_+$ single-generator residuals decay towards $\theta_d\to-\infty$ for $f$/id and $+\infty$ for $e$/$\tau$, and stay at 1.00 in the opposite directions). (b) No full twisted diagonal in the limits: in the diagonal limit the $y$ and $y^{-1}$ coefficients force $\chi(\Delta k_1)=t^{-1}AD^{-1}$ and $\chi(\Delta k_1^{-1})=t^{-1}D^{-1}A$ (the same operator, product $t^{-2}q^{\mp(4n+2s)}$); in the off-diagonal limit $t^{-1}CB^{-1}$ and $t^{-1}B^{-1}C$ both raise $n$ by 2. Added beyond the report: in each limit the twist of the other kind fails already for single generators. Each charge sector is the twist line $t_1=1/t_0=q^{-2n-s}$ (I$_-$). (4) for I$_+$ ($1\otimes f_0$, $y^{-2}$ block $-C$) verified and added. The I$_-$ second module: all four generator equations and all relations verified exactly for symbolic $n$, and the free constant $c(1-q^4)=0$; the proof block now contains the complete argument. *Modification:* the replacement (5) is worded as a matrix-limit statement plus the no-module statement, labelled Theorem (proof in the block); the report's sentence about relative corrections $e^{-\lambda|\theta_d|}$ is kept as a numerical remark. Changes: label **[Theorem (1), (3)–(5), and (2) for type I$_-$; Numerical (2) for type I$_+$]** (was [Theorem (1), (4); Numerical (2), (3), (5)]); items 3 and 5; the whole proof rewritten in the notation $\mathbb L=\mathcal D+y\mathcal O$; the numbers paragraph; the chapter opening; "Reading the result"; the figure lead-in, caption and alt text; the figure script's labels ("matrix limit", "$C$ with a charge shift", "line $C$", "each $n$") and docstring (5 PASS, only its own SVG regenerated); `tbl-dbd-symmetries` rows (type I$_\mp$ jump [Theorem], last row rewritten); "Beyond sine-Gordon" (the $a_r$ limit is the rotation dressed by a charge shift). Ch. 43 lines 152, 208 and 218 do not repeat the "restored" claim and were not changed for it.
- **Item 3. Verified.** The UV–IR point: a classical limit fixes only $\theta_d/\theta_d^{\rm cl}\to1$. CZ (2.19): confirmed by direct fusion in CZ's own conventions (PDF read): (2.17) is right for the general (2.14) (exact, sympy); at the type-II parameters it gives $h=-(|b_-|^2+|b_+|^2)Q^{-2}$, and (2.19) matches the fusion entry by entry only with this $h$ *and* $s_-$ of the opposite sign; either correction alone fails (gauge-invariant deviations 2.0 and 0.10); since (2.17) is the paper's own formula, (2.19) is the side that is wrong (*Misprint*). (2.18) also checked exactly. The boxed scalar identity was derived by hand from the sector blocks ($k_++y_d^{-2}k_-=4y_d^{-1}P_0^+$) and checked to 7e-9 at new points. Changes: `prp-defect-fusion-k` label **[Theorem (1); Numerical (2)–(3); (4) by hand, given the quantum UV–IR matching]** (was "by hand (4)"), (4) and its proof, the "new points" sentence, the pointer after "(2.19)", the scalar identity in the check paragraph, `exr-dbd-fusion-gz`, Notes bullet, Appendix D new entry *Corrigan and Zambon 2012*.
- **Item 4. Verified.** $V(1/y)\circ\omega_q=\bar V(y)$ (chapter right); ch. 33 l. 72 and ch. 34 l. 72 are wrong as literally written. Not edited (Part VII is on another branch): exact replacement texts recorded under "Global pending edits" in the status file. Maximality: the coideal-subspace restatement and the $c^{\pm1}\otimes1$ counterexample verified by hand. Changes: `prp-diagonal-coideal` (4) restated; the paragraph after it (finite space, not subalgebras; the counterexample; Appel–Vlaar 4.1.3 "mirrors"); `tbl-dbd-symmetries` row 1.
- **Item 5. Verified.** The mixed relations follow from RTT by multiplication with inverses; setting generators to zero is a quotient; $\Delta_A(\Delta(e))$ computed by hand. Changes: the paragraph, retitled "The RTT presentation of the off-diagonal doubled boundary"; ch. 40 table row "Quadratic algebra".
- **Item 6. Verified; no chapter edit.** No Part VIII text calls the fused higher-rank system integrable. Part VII ch. 35 problem 9 does ("integrable dynamical $a_n^{(1)}$ boundaries"): replacement recorded as a pending Part VII edit. Conventions §12 gains the wording rule.
- **Item 7. Verified.** Appel–Vlaar Thm. 5.6.1 acts on $M\otimes V$ for $M\in\mathcal W_\theta$, modules of $U_q(\mathfrak k)$; the restriction hypothesis is in Thm. 6.5.2, which identifies the tensor K-matrix only up to a formal scalar. Shimamori–Yamaguchi: vacua $|a\rangle=L_a|\tilde0\rangle$ carry the fusion action; for $p=4$ scenario 2 only the outer-vacuum-2 block ($X^2_{13}$) is off-diagonal, $U^{1}_{22}=U^3_{22}=R\ne0$ stays diagonal; (4.29) is the $p=4$ solution, (4.35) the $p\ge5$ one. Changes: "proportional to the trigonometric tensor K-matrix", the jump-defect paragraph of `sec-tensor-k`, "carry the fusion action", the $p=4$ consequences.
- **Item 8. Verified.** Attributions read in the four abstracts/introductions (Bao–Wang introduced $\Upsilon$; Balagović–Kolb: general construction and universal K in finite type; Kolb: Kac–Moody coideals; Appel–Vlaar: universal K of Kac–Moody type). Changes: Notes paragraph, "What is proved, and where" and the Notes bullet ("their constructions of the coideals and of the universal K-matrices, respectively, do not directly cover"), ch. 43 problem 12, Appendix D (Bao–Wang entry), `exr-dbd-fusion-gz` (item 3).

**Other observations of the report.** The higher-rank rotation-twisted second module belongs to the chapter-43 report and was not touched. The figure-language remark is covered by item 2.

**Book script after integration:** `defect_doubling.py` 20 PASS (was 17), 15 s. New or changed checks: the exact obstruction (all four failing generators, (4) for both types, both limits, forced Cartan operators); $\tilde\rho$ exactly for symbolic $n$; (4) numerically for both types with the algebra step; the limits for both types and the opposite directions; the sector twist lines; CZ (2.17)/(2.19)/(2.18) against direct fusion. Figure script 5 PASS. `defect_common.py` unchanged. Reference check over all 65 book `.qmd` files: 1162 ids, no duplicates, every `@sec/eq/prp/fig/exr/thm/cnj/tbl` reference resolves, every citation key is in `references.bib`.

**Labels changed.** `cnj-twisted-diagonal` (wording of the bracket); `prp-diagonal-borel` (3), (5) and (2) for type I$_-$ upgraded to Theorem; `prp-defect-fusion-k` (4) made conditional; `tbl-dbd-symmetries` status cells (rows 1, 5, 6, 7); ch. 40 `tbl-boundary-defect-algebra` status cell. The upgrades have had no independent referee pass beyond this report, which supplied the arguments for (3), (4)-I$_+$ and (5).

**Pending.** The three Part VII edits in the status file (ch. 33 l. 72, ch. 34 l. 72, ch. 35 problem 9).
