# Independent referee report: chapters 38 and 40

Date: 2026-10-12. Model: GPT-6 (Codex). Initial review commit: d855e742e0929a7a85a24a66a3c13eee73a1ba6d. Commit on resuming the review: c2dd82b0c4584c9b4ec08964cc5178b863a7fa6b. Commit at final content validation: b3bbc3572652cede5a25441681627b7192ad5fae.

The requested report path did not exist at the start. I followed the existing brief, briefs/referee_defects_38_40.md, which specifies this path as its output. No pre-existing book, brief, or textbook-code file was edited; no index-changing git command or Quarto render was run. Other work advanced the repository during this review. Chapters 38 and 40 remained identical to the saved review copies. Their SHA-256 hashes are respectively 777dcfdf78049c189f005c18f32a5ec62f21a06971ee3bafc01df624f2fe985e and f3779a45b07d2c6775534c7a5b6bc2d0c65878c81a2656f4d1f5177261a6570f. Line references below refer to these copies. References to chapters 30 and 36 refer to the files inspected on resumption.

All independent scripts, outputs, and four explicit spin-two corrections are under [referee_scratch/defects_38_40](../referee_scratch/defects_38_40/). Run scripts using:

~~~sh
PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_38_40/<script>.py
~~~

They use Python, NumPy, SymPy and SciPy. None imports the book's computations except the explicitly separated comparison with boundary_common.rcheck in final_checks.py; that comparison is not the derivation of the R-matrix. Exact symbolic zeros, numerical residuals, and formal arguments are distinguished below.

The principal classical correction to BCZ is supported: I independently constructed both local spin-two corrections for both defect types in \(a_2\), and checked them on a transmitted soliton. Several subsidiary claims need repair. The most substantial are the insufficient root-of-unity exclusion in the uniqueness theorem, the distinction between an independent defect Poisson algebra and the full equal-time field theory, the scope and meaning of Weyl traces, and the Ising disorder-field exercise. The all-rank local-charge construction and all-length TQ interpretation retain qualifications described below.

## 1. Both spin-two charges, the cubic obstruction, and the structural argument

**Verdict: correct for the independently constructed \(a_2\) charges and the all-\(a_r\) cubic cancellation; wrong tensor normalization in the text; unclear as an unconditional theorem about all local charges at both spectral ends.** The chapter appropriately labels the latter argument a sketch.

**Independent construction.** [spin_check.py](../referee_scratch/defects_38_40/spin_check.py), with full output in [spin_output.txt](../referee_scratch/defects_38_40/spin_output.txt), derives the bulk current before imposing the defect. It uses an orthonormal basis of the sum-zero subspace of \(\mathbb R^3\), whose defining weights are the rows of

\[
U=\begin{pmatrix}1/\sqrt2&1/\sqrt6\\-1/\sqrt2&1/\sqrt6\\0&-2/\sqrt6\end{pmatrix},\quad
\alpha_j=U^T(e_j-e_{j+1}),\quad
D=\begin{pmatrix}1&1/\sqrt3\\-1/\sqrt3&1\end{pmatrix},\quad B=1-D.
\]

For this calculation \(b=m=1\), \(\sigma\ne0\), and \(d_e=\partial_t+e\partial_x\), \(e=\pm1\), without a \(1/\sqrt2\). Solving bulk conservation gives

\[
C(u)=\frac{\sqrt6}{9}u_1(3u_0^2-u_1^2)
=\frac23\sum_j(\mu_j\cdot u)^3,\quad
T_e=C(d_e\phi)+(d_e^2\phi)^TBd_e\phi,\quad
\Theta_e=-(d_e\phi)^TB\nabla V.
\]

Thus \(\partial_{-e}T_e=\partial_e\Theta_e\), and \(Q_{-2e}=\int(T_e-\Theta_e)\,dx\) has derivative \(eG_e\), where \(G_e=(T_e+\Theta_e)_{\rm left}-(T_e+\Theta_e)_{\rm right}\). After substituting the sewing conditions and bulk equations, I integrate the acceleration coefficient of \(G_e\); the remaining one-form in the four boundary fields is closed. The resulting \(F_e\) obeys **exactly** \(G_e-dF_e/dt=0\), so \(\mathcal B_{-2e}=-eF_e\).

The four explicit answers, in both SymPy and pasteable LaTeX, are:

| Type | Charge | Explicit correction |
|---|---|---|
| I\(_-\) | \(+2\) | [type_+1_spin_+2.txt](../referee_scratch/defects_38_40/type_+1_spin_+2.txt) |
| I\(_-\) | \(-2\) | [type_+1_spin_-2.txt](../referee_scratch/defects_38_40/type_+1_spin_-2.txt) |
| I\(_+\) | \(+2\) | [type_-1_spin_+2.txt](../referee_scratch/defects_38_40/type_-1_spin_+2.txt) |
| I\(_+\) | \(-2\) | [type_-1_spin_-2.txt](../referee_scratch/defects_38_40/type_-1_spin_-2.txt) |

In those files \(\phi=(p,r)\), \(\psi=(v,w)\), and p1,r1,v1,w1 are their time derivatives. The filename's type_+1 denotes the original sewing sign, hence book I\(_-\). The answers are linear in velocities plus exponentials, with arbitrary symbolic \(\sigma\). For nonzero \(b,m\), use dimensionless fields \(b\phi,b\psi\) and coordinates \(mx,mt\); multiplying the dimensionless charges and corrections by \(m^2/b^2\) gives the current normalized with cubic term \(bC\). Holomorphic continuation gives imaginary coupling. No numerical fit enters these corrections.

For an independent dynamical check, I used the exact complex \(a_2\) soliton

\[
\phi=-\sum_{j=0}^2\alpha_j\log(1+e^{2\pi ij/3}E),\quad
E=e^{\sqrt3\cosh\theta\,x-\sqrt3\sinh\theta\,t+\xi},\quad
\psi=\phi\big|_{E\mapsto zE},
\]

with \(\theta=0.9,\xi=0.2+0.45i,\sigma=0.71\). I integrated the bulk densities on \([-24,0]\cup[0,24]\) at five times. Time reversal, with the reciprocal parameter in the starting solution, gives I\(_+\).

| Type | Spin | Maximum sewing residual | Bulk charge variation | Corrected charge variation |
|---|---:|---:|---:|---:|
| I\(_-\) | \(+2\) | \(3.16\times10^{-15}\) | 26.198 | \(4.27\times10^{-14}\) |
| I\(_-\) | \(-2\) | \(3.16\times10^{-15}\) | 0.71582 | \(4.87\times10^{-15}\) |
| I\(_+\) | \(+2\) | \(2.84\times10^{-15}\) | 0.56566 | \(8.62\times10^{-15}\) |
| I\(_+\) | \(-2\) | \(2.84\times10^{-15}\) | 20.702 | \(5.48\times10^{-14}\) |

The nonintegrable control \(D_{\rm bad}=1-0.6B\) retains \(D_{\rm bad}+D_{\rm bad}^T=2\). Its two flux Euler derivatives have residuals 1.16358 and 1.16027, excluding local primitives. Substituting the same transmitted soliton into its sewing conditions fails by 4.51571. This is a control of integrability and sewing, not a claim that this soliton solves the modified theory.

**All-rank cancellation and the precise source error.** For nonzero Fourier modes, \(D\) has eigenvalues \(d_k=1+i\cot(\pi k/h)\). The cubic invariant only couples \(k_1+k_2+k_3=0\bmod h\). Writing \(c_i=\cot(\pi k_i/h)\), each disputed coefficient is

\[
\prod_i(1+ic_i)+\prod_i(1-ic_i)
=2[1-c_1c_2-c_2c_3-c_3c_1]=0.
\]

The last equality follows by expanding \(\sin(a+b+c)=0\) and dividing by \(\sin a\sin b\sin c\). [classical_checks.py](../referee_scratch/defects_38_40/classical_checks.py) independently checks ranks 2–9. At rank 2, SymPy gives \(C(Du)+C(D^Tu)=0\) exactly, while

\[
C(Du)=-\frac{8\sqrt6}{81}u_0(u_0+\sqrt3u_1)(\sqrt3u_0-3u_1)
\]

is nonzero. BCZ's error is the sentence immediately after (5.11), inferring nonvanishing from \(D=2-D^T\). That implication is false for this cubic invariant. In Penati–Polvara App. B the same tensor occurs in (B.4); polarization cancels its mixed cubic terms. With \(E=F=1-D\), \(D,D^T\) commute with \(B\), so the acceleration terms cancel too. Their assertion following (B.4) that it remains nonzero is false. Vanishing of highest terms alone would not establish a conserved charge; the exact primitive above supplies the missing lower-order check in \(a_2\). This does not settle their quantum anomaly calculations. Sources: [BCZ 2004](https://arxiv.org/abs/hep-th/0401020), [Penati–Polvara 2019](https://arxiv.org/abs/1902.10690).

There is a factor-of-three error in the book's tensor identification. BCZ (5.2) has \(T=(1/3)A_{abc}u_au_bu_c+B_{ab}v_au_b\). Their (5.5), with \(B=1-D\), therefore requires

\[
A(u,u,u)=2\sum_j(\mu_j\cdot u)^3=3C(u),
\]

not \(A=C\). This error does not change the zero obstruction, but makes the cited tensor equation and the exercise false as written.

**Structural argument.** In the defining representation \(S\) is a cyclic shift, \(S^h=1\), and

\[
\widetilde L_-=1+\frac\sigma\lambda S^T
=\frac\sigma\lambda S^T\left(1+\frac\lambda\sigma S\right).
\]

This proves invertible leading Laurent coefficients at both ends for \(\sigma\ne0\), not holomorphy of \(L\) at both ends. The text's careful Laurent formulation is correct. An invertible matrix can have zero diagonal entries, so generic nonvanishing of \(N_{kk}\) needs an extra argument: at the vacuum \(S\) is diagonal in the principal eigenbasis and those leading entries are nonzero; analyticity then proves generic nonvanishing. The passage to a regularized infinite-line generating function, and exceptional jumps, remains a sketch. “All spins” must mean spins present in the bulk hierarchy; \(a_1\) has no spin-two charge.

**Paste-ready corrections.** At part8-defects/38-defect-conditions.qmd:163, replace the opening tensor identification, and make the same normalization replacement at line 328:

~~~markdown
With the normalization $T_{\pm3}=\frac13A_{abc}\partial_\pm\phi_a\partial_\pm\phi_b\partial_\pm\phi_c+B_{ab}\partial_\pm^2\phi_a\partial_\pm\phi_b$ of [@bcz2004b, eq. (5.2)], the tensor is $A(\vec u,\vec u,\vec u)=2\sum_{k=1}^h(\vec\mu_k\cdot\vec u)^3$. Its contribution to the density is therefore $\frac23\sum_k(\vec\mu_k\cdot\vec u)^3$. With $B=1-\mathsf D$, this tensor satisfies their (5.5).
~~~

At line 165, replace the inference from invertibility to nonzero diagonal entries with:

~~~markdown
At the vacuum, the leading defect matrix is diagonal with nonzero entries in the principal eigenbasis. Its leading diagonal entries remain nonzero for generic defect values by analyticity. This gives a formal local expansion on that open set; control of the regularized infinite-line limit and of exceptional jumps remains part of the sketch.
~~~

At line 353, qualify the all-rank conclusion:

~~~markdown
The purported cubic obstruction to conserving both spin-$\pm2$ charges vanishes for every $a_r^{(1)}$ with a nonzero cubic invariant. Both local defect corrections have been constructed explicitly for $a_2^{(1)}$, with numerical evidence for $a_3^{(1)}$; the general two-chirality hierarchy is supported by the formal argument in @prp-spin-two.
~~~

## 2. The type-I classification and the affine-diagram scan

**Verdict: correct under the stated ansatz.** It is not a classification of defects with auxiliary fields, singular kinetic matrices, different theories on the two sides, or arbitrary nonexponential potentials.

**Independent computation.** classical_checks.py constructs finite simple-root systems, untwisted extensions, dual systems and \(BC_r\) twisted extensions without importing book data. It enumerates Gram-matrix-preserving permutations, solves \(D(1-P)=2\), and tests forbidden products. [classical_output.txt](../referee_scratch/defects_38_40/classical_output.txt) gives:

~~~text
23 diagrams, 78 automorphisms
fixed node: 64; singular: 5; forbidden cross terms: 2; pass: 7
passes: a1 exchange; a2, a3, a4 rotations by +1 and -1
~~~

For \(a_4^{(1)}\), \(p:j\mapsto j+2\), the computed \((\alpha_i\cdot D\alpha_j)\) is

\[
\begin{pmatrix}
2&0&-2&2&-2\\-2&2&0&-2&2\\2&-2&2&0&-2\\-2&2&-2&2&0\\0&-2&2&-2&2
\end{pmatrix}.
\]

It contains ten forbidden entries. A fixed-point-free automorphism is insufficient.

The finite exponential ansatz does not hide cancellations: the map from \((x,y)\) to the joint exponent of \(f_xg_y\) is injective. Adding its two field coefficients recovers \(x\), since \(D+D^T=2\), and then recovers \(y\). Bulk terms require \(x=\alpha_i\). An extra \(f\) exponent would have to be orthogonal to all required \(g\) exponents, which span the space, so it is zero (an irrelevant constant). Extra \(g\) exponents are excluded similarly. Thus unwanted products can be tested separately.

If a permutation has multiple orbits, a proper orbit sum of roots is a nonzero fixed vector: the unique affine dependence has positive coefficients on every node. Hence \(1-P\) is singular. A finite tree with more than two vertices has an invariant centre and no transitive automorphism. This leaves the \(a_r\) cycle and the equal two-node \(a_1\) diagram. The support condition and \(D+D^T=2\) force the partner to be adjacent, leaving two orientations. This agrees with CZ (2.21)–(2.35), rather than relying on a rank-four scan for the all-rank result. [CZ 2009a](https://arxiv.org/abs/0902.1307).

For an order-\(k\) twisted presentation with literature long-root length squared \(2k\), set \(\alpha_{\rm lit}=\sqrt{k}\alpha\), \(\beta_{\rm lit}=\beta/\sqrt{k}\). Then \(\beta_{\rm lit}\alpha_{\rm lit}=\beta\alpha\); the exponentials are unchanged and the common Gram-matrix factor cannot alter the classification. The \(a_3^{(2)}=d_3^{(2)}\) identification and omission of duplicate \(b_2^{(1)}=c_2^{(1)}\) presentations account for the count.

No correction to prp-ar-only is required. The nearby exclusion at 38-defect-conditions.qmd:145 overlooks its mirror. Replace that sentence with:

~~~markdown
Within this field-independent Lax ansatz, the allowed matrices are (-@eq-d-matrix) and its transpose, up to the diagram relabelling that exchanges the two orientations. Other choices with $\mathsf D+\mathsf D^{\mathsf T}=2$ do not satisfy the corresponding compatibility equations under the same assumptions.
~~~

## 3. Dirac brackets, the transport sign, and involution

**Verdict: correct as a conditional Poisson-algebra construction; wrong transport sign in chapter 30's proof; unclear if promoted to the full equal-time field-theory bracket.** Chapter 38's explicit caveat is necessary.

**Independent derivation.** Write the defect kinetic term as \(\frac12z^TK\dot z\), \(z=(\phi,\psi)\). Its canonical one-form is \(a_i=-\frac12(Kz)_i\). For constraints \(\chi_i=p_i-a_i\),

\[
\{\chi_i,\chi_j\}=\partial_i a_j-\partial_j a_i=K_{ij},
\qquad\{z_i,z_j\}_D=(K^{-1})_{ij}.
\]

At rank one, I\(_-\) has \(K=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)\), hence \(\{\phi,\psi\}=-1\). Reversing the kinetic term reverses this bracket. This is valid for independent defect coordinates; it does not define a continuum bulk bracket restricted to a point.

classical_checks.py differentiates the entries of \(L_d\), contracts with \(K^{-1}\), and independently builds chapter 30's \(r\) from matrix units. For ranks 1–4, both types and both couplings, the residual in

\[
\{L_{d1},L_{d2}\}=+[r,L_{d1}L_{d2}]
\]

is at most \(6.48\times10^{-16}\). The independent canonical calculation for \(a_x\) gives

\[
\{a_{x1}(x),a_{x2}(y)\}=[r,a_{x1}(x)+a_{x2}(y)]\delta(x-y)
\]

with residual below \(1.61\times10^{-16}\). These establish the coefficient and sign without fitting them. They are finite-rank checks, not a supplied all-rank component proof.

The transport sign follows exactly. Put \(A=-a_x\), so \(\partial_xT=AT\). Then \(\{A_1,A_2\}=-[r,A_1+A_2]\delta\). Differentiating

\[
T_1(b,x)T_2(b,x)\,r\,T_1(x,a)T_2(x,a)
\]

gives the insertion \([r,A_1+A_2]\). The negative sign in the local bracket therefore yields \(\{T_1,T_2\}=-[r,T_1T_2]\). Inverting \(L_d\) gives precisely this negative quadratic bracket. If the three factors Poisson-commute, \(T_+L_d^{-1}T_-\) also has it, and tracing the commutator proves involution. This last implication is exact.

**Paste-ready corrections.** At part7-half-line/30-boundary-conditions.qmd:120, replace the final sentence with:

~~~markdown
For the convention $\partial_xT=-a_xT$ and the stated $r$, the canonical local bracket gives $\{T_1,T_2\}=-[r,T_1T_2]$. The Leibniz derivation of involution uses this sign consistently. Reversing the overall sign of $r$ leaves the homogeneous classical reflection equation and the conclusion of involution unchanged.
~~~

This does not alone invalidate the boundary classification or commuting traces. It matters when bulk and defect algebras are combined: the sign must not be changed in only one.

At 38-defect-conditions.qmd:15, replace the claim that the field-theory action has induced the required bracket by:

~~~markdown
The defect kinetic term defines a finite-dimensional symplectic form on independent defect coordinates. With its inverse bracket, the defect matrix satisfies the quadratic algebra with exactly the bulk $r$-matrix. The relation of this algebra to the full equal-time field theory requires the qualifications in @sec-defect-r-matrix.
~~~

At line 191, use **[Conditional proposition; numerically checked]**, retaining the independence assumption. The time-like canonical-transformation argument and cited sine-Gordon results are distinct from completing an \(a_r\) equal-time reduction. [Caudrelier 2015](https://arxiv.org/abs/1411.5171).

## 4. The \(a_r\) delay and twisted Bianchi identities

**Verdict: correct, including the minus sign in the self-conjugate delay.**

**Independent computation.** classical_checks.py substitutes independently built soliton fields and analytic derivatives into both sewing equations for every species at ranks 1–4. Maximum residual: \(2.14\times10^{-14}\). Delay identities have residual at most \(1.24\times10^{-14}\); the nonintegrable-\(D\) control fails by 4.51571.

Set \(t=e^{-(\theta-\theta_d)}\), \(a=\pi a_{\rm species}/h\). Direct multiplication gives

\[
z_a=\frac{t-ie^{ia}}{t-ie^{-ia}},\quad
z_a z_{h-a}=\frac{t^2+2t\sin a+1}{t^2-2t\sin a+1}
=\frac{\cosh(\theta-\theta_d)+\sin a}{\cosh(\theta-\theta_d)-\sin a}.
\]

Conjugation and \(e^{i\pi(h-a)/h}=-e^{-ia}\) prove the conjugate-species identity. With \(\varepsilon_a=-ie^{ia}\), \(\mathsf a(\theta)=-e^{-\theta}\), direct division proves the twisted Bianchi form. For \(a=h/2\), \(z=(t+1)/(t-1)=-\coth((\theta-\theta_d)/2)\).

BCZ (6.2) uses the opposite rapidity sign; CZ 2007 (2.20) requires \(a\mapsto h-a\), consistent with transposing \(D\); CZ 2009a (4.10) agrees in the dictionary. Substitution into CZ 2007 (2.20) verifies that their (2.21) should contain \(-\coth\). Sources: [BCZ 2004](https://arxiv.org/abs/hep-th/0401020), [CZ 2007](https://arxiv.org/abs/0705.1066), [CZ 2009a](https://arxiv.org/abs/0902.1307).

No replacement of eq-dc-ar-delay or its identities is needed. They are meromorphic identities away from delay poles and singular soliton profiles.

## 5. Type II: sectors, normalization, and corrections to CZ 2009b

**Verdict: correct for the SG and Tzitzéica formulas and source corrections; correct identification of the Bristow–Bowcock series.** This does not establish every generalized type-II construction or its Liouville integrability.

**Independent computation.** [type_two_checks.py](../referee_scratch/defects_38_40/type_two_checks.py), [type_two_output.txt](../referee_scratch/defects_38_40/type_two_output.txt), differentiates \(f,g\) and verifies exactly

\[
\tfrac12(f_qg_\lambda-f_\lambda g_q)=V(p+q)-V(p-q)
\]

in dimensionless source variables for both models. It then tests all three sewing equations on solitons, using analytic bulk derivatives and a fourth-order derivative for the auxiliary field.

~~~text
SG momentum residual: 0
Tz momentum residual: 0
SG tanh sector: corrected sewing 2.38e-11; printed-source control 7.23914
SG coth sector: sewing 9.18e-11
Tz three sectors: 7.42e-11, 9.03e-11, 6.56e-11
Tz z1*z2*z3 - 1: 6.68e-16
Tz printed-f momentum control: 0.146580
~~~

For the no-jump SG vacuum, \(e^{\hat\lambda}=+1/(\sqrt2\cosh\tau)\) gives coth–coth; the negative root gives tanh–tanh. The latter is CZ's \(z_1\). At \(\theta=1.13,\eta=0.27,\tau=0.36\), the vacuum values are \(\pm0.6636366742\). In the tanh sector,

\[
e^{\lambda_1}=-\frac{(1+E_0)(1+z_1E_0)}{\sqrt2\cosh\tau(1-\rho E_0)(1-\widetilde\rho E_0)},\quad
\rho=\tanh\frac{\theta-\eta+\tau}{2},\quad
\widetilde\rho=\tanh\frac{\theta-\eta-\tau}{2}.
\]

This confirms both the overall sign and denominator signs in the chapter's correction to source (4.4). For imaginary \(E_0\) and real parameters the corrected phase formula (4.5) is

\[
\frac{e^{\lambda_1}}{e^{\bar\lambda_1}}=
\frac{1+E_0}{1-E_0}\frac{1+z_1E_0}{1-z_1E_0}
\frac{1+\rho E_0}{1-\rho E_0}\frac{1+\widetilde\rho E_0}{1-\widetilde\rho E_0}.
\]

The final two factors reverse those obtained from the erroneous printed (4.4). Source: [CZ 2009b](https://arxiv.org/abs/0908.3126).

For SG, \(\hat x=m_0x/\sqrt2\) and the prefactor \(m_0/(\sqrt2\beta_{\rm SG}^2)\) consistently convert the source equations. For Tzitzéica, differentiating the bulk potential gives \(V''(0)=3m^2\). With \(u=\beta\phi/\sqrt2\), \(\hat x=mx/\sqrt2\), its equation becomes \(u_{\hat t\hat t}-u_{\hat x\hat x}=-2e^{2u}+2e^{-u}\). The defect prefactor is consequently \(\sqrt2m/\beta^2\), as printed. The factor \(1/\hat\sigma\) must multiply both terms of \(f\); otherwise the momentum identity fails.

CZ (4.9)'s \(z_2\) equals the chapter's \(z_2\), whereas printed (4.11) is exactly the chapter's \(z_3\). The set of delays is unchanged. The asymptotic auxiliary-field values are the three cube roots stated in the chapter.

I checked Bristow–Bowcock's bulk exponent vectors in (4.19), (4.23), (4.27), (4.30) against the \(e_i-e_{i+1}\), \(e_i+e_{i+1}\), \(e_i\), and \(2e_i\) root presentations. In particular \(D_n\) has affine labels \(1,1,2,\ldots,2,1,1\), and \(C_n\) has \(1,2,\ldots,2,1\). These identify untwisted \(a_n^{(1)},b_n^{(1)},c_n^{(1)},d_n^{(1)}\), not their affine duals. [Bristow–Bowcock](https://arxiv.org/abs/1612.03002). The chapter's qualification is appropriate for their full defect constructions.

No substantive replacement is required. The corrected phase formula above can be inserted after 38-defect-conditions.qmd:356 to make the correction to (4.5) explicit.

## 6. The antilinear symmetry and the BCZ transpose

**Verdict: correct for fixed nondegenerate parameters, \(\sigma\in\mathbb C^\times\), and the stated action of \(\Theta\).**

**Independent derivation.** At imaginary coupling, all exponent vectors and \(D,E\) are real, while \(i\beta\) changes sign under conjugation. Thus

\[
\overline{f_j(\phi,\psi)}=f_j(-\bar\phi,-\bar\psi),\qquad
\overline{g_j(\phi,\psi)}=g_j(-\bar\phi,-\bar\psi).
\]

Conjugating \(\sigma\) proves the potential identity; conjugating its gradients proves the sewing identity. The \(f\)- and \(g\)-families have distinct joint exponent vectors, so equality for all fields forces \(\sigma=\bar\sigma\). This is a symmetry of the fixed theory, not merely a map to the theory with conjugated parameters. At \(m=0\) the “only if” inference is unavailable; the theorem concerns the massive model. classical_checks.py tests complex \(\sigma\) at ranks 1–4, both types, with residual below \(10^{-12}\).

For the BCZ transpose, let \(a=\lambda_1,b=\lambda_2\), with \(\lambda_0=\lambda_3=0\), in \(a_2\). Their (3.28) is \(2[(a-b)a^T+bb^T]\); (5.43) is \(2[a(a-b)^T+bb^T]\). These are unequal transposes because \(a,b\) are independent. Substituting their antisymmetric \(E\) in (5.42) gives the first expression for \(1-E\). The index argument extends to all ranks and confirms the chapter's correction. [BCZ 2004](https://arxiv.org/abs/hep-th/0401020).

For explicit scope, the first sentence at 38-defect-conditions.qmd:300 can read:

~~~markdown
**[Theorem]** For both types of massive $a_r^{(1)}$ defect at imaginary coupling, with $\beta\ne0$ and $\sigma\in\mathbb C^\times$, the following identity holds.
~~~

There is a dictionary sign error at 38-defect-conditions.qmd:65. Replace its last cell by:

~~~markdown
their $\gamma/2=b_{\rm SG}^2/(8\pi)$ at real coupling, continued to $-\beta_{\rm SG}^2/(8\pi)$ when $b_{\rm SG}=i\beta_{\rm SG}$
~~~

This follows by continuing the squared coupling in the real Euclidean normalization.

## 7. Universal R through cubic order, the two Borels, and BGKNR

**Verdict: correct for the universal-R direction and the stated expansion through \(O(x^3)\); wrong unqualified unitarity claim for complex \(s\).** A finite-order evaluation is not an independent proof of the complete universal product.

**Independent computation.** [quantum_checks.py](../referee_scratch/defects_38_40/quantum_checks.py), [quantum_output.txt](../referee_scratch/defects_38_40/quantum_output.txt), first solves the spin-half coproduct equations, independently of the KT formula. Put \(x=y_2/y_1\), normalize the highest-weight entry to one, and obtain

\[
R(x)=\begin{pmatrix}
1&0&0&0\\
0&\dfrac{q(1-x^2)}{q^2-x^2}&\dfrac{(q^2-1)x}{q^2-x^2}&0\\
0&\dfrac{(q^2-1)x}{q^2-x^2}&\dfrac{q(1-x^2)}{q^2-x^2}&0\\
0&0&0&1
\end{pmatrix}.
\]

Writing \(E=E_{12},F=E_{21}\), the separately constructed root vectors, with spectral powers stripped off, are

\[
e_\delta=EF-q^2FE=\operatorname{diag}(1,-q^2),\quad
f_\delta=EF-q^{-2}FE=\operatorname{diag}(1,-q^{-2}),
\]
\[
e_{\alpha+\delta}=-qE,\quad e_{2\delta-\alpha}=-qF,\quad
f_{\alpha+\delta}=-q^{-1}F,\quad f_{2\delta-\alpha}=-q^{-1}E.
\]

With \(K_c/q^{1/2}=\operatorname{diag}(1,q^{-1},q^{-1},1)\), multiply, in precisely the displayed order,

\[
K_c(1+\kappa xE\otimes F)
(1+\kappa x^3f_{2\delta-\alpha}\otimes e_{2\delta-\alpha})
\left(1+\frac{\kappa x^2}{[2]}f_\delta\otimes e_\delta\right)
(1+\kappa x^3f_{\alpha+\delta}\otimes e_{\alpha+\delta})
(1+\kappa xF\otimes E).
\]

Normalize the product's 00 entry and discard powers \(x^4\) and higher. Every coefficient of the difference from \(R(x)\) is **exactly zero**. Omitting the imaginary root produces an \(x^2\) error; omitting the height-three real roots produces an \(x^3\) error. These controls exclude a match from the first-order terms alone.

The final comparison in [final_checks.py](../referee_scratch/defects_38_40/final_checks.py), [final_output.txt](../referee_scratch/defects_38_40/final_output.txt), gives

~~~text
independent R versus P*boundary_common.rcheck:
q=exp(0.71i): relative error 1.52e-16, nullity 1
q=1.37:      relative error 1.04e-15, nullity 1
~~~

The Borel direction also follows without identifying defect types by name in an external source. The second leg of \(\mathcal R\) is positive Borel, and its Cartan factor on \(v_\pm\otimes|n\rangle\) is \(q^{\pm(n+s/2)}\). The first off-diagonal terms have \(+i\,y_d/y\). This is the book's I\(_+\) matrix. Conversely \(\mathcal R_{21}^{-1}\) has negative Borel in its second leg, inverse Cartan factor, and \(-i\,y/y_d\); it produces I\(_-\). The universal intertwining identity followed by \(P\) has exactly the spatially ordered coproduct in eq-borel-intertwiner. Thus the direction is correct in these fixed conventions.

BGKNR use \(\Delta e=e\otimes1+q^{-h}\otimes e\). Replacing their deformation parameter by its inverse, while consistently identifying the Cartan generators, gives the book's coproduct. Their positive-negative ordered R corresponds to the book's flipped inverse. Inverting a product reverses its order; flipping reverses its Borel legs; \(\kappa\) changes sign; and
\(\exp_{q^2}(-z)^{-1}=\exp_{q^{-2}}(z)\).
These changes explain the Cartan factor, ordering and \(q\)-exponential parameter in eq-dq-kt. A bare substitution in coefficients without these operations is insufficient. The independent cubic calculation confirms the translated formula, including its imaginary factor. [BGKNR, §§4–5](https://arxiv.org/abs/1205.1631).

Part (4) of prp-two-borels says “unitary ... for every \(s\),” although the preceding theorem allows \(s\in\mathbb C\). For \(|q|=1\), diagonal amplitudes \(q^{\pm(n+s/2)}\) need not have unit modulus for complex \(s\). [supplemental_checks.py](../referee_scratch/defects_38_40/supplemental_checks.py) normalizes the charge block by its scalar modulus and obtains:

~~~text
s=0.3:       ||L†L-1|| = 1.36e-16
s=0.3+0.2i:  ||L†L-1|| = 0.201999
~~~

**Paste-ready correction at part8-defects/40-defect-quantum-groups.qmd:132:** replace part (4) by

~~~markdown
(4) In the normalization of (1), the scalar factor is independent of the real offset $s$: $\Phi(\theta-\theta_d)\,q^{-s/2}$ times the matrix of (-@eq-sg-transmission-matrix) obeys physical unitarity for real rapidities, real $s$, and the physical $|q|=1$ parameters. The crossing identity extends algebraically to complex $s$; Hilbert-space unitarity does not.
~~~

## 8. The triangle theorem, its Hom hypothesis, and uniqueness

**Verdict: correct conditional theorem after requiring a nonzero normalization element; wrong uniqueness theorem with only \(q^4\ne1\); generic Hom dimension independently supported by exact calculations on the infinite modules.**

**Proof audit.** Coassociativity and \(\Delta(B)\subset B\otimes B\) make both triangle products intertwiners. If the relevant Hom space has dimension one, they are proportional. Agreement on a **nonzero** matrix element fixes the proportionality; agreement on a zero element does not. Multiplying \(T(\theta_1),T(\theta_2)\), or \(\check R\) by scalar functions multiplies both products by the same scalar. Thus condition (iii) is a compatibility check, not an adjustable normalization that can fix a failed triangle identity.

The proof needs closure of both coproduct legs. An antipode is not used in this argument. A single left coideal is insufficient for the displayed middle-leg argument; this does not prove that every conceivable transmitting-defect symmetry must be a Hopf subalgebra. More general actions or coactions are outside the theorem.

**Independent infinite-module computation.** A finite truncation of a bilateral module can create or remove intertwiners. I instead used finite necessary equations which every endomorphism of the *infinite* module must satisfy. In [hom_checks.py](../referee_scratch/defects_38_40/hom_checks.py), let the total charge be \(2l\). Each block of \(V_1\otimes V_2\otimes\mathcal W\) has four states, with defect index \(n=l-(a+b)/2\). Let \(A_l\) raise \(l\) by one and \(B_l\) lower it by one. A weight-preserving endomorphism satisfies

\[
X_{l+1}=A_lX_lA_l^{-1}.
\]

Each \(A_l\) is invertible: in particle-weight order its Weyl-shift contribution is nonzero diagonal and its particle-flip contribution is triangular. Transporting \(B_{l+1}A_l\) back to \(l=0\) gives finite commutator constraints on \(X_0\). At the exact rational point
\(q=2,y_1=3/2,y_2=5/4,y_d=7/6,s=0\),
three such loops give:

~~~text
Borel +: rank 15 of 16, commutant dimension 1; det A0=2401/6561
Borel -: rank 15 of 16, commutant dimension 1; det A0=256/2401
~~~

A nonzero 15-by-15 minor establishes the same rank on a generic parameter set. The scalar solution always exists. Propagation makes the whole infinite-module endomorphism scalar. Composing with an invertible generic product of local intertwiners identifies the desired Hom space with this endomorphism space, giving dimension one. No artificial upper or lower Weyl boundary was introduced.

For \(\mathcal F^+\), final_checks.py uses its genuine lowest total-charge block, of dimension one. Acting upward and with charge-preserving loops generates the next block with rank 3 and the following block with rank 4, exactly at the same rational point. All higher raising blocks are triangular with nonzero diagonal. Thus the lowest vector is generically cyclic; an endomorphism is scalar on it and hence everywhere. Output:

~~~text
Fock+ lowest-weight cyclicity: block0 rank 3; block1 rank 4
det A1 = 2401/6561
~~~

This verifies the generic *Hom-dimension* hypothesis more strongly than a window test. It does not assert irreducibility: a scalar endomorphism ring is not, in general, a proof that a module has no proper submodule.

**Counterexample to prp-sg-t-unique.** Take \(q=e^{i\pi/3}\), so \(q^4\ne1\) but \(q^6=1\). On a bilateral Weyl module let \(S|n\rangle=|n+1\rangle\). Then \(S^3\) commutes with \(q^{\mathsf h}\) and with both constant shifts \(f_0,f_1\). Therefore

\[
(S^3\otimes1)\mathbb T
\]

is another intertwiner, linearly independent of \(\mathbb T\). Arbitrarily many powers give further ones. The failure occurs at the start of the proof: the \(k\)-equations no longer force the four charge channels written there, because they distinguish weights only modulo the root-of-unity period. The homogeneous calculation \(\alpha(q^{-3}-q)=0\) is correct *after* imposing those four channels, but does not exclude the extra channels. quantum_checks.py explicitly checks the commuting shift; its \([S^3,k]\) numerical residual is below \(10^{-15}\).

For \(q\) not a root of unity, the four-channel argument is valid and the printed recurrence proves uniqueness for \(y,y_d\ne0\). The “for all \(y,y_d\)” wording is also too broad: the representation uses their inverses.

The level-zero conclusion can be checked directly, without a numerical window. If \(k_0k_1=c\), use charge block \((v_+\otimes|n\rangle,v_-\otimes|n+1\rangle)\). The \(f_1\) equations give \(B_{n+1}=qB_n\), \(C_{n+1}=q^{-1}C_n\), and \(D_n=q^{2n+s}A_n\). The \(f_0\) equations give instead \(D_n=c^{-1}q^{2n+s}A_n\). For \(c\ne1\), \(A=D=0\), and the remaining equations force \(B=C=0\). This assumes the nonzero constant Weyl shifts and the generic charge-channel condition. final_checks.py independently solves necessary finite equations at \(c=1,-1,4,1/4\), obtaining nullities \(1,0,0,0\).

**Paste-ready corrections.** At 40-defect-quantum-groups.qmd:103 replace (iii) by:

~~~markdown
(iii) the two composite intertwiners agree on a matrix element that is nonzero for generic parameters.
~~~

After the proof at line 107 insert:

~~~markdown
Condition (iii) is an actual check: scalar rescalings of the individual intertwiners multiply the two sides by the same factor. The coproduct closure used here is sufficient for this proof; no necessity theorem for all possible defect symmetry structures is being asserted.
~~~

At line 114 replace the opening of prp-sg-t-unique by:

~~~markdown
**[Theorem]** Suppose $q\in\mathbb C^\times$ is not a root of unity and $y,y_d\in\mathbb C^\times$. For each offset $s$, with a fixed choice of $q^s$, the algebraic $U_q(\mathfrak b^-)$-intertwiners $V(y)\otimes\mathcal W^-(y_d)\to\mathcal W^-(y_d)\otimes V(y)$ form a one-dimensional space.
~~~

Retain the matrix formulas and qualify the level-zero statement by the same assumptions. Alternatively retain roots of unity only by explicitly restricting to intertwiners preserving the unexponentiated integer total charge; that is a different theorem.

At 40-defect-quantum-groups.qmd:65 replace the unqualified necessity claim by:

~~~markdown
The Borel Hopf subalgebras provide a sufficient setting: both coproduct legs remain in the algebra represented on the defect, so the same action is defined when the defect occupies a middle tensor position. A single left-coideal action does not provide this step of the triangle proof. This argument establishes the setting used here, not a classification of all possible transmitting-defect symmetries.
~~~

At line 71 add a scope sentence before introducing the oscillator modules:

~~~markdown
Unless a root-of-unity specialization is stated explicitly, assume $q\in\mathbb C^\times$ is not a root of unity and fix the powers $q^s$. Reducibility and additional intertwining channels at roots of unity require separate treatment.
~~~

## 9. Inverse L-operators, Weyl and Fock traces, and chapter 36

**Verdict: correct inverse orientation and \(\mathcal W^-\) factorization; wrong if that factorization is stated for both Weyl modules; correct absence of ordinary bilateral convergence; unclear analytic justification for every claimed infinite-auxiliary commutation.**

**Independent computation and derivation.** The local relation
\(R_{aj}X L_j=L_jXR_{aj}\) implies
\[
XR_{aj}L_j^{-1}=L_j^{-1}R_{aj}X.
\]
This is the train relation required by the chosen order \(R_{aN}\cdots R_{a1}\). It fixes the inverse, rather than allowing its choice by analogy. The scalar-free inverse represents the matrix part of right-to-left transmission; physical transmission also includes its scalar factor.

quantum_checks.py constructs \(\mathcal D\) by enumerating the finitely many exact auxiliary paths at each site and selecting paths returning to the starting \(n\). There is no Weyl cutoff. It checks every charge sector at \(N=2,\ldots,6\), two generic unit-circle values of \(q\), nonzero offset, and unequal rapidities. Maximum relative residuals are:

| Check | Maximum |
|---|---:|
| \(\langle2|M|2\rangle=q^{2m}\langle0|M|0\rangle\) | \(8.3\times10^{-17}\) |
| \([t_{q^{-m/2}},\mathcal D]\) | \(6.6\times10^{-17}\) |
| \([\mathcal D(y_d),\mathcal D(y'_d)]\) | \(1.9\times10^{-17}\) |
| Wrong-twist control in nontrivial sectors | 0.057–0.137 |
| Replacing inverse L by L, \(N\ge3\) | 0.317–0.612 |

For I\(_-\), shifting every auxiliary index by one multiplies a column of \(L^-\) by \(q^{-a}\), \(a=\pm1\); it multiplies a row of \((L^-)^{-1}\) by \(q^a\). On a closed path their product is \(q^m\). Thus, exactly,

\[
\langle n|M^-|n\rangle=q^{mn}\langle0|M^-|0\rangle,\quad
\operatorname{Tr}_{\mathcal W^-}(\omega^{\mathsf h}M^-)
=\omega^s\sum_{n\in\mathbb Z}(\omega^2q^m)^n\,\mathcal D.
\]

No nonzero complex value of \(\omega^2q^m\) makes both tails of this geometric series converge; at modulus one its terms do not tend to zero. The zero operator in a special sector is a trivial exception, not a definition of the generic trace.

For \(\mathcal W^+\), the inverse instead gains the column factor \(q^{-a}\), giving

\[
\langle n|M^+|n\rangle=q^{-mn}\langle0|M^+|0\rangle.
\]

final_checks.py verifies this already at \(N=1\). For \(q=e^{0.71i}\), the \(v_+\) diagonal ratio is \(0.758361876-0.651833771i=q^{-1}\), whereas the unqualified proposition would give \(q\). Thus “for a Weyl module” must be restricted to \(\mathcal W^-\), or both signs supplied. The formal support for \(\mathcal W^+\) is \(\omega^2=q^{m}\).

The series \(\sum_{n\in\mathbb Z}r^n\) can be treated as a formal multiplicative delta, or as a periodic distribution when \(r\) lies on the unit circle. It is not an analytic function supported at a point of \(\mathbb C^\times\). For \(r=e^{i\chi}\) its symmetric partial sum is
\[
\sum_{n=-K}^K r^n=\frac{\sin((K+\tfrac12)\chi)}{\sin(\chi/2)}.
\]
It converges distributionally to the periodic delta (with the usual normalization), not pointwise to zero off its support. The formal trace argument needs this interpretation and a compatible twisted cyclic functional; ordinary cyclicity cannot be assigned to a divergent trace.

The brief asks whether commutation holds “only” at the given twist. That is not literally true, and the chapter itself more carefully says “for.” Since \(t_{-\omega}=-t_\omega\), both square roots give commutation. For \(N=1\) every sector is one-dimensional and every twist commutes; extremal sectors at any \(N\) have the same exception. The independent \(N=1\) commutator is exactly zero for an arbitrary test twist. Keep a sufficient sector-dependent twist statement, not an iff claim.

For \(\mathcal F^+\), at \(|q|=1\), \(|\omega|<1\), and away from local spectral poles, the explicit coefficients are bounded in \(n\), so the weighted trace converges absolutely. At \(\omega=0.63,N=3\), increasing the genuine Fock cutoff from 45 to 65 changes the trace by less than numerical resolution; \([t_\omega,Q]\) has residual at most \(1.86\times10^{-16}\). Such transfer operators are Q-operators. The arbitrary-complex-\(q\) claim does not follow from \(|\omega|<1\) alone. Also, a formal auxiliary-auxiliary intertwiner does not automatically supply the boundedness and trace manipulations needed for an analytic proof of \([Q,Q']=0\); the chapter's sketch label should remain.

**Chapter 36 consistency.** The rotated-channel scalar transmission at \(i\pi/2-\theta\) and the finite-chain algebraic operator are related constructions, not literally identical formulas before orientation, crossing and scalar normalizations are specified. The unrestricted assertion at 36-defect-scattering.qmd:225 that infinite Borel-module traces are Q-operators is too broad for the bilateral Weyl modules now discussed.

**Paste-ready corrections.** At 40-defect-quantum-groups.qmd:205 change part (1)'s domain to:

~~~markdown
(1) For the Fock modules in the physical regime $|q|=1$, with $|\omega|<1$ and away from spectral poles, the weighted trace converges absolutely. The train argument gives commutation with $t_\omega$ when the required auxiliary intertwining and trace manipulations are defined. Mutual commutation uses the auxiliary-auxiliary intertwiner; we retain the sketch qualification for its analytic implementation.
~~~

Change the beginning of part (2) to:

~~~markdown
(2) For the bilateral Weyl module $\mathcal W^-$, the generic trace has no ordinary convergent twist. Its sector factorization is $\langle n|M|n\rangle=q^{mn}\langle0|M|0\rangle$. The series $\sum_{n\in\mathbb Z}(\omega^2q^m)^n$ is a formal multiplicative delta, or a distribution on the unit circle. For $\mathcal W^+$ the exponent is $-mn$ and the formal support is instead $\omega^2=q^m$.
~~~

At line 213 replace “the average over \(n\) is then a cyclic trace” by:

~~~markdown
The bilateral expression must be interpreted as a formal twisted trace. On its sector support the auxiliary shift factors cancel; this is the algebraic content of the train argument. It does not define an ordinary convergent trace. The finite matrix element $\mathcal D$ is well defined separately.
~~~

At line 201 replace the physical-inverse sentence by:

~~~markdown
With scalar factors restored, defect unitarity identifies this orientation with right-to-left transmission. In the scalar-free normalization used here, $\mathbb L^{-1}$ is its matrix part.
~~~

At part8-defects/36-defect-scattering.qmd:225 replace the general trace/Q assertion with:

~~~markdown
For suitable auxiliary modules and a convergent twisted trace, this construction gives a commuting transfer operator. In the conventions of @sec-q-operators the product uses inverse L-operators. Fock auxiliary modules give Q-operators; bilateral Weyl modules require the separate charge-sector construction of @prp-defect-transfer, since their ordinary traces do not converge.
~~~

## 10. The Weyl TQ relation and its bootstrap interpretation

**Verdict: correct in the independent finite-chain tests through \(N=6\); unclear as an all-\(N\) theorem or a fully normalized physical bootstrap identity.** The conjecture label should remain.

**Independent computation.** quantum_checks.py uses the independently solved R and exact closed Weyl paths, without fitted coefficients. At \(y=i y_d\), it tests precisely

\[
t_{q^{-m/2}}(i y_d)\mathcal D(y_d)
=q^{-N/2}\mathcal D(qy_d)
+\left[\prod_j\frac{q^{1/2}(qz_j^2+q^{-1})}{1+z_j^2}\right]\mathcal D(y_d/q)
\]

in every sector. Maximum relative errors at \(N=2,3,4,5,6\) for \(q=e^{0.71i}\) are respectively
\(1.1\times10^{-16},7.3\times10^{-17},7.6\times10^{-17},1.1\times10^{-16},1.3\times10^{-16}\).
For the physical choice \(q=e^{-0.37\pi i}\), corresponding to \(\lambda=1.37\) on the chosen branch, they are
\(2.5\times10^{-16},1.5\times10^{-16},2.6\times10^{-16},3.0\times10^{-16},6.6\times10^{-16}\).
The offset is \(s=0.3\) and the site rapidities are unequal.

The same independent code verifies the printed Fock relation, including the \(s\mapsto s\pm1\) shifts and its different prefactors, using a convergent trace. The residual is at most \(5.30\times10^{-16}\). This also checks that the Weyl formula was not obtained by simply forgetting the Fock offsets.

There is an elementary independent check of the bootstrap parameter reading. Every closed Weyl path has an even number of off-diagonal steps, so \(\mathcal D(-y_d)=\mathcal D(y_d)\). Changing \(y\to-y\) in \(R\) is conjugation by \(\sigma_z\) on the auxiliary particle; it leaves the diagonally twisted trace invariant. Finally
\[
q^{\pm1}y_d=-e^{\lambda(\theta_d\mp i\pi)},\qquad
e^{\lambda(\theta_d\pm i\pi/(2\lambda))}=\pm i y_d.
\]
These identities justify the two odd-sector parameters and either resonance sign used in the chapter.

They do not alone identify an image module, a kernel module, or the physical scalar residue. A short exact sequence with explicit embedding/projection, followed by the local fusion identity in this inverse-L normalization, would prove the all-length identity and settle those assignments. Weston's §2.4.1 supplies the corresponding mechanism, parameter shifts and explicit maps in his conventions; this is strong prior support, not a proof that every normalization in this chapter has already been translated. [Weston 2010](https://arxiv.org/abs/1006.1555).

**Paste-ready correction at 40-defect-quantum-groups.qmd:235:** replace the final bootstrap conclusion with:

~~~markdown
The spectral shifts agree with the odd-sector shifts in the defect bootstrap. This supports interpreting the TQ relation as its scalar-free transfer-operator counterpart. To identify the two terms with the image and kernel of the physical resonance fusion, one must translate the exact-sequence maps and restore the scalar transmission factors; the parameter identities alone do not establish those identifications.
~~~

At line 226 the numerical range can be updated to \(N\le6\) if these independent tests are adopted, but it must not be relabelled a theorem on that basis.

## 11. The classical limit and the two different deformation parameters

**Verdict: correct, conditional on the classical defect-coordinate Poisson algebra of item 3.**

**Independent derivation.** In the \(a_1\) defining representation, write
\[
u=e^{-b(\phi-\psi)/(2\sqrt2)},\quad w_\pm=e^{\pm b\phi/\sqrt2}.
\]
From \(\{\phi,\psi\}=-1\),
\[
\{\log u,\log w_\pm\}=\pm b^2/4,\qquad
\{u,w_\pm\}=\pm(b^2/4)uw_\pm.
\]
supplemental_checks.py calculates these brackets symbolically, returning exactly \(+b^2/4,-b^2/4\). Multiplying the three classical Lax factors gives
\[
L_d=\operatorname{diag}(u^{-1},u)
\left[1+\frac{\sigma}{\lambda}(E_{21}w_++E_{12}w_-)\right].
\]
Indeed its upper off-diagonal exponential is \(e^{-b(\phi+\psi)/(2\sqrt2)}\), and its lower one is the reciprocal, as required.

On the Weyl module take \(U|n\rangle=q^{n+s/2}|n\rangle\), and
\[
W_+|n\rangle=q^{-2n-1-s}|n+1\rangle,\qquad
W_-|n\rangle=q^{2n-1+s}|n-1\rangle.
\]
Then \(W_+W_-=W_-W_+=1\), \(UW_\pm=q^{\pm1}W_\pm U\), and multiplying the factors in eq-dq-weyl-l reproduces the entries of \(L^-\). The independent numerical residual is below \(3.34\times10^{-16}\). Operator ordering matters in the off-diagonal terms.

Quantization with \(q_{\rm local}=e^{i\hbar b^2/4}\) gives the displayed Poisson bracket on dividing the commutator by \(i\hbar\). Since \(b=i\beta=i\beta_{\rm SG}/\sqrt2\) in the \(a_1\) convention, this is \(\hat q=e^{-i\hbar\beta_{\rm SG}^2/8}\). This local deformation parameter tends to one with \(\hbar\), whereas the physical S-matrix parameter \(q=-e^{-i\pi\lambda}=e^{-8\pi^2i/\beta_{\rm SG}^2}\) has inverse-coupling dependence. They cannot be identified by a naive \(q\to1\) substitution at fixed physical coupling.

The classical dictionary and the warning that \(\hat q\) is not the S-matrix \(q\) are correct. No replacement is needed apart from preserving the conditional status of the underlying Poisson bracket.

## 12. CFT lines, the massless limit, and chiral perturbations

**Verdict: correct local gluing, chiral decomposition, and Ising/Lee–Yang eigenvalues; wrong unrestricted compact-boson interpretation and Schur-lemma hypothesis; unclear as an all-orders derivation of the massive Borel charges from first-order perturbation theory.**

**Independent computation.** In the massless equations write
\(\phi_t=A+B,\phi_x=A-B,\psi_t=C+D,\psi_x=C-D\).
The two sewing equations give \(C=A,D=-B\). Integrating yields \(\psi_L=\phi_L+c_L,\psi_R=-\phi_R+c_R\); one constant can be absorbed in the decomposition. Consequently \(\phi+\psi=2\phi_L+\text{const}\), \(\phi-\psi=2\phi_R+\text{const}\). The two cosine terms have opposite chirality at this conformal gluing, with classical coefficients \(m_0\sigma\), \(m_0\sigma^{-1}\). Their product is proportional to \(m_0^2\); their opposite boost dependence agrees with the classical rapidity parameter. This does not determine a renormalized quantum relation between a Lagrangian parameter and \(\theta_d\).

For a compact boson, local stress-tensor gluing is insufficient to define an unrestricted defect of the same theory. In units where
\[
p_L=n/R+wR,\quad p_R=n/R-wR,
\]
the reflection \(p_R\mapsto-p_R,p_L\mapsto p_L\) would send \((n,w)\) to \((R^2w,n/R^2)\). At \(R^2=2\), the state \((1,0)\) would have \(w'=1/2\), outside the lattice. The gluing naturally defines the T-duality interface between radii \(R\) and \(1/R\); a same-theory realization requires self-duality or suitable lattice restrictions and possibly a noninvertible interface. supplemental_checks.py records this counterexample.

For the Ising theory the independently evaluated matrix is
\[
S=\frac12\begin{pmatrix}1&1&\sqrt2\\1&1&-\sqrt2\\\sqrt2&-\sqrt2&0\end{pmatrix}.
\]
The ratios \(S_{ai}/S_{0i}\) give
\[
D_\varepsilon=(1,1,-1),\qquad D_\sigma=(\sqrt2,-\sqrt2,0),
\]
and \(D_\sigma^2=1+D_\varepsilon\) exactly. For Lee–Yang I used
\[
S=\frac2{\sqrt5}\begin{pmatrix}-\sin(2\pi/5)&\sin(\pi/5)\\
\sin(\pi/5)&\sin(2\pi/5)\end{pmatrix}.
\]
The two ratios are \((1-\sqrt5)/2,(1+\sqrt5)/2\), and each satisfies \(d^2-d-1=0\) exactly. Outputs are in [supplemental_output.txt](../referee_scratch/defects_38_40/supplemental_output.txt). The eigenvalues and line assignments in the body text are correct.

The Schur argument at line 23 requires commuting with the **full chiral algebra** whose irreducible modules are \(R_i\). Commutation with Virasoro alone does not suffice when those modules are reducible under Virasoro. For diagonal Virasoro minimal models, including these two examples, the stated argument works.

The first-order current discontinuity does show why one chirality can remain continuous. It is not by itself an all-orders construction of defect contributions to the massive nonlocal quantum charges. The exact algebraic statement about the chosen Weyl module is stronger and can be proved without a least-squares residual: for generic \(q\), an extension \(e_1\) must have the form \(e_1|n\rangle=c_n|n+1\rangle\). The required relation \([e_1,f_0]=0\), with \(f_0\) a nonzero constant upward shift, forces \(c_{n+1}=c_n\). Then \(e_1\) commutes with the constant downward shift \(f_1\), contradicting
\[
[e_1,f_1]=(k_1-k_1^{-1})/(q-q^{-1}).
\]
Thus this generic Weyl Borel representation does not extend to the full quantum affine algebra. That is not a classification excluding every other conceivable defect symmetry.

**Paste-ready corrections.** At part8-defects/40-defect-quantum-groups.qmd:23 replace “commutes with both Virasoro algebras” by “commutes with both chiral algebras defining the irreducible modules \(R_i,\bar R_i\).”

After line 37 insert:

~~~markdown
These are local gluing equations. For a compact boson, a defect must also respect the momentum–winding lattice. Reflection of one chirality is naturally a T-duality interface; realizing it as a line within the same compactified theory requires compatible radii and lattice data. We use the local gluing here to identify the chiral fields.
~~~

At line 61 replace the numerical no-extension claim and final universal-symmetry sentence by:

~~~markdown
For generic $q$, this Weyl representation cannot extend to all of $U_q(\widehat{sl}_2)$: the weight relation forces $e_1|n\rangle=c_n|n+1\rangle$, and $[e_1,f_0]=0$ forces $c_n$ to be constant, contradicting $[e_1,f_1]=(k_1-k_1^{-1})/\kappa$. The exact transmission matrices therefore realize the indicated Borel symmetry. The first-order conformal argument motivates this result but does not itself construct the all-orders massive conserved charges.
~~~

## 13. Weston, transposition, bilateral traces, and BCZ's labels

**Verdict: correct reconciliation once reversal of intertwining direction and tensor order are made explicit; correct bilateral-trace caveat; correct need to relabel the defect when conjugating the transmitted particle.** The trace issue is a qualification of a formal construction, not evidence that Weston asserted an ordinarily convergent bilateral trace.

**Independent derivation.** For \(\tau(f)=k^{-1}e,\tau(k)=k\),
\[
\Delta(k^{-1}e)
=(k^{-1}\otimes k^{-1})(e\otimes1+k\otimes e)
=k^{-1}e\otimes k^{-1}+1\otimes k^{-1}e
=(\tau\otimes\tau)\Delta(f).
\]
The anti-homomorphism property turns transposes of positive-Borel representations into negative-Borel representations. supplemental_checks.py verifies this coproduct identity exactly. But transposition gives
\[
L\Delta'=\Delta'_{\rm op}L
\quad\Longrightarrow\quad
L^T\Delta_{\rm op}=\Delta L^T,
\]
which reverses the intertwining direction. To return to a convention \(L_{\rm new}\Delta=\Delta_{\rm op}L_{\rm new}\), one must also reverse the tensor order (or take the corresponding inverse). Weston's operator is on \(\mathcal W\otimes V\), whereas the book evaluates on \(V\otimes\mathcal W\). This is exactly the relevant reversal, not an arbitrary identification of a transpose with the same oriented matrix.

Weston's (4.15) and footnote 10 explicitly flag transposed component conventions. With \(j=-n,r_0=q^s,\zeta=i\,y/y_d\), the charge-block ratios become those of I\(_-\): \(A_{n+1}/A_n=q^{-1}\), \(D_n/A_n=q^{2n+s}\), and \(C_nB_{n+1}/(A_nD_{n+1})=-(y/y_d)^2\), after the indicated diagonal gauge. These ratios also follow directly from the independently determined negative-Borel matrix in item 8. [Weston 2010](https://arxiv.org/abs/1006.1555).

For \(r_1=r_2=0\), the absence of a highest or lowest auxiliary weight produces the bilateral geometric factor derived in item 9. At \(|q|=1\) there is no ordinary twist that makes it converge. Weston's language is a *regularized* trace; it leaves the regularization to be specified. Rossi–Weston had already isolated this delta factor explicitly, including its singular neutral block; see item 17. It would be inaccurate to present that issue as a newly discovered error in the existence of formal Q-operators.

**Independent two-label calculation.** Let the old labels be \(n=(b-a)/2,m=(a+b)/2\), and put \(Y=y/y_d\). Lift the single-jump amplitudes to their vacuum transitions:
\[
\begin{array}{c|c|c}
\text{input}&(\text{new }n,\text{new }m;\text{output})&\text{amplitude}\\ \hline
+&(n,m-1;+)&q^{-n}\\
+&(n+1,m;-)&-iYq^{-n}\\
-&(n-1,m;+)&-iYq^{n+s}\\
-&(n,m+1;-)&q^{n+s}
\end{array}
\]
The map \(b\mapsto-b\) sends \((n,m)\mapsto(-m,-n)\); conjugate the outgoing particle as well. Now apply the diagonal basis factor \(g(n,m)=q^{m(m+1)/2-nm}\). For example the new conversion amplitude from the old \(q^m\) is multiplied by \(g(n+1,m)/g(n,m)=q^{-m}\). The four amplitudes become, exactly,
\[
A'_n=-iYq^n,\quad C'_n=1,\quad B'_n=q^s,\quad
D'_n=-iYq^{s+1-n}.
\]
They are independent of the second label. Therefore
\[
A'_{n+1}/A'_n=q,\quad D'_n/A'_n=q^{-2n-s_+},\quad
\frac{C'_nB'_{n+1}}{A'_nD'_{n+1}}=-Y^{-2},\quad s_+=-1-s.
\]
These are precisely the positive-Borel gauge invariants. This calculation verifies the offset change and the switch to \(y_d/y\). Merely conjugating the outgoing soliton while retaining the old single defect label violates total-charge conservation. BCZ (3.57)–(3.59) describe the necessary two-label freedom; the book's clarification is justified. [BCZ 2005](https://arxiv.org/abs/hep-th/0506169).

**Paste-ready correction to the exercise at 40-defect-quantum-groups.qmd:306:**

~~~markdown
Show that $\tau(f_j)=k_j^{-1}e_j$, $\tau(k_j)=k_j$ is an algebra anti-isomorphism with $(\tau\otimes\tau)\Delta=\Delta\tau$. Show that transposing $L\Delta=\Delta^{\rm op}L$ gives $L^{\mathsf T}\Delta^{\rm op}=\Delta L^{\mathsf T}$ on the transposed modules. Account explicitly for the reversal from Weston's $\mathcal W\otimes V$ to the book's $V\otimes\mathcal W$ before identifying the negative-Borel L-operator.
~~~

The same clarification can be appended to line 136's proof. Its displayed reversed equation is already correct.

## 14. The rapidity limits and the scalar phase

**Verdict: correct limiting transmission matrix and twist angle for real physical parameters; wrong statement that the Gamma ratio itself tends to one.** The large-positive-rapidity limit includes the defect-charge shift and is not a particle-only operator on a fixed-\(n\) subspace.

**Independent scalar computation.** Instead of truncating the book's Gamma product, I summed the standard integral representation of its logarithm. For \(\delta=\theta-\theta_d\), \(\lambda>0\), the scalar in the book's minimal normalization can be written
\[
\Phi(\delta)=\frac{\exp\left\{i\pi(1+\lambda)/4+
2i\int_0^\infty\frac{\sin(\lambda\delta v/2\pi)}{v}H(v)\,dv\right\}}
{1+i e^{\lambda\delta}},
\]
\[
H(v)=\frac{e^{-v/4}-e^{-(\lambda/2+3/4)v}}
{(1+e^{-\lambda v/2})(1-e^{-v})},\qquad H(0)=\frac{1+\lambda}{4}.
\]
The differences of the linear Gamma terms cancel before summation. The small-\(v\) limit and the Dirichlet integral give
\[
\Phi(-\infty)=1,\qquad
e^{\lambda\delta}\Phi(\delta)\longrightarrow e^{i\pi\lambda/2}
\quad(\delta\to+\infty).
\]
In particular the Gamma ratio, excluding the explicit prefactor, tends to \(e^{-i\pi(1+\lambda)/4}\), not one, at the negative end.

supplemental_checks.py integrates this expression independently. For \(\lambda=0.73,1.37,2.6\), the errors \(|\Phi(-10)-1|\) are \(5.41\times10^{-5},4.98\times10^{-5},7.99\times10^{-5}\). The correspondingly normalized positive limits have the same errors. For \(\lambda=1.37\), a full charge-block comparison with the proposed diagonal limit falls from 0.0038653 at \(\delta=-6\) to \(7.0464\times10^{-5}\) at \(\delta=-10\). This tests the phase as well as the modulus.

The limiting diagonal entries are \(q^{-n},q^{n+s}\), hence
\[
q^{s/2}\operatorname{diag}\big(q^{-(2n+s)/2},q^{(2n+s)/2}\big).
\]
With \(q^{1/2}=e^{i\pi(1-\lambda)/2}\), the twist angle is
\(\alpha_n=-\pi(1-\lambda)(2n+s)/2\bmod2\pi\).
This is real for real \(s\). The common \(q^{s/2}\) reflects the chapter-37 normalization; in the normalized unitary matrix of item 7 it is removed.

For \(\mathcal D\), the auxiliary state is fixed at \(n=0\), and the inverse normalized diagonal matrix contributes \(q^{as/2}\) per particle. Therefore \(\mathcal D\to q^{ms/2}\), exactly as stated. At the other end, the surviving amplitudes are proportional to \(-iq^{-n}|n+1\rangle\otimes v_-\) and \(-iq^{n+s}|n-1\rangle\otimes v_+\). They exchange particle charge and shift the defect simultaneously. Projecting to a fixed defect-charge state would remove them, not give a charge-conjugation matrix there.

**Paste-ready corrections.** At 40-defect-quantum-groups.qmd:239 specify \(\lambda>0\), real rapidities and real \(s\). At line 247 replace the Gamma-ratio sentence with:

~~~markdown
The scalar satisfies $\Phi(\delta)\to1$ as $\delta\to-\infty$. The Gamma-function ratio tends to $e^{-i\pi(1+\lambda)/4}$, cancelling the explicit prefactor $e^{i\pi(1+\lambda)/4}$; its modulus alone would not determine this limit.
~~~

At line 252 replace “Far below, it is the charge-conjugation line” with:

~~~markdown
Far below, it converts the particle into its conjugate while shifting the defect charge, as in part (3); this limit acts on the combined particle–defect space.
~~~

## 15. Status labels and the comparison table

**Verdict: mixed. The conditional triangle theorem and fixed-parameter \(\Theta\) theorem are justified; the unrestricted uniqueness theorem is wrong; several “Established” or universal descriptions exceed the evidence.**

This status audit uses the exact counterexample and infinite-module rank computations in item 8, the exact \(a_2\) flux identities in item 1, and the direct bracket calculation in item 3. For example \(q=e^{i\pi/3}\) satisfies the printed theorem's hypothesis while admitting \((S^{3k}\otimes1)T\) for every integer \(k\); no numerical generic check can repair that quantified theorem. Conversely, exact vanishing of the \(a_2\) flux minus a total derivative establishes a symbolic identity and deserves a stronger label than a floating-point check.

| Object | Appropriate status and scope |
|---|---|
| prp-defect-momentum | Established within the stated first-order ansatz; energy/momentum algebra is reproducible. |
| prp-ar-only | Established under the restrictive ansatz; rank-four enumeration is corroboration, not the all-rank proof. |
| prp-defect-lax | Established for the known family, with mirror/relabel equivalence included. |
| prp-spin-two | Exact obstruction cancellation for all ranks; exact local \(a_2\) result; numerical \(a_3\); formal/generic sketch for the full hierarchy. |
| prp-defect-involution | Conditional Poisson-algebra proposition, with explicit independence assumption; full equal-time reduction not established here. |
| prp-defect-theta | Theorem for nonzero massive parameters and the fixed action of \(\Theta\). |
| thm-t-triangle | Theorem under all three hypotheses, including a nonzero comparison element. |
| prp-sg-t-unique | Theorem after generic-\(q\), nonzero-spectral-parameter correction; false as currently quantified. |
| prp-two-borels | Algebraic direction established generically; physical unitarity requires the real parameter domain. |
| prp-l-operator | Exact formulas supported by intertwining and KT checks; window uniqueness is not an infinite-module proof by itself. |
| prp-defect-transfer | Exact sector factorization with a Borel sign; formal trace/sketch qualifications and analytic domains required. |
| cnj-defect-tq | Independently checked through \(N=6\); keep the all-\(N\) conjecture in these conventions unless a fusion proof is supplied. |
| prp-defect-twist-limit | Correct after fixing the scalar-phase explanation and real-parameter domain. |

“Symbolically to \(10^{-15}\)” at 38-defect-conditions.qmd:167 conflates exact simplification with numerical evaluation of a symbolic expression. The supplied independent calculation gives exact zeros for \(a_2\); the existing book script's finite residual should be described as a numerical evaluation unless its expression has also been proved zero.

**Paste-ready replacement at 38-defect-conditions.qmd:155:**

~~~markdown
**[Symbolic for $a_2$; Numerical for $a_3$; Proposition, sketch for general rank]** For $a_r^{(1)}$ with $r\ge2$, the claimed obstruction to conserving both spin-$\pm2$ charges vanishes. Both local defect corrections exist explicitly for $a_2^{(1)}$, with numerical evidence for $a_3^{(1)}$. More generally, the formal expansion at the two spectral ends supports separate defect-corrected charges of both chiralities for each spin present in the bulk hierarchy.
~~~

For tbl-boundary-defect-algebra at 40-defect-quantum-groups.qmd:260–269, the following cell replacements avoid overstating the conclusions:

~~~markdown
Symmetry status: Borel intertwining established for the stated generic modules; no extension of this Weyl module to the full quantum affine algebra. A classification of all defect symmetries is not asserted.

Quadratic equation status: [Theorem] under hypotheses (i)–(iii).

Quadratic algebra, defect cell: Yang–Baxter (RTT) algebra, with $\mathbb L^+=(\pi\otimes\rho^+)(\mathcal R)$ and $\mathbb L^-=(\pi\otimes\rho^-)(\mathcal R_{21}^{-1})$, up to scalar factors.

Operator on the circle, defect cell: Fock trace Q-operator; for bilateral Weyl modules, the separate sector matrix element of @prp-defect-transfer.

Classical-limit status: Conditional Poisson-algebra construction; local bracket checked directly. Identification with the full equal-time field theory requires a separate reduction.
~~~

The introductory learning goal that “every intertwiner” solves the triangle relation should likewise say “under the Hom-dimension and normalization hypotheses.” An intertwiner alone does not imply those hypotheses.

## 16. Exercises: independent spot-checks

**Verdict: at least three exercises in each chapter are solvable and correct after the stated qualifications; the cubic-tensor normalization and Ising disorder-field question require correction.** I checked at least six from each chapter.

**Chapter 38.**

1. **exr-dc-energy, line 316 — correct with the stated local modified-energy ansatz.** For velocities \(v=\phi_t,w=\psi_t\), the trial conditions give
   \[
   \dot E=v^TAv+v^T(B-C^T)w-w^TDw+v\cdot X-w\cdot Y.
   \]
   Cancellation with \(d\mathcal B/dt\) for arbitrary velocities forces antisymmetric \(A,D\), \(C=B^T\), \(X=-\nabla_\phi\mathcal B\), \(Y=\nabla_\psi\mathcal B\). This independently solves the converse and fixes the gradient signs.
2. **exr-dc-d-matrix, line 320 — correct.** On the sum-zero space \(D=2(1-\pi)^{-1}\), so \(D+D^T=2\), and on a \(\pi\)-eigenvector \(D=2/(1-e^{2\pi ik/h})=1+i\cot(\pi k/h)\). Conjugating \(\pi\) to \(\pi^{-1}\) gives \(D^T\). The root/weight formula follows by multiplying by \(1-\pi\), and the mirror-species delay then follows from item 4.
3. **exr-dc-scan, line 324 — correct for the indicated diagrams.** A proper orbit fixes its nonzero root sum; \(c_3^{(1)}\)'s endpoint pair is an example. The \(a_4\) matrix and ten forbidden entries are given in item 2; \(a_1\)'s single spatial direction gives \(D=1\). These checks are in classical_checks.py.
4. **exr-dc-cot, line 328 — wrong tensor coefficient, otherwise correct.** The sine-addition proof, exact nonzero \(C(Du)\), and cancellation are in item 1. Replace \(2/3\) by \(2\) when naming BCZ's tensor \(A\), or retain \(2/3\) and name the density polynomial \(C=A/3\).
5. **exr-dc-bracket, line 336 — correct as the finite-dimensional bracket exercise.** \(K^{-1}_{12}=-1\); both signs and the full rank-one matrix identity were independently computed in item 3. I\(_+\) reverses the symplectic sign and uses its corresponding Lax matrix; changing just the sign of the bracket while retaining the old Lax matrix is not the exercise's solution.
6. **exr-dc-delay and exr-dc-type-two, lines 340,344 — correct.** The explicit delay identities and symbolic type-II momentum identities in items 4–5 supply independent solutions. The two no-jump auxiliary roots give inverse tanh/coth products. The auxiliary vacuum sector, not just a choice of algebraic expression for \(z\), must be retained.

**Chapter 40.**

1. **exr-dq-verlinde, line 282 — wrong final question.** The fusion coefficients computed from the independent Ising \(S\) give
   \[
   \mathcal H_\sigma=(R_1\otimes\bar R_\sigma)\oplus(R_\varepsilon\otimes\bar R_\sigma)
   \oplus(R_\sigma\otimes\bar R_1)\oplus(R_\sigma\otimes\bar R_\varepsilon).
   \]
   None has \((h,\bar h)=(1/16,1/16)\). The disorder field \(\mu\) instead belongs to the \((\sigma,\sigma)\) summand in
   \[
   \mathcal H_\varepsilon=(R_1\otimes\bar R_\varepsilon)\oplus
   (R_\varepsilon\otimes\bar R_1)\oplus(R_\sigma\otimes\bar R_\sigma).
   \]
   The body text correctly attaches the disorder field to an \(\varepsilon\)-line; the exercise contradicts it.
2. **exr-dq-chiral-defect, line 286 — locally correct.** Solving \(A-B=C+D,C-D=A+B\) gives the gluing in item 12. Its use for a compactified theory needs the separate momentum–winding compatibility condition.
3. **exr-dq-oscillator, line 290 — correct for a nondegenerate shift module at generic \(q\).** The recurrence's characteristic polynomial is
   \[
   1-[3]r+[3]r^2-r^3=(1-r)(1-q^2r)(1-q^{-2}r),
   \]
   giving \(u_n=A+Bq^{2n}+Cq^{-2n}\). At repeated roots this answer needs polynomial-in-\(n\) factors; division by a vanishing shift amplitude also requires separate treatment. Substitution of \(u_n=\kappa^{-2}\) and \(u_n=\kappa^{-2}(1-q^{2n})\), up to the spectral homomorphism factors, verifies Weyl and Fock cases.
4. **exr-dq-hopf, line 294 — correct.** A left coideal has \(\Delta b=\sum a_i\otimes b_i\) with \(a_i\) in the full algebra. The middle-leg decomposition would require the defect to carry the unrestricted \(a_i\) action; it need not. A boundary remains on the one permitted end, so its proof does not require this step.
5. **exr-dq-unique, line 298 — inherits the false hypothesis.** Item 8 solves the recurrence under generic-\(q\) assumptions and supplies the additional sixth-root counterexample. Asking only about \(q^4=1\) misses other roots of unity.
6. **exr-dq-transpose and exr-dq-weyl-trace, lines 306,310 — correct after stating orientation and formal-trace qualifications.** The exact coproduct calculation is in item 13; the finite Dirichlet kernel and limiting \(q^{ms/2}\) operator are in items 9 and 14.

**Paste-ready replacement for the final sentence of exr-dq-verlinde at line 282:**

~~~markdown
Repeat the calculation for a circle cut by the $\varepsilon$-line. In which twisted Hilbert space does a primary with $(h,\bar h)=(1/16,1/16)$ occur, and how does it represent the disorder field $\mu$?
~~~

At line 290 prepend:

~~~markdown
Assume $q$ is generic and that the upward shift is nonzero along the module; treat the Fock endpoint separately.
~~~

At line 298 replace the root-of-unity question by:

~~~markdown
Work first with $q$ not a root of unity. Explain why the four-channel form need not follow from the $k$-equations at a root of unity, and construct an additional intertwiner when $q=e^{i\pi/3}$. Separately examine the degeneracy of the displayed recurrence at $q^4=1$.
~~~

The transpose replacement is supplied in item 13. The Weyl-trace exercise should explicitly ask for a formal or distributional interpretation, rather than pointwise convergence.

## 17. Bibliography and prior art

**Verdict: wrong to treat Weyl-trace factorization or the general TQ/defect mechanism as new; correct to present the universal-R direction as a convention translation; unclear whether the explicit correction of the spin-two obstruction has prior publication.** This is a bounded primary-source audit, not proof of priority.

I checked the named BCZ/CZ sources, Weston, the relevant BLZ formulas, Robertson's folding, Bristow–Bowcock, the cited Caudrelier work, and additional primary material on quantum-group defects. The closest missing comparison is **Rossi–Weston 2002, §§4.3–5**, which is already in the bibliography but is not used where it matters most. I downloaded and read its primary text, together with BLZ III; copies are retained in the scratch directory.

**An independent calculation establishing the prior-art match.** For input particle signs \(a_j\) and output signs \(b_j\), put \(\sum a_j=\sum b_j=m\) and \(z_j=y_j/y_d\). Along a closed inverse-\(L^-\) path,
\[
n_j=n_0+\frac12\sum_{k\le j}(a_k-b_k).
\]
Its local weight is
\[
\frac{(iz_j)^{(1-a_jb_j)/2}}{1+z_j^2}\,
q^{b_j(n_j+s/2)}.
\]
There is at most one auxiliary-index path for a fixed pair of particle strings. Multiplying and simplifying gives the explicit matrix element
\[
\mathcal D_{\mathbf b,\mathbf a}
=q^{sm/2+\Omega(\mathbf a,\mathbf b)}
\prod_{j=1}^N
\frac{(iz_jq^{-1/2})^{(1-a_jb_j)/2}}{1+z_j^2},
\quad
\Omega=\frac14\sum_{i<j}(a_ib_j-b_ia_j).
\]
Different total charges give zero. The starting index contributes the additional factor \(q^{mn_0}\) before taking the \(n_0=0\) element.

[prior_art_check.py](../referee_scratch/defects_38_40/prior_art_check.py), [prior_art_output.txt](../referee_scratch/defects_38_40/prior_art_output.txt), compares this independently derived closed formula with the independent path-enumeration implementation for \(N=2,\ldots,6\), homogeneous and inhomogeneous rapidities. All relative errors are below \(4.34\times10^{-17}\).

In the homogeneous case, Rossi–Weston (5.5) is exactly this matrix, multiplied by its formal delta factor, under
\[
q_{\rm RW}=q^{-1},\quad s_{0,\rm RW}=q^{-s},\quad
\frac{z_{\rm RW}}{w_{\rm RW}}=-z_j^2,\quad
\widetilde\rho_{\rm RW}=\frac1{1+z_j^2},
\quad
\sqrt{\frac{z_{\rm RW}q_{\rm RW}}{w_{\rm RW}}}=iz_jq^{-1/2}.
\]
Their antisymmetric exponent has the opposite sign, which is exactly cancelled by \(q_{\rm RW}=q^{-1}\). Their (5.5) contains \(\delta(q_{\rm RW}^m)=\sum_{n\in\mathbb Z}q_{\rm RW}^{mn}\). Equations (5.9)–(5.10) explicitly identify and remove this factor to obtain Baxter's matrix. Their §4.3 discusses the divergent neutral block; (5.11) records the associated TQ statement and numerical checks through \(N=10\); Appendix B proves mutual commutation for the specialized matrices. Thus the basic factorization, divergent trace and stripped matrix are definite prior art, not merely a similar oscillator construction. [Rossi–Weston 2002](https://arxiv.org/abs/math-ph/0207004).

The chapter's inhomogeneous inverse-L normalization and nonneutral sector-dependent twist require a further convention translation before claiming that the *exact printed TQ identity* is already proved there. My independent checks verify the current coefficients through \(N=6\), but I have not supplied that full all-\(N\) translation. The homogeneous stripped operator itself is already identified exactly.

The resulting source assessment is:

| Claim | Primary-source comparison | Assessment |
|---|---|---|
| Both classical spin-two corrections | BCZ §5 and Penati–Polvara App. B contain the contradicted obstruction; later Lax/generating-function work supplies broader integrability context. | The explicit all-rank cancellation and independently integrated \(a_2\) correction are substantive. I did not locate the same correction in the sources inspected; priority remains unestablished. |
| I\(_+\)\(\leftrightarrow\mathcal R\), I\(_-\)\(\leftrightarrow\mathcal R_{21}^{-1}\) | KT/BGKNR fix universal-R conventions; Weston (4.15) and footnote 10 explain the transpose in the defect identification. | A useful convention reconciliation, not a new universal R or oscillator realization. |
| Weyl trace factorization and no ordinary bilateral convergence | Rossi–Weston §§4.3–5 explicitly give the delta factor and the matrix after removing it. | Already present, with the exact homogeneous matrix match above. |
| TQ as fusion of oscillator and particle modules | BLZ, Rossi–Weston, and Weston §2.4.1 provide functional relations and exact-sequence fusion. | The mechanism is established prior art. Current orientation, scalar factors, and the nonneutral twist need explicit translation. |
| Borel transmission matrices | Weston and CZ 2011 construct defect matrices from infinite-dimensional Borel representations. | Established; the positive/negative naming is convention-dependent. |
| Defects, transfer operators and moving defects through chains | Bajnok's spin-chain construction and Weston's construction already use these algebraic mechanisms. | Established background, not evidence for ordinary convergence of every infinite auxiliary trace. |
| Type-II and folding | CZ, Robertson, and Bristow–Bowcock give the constructions cited in the chapter. | The current text mostly attributes these appropriately. |

Relevant primary sources for these comparisons are [Weston 2010](https://arxiv.org/abs/1006.1555), [CZ 2011](https://arxiv.org/abs/1012.4186), [BLZ II](https://arxiv.org/abs/hep-th/9604044), [BLZ III](https://arxiv.org/abs/hep-th/9805008), [BGKNR](https://arxiv.org/abs/1205.1631), [Bajnok 2006](https://arxiv.org/abs/hep-th/0601107), [Robertson 2014](https://arxiv.org/abs/1304.3129), [Bristow–Bowcock](https://arxiv.org/abs/1612.03002), and [Bristow 2017](https://arxiv.org/abs/1709.03927). Bristow's abstract explicitly distinguishes an infinite set of charges from a completed Liouville-integrability proof; these notions should not be conflated in the book.

**Two specific citation qualifications.** First, BLZ III (2.6) explicitly uses the coproduct with \(q^{-h_i}\), confirming the book's \(q\mapsto q^{-1}\) dictionary. Second, BLZ III (3.11)–(3.12) defines its oscillator traces in a complex-momentum convergence domain and then analytically continues. The phrase at 40-defect-quantum-groups.qmd:187 that the twist “ensur[es] convergence” needs that domain; for physical real momentum it is not an ordinary decaying Fock weight.

**Paste-ready corrections.** Add after the Weyl-trace discussion at 40-defect-quantum-groups.qmd:215:

~~~markdown
The formal delta factor and its removal have a direct antecedent in Rossi–Weston [@rossi_weston2002, Secs. 4.3–5, especially eqs. (5.5), (5.9)–(5.11)]. In the homogeneous case, the matrix after removing that factor is Baxter's Q-matrix up to the present scalar and gauge conventions. Our sector-dependent twist and inverse-L normalization should be compared with that construction before assigning novelty to the functional relation.
~~~

At line 187 replace the convergence phrase with:

~~~markdown
with the oscillator traces first defined in a complex-momentum domain where they converge, and the resulting operators then continued analytically to other momenta [@blz1999, eqs. (3.11)–(3.12)]
~~~

At 38-defect-conditions.qmd:15 replace “Two results are new” by:

~~~markdown
Two calculations require particular attention: the cancellation of the purported spin-two obstruction, and the quadratic Poisson algebra of the defect matrix in the stated defect-coordinate construction.
~~~

This does not deny possible novelty of the explicit correction; it avoids claiming priority or a completed field-theory reduction on the basis of this limited audit. All named bibliography entries used in the principal comparisons exist. The main bibliographic defect is incomplete use of the already-cited Rossi–Weston paper, rather than a missing bibliography key.

## Other observations — unverified

- **Unverified:** I have not reconstructed the general-rank Drinfeld–Sokolov expansion with its infinite-line subtraction, Stokes sectors and exceptional defect jumps. The report establishes the cubic cancellation generally and the full local spin-two correction explicitly at \(a_2\).
- **Unverified:** I have not independently reproduced the separately reported \(a_3\) spin-two primitive. The book says that run is not part of its retained script; preserving a reproducible version would improve the evidence.
- **Unverified:** I have not completed an equal-time constrained-field-theory reduction proving that independent defect coordinates may be identified with both bulk endpoint values while retaining the assumed vanishing mixed brackets.
- **Unverified:** I have not established the boundedness/completion of every auxiliary-auxiliary intertwiner needed for analytic mutual commutation of the Fock operators at physical \(q\), or a unique regularization of the bilateral trace.
- **Unverified:** I have not translated Weston's full exact sequence, with scalar fusion factors and physical pole residues, into the book's inverse-L, inhomogeneous, nonneutral-sector TQ identity for arbitrary \(N\). The report's direct checks extend through \(N=6\).
- **Unverified:** The quantum higher-spin anomaly calculations of Penati–Polvara have not been re-audited. The classical cubic cancellation is not a result about those anomalies.
- **Unverified:** The full Robertson folding potentials, every Bristow–Bowcock generalized defect, and the exceptional-algebra existence question were not rederived. The affine-series identification and the explicitly tested SG/Tzitzéica formulas are the verified subset.
- **Unverified:** The compact-boson line requires a global lattice construction at the intended SG compactification radius; only its local gluing and the generic same-radius obstruction were checked here.
- **Unverified:** I have not audited all downstream Part VII formulas for applications of the transport sign. The homogeneous reflection equation and trace involution are unaffected by a consistent global sign, but mixed bulk/defect constructions must be checked individually.
- **Unverified:** The source search cannot establish that no earlier or later publication contains the explicit BCZ spin-two correction. It does establish a concrete earlier formula for the Weyl trace factorization and homogeneous stripped Q-matrix.

## Integration (2026-10-12)

The referee's scripts and outputs were copied to `briefs/referee_defects_38_40_scripts/` (no path changes were needed: `final_checks.py` finds `textbook_code/` at the same depth). All eight were rerun from there and reproduce the saved outputs exactly: `classical_checks.py` 0.8 s, `quantum_checks.py` 4.3 s, `hom_checks.py` 6.5 s, `supplemental_checks.py` 1.9 s, `type_two_checks.py` 1.3 s, `final_checks.py` 1.1 s, `prior_art_check.py` 0.9 s, `spin_check.py` 50 s (it regenerates the four `type_*_spin_*.txt` files, byte-identical). The book cites only `textbook_code/`. Independent checks were made in scratch scripts of the integrating session (no referee imports): the ch. 30 transport sign, the BCZ tensor normalization, an exact sympy version of the $a_2$ spin-two identity with the book's own code, the mirror defect's Lax matrix, the Gamma-ratio limit of `eq-sgd-scalar`, the $\mathcal W^\pm$ factorization exponents, the TQ relation at $N=5,6$, the closed form of the book's $\mathcal D$, the level-zero recursion and the infinite-module Hom argument; source texts of BCZ 2004 §5, CZ 2009b §4, Penati–Polvara §4, BLZ II (2.29) and BLZ III (3.11)–(3.12), and Rossi–Weston §§4.3–6 were read.

- **Item 1. Accepted.** Independent: solving BCZ (5.5) with $B=1-\mathsf D$ for $A_{abc}=\kappa\sum_k\mu_{k,a}\mu_{k,b}\mu_{k,c}$ by least squares gives $\kappa=2.000000$ for $r=2,3,4$ (residual $10^{-15}$; $B=1-\mathsf D^{\mathsf T}$ gives $-2$), so the book's "$A=\frac23\sum(\mu\cdot u)^3$" was off by the factor 3 of BCZ (5.2); the book script's numerical test was already (5.5) with $A=2\sum\mu^3$, only its label was wrong. The cot-identity expansion was rechecked by hand. With the book's own sympy code, the $a_2$ flux minus $d\mathcal B_{\pm2}/dt$ expands to exactly 0 with $\sigma$ symbolic, for both types and both spins (about 6 s per case). The vacuum-plus-analyticity argument for nonzero diagonal entries is correct (at the vacuum the leading coefficient of $N$ is $1$ at $\lambda\to\infty$ and $S$, diagonal in its eigenbasis, at $\lambda\to0$). Changes: ch. 38 `prp-spin-two` label (Theorem for (1) and for $a_2$ in (3); Numerical $a_3$; sketch (2)), statement ("every spin present in the bulk hierarchy", $r\ge2$, exactness for $a_2$), proof (1) tensor normalization, (2) vacuum argument and the remark that invertibility does not suffice, regularized generating function and exceptional values, (3) "symbolically to $10^{-15}$" replaced by random-point evaluation plus the exact identity, and the $a_3$ run marked as not independently reproduced; the BCZ paragraph ("symmetric tensor of the cubic term"); `exr-dc-cot` ($A=2\sum$, density $\frac13A$). Script `defect_classical_an.py`: check 11 docstring/label corrected; the $\mathcal B_{\pm2}$ construction is now exact (sympy, $\sigma$ symbolic) for both types.
- **Item 2. Accepted.** Independent: with the book's `ltilde_space`, the mirror defect $(\mathsf D^{\mathsf T},-E)$, with the relabelled potential, has a one-dimensional space of field-independent $\tilde L$ for $a_2$, $a_3$, both couplings and types, so the sentence in `prp-defect-lax` was false. Changes: ch. 38 `prp-defect-lax` (allowed matrices are (-@eq-d-matrix) and its mirror) and its proof; `defect_classical_an.py` check 4 now includes the mirror (dimension 1) alongside the $E\to0.6E$ control. `prp-ar-only` unchanged.
- **Item 3. Accepted.** Independent (scratch `sklyanin_sign.py`): with ch. 30's $a_x$ and $r$, the canonical local bracket is $+[r,a_{x1}+a_{x2}]\delta$ (coefficient $+1.000000$, residual $10^{-15}$), and a four-site lattice transport with $\partial_xT=-a_xT$ gives $c_T=-1.0000$ (residual $2\times10^{-7}$; control $+a_x$ gives $+1$), for $a_1$–$a_3$ at both couplings. The Dirac bracket $K^{-1}$ was rederived by hand. Changes: ch. 30 `prp-sklyanin` proof (bracket of $T$ is $-[r,TT]$; the reflection equations and the involution are homogeneous in $r$; the sign matters when combined with independently fixed brackets; cites `defect_lax_rmatrix.py` for the sign and notes that `boundary_classical.py` uses the opposite overall sign harmlessly). Ch. 38: learning goal; introduction ("two calculations need particular attention", finite-dimensional symplectic form on independent defect values, relation to the equal-time theory deferred); "regarded as independent coordinates"; `prp-defect-involution` label **[Numerical for $a_1^{(1)}$–$a_4^{(1)}$]** (the referee's "Conditional proposition" is not one of the book's labels; the condition is now in the statement), reference to `prp-sklyanin` for the bulk sign, and the remark that the step to involution is exact; "What the proposition does not settle" now says that the proposition concerns a finite-dimensional algebra. `eq-dc-defect-bracket` and `prp-defect-involution` are consistent with the corrected ch. 30. `defect_lax_rmatrix.py` check 3 extended to $a_4$; the control residual is now 0.2 (text updated).
- **Item 4. Correct; no change.** Independent: the $z_az_{h-a}$, conjugation and self-conjugate identities hold to $6\times10^{-14}$ for $h=2,\dots,6$ (scratch); the book script already checks the delays.
- **Item 5. Accepted (the optional insertion).** The referee's sewing checks of the corrected (4.4) are convincing (exact momentum identity, sewing $10^{-11}$, printed control 7.2); (4.5) follows algebraically. Independent: tracking the phase of the corrected $e^{\lambda_1}$ along imaginary $E_0$ gives a winding of $2\pi$ when $\rho,\tilde\rho>0$ and none when both are negative, the reverse of the printed formula and of CZ's discussion after (4.5). Change: ch. 38 Notes bullet on [@cz2009b, eq. (4.4)] now gives the corrected (4.5) and this consequence.
- **Item 6. Accepted.** The scope wording is right (the "only if" uses the independence of the $f$ and $g$ families, which needs $m\ne0$). The Penati–Polvara dictionary cell had the sign wrong: PP define $\gamma\equiv\beta^2/2\pi$ for the sinh-Gordon coupling (their §4; ch. 37 already uses $\gamma/2=b^2/8\pi$). Changes: ch. 38 `prp-defect-theta` first sentence (massive, $\beta\ne0$, $\sigma\in\mathbb C^\times$); `tbl-defect-dictionary` PP row ($b_{\rm SG}^2/8\pi$, continued to $-\beta_{\rm SG}^2/8\pi$).
- **Item 7. Accepted.** Trivially $|q^{s/2}|\ne1$ for complex $s$ at $|q|=1$ (scratch: 1.12 at $s=0.3+0.2i$). Changes: ch. 40 `prp-two-borels` (4) (crossing for all $s$; unitarity for real rapidities, real $s$, real $\theta_d$, $|q|=1$) and its proof. The universal-R direction and the BGKNR translation need no change.
- **Item 8. Accepted.** The counterexample at $q=e^{i\pi/3}$ was checked by hand and numerically ($[S^3,k]=0$ to $10^{-15}$); the level-zero recursion was rederived by hand (the $B$, $C$ recursions do not involve $c$; the two $D/A$ relations differ by $c^{-1}$). The infinite-module Hom argument was reimplemented independently in book conventions (exact sympy, $q=2$, $y_1=\frac32$, $y_2=\frac54$, $y_W=\frac76$): commutant dimensions 1, 1 for the Weyl modules of $\mathfrak b^\pm$, Fock cyclicity ranks 3, 4, raising generators invertible. Changes: ch. 40 learning goal; introduction ("what the argument uses"); "Why a Hopf subalgebra" (sufficient setting, no classification) and a global scope sentence ($q$ not a root of unity, fixed $q^s$); `thm-t-triangle` (iii) nonzero element, plus a remark that (iii) is a genuine check and closure is sufficient, not necessary; the hypothesis-(ii) paragraph (windows plus the exact infinite-module argument; not an irreducibility proof); `prp-sg-t-unique` hypotheses ($q$ not a root of unity, $y,y_d\ne0$, fixed $q^s$), proof (why the $k$-equations force the form; level-zero argument written out) and a root-of-unity remark; `exr-dq-unique`, `exr-dq-oscillator`. Ch. 41 line 11 (minimal): "once a uniqueness hypothesis holds and the two sides agree on one nonzero matrix element". `defect_universal_r.py`: new checks 11 (exact Hom argument) and 12 (root of unity), 16 checks.
- **Item 9. Accepted.** Independent, with the book's code: $\langle n|M^\pm|n\rangle=q^{\mp mn}\langle0|M^\pm|0\rangle$ for $\mathcal W^\pm$ (residual $10^{-16}$; the wrong exponent fails by 0.17–0.19). Changes: ch. 40 `prp-defect-transfer` (1) domain ($|q|=1$, $|\omega|<1$, away from poles; commutation when the manipulations are defined), (2) restricted to $\mathcal W^-$ with formal-delta/distribution reading, $\omega=\pm q^{-m/2}$, and the $\mathcal W^+$ exponent and support; proof (row factors of the inverse, formal twisted trace, boundedness at $|q|=1$); the physical-inverse sentence; ch. 36 `sec-defect-circle` paragraph "The defect as an operator" (convergent twisted trace, inverse L-operators, Fock = Q-operators, Weyl modules deferred to `prp-defect-transfer`). `defect_transfer.py` check 4 now also checks $\mathcal W^+$.
- **Item 10. Accepted.** With the book's code the Weyl TQ relation holds at $N=5$ and $6$ ($4\times10^{-15}$), so the numerical range was extended. Changes: ch. 40 `cnj-defect-tq` (**[Numerical]** for $N\le6$; conjecture kept), the paragraph after it, and the bootstrap reading (spectral shifts agree; the image/kernel identification needs Weston's exact sequence in these conventions and the scalar factors). `defect_transfer.py` check 6 runs $N=2,\dots,6$.
- **Item 11. Correct; no change.**
- **Item 12. Accepted.** The compact-boson lattice map $(n,w)\to(R^2w,n/R^2)$, the Ising fusion rules ($\mu$ in the $\varepsilon$-twisted sector) and the no-extension argument ($e_1|n\rangle=c_n|n+1\rangle$, $[e_1,f_0]=0$, contradiction with $[e_1,f_1]$) were rederived by hand. Changes: ch. 40 Schur sentence (chiral algebras), compact-boson caveat after (-@eq-dq-massless-defect), "Which non-local charges survive" (exact no-extension argument; the numerical residual kept as corroboration; no classification claimed); `exr-dq-verlinde` final question. `defect_universal_r.py` check 9 wording.
- **Item 13. Accepted (clarification).** The coproduct identity was checked by hand. Changes: ch. 40 proof of `prp-two-borels` (2) (reversal of the intertwining direction and of the tensor order) and `exr-dq-transpose`.
- **Item 14. Accepted.** Independent (scratch `item14_scalar.py`, log-Gamma products to $N=4\times10^5$): the Gamma ratio of `eq-sgd-scalar` tends to $e^{-i\pi(1+\lambda)/4}$ and $|\Phi-1|=1.2\times10^{-5}$ at $\theta=-12$ for $\lambda=0.73,1.37,2.6$. Changes: ch. 40 `prp-defect-twist-limit` hypotheses ($\lambda>0$, real rapidities, real $s$), the Gamma-ratio sentence of its proof, and "far below" (conversion with a shift of the defect charge, on the combined space). The book script already used the full $\Phi$ and passed, so only the explanation was wrong.
- **Item 15. Accepted, with the book's label vocabulary.** "Symbolic" and "Conditional proposition" are not book labels; the corresponding scope is written into the statements. Changes as listed under items 1, 3, 7, 8, 9, 14, plus `tbl-boundary-defect-algebra` (symmetry, quadratic equation, quadratic algebra with both $\mathbb L^\pm$, operator on the circle, functional relations $N\le6$, classical limit). `prp-l-operator` keeps its label; its proof already says "numerically on windows" for the Fock modules.
- **Item 16. Accepted.** Changes: `exr-dc-cot` (item 1), `exr-dq-verlinde` (item 12), `exr-dq-oscillator` (generic $q$, nonzero shift, Fock endpoint), `exr-dq-unique` (item 8), `exr-dq-transpose` (item 13), `exr-dq-weyl-trace` (formal or distributional support; $\mathcal W^+$).
- **Item 17. Accepted.** Rossi–Weston was read: §4.3 (matrix elements $\sum a_k\delta(q^k)$, "some care is needed" at $k=0$), (5.5), (5.9)–(5.10) (stripped $\tilde Q=s_0^{-n/2}Q_{\rm Bax}$), (5.11) (TQ for $n=0$ checked to $N=10$; for $n\ne0$ only at $q^n=1$), (5.12) and App. B. With the book's code, $\mathcal D$ has exactly the closed form $q^{sm/2+\Omega}\prod_j(iz_jq^{-1/2})^{(1-a_jb_j)/2}/(1+z_j^2)$, $\Omega=\frac14\sum_{i<j}(a_ib_j-b_ia_j)$ (residual $5\times10^{-16}$; $-\Omega$ fails), which has the structure of RW (5.5). BLZ II states that (2.29) makes sense only for $\operatorname{Im}p\gtrless0$ and is continued analytically. Changes: ch. 38 introduction ("two calculations need particular attention") and Notes (no priority claimed for the spin-two cancellation); ch. 40 new paragraph after `prp-defect-transfer` crediting Rossi–Weston (delta factor, stripped matrix, Baxter identification, TQ to $N=10$, commutation; our twist plays the role of their $q^n=1$; no novelty assigned before a detailed comparison); the BLZ sentence (convergence domain and continuation, cited to BLZ II after (2.29) rather than BLZ III, since the book cites BLZ II there); the Notes paragraph (RW credited; "stated here" list revised) and the Weston bullet (now a clarification crediting RW). No new bibliography entries were needed.

**Other observations — unverified.**
- General-rank Drinfeld–Sokolov expansion with infinite-line subtraction and exceptional jumps: not done; `prp-spin-two` (2) stays a sketch, and its proof now names the regularized generating function and exceptional values.
- The $a_3$ spin-two primitive: not reproduced; the text now says so. Retaining a reproducible run is listed as a pending item in the status file.
- Equal-time constrained reduction: not done; ch. 38 states the conditional nature in the proposition, the introduction and "What the proposition does not settle"; listed as an open problem for ch. 43.
- Boundedness of auxiliary–auxiliary intertwiners and a unique regularization of the bilateral trace: not done; `prp-defect-transfer` keeps the sketch label and now states the conditions.
- Weston's exact sequence in book conventions: not done; the bootstrap reading is now qualified; listed for ch. 43.
- Penati–Polvara quantum anomalies: not re-audited; nothing in ch. 38 relies on them beyond the existing citation.
- Robertson folding, Bristow–Bowcock generalized defects, exceptional algebras: not rederived; ch. 38 already says "We have not checked these results".
- Compact-boson line at the SG radius: not resolved; ch. 40 now carries the lattice caveat.
- Downstream uses of the transport sign: Part VII and chapters 36–42 were searched; only ch. 30's proof and ch. 38's `prp-defect-involution` use the sign of the bracket of $T$, and both are now consistent. `boundary_classical.py`'s `tau_bracket` uses $+[r,TT]$, harmless for involution (not edited; outside this integration's file list).
- Literature search for the spin-two correction: ch. 38 Notes now say that the sources checked do not contain it, without claiming priority.

**Book scripts after integration:** `defect_classical_an.py` 24 PASS (114 s), `defect_lax_rmatrix.py` 6 PASS (1.8 s), `defect_universal_r.py` 16 PASS (92 s), `defect_transfer.py` 9 PASS (29 s); no FAIL. (Runtimes on a shared machine.) `defect_doubling.py` and `defect_an_transmission.py` were not rerun: they import only `defect_common.py` and `boundary_common.py`, which were not changed.
