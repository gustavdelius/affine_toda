# Independent referee report: chapter 41

Date: 2026-10-12–13 (Europe/London). Model: GPT-6 (Codex). Commit at first inspection: ac0410f34631fc8619cc50ece360a0b1bc035d44. Commit when review copies were saved: 54dc21087e351dcf60f1fce29e2ee93674861187. Other agents advanced the repository during the review.

The reviewed chapter is part8-defects/41-soliton-transmission.qmd, SHA-256 57c5a7a3a5e2e713093bbaf7aa3b372204a57631c000c9d2d46f0e74a6d7d544. Its saved copy and manifest are in [referee_scratch/defects_41](../referee_scratch/defects_41/). Line numbers below refer to that copy. No pre-existing file was edited, no index-changing git command was run, and Quarto was not run. The only new files are this report and files under the permitted scratch directory.

The main matrix formulas survive independent checks. The principal corrections concern the scope and evidence for uniqueness, normalization of the reversed orientation, the closed endpoint of the instability window, the identification of a kernel restriction with a physical excited defect, and the ensuing unconditional charged-class resonance claims. There is also an inconsistency in the cited Tzitzéica pole formulas and a current/charge spin mismatch in the Penati–Polvara citation.

All scripts are independent of the chapter's verification script and its helpers. The final part of supplement_checks.py imports boundary_common **only to compare** its R-matrix with one already independently constructed; it is not the derivation of that R-matrix. To reproduce a calculation, run:

    OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_41/<script>.py

The scripts use NumPy and SciPy. exact_checks.py additionally implements exact arithmetic modulo 1009. The corresponding output files preserve the numerical results. There is no finite cutoff of the defect lattice: all matrices act between complete sectors of fixed total charge. Numerical checks, exact rank certificates, and analytic derivations are distinguished below.

## 1. Weyl module, intertwiners, gauge, triangle relation, and orientation

**Verdict: correct displayed intertwiners at generic parameters; unclear all-rank uniqueness; correct reversed-orientation family with a missing normalization qualification.**

**Independent evidence:** [independent_ar.py](../referee_scratch/defects_41/independent_ar.py), [ar_output.txt](../referee_scratch/defects_41/ar_output.txt), [exact_checks.py](../referee_scratch/defects_41/exact_checks.py), [exact_output.txt](../referee_scratch/defects_41/exact_output.txt), [supplement_checks.py](../referee_scratch/defects_41/supplement_checks.py), [supplement_output.txt](../referee_scratch/defects_41/supplement_output.txt).

### Module and Serre relations

Write \(n=N+s\) and use projected weights \(\mu_i=e_i-h^{-1}\mathbf1\). Projections are implicit in the chapter whenever \(e_i\) denotes a soliton weight. Put \(m=\mu/y_d\). Direct composition gives, for \(r\ge2\),

\[
f_jf_{j+1}|n\rangle
=m^2q^{n_{j+1}+n_{j+2}-1}|n-\alpha_j-\alpha_{j+1}\rangle,
\]
\[
f_{j+1}f_j|n\rangle
=m^2q^{n_{j+1}+n_{j+2}}|n-\alpha_j-\alpha_{j+1}\rangle.
\]

This includes the affine indices. Nonadjacent shifts leave each other's coefficient unchanged. Also

\[
k_if_jk_i^{-1}=q^{-\alpha_i\cdot\alpha_j}f_j,\qquad \prod_jk_j=1.
\]

If \(f_jf_k=cf_kf_j\), the quadratic Serre expression is

\[
[c^2-(q+q^{-1})c+1]f_kf_j^2
=(c-q)(c-q^{-1})f_kf_j^2.
\]

Thus the “if and only if” is correct **for adjacent nodes with nonzero shift operators**. The operator being multiplied is invertible in this bilateral module when \(m\ne0\), making cancellation legitimate. It is not an assertion about arbitrary modules where that operator could vanish. Reversing the exponent exchanges \(q^{-1}\) and \(q\); changing to the positive Borel reverses the shift and again gives Serre. For \(h=3,4,5\), both Borels and orientations have maximum residual \(2.01\times10^{-16}\). Plain shifts fail with residual 0.431535.

For \(r=1\), the affine Cartan entry is \(-2\), so the cubic rather than quadratic relation applies. Commuting inverse Weyl shifts satisfy it because \(1-[3]_q+[3]_q-1=0\). This verifies the rank-one portion of the Serre exercise.

### Solving without a lattice cutoff

Let \(A_j(Q)\) be the matrix of \(\Delta(f_j)\) between sectors \(Q\) and \(Q-\alpha_j\) in \(V\otimes\mathcal W\), and \(B_j(Q)\) its counterpart in \(\mathcal W\otimes V\). Each is an invertible diagonal matrix plus a particle lowering term. The intertwiner must obey

\[
T_{Q-\alpha_j}A_j(Q)=B_j(Q)T_Q.
\]

I form the two paths lowering by \(\alpha_j+\alpha_k\), then return by the inverse of one path. The resulting loops satisfy \(T_QL_A=L_BT_Q\). Solving this finite system determines an initial block, and transport determines other blocks. The positive Borel uses \(Q+\alpha_j\). Substitution then checks all original generator equations, including the affine one.

At \(\omega=1.37\), generic offset and \(y=1.19\), the solutions agree with the three displayed closed forms:

| \(h\) | \(V,I_-\) error | \(\bar V,I_-\) error | \(\bar V,I_+\) error | Loop-system nullity |
|---|---:|---:|---:|---:|
| 3 | \(1.34\times10^{-15}\) | \(2.17\times10^{-15}\) | \(9.50\times10^{-16}\) | 1 each |
| 4 | \(1.27\times10^{-15}\) | \(4.85\times10^{-15}\) | \(7.18\times10^{-16}\) | 1 each |

All three formulas also satisfy the original generator equations for \(h=5\), to \(1.64\times10^{-15}\). I separately construct the undisplayed \(I_+\) matrix on \(V\): it has respectively 3 and 4 nonzero entries per column, generator residuals \(2.94\times10^{-14}\) and \(6.34\times10^{-15}\), and triangle residuals \(1.52\times10^{-13}\) and \(2.31\times10^{-14}\).

There is an exact generic-uniqueness certificate. Composing another intertwiner with an invertible known one produces an endomorphism determined by one block and commuting with every transported loop. At the rational point \(q=2,\mu=3/7,y=5/4,s=0\), exact_checks.py finds modulo 1009:

    h=3, both multiplets, both Borels: block dimension 3, commutator rank 8, nullity 1
    h=4, both multiplets, both Borels: block dimension 4, commutator rank 15, nullity 1
    h=3, V(y1) tensor V(y2) tensor W: block dimension 9, rank 80, nullity 1

Every inverse exists in that finite field. These ranks certify nonzero rational minors, hence maximal rank at generic parameters. The identity supplies the matching upper bound. Thus, given existence of the displayed invertible intertwiners, this is a generic uniqueness proof without truncating \(\mathcal W\). The last line verifies hypothesis (ii) generically for \(a_2\), by the same endomorphism argument. It does not establish an all-rank theorem or uniqueness at every exceptional spectral point.

The proposition does not repeat its genericity assumptions. The conventions document already says that quantum statements assume generic \(q\), not a root of unity, so the following is a boundary counterexample to an unrestricted reading, not a contradiction of that intended convention. It makes the exclusion concrete and explains why it should accompany the standalone proposition and exercise. At \(q^5=1\), define

\[
S|N\rangle=|N+5\alpha_1\rangle.
\]

This non-scalar operator commutes with every Weyl-module generator, for either Borel. Hence \((S\otimes1)T\) is another intertwiner. This includes the physical coupling \(\omega=0.6\), where \(q=e^{2\pi i/5}\); it is not just the singular case \(q=\pm1\). Numerical \(k,f\) commutators are \(3.51\times10^{-16}\) and \(1.67\times10^{-16}\), with nonzero displacement \((5,-5,0)\). The exact reason is \(q^{5n}=1\). At such couplings, the \(k_j\)-equations alone do not enforce conservation of full topological charge.

Even under the existing generic-q convention, the cited matrices and finite-window checks do not establish the proposition's all-rank uniqueness assertion. The exact certificates here cover \(a_2,a_3\) generically; all ranks and all offsets would need an additional argument.

### Gauge and Part VI R

In [CZ2009a](https://arxiv.org/abs/0902.1307), equations (5.8), (5.11), set

\[
\alpha=n,\quad l_i=\mu_i,\quad Q_{\rm CZ}=q^{-1},
\quad \hat x=X,\quad \Delta=\theta_d+i\pi/2.
\]

Then \(Q_{\rm CZ}^{\alpha\cdot l_i}=q^{-n_i}\); the first charge delta is \(n'=n+\mu_i-\mu_{i-1}\). Also

\[
Q_{\rm CZ}^{-\alpha\cdot(l_i+\cdots+l_{i+j})}
=q^{n_i+\cdots+n_{i+j}},
\]

and the conjugate delta is \(n'=n-\mu_i+\mu_{i+j}\). This is an entry-by-entry match in the same gauge. The denominator in source (5.11) belongs to the scalar.

I independently solve the full finite-particle \(e,f,k\) coproduct equations for \(\check R\), retaining all equal-weight entries. Composing maps between the actual three tensor orders then gives triangle residuals below \(2.04\times10^{-15}\) for \(VV,\bar V\bar V,V\bar V\), \(h=3,4\), for \(I_-\). The displayed \(I_+\) matrix gives residuals below \(1.88\times10^{-15}\) for \(\bar V\bar V\). The independently obtained R-matrices agree up to normalization with boundary_common.rcheck within \(8.08\times10^{-16}\). These checks therefore use the Part VI R in the stated convention.

### Reversed orientation

The exponent reversal gives a full \(V\) matrix and a two-entry \(\bar V\) matrix. Direct solution gives, in cyclic indices,

\[
M^-_{\bar V}(\bar v_i\otimes|n\rangle)
=q^{n_i}|n\rangle\otimes\bar v_i
+Xq^{n_i-n_{i+1}}|n-\alpha_i\rangle\otimes\bar v_{i+1}.
\]

Its residual against the independent \(a_2\) solution is \(1.53\times10^{-15}\). The barred index moves forward, giving the reversed rotation.

But the normalization matters. In the basis \(|n\rangle'=q^{-|n|^2/2}|n\rangle\), the converting coefficient is **\(qX\)**. CZ2007 (5.10) has coefficient \(\hat x\). Thus unchanged \(\mu\) gives a member of the mirror family with \(\hat x_{\rm mirror}=qX\), not the same physical parameter \(\hat x=X\). To keep that parameter, use \(\mu_-=\mu_+/q\). The same result follows by applying \(n_i\mapsto-n_{-i}\) to the original module and making a quadratic basis change.

The classical source dictionary is correct: CZ2007 (2.12) is
\(2\sum\lambda_a(\lambda_a-\lambda_{a+1})^T=D_{\rm book}^T\), and (2.20) is \(z^{\rm book}_{h-a}\). supplement_checks.py verifies these for \(h=3,4,5\), with errors 0 and at most \(2.69\times10^{-15}\). This establishes the transposed-defect/species identification, while retaining the quantum normalization qualification.

**Ready-to-paste corrections:**

- At line 62, replace the blanket status and scope with:

  > **[Established type-I$_-$ matrix solutions; generic uniqueness checked for $a_2,a_3$; other statements numerically checked as specified below]** Let $r\ge2$, take $q$ not a root of unity and generic spectral parameters, and let $\mathcal W$ be (-@eq-dt-weyl). The uniqueness statements below are restricted to the ranks checked here; no all-rank uniqueness theorem is claimed.

  Qualify items 1–2 themselves accordingly. If the exact loop argument above is incorporated, “checked” can be strengthened to “proved generically” for \(a_2,a_3\).

- Append to line 36:

  > This identifies a family of mirror defects. With the same $\mu$ as in (-@eq-dt-weyl), a quadratic change of the defect basis gives the CZ2007 gauge with $\hat x=qX$. To retain $\hat x=X$, use $\mu_-=\mu/q$ for the reversed-orientation module.

- At line 238, start the exercise with “Assume \(q\) is not a root of unity and the parameters are generic.”

## 2. Crossing, the q-sign, and unitarity

**Verdict: correct within the fixed conventions at generic parameters.** The broader Klein-factor question for odd \(h\) remains unclear, as the chapter explicitly acknowledges.

**Independent evidence:** independent_ar.py and ar_output.txt, including both crossed channels and complete charge blocks at three offsets.

The antipode gives

\[
\pi_y^*(e_j)=(-k_j^{-1}e_j)^T=-\frac yq E_{j+1,j},
\qquad \pi_y^*(k_j)=q^{-H_j}.
\]

Hence \(V(y)^*=\bar V(-y/q)\) with \(C=1\). At physical \(q\), this is exactly \(\theta\mapsto\theta+i\pi\).

In a charge block \(M_V=D_Q+XS_-\),

\[
(D_Q)_{ii}=q^{(h-1)/h-Q_i},\qquad \prod_i(D_Q)_{ii}=q^{h-1}.
\]

The finite geometric inverse has denominator \(1+(-q^{-1})^{h-1}X^h\). At \(X'=-qX\) this becomes \(1-qX^h\), explaining \(c_{\rm II}\). Crossing requires inversion of the conjugate block of charge \(n-\mu_j\) for incoming particle \(i\), outgoing \(j\), incoming defect \(n\); inverting the original block instead would use the wrong charge.

Implementing those indices for both channels gives

\[
c_{\rm I}=1+(-q)^{h-1}X^h,\qquad c_{\rm II}=1-qX^h.
\]

The product of diagonal factors is independent of all offsets with \(\sum s_i=0\). For \(s=0,\lambda_1\), and a generic \(s\), the maximum crossing residuals for \(h=3,4,5\) are \(3.07,5.08,5.51\times10^{-16}\); the maximum scalar error is \(1.67\times10^{-15}\).

Eliminating \(\bar g(\theta+i\pi)\) from crossing II at \(\theta+i\pi\) and crossing I at \(\theta+2\pi i\) gives

\[
\frac{g(\theta+2\pi i)}{g(\theta)}
=\frac{c_{\rm II}(\theta+i\pi)}{c_{\rm I}(\theta+2\pi i)}
=\frac{1+(-q^{-1})^{h-1}X^h}{1+(-q^{-1})^{h+1}X^h}.
\]

This is CZ2009a (6.8).

At the wrong sign of \(q\), the crossed representation differs by \(y\mapsto-y\). A Cartan-preserving conjugation is diagonal and must satisfy \(c_{i+1}=-c_i\). Closing the affine cycle requires \((-1)^h=1\). More invariantly, \((\sum e_j)^h=y^h1\), so for odd \(h\) changing the sign of \(y\) changes a similarity invariant: no constant conjugation works. For even \(h\), \(C'=\operatorname{diag}((-1)^i)\) works.

The independent wrong-q crossing residuals, after fitting an arbitrary common scalar, are:

| \(h\) | \(C=1\) | \(C'\) |
|---|---:|---:|
| 3 | 0.9851 | 0.9833 |
| 4 | 0.9851 | \(4.61\times10^{-16}\) |
| 5 | 0.9851 | 0.9979 |

This proves the fixed-representation claim, not impossibility of a changed exchange algebra. A Klein-type transformation could also change mixed amplitudes and the identification of fused antiparticles. Choosing \(\bar V(-y)\) at the wrong \(q\) removes the representation mismatch for odd \(h\), but changes that dictionary; no consistent full Klein transformation has been constructed here. The reservation at line 131 must remain.

For \(h\ge3\), \(S_-\) and \(S_-^\dagger\) connect distinct entries, so

\[
M^\dagger M=(1+|X|^2)1+XD^\dagger S_-+X^*S_-^\dagger D
\]

has nonzero off-diagonal terms, which no scalar factor removes. The largest off-diagonal entry at the test point is exactly \(|X|=1.19\) for \(h=3,4,5\). This is failure of Hilbert-space unitarity in the displayed basis. Algebraic defect unitarity instead defines the reverse transmission by \(\bar{\mathbb T}(-\theta)\mathbb T(\theta)=1\). The chapter's distinction is correct.

**Correction required:** none to the displayed crossing formulas or the convention-qualified q-sign conclusion.

## 3. Scalar factor, printed-product correction, and bulk bootstrap

**Verdict: correct convergent scalar and bootstrap; correct \(K^{-z}\) drift after normalization at \(z=0\).** That normalization should be explicit in the product comparison and Gamma exercise.

**Independent evidence:** [scalar_checks.py](../referee_scratch/defects_41/scalar_checks.py), [scalar_output.txt](../referee_scratch/defects_41/scalar_output.txt), [fusion_checks.py](../referee_scratch/defects_41/fusion_checks.py), [fusion_output.txt](../referee_scratch/defects_41/fusion_output.txt).

Put \(H=h\omega\), \(A=\frac12+(h-1)\omega/2\), \(B=A+\omega\), and call the chapter's normalized product \(\mathcal H(z)\). Cancelling constant and linear terms in the integral for \(\log\Gamma\) gives independently

\[
\log\mathcal H(z)=
\int_0^\infty
\frac{2e^{-At}(1-e^{-\omega t})(\cosh zt-1)}
{t(1-e^{-t})(1-e^{-Ht})}\,dt,\qquad |\operatorname{Re}z|<A.
\]

Expanding the last denominator produces exactly the equal-\(k\) Gamma product; meromorphic continuation extends it. Its summand satisfies

\[
\log\frac{F_k(z)}{F_k(0)}
=z^2[\psi'(A+kH)-\psi'(B+kH)]+O(k^{-4})
=\frac{\omega z^2}{(A+kH)^2}+O(k^{-3}).
\]

This proves local convergence away from poles, independently of numerical extrapolation.

Telescoping the Gamma factors gives the functional identities

\[
\frac{\mathcal H(z-H)}{\mathcal H(z)}
=\frac{\sin\pi(A-z)}{\sin\pi(B-z)}
\]

and

\[
\frac{\prod_{j=0}^{h-1}\mathcal H(z+\omega(h/2-j))}
{\prod_{j=0}^{h-1}\mathcal H(\omega((h-1)/2-j))}
=\frac1{\cos\pi(z+\omega/2)}.
\]

For the second identity, telescope first in \(j\), then in the cutoff index \(k\). The low-end pair is \(\Gamma(\frac12+u)\Gamma(\frac12-u)\), \(u=z+\omega/2\). After comparison with \(u=0\), the far-end pair tends to 1; Euler reflection gives the cosine.

Now \(X=e^{-2\pi iz/h}\) and \(g=c_he^{\pi iz/h}\mathcal H(z)\). The constituent shifts in \(\bar g\), together with \(g(\theta-i\pi)\), give precisely the \(h\) arguments in the second identity, while

\[
1-qX^h=2e^{-i\pi(z+\omega/2)}\cos\pi(z+\omega/2).
\]

Thus crossing II fixes

\[
c_h=\left[
2\prod_{j=0}^{h-1}\mathcal H\!\left(\omega\left(\frac{h-1}{2}-j\right)\right)
\right]^{-1/h}.
\]

All factors defining this constant are positive and inside the convergence strip. The first identity, including the factor \(e^{-i\pi\omega}\) from \(X^{-1/2}\), gives (6.8); crossing I follows. This constructs the stated minimal representative, with homogeneous/CDD freedom still present. Functional equations alone do not remove that freedom.

Independent output:

| \(h\) | \(\omega\) | \(c_h\) | Maximum error: both crossings and (6.8) |
|---|---:|---:|---:|
| 3 | 1.37 | 0.398306002479 | \(6.63\times10^{-10}\) |
| 4 | 1.37 | 0.393286370853 | \(5.46\times10^{-10}\) |
| 5 | 1.37 | 0.390976898921 | \(1.35\times10^{-9}\) |
| 3 | 0.6 | 0.631547074024 | \(4.01\times10^{-10}\) |
| 4 | 0.6 | 0.643078061887 | \(2.11\times10^{-10}\) |
| 5 | 0.6 | 0.648797400259 | \(1.54\times10^{-9}\) |

The implementation uses cutoffs 100, 200, 400, 800 and extrapolates the inverse-cutoff expansion just derived. The integral is used separately for the Wigner test.

### The exact printed-product discrepancy

Let \(P_K(z)\) be the Gamma part of printed CZ2009a (6.10), including the Gamma outside its product, cut off at \(k=K\). Set

\[
\mathcal H_K(z)=\prod_{k=0}^{K-1}\frac{F_k(z)}{F_k(0)},\qquad A_K=A+KH.
\]

Reading the indices in the source gives the exact finite-product equality

\[
\boxed{\frac{P_K(z)}{P_K(0)}
=\mathcal H_K(z)\frac{\Gamma(A_K-z)}{\Gamma(A_K)}}.
\]

The normalized printed product is therefore asymptotic to \(\mathcal H(z)(h\omega K)^{-z}\). The raw product also contains \(P_K(0)\), a cutoff-dependent factor; it is not literally just \(K^{-z}\) times a convergent product. Normalizing at zero and removing the endpoint factor yields the book's product.

For \(h=3,\omega=1.37,z=0.43+0.72i\), doubling the cutoff from 2000 to 4000 changes the normalized printed logarithm by \(-0.298014-0.499012i\), versus \(-z\log2=-0.298053-0.499066i\). The equal-\(k\) logarithm changes only by \(-6.76\times10^{-6}+1.26\times10^{-5}i\). The rapidity-dependent drift is genuine.

### Independent matrix bootstrap

I solve the finite-dimensional coproduct equations for

\[
\bar V(\theta)\longrightarrow
\bigotimes_{m=0}^{r-1}V\!\left(\theta+\frac{i\pi(2m-r+1)}h\right),
\]

then transmit the constituents from right to left. The result agrees **without a proportionality fit** with the embedded \(M_{\bar V}\), with errors \(7.39\times10^{-16}\) for \(a_2\) and \(1.01\times10^{-14}\) for \(a_3\). The embedding nullity is one. This fixes the matrix ratio to 1 and hence gives \(\bar g=\prod g\).

For even \(r\), the shifts are the pairs \(\pm(2a+1)i\pi/h\); for odd \(r\), they are 0 and \(\pm2ai\pi/h\). These are CZ2009a (6.6)–(6.7).

**Ready-to-paste corrections:**

- At line 133, replace “The product is the same, but its partial products…” with:

  > After division by their value at $z=0$, the printed partial products equal $\prod_{k=0}^{K-1}F_k(z)/F_k(0)$ times $\Gamma(A_K-z)/\Gamma(A_K)$, where $A_K=\frac12+(hK+\frac r2)\omega$. The endpoint factor behaves as $(h\omega K)^{-z}$, producing a rapidity-dependent drift. Removing it yields (-@eq-dt-g).

- At line 250, replace the product-comparison sentence with:

  > Let $P_K(z)$ be the Gamma part of [@cz2009a, eq. (6.10)], including its prefactor Gamma function, truncated at $k=K$. Show that $P_K(z)/P_K(0)=\left[\prod_{k=0}^{K-1}F_k(z)/F_k(0)\right]\Gamma(A_K-z)/\Gamma(A_K)$, and derive its $(h\omega K)^{-z}$ drift.

## 4. Poles, instability windows, and the kernel bootstrap

**Verdict: correct generic pole, determinant, rank, species dictionary, and algebraic kernel restriction; wrong inclusion of the zero-width endpoint; unclear physical excited-defect identification, with an invalid argument for selecting its scalar-factor shift.**

**Independent evidence:** [independent_ar.py](../referee_scratch/defects_41/independent_ar.py), [ar_output.txt](../referee_scratch/defects_41/ar_output.txt), [scalar_checks.py](../referee_scratch/defects_41/scalar_checks.py), [scalar_output.txt](../referee_scratch/defects_41/scalar_output.txt), [fusion_checks.py](../referee_scratch/defects_41/fusion_checks.py), [fusion_output.txt](../referee_scratch/defects_41/fusion_output.txt), [exact_checks.py](../referee_scratch/defects_41/exact_checks.py), [exact_output.txt](../referee_scratch/defects_41/exact_output.txt).

### Pole and rank

Solving \(A-z=0\), with \(z=ih\omega(\theta-\theta_d-i\pi/2)/(2\pi)\), gives

\[
\theta_1=\theta_d+\frac{i\pi}2-i\left(\frac{\pi(h-1)}h+\frac{\pi}{h\omega}\right).
\]

On the complete sector with total projected charge \(Q\), the diagonal entry is
\(D_{ii}=q^{(h-1)/h-Q_i}\), and the cyclic off-diagonal entry is \(X\).
Only the identity permutation and the full cycle contribute to the determinant. Since \(\sum_iQ_i=0\),

\[
\det M_Q=q^{h-1}+(-1)^{h-1}X^h
=q^{h-1}\left[1+(-q^{-1})^{h-1}X^h\right].
\]

At the pole, the bracket vanishes. The kernel recurrence is
\(u_{i+1}=-D_{ii}u_i/X\); every component is nonzero. A proper diagonal minor is nonzero, so the rank is exactly \(h-1\), not merely at most \(h-1\). For generic coupling the scalar pole is simple and nonzero, giving the same rank for the residue. Exceptional coincidences of Gamma poles or zeros require separate treatment.

The singular values from independent blocks at the pole are:

| \(h\) | Singular values |
|---|---|
| 3 | \(1.73205,\ 1.73205,\ 9.61\times10^{-17}\) |
| 4 | \(2,\ 1.41421,\ 1.41421,\ 4.52\times10^{-16}\) |
| 5 | \(1.90211,\ 1.90211,\ 1.17557,\ 1.17557,\ 5.14\times10^{-16}\) |

The \(m=0\) factor in the conjugate bootstrap has its pole at
\(\theta_d+i\pi/2-i(\pi/h+\pi/(h\omega))\), as stated. At \(h=3,\omega=1.37\), the estimates \(\epsilon g(\theta_1+\epsilon)\) and \(\epsilon\bar g(\theta_r+\epsilon)\), with \(\epsilon=10^{-5}\), are respectively
\(-0.902644-0.368365i\) and \(-0.264932-0.341552i\). The output records convergence from \(\epsilon=10^{-3}\), and analogous checks for \(h=4,5\).

### The window needs both a strict endpoint and a specified strip

Direct substitution gives

\[
\frac{E_a}{M_a}=\cosh\theta_d\sin\vartheta_a
+i\sinh\theta_d\cos\vartheta_a,\qquad
\vartheta_a=\frac{\pi(h-a)}h+\frac{\pi}{h\omega}.
\]

For \(\theta_d>0\), in the nearest lower strip
\(-\pi/2<\operatorname{Im}\theta_a<0\), the two strict energy inequalities are equivalent to

\[
\frac\pi2<\vartheta_a<\pi,\qquad
\boxed{\max(0,a-h/2)<1/\omega<a,\quad \omega>0.}
\]

The lower bound \(1/\omega>0\) is automatic at finite positive \(\omega\). It is appropriate to discuss \(1/\omega=0\) as a classical limit, not as an additional finite-coupling state.

At \(\vartheta_a=\pi/2\), the pole is real and \(\operatorname{Im}E_a=0\). For example \(h=3,a=2,\omega=2,\theta_d=0.7\) gives \(E_a/M_a=1.25516900563+0i\). It is a zero-width endpoint, not a decaying state. CZ2009a (6.15) already has this endpoint inconsistency; reproducing it is not independent validation. CZ2007's strict window \(1/2<\gamma<2\) has the correct endpoint.

The energy-sign conditions alone are periodic. Without restricting the pole to the nearest strip, they also hold when

\[
a-h/2+2mh<1/\omega<a+2mh,\qquad m\in\mathbb Z,
\]

intersected with \(1/\omega>0\). For instance, \(h=3,a=1,\omega=1/6.5,\theta_d=0.7\) gives a pole \(0.7-7.33038i\) and \(E/M=0.627585-0.656953i\), outside the stated first window. This is not a proposal to interpret remote-sheet poles as additional unstable particles; it shows why the strip restriction is needed in the claimed equivalence.

The first-strip windows, with \(t=1/\omega\), are:

| Theory | Species \(a\) | Strict interval for \(t\) |
|---|---|---|
| \(a_2\) | 1, 2 | \((0,1),\ (1/2,2)\) |
| \(a_3\) | 1, 2, 3 | \((0,1),\ (0,2),\ (1,3)\) |
| \(a_4\) | 1, 2, 3, 4 | \((0,1),\ (0,2),\ (1/2,3),\ (3/2,4)\) |

For \(a_2\), species 1 has \(\omega>1\), or \(0<\beta^2<2\pi\). Species 2 has \(1/2<\omega<2\), or \(4\pi/3<\beta^2<8\pi/3\).

The dictionary with CZ2007 is correct: direct construction of their \(D=2\sum_a\lambda_a(\lambda_a-\lambda_{a+1})^T\) gives the transpose of the book's \(D\), and their delay for species \(a\) is the book's delay for \(h-a\). The independent residual is zero for the transpose and at most \(2.69\times10^{-15}\) for the delay identity at \(h=3,4,5\), as recorded in supplement_output.txt. Thus their species 1 is the book's species 2 at \(h=3\).

The correction to CZ2009a's species-independent bound is also correct. Their bound agrees with \(t<a\) only for the largest integer \(a<h/2\). For example \(h=5,a=1,t=1.5\) gives \(\vartheta_a=1.1\pi\), hence negative real energy despite satisfying their quoted bound \(t<2\).

### What the kernel calculation actually establishes

I construct a null vector of \(M_V(\theta_1)\), embed the corresponding one-dimensional sector, and apply the independently constructed particle \(\check R_{12}\) followed by \(M_{23}\). For \(h=3,4\), at offsets \(s=0,\lambda_1\), the output gives:

| Quantity | Largest error |
|---|---:|
| Invariance of kernel subspace | \(2.39\times10^{-15}\) |
| Diagonal ratios of the restricted matrix | \(6.67\times10^{-16}\) |
| Entries forbidden in a two-entry cyclic matrix | \(6.04\times10^{-16}\) |
| Cycle ratio compared with \(q^2\) | \(1.66\times10^{-15}\) |

Here \(q^2=-0.684547105929-0.728968627421i\). The offset becomes \(s+\lambda_1\), and the product of the cyclic conversion entries relative to the diagonal entries acquires \(q^2\). This independently verifies the **matrix equivalence up to diagonal gauge and overall scalar**, \(X^h\mapsto q^2X^h\).

It does not identify the kernel with the physical excited defect defined earlier in the book. Let \(\mathcal R=\operatorname{Res}\mathbb T\). The applicable exact sequence is

\[
0\longrightarrow\ker\mathcal R
\longrightarrow V(\theta_1)\otimes\mathcal W
\longrightarrow\operatorname{im}\mathcal R
\longrightarrow0,
\]

with the last map understood as the surjection onto the image in \(\mathcal W\otimes V(\theta_1)\). The residue-image prescription uses the quotient by the kernel, with dimension \(h-1\) per sector. The kernel has dimension 1. They are not automatically isomorphic, and the kernel is not “its complement”: they initially live in different ordered tensor products, and no splitting has been established.

Moreover, the cycle invariant permits every representative

\[
\delta\theta_d=\frac{2\pi ip}h+\frac{2\pi ik}{h\omega},\qquad k\in\mathbb Z.
\]

Every fixed \(k\) has the same classical limit. Therefore the sentence “the classical limit selects \(2\pi i/h\)” is false. Neither the matrix invariant nor that limit selects \(k=0\). An excited-defect scalar bootstrap, and possible CDD factors, remain necessary.

There is a concrete counterexample to the proposed inference about quantum widths. Suppose provisionally that the scalar is the shifted minimal scalar. For a real classical absorption with \(a+2p=h/2\), the representatives give

\[
\theta_{\rm pole}=\theta_d+\frac{i\pi(2k-1)}{h\omega}.
\]

With \(h=6,a=p=1,\omega=1.37,\theta_d=0.7\), the independent output is:

| \(k\) | Pole | \(E/M\) | Cycle invariant |
|---|---|---|---|
| 0 | \(0.7-0.382189i\) | \(1.164609-0.282916i\) | \(q^2\) |
| 1 | \(0.7+0.382189i\) | \(1.164609+0.282916i\) | \(q^2\) |

The same verified matrix invariant and the same classical limit allow opposite signs of the imaginary energy under these scalar assignments. This does not prove that both choices define physical defects. It proves that the given evidence cannot select the claimed resonance. The formula at \(k=0\) is algebraically correct under the stated scalar conjecture, including the \(a_5\), species-1, class-1 example; its unconditional physical interpretation is not established.

**Ready-to-paste corrections** in part8-defects/41-soliton-transmission.qmd:

- Replace lines 173–177, retaining the equation label:

  > **Unstable states.** For $\theta_d>0$ we consider the pole in the nearest lower strip $-\pi/2<\operatorname{Im}\theta_a<0$. The conditions $\operatorname{Re}E_a>0$ and $\operatorname{Im}E_a<0$ are then equivalent to $\pi/2<\vartheta_a<\pi$. Thus, at finite positive $\omega$, the window is
  >
  > $$\max(0,a-h/2)<1/\omega<a.$$
  >
  > The endpoint $\vartheta_a=\pi/2$ has zero imaginary energy and is excluded. The non-strict endpoint in [@cz2009a, eq. (6.15)] does not satisfy the stated decay condition.

- In line 185 replace the species-2 intervals by:

  > $\frac12<\omega<2$, that is $\frac{4\pi}3<\beta^2<\frac{8\pi}3$.

- Replace line 196:

  > **Possible quantum class assignment.** The kernel restriction in item 3 has offset $s+\vec\lambda_1$ and the cyclic matrix invariant of a class-1 defect. This is compatible with the classical class shift, but does not identify the kernel with the excited defect defined by the residue image. The matrix invariant fixes the rapidity shift only modulo $2\pi i/(h\omega)$, and every fixed representative has the same classical limit. Assigning the scalar $g(\theta-2\pi ip/h)$ to class $p$ is therefore an additional conjecture; it has not been derived from the physical excited-defect bootstrap.

- Replace the first three sentences of line 198:

  > **The image of the residue.** The excited-defect prescription of @sec-defect-bound-states uses the image of the residue in $\mathcal W\otimes V(\theta_1)$. It has dimension $h-1$ per sector and is isomorphic to the quotient of $V(\theta_1)\otimes\mathcal W$ by the kernel. The kernel has dimension 1 per sector and carries the module described above; no identification with the image or direct-sum splitting has been established.

- Replace the opening of line 200:

  > **Conditional poles in charged classes.** If the class-$p$ scalar is chosen to be the shifted minimal scalar $g(\theta-2\pi ip/h)$, with this representative and no additional CDD factor, its species-$a$ pole lies at

- After the pole formula, replace the physical conclusion in line 206:

  > Under this scalar assignment, the first-strip window is obtained by replacing $a$ by $a+2p$, for $0<a+2p\le h$. Real classical absorptions have candidate poles at $\theta_d-i\pi/(h\omega)$, which have decaying energies for $\theta_d>0$ and $0<1/\omega<h/2$. This includes the species-1, class-1 example in $a_5^{(1)}$. Their identification as resonances of a quantum Lagrangian defect remains conditional on the scalar assignment and the excited-defect bootstrap.

- In line 254 explicitly restrict the first question to the nearest lower strip, use the strict interval, and replace “For a defect of class \(p\), show…” by:

  > Assuming the shifted minimal scalar assigned in the charged-class discussion, derive the corresponding first-strip window with $a\to a+2p$.

**Concurrent addition checked:** at final inspection, commit 70b57c11a86309bfa62c409c59d707dd25c094fa, the chapter also contains a new paragraph and figure caption between saved lines 192 and 194. The caption repeats the closed endpoint. In that caption replace the green-band condition by $-\frac12<\operatorname{Im}(\theta_a-\theta_d)/\pi<0$, use the strict species-2 intervals above, and identify the top boundary as a zero-width endpoint. The figure's endpoint markers should follow the same convention. This addition moves subsequent live-file line numbers by four; the report's original line references remain tied to the saved copy. I have not audited the figure-generation script.

## 5. Selection rule and semiclassical test

**Verdict: correct suppression and Wigner prefactor; correct numerical semiclassical test at the stated sample points; wrong literal “at the edge” formulation in the exercise.**

**Independent evidence:** [scalar_checks.py](../referee_scratch/defects_41/scalar_checks.py), [scalar_output.txt](../referee_scratch/defects_41/scalar_output.txt), [exact_checks.py](../referee_scratch/defects_41/exact_checks.py), [exact_output.txt](../referee_scratch/defects_41/exact_output.txt).

### Suppression

For real charges and physical \(q\), all the explicit \(q\)-powers have modulus 1. Put \(x=|X|=e^{\omega(\theta-\theta_d)}\). The largest of the \(h\) entries in an incoming column of \(M_{\bar V}\) has modulus \(\max(1,x^{h-1})\); hence the ratio of entry \(j\) to that maximum is exactly

\[
\frac{x^j}{\max(1,x^{h-1})}
=\min(x^j,x^{j-h+1}).
\]

For \(1\le j\le h-2\), this tends to zero for fixed \(\theta-\theta_d\ne0\), on either side. At zero every ratio is 1; the transition has width \(O(\omega^{-1})\). At \(|\theta-\theta_d|=0.3\), the largest forbidden ratio for \(h=3,4\) is:

| \(\omega\) | Ratio |
|---|---:|
| 1.37 | 0.662986932 |
| 10 | 0.0497870684 |
| 40 | \(6.14421235\times10^{-6}\) |

The allowed endpoint transitions are \(\mu_i\to\mu_i,\mu_{i-1}\) for \(V\), and \(-\mu_i\to-\mu_i,-\mu_{i-1}\) for \(\bar V\). The latter conversion is \(j=h-1\). These agree with the stored-charge rule of chapter 39. The intermediate weights are forbidden **for that incoming soliton and defect orientation**; they are not weights absent from the classical species-\(r\) multiplet. All its \(h\) weights occur classically.

The relative quantity \(|X|^2/(1+|X|^2)\) tends to the stated step function. Since the matrices are not unitary, its use here is a normalized comparison of two entries, not a derivation of a physical transition probability.

### Wigner normalization and an independent derivative

Combining the classical displacement and Wigner relation directly gives

\[
\frac{d\delta}{d\theta}
=-M_1\cosh\theta\,\Delta x
=\frac{M_1}{m_1}\log|z_1|
=\frac{2h}{\beta^2}\log|z_1|.
\]

Thus the multiplier of the quantum phase derivative is exactly

\[
\frac{\beta^2}{2h}=\frac{2\pi}{h(\omega+1)}.
\]

The dominant matrix entries have rapidity-independent phase, so their scalar phase derivative is the derivative of \(\arg g\). The real exponential prefactor in \(g\) contributes no phase derivative on the real rapidity axis.

I differentiate the convergent integral in section 3, rather than finite-differencing the chapter's product. In the implementation the \(k=0\) term is differentiated with digamma functions, and the remaining tail is integrated. This separation keeps the integral in its convergence domain for the shifted arguments needed for \(\bar g\). At \(\theta-\theta_d=-0.7,0.5,1.1\), the maximum absolute errors against the independently evaluated classical \(\log|z_1|\) are:

| \(h\) | \(\omega=10\) | \(\omega=40\) | \(\omega=160\) |
|---|---:|---:|---:|
| 3 | 0.0922313419 | 0.0247587021 | 0.00630660838 |
| 4 | 0.0672514677 | 0.0179846130 | 0.00457626487 |

The corresponding \(\bar g\) values agree to the shown digits. This supports the chapter's numerical proposition and its \(O(1/\omega)\) description on the tested set; it is not a proof of uniform convergence at all rapidities.

### The phase-sector exercise

Let \(u=\theta-\theta_d\) and \(\gamma=\pi/h\). Rationalizing the classical delay gives

\[
\arg z_1(u)=\arg\!\left(e^{-2u}+e^{2i\gamma}\right).
\]

At \(u=0\) this is \(\gamma\), as stated. Let the regular incoming soliton be a phase distance \(\epsilon\), with \(0<\epsilon<\gamma\), inside the boundary approached as \(\arg z_1\) increases. Its conversion threshold solves \(\arg z_1(u_*)=\epsilon\):

\[
u_*=-\frac12\log\left(\frac{\sin2\gamma}{\tan\epsilon}-\cos2\gamma\right)<0.
\]

For \(h=3\), \(\epsilon=0.03,0.2,\pi/3\) gives respectively \(u_*=-1.68979706,-0.78140763,0\). As \(\epsilon\downarrow0\), \(u_*\to-\infty\). An exact boundary value makes a Hirota tau function vanish at some spatial point, so it is a singular soliton, not a regular example with a finite threshold. The exercise should use a point inside the edge and discuss the limit. The quantum matrix does not carry this continuous incoming phase parameter.

**Ready-to-paste corrections** in part8-defects/41-soliton-transmission.qmd:

- In line 227 replace “the classically missing weights” by:

  > outgoing charges forbidden by the classical rule for that incoming soliton

- Replace the second sentence of line 262:

  > Using (-@eq-ds-delay-phase), show that $\arg z_1(\theta_d)=\pi/h$. Take a regular incoming soliton whose phase is a distance $0<\epsilon<\pi/h$ inside the charge-sector boundary approached as $\arg z_1$ increases. Show that its classical conversion threshold is $\theta_*-\theta_d=-\frac12\log[\sin(2\pi/h)/\tan\epsilon-\cos(2\pi/h)]<0$. Discuss the limit $\epsilon\to0$ and why the exact boundary is singular. Why can the quantum matrix not see this continuous phase parameter?

## 6. The \(a_2^{(2)}\) construction

**Verdict: correct algebra conversion and numerical intertwiner/triangle/crossing statements at the tested generic parameters, up to an overall scalar; unclear physical scalar dictionary. The cited source's pole locations are internally inconsistent with its displayed Gamma factors.**

**Independent evidence:** [a22_checks.py](../referee_scratch/defects_41/a22_checks.py), [a22_output.txt](../referee_scratch/defects_41/a22_output.txt), [exact_checks.py](../referee_scratch/defects_41/exact_checks.py), [exact_output.txt](../referee_scratch/defects_41/exact_output.txt). The source formula was also inspected visually in the locally rendered [CZ2011 page](../referee_scratch/defects_41/cz2011_p7.png), to exclude a text-extraction error.

### Algebra and module

Starting with the source coproduct and multiplying in the stated order gives

\[
\Delta(K_iX_i^+)
=K_iX_i^+\otimes1+K_i^2\otimes K_iX_i^+,
\]
\[
\Delta(X_i^-K_i^{-1})
=X_i^-K_i^{-1}\otimes K_i^{-2}+1\otimes X_i^-K_i^{-1}.
\]

The Cartan and commutator relations become the standard ones with \(k_i=K_i^2\), symmetrizers \(d_0=4,d_1=1\), and \(q_i=q^{d_i}\). In particular \([e_i,f_i]=(k_i-k_i^{-1})/(q_i-q_i^{-1})\); for unequal nodes the two ordering factors cancel. This checks the conversion, including its ordering.

On the bilateral defect module the converted operators are explicitly

\[
e_1|n\rangle=q^{n+1}|n+1\rangle,\qquad
e_0|n\rangle=q^{-2(n-2)}F(n)F(n-1)|n-2\rangle,
\]
\[
k_1|n\rangle=q^{2n}|n\rangle,\qquad
k_0|n\rangle=q^{-4n}|n\rangle.
\]

I evaluate each word of both Serre polynomials on actual integer labels. There is no cutoff of the operators. For \(q=e^{-0.77i}\), \(n=-3,\ldots,3\), and
\((f_1,f_2)=(1,0.7),(1,0),(0,1),(0.4+0.3i,-1.1)\), the largest normalized residual is \(3.63\times10^{-16}\). The controls \(F(n)=1\) and \(F(n)=(-1)^n\) each have residual 1. This verifies the stated family and controls, not a classification of all possible oscillator functions.

The particle representation has Cartan/commutator residual \(1.81\times10^{-15}\). Its Serre identities also follow directly from \(e_1^3=e_0^2=e_0e_1e_0=0\), and the analogous identities for \(f\).

### Uniqueness and triangle without an artificial charge boundary

Use total charge \(Q=n+m\), where the three particle charges are \(m=1,0,-1\). Each complete sector has dimension 3. The three words with two \(e_1\)'s and one \(e_0\) return \(Q\) to itself. Solving their intertwining equations gives a one-dimensional space at the seed sector; the next singular value is at least 1.4486 in all four examples. Since the \(e_1\) sector map is invertible, its equations transport the seed uniquely to every integer \(Q\).

Testing both generators on the transported sectors gives maximum residual \(3.08\times10^{-15}\). Independently solving the full finite-dimensional particle intertwiner and evaluating the triangle gives at most \(1.06\times10^{-14}\). This is stronger evidence than an arbitrary finite window, while remaining a numerical generic-parameter check. “Unique” here and in the proposition means **unique up to a scalar**, and is not a claim about exceptional spectral parameters.

### Explicit conjugation and the sign test

In the source's three-state basis one convenient normalization is

\[
C(q)=
\begin{pmatrix}
0&0&1\\
0&-q^{-1}&0\\
q^{-2}&0&0
\end{pmatrix},
\qquad x_c=-q^{-6}.
\]

Applying the antipode and transposing gives
\(C\rho_x^*(a)=\rho_{xx_c}(a)C\), with residual \(5.84\times10^{-16}\). The node-1 equation fixes the successive antidiagonal ratios; node 0 then fixes \(x_c\).

The left dual of the crossed representation \(V(x/x_c)\) is \(V(x)\). Thus, in this choice of matrix identifications, the contraction is \(C\,U^{-t_V}C^{-1}\), preserving the defect labels. For output particle index \(j\) and input \(i\), its coefficient is the appropriate inverse-sector entry times \(C_{j,2-j}/C_{i,2-i}\), using indices \(0,1,2\). This also explains why replacing a symmetric conjugation by this antidiagonal one requires care with the transpose in the chapter's index convention.

Allowing the unspecified overall crossing scalar, the residuals in the four examples are \(1.84,2.69,2.26,2.43\), in units of \(10^{-15}\). At \(-q\), \(x_c\) is unchanged and only the middle entry of \(C(q)\) changes sign. Using \(C(-q)\) gives residual at most \(2.69\times10^{-15}\); retaining \(C(q)\) gives errors between 0.9889 and 0.9914. The convention-qualified sign statement is correct. These are algebraic checks and do not determine the physical map from \(x,q\) to rapidity and coupling.

### A source inconsistency which the chapter currently repeats

The two Gamma factors outside the infinite product in [CZ2011 (3.12)](https://arxiv.org/abs/1012.4186) are

\[
\Gamma\!\left(\frac12-z+\frac{\pi}{3\xi}\right)
\Gamma\!\left(\frac12-z+\frac{2\pi}{3\xi}\right),
\qquad z=\frac{i(\theta-\nu)}{\xi}.
\]

Their first poles are therefore, at generic \(\xi\),

\[
\theta_1=\nu-i\left(\frac{\xi}{2}+\frac\pi3\right),\qquad
\theta_2=\nu-i\left(\frac{\xi}{2}+\frac{2\pi}3\right).
\]

The following source paragraph prints \(\xi\) in place of \(\xi/2\). At its printed pole locations the respective outside Gamma argument is \(-1/2\), which is finite; at the locations just derived it is zero. exact_output.txt records both substitutions. The discrepancy is present in the PDF, not just its extracted text.

If the displayed Gamma factors and \(z\) definition are retained, substitution of \(\nu=\eta+i\pi/2\) gives the formal first-strip energy windows \(\pi/3<\xi<4\pi/3\) and \(-\pi/3<\xi<2\pi/3\), respectively, before intersection with the physical coupling domain. These are twice the intervals quoted after (3.13). This local inconsistency does **not** establish the correct complete physical scalar factor: resolving it requires checking that factor against crossing/bootstrap and the physical dictionary. The chapter honestly says it has not checked the quoted windows, but should flag the inconsistency instead of repeating them as a physical result.

The lightest-breather comparison immediately following it does pass a direct substitution: (3.14) becomes

\[
\tanh\left(\frac{\theta-\eta}{2}-\frac{i\pi}{3}\right)
\tanh\left(\frac{\theta-\eta}{2}-\frac{i\pi}{6}\right),
\]

which is exactly the factored form of (3.15).

Source PDF: :codex-file-citation{path="/home/gustav/referee-copy/papers/1012.4186.pdf" purpose="source"}

**Ready-to-paste corrections** in part8-defects/41-soliton-transmission.qmd:

- In line 152 replace “is unique” by “is unique up to a scalar at the tested generic spectral parameters”.

- In line 160 replace the sentence quoting the two instability windows:

  > Corrigan and Zambon discuss two scalar-factor poles approaching the classical absorption poles when $\nu=\eta+i\pi/2$ [@cz2011, eqs. (3.12)–(3.13)]. Their printed pole positions and instability windows are inconsistent with the Gamma factors and the definition of $z$ in (3.12): those factors give $\nu-i(\xi/2+\pi/3)$ and $\nu-i(\xi/2+2\pi/3)$, whereas the subsequent text uses $\xi$ instead of $\xi/2$. We leave the physical windows unresolved pending a check of the complete scalar factor and physical dictionary.

## 7. The square root in CZ2007 (6.14)

**Verdict: correct correction.**

**Independent evidence:** the direct substitution below and the last check in [scalar_checks.py](../referee_scratch/defects_41/scalar_checks.py), with output in [scalar_output.txt](../referee_scratch/defects_41/scalar_output.txt).

With \(q=-e^{-i\pi\gamma}\) and \(\Delta=\eta+i\pi/2\),

\[
\hat x=e^{\gamma(\theta-\Delta)}
=e^{-i\pi\gamma/2}e^{\gamma(\theta-\eta)}
=(-q)^{+1/2}e^{\gamma(\theta-\eta)},
\]

where the square-root branch is explicitly \((-q)^{1/2}=e^{-i\pi\gamma/2}\). It need not be the principal square root at every real \(\gamma\).

Substituting this in both converting entries gives the sine-Gordon matrix immediately preceding [CZ2007 (6.14)](https://arxiv.org/abs/0705.1066). The diagonal entries do not use \(\hat x\) and require no change. The printed inverse square root changes the converting phase to \(e^{+i\pi\gamma/2}\). At \(\gamma=1.37\), the corrected phase residual is zero and the printed phase differs by 1.67161472, before the common real exponential is included.

No correction is needed to chapter 41, line 274.

## 8. Lagrangian caveat, status labels, exercises, and citations

**Verdict: correct and necessary separation between algebraic scattering solutions and a Lagrangian identification; wrong current/charge spin in the anomaly citation; unclear all-rank uniqueness and charged-class scalar claims under the blanket “Established” language.**

**Independent evidence:** the spin-counting calculation below; the exact finite-field rank certificates and root-of-unity boundary example in section 1; the independent class-representative counterexample in section 4; the eight exercise checks below. These use [exact_output.txt](../referee_scratch/defects_41/exact_output.txt), [fusion_output.txt](../referee_scratch/defects_41/fusion_output.txt), and [scalar_output.txt](../referee_scratch/defects_41/scalar_output.txt), rather than the book's numerical pass labels. Source attributions were checked against the named local texts, supplemented by the primary-source page for Robertson.

### What the anomaly source says

The cubic chiral density is a spin-3 **current**. Integrating the corresponding density over a spatial slice gives a charge of spin 2; more generally a conserved current with leading spin \(s+1\) supplies \(Q_s\). The stress tensor/current of spin 2 and energy-momentum charges of spin 1 give the familiar lower-spin check. This is also the book's convention and the source's notation \(Q_{\pm2}\).

[Penati–Polvara, section 5.2](https://arxiv.org/abs/1902.10690) study this spin-3 current for \(A_r\), \(r>1\), and report explicit anomaly-cancellation calculations through \(r=5\). With the finite renormalizations of the derivative coupling, improvement coefficients, and defect potential allowed in their analysis, they find no cancellation. The chapter's “spin-3 charge” should therefore be “spin-3 current, associated with spin-\(\pm2\) charges.” The finite range of their explicit calculations should also be stated.

The source questions quantum integrability; it does not prove that every possible quantum completion or every algebraic transmission matrix is impossible. Conversely, finding a triangle solution does not establish that it is the amplitude of this Lagrangian. The book's independent classical correction concerning the other chirality does not by itself remove a quantum anomaly. These distinctions make the caveat appropriate.

**Ready-to-paste replacement for line 40** of part8-defects/41-soliton-transmission.qmd:

> **The Lagrangian caveat.** This chapter constructs algebraic transmission matrices. Where crossing, bootstrap and scalar factors are supplied, they solve the specified scattering axioms to the stated level of proof or numerical verification. Their identification with the Lagrangian defects of @sec-defect-conditions is a separate conjecture, supported by classical limits but not derived here. Penati and Polvara found quantum anomalies in the spin-3 current, associated with spin-$\pm2$ charges, and reported that their explicit calculations for $2\le r\le5$ admitted no cancellation by the finite renormalizations and improvements considered [@penati_polvara2019, Sec. 5.2]. They consequently questioned the quantum integrability of the Lagrangian defects. This does not invalidate the algebraic matrix solutions, but it prevents treating their Lagrangian interpretation as established.

### Status should attach to individual claims

| Claim | Supported status |
|---|---|
| Type-I\(_-\) vector and conjugate matrix formulas | Established source formulas, independently matched in the same gauge |
| Weyl-module relations and displayed intertwining identities | Direct algebra plus independent checks in the ranks stated |
| Generic uniqueness for \(a_2,a_3\) | Exact generic rank certificates supplied here, conditional on existence; stronger than the chapter's window evidence |
| Uniqueness for every \(r\), every offset, or exceptional parameters | Not established by the supplied arguments; unrestricted root-of-unity version is false |
| Type-I\(_+\) displayed matrix and the undisplayed \(V\) matrix | Algebraic intertwiners and numerical triangle checks as specified; a complete physical scalar construction has not been supplied |
| Crossing and equal-\(k\) scalar factor for type I\(_-\) | Independent algebraic/analytic derivations in sections 2–3, with numerical checks |
| Determinant, generic residue rank, strict first-strip window | Direct derivations; endpoint correction required |
| Kernel restriction and \(X^h\mapsto q^2X^h\) | Independently verified matrix equivalence at the ranks and offsets tested |
| Kernel as the physical excited defect, its chosen scalar, and charged-class resonance assignment | Conjectural; not consequences of the verified kernel invariant |
| \(a_2^{(2)}\) intertwiner, triangle, algebraic crossing | Numerical generic-parameter evidence; physical scalar/dictionary remains unverified |
| Identification with the quantum Lagrangian defects | Conjectural |

Thus “[Established; numerical re-derivation]” is justified for the cited type-I\(_-\) formulas, not for every clause of prp-an-transmission. Use the scoped replacement in section 1. The conventions document already supplies generic \(q\) and fixed branches for \(q^s\); retain that convention and make it visible in the proposition. Even with that assumption, an all-rank uniqueness theorem is missing.

The blanket statement that “everything” is a solution of all scattering axioms should be narrowed: a matrix-only triangle check does not supply its scalar, crossing, or physical bootstrap. The final notes at line 267 likewise need to call the charged-class scalar and resonance interpretation conjectures. A suitable replacement for the ending of that paragraph is:

> The single-label module, orientation analysis, numerical type-I$_+$ construction, offset-independent crossing, sign test, convergent scalar product, residue rank, and algebraic kernel restriction are developed here, with the qualifications stated above. The all-rank uniqueness statement requires a proof; the charged-class scalar assignment and its proposed resonance interpretation remain conjectural.

The drafting report, status document, and conventions document repeat some of the chapter's window and class-assignment claims. They should be updated when corrections are integrated. They were read but have not been edited in this referee pass.

### Exercise checks

All eight exercises were inspected; the following are independent answers or counterexamples, not just consistency checks against the chapter's script.

| Exercise, saved line | Result |
|---|---|
| exr-dt-serre, 234 | Correct for generic \(q\): the factorization \((c-q)(c-q^{-1})\) proves the adjacent relation. At \(r=1\) use the cubic relation, which vanishes for commuting inverse shifts. |
| exr-dt-intertwiner, 238 | Correct intended generic scope; the \(k\)-equations distinguish full charge only away from roots of unity. The exact \(q^5=1\) shift example shows why that hypothesis matters. Generic \(a_2,a_3\) uniqueness is certified in section 1. |
| exr-dt-det, 242 | Correct, using projected weights. The two determinant permutations and the nonzero-component kernel recurrence in section 4 give the requested result, including rank 1 at \(h=2\). |
| exr-dt-crossing, 246 | The difference equation is correct. The rescaling multiplies both crossing scalars by \(\varphi(s)\bar\varphi(s)\). Independence from \(s\) alone fixes this product to an arbitrary constant; it equals 1 only when retaining the displayed scalar normalization. |
| exr-dt-gamma, 250 | Convergence and fusion shifts are correct. The printed-product comparison requires normalization by \(P_K(0)\); section 3 gives the exact finite identity. |
| exr-dt-window, 254 | Needs the strict endpoint, specified lower strip, and conditional charged-class scalar. The corrected window table is in section 4. Solving \(a+2p\equiv h/2\pmod{2h}\), with \(0\le p<h\), gives \((a,p)=(2,0)\) for \(a_3\), and \((1,1),(3,0),(5,5)\) for \(a_5\). The last pair uses a different representative from \(0<a+2p\le h\) and must not be lost when listing classical real absorptions. |
| exr-dt-qsign, 258 | Correct in the fixed conventions. The recurrence \(c_{i+1}=-c_i\) closes iff \(h\) is even. For \(h=2\), composing the diagonal grading with the identification of the conjugate doublet by \(\sigma_x\) produces an antisymmetric conjugation, up to normalization. |
| exr-dt-semiclassical, 262 | Prefactor and central phase are correct. A regular soliton must be inside the sector boundary; the exact edge is singular. Section 5 supplies the corrected threshold calculation. |

**Additional ready-to-paste correction at line 246:**

> Show that $s$-independent crossing scalars require $\varphi(s)\bar\varphi(s)$ to be constant. If the displayed crossing scalars are retained with their given normalization, this constant is 1.

### Citation audit

The table covers every distinct bibliography key cited in the saved chapter and the claims attached to it. Line references include the notes and source-correction list; repeated uses of a source are grouped.

| Source and chapter use | What was checked; verdict |
|---|---|
| [CZ2009a, 0902.1307](https://arxiv.org/abs/0902.1307), lines 64, 71, 74, 80–133, 171–194, 267, 271–272, 275 | Equations (5.8), (5.11) give the displayed type-I\(_-\) matrices under the explicit dictionary. The non-unitarity remark, selection-rule asymmetry, crossing (6.8), fusion (6.6)–(6.7), and species pole (6.13) support their cited uses. The correction to the text after (6.15) is valid, but its zero-width endpoint also needs correcting. The product drift is valid after normalization as detailed above. The source does not prove the chapter's new all-rank Weyl-module uniqueness assertion. |
| [CZ2007, 0705.1066](https://arxiv.org/abs/0705.1066), lines 36, 74, 185, 229, 267, 273–274 | The transpose/species dictionary follows from (2.12), (2.20). The tunnelling interpretation is in section 6. Their instability interval is strictly \(1/2<\gamma<2\), so the chapter's closed endpoint is not a faithful restatement. Equation (6.14) has the corrected square-root sign verified in section 7. Section 7 supports the lightest-breather/classical-factor comparison with its parameter identification; it is not independent verification of the later chapter's all-rank breather result. |
| [CZ2011, 1012.4186](https://arxiv.org/abs/1012.4186), lines 34, 137–160, 267 | Equations (5.7)–(5.9) support the \(a_2^{(1)}\) oscillator discussion; (4.1)–(4.10) support the Tzitzéica algebra/module conversion. The stated family without \((-1)^N\) is a subset of their more general family. Sections 3–4 support the two construction routes. The source really prints the quoted windows, but its Gamma/pole inconsistency must be flagged. The substitution (3.14)\(\to\)(3.15) is correct. |
| [Penati–Polvara, 1902.10690](https://arxiv.org/abs/1902.10690), lines 40, 267 | The quantum-integrability reservation is supported. Section 5.2 concerns a spin-3 current and spin-2 charge, with explicit calculations through \(r=5\); the present wording confuses current and charge and conceals that extent. |
| [CZ2010, 1006.0939](https://arxiv.org/abs/1006.0939), line 76 | Section 4.1 exhibits the sine-Gordon soliton S-matrix as a finite block of a specially parametrized type-II transmission matrix, including the matching scalar. This supports the analogy with fused defects. It does not settle the chapter's open \(a_r\) type-I\(_-\)/type-I\(_+\) question. |
| [Robertson, 1408.1960](https://arxiv.org/abs/1408.1960), lines 76, 267 | The primary-source abstract explicitly describes classical and quantum fusing of \(a_r^{(1)}\) defects and a new species-two transmission matrix in \(a_3^{(1)}\). This supports the broad attribution. Neither this report nor the chapter verifies that additional matrix. |
| [Weston, 1006.1555](https://arxiv.org/abs/1006.1555), line 267 | Equations (4.15), (4.16) explicitly connect type-I and type-II sine-Gordon transmission matrices to the oscillator \(L\)-operator. The historical attribution is supported; it is not an all-rank uniqueness theorem for this chapter. |
| [Bowcock–Corrigan–Zambon, hep-th/0506169](https://arxiv.org/abs/hep-th/0506169), lines 133, 267 | Equation (3.51) gives the sine-Gordon matrix, supporting the attribution. The differently indexed Gamma factors in (3.42), (3.47) support the analogous product-normalization issue; the exact calculation required for the present chapter is supplied independently in section 3. |
| [Konik–LeClair, hep-th/9703085](https://arxiv.org/abs/hep-th/9703085), line 267 | Equations (5.21), (5.25) give the sine-Gordon transmission-amplitude ansatz and minimal solution. The historical attribution is supported. |

Internal cross-references were checked against the specified prerequisite conventions, crossing/triangle/bootstrap equations, classical delay, and class/charge statements. This is not a full independent review of every theorem in those other chapters. In particular the final lightest-breather generalization points forward to a separate chapter and is not established by the cited \(a_2\) source alone.

## Other observations — unverified

- A uniform all-rank generic uniqueness proof may be obtainable from a universal-R or representation-theoretic construction. The finite-field certificates here cover only the explicitly stated ranks and do not supply that theorem.
- A consistent Klein-type equivalence between the two signs of \(q\) for odd \(h\), including mixed amplitudes, fusion, and scalar factors, has not been constructed or ruled out beyond the fixed-convention obstruction.
- The physical bootstrap on the \((h-1)\)-dimensional residue image, its scalar factor, and a possible fused-defect interpretation remain open in this review.
- Robertson's additional \(a_3\) matrix, the full physical Tzitzéica scalar/dictionary, and the later chapter's universal lightest-breather claims have not been independently checked here.
- This report verifies the cited attributions and computations needed for the brief; it is not an exhaustive priority search for every construction described as new.
