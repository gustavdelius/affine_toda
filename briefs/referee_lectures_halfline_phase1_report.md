# Half-line lectures: blind phase-1 referee report

- Date: 2026-10-25.
- Referee: Codex (GPT-6; exact runtime model identifier not exposed in this session).
- Review copy: `/home/gustav/referee-copy/referee_lectures_halfline-codex`.
- Commit: the brief identifies the lecture snapshot as `035351d`. The required `git rev-parse HEAD` was attempted and returned `fatal: not a git repository (or any parent up to mount point /home)`. This exported copy has no `.git`; I cannot independently certify its commit hash.
- Scope: all six lectures and preface read in full before this risk list was written. No forbidden earlier reports or phase-2 list read. No subagents, changes to existing source files, Git mutations, rendering, or author figure-script execution.
- Locations below use the line numbers in this copy. L1–L6 abbreviate the six numbered `.qmd` files in `lecture-notes/half-line/`; Index means `index.qmd`.

## 1. Risk list (fixed before checking)

The following order was recorded before independent calculations, source comparisons, or examination of author scripts. Entries will not be reordered or deleted.

1. **L3:224–243, `eq-lec3-cj-from-chi`, figure caption at 175.** Claim: “$C_j = \epsilon_{j+1}\,\epsilon_{j-1}$”, with $\chi\in(-\pi/3,\pi/3)$ giving $(-1,+1,-1)$. Risk: square-root branches and cyclic versus antiperiodic cosine indices can change boundary signs. Check the real field, its derivative and the boundary potential directly; enumerate all eight sign patterns.
2. **L6:195–209, rotor callout.** Claim: “Computing the classical kinetic moment of inertia ... gives $I = \frac{6}{\beta^2 m_1}$”, and a “full $2\pi$ circle” connects the three patterns. Risk: the classical family has singular endpoints and fixed-boundary parameter intervals; a coordinate-dependent metric need not describe a compact rotor. Integrate $|\partial_\chi\phi|^2$ and inspect endpoint distances and boundary signs; compare the book's status.
3. **L1:464–485, `eq-lec1-ly-shift`.** Claim: “$K_b^{(1)}(\theta)=K_{b+2}(\theta)$”. Risk: the last block cancellation does not visibly follow. Multiply the defined blocks, find the actual parameter transformation, and read the source.
4. **L1:369–404 and Exercise 1.3:535.** Claim: “$S_{\rm LY}(2\theta)=-(2)(4)$” and “$f_x(\theta)=(2x)(6-2x)$”. Risk: doubled rapidity, block scale and overall signs are mixed in the advertised step-by-step proof. Derive the doubling identity and independently test the corrected bootstrap.
5. **L2:323–327, full Dirichlet obstruction.** Claim: the Cartan part “forces $[K(\lambda),\vec H]=0$”, hence $K$ diagonal. Risk: the compatibility equation has a plus sign and may force an anticommutator instead. Derive the condition and solve it on the defining representation; delimit the no-go scope.
6. **L6:214–227, Coleman–Thun example.** Claim: “$5-2(2)-1=1$”, residue $-1.1875$, and an identity verified “to machine precision”. Risk: topological pole counting, rapidity/block-label derivatives and residue normalizations may disagree. Count pole orders from the printed factors, compute residues with epsilon scaling, and inspect diagram/source and script scope.
7. **L6:91–99, imaginary-to-real continuation.** Claim: “$(\frac12-\frac h2\omega)_I\to(-\frac B2+h)$”. Risk: direct substitution into $2a/\omega$ seems sensitive to a lost minus sign. Perform the substitution in all three factors and compare full analytic functions.
8. **L5:155–196 and Exercise `exr-lec5-serre-ej`.** Claims: cross terms cancel and “$\kappa=k_i^2k_j=F_{ij}(b_i,b_j)-b_j$” can be used for conjugation inside the coideal. Risk: coefficient omission, unauthorized inversion in a subalgebra, and hidden noncommutative cancellations. Independently expand the Serre polynomial and distinguish an algebra identity from an intertwiner obstruction proof.
9. **L5:124–131 and Exercise `exr-lec5-cyclic-unitarity`.** Claims: $K_s(-\theta)^T K_s(\theta)=U_s1$ and “$Z^T=sZ^{-1}$”. Risk: twisted cyclic shifts, branches of fractional powers and scalar normalization. Directly construct matrices at generic irrational coupling and derive the eigenvalue cancellation.
10. **L3:116–127 and Exercise `exr-lec3-mirror-square`.** Claims: “$k=m_a\cosh\theta=2ms$”, displayed derivative at 298, and the perfect-square proof of the uniform boundary condition. Risk: missing powers of $\cosh\theta$, incorrect common denominator and vacuum shift. Differentiate the tau functions and prove the boundary identity including its square-root branch.
11. **L1:157–172 and Exercise `exr-lec1-twisted`.** Claim: both twisted reflection-equation sides “cleanly map” the specified tensor products. Risk: one mixed braid's subscripts may be reversed; prose calls the second tensor slot the first. Trace each domain and codomain in spatial order, and derive the fixed-slot/braid conversion.
12. **L4:295–310 and Exercise `exr-lec4-1d-submodule`.** Claim: the common eigenvector exists “precisely when” $y=i(-1)^n e^{-i\eta}$, with printed $\beta_j$. Risk: generalized versus ordinary eigenvalue equations and additional reducibility branches. Solve the two quadratic eigenvector equations independently.
13. **L6:254–268, Exercise `exr-lec6-kim-oneloop`.** Claims: the block derivative, alternative formula for $K_1^N$, and the quoted one-loop coefficient. Risk: a wrong block identity can make the exercise unsolvable. Differentiate the blocks symbolically and numerically, compare the resulting coefficient with the source.
14. **L6:237–249, Exercise `exr-lec6-bootstrap-a2`.** Claim: $u_{11}^2=\pi/3$ for equal-mass $a_2$ particles. Risk: fusion angle confused with constituent shifts. Derive from the mass triangle and independently check the bootstrap identity.
15. **L3:266–279, `eq-lec3-linear-reflection`.** Claim: for $C\in\{0,+1,-1\}$, $K^{\rm cl}=-[(a)(h-a)]^{-C}$. Risk: the block expression may fail at $C=0$. Solve the Robin boundary equation and check each value.
16. **L4:257–263, sign of $q$ and figure caption.** Claim: $q\mapsto-q$ conjugates the bulk braid matrix by $\Pi\otimes1$; wrong-sign crossing has an $O(1)$ error checked by the script. Risk: spectral reparametrization and Klein twist are not ordinary conjugation; the script may not test the stated matrix crossing equation. Construct a six-vertex intertwiner from the coproduct and inspect exactly what is tested.
17. **Index:45–49; L5:230–251; L6:117–151.** Claim: the Neumann/uniform puzzle is resolved and the amplitudes are exact. Risk: a conjectural physical identification or restricted fusion check may have been silently promoted. Compare book labels and source claims; independently check the SG UV–IR Taylor series and formula-level duality, separating them from identification of boundaries.
18. **L6:48–83,119–125, singlet fusion.** Claims: an affine “singlet” and Gamma products telescope to the stated three blocks, with $c_+(\mu)$ verified by the lecture script. Risk: finite versus affine singlet, reversed fusion rapidities, unprovided scalar factors and branch normalization. Check representation domains and scalar formula against the book/source; inspect script coverage and independently evaluate any fully specified identity.
19. **L2:227–311, classical matrices and rigidity proof.** Claim: explicit $a_1$ and $a_2$ matrices solve all node equations; one bad spectator determinant forces the Neumann branch. Risk: mixed parameter branches are not explicitly eliminated. Multiply the matrices and prove propagation, including generic spectral invertibility.
20. **L1:20–35,188–202 and figure caption.** Claims: continuity equations convert to the printed sum/difference law; $(x,t_E)\mapsto(-t_E,-x)$ is a rotation; $\tilde E=-ip,\tilde p=-iE$ implies $\theta=i\pi/2-\tilde\theta$. Risk: independent sign/orientation errors at foundational steps. Apply coordinate differentiation and elementary hyperbolic identities.
21. **Index:27; L1:206–215,278–293.** Claims: the boundary state “can only emit” opposite-rapidity pairs; fixed versus moving poles diagnose one-point versus bound-state origins. Risk: zero-momentum one-particle terms and Coleman–Thun alternatives are omitted. Apply charge constraints to a zero-rapidity particle and read the boundary-state/source qualifications.
22. **L1:353–357,414–452.** Claims: GZ printed the opposite Lee–Yang sign and $\Gamma^2=-2\sqrt3(2+\sqrt3)$. Risk: corrections to literature and residue normalizations require source-level evidence. Read the exact GZ line, compute bulk and boundary residues directly and compare the book.
23. **L3:211–222 and Exercises `exr-lec3-breathers`, `exr-lec3-self-conjugate`.** Claims: both quadratic roots lie at $u\ge1$; regular breathers and a self-conjugate zero-energy family follow. Risk: negative/complex roots, endpoint singularities and $s=1$ degeneracy. Factor the lower bound, check boundary regularity and directly calculate the self-conjugate energy.
24. **L6:175–203,273–284 and spectrum figure.** Claim: the excited factor formula describes $b_{n,m}$ for $n+m\ge0$, including the rotor and exact degeneracy at $B=1/2$. Risk: state-existence ranges and cancellations are not given. Verify energy sums, ground-state specialization and degeneracy; distinguish analytic continuation from physical states.
25. **Index:69 and all six script/figure claims.** Claim: “Every non-trivial calculation and figure ... verifies the algebra ... from first principles.” Risk: scripts may test end results or tautologies while missing the printed steps; figures can depict canceled poles or inadmissible states. Read (do not execute) the six scripts, inspect the six committed figures, and map actual assertions to lecture claims using the independent checks above.

## 2. Checks (in rank order)

### 2.1. Risk 1 — vacuum sign assignments

**Verdict: wrong (two interval assignments and the figure caption); the displayed cosine-product formula is correct with extended integer indices.**

Independently used roots $\alpha_0=(-1,0,1)$, $\alpha_1=(1,-1,0)$, $\alpha_2=(0,1,-1)$, $m=\beta=1$, and the printed tau functions. For each of the eight patterns I evaluated
\[
 r_C=\phi'(0)+\sum_{j=0}^2C_j\alpha_j\exp(\alpha_j\cdot\phi(0)/2),\quad
 \phi=-\sum_j\alpha_j\log\tau_j,\quad
 \phi'=-\sum_j\alpha_j\tau'_j/\tau_j.
\]
`referee_scratch/lectures_halfline/check01_vacuum.py` gives:

| $\chi$ | unique pattern with $r_C=0$ | residual of the printed pattern |
|---|---|---|
| 0.25 | $(+,-,-)$ | 11.6048381456 |
| 2 | $(-,-,+)$ | $2.23\,10^{-15}$ |
| −2 | $(-,+,-)$ | 4.82785507479 |

The cosine signs are **antiperiodic**, $\epsilon_{j+h}=(-1)^a\epsilon_j$, not periodic. At small $\chi$, $\epsilon_{-1}=+$ whereas $\epsilon_2=-$. Substituting $\epsilon_2$ for $\epsilon_{-1}$ produces the erroneous first row. The product over all $C_j$ remains $+1$: the two wrapping factors each contribute $(-1)^a$.

**Replacement, L3:229–243:** retain the cosine-product equation but add:

> In this formula the cosine is evaluated at the integer indices $j-1$ and $j+1$, before reducing any index modulo $h$. Thus $\epsilon_{j+h}=(-1)^a\epsilon_j$. The two wraparound signs cancel in $\prod_j C_j$, giving $+1$. For $h=3$, $a=1$, the three intervals give
> \[
> (C_0,C_1,C_2)=
> \begin{cases}
> (-1,+1,-1),&-\pi<\chi<-\pi/3,\\
> (+1,-1,-1),&-\pi/3<\chi<\pi/3,\\
> (-1,-1,+1),&\pi/3<\chi<\pi.
> \end{cases}
> \]

At L3:175 replace `$(C_0,C_1,C_2)=(-1,+1,-1)$` by `$(C_0,C_1,C_2)=(+1,-1,-1)$`. The book's `prp-solitonic-vacua` does not give the erroneous interval table.

### 2.2. Risk 2 — moment of inertia and rotor interpretation

**Verdict: correct for the moment of inertia and leading energy match; unclear for the asserted global quantum rotor.**

From the lecture definitions, putting $z=e^{m_1x}$, $s=\sqrt3/2$ and $v_j=\chi+2\pi j/3$ gives
\[
 \partial_\chi\phi=-\beta^{-1}\sum_j\alpha_j
 \frac{-2z\sin v_j}{(1+s)\tau_j},\qquad
 I=\frac1{m_1}\int_0^1\frac{|\partial_\chi\phi|^2}{z}\,dz.
\]
Independent quadrature in `check02_metric.py` gives $\beta^2m_1I=6.00000000000$ at $\chi=0,0.25,0.6,1$, and $5.99999999978$ at $\pi/3-0.001$. Near an endpoint the integrand has a very narrow boundary layer; unpartitioned quadrature initially missed it. The retained script uses a positive, cancellation-resistant tau expression and split integration. Thus this check **does not** refute the inertia claim.

Also $F(n)+F(-n)=\sin^2(n\delta)/\sin\delta=n^2\delta+O(\delta^3)$, and $\delta=\beta^2/12+O(\beta^4)$, so $I=6/(\beta^2m_1)$ is precisely the required coefficient.

However, at fixed boundary potential the regular coordinate interval has length $2\pi/3$. At its endpoints a tau function vanishes at the wall; continuing changes the boundary potential (check 2.1). Its metric length is finite, $2\pi\sqrt I/3$. The same local kinetic operator on this interval has different spectra for different self-adjoint endpoint conditions; imposing Dirichlet endpoints, for example, gives $E_k=(\pi k/(2\pi/3))^2/(2I)$, not the proposed circle spectrum. Equal reflection amplitudes for different boundary theories do not themselves prescribe gluing their configuration spaces. A quantum identification and endpoint conditions remain to be supplied.

**Replacement, L6:206–210:**

> The classical modulus has metric $I=6/(\beta^2m_1)$, and the small-$B$ diagonal tower has energies $n^2/(2I)$. This is evidence for a rotor interpretation. For a fixed sign pattern, however, the regular classical family occupies only an interval of length $2\pi/3$, with singular endpoints. A rotor of period $2\pi$ would require a quantum identification of the three intervals and appropriate endpoint conditions. We have not derived that identification.

**Book finding:** the inertia statement at `part7-half-line/36-particle-reflection.qmd:282` survives; its inference “So the rotor must run through all three quantum-equivalent patterns” needs the same qualification. The book does call this an interpretation; the lecture incorrectly makes it compulsory.

### 2.3. Risk 3 — Lee–Yang excited-state shift

**Verdict: correct.**

Let $k_b=K_b/K_{(\mathbb1)}$. Applying the shift rule to $S(\theta\pm iv_1)$ gives
\[
 S(\theta-iv_1)S(\theta+iv_1)
 =(3+b)(1-b)(5+b)(3-b).
\]
Multiplying by $k_b$ cancels $(1-b)$ and $(5+b)$ and leaves
\[
 (1+b)(3+b)(3-b)(5-b)=k_{b+2}.
\]
Here $(-x)=(x)^{-1}$ and $(12-x)=(x)^{-1}$ have been used in $k_{b+2}$. It also equals $k_{4-b}$, since $K_p=K_{6-p}$. Thus the alternative parametrization $4-b$ is **not** a counterexample to the printed $b+2$.

`check03_04_ly.py`, at $\theta=0.71$, $b=0.37$, gives residuals $6.75\,10^{-16}$ for $b+2$ and $4.04\,10^{-16}$ for $4-b$. [Dorey–Pocklington–Tateo–Watts](https://arxiv.org/abs/hep-th/9712197), source `sources/9712197/Paper4.tex`, equation labelled `eq:rfb`, prints exactly the family
`\usbl{\fract{1-b}{2}}^{-1}\!\usbl{\fract{1+b}{2}}\usbl{\fract{5-b}{2}}\usbl{\fract{5+b}{2}}^{-1}`
relative to $R_{(1)}$ (their block arguments are half ours). No correction to the shift is needed.

### 2.4. Risk 4 — Lee–Yang doubling proof and block exercise

**Verdict: wrong, despite the correct final identity-boundary amplitude.**

Direct factorization of $\sinh(2z)$ gives, with blocks of denominator $h$,
\[
 (x)(2\theta)=-(x/2)(\theta)/(h-x/2)(\theta),\qquad
 f_z(\theta)=-(hz)(h-hz).
\]
Consequently, at $h=6$,
\[
 S_{\rm LY}(2\theta)=-\frac{(1)(2)}{(4)(5)},
\]
not $-(2)(4)$ at argument $\theta$. At $\theta=0.71$ the printed equality has absolute error $0.824426363164$; the corrected one has error $4.58\,10^{-16}$ (`check03_04_ly.py`). The later display identifying $(y)(2\theta)$ with $f_{y/6}(\theta)$ is also false. The displayed factorization at L1:385 has RHS $\sinh(2\theta)+i\cos\alpha$, so is valid, but does not repair the preceding false identities. The abandoned ellipses and self-corrections are not a proof.

**Replace L1:369–403 by:**

> The doubling identity for our $h=6$ blocks is
> \[
> (x)(2\theta)=-\frac{(x/2)(\theta)}{(6-x/2)(\theta)}.
> \]
> Hence $S_{\rm LY}(2\theta)=-(1)(2)/((4)(5))$. Since $(6)=-1$ and $(-1)(1)=1$, the shifted product is
> \[
> K(\theta+i\pi/3)K(\theta-i\pi/3)=-\epsilon^2\frac{(3)(5)}{(2)}.
> \]
> Multiplying by $S_{\rm LY}(2\theta)$ gives $\epsilon^2(1)(3)/(4)$. The bootstrap therefore requires $\epsilon=\epsilon^2$, and the nonzero solution is $\epsilon=+1$.

**L1:535:** replace `$f_x(\theta)=(2x)(6-2x)$` by `$f_x(\theta)=-(6x)(6-6x)$`. The ratio requested in the exercise is correct because the two minus signs cancel. The positive identity amplitude independently satisfies its bootstrap to $7.02\,10^{-16}$.

### 2.5. Risk 5 — Dirichlet argument

**Verdict: wrong derivation and overstated scope; correct exclusion for the stated $a_n$ defining-representation laboratory.**

For $\sigma=T$, the free Cartan term in $Ka_t+a_t(1/\lambda)^TK$ is proportional to
\[
 KH+HK=\{K,H\},
\]
not $[K,H]$. On the defining $sl_N$ module, write $H=\operatorname{diag}(h_1,\ldots,h_N)$ with $\sum h_i=0$. Then $(h_a+h_b)K_{ab}=0$ for every such $H$. If $N\ge3$, for each $(a,b)$ there is an $H$ with $h_a+h_b\ne0$, so every $K_{ab}=0$. For $N=2$ the off-diagonal entries survive; this is exactly the antidiagonal Dirichlet limit in the preceding paragraph. This independent derivation establishes the exclusion without invoking the spectral-cycle argument.

**Replace L2:327 by:**

> With $\sigma=T$, the unconstrained Cartan term requires $\{K,H\}=0$ for every Cartan generator. In the defining representation of $a_n$, $n\ge2$, this implies $K=0$: its weights contain no opposite pair. For $a_1$ the opposite weights allow an antidiagonal $K$, as in the Dirichlet limit above. Thus full Dirichlet has no nonzero field-independent $K$ in our $a_n$ laboratory. The book's more general no-go statement specifies untwisted algebras, finite-dimensional evaluation representations without trivial summands, and the allowed anti-automorphisms; it does not rule out every conceivable construction.

The book's `prp-no-dirichlet` does state these qualifications and handles the dual representation correctly. Do not replace its proof with the lecture argument.

### 2.6. Risk 6 — Coleman–Thun counting and residue

**Verdict: wrong pole-counting arithmetic; correct residue and residue identity; the separate propagator/loop counts are unclear.**

The printed eight-block amplitude has exactly one uncanceled denominator zero at $p=47/40$ for $B=9/20$. The intermediate factor has one numerator zero at $p_R=7/40$. My independent 65-digit calculation `check06_ct.py` evaluates residues from $\varepsilon K(\theta_p+\varepsilon)$ and the zero derivative from $K(\theta_R+\varepsilon)/\varepsilon$, using $\varepsilon=10^{-10},10^{-20}$:

- $-i\operatorname{Res}K_1^{b_{-2,2}}=-1.1874897238607920702$;
- $(K_2^{b_{-1,1}})'(\theta_R)=-2.2114724331825i$;
- discrepancy in the stated residue identity: $1.45\,10^{-9}$ then $1.45\,10^{-19}$.

Here $g^2=2\operatorname{Res}K_1^{b_{-1,1}}(33i\pi/120)/i$ and $\Gamma^2=\operatorname{Res}S_{11}(2i\pi/3)/i$, computed independently from the blocks. The convergence establishes the orders and normalization, rather than detecting a pole by a residue threshold.

But $5-4-1=0$, not $1$. If the loop count $L=2$ and formula $P-2L-Z$ are retained, a naive double pole would require $P=6$. That is an algebraic consistency requirement, not an independent verification of the topology. The lecture supplies no drawing with the internal segments and loop variables identified; I cannot certify whether its $P$, $L$, or the application of the counting rule needs changing. The independently checked scalar factors do establish a simple pole and the stated residue identity.

**L6:223 replacement sentence:**

> The reflection factor has a simple pole, and the intermediate reflection factor has a simple zero. The residue identity below supports the proposed Coleman–Thun explanation. A diagram with its internal propagators and integration variables identified is needed to justify a separate topological pole count.

Retain the residue identity and the rounded value $-1.1875$. This numerical identity is evidence for the proposed explanation, not a general theorem assigning every negative-residue pole to this diagram.

### 2.7. Risk 7 — continuation of imaginary blocks

**Verdict: correct, with a missing explanation of periodic reduction.**

The three labels $2a/\omega$ after continuation are
\[
 -n,\quad -h-B/2,\quad -1+B/2.
\]
The middle label is not literally $h-B/2$, but it is equivalent modulo the block period $2h$. Thus all three printed final blocks and $K_1$ are correct. `check07_09_algebra.py`, $h=3,B=0.47,\theta=0.63+0.11i$, obtains labels $(-2,-3.235,-0.765)$ and a full-product residual $1.27\,10^{-16}$.

**Clarification at L6:94:** write `$(-h-B/2)=(h-B/2)$ by the period $2h$ of the block`. No change to the amplitude is warranted.

### 2.8. Risk 8 — Serre identity versus obstruction argument

**Verdict: correct Serre identity; wrong equality and an unproved subalgebra step in the following argument.**

Independently separating the weight-zero terms gives
\[
 (1-q^2)e_i(c_if_ik_i)+(1-q^{-2})(c_if_ik_i)e_i
 =c_i(1-q^2)[e_i,f_i]k_i=-qc_i(k_i^2-1).
\]
The $e_j$ mixed terms reduce to $-q^2D_ie_j+e_jD_i=qc_ie_j$, since $k_i^2e_j=q^{-2}e_jk_i^2$; the $\varphi_j$ terms reduce to $D_i\varphi_j-q^2\varphi_jD_i=qc_i\varphi_j$. Together with the $k_i^2$ coefficient $2-q-q^{-1}$ this proves the displayed identity. An independent coproduct implementation on dimensions 3, 9 and 27, with generic irrational physical $q$ and arbitrary complex $c_i,\hat\epsilon_i$, gives maximum errors $1.24,2.29,3.27$ times $10^{-16}$ (`check07_09_algebra.py`).

At $c_i=q^{-1}$, define $d_{ij}=\hat\epsilon_j((1-q)(1-q^{-1})\hat\epsilon_i^2-1)$. The actual consequence is
\[
 k_i^2k_j=(F_{ij}(b_i,b_j)-b_j)/d_{ij}\in B\quad(d_{ij}\ne0).
\]
The lecture omits $d_{ij}^{-1}$. More seriously, membership of an invertible ambient element $\kappa$ in a unital subalgebra does **not** imply $\kappa^{-1}\in B$. The claimed conjugations inside $B$ need justification. The elementary counterexample is $\mathbb C[t]\subset\mathbb C[t,t^{-1}]$: $t$ belongs to the smaller algebra, $t^{-1}$ does not.

The intertwiner conclusion can be proved without claiming that ambient inverse belongs to $B$. In the $a_2$ laboratory, for $(i,j)=(1,2)$ the two represented Cartan elements have spectra
\[
 \{q^2,q^{-1},q^{-1}\},\qquad\{q^{-2},q,q\}.
\]
They are disjoint at generic non-root-of-unity $q$, so an intertwiner between them is zero. Alternatively, on any fixed finite-dimensional invertible representation, its inverse is a polynomial in its matrix by Cayley–Hamilton; the conjugation argument works **in the represented algebra**, without establishing the universal membership claim.

**Replace L5:194–196 by the formula for $d_{ij}$ and the preceding $a_2$ spectral argument.** For an all-rank statement, explicitly use the represented-algebra Cayley–Hamilton argument before the Vandermonde separation. This is ready to paste as:

> If $d_{ij}\ne0$, then $\kappa=(F_{ij}(b_i,b_j)-b_j)/d_{ij}=k_i^2k_j$ belongs to $B$. On a fixed finite-dimensional representation, $\pi(\kappa)^{-1}$ is a polynomial in the invertible matrix $\pi(\kappa)$. Thus the conjugation and separation argument applies to the represented algebra. It does not by itself prove that $\kappa^{-1}$ belongs to the universal coideal.

**Book finding:** the same unjustified universal inverse step occurs at `part7-half-line/35-quantum-symmetric-pairs.qmd:158`. The representation-level repair preserves the proposed no-intertwiner result.

### 2.9. Risk 9 — twisted cyclic unitarity exercise

**Verdict: correct matrix and unitarity factor; wrong transpose identity in the exercise.**

For $Z=\left(\begin{smallmatrix}0&1&0\\0&0&1\\s&0&0\end{smallmatrix}\right)$ and $s=\pm1$, direct multiplication gives $Z^TZ=1$, so $Z^T=Z^{-1}$, not $sZ^{-1}$. The latter fails by entries of magnitude 2 when $s=-1$.

With $t=X\zeta$, $t^3=sX^3$, the polynomial eigenvalue is
\[
 \kappa_s(X)+X^{-1/2}\zeta+X^{1/2}\zeta^2.
\]
Combining it over the common denominator and using $1+t+t^2=(t^3-1)/(t-1)$ gives the rational eigenvalue printed at L5:296. For the reflected matrix, $X'=w^2/X$, $\zeta'=\zeta^{-1}$, and $t'=w^2/t=-q/t$. Therefore $t'+q=q(t-1)/t$ and $t'-1=-(t+q)/t$; these cancel the eigenvector dependence. The result is the printed $U_s$.

Using continuous powers $X^p=\exp[p\omega(\theta-i\pi/2)]$, `check07_09_algebra.py` gives unitarity errors $4.58\,10^{-16}$ and $9.16\,10^{-16}$ at $\omega=\sqrt2,\theta=0.43$ for $s=+1,-1$ respectively.

**L5:298:** replace `Since $Z^T=sZ^{-1}$` by `Since $Z^T=Z^{-1}$`. The subsequent substitution $\zeta\mapsto\zeta^{-1}$ was already correct.

### 2.10. Risk 10 — perfect-square mirror calculation

**Verdict: wrong intermediate equations; correct mirror shift after repair.**

Let $a_* =\sin(\pi a/h)$, $D=1+Cs$, $z_j=\zeta^{ja}U$, and $c_j=\cos(\psi+\pi ja/h)$ with $U=e^{2i\psi}$. Direct differentiation, before making any boundary simplification, gives
\[
 \tau_j(0)=\frac{2+z_j+z_j^{-1}}D=\frac{4c_j^2}D,\quad
 \tau'_j(0)=\frac{k}{D}\big[z_j+z_j^{-1}+2(1-Cs)\big].
\]
There is no extra $D$ in the last term. Also
\[
 k=m_a\cosh\theta=2m a_*\cosh\theta=\frac{2m a_*^2}{s},
\]
not either $2ms$ (L3:127) or $2ma_*/s$ (L3:305).

Taking the algebraic square roots consistently gives $\sqrt{\tau_{j+1}\tau_{j-1}}=4(c_j^2-a_*^2)/D$. Consequently
\[
 \frac{\tau'_j-mC\sqrt{\tau_{j+1}\tau_{j-1}}}{\tau_j}
 =k-mC+\frac{C}{c_j^2}\left(ma_*^2-\frac{ks}{2}\right)=k-mC.
\]
This is the missing three-line proof, with its necessary branch qualification. At imaginary coupling it is relative to the initial boundary vacuum, as L3:116 says; using arbitrary principal square roots is not equivalent.

`check10_mirror.py` at $h=3,a=1,\theta=0.6,U=e^{0.2+0.7i}$ finds errors $0.26994$ and $71.4993$ in the printed derivative for $C=+1,-1$; the corrected derivative and constant ratios agree within $10^{-14}$.

**Ready-to-paste replacements:** L3:127 and 305 use `$k=m_a\cosh\theta=2m\sin^2(\pi a/h)/s$`. At L3:298 use
\[
 \partial_x\tau_j(0,t)=\frac{k}{1+Cs}\left(\zeta^{ja}U+\zeta^{-ja}U^{-1}+2(1-Cs)\right).
\]
Use the constant-ratio computation above after that equation. Also, L3:161 should factor the dominant terms as `$\zeta^{ja}E_1(1+\zeta^{-ja}AE_2)$`, not with $\zeta^{-2ja}$; this does not change the location $AE_2=1$.

### 2.11. Risk 11 — reflection-equation tensor domains

**Verdict: wrong.**

The RHS of `eq-lec1-reflection-twisted`, read right to left, must have domains
\[
 V_\mu\otimes V_\nu
 \xrightarrow{\check R_{\mu\nu}}V_\nu\otimes V_\mu
 \xrightarrow{1_\nu\otimes K_\mu}V_\nu\otimes V_{\bar\mu}
 \xrightarrow{\check R_{\nu\bar\mu}}V_{\bar\mu}\otimes V_\nu
 \xrightarrow{1_{\bar\mu}\otimes K_\nu}V_{\bar\mu}\otimes V_{\bar\nu}.
\]
The printed third operator $\check R_{\bar\mu,\nu}$ instead requires $V_{\bar\mu}\otimes V_\nu$ as its input and cannot be applied. The LHS is typed correctly. L4:148 and the book's twisted equation have the correct RHS mixed braid.

There is a separate conversion error. Substituting $\check S=S_{12}P=PS_{21}$ in the second-slot braid equation gives the fixed-slot equation term by term, because
\[
 (S_{12}(d)P)K_2(\theta_1)(S_{12}(s)P)K_2(\theta_2)
 =S_{12}(d)K_1(\theta_1)S_{21}(s)K_2(\theta_2).
\]
Substitution of the lecture's $PS_{12}$ instead exchanges every $12$ with $21$. This distinction is invisible for the permutation-symmetric six-vertex matrix used by its check script.

**Replacements:** at L1:168 and 513 change `$\check R_{\bar\mu,\nu}$` to `$\check R_{\nu,\bar\mu}$`. Replace L1:157's conversion prose by:

> For the second-slot braid convention used below, set $\check S=S_{12}P=PS_{21}$. Using $PK_2=K_1P$ converts the braid equation into the fixed-slot equation term by term. Both reflection matrices act in the second tensor slot, adjacent to the wall.

Delete the claim that multiplying once on the left by $P$ proves the displayed conversion. The book gives the correct convention.

### 2.12. Risk 12 — sine-Gordon boundary submodule

**Verdict: correct eigenvector/eigenvalues at the physical poles; wrong determinant detour and overly restrictive “precisely when”.**

For a common eigenvector $(1,\alpha)^T$ the two equations reduce independently to
\[
 y\alpha^2+\hat\epsilon_1(q-q^{-1})\alpha-y^{-1}=0,\qquad
 y^{-1}\alpha^2-\hat\epsilon_0(q-q^{-1})\alpha-y=0.
\]
At $y_n=i(-1)^ne^{-i\eta}$, $\alpha_n=(-1)^ne^{-\vartheta}$ satisfies both. Substitution in the first rows gives exactly the printed $\beta_1,\beta_0$. `check12_submodule.py` obtains residuals $5.16\,10^{-16}$ and $5.83\,10^{-16}$ at $\lambda=\sqrt2,\eta=2.1,\vartheta=0.43$.

But the determinant for an ordinary eigenvalue is
\[
 \det(\pi_y(b_1)-\beta_1 1)
 =\beta_1^2-\hat\epsilon_1(q+q^{-1})\beta_1+\hat\epsilon_1^2-1,
\]
not the generalized eigenvalue determinant with $\beta_1k_1$. The latter is independent of $y$ and cannot select the bound-state rapidities.

Factoring the two quadratic equations gives four generic common-eigenvector branches:
\[
 y^2\in\{-e^{2i\eta},-e^{-2i\eta},-e^{2\vartheta},-e^{-2\vartheta}\}.
\]
The independent script checks all four, with residuals below $7\,10^{-16}$. Selecting actual boundary states additionally needs scalar poles and the physical strip.

**Replace L4:299–304 by:**

> At the reflected physical-pole rapidity $y_n=i(-1)^ne^{-i\eta}$, a simultaneous eigenvector is $v_n=(1,(-1)^ne^{-\vartheta})^T$. Substituting it in both coideal matrices gives the eigenvalues below. There are also other algebraic reducibility branches; a common invariant line alone does not identify a physical boundary bound state.

The exercise is solvable with this vector and needs no change to its target eigenvalues.

### 2.13. Risk 13 — one-loop exercise

**Verdict: correct derivative and final one-loop coefficient; wrong intermediate equality and proposed self-dual product.**

Direct differentiation yields
\[
 \partial_x\log(x)=\frac{i\pi}{h}\frac{\sinh\theta}{\cosh\theta-\cos(\pi x/h)}.
\]
Thus the exercise's first-order block expansion is correct. For $h=3$,
\[
 [\beta^2]K_1^N
 =-\frac{i\sinh\theta}{8(\cosh\theta+1/2)(\cosh\theta-1)}
 =-\frac14\frac{\cos(\theta/2i)}{\sin(3\theta/2i)}.
\]
`check13_15_particles.py` independently differentiates the amplitude: at $\theta=0.71$ both expressions equal $-0.20806549381901771696i$.

However, $(B/2)$ is not $-(3-B/2)$. The second expression for $K_1^N$ at L6:262 differs from the first by $0.702318830002$ at $B=0.43$. Also the proposed $K^{\rm sd}=K^N K^{\rm cl}/K$ is not self-dual: duality exchanges its ratio $K^N/K$ with its inverse. At the same point its values at $B$ and $2-B$ differ by $1.63360441653$.

**L6:262:** delete the second equality, retaining only $K_1^N=(2+B/2)/((2)(B/2))$.

**Replace part 3 (L6:268) by:**

> Define $\widetilde K(\theta;B)=K_1^N(\theta;B(2-B)/2)$. Show that $\widetilde K(\theta;2-B)=\widetilde K(\theta;B)$ and that $\widetilde K-K_1^N=O(B^2)$. This illustrates why one-loop agreement alone cannot establish self-duality. It does not assert that $\widetilde K$ solves crossing and bootstrap with the original bulk coupling $B$.

The independent expansion validates the coefficient, not its attribution to a particular historical perturbative calculation. I have not read Kim's source for that attribution; it remains **unclear** in this pass.

### 2.14. Risk 14 — $a_2$ fusion angle

**Verdict: wrong introductory angle; correct bootstrap exercise after correction.**

For equal masses $m_1=m_2=m$, a bound state of two particles at rapidity difference $iu$ has $m_c^2=2m^2(1+\cos u)$. Requiring $m_c=m$ gives $\cos u=-1/2$, hence $u=2\pi/3$. The shifts of each constituent relative to the center are $\pm i\pi/3$. The printed $u=\pi/3$ instead gives mass $\sqrt3m$.

With the printed blocks, the actual bootstrap $K_2(\theta)=K_1(\theta+i\pi/3)S_{11}(2\theta)K_1(\theta-i\pi/3)$ holds to $2.34\,10^{-41}$ at $B=\sqrt2/3,\theta=0.71$ in `check13_15_particles.py`. The doubling identity in this exercise is correct: $(x)(2\theta)=-(x/2)(x/2+3)$ at $h=3$.

**L6:237:** replace `$u_{11}^2=\pi/3$ ($\bar u_{11}^2=2\pi/3$)` by `$u_{11}^2=2\pi/3$ ($\bar u_{11}^2=\pi/3$)`. No change to the constituent shifts.

### 2.15. Risk 15 — Neumann limit of classical reflection

**Verdict: wrong block expression at $C=0$; correct rational expression.**

Writing $\gamma=Cm_a^2/(2m)$, the Robin equation gives $K=(ip+\gamma)/(ip-\gamma)$, equal to the printed rational expression. Since $(a)(h-a)=-f_{a/h}(\theta)$, the block form for all three permitted $C$ values is
\[
 K_a^{\rm cl}(\theta)=\big[-(a)(h-a)\big]^{-C},\qquad C\in\{-1,0,1\}.
\]
The lecture puts the minus sign outside the power, giving $-1$ at Neumann. `check13_15_particles.py` finds difference exactly 2 at $C=0$ and less than $3\,10^{-41}$ at $C=\pm1$.

**L3:278:** replace `$-[(a)(h-a)]^{-C}$` by `$[-(a)(h-a)]^{-C}$`, or give $K=1$ separately for $C=0$. The book's `eq-linear-reflection` separates the Neumann case correctly.

### 2.16. Risk 16 — sign of $q$

**Verdict: correct for the normalized six-vertex matrix and matrix crossing test, with scalar-factor scope to be stated.**

Rather than importing an R-matrix, I solved its two central coefficients from the coproduct acting on $|22\rangle$. With $z=y_1/y_2$, diagonal central entry $c$ and off-diagonal entry $b$, the equations are $zc+q^{-1}b=1$, $q^{-1}c+zb=z/q$. They give
\[
 c=\frac{q-q^{-1}}{qz-q^{-1}z^{-1}},\qquad
 b=\frac{z-z^{-1}}{qz-q^{-1}z^{-1}}.
\]
Changing $q$ to $-q$ preserves $c$ and negates $b$, exactly the stated conjugation by $\Pi\otimes1$ of this two-particle braid matrix. This is consistent with a Klein twist; it is not a one-particle gauge equivalence of the entire physical theory.

`check16_qsign.py` computes boundary-state vectors from the independently specified coideal K-matrix and tests proportionality under the braid matrix. At $\lambda=\sqrt2,\theta=0.37-0.12i$:

| bulk sign | ordinary $C$ | antisymmetric $C$ |
|---|---:|---:|
| physical $q$ | $1.86\,10^{-16}$ | 0.772517 |
| $-q$ | 0.772517 | $1.86\,10^{-16}$ |

The lecture script really does check this proportionality, not an unrelated identity. It does **not** check the full scalar crossing equation, nor prove an “if and only if” statement at exceptional boundary parameters. Replace “boundary crossing-unitarity holds if and only if” at L4:263 with “for generic boundary parameters, the matrix part of boundary crossing-unitarity is compatible with ordinary charge conjugation at the physical sign”. The existing mention of the matrix part later in that sentence is appropriate.

### 2.17. Risk 17 — conjectural dictionary presented as a resolution

**Verdict: wrong epistemic status; correct conditional SG expansion and formula-level duality.**

For the proposed SG UV–IR map $f(\mu)=-\cos((\lambda+1)\arccos\mu)/\sin\pi\lambda$, differentiation gives
\[
 f(0)=\epsilon_*,\qquad f'(0)=-\frac{\lambda+1}{2\sin(\pi\lambda/2)}.
\]
At $\lambda=\sqrt2$, independent finite differentiation in `check17_18_fusion.py` gives intercept $-0.825491348399$ and slope $-1.517050514988$, versus the analytic $-1.517050514958$. Substituting $B\mapsto2-B$ in each uniform block proves $K_a^N(B)=K_a(2-B)$ exactly as a formula identity. Neither computation identifies the higher-rank physical boundary to all orders.

The book explicitly labels `cnj-neumann-coideal`, `prp-dm-first-order` (despite its `prp-` id), and `cnj-neumann-amplitude` **[Conjecture]**. In particular `prp-dm-first-order` is conjectural for $a_n$, $n\ge2$, and conditional on the UV–IR relation in rank one. The lecture's assertive callouts erase those qualifications even though some links have “Conjecture” in their captions.

**Replace L5:232–235 with:**

> The proposed reconciliation is the book's conjectural dictionary: $\hat\epsilon=0$ is assigned to the uniform boundary and $s=+1$ to Neumann. The sine-Gordon UV–IR relation gives this assignment in rank one; higher-rank support comes from matrix fusion, classical limits and semiclassical checks. A general field-theoretic derivation remains open.

**Replace L5:239–241's categorical explanation with:**

> A possible explanation is a term non-analytic in the bulk coupling, invisible in the stated perturbative expansion. The book presents this as a supporting argument for a conjecture, not a derivation of the higher-rank Neumann charges.

Use “proposed amplitudes” and “conjectural physical assignment” in the corresponding Index, L5:251 and L6 callouts. Do not weaken the exact identity between the two block formulas.

**Figure correction:** for fixed finite $\lambda$, the zero used in L5's figure is at $\mu_0=\cos[\pi/(2(\lambda+1))]$, not $\mu=1$. At $\mu=1$, $\hat\epsilon=-1/\sin\pi\lambda\ne0$. Change the caption and annotation “Uniform limit ($\mu\to1$)” to “$\hat\epsilon=0$ at $\mu=\mu_0$; $\mu_0\to1$ only as $\lambda\to\infty$”. The plotted zero is correctly located. The historical normalization and attribution to Delius–MacKay were not source-checked here and remain **unclear**.

### 2.18. Risk 18 — affine singlet and scalar fusion

**Verdict: correct in the independently checked $a_2$ case; all-rank analytical telescoping not established by this pass.**

For $V(y_a)\otimes\bar V(y_b)$, write a candidate invariant vector as $v=\sum_jr^j|j,j\rangle$. The finite-node coproduct equations give $r=-qy_b/y_a$. The affine node imposes $r^h=1$. At the lecture's fusion shifts,
\[
 y_a/y_b=\exp(-i\pi\omega+2\pi i/h),\qquad r=e^{-2\pi i/h},
\]
so it really is an **affine** invariant vector at this special spectral ratio, not just a finite-algebra singlet. Reversing the constituent ordering would not give this equation generically.

In `check17_18_fusion.py`, I built the representations from elementary matrices and solved the bulk braid intertwiner in its 15-dimensional weight-preserving subspace, normalized on $|1,1\rangle$. I used no author functions. At $\omega=\sqrt2,\mu=0.29-0.08i$, the affine-singlet residual is $1.39\,10^{-15}$; the non-diagonal fused vector differs from the printed $c_+(\mu)v$ by $1.04\,10^{-14}$ relative norm. This verifies the matrix factor at a different coupling from the lecture script.

I also independently evaluated the **printed** diagonal Gamma product and replaced $F_{11}$ using its crossing equation. At 45-digit precision, `check18_scalar.py` gives relative differences from the three-block formula of $1.44110\,10^{-7}$, $1.80194\,10^{-8}$, $2.25261\,10^{-9}$ and $1.44170\,10^{-10}$ for cutoffs 30, 60, 120 and 300. This is consistent with convergence as the cutoff cubed. Large-cutoff double-precision subtraction is unreliable; the retained calculation uses arbitrary precision.

These computations support the stated diagonal fusion and the $a_2$ matrix projection. They do not provide the advertised “zero skipped algebra” derivation of the infinite-product telescoping, or determine CDD factors uniquely from unitarity and crossing. The non-diagonal scalar $A^{(+)}$ and $F_{11}$ are not fully specified inside this lecture. **L6:80 replacement introductory sentence:**

> Substituting the specified minimal scalar factor gives the following three-block result; the Gamma-product reduction is supplied in the book. Unitarity and crossing alone still allow boundary CDD factors.

The general-rank non-diagonal scalar fusion and the literature attributions to Gandenberger were not independently completed; those parts are **unclear**, not certified by the $a_2$ test.

### 2.19. Risk 19 — classical K-matrices and mixed branches

**Verdict: correct matrices; incomplete branching argument and missing generic-spectral qualifier.**

`check19_classical.py` constructs $E_{ab},F_{ab},H_{ab}=[E,F]$ from scratch and simplifies every entry of $K(-C_jH_j+\lambda E_j-\lambda^{-1}F_j)-(C_jH_j+\lambda E_j-\lambda^{-1}F_j)K$ symbolically. The result is identically zero for arbitrary $C_0,C_1$ in rank one and for all eight rank-two sign patterns. Every rank-two determinant is a nonzero rational function of $\lambda$.

The mixed case can be excluded without assuming that two spectator determinants are nonzero. If $C_1^2\ne1$, node 1 kills $k_{13},k_{23},k_{31},k_{32}$. The active-block equations for all nodes set the diagonal entries to a common $d$. Invertibility then requires $d\ne0$, and node 2's active equation forces $C_2=0$ because $k_{23}=k_{32}=0$. Its determinant is now nonzero, killing $k_{12},k_{21}$ as well. The remaining active equations force all $C_j=0$.

**Replace the first two sentences of L2:286 with this propagation argument.** At L2:238 replace “non-zero and invertible for arbitrary complex numbers” with “generically invertible as a function of $\lambda$, for every fixed pair of complex parameters”. For example $C_0=C_1=0,\lambda=1$ makes the displayed normalization zero; special spectral values must not be included in the assertion.

### 2.20. Risk 20 — continuity and Euclidean coordinate algebra

**Verdict: wrong at both foundational derivations.**

Expanding the first printed continuity equation gives
\[
 (\partial_t-\partial_x)T=(\partial_t+\partial_x)\Theta
 \quad\Longrightarrow\quad
 \partial_t(T-\Theta)=\partial_x(T+\Theta).
\]
This is the opposite of the printed sum/difference placement. Likewise the negative-spin equation gives $\partial_t(T_--\Theta_-)=-\partial_x(T_-+\Theta_-)$. A concrete counterexample is $T=t+x$, $\Theta=t-x$: the initial equation holds ($0=0$) but the claimed transformed equation says $2=0$.

**Replace L1:27–48's densities and fluxes consistently by:**
\[
 \partial_t(T_{s+1}-\Theta_{s-1})=\partial_x(T_{s+1}+\Theta_{s-1}),\quad
 \partial_t(T_{-s-1}-\Theta_{-s+1})=-\partial_x(T_{-s-1}+\Theta_{-s+1}),
\]
\[
 Q_{\pm s}=\int_{-\infty}^0(T_{\pm s\pm1}-\Theta_{\pm s\mp1})\,dx,
\quad
 \frac d{dt}(Q_s+Q_{-s})=
 [T_{s+1}+\Theta_{s-1}-T_{-s-1}-\Theta_{-s+1}]_0=-\dot\Sigma_s.
\]
Here the compact $\pm$ notation means the two explicit densities just displayed. The final $\widehat Q_s=Q_s+Q_{-s}+\Sigma_s$ remains valid. Alternatively one may redefine the symbol $\Theta$ with the opposite sign and change the initial light-cone equations; both conventions work, but mixing them does not.

For the Euclidean map at L1:190, the matrix $\left(\begin{smallmatrix}0&-1\\-1&0\end{smallmatrix}\right)$ has determinant $-1$: it is a reflection, not a $90^\circ$ rotation. A proper rotation taking $x\le0$ to $\tilde t_E\ge0$ is $(\tilde x,\tilde t_E)=(t_E,-x)$.

The stated $\tilde E=-ip,\tilde p=-iE$ implies $\tilde\theta=\theta-i\pi/2$, not $\theta=i\pi/2-\tilde\theta$. Indeed substituting the latter gives $-iE=-m\sinh\tilde\theta$, the wrong sign. The standard pair-amplitude relation can be retained, but it includes crossing/orientation reversal of a leg, not just the asserted coordinate rotation.

**Replace L1:188–202 by:**

> Continue to Euclidean time and rotate the half-plane using $(\tilde x,\tilde t_E)=(t_E,-x)$. The wall becomes an initial surface at $\tilde t_E=0$, and the bulk occupies $\tilde t_E\ge0$. Analytic continuation of a worldline shifts its rapidity by a quarter turn. With one leg crossed and the pair ordered as $(-\tilde\theta,+\tilde\theta)$, the boundary pair-emission convention is $\mathcal K^{ab}(\tilde\theta)=K_{\bar a}^{b}(i\pi/2-\tilde\theta)$.

Remove the inconsistent energy–momentum intermediate equation. Make the same coordinate-map correction in L1's figure caption and the preface's geometry description. This does not invalidate the subsequent scalar crossing equation.

### 2.21. Risk 21 — boundary state and pole interpretation

**Verdict: wrong as an unconditional pair-only statement.**

A single zero-rapidity particle obeys every stated odd-charge constraint: $(Q_s-Q_{-s})|a(0)\rangle=2\gamma_a^{(s)}\sinh0\,|a(0)\rangle=0$. Thus those constraints do not forbid a one-particle component. The lecture later computes a nonzero Lee–Yang coupling $g$, directly contradicting its earlier pair-only formula if that formula is asserted for this boundary.

[Ghoshal–Zamolodchikov](https://arxiv.org/abs/hep-th/9306002), source `sources/9306002/9306002.tex` around (3.42), explicitly describes a “contribution of the zero-momentum particle $A_c$, i.e.” and writes a separate one-particle term in $|B\rangle$. The coefficient convention for this term should not be inferred from the pole without specifying normalization; the pair-only assertion already fails regardless of that convention.

**Add before L1:213 and qualify Index's pair-only sentence:**

> At nonzero rapidity the integrable boundary state creates opposite-rapidity pairs. It may also contain zero-momentum one-particle components when the corresponding boundary couplings are nonzero. The following exponential gives the pair sector, or the full state when these one-particle components are absent.

At L1:279 use “A physical-strip pole may represent a boundary bound state” rather than “a pole ... signals”. At L1:293 replace the absolute fixed/moving-pole diagnostic with:

> A one-point pole has a location fixed by bulk kinematics. Motion with a boundary parameter can help distinguish it from a bound-state pole, but neither this observation nor a simple pole alone excludes Coleman–Thun mechanisms or establishes a new state.

The latter limitation is already demonstrated by the lecture's own independently verified Coleman–Thun example (check 2.6).

### 2.22. Risk 22 — Lee–Yang literature correction and bulk coupling

**Verdict: correct GZ-sign criticism and boundary coupling; wrong bulk three-point coupling.**

I read the actual GZ source. Its (3.51) begins exactly:

```tex
$$R_{(1)}(\theta) = -{{\sinh({\theta\over 2} +{{i\pi}\over
```

The continued source line contains the $\pi/4,\pi/12,-\pi/3$ factors, precisely $-(3)(1)/(4)$ in the lecture notation. This is not a reconstruction from the book. The independently established quadratic bootstrap in check 2.4 selects the opposite overall sign. The GZ correction is therefore supported.

For the explicitly defined bulk amplitude,
\[
 \operatorname{Res}_{\theta=2\pi i/3}S_{\rm LY}(\theta)
 =\frac{2i\sin(2\pi/3)}{\cos(2\pi/3)}=-2i\sqrt3,
 \qquad\Gamma^2=-2\sqrt3.
\]
There is no factor $2+\sqrt3$. Independently, the boundary residue is $-2i(2\sqrt3-3)$, giving $g^2=-4(2\sqrt3-3)$ as printed. `check22_24_residues_energy.py` obtains $\Gamma^2=-3.4641016151377545871$ and $g^2=-1.8564064605510183482$, with imaginary errors below $5\,10^{-20}$ from the limiting procedure.

**L1:452:** replace `$\Gamma^2=-2\sqrt3(2+\sqrt3)$` by `$\Gamma^2=-2\sqrt3$`. The book's Lee–Yang discussion and the lecture script compute the latter residue; the extra factor is a lecture error. The sign of $g$ itself is a field-normalization choice; the residue fixes $g^2$.

### 2.23. Risk 23 — tau positivity, breathers and self-conjugate family

**Verdict: wrong root-location proof; correct open-half-line positivity and self-conjugate family, with endpoint qualifications.**

For $u=e^{m_ax}\in(0,1)$ the exact lower bound is
\[
 \tau_j\ge1-\frac{2u}{1+s}+\frac{1-s}{1+s}u^2
 =\frac{(1-u)^2+s(1-u^2)}{1+s}>0.
\]
It proves the required open-half-line positivity without asserting where every root lies. At $\chi+2\pi ja/h=0$, the roots are $-1$ and $-(1+s)/(1-s)$: both negative, refuting “both roots ... lie at $u\ge1$”. Complex roots are also possible. At the singular endpoint with cosine $-1$, the derivative with respect to $u$ at $u=1$ is $-2s/(1+s)<0$, contrary to the derivative sign in the attempted proof.

**Replace L3:211–222 by the lower-bound calculation above, followed by:**

> At the wall, $\tau_j(0)=4\cos^2(\chi/2+\pi ja/h)/(1+s)$. Thus the solution is regular on the closed half-line only away from the isolated values where one of these cosines vanishes.

The same lower bound proves Exercise 3.2(2) for $y<0$. It does **not** make those periodic solutions regular at $x=0$: for $\rho_-=0$, the phase periodically reaches $\psi=\pi$, making a boundary tau vanish. The book explicitly distinguishes this open-half-line regularity from a fully regular boundary breather. Add that warning to the exercise.

For Exercise 3.2(3), when $0<s<1$ and $\psi=\pi$,
\[
 \tau=(1-e^y)\left(1-\frac{1+s}{1-s}e^y\right).
\]
It has an interior zero at $e^y=(1-s)/(1+s)$ and is negative just before the wall, where it returns to zero. At $s>1,\psi=0$, $\tau(0)=4/(1-s)<0$, so there is again an interior zero. Replace “changes sign between ... and $y=0$” by these precise cases; $s=1$ needs a separate limiting prescription since the displayed normalization diverges.

For the self-conjugate family, $u_j=(-1)^je^{2mx+\xi}$ with $\xi<0$ gives
\[
 \phi'=-\frac{2m}{\beta}\sum_j\alpha_j\frac{u_j}{1+u_j}
 =\frac m\beta\sum_j\alpha_j\frac{1-u_j}{1+u_j},
\]
using $\sum\alpha_j=0$. This is exactly the $C=-1$ boundary condition. The static Hirota solution has a continuous position modulus, so differentiating its energy gives zero by the bulk equation and boundary condition; as $\xi\to-\infty$ it tends to the zero-energy vacuum, establishing $E=0$. Its linearized mode has $p=2mi$, $\Omega=0$, consistent with the pole at $\theta=i\pi/2$. That pole is at the **edge**, not the interior, of the boundary physical strip. This part of Exercise 3.3 is correct with that wording.

### 2.24. Risk 24 — excited spectrum and allowed ranges

**Verdict: correct algebraic energy formulas and degeneracy; incomplete physical-state specification.**

The elementary identity
\[
 2\sin\delta\cos(\pi/6-(2k+1)\delta)
 =\sin(\pi/6-2k\delta)-\sin(\pi/6-2(k+1)\delta)
\]
proves the finite sum for $F(n)$ by telescoping. Its algebraic extension gives the stated $F(n)+F(-n)$, but negative labels are reached along different bootstrap paths, not by binding a negative number of particles.

At $n=m=0$, the second and fourth pairs of the excited amplitude cancel, leaving the printed four-block ground-state factor. At $B=1/2$, $\delta=\pi/12$ and
\[
 F(2)+F(-2)=\frac{\sin^2(\pi/6)}{\sin(\pi/12)}
 =\cos(\pi/12)=F(1).
\]
`check22_24_residues_energy.py` gives both sides $0.96592582628906828675$ and verifies the telescoping sums independently for $n=1,\ldots,7$.

The condition $n+m\ge0$ alone does not define the spectrum at fixed $B$. [Delius–Gandenberger](https://arxiv.org/abs/hep-th/9904002), source `sources/9904002/final2.tex`, equation labelled `boundsa2`, explicitly prints
`-\frac 1{2B} - \frac 12 < n,m < \frac 1{2B}+\frac 12`
for $n+m>0$, and `$0\leq n < \frac 3{2B}+\frac 12$` for the diagonal towers. The independent pole condition for their last binding step, $3/2-(2n-1)B/2>0$, gives the second bound. The corresponding ordinary binding condition gives the positive-label part of the first bound. I read the source for its full negative-label range rather than deriving it from the energy formula alone.

**Add after L6:181:**

> The factors are algebraically defined for integer labels, but the physical spectrum is finite. For $n+m>0$ the allowed states satisfy $|n|,|m|<1/(2B)+1/2$. The diagonal towers satisfy $0\le n<3/(2B)+1/2$. These restrictions come from the physical-strip bootstrap, not just positivity of the energy. Negative labels require the diagonal bootstrap paths.

The figure script does restrict the diagonal curves to their allowed ranges; the plotted $b_{1,0}$ and $b_{1,1}$ are valid over its $B<0.98$ range. The word “rotor” in that figure still needs the qualification of check 2.2.

### 2.25. Risk 25 — what the scripts and figures actually certify

**Verdict: wrong blanket assurance and several concrete caption/index claims.**

I read the six scripts without executing them. Their four/five PASS summaries aggregate particular samples; they do not test every printed blackboard step, exercise, attribution or physical interpretation. Independent calculations above expose errors despite the advertised checks. `check25_coverage.py` additionally parses the report-call inventory, computes the stationary-vacuum energy from the field and potential, and tests the braid conversion on arbitrary matrices (an algebraic identity must hold before imposing special symmetry or the reflection equation).

| Script | What the inspected assertions establish | Important limitation relevant to this report |
|---|---|---|
| L1 | Six-vertex reflection equations at one sample; three scalar Lee–Yang crossing forms; bootstrap/sign; residues; excited-state shift | Uses permutation-symmetric six-vertex R and **first-slot** braid K; its conversion comparison subtracts two already-zero equation residuals. It does not test the typed mixed-representation equation or the false doubling steps. It computes the correct bulk residue numerically, never comparing it to the erroneous $\Gamma^2$ printed in the lecture. |
| L2 | Explicit classical matrices, all eight $a_2$ signs, one off-branch sample, classical reflection equations for two examples, one $\dagger$ boundary | No Dirichlet Cartan derivation; an off-branch sample is not an iff proof. The two-algebra reflection loop accounts for five PASS lines from four report-call sites. |
| L3 | Mirror boundary condition with vacuum shift; delay/position formula consistency; static EOM/BC/energy at three modulus centers; Robin rational expression for $C=\pm1$ | Chooses the **best** sign pattern independently at each modulus sample; never asserts the incorrect interval assignment or legend. Does not check the printed exercise derivative, the false $k$ identities, or the block expression at $C=0$. |
| L4 | Coideal coproduct; SG intertwiner and braid RE; **matrix proportionality** under crossing; common eigenvalue tests at selected $y_n$ | Does not compute full scalar crossing; does not test the generalized determinant detour or the exclusivity of the physical-pole reducibility branch. The header promises a test of general $c_0c_1=q^{-2}$, but the assertions use only $c_j=q^{-1}$. |
| L5 | Sampled existence/nonexistence, K formula, Serre identities on a 27-dimensional module, unitarity and matrix crossing | Mixed $V,\bar V$ reflection equations are tested **up to a fitted scalar**; this is weaker than normalized equality. It does not test the universal inverse-membership argument or the erroneous $Z^T$ statement. |
| L6 | Particle-factor axioms/classical limits for ranks 1–4, one $a_2$ matrix fusion, excited bootstrap/energy steps, finite-difference CT residue | Does not integrate the modulus metric despite its header's inertia promise; does not determine a global rotor Hilbert space; no full Gamma-product telescoping or WKB computation. CT tolerance is $10^{-8}$, not machine precision. The formulas are correct in my independent high-precision residue check. |

The independent arbitrary-matrix braid conversion test has residual $4.00\,10^{-16}$ with $\check S=S_{12}P$ and $10.0376$ with the lecture's convention. This demonstrates why the symmetric special case cannot certify the claimed conversion.

**Script-index errors (existing files left untouched):** `code/lectures/README.md:63–65` contains nonexistent lecture/exercise ids and stale formulas. In particular, row L1 calls the excited factor $K_{-b}$; the checked result is $K_{b+2}=K_{4-b}$. Row L3 interchanges attractive/repulsive $C$ signs and gives the vacuum energy $2m(3-\sqrt3)/\beta^2$, whereas direct quadrature in `check25_coverage.py` gives $-2.3923048454132623$ at $m=\beta=1$, equal to $2(4-3\sqrt3)$; the README value is $+2.5358983848622456$. Correct these rows to the current ids and the actual tested formulas. Row L1 should also describe its tested **identity-boundary** residues rather than claiming the quoted $g^2$ for every $b$.

**Figures:** all six committed SVGs were converted by `rsvg-convert` into PNGs under `referee_scratch/lectures_halfline/figures/` and visually inspected. No source figure-generation script was run, and no SVG was overwritten.

- **Figure 1:** the annotation “$b=0$: identity wall” confuses equal reflection factors with equal boundary conditions. Independently, the ratio of its block factors is exactly 1 at $b=0$; this establishes equality of amplitudes only. DPTW's conclusion explicitly says that reflection-factor data “is not enough to characterise a boundary condition completely.” It distinguishes the identity boundary from every finite-$h$ member of $\Phi(h)$ even though $R_{(1)}=R_{b=0}$. Replace the annotation by “$b=0$: $K_b=K_{(\mathbb1)}$ (equal amplitudes)”. The vertical axis also mixes a normalized label $\operatorname{Im}\theta/(\pi/6)$ with unnormalized tick labels $\pi/6,\pi/3,\pi/2$. Use ticks $1,2,3$, or change the axis title to $\operatorname{Im}\theta$.
- **Figure 2:** the green equal-coupling slice vanishes at $0,\pm1$, but the orange slice with $(C_1,C_2)=(+1,-1)$ vanishes only at $C_0=\pm1$. At $C_0=0$ it is not Neumann. The mixed-branch argument in check 2.19 proves this independently. Clarify this distinction in L2:172's caption.
- **Figure 3:** both caption and plotted legend have the wrong $C_0,C_1$ at $\chi=0.25$. Replace the legend entries by $C_0=+1,C_1=-1,C_2=-1$; the actual tau curves and energy are consistent with those corrected signs (checks 2.1 and 2.25).
- **Figure 4:** L4:263/312 says the figure shows crossing residuals and excitation energies. The actual SVG has a coideal schematic and pole-position curves, with neither a residual plot nor an energy plot. Delete “and shown in @fig-lec4-coideal” at L4:263; remove “together with the residual norm ...” and “and their excitation energies ...” from its caption, or add genuinely computed panels. The pole kinematics themselves agree with $\nu_n=\eta/\lambda-(2n+1)\pi/(2\lambda)$.
- **Figure 5:** the fixed-$\lambda$ zero is correctly plotted at $\mu_0$, but caption/annotation misidentify it with $\mu\to1$; use the replacement in check 2.17. Only as $\lambda\to\infty$ does $\mu_0\to1$.
- **Figure 6:** the plotted formula-level duality and spectral truncation are consistent with checks 2.17 and 2.24. Rename the “Quantized Vacuum Rotor” interpretation as “Rotor-like diagonal tower” or make its conjectural status explicit.

**Replace Index:69 by:**

> The scripts in `code/lectures/` check selected formulas and numerical examples and generate the figures. Their PASS results apply to the assertions actually implemented; they do not verify every displayed derivation, exercise, source attribution or physical interpretation.

Replace L6:227's “to machine precision” by “numerically, with relative tolerance $10^{-8}$ in the lecture script”. A future tightened assertion should compare against a controlled limit, such as check 2.6, before claiming machine precision.

## 3. Entries not checked, and remaining limits

All **25 ranked entries received an independent calculation or explicit derivation**, in the recorded order. No ranked entry was dropped. A mixed verdict above applies to the checked components, not to an entire surrounding chapter.

The following components remain **unclear / unchecked**:

- **Risk 2:** quantum gluing of the three modulus intervals and the self-adjoint boundary conditions for the quantum collective coordinate. The metric and leading energy coefficient were checked; a physical rotor quantization was not established.
- **Risk 6:** the separate internal-propagator and loop counts of the proposed boundary diagram. The printed arithmetic is false; the residue calculation does not by itself fix its topology.
- **Risk 8:** a universal-coideal inverse-membership theorem. The repair proves a statement in represented finite-dimensional algebras; it does not assert the unsupported universal statement. Letzter–Kolb parameter-set and generalized q-Onsager citation scopes were not audited source by source.
- **Risk 13:** attribution of the verified one-loop coefficient to Kim, and historical claims about every earlier self-dual proposal. No source-level verification of those claims was completed.
- **Risk 16:** full scalar SG crossing-unitarity and exceptional-parameter cases. The normalized matrix obstruction and Klein-twisted version were checked independently.
- **Risk 17:** the historical field-normalization argument in Delius–MacKay and an all-rank derivation of the physical boundary/coideal identification. Conditional SG Taylor coefficients, book status labels, and the algebraic duality identity were checked.
- **Risk 18:** an all-rank analytical proof of Gamma-product telescoping, non-diagonal scalar factors/CDD selection, and the complete historical attribution. The diagonal product was checked by controlled numerical convergence for $a_2$, and the $a_2$ matrix fusion was independently constructed.
- **Risks 23–24:** full classification of regular time-periodic solutions beyond the displayed family; completeness of the physical bound-state spectrum at all exceptional couplings. The particular lower bounds, static family, energy identities and quoted generic state ranges were checked.
- **Risk 25:** this was a source audit of the author scripts, not a rerun. No assertion is made that all their PASS tests still execute successfully in this environment. Their scientific coverage was compared with independent calculations. No full HTML/PDF rendering or site-layout review was performed.

The other course exercises were read, but not every subpart received an independent solution. In particular Exercise 2.1's $d_4$ extension, Exercise 2.3's full equivalence-class wording, and the general coefficient constraint in Exercise 4.1 were not completed. No verdict on an unread literature source is implied by a correct lecture formula.

Validation command: `python3 check-refs.py lecture-notes/half-line` passed all seven checks (8 files, 102 ids, 351 bibliography keys). This checks reference structure, not mathematical correctness, and does not detect the stale README ids.

Independent scripts and all downloaded source material are retained under `referee_scratch/lectures_halfline/`. Downloaded TeX sources came directly from arXiv for hep-th/9306002, hep-th/9712197, hep-th/9810098 and hep-th/9904002. The last additional source, hep-th/9810098, was consulted for its discussion of zero-induced pole reduction; its Lee–Yang example does not settle the separate $a_2$ topology count in risk 6. No prior half-line report was consulted.

## 4. Other observations — unverified

These are leads, not additional established findings:

- Index's “solitons are their own antiparticles ($V\simeq\bar V$)” may conflate equivalence of the sine-Gordon multiplet with equality of its soliton and antisoliton states. The opening “all higher-rank affine Toda theories” rigidity language also appears broader than the non-simply-laced classification permits.
- L1's equal coefficients $\gamma_a^{(s)}$ for both chiralities and its claim that linear independence alone proves the multi-species elasticity result may need the extra parity/species hypotheses stated in the book. The lecture's chirality and spin naming merits a separate systematic audit after fixing its continuity equations.
- The “half the Lorentz boosts” wording in L1 and the opening story about translation invariance are physically loose. I have not turned these into referee findings because the consequential mathematical issues above take priority.
- L4's non-local currents are called “spin-1” for generic coupling. This deserves checking against the actual anomalous spin and scaling conventions.
- L5's Schur-lemma proof refers to “matching highest weights” for the proportionality constant without spelling out normalization across mixed modules. The author script checks mixed reflection equations only up to scalar. A normalized all-module proof would be useful.
- Several theorem/chapter numbers are typed into link text despite repository rules. The reference checker cannot make those link captions track renumbering.

Phase-1 work stops here. This report is to remain unchanged after delivery; any directed follow-up belongs in the phase-2 report.

---

*Integrated 2026-10-25: see the `## Integration` section of `referee_lectures_halfline_combined.md`.*
