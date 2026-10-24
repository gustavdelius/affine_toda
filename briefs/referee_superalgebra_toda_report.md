# Referee report: affine Toda theories based on Lie superalgebras

Date: 2026-10-23. Referee: Codex (GPT-6).

Reviewed commit: `9820f47dfd28f50fbf28b4707bc76bc6f78881ad`, “Add chapter 16, affine Toda theories based on Lie superalgebras”, branch `claude/superalgebra-toda`, worktree `.claude/worktrees/superalgebra-toda`.

Assignment: `briefs/referee_superalgebra_toda.md`. All line numbers below refer to this commit. No book files were edited, no subagents were used, and no rendering or index-changing git commands were run.

## Recommendation and evidence

Revise before accepting. The component action, Lax curvature, two super sinh-Gordon currents, DGPZ free-fermion correction, negative sign of gamma, and the tested soliton masses survive independent checks. The three proposed Appendix D corrections have substantive support. However, the surrounding exposition contains errors that the author's successful checks do not catch:

1. The heaviest boson of the C family is **not** obtained as a fermion–fermion bound-state pole. The source explicitly excludes it from the bootstrap derivation.
2. The exceptional-family cubic Yukawa term has the wrong sign in the stated source normalization.
3. Gandenberger's sign is `(-1)^(ab)` after fusion, not a universal minus sign for all scalar factors.
4. The displayed central-charge parametrization is undefined for the central example m=2, a=1; its finite limit should be written explicitly.
5. Changing from projected psl root data to a gl theory is not generally the addition of decoupled free superfields.
6. The pure-massiveness criterion needs the hypothesis “affine”. The boundary discussion overstates both necessity of a boundary fermion and what the cited conservation tests establish.

The complete physical charge-conjugation derivation of crossing remains unclear. The all-loop current calculations are verified here as accurately scoped source claims, not independently reproduced all-loop calculations. Some historical citations and bibliographic fields remain unchecked, listed at the end.

Independent evidence is under `referee_scratch/superalgebra_toda/`:

- `independent_checks.py` and `.log`: exterior-algebra expansion, superjet current reduction, matrix relations, root Hessians, tree scattering, braid YBE, C-family amplitudes and exceptional-family checks.
- `soliton_scalar_checks.py` and `.log`: projected roots, potentials, tau-function equations, first-order equations, integrated energies, central terms, Gamma products, residues and bootstrap.
- `scope_checks.py` and `.log`: ODE elimination, boundary supercharge equations, reality signatures, threshold kinematics and N=1 scalar factors.

These scripts import no book check code. They use independently implemented Grassmann multiplication and differentiation, direct differentiation of root potentials, and the formulas displayed in the chapter or primary papers. The soliton logarithms are differences of logarithms of tau functions, so that the field is traceless; taking principal logarithms of their ratios independently can introduce a spurious trace and spoil the supersymmetry check.

Run from the reviewed worktree root:

```sh
python3 -u referee_scratch/superalgebra_toda/independent_checks.py
python3 -u referee_scratch/superalgebra_toda/soliton_scalar_checks.py
python3 -u referee_scratch/superalgebra_toda/scope_checks.py
```

Final results: 34, 44 and 8 PASS lines respectively, no FAIL, all exit status 0. These 86 checks include negative controls (detecting the wrong printed identity), not 86 independent mathematical claims. The numerical Gamma-product tolerance is 3e-6; observed bootstrap errors are below 1.4e-9. The algebraic checks are symbolic or use the tighter tolerances recorded in the scripts.

As a secondary check, I ran all six author scripts with bytecode writing disabled: 9, 31, 11, 27, 7 and 2 PASS lines, respectively, no FAIL, and exit status 0. Approximate runtimes were 1.1, 13.8, 4.6, 25.3, 40.5 and 3.2 seconds. Their full output is in `author_*.log`. `python3 check-refs.py` also passes, with only the existing unnumbered-page warning for `citation.qmd`. These runs are corroboration, not the independent evidence.

Primary texts were found in the drafter's temporary source directory, including both DGPZ papers and Evans–Hollowood 1991. Copies of the extracted texts are retained in the referee scratch `sources/` directory. Short quotations below preserve source notation; mathematical derivations that follow them are mine. ArXiv identifiers and equation numbers identify the sources independently of the temporary directory.

## 1. Component Lagrangian and signs

**Verdict: correct.**

Let theta-plus precede theta-minus, and place both before the component fermions. All odd generators anticommute; differentiation with respect to theta is from the left. For one component, direct differentiation gives

\[
D_+\Phi=i\psi_+-i\theta^+\partial_+\phi+i\theta^-\sigma
 +\theta^+\theta^-\partial_+\psi_-,
\qquad
D_-\Phi=i\psi_--i\theta^+\sigma-i\theta^-\partial_-\phi
 -\theta^+\theta^-\partial_-\psi_+.
\]

Taking i times the theta-plus theta-minus coefficient of the kinetic term gives

\[
\tfrac12\partial_+\phi\partial_-\phi+\tfrac12\sigma^2
+\tfrac i2(\psi_+\partial_-\psi_++\psi_-\partial_+\psi_-).
\]

For an odd-root exponential its integrated contribution is

\[
-\frac{\mu}{\beta}(\alpha\cdot\sigma)e^{\beta\alpha\cdot\phi}
+i\mu(\alpha\cdot\psi_+)(\alpha\cdot\psi_-)e^{\beta\alpha\cdot\phi}.
\]

The explicit even-root spurion integrates to `-mu exp(beta alpha.phi)/beta^2`, with no auxiliary or Yukawa term. Completing the auxiliary square gives exactly `eq-super-components`. The generic calculation extends componentwise with the metric G; no positivity assumption is used.

Writing W for the odd-root superpotential, the component equations are

\[
\sigma=\nabla W,\quad
\partial_-\psi_+=-W''\psi_-,\quad
\partial_+\psi_-=W''\psi_+,\quad
\Box\phi=-W''\nabla W+iW'''(\psi_+,\psi_-)-\nabla V_{\rm even}.
\]

The theta-zero, theta-minus, theta-plus and theta-plus theta-minus coefficients of `eq-super-eom` reproduce them. In particular the exercise's real fermion mass equations have the correct signs; there is no additional i in the mass parameter M of that exercise.

For SSG, `W=4m cosh(beta phi/sqrt(2))/beta^2`. Therefore

\[
V=4m^2\sinh^2(\beta\phi/\sqrt2)/\beta^2,
\quad V''(0)=4m^2,
\quad W''(0)=2m.
\]

The bosonic truncation is the book's a1 affine potential including the vacuum subtraction. Continuing beta to i beta gives `4m^2 sin^2(beta phi/sqrt(2))/beta^2` and the stated cosine Yukawa term. No extra coupling rescaling is needed.

For the mass exercise, put `w=dW` in coordinates. At a **supersymmetric** vacuum w=0,
`Hess(V)=Hess(W) G^-1 Hess(W)`, hence the linearized bosonic operator is `(G^-1 Hess(W))^2`. Stationarity of V alone would not suffice, but the exercise explicitly imposes `sum mu_j alpha_j=0` and is sound. The source formula is Evans–Madsen 1997 (3.4); the derivation above checks it rather than assuming it.

Output: kinetic, odd exponential, even spurion, SSG mass and truncation checks all PASS in `independent_checks.log`.

## 2. Graded Lax pair and algebra identification

**Verdict: correct for the stated curvature and root relations.** The four-relation matrix test is not by itself a presentation theorem for the entire affine superalgebra.

Use the Grassmann envelope convention
`(c X)(d Y)=(-1)^(|X||d|) cd XY`, with coefficients on the left. In the anticommutator of the two odd connections, the Cartan/root term is

\[
\{\beta D_+\Phi\cdot H,\eta_j e^{\beta\alpha_j\cdot\Phi}f_j\}
=-\beta(\alpha_j\cdot D_+\Phi)\eta_j e^{\beta\alpha_j\cdot\Phi}f_j.
\]

This cancels the differentiated exponential at order lambda^-1. There is no order-lambda term. At order zero the curvature is

\[
-\beta D_+D_-\Phi\cdot H+\sum_j\eta_j e^{\beta\alpha_j\cdot\Phi}\alpha_j\cdot H.
\]

Taking `eta_j=-i mu_j` gives precisely the multiple of the field equation printed at line 138. The one-field realization gives the extra sqrt(2) because `H=sqrt(2)h`.

My direct 3-by-3 products verify all four anticommutators, both raising-generator weights and zero supertrace. The corresponding lowering weights also follow directly from the displayed matrices. For the suggested osp ansatz, `(E31-E23)^2=-E21`, so using the same odd lowering generator in e0 and f1 leaves a nonzero cross anticommutator.

The source identifications agree: Evans–Madsen 1997 Table 4 has `A(1,0)^(2)=C(2)^(2)`, and Ito–Zhu Sec. 3 treats `C(2)^(2)=osp(2|2)^(2)`. The latter's purely odd rank-one system is the relevant one. B(0,1) untwisted has an even affine simple root and does not give this manifestly supersymmetric action. The failed osp ansatz is a useful obstruction, but the source classification, not that single failed ansatz, establishes the algebra identification.

Output: both matrix checks PASS. The order-by-order calculation above is the independent general Lax check.

## 3. Supercurrents and Evans–Madsen correction

**Verdict: correct. Keep the Appendix D entry.**

In the source normalization put `u_k=D_+^k Phi`. Then `i D_+D_-Phi=-M sinh(Phi)` and

\[
D_-u_1=-iM\sinh\Phi,\qquad D_-u_{k+1}=-D_+(D_-u_k).
\]

The independent superjet algebra uses commuting even u's and exterior odd u's. For `T=u1 u2-u3`, direct reduction gives

\[
D_-T=iM e^{-\Phi}u_2
=D_+\big(iM e^{-\Phi}u_1\big)
=D_+\big(-iM D_+e^{-\Phi}\big).
\]

Applying the same odd Leibniz rule to `T D_+T` gives

\[
D_-(T D_+T)=D_+\left[iM(u_1u_2^2-u_2u_3)e^{-\Phi}\right].
\]

Both are exact PASS checks, not tests at particular component values. The source (5.44) prints `D_+J=D_-bar J`; (5.45) prints `bar J=(m/2)iD_+(exp(-Phi))`. These disagree with the result from its own (1.5): the first partner needs multiplication by -2, the second by +2.

There is a separate convention-independent check of the derivative ordering. Substitute the source's (5.47), `J=alpha+theta+ j+theta- h+i theta+ theta- beta`, and its barred counterpart into `D_-J=D_+bar J`. The four coefficients give

\[
h=\bar h,\quad \partial_-\alpha=\bar\beta,\quad
\partial_+\bar\alpha=\beta,\quad \partial_-j+\partial_+\bar j=0.
\]

These are exactly the source's (5.48)–(5.51). The printed ordering instead relates the other derivatives and components. It cannot be rescued by simply exchanging the derivative labels while keeping the source's field expansion and action.

At zero fermions, the theta-plus coefficient of `T D_+T` is `(p^2-p')^2`, where `p=partial_+ hat phi`. Its difference from `p'^2+p^4` is `-(2/3) partial_+(p^3)`. With `hat phi=beta phi/sqrt(2)`, this is the density of `eq-spin3-current` up to an overall normalization and a total derivative, as stated.

The source uses m, whereas the chapter uses M=-2m when comparing actions. The Appendix D entry correctly keeps the source's m separate from the book's physical mass parameter.

## 4. Boson–fermion S-matrix

**Verdict: unclear only for the claimed complete physical crossing justification; correct for the independently checked algebraic formulas, tree sign, poles and bootstrap.**

DGPZ (2.4) gives the conversion amplitudes with opposite signs. In the ordered basis `(phi phi, phi psi, psi phi, psi psi)` the chapter's matrix is

\[
\check S=\begin{pmatrix}B&0&0&-D\\0&R&T&0\\0&T&R&0\\D&0&0&F\end{pmatrix}.
\]

Its free limit is the graded permutation. All entries preserve total parity, so the embeddings in three tensor factors are ordinary Kronecker products. Expanding the braid equation gives the component relations of DGPZ (2.7), with the conversion sign retained; the matrix residuals in my test at n=2,3 are below `2e-16`. Changing only F to -F changes this equation and is not an equivalent convention change. A consistent change of convention would have to transform the entire operator and tensor prescription.

Write `r=R/T`, `d=D/T`. Direct multiplication gives

\[
(\check S/T)(\theta)(\check S/T)(-\theta)
=\left(1+\gamma^2/\sinh^2(h\theta/2)\right)1.
\]

Thus the scalar unitarity equation in the chapter is exactly the needed one. My separate summation of the Gamma product, with 2000/4000-term tail extrapolation, gives scalar unitarity errors `2.5e-9` and `4.4e-9`. Physical unitarity also passes. This is numerical corroboration; it is not an exact convergence proof. Analytically, in the open physical strip every Gamma argument in the numerator has positive real part: the smallest shifts are `eta+Omega` and `eta+1-Omega`. The same holds for the crossed product. The paired log-Gamma terms have an O(l^-2) tail, uniformly on compact subsets. Consequently there are no scalar poles in that strip; a grid scan is unnecessary for this claim.

**Independent tree sign.** Set beta_DGPZ=1 when differentiating the potential, then restore its power. For A, the even-root exponential coefficients are `1,1,2,...,2` on `-e1-e2,e1-e2,e2-e3,...`, and the odd-root square gives coefficient 1 on `2en`. The Hessian's eigenvector of eigenvalue 2 is the external particle. Let `g_k=V'''(v,v,v_k)` and `c4=V''''(v,v,v,v)`. I obtain `c4=2` for n=2,3 and evaluate

\[
\mathcal M/\beta_D^2=-c_4-\sum_k g_k^2
\left[(s-m_k^2)^{-1}-m_k^{-2}+(u-m_k^2)^{-1}\right].
\]

This includes the contact term and all s,t,u exchanges, with t=0. The S coefficient is `i M/(8 sinh theta)`. Independently differentiating the displayed product at Omega=0 gives

\[
\left.\partial_\Omega\log\Pi(\theta)\right|_0
=h^{-1}\left[\psi(\eta/h)-\psi(\eta/h+1/2)\right].
\]

Adding the crossed term and the matrix contribution `-2i pi/sinh(h theta)`, then dividing by 4 pi, reproduces the tree coefficient to less than `4e-17` at theta=.71 and 1.19, for both ranks. The positive-gamma branch reverses the matrix contribution and fails. This fixes the sign locally on the perturbative branch. Extending it as `-sin(pi Omega)` throughout 0<Omega<1 assumes the continuous branch specified in the conjecture.

At a direct-channel pole `theta=i pi(h-2a)/h`, the scalar is finite and positive: all Gamma arguments are positive real numbers there. The residue of B and F is `-2i gamma T/h`, hence i times a positive number. The conversion residue has equal magnitude and opposite off-diagonal signs; the 2-by-2 determinant vanishes. My numerical eigenvalues of `Res(S)/i` are `(0,.6578743)` for n=2 and `(0,.4784043)`, `(0,.4161382)` for n=3. The scripts also check 1/epsilon scaling. Kinematics gives `s=4m_n^2 cos^2(u/2)=8 sin^2(pi a/h)`. Substitution in DGPZ's residue bootstrap (2.44) reproduces the displayed block product for S_an in the tested ranks. This checks a bootstrap identity, not the identification of the solution with the full quantum field theory.

**Crossing limitation and replacement.** The hyperbolic formulas give `D(i pi-theta)=i R(theta)` directly. In a chiral gamma-matrix basis, `u=(e^(theta/2),e^(-theta/2))` and `v=(e^(theta/2),-e^(-theta/2))` satisfy `u(theta+i pi)=i v(theta)`; one may use `C=[[0,1],[-1,0]]` with gamma-zero off diagonal. This is a check on spinor conventions. It does not by itself derive crossing of ordered many-particle amplitudes, including LSZ phases and fermionic reordering. The book gives no such derivation. Replace the second sentence of chapter 16 line 252 by:

> This phase agrees with the analytic continuation of the Majorana spinors used in Delius et al., App. A: $u(\theta+i\pi)=i\,v(\theta)$. We use their amplitude crossing convention; a derivation from charge conjugation and ordered asymptotic states is not given here.

The conjecture label is appropriate. Its evidence should not be strengthened to an all-orders identification, and the C-family exception in item 6 must be added to its scope. Lines 15 and 281, and chapter 14 line 96, should say that mass ratios are preserved **at one loop**, while a fixed Coxeter number belongs to the proposed exact amplitudes.

## 5. DGPZ free-fermion condition

**Verdict: correct correction. Keep Appendix D line 19.**

DGPZ Appendix B explicitly prints `BF+T^2-D^2-R^2=0` in (B.1), and immediately continues with `d(theta)=D(theta)/(iR(theta))` in (B.2). Thus Appendix B has not silently replaced the amplitude D by iD.

Put `z=h theta/2`. The displayed amplitudes give

\[
BF/T^2=-1-\gamma^2/(\sinh^2z\cosh^2z),\quad
D^2/T^2=-\gamma^2/\cosh^2z,\quad R^2/T^2=-\gamma^2/\sinh^2z.
\]

It follows identically that `BF+T^2+D^2-R^2=0`; the printed combination instead equals `2 gamma^2 T^2/cosh^2 z`. My n=2,3 checks verify both the correct zero and the nonzero counterexample. The eight-vertex interpretation agrees, since the two conversion weights multiply to `D(-D)`. This is an internal inconsistency of the source formulas, not a claim that every convention must write the free-fermion relation with this sign.

## 6. C family, root reconstruction and masses

**Verdict: wrong in the bound-state scope and one root sentence; masses and fermion block identities are correct.**

Starting from the chapter's roots and eliminating the odd-root auxiliary contribution gives the following bosonic exponent sets:

- A: `-e1-e2`, `e1-e2`, `e2-e3`, ..., `e(n-1)-en`, `2en`, with weights `1,1,2,...,2,1`.
- C: `-2e1`, `e1-e2`, ..., `e(n-1)-en`, `2en`, with weights `1,2,...,2,1`.

Their weighted root sums vanish. The Hessians are sums of weighted outer products. Independently obtained eigenvalues are A n=2 `(2,6)`, A n=3 `(2,5-sqrt(5),5+sqrt(5))`, C n=2 `(4,8)`, C n=3 `(2,6,8)`. The odd superpotential coefficients sqrt(2) give fermion eigenvalues sqrt(2) on en in A, and sqrt(2) on e1,en in C. This verifies `eq-dgpz-masses` and the source diagram reconstruction for the requested ranks.

Chapter 16 line 195 should begin:

> The bosonic potential contains the even root $2\epsilon_n$ in both families, and also $-2\epsilon_1$ in the $C^{(2)}(n+1)$ family. Its other simple roots are the even simple roots listed above.

The old parenthesis can be read as assigning -2e1 to both potentials; it is absent from the A potential.

The block definitions imply unitarity of every factor. Under crossing the C blocks pair with their h-x partners. Direct tests at n=2,3,4,5 confirm that crossing interchanges the two even-n fermion amplitudes and preserves each odd-n amplitude. The order of the two odd-n products agrees with DGPZ (3.13).

**Heaviest-particle error.** After (3.17), DGPZ explicitly says the nth-boson pole is “actually absent in the fermion-fermion amplitude”. The letter states the same exception after its (26). Independently, `m_n^2=8=4m_psi^2` places this particle at threshold u=0, not at a physical-strip bound-state pole. Replace the end of chapter 16 line 269, from “The bosons are bound states”, by:

> The bosons $\phi_a$, $a=1,\ldots,n-1$, occur as fermion–fermion bound states. Their bootstrap gives $S_{\chi a}=\{n-a+1\}\{n-a+3\}\cdots\{n+a-1\}$ and the boson amplitudes below, initially for $a,b<n$. The heaviest boson has $m_n=2m_\psi$, and its putative fermion–fermion pole is absent. Delius et al. propose extending these formulas to $a=n$ or $b=n$ because the extension agrees with tree amplitudes; those cases are not derived by this bootstrap [@dgpz1991b, after eq. (3.17)].

Add to `cnj-dgpz-smatrix`, line 278:

> In the $C^{(2)}(n+1)$ family, amplitudes involving the heaviest boson are proposed extensions of the bootstrap formulas. Their one-loop threshold analysis requires separate care because $m_n=2m_\psi$.

The description of assuming no fermion-pair annihilation is faithful, but the nontrivial kinematic case is n divisible by 3: then `m_(n/3)=sqrt(2)`. Otherwise a boson has no degeneracy allowing that process in an elastic two-particle S-matrix. The source makes this distinction. The folding history and need for a fermion box contribution are accurately attributed to the letter, Sec. 4. I checked the attribution and threshold logic, not the full one-loop box coefficient independently.

## 7. N=2 family and solitons

**Verdict: wrong in the unrestricted gl-decoupling assertion and singular parametrization; the main root, potential, soliton and mass formulas are correct.**

Taking i=j in `(lambda_i+lambda_j)a_ij=0` gives `a_jj=0`. A nonzero edge forces opposite eigenvalues. A four-cycle is bipartite; a complete four-vertex graph is not. D(2,1;alpha) is exceptional here only for admissible alpha, excluding 0 and -1 where the algebra degenerates. This proves the stated N=2 obstruction without using a check script.

With `P=1-11^T/m`, realize the roots as
`alpha_(2j-2)=(-Pe_j,iPe_(j-1))`,
`alpha_(2j-1)=(Pe_j,-iPe_j)`.
Their ordinary bilinear Gram matrix has rank 2m-2 and the two indicated null vectors. For `z=(p,q)`, summing the squared superpotential gradient cancels all non-neighbouring terms and gives

\[
V=\frac{\mu^2}{\beta^2}\sum_j
\left[e^{\beta(p_j-p_{j+1})}-e^{i\beta(q_j-q_{j+1})}\right].
\]

The vacuum constants cancel between the two Toda potentials. My script independently checks this decomposition and the complex structure for m=2,3,4.

For `tau_j=1+exp(rho(x-vt)+2 pi i a j/m+xi)`, use
`q_j=(i/beta)(log tau_j-log tau_(j-1))`. Differentiating the rational tau expressions verifies the field equation, with `rho^2(1-v^2)=4mu^2 sin^2(pi a/m)`. The residual is below `1.5e-14` in all tested sectors of m=2,3,4. Direct rest-energy integration gives `8`, `6sqrt(3)`, and `(8sqrt(2),16,8sqrt(2))` for mu=beta=1, with errors below `1e-12`. The asymptotic root exponentials are w on odd roots and w^-1 on even roots, where `w=exp(2 pi i a/m)`; substituting gives the same central term.

**First-order sign.** From the stated supercharges and the scalar superfield expansion one obtains

\[
\delta\psi_+=-\epsilon^+\partial_+\phi+\epsilon^{+'}J\partial_+\phi
+\epsilon^-F+\epsilon^{-'}JF,
\]
\[
\delta\psi_-=-\epsilon^-\partial_-\phi+\epsilon^{-'}J\partial_-\phi
-\epsilon^+F-\epsilon^{+'}JF.
\]

Substitute `epsilon-=lambda epsilon+`, `epsilon+'=sigma epsilon+`, `epsilon-'=lambda sigma epsilon+`. This derives the first-order equation, not merely a test of it. For positive rho and the chapter's `J alpha_j=i(-1)^j alpha_j`, the solution has negative lambda. My residuals are below `4e-15`. Crucially, Evans–Hollowood (5.3) specifies **rho squared**, not its sign. Choosing negative rho changes the required lambda sign. Their positive sign need not be a physical error. The Notes should explain the orientation convention rather than leave an unspecified suspicion about the supersymmetry variations. Their printed second line of (5.5) also lacks J in its primed derivative term; the transformation law above supplies it directly.

Suggested replacement at chapter 16 line 369:

> For the soliton orientation $\rho>0$, with $J\alpha_j=i(-1)^j\alpha_j$, the first-order equation uses $\lambda=-\sqrt{(1-v)/(1+v)}$. Reversing the sign of $\rho$ reverses this sign. Evans and Hollowood specify $\rho^2$, so this is an orientation convention, not a correction to their mass formula.

**Removable singularity.** In `eq-n2-central`, `sigma=tan(pi a/(2m)+pi/4)` is undefined at a=m/2, including m=2,a=1. The mass remains finite, but the displayed definition needs an explicit limiting prescription. Replace line 318 by the equivalent regular formula

\[
M=\frac\mu{\beta^2}\left[c_a\sum_{j\ {\rm even}}e^{\beta\alpha_j\cdot\phi}
+c_a^{-1}\sum_{j\ {\rm odd}}e^{\beta\alpha_j\cdot\phi}\right]_{-\infty}^{+\infty},
\qquad c_a=i e^{i\pi a/m}.
\]

Explain that this equals the sigma expression wherever sigma is finite and extends it continuously. For the first-order equation multiply through by `cos(pi a/(2m)+pi/4)` before taking that limit. My checks include the otherwise omitted m=2 and m=4,a=2 cases.

**Notation error.** The sl/psl quotient statement is right. The last sentence of line 299 is too strong. In gl(m|m), I is null under supertrace and pairs with `K=diag(1_m,-1_m)`: `str(IK)=2m`. In a field `Y+uI+vK`, the kinetic term couples u and v, while the odd-root exponentials contain `exp(plus or minus 2 beta v)`. Variation with respect to u makes v free, but nonconstant v still enters the interacting field equations. Free equation does not mean a decoupled direct summand. Evans–Hollowood's actual reduction says “We may therefore consistently set Φ′ = 0”; that is a truncation. Replace the last sentence of line 299 by:

> Here we work with the projected root space after fixing the null superfield combinations to zero. This gives the interacting system on the $2m-2$ dimensional space with non-degenerate form. An extension to a larger Cartan space requires specifying its additional fields and their couplings; it is not in general the addition of decoupled free superfields.

The classical-only qualification at line 321 is appropriate. A positive-norm BPS inequality for every complex or indefinite theory is not established by this boundary identity, and the text should continue to avoid claiming one. The N=2 S-matrix remains a conjecture. The optional m>2 extrapolation should retain Evans–Hollowood's caveat “without the factor Smin” after (5.10), because the minimal supersymmetric factor supplies that pole structure. I have not independently rederived the quantum coupling dictionary of Kobayashi–Uematsu.

## 8. Exceptional family

**Verdict: wrong sign for the cubic Yukawa term; the root data and mass/sign calculations are correct in their domains.**

Using the diagonal metric `G=diag(-(1+alpha)/2,1/2,alpha/2)` and the four sign-vector roots in the chapter gives zero diagonal Gram entries and the three products alpha,1,-1-alpha. The fermion operator, divided by mu, is
`diag(-2(1+alpha),2,2alpha)` in this root-coordinate basis. Squaring gives the stated boson masses. For unrestricted real alpha the positive mass magnitudes are `2|mu| |1+alpha|`, `2|mu|`, `2|mu| |alpha|`.

The real-field kinetic coefficients are proportional to `(-1/(1+alpha),1,1/alpha)`. Continuing fields 2,3 for alpha>0 makes all three negative; reversing the overall action gives positive kinetic energy. Continuing 1,3 for -1<alpha<0 makes them positive. Continuing 1,2 for alpha<-1 makes them negative. These are direct sign calculations, independently implemented in `scope_checks.py`.

The source's cubic formula in Sec. 3.3 is

\[
L_{\rm yuk}=-4\mu\sqrt{\alpha(1+\alpha)}
(\bar\psi_1\psi_2\phi_3+\bar\psi_1\psi_3\phi_2+\bar\psi_2\psi_3\phi_1).
\]

This is **negative**, not positive as in chapter 16 line 327. It can also be recovered without importing that final formula: expand source (3.32) to cubic order, reverse the overall negative kinetic sign, and rescale the three fields by `sqrt(1+alpha),1,sqrt(alpha)`. The three mixed coefficients become `-2m sqrt(alpha(1+alpha))=-4mu sqrt(alpha(1+alpha))`. A different fermion phase convention would require a stated dictionary; the chapter explicitly claims the source normalization.

Replace the Yukawa phrase in line 327 by:

> In the canonically normalized, positive-kinetic convention of Evans and Madsen, with their coupling set to one, the cubic Yukawa interaction is $\mathcal L_Y^{(3)}=-4\mu\sqrt{\alpha(1+\alpha)}(\bar\psi_1\psi_2\phi_3+\bar\psi_1\psi_3\phi_2+\bar\psi_2\psi_3\phi_1)$. The full Yukawa interaction is non-polynomial, as in their (3.32).

Also state `alpha !=0,-1` at the beginning, and either use absolute mass magnitudes or explicitly restrict to alpha>0 before listing them. Restoring beta in canonical physical fields multiplies this cubic coefficient by beta.

The folding/limit summary is consistent with the source's Sec. 4, but it is an existence statement at the appropriate parameters: at alpha=1 the interchange of sectors 2 and 3 gives the B(1,1) folding; further truncation gives the rank-one model. In the alpha->0 limit a rescaling is needed to separate the degenerating/free sector. This is not literal substitution into `1/alpha`, and not a folding at every alpha. Add these qualifications if retaining the sentence. I did not independently prove every Lax-algebra folding map.

## 9. Quantum-integrability claims

**Verdict: correct as source reports, with a needed one-loop qualifier on the mass statements outside this subsection.**

Penati–Zanon's conclusion says “renormalized expression to all-loop order”. Their equations (10)–(12) specify normal ordering and massless superspace perturbation theory, extracting local anomalies from one interaction insertion. The fact that this is first order in the interaction does not mean first order in beta: repeated contractions with that insertion generate the higher powers. Independently counting fields in their spin-3 ansatz explains why a finite polynomial in beta can express an all-orders result: only finitely many current fields can contract with the single normal-ordered exponential. This is a structural consistency check, not a recalculation of the coefficient system.

Their (17) gives mass squares 2 and 8. The one-loop self-energy expressions following (29) cancel on shell, and the conclusion explicitly limits the mass statement to “lowest order in perturbation theory”. Their charge formulas (31)–(32) preserve omega2/omega1=1/4. The elementary consistency ratio `m1^2/m2^2=1/4` and threshold `m2=2m1` match the conservation-law discussion. The vanishing three-point function and lack of one-loop mass shifts are accurately reported at chapter 16 line 285. I have not independently evaluated those loop integrals.

Evans–Madsen Sec. 6 similarly normal-orders the exponential interaction and uses massless propagators. Their counting rules eliminate contractions of two odd-derivative fields; the highest nontrivial contractions in their spin-7/2 ansatz are triple contractions. Solving for two independent current coefficients gives the two currents, with a finite beta expansion; it is not merely a two-loop truncation with unexamined higher contractions. They state that the calculation does not explicitly use the reality conditions. Thus the chapter correctly describes perturbative quantum integrability, including the twisted versions. Exclude alpha=0,-1 from “every alpha”, as the algebra and propagator setup degenerate there.

For GPZ's one-loop mass results I verified the attribution against DGPZ's introduction and the letter's bibliography, not against a newly obtained full GPZ paper. The letter points to PLB 253 (1991) 357. DGPZ NPB reference [16], PLB 234 (1990) 88, is not the mass paper. The book's `gpz1991` is the correct reference. This bibliographic mismatch should not be used as evidence for a physical correction.

Replacement for the relevant part of chapter 14 line 96:

> Their classical mass ratios survive the one-loop calculation. The proposed exact S-matrices keep the Coxeter number fixed at $2n-1$ and $2n$, respectively.

Use the same distinction in chapter 16 lines 15 and 281. The conjecture's status must govern these unqualified summaries as well.

## 10. Gandenberger scalar factors

**Verdict: correct Appendix D correction for (1,1), but wrong universal-sign wording in chapter 23.**

The source states “they are identical” before (5.10). Substituting lambda=1/h and a=b=1 in (3.28), then dividing (5.6) by it, uses only `Gamma(z+1)=z Gamma(z)`. The ratio through M factors is

\[
\prod_{j=1}^M\frac{(j-u)(u+j-1)}{(j-1-u)(u+j)}
=-\frac{M-u}{M+u}\longrightarrow-1.
\]

The telescoping calculation is exact; my script checks it symbolically at M=1,2,7. This confirms the (1,1) correction without assuming a regulator-dependent phase.

However, the source's (3.12) and (5.7) fuse a*b copies of the respective (1,1) factor with identical shifts at lambda=1/h. Therefore

\[
G^{(N=2)}_{a,b}=(-1)^{ab}F_{a,b}|_{\lambda=1/h}.
\]

Already a=2,b=1 is a counterexample to a universal -1. Replace the second sentence of chapter 23 line 125 by:

> At this coupling the fundamental scalar factors obey $G^{(N=2)}_{1,1}=-F_{1,1}|_{\lambda=1/h}$. Their fusion formulas give $G^{(N=2)}_{a,b}=(-1)^{ab}F_{a,b}|_{\lambda=1/h}$ in the same normalization.

Appendix D line 63 is correct as it explicitly concludes only about (1,1); it may add the fused formula. The author's script checks only (1,1), which explains why its PASS did not catch the broader prose claim.

For N=1, substituting h=2n+1 and lambda=1/h in (4.4) matches all eight numerator Gamma arguments of the product in (5.14); the denominator follows by u->-u. This is independently checked symbolically in `scope_checks.py`. The remaining sine prefactor is exactly (5.17). The caution that scalar-factor identities do not establish a field-theory equivalence is appropriate. I did not independently verify the further N=2/N=1 folding identity (5.22); do not count the two author-script checks as evidence for that third identity.

## 11. Boundary and defect integration

**Verdict: wrong in boundary scope and one typesetting token; defect account is correct at the level claimed.**

The boundary source reports conservation of P_(1/2) and P_3 and calls the latter “strong evidence” for integrability. It then says “We conjecture” in introducing the quantum integrable model and reflection factor. These statements should not become a proof of a whole boundary hierarchy.

There is an independent check on the claimed construction. From Nepomechie (19), let
`N=alpha sin((phi-phi0)/2)-4 sin(phi/2)` and
`g^2=8[alpha cos((phi-phi0)/2)-4 cos(phi/2)+C]`.
Then `f=N/g` satisfies both `gf=N` and `g'=-2f`. This reproduces the supercharge constraint and the remaining integration constant. The additional P_3 constraint is what fixes it. This calculation is in `scope_checks.py`.

The two-parameter family uses a boundary fermion. That does not prove every integrable supersymmetric boundary needs one: the parameter-free cases without it are the starting point of the cited paper. Thus the universal sentence in chapter 16 line 335 and “forced by supersymmetry” in chapter 30 line 303 are too strong.

Replacements:

- Chapter 16 line 335, first sentence: “Nepomechie's two-parameter family of supersymmetric sine-Gordon boundary interactions uses a fermionic boundary degree of freedom (@sec-bsg-super).”
- Chapter 30 line 299, first two sentences: “The supersymmetric sine-Gordon theory adds a Majorana fermion to the boson, with the Yukawa coupling fixed by supersymmetry (@sec-super-shg). On the half-line, Nepomechie tests conservation of a boundary supercharge and a spin-3 charge; the latter is evidence for integrability.”
- Chapter 30 line 301, last sentence: “The result is a two-parameter family preserving the boundary supercharge and the tested spin-3 charge [@nepomechie2001, eqs. (16), (29)–(30)].”
- Chapter 30 line 303, replace the two sentences beginning “Whether a boundary degree”: “The extra boundary fermion allows the two-parameter family. Supersymmetry alone does not require it for every special boundary condition.”
- At the same line, replace the literal tab in `$2<TAB>imes2$` by `$2\times2$`.

The proposed lowest-breather reflection factor, its three parameters, the expected relation between them, and the cited algebraic properties accurately reflect the 2001 paper. Phrase the unknown relation historically (“was not determined there”), rather than imply a survey of all subsequent work.

For defects, the source 2015 introduction explicitly describes a “partially” type-II defect, and the older action contains `2 f1 partial_t f1`. The distinction from an auxiliary bosonic defect coordinate is correctly retained. The Lax pair, gauge matrix, modified energy/momentum, generating function and invariance of the total action are all present in the cited sections. An independent local check of the generating-function logic is: if `Psi_2=K Psi_1`, then the difference of logarithmic spatial derivatives of a chosen nonzero component is `partial_x log((K Psi_1)_i/(Psi_1)_i)`. Its time counterpart from the gauge equation is the defect boundary term in the integrated conservation law. This supports the stated mechanism, not a new proof of involutivity of the charges. No quantum transmission matrix is asserted, correctly.

## 12. Citation fidelity, classification and prior art

**Verdict: wrong in the unrestricted massiveness criterion and one integration sentence; unclear for an exhaustive all-citation audit.**

Evans–Madsen 1997 Proposition 6 explicitly begins its two cases with “Affine diagrams with no grey nodes” and “Affine diagrams which are entirely fermionic”. Chapter 16 line 107 drops “affine”. The surrounding construction includes finite systems, so this is a substantive missing hypothesis. Independent counterexample: one black simple root gives `W=mu exp(beta a phi)/beta^2`, whose bosonic potential is a single exponential with no finite stationary vacuum. It is super-Liouville, not a purely massive theory.

Replace line 107's first criterion by:

> Within the affine root systems considered here, the theory is purely massive if and only if the Dynkin diagram has no grey node or all its simple roots are odd [@evans_madsen1997, Prop. 6].

The list in `prp-super-susy` agrees with Table 1(c), including the low-rank convention `D(1,m)=C(m+1)` used in that table. The rank-one C(2) example is covered through D(1,1)^(2), rather than requiring an m=0 entry. The claimed uniqueness is uniqueness of the purely fermionic simple-root system in that classification, not uniqueness of all root bases.

For `prp-super-bosonic`, summing the odd-root products in V gives precisely the B-plus/B-minus exponent sets. The source's sign conditions are more useful than “suitable signs”: take the even-root mu with the sign of its root square, and all odd-root mu of one sign. Its Proposition 1 supplies orthogonality and definiteness, and Proposition 2 supplies the free orthogonal complement. The chapter accurately reports them; it does not provide a new general classification proof. The explicit A, C and D computations above independently test the mechanism.

The Evans–Madsen 1996 attribution concerns (1,0) extensions and agrees with their abstract and argument for the classical series. To avoid including the surviving stress tensor in “none of the W-currents”, say “the additional higher-spin W-generators examined there”. Their conclusions are not a theorem about every possible fermionic extension. Chapter 16 line 11 is much broader than the specific evidence; preferably make it a description of the cited construction, not an unqualified current-literature assertion.

Ito–Zhu's squared-potential claim is independently verified by eliminating the second component of the lower 2-by-2 block in their (3.7):

\[
(\partial_x+\ell/x)y=z,\quad(\partial_x-\ell/x)z=p^2y
\quad\Rightarrow\quad
[-\partial_x^2+\ell(\ell+1)/x^2+p^2]y=0.
\]

This gives the squared potential and the centrifugal term of (3.10), with `p=x^(2M)-E`. Their conclusion does not establish a general massive quantum correspondence. Both chapter 16 and the chapter 15 integration are appropriately limited. See [Ito–Zhu](https://arxiv.org/abs/2206.08024).

Chapter 6 line 185 says everything in Parts III–IX concerns ordinary affine Lie algebras while pointing to a superalgebra chapter in Part IV. Replace that sentence by:

> The root-system formulas of this chapter concern ordinary affine Lie algebras. The graded construction and Toda theories based on affine Lie superalgebras are treated in @sec-super-toda.

The new material is largely a literature review and has no new theorem labels. However, the statement “no Appendix E impact” is not fully justified under the supplied AGENTS.md, which explicitly requests updating the original-results list when corrections to papers are added. Appendix E's current introduction says corrections are listed separately in Appendix D. This is a repository-policy inconsistency to resolve editorially; it is not a new mathematical result. Do not silently treat the prose in Appendix E as overriding the explicit working rule. The corrections should remain labelled identities or misprints, not promoted to theorems just because this pass checked them.

Source-fidelity limits: I have checked the relevant passages in Evans–Hollowood 1991/1992, Evans–Madsen 1996/1997/1999, both DGPZ papers, Penati–Zanon, Gandenberger, Ito–Zhu, Nepomechie and the two defect papers. I have not independently read every cited foundational or scattering paper (Olshanetsky, Kac, van de Leur, Shankar–Witten, Schoutens, Ahn, Witten–Olive and Hollowood–Mavrikis), nor audited every old citation in the edited paragraphs. Their attribution must not be counted as fully independently certified by this report. No unsupported priority or novelty claim should be added on its strength.

## 13. Bibliography

**Verdict: unclear for the complete set; the spot-checked entries agree.**

The following independent metadata comparisons succeeded:

| Entry | Evidence and result |
|---|---|
| `evans_hollowood1991` | Original paper title page: NPB 352 (1991), 723–768. Erratum location NPB 382 (1992), 662 confirmed bibliographically; its content remains unread. |
| `dgpz1991a` | Original letter title page: PLB 256 (1991), 164–172, four authors in the stated order. |
| `dgpz1991b` | Original NPB title page: 359 (1991), 125–167, stated title and authors. |
| `evans_madsen1997` | [ArXiv record](https://arxiv.org/abs/hep-th/9703065): NPB 503 (1997), 715–746; DOI matches; v2 is dated 1 March 2000. |
| `ito_zhu2022` | [ArXiv record](https://arxiv.org/abs/2206.08024): authors, title and DOI match; revised published version is v2. |
| `nepomechie2001` | [ArXiv record](https://arxiv.org/abs/hep-th/0103029): PLB 509 (2001), 183–188; DOI matches. |
| `gomes2006` | [ArXiv record](https://arxiv.org/abs/hep-th/0601014): J. Phys. A 39 (2006), 7471–7483; DOI matches. Published-title difference is corroborated by the 2015 paper's reference list and is already explained in the entry. |
| `aguirre2015` | [ArXiv record](https://arxiv.org/abs/1412.2579): JHEP 02 (2015) 175, four authors and DOI match. |
| `gpz1991` | The original DGPZ letter's reference list supports PLB 253 (1991), 357, unlike the NPB paper's incorrect reference [16]. Full GPZ paper not independently read. |

For van de Leur, searches consistently locate DOI `10.1080/00927878908823823` and pages 1815–1841, volume 17(8), 1989. The publisher page returned 403 and the Crossref endpoint was unavailable through the browsing tool. I therefore do not mark all its metadata independently verified from the publisher. The existing start-page-only entry is not demonstrably wrong; adding the page range and DOI is a proposed bibliographic completion subject to checking.

The other new entries were compared where they appear in the primary papers' bibliographies, but that is weaker than a publisher/INSPIRE record check of every field. I did not reproduce the drafter's claim of a full 23-entry metadata audit. No fabricated metadata correction is proposed.

## Other observations — unverified or editorial

- The two introductory labels “not known beyond rank one” should be dated or scoped to the specific extensions discussed, unless a current literature survey is supplied. I have not established a current universal absence result.
- The purported exact fixed mass ratios deserve a source or a conjectural qualifier even though the one-loop and proposed S-matrix claims are well supported.
- The source's rho orientation resolves the first-order sign test; a check script which fixes rho>0 but calls the opposite sign intrinsically wrong overstates its result.
- The drafting report contains stale statements that the PDF render was not run alongside a later detailed PDF-render account. This report does not resolve that process-history inconsistency and did not render the book.

## Claims not independently completed

1. Full fermionic LSZ/crossing derivation from the physical charge-conjugation convention and ordered asymptotic states.
2. Independent loop-integral calculations for GPZ mass shifts, the DGPZ anomalous-threshold coefficients, and the complete Penati–Zanon/Evans–Madsen quantum-current coefficient systems.
3. All-parameter proof of the C-family proposed heaviest-particle amplitudes; the source itself excludes them from its direct bound-state bootstrap.
4. A presentation-level isomorphism between the entire generated loop superalgebra and the named twisted affine algebra; the needed relations and source identification are checked.
5. Full general-rank classification proofs and every exceptional-family folding map, beyond the explicit root, mass and signature calculations.
6. The Evans–Hollowood 1991 erratum, the original GPZ mass paper, the full quantum N=2 coupling dictionary and the further Gandenberger N=2/N=1 folding identity.
7. Exhaustive verification of every historical citation and every bibliographic field of all new entries.

These limitations do not invalidate the independently checked calculations. They prevent interpreting this report as certification of stronger claims than the evidence above supports.

Final preservation check: the chapter, Appendix D and bibliography retain their initial SHA-256 hashes. Another session advanced the branch to `d119e26633234d6679b0a9995e452cd1ce7c56ec` (“Chapter 9: point to the superalgebra chapter”) during this review. Its changes are confined to chapter 9 and three drafting/brief records; all book files in this report's stated review scope are identical to `9820f47`. The report remains against that reviewed commit. I made no tracked-file edits. `git status --short` reports only this new report; the referee scratch directory is ignored.

## Integration

Integrated 2026-10-24 by Claude (Opus 5.5) on branch `claude/superalgebra-toda`, on top of the merge of main (61a4626); the report and its scripts were first committed as received (b5bb638), with the scripts in `briefs/referee_superalgebra_toda_scripts/`. Rerun from a scratch directory, the three scripts reproduce their recorded PASS lines exactly (34, 44, 8; no FAIL). Following Rule 1, every item below was checked again by reading the source passage (the extracted texts the referee kept in `referee_scratch/superalgebra_toda/sources/`, not committed), by a derivation written out, or by a new check in the chapter's scripts. The referee's scripts were not used as evidence. No status label changed.

| Item | Verdict | How it was checked | Files changed |
|---|---|---|---|
| 1. Component Lagrangian, signs, SSG truncation, mass exercise | verified, no change | agrees with the drafting derivation and `components.py` | — |
| 2. Lax pair; algebra identification | verified; one sentence modified | the osp(1\|2) obstruction shows only that the obvious loop algebra fails; the identification rests on Evans–Madsen 1997, Table 4 ("A(1,0)^(2) ≅ C(2)^(2)") | ch. 16 (`sec-super-shg`) |
| 3. Supercurrents and the Evans–Madsen 1999 correction | verified, no change | as in the drafting report; Appendix D entry kept | — |
| 4a. DGPZ operator, braid YBE, unitarity, tree sign of γ, poles, bootstrap | verified, no change | agrees with `dgpz_smatrix.py` | — |
| 4b. Crossing phase | modified | DGPZ App. A: $u=(e^{\theta/2},e^{-\theta/2})$, $v=C\gamma^0u=(e^{\theta/2},-e^{-\theta/2})$, so $u(\theta+i\pi)=iv(\theta)$; this is a spinor identity, not a derivation of crossing for ordered amplitudes. The text now says the amplitudes use DGPZ's crossing convention | ch. 16 (`sec-super-dgpz`) |
| 4c / 9. "At one loop" for the mass ratios; fixed Coxeter number belongs to the proposed S-matrices | verified | GPZ and DGPZ claim one-loop mass ratios only; the exact S-matrices are conjectures | ch. 16 (intro, `sec-super-fixed-masses`, comparison paragraph), ch. 14 Notes, ch. 12 Notes |
| 5. DGPZ (B.1) free-fermion condition | verified, no change | as in the drafting report; Appendix D entry kept | — |
| 6a. "$-2\epsilon_1$" read as in both potentials | verified | eliminating the odd roots: $A$ has $-\epsilon_1-\epsilon_2$ and $2\epsilon_n$, $C$ has $-2\epsilon_1$ and $2\epsilon_n$; the $C$ cross term has $(-\epsilon_1)\cdot\epsilon_n=0$. New check of the exponent sets in `components.py` (n = 2, 3, 4) | ch. 16 (`sec-super-fixed-masses`) |
| 6b. Heaviest $C^{(2)}(n+1)$ boson is not a fermion–fermion bound state | verified | DGPZ NPB after (3.17): "These results hold for a,b = 1, 2, …, n−1 while the amplitudes involving the boson φn cannot be determined by the bootstrap since the corresponding pole is actually absent in the fermion-fermion amplitude"; the letter says the same after its (26). Kinematics: $m_n^2=8=4m_\psi^2$ puts the pole at $s=4m_\psi^2$, i.e. $\theta=0$. New check in `components.py` | ch. 16 (`sec-super-dgpz`, `cnj-dgpz-smatrix`) |
| 6c. Annihilation matters only for $n=3k$ | verified | DGPZ NPB Sec. 3: "if n = 3k then m_k = √2 … we will assume that there is no amplitude for a pair of fermions to annihilate into a pair of bosons"; $8\sin^2(\pi a/2n)=2$ iff $a=n/3$ (new check, n = 2..12) | ch. 16 (`sec-super-dgpz`) |
| 7a. N = 2 obstruction; admissible α | verified | the bipartite argument; α ∈ ℂ∖{0, −1} (Evans–Madsen 1999, Sec. 4) | ch. 16 (`sec-super-n2`) |
| 7b. sl/psl/gl decoupling sentence | verified, replaced | own computation: on $K=\operatorname{diag}(1_m,-1_m)$ every simple root $\epsilon_i-\delta_j$ takes the value ±2, and $\operatorname{str}(IK)=2m$, $\operatorname{str}I^2=\operatorname{str}K^2=0$, so the extra field enters the exponentials; Evans–Hollowood 1992, Sec. 1: "We may therefore consistently set Φ′ = 0" (a truncation) | ch. 16 (`sec-super-n2`) |
| 7c. Central charge singular at $a=m/2$ | verified, formula replaced | $(1+i\sigma)/(1-i\sigma)=e^{2ix}=ie^{i\pi a/m}$ for $\sigma=\tan x$, $x=\frac{\pi a}{2m}+\frac\pi4$. `eq-n2-central` now uses $c_a=ie^{i\pi a/m}$; `n2_solitons.py` checks the identity, adds $m=2$ and $a=m/2$, and multiplies (5.7) by $\cos x$ | ch. 16, `n2_solitons.py` |
| 7d. Sign of λ is an orientation convention | verified, Notes rewritten | Evans–Hollowood (5.3) fix only "ρ²(1 − v²) = 4µ² sin²(πa/m)". At fixed $z=\rho(x-vt)$, reversing ρ flips $\partial_+\phi$ and not $F$, so λ flips. `n2_solitons.py` now checks both orientations | ch. 16 (`sec-super-notes`, `sec-super-n2`), `n2_solitons.py`, conventions brief |
| 7e. "without the factor S_min" for m > 2 | verified | Evans–Hollowood after (5.10): "(but without the factor Smin, since in this case S N=2 [m] would provide the necessary pole structure)" | ch. 16 (`cnj-n2-sine-gordon`) |
| 7f. Missing $J$ in the second line of Evans–Hollowood (5.5) | verified, not integrated | the printed $\varepsilon_-'\partial_-\phi$ should be $\varepsilon_-'J\partial_-\phi$ by symmetry with the first line. The book does not use (5.5), so no correction is added | — |
| 8a. Sign of the $D(2,1;\alpha)$ cubic Yukawa term | verified, corrected | Evans–Madsen 1999, Sec. 3.3 prints "L_yuk = −4µ√(α(1+α))(ψ̄1ψ2φ3 + ψ̄1ψ3φ2 + ψ̄2ψ3φ1)". Own derivation: (3.32) to first order gives $+2m\bar\psi_1\psi_2\phi_3$ etc.; reversing the overall sign (negative kinetic terms) and normalizing by $\sqrt{1+\alpha},1,\sqrt\alpha$ gives $-2m\sqrt{\alpha(1+\alpha)}$. New check in `components.py` (the old one tested only the modulus) | ch. 16 (`sec-super-n2`), `components.py` |
| 8b. α ≠ 0, −1; α > 0 before the masses; foldings and the α → 0 limit | verified | Evans–Madsen 1999: α ∈ ℂ∖{−1, 0}; "we can without loss of generality assume that α is positive"; Sec. 4: $C(2)^{(2)}$ by folding, $B(1,1)^{(1)}$ by folding $D(2,1)^{(1)}$ with $\alpha_2\leftrightarrow\alpha_3$, $A(1,1)^{(1)}$ as the singular limit with $\alpha_2=-\alpha_1$, $\Phi_3=0$ | ch. 16 (`sec-super-n2`, `sec-super-quantum-currents`) |
| 9. Quantum-integrability claims | verified, "every admissible α" | as reported in the sources; the book's `gpz1991` reference is the right one | ch. 16 |
| 10. Gandenberger: $(-1)^{ab}$ after fusion | verified, corrected | his (3.12) and (5.7) fuse the (1,1) factors with the same shifts $(a+b-2j-2k+2)/2h$ at λ = 1/h, so $G_{a,b}=(-1)^{ab}F_{a,b}$; new fused check in `gandenberger_scalar.py` ((a, b) = (1,2), (2,1), (2,2), (1,3)) | ch. 23 Notes, Appendix D, `gandenberger_scalar.py` |
| 11a. Boundary scope | verified, modified | Nepomechie: "The fact that P3 is an integral of motion is strong evidence that the model is integrable" and "We conjecture that the boundary SSG model is integrable"; the parameter-free boundaries of Inami et al. have no boundary fermion | ch. 16 (`sec-super-connections`), ch. 30 (`sec-bsg-super`) |
| 11b. Literal tab in `$2\times2$` | verified, fixed | the file contained a tab before "imes2" | ch. 30 |
| 11c. Defect account | verified, no change | — | — |
| 12a. Purely massive criterion needs "affine" | verified | Prop. 6: "(a) Affine diagrams with no grey nodes … (b) Affine diagrams which are entirely fermionic"; one black root gives super-Liouville | ch. 16 (`sec-super-superspace`) |
| 12b. Intro and "beyond rank one" | verified, rewritten | Evans–Hollowood 1992, Sec. 1: "there is no obvious way to supersymmetrize a given bosonic Toda theory whilst maintaining integrability"; "rank one" was the drafter's gloss | ch. 16 (intro, `sec-super-superspace`) |
| 12c. "none of the W-currents" | verified | Evans–Madsen 1996: T and G survive; the higher-spin generators do not | ch. 16 |
| 12d. Sign conditions in `prp-super-bosonic` | verified | Evans–Madsen 1997, Prop. 2: positive/negative definite if μ_i ≷ 0 for α_i² ≷ 0 on even roots and all odd-root μ_i of one sign | ch. 16 |
| 12e. Chapter 6 sentence about Parts III–IX | verified | it contradicted the Part IV chapter | ch. 6 |
| 12f. Appendix E impact | rejected | Appendix E says "Corrections to published papers are listed separately, in @sec-app-corrections", and no part lists a correction as an item; the AGENTS.md clause concerns items affected by a correction. Left to the editor if the policy should change | — |
| 13. Bibliography: van de Leur | verified, completed | Crossref record of DOI 10.1080/00927878908823823: Commun. Algebra 17(8), 1815–1841 (1989) | `references.bib` |
| Other: drafting report's stale "PDF render was not run" | verified, fixed | the report also records the PDF render | `briefs/draft_superalgebra_toda_report.md` |

**Scripts.** `components.py` 31 → 39 checks (Yukawa sign from (3.32); exponent sets of the A and C potentials; heaviest-boson threshold and the n = 3k degeneracy); `n2_solitons.py` 7 → 9 (c_a identity; m = 2 and a = m/2; both orientations of ρ); `gandenberger_scalar.py` 2 → 3 (fusion). All six scripts pass (9, 39, 11, 27, 9, 3), and `code/super_toda/README.md` is updated.

**Not integrated, and why.**
- Evans–Hollowood 1992 (5.5): a missing $J$ in a formula the book does not use (7f).
- Appendix E: rejected as above (12f).
- The open points the report lists under "Claims not independently completed" (crossing from charge conjugation, loop integrals, the presentation of the loop superalgebra, the Evans–Hollowood 1991 erratum, the further Gandenberger identity (5.22)) remain open. The chapter's labels ([Conjecture] for both S-matrices, [Established] for the two literature propositions) already reflect this, and the text now states the crossing convention as a convention.
