# Referee report: the boundary sine-Gordon Dirichlet sign

- Date: 2026-10-25.
- Model: GPT-6 (Codex).
- Commit reviewed: `ed0585ae78572d7ada8986337f89ec685c0b6462`.
- Worktree: `.claude/worktrees/eager-bun-c1cfa6`. Its HEAD was already `b1d214999bf5c6a06a0195ba20c976e2892f7817`; this report reviews the requested commit, using files extracted with `git show`, not the later HEAD. All book line numbers below refer to **ed0585a**.
- Independent evidence: `briefs/referee_bsg_dirichlet_sign_scripts/{independent_checks.py,algebra_checks.py,analytic_flux_check.py}`, with matching `.out` files preserved alongside them. See that directory’s README for reproduction commands and source provenance. No book module is imported by these scripts. Downloaded arXiv source archives and extracted TeX are under `referee_scratch/bsg_dirichlet_sign/src/`; the reviewed snapshot is under `snapshot/`.

**Recommendation: retain the negative Dirichlet sign and the same-node pairing, with revisions.** The charge argument, an independently constructed exact classical scattering solution, and the current OPE agree. The exact solution makes the classical phase check substantially stronger than the supplied leapfrog calculation. This does not prove the finite-coupling UV–IR conjecture or establish priority. The main required revisions are the direct/crossed-channel wording, the missing ground-state subtraction in the phase argument, the Neumann branch qualification, and the Appendix E omission. There is also a stale sign in a script's printed description. The prior-art verdict remains **unclear** in the specific sense explained in item 6.

During the original referee pass, no existing book, code, lecture, or index file was edited, and no rendering or mutating git command was run. The subsequent editor-requested report PR archives this report and its evidence and updates the briefs index; book integration remains pending.

## Follow-up for Claude after PR #63

**Status: report received; integration pending.** This report reviews ed0585a. A subsequent read-only audit checked PR #63's final head `af5cc340c63fcba654aafee43b01fd2b57fee8cb` and merge commit `9f0be486803d34be628e56c01e4958d802fdc122`; their trees are identical. The merged chapter 30 is unchanged from the reviewed version. The report was not included in that PR.

Claude: do not spawn subagents. Independently verify each finding under Rule 1, apply accepted changes, and add an `## Integration` section recording verified, modified and rejected items. The editor requested this PR as a handoff for later work; this PR does not integrate the findings into the book.

Outstanding at the merge of #63:

- Item 1: distinguish zero direct-channel energy from zero crossed-channel momentum (chapter 30, line 108).
- Item 2: state the ground-state subtraction in the scattering-action derivative (chapter 30, line 116). The simulation already subtracts it. The independent exact-flux derivation and archived scripts provide additional evidence to verify.
- Item 3: optionally state the analytic branch of the arccos expression (chapter 33, line 273).
- Item 4: explain the negative Neumann representative approached for positive phi_0 and nonnegative vartheta, and restrict the spectral invariance statement to unlabelled energies (chapter 30, lines 131 and 142).
- Item 5: qualify the Mattsson–Dorey inheritance claim and the scope of the unchanged spectrum; make the UV–IR correction category precise (Appendix D, lines 90 and 109 in the merged version).
- Item 6: retain the unresolved priority and lattice-dictionary qualifications. The prior-art investigation does not establish novelty.
- Item 7: fix the printed `-+` to `+-` in item (e) of `code/lectures/defects_l9_fusion.py`, line 688. The updated item (h) is correct.
- Item 8: add the conservative Appendix E entry or explicitly reconcile the repository's correction-list rule with the appendix policy. Keep the quantum relation labelled [Conjecture].

The line numbers in this checklist refer to merge commit 9f0be486. The numbered report sections below retain their original ed0585a line numbers. PR #63 also corrected the description of Delius–MacKay's inverse q as a convention; that separate change does not resolve the outstanding items above.

## 1. Charge argument and the Dirichlet sign

**Verdict: correct for the sign; wrong for the phrase “absorbed at zero momentum” in the direct channel. The exact linear relation remains conjectural at finite coupling.**

### Source check

[Ghoshal–Zamolodchikov](https://arxiv.org/abs/hep-th/9306002), (2.8), integrates over the negative half-line. The text around (3.16) explicitly has Minkowski space with `x<0` and incoming positive rapidities. Their statement immediately before (5.3) is:

> The soliton (antisoliton) carries a positive (negative) unit of “topological charge”

The integral in (5.3) is the increasing-field charge used in the brief. Equations (5.12) and (5.26) give the diagonal prefactors with the same signs as the book. Equation (5.27) has the positive sign. The emission wording quoted in Appendix D is present in the source (TeX lines 1633–1634); it is not a quotation error. The source's (3.41)–(3.42) concern a zero-momentum particle in the **rotated boundary-state channel**.

### Independent static computation

Write `a = beta*phi_0`, with `0<a<pi`, and `u=beta*phi`. Finite energy and the static equation give

\[
\tfrac12u_x^2=m_0^2(1-\cos u).
\]

A nonconstant solution cannot cross another vacuum at finite x: at a vacuum both u_x and u_xx vanish, and uniqueness of the ODE gives the constant solution. Consequently there are exactly two real, regular, finite-energy static solutions with this boundary value:

\[
 u_g(x)=4\arctan[\tan(a/4)e^{m_0x}],\qquad
 u_e(x)=2\pi-4\arctan[\cot(a/4)e^{m_0x}].
\]

They approach 0 and 2pi respectively at negative infinity. Their energies, excluding the common Dirichlet boundary constant, are

\[
 E_g=\frac{4m_0}{\beta^2}(1-\cos(a/2)),\quad
 E_e=\frac{4m_0}{\beta^2}(1+\cos(a/2)),\quad
 E_e-E_g=M_{\rm cl}\cos(a/2),\quad M_{\rm cl}=\frac{8m_0}{\beta^2}.
\]

Their charges are `a/(2pi)` and `a/(2pi)-1`. At a=pi their energies coincide. A soliton approaching from the left has its left vacuum one step **below** its right vacuum. An antisoliton approaching a ground-state boundary instead has vacuum 2pi on its left and 0 on its right. Absorbing it leaves precisely u_e. Thus the first excited state is in the antisoliton channel.

Independent integration in `algebra_checks.py`, with m0=beta=1:

```text
a=0.6: E_g=0.1786540434975759, E_e=7.821345956502425
        gap=7.642691913004849, 8*cos(a/2)=7.642691913004848
a=1.4: E_g=0.9406312508620462, E_e=7.0593687491379535
        gap=6.118737498275907, 8*cos(a/2)=6.118737498275908
a=pi:  E_g=3.9999999999999973, E_e=3.9999999999999996
```

### Pole and scalar-factor computation

For generic positive eta, at `nu_0=(eta-pi/2)/lambda`,

\[
\cos(\lambda\nu_0-\eta)=0,\qquad
\cos(\lambda\nu_0+\eta)=\sin(2\eta).
\]

Thus P-minus cancels the first scalar pole and P-plus retains it, barring exceptional coincident zeros. My direct Gamma-product computation at lambda=2+sqrt(2), eta=2.2 gives:

| u−nu_0 | absolute P-plus, without the common R0/cos eta | absolute P-minus, same normalization |
|---|---:|---:|
| 0.001 | 107.7862601 | 0.3862974 |
| 0.0001 | 1079.0408544 | 0.3871019 |
| 0.00001 | 10791.5823751 | 0.3871824 |

This verifies pole order by scaling, not a residue threshold.

Evenness of sigma is an exact identity, despite the asymmetric appearance of its product. Set z=eta/pi and t=lambda*u/pi. On changing z to −z the symmetric Gamma numerators cancel. The ratio of the first N product terms telescopes to

\[
 \frac{\Gamma(\frac12-z-t)}{\Gamma(\frac12+z-t)}
 \frac{\Gamma(\frac12+z+2N\lambda-t)}{\Gamma(\frac12-z+2N\lambda-t)},
\]

divided by the same expression at −t. The large-N endpoint tends to 1 in this quotient. Gamma reflection turns the remaining factor into `cos(eta-lambda*u)/cos(eta+lambda*u)`, cancelling the ratio of the explicit cosine prefactors. Hence sigma(−eta,u)=sigma(eta,u). A finite cutoff is not exactly even; the numerical ratio at N=20000 is `0.999999999841 + 0.0000178325283 i`, consistent with the endpoint correction.

The physical pole must therefore use negative eta when phi_0>0 on x<=0. Its edge value fixes `eta=-pi(lambda+1)/2` at `phi_0=pi/beta`. Assuming linearity gives the proposed negative relation. On x>=0 the charges of the corresponding two static states are `-a/(2pi)` and `1-a/(2pi)`; the absorbed particle has charge +1. The positive relation is appropriate there, with rapidity measured toward the boundary.

[Saleur–Skorik–Warner](https://arxiv.org/abs/hep-th/9408004) explicitly use `[0,infinity)` in section 2. Near the end of section 3 they write:

> the kink sitting in the middle and the incoming one are of opposite topological charges

This is consistent with the parity argument. Their time-delay calculation alone does not fix the rapidity-independent relative phase; see item 2.

The crossed-channel emission language can explain how a sign confusion arose, but does not repair GZ's simultaneous direct-channel definitions. For 0<a<pi the pole is strictly inside the physical strip and its charge is fixed by direct-channel absorption. At the edge the direct-channel energy tends to zero, whereas the momentum tends to iM. Calling a particle a soliton after crossing cannot turn this direct-channel antisoliton pole into P-plus without changing the charge or amplitude dictionary.

### Replacement text

At `part7-half-line/30-boundary-sine-gordon.qmd:108`, replace the clause beginning “so one amplitude must have” by:

> so one amplitude must have a boundary-bound-state pole whose direct-channel energy tends to zero, at $\theta=i\pi/2$. After the channel rotation this corresponds to a zero-momentum particle in the boundary state, as in (-@eq-one-point-residue).

Retain the negative equation at line 111 and the explicit linearity assumption.

## 2. Semiclassical check

**Verdict: correct sign, integral, and phase convention; the displayed action variation needs a ground-state subtraction.**

### Action and scattering-phase convention

For a finite time interval, the on-shell variation of the bulk action has the spatial endpoint term `-integral phi_x(0,t) delta phi_0 dt`. Equivalently, at finite boundary mass, stationarity removes the implicit field variation and `dS/dphi_0=-integral B_phi0 dt`; the boundary condition gives `B_phi0=phi_x(0,t)`. The two calculations agree in the Dirichlet limit.

However, the static ground-state slope is `phi_x^g(0)=2m0 sin(a/2)/beta`. The unsubtracted integral over all time diverges. For the scattering action with the static background and free asymptotic propagation subtracted, at fixed incoming/outgoing scattering data,

\[
 \frac{\partial S_{\rm sc}}{\partial\phi_0}
 =-\int_{-\infty}^{\infty}
 [\phi_x(0,t)-\phi_x^g(0)]\,dt.
\]

Forward time evolution is exp(−iHt). Time slicing its configuration-space kernel gives exp(i integral(p*dot(q)−H)dt), hence exp(+iS_sc), not exp(−iS_sc). This fixes the sign independently of naming a soliton. The difference of two actions is finite without explicitly displaying the identical background subtraction, but their sum is not.

### Independent exact scattering solution

I used the eight-term Hirota polynomial, independently imposing the boundary value and checking the sine-Gordon equation, rather than rerunning or adapting the author's leapfrog solver. In units m0=1, let

\[
 g=\cosh\theta,\quad h=\sinh\theta,\quad v=\tanh\theta,\quad
 w=\tanh^2(\theta/2),\quad r=\tan(a/4),
\]
\[
 D^2=\frac{w+r^2}{v^2w(1+r^2w)},\qquad
 z_1=sD e^{gx-ht},\quad z_2=sD e^{gx+ht},\quad z_3=(r/w)e^x,
\]

where s=+1 is an incoming increasing-field soliton and s=−1 an antisoliton. Define

\[
 \tau=1-v^2z_1z_2-wz_3(z_1+z_2)
       +i[z_1+z_2+z_3-v^2w^2z_1z_2z_3].
\]

The continuous field is `u=4 arg(tau)-2pi*s`, with the argument branch chosen continuously from the incoming configuration. It has u(0,t)=a, the ground-state wall profile before and after reflection, and a particle of charge s. Substitution gives the PDE residual reported below. This is the standard three-soliton construction of SSW, but the parameters, boundary constraint, derivatives, and integrations here were computed independently.

Direct differentiation simplifies the boundary slope to

\[
 u_x(0,t)-2\sin(a/2)
 =\frac{4h\sqrt{1-k^2}}{s\cosh(ht)+k},\qquad
 k=\frac{\sin(a/2)}{\sqrt{h^2+\sin^2(a/2)}}.
\]

The elementary substitution z=tanh(ht/2) now gives, for |a|<pi and theta>0,

\[
 I_s:=\int [u_x(0,t)-2\sin(a/2)]dt
 =4\pi s-8\arcsin k
 =4\pi s-8\arctan\frac{\sin(a/2)}{\sinh\theta}.
\]

In particular `I_+−I_−=8pi`, exactly. Restoring beta gives the integral 8pi/beta in the brief; restoring m0 rescales the time and slope inversely and leaves the integral unchanged. Since `d(delta_s−delta_antis)/da=−8pi/beta^2` and `P-plus/P-minus ~ exp(2i eta)`, the semiclassical relation is `eta=−4pi*a/beta^2`, with the constant fixed by charge conjugation at a=0. This is leading semiclassical information, not a proof excluding finite-coupling corrections.

My independent numerical integration checks the analytic result as follows:

| theta | a | I-plus | I-minus |
|---:|---:|---:|---:|
| 0.4 | 0.6 | 7.5770276771 | −17.5557135516 |
| 0.4 | 1.4 | 4.5408621031 | −20.5918791256 |
| 0.4 | −0.9 | 19.0785362647 | −6.0542049640 |
| 0.8 | 0.6 | 9.9965430618 | −15.1361981669 |
| 0.8 | 1.4 | 7.5458929200 | −17.5868483087 |
| 0.8 | −0.9 | 16.2097964588 | −8.9229447699 |
| 1.3 | 0.6 | 11.1881612775 | −13.9445799512 |
| 1.3 | 1.4 | 9.6660016476 | −15.4667395811 |
| 1.3 | −0.9 | 14.5721062855 | −10.5606349432 |

The velocities are approximately 0.380, 0.664 and 0.862. Maximum errors: PDE `5.11e-15`; boundary value modulo 2pi `1.54e-15`; closed boundary-slope formula `5.33e-15`; integral difference from 8pi `3.94e-13`. The analytic branch choice above fixes the actual boundary value, not merely its value modulo 2pi.

### Sign-blind common phase and time delays

Differentiating the displayed Gamma product with digamma functions gives, at eta=(lambda+1)a/2>0,

\[
 \partial_\eta\arg\sigma(\eta,-i\theta)
 \longrightarrow -\frac{I_++I_-}{8\pi}
 =\frac2\pi\arctan\frac{\sin(a/2)}{\sinh\theta}.
\]

At theta=0.8, a=1.4, the target is 0.399516920868. Removing the leading 1/N cutoff error with N=4000 and 16000 gives:

```text
lambda= 301.414213562: 0.400352294837
lambda=1001.414213562: 0.399768130359
lambda=3001.414213562: 0.399600714321
```

The error decreases as 1/lambda. At the first coupling, raw N=1000,4000,16000 values are 0.400271238, 0.400332031, 0.400347229. Thus the comparison is not an unchecked finite-product artifact. The sign would fail for exp(−iS_sc). For negative a, evenness of sigma and charge conjugation give the corresponding odd derivative.

Both signs s use the same D and the same magnitudes of all asymptotic position coefficients. Hence their classical time delays coincide. Equivalently, their leading phase difference is independent of theta, so its energy derivative vanishes. Classical delays cannot fix this sign.

### Replacement text

At `part7-half-line/30-boundary-sine-gordon.qmd:116`, replace the sentence beginning “Semiclassically the amplitude” by:

> Semiclassically the amplitude is $e^{iS_{\rm sc}}$, where $S_{\rm sc}$ is the on-shell action with the static ground-state background and free asymptotic propagation subtracted. At fixed scattering data, $\partial S_{\rm sc}/\partial\phi_0=-\int dt\,[\partial_x\phi(0,t)-\partial_x\phi_g(0)]$, with $\partial_x\phi_g(0)=2m_0\sin(\beta_{\rm SG}\phi_0/2)/\beta_{\rm SG}$. The background cancels in the difference of the soliton and antisoliton integrals and must be subtracted separately in their sum.

The existing check script already performs this subtraction; the defect is in the prose.

## 3. Node pairing and first-order expansion

**Verdict: correct, on the stated analytic branch and with the book's chirality and q conventions. It is an independent consistency check, not a proof of the full UV–IR relation.**

Let A=eta/(lambda+1), B=vartheta/(lambda+1), a=beta*phi_0 and rho=tilde(mu)/tilde(mu_crit). Expanding the complex cosine gives

\[
 \cos(A-iB)=\cos A\cosh B+i\sin A\sinh B
 =\rho e^{-ia/2}=\mu_-,
\]

and its conjugate is mu-plus. Consequently, locally on the branch continued from Neumann,

\[
 \hat\epsilon_1=f(\mu_-),\quad\hat\epsilon_0=f(\mu_+),\qquad
 f(z)=-\frac{\cos[(\lambda+1)\arccos z]}{\sin\pi\lambda}.
\]

This is a branch-dependent analytic function for noninteger lambda+1, not a globally single-valued assertion about arbitrary independent complex couplings. The chapter's physical branch is sufficient.

### Independent current computation

With the canonically normalized chiral covariance `-log(z-w)/(4pi)`, the OPE of charges A and B has power `AB/(4pi)`. The power for `exp(-ia_current phi_R)` and the right-moving part of `exp(+i beta phi)` is

\[
 \frac{(-8\pi/\beta)\beta}{4\pi}=-2.
\]

It is +2 for the opposite bulk exponential. At a Neumann boundary phi=2phi_R=2phi_L, so `exp(i beta phi/2)` has the same double pole with either negative-exponent current. Its coefficient in the cosine boundary interaction is proportional to `exp(-i beta phi_0/2)`, namely C1. Thus b1, containing the two negative-exponent currents, meets C1. This calculation does not use the Dirichlet pole or the scattering solution.

The canonical string definitions also fix the charge, rather than relying on the current's subscript. Using `phi_R=(phi-integral^x pi)/2` and `[phi(x),pi(y)]=i delta(x-y)` gives the disorder commutators

\[
 [T,\phi_R]=i\beta/(4\pi),\qquad
 [T,\phi_L]=-i\beta/(4\pi).
\]

Therefore `[T,exp(-ia_current phi_R)]=2 exp(-ia_current phi_R)`. The right-moving charge raises T by two units on a one-particle doublet, exactly as E maps the antisoliton into the soliton. Its spin is `a_current^2/(8pi)-1=lambda`, so its rapidity factor is exp(lambda*theta). This independently confirms the particular chirality and representation identifications needed here.

[Delius–MacKay](https://arxiv.org/abs/hep-th/0112023), source `paper.tex`, lines 564–566, use z=(t+ix)/2 and Euclidean t=i*t_M. Their holomorphic field is therefore the book's left-mover. Their source current includes the exact fragment

```tex
J_\pm = : e^{\pm \frac{2i}{\hb}\varphi}:
```

and their boundary perturbation, lines 713–715, includes

```tex
\epsilon_+ e^{-i\hb
  \tilde{\phi}(0,t)/2}
```

Their Q-hat-plus combines the positive-exponent left and right currents. Moreover their q is the inverse of the book's q, so their q^T becomes the book's k0. Thus their Q-hat-plus corresponds to b0 after the stated normalization changes, and their epsilon-plus couples to C0. Merely comparing the plus subscript to the book's Q-plus would give the wrong dictionary.

### Taylor coefficient

Direct differentiation gives

\[
 f(0)=\frac1{2\cos(\pi\lambda/2)}=\epsilon_*,\qquad
 f'(0)=-\frac{\lambda+1}{2\sin(\pi\lambda/2)}.
\]

My finite-difference check at lambda=sqrt(2), e, 2+sqrt(2) has intercept errors at most 1.12e-16 and derivative errors at most 1.29e-10. This verifies the expansion in `prp-dm-first-order` and the changed lecture formula. It does not independently establish the proposed nonperturbative Neumann offset mechanism for higher rank.

No replacement is needed for the changed pairing or expansion. Optionally insert “on the branch continued from Neumann” before the arccos formulas at `part7-half-line/33-boundary-quantum-groups.qmd:273`.

## 4. General UV–IR relation and its tests

**Verdict: correct as an explicitly sign-adjusted conjecture; the Neumann bullet needs a branch qualification. The cited TCSA energy tests are blind to the global charge-conjugation choice.**

### Independent limits and symmetry argument

The two equations in item 3 are exchanged by simultaneously swapping C0,C1 and epsilon0,epsilon1. Changing only the proposed sign of the UV–IR dictionary, at fixed named Lagrangian couplings, exchanges the charged reflection channels. This is a concrete conjectural dictionary change. Symmetry and the Dirichlet limit do not uniquely determine all finite-mass corrections, so the existing [Conjecture] label must remain.

The four advertised limits follow without numerical fitting:

1. At rho=0, cos A=0 and sinh B=0 on the physical branch, so B=0 and A=+pi/2 or −pi/2. These give the identical Neumann reflection matrix and the same epsilon-star.
2. For rho tending to infinity with B tending to positive infinity, dividing the equations gives `tan A=-tan(a/2)`. In `|A|<=pi/2`, `|a|<pi`, this gives A=−a/2, and hence the negative Dirichlet relation.
3. At a=0, eta=pi/2 and vartheta=0, both coideal parameters vanish and rho=cos[pi/(2(lambda+1))]. The negative-eta representative is equivalent at vartheta=0.
4. With mu=m0^2/beta^2, `sin(beta^2/8)~beta^2/8`, so tilde(mu_crit) tends to 4m0/beta^2, giving C0=C1=1.

There is a minor parameter-chart issue: for a>0, rho>0 and vartheta>=0, the new second equation requires eta<0 in the fundamental range. Thus this branch approaches **negative** eta_N at Neumann. The positive Neumann representative quoted in the book is physically equivalent, but is not that branch's signed limit.

For the TCSA point, charge conjugation U_C sends `H(phi_0)` to `H(-phi_0)`. In a free-boson truncated basis respecting charge conjugation, the energy truncation commutes with U_C, so the truncated Hamiltonians are related by the same unitary similarity. Neumann at the other end is invariant. Equal Dirichlet ends map to equal opposite Dirichlet ends; the unlabelled energy sets are again the same. For scattering, let C be the 2x2 swap matrix. The displayed matrix obeys

\[
 K(-\eta,\vartheta)=C K(\eta,\vartheta) C,
\]

including its even scalar factor. The bulk S-matrix respects simultaneous charge conjugation, and Neumann K commutes with C. Hence the double-row operators are conjugate and their unlabelled eigenvalues agree. `algebra_checks.py` verifies the matrix identities exactly to floating precision and the sample double-row determinant to 4.56e-16.

[Bajnok–Palla–Takács](https://arxiv.org/abs/hep-th/0108157), `uv_ir3.tex:890`, state:

> only at one end of the strip.

At lines 1229–1230 they also state:

> identical boundary conditions are imposed.

These are the two numerical setups claimed in the brief. The source check and the similarity derivation establish blindness of these **energy comparisons**. I did not rerun their TCSA program. This argument does not assert blindness of charge-labelled sectors or signed matrix elements of a fixed charged operator.

### Replacement text

At `part7-half-line/30-boundary-sine-gordon.qmd:131`, replace the Neumann bullet by:

> - $\tilde\mu\to0$ gives $\vartheta=0$ and $\eta=\pm\frac\pi2(\lambda+1)$, two equivalent representatives of the Neumann reflection matrix. For $\phi_0>0$ and $\vartheta\ge0$, the corrected relation approaches the negative representative.

At line 142, replace the sentence beginning “Results at” by:

> Results at $\phi_0=0$, where $C_0=C_1$, and charge-conjugation-invariant observables cannot see the sign. These include the unlabelled half-line bound-state energies and the energy spectra in the TCSA comparisons of [@bpt2002]. The assignment of a level or a pole to a fixed topological-charge sector does change.

## 5. Corrections to the literature

**Verdict: the quoted signs and the charge-sign corrections are correct. Qualify “inherited” for the Mattsson–Dorey arithmetic error and specify which spectrum is unchanged.**

The evidence is the independent charge/pole computation in item 1 and the conjugation calculation in item 4, not the author's correction list. Here are the source-line checks for each new entry.

- **Ghoshal–Zamolodchikov (5.27).** `src/9306002/9306002.tex:1633–1643` contains the stated emission argument and the positive relation. In the source notation the exact equation is `\xi = {{4\pi}\over \beta}\varphi_0.` The inconsistency is internal to their x<0, positive-charge soliton, and P-plus conventions. *Error* is justified with those conventions held fixed. The solution of the reflection equation is not being refuted; its charge-resolved Lagrangian dictionary is being corrected. Retain the entry, with the direct/crossed-channel distinction of item 1.
- **Mattsson–Dorey section 2.2.** [Source](https://arxiv.org/abs/hep-th/0008071), `bsgdht.tex:378–383`, first defines the half-line charge as the boundary value minus the value at negative infinity and then prints:

  ```tex
  with the charge of the first excited state being \m{1-\frac{\beta
  \varphi_{0}}{2\pi}}
  ```

  Substituting their own asymptotic vacuum 2pi/beta in their charge integral gives a/(2pi)−1. This is a definite arithmetic sign error. The printed value is indeed the x>=0 answer. Their subsequent soliton/antisoliton assignment follows the positive GZ dictionary, but the source does not establish the historical cause of this particular arithmetic error. Calling the whole sentence simply “inherited” overstates that causal claim.
- **BPTT section 4.** [Source](https://arxiv.org/abs/hep-th/0106070), `bsgrevision.tex:833–857`, explicitly invokes the GZ Dirichlet limit and proposes the general relation from it and symmetry. The exact positive right-hand side of its second equation is:

  ```tex
  = & \frac{M_{0}}{M_{crit}}\sin \left( \frac{\beta \phi _{0}}{2}\right)
  ```

  Thus the book has not misquoted the sign, and “inherited from the Dirichlet limit” is supported by the source. The proposed replacement is correct for the book's charge-resolved dictionary, subject to the conjectural status of the general relation.
- **BPT (3.2).** `src/0108157/uv_ir3.tex:352` prints the same positive sine right-hand side, in its normal-ordered coupling normalization:

  ```tex
  = & \frac{\tilde{\mu }}{\tilde{\mu }_{\mathrm{crit}}}\sin
  ```

  Its first equation and normalization need no sign change. Their numerical energy tests do not choose the charged-channel sign, for the reason established in item 4.

The half-line energies depend on |eta| through the admissible pole angles; changing eta's sign exchanges P-plus and P-minus. In particular the first energy is `M cos[(|eta|-pi/2)/lambda]`, not a formula obtained by inserting negative eta into a pole list restricted to positive eta. Therefore “the spectrum is unchanged” is true for the **unlabelled half-line energy spectrum**, not the charge-resolved spectrum. The chapter's positive-eta spectrum discussion is usable by charge conjugation, but that convention must be applied when returning to phi_0>0 Dirichlet data.

### Replacement text

At `appendices/D-corrections.qmd:89`, replace the last two sentences of the new Mattsson–Dorey item by:

> They use the positive Dirichlet relation of Ghoshal and Zamolodchikov (5.27), so for $\varphi_0>0$ the soliton and antisoliton assignments to the boundary states must be exchanged. The unlabelled half-line bound-state energies are unchanged. *Charge-sign error; the channel assignment follows the inherited Dirichlet sign* (@sec-bsg-parameters).

At `appendices/D-corrections.qmd:109`, change the category to:

> *Inherited sign error in the charge-resolved UV–IR dictionary* (@sec-bsg-parameters; `code/boundary_sine_gordon.py`).

This more precise category preserves the distinction between a wrong sign assignment and a claimed proof or disproof of the entire conjecture. The new correction at `part7-half-line/30-boundary-sine-gordon.qmd:347` is supported by the same computations; it should use the same qualifications wherever it summarizes these entries.

## 6. Prior art and lattice/NLIE checks

**Verdict: unclear on priority and on a complete independent lattice confirmation. The specified sources were checked; none supplies an unambiguous counterexample to the negative sign in the book's fully fixed dictionary. Some use the positive sign on the other half-line, and one NLIE formula gives useful evidence for a particle-label reversal.**

This is not a finding that nobody noticed the sign previously. Targeted searches for a GZ Dirichlet sign correction did not identify an explicit erratum. Absence from these searches is not evidence of novelty.

### Sources checked

1. **Skorik–Saleur, hep-th/9502011.** Their model is explicitly on x>=0: source `9502011.tex:125` contains `\int_0^{\infty}`. At line 184 they write `\xi={4\pi\over\beta}\varphi_0`. Their boundary strings and positive-channel poles therefore use the sign appropriate to their half-line. The boundary-charge change is independently +1 by the static calculation in item 1. This supports the orientation explanation, but is not a previous explicit correction of GZ's x<=0 text. [Source](https://arxiv.org/abs/hep-th/9502011).
2. **Fendley–Saleur, hep-th/9402045.** The source has `\phi(0)=\phi(L)=\phi_0` and, after comparing the reflection kernels, `\xi={\pi \over 2}(t-1-H)`. It ties a lattice boundary parameter to xi. It does not, in the passages checked, independently fix the absolute increasing-field charge of the lattice excitation together with a signed phi_0(H) at each endpoint. Global reversal of the arrow/field dictionary remains relevant. It cannot by itself disprove the proposed sign. [Source](https://arxiv.org/abs/hep-th/9402045), `paper.tex:179–186,540–551`.
3. **LeClair–Mussardo–Saleur–Skorik, hep-th/9503227.** Their continuum boundary energy is expressed using the exact source terms `\cos(t-1)\xi` and `\cos(t-1)\xi'` (`9503227.tex:1529–1530`). Each term is even in its boundary argument. Independently applying xi->−xi to that expression leaves it unchanged, so this boundary-energy comparison cannot select the charge-resolved sign. I checked the boundary/TBA and continuum-energy passages; I did not redo their lattice construction. [Source](https://arxiv.org/abs/hep-th/9503227).
4. **Ahn–Bellacosa–Ravanini, hep-th/0312176.** They use a strip with phi-minus at x=0 and phi-plus at x=L, and print `H_{\pm}=p(1\mp\frac{8}{\beta}\phi_{\pm})`. The text immediately preceding this formula explicitly identifies xi with phi through GZ. Thus their printed dictionary is not wholly independent input for testing GZ. Their analysis also restricts the boundary fields to a regime without boundary bound states. The differing endpoint signs must be retained when comparing with a single right boundary. [Source](https://arxiv.org/abs/hep-th/0312176), `letter.tex:405–435`.
5. **Ahn–Bajnok–Palla–Ravanini, 0711.4101.** This is the most relevant additional check; the calculation below isolates what its ultraviolet formula does and does not establish. It contains the positive xi-plus/phi-plus assignment and calls the corresponding factor a soliton amplitude, but its winding formula is not consistent with simply identifying a positive-counting-function hole with increasing-field topological charge. [Source](https://arxiv.org/abs/0711.4101).
6. **Caux–Saleur–Siano, cond-mat/0306328.** Their two-boundary massless problem places a relative phase chi at the right boundary. They explicitly say the ground-state energy is an “even function” of chi (`NPBsub.tex:1041–1042`). Independently, phi->−phi sends chi->−chi in their Hamiltonian while keeping the left cosine invariant. Their ground-state-energy test therefore cannot select a global charge-conjugation convention. This observation does not make arbitrary charge-labelled two-boundary observables sign blind. [Source](https://arxiv.org/abs/cond-mat/0306328).
7. **Kormos–Takács, 0712.1886.** The requested authors' paper uses the “boundary Lee-Yang model as testing ground” (`bfftcsa_revised.tex:59`), not a charge-resolved sine-Gordon soliton test. I also checked the more directly relevant **Lencsés–Takács, 1106.1902**, on breather boundary form factors. Their non-diagonal matrix-element comparison uses absolute values; the exact source phrase is “only the absolute values can be compared” (`bsgfftcsa.tex:354–355`). This paper also tests a vacuum expectation value, so it must **not** be dismissed wholesale as modulus-only. Its interacting boundary is at x=0, with Neumann at x=L (lines 433 onward), so the positive GZ sign there is consistent with a bulk to the right. I found no fixed-charge right-boundary counterexample. [Kormos–Takács](https://arxiv.org/abs/0712.1886), [Lencsés–Takács](https://arxiv.org/abs/1106.1902).

### Independent winding computation using the NLIE formula

In 0711.4101, the source explicitly defines the physical strip charge as Q=m-plus−m-minus. Its ultraviolet conformal dimension, after eliminating the momentum term, is

\[
 \Delta_0=\tfrac12[(\phi_+-\phi_-)/\sqrt\pi+mR]^2,\qquad
 R=\sqrt{4\pi}/\beta.
\]

The exact source line for its winding/counting relation (`bbs.tex:1087`) is:

```tex
m=\frac{1}{2}({\rm sign}(H_{+})+{\rm sign}(H_{+}))-1-S.
```

The repeated H-plus in that line is left as printed here. I only use the regime H-plus>0 and H-minus>0, so that possible typo has no effect on the following calculation. Their counting variable S is N_H−2N_S−M_C−2M_W step(p−1).

Independently, minimize the ultraviolet free-boson energy with endpoint difference

\[
 \Delta\phi=\phi_+-\phi_-+2\pi Q/\beta.
\]

The minimizer is linear, `E=(Delta phi)^2/(2L)`, and thus `Delta=EL/pi=(Delta phi)^2/(2pi)`. This is exactly the quoted NLIE dimension with m=Q. For a one-hole state without the other listed sources, S=1 and their relation gives **Q=m=−1**. A zero-hole state has m=0. Therefore, on these formulas, the elementary hole lowers rather than raises the increasing-field charge. My independent numerical substitution, beta=1.3 and phi-plus−phi-minus=0.4, gives:

```text
Q=-1: direct free-field dimension 3.127941304018718
Q= 0: direct free-field dimension 0.0254647908947033
Q=+1: direct free-field dimension 4.358710534787949
```

All agree with the NLIE zero-mode expression to 9e-16. The unequal endpoint values distinguish the two charged sectors. This is useful prior evidence that their hole/reflection-factor labels need translation; a printed P-plus is not by itself a counterexample to item 1. It is also compatible with the antisoliton assignment at a right boundary.

**Limit of this check:** I have not reconstructed the complete microscopic arrow-to-field map, analytically continued every boundary string, or independently solved the NLIE with charge-resolved poles. Consequently I do not certify that the whole two-boundary lattice dictionary has been repaired, nor that this report establishes historical priority. An exact lattice comparison with a fixed physical charge would be a valuable further check. No “first discovery” claim should be added on the strength of this search.

## 7. Downstream text and code

**Verdict: the changed mathematical formulas and item (h) are correct; one printed code description still has the wrong ordering of signs. The sweep found no remaining operative positive Dirichlet dictionary in the current book/lecture/code formulas.**

### Independent continuation and fusion computations

For the specified continuation beta_SG=ib, lambda=−2/B, the corrected Dirichlet relation gives

\[
 E_D=\frac B\pi\left[-\frac{\lambda+1}{2}ib\phi_D\right]
 =\frac{i}{\pi}(1-B/2)b\phi_D.
\]

This is the new sign at `part8-defects/38-defect-sine-gordon.qmd:270`. Choosing beta_SG=−ib reverses it. The reflection factor is even in E_D, so neither this agreement nor the sign of F constitutes an independent physical sign test. My numerical substitution at B=0.73 and phi_D=0.37 gives zero residual to displayed precision.

For classical fusion write a=beta_SG*phi_D and t=theta_d^cl. The classical potential coefficients give

\[
 C_0=\cosh(t+ia/2),\qquad C_1=\cosh(t-ia/2).
\]

The corrected complex UV–IR equation then yields

\[
 \cos\frac{\eta-i\vartheta}{\lambda+1}=C_1
 =\cos(-a/2-it),
\]

and the conjugate equation gives C0. On the chosen continuation branch,

\[
 \eta=-(\lambda+1)a/2,\qquad
 \vartheta=(\lambda+1)t.
\]

Thus the sign in `part8-defects/43-defects-and-boundaries.qmd:152` agrees with the Dirichlet seed used in item (3). Comparing vartheta=lambda*theta_d gives the proposed defect dictionary only when the extrapolation beyond the classical limit is assumed. The existing paragraph correctly says that classical agreement fixes only theta_d/theta_d^cl->1 and does not fix all subleading terms.

My independent cosine/cosh residuals at `(t,a)=(0.6,0.8),(-0.35,2.3)` are zero to floating precision. These calculations also verify the actual assertions of item (h) in `code/lectures/defects_l9_fusion.py:823–832`. The whole supplied lecture script subsequently passed all 10 checks, including (h).

The sinh-Gordon sign match can be checked directly too. Put eta=pi*E/B, vartheta=−i*pi*F/B. Since lambda+1=−(2−B)/B, the corrected complex equations become

\[
 C_1=\cos\frac{\pi(E-F)}{2-B},\qquad
 C_0=\cos\frac{\pi(E+F)}{2-B}.
\]

Choosing E=(a0+a1)(1−B/2), F=(a0−a1)(1−B/2) indeed yields Cj=cos(pi*aj). Flipping F interchanges the nodes. This is algebraic consistency, not a sign-sensitive measurement.

### Search scope and actual leftover

I searched the snapshot's book chapters, appendices, lecture projects and code for the requested reference ids, mu variants, Dirichlet 4pi/beta formulas, psi_0 and phi_D, together with node-crossing phrases. Historical briefs were not treated as active formulas needing replacement. Broad search outputs are saved in `downstream_scan.txt` and `sign_scan.txt`; the latter also includes historical material and must be filtered accordingly.

The operative equations in the changed three lecture files have the negative sign and the mu-minus/mu-plus pairing. Positive eta used as the spectrum parameter, quoted positive formulas attributed to the sources, and phase-even sinh-Gordon tests are not stale Dirichlet dictionaries. No additional operative old-sign formula was identified in this sweep. Text search is not a proof that every equivalent spelling has been found.

There is nevertheless a concrete misleading output at `code/lectures/defects_l9_fusion.py:688`: item (e) prints `C0, C1 = cosh(theta_d -+ i psi0/2)`. Its own header, potential, and item (h) require `+-`, not `-+`. This is an existing description error exposed by the audit, not a failure of the newly changed item (h).

**Replacement at line 688:**

```python
report("(e) freezing is fusing: B(Phi, psi0) is the GZ potential with C0, C1 = cosh(theta_d +- i psi0/2), every "
```

Leave the subsequent arguments/continuation of the existing report call intact.

### Supplied-code comparison, separate from independent evidence

All executions used the ed0585a snapshot, with bytecode generation disabled. Results:

| Script | PASS lines | FAIL lines |
|---|---:|---:|
| `code/bsg_dirichlet_sign.py` | 5 | 0 |
| `code/boundary_sine_gordon.py` | 32 | 0 |
| `code/boundary_coideal.py` | 59 | 0 |
| `code/nonlocal_charge_chirality.py` | 17 | 0 |
| `code/shg_boundary_breathers.py` | 9 | 0 |
| `code/lectures/defects_l9_fusion.py` | 10 | 0 |

All processes exited successfully. Their logs are the six `book_*.out` files in `briefs/referee_bsg_dirichlet_sign_scripts/`. These are secondary comparisons; the evidence for this report is the preceding independent derivations and scripts.

## 8. Labels and citations

**Verdict: correct to leave the existing conjectural/numerical labels unchanged; wrong to leave Appendix E entirely unchanged under the repository's stated rule.**

### Independent audit

Inspection of the commit diff shows no promotion to [Theorem]. The changed first-order statement remains [Conjecture], the UV–IR relation remains [Conjecture], and the fusion claims retain their existing mixed status. Identities obtained by inserting one conjectured relation into another do not justify upgrades. The new sign argument should not be presented as an all-coupling theorem merely because its classical limit can be checked exactly.

The source checks support the citations in the changed paragraphs: GZ has the positive printed relation and calls its linearity conjectural; BPTT invokes that limit; BPT gives the quoted normal-ordering normalization; Delius–MacKay has the stated exponent/coupling pairing after the dictionary in item 3. The SSW citation must retain its opposite-half-line qualification.

The sinh-Gordon citations were checked against their own TeX. [Corrigan–Taormina](https://arxiv.org/abs/hep-th/0008237), source `ct.tex:268`, prints:

```tex
E=(a_0+a_1)(1-B/2),\qquad F=(a_0-a_1)(1-B/2).
```

[Chenaghlou–Corrigan](https://arxiv.org/abs/hep-th/0002065), `acec.tex:1138`, prints the same two expressions with `\qquad` between them and no intervening comma. These citations support the quoted sign convention, but their reflection-factor tests cannot distinguish F from −F.

I independently compared the newly added correction/phase-check claims with the entire half-line list in Appendix E: there is no corresponding entry. The applicable AGENTS.md explicitly requires updating that list “when a correction to a paper is added or dropped.” Appendix E's introductory sentence says corrections are listed separately in Appendix D; that explains the omission but does not override the explicit repository instruction. The phase check and charge-resolved dictionary correction should therefore receive a conservative entry, without claiming established priority.

The reference checker was run against the pinned snapshot after including the configuration files and code README files needed for a faithful check. It passed all six failure-producing categories, including the 42 existing Appendix E items. Its sole warning is the pre-existing title-style reference `citation.qmd:11 @sec-about`. An initial incomplete snapshot lacked four code README targets; after extracting those unmodified files from ed0585a, those artificial failures disappeared. The checker cannot detect a missing Appendix E item, as the repository rules correctly warn.

### Replacement text

Insert under the half-line heading at `appendices/E-original-results.qmd:30`, inside the existing original-results div:

> - **[Numerical]** A charge-resolved check of the sine-Gordon Dirichlet parameter identifies the negative sign on $x\le0$ with the soliton of positive topological charge. The classical scattering phase supports this correction to the sign printed by Ghoshal and Zamolodchikov. With this sign the conjectured boundary UV–IR relation pairs each coideal node with its own boundary coupling (@sec-bsg-parameters).

This entry reports what has been checked and uses the section id, as permitted for entries without a single labelled statement. It does not promote the exact quantum relation to a theorem. If the editor instead decides that Appendix E should categorically exclude corrections, that repository rule and the appendix policy need to be reconciled explicitly; this brief contains no exception authorizing silent omission.

## Claims not fully checked

- No proof of the exact finite-coupling Dirichlet linearity or the full two-parameter UV–IR conjecture is supplied or claimed here. The charge argument fixes the channel and edge value; the exact scattering calculation fixes the leading semiclassical slope.
- The prior-art search does not establish novelty. I have not reconstructed the full lattice arrow/field dictionary or solved the NLIE independently in charged boundary-string sectors. Item 6 records the concrete ultraviolet winding evidence and its limit.
- I did not rerun the original TCSA implementations. Blindness of the cited energy tests follows from the stated setups and the exact charge-conjugation similarity argument, not from a fresh TCSA dataset.
- I did not audit every signed boundary form factor or vacuum expectation value in the later literature. In particular, the Lencsés–Takács paper must not be described as entirely modulus-only.
- Higher-rank claims in `prp-dm-first-order`, the proposed nonanalytic Neumann offset, and the full fusion construction beyond the sign-dependent part of (4) were not re-refereed. They are outside this sign check.
- The scan did not establish absence of every possible algebraically disguised stale sign. The concrete leftover found is listed in item 7.

## Other observations — unverified

- It may be useful to trace whether the minus sign has previously been absorbed into endpoint orientation, the lattice-hole charge label, or the definition of phi_0 in less commonly cited treatments. This report does not claim such a source exists.
- A microscopic lattice calculation with a fixed increasing-field charge and unequal endpoint values could settle the remaining convention translation more cleanly than additional unlabelled energy tests. No such new calculation was performed here.
