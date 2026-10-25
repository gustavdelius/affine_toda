# Referee report: Lecture 9 additions to the defect chapters

Date: 2026-10-25. Referee: Codex, GPT-6 (the runtime does not expose a more specific model identifier).

Reviewed commit: `c444a86ece9c0b19e84791f33544901c9b669249`, in `.claude/worktrees/jovial-antonelli-69b2cf`. HEAD when reviewed was `3fc403d6702f6f1a9a7a28ebcb4301884b492d7f`; its changes after the reviewed commit concern briefs only. Line numbers below refer to the reviewed book files. No book file was changed and no subagent was used.

The two main formulas survive independent calculation. The scalar limit needs an explicit analytic step before the Riemann–Lebesgue lemma applies. The new chirality reading in chapter 43 is correct, but chapters 21 and 41 contain the opposite assignment. The novelty language needs narrowing: Bajnok–Simon explicitly say that they checked fusion with solitonic transmission factors. They do not print the integral proof under review.

Independent scripts and complete outputs are in `referee_scratch/defects_l9/`:

- `independent_symbolic.py`, `symbolic_output.txt`: the fused spin-3 charge, endpoint flux, sector algebra, folding identities and parameter conversions.
- `independent_followup.py`, `followup_output.txt`: all four matrix entries checked separately, a direct boundary primitive, and the Tzitzéica identity. This imports only my own first script.
- `independent_scalars.py`, `scalar_output.txt`: separately transcribed Gamma products, integral evaluation and cutoff extrapolation; no book code imported.
- `book_comparison_output.txt`: the author's nine checks also pass. This is corroboration, not my evidence. I did not rerun the older window calculation in `defect_doubling.py`.

## 1. Scalar identity and defect–Dirichlet fusion

**Verdict: correct identity and sector equivalence; the printed proof of the Dirichlet limit is incomplete as stated.** Nothing in the equality of the specified matrices needs to remain merely numerical after the short repair below. This is a mathematical identity for the specified scalar factors, not a proof that the conjectured UV–IR relation or a particular quantum defect action is exact.

### Independent derivation of the integral

Take real `lambda > 0`, put `w=lambda u/pi`, `X=xt/pi`, and let the four positive-real-part Gamma base arguments be

\[
a_1=\tfrac12+x/\pi+(2l-1)\lambda,\quad
a_2=\tfrac12-x/\pi+(2l-1)\lambda,\quad
a_3=\tfrac12-x/\pi+(2l-2)\lambda,\quad
a_4=\tfrac12+x/\pi+2l\lambda.
\]

Apply Malmstén separately and then form the signed combination `(+,+,-,-)`. The linear terms cancel. Directly summing the exponentials gives

\[
\sum_{r=1}^4\epsilon_r e^{-a_rt}
=e^{-t/2-(2l-1)\lambda t}
 (e^{-X}+e^X-e^{X+\lambda t}-e^{-X-\lambda t})
=-4e^{-t/2-(2l-1)\lambda t}
 \sinh(\lambda t/2)\sinh(X+\lambda t/2).
\]

This verifies the sign and the shift in the book. The grouped `l`-integrand is bounded near zero by `C t exp(-(2l-1)lambda t)`. Summing these bounds is `O(1)`, hence integrable at zero. At infinity the sum of grouped absolute values is bounded by

\[
\frac C t\exp\{- [\tfrac12-|\operatorname{Re}x|/\pi
                          -|\operatorname{Re}(\lambda u)|/\pi]t\}.
\]

The stated strict strip therefore suffices for absolute convergence of the grouped sum-integral and for all the Gamma arguments to have positive real parts. This also gives local uniform convergence on compact subsets. It is a sufficient strip, not a claim about the maximal continuation domain. Summing the geometric progression gives exactly

\[
\log\sigma(x,u)=\log\frac{\cos x}{\cos(x+\lambda u)}
-\int_0^\infty\frac{\sinh(\lambda ut/\pi)
 \sinh(xt/\pi+\lambda t/2)}{t\sinh(t/2)\cosh(\lambda t/2)}dt.
\]

Use the analytic logarithm normalized at `u=0`; thereafter the identity is meromorphically continued. No ungrouped sum of individual Malmstén integrands is needed.

### The Dirichlet limit

The result `sigma(i vartheta,-i theta) -> 1` is correct. The raw integrand in the printed argument does **not** tend to zero by the Riemann–Lebesgue lemma. Write `d=vartheta/lambda`, `a=pi/2`, `b=pi/(2lambda)`. The real oscillatory term has amplitude

\[
\frac{\sin(k\theta)}{k\sinh(bk)}\sim\frac{\theta}{bk},
\]

which is not in `L^1` at zero. Its contribution tends to a nonzero constant that cancels the elementary logarithm. The exact cosine-difference integral instead gives

\[
\log\sigma(i\lambda d,-i\theta)
=\log\cosh(\lambda d)
-\tfrac12\log\cosh\lambda(\theta-d)
-\tfrac12\log\cosh\lambda(\theta+d)
+i\int_0^\infty
 \frac{\sin(k\theta)\tanh(\pi k/2)}{k\sinh(\pi k/(2\lambda))}
 \cos(kd)\,dk. \tag{R1}
\]

For real `theta`, the remaining amplitude is bounded at zero and exponentially decreasing at infinity. Now Riemann–Lebesgue applies, and the three elementary terms tend to zero. The same argument works locally uniformly for complex `theta` in `|Im theta|<pi/(2lambda)`; continuation covers the meromorphic identity elsewhere. Thus the limit used for `K_D` is justified. The normalized diagonal limits are `cos(lambda u +/- eta_D)` and the off-diagonal entry vanishes, as claimed.

**Replacement, chapter 43 line 151:** replace “by (-@eq-dbd-sigma-integral) below and the Riemann–Lebesgue lemma” with “by the limit derived below”. Insert after `eq-dbd-sigma-integral`:

> To take the Dirichlet limit, set $d=\vartheta/\lambda$ and first take $\theta$ real. The cosine-difference integral gives
> $$
> \log\sigma(i\lambda d,-i\theta)=\log\cosh\lambda d
> -\frac12\log\cosh\lambda(\theta-d)-\frac12\log\cosh\lambda(\theta+d)
> +i\int_0^\infty\frac{\sin(k\theta)\tanh(\pi k/2)\cos(kd)}{k\sinh(\pi k/(2\lambda))}\,dk.
> $$
> The amplitude of $\cos(kd)$ is integrable: it is bounded at zero and decreases exponentially at infinity. The Riemann–Lebesgue lemma makes this integral vanish as $d\to+\infty$. The elementary terms also tend to zero. Thus $\sigma(i\vartheta,-i\theta)\to1$.

### Sector algebra, scalar matching and offset

Direct multiplication gives, with `Y=y/yd`, `Y'=1/(y yd)`,

\[
A(Y')^{-1}\operatorname{diag}(k_+,k_-)A(Y)
=\frac{4/y_d}{1+y^{-2}y_d^{-2}}\Lambda^{-1}
 \begin{pmatrix}P_0^+&Q_0\\Q_0&P_0^-\end{pmatrix}\Lambda,
\quad\Lambda=\operatorname{diag}(1,-e^{i\eta_D}).
\]

My symbolic computation gives zero for each of the four entries separately. In particular,
`k_+ + yd^-2 k_- = 4 P_0^+/yd`, with the analogous minus equation, and
`i(y^-2-y^2)/yd = -4 Q_0/yd`. Including the sector matrix `D_n` gives the conjugation `(Lambda D_n)^-1 K_GZ (Lambda D_n)`. The factor `q^-s/2` cancels between a transmission matrix and its inverse. The remaining `n,s` occur only in this constant similarity. No finite charge window is needed for this argument.

For completeness, subtracting the independently derived logarithm of the Phi ratio from the sigma integral leaves

\[
\lambda\theta+\log\cosh\lambda\theta_d
-\log\cosh\lambda(\theta+\theta_d)
=\log\frac{1+y_d^{-2}}{1+y^{-2}y_d^{-2}}.
\]

The difference of kernels is
`[cos k(theta+theta_d)-cos k(theta-theta_d)]/[k sinh(pi k/(2lambda))]`.
At `theta=0` both sides of the proposed scalar identity equal 1, fixing its multiplicative constant. Continuation then proves `eq-dbd-scalar-identity`; there is no residual CDD factor for these specified inputs. This does not prohibit multiplying an input transmission amplitude by an additional admissible CDD factor and obtaining a different dressed boundary.

The independent Gamma calculation used 38-digit arithmetic, common cutoffs, and degree-four extrapolation in `1/N`, comparing cutoffs `32,...,512` with `64,...,1024`. At `lambda=sqrt(2), e/2, sqrt(5)`, `x=.21+.31i`, `u=.13-.19i`, the sigma integral/product errors are respectively `5.87e-16`, `5.84e-16`, `6.42e-16`; cutoff changes are below `2e-14`. At `theta_d=.61` and `theta=.43, .37+.16i`, all six direct Gamma-product fusion residuals are below `4.7e-16`. These are measured convergence checks, not certified error bounds. For `lambda=sqrt(2), theta=.43`, the errors `|sigma(i vartheta,-i theta)-1|` at `vartheta=2,4,8,12` are `.254166`, `.0591254`, `.00345772`, `.000204349`.

The numerical paragraph at line 179 is compatible with the algebra. A sign change of `vartheta` generically changes the diagonal ratio; it is not an innocuous diagonal gauge. Special parameter choices can make a negative control degenerate, so “fail” there should be understood as referring to the sampled generic cases. I independently confirm the identities, not every historical residual or the claimed historical order of discovery. The Notes sentence at line 272 needs the prior-art qualification in item 6.

**Label:** retain `[Theorem]` for (3) after making the Dirichlet-limit argument explicit. Item (4) remains conditional on the UV–IR hypothesis.

## 2. The fused spin-3 charge

**Verdict: correct**, for the usual nonzero `sigma`, `m_0`, `beta` and smooth solutions for which the charge and endpoint fluxes exist. “For all sigma” must not literally include zero, since the potential itself contains `sigma^-1`.

Here is a direct check on the fused boundary, independent of a shrinking-interval argument. Put `f=phi(0,t)`, `v=f_t`, `a=f_tt`, `m=m_0`,

\[
c_\pm=\cos\tfrac\beta2(f\pm\phi_D),\quad
s_\pm=\sin\tfrac\beta2(f\pm\phi_D),\quad
F(f)=-\frac m\beta(\sigma s_++\sigma^{-1}s_-),\quad
W(f,v)=v+\frac m\beta(\sigma s_+-\sigma^{-1}s_-).
\]

Then `phi_x=F`, `psi_x=W`, and the field equation gives

\[
\partial_\pm\phi=v\pm F,\qquad
\partial_\pm^2\phi=2a\pm2F'v+(m^2/\beta)\sin\beta f.
\]

Insert these into the two printed currents and into the printed defect terms with `partial_+ psi=W`, `partial_- psi=-W`. The script treats `exp(i beta f/2)` and `exp(i beta phi_D/2)` as independent nonzero variables, so the calculation is a rational identity for symbolic `beta,m,sigma,phi_D`, not a sampling of angles. It obtains

\[
\frac d{dt}\left(\mathcal D_{-3}+\mathcal D_{+3}
-\frac\beta m\sin\beta\phi_D\,W\right)=j_{-3}(\phi)+j_{+3}(\phi).
\]

Since the time derivative of the half-line bulk integral is minus this current, the charge is conserved. The output is `PASS fused spin3, symbolic beta,m,sigma,phi_D: 0`.

There is also a concise direct proof suitable for the book. Write the displayed boundary term as `H(f,v)`, and define `U(f)=H(f,0)`. Independent simplification gives

\[
H=-\frac{\beta^2}{m^3}F'v^2+U(f),\qquad F''=-\frac{\beta^2}{4}F,
\]
\[
U'=\frac{\beta^2}{2m}\cos(\beta f)F
-\frac\beta m\sin(\beta f)F'
+\frac{\beta^4}{4m^3}F^3.
\]

Expanding the two currents gives exactly the derivative of this expression. The follow-up script checks both identities. This removes any need to assume a convergent family of interval solutions as `a -> 0`.

If the Dirichlet term is omitted, the residual is

\[
\frac\beta m\sin(\beta\phi_D)
\left[a+\frac m2(\sigma c_+-\sigma^{-1}c_-)v\right].
\]

It is generically nonzero; it vanishes at special Dirichlet values `sin(beta phi_D)=0`, so the negative control is not a universal failure statement. Both the coefficient and the minus sign of the correction are necessary.

At the far Dirichlet endpoint, `psi_t=psi_tt=0` and
`partial_+^2 psi=M+2 psi_xt`, `partial_-^2 psi=M-2 psi_xt`, where `M=m^2 sin(beta phi_D)/beta`. Thus their difference of squares is `-8 M psi_xt`, and

\[
j_{-3}+j_{+3}=-\frac\beta m\sin(\beta\phi_D)\psi_{xt}.
\]

The two-half-line-plus-interval bulk derivative has the additional contribution **minus** this endpoint current. Subtracting `(beta/m) sin(beta phi_D) psi_x` therefore gives the sign printed in the theorem. The second defect equation indeed gives `W`; setting `psi_x=0` would be wrong.

For a static elastic reflecting boundary, a combination with one-particle eigenvalue `A exp(3 theta)+B exp(-3 theta)` can survive reflection at all rapidities only if `A=B`. This verifies the sum-versus-difference statement in its intended setting. It is not a classification of every possible dynamical or species-dependent boundary charge.

Ghoshal–Zamolodchikov's [Appendix](https://arxiv.org/abs/hep-th/9306002), equations (A.3)–(A.9), computes this spin-3 boundary-current condition. Their exact line after (A.7) is: “This reduces to a total d/dy derivative if the function V (φ) satisfies”. They then give `4 V''+V=0` and `V=-Lambda cos((phi-phi_0)/2)`. Independently, the fused potential is a linear combination of `cos(beta f/2)` and `sin(beta f/2)` and obeys precisely this differential equation after rescaling the field and dropping its constant. Thus conservation for this boundary family is established prior art; the particular decomposition into two defect contributions and a Dirichlet correction is the narrower contribution here.

**Suggested precision edits, chapter 38:** line 375, say “for all nonzero $\sigma$ and all $\phi_D$”; line 385, begin “A static elastic reflecting boundary maps $\theta$ to $-\theta$. Among these two equally normalized local charges it can preserve only their sum, with a boundary correction.” The theorem label is justified by the exact direct calculation above. The interval discussion is a useful explanation but should either assume a smooth limit or include the direct proof.

## 3. The boost heuristic

**Verdict: correct as a heuristic with the assumptions explicitly printed.** It does not establish a renormalization prescription or the additive rapidity constant.

The Gaussian two-point function in this normalization is `-log(z-w)/(4pi)` for one chiral field. Wick contraction gives weight `a^2/(8pi)` for `:exp(i a phi_L):`. Setting `a=beta_SG` gives
`h=beta_SG^2/(8pi)=1/(lambda+1)`.
Solving the massless type-I-minus gluing gives `psi_L=phi_L+constant`, `psi_R=-phi_R+constant`; hence the two cosines reduce to the two stated chiral fields, up to constant phases. The factor of two in `phi +/- psi` is essential and is correctly included.

For a chiral coordinate rescaling `z'=c z`, the integrated perturbation transforms as
`g integral dz O(z) = g c^(h-1) integral dz' O'(z')`.
Assigning `c=exp(+/- zeta)` to the two chiralities gives the stated opposite factors. If the bare parameter denotes precisely the ratio of these normal-ordered couplings, and the scattering parameter shifts by `zeta`, this gives
`theta_d^cl=(1-h)theta_d+constant`.
Algebraically, `1/(1-h)=(lambda+1)/lambda`; setting `lambda=-2/B` gives `1-B/2`. These identities pass the independent symbolic check. Relating this coordinate argument to the family of static massive defects is part of the heuristic assumption, not an additional theorem about boosts of a fixed defect worldline.

Konik–LeClair, [hep-th/9703085](https://arxiv.org/abs/hep-th/9703085), section 5.1, p. 11, states: “The scaling dimension of g is d(g) = 1 − β̂²/2”. Their immediately following definition is `beta_hat=beta/sqrt(4pi)`, and their action (5.1) has kinetic coefficient `1/2`. Consequently their dimension is exactly `1-h` here.

Bajnok–Holló–Watts, [1307.4536](https://arxiv.org/abs/1307.4536), equations (101)–(102), give

\[
\hat\mu=\frac{|h_c|m^{6/5}}2e^{-i(b+3)\pi/5},\qquad
\hat{\bar\mu}=\frac{|h_c|m^{6/5}}2e^{+i(b+3)\pi/5}.
\]

Their exact introduction to these equations is “Using this, we see that the correct identification is”. Equation (97) is the boundary UV–IR relation, (98) defines `h_c`, (99) gives the two possible sign assignments, (100) checks the product against the bulk integrability condition, and the intervening text fixes the signs using `b=-3+i alpha`, `alpha -> +/- infinity`. Equation (96) belongs to the preceding energy calculation, not the coupling identification. Thus the book's range (97)–(102) is right.

Their defect-field decomposition (22), together with the highest weight `h=-1/5`, supplies fields of weights `(-1/5,0)` and `(0,-1/5)`. Independently, substituting `b=3+6i theta_b/pi` gives
`-i(b+3)pi/5=(6/5)(theta_b-i pi)` and its negative for the other coupling. The exponent `6/5=1-h` and the signs in the bullet are correct. The final two sentences explicitly state the needed assumption and disclaim a proof. No substantive change recommended.

## 4. Folding and linear waves

**Verdict: correct formulas and attributions.** Distinguish the criterion for the linear scattering map to preserve the fold from the full nonlinear folding construction.

Let `rho` be the diagram reflection and `pi` the cyclic rotation. Since `rho pi rho=pi^-1`, a `pi`-eigenvector of eigenvalue `omega^k` is taken to one of eigenvalue `omega^-k`. On these modes `D=2/(1-omega^k)`; replacing `D` by `D^T` sends this to `2/(1-omega^-k)=2-D`. Thus reflection interchanges `k` and `h-k` in the linearized defect equations, including their potential Hessians, and the mirror factor on mode `k` is `T_(h-k)`, with the same parameter. This is an internal Dynkin reflection, not a reversal of the direction of incidence.

Writing `r=exp(i gamma)`, direct rational algebra gives

\[
T_k=\frac{1+xr}{1+x/r},\quad T_{h-k}=\frac{1-x/r}{1-xr},\quad
\frac{T_k}{T_{h-k}}=\frac{1-x^2r^2}{1-x^2/r^2}.
\]

Cross-multiplying at `x_1` and `x_2` gives
`(x_1^2-x_2^2)(r^-2-r^2)=0`. For `1<=k<=h-1`, the coefficient vanishes only at `2k=h`. Otherwise equality for the two reflected modes requires and suffices for `sigma_1^2=sigma_2^2`. The self-conjugate mode imposes no restriction by itself. The independent script verifies both transparency at opposite signs and

\[
T_kT_{h-k}=\frac{\sinh u-i\sin\gamma}{\sinh u+i\sin\gamma}
=f_{k/h}(u)^{-1}.
\]

Replacing `exp u` by `-i exp u` gives `(cosh u+sin gamma)/(cosh u-sin gamma)`. For positive sigma and real `u` this factor is positive, so it causes no change of topological charge. For a non-self-conjugate mode it is finite on the real axis. The exceptional self-conjugate sine-Gordon mode has the familiar singularity at `u=0` and should not be swept into a blanket claim of finite delay.

At `h=3`, the exact reduction is

\[
\frac{\cosh u+\sqrt3/2}{\cosh u-\sqrt3/2}
=\frac{(e^{-u}+e^{i\pi/6})(e^{-u}+e^{-i\pi/6})}
       {(e^{-u}-e^{i\pi/6})(e^{-u}-e^{-i\pi/6})}
=\coth(u/2-i\pi/12)\coth(u/2+i\pi/12).
\]

This is Corrigan–Zambon's `z_1` in (4.9)–(4.10), in the sector `exp(hat lambda)->1`: their (4.12) tends to 1 as the incoming exponential tends to zero. The other two cube-root sectors do not have this delay. The book already prints `exp(theta_d)=sqrt(2) hat sigma` in the preceding section; matching the delay as a function of relative rapidity does not by itself prove equality of the full folded actions.

Robertson, [1304.3129](https://arxiv.org/abs/1304.3129), section 3.1, derives `d_tilde=+/-d` in (3.16) and again from the soliton delays on pp. 14–15. His exact concluding phrase is “there is a non-trivial time delay or advance and no change of topological charge.” His opposite-sign case has `z_p z_tilde_p=1` and vanishing defect Lagrangian when the outer fields agree. Corrigan–Zambon, [0908.3126](https://arxiv.org/abs/0908.3126), section 5, explicitly send `sigma -> 1/(sqrt(2) sigma)` before (5.8). Thus “rescale” is a loose description of a substitution that includes inversion, but the attribution is right.

**Suggested replacement, chapter 39 line 422:** replace “so the pair folds exactly when” with “so its linear scattering preserves the folding constraint exactly when”. At line 431, optionally specify “the sector $e^{\hat\lambda}\to1$” after the identification with Tzitzéica `z_1`. These sharpen the scope without changing the result.

## 5. Chirality and the Borel halves

**Verdict: the new chapter 43 text is correct; the cited older assignments in chapters 21 and 41 are wrong in the book's stated spin convention.** This is a convention inconsistency in the exposition, not evidence against the computed intertwiners.

Independently, for `phi=L(t+x)+R(t-x)` and `P=-integral phi_t phi_x`,

\[
E+P=2\int(R')^2dx,\qquad E-P=2\int(L')^2dx.
\]

The left-moving sector therefore has negative Lorentz spin. A weight-`h_J` integrated left-moving current has spin `-(h_J-1)`, whereas its right-moving counterpart has spin `+(h_J-1)`. On the evaluation module `e_j` scales as `exp(lambda theta)` and `f_j` as `exp(-lambda theta)`. Consequently `e_j` belongs to the positive-spin/right-moving half and `f_j` to the negative-spin/left-moving half.

Solving `phi_x=psi_t`, `psi_x=phi_t` gives `psi_L'=phi_L'`, `psi_R'=-phi_R'`. Changing the kinetic sign gives `psi_L'=-phi_L'`, `psi_R'=phi_R'`. These are precisely the type-I-minus and type-I-plus assignments of chapter 43. The Weyl action `rho(f_j)=mu yd^-1 times shift` scales by `exp(-lambda zeta)`. The two coefficients of the second module in `eq-dbd-second-module` are proportional to `yd` and scale by `exp(+lambda zeta)`. Thus the boost statement is correct too. “Dominate” refers to the rapidity dependence; a rigorous massless limit also rescales the dimensionful charges.

### A local repair that preserves the physical q and the evaluation matrices

It is insufficient just to relabel the word “left” in chapter 21 while leaving its vertex exponent and continuity equation unchanged. With the canonical bracket `[phi(x),pi(y)]=i delta(x-y)` and a common left-endpoint convention, set

\[
\phi_L=\tfrac12\left(\phi+\int_{-\infty}^x\pi\,dy\right),\qquad
\phi_R=\tfrac12\left(\phi-\int_{-\infty}^x\pi\,dy\right).
\]

Then

\[
[\phi_L(x),\phi_L(y)]=-\tfrac i4\operatorname{sgn}(x-y),\qquad
[\phi_R(x),\phi_R(y)]=+\tfrac i4\operatorname{sgn}(x-y),
\]
\[
[T,\phi_L]=-i\beta/(4\pi),\qquad [T,\phi_R]=+i\beta/(4\pi),
\quad T=\frac\beta{2\pi}\int\phi_xdx.
\]

Thus a consistent leading-current convention retaining `[T,Q_+]=2Q_+`, `Q_+ -> e_1` is

\[
J_\pm=:e^{\mp ia\phi_R}:,\qquad
\bar J_\pm=:e^{\pm ia\phi_L}:,\qquad a=8\pi/\beta.
\tag{R2}
\]

Both the chirality and the sign of the unbarred exponent change; the barred labels here also have topological charge `+/-2`, so `bar Q_-` can represent `f_1` as the book requires. The continuity equation becomes `partial_+ J=partial_- H`. In the sentence about the double pole, the unbarred `J_+` now pairs with `exp(+i beta phi)`, not `exp(-i beta phi)`.

This also fixes the braiding explanation rather than changing the physical q. For right currents the equal-time exchange factor of two same-sign exponentials is `exp(-i a^2 sgn(x-y)/4)`, the inverse of the factor for left currents. More directly, moving a charge density past a topological step `Delta phi=2pi T/beta` gives the right current in (R2) the phase

\[
\exp(-ia\,\Delta\phi/2)=\exp(-8\pi^2 i T/\beta^2)=q^T.
\]

The original `exp(+ia phi_L)` instead gives the inverse phase with the same ordering and endpoint convention. Changing chirality without tracking contour/order and exponent signs would therefore confuse `q` with `q^-1`. Nothing here gives permission to replace `q` by `-q`. With (R2) and the spatial ordering used for the book's coproduct, `q=exp(-i pi h_J)=-exp(-i pi lambda)` stays exactly as printed. This is an independent canonical derivation; it is not a claim that all source papers use (R2).

For a source check, Delius–MacKay [hep-th/0112023](https://arxiv.org/abs/hep-th/0112023), footnote 2 and equations (3.5), (3.10), explicitly specify their Euclidean coordinates and their q. Their footnote says “the Euclidean time t is related to Minkowski time tM by t = itM.” Their q is inverse to the book's value. This illustrates why copying a holomorphic-current label without transporting the ordering convention is unsafe; it does not invalidate the physical q selected by the book's scattering conventions.

**Replacement text and locations:**

- Chapter 21 lines 16–25: replace the leading-current formulas and conservation equation by (R2) and `partial_+J_\pm=partial_-H_\pm`, and change the double-pole partner to `e^{\pm i beta_SG phi}`. Introduce them as “In our rapidity, topological-charge and coproduct conventions, choose the leading currents ...”, rather than attributing these precise signs verbatim to Bernard–LeClair.
- Chapter 21 line 34: “The unbarred charges are built from right-moving currents and have spin $h_J-1$; the barred charges are built from left-moving currents and have spin $-(h_J-1)$. Both chiral fields are non-local functionals of $\phi$ and its conjugate momentum.”
- Chapter 21 lines 42, 46, 60 and 66: the chiral field used to construct the **unbarred** charges/e-generators must be `phi_R`; keep `phi_L(x^+)` as the actual left-moving field in the decomposition. Insert the phase derivation above near line 52 to fix the orientation explicitly. Keep `eq-q-physical`, the coproduct and the evaluation representations unchanged.
- Chapter 41 line 61, first sentence: “The non-local charges $Q_j$, $\bar Q_j$ of @sec-nonlocal-charges are built from right- and left-moving vertex operators, respectively, and $Q_j\to e_j$, $\bar Q_jq^{-T_j}\to f_j$, $q^{T_j}\to k_j$.” Replace the Borel sentence by: “These generate $U_q(\mathfrak b^-)=\langle f_j,k_j^{\pm1}\rangle$ when the perturbation is right-moving, and $U_q(\mathfrak b^+)=\langle e_j,k_j^{\pm1}\rangle$ when it is left-moving.” This describes the untouched opposite chirality; it does not derive the full two-perturbation defect from first-order conformal perturbation theory.

No reversal of the new chapter 43 paragraph is recommended.

## 6. Prior art, novelty and citations

**Verdict: the Section 3 citation correction is correct; novelty of the specific explicit derivations is unclear, and novelty of fusion or boundary spin-3 conservation as such is not justified.** The following distinction should be made in Appendix E and in the Notes.

I searched the named papers, the local primary-source copies, and online for defect/boundary fusion and work by Caudrelier, Avan, Doikou and Delius. The mathematical comparisons are as follows.

1. Bajnok–George, [hep-th/0404199](https://arxiv.org/abs/hep-th/0404199), equation (19), already gives the transmission–reflection–transmission dressing, with the reflection equation inherited from the transmitting and boundary equations. Its index contraction becomes the book's formula after translating right-incidence rapidities. Thus the fusion mechanism is not new.
2. Bajnok–Simon, [0712.4292](https://arxiv.org/abs/0712.4292), equations (10)–(17), explicitly match Dirichlet and general sinh-Gordon reflection factors and obtain the parameter renormalization. More significantly, their conclusion, p. 21, says: **“We have checked this solution by performing the fusing procedure on the solitonic transmission factors.”** This is a direct earlier claim of the sine-Gordon soliton-fusion check. Their PDF's extracted text drops some letter `c` glyphs; the quotation restores those glyphs, not the wording. I did not find the sigma/Phi identity or the Malmstén proof printed there. Accordingly, this source is evidence against claiming the whole soliton-fusion result or parameter map as a new discovery; it does not establish that the particular proof here was already published.
3. Corrigan–Zambon, [1202.6016](https://arxiv.org/abs/1202.6016), (2.16)–(2.18), compute the matrix coefficients. Their exact statement after (2.18) is “which is equivalent to the full Ghoshal-Zamolodchikov solution.” The displayed calculation suppresses an overall factor and does not supply this integral evaluation. Section **3**, not 2, contains the Habibullin discussion, including “using a combination of two carefully designed Bäcklund transformations.” The citation change is correct. Their (2.20)–(2.23) uses an energy-like spin-three condition for a dynamical/type-II boundary; that is not the explicit type-I-minus fused term printed here.
4. Ghoshal–Zamolodchikov's Appendix already proves the spin-three boundary condition for the potential to which the fused potential reduces. My computation in item 2 demonstrates the reduction and the equality of its flux primitive with the new assembled term. The existence of this conserved charge is therefore old. A presentation of its boundary correction as `D_-3+D_+3-(beta/m) sin(beta phi_D) psi_x` is a narrower result.
5. Caudrelier, [0704.2326](https://arxiv.org/abs/0704.2326), (3.28)–(3.30), gives an explicit third-order defect charge in light-cone coordinates. The text calls (3.28) “the third order conserved quantity”. It is a charge across an internal defect, not the Dirichlet-fused boundary expression above. [1411.5171](https://arxiv.org/abs/1411.5171) treats the sine-Gordon defect's multisymplectic and canonical structure; I did not locate this fused term there.
6. Avan–Doikou, [1205.1661](https://arxiv.org/abs/1205.1661), says in the charge construction that “Schwartz boundary conditions are imposed at the end point of the system ±L.” Its defect/sewing calculation is not the finite Dirichlet boundary in this brief. Doikou's [1404.7329](https://arxiv.org/abs/1404.7329) treats simultaneous boundaries and defects and explicitly says, “Based on this construction the vector NLS model was investigated.” Its reflected defect contributions supply relevant general methodology, not the sine-Gordon scalar-factor identity. These searches did not locate the precise two new formulas in those works.

This is a bounded literature search, not a proof of absence. In particular, “Corrigan and Zambon found the matrix part” should not be the only acknowledgement attached to the soliton-fusion result: it omits Bajnok–Simon's earlier explicit statement. Appendix E currently bundles established physical content with a possibly new explicit proof and decomposition.

**Replacement, Appendix E line 50:**

> - **[Theorem]** An explicit integral proof matches the scalar factors in the previously studied fusion of a type-I$_-$ sine-Gordon defect with a Dirichlet boundary, at $\vartheta=\lambda\theta_d$, with no additional CDD factor for the specified inputs. An explicit formula assembles the classical fused boundary's spin-3 correction from the two defect terms and a Dirichlet term (@prp-defect-fusion-k; @prp-sg-fused-spin-three).

**Replacement, chapter 43 line 181, the sentence beginning “The new points are”:**

> Bajnok and Simon also report checking fusion with the solitonic transmission factors [@bajnok_simon2008, Sec. 7]. Here the parameter map is written in our conventions and the scalar-factor identity is proved explicitly through integral representations. Part (4) states separately the UV–IR hypothesis needed to identify the quantum and classical defect parameters. We have not located this explicit scalar-factor proof in the sources cited here.

**Replacement within the long Notes sentence, chapter 43 line 272:** replace the segment beginning “the parameter map of @prp-defect-fusion-k” and ending at “(-@eq-dbd-sigma-integral)” with:

> the parameter dictionary and conditional UV–IR comparison of @prp-defect-fusion-k (3)–(4), following the soliton-fusion check reported by Bajnok and Simon [@bajnok_simon2008, Sec. 7], and the explicit integral proof of the scalar identity (-@eq-dbd-scalar-identity) through (-@eq-dbd-sigma-integral)

The Notes may say the derivation is “stated here”; they should not imply that neither fusion nor its solitonic check had appeared before. The Appendix E entry should be retained as an original **explicit proof/formula** only if that narrow scope is what the editor intends. The broader claim should not be retained as written.

## Other observations — unverified

- I did not audit all downstream uses of the chapter 21 current labels, including other algebras, boundary non-local charges, or exercise solutions. The local repair above is derived for sine-Gordon and fixes the displayed spin/coproduct mismatch; a repository-wide convention audit is separate work.
- I did not establish equality of the complete nonlinear folded actions in the Robertson and Corrigan–Zambon parameterizations. I verified the linear condition and the specified soliton-delay sector, as requested.
- I did not establish that no later paper or unpublished computation contains the same scalar-factor integral proof or the same assembled boundary term. The novelty recommendation is deliberately limited to the sources actually inspected.

## Integration

Integrated by Claude (Opus 5.5) on 2026-10-25, on branch `claude/defects-l9-book`, after the report was committed as received (759efaf). The referee's three scripts and their outputs are in `briefs/referee_defects_l9_report_scripts/`; rerun from there, all three reproduce the recorded outputs exactly. That rerun is corroboration only. Each item below was checked again by hand and with new checks in `code/defect_fused_boundary.py`, which uses the book's definitions and does not import the referee's scripts. The source lines were read in the arXiv texts.

| Item | Verdict | How it was checked | Files changed |
|---|---|---|---|
| 1. Dirichlet limit in the proof of `prp-defect-fusion-k` (3) | **Verified.** The raw integrand at $x=i\lambda\theta_d$ contains $-\sin k\theta\sin k\theta_d/(k\sinh\frac{\pi k}{2\lambda})$, whose amplitude grows like $1/k$ at $k=0$; that part tends to $-\lambda\theta$, not 0, and cancels the limit $\lambda\theta$ of the elementary logarithm. (R1) follows from the book's display by splitting $\sinh(\frac{\pi k}2+ik\theta_d)$ and using $\int_0^\infty(1-\cos kv)\frac{dk}{k\sinh(\pi k/2\lambda)}=\ln\cosh\lambda v$. **Modified** in presentation: the step is placed after the existing display at $x=i\lambda\theta_d$, not after `eq-dbd-sigma-integral`, and stated for real $\theta$ | By hand. Check 1 of the script, extended: (R1) against the Gamma product of $\sigma$ at $\lambda=\sqrt2$, $\theta=0.43$, $d=1,2,4,8$ (max rel. err. over check 1: $2.3\times10^{-25}$); $\lvert\ln\sigma\rvert=3\times10^{-4}$ at $d=8$; the raw $\sin k\theta\sin kd$ integral tends to $\lambda\theta$ | ch. 43 (proof of (3); the checking paragraph) |
| 1. Label of (3) | **Verified**, [Theorem] kept: the proof is now complete | — | — |
| 2. `prp-sg-fused-spin-three` | **Verified.** "For all $\sigma$" → $\sigma\ne0$. The direct proof is correct: the flux at $x=0$ is $-\frac{2\beta^2}{m_0^3}F'va+\frac{\beta^4}{4m_0^3}Fv^3+U'v$, which is $\frac{d}{dt}\big[-\frac{\beta^2}{m_0^3}F'v^2+U\big]$ exactly because $F''=-\frac{\beta^2}4F$. By hand, the $v^2$ terms of $\mathcal D_{-3}+\mathcal D_{+3}$ are $\frac{\beta^2}{2m_0^2}(\sigma c_++\sigma^{-1}c_-)v^2=-\frac{\beta^2}{m_0^3}F'v^2$, and their $v$-linear terms are $\frac\beta{m_0}\sin(\beta\phi_D)v$, which the Dirichlet term cancels. That is the reason for the Dirichlet term. **Modified:** the direct proof is now the proof; the interval argument is kept as the origin of the Dirichlet term. The sentence after the proof says "static boundary" and gives the $e^{\mp3\theta}$ reason | By hand (flux expansion, $v^2$ and $v$ coefficients). Check 6, extended with exact sympy tests of each step (flux formula, $F''$, $v^2$, $v^1$ and $v^0$ parts of the boundary term; the $v$-linear term of the defect terms alone). Controls for non-vacuity fail as expected | ch. 38 (statement, proof, following paragraph) |
| 2. Prior art: Ghoshal–Zamolodchikov App. A | **Verified.** Read in hep-th/9306002: after (A.7), "This reduces to a total d/dy derivative if the function V(φ) satisfies" (A.8) $4V''+V=0$, solved by (A.9). The fused potential is of this form ($F''=-\frac{\beta^2}4F$). Conservation is therefore theirs; the form of the boundary term is what the proposition adds | Source read; the reduction is the $F''$ identity above | ch. 38 (proof and following paragraph cite [@gz1994, App. A]); Appendix E |
| 3. Boost heuristic | **Verified**, no change, as the referee recommends | (Drafting checks, check 7) | — |
| 4. Folding | **Verified.** "folds exactly when" overstates a linear criterion; now "the linear scattering of the pair preserves the folding exactly when". For $2k=h$ the ratio is 1 (no condition) and the delay is $\coth^2\frac u2$, with a pole at $u=0$; the delay is finite and positive only for $2k\ne h$. The Tzitzéica $z_1$ is the sector $e^{\hat\lambda}\to1$, as the book's Tzitzéica paragraph already says. Corrigan–Zambon "send σ → 1/(√2σ)" before their (5.8) (read in 0908.3126): "rescale" replaced by the substitution | By hand; check 8 extended with $\coth^2\frac u2$ at $\gamma=\frac\pi2$ | ch. 39 |
| 5. Chirality, ch. 43 | **Verified**, no change | (Drafting check 9) | — |
| 5. Chirality, ch. 21 and ch. 41 | **Not integrated here.** These files belong to the session that fixes the chirality labels (task "Fix chirality labels of non-local charges (ch. 21, 41)"); the report's (R2) and its replacement text go to that session | — | — |
| 6. Bajnok–Simon | **Verified.** Read in 0712.4292, Sec. 7 (Conclusions), p. 21: "We have checked this solution by performing the fusing procedure on the solitonic transmission factors." "This solution" is their soliton transmission factors with the parameter relation continued from sinh-Gordon theory, so their check covers the matrix part and (4); no computation is printed. Ch. 43 and its Notes now credit it, and say that the explicit form of that check is what is stated here | Source read | ch. 43 (paragraph after the checks; Notes, two sentences); Appendix E |
| 6. Corrigan–Zambon 2012, Sec. 3 | **Verified** (already corrected in c444a86): the Habibullin remark is in Sec. 3, "Discussion", of 1202.6016 | Source read | — |
| 6. Appendix E item | **Verified; modified.** Narrowed to an explicit proof of the fusion with its scalar factor, crediting Corrigan–Zambon and Bajnok–Simon, and to the form of the boundary term of Ghoshal–Zamolodchikov's spin-3 charge. The referee's replacement text was adapted, not pasted | `check-refs.py` (labels match) | Appendix E |
| 6. Caudrelier, Avan–Doikou, Doikou | No change needed: the book makes no claim that these sources contradict | — | — |

**Not changed, with reasons.**

- Item 1: the remark that $\vartheta\to-\vartheta$ is not a diagonal gauge agrees with the text, which lists $\vartheta=-\lambda\theta_d$ as a failing control. The remark on extra CDD factors is covered by the statement, which fixes $\mathbb T$ with its scalar factor.
- Item 1: the referee's local uniformity for complex $\theta$ is not needed. The Dirichlet matrix is the limit at real $\theta$, continued.
- **Lecture notes.** The same issues in the lecture notes are handled on branch `claude/defects-lecture-9` (PR #56), by the session writing them: the reversed $e_j/f_j$ chirality in Lecture 7 (l. 43), the Corrigan–Zambon Sec. 3 citation, the Bajnok–Simon credit and the folding scope in Lecture 9, and the dropped sentence that the book states the scalar identity only numerically.

**Checks.** `code/defect_fused_boundary.py`: 9 PASS (223 s, with a render running alongside). `python3 check-refs.py`: all PASS (the WARN on `citation.qmd` is from before this branch). `./render-scratch.sh`: HTML renders (exit 0), no unresolved cross-references; the one TeX-math warning ($U_q(\widehat{sl}_2)$) is from before this branch.
