# Independent Codex referee report: Part III, phase 1

Date: 2026-10-15. Model: GPT-6 (Codex; no more specific model identifier exposed).
Reviewed commit: `bcd82b9ee9667975eed02109f08a257d50ba38b1`, branch `claude/referee-briefs-part-iii`, worktree `.claude/worktrees/referee-briefs-part-iii`.

The editor requested an independent parallel review with Codex-specific filenames. This overrides the brief's default report filename. Scratch work is confined to `referee_scratch/classical_09_11/codex/`; every created file has `codex` in its name. No book, author script, existing report, index or git state is edited. No subagents or figure scripts are run. Earlier reports, the reading brief, Gemini's work and the phase-2 list are not consulted. Chapter snapshots and SHA-256 hashes preserve the reviewed text.

## 1. Risk list (frozen before checking)

I read all three chapters, the part page, captions, exercises and notes, and the specified background before writing this list. Rank is fixed. Locations below refer to the reviewed worktree, not the later main checkout. Quoted snippets are from the book.

1. **Singular energy integrals.** 11:60,91, `sec-singular-solutions`, `exr-singular-points`: “double poles with zero residue, so integrated charges are unaffected.” Zero residue need not make a real integral converge. Check a tuned sine-Gordon tau zero, Laurent coefficients and cutoff integration.
2. **Energy formula at complex rapidity.** 11:32–39, `prp-hirota-energy`: “For a Hirota solution” and “(possibly complex) rapidities”. The boundary evaluation depends on which exponent grows at each spatial end and on regularity. Derive the boundary formula and test a rapidity outside the usual positive-real-decay strip; inspect Takács–Watts eq. 4.6.
3. **Scope of the selection rule and label.** 11:45–49 and Appendix E: “every conserved charge” and “[Theorem, given additivity of the higher charges]”. The data omit internal phases and topological charges, and symmetry alone is only sufficient. Derive the rule under precise hypotheses, test the kink's topological charge, and distinguish invariant fields from invariant reduced data.
4. **Breather existence and parameter band.** 10:130–136, 11:26–28,87: “suitably related” phases give “a localized, time-periodic solution of zero topological charge”; the figure uses all `u` in `(0,pi/2)`. Derive the two-tau zero criterion for conjugate species; compare Harder–Iskandar–McGhee and reconstruct sine-Gordon phases.
5. **Lax pair and twisted scope.** 09:57–77, `eq-lax-pair`, `exr-lax`: “For j=0 ... root vectors E^{-theta} ...” for every affine theory. Independently calculate curvature and boost covariance, then test the description of the affine generator in a twisted rank-one example.
6. **Fusion range and endpoint.** 10:126,163: “single soliton of species a+b”. Species only range from 1 to h-1; at a+b=h the asserted particle disappears and the pole changes order. Compute the limit for sums below, equal to and above h.
7. **Area rule and nonzero-coupling qualification.** 09:128–140, `eq-area-rule`: `|C_abc|=4 beta Delta_abc/sqrt(h)`. This equality is displayed without the selection rule even though many triangles have zero coupling; mass units also merit checking. Compute A2 cubic couplings in a charge basis, at general m, against Heron's formula.
8. **Time delay and semiclassical sign.** 10:122–124: “Delta t=-Delta x/v_1<0” from theta_1>theta_2; the claimed sine-Gordon phase derivative. Check asymptotic trajectories in boosted frames, the equal-rapidity endpoint of A, and the small-coupling limit of chapter 3's integral.
9. **Spin-3 current, charge and rapidity convention.** 09:90–96,167: the displayed current and “chi m^3 e^{3 theta}”. Light-cone signs and integration conventions can reverse the spin. Differentiate using the equation of motion and evaluate the current on a boosted sine-Gordon kink.
10. **Topological-charge cycle.** 10:82–92,155; 11:66: “h/gcd(a,h)” weights of `Lambda^a C^h`. Logarithm orientation may conjugate the representation or alter the count. Derive coordinates in the standard A root basis and enumerate phase chambers independently.
11. **Lie duality and folded masses.** 10:104,159: c_n soliton masses proportional to B_n particle masses. Folding may rescale roots, kinetic energy or m. Construct B_n and C_n mass matrices directly and compare the coincident-pair spectrum, including c_2's explicit polynomial.
12. **Hirota transpose and factors.** 10:33–59,143–147: the bilinear equations for all affine Cartan matrices and the moving one-soliton. Independently substitute the coroot ansatz; evaluate two-soliton residuals with unequal masses and complex phases.
13. **Soliton masses and stress tensor.** 10:96–102,151,174: `M_a=2 h m_a/beta^2` and the total-derivative stress tensor. Derive the static boundary term and independently integrate the canonical energy of regular one- and two-solitons.
14. **Dorey's rule, degeneracies and D4 figure.** 09:128,142–158: centre charge eigenvectors identify particles and the three outer orbits project to the same hexagon. Independently construct D4 reflections/orbits and cubic tensor; separate finite examples from the general established claim.
15. **Twisted charge spins and Cartan eigenvectors.** 09:88,98: the explicit spin lists and unqualified Cartan-eigenvector claim. Compare twisted Coxeter spectra, and test the spin-one mass vector for a non-simply-laced theory.
16. **Vertex-operator formula and polynomial degree.** 10:69–77: displayed tau matrix element, degree bounds and orbit degrees. The formula appears to omit spacetime evolution; nilpotency gives an upper bound only. Compare chapter 6 and Olive–Turok–Underwood, then test the explicit c_2 polynomial.
17. **Vacua and “minimum”.** 09:28–30,51: a complex potential has “a minimum” and its zero-gradient zero set is the coweight lattice for every affine theory. Check the real-field argument and a twisted affine-root constraint.
18. **Originality and source support for Theta.** 11:96, Appendix E: “The formulation as a selection rule from Theta is new here.” Inspect the cited reality papers and search relevant literature, including Delius; distinguish an elementary conjugation observation from a demonstrated new result.
19. **Missing static saddles.** 11:72: “no static saddle of the right energy carries the missing weights.” A later chapter's classification could be conditional or narrower than this wording. Read the referenced result and inspect its hypotheses before endorsing this summary.
20. **Exceptional masses.** 09:114–116: the D_n and E8 mass table. A plausible trigonometric entry may be misindexed. Build mass matrices from explicit roots and marks; compare eigenvalues and the trace rule.

## 2. Checks, in the original rank order

All calculations below are independent of the book's scripts. Executable evidence is in `referee_scratch/classical_09_11/codex/codex_checks_v2.py` and its `codex_checks_v2_output.txt`. “PASS” in that output means that the independent assertion passed, including assertions demonstrating an error in the book. The first version's curvature check produced a false FAIL because SymPy did not rewrite `sinh` as exponentials; v2 performs that rewrite and gives the zero matrix. Both versions are retained. No author's checks were run.

### 1. Singular energy integrals — wrong (major)

For sine-Gordon take the book's `m=beta=1`, hence `m_0=2`, and the static Hirota solution `tau_0=1+exp(2x)`, `tau_1=1-exp(2x)`. Direct substitution in the canonical energy gives

\[
\mathcal E(x)=-8\operatorname{csch}^2(2x)
=-\frac2{x^2}+\frac83+O(x^2).
\]

Its residue is indeed zero. But the integral excluding `(-epsilon,epsilon)` is

\[
E_\epsilon=8(1-\coth 2\epsilon)=-4/\epsilon+8+O(\epsilon),
\]

not a convergent energy. Independent outputs: `-392.053331911165289` at epsilon=.01 and `-3992.00533333191111` at .001. Even the symmetric Cauchy principal value diverges. A Hadamard finite part would give 8, but that is an extra definition, not the ordinary conserved energy. Near any simple tau zero on a time slice, `d_x^2 log tau=-1/(x-x_*)^2+O(1)` gives the same distinction. The implicit function theorem does support isolated spacetime zeros when the real/imaginary Jacobian has rank two.

**Replacement, `part3-classical/11-reality.qmd:60`:**

> Zeros of a complex tau-function are generically isolated spacetime points when its real and imaginary parts have independent gradients. Special data can instead give curves of zeros. At a simple spatial zero, the total-derivative energy density has a double pole with zero residue. The ordinary real-line energy integral nevertheless diverges. A finite-part or contour prescription may continue the formal boundary charges through the singularity, but this requires a separate definition and does not make the singular field a regular finite-energy solution.

**Replacement exercise, 11:91:**

> Show that a transverse zero of a complex tau-function is isolated in the real spacetime plane. For the static sine-Gordon data `tau_0=1+exp(2mx)`, `tau_1=1-exp(2mx)`, compute the Laurent expansion of the energy density and the energy with a symmetric interval about the zero removed. Explain why zero residue does not imply a finite energy integral.

### 2. Energy formula at arbitrary complex rapidity — wrong as stated (major)

There is a regular counterexample, not merely the singular case in item 1. In sine-Gordon set `theta=i*pi`, `xi=i*pi/2`, `m=beta=1`. Then

\[
\tau_0=1+i e^{-2x},\quad \tau_1=1-i e^{-2x},\quad
\phi=2\sqrt2\arctan e^{-2x}.
\]

This is an ordinary real static antikink, with energy `8` and momentum `0`. The proposition's unrestricted formula gives `8*cosh(i*pi)=-8`. Independent quadrature gives `8.0`.

The missing hypothesis is visible directly in the boundary evaluation. Write `k_i=m_ai cosh theta_i`. For a generic N-soliton polynomial with all `Re k_i>0` and nonzero leading interaction product, all exponentials tend to zero at the left end and the full product dominates at the right. Thus

\[
E=\frac2{\beta^2}\sum_j[\partial_x\log\tau_j]_{-\infty}^{\infty}
=\frac{2h}{\beta^2}\sum_i k_i.
\]

If `Re k_i<0`, the assignment of the two spatial ends changes. For the counterexample the derivative is `-2` at the left and 0 at the right, giving positive energy. At `Re k_i=0` even localization is not automatic. Fusion limits also change the leading polynomial.

**Source check.** Takács–Watts, [hep-th/9810006](https://arxiv.org/abs/hep-th/9810006), p.11, eq.4.6, really displays the sums of coshes and sinhes. This confirms the citation's formula, not the book's unrestricted domain. The following subsection explicitly calls `[-pi/2,pi/2]` the “fundamental region of α”. The source's formula uses its own mass normalization; I do not transplant that normalization into this book.

**Replacement, 11:32–39 (retain the displayed sums):**

> **[Established; reformulated]** For regular localized Hirota solitons in the standard spatial-decay sector, the energy and momentum densities are total derivatives. For the `a_n^{(1)}` N-soliton expression, suppose `Re(m_{a_i} cosh theta_i)>0`, the leading interaction product is finite and nonzero, and no tau-function vanishes on the integration slice. Boundary evaluation then gives the displayed additive energy and momentum. Extension to other complex data requires checking the spatial asymptotics and convergence; the sums cannot be applied to arbitrary complex rapidities without this check. Higher local charges require the additivity and phase-independence assumptions stated below.

This is a corrected sufficient domain for the explicit A-series formulas, not a new claim covering all exceptional/folded solutions.

### 3. Selection rule and Appendix E label — wrong in scope; conditional algebra correct (major)

Let `D={(a_i,theta_i)}` and suppose a particular charge has the well-defined value `F_Q(D)`, independent of all omitted phases and positions, on the regular solutions being discussed. If `Theta D=D` and `Q[Theta phi]=epsilon_Q conjugate(Q[phi])`, then

\[
Q[\phi]=F_Q(D)=F_Q(\Theta D)=Q[\Theta\phi]
=\epsilon_Q\overline{Q[\phi]}.
\]

That is a complete elementary conditional derivation. Additivity plus phase-independent one-soliton eigenvalues is one way to obtain `F_Q`; additivity alone does not imply phase independence.

“Every conserved charge” is false. The real sine-Gordon kink has invariant reduced data `(a=1,theta=0)` and topological charge `Q_top=1`. Under `phi -> -conjugate(phi)`, `Q_top -> -conjugate(Q_top)`, so its parity is -1. The claimed universal rule would require `1=-1`. The independent integral in the script gives `Q_top=1.0`. The topological charge depends precisely on the omitted phase. More generally one must choose a charge basis with definite antilinear parity; an arbitrary complex linear combination need not obey the same assigned sign.

The text correctly notices that a field invariant under this particular Theta is imaginary. A nonzero real sine-Gordon breather is consequently not an invariant *field*. Its reduced species/rapidity data are invariant. The final sentence of 11:26 conflates these notions. The same distinction must be made in 11:54. Real energy is a sufficient consequence under the stated conditions, not a necessary characterization of a physical sector; regular nonconjugate-species single solitons already illustrate the difference.

**Replacement, 11:45:**

> **Conditional selection rule.** Consider a local hierarchy charge of definite Theta parity whose value is determined by the species and rapidities alone. If this reduced data multiset is invariant under `(a,theta) -> (bar a,bar theta)`, then `Q=epsilon_Q conjugate(Q)`. In particular, the regular energy and momentum of the preceding proposition are real. For higher charges we assume additivity and independence of the internal phases in the domain under consideration. The assertion does not cover topological charges, which depend on those phases.

**Replacement, final sentence of 11:26:**

> A breather has Theta-invariant reduced species/rapidity data; this does not assert invariance of its field.

**Replacement, 11:49:**

> This gives a sufficient reality criterion for the stated charges on regular solutions. It does not by itself select a quantum physical sector or resolve the singularity problem.

**Appendix E, classical item:** recommend **[Established; reformulated, conditional]**, and describe only this restricted implication. The brief's referee rule excludes promoting convention/symmetry reformulations to a theorem. The purported originality is separately assessed in item 18. Chapter 11's learning goal promising real *higher* charges also needs “real or imaginary according to their parity”.

### 4. Breather range and the figure — wrong without a regularity restriction (major)

For `b=h-a`, `theta_1=-theta_2=iu`, the printed interaction coefficient becomes

\[
A=1-\frac{\sin^2(\pi a/h)}{\cos^2 u}.
\]

The nonsingular conjugate-species (type-B) family requires `A<0`. For `0<u<pi/2`, this means

\[
u>\left|\frac\pi2-\frac{\pi a}{h}\right|,
\qquad 0<M_B<2M_a\sin(\pi a/h).
\]

Here is an independent zero argument, allowing unequal constituent magnitudes. After translations write

\[
\tau_j=1+r e^{i\delta}(s e^{i\psi_j}+s^{-1}e^{-i\psi_j})
       +A r^2 e^{2i\delta},\quad r>0,\ s>0.
\]

As x and t vary, r ranges over positive numbers and each `psi_j` over a circle. When `0<A<1`, the origin lies on one of these ellipses for some r. For `s!=1`, divide by `r exp(i delta)` and test the centre against the ellipse: the squared normalized distance is

\[
f(r)=\frac{(r^{-1}+Ar)^2\cos^2\delta}{(s+s^{-1})^2}
+\frac{(Ar-r^{-1})^2\sin^2\delta}{(s-s^{-1})^2}.
\]

It tends to infinity at both ends, but `f(1/sqrt(A))=4A cos²(delta)/(s+s^-1)²<1`. Continuity gives a tau zero. For `s=1`, take `r=1/sqrt(A)` and `cos psi_j=-sqrt(A) cos delta`. Thus no choice of phases removes the obstruction. Conversely, for `A<0`, equal magnitudes and `delta=pi/2` give `Re tau_j=1+|A|r²>0`, an explicit regular family.

For A2, `u=.2` gives `A=.2191814811280546`; at `r=2.135984374348023`, `psi=pi/2`, `delta=pi/2`, the script finds a zero to `2e-36`. For `u=pi/3`, `A=-2`, the positive-real-part construction is regular. Consequently the A2 band is `0<M_B<sqrt(3) M_a`, not all of `(0,2M_a)`.

**Source support.** Harder–Iskandar–McGhee, [hep-th/9409035](https://arxiv.org/abs/hep-th/9409035), §3.2, states exactly: “it is not possible to have type B breathers with positive interaction coefficient”. Their paper also treats charged same-species breathers; the book describes only the neutral type-B family. Takács–Watts §4.2 independently gives the A2 bound. These are restrictions already in the cited sources, not corrections to them.

The exercise with `xi_2=conjugate(xi_1)` is sound in the nondegenerate range: the tau-functions are real; for `A<0` their large-r sign is negative, forcing zeros. For sine-Gordon the regular choice `xi_1=-log(tan u)`, `xi_2=xi_1+i*pi` gives, with `y=m_0 x cos u`, `z=m_0 t sin u`,

\[
\tau_0=1+e^{2y}-2i\cot u\,e^y\sin z,\quad
\tau_1=\overline{\tau_0}.
\]

The field is the sine-Gordon breather, up to overall sign/time origin. It is time-periodic in its rest frame; a boosted localized breather is not periodic at fixed laboratory x.

**Replacement, 10:130 and the range explanation at 10:136:**

> A soliton and its conjugate species can form a regular neutral breather with rapidities `theta_0 +/- iu`, provided `|pi/2-pi a/h|<u<pi/2` and the internal constants satisfy the regularity conditions. The solution is periodic in its rest frame. Its rest energy is `2M_a cos u`, so this regular family fills `0<M_B<2M_a sin(pi a/h)`. For a self-conjugate species the upper bound is `2M_a`. Other breather families can carry topological charge.

**Figure caption, 11:28:** retain the cosine as a formal energy curve, but indicate the allowed species-dependent interval and mark the excluded singular part. The existing SVG was not regenerated.

### 5. Lax pair — correct curvature and boost; wrong description for twisted algebras

Put `U_j=exp(beta alpha_j.phi/2)`. The coefficient of `e_j` in the curvature is

\[
-m\lambda c_j U_j\frac\beta2\alpha_j.\partial_-\phi
+m\lambda c_j U_j\frac\beta2\alpha_j.\partial_-\phi=0,
\]

and similarly for `f_j`. The Cartan term is

\[
-\beta\Box\phi.H-m^2\sum_j n_j e^{\beta\alpha_j.\phi}\alpha_j.H.
\]

Thus the normalization and sign reproduce the printed real field equation. Independent 2x2 matrix differentiation for the stated A1 exercise gives the zero matrix after substitution of `Box phi=-(2sqrt(2)m²/beta)sinh(sqrt(2)beta phi)`. Under the printed active boost, `partial_+ phi` scales by `exp(-vartheta)`, `partial_- phi` by `exp(vartheta)`, and the spectral terms transform exactly with `lambda -> exp(vartheta)lambda`. No factor of two is missing.

The generator description is restricted to the untwisted case. In A2^(2), the projected affine root is `alpha_0=-2alpha_1`; the fixed rank-one finite algebra has only roots `+/-alpha_1`. Hence `E^{-theta}` in that algebra cannot be the affine generator. Abstract twisted affine Chevalley generators still give the displayed curvature calculation.

**Replacement, 09:57, the two sentences specifying generators:**

> For an untwisted affine algebra, the generators with `j>=1` belong to the finite algebra and `e_0`, `f_0` are its negative/highest-root vectors multiplied by the loop variable and its inverse. For twisted affine algebras use the Chevalley generators in the appropriate twisted loop algebra; the affine generator need not be a root vector of the fixed finite subalgebra.

### 6. Fusion — wrong without ranges; correct for a+b<h

Write `A=pi a/h`, `B=pi b/h`. Independently,

\[
\sin A e^{iB}+\sin B e^{-iA}=\sin(A+B).
\]

The same identity with both exponents conjugated supplies the other light-cone component. For `a+b<h` the exponent of the product is therefore exactly the desired soliton exponent, and scaling both individual amplitudes to zero with their product times the simple interaction pole fixed gives the stated tau limit.

At `a+b=h`, the sum is zero and `omega^{j(a+b)}=1`: the tau limit is a common constant, giving a vacuum, not a new species. Also `cosh(i*pi+epsilon)+1=-epsilon²/2+...`, so the interaction pole is double, not a generic simple fusion pole. At `a+b>h` the displayed sine mass is negative. For A2 the script obtains `sqrt(3),0,-sqrt(3)` for `(a,b)=(1,1),(1,2),(2,2)` respectively.

For the upper-sum channel use `c=a+b-h`, `theta_1=theta+i*pi(h-b)/h`, `theta_2=theta-i*pi(h-a)/h`; their mass-weighted exponent sums equal the positive `m_c`, and the rapidity difference is `i*pi(2h-a-b)/h`.

**Replacement, 10:126:** add “For `a+b<h`” to the first sentence, then append:

> For `a+b>h`, the physical fusion channel has species `a+b-h` and rapidity difference `i*pi(2h-a-b)/h`; use `theta_1=theta+i*pi(h-b)/h` and `theta_2=theta-i*pi(h-a)/h`. The case `a+b=h` is not fusion to a soliton species.

In 10:163 specify species `1+1 -> 2` (or also treat `2+2 -> 1` using the second prescription).

### 7. Area rule — wrong as an unconditional displayed identity

In A2 take the three cyclic differences in R³, with charge eigenvectors `(1,omega^a,omega^(2a))/sqrt(3)`. Then

\[
C_{abc}=m^2\beta\sum_j(\alpha_j.e_a)(\alpha_j.e_b)(\alpha_j.e_c).
\]

The geometric sum vanishes unless `a+b+c=0 mod 3`. Independently, `|C_111|=|C_222|=3m² beta`, while `C_112=C_122=0`. All masses are `sqrt(3)m`, so *all* these triangles have area `3sqrt(3)m²/4`; the displayed area-rule right side is `3m² beta`, even for `112`. The numerical output has `|C_112|=6.1e-16`, with right side 3 at unit parameters. The missing selection-rule factor is substantive, even though the following section informally acknowledges it.

The dimensionful version works directly if Delta uses the actual masses. If Delta is evaluated with m=1 while the coupling retains general m, an additional m² is needed.

**Replacement, 09:130–136:**

> For simply-laced theories, every **nonzero** cubic coupling in this charge basis obeys the following area rule. Here `Delta_abc` is the area of the triangle of the actual masses `m_a,m_b,m_c`. Couplings forbidden by Dorey's rule vanish even when that triangle exists.

Retain the equation, remove “(with m=1)”, and retain the discussion of the fusing angles. This calculation checks the normalization and the counterexample, not a proof of the rule for all ADE types.

### 8. Time delays — correct shift and semiclassical sign; wrong unrestricted inequality

On the trajectory of particle 1, the second exponent contains `m_b cosh(theta_2)(v_1-v_2)t`. For `theta_1>theta_2`, it goes from zero to infinity. Thus `delta xi_1=log A`, `delta x_1=-log A/(m_a cosh theta_1)`. The slower particle instead has `delta xi_2=-log A` and the opposite displacement. Further,

\[
\cos\frac{\pi(a-b)}h-\cos\frac{\pi(a+b)}h
=2\sin\frac{\pi a}h\sin\frac{\pi b}h>0.
\]

Hence `0<A<1` for distinct real rapidities. At equal rapidities and equal species, `A=0`; this degenerate case is not a scattering collision.

The sign of `delta t=-delta x/v` is frame dependent. With A2 species 1, `theta_1=-.4>theta_2=-1`, the independent output is `delta x_1=1.17861800503`, `delta t_1=3.10204296397>0`. The assertion `<0` requires positive `v_1`, for example the rest frame of particle 2. At `v_1=0`, the passage-time definition is unavailable.

The sine-Gordon semiclassical comparison is **correct**. With chapter 3's small `xi`, its kernel is `2 tanh(pi t/2)/(xi t)` to leading order. The Abel-regularized cosine transform obeys

\[
\int_0^\infty\frac{dt}{t}\tanh(\pi t/2)\cos(\theta t)
=-\log\tanh(\theta/2),\quad\theta>0.
\]

One may derive it by expanding `tanh=1-2/(e^{pi t}+1)`, differentiating with respect to theta, using the partial-fraction expansion of `csch(theta)`, then fixing the constant at theta=infinity. Therefore `delta'(theta)=(2/xi)log tanh(theta/2)`, and `xi~beta_SG²/8` gives the printed `16/beta_SG²`. Also `(2h/beta²)log A` gives the same coefficient for h=2 and `A=tanh²(theta/2)`.

**Replacement, 10:122, inequality sentence:**

> For distinct real rapidities, `0<A_ab<1`, so soliton 1 has a positive spatial displacement. If `v_1!=0`, its passage-time shift is `Delta t=-Delta x/v_1`. In the rest frame of soliton 2, `v_1>0` and this is a time advance.

### 9. Spin-3 current — correct identity, opposite rapidity dependence

Let `p=partial_+phi`, `r=partial_+²phi`. From the sinh-Gordon equation,

\[
\partial_-p=-m²\sinh(\beta\phi)/\beta,\qquad
\partial_-r=-m²\cosh(\beta\phi)p.
\]

Differentiating the printed T and X gives respectively `-2m² cosh(beta phi)pr-beta m² sinh(beta phi)p³` and its negative. The symbolic check is exactly zero. The equation is `partial_t(T+X)=partial_x(T-X)`, so the *sum* is the right conserved density.

But the book's moving kink has `z=m_0(x cosh theta-t sinh theta)` and `partial_+ z=m_0 exp(-theta)`. With `b=beta_SG`, its current is

\[
T=(\partial_+²\phi)^2-\frac{b²}4(\partial_+\phi)^4,
\quad X=m_0²\cos(b\phi)(\partial_+\phi)^2.
\]

Using `phi_z=2 sech z/b` and the integrals of `sech² z` and `sech⁴ z` gives

\[
\int T\,dx=-\frac{8m_0³}{3b²}\frac{e^{-4\theta}}{\cosh\theta},\quad
\int X\,dx=-\frac{8m_0³}{3b²}\frac{e^{-2\theta}}{\cosh\theta},
\quad Q=-\frac{16m_0³}{3b²}e^{-3\theta}.
\]

For `m_0=b=1`, independent quadrature gives `Q(0)=-5.3333333333333333`, `Q(.6)=-.8815940705151282`, ratio `.16529888822158654=exp(-1.8)`. With the rapidity convention in chapter 2, this is the spin **-3** member. The parity-reflected current has spin +3. The general potential parameters in the sinh-Gordon example are conventional scalar parameters, not automatically the root-normalized `m,beta` used earlier; this dictionary should be explicit in the exercise.

**Replacement, 09:96:**

> With our rapidity convention this current gives `Q_{-3}=integral dx (T_4+X_2)`, whose one-particle value is proportional to `e^{-3 theta}`. Exchanging `partial_+` and `partial_-` gives `Q_3`. For sine-Gordon replace the quartic coefficient by `-beta_SG²/4` and `cosh(beta phi)` by `cos(beta_SG phi)`, using the sine-Gordon parameters `m_0,beta_SG`.

### 10. Topological-charge cycle — correct

Use `alpha_j=e_j-e_{j+1}` for `1<=j<h`, `alpha_0=e_h-e_1`. Put `r_j=[delta+2pi ja/h]_{(-pi,pi)}`. The k-th Cartesian coordinate of the charge is

\[
Q_k=\frac{r_{k-1}-r_k}{2\pi}\in\{-a/h,\ 1-a/h\},
\]

where the final cyclic difference uses r_0. Their sum is zero, so precisely a entries equal `1-a/h`. These are exactly weights of `Lambda^a C^h`, with the representation orientation asserted in the text. There are `h/gcd(a,h)` distinct phases of `omega^{ja}`; crossing each forbidden ray changes a nonempty set of coordinates and cycles through that many distinct weights. No change occurs inside a chamber. The excluded endpoints are exactly singular rays.

The independent enumeration checks every species for `h=2,...,14` (91 cases): the predicted number of distinct charges and the weight condition hold, with maximum floating residual `4.11e-15`. This also confirms the A3 figure's four-versus-two count and exercise. It does not prove that no other classical solutions carry the missing weights, nor the separate assertion about defects.
### 11. Lie duality and folded masses — correct for the stated C_n example

Direct Cartesian constructions avoid borrowing the author's folding code. For B_n use `alpha_0=-e_1-e_2`, `alpha_i=e_i-e_{i+1}` for `i<n`, `alpha_n=e_n`, with marks `(1,1,2,...,2)`. For C_n use `alpha_0=-sqrt(2)e_1`, `alpha_i=(e_i-e_{i+1})/sqrt(2)`, `alpha_n=sqrt(2)e_n`, with marks `(1,2,...,2,1)`. Form `sum n_j alpha_j alpha_j^T` in each case.

For B_n, the first coordinate has eigenvalue 2 and decouples; the remaining tridiagonal chain has diagonal 4 and off-diagonal -2. Its eigenvalues are `4-4cos(pi a/n)`, `a=1,...,n-1`. For C_n, the tridiagonal matrix has end diagonals 3, interior diagonals 2, off-diagonals -1; its eigenvalues are `4sin²(pi a/(2n))`, `a=1,...,n`. Consequently the B_n masses are `sqrt(2)` times `{1,2sin(pi a/(2n))}`, as required.

Independent matrices at `n=2,3,4,7` agree to at most `7.2e-16`. Orthogonal restriction to the diagram-fixed subspace preserves the kinetic norm, so a regular coincident parent pair has the sum of the two parent energies; the singleton species n has one. This yields exactly the claimed C_n soliton ratios. The explicit c_2 polynomial is checked in item 16. The claim for every other untwisted non-simply-laced algebra is not independently verified here; no general Lie-duality proof is inferred from these examples.

### 12. Hirota equations — correct

From `phi=(i/beta)sum_j alpha_j^vee log tau_j`, the exponential is `prod_l tau_l^{-a_lj}`. Multiplication of the proposed log-tau equation by `alpha_j^vee` gives `m² n_j alpha_j(prod-1)`. The constant sums to zero, and multiplication by `i/beta` gives exactly the imaginary field equation. Thus the transpose in the printed product is necessary and correct. This is a sufficient choice of tau equations; the redundant coroot combination means it should not be called a unique gauge-free equation for tau.

Finally,

\[
\frac12(D_x²-D_t²)\tau.\tau
=\tau(\tau_{xx}-\tau_{tt})-\tau_x²+\tau_t²
=-\tau²\Box\log\tau.
\]

Using `a_jj=2` gives the displayed right side. For one soliton the left side is `(k²-nu²)omega^{ja}E=m_a² omega^{ja}E`, proving the moving formula with no derivative convention ambiguity.

The independent finite-exponential evaluator differentiates each monomial analytically and checks the bilinear residuals. Tests use A2 and A3 two-solitons with distinct complex phases, unequal masses in A3, and an A2 conjugate-rapidity breather. At three spatial points the largest normalized residual is `3.87e-35` (35-digit arithmetic). No finite-difference or author's helper was used.

### 13. Soliton energy and stress tensor — correct in the regular sector checked

For a static solution the field equation gives

\[
\frac{d}{dx}\left(\frac12\phi'^2-V\right)=0.
\]

Vacuum boundary conditions make this constant zero. From the independently derived bilinear equation,

\[
V=\frac1{\beta²}\sum_j\frac2{\alpha_j²}\partial_x²\log\tau_j,
\qquad E=\frac2{\beta²}\sum_j\frac2{\alpha_j²}
[\partial_x\log\tau_j]_{-\infty}^{\infty}.
\]

For A-series single solitons every tau has degree 1, so `E=2hm_a/beta²`. In a simply-laced case with degree exactly n_j, the same computation gives `sum n_j=h` and the same formula; that general extension depends on the asserted degrees.

Independent canonical-energy quadrature, without complex conjugating field derivatives, gives `10.3923048454132638` for A2 species 1 and `16.0` for A3 species 2 at `m=beta=1`, at phases .3 and 1.1 respectively. These equal the boundary answers. The finite-exponential calculation also compares the canonical density to `2sum_j partial_x² log tau_j` for the two-solitons and the regular breather in item 12: maximum discrepancy `1.47e-34`.

For the trace, the proposed tensor has `T^mu_mu=Box C=2V`, which independently follows from the log-tau equations. This confirms its sign and prefactor. The tested full density identity is stronger than the trace check, but it is not a proof for every exceptional algebra. Neither the boundary identity nor this computation licenses integration through poles (item 1) or arbitrary complex rapidities (item 2). Replace “and the energies are real” at 10:174 by “and their energies follow from the spatial boundary terms; the regular real-rapidity solitons have real energies.”

### 14. Dorey rule and D4 figure — correct in independently tested cases; established source supports the rule

The extra script `codex_algebra_checks.py` constructs D4 from `e_1-e_2,e_2-e_3,e_3-e_4,e_3+e_4`, with outer nodes black and the central node white. It computes `w=s_1 s_3 s_4 s_2`, fundamental weights by inversion, and `(1-w^-1)lambda_a`. Results:

- Four orbits, each of length 6, partition all 24 roots.
- Projection to the eigenplane of `exp(i*pi/3)` gives radius ratios `(1,sqrt(3),1,1)`.
- The three outer hexagons coincide as sets.
- An explicit closing triple in orbits `(1,1,2)` is `(1,0,-1,0)+(-1,1,0,0)+(0,-1,1,0)=0`.

The triangle's sides are `1,1,sqrt(3)`, giving angles `pi/6,pi/6,2pi/3` by the cosine rule. This confirms the caption independently.

`codex_d4_couplings.py` constructs the two commuting centre transformations as orthogonal maps on the roots, diagonalizes their joint projectors, and computes the cubic tensor in the resulting charge basis. All 20 unordered cubic triples agree with the root-orbit test. The nonzero triples are `(112),(134),(222),(233),(244)`; their area-rule errors are below `8.89e-16`. Independently in A2, the orbit condition holds only for `(111),(222)`, agreeing with item 7. These checks substantiate the need for a charge basis in a degenerate mass space.

**Source check.** Dorey's sequel [hep-th/9110058](https://arxiv.org/pdf/hep-th/9110058), eq.2.8, uses precisely the book's white/black orbit representatives, and its p.16 summarizes the same root-triangle criterion. Its area-rule discussion on p.12 explicitly says “which, if non-vanishing, obey the ‘area rule’”. This independently supports the correction in item 7. I did not obtain the complete pre-arXiv 1991 paper cited by `dorey1991`; attribution to that particular source remains unchecked, while the sequel is directly verified. Retain **[Established]**, not a new theorem label. Full exceptional cubic tensors were not recomputed.

### 15. Twisted spin lists and the Cartan-eigenvector statement — unclear in full scope

An independent arithmetic consistency check of the printed residue classes finds no contradiction in the dual-pair statement: `{1,5,7,11}+12Z` and `{1,5}+6Z` both equal integers coprime to 6; `{1,5,7,11,13,17}+18Z` is the same set. Thus the listed D4^(3)/G2 and E6^(2)/F4 spin sets agree; the B/C pairs both have odd spins. This checks equality of the asserted sets, not that the Heisenberg algebra actually has those degrees or multiplicities. Kac's source and an independent twisted Heisenberg construction were not obtained, so the complete classification is **unclear**, not certified by agreement between chapters.

The Cartan statement needs left/right and classical/quantum qualifications outside the simply-laced case. For example the directly computed B3 mass vector in simple-node order is proportional to `v=(1,sqrt(3),1)`. With the book's `a_ij=alpha_i^vee.alpha_j`,

\[
A=\begin{pmatrix}2&-1&0\\-1&2&-1\\0&-2&2\end{pmatrix},\quad
A^T v=(2-\sqrt3)v,
\]

whereas `Av` is not proportional to v. This does not disprove a statement allowing left eigenvectors. It does prevent an unqualified transfer of the symmetric simply-laced formulation. The higher quantum charges of non-simply-laced theories have not been calculated in this pass.

**Suggested clarification, 09:98:** “For simply-laced theories, their one-particle eigenvalues form eigenvectors of the Cartan matrix.” This is the scope actually used immediately below; broader claims should specify conventions and quantum corrections.

### 16. Vertex formula — wrong as displayed; explicit polynomial and nilpotency bound supported

The matrix element in 10:72 contains only constant `Q_k`, `z_k` and operators. Its derivatives with respect to x and t are both zero. It therefore cannot equal the nonconstant tau-function `1+omega^{ja}exp(m_a x+xi)` even in the A-series one-soliton case. Chapter 6 includes the missing factors `W_k(x,t)`.

**Replacement, 10:72 and 10:75:**

\[
\tau_j(x,t)=\langle\Lambda_j|\prod_k
\exp\!\big(Q_k W_k(x,t)F^{a_k}(z_k)\big)|\Lambda_j\rangle,
\]

> Here `W_k(x,t)` is the spacetime exponential defined in `sec-vertex-operators`; `Q_k` is its constant position and phase factor.

Nilpotency of order `n_j+1` implies degree **at most** n_j by expanding the exponential; it does not prove a nonzero top coefficient. The book makes that distinction and records its stronger degree claims as numerical/conditional. Olive–Turok–Underwood [hep-th/9305160](https://arxiv.org/abs/hep-th/9305160), p.2, states that powers above `2x/gamma(i)^2` vanish at level x and adds exactly “while lower powers usually have finite matrix elements.” For simply-laced roots this gives the printed bound.

The independent symbolic c_2 check sets `d_x=sqrt(2)E d_E` and substitutes `(1+2E+E²/2,1+E²/2,1-2E+E²/2,1+E²/2)` into all four parent Hirota equations. Residuals: `[0,0,0,0]`. This confirms the example, reflection invariance and its degree 2. As a real-E polynomial it can have zeros, so its formal validity is not a claim of regularity for every phase. The universal folded degree assertion remains conditional, not independently proved here.

### 17. Vacua — correct; “minimum” should mean the real part on real fields

Let `z_j=exp(i beta alpha_j.phi)`. The r by (r+1) root matrix has one-dimensional nullspace generated by the positive marks. If the gradient vanishes, `sum_j n_j z_j alpha_j=0` implies `z_j=z` for all j. If also `V=0`, then `h(z-1)=0`, so every `z_j=1`. This proof works even for complex fields. The equations force every `alpha_j.Im(phi)=0`, hence `Im(phi)=0` since the roots span. For `n_0=1`, the affine condition is redundant once the finite simple-root phases are integral multiples of `2pi`; one gets the stated coweight lattice.

For example the rank-one twisted data `alpha_0=-2alpha_1`, marks `(1,2)`, impose no extra fractional-period vacuum: `alpha_1.phi=2pi k/beta` is sufficient and necessary. Thus the root-length difference does not defeat the lattice argument.

For real fields, `Re V_imag=(m²/beta²)sum_j n_j[1-cos(beta alpha_j.phi)]>=0`. Its Hessian at zero is the stated positive-definite matrix. A complex-valued V itself has no ordered minimum. **Replacement, 09:30:** “The sign makes the quadratic term at zero positive definite; on real fields, the real part of the potential is minimized there.” The positive-fugacity assertion follows directly by expanding `exp(-S_E)` with this sign and is correct.

### 18. Originality of the Theta formulation — unclear; broad novelty is not established

The independent proof in item 3 is a general identity: an antilinear involution acting on a charge, together with invariance of the data determining that charge, gives its reality/parity property. It does not require special Toda structure. Its Toda-specific content lies in proving which charges depend on the reduced data and on which regularity domain. Those are exactly the hypotheses the chapter partly assumes.

Relevant prior art found and read:

- Zhu–Caldi, [hep-th/9404097](https://arxiv.org/abs/hep-th/9404097), eq.2.5, pairs tau-functions with their conjugates to obtain a real **density**. The source describes this by the exact phrase “The above equation will be referred to as the ‘pairing condition.’” Its density restriction is stronger and different from phase-independent integrated energy. This paper should not be paraphrased as the same solution-selection criterion.
- Cen–Fring, [Complex solitons with real energies (2016)](https://openaccess.city.ac.uk/id/eprint/15266/1/complexsolitons.pdf), pp.1–2, already discusses antilinear symmetry, real integrated energy and decomposition into asymptotically separated constituents, including complex sine-Gordon. Their symmetry acts on spacetime too, unlike the book's field-only Theta.
- Cen–Correa–Fring, [arXiv:1608.01691](https://arxiv.org/html/1608.01691), §3, states exactly: “The reality of all charges build on one-soliton solutions is guaranteed by” PT symmetry. The paper then uses additive separated-soliton contributions. Its model is KdV, so this is prior art for the reasoning, not proof that the precise affine-Toda reduced-data criterion has appeared before.
- A targeted search including Delius found his work on antilinear boundary constructions, but no verified source with this exact bulk criterion. I do not infer novelty from that negative search. Fring–Taira [arXiv:2007.15425](https://arxiv.org/abs/2007.15425) is another relevant antilinear/BPS reality treatment; only its abstract was inspected here.

Therefore “new here” is **unclear**, and the stronger suggestion that symmetry-based reality arguments originate here is unsupported. A complete novelty search is not established by this pass.

**Replacement, 11:96 final sentence:**

> Here we organize the integrated-charge reality conditions using the field antilinear symmetry Theta and the reduced Hirota data. This is distinct from conditions for pointwise real energy density. Related antilinear-symmetry arguments for complex solitons occur in the literature; the higher-charge statement here is conditional on additivity and phase independence.

The editor should add the relevant bibliography entries if using this replacement and revise Appendix E's originality claim accordingly. No entries were edited by this referee.

### 19. Missing static saddles in A3 — correct for regular solutions approaching vacua; stronger surrounding numerical assertions unchecked

I read the referenced proposition in chapter 17 and independently checked the mechanism rather than rerunning its script. At `m=beta=1`, the static linearized field equation has three positive rates `(sqrt(2),2,sqrt(2))`. No nonlinear sum of two or more positive rates equals one of these rates, since `2sqrt(2)>2`; the repeated rate is a degeneracy, not such a resonance.

There is also a direct local-completeness argument. Construct the static three-species Hirota family with free amplitudes `(p,q,r)` and pair factors `A_12=A_23=3-2sqrt(2)`, `A_13=1/2`. The new symbolic script `codex_static_a3_check.py` verifies all four bilinear equations identically as polynomials: residuals `[0,0,0,0]`. The field's three linear tangent vectors are Fourier combinations `i sum_j alpha_j i^{ja}`, `a=1,2,3`, with rank 3. Thus at a sufficiently negative x this analytic family gives a local coordinate chart of the three-complex-dimensional unstable manifold of the vacuum. Uniqueness of the analytic ODE extends the identification along any regular continuation. Any regular static solution tending to that vacuum as x goes to minus infinity is therefore in this family; a vacuum shift handles other starting vacua.

For any subset of nonzero amplitudes, the unique leading full-subset term at plus infinity has a nonzero coefficient (all three pair factors are nonzero). The boundary energy is the corresponding subset sum of `(8sqrt(2),16,8sqrt(2))`. The independent enumeration is

`11.3137084989848, 16, 11.3137084989848, 22.6274169979695, 27.3137084989848, 27.3137084989848, 38.6274169979695`.

Only the singleton species 2 has energy 16. Item 10 proves it reaches only two of the six weights. Hence the four missing weights cannot be carried by a regular static vacuum-to-vacuum solution of this energy. This supports the chapter 11 summary with its usual soliton boundary conditions. It does not cover singular fields with assigned finite-part energies, lattice critical points or time-dependent complex saddles. I did not independently verify chapter 17's additional enumeration that only the triple-species subset reaches the missing weights, or its bound-state charge scan and thimble conjecture.

**Scope clarification, 11:72:** replace “no static saddle” by “no regular static solution approaching vacua at both spatial ends”. This makes explicit what the argument proves.

### 20. D_n and E8 mass tables — correct

The independent root construction uses standard D_n simple roots `e_i-e_{i+1}`, `e_{n-1}+e_n`, highest-root marks `(1,2,...,2,1,1)`, and the negative highest root as alpha_0. The mass matrix yields the table for `n=4,5,6,7,8,10`, with eigenvalue-derived mass errors at most `6.67e-16`; the traces are respectively `12,16,20,24,28,36`, equal to `2h`.

For E8 the script specifies all eight Cartesian simple roots (including the half-integer root) and highest-root marks `(2,3,4,6,5,4,3,2)`. All nine affine roots have squared length 2. Diagonalizing their marked outer-product sum gives trace 60; every displayed trigonometric mass ratio agrees to `2.67e-15`. This verifies the exceptional table independently of the book's root-data helper.

## 3. Items not fully checked and limitations

All twenty ranked entries received an independent derivation, calculation or a clearly delimited partial check. The top ten have verdicts and evidence. The following parts remain **unclear / not checked**, rather than silently endorsed:

- **Item 11:** general Lie duality beyond the explicitly derived C_n/B_n case; Kneipp–Olive's full all-algebra proof was not audited.
- **Item 13:** a general proof of the full stress tensor identity for all exceptional/folded solitons and the inaccessible pre-arXiv `olive1993a` source. The static derivation, trace and explicit A-series tests are the evidence here.
- **Item 14:** the full text of `dorey1991` was not available; its sequel was checked directly. No exceptional cubic tensors were recomputed. The full centre-basis classification in every accidental degeneracy remains outside the finite A2/D4 tests.
- **Item 15:** the twisted Heisenberg spectrum and multiplicities, and the non-simply-laced higher quantum-charge assertion. No independent construction or complete source verification was obtained. This is the principal ranked entry with an overall unclear verdict.
- **Item 16:** universal exact polynomial degrees and all coincident-orbit bounds. The cited nilpotency bound and explicit c_2 example do not prove those stronger statements.
- **Item 18:** exact novelty of the particular Toda reduced-data formulation. Relevant prior art was found, but a universal negative literature claim cannot be certified by this search.
- **Item 19:** the later numerical bound-state charge enumeration, intersection numbers and all quantum interpretations; only the static continuum obstruction under regular vacuum boundary conditions was checked.

The reviews' purpose is comparison. This report does not import any verdict from the earlier Claude reports or Gemini, and phase 2 has not been read or started.

## 4. Other observations (unverified unless explicitly stated above)

- The sentence at 10:59 gives position `-Re xi/m_a` without saying “in the rest frame”. Direct inspection of the printed exponent suggests laboratory intercept `-Re xi/(m_a cosh theta)`. This minor wording issue was not developed as a ranked check.
- Several broad quantum statements in Part III, including completeness of quantum multiplets and the restriction to a “physical sector”, need assumptions supplied by later parts. This pass does not validate them from classical energy reality.
- The general claim that a real-coupling theory's complete spectrum consists of the classical mass-matrix particles “and their bound states” is not a theorem supplied by convexity. I did not audit its quantum meaning.
- Figure inspection used read-only conversion of the four existing SVGs to Codex PNG previews with `rsvg-convert`; the author's figure scripts were not run. The geometric contents agree with the checks above, apart from the breather-domain issue in item 4. These previews are not rendered-book layout QA.

### Evidence and handoff

The final scripts and outputs are also copied to `briefs/referee_classical_09_11_codex_scripts/` so the evidence survives scratch cleanup. Every filename there includes `codex`. The snapshot manifest identifies the chapter text reviewed. The only false FAIL in a preliminary script was the symbolic simplification issue documented at the start of section 2; all final script assertions pass. No git mutation, book edit, author's script execution or Quarto render was performed.

The most consequential corrections are the treatment of singular integrals, the domain of the energy proposition, the scope of the Theta selection rule, and the allowed breather band. Other concrete corrections concern twisted Lax-generator wording, fusion ranges, the nonzero-coupling condition in the area rule, the frame restriction for time advances, the spin sign of the displayed current, and the missing spacetime factor in the vertex formula.

Phase 1 complete. This report is frozen on completion; phase 2 requires the editor's separate list.

## Integration

Integrated on 2026-10-24 on branch `part3-referee-integration`, together with the other three Part III reports, after independent checks (AGENTS.md, Rule 1). The verdict on each finding, the coordinator's checks, the decisions left to the editor and the overlap table of the two-phase experiment are in the `## Integration` section of `briefs/referee_classical_09_11_phase1_combined.md`.

- **Verified and integrated:** C1–C9, C13 (scope of "the energies are real"), C16–C19, and the remark on the soliton position. C15 (twisted spins, "unclear") was spot-checked for $a_4^{(2)}$: correct, no edit.
- **Confirmed correct:** C10–C12, C14, C20.
- **Rejected:** none. No false claim was found in this report.
