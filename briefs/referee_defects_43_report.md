# Independent referee report: chapter 43

Date: 2026-10-13. Audited commit: `68a9bea6b0556869751b13c9deb0a6c4405900a3`, branch `claude/integrable-defects-framework-a54a84`, worktree `.claude/worktrees/integrable-defects-framework-a54a84`.

Models: **GPT-6 (Codex)** wrote items 1–4 and the four original scripts, and did the source audit. When Codex ran out of quota, **Claude Opus 5.5** took over from its saved handoff. Claude re-derived every handoff conclusion it uses, added `second_module_checks.py`, wrote items 5–8 and the closing observations, and assembled the report. Where a section says "second referee", the derivation is Claude's; numbers marked "Codex" come from Codex's scripts.

The chapter reviewed is the current, partly integrated version, not only its original draft at `f677c84`. Its SHA-256 is `a3c56ae44321eabeab1b56a70e8e6c6a5024614d67b56d233a33a9083f715c66`. Line numbers refer to [the saved review copy](../referee_scratch/defects_43/chapter43_snapshot.qmd), which is byte-identical to the current file. At assembly, all nine Part VIII chapter files, `defect_common.py` and `defect_particle_transmission.py` still matched the hashes recorded at the start of the review ([initial_hashes.json](../referee_scratch/defects_43/initial_hashes.json)). All file/line corrections refer to `part8-defects/43-particle-transmission.qmd` unless another file is specified.

No pre-existing file was edited, no index/working-tree-changing git command was run, and no Quarto render was run. All scripts and their outputs are under [referee_scratch/defects_43](../referee_scratch/defects_43/). The named local source archive was absent; primary papers were downloaded from arXiv into that directory. The [citation manifest](../referee_scratch/defects_43/citation_manifest.json) records all 29 citation keys, source URLs, cited lines, and text hashes. The relevant Bajnok–Simon and CZ equations were also inspected in rendered PDF pages.

The computations below are independent: none imports a textbook computation or a previous referee's code, and no book verification script was used as evidence. Run them from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_43/independent_breathers.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_43/independent_doubled.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_43/independent_particles.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_43/supplemental_checks.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 referee_scratch/defects_43/second_module_checks.py
```

The main breather identity survives the independent computation, including every species of ranks 2–4. The particle identities and the vacuum theorem also survive. Corrections are needed to the residue ranks in the spectrum table, the proposed soliton CDD condition, the imaginary-parameter exercise, and several physical or bibliographic conclusions. The missing second Borel module can be supplied explicitly.

## Summary of verdicts

| Item | Verdict | Lines to change |
|---|---|---|
| 1 `prp-lightest-breather-universal` | Correct as `[Numerical]` for \(r\le4\); "for every \(r\)" stated only as an expectation. Line 18 is wrong for higher sine-Gordon breathers. | 18, 20, after 24 |
| 2 `prp-an-particle-transmission` | Correct, both types; "no S-factor" correct for bulk fusion. "No internal states" and "no bound states" need scoping. | 69, 102 |
| 3 Bajnok–Simon, continuation | One-loop identity exact; assumptions listed; prediction falsifiable; duality conditional. | 12, 134 |
| 4 `prp-dpt-vacuum` | Correct; the `[Theorem]` is proved. | none |
| 5 `prp-dpt-doubled` | Correct; \(Y,\kappa\) solved for. The second module is explicit and unique, and equals ch. 42's for \(h=2\). | 202, 218 |
| 6 Spectrum | Positions and windows correct in all rows. The rank for species \(r\) and the Lee–Yang direction are wrong; \(\Theta\) at \(v=\pi\) and the \(\theta_d=0\) sentence are overstated. | 163, 165, 168, 170, 172, 216 |
| 7 Finite volume, problems | The CDD condition is wrong for \(h>2\), though its conclusion is right and now computed. The Jiang paraphrase overstates. TBA reality needs assumptions. Absence and novelty claims need bounding. | 97, 100, 134, 180, 183, 208, 211, 221, 259, 263 |
| 8 Labels, exercises, bibliography | Labels justified (one could be upgraded); `exr-dpt-imaginary` is wrong as posed; four citation fixes. | 20, 77 (optional), 180, 246, 250 |


## 1. `prp-lightest-breather-universal`

**Verdict: correct as a numerical statement for the tested ranks at generic parameters; the physical breather domain needs stating. The introductory sine-Gordon sentence is wrong for higher breathers.** The all-rank claim is explicitly an expectation at line 50, appropriately so. The `[Numerical]` label is justified; the displayed proof is a description of numerical evidence, not a theorem for every coupling.

**Independent construction.** [independent_breathers.py](../referee_scratch/defects_43/independent_breathers.py), [full output](../referee_scratch/defects_43/breathers_output.txt), builds the tensor generators from

\[
e_j=yE_{j,j+1},\quad f_j=y^{-1}E_{j+1,j},\quad k_j=q^{E_{jj}-E_{j+1,j+1}},
\]

with cyclic indices and the Chevalley-conjugate matrices for \(\bar V\). It applies the iterated coproduct directly to basis words. For \(h\) vectors, the zero-weight basis consists of all \(h!\) permutations; for \(\bar V\otimes V\) or \(V\otimes\bar V\), it is the \(h\) pairs with opposite weights. The common kernel of all \(e_j,f_j\) is computed there; the Cartan equations hold identically on this basis.

For species \(a\), I use the two constituent strings at exactly the shifts on line 26. Every constituent then propagates, rightmost first, through a separately implemented `eq-dt-tv`, retaining every reached defect label. There is no finite defect window. I also independently propagate the two-module singlets through `eq-dt-tv` and `eq-dt-tvbar`.

The scalar implementation is different from the book's. With \(a_k=a_0+kh\omega\), expand each logarithmic Gamma ratio using

\[
\log\Gamma(A+c)\sim(A+c-\tfrac12)\log A-A+\tfrac12\log(2\pi)
+\sum_{n\ge1}\frac{(-1)^{n+1}B_{n+1}(c)}{n(n+1)A^n}.
\]

For the shifts \((z,-z,\omega+z,\omega-z,0,\omega)\) with multiplicities \((1,1,-1,-1,-2,2)\), the logarithmic and \(A^{-1}\) terms cancel. I sum 32 Gamma terms at 38-digit precision and the tail through \(A^{-18}\), using Hurwitz zeta functions. Crossing fixes the positive scalar normalization at \(D=0.193\), and is then tested at another complex rapidity. Increasing the cutoff to 64 changes the checked scalar values by at most \(1.84\times10^{-31}\); the independent crossing residual is at most \(3.65\times10^{-16}\). This is a numerical convergence check, not a rigorous tail bound.

The grid uses \(\omega=\sqrt2,\sqrt5,0.6\), every species, a complex breather rapidity, two defect lattice labels, zero and generic offsets, and additionally \(s=\lambda_1,\lambda_2\) for \(a_3\). Irrational couplings provide tests beyond the rational/root-of-unity examples in the drafting script.

| Rank | Largest singlet residual | Smallest next relative singular value | Largest leakage from singlet × original defect label | Largest error versus \(T_a^{\rm cl}\) |
|---|---:|---:|---:|---:|
| 2 | \(3.36\cdot10^{-16}\) | 0.7071 | \(1.73\cdot10^{-15}\) | \(1.20\cdot10^{-15}\) |
| 3 | \(2.65\cdot10^{-16}\) | 0.5000 | \(2.15\cdot10^{-14}\) | \(1.55\cdot10^{-15}\) |
| 4 | \(8.89\cdot10^{-16}\) | 0.3717 | \(4.43\cdot10^{-15}\) | \(2.67\cdot10^{-15}\) |

The independent two-module construction at \(\omega=\sqrt2\) gives \(T_1^{\rm cl}\) from \(\bar V\otimes V\) and \(T_r^{\rm cl}\) from \(V\otimes\bar V\), to \(4.01\times10^{-16}\). Reversing the shifts, with the species order fixed, gives smallest relative singular values 0.2608, 0.2762, 0.04631 for \(h=3,4,5\): no singlet at these generic points. This does not assert nonexistence at exceptional root-of-unity degeneracies.

**Labels and CZ.** Writing \(x=e^{-D}\), direct simplification of CZ (7.5) gives

\[
e^{-i\pi/3}\frac{\sinh(D/2-i\pi/6)}{\sinh(D/2+i\pi/6)}
=\frac{x+e^{2\pi i/3}}{x+e^{-2\pi i/3}},
\]

and the other factor gives the same expression with \(\gamma=\pi/3\). The independent numerical error is \(4.48\times10^{-16}\). CZ (7.3) places their species 2 at the lower shift and species 1 at the upper shift. With the transposed defect orientation, these are book species 1 and 2, respectively, hence book \(V\otimes\bar V\) and particle 2. The label swap is justified by the Fourier eigenmode convention of chapter 39, as well as by the computed factors. It is a convention dictionary, not an error in CZ. [CZ 2007, §7](https://arxiv.org/abs/0705.1066).

**Corrections ready to paste.**

At **line 18**, replace “the result is a scalar on … independent … of the coupling” with:

```markdown
By @prp-sg-breather-transmission the result is a scalar on $\mathcal W$, independent of the defect state. For $B_1$ it is also independent of the coupling and equals the classical factor (-@eq-shg-linear-transmission). Higher-breather factors depend on the coupling through their fusion shifts.
```

This is a real error, not just missing qualification: at \(D=0.13\), the chapter-37 formula for \(B_2\) gives \(-0.951042862-0.309059015i\) at \(\lambda=2.6\) and \(-0.956504391-0.291717928i\) at \(\lambda=3.1\). Both breathers exist. BCZ themselves make this distinction in (5.7) and (5.10). [BCZ 2005, §5](https://arxiv.org/abs/hep-th/0506169).

After **line 24**, add:

```markdown
A physical lowest breather requires $h\omega>2$, so that its soliton–antisoliton pole has $0<\operatorname{Im}(2i\bar u)<\pi$. Below this threshold the singlet calculation is an algebraic continuation, not a physical bound-state construction. The absence of a singlet with reversed shifts is a statement at generic parameters.
```

Indeed \(\omega=0.6,h=3\), used in both verification grids, lies below this threshold. Its algebraic identity remains valid. At **line 20**, cite `[@gandenberger1995; @gandenberger1998, Sec. 3.6]`: the 1995 paper treats \(a_2\), while the 1998 paper establishes the general \(a_n\) breather–particle comparison in (3.48)–(3.52).

## 2. `prp-an-particle-transmission`

**Verdict: correct for the stated scalar solutions, both types; wrong to infer from a single bulk vacuum alone that a defect has no internal states. The no-bound-state conclusion must be scoped to the proposed minimal amplitudes.**

**Independent algebra.** Set \(x=e^{-D},u=e^{i\gamma_a},v=e^{i\gamma_b}\). For \(a+b<h\),

\[
T_a(D+i\gamma_b)T_b(D-i\gamma_a)
=\frac{x+uv}{x+v/u}\frac{x+v/u}{x+(uv)^{-1}}=T_{a+b}(D).
\]

For \(a+b>h\),

\[
T_a(D+i(\pi-\gamma_b))T_b(D-i(\pi-\gamma_a))
=\frac{x-u/v}{x-(uv)^{-1}}\frac{x-uv}{x-u/v}=T_{a+b-h}(D).
\]

Repeating with \(x=e^D\) and the minus signs in \(T_a^+\) gives the two type-I\(_+\) identities. [supplemental_checks.py](../referee_scratch/defects_43/supplemental_checks.py) simplifies all four rational residuals to **exactly zero**.

Unitarity follows from \(e^{i\gamma_{h-a}}=-e^{-i\gamma_a}\): \(T_a(D)T_{h-a}(D+i\pi)=1\). Defining \(\bar T_a(\theta)=T_{h-a}(i\pi-\theta)\), with the same fixed defect parameter, then gives crossing and \(T_a(\theta)\bar T_a(-\theta)=1\). For real \(D\), numerator and denominator are complex conjugates, so \(|T_a|=1\). Consequently \(\bar T_a(\theta)=T_a(-\theta)^*\). The same reasoning works for \(T^+\).

The only poles modulo \(2\pi i\) are \(D=-i(\pi-\gamma_a)\) and \(D=-i\gamma_a\) for the two types. There are none in \(0<\operatorname{Im}\theta<\pi/2\) for real \(\theta_d\), or even in the larger strip up to \(\pi\). The conjugate-direction factors have the corresponding lower-half-plane representatives. The quoted bulk fusing poles follow from the endpoint blocks in `eq-an-smatrix`, at interacting coupling \(0<B<2\); the poles disappear at the free endpoint \(B=0\).

**Independent classical calculation.** [independent_particles.py](../referee_scratch/defects_43/independent_particles.py), [output](../referee_scratch/defects_43/particles_output.txt), constructs the affine roots in \(\mathbb R^h\), the rotation \(\pi\), and \(D=2(1-\pi)^{-1}\) on the sum-zero subspace. It differentiates the exponential potential at zero by summing outer products of exponent vectors. Substitution of \(\phi=(1+R)\rho_a,\psi=T\rho_a\) into both full vector sewing conditions determines \(R,T\), without importing the chapter-39 reduction. Both kinetic signs and the reciprocal potential parameter for I\(_+\) are included.

```text
scalar checks h=2..8, both types, fusion cases=224
unitarity: 8.006e-16
Hermitian analyticity: 7.301e-16
phase: 2.220e-16
bootstrap: 1.977e-15
block identities: 5.118e-16
bulk-S CDD bootstrap: 5.447e-15
independent Hessian/vector sewing h=2..8:
  reflection=1.555e-15, transmission-factor error=1.633e-15
```

**No S-factor.** Correct for the **bulk-particle fusion** bootstrap being discussed. If \(C:V_c\to V_a\otimes V_b\) is the fusion map, then \((T_a\otimes1)(1\otimes T_b)(C\otimes1)=(1\otimes C)T_c\). The constituents cross successively without exchanging their spatial order, so no bulk S-matrix is inserted. A different equation, excitation of the **defect** by absorbing a particle, does contain an S-factor; the chapter uses it correctly in the last part of line 172.

**Paste-ready scope corrections.** At **line 69**, replace the sentence beginning “At real coupling there is a single vacuum” by:

```markdown
At real coupling there is a single bulk vacuum and no soliton topological charge. We consider transmission on the ground-state defect, with no additional internal multiplicity. For this scalar ansatz and the diagonal bulk S-matrix, the triangle relations are automatic.
```

At **line 102**, replace the first two sentences by:

```markdown
**No particle bound-state poles in the proposed ground-state solution.** For real $\theta_d$, the minimal factors of @prp-an-particle-transmission have no pole in the physical strip, so this bootstrap solution supplies no particle bound state. This is not a classification of all possible defect degrees of freedom or of nonminimal solutions.
```

## 3. Bajnok–Simon, continuation, and duality

**Verdict: correct for the one-loop identification and conditional parameter/duality statements. The prediction is falsifiable; a failed prediction would not by itself identify the source of failure as the Penati–Polvara anomaly.**

**Independent derivation.** With \(x=e^{\theta+\mu}\),

\[
T^{\rm cl}=\frac{1+ix}{1-ix},\qquad
\partial_\mu T^{\rm cl}=\frac{2ix}{(1-ix)^2}
=-\frac1{1-i\sinh(\theta+\mu)}.
\]

The rational identity simplifies to zero in `independent_particles.py`. Write \(t=b^2/(8\pi)\). Since \(B/2=t/(1+t)\), \(\mu_q=\mu(1-B/2)=\mu-\mu t+O(t^2)\), whence

\[
T(\theta,\mu_q)=T(\theta,\mu)+\frac{\mu b^2}{8\pi}\frac1{1-i\sinh(\theta+\mu)}+O(b^4).
\]

This is exactly the source's §4.2 result, with \(B_{\rm BS}\vartheta=-\theta_d\), \(\mu=-\theta_d^{\rm cl}\), and \(B_{\rm BS}=B/2\). Equations (16), (17), and (23) agree with that dictionary. The book uses \(b=\beta_{\rm SG}=\sqrt2\beta\), so the denominator \(8\pi+b^2\) is also correct. The source checks absence of reflection at this order and explicitly leaves the independent mass-scale renormalization for a further calculation. [Bajnok–Simon](https://arxiv.org/abs/0712.4292).

The assumptions at lines 108–111 explicitly cover minimality, the vacuum, parameter renormalization and quantum integrability. The conjecture repeats the relevant ground-state/minimal restrictions. The higher-rank one-loop prediction is the testable equality

\[
\frac{\delta T_a(\theta)}{\partial_{\theta_d}T_a^{\rm cl}(\theta)}
=\delta\theta_d\quad\text{for every }a,\theta,
\]

together with zero reflected and species-changing amplitudes. This follows by a Taylor expansion of the conjectured family; it is not claimed as an already performed perturbative calculation.

For duality, fixing the quantum parameter gives

\[
\frac{2-B}{2}\theta_d^{\rm cl}=\frac B2\theta_d'^{\rm cl},\qquad
\theta_d'^{\rm cl}=\frac{2-B}{B}\theta_d^{\rm cl}.
\]

Thus line 150 is correct conditional on the parameter conjecture and the identification with the quantum Lagrangian. It establishes a scattering-parameter map, not a separate nonperturbative proof of duality of defect QFTs.

The Penati–Polvara description agrees with §5.2: explicit ranks through 5, a spin-3 current, and allowed renormalizations of the potential, derivative matrix and improvement terms. Their spin-4 sinh-Gordon current has a spin-3 charge. Their conclusions leave the integrability question open; the book correctly does not treat the anomaly as a theorem excluding every quantum completion. [Penati–Polvara](https://arxiv.org/abs/1902.10690).

**Paste-ready correction at line 134**, last sentence:

```markdown
A species-dependent parameter shift, reflection, or species conversion would falsify this prediction of the continuation conjecture. Relating such a discrepancy specifically to the current anomaly of Penati and Polvara would require a calculation in a common renormalization scheme, also checking the minimality and parameter assumptions.
```

For the opening overview at **line 12**, replace “the only quantum effect left in particle transmission is” by “under the continuation conjecture, the only quantum effect in these ground-state particle factors is”. The body already makes this qualification, but the overview should too.

## 4. `prp-dpt-vacuum`

**Verdict: correct; the `[Theorem]` label and the proof in the text are justified for real root-space fields, real nonzero coupling, and positive mass scale.**

**Independent proof.** Set

\[
u_j=\frac b2\alpha_j\cdot(D^T\phi+D\psi),\qquad
v_j=\frac b2\alpha_j\cdot D(\phi-\psi).
\]

Both vectors have zero sum. The arithmetic–geometric mean inequality gives \(\sum_j e^{u_j}\ge h\), since \(\prod_j e^{u_j}=1\), and equality requires every \(u_j=0\); the same holds for \(v\). Positive \(\sigma\) makes both contributions to the book-normalized potential nonnegative. Equality requires \(u=v=0\). The affine roots span the root space, so \(v=0\) implies \(D(\phi-\psi)=0\). The rotation has no eigenvalue one on that space, and therefore \(D=2(1-\pi)^{-1}\) is invertible. Hence \(\phi=\psi\). Now \(u=0\) gives \((D+D^T)\phi=2\phi=0\). This proves uniqueness for every rank, including rank one. Interchanging \(\sigma\) and \(\sigma^{-1}\) proves the other type.

For \(\sigma<0\), choose a nonzero real root-space vector \(v\), take \(\phi=\psi=tv\), and choose its sign so some \(b\alpha_j\cdot v>0\). Then every \(g_j=1\), whereas

\[
\mathcal B_d(tv,tv)=\frac{m\sigma}{b^2}
\left(\sum_j e^{bt\alpha_j\cdot v}-h\right)\longrightarrow-\infty.
\]

For I\(_+\), replace the negative coefficient by \(\sigma^{-1}\). This is an explicit unbounded direction, not an inference from a finite random sample.

The extension to static field energy used elsewhere in the chapter also has an independent justification. Extend the boundary values \(tv\) to zero using linear layers of width \(\ell=e^{-At/2}\), where \(A=\max_j b\alpha_j\cdot v>0\). Their positive bulk gradient and potential energies are \(O(t^2e^{At/2})\), while the negative defect term is of order \(-e^{At}\). The total energy therefore tends to \(-\infty\). For positive \(\sigma\), all bulk static terms and the defect term are nonnegative and vanish only for the constant zero fields.

No correction to this proposition is needed. At imaginary coupling the exponents are not real, so the real convexity/AM–GM inequality is unavailable; periodic sine-Gordon vacua do not contradict it.

## 5. `prp-dpt-doubled`

**Verdict: correct as a `[Numerical]` statement, with the right logical scope: \(Y\) and \(\kappa\) are solved for, not fixed. The open end of the proof, at line 202, can be closed. The \(e\)-type solution is explicit, the module relations select it uniquely, and for \(h=2\) it is chapter 42's second module (`eq-dbd-second-module`). The two Borel halves of the \(a_r\) defect also do not combine into one coideal, by the argument of `prp-diagonal-borel` (4).**

**Scope of the book's evidence.** Check (f) (`textbook_code/defect_particle_transmission.py`, lines 407–455) has one unknown coefficient for \(Y\) and one for \(\kappa\) at every window label. These are \(\theta\)-independent and shared between \(V\) and \(\bar V\) and between rapidities; only the weights are imposed. So line 202 describes the computation correctly. There are two gaps. Equations are imposed only on interior labels of a finite window, and for \(a_3\) the script omits the \(f\)-type cases \(\psi=\pi^{\pm1}\), although the proposition asserts their failure for \(a_3\). The independent computation below closes both gaps.

**Independent computation 1** ([independent_doubled.py](../referee_scratch/defects_43/independent_doubled.py), [output](../referee_scratch/defects_43/doubled_output.txt)). Every input \(v_i\otimes|N\rangle\) is propagated exactly through `eq-dt-tv` or `eq-dt-tvbar`, and every reached output equation is kept. There is no window truncation. Both \(Y\) and \(\kappa\) are unknown at every reached label, common to \(\theta=-0.32,\ 0.41,\ 0.82+0.11i\) and to \(V,\bar V\), with \(\omega=\sqrt2\), \(\theta_d=0.27\) and a generic offset. All nodes are solved, for \(h=3\) and \(h=4\). Relative least-squares residuals (ranges over nodes):

| \(h\) | \(f\), \(\psi=\mathrm{id}\) | \(e\), \(\psi=\pi\) | \(e\), id | \(e\), \(\pi^{-1}\) | \(f\), \(\pi\) | \(f\), \(\pi^{-1}\) |
|---|---:|---:|---:|---:|---:|---:|
| 3 | \(\le1.7\cdot10^{-13}\) | \(\le1.6\cdot10^{-13}\) | 0.918–0.928 | 0.918–0.927 | 0.826–0.827 | 0.845–0.846 |
| 4 | \(\le1.8\cdot10^{-12}\) | \(\le5.0\cdot10^{-13}\) | 0.916–0.929 | 0.915–0.929 | 0.693–0.703 | 0.743–0.758 |

Here \(\pi(e_j)=e_{j-1}\). The smallest failing residual, 0.693, is consistent with the book's "at least 0.67". The surviving cases are iterative (`lsmr`) solutions; substituting the closed forms below instead leaves at most \(6.5\times10^{-16}\).

**Independent computation 2: the \(e\)-type solution by hand** (second referee). Take \(\psi=\pi\), \(\pi(e_j)=yE_{j,j+1}\) on \(V\), \(Y|N\rangle=Y_N|N+\alpha_j\rangle\), \(\kappa|N\rangle=\kappa_N|N+\alpha_j-\alpha_{j-1}\rangle\), \(n=N+s\), and \(X=y/t\) with \(t=e^{\omega(\theta_d+i\pi/2)}\). Apply both sides of the \(e\)-type equation at line 191 to \(v_i\otimes|N\rangle\) with `eq-dt-tv`, and sort by powers of \(y\):

- \(y^0\): \(Y_Nq^{-n_i}\) on both sides, identically.
- \(y^2\): only for \(i=j+1\), output \(|N+e_j-e_{j-1}\rangle\otimes v_{j-1}\). The left side gives \(yX\), from \(\pi(e_j)\) followed by the converting entry of \(M_V\). The right side gives \(Xy\,\kappa_{N-\alpha_j}\), from the converting entry followed by \(\kappa\otimes\pi(e_{j-1})\). Hence \(\kappa\equiv1\).
- \(y^1\), \(i=j+1\): \(Y_{N-\alpha_j}-q^{-1}Y_N=t\,q^{-n_j}\).
- \(y^1\), \(i=j\): \(qY_N-Y_{N-\alpha_{j-1}}=t\,q^{-n_j}\).
- \(y^1\), \(i\ne j,j+1\): \(Y_{N+e_i-e_{i-1}}=Y_N\), invariance under the other simple shifts.

The general solution is \(Y_N=C\,q^{-n_j}+c_j\,q^{n_j}\), with \(C=t/(q-q^{-1})\) and one free constant \(c_j\) per node. This is the "free constant" of line 202. Now impose the relations of \(U_q(\mathfrak b^+)\). Write \(\tilde\rho(k_i)=\kappa_i\), the shift by \(w_i=\alpha_i-\alpha_{i-1}\). The \(j\)-th component of \(w_i\) is \(2\delta_{ij}-\delta_{j,i+1}-\delta_{j,i-1}=a_{ij}\), the affine Cartan matrix for \(h\ge3\); for \(h=2\) it is \(2\) or \(-2\), again \(a_{ij}\). Therefore

\[
\tilde\rho(k_i)\,Y_j\,\tilde\rho(k_i)^{-1}=q^{a_{ij}}Y_j
\iff c_j\,q^{n_j}\big(q^{-a_{ij}}-q^{a_{ij}}\big)=0 ,
\]

and \(i=j\) forces \(c_j=0\) for generic \(q\). With \(c_j=0\), the \(\tilde\rho(k_i)\) commute and have product 1, since \(\sum_iw_i=0\). Moreover \(Y_{j+1}Y_j=q\,Y_jY_{j+1}\), and non-adjacent \(Y\) commute. By the lemma after `eq-dt-weyl` (\(c=q^{\pm1}\) gives the quadratic \(q\)-Serre relations), the \(q\)-Serre relations hold for \(h\ge3\). For \(h=2\), \(Y_1Y_0=Y_0Y_1=C^2q\) is central, so the cubic relations hold too. Thus \(\tilde\rho\) is a \(U_q(\mathfrak b^+)\)-module, and it is the unique choice among the \(\theta\)-independent solutions.

The script [second_module_checks.py](../referee_scratch/defects_43/second_module_checks.py) ([output](../referee_scratch/defects_43/second_module_output.txt)) checks (C): with \(c_j=0.3\) the relation fails by 0.23–0.25 for \(h=3,4,5\), and with \(c_j=0\) it holds to \(6.3\times10^{-16}\). Codex's `module_relations` checks all Cartan and Serre relations of \(\tilde\rho\) for \(h=3,\dots,6\), to \(5.5\times10^{-16}\). The same \(\tilde\rho\) jointly solves the \(V\) and \(\bar V\) equations to \(6.5\times10^{-16}\) (\(h=3,4\), all nodes, three rapidities). The \(\bar V\) part is numerical only.

**Consistency with chapter 42.** For \(h=2\), take \(N=(n,-n)\) and \(s=(\frac s2,-\frac s2)\), and set \(\tau=-ie^{i\pi\omega/2}\), so that \(\tau^2=q^{-1}\). In the basis \(d_n|n\rangle\), \(d_n=\tau^{-n}q^{n(n-1)/2+ns/2}\), the matrix \(q^{s/2}M_V\) of `eq-dt-tv` becomes the matrix part of `eq-sg-transmission-matrix`. Its invariant \(C_nB_{n+1}=-Y^2q^{1+s}\) agrees, and the entries agree to \(4.1\times10^{-15}\) (check (A)). In the same basis, \(\tilde\rho\) becomes exactly `eq-dbd-second-module`:
\(\tilde\rho(e_1)\to\frac{iy_d}{q^2-1}q^{-2n-s}\), \(\tilde\rho(k_1)\to q^{-2n-2-s}\) with shift \(+2\), \(\tilde\rho(e_0)\to\frac{iy_d}{q^2-1}q^{2n+s}\), and \(\tilde\rho(k_0)\to q^{2n-2+s}\) with shift \(-2\). The coefficient computation is \(C\tau=iy_d/(q^2-1)\), and the entries agree to \(5.5\times10^{-15}\) (check (B)). For \(h=2\), \(\pi=\tau\). So chapter 42's part (2) is the rank-one member of the same family, and the \(f\)-type identity case is its part (1), as line 202 says. In fact the \(f\)-type solution is \(Y=\rho(f_j)\), \(\kappa=\rho(k_j)^{-1}\), so that the equation is `eq-borel-intertwiner` itself.

**The two halves do not form one coideal (\(a_r\)).** In \(U\otimes U^{\rm cop}\),
\[
\Delta(k_j)\big(e_j\otimes1+k_j\otimes e_{j-1}\big)\Delta(k_j)^{-1}=q^2\,e_j\otimes1+q^{-1}\,k_j\otimes e_{j-1}\qquad(h\ge3),
\]
so for \(q^3\ne1\) any subalgebra containing both halves contains \(e_j\otimes1\). That element acts on \(V\otimes\mathcal W\) as \(\pi(e_j)\otimes1+\pi(k_j)\otimes Z\) and on the reflected copy as \(1\otimes Z\). The \(y^2\) term of \(\mathbb T(\pi(e_j)\otimes1)\) on \(v_{j+1}\otimes|N\rangle\) is \(yX\,|N+e_j-e_{j-1}\rangle\otimes v_{j-1}\), and no other term reaches order \(y^2\). So no \(Z\) exists. Problem 12 currently cites the sine-Gordon proposition for this; the \(a_r\) case follows by the same argument.

**Corrections ready to paste.**

At **line 202**, replace the last two sentences ("As for sine-Gordon theory, the $e$-type solution leaves free constants. We have not checked …") with:

```markdown
The $f$-type solution with $\psi=\mathrm{id}$ is $Y=\rho(f_j)$, $\kappa=\rho(k_j)^{-1}$, so that case is the intertwining property (-@eq-borel-intertwiner) itself. The $e$-type solution with $\psi=\pi$ is explicit. On $V$ the terms of order $y^2$ force $\kappa$ to be the plain shift by $\alpha_j-\alpha_{j-1}$, and the terms of order $y$ leave $Y=\tilde\rho(e_j)$ plus $c_j\,q^{(N+s)_j}$ times the shift by $\alpha_j$, one free constant $c_j$ for each node, where
$$
\tilde\rho(e_j)\lvert N\rangle=\frac{e^{\omega(\theta_d+i\pi/2)}}{q-q^{-1}}\,q^{-(N+s)_j}\,\lvert N+\alpha_j\rangle,\qquad
\tilde\rho(k_j)\lvert N\rangle=\lvert N+\alpha_j-\alpha_{j-1}\rangle .
$$
The $j$-th component of $\alpha_i-\alpha_{i-1}$ is the Cartan entry $a_{ij}$, so $\tilde\rho(k_i)\tilde\rho(e_j)\tilde\rho(k_i)^{-1}=q^{a_{ij}}\tilde\rho(e_j)$, and the same relation for $Y$ with $i=j$ forces $c_j=0$. The $\tilde\rho(k_i)$ commute with product 1; adjacent $\tilde\rho(e_j)$ commute up to $q^{\pm1}$ and the others commute, so the $q$-Serre relations hold as for (-@eq-dt-weyl). Thus $\tilde\rho$ is a module of $U_q(\mathfrak b^+)$. It also solves the equations for $\bar V$ (checked numerically for $a_2$ and $a_3$). For $h=2$, in the basis $\tau^{-n}q^{n(n-1)/2+ns/2}\lvert n\rangle$ with $\tau=-ie^{i\pi\omega/2}$, and with $\mathbb T$ multiplied by $q^{s/2}$, it is (-@eq-dbd-second-module).
```

At **line 218** (problem 12), replace the sentence "Still open are the module relations of the second structure $\tilde\rho$ for $a_r$, the type-I$_+$ case, and the algebra generated by the two Borel halves, which do not combine into one coideal (@prp-diagonal-borel)." with:

```markdown
Still open are the type-I$_+$ case and the algebra generated by the two Borel halves. They do not combine into one coideal: for $q^3\ne1$, conjugation by $\Delta(k_j)$ separates $e_j\otimes1$ from $k_j\otimes e_{j-1}$, and $e_j\otimes1$ cannot act, as in @prp-diagonal-borel (4).
```

The proposition's `[Numerical]` label should stay. The \(V\) half of the \(e\)-type statement is now a derivation, but the \(\bar V\) half and the failures of the other twists are numerical.

## 6. `tbl-dpt-spectrum` and `sec-defect-spectrum`

**Verdict: the pole positions and coupling windows in all nine rows are correct. Two statements are wrong: the residue rank "of species \(a\)" in row 6 (and problem 10), and the table caption's direction of the Lee–Yang poles. Two statements are overstated: the first sentence of "Reading the table", and "\(\Theta\) broken" for imaginary \(\theta_d\), which fails at \(v\in\pi\mathbb Z\). At \(v=\pi\) the particle pole of row 9 sits on the \(\sigma=-1\) defect, whose classical energy is unbounded below.** The brief's example is handled correctly. Row 6 uses the general window \(\max(0,a+2p-\frac h2)<\frac1\omega<a+2p\), which for \(a_2\), species 1, \(p=0\) gives all \(\omega>1\), as in `prp-a2-unstable` (1). A stale window is not repeated.

**Residue ranks (derivation, second referee).** At the species-\(r\) pole of \(\bar g\), \(\theta_r=\theta_d+\frac{i\pi}2-i(\frac\pi h+\frac\pi{h\omega})\), we have \(X=e^{-i\pi/h-i\pi\omega/h}\), so \(X^h=-e^{-i\pi\omega}=q\), not the value \((-q^{-1})^rX^h=-1\) of the vector pole. On the sector of total charge \(Q\) (\(\sum_iQ_i=-1\)), incoming \(\bar v_i\otimes|Q+e_i\rangle\), `eq-dt-tvbar` gives

\[
(M_{\bar V})_{ki}=q\,u_k\,w_i\times\begin{cases}1,&k\ge i,\\ X^hq^{-1},&k<i,\end{cases}\qquad
u_k=X^kq^{P_k},\quad w_i=X^{-i}q^{-P_{i-1}},\quad P_k=\sum_{m\le k}(Q+s)_m ,
\]

using \(P_h=-1\) for the entries that wrap around. So \(M_{\bar V}=q\,u\,w^{\mathsf T}\) has rank 1 exactly when \(X^h=q\). By contrast, \((M_V)_Q\) is cyclic bidiagonal with non-zero entries, so where \(\det(M_V)_Q=0\) it has rank exactly \(h-1\). The pole of \(\bar g\) is not cancelled. [supplemental_output.txt](../referee_scratch/defects_43/supplemental_output.txt) gives \(\epsilon\,\bar g(\theta_r+\epsilon)\to0.0617-0.2446i\), \(0.0691-0.0870i\) and \(0.0453-0.0342i\) for \(h=3,4,5\), stable between \(\epsilon=10^{-5}\) and \(5\cdot10^{-6}\). Hence **the residue of \(\mathbb T_{\bar V}\) at the species-\(r\) pole has rank 1, and its kernel has dimension \(h-1\).** Independently, Codex's [particles_output.txt](../referee_scratch/defects_43/particles_output.txt) gives normalized singular values \([1,1,8\cdot10^{-17}]\) for \(h=3\), \(a=1\) (rank 2) and \([1,10^{-16},3\cdot10^{-17}]\) for \(a=2\) (rank 1). For \(h=4,5\) it gives rank \(h-1\) for \(a=1\) and rank 1 for \(a=r\). Chapter 41 (lines 200, 215) claims rank \(h-1\) only for \(\mathbb T_V\), and correctly; chapter 43 extends it to every species.

**Lee–Yang direction** ([supplemental_checks.py](../referee_scratch/defects_43/supplemental_checks.py)). At \(b=0.5\), near \(\theta=i\pi/4\): \(\epsilon\,\mathbb T_b\to0\) and \(\epsilon\,\bar{\mathbb T}_b\to-3.4641i\). The pole is in \(\bar{\mathbb T}_b\), as `prp-ly-transmission` and Bajnok–Simon §6.1 say (their \(T_+=[5-b][-5-b]\) has the pole at \(i\pi(b+1)/6\) for \(b\in[-1,2]\)). The caption at line 170 says that all rows list poles of \(\mathbb T\).

**Imaginary \(\theta_d\) (derivation).** For \(\theta_d=iv\), \(|\mathbb T_a(\theta)|^2=\frac{1+e^{-2\theta}+2e^{-\theta}\cos(v-\gamma_a)}{1+e^{-2\theta}+2e^{-\theta}\cos(v+\gamma_a)}\) for real \(\theta\), which equals 1 iff \(\sin v\sin\gamma_a=0\), that is \(v\in\pi\mathbb Z\). The window \(\pi-\gamma_a<v<\frac{3\pi}2-\gamma_a\) contains \(v=\pi\) exactly when \(a<\frac h2\). There \(\sigma=e^{-i\pi}=-1\) is real, so \(\Theta\) (complex conjugation at real coupling, `prp-defect-theta`) is unbroken and \(\mathbb T_a\) is a phase. But by `prp-dpt-vacuum` this defect's classical energy is unbounded below. Codex's example: \(h=3\), \(a=1\), \(v=\pi\) gives \(u=\pi/3\) and \(|\mathbb T(0.4)|=1\) to machine precision. The same exception applies to the sine-Gordon row 4, where \(v=\pi\) is in the window for \(\lambda<1\).

**Row-by-row audit.**

| Row | Claim | Verdict and evidence |
|---|---|---|
| Lee–Yang | \(\lvert1^+\rangle\), \(-1<b<2\); \(\lvert2^+\rangle\), \(1<b<2\) | Ranges correct (BS §6.1, Fig. 9). The poles are in \(\bar{\mathbb T}_b\), not \(\mathbb T_b\); the caption must say so. |
| shG, real \(\theta_d\) | none | Correct: \(\gamma=\frac\pi2\), pole at \(\theta_d-\frac{i\pi}2\), outside the strip. |
| SG, real \(\theta_d\) | resonance \(\theta_d-\frac{i\pi}{2\lambda}\); odd defect \(\theta_d-i\pi\) | Correct (ch. 37, BCZ §4); \(\theta_d-i\pi\) is from the image, not the kernel. |
| SG, \(\theta_d=iv\) | strip for \(\frac\pi{2\lambda}<v<\frac\pi2+\frac\pi{2\lambda}\); \(\Theta\) broken | Window correct. \(\Theta\) is broken only for \(v\notin\pi\mathbb Z\); \(v=\pi\) is in the window for \(\lambda<1\). |
| SG type II | \(\theta_d\pm\tau-\frac{i\pi}{2\lambda}\) | Correct as quoted, with ch. 37's correction of CZ2010 (3.17). |
| \(a_r\) solitons | positions, windows, \(\theta_d=0\) candidates; "residue of rank \(h-1\)" | Positions and windows are correct. At \(\theta_d=0\), \(\operatorname{Im}\theta_a>0\) iff \(\frac1\omega<a+2p-\frac h2\), and \(<\frac\pi2\) always. **Rank wrong for \(a=r\)** (rank 1); intermediate species not computed. |
| \(a_2^{(2)}\) | windows \((\frac\pi3,\frac{4\pi}3)\), \((-\frac\pi3,\frac{2\pi}3)\) | Correct: poles at \(\eta+i(\pm\frac\pi6-\frac\xi2)\) with \(-\frac\pi2<\operatorname{Im}<0\). Consistent with ch. 41 line 189. |
| \(a_r\) particles, real \(\theta_d\) | none | Correct for the minimal factors (§2). |
| \(a_r\) particles, \(\theta_d=iv\) | pole \(i(v-\pi+\gamma_a)\), window; \(\Theta\) broken | Position, window and residue \(e^{2i\gamma_a}-1\) are correct. \(\Theta\) is broken only for \(v\ne\pi\); at \(v=\pi\) (\(a<\frac h2\)) the vacuum is unstable. |

**Corrections ready to paste.**

**Line 165**, third column of the \(a_r\)-soliton row, replace "residue of rank $h-1$, whose kernel carries the class-$(p+1)$ matrix" by:

```markdown
at the species-1 pole the residue of $\mathbb T_V$ has rank $h-1$, and its kernel carries the class-$(p+1)$ matrix; at the species-$r$ pole the residue of $\mathbb T_{\bar V}$ has rank 1; other species not computed
```

**Line 163** and **line 168**, replace "$\Theta$ broken" by "$\Theta$ broken unless $v\in\pi\mathbb Z$". In line 168 append "; at $v=\pi$ ($a<\frac h2$) $\sigma=-1$ and the vacuum is unstable (@prp-dpt-vacuum)".

**Line 170**, caption, replace the second sentence by:

```markdown
Rapidities are pole positions of $\mathbb T$ for a particle or soliton arriving from the left, except in the Lee–Yang row, whose poles are those of $\bar{\mathbb T}_b$ (@prp-ly-transmission); $\bar{\mathbb T}$ has the mirror poles.
```

**Line 172**, replace the first sentence ("At real $\theta_d$ no affine Toda defect has a stable excited state.") by:

```markdown
For real $\theta_d\ne0$ none of these amplitudes has a stable excited state.
```

In the same paragraph, after "the energy of the bound state is $m_a\cos(v-\pi+\gamma_a)$", insert:

```markdown
At $v=\pi$, inside the window for $a<\frac h2$, $\sigma=-1$ is real and $\mathbb T_a$ is a phase, but the classical energy is unbounded below (@prp-dpt-vacuum), so the pole alone does not establish a stable excited defect.
```

**Line 216** (problem 10), replace the first sentence by:

```markdown
At the species-1 pole the residue of $\mathbb T_V$ has rank $h-1$; at the species-$r$ pole the residue of $\mathbb T_{\bar V}$ has rank 1, with an $(h-1)$-dimensional kernel; the intermediate species have not been computed.
```

In the same problem, replace "For imaginary $\theta_d$, stable excited defects exist for sine-Gordon solitons and for $a_r$ particles (@tbl-dpt-spectrum)." by:

```markdown
For imaginary $\theta_d$ the sine-Gordon soliton factors and the $a_r$ particle factors have poles in the physical strip (@tbl-dpt-spectrum); whether they describe stable excited defects is open.
```

## 7. Finite volume, open problems, novelty, Jiang, and the CDD claim

**Verdict: the finite-volume section and the list of open problems are accurate in substance, and the Rossi–Weston/BLZ/Weston credit in problem 11 is now correct. Corrections are needed in five places. (i) The soliton CDD condition at line 100 is wrong for \(h>2\), although the conclusion drawn from it is right; it is now computed, not only argued. (ii) The Jiang paraphrase overstates. (iii) The TBA reality sentence needs its assumptions. (iv) Absolute absence claims should be bounded. (v) The novelty list omits Corrigan–Zambon's (7.11).**

### 7a. The soliton CDD factor (line 100)

**Derivation (checked by both referees).** Multiply \(g\) by \(f(\theta-\theta_d)\). By `eq-dt-bootstrap`, \(\bar g\) is multiplied by \(\prod_{m=0}^{r-1}f(\theta+i\pi(2m-r+1)/h)\). The second crossing relation of chapter 41, \(\bar g(\theta)g(\theta-i\pi)c_{\rm II}=1\), then requires

\[
\prod_{k=0}^{h-1}f\Big(\theta-i\pi+\frac{2\pi ik}h\Big)=1 ,
\]

the product over a string of \(h\) rapidities spaced by \(\frac{2\pi}h\), the trivial-representation string of line 28. The first crossing relation is the same condition shifted by \(-i\pi+\frac{2\pi i}h\). For \(h=2\) it is \(f(\theta)f(\theta+i\pi)=1\), the condition printed at line 100. For \(h>2\) the printed condition is neither necessary nor sufficient:

- \(f=(\sinh D+0.3i)/(\sinh D-0.3i)\) satisfies the printed condition to \(2.2\times10^{-16}\), but at \(h=3\) the actual condition fails by 0.903 (`independent_particles.py`).
- \(f=\mathbb T^{\rm cl}_1\) satisfies the actual condition for every \(h\). With \(x=e^{-D}\) and \(u=e^{i\pi/h}\), the shifted arguments give \(e^{-D}\to-x\,e^{-2\pi ik/h}\), and \(\prod_k(u-x\,e^{-2\pi ik/h})/(u^{-1}-x\,e^{-2\pi ik/h})=(u^h-x^h)/(u^{-h}-x^h)=1\), since \(u^{\pm h}=-1\). But \(\mathbb T^{\rm cl}_1(D)\mathbb T^{\rm cl}_1(D+i\pi)=(u^2-x^2)/(u^{-2}-x^2)\ne1\) for \(h>2\).

**Coupling dependence, computed.** The breather string is the trivial string with its halves pulled together by \(\frac{2\pi}{h\omega}\), so the CDD product over it is generically not 1. Two admissible families:

- \(f=e^{ic\sinh D}\) satisfies the condition because \(\sum_k e^{\pm 2\pi ik/h}=0\). On species \(a\) it fuses to \(\exp\big(ic\frac{\sin(\pi a/h)}{\sin(\pi/h)}\sinh D\big)\). On \(B^{(a)}_1\) it gives \(\exp\big(2ic\frac{\sin(\pi a/h)}{\sin(\pi/h)}\sin\frac\pi{h\omega}\sinh D\big)\), explicitly \(\omega\)-dependent. At \(h=4\), \(a=2\), \(c=0.19\), \(D=0.41\): \(0.99287+0.11917i\) at \(\omega=\sqrt2\) and \(0.99696+0.07787i\) at \(\omega=\sqrt5\). The closed formula agrees to \(2.3\times10^{-16}\).
- The bounded meromorphic \(f=\mathbb T^{\rm cl}_1\), at the same point: \(-0.21999-0.97550i\) at \(\omega=\sqrt2\) and \(0.20044-0.97971i\) at \(\omega=\sqrt5\). The crossing residual is \(<7\times10^{-16}\).

So the claim "in general that product depends on the coupling" is **correct**, and minimality is a genuine assumption. BCZ already state the sine-Gordon analogue, just after their (5.7); the point is not new, and the chapter does not claim it is.

**Paste-ready, line 100**, replace "A soliton defect CDD factor $f(\theta-\theta_d)$ with $f(\theta)f(\theta+i\pi)=1$ multiplies the breather factor by a product of $h$ shifted copies of $f$. In general that product depends on the coupling, and @prp-lightest-breather-universal would fail." by:

```markdown
A soliton defect CDD factor $f(\theta-\theta_d)$ multiplies $g$, and by (-@eq-dt-bootstrap) it multiplies $\bar g$ by $r$ shifted copies of $f$. The crossing relations of @sec-transmission-unitarity-crossing then require
$$
\prod_{k=0}^{h-1}f\Big(\theta-i\pi+\frac{2\pi ik}h\Big)=1,
$$
the product over the string of rapidities that carries $\Lambda^hV$; for $h=2$ this is $f(\theta)f(\theta+i\pi)=1$. The breather factor is multiplied by the product of $f$ over the same string with its halves pulled together by $\frac{2\pi}{h\omega}$, and that product depends on the coupling. For example $f=\mathbb T^{\rm cl}_1$ satisfies the condition for every $h$; for $a_3$, species 2, at $\theta-\theta_d=0.41$, it multiplies the breather factor by $-0.2200-0.9755i$ at $\omega=\sqrt2$ and by $0.2004-0.9797i$ at $\omega=\sqrt5$. With such a factor @prp-lightest-breather-universal would fail.
```

At **line 97**, replace "a standing particle of species $e$ at rapidity $\theta_0$ [@bajnok_george2006]" by "a defect that scatters like a particle of species $e$ at rapidity $\theta_0$, a standing particle for $\theta_0=0$ [@bajnok_george2006]". Bajnok and George (hep-th/0404199, after (18)): "If u = 0 the solution describes a standing particle."

### 7b. Finite volume

**TBA (line 180).** For real \(\theta_d\), set \(F_a(\theta)=\mathbb T_{\bar a}(\frac{i\pi}2+\theta)\) and \(x=e^{\theta_d-\theta}>0\). By hand, \(F_a=(ix+e^{-i\gamma_a})/(ix+e^{i\gamma_a})\), so \(F_{\bar a}=\overline{F_a}\) and

\[
|F_a|^2=\frac{1+x^2-2x\sin\gamma_a}{1+x^2+2x\sin\gamma_a}<1 .
\]

Codex's check gives \(4.0\times10^{-16}\) for \(h=2,\dots,8\). For the self-conjugate species of even \(h\), \(F_{h/2}=-\tanh\frac{\theta-\theta_d}2\) is real and changes sign (\(-0.2260\) at \(\theta=0.73\), \(\theta_d=0.27\)). So "driving term \(-\log\)" needs a branch. The form that Bajnok and Simon use in §5.1 (after (25)), \(\log(1+F\,e^{-\epsilon})\), avoids the issue. Conjugation together with \(a\leftrightarrow\bar a\) is a symmetry of the system, and reality of \(E_0\) follows for a symmetric solution, which is what one gets at large \(L\). Bajnok and Simon note that their DTBA derivation "generalizes to any diagonal scattering theory" (§7). So the \(a_r\) DTBA is written by that recipe; what is missing is its analysis. The quoted energy \(m\cosh\theta_d/(2\sin\frac{\pi B}2)\) is BS (18), with \(B_{\rm BS}\vartheta=-\theta_d\) and \(B_{\rm BS}=\frac B2\). It is derived in their §3.1 and recovered from the DTBA in §5.3, in their additive normalization.

Paste-ready, **line 180**, replace the last three sentences, from "With (-@eq-dpt-particle), …", by:

```markdown
With (-@eq-dpt-particle) the same derivation, which Bajnok and Simon note applies to any diagonal scattering theory [@bajnok_simon2008, Sec. 7], gives the $a_r$ defect TBA with fugacities $F_a(\theta)=\mathbb T_{\bar a}(\frac{i\pi}2+\theta)$ in $\log\big(1+F_a\,e^{-\epsilon_a}\big)$. For real $\theta_d$, $F_{\bar a}=\overline{F_a}$ and $|F_a|<1$ on the real line; for $h$ even, $F_{h/2}=-\tanh\frac{\theta-\theta_d}2$ is real but changes sign, so a driving term $-\log F_{h/2}$ needs a branch. The equations are invariant under complex conjugation combined with $a\leftrightarrow\bar a$, so the ground-state energy is real for a solution with this symmetry. Its ultraviolet limit, and a numerical solution, have not been worked out.
```

**Jiang (line 183).** Full text inspected. The abstract and the conclusions say that topological, purely transmitting defects do not affect the universal UV scaling of the entanglement entropy. They also say the defect "modifies slightly the leading exponential correction in the IR". The leading correction (5.5) depends only on the spectrum and the transmission amplitude, and the defect contribution is (5.6). The examples are sinh-Gordon (§5) and scaling Lee–Yang (§6). Jiang computes only the leading IR correction, so "changes **only** the leading infrared corrections" is not supported. In addition, "(from the abstract)" should go, since the paraphrase is now checked against the paper.

Paste-ready, **line 183**:

```markdown
- **Entanglement.** Using the form-factor bootstrap, Jiang showed that a topological, purely transmitting defect does not change the universal ultraviolet scaling of the entanglement entropy, and computed its contribution to the leading exponential infrared correction, which depends only on the spectrum and the transmission amplitude. His examples are sinh-Gordon theory and the scaling Lee–Yang model, one unitary and one non-unitary [@jiang2017, Secs. 4–6].
```

**BHW 2014 and Kormos–Runkel–Watts 2009** (line 181–182): accurate. BHW compare defect TCSA with the DTBA, Lüscher corrections and the defect bound states; KRW study RG flows of perturbed topological defects by TCSA.

**Absence claims.** "No one-loop, TBA or truncated-space calculation exists for an $a_r$ defect" (line 134), "This TBA has not been written down" (line 180) and "No defect TBA, truncated-space or Lüscher calculation exists for an affine Toda defect beyond sinh-Gordon theory" (line 221) cannot be established by a bounded audit. Replace "No … exists" by "We know of no …" at lines 134 and 221. Line 180 is covered above.

### 7c. The fifteen open problems

| # | Verdict |
|---|---|
| 1 Anomalies | Accurate (PP §§4–5: spin-3 \(a_r\), spin-4 shG, \(r\le5\), the \(\sqrt{1+\gamma/2}\) momentum fix). See §3 for the one-loop qualification. |
| 2 Second chirality | Accurate. "In the quantum theory the type-I\(_-\) defect keeps…" refers to the axiomatic matrices, not to the Lagrangian; write "The axiomatic type-I\(_-\) transmission matrices keep …". `prp-spin-two` not re-derived here. |
| 3 UV–IR | Accurate, conditional as stated. BS leave the mass-scale renormalization to a further calculation (§3 above). |
| 4 Classical integrability | Accurate. Caudrelier's time-like bracket and Avan–Doikou's reported discrepancies are cited for what they say. |
| 5 Beyond \(a_r\) | Needs bounding. Bristow-Bowcock (1612.03002) finds **momentum-conserving** type-II defects for \(A_n,B_n,C_n,D_n\), and says the Bäcklund structure "strongly suggests" integrability. Their conclusion lists \(E_6,E_7,E_8\) as the remaining cases. Replace "Type-II defects are known for …" by "Momentum-conserving type-II defects, with Bäcklund transformations, are known for …", and "No defect is known for $e_6^{(1)}$, …" by "No Lagrangian jump defect has been found for $e_6^{(1)}$, $e_7^{(1)}$, $e_8^{(1)}$ [@bristow_bowcock2017]". "Quantum transmission matrices exist only for \(a_2^{(2)}\)" is correct in context (beyond \(a_r^{(1)}\)). The closing "suggests that @prp-lightest-breather-universal holds for every integrable jump defect" is suitably tentative. |
| 6 Fusing | Accurate. Robertson 2014b's fused \(a_3\) matrix is acknowledged as unchecked. |
| 7 Type I\(_+\) | Accurate. Hypothesis (ii) is proved only for \(a_2\), consistent with ch. 41 line 44. |
| 8 Sign of \(q\) | Accurate as a convention question. |
| 9 All-\(r\) breathers | Accurate. A finite grid is not a proof. The last sentence agrees with the corrected line 18 (only \(B_1\) is coupling-independent). |
| 10 Spectrum | **Wrong** rank sentence and overstated existence claim; see §6. |
| 11 TQ | Accurate. RW (5.9) delta factor, (5.10) stripped matrix, checks "up to N = 10"; Weston §2.4.1; BLZ. Runkel 2010 says that convergence "will not be addressed here" and treats the perturbed defect operators as formal power series. No novelty is claimed. |
| 12 Doubled picture | First item now closed (§5). Kolb 2014 §2 and Appel–Vlaar 2022 §3.6 assume an indecomposable Cartan matrix (checked in the sources); the diagonal \(A\oplus A\) is decomposable, so the statement is accurate as a statement of those theorems' scope. |
| 13 Boundaries | Accurate as a list. The ch. 42 dynamical-boundary results were not re-refereed here. |
| 14 Reality | Accurate. Vacuum proved (§4); no Krein or Bethe–Yang reality result is claimed. |
| 15 Finite volume | Bound the absence claim (above). |

### 7d. Novelty (line 259)

"New here: … the particle factors (-@eq-dpt-particle) as solutions of the axioms" must credit Corrigan and Zambon. Their (7.11) gives the classical linear factors for **every** \(r\) (with the denominator misprint corrected in ch. 39 and Appendix D). Right after, they observe that for \(a_2\) these equal the lightest-breather factors (7.5), and they take this as support for their matrices (CZ 2007, §7). What is new here is the proposal of these functions as exact real-coupling factors, with the scalar bootstrap and limits, the all-species breather check for \(r\le4\), the vacuum theorem, and the doubled-symmetry test.

Paste-ready, **line 259**, replace the last sentence by:

```markdown
New here: the all-species check of @prp-lightest-breather-universal for $r\le4$ and the ordering of the constituents; the proposal of the classical factors of [@cz2007, eq. (7.11)] (@prp-an-linear-transmission) as exact particle transmission factors, (-@eq-dpt-particle), with their bootstrap and limits; @cnj-transmission-continuation with its assumptions; @prp-dpt-vacuum; the duality remark; and the test of the $a_r$ expectation of @sec-jump-coideal, with the second module structure (@prp-dpt-doubled).
```

At **line 263**, prefix the first bullet with "Clarification (conventions):". The label exchange is a dictionary between transposed conventions, not an error in CZ. Appendix D has no entry for it, which is correct.

## 8. Status labels, exercises, bibliography

**Verdict: the single `[Theorem]` is proved. The `[Numerical]` labels are justified for what was computed. Seven of the eight exercises are correct; `exr-dpt-imaginary` asks the reader to prove a false statement. The citations are accurate apart from the items listed.**

**Status labels.** The only `[Theorem]` is `prp-dpt-vacuum` (line 115), and §4 proves it independently. `prp-lightest-breather-universal` and `prp-dpt-doubled` are correctly `[Numerical]`. `prp-an-particle-transmission` (line 77) is labelled `[Numerical]` for \(h=2,\dots,8\). Its items 1–4 are exact rational identities: [supplemental_checks.py](../referee_scratch/defects_43/supplemental_checks.py) reduces all four bootstrap residuals to 0 with sympy, and unitarity, Hermitian analyticity and the pole positions are one-line computations. Item 5 is `prp-an-linear-transmission`, proved in ch. 39 for every \(\gamma_k\). So the proposition could be `[Theorem]` for every \(h\), with the hand proof of `exr-dpt-telescope` moved into the proof. This is optional; the current label understates, it does not overstate. The \(a_5\) run at line 50 is described as a scratch run "not in the script". Neither referee reproduced it (see the unverified observations).

**Exercises (all eight checked).**

| Exercise | Verdict | Evidence |
|---|---|---|
| `exr-dpt-telescope` | correct | \(\mathbb T^{\rm cl}_a=e^{i\gamma_a}\cosh\frac12(D+i\gamma_a)/\cosh\frac12(D-i\gamma_a)\) by multiplying through by \(e^{D/2}\); the cancellation in the hint is exact; four exact sympy zeros; block forms to \(5.1\times10^{-16}\). |
| `exr-dpt-cz` | correct | §1: (7.5) simplifies to \(\mathbb T^{\rm cl}_{2},\mathbb T^{\rm cl}_1\) to \(4.5\times10^{-16}\); (7.3) ordering checked on the rendered page ([cz75.png](../referee_scratch/defects_43/sources/cz75.png)). |
| `exr-dpt-vacuum` | correct | §4 (AM–GM proof). At imaginary coupling the exponents are imaginary and the inequality has no meaning, which is the intended answer. |
| `exr-dpt-bs` | correct | §3: \(\partial_\mu\) identity exact; \(\theta_d^{\rm cl}(1-\frac{b^2}{8\pi})\) is the \(O(b^2)\) truncation; \(1-\frac B2\) and \(e^{-B/2}\) differ at \(O(B^2)=O(b^4)\). |
| `exr-dpt-limits` | correct | \(\theta_d\to+\infty\): \(e^{-(\theta-\theta_d)}\to\infty\), \(\mathbb T_a\to1\); \(\theta_d\to-\infty\): \(\mathbb T_a\to e^{2i\gamma_a}\); exchanged for \(\mathbb T^+_a\). On the weight-zero singlet the stored-charge twist phases cancel. |
| `exr-dpt-imaginary` | **wrong as posed** | "Show that \(\mathbb T_a\) is no longer a phase for real \(\theta\)" is false at \(v=\pi\), which is inside the window for \(a<\frac h2\) (§6). Residue \(e^{2i\gamma_a}-1\): the numerator at the pole is \(2i\sin\gamma_a\) and the derivative of the denominator is \(e^{-i\gamma_a}\). |
| `exr-dpt-cdd` | correct; wording | The bulk-\(S\) factor satisfies all three conditions (Codex: bootstrap \(5.4\times10^{-15}\), 224 fusings). "Standing at rapidity \(\theta_0\)" only for \(\theta_0=0\) (BG). The question about \(\omega\)-dependence is answered in §7a. |
| `exr-dpt-doubled` | correct | Derivation below. |

`exr-dpt-doubled`, by hand with `eq-dbd-coproduct`. Write \(p=\psi(j)\), \(X=e_j\otimes1+k_j\otimes e_p\), \(K=k_j\otimes k_p\). Then

\[
\Delta_AX=(e_j\otimes1)\otimes1+(k_j\otimes1)\otimes X+(k_j\otimes e_p)\otimes K .
\]

On \(V\otimes\mathcal W\) this is \(\pi(e_j)\otimes1+\pi(k_j)\otimes Y\), since \(\varepsilon(e_p)=0\). On the reflected copy it is \(1\otimes Y+\pi'(e_p)\otimes\kappa\), with \(\kappa=\chi(K)\). For \(F=f_j\otimes k_p^{-1}+1\otimes f_p\),

\[
\Delta_AF=(f_j\otimes k_p^{-1})\otimes K^{-1}+(1\otimes k_p^{-1})\otimes F+(1\otimes f_p)\otimes1 ,
\]

which gives \(\pi(f_j)\otimes\kappa+1\otimes Y\) and \(\pi'(k_p)^{-1}\otimes Y+\pi'(f_p)\otimes1\), now with \(\kappa=\chi(K^{-1})\). These are the two equations at line 191, and charge conservation gives exactly the weights stated at line 202. The \(h=2\) case is §5.

Paste-ready, **line 246** (`exr-dpt-imaginary`):

```markdown
Let $\theta_d=iv$. Show that $\mathbb T_a$ has a pole at $\theta=i(v-\pi+\gamma_a)$, with residue $e^{2i\gamma_a}-1$, which lies in the physical strip for $\pi-\gamma_a<v<\frac{3\pi}2-\gamma_a$, and compute the energy $m_a\cos(v-\pi+\gamma_a)$ of the would-be bound state. Show that $\mathbb T_a$ is a phase for real $\theta$ only if $v\in\pi\mathbb Z$. For $a<\frac h2$ the window contains $v=\pi$, where $\sigma=-1$: what does @prp-dpt-vacuum say about that defect, and what does this mean for interpreting the pole? For sinh-Gordon theory, compare with the Lee–Yang defect of @sec-defect-lee-yang, whose defect rapidity is also imaginary.
```

At **line 250** (`exr-dpt-cdd`), replace "Interpret it as a particle of species $e$ standing at rapidity $\theta_0$" by "Interpret it as a defect that scatters like a particle of species $e$ at rapidity $\theta_0$, standing for $\theta_0=0$", and add after the second sentence "Which condition must $\varphi$ satisfy for $g\varphi$ to solve the crossing relations of @sec-transmission-unitarity-crossing?".

**Bibliography.** All 29 cited keys resolve in the bibliography, and all 29 primary texts were downloaded from arXiv and searched for the cited passages ([manifest](../referee_scratch/defects_43/citation_manifest.json), with the chapter line of every citation). The second referee re-checked a subset directly against the source text: BS §§3.1, 5.1, 5.3, 6.1, 7; Jiang §§4–7; BG after (18); CZ 2007 §7, (7.10)–(7.11); KL (5.30)–(5.37); Gandenberger 1998 §3.6; Kolb §2; AV §3.6; Runkel 2010; RW §5; Bristow-Bowcock 2017.

| Key | Status |
|---|---|
| `bcz2005`, `konik_leclair1999`, `cz2010`, `cz2011`, `cz2012`, `penati_polvara2019`, `bajnok_simon2008` (one loop, DTBA, LY), `bhw2014`, `kormos_runkel_watts2009`, `akr2002`, `delius_gandenberger1999`, `caudrelier2015`, `avan_doikou2012b`, `robertson2014a`, `robertson2014b`, `cz2009b`, `blz1997`, `rossi_weston2002`, `weston2010`, `runkel2008`, `runkel2010`, `kolb2014`, `appel_vlaar2022`, `gandenberger1998` | accurate as cited |
| `gandenberger1995` (line 20) | covers \(a_2\) only; add `@gandenberger1998, Sec. 3.6` for general \(a_n\) (§1) |
| `cz2007` | (7.5), (7.3) accurate; credit (7.11) in the notes (§7d); label swap is a convention (§7d) |
| `jiang2017` | paraphrase overstates (§7b) |
| `bajnok_george2006` | "standing" only at \(\theta_0=0\) (§7a, §8) |
| `bristow_bowcock2017` | momentum-conserving, integrability suggested (§7c, problem 5) |
| `bajnok_simon2008`, line 180 | "Sec. 5, eq. (18)": (18) is in §3.1 and is recovered in §5.3; write "eq. (18), Sec. 5.3" |

## Other observations (unverified)

The following were noticed during the review but are **not** established by an independent computation here. They are kept separate from the verified corrections above.

1. **All-rank breather statement.** Neither referee has a proof of `prp-lightest-breather-universal` for general \(r\). The evidence covers ranks 2–4 at generic points: Codex's grid, plus the book's grid. Root-of-unity couplings, where the singlet might degenerate, were not examined.
2. **The \(a_5\) run (line 50).** It is described as a ten-minute scratch run that is not in the script. Neither referee reproduced it. Either add it to the script as an optional check or describe it as unrecorded.
3. **The species-\(r\) excited defect may be the easier case of problem 10.** At the species-\(r\) pole the residue of \(\mathbb T_{\bar V}\) has rank 1 (§6, verified). So its image, which the book takes to be the excited defect, is one-dimensional in each sector, as for sine-Gordon theory. The bootstrap on that image may then be computable in closed form. This was not attempted.
4. **Quantum integrability and the Lagrangian.** No check in this review bears on whether the Lagrangian defect is quantum integrable, or on its identification with the axiomatic amplitudes. The chapter flags this correctly as a conjecture.
5. **The candidate stable states.** The residues, bootstrap, energies and reality of the \(\theta_d=0\) soliton candidates and of the imaginary-\(\theta_d\) particle poles were not analysed beyond the \(v=\pi\) observation of §6.
6. **TBA reality for all \(L\).** The symmetry argument of §7b gives a real \(E_0\) for a conjugation-symmetric solution. That the physical solution has this symmetry for every \(L\), rather than only at large \(L\), was not shown.
7. **Literature absence and priority.** The source audit covers the 29 cited works. There was no systematic search of later literature, which is why the corrections say "We know of no …" rather than "no … exists".
8. **Knock-on outside chapter 43.** The \(v\in\pi\mathbb Z\) exception of §6 applies verbatim to the sentence "Then $\sigma$ is a phase and $\Theta$ is broken" in `sec-sg-defect-spectrum` (ch. 37): \(v=\pi\) lies in its window when \(\lambda<1\). The argument is the one verified in §6, but chapter 37 was not audited.
9. **Not re-refereed.** The classical spin-\(\pm2\) claims (`prp-spin-two`) behind problem 2, the chapter-42 dynamical-boundary results behind problem 13, and the claims of Robertson's fused \(a_3\) matrix were taken as stated.

## Integration (2026-10-13)

The referee's scripts are in `briefs/referee_defects_43_scripts/` with their outputs. They were not used as verification. The integrating session (Claude Opus 5.5) checked every item with its own scripts, written from the chapter's definitions and importing only the book helpers (`defect_common`, `boundary_common`): `integration_breathers.py` (higher sine-Gordon breathers bootstrapped from `sg_tmatrix` at two couplings; the threshold $h\omega>2$), `integration_particle_identities.py` (sympy: all four bootstrap identities, unitarity/crossing and Hermitian analyticity for symbolic $u$, $v$, both types; pole representatives), `integration_bs_oneloop.py` (sympy: Bajnok–Simon §4.2 and (17), (23)), `integration_vacuum.py` (the direction $\phi=\psi=tv$ against the closed form, and the total static energy with linear layers), `integration_second_module.py` (exact propagation on dicts of labels, no window: $\tilde\rho$ for $V$ and $\bar V$, $h=3,4,5$, all nodes; the Cartan equation; the free constant; all $U_q(\mathfrak b^+)$ relations; the $h=2$ reduction to ch. 37 and ch. 42), `integration_spectrum.py` (sector blocks of `eq-dt-tv`/`eq-dt-tvbar` at the species-1 and species-$r$ poles, $h=3,\dots,6$, random offsets; the pole of $\bar g$; Lee–Yang; $\theta_d=iv$), `integration_cdd_tba.py` (both crossing relations with $g\to gf$; the examples; the breather multipliers; the TBA fugacities) and `integration_refcheck.py`; all with `*_output.txt` in the same directory. `integration_book_b3_a5_output.txt` is the record of the new optional $a_5$ check. Source passages read: Gandenberger hep-th/9703158 §3.3.1 ((3.33)–(3.34), $p\le\frac{h\lambda}2$, $\lambda$ = book $\omega$) and §3.6 ((3.48)–(3.52)); Bajnok–Simon 0712.4292 §3.1 ((16)–(18)), §4.1 ((23)), §4.2, §5.1 ((25) and the fugacity form), §5.3, §7; Bajnok–George hep-th/0404199 (18) and the sentence after it; BCZ hep-th/0506169 (5.7)–(5.10); Jiang 1703.03562 abstract, §1, (5.5)–(5.6), §7; Bristow–Bowcock 1612.03002 abstract and §5; CZ 0705.1066 §7, (7.3)–(7.11) and the sentence after (7.11).

- **Item 1. Verified.** $B_2$ bootstrapped from the soliton matrix at $D=0.13$: $-0.951042862-0.309059015i$ at $\lambda=2.6$ and $-0.956504391-0.291717928i$ at $\lambda=3.1$, equal to the product formula to $4\times10^{-14}$; $B_1$ agrees between the couplings to $2\times10^{-13}$. Gandenberger §3.6 does compare the general $a_n$ lowest breathers with the particles; his §3.3.1 gives $p\le\frac{h\omega}2$, so the physical threshold is $h\omega>2$, now cited for it. Changes: line 18 (only $B_1$ is coupling-independent; BCZ (5.7), (5.10)); line 20 (`@gandenberger1998, Sec. 3.6`); two sentences after the mass (physical only for $h\omega>2$, otherwise an algebraic continuation); "generic $\omega$" in line 24 and in item 1 of the proposition (beyond the report: at $h\omega=2$ the reversed shifts coincide with the original ones); the proof notes that $h=3$, $\omega=0.6$ tests only the algebraic identity. No label change.
- **Item 2. Verified; optional upgrade done.** All identities are exact (sympy zeros for symbolic $u,v$, both types). Line 69 now states the ground-state scalar ansatz explicitly; line 102 is retitled "No bound-state poles" and scoped to the minimal ground-state factors. **Label change:** `prp-an-particle-transmission` [Numerical, $h=2,\dots,8$] → **[Theorem]**, every $h\ge2$; the proof block now contains the full hand proof (items 1–4 for both types, both fusings; item 5 is `prp-an-linear-transmission`); check (c) contains the exact identities. The table cell of the row "$a_r$ particles, real $\theta_d$" says "[Theorem] for the factors".
- **Item 3. Verified.** The one-loop identity, (17) $\Rightarrow B_{\rm BS}\vartheta=(1-B_{\rm BS})\mu$ and (23) are exact in sympy. Changes: line 12 conditional on the continuation conjecture; line 134 ("We know of no …"; a failed prediction falsifies the conjecture but does not by itself identify the Penati–Polvara anomaly).
- **Item 4. Verified; no label change.** The proof's last sentence is replaced by the explicit direction $\vec\phi=\vec\psi=t\vec v$ and the linear-layer estimate, which supports the chapter's "energy unbounded below" (lines 109, problem 14). Check (g) also tests the closed form along that direction.
- **Item 5. Verified.** The $y^0$, $y^1$, $y^2$ equations were re-derived by hand; $\kappa$ is the plain shift, $Y_N=\frac t{q-q^{-1}}q^{-n_j}+c_jq^{n_j}$, and the Cartan relation forces $c_j=0$ ($q^4\ne1$). The $\bar V$ equations, the Cartan-element equation, all relations and the $h=2$ reduction (entries of $q^{s/2}M_V$ against `eq-sg-transmission-matrix` to $9\times10^{-16}$, $\tilde\rho$ against `eq-dbd-second-module` to $4\times10^{-15}$, shifts $+1,+2,-1,-2$) were checked on exact labels. Beyond the report: the failing twists fail for any $\theta$-independent $Y,\kappa$ by a $y$-power argument on $V$ alone (e-type $\psi(j)\ne j-1$ at $y^2$ on $v_{j+1}$; f-type $\psi(j)\ne j$ at $y^{-1}$ on $v_{\psi(j)}$), for every $r$, which closes the report's remark that the book's check (f) omits the $a_3$ f-type twists. The coideal obstruction for $a_r$ ($q^3\ne1$) was re-derived by hand. Changes: statement (uniqueness sentence), proof (new display `eq-dpt-second-module`, the derivation, the $h=2$ reduction, the obstruction paragraph), problem 12 (module relations removed; the coideal argument; the ch. 42 Kolb/Appel–Vlaar clause kept), ch. 42 "Beyond sine-Gordon" last sentence. Label stays [Numerical] (the $\bar V$ half is numerical).
- **Item 6. Verified, with one modification.** Ranks: rank $h-1$ at the species-1 pole of $\mathbb T_V$ and rank 1 at the species-$r$ pole of $\mathbb T_{\bar V}$ ($X^h=q$) in every sector, $h=3,\dots,6$; the pole of $\bar g$ is simple and uncancelled. Lee–Yang: poles in $\bar{\mathbb T}_b$ only. $|\mathbb T_a|=1$ on the real line for $\theta_d=iv$ iff $v\in\pi\mathbb Z$; $v=\pi$ is in the particle window iff $a<\frac h2$. *Modification:* the sine-Gordon window contains $v=\pi$ for $\frac12<\lambda<1$, not for all $\lambda<1$ (§6 and observation 8). Changes: table rows 4, 6, 9 and the caption ("bound-state poles … of $\bar{\mathbb T}_b$"), the first sentence of "Reading the table" (real $\theta_d\ne0$; the next sentence also scoped to the factors (-@eq-dpt-particle)), the $v=\pi$ sentence, problem 10 (both ranks; "have poles in the physical strip"), ch. 37 `sec-sg-defect-spectrum` ("$\Theta$ is broken unless $v\in\pi\mathbb Z$", with the window $\frac12<\lambda<1$). Ch. 41 not changed: its rank claim is for $\mathbb T_V$ only.
- **Item 7. Verified, with one modification.** (a) Both crossing relations with $g\to gf$ reduce to the string condition (the first at a shifted $D$), now `eq-dpt-soliton-cdd`; $\mathbb T^{\rm cl}_1$ satisfies it for every $h$ (by hand, $\prod_k(u-x\omega^{-k})/(u^{-1}-x\omega^{-k})=(u^h-x^h)/(u^{-h}-x^h)=1$ with $u^{\pm h}=-1$; numerically $h=2,\dots,12$), the sinh-type factor fails it at $h=3$ (0.903) and $h=5$ (it passes for even $h$, where the string pairs $D$ with $D+i\pi$); the breather multipliers agree with the report. (b) TBA: fugacity form, conjugacy, $|F_a|<1$, the sign change of $F_{h/2}$, conditional reality, BS §7, and the location of (18) all confirmed in the source. Jiang: rewritten from the paper; "(from the abstract)" dropped. Lines 134, 221: "We know of no …". Problem 2: "The axiomatic type-I$_-$ transmission matrices keep …". Problem 5: momentum-conserving type-II defects, Bäcklund transformations suggest integrability, $e_6$–$e_8$ listed by Bristow and Bowcock as remaining ("We know of no defect for …"). Line 97 and `exr-dpt-cdd`: "standing" only for $\theta_0=0$. *Modification* (7d): CZ take the agreement of (7.11) with (7.5) as "a further justification for the choice made in (6.12) for the constant $\Delta$", so the Notes say that, not "support for their matrices"; the novelty list credits CZ (7.11) and adds `eq-dpt-soliton-cdd` and the second module. Line 263 is labelled "Clarification (conventions)"; Appendix D correctly has no entry for it (and already has the (7.11) misprint).
- **Item 8. Verified.** `exr-dpt-imaginary` rewritten (residue $e^{2i\gamma_a}-1$ checked by hand; $v\in\pi\mathbb Z$; the $v=\pi$ question); `exr-dpt-cdd` asks for `eq-dpt-soliton-cdd`; the BS citation of line 180 fixed. The $a_5$ run is now the optional check (b3) (`DPT_A5=1`); run once: PASS, all five species, $4\times10^{-9}$, about 2.3 minutes, not ten. The chapter says so.

**Rejected:** nothing. **Modified:** item 6 (the sine-Gordon window statement) and item 7d (CZ's own claim concerns the constant $\Delta$).

**Other observations of the report.** (1) All-rank proof: still open (problem 9). (2) $a_5$: now check (b3). (3) The rank-1 residue at the species-$r$ pole is stated in problem 10 (one-dimensional image, as for sine-Gordon); its bootstrap was not attempted. (4), (5), (9): unchanged, no claims added. (6) TBA reality: the text says "for a solution with this symmetry". (7) Absence claims bounded. (8) Done (ch. 37), with the corrected window.

**Book script after integration:** `defect_particle_transmission.py` 11 PASS (was 8), 35 s alone, 49 s under load: (c) gains the exact sympy identities; (g) the explicit $\sigma<0$ direction; new (h) the second module and the $y$-power obstructions, (i) the soliton CDD condition, (j) residue ranks, the Lee–Yang direction and $v=\pi$; optional (b3) $a_5$. `defect_common.py` unchanged, so no other script was rerun. Reference check over all 65 book `.qmd` files: 1164 ids (two new: `eq-dpt-second-module`, `eq-dpt-soliton-cdd`), no duplicates, every `@sec/eq/prp/fig/exr/thm/cnj/tbl` reference resolves, every citation key is in `references.bib`.

**Labels changed.** `prp-an-particle-transmission`: [Numerical] ($h=2,\dots,8$) → [Theorem] (every $h\ge2$); `tbl-dpt-spectrum`, row "$a_r$ particles, real $\theta_d$": [Numerical] → [Theorem] for the factors. The upgrade has had no independent referee pass beyond this report, which suggested it.

**Pending.** A written proof of the $\bar V$ half of `eq-dpt-second-module` (sketch, not in the book: the $e$-type and Cartan equations for $V(\theta_1)$, $V(\theta_2)$ imply the $e$-type equation for $(\mathbb T\otimes1)(1\otimes\mathbb T)$ with the coproduct actions, checked by hand; the remaining step is the bootstrap embedding of $\bar V$). No Part VII or `index.qmd` edits are needed.
