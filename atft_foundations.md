# Affine Toda field theories as quantum field theories: a review, a critique, and a foundation

*Draft, October 2026*

## Abstract

Affine Toda field theories (ATFTs) are a family of two-dimensional integrable quantum field theories, one for each affine Kac–Moody algebra $\hat g$. With real coupling they are unitary theories of $\operatorname{rank}(g)$ scalar particles. Their exact S-matrices are known for every $\hat g$ and have been checked perturbatively and numerically.

With imaginary coupling they have a complex action, solitons that are complex classical saddles, and S-matrices that are generally not unitary. These features have long been taken to cast doubt on whether the imaginary-coupling theories exist as quantum field theories at all.

We review the field and list twelve gaps between what is used in practice and what has been established. We then build a foundation that closes as many as possible. Our results are:

- **Positivity.** On a lattice, every imaginary-coupling ATFT is a well-defined statistical model. After integrating out the free field, its partition function is a convergent sum of *positive* Coulomb-gas weights. The "complex action" therefore never produces a sign problem.
- **Krein structure.** The lattice transfer matrix is self-adjoint with respect to an indefinite (Krein) inner product whose fundamental symmetry is field parity $\vec\phi\mapsto-\vec\phi$. It commutes with an antilinear PT-type symmetry $\Theta$. So the spectrum is closed under complex conjugation.
- **Real vacuum energy.** Positivity of $\operatorname{Tr}T^N$ together with Pringsheim's theorem proves that the finite-volume vacuum energy is real and has the smallest real part in the spectrum. Excited levels can be complex, and an exactly solvable reduction shows they are.
- **Continuum limit.** In finite volume, with Wick ordering done site by site, the lattice theory has a continuum limit for $\beta^2<2\pi$ and every algebra, with every term of the Coulomb-gas expansion positive (Theorem E). Without further counterterms the limit fails at and above the first collapse threshold $\beta^2_{\rm UV}$, and with site-independent Wick ordering it already fails at $2\pi$. In the window $2\pi\le\beta^2<\beta^2_{\rm UV}$ every term of the expansion converges (at sketch level), and summability is reduced to a version with cluster factors of a bound on products of nearest-neighbour distances, which we prove (Theorem 6.9). Summability is equivalent to a stability bound on the free energy at real fugacity (Lemma 6.15), whose expected exponent is set by the bulk free energy.
- **Classical reality.** Classical energies are real precisely for $\Theta$-invariant asymptotic data, which turns the old ad hoc restriction to "good" solutions into a symmetry statement.
- **Semiclassics.** The Picard–Lefschetz form of the semiclassical expansion requires eigenvalues to be counted with algebraic multiplicity. This resolves the long-standing disagreements over one-loop soliton masses, which we check against a regularization-independent lattice computation and exact S-matrices. We also derive the counting rule from the transmission factors. The size of each Jordan block equals the order of a zero of a transmission factor, also for eigenvalues embedded in the continuum (Theorem F). The weight with which the one-loop trace counts a given energy is the net order of the product of all transmission factors; this factorization is proved whenever the fluctuation solutions have Hirota form (Theorem G). A zero in one channel can be cancelled by a pole in another (Proposition 8.6). At the $c_n^{(1)}$ blocks, the second unit of weight is carried by a bounded threshold resonance rather than an $L^2$ vector (Proposition 8.7).
- **Scattering.** Braiding (R-matrix) unitarity, crossing and the bootstrap follow from the universal R-matrix [Delius 1995], but they are not QFT unitarity, which in addition requires Hermitian analyticity [Miramontes 1999]. What the Krein structure does imply is that eigenvalues of S-matrices and multi-particle transfer matrices are pure phases or come in pairs $(s,1/\bar s)$. This turns the Takács–Watts criterion for a real finite-volume spectrum into a criterion for unbroken $\Theta$.
- **What cannot be rescued.** We list the physical consequences of unitarity that fail in the unrestricted imaginary-coupling theories of rank $r\ge2$ and that no choice of inner product, basis, gradation or quantum-group restriction can restore, because their failure is detected by basis-independent quantities (Section 10).

The resulting picture is that imaginary-coupling ATFTs are valid **non-unitary** quantum field theories, in the same sense as the scaling Lee–Yang model. They have a positive Euclidean measure, a Krein-space Hamiltonian formulation, a real vacuum energy, real masses for the asymptotic particles in every known exact solution, and a consistent S-matrix theory. Unitarity is not one of their properties: it fails in identifiable sectors, and requiring it was a category error. Section 10 states exactly which of its usual consequences are lost and cannot be recovered. We also correct a premise common in the literature: exact soliton S-matrices are *not* known for every algebra. They are missing, to our knowledge, for $c_n^{(1)}$ with $n\ge3$, $f_4^{(1)}$, $e_6^{(2)}$ and $a_{2n-1}^{(2)}$.

**Status labels.** Every claim is labelled:

- **[Theorem]**, with a complete proof given here;
- **[Proposition, sketch]**, with a proof outline that a specialist could complete;
- **[Established]**, proved or computed in the cited literature;
- **[Conjecture]**, with the supporting evidence stated.

---

## 1. Introduction

### 1.1 Why affine Toda theory matters

For each affine Kac–Moody algebra $\hat g$ there is a Lagrangian field theory of $r=\operatorname{rank}g$ scalars with exponential interactions, the affine Toda field theory. These theories have been a laboratory for several central ideas in two-dimensional quantum field theory:

- exact factorized S-matrices and the bootstrap [Zamolodchikov–Zamolodchikov 1979; Braden–Corrigan–Dorey–Sasaki 1990];
- the geometry of root systems in scattering (Dorey's fusing rule [Dorey 1991]);
- strong–weak coupling duality [Delius–Grisaru–Zanon 1992];
- solitons with quantum-group symmetry [Bernard–LeClair 1991; Hollowood 1993];
- integrable perturbations of conformal field theories [Zamolodchikov 1989; Eguchi–Yang 1989; Hollowood–Mansfield 1989];
- the ODE/IM correspondence [Lukyanov–Zamolodchikov 2010; Dorey–Faldella–Negro–Tateo 2013].

The theories come in two regimes.

- **Real coupling.** The theory is an ordinary unitary QFT with a unique vacuum and $r$ massive particles.
- **Imaginary coupling.** The potential is periodic and complex. There are infinitely many degenerate vacua, solitons interpolating between them, and breathers. These are the theories that describe perturbed (W-algebra) minimal models after quantum-group restriction, and they carry the richest physics.

The theories are not only of historical interest. They are rich enough to show generic phenomena of quantum field theory, yet solvable enough that claims about them can be checked exactly. Five features keep them current.

- **A family, not a single model.** There is one theory for each affine algebra, and results come with Lie-theoretic structure: Dorey's fusing rule, Perron–Frobenius masses, Lie duality. A general claim about integrable or massive two-dimensional QFT can therefore be tested across a whole family, not in one example.
- **Exact benchmarks.** At real coupling the S-matrices are known for every algebra, together with the thermodynamic Bethe ansatz, form factors and the exact mass–coupling relation (Section 3.1). New methods are calibrated against them. Recent examples are tree-level integrability for every ATFT, derived from properties of root systems [Dorey–Polvara 2022], one-loop integrability of the simply-laced theories [Polvara 2023], and the comparison of one-loop amplitudes with the bootstrap S-matrices for the whole simply-laced class [Fabri–Polvara 2024].
- **A non-unitary QFT with structure.** The imaginary-coupling theories are the field-theoretic parents of non-unitary perturbed minimal and W-algebra models. We argue that their non-unitarity is organised by a Krein structure and an antilinear $PT$-type symmetry $\Theta$, and is not a pathology. That is the question non-Hermitian physics now asks of many systems (Section 1.5).
- **Unfinished business.** Soliton S-matrices are still missing for four families (Section 3.3). The disagreements over one-loop soliton masses needed a new counting rule to resolve (Section 8). The missing-topological-charge problem is open. Open problems in a theory whose surroundings are this well mapped are unusually well posed.
- **A meeting point of mathematics and physics.** Quantum affine algebras and their R-matrices, Hirota tau-functions, Coulomb gases and Picard–Lefschetz theory all act on the same model. Progress in any one of them tends to have consequences here.

### 1.2 Three objections

Three objections are routinely raised against the imaginary-coupling theories.

1. **Complex action.** For real fields the Euclidean weight $e^{-S_E}$ is complex, so it is unclear what measure is being integrated and whether the theory has a sign problem.
2. **Complex saddles.** The solitons are complex solutions of the field equations. A path integral over real fields gives no obvious reason to include them, and no criterion for which ones to include [MacKay–Watts 1995].
3. **Non-unitarity.** The Hamiltonian is not Hermitian and the S-matrices are not unitary. In some sectors the finite-volume spectrum is complex [Takács–Watts 1999, 2002].

### 1.3 What we mean by a valid QFT

We adopt operational criteria that apply equally to unitary and non-unitary theories, the scaling Lee–Yang model being the standard example of the latter [Fisher 1978; Cardy–Mussardo 1989].

- **(C1) Definition.** A regularized theory exists, its correlation functions are finite, and there is a continuum limit.
- **(C2) Euclidean consistency.** The Euclidean correlation functions of a generating set of local observables satisfy the Osterwalder–Schrader axioms. Reflection positivity may be replaced by a twisted reflection positivity that yields a Krein space.
- **(C3) Hamiltonian structure.** There is a Hamiltonian, self-adjoint with respect to a possibly indefinite inner product, with an antilinear symmetry. Its vacuum energy is real.
- **(C4) Asymptotic spectrum.** The one-particle states (particles, solitons, breathers, excited solitons) have real masses.
- **(C5) Scattering.** The S-matrix satisfies factorization, the Yang–Baxter equation, crossing, the bootstrap, and pseudo-unitarity with respect to the Krein structure.
- **(C6) Predictivity.** Finite-volume spectra, form factors and the mass–coupling relation are computable and agree with independent numerical or nonperturbative checks.

Unitarity (a positive inner product preserved by the S-matrix) is **not** on this list. Section 9 shows that it fails in identifiable sectors of most imaginary-coupling ATFTs, and Section 10 lists which of its usual consequences cannot be recovered by any means. That failure is a physical property of these theories, as it is of Lee–Yang, and not an inconsistency.

### 1.4 Summary of results

| Gap | Content | Status after this paper |
|---|---|---|
| H1 | existence of the regularized theory | closed: Theorems A and B (lattice) |
| H2 | sign problem from the complex action | closed: Theorem B (positive Coulomb gas) |
| H3 | continuum limit | closed in finite volume for $\beta^2<2\pi$ (Theorem E); for $2\pi\le\beta^2<\beta^2_{\rm UV}$ every term converges (Proposition 6.13, sketch) and summability is reduced to Conjecture 6.14, a version with cluster factors of a geometric bound proved as Theorem 6.9; further counterterms needed at and above $\beta^2_{\rm UV}$ (Proposition 6.7); infinite volume open (Conjecture 6.2) |
| H4 | Hilbert space and Hamiltonian | closed: Theorem C (Krein space, $\Theta$-symmetry) |
| H5 | reality of the vacuum energy | closed: Theorem D |
| H6 | reality of classical energies; singular and complex-energy solutions | closed: Proposition 7.1 and the $\Theta$ selection rule |
| H7 | which complex saddles contribute | reduced: thimble formulation; masses shown independent of the choice |
| H8 | semiclassical mass corrections | closed: algebraic-multiplicity rule, with lattice and exact checks; counting rule derived from the transmission factors (Theorems F and G, Proposition 8.5) and checked for $a_n^{(1)}$ and $c_n^{(1)}$ (Propositions 8.6 and 8.7) |
| H9 | non-unitarity of S-matrices and complex finite-volume spectra | reclassified: Krein-unitarity (Proposition 9.1), exact on an integrable lattice regularization (Proposition 9.2); sector-by-sector classification; consequences of unitarity that cannot be rescued identified (Section 10) |
| H10 | completeness and correctness of exact S-matrices | partly open: soliton S-matrices missing for four families |
| H11 | missing classical topological charges and solutions | open, reformulated via thimbles |
| H12 | unitary subtheories and RSOS restriction | established in specific cases; general criterion conjectured |

### 1.5 Connections to current research

Read as statements about a non-unitary quantum field theory, these results connect to several active fields. The links are strongest for non-Hermitian physics and for resurgence. For both, the theorems above give exact field-theoretic examples of ideas that are mostly studied in quantum mechanics or in toy models. The connections below are our assessment, and the citations sample each field rather than survey it.

- **Non-Hermitian and $PT$-symmetric physics.** The transfer matrix is self-adjoint for a Krein form and commutes with an antilinear symmetry $\Theta$ (Theorem C). So its spectrum is closed under complex conjugation, with $\Theta$-unbroken and $\Theta$-broken sectors (Section 9.2). The Jordan blocks of the soliton fluctuation operators (Section 8) are exceptional points, and Section 10 (U4, U8) spells out what they imply. Pseudo-Hermiticity, $PT$ symmetry, Jordan forms and exceptional points are the organising concepts of non-Hermitian physics [Bender–Boettcher 1998; Mostafazadeh 2002; Ashida–Gong–Ueda 2020]. That field now spans photonics, open quantum systems and topological phases.
- **Resurgence, Lefschetz thimbles and complex saddles.** The complex solitons are critical points of the thimble decomposition of the lattice integral (5.1). Their one-loop determinants must be counted with algebraic multiplicity (Theorems F and G). Intersection numbers decide which topological sectors contribute (Conjecture 8.3). Resurgence relates perturbative data to the global structure of saddles in quantum mechanics and QFT [Witten 2010; Dunne–Ünsal 2015; Dunne 2025].
- **The sign problem.** Theorem B turns a complex action into a positive Coulomb gas. This is a rare case where a complex weight is provably harmless. Complexifying field space by Picard–Lefschetz theory is a leading approach to sign problems such as QCD at finite density and the Hubbard model [Alexandru–Başar–Bedaque–Warrington 2022]. Here the answer is known exactly, which makes these theories a test case for those methods.
- **Constructive QFT and probability.** The continuum limit of Section 6 is a statement about multi-species Coulomb gases and imaginary multiplicative chaos [Junnila–Saksman–Webb 2020]. The conformal endpoint is now being built rigorously. Compactified imaginary Liouville CFT has been constructed and shown to satisfy the CFT axioms [Guillarmou–Kupiainen–Rhodes 2025], and the construction has been extended to higher-rank compactified imaginary Toda CFT [Yao 2026]. Imaginary-coupling ATFT is this Toda theory perturbed by the vertex operator of the affine root.
- **S-matrix bootstrap and perturbative integrability.** Crossing and the bootstrap have a Hopf-algebraic origin [Delius 1995]. Krein-unitarity replaces unitarity (Propositions 9.1 and 9.2), and Hermitian analyticity is the missing condition for QFT unitarity [Miramontes 1999]. The numerical S-matrix bootstrap explores the space of unitary two-dimensional S-matrices, with integrable models on its boundary [Paulos et al. 2017; Bercini et al. 2019]. The ladder of conditions in Section 9.2 raises the analogous question for Krein-unitary S-matrices.
- **Non-unitary and complex CFTs.** RSOS restrictions give non-unitary minimal models of Lee–Yang type. These models distinguish the effective central charge from the central charge, and their complex finite-volume spectra are seen in the truncated conformal space approach (Sections 9.4 and 10, U5–U7). Complex CFTs control walking renormalization-group flows and weakly first-order transitions such as the $Q>4$ Potts model [Gorbenko–Rychkov–Zan 2018]. Conjecture 9.3 would give a criterion, checkable from the S-matrix, for when a restricted model is unitary.
- **ODE/IM, exact WKB and gauge/string theory.** The linear problems of ATFT underlie the massive ODE/IM correspondence [Lukyanov–Zamolodchikov 2010; Dorey–Faldella–Negro–Tateo 2013]. The ultraviolet limit of the $A_r$ massive ODE/IM correspondence gives non-unitary $WA_r$ minimal models [Ito–Shu 2018]. A recent monograph reviews these correspondences, with applications to exact WKB, minimal surfaces and gluon amplitudes in AdS/CFT [Ito–Shu 2025].
- **$q$-deformed world-sheet theories.** IRF bases restore unitarity of $q$-deformed S-matrices [Hoare–Hollowood–Miramontes 2013]. This is the mechanism we propose for unitary restrictions of ATFT (Section 9.4).
- **Integrable circuits and analogue simulation.** The light-cone regularization of Proposition 9.2 is a brickwork circuit of R-matrix gates. It is unitary for a Krein form rather than for the Hilbert-space inner product. The zero-mode model of Section 5.6 and few-site versions of this circuit are small enough that complex-conjugate level pairs and exceptional points could be realised in $PT$-symmetric photonic or cold-atom systems [Ashida–Gong–Ueda 2020]. This link is more speculative than the others.

These connections suggest where the results can be used. First, Theorems A–D use neither integrability nor the Lie-algebraic data, beyond the non-negativity of the fugacities (Section 5.1). So they give a template for other theories with periodic complex potentials, such as complex Liouville and Toda CFTs and their deformations. Second, the open problems of Section 12 are natural targets for these communities. The thimble intersection numbers (problem 4) can be attacked with the algorithms reviewed in [Alexandru et al. 2022]. The missing soliton S-matrices (problem 2) could be constrained by extending the on-shell loop methods of [Polvara 2023; Fabri–Polvara 2024] to non-simply-laced algebras and imaginary coupling. The $\Theta$-broken sectors (problems 6 and 8) can be mapped numerically with the truncated conformal space approach, tensor networks and the lattice of Proposition 9.2. The infinite-volume limit and mass gap (Conjecture 6.2) might be reached by adding the affine perturbation to the probabilistic constructions above, using Theorems B and E as input.

---

## 2. The theories

### 2.1 Lie-algebraic data

Let $\hat g$ be an affine Kac–Moody algebra with simple roots $\vec\alpha_0,\dots,\vec\alpha_r\in\mathbb R^r$, normalized so that the longest roots have $\lvert\vec\alpha\rvert^2=2$. The positive integers $n_j$ (Kac labels, with the normalization of [MacKay–Watts 1995]) satisfy $\sum_jn_j\vec\alpha_j=0$.

The algebras fall into three classes:

- the self-dual algebras $a_n^{(1)},d_n^{(1)},e_n^{(1)},a_{2n}^{(2)}$;
- the dual pairs of non-simply-laced algebras $(b_n^{(1)},a_{2n-1}^{(2)})$, $(c_n^{(1)},d_{n+1}^{(2)})$, $(g_2^{(1)},d_4^{(3)})$ and $(f_4^{(1)},e_6^{(2)})$.

Duality reverses the arrows of the affine Dynkin diagram. Folding a simply-laced diagram by a diagram automorphism gives every non-simply-laced algebra [Olive–Turok 1983].

### 2.2 Lagrangians

In Minkowski signature,

$$
\mathcal L=\tfrac12\partial_\mu\vec\phi\cdot\partial^\mu\vec\phi-V(\vec\phi),
$$

with

$$
V_{\rm real}=\frac{m^2}{\beta^2}\sum_jn_j\left(e^{\beta\vec\alpha_j\cdot\vec\phi}-1\right),\qquad
V_{\rm imag}=-\frac{m^2}{\beta^2}\sum_jn_j\left(e^{i\beta\vec\alpha_j\cdot\vec\phi}-1\right). \tag{2.1}
$$

Here $\beta>0$ in both cases. Expanding about $\vec\phi=0$ gives the classical mass matrix $M^2=m^2\sum_jn_j\vec\alpha_j\vec\alpha_j^{\mathsf T}$. Its eigenvalues are the classical particle masses, which form the Perron–Frobenius eigenvector of the Cartan matrix [Braden et al. 1990; Freeman 1991].

The sign in $V_{\rm imag}$ matters. It makes $\vec\phi=0$ a minimum ($V\simeq\tfrac12\vec\phi^{\mathsf T}M^2\vec\phi$). In the Euclidean theory it gives positive fugacities (Section 5).

The theories are classically integrable. They have a Lax pair and conserved charges of spins equal to the exponents of $\hat g$ modulo the Coxeter number [Olive–Turok 1985]. The higher-spin charges survive quantization [Feigin–Frenkel 1993; Delius–Grisaru–Zanon 1992; Niedermaier 1994].

### 2.3 Imaginary coupling: vacua, solitons, symmetry

$V_{\rm imag}$ is periodic under $\vec\phi\to\vec\phi+\tfrac{2\pi}{\beta}\vec\lambda$ for $\vec\lambda$ in the co-weight lattice, and vanishes at every lattice point. Solitons interpolate between these vacua.

Hirota's method gives explicit solutions,

$$
\vec\phi=\frac i\beta\sum_j\frac{2\vec\alpha_j}{\alpha_j^2}\ln\tau_j,
$$

with each $\tau_j$ a finite sum of products of exponentials of linear functions of $x$ and $t$ [Hollowood 1992; Olive–Turok–Underwood 1993; Kneipp–Olive 1996]. The soliton masses are proportional to particle masses:

- in the simply-laced and twisted cases, to the particle masses of the same theory;
- in the untwisted non-simply-laced cases, to those of the theory based on $(g^\vee)^{(1)}$ ("Lie duality").

The quantum theory has a non-local $U_q(\hat g^\vee)$ symmetry with $q=\pm e^{4\pi^2i/\beta^2}$ [Bernard–LeClair 1991; Felder–LeClair 1992]. Soliton S-matrices are therefore proportional to $U_q(\hat g^\vee)$ R-matrices.

---

## 3. Literature overview

### 3.1 Real coupling: the particle S-matrices are known for every algebra

**Simply-laced algebras.** The S-matrices are products of building blocks

$$
\{x\}=\frac{(x-1)(x+1)}{(x-1+B)(x+1-B)},\qquad (x)=\frac{\sinh\left(\frac\theta2+\frac{i\pi x}{2h}\right)}{\sinh\left(\frac\theta2-\frac{i\pi x}{2h}\right)},\qquad B(\beta)=\frac{1}{2\pi}\frac{\beta^2}{1+\beta^2/4\pi},
$$

constructed by Braden, Corrigan, Dorey and Sasaki, Christe and Mussardo, and Klassen and Melzer (all 1990), with a universal root-system formula due to Dorey [Dorey 1991].

**Non-simply-laced algebras.** The S-matrices were found by Delius, Grisaru and Zanon [1992] and by Corrigan, Dorey and Sasaki [1993]. Here the Coxeter number "floats" with the coupling and interpolates between $h$ and $h^\vee$ of the dual algebra. The S-matrices of a dual pair are exchanged under $\beta\to4\pi/\beta$.

**Checks.** These S-matrices have been checked at one loop [Braden et al. 1990; Delius–Grisaru–Zanon 1992], by Monte Carlo simulation (hep-th/9206112, hep-th/9508007), and more recently by tree-level and loop integrability arguments [Dorey–Polvara 2022].

**Further results.** Further [Established] results include:

- thermodynamic Bethe ansatz [Al. Zamolodchikov 1990; Klassen–Melzer 1990];
- form factors [Fring–Mussardo–Simonetti 1993];
- the exact mass–coupling relation [Fateev 1994];
- ODE/IM descriptions [Lukyanov–Zamolodchikov 2010; Dorey et al. 2013; Ito–Locke 2014].

**Rigour.** The real-coupling potentials are bounded below and real, so these theories fall into the class of exponential interactions constructed rigorously by Albeverio and Høegh-Krohn [1974]. That construction gives the Wightman axioms and a mass gap for a single field in the super-renormalizable range. Extending it to several fields is expected to be routine, but to our knowledge has not been written down.

### 3.2 Imaginary coupling: classical theory

**Solutions.** Single and multi-soliton solutions are known for all algebras [Hollowood 1992; Olive–Turok–Underwood 1993; Kneipp–Olive 1996; MacKay–McGhee 1993; Aratyn–Constantinidis–Ferreira 1993; Zhu–Caldi 1995]. So are breathers [Harder–Iskandar–McGhee 1995] and soliton time delays [Fring–Johnson–Kneipp–Olive 1994].

**Reality.** Although the fields are complex, the energy–momentum tensor is a total derivative up to an improvement term, and the energies are real [Olive–Turok–Underwood 1993; Zhu–Caldi 1995]. So are all the higher conserved charges [Freeman 1995].

**Known problems.**

- The topological charges of single solitons do not fill the expected representations [McGhee 1994].
- Some solutions are singular at isolated spacetime points.
- There are solutions with complex energy that become singular in finite time, and solutions with negative energy [Khastgir–Sasaki 1996].
- In $a_2^{(1)}$, soliton–breather scattering produces *complex* classical time delays and can drive regular initial data into singular final states [Takács–Watts 1999].

### 3.3 Imaginary coupling: quantum theory

**Soliton S-matrices.** Their status is summarized below. They are built from $U_q(\hat g^\vee)$ R-matrices and fixed by crossing, unitarity in the R-matrix sense, and the bootstrap. In every case examined, the lowest breathers reproduce the real-coupling particle S-matrices continued to imaginary coupling (breather–particle identification [Smirnov 1991; Gandenberger 1995; Gandenberger–MacKay 1995; GMW 1996; Takács 1997]).

| $\hat g$ | soliton S-matrix | charge algebra |
|---|---|---|
| $a_n^{(1)}$ | Hollowood 1993; bound states Gandenberger 1995 | $U_q(a_n^{(1)})$ |
| $a_2^{(2)}$ | Smirnov 1991; Efthimiou 1993 | $U_q(a_2^{(2)})$ |
| $a_{2n}^{(2)}$ | Gandenberger 1997 | $U_q(a_{2n}^{(2)})$ |
| $d_{n+1}^{(2)}$ | Gandenberger–MacKay 1995; Delius 1995 ($U_q(c_n^{(1)})$-symmetric S-matrices with non-rigid poles) | $U_q(c_n^{(1)})$ |
| $b_n^{(1)}$ | Gandenberger–MacKay–Watts 1996 | $U_q(a_{2n-1}^{(2)})$ |
| $g_2^{(1)}$ | Takács 1997 | $U_q(d_4^{(3)})$ |
| $d_4^{(3)}$ | Takács 1997 | $U_q(g_2^{(1)})$ |
| $d_n^{(1)},e_n^{(1)}$ | Johnson 1997 (some scalar factors undetermined) | self-dual |
| $c_2^{(1)}\simeq b_2^{(1)}$ | via GMW 1996 | $U_q(a_3^{(2)})$ |
| $c_n^{(1)}$, $n\ge3$ | **not constructed** | $U_q(d_{n+1}^{(2)})$ |
| $f_4^{(1)}$ | **not constructed** (only the $\mathbf{27}$ R-matrix is known [GMW 1996]) | $U_q(e_6^{(2)})$ |
| $e_6^{(2)}$ | **not constructed** | $U_q(f_4^{(1)})$ |
| $a_{2n-1}^{(2)}$ | **not constructed** (spinor R-matrices missing) | $U_q(b_n^{(1)})$ |

The bold entries are, to our knowledge, still open. The obstacle in each case is constructing R-matrices in representations whose tensor products are not multiplicity-free [GMW 1996; Delius–Gould–Zhang 1996].

For $a_{2n-1}^{(2)}$ there is lattice information. Vernier, Jacobsen and Saleur [2016] identify the scaling limit of the staggered $U_q(a_{2n-1}^{(2)})$ vertex model, in two of its three regimes, with imaginary-coupling $a_{2n-1}^{(2)}$ Toda theory, and find $n$ soliton types with coupling-dependent mass ratios [Established; computed, not proved]. They note that no systematic S-matrix study exists for this family beyond $a_3^{(2)}=d_3^{(2)}$, where the lattice reproduces [Gandenberger–MacKay 1995]. In their construction the lattice carries $U_q(a_{2n-1}^{(2)})$, while the soliton S-matrix carries $U_q(b_n^{(1)})$ (table above). In the third regime the same vertex model has a non-compact continuum limit.

**Semiclassical masses.** One-loop soliton masses were computed by the DHN method for $a_n^{(1)}$ [Hollowood 1993], for $c_2^{(1)}$ [Watts 1994], and for all algebras [MacKay–Watts 1995], with a competing result for $c_n^{(1)}$ [Delius–Grisaru 1995]. The exact $b_n^{(1)}$ masses disagree with MacKay and Watts' tabulated semiclassics [GMW 1996], and the issue was left unresolved [Gandenberger–MacKay 1996]. Section 8 resolves it.

**Restrictions.** At $q$ a root of unity, quantum-group (RSOS) restriction yields perturbed minimal models [Reshetikhin–Smirnov 1990; Smirnov 1991; Efthimiou 1993; de Vega–Fateev 1991; Takács–Kausch–Watts 1997; Takács 1997]. Some of these are unitary. Others have complex finite-volume spectra, confirmed by the truncated conformal space approach [Takács–Kausch–Watts 1997; Takács–Watts 2002].

**Non-unitarity.** Takács and Watts [1999, 2002] showed the following.

- In $a_2^{(1)}$, soliton–antisoliton and soliton–breather amplitudes do not have phase eigenvalues for real rapidity. The soliton–breather case is mirrored exactly by complex classical time delays, and it survives RSOS restriction because breathers are quantum-group singlets.
- In $a_2^{(2)}$, the fundamental soliton–soliton amplitudes are well behaved, but soliton–excited-soliton amplitudes are not, so the spectrum is real only in the repulsive regime.
- All theories based on $a_1^{(1)}$ have real finite-volume spectra.
- A formal S-matrix description remains valid, and agrees with numerics, even when the spectrum is complex.

**Two notions of unitarity.** Delius [1995] showed that S-matrices built from the universal R-matrix of $U_q(\hat g)$ satisfy crossing, the bootstrap and the relation $S_{ba}(-\theta)S_{ab}(\theta)=1$ as consequences of the Hopf-algebra axioms. He also showed that the Lorentz spins of the charges fix the gradation, which makes the pole positions coupling-dependent through a "quantum" dual Coxeter number. He stressed that the relation usually called unitarity is unrelated to whether the field theory is unitary.

Miramontes [1999] identified the missing condition. Braiding unitarity is equivalent to QFT unitarity when the S-matrix is also Hermitian analytic, a basis-dependent property that quantum-group S-matrices generally violate in the particle basis. Hoare, Hollowood and Miramontes [2013] later found that, for $q$ a root of unity, the IRF (RSOS) form of the $U_q(su(2)^{(1)})$-based S-matrices they studied (a $q$-deformed world-sheet theory and a generalized sine-Gordon theory) is manifestly Hermitian analytic, and hence unitary.

---

## 4. Critique: where the rigorous treatment has holes

We list the gaps between what the literature uses and what it establishes. In each case we state precisely what is missing.

**H1. Existence of the regularized theory.** For imaginary coupling the Euclidean weight $e^{-S_E}$ is complex. No paper we know of shows that a regularized imaginary-coupling ATFT defines finite, well-defined correlation functions, still less a probability-like measure.

**H2. Sign problem.** A complex weight usually means cancellations between configurations. If those cancellations were severe, the theory would be ill-defined in practice even on a lattice, and expectation values could be dominated by numerical noise.

**H3. Continuum limit.** The interaction is a sum of complex vertex operators. One must ask in which range of $\beta$ Wick ordering suffices, and whether the infinite-volume limit exists and has a mass gap. For real coupling, Albeverio–Høegh-Krohn methods apply. For imaginary coupling there is no analogous result.

**H4. Hilbert space.** The Hamiltonian is not Hermitian in the Schrödinger representation. The literature uses "non-unitary" without specifying an inner product, an adjoint, or the structure that replaces Hermiticity.

**H5. Vacuum.** With a non-Hermitian Hamiltonian, nothing a priori guarantees a real ground-state energy, or even that the vacuum is the state of lowest real energy.

**H6. Classical reality.** Reality of soliton energies has been shown for specific solution classes. There are also complex-energy and negative-energy solutions [Khastgir–Sasaki 1996], and singular ones. No general criterion distinguishes the solutions that should count as physical.

**H7. Complex saddles.** DHN semiclassics expands around complex solitons. The justification usually offered (steepest descent may leave the real contour) does not say *which* complex saddles contribute, nor with what weight [MacKay–Watts 1995].

**H8. Semiclassical quantization.** The fluctuation operator around a complex soliton is not Hermitian. Mode sums with orthonormal eigenfunctions assume it is diagonalizable, and independent calculations disagree ($c_n^{(1)}$ [Delius–Grisaru 1995; MacKay–Watts 1995]; $b_n^{(1)}$ semiclassics versus exact [GMW 1996]).

**H9. Non-unitarity.** Most soliton S-matrices are not unitary. Some amplitudes do not have phase eigenvalues, and finite-volume spectra can be complex [Takács–Watts 1999, 2002]. It is unclear what replaces unitarity and what consistency conditions survive.

**H10. Exact S-matrices.** They are conjectures, built on quantum-group symmetry, R-matrix unitarity, crossing and the bootstrap. They are incomplete for four families (Section 3.3). Even where they exist, identifying the full spectrum is subtle: a closed bootstrap can contain particles that are neither solitons nor breathers [Takács–Watts 1999].

**H11. Missing classical solutions.** The topological charges of classical solitons do not fill the quantum multiplets [McGhee 1994]. Classical counterparts of some quantum bound states are not known [Takács–Watts 1999].

**H12. Unitary subtheories.** RSOS restriction sometimes yields unitary theories and sometimes not, depending on the coupling, the gradation, and whether breathers survive. There is no general criterion.

---

## 5. Foundations I: the lattice theory

### 5.1 Definition

Fix a lattice spacing $a$ and integers $N_t$ and $N_x$. Let

$$
\Lambda=\{(t,x):t\in a\mathbb Z_{N_t},\ x\in a\{1,\dots,N_x-1\}\},
$$

periodic in Euclidean time, with Dirichlet conditions $\vec\phi=0$ at $x=0$ and $x=aN_x$. The lattice action is

$$
S_\Lambda[\vec\phi]=\underbrace{\tfrac12\sum_{\langle uv\rangle}\lvert\vec\phi_u-\vec\phi_v\rvert^2}_{S_0}+a^2\sum_{u\in\Lambda}V(\vec\phi_u),
\qquad
V(\vec\phi)=-z\sum_jn_j\left(e^{i\beta\vec\alpha_j\cdot\vec\phi}-1\right),\quad z>0. \tag{5.1}
$$

Let $C=C_a$ denote the covariance of the Gaussian measure $d\mu_0\propto e^{-S_0}\prod d\vec\phi_u$, which is positive definite because of the Dirichlet conditions.

The bare fugacity may depend on $a$, on the root and on the site. Wick ordering replaces $z$ by $z_{a,j}(u)=z\,e^{\beta^2\lvert\vec\alpha_j\rvert^2C_a(u,u)/2}$, root by root and site by site. Because of the Dirichlet conditions $C_a(u,u)$ is not constant, and a site-independent Wick constant has no continuum limit for $\beta^2\ge2\pi$ (Proposition 6.3). Theorems A–D below use only that each fugacity is non-negative and independent of Euclidean time. They therefore hold verbatim for fugacities $z_{x,j}\ge0$ that depend on the root and on the spatial coordinate $x$, with $z$ replaced by $z_{x,j}$ in (5.1)–(5.3); in particular they hold for Wick-ordered fugacities.

### 5.2 Existence

**Theorem A [Theorem].** For all $a$, $\Lambda$, $\beta>0$ and $z\ge0$, the weight $e^{-a^2\sum V}$ is bounded by 1 in modulus. Hence

$$
Z_\Lambda=\int e^{-a^2\sum_uV(\vec\phi_u)}\,d\mu_0
$$

and $\int F\,e^{-a^2\sum V}d\mu_0$ exist for every $F\in L^1(d\mu_0)$.

*Proof.* $\operatorname{Re}V(\vec\phi)=-z\sum_jn_j(\cos\beta\vec\alpha_j\cdot\vec\phi-1)\ge0$. So $\lvert e^{-a^2\sum V}\rvert=e^{-a^2\sum\operatorname{Re}V}\le1$. $\square$

Theorem A closes the existence half of H1. The complex weight is bounded, so nothing diverges.

### 5.3 Positivity: there is no sign problem

**Theorem B [Theorem].** Write $q_j=\beta\vec\alpha_j$. Then

$$
Z_\Lambda=e^{-a^2\lvert\Lambda\rvert z\sum_jn_j}\sum_{N\ge0}\frac{(a^2z)^N}{N!}\sum_{j_1,\dots,j_N}\Big(\prod_kn_{j_k}\Big)\sum_{u_1,\dots,u_N\in\Lambda}\exp\Big(-\tfrac12\sum_{k,l=1}^Nq_{j_k}\cdot q_{j_l}\,C(u_k,u_l)\Big). \tag{5.2}
$$

Every term is positive and at most 1, and the series converges absolutely. In particular $Z_\Lambda>0$.

More generally, for any real vectors $\vec w_1,\dots,\vec w_M$ and sites $v_1,\dots,v_M$, the quantity

$$
Z_\Lambda\Big\langle\prod_me^{i\beta\vec w_m\cdot\vec\phi_{v_m}}\Big\rangle
$$

is given by the same expansion with the charges $\beta\vec w_m$ added at $v_m$. It is therefore real and non-negative.

*Proof.* Expand $e^{-a^2\sum V}=e^{-a^2\lvert\Lambda\rvert z\sum n_j}\exp\big(a^2z\sum_{u,j}n_je^{iq_j\cdot\phi_u}\big)$ in powers. Use

$$
\int\prod_ke^{iq_k\cdot\phi_{u_k}}d\mu_0=\exp\Big(-\tfrac12\sum_{k,l}q_k\cdot q_l\,C(u_k,u_l)\Big).
$$

The exponent is minus a positive-semidefinite quadratic form (the covariance $C\otimes\mathbb 1$ evaluated on the "charge distribution" $\sum_kq_k\delta_{u_k}$), so each term lies in $(0,1]$. Absolute convergence follows by comparison with $\exp(a^2z\lvert\Lambda\rvert\sum n_j)$. $\square$

**Interpretation.** Integrating out the free field turns the imaginary-coupling ATFT into the grand-canonical ensemble of a **multi-component lattice Coulomb gas**. The charges are the affine roots $\beta\vec\alpha_j$, the fugacities $a^2zn_j$ are positive, and the Boltzmann weights are real. The "complex action" is a property of the field representation, not of the statistical model: the natural observables (vertex operators) have positive, real correlation sums.

This is the multi-component analogue of the well-known equivalence of sine-Gordon theory with the two-component Coulomb gas. Here, however, it holds without any reality of the potential. It closes H2.

The sign of $V$ in (2.1) is essential: the opposite sign gives fugacities of alternating sign. Appendix A shows that this sign choice decides whether the partition function is positive.

### 5.4 Hamiltonian structure: a Krein space

Let $\mathcal H=L^2(\mathbb R^{r(N_x-1)})$ be the space of wavefunctions of a time slice, $P$ the field parity $(P\psi)(\vec\phi)=\psi(-\vec\phi)$, and $K$ complex conjugation in the field representation. The transfer matrix is

$$
T=e^{-\frac a2\mathcal V}\,\mathcal K\,e^{-\frac a2\mathcal V},\qquad (\mathcal V\psi)(\vec\phi)=\sum_xa\,V(\vec\phi_x)\psi(\vec\phi), \tag{5.3}
$$

where $\mathcal K$ is the (real, symmetric, positive) Gaussian kernel of the kinetic and spatial-gradient terms. With this definition, $Z_\Lambda=\operatorname{Tr}T^{N_t}$.

**Theorem C [Theorem].**

1. $T$ is trace class.
2. $T^\dagger=PTP$. That is, $T$ is self-adjoint with respect to the indefinite inner product $[\psi,\chi]=\langle\psi,P\chi\rangle$, so $(\mathcal H,[\cdot,\cdot])$ is a Krein space with fundamental symmetry $P$.
3. The antilinear operator $\Theta=PK$ commutes with $T$.
4. The spectrum of $T$, with algebraic multiplicities, is invariant under complex conjugation.

*Proof.*

1. $\mathcal K$ is the heat kernel of a harmonic-oscillator Hamiltonian, because the Dirichlet gradient term is a positive-definite quadratic form, so it is trace class. The factors $e^{-\frac a2\mathcal V}$ are bounded by Theorem A. Hence $T$ is trace class.
2. $\mathcal K$ is real and $P$-invariant. Also $\overline{V(\vec\phi)}=V(-\vec\phi)$, so $(e^{-\frac a2\mathcal V})^\dagger=e^{-\frac a2\overline{\mathcal V}}=Pe^{-\frac a2\mathcal V}P$. Therefore $T^\dagger=PTP$.
3. $KTK=\overline T=P T P$ by the same identity, so $\Theta T\Theta^{-1}=T$.
4. $\Theta$ maps a Jordan chain for $\lambda$ to a Jordan chain for $\bar\lambda$. $\square$

The Euclidean counterpart is a twisted reflection positivity. Write $\vartheta$ for time reflection. Since $\overline{V(-\vec\phi)}=V(\vec\phi)$, the interacting measure satisfies

$$
\int\overline{(\vartheta PF)}\,F\,e^{-S}=\langle G,PG\rangle_{\rm OS,0},\qquad G=Fe^{-S_+},
$$

for $F$ supported at positive times, where $\langle\cdot,\cdot\rangle_{\rm OS,0}$ is the (positive) Osterwalder–Schrader form of the free field. The OS reconstruction therefore yields a Krein space with fundamental symmetry $P$, not a Hilbert space. This is the precise form of the statement "the theory is non-unitary", and it closes H4.

### 5.5 The vacuum energy is real

**Theorem D [Theorem].** The spectral radius $\rho(T)$ is positive and is an eigenvalue of $T$. Consequently the lattice vacuum energy $E_0=-a^{-1}\ln\rho(T)$ is real, and every energy $E=-a^{-1}\ln\lambda$ with $\lambda\in\sigma(T)\setminus\{0\}$ satisfies $\operatorname{Re}E\ge E_0$.

*Proof.* By Lidskii's theorem,

$$
\operatorname{Tr}T^N=\sum_\lambda m_\lambda\lambda^N,
$$

where $m_\lambda$ are the algebraic multiplicities. By Theorem B this equals $Z_\Lambda(N_t=N)>0$ for every $N\ge1$. In particular $T$ is not quasi-nilpotent, so $\rho(T)>0$.

Consider

$$
f(w)=\sum_{N\ge1}\operatorname{Tr}(T^N)\,w^N=\sum_\lambda m_\lambda\frac{\lambda w}{1-\lambda w}.
$$

It has simple poles at $w=1/\lambda$ with nonzero residues, and no cancellation is possible because distinct $\lambda$ give distinct poles. So its radius of convergence is exactly $1/\rho(T)$. The coefficients of $f$ are positive, so by Pringsheim's theorem the point $w=1/\rho(T)$ on the circle of convergence is a singularity. The only singularities are at $w=1/\lambda$, hence $\rho(T)\in\sigma(T)$. $\square$

Theorem D closes H5 at the lattice level. It is sharp: it says nothing about excited levels, and these can indeed be complex.

### 5.6 An exactly solvable check

The zero-mode reduction of (5.1) is quantum mechanics on the field-space torus,

$$
H=\tfrac12\lvert\vec p\rvert^2-g\sum_jn_je^{i\vec\alpha_j\cdot\vec x},\qquad \vec p\in\text{root lattice},\quad g>0.
$$

This $H$ has the same $\Theta$-symmetry, and its partition function has the same positive Coulomb-gas expansion. Diagonalizing in a momentum basis (cutoff $P=16$, converged) gives:

| algebra | $g$ | lowest level | complex levels among the lowest 40 | $Z(T)$ for $T=0.5,2,8$ |
|---|---|---|---|---|
| $a_2$ | 0.5 | $-0.1625$ (real) | 12 | 3.66, 1.48, 3.67 |
| $a_2$ | 1.0 | $-0.8077$ (real) | 14 | 3.91, 5.12, 640.2 |
| $c_2$ | 0.5 | $-0.4823$ (real) | 20 | 6.32, 2.88, 47.4 |
| $c_2$ | 1.0 | $-1.7542$ (real) | 20 | 6.90, 34.3, $1.24\times10^6$ |

As Theorems B and D require, $Z(T)>0$ and the lowest level is real. Excited levels include complex-conjugate pairs: $\Theta$ is broken in parts of the finite-volume spectrum, while the vacuum is protected.

With the *opposite* sign of the potential, $Z(T)$ for $a_2$ at $g=1$ is negative at $T=2$ ($-0.72$), and the lowest level is the complex pair $0.410\pm1.327\,i$. For $c_2$ the two signs are equivalent, because a shift of $\vec x$ flips the sign of every $e^{i\vec\alpha_j\cdot\vec x}$ when $\sum_jn_j$ is even. Positivity of the fugacity is the mechanism that protects the vacuum.

---

## 6. Foundations II: the continuum limit

### 6.1 Collapse thresholds

In the Coulomb-gas picture, ultraviolet divergences come from clusters of charges collapsing to a point. Consider a cluster $S$ (a finite multiset of affine roots, possibly with repeats) shrinking uniformly by a factor $s\to0$, with continuum covariance $C(x)\simeq-\frac{1}{2\pi}\ln\lvert x\rvert$. Its weight scales as

$$
s^{\,2(\lvert S\rvert-1)}\cdot s^{-\frac{\beta^2}{4\pi}\left(\sum_{k\in S}\lvert\vec\alpha_k\rvert^2-\lvert\sum_{k\in S}\vec\alpha_k\rvert^2\right)} .
$$

Define the **first collapse threshold**

$$
\beta^2_{\rm UV}(\hat g)=\min_S\ \frac{8\pi(\lvert S\rvert-1)}{\sum_{k\in S}\lvert\vec\alpha_k\rvert^2-\lvert\sum_{k\in S}\vec\alpha_k\rvert^2},
$$

where the minimum runs over clusters with a positive denominator. Some values:

- For a pair $(\vec\alpha_i,\vec\alpha_j)$ with $\vec\alpha_i\cdot\vec\alpha_j<0$ the bound is $4\pi/\lvert\vec\alpha_i\cdot\vec\alpha_j\rvert$.
- For $a_1^{(1)}$ the pair $(\vec\alpha,-\vec\alpha)$ gives $2\pi$. This is the familiar sine-Gordon bound $\beta_{\rm SG}^2<4\pi$, since $\beta_{\rm SG}^2=2\beta^2$.
- For simply-laced $\hat g$ the full neutral cluster $\sum_jn_j\vec\alpha_j=0$ gives $4\pi(h-1)/h$, and this is the minimum. A cluster whose net charge is not a multiple of $\sum_jn_j\vec\alpha_j$ has $\lvert\sum_k\vec\alpha_k\rvert^2\ge2$, which gives a bound of at least $4\pi$; a cluster with net charge zero is $k$ copies of the full neutral cluster and gives $4\pi(1-1/kh)$.
- For $b_n^{(1)}$, $c_n^{(1)}$, $g_2^{(1)}$ and $a_2^{(2)}$ the minimum is $4\pi$, attained by a pair of adjacent roots with $\vec\alpha_i\cdot\vec\alpha_j=-1$.
- For every $\hat g\ne a_1^{(1)}$, $8\pi/3\le\beta^2_{\rm UV}\le4\pi$. The upper bound comes from a longest root and any neighbour, whose inner product is $-1$. For the lower bound, the denominator is at most $\sum_k\lvert\vec\alpha_k\rvert^2\le2\lvert S\rvert$, so clusters with $\lvert S\rvert\ge3$ give at least $4\pi(1-1/\lvert S\rvert)\ge8\pi/3$. Two distinct simple affine roots have $\lvert\vec\alpha_i\cdot\vec\alpha_j\rvert\le1$ unless $\vec\alpha_i=-\vec\alpha_j$, which happens only for $a_1^{(1)}$; two copies of one root give a negative denominator. In particular $\beta^2_{\rm UV}>2\pi$ for every algebra except $a_1^{(1)}$, where $\beta^2_{\rm UV}=2\pi$.

With pointwise Wick ordering (Section 6.2), Hepp-sector power counting shows that the continuum integrand of each term of (5.2) is integrable near every *bulk* collapse exactly when $\beta^2<\beta^2_{\rm UV}$. In a Hepp sector the iterated integral at a node $S$ converges if and only if $2(\lvert S\rvert-1)>\frac{\beta^2}{4\pi}\big(\sum_{k\in S}\lvert\vec\alpha_k\rvert^2-\lvert\sum_{k\in S}\vec\alpha_k\rvert^2\big)$. Proposition 6.7 gives the "only if" direction rigorously. Proposition 6.13 gives the "if" direction at sketch level, including the boundary layer for clusters near a Dirichlet wall. With a site-independent Wick constant even the one-charge term diverges for $\beta^2\ge2\pi$ (Proposition 6.3).

**The thresholds as resonances.** Let $\hat g$ be simply laced, write $x=\beta^2/2\pi$ for the scaling dimension of $e^{i\beta\vec\alpha_j\cdot\vec\phi}$, and let $g$ be the physical coupling (Section 6.2). By dimensional analysis the bulk vacuum-energy density is proportional to $g^{2/(2-x)}$, while a cluster of $k$ copies of the neutral cluster first contributes to it at order $g^{kh}$. The two have the same dimension exactly when $x=2-2/(kh)$, that is $\beta^2=4\pi(1-1/kh)$, which are the neutral thresholds above. For sine-Gordon ($h=2$) these are the poles of the exact bulk energy $-\frac{m^2}4\tan\frac{\pi\xi}2$, with $\xi=\beta^2/(4\pi-\beta^2)$, at odd $\xi$ [Established; Al.B. Zamolodchikov 1995]. By Theorem B the coefficient of $g^{kh}$ is positive and cannot vanish accidentally. For the other simply-laced algebras we therefore expect the exact bulk energy, written in terms of the bare coupling, to be singular at these couplings; no such resonance falls at $\beta^2=2\pi$ when $h\ge3$ [heuristic]. This can be tested with the integrable lattice regularizations of open problem 8, in which the thresholds of $a_{N-1}^{(1)}$ sit at $q^{kN}=-1$.

### 6.2 Wick ordering and the statement

Fix $L_t=aN_t$ and $L_x=aN_x$, and let $\Omega=(\mathbb R/L_t\mathbb Z)\times(0,L_x)$ be the continuum cylinder, with distances periodic in $t$. Let $G_\Omega$ be its Dirichlet Green function, $-\Delta G_\Omega=\delta$. Write $q_j=\beta\vec\alpha_j$ and $\kappa=\beta^2/4\pi$.

**Pointwise Wick ordering.** Because of the Dirichlet conditions, $C_a(u,u)$ depends on the distance of $u$ from the boundary. We therefore Wick order site by site,

$$
:e^{iq\cdot\vec\phi_u}:\;=e^{\frac12\lvert q\rvert^2C_a(u,u)}\,e^{iq\cdot\vec\phi_u},
$$

which replaces $z$ by $z_{a,j}(u)=z\,e^{\frac12\beta^2\lvert\vec\alpha_j\rvert^2C_a(u,u)}$ in (5.1). These fugacities depend on the root and on the spatial coordinate of $u$ only, so Theorems A–D apply to them (Section 5.1). External vertex operators are Wick ordered in the same way.

**What pointwise Wick ordering means.** For fixed interior $u$, $C_a(u,u)=\frac1{2\pi}\ln\frac{R_\Omega(u)}a+c_0+o(1)$ as $a\to0$. Here $R_\Omega(u)$ is the conformal radius of $\Omega$ at $u$, defined by $G_\Omega(u,v)=\frac1{2\pi}\ln\frac{R_\Omega(u)}{\lvert u-v\rvert}+o(1)$ as $v\to u$, and $c_0$ is a lattice constant [Established; Lawler–Limic 2010]. Normalize vertex operators instead at the lattice scale, by the $u$-independent factor $e^{\frac12\lvert q_j\rvert^2(\frac1{2\pi}\ln\frac1a+c_0)}$. Relative to that normalization, pointwise Wick ordering is the position-dependent coupling

$$
g_j(u)=z\,R_\Omega(u)^{x_j},\qquad x_j=\kappa\lvert\vec\alpha_j\rvert^2,
$$

where $x_j$ is the scaling dimension of $e^{i\beta\vec\alpha_j\cdot\vec\phi}$. Near a wall $R_\Omega(u)\simeq2d$, so the coupling vanishes like $(2d)^{x_j}$, consistently with (G3). Three consequences follow.

- The walls are harmless in Theorem E and Proposition 6.13 because the interaction is switched off there.
- A site-independent Wick constant keeps the coupling constant up to the wall. Heuristically, the boundary free energy then scales as $g^{1/(2-x)}$, which has the dimension of the first-order term exactly at $x=1$, that is at $\beta^2=2\pi$ for a long root, whatever the algebra. This is the threshold of Proposition 6.3.
- The finite-volume theory of Theorem E is the ATFT with this coupling profile, not with a constant coupling. Comparisons with exact S-matrices, TCSA or integrable lattice data must use the profile. In Conjecture 6.2 the infinite-volume limit is taken at fixed physical coupling, that is, with the fugacity of root $j$ scaled like $L^{-x_j}$, or equivalently with Wick ordering at a fixed reference scale.

**Renormalized partition function.** With these fugacities in both terms of (5.1), $Z_\Lambda=e^{-a^2\sum_u\sum_jn_jz_{a,j}(u)}\,\Xi_a$, where

$$
\Xi_a=\sum_{N\ge0}\frac{(a^2z)^N}{N!}\sum_{j_1,\dots,j_N}\Big(\prod_kn_{j_k}\Big)\sum_{u_1,\dots,u_N\in\Lambda}W_a(u),\qquad
W_a(u)=\exp\Big(-\tfrac12\sum_{k\ne l}q_{j_k}\cdot q_{j_l}\,C_a(u_k,u_l)\Big). \tag{6.1}
$$

Pointwise Wick ordering removes exactly the diagonal terms $k=l$; pairs $k\ne l$ at the same site are kept. The prefactor is field-independent (a vacuum-energy counterterm) and cancels from all correlation functions. For external charges $\beta\vec w_m$ at sites $v_m$, let $\Xi_a(w,v)$ be the same expansion with these charges added (not summed over) and with $W_a$ including all pairs. Then $\big\langle\prod_m:e^{i\beta\vec w_m\cdot\vec\phi_{v_m}}:\big\rangle=\Xi_a(w,v)/\Xi_a$. The continuum versions replace $C_a$ by $G_\Omega$ and $a^2\sum_u$ by $\int_\Omega$.

**Proposition 6.1 [Theorem for $\beta^2<2\pi$ (Theorem E); Conjecture for $2\pi\le\beta^2<\beta^2_{\rm UV}$, reduced to Conjecture 6.14].** Let $0<\beta^2<\beta^2_{\rm UV}(\hat g)$, fix the finite volume, and use pointwise Wick ordering.

1. $\Xi_a$ converges as $a\to0$ to the continuum version of (6.1). Every term is positive and finite, and the series converges.
2. The same holds for $\Xi_a(w,v)$, with fixed distinct interior points $v_m$ (approximated by lattice sites), provided that for every $m$ and every non-empty multiset $S$ of affine roots
   $$
   2\lvert S\rvert>\kappa\Big(\sum_{k\in S}\lvert\vec\alpha_k\rvert^2+\lvert\vec w_m\rvert^2-\big\lvert\vec w_m+\textstyle\sum_{k\in S}\vec\alpha_k\big\rvert^2\Big). \tag{6.2}
   $$
   Hence the normalized vertex-operator correlation functions converge.

Condition (6.2) says that no cluster collapsing onto an external charge is divergent. Each ingredient of the proposition is necessary: Proposition 6.3 shows that part 1 fails for a site-independent Wick constant when $\beta^2\ge2\pi$, that $Z_\Lambda$ itself tends to zero, and that part 2 fails whenever (6.2) is violated with $\lvert S\rvert=1$. Proposition 6.7 shows that part 1 fails for $\beta^2\ge\beta^2_{\rm UV}$.

**Conjecture 6.2 [Conjecture].** For $0<\beta^2<4\pi$ (in our normalization, the range in which the vertex operators $e^{i\beta\vec\alpha_j\cdot\vec\phi}$ of the long roots are relevant), the infinite-volume Schwinger functions of vertex operators exist and decay exponentially, so there is a mass gap. For $\beta^2_{\rm UV}\le\beta^2<4\pi$ this requires additional counterterms generated by the collapsing clusters (Proposition 6.7 shows that some are needed). Only finitely many cluster types diverge in this range, since a divergent cluster has $\lvert S\rvert\le1/(1-\kappa)$, and we conjecture that the corresponding counterterms suffice, as in the Benfatto–Gallavotti–Nicolò treatment of sine-Gordon [BGN 1982]. For simply-laced $\hat g$ the counterterms can be specified. By Section 6.1, a bulk cluster with non-zero net charge diverges only at $\beta^2\ge4\pi$. So in this range every divergent bulk cluster consists of $k$ copies of the neutral cluster, and the bulk counterterms are field-independent: vacuum-energy terms at orders $z^{kh}$, at the resonances described in Section 6.1. Boundary terms near the Dirichlet walls are not covered by this count. For sine-Gordon the same count puts the first non-neutral divergence at $\beta_{\rm SG}^2=8\pi$.

The evidence for this conjecture is the sine-Gordon case ($a_1^{(1)}$), where it is established, together with the existence and internal consistency of the exact S-matrices throughout this range.

**Lattice Green function.** We use five standard estimates. Distances are physical, the constants $D,K$ depend only on $\Omega$, and $m$ denotes a lattice distance.

- **(G1)** $C_a(u,v)\le\frac1{2\pi}\ln\frac{D}{\max(\lvert u-v\rvert,a)}+K$ for all sites $u,v$, including $u=v$.
- **(G2)** $C_a(u_a,v_a)\to G_\Omega(x,y)$ whenever $u_a\to x$, $v_a\to y$ and $x\ne y$.
- **(G3)** If $u$ is at lattice distance $m$ from the nearer Dirichlet row, then $C_a(u,u)\le\frac1{2\pi}\ln(2m)+K$. If moreover $m\le\min(N_x,N_t)/4$, then $C_a(u,u)\ge\frac1{2\pi}\ln m-K$.
- **(G4)** $C_a(u,v)\le\frac1{4\pi}\ln\big(1+4d_ud_v/\max(\lvert u-v\rvert,a)^2\big)+K$, where $d_u$ is the distance from $u$ to the nearer Dirichlet row.
- **(G5)** For the lattice ball $B$ of radius $\rho$ about $x$, $G_B(x,y)=\frac1{2\pi}\ln\frac{\rho}{\max(\lvert x-y\rvert,a)}+O(1)$ for $\lvert x-y\rvert\le\rho/2$, and $G_B(x,y)\le c$ for $\lvert x-y\rvert\ge\rho/2$.

These follow from domain monotonicity, from the image formula on the half-cylinder (the half-plane Green function $\mathfrak a(u-\bar v)-\mathfrak a(u-v)$, with $\bar v$ the reflection of $v$ in the wall, plus summable periodic images), from comparison with a disc of radius $m$, and from the potential-kernel asymptotics $\mathfrak a(x)=\frac1{2\pi}\ln\lvert x\rvert+O(1)$ [Lawler–Limic 2010, Thm 4.4.4, rescaled to $-\Delta$]. (G4) is the image bound $G(x,y)=\frac1{4\pi}\ln(1+4d_xd_y/\lvert x-y\rvert^2)$ for the half-plane, and (G5) is the standard asymptotics of Green functions of lattice balls [Lawler–Limic 2010, Ch. 6].

### 6.3 The corrections are necessary

**Proposition 6.3 [Theorem].**

1. **The constant term.** If the constant term of (5.1) carries the Wick-ordered fugacities, the prefactor $e^{-a^2\sum_u\sum_jn_jz_{a,j}(u)}$ tends to $0$. So $Z_\Lambda$ has no non-trivial continuum limit; the object with a limit is $\Xi_a$.
2. **Site-independent Wick ordering diverges at the walls for $\beta^2\ge2\pi$.** Let $z_{a,j}=z\,e^{\frac12\beta^2\lvert\vec\alpha_j\rvert^2c_a}$ with $c_a$ independent of $u$ and $c_a\ge\frac1{2\pi}\ln(L/a)-K$ for some $L,K>0$. Any constant that makes bulk one-point functions converge has this form, for example $c_a=C_a(u_*,u_*)$ at a fixed interior point, by (G3). Then for every affine algebra and every $\beta^2\ge2\pi$, the series (5.2) with these fugacities tends to $+\infty$ as $a\to0$.
3. **External charges.** Suppose $\beta^2\vec\alpha_j\cdot\vec w\le-4\pi$ for some $j$; this violates (6.2) with $S=\{\vec\alpha_j\}$. Then, in any range where $\Xi_a$ converges, the normalized correlation function containing $:e^{i\beta\vec w\cdot\vec\phi(v)}:$ diverges. For every $\beta>0$ this happens for some $\vec w$, for example $\vec w=-m\vec\alpha_j$ with $m$ large.

*Proof.*

1. By (G3), sites with $m\ge\min(N_x,N_t)/4$ have $C_a(u,u)\ge\frac1{2\pi}\ln\frac{\min(N_x,N_t)}4-K$. They occupy a fixed fraction of $\Omega$, so $a^2\sum_uz_{a,j}(u)\ge c\,z\,a^{-\kappa\lvert\vec\alpha_j\rvert^2}\to\infty$.
2. By Theorem B every term of (5.2) is positive, so the series is at least its $N=1$ term. Keep only a long simple root $\vec\alpha_{j_0}$, which every affine algebra has; then $\lvert q_{j_0}\rvert^2=2\beta^2$. By (G3), a site at lattice distance $m\le N_x/2$ from a wall has $C_a(u,u)\le\frac1{2\pi}\ln(2m)+K$. There are $L_t/a$ such sites for each $m$, so the $N=1$ term is at least
   $$
   a^2z\,n_{j_0}\sum_ue^{\beta^2(c_a-C_a(u,u))}\ge z\,n_{j_0}e^{-2K\beta^2}L_t\,a\Big(\frac L{2a}\Big)^{\beta^2/2\pi}\sum_{m\le M}m^{-\beta^2/2\pi}=c\,a^{1-\beta^2/2\pi}\sum_{m\le M}m^{-\beta^2/2\pi},
   $$
   with $M\asymp L_x/a$. This tends to $\infty$ for $\beta^2>2\pi$, and grows like $\ln(1/a)$ at $\beta^2=2\pi$.
3. The $N=1$ term of $\Xi_a(w,v)$ contains $a^2zn_j\sum_ue^{-\beta^2\vec\alpha_j\cdot\vec w\,C_a(u,v)}$. By Fatou and (G2) its limit inferior is at least $zn_j\int_\Omega e^{\beta^2\lvert\vec\alpha_j\cdot\vec w\rvert G_\Omega(x,v)}dx\ge c\int_{\lvert x-v\rvert<r}\lvert x-v\rvert^{-\beta^2\lvert\vec\alpha_j\cdot\vec w\rvert/2\pi}dx=\infty$. $\square$

The mechanism in part 2 is the image charge. A single charge $q$ at distance $d$ from the wall carries the weight $(2d)^{-\lvert q\rvert^2/4\pi}$ relative to the bulk. For a cluster collapsing onto the wall, the interior and image pair factors cancel, so a single long root is the worst case. The resulting threshold $2\pi$ is algebra-independent, and by Section 6.1 it lies strictly below $\beta^2_{\rm UV}$ for every algebra except $a_1^{(1)}$, where the two coincide.

### 6.4 Proof below $2\pi$

Write $\gamma_k=\lvert q_k\rvert^2/4\pi=\kappa\lvert\vec\alpha_{j_k}\rvert^2\le2\kappa=\beta^2/2\pi$.

**Theorem E [Theorem].** Use pointwise Wick ordering, let $0<\beta^2<2\pi$, and let the external charges $\beta\vec w_m$ sit at fixed distinct interior points $v_m$ with $\beta^2\lvert\vec w_m\rvert^2<4\pi$. Set
$$
\gamma=\max\Big(\frac{\beta^2}{2\pi},\ \max_m\frac{\beta^2\lvert\vec w_m\rvert^2}{4\pi}\Big)<1.
$$
Then parts 1 and 2 of Proposition 6.1 hold, for every affine algebra. Uniformly in $a$, the $N$-th term of $\Xi_a(w,v)$ is at most $K_0(Kz)^N(N+M)!^{\gamma/2}/N!$.

For sine-Gordon, $\beta^2_{\rm UV}=2\pi$, so Theorem E covers the whole range of Proposition 6.1, consistent with [Fröhlich 1976].

**Lemma 6.4 (lattice Newton theorem) [Theorem].** Let $s_1,\dots,s_n$ be distinct interior sites carrying charges $Q_i\in\mathbb R^r$. Let $d_0=\min(L_t,L_x)/4$, let $\delta_i=\min\big(d_0,\min_{i'\ne i}\lvert s_i-s_{i'}\rvert\big)$, and let $\rho_i=(\delta_i-a)/2-a/8$. Let $B_i$ be the set of interior sites $y$ with $\lvert y-s_i\rvert\le\rho_i$ (empty if $\rho_i<0$), and let $\mu_i$ be the exit law from $B_i$ of simple random walk started at $s_i$ and killed on the Dirichlet rows (so $\mu_i=\delta_{s_i}$ if $B_i$ is empty). Then
$$
\sum_{i\ne i'}Q_i\cdot Q_{i'}\,C_a(s_i,s_{i'})\ \ge\ -\sum_i\lvert Q_i\rvert^2\,C_a(\mu_i,\mu_i),\qquad
C_a(\mu_i,\mu_i)\le\min\Big(C_a(s_i,s_i),\ \max_{y\notin B_i}C_a(s_i,y)\Big).
$$

*Proof.*

1. For $y\notin B_i$, the function $C_a(\cdot,y)$ is discrete harmonic on $B_i$ and vanishes on the Dirichlet rows. By optional stopping, $\sum_{y'}\mu_i(y')C_a(y',y)=C_a(s_i,y)$.
2. The support of $\mu_i$ lies outside $B_i$, within distance $\rho_i+a$ of $s_i$. Since $\rho_i+\rho_{i'}+a<\lvert s_i-s_{i'}\rvert$, each point of $\operatorname{supp}\mu_{i'}$ lies outside $B_i$. Applying step 1 twice gives $C_a(\mu_i,\mu_{i'})=C_a(s_i,s_{i'})$ for $i\ne i'$.
3. $C_a$ is positive definite, so $\sum_{i,i'}Q_i\cdot Q_{i'}\,C_a(\mu_i,\mu_{i'})\ge0$. With step 2 this gives the inequality.
4. By step 1, $C_a(\mu_i,\mu_i)=\sum_y\mu_i(y)\,C_a(s_i,y)$, with $\mu_i$ a sub-probability measure on sites outside $B_i$. This is at most $\max_{y\notin B_i}C_a(s_i,y)$. It is also at most $C_a(s_i,s_i)$, because $C_a(s,y)=P_y(\text{hit }s)\,C_a(s,s)\le C_a(s,s)$. $\square$

**Lemma 6.5 (multiply occupied sites) [Theorem].** For a configuration $u$ of all $N+M$ charges, group the charges by occupied site. Then
$$
W_a(u)\le K_1^{N+M}\prod_{s\ \text{singly occupied}}\Big(\frac D{\delta_s}\Big)^{\gamma_{k(s)}}\ \prod_{s\ \text{multiply occupied}}\ \prod_{k\ \text{at}\ s}\Big(\frac Da\Big)^{\gamma_k},
$$
where $\delta_s$ is as in Lemma 6.4 and $k(s)$ is the charge at a singly occupied site.

*Proof.*

1. Let $Q_s$ be the total charge at site $s$. Split the exponent of $W_a$ into on-site pairs and inter-site pairs, and apply Lemma 6.4 to the inter-site part. The contribution of site $s$ is then at most
   $$
   \tfrac12\sum_{k\ \text{at}\ s}\lvert q_k\rvert^2C_a(s,s)-\tfrac12\lvert Q_s\rvert^2\big(C_a(s,s)-C_a(\mu_s,\mu_s)\big).
   $$
2. At a singly occupied site this equals $\tfrac12\lvert q_k\rvert^2C_a(\mu_s,\mu_s)$. By Lemma 6.4 and (G1) it is at most $\gamma_k\big(\ln\frac D{\max(\rho_s,a)}+2\pi K\big)$. Finally, $\max(\rho_s,a)\ge\delta_s/4$.
3. At a multiply occupied site, $C_a(\mu_s,\mu_s)\le C_a(s,s)$, so the contribution is at most $\tfrac12\sum_k\lvert q_k\rvert^2C_a(s,s)$. Apply (G1). $\square$

**The continuum majorant.** Embed a lattice configuration as $x\in\Omega^N$, with $x_k$ in the cell (square of side $a$) of $u_k$. Let $\delta_k(x)$ be the distance from $x_k$ (or from $v_m$, for an external charge) to the nearest other point among $x_1,\dots,x_N,v_1,\dots,v_M$, capped at $d_0$. Then:

- a lone charge has $\delta_k(x)\le\delta_s+\sqrt2a\le3\delta_s$;
- a charge at a multiply occupied site, or an internal charge at the site of some $v_m$, has $\delta_k(x)\le2a$.

Since $D/\delta_k(x)\ge1$, every $\gamma_k$ may be replaced by $\gamma$, and Lemma 6.5 gives

$$
W_a(u(x))\le K_2^{N+M}\prod_{k=1}^{N+M}\Big(\frac D{\delta_k(x)}\Big)^{\gamma}. \tag{6.3}
$$

**Lemma 6.6 (nearest-neighbour integral) [Theorem].** For $0\le\gamma<1$ and fixed distinct $v_1,\dots,v_M\in\Omega$,
$$
\int_{\Omega^N}\prod_{k=1}^{N+M}\delta_k(x)^{-\gamma}\,dx\le c(v)\,C^N(N+M)!^{\gamma/2},
$$
with $C$ depending only on $\gamma$ and $\Omega$.

*Proof.* Induct by adding one integrated point $x$ to a set $X'$ of $n$ points (integrated or fixed). Write $\delta'_k$ for the nearest-neighbour distances within $X'$. Adding $x$ changes $\delta_k$ to $\min(\delta'_k,\lvert x-x_k\rvert)$ and adds the factor $\operatorname{dist}(x,X')^{-\gamma}$.

1. Let $A=\{k:\lvert x-x_k\rvert<\delta'_k\}$. For distinct $k,k'\in A$, $\lvert x_k-x_{k'}\rvert\ge\max(\delta'_k,\delta'_{k'})$ is the strictly longest side of the triangle $(x,x_k,x_{k'})$. So the angle at $x$ exceeds $60^\circ$, and $\lvert A\rvert\le5$.
2. Let $k_1\in A$ be the point of $A$ closest to $x$. For the other $k\in A$, $\delta'_k\le\lvert x_k-x_{k_1}\rvert\le2\lvert x_k-x\rvert$.
3. Hence $\prod_{k\in X'}\delta_k^{-\gamma}\le2^{4\gamma}\prod_k\delta_k'^{-\gamma}\Big(1+\sum_k\mathbb 1[\lvert x-x_k\rvert<\delta'_k]\big(\delta'_k/\lvert x-x_k\rvert\big)^\gamma\Big)$.
4. Multiply by $\operatorname{dist}(x,X')^{-\gamma}$ and integrate over $x$. On $\lvert x-x_k\rvert<\delta'_k/2$ the distance to $X'$ is $\lvert x-x_k\rvert$, which gives $c\,\delta_k'^{2-\gamma}$; this uses $2\gamma<2$. On the annulus $\delta'_k/2\le\lvert x-x_k\rvert<\delta'_k$ the ratio factor is at most $2^\gamma$, and each $x$ lies in at most 5 such annuli by step 1.
5. $\int_\Omega\operatorname{dist}(x,X')^{-\gamma}dx\le c\,n^{\gamma/2}\lvert\Omega\rvert^{1-\gamma/2}$. To see this, split $\Omega$ into Voronoi cells, replace each cell by the centred disc of equal area (rearrangement), and use concavity.
6. The discs of radius $\delta'_k/2$ are disjoint and lie in the $d_0/2$-neighbourhood of $\Omega$. So $\sum_k\delta_k'^{2-\gamma}\le n^{\gamma/2}\big(\sum_k\delta_k'^2\big)^{1-\gamma/2}\le c\,n^{\gamma/2}$.
7. Each step therefore costs at most $C\,n^{\gamma/2}$. $\square$

*Proof of Theorem E.*

1. Write the $N$-th term of $\Xi_a(w,v)$ as $\frac{z^N}{N!}\sum_{j}\prod_kn_{j_k}\int W_a(u(x))\,dx$, integrating over the union of cells.
2. By (6.3) and Lemma 6.6, it is at most $K_0(Kz)^N(N+M)!^{\gamma/2}/N!$ uniformly in $a$. Since $\gamma<1$, this is summable in $N$.
3. For almost every $x\in\Omega^N$ (distinct points, distinct from the $v_m$), $u(x)\to x$ and $W_a(u(x))\to W(x)$ by (G2). The majorant in (6.3) is independent of $a$ and integrable by Lemma 6.6. Dominated convergence gives convergence of each term, and Tannery's theorem, using step 2, gives convergence of the series.
4. The limiting terms are positive and finite, and $\lim\Xi_a\ge1$, so the normalized correlation functions converge. $\square$

### 6.5 The threshold is sharp

**Proposition 6.7 [Theorem].** With pointwise Wick ordering and $\beta^2\ge\beta^2_{\rm UV}(\hat g)$, $\Xi_a\to\infty$ as $a\to0$. So additional counterterms are necessary at and above the first collapse threshold.

*Proof.*

1. Pick a cluster $S=(j_1,\dots,j_n)$ with $\kappa D(S)\ge2(n-1)$, where $D(S)=\sum_k\lvert\vec\alpha_{j_k}\rvert^2-\lvert\sum_k\vec\alpha_{j_k}\rvert^2>0$. Every term of (6.1) is positive, so $\Xi_a$ is at least the $n$-th term restricted to the species sequence of $S$ and to configurations in $B^n$, for a small closed disc $B\subset\Omega$.
2. By Fatou and (G2), the limit inferior of this restricted term is at least $\frac{z^n}{n!}\prod_kn_{j_k}\int_{B^n}W$.
3. On $B$, $G_\Omega(x,y)=\frac1{2\pi}\ln\frac1{\lvert x-y\rvert}+O(1)$. So $W\ge c\prod_{k<l}\lvert x_k-x_l\rvert^{-\kappa D_{kl}}$, with $D_{kl}=-2\vec\alpha_{j_k}\cdot\vec\alpha_{j_l}$ and $\sum_{k<l}D_{kl}=D(S)$.
4. Restrict to uniform collapse: $x_k=x_0+s\,y_k$, with $y$ in a small open set of unit-norm relative configurations. The measure is $c\,s^{2n-3}\,ds$ and the integrand is at least $c\,s^{-\kappa D(S)}$. Since $2n-3-\kappa D(S)\le-1$, the integral is infinite. $\square$

### 6.6 The window $2\pi\le\beta^2<\beta^2_{\rm UV}$

By Section 6.1 this window is non-empty for every algebra except $a_1^{(1)}$. The natural extension of Fröhlich's method does not reach it.

**Proposition 6.8 [Theorem].** Let $\vec\alpha_j$ be a long root and $\beta^2\ge2\pi$.

1. **Per-charge (Onsager) bounds.** Any bound on $W_a$ that depends on the charges only through their moduli $\lvert q_k\rvert$ also bounds the gas with charges $\pm q_j$. For that gas the two-charge term satisfies $\liminf_{a\to0}a^4\sum_{u_1,u_2}e^{\lvert q_j\rvert^2C_a(u_1,u_2)}\ge\iint_{\Omega^2}e^{2\beta^2G_\Omega(x,y)}dx\,dy=\infty$. So no such bound is integrable uniformly in $a$.
2. **Moment bounds.** Let $M_j=\int_\Omega:e^{iq_j\cdot\vec\phi}:\,dx$. Then $\mathbb E\lvert M_j\rvert^2=\iint e^{2\beta^2G_\Omega(x,y)}dx\,dy=\infty$, and the lattice second moments diverge as $a\to0$ by Fatou and (G2). So any bound on the $N$-th term through moments of the $M_j$ (Hölder in the species, or chaos-based arguments) fails.

*Proof.* Near the diagonal, $e^{2\beta^2G_\Omega(x,y)}\ge c\lvert x-y\rvert^{-\beta^2/\pi}$, which is not integrable in two dimensions for $\beta^2\ge2\pi$. The lattice statements follow by Fatou and (G2). $\square$

Consistently with part 2, imaginary multiplicative chaos for a single root is constructed for $\lvert q\rvert^2<4\pi$ [Junnila–Saksman–Webb 2020; see also Lacoin–Rhodes–Vargas 2015 for complex chaos], that is for $\beta^2<2\pi$ on long roots. So the probabilistic route reaches the same range as Theorem E and no further.

Below $\beta^2_{\rm UV}$ neither failure reflects the true weights. Both arise because the bounds introduce the conjugate charge $-\vec\alpha_j$, which is not in the gas for rank $\ge2$: no simple affine root is the negative of another except in $a_1^{(1)}$. Closing the window therefore needs a genuinely vector-valued stability estimate that uses $\vec\alpha_j+\vec\alpha_k\ne0$.

**Estimate (S).** Let $0<\beta^2<\beta^2_{\rm UV}(\hat g)$ and use pointwise Wick ordering. Estimate (S) consists of two statements:

- **(S1)** There are $C<\infty$ and $\theta<1$ such that, uniformly in $a$, $\sum_{j\in J^N}\prod_kn_{j_k}\int W_a(u(x))\,dx\le C^N N!^{\theta}$.
- **(S2)** For each $N$, the functions $x\mapsto W_a(u(x))$ are uniformly integrable on $\Omega^N$ as $a\to0$.

The same statements are made with external charges satisfying (6.2). Given (S), the proof of Theorem E goes through with Vitali's theorem in place of dominated convergence, and proves Proposition 6.1 in the window.

**Status of (S).** The rest of this section establishes two things, the first at sketch level. First, (S2) holds throughout the window (Proposition 6.13, sketch). This means every term of the expansion has a finite continuum limit, including near the Dirichlet walls. Second, (S1) follows from a bound on products of nearest-neighbour distances of points in $\Omega$. That bound is proved below as Theorem 6.9. The version needed here also carries the intra-unit factors of (6.5); this is Conjecture 6.14. The mechanism is to smear adaptively chosen clusters of charges, called units, with a common radius. This keeps each unit's internal interactions exact and charges self-energy only at the scale that separates the unit from everything else.

**Theorem 6.9 (products of nearest-neighbour distances) [Theorem].** For a configuration $x\in\Omega^N$, let $\delta_k^{(m)}$ be the distance from $x_k$ to its $m$-th nearest other point, capped at $D$, with $\delta_k^{(m)}=D$ when $N-1<m$. For every $m\ge1$ and $0\le\gamma<2m/(m+1)$,

$$
\int_{\Omega^N}\prod_{k=1}^N\Big(\frac{D}{\delta_k^{(m)}}\Big)^{\gamma}dx\le C^N\,N!^{\gamma/2},
$$

with $C$ depending only on $m$, $\gamma$ and $\Omega$.

For $m=1$ this is Lemma 6.6. The exponent $\gamma/2$ is the natural one, since typical $m$-th-neighbour distances are of order $N^{-1/2}$. The range of $\gamma$ is the local-integrability range: a collapse of $m+1$ points costs $r^{-(m+1)\gamma}$ against a measure $r^{2m}$. The one-point insertion argument of Lemma 6.6 cannot reach $m\ge2$, because the singular cost belongs to a group of $m+1$ points rather than to a single point. The proof below therefore works with clusters.

*Proof.* Fix $\varepsilon\in(0,1]$, to be chosen at the end, and let $\ell=\varepsilon N^{-1/2}$ and $\alpha=2m-(m+1)\gamma>0$. On the cylinder every step uses the periodic metric, and for small $\varepsilon$ the rescaled circumference below exceeds $2$.

**Step 1: coarse scales.** Write $g_k=\big((\ell/\delta_k^{(m)})^{\gamma}-1\big)_+$. Then

$$
\Big(\frac{D}{\delta_k^{(m)}}\Big)^{\gamma}=\Big(\frac{D}{\max(\ell,\delta_k^{(m)})}\Big)^{\gamma}(1+g_k)\le\Big(\frac{D\sqrt N}{\varepsilon}\Big)^{\gamma}(1+g_k).
$$

Using $N^N\le e^NN!$, the product over $k$ is at most $\big(D^\gamma e^{\gamma/2}\varepsilon^{-\gamma}\big)^NN!^{\gamma/2}\prod_k(1+g_k)$. It therefore suffices to show $\int_{\Omega^N}\prod_k(1+g_k)\,dx\le e^N\lvert\Omega\rvert^N$.

**Step 2: factorization over clusters.** Rescale to $u=x/\ell$, so that $\Lambda=\Omega/\ell$ has area $N\lvert\Omega\rvert/\varepsilon^2$, and write $\tilde\delta_k=\delta_k^{(m)}/\ell$. Join two points when $\lvert u_i-u_j\rvert\le1$.

- On a connected component with at most $m$ points, every factor $1+g_k$ equals $1$, because each point has at most $m-1$ other points within distance $1$.
- On a component $C$ with at least $m+1$ points, each factor depends only on the points of $C$. If $\tilde\delta_k<1$, all $m$ nearest neighbours of $u_k$ lie in $C$. Otherwise the factor is $1$, whether $\tilde\delta_k$ is computed from all points or from $C$ alone.

So $\prod_k(1+g_k)=\prod_CW(C)$, the product running over components with at least $m+1$ points, where $W(C)=\prod_{k\in C}\max(1,\tilde\delta_k^{-\gamma})$ with $\tilde\delta_k$ computed within $C$.

**Step 3: the cluster bound.** For $n\ge m+1$ let

$$
Z_n=\int_{\{u_1=0,\ u\ \text{1-connected}\}}\ \prod_{k=1}^{n}\max\big(1,\tilde\delta_k^{-\gamma}\big)\,du_2\cdots du_n ,
$$

where "1-connected" means the graph joining points at distance at most $1$ is connected. Then $Z_n\le A^n\,n!$, with $A=e\max(1,2\pi m/\alpha)$.

1. *Split by spanning tree.* For almost every configuration the Euclidean minimum spanning tree $\tau$ is unique and all distances are distinct. By Cayley's formula there are $n^{n-2}\le e^nn!$ labelled trees. On $\{\text{EMST}=\tau\}$, 1-connectedness forces every edge of $\tau$ to have length at most $1$, because a minimum spanning tree also minimizes its longest edge.
2. *Bound the weights.* Single-linkage clustering is Kruskal's algorithm on $\tau$. Its dendrogram $T$ has the edges of $\tau$ as internal nodes, and $s_j$ denotes the length of edge $j$. For each point $k$, let $i(k)$ be the lowest ancestor of $k$ whose subtree has at least $m+1$ leaves.
   - Just before merger $i(k)$, the cluster of $k$ has at most $m$ points.
   - Every point outside that cluster is at distance at least $s_{i(k)}$ from $u_k$, because the tree path between two points has no edge longer than their distance.
   - Hence $\tilde\delta_k\ge s_{i(k)}$. This holds whether or not the cluster that $k$'s cluster merges with already has $m+1$ points.
   - Since all $s_j\le1$, the integrand is at most $\prod_js_j^{-\gamma c_j}$, where $c_j=\#\{k:i(k)=j\}$.
3. *Change variables.* Pass to the edge vectors of $\tau$; the map is linear with determinant $\pm1$. Then use polar coordinates. The orderings of $\tau$'s edges for which Kruskal's algorithm produces $T$ are exactly the linear extensions of $T$'s tree order. Dropping the spanning-tree constraint, the region for $T$ is $R(T)=\{s_j\le s_{{\rm parent}(j)},\ s\le1\}$.
4. *Integrate the nested scales.* Let $h_j$ be the number of internal nodes in the subtree of $j$, and let $C_j=\sum_{i\le j}c_i$ be the sum over that subtree. Integrating from the leaves up gives
   $$\int_{R(T)}\prod_js_j^{\,1-\gamma c_j}\,ds=\prod_j\frac1{E_j},\qquad E_j=2h_j-\gamma C_j .$$
   - If the subtree of $j$ has at most $m$ leaves, then $C_j=0$ and $E_j=2h_j$.
   - Otherwise every leaf below $j$ is already critical, so $C_j=h_j+1$ and $E_j=(2-\gamma)h_j-\gamma\ge\alpha>0$.
   - The ratio $2h/((2-\gamma)h-\gamma)$ decreases in $h$, and $h\ge m$ here, so $2h_j/E_j\le2m/\alpha$.
   - Hence the integral is at most $(2m/\alpha)^{n-1}\prod_j1/(2h_j)$. This is the only place where $\gamma<2m/(m+1)$ is used.
5. *Sum over dendrograms.* By the hook-length formula for rooted trees, $T$ has $(n-1)!/\prod_jh_j$ linear extensions. Each of the $(n-1)!$ orderings of $\tau$'s edges yields exactly one $T$. So $\sum_{T\sim\tau}\prod_j1/h_j=1$.
6. *Total.* Each $\tau$ contributes at most $(2\pi)^{n-1}(2m/\alpha)^{n-1}2^{-(n-1)}=(2\pi m/\alpha)^{n-1}$. Summing over $\tau$ gives $Z_n\le n^{n-2}(2\pi m/\alpha)^{n-1}\le A^nn!$.

**Step 4: summing over clusters.** Write $\prod_CW(C)=\prod_C\big(1+(W(C)-1)\big)$ and expand. Since every term is non-negative, the expansion may be enlarged to all families of disjoint 1-connected sets with at least $m+1$ points each, which gives a pointwise upper bound.

- A set of $n$ points integrates to at most $\lvert\Lambda\rvert\,Z_n$.
- There are at most $N^n/n!$ ways to choose it.
- Hence each set contributes at most $N A^n(\varepsilon^2/\lvert\Omega\rvert)^{n-1}$ relative to $\lvert\Lambda\rvert^n$.
- Summing over families gives
  $$\int_{\Lambda^N}\prod_k(1+g_k)\,du\le\lvert\Lambda\rvert^N\exp\big(N\zeta(\varepsilon)\big),\qquad\zeta(\varepsilon)=\sum_{n\ge m+1}A^n\Big(\frac{\varepsilon^2}{\lvert\Omega\rvert}\Big)^{n-1}.$$

Choose $\varepsilon$ small enough that $\zeta(\varepsilon)\le1$. Undoing the rescaling ($\ell^{2N}\lvert\Lambda\rvert^N=\lvert\Omega\rvert^N$) gives $\int_{\Omega^N}\prod_k(1+g_k)\,dx\le e^N\lvert\Omega\rvert^N$, and Step 1 completes the proof. $\square$

*Checks.* For random labelled trees with up to 8 vertices, all edge orderings were enumerated. The grouping by Kruskal dendrogram matches the hook-length count, and $\sum_{T\sim\tau}\prod_j1/h_j=1$ exactly. The nested-integral formula of Step 3.4 was checked symbolically on a 4-leaf dendrogram with $m=2$ and $\gamma=6/5$: both sides equal $25/24$. The argument has been reviewed only by AI-assisted checks so far, and a referee's reading of Step 3 would be valuable.

**Units.** Fix $K\ge10$. A set $S$ of charges is *$K$-tight* if $\operatorname{diam}S\le\operatorname{dist}(S,\text{rest})/K$, with distances capped at $D$; singletons are tight.

- *Laminarity.* Two tight sets are either nested or disjoint. Suppose instead that $S$ and $T$ are tight and cross, and pick $s\in S\setminus T$, $t\in T\setminus S$ and $u\in S\cap T$. Since $s\notin T$, $\lvert s-u\rvert\ge K\operatorname{diam}T$; since $s,u\in S$, $\lvert s-u\rvert\le\operatorname{diam}S$. So $\operatorname{diam}S\ge K\operatorname{diam}T$, and by symmetry $\operatorname{diam}T\ge K\operatorname{diam}S$. Hence $\operatorname{diam}S=0$, which contradicts $s\ne u$.
- *Units.* For $m\ge1$, the *units* $\mathcal U_m(x)$ are the maximal tight sets with at most $m$ charges. By laminarity they partition the configuration.
- *Notation.* For a unit $B$, write $R_B=\operatorname{dist}(B,\text{rest})$, $\rho_B=R_B/3$, and $d_k$ for the distance from $x_k$ to the nearer Dirichlet row.

**Lemma 6.10 (overlap identity) [Theorem].** Let $x_1,\dots,x_n$ be interior sites with charges $q_k$, and choose radii $\rho_k\ge0$. Let $B_k$ be the set of interior sites within distance $\rho_k$ of $x_k$. Let $\mu_k$ be the exit law from $B_k$ of simple random walk started at $x_k$ and killed on the Dirichlet rows, with $\mu_k=\delta_{x_k}$ if $B_k$ is empty. Let $G_B$ be the Green function of $B$ with Dirichlet conditions outside $B$. Then for $k\ne l$,

$$
\Delta_{kl}:=C_a(x_k,x_l)-C_a(\mu_k,\mu_l)=\mathbb 1[x_l\in B_k]\,G_{B_k}(x_k,x_l)+\sum_{y\in B_l}\mu_k(y)\,G_{B_l}(x_l,y)\ \ge0 .
$$

With $\nu=\sum_kq_k\mu_k$,

$$
-\tfrac12\sum_{k\ne l}q_k\cdot q_l\,C_a(x_k,x_l)=-\tfrac12\langle\nu,C_a\nu\rangle+\tfrac12\sum_k\lvert q_k\rvert^2C_a(\mu_k,\mu_k)-\sum_{k<l}q_k\cdot q_l\,\Delta_{kl}. \tag{6.4}
$$

*Proof.*

1. By the strong Markov property at the exit time from $B_k$, $C_a(x_k,y)-C_a(\mu_k,y)=G_{B_k}(x_k,y)$ for $y\in B_k$, and $0$ for $y\notin B_k$.
2. Apply this with $y=x_l$. Then apply it with the roles of $k$ and $l$ exchanged, inside the $\mu_k$-average, using the symmetry of $C_a$. This gives the formula for $\Delta_{kl}$.
3. Identity (6.4) follows by expanding $\langle\nu,C_a\nu\rangle$. $\square$

When the balls are disjoint, every $\Delta_{kl}$ vanishes and (6.4) reduces to Lemma 6.4. When balls overlap, the overlapping pairs keep their exact short-distance interaction.

**Lemma 6.11 (unit majorant) [Theorem].** Let $\mathcal P$ be any partition of the charges into $K$-tight sets of at most $m$ charges each. For $k,l$ in the same set $B$, write $s_{kl}=\min(\rho_B,\sqrt{d_kd_l})$ and $D_{kl}=-2\vec\alpha_{j_k}\cdot\vec\alpha_{j_l}$. Then

$$
W_a(u)\le e^{c(N+M)}\prod_{B\in\mathcal P}\Big[\Big(\frac{D}{\max(\rho_B,a)}\Big)^{\sum_{k\in B}\gamma_k}\prod_{k<l\in B}\max\Big(\frac{s_{kl}}{\max(r_{kl},a)},1\Big)^{\kappa D_{kl}}\Big],
$$

with $c$ depending only on $K$, $m$, $\beta$ and $\Omega$. The bound keeps both the attractive pairs ($D_{kl}>0$) and the repulsive pairs ($D_{kl}<0$) inside each set.

*Proof.*

1. Apply Lemma 6.10, giving every charge of $B$ the radius $\rho_B$ centred at its own position, and drop $-\frac12\langle\nu,C_a\nu\rangle\le0$.
2. *Different sets.* For $k\in B$ and $l\in B'\ne B$, tightness gives $r_{kl}\ge\max(R_B,R_{B'})>\rho_B+\rho_{B'}+a$. So $x_l\notin B_k$ and $\operatorname{supp}\mu_k\cap B_l=\emptyset$, and therefore $\Delta_{kl}=0$.
3. *Self-energies.* $C_a(\mu_k,\mu_k)\le\max_{y\notin B_k}C_a(x_k,y)\le\frac1{2\pi}\ln\frac D{\max(\rho_B,a)}+K$, by (G1).
4. *Same set: upper bound.* Here $r_{kl}\le\operatorname{diam}B\le R_B/K<\rho_B/3$. The first term of $\Delta_{kl}$ is at most $\frac1{2\pi}\ln\frac{\rho_B}{r_{kl}}+c$ by (G5). The second term is at most $c$, because $\operatorname{supp}\mu_k$ lies at distance at least $\rho_B/2$ from $x_l$. Also $\Delta_{kl}\le C_a(x_k,x_l)\le\frac1{4\pi}\ln\big(1+4d_kd_l/r_{kl}^2\big)+K$ by (G4). Together, $\Delta_{kl}\le\frac1{2\pi}\ln_+(s_{kl}/r_{kl})+c$.
5. *Same set: lower bound.* Suppose $r_{kl}\le s_{kl}/8$. Then $d_k\ge7r_{kl}$, and the ball of radius $t=\min(\rho_B,d_k)/2\ge s_{kl}/2.2$ about $x_k$ lies in $B_k$ and avoids the walls. So $\Delta_{kl}\ge G_{B_k}(x_k,x_l)\ge\frac1{2\pi}\ln(t/r_{kl})-c$ by (G5). If instead $r_{kl}>s_{kl}/8$, simply use $\Delta_{kl}\ge0$. In both cases $\Delta_{kl}\ge\frac1{2\pi}\ln_+(s_{kl}/r_{kl})-c$.
6. Use the upper bound for attractive pairs and the lower bound for repulsive pairs, and note $-q_k\cdot q_l/2\pi=\kappa D_{kl}$. A set has at most $m^2/2$ pairs, so the $O(1)$ errors total at most $e^{c(N+M)}$. Distances below $a$ (multiply occupied sites) are handled as in Lemma 6.5. $\square$

**Lemma 6.12 (tightness forces $m$-th neighbours) [Theorem].** Let $\mathcal P=\mathcal U_m(x)$ be the partition into units. Then for every $B\in\mathcal U_m$ and every $k\in B$,

$$
R_B\ \ge\ c_{K,m}\,\delta_k^{(m)},\qquad c_{K,m}=\tfrac12(3K+3)^{1-m}.
$$

*Proof.*

1. Order the other units $B_1,B_2,\dots$ by $r_i=\operatorname{dist}(B,B_i)$, so that $r_1=R_B$. Let $S_j=B\cup B_1\cup\dots\cup B_j$.
2. Suppose $S_j$ has at most $m$ charges. Then $S_j$ is not tight, by maximality of $B$. So some unit outside $S_j$ lies within $K\operatorname{diam}S_j$ of it, and hence $r_{j+1}\le(K+1)\operatorname{diam}S_j$.
3. Units are tight, so $\operatorname{diam}B_i\le r_i/K$. Hence $\operatorname{diam}S_j\le3r_j$, and $r_{j+1}\le3(K+1)r_j$.
4. Let $\ell$ be the first index with $\lvert S_\ell\rvert>m$; then $\ell\le m$. The $m$-th nearest other charge of $x_k$ lies in $S_\ell$, so $\delta_k^{(m)}\le r_\ell+\operatorname{diam}B_\ell+\operatorname{diam}B\le2r_\ell\le2(3K+3)^{m-1}R_B$. $\square$

**The combined majorant.** Lemmas 6.11 and 6.12, with $\mathcal P=\mathcal U_m(x)$, give

$$
W_a(u)\le C^{N+M}\prod_k\Big(\frac{D}{\max(\delta_k^{(m)},a)}\Big)^{\gamma_k}\prod_{B\in\mathcal U_m(x)}\ \prod_{k<l\in B}\max\Big(\frac{s_{kl}}{\max(r_{kl},a)},1\Big)^{\kappa D_{kl}}. \tag{6.5}
$$

This replaces per-charge smearing at the nearest-neighbour distance (Lemma 6.5) with three ingredients: per-charge smearing at the $m$-th-neighbour distance, the exact interaction inside each unit, and wall-aware pair scales.

**Proposition 6.13 [Proposition, sketch].** Let $0<\beta^2<\beta^2_{\rm UV}(\hat g)$ and choose $m>\kappa/(1-\kappa)$, so that $\gamma_k\le2\kappa<2m/(m+1)$. Then for every $N$ the continuum version of the majorant (6.5) is integrable on $\Omega^N$, and the lattice weights are uniformly integrable. Hence (S2) holds in the whole window, and every term of $\Xi_a$ converges to a finite continuum limit, including near the Dirichlet walls. External charges satisfying (6.2) are handled in the same way, with the units containing an external point.

*Sketch.* For fixed $N$, split $\Omega^N$ into the finitely many regions on which $\mathcal U_m(x)$ is constant. Within a region, there are three kinds of collapse.

1. **Inside a unit.** The intra-unit factors in (6.5) reproduce the true power counting, because both attraction and repulsion are kept. So these collapses are integrable for $\beta^2<\beta^2_{\rm UV}$ (Section 6.1). Near a wall the pair scales $s_{kl}$ only weaken the factors.
2. **A group of at most $m$ charges spanning several units.** There is no singular factor, because each member's $m$-th neighbour lies outside the group.
3. **A group of $P>m$ charges.** The $m$-th-neighbour factors give at most $r^{-P\gamma}$ against a measure $r^{2(P-1)}$. This is integrable when $P\gamma<2(P-1)$, which follows from $\gamma<2m/(m+1)$.

A full proof bounds each $m$-th-neighbour factor by the geometric mean of the distances to the $m$ nearest charges. It then integrates each unit's internal coordinates first, at scale $\rho_B$, and applies the power-counting criterion: near the diagonal, $\int\prod_{k<l}\lvert x_k-x_l\rvert^{-a_{kl}}$ with real $a_{kl}$ converges if and only if $\sum_{k<l\in S}a_{kl}<2(\lvert S\rvert-1)$ for every subset $S$ with $\lvert S\rvert\ge2$.

The criterion must be used with signed exponents, because the repulsive intra-unit factors cannot be dropped. For example, in $a_4^{(1)}$ a unit $\{\vec\alpha_1,\vec\alpha_1,\vec\alpha_2,\vec\alpha_2\}$ is allowed once $m\ge4$, which the condition $m>\kappa/(1-\kappa)$ forces for $\kappa\ge3/4$. Its attractive pairs alone give $\sum a_{kl}=8\kappa$ against $2(\lvert S\rvert-1)=6$, which fails for $3\pi\le\beta^2<16\pi/5$, inside the window; with the repulsive pairs included the sum is $0$. For signed exponents the criterion follows from the single-linkage sectors used in the proof of Theorem 6.9. On the sector of a dendrogram, $s_{{\rm lca}(k,l)}\le\lvert x_k-x_l\rvert\le(n-1)\,s_{{\rm lca}(k,l)}$, and integrating the nested scales from the leaves up gives exactly the condition $\sum_{k<l\in S}a_{kl}<2(\lvert S\rvert-1)$ at every node $S$; necessity follows by uniform collapse, as in Proposition 6.7. We have checked the bookkeeping in the cases above but have not written out the general case. Structures below the lattice scale contribute at most $a^{2(\lvert S\rvert-1)-\kappa D(S)}\to0$. Pointwise convergence is (G2), and Vitali's theorem then gives convergence of each term.

So Proposition 6.13 supplies, at sketch level, the boundary layer that Section 6.1 left open.

**Conjecture 6.14 (Theorem 6.9 with unit factors) [Conjecture].** Let $0<\beta^2<\beta^2_{\rm UV}(\hat g)$ and $m>\kappa/(1-\kappa)$. Then the continuum version of the right-hand side of (6.5) satisfies

$$
\int_{\Omega^N}\prod_k\Big(\frac{D}{\delta_k^{(m)}}\Big)^{\gamma_k}\prod_{B\in\mathcal U_m(x)}\ \prod_{k<l\in B}\max\Big(\frac{s_{kl}}{r_{kl}},1\Big)^{\kappa D_{kl}}dx\le C^N\,N!^{\kappa}.
$$

By (6.5), Conjecture 6.14 implies (S1) with $\theta=\kappa<1$. Together with Proposition 6.13, it therefore implies Proposition 6.1 in the window.

Evidence:

- Without the unit factors this is Theorem 6.9, with $\gamma=2\kappa$.
- The unit factors involve only groups of at most $m$ charges. They are scale-invariant, and their singularities are integrable for $\beta^2<\beta^2_{\rm UV}$.
- The proof of Theorem 6.9 should extend by integrating each unit's internal coordinates inside the cluster bound of Step 3. This has not been checked.
- The exponent $\theta=\kappa$ is the one dictated by the scaling of the bulk free energy, and heuristically it cannot be improved (Lemma 6.15 and the remark after it).

**Lemma 6.15 (stability form of (S1)) [Theorem].** Write $c_N(a)=\sum_{j\in J^N}\prod_kn_{j_k}\int W_a(u(x))\,dx$, so that $\Xi_a(z)=\sum_{N\ge0}z^Nc_N(a)/N!$ with $c_N(a)\ge0$ (Theorem B).

1. If $c_N(a)\le C^NN!^{\theta}$ for all $N$ and $a$, with $0\le\theta<1$ (no loss of generality), then $\log\Xi_a(z)\le\ln2+(1-\theta)\,2^{\theta/(1-\theta)}(Cz)^{1/(1-\theta)}$ for all $z>0$ and all $a$.
2. Conversely, if $\log\Xi_a(z)\le A(1+z^\rho)$ for all $z>0$ and all $a$, then $c_N(a)\le e^A(eA\rho)^{N/\rho}N!^{\,1-1/\rho}$ for $N\ge1$, which is (S1) with $\theta=1-1/\rho$.

So (S1) holds for some $\theta<1$ if and only if the renormalized free energy at real fugacity grows at most polynomially in $z$, uniformly in $a$. In other words, (S1) is a stability bound.

*Proof.*

1. Put $s=1-\theta$ and $y=(Cz)^{1/s}$, so that $\Xi_a(z)\le\sum_N(y^N/N!)^s$. For $s=1$ this is $e^y$. For $0<s<1$, write $(y^N/N!)^s=\big((2^{(1-s)/s}y)^N/N!\big)^s\cdot2^{-(1-s)N}$ and apply Hölder's inequality with exponents $1/s$ and $1/(1-s)$. This gives $\sum_N(y^N/N!)^s\le2^{1-s}\exp\big(s\,2^{(1-s)/s}y\big)$.
2. All terms are non-negative, so $c_N(a)r^N/N!\le\Xi_a(r)\le e^{A(1+r^\rho)}$ for every $r>0$. Choose $r^\rho=N/(A\rho)$ and use $N^N\ge N!$. $\square$

**Remark (the exponent $\theta=\kappa$ is sharp) [heuristic].** For a long root $x=2\kappa$, so the growth exponent corresponding to $\theta=\kappa$ is $\rho=1/(1-\kappa)=2/(2-x)$. This is exactly how the bulk free energy scales with the physical coupling. With the coupling profile of Section 6.2 one expects $\log\Xi(z)\simeq K(\beta)\,z^{2/(2-x)}\int_\Omega R_\Omega^{2x/(2-x)}$ for large $z$, with $K(\beta)>0$ fixed by the exact bulk energy. If $K(\beta)\ne0$, then $\Xi$ has order exactly $1/(1-\kappa)$ as an entire function of $z$, and by part 1 no bound (S1) with $\theta<\kappa$ can hold. So the exponent of Conjecture 6.14 is optimal, and Theorem E already attains it (without external charges, $\gamma/2=\kappa$ there). What the window requires is a finite constant. That constant must diverge as $\beta^2\uparrow\beta^2_{\rm UV}$, where $\Xi_a$ itself diverges (Proposition 6.7) and the bulk energy has the resonance described in Section 6.1.

Given Proposition 6.13, Estimate (S), and hence Proposition 6.1 in the window, is reduced to Conjecture 6.14.

**Status of H1–H5.** With Theorems A–D, Theorem E closes H1–H5 in finite volume for $0<\beta^2<2\pi$, for every affine algebra. For sine-Gordon this is the whole range below the first collapse threshold. The window $2\pi\le\beta^2<\beta^2_{\rm UV}$ is non-empty for every other algebra. There, every term of the expansion converges (Proposition 6.13, sketch). Summability is reduced to Conjecture 6.14, the version of Theorem 6.9 that carries the unit factors of (6.5). At and above $\beta^2_{\rm UV}$ further counterterms are needed (Proposition 6.7). The infinite-volume limit remains Conjecture 6.2.

---

## 7. The classical theory: a selection rule from $\Theta$

### 7.1 Reality of energy and momentum

**Proposition 7.1 [Established; reformulated].** For a Hirota solution, the energy and momentum densities are total derivatives [Olive–Turok–Underwood 1993]. They integrate to

$$
E=\sum_iM_{a_i}\cosh\theta_i,\qquad P=\sum_iM_{a_i}\sinh\theta_i,
$$

where $(a_i,\theta_i)$ are the species and (possibly complex) rapidities of the constituent solitons [Takács–Watts 1999, eq. 4.6]. The masses $M_a$ are real and positive. The same holds for all higher conserved charges [Freeman 1995].

On Hirota data, $\Theta$ acts as $\tau_j\mapsto\bar\tau_j$. That is, $(a,\theta,p)\mapsto(\bar a,\bar\theta,\bar p)$, where $\bar a$ is the conjugate species ($\omega^{ja}\mapsto\omega^{-ja}$). Consequently:

- **Selection rule.** If the asymptotic data $\{(a_i,\theta_i)\}$ are $\Theta$-invariant, meaning the rapidities are real or occur in conjugate pairs with conjugate species, then $E$, $P$ and all conserved charges are real. Solitons and breathers are of this kind.
- **Complex energies.** The complex-energy solutions of Khastgir and Sasaki have asymptotic data that are not $\Theta$-invariant. They are not in the physical sector, for the same reason that complex eigenvalues of a $\Theta$-symmetric operator are not $\Theta$-invariant states.

This closes H6, apart from the question of singularities. The restriction to "good" solutions, previously ad hoc, is the classical shadow of Theorem C.

### 7.2 Singular solutions

Zeros of $\tau_j$ in the real $(x,t)$ plane are isolated points: one complex condition, two real conditions. There the energy density has double poles with zero residue, so integrated charges are unaffected, and the Hirota data continue analytically through the singularity.

The genuinely problematic phenomenon is different. In theories with non-self-conjugate species, soliton–breather scattering produces complex time delays and can carry regular data into singular regions [Takács–Watts 1999]. This is the classical image of a quantum sector in which $\Theta$ is broken (Section 9.3). It is not a failure of the formalism.

---

## 8. Semiclassics: complex saddles, thimbles, and the counting rule

### 8.1 Thimbles

The lattice integral (5.1) is a finite-dimensional integral of a holomorphic function of the complexified fields. Picard–Lefschetz theory therefore applies to it directly [Pham 1983; Witten 2010]. The integral decomposes as

$$
\int_{\mathbb R^n}=\sum_\sigma n_\sigma\int_{\mathcal J_\sigma},
$$

a sum over steepest-descent thimbles $\mathcal J_\sigma$ attached to complex critical points, with integer intersection numbers $n_\sigma$. Complex solitons are critical points of the continuum action, and the semiclassical expansion of a given topological sector is the expansion around the contributing thimbles. This answers the "why complex saddles" half of H7.

$\Theta$ maps thimbles to thimbles and conjugates their contributions, so the total is real (consistent with Theorem B). In the continuum the soliton critical points form complexified translation orbits $\{\vec\phi(x;\xi):\xi\in\mathbb C\}$. Removing the values of $\xi$ at which some $\tau_j$ vanishes on the real line cuts these orbits into strips in $\operatorname{Im}\xi$, which are the topological sectors.

### 8.2 The counting rule

**Proposition 8.1 [Theorem].** Let $B$ be a complex symmetric matrix with spectrum off $(-\infty,0]$, and let $\mathcal J$ be the thimble of $\tfrac12q^{\mathsf T}Bq$. Then

$$
\int_{\mathcal J}e^{-\frac12q^{\mathsf T}Bq}\,d^nq=(2\pi)^{n/2}\prod_\lambda\lambda^{-m_\lambda/2},
$$

where $m_\lambda$ is the **algebraic** multiplicity, whether or not $B$ is diagonalizable.

*Proof.* Triangularize $B$ (Schur form), deform the contour accordingly, and integrate the triangular quadratic form iteratively. $\square$

Applied to the Euclidean fluctuation operator $-\partial_\tau^2+\mathcal A$ about a static saddle, this gives the one-loop energy

$$
\Delta M=\tfrac12\operatorname{Tr}\big(\sqrt{\mathcal A}-\sqrt{\mathcal A_0}\big)-\tfrac14\operatorname{Tr}\big[(\mathcal A-\mathcal A_0)\mathcal A_0^{-1/2}\big], \tag{8.1}
$$

which counts eigenvalues by algebraic multiplicity. The second term is the normal-ordering counterterm.

The traditional DHN mode sum expands fluctuations in orthonormal eigenfunctions of $\mathcal A$. That presupposes $\mathcal A$ is diagonalizable with non-self-orthogonal eigenvectors, which can fail for a complex symmetric $\mathcal A$.

**An exactly solvable example [Theorem].** On $L^2(S^1)$, the operator $H_g=-\partial_x^2+g\,e^{ix}$ (with $g\ne0$) is lower triangular in the Fourier basis. Its spectrum is $\{k^2:k\in\mathbb Z\}$, and each $n^2$ with $n\ge1$ is a $2\times2$ Jordan block.

The path integral for $\operatorname{Tr}e^{-TH_g}$ is independent of $g$: every order in $g$ contains a zero-mode integral $\int_0^{2\pi}e^{imx_0}dx_0=0$. So it equals $\sum_ke^{-Tk^2}$, which counts each block twice. Counting eigenvectors would give $1+\sum_{n\ge1}e^{-Tn^2}$, which is wrong.

### 8.3 Resolution of the semiclassical mass disagreements

In $c_n^{(1)}$, the fluctuation operator around soliton $a<n$ has a Jordan block at frequency $\nu_0=m_{n-a}$. Its eigenvector comes from a double zero of the transmission factor of the perturbing soliton $n$. The eigenvector is self-orthogonal, and its generalized eigenvector sits at the $(n-a)$ threshold. Counting this block once gives the MacKay–Watts result; counting it with its algebraic multiplicity gives Delius–Grisaru's. Proposition 8.7 (Section 8.7) makes this precise. The eigenvector $\psi_n(\lambda_0)$ is in $L^2$ and self-orthogonal. The generalized eigenvector $\partial_\lambda\psi_n(\lambda_0)$ tends to a nonzero constant in channel $n-a$, so it is a bounded threshold resonance rather than an $L^2$ vector. The one-loop trace nevertheless counts $\lambda_0$ twice.

Three independent checks favour the algebraic-multiplicity count:

- **Lattice.** Evaluating (8.1) on a lattice reproduces the algebraic-multiplicity result to about $10^{-4}$ for all ten cases with $n\le5$. The same method reproduces Hollowood's simply-laced $a_{N-1}^{(1)}$ masses for $N\le6$.
- **Exact S-matrix.** The exact S-matrix of $c_2^{(1)}\simeq b_2^{(1)}$ agrees with the algebraic-multiplicity result to $O(\beta^4)$.
- **Fate of the block.** That exact S-matrix resolves the Jordan block into two real bound states, bound below threshold by amounts proportional to $\beta^2$ and $\beta^4$.

The details are in our companion analysis. For $b_n^{(1)}$, the exact masses of [GMW 1996] agree with the multiplicity-counted semiclassics, as those authors observed.

**Conjecture 8.2 [Conjecture].** For every $\hat g$ and every soliton obtained by folding Hirota multisolitons, the fluctuation solutions have the Hirota form (H) of Section 8.5. Then the factorization (F) holds by Theorem G, so the algebraic multiplicity of each isolated eigenvalue is the net order (zeros minus poles) of $\prod_bX_{ab}$, and the same net order is the weight in the one-loop trace. Wherever Hypothesis (R) holds, the Jordan blocks have the sizes of the individual zeros (Theorem F). Counting accordingly, the one-loop masses of all solitons obtained by folding multisolitons ($b_n^{(1)}$, $c_n^{(1)}$, $g_2^{(1)}$, $f_4^{(1)}$) satisfy

$$
\frac{M_a^{\rm qu}}{m_a^{\vee,\rm qu}}\ \text{independent of }a\ \text{to }O(\beta^2),
$$

where $m^\vee$ are the particle masses of the Lie-dual real-coupling theory.

Evidence:

- The $c_n^{(1)}$ blocks were constructed explicitly as Jordan chains.
- The lattice evaluation of (8.1) agrees with the algebraic-multiplicity count.
- The exact $b_n^{(1)}$ and $c_2^{(1)}$ S-matrices agree with it.
- (H) holds for $a_n^{(1)}$ and $c_n^{(1)}$, so for these the counting part is a theorem (Propositions 8.6 and 8.7). For $a_n^{(1)}$ with $n\le6$ the net count was also checked against the numerical spectrum.
- For $c_n^{(1)}$ with $n\le4$, the net count reproduces the full heat trace, including weight 2 at the threshold blocks (Proposition 8.7).

The counting part is a theorem once (H) or (R) is known (Theorems F and G, Proposition 8.5). The weight relevant to (8.1) is the net order from (F). It agrees with the dimension of the $L^2$ root space where (R) holds, but not at the threshold blocks of $c_n^{(1)}$ (Proposition 8.7). The simpler statement, that every zero of $X_{ab}$ in the bound-state half-plane is an eigenvalue of multiplicity equal to its order, is false: in $a_4^{(1)}$ a zero of $X_{12}$ is cancelled by a pole of $X_{11}$ (Proposition 8.6).

This closes H8 for $c_n^{(1)}$ and $b_n^{(1)}$. Elsewhere it reduces H8 to checking the Hirota form (H) of the explicit fluctuation solutions, for the $d_n^{(1)}$ and $e_n^{(1)}$ solitons and their foldings, together with the mass-ratio statement.

### 8.4 Which saddles contribute

The one-loop mass (8.1) is independent of $\operatorname{Im}\xi$. This was checked numerically across topological sectors, including nearly singular ones. Classically the same holds by Proposition 7.1. So the intersection numbers $n_\sigma$ determine which topological sectors appear, but not the masses.

**Conjecture 8.3 [Conjecture].** The intersection numbers are nonzero precisely for the sectors needed to fill the $U_q(\hat g^\vee)$ multiplets, and this resolves the missing-topological-charge problem (H11) at the quantum level.

This is the least developed part of the foundation. Computing $n_\sigma$ even in a reduced model would be valuable.

### 8.5 The counting rule from the transmission factors

**Setting.** The fluctuation operator about a static soliton is $\mathcal A=-\partial_x^2+M(x)$ on $L^2(\mathbb R,\mathbb C^r)$, with $M(x)=m^2\sum_jn_j\vec\alpha_j\vec\alpha_j^{\mathsf T}e^{i\beta\vec\alpha_j\cdot\vec\phi_s(x)}$. The matrix $M$ is complex symmetric, and it tends exponentially to $M_0=m^2\sum_jn_j\vec\alpha_j\vec\alpha_j^{\mathsf T}$ at both ends, because both vacua lie on the co-weight lattice.

- **Channels.** Let $e_b$ be eigenvectors of $M_0$ with eigenvalues $m_b^2$, adapted to the species and normalized by $e_b^{\mathsf T}e_c=\delta_{c\bar b}$; so $e_b$ is isotropic when $b\ne\bar b$. Let $\kappa_b(\lambda)=\sqrt{m_b^2-\lambda}$ with $\operatorname{Re}\kappa_b\ge0$. At a point $\lambda_0$, channel $c$ is *closed* if $\operatorname{Re}\kappa_c(\lambda_0)>0$ (always, for non-real $\lambda_0$), *open* if $\lambda_0>m_c^2$, and *at threshold* if $\lambda_0=m_c^2$. Write $C$ for the set of closed channels.
- **Bilinear structure.** The form $\langle u,v\rangle=\int u^{\mathsf T}v$ satisfies $\langle\mathcal Au,v\rangle=\langle u,\mathcal Av\rangle$. The Wronskian $W(u,v)=u^{\mathsf T}v'-u'^{\mathsf T}v$ is constant for two solutions at the same $\lambda$, because $M^{\mathsf T}=M$.
- **Root space.** Let $\mathcal R(\lambda_0)=\bigcup_j\ker(\mathcal A-\lambda_0)^j$, taken in $L^2$. For an isolated eigenvalue ($\lambda_0\notin[m_1^2,\infty)$), $\mathcal R(\lambda_0)$ is the range of the Riesz projection, so $\dim\mathcal R(\lambda_0)$ is the algebraic multiplicity. For an eigenvalue embedded in the continuum there is no Riesz projection, and $\dim\mathcal R(\lambda_0)$ is the natural replacement.
- **Transmission factors.** In the convention of Section 8.3, the channel-$c$ fluctuation is $e^{ikx}e_c$ as $x\to-\infty$ and $X_c(k)e^{ikx}e_c$ as $x\to+\infty$, and bound states sit at zeros with $\operatorname{Im}k<0$. We write $X_c(\lambda):=X_c(-i\kappa_c(\lambda))$. Away from threshold $dk/d\lambda\ne0$, so the order of a zero is the same in $\lambda$ and in $k$.

**Hypothesis (R) at $\lambda_0$ (reflectionless, species-diagonal).** For each $c\in C$ and each $\lambda$ in a neighbourhood $U$ of $\lambda_0$, there are solutions $\psi_c(\cdot,\lambda)$ and $G_c(\cdot,\lambda)$ of $(\mathcal A-\lambda)u=0$, analytic in $\lambda$, and an analytic function $X_c$ on $U$, such that:

- **(R1)** $\psi_c=e^{\kappa_cx}(e_c+o(1))$ as $x\to-\infty$, and $\psi_c|_{\mathbb R_-}$ is analytic as an $L^2(\mathbb R_-)$-valued function;
- **(R2)** $(\psi_c-X_cG_c)|_{\mathbb R_+}$ is analytic as an $L^2(\mathbb R_+)$-valued function;
- **(R3)** $G_c=e^{\kappa_cx}(e_c+o(1))$ as $x\to+\infty$.

In words: $\psi_c$ is the channel-$c$ fluctuation solution, it is reflectionless, and it leaves through channel $c$ only, with nothing in open or threshold channels. This is the structure the semiclassical literature works with [Hollowood 1993; MacKay–Watts 1995]. It is not automatic. Proposition 8.6 shows that it fails where a zero of one transmission factor is cancelled by a pole of another. Proposition 8.7 shows that it fails at the $c_n^{(1)}$ blocks, where a forced tail of $\psi_c$ resonates with the threshold of another channel.

**Theorem F [Theorem].** Assume (R) at $\lambda_0$, and let $p_c=\operatorname{ord}_{\lambda_0}X_c$ for $c\in C$. Then the vectors

$$
\chi_{c,j}=\frac1{j!}\,\partial_\lambda^j\psi_c\Big|_{\lambda_0},\qquad0\le j<p_c,\quad c\in C,
$$

form a basis of $\mathcal R(\lambda_0)$ made of Jordan chains: $(\mathcal A-\lambda_0)\chi_{c,0}=0$ and $(\mathcal A-\lambda_0)\chi_{c,j}=\chi_{c,j-1}$. Consequently:

1. $\dim\mathcal R(\lambda_0)=\sum_{c\in C}p_c$. For isolated $\lambda_0$ this is the algebraic multiplicity.
2. The Jordan blocks have sizes $p_c$, one block for each channel with a zero. So $\mathcal A$ is diagonalizable at $\lambda_0$ if and only if all these zeros are simple.
3. Items 1 and 2 hold whether $\lambda_0$ is isolated or embedded in the continua of open channels, including when $\lambda_0$ is the threshold of another channel, provided (R) holds there. Proposition 8.7 gives a case where it does not.
4. The eigenvector heading a block of size at least 2 is self-orthogonal.

*Proof.*

**Step 1: the chains exist and lie in $L^2$.**

1. Differentiating $(\mathcal A-\lambda)\psi_c=0$ $j$ times at $\lambda_0$ gives the chain relations.
2. By (R1), each $\chi_{c,j}$ lies in $L^2(\mathbb R_-)$.
3. On $\mathbb R_+$, (R2) gives $\chi_{c,j}=\sum_{i\le j}\frac{X_c^{(i)}(\lambda_0)}{i!}\,\frac{\partial^{j-i}G_c}{(j-i)!}+L^2(\mathbb R_+)$. Every coefficient $X_c^{(i)}(\lambda_0)$ with $i\le j<p_c$ vanishes, so $\chi_{c,j}\in L^2(\mathbb R_+)$.
4. Each $\chi_{c,j}$ is smooth with $\chi$ and $\mathcal A\chi$ in $L^2$, so it lies in $H^2$.

**Step 2: independence.** The eigenvectors $\psi_c(\lambda_0)$ are independent, because their leading behaviours at $-\infty$ lie in the independent directions $e_c$. Jordan chains headed by independent eigenvectors are jointly independent.

**Step 3: the chains span $\mathcal R(\lambda_0)$.**

1. Consider the solutions of $(\mathcal A-\lambda_0)u=0$ that lie in $L^2(\mathbb R_-)$. Because the perturbation decays exponentially, asymptotic integration applies [Coddington–Levinson 1955, Ch. 3]. Open channels give oscillating solutions and threshold channels give $1$ and $x$, none of which are in $L^2$. So this space has dimension $\lvert C\rvert$, and the $\psi_c(\lambda_0)$ span it.
2. Let $u_0,\dots,u_{q-1}$ be any Jordan chain in $L^2$. We show by induction on $j$ that there is a polynomial vector $\alpha(\lambda)=(\alpha_c(\lambda))_{c\in C}$ such that $v(\lambda)=\sum_c\alpha_c(\lambda)\psi_c(\lambda)$ satisfies $u_i=\frac1{i!}\partial^iv(\lambda_0)$ for $i\le j$, and $(\alpha_cX_c)^{(i)}(\lambda_0)=0$ for all $i\le j$ and all $c$.
3. *Base case $j=0$.* By item 1, $u_0=\sum a_c\psi_c(\lambda_0)$. On $\mathbb R_+$, $u_0=\sum a_cX_c(\lambda_0)G_c(\lambda_0)+L^2$. The $G_c(\lambda_0)$ are independent modulo $L^2(\mathbb R_+)$ by (R3), so $a_cX_c(\lambda_0)=0$ for every $c$. Take $\alpha=a$.
4. *Step $j\to j+1$.* Since $(\mathcal A-\lambda_0)\big[\frac1{(j+1)!}\partial^{j+1}v(\lambda_0)\big]=u_j$, the difference between $u_{j+1}$ and that term is an $L^2(\mathbb R_-)$ solution, hence equal to $\sum b_c\psi_c(\lambda_0)$. Add $b_c(\lambda-\lambda_0)^{j+1}$ to $\alpha_c$; this does not change the lower derivatives. By (R2) and the induction hypothesis, $u_{j+1}=\sum_c\frac{(\alpha_cX_c)^{(j+1)}(\lambda_0)}{(j+1)!}G_c(\lambda_0)+L^2(\mathbb R_+)$. Since $u_{j+1}\in L^2$, each coefficient $(\alpha_cX_c)^{(j+1)}(\lambda_0)$ vanishes.
5. So $\alpha_cX_c=O((\lambda-\lambda_0)^q)$, which forces $\alpha_c=O((\lambda-\lambda_0)^{q-p_c})$. Then $u_j=\sum_c\sum_i\frac{\alpha_c^{(i)}(\lambda_0)}{i!}\chi_{c,j-i}$, and only terms with $j-i\le p_c-1$ survive. Hence every $u_j$ lies in $\operatorname{span}\{\chi_{c,m}:m<p_c\}$.

**Step 4: self-orthogonality.** $\langle\chi_0,\chi_0\rangle=\langle(\mathcal A-\lambda_0)\chi_1,\chi_0\rangle=\langle\chi_1,(\mathcal A-\lambda_0)\chi_0\rangle=0$. $\square$

**Remark.** Without the species-diagonal assumption, the same proof applies with $\psi_c-\sum_dT_{cd}G_d\in L^2(\mathbb R_+)$. It then gives $\dim\mathcal R(\lambda_0)=\operatorname{ord}_{\lambda_0}\det T$, with block sizes equal to the partial multiplicities of the matrix function $T$ at $\lambda_0$ [Gohberg–Sigal 1971]. Species-diagonality is what turns these into the individual orders. For isolated eigenvalues this recovers the Evans-function multiplicity theorem [Alexander–Gardner–Jones 1990] in this setting.

**Lemma 8.4 [Theorem].** Suppose $\psi_b(\cdot,k)$ is an exact single-channel solution at both $k$ and $-k$, as the Hirota fluctuation solutions are. Then $X_b(k)\,X_{\bar b}(-k)=1$. In particular $X_b(0)\ne0$, so no transmission factor vanishes at its own threshold.

*Proof.* The Wronskian $W(\psi_b(k),\psi_{\bar b}(-k))$ is constant. It equals $-2ik$ at $-\infty$ and $-2ik\,X_b(k)X_{\bar b}(-k)$ at $+\infty$. $\square$

**Proposition 8.5 (the net count).** Write $\mathcal A-\mathcal A_0=V_1V_2$ with exponentially decaying factors. Let $K(\lambda)=V_2(\mathcal A_0-\lambda)^{-1}V_1$ be the Birman–Schwinger operator, which is trace class, and let $d(\lambda)=\det(I+K(\lambda))$. Assume the factorization

$$
d(\lambda)=\prod_bX_b(\lambda) \tag{F}
$$

on the physical sheet $\mathbb C\setminus[m_1^2,\infty)$. Theorem G below proves (F) whenever the fluctuation solutions have the Hirota form (H).

1. **[Theorem]** For every isolated eigenvalue $\lambda_0$, $m_{\rm alg}(\lambda_0)=\sum_b\operatorname{ord}_{\lambda_0}X_b$, where zeros count positively and poles negatively.
2. **[Proposition, sketch; checked numerically for the $c_n^{(1)}$ solitons with $n\le4$ (Proposition 8.7)]** Suppose in addition that no $X_b(k)$ vanishes at real $k$. Then for $t>0$,
   $$
   \operatorname{Tr}\big(e^{-t\mathcal A}-e^{-t\mathcal A_0}\big)=\sum_b\Big[\sum_\zeta\operatorname{ord}_\zeta X_b\;e^{-t(m_b^2+\zeta^2)}+\frac1{\pi i}\int_0^\infty e^{-t(m_b^2+k^2)}\,\partial_k\log X_b(k)\,dk\Big],
   $$
   where the inner sum runs over the zeros and poles of $X_b$ with $\operatorname{Im}\zeta<0$, signed as in part 1. Each zero is counted with its order whether its energy is isolated or embedded in other channels' continua.

*Proof of part 1.*

1. The function $d$ is analytic off $[m_1^2,\infty)$. The Konno–Kuroda resolvent formula $R=R_0-R_0V_1(I+K)^{-1}V_2R_0$ and cyclicity of the trace give $d'/d=\operatorname{tr}(R_0-R)$.
2. Integrate over a small circle around $\lambda_0$. The term $R_0$ is analytic there, and $-\frac1{2\pi i}\oint R\,d\lambda$ is the Riesz projection $P$. So $\operatorname{ord}_{\lambda_0}d=\operatorname{tr}P=m_{\rm alg}(\lambda_0)$.
3. By (F), $\operatorname{ord}_{\lambda_0}d=\sum_b\operatorname{ord}_{\lambda_0}X_b$. $\square$

*Sketch of part 2.*

1. Write the heat trace as $\frac1{2\pi i}\int_\Gamma e^{-t\lambda}\,\partial_\lambda\log d\,d\lambda$, with $\Gamma$ a contour around the spectrum.
2. By (F) the integrand splits into one term per channel. The $b$-th factor $X_b$ depends on $\lambda$ only through $k_b$, so it is meromorphic off $[m_b^2,\infty)$. Deform each channel's contour to wrap $[m_b^2,\infty)$ and encircle its zeros and poles.
3. The jumps across the cuts, summed over $b$ using Lemma 8.4, give the integrals. Since $X_b(0)\ne0$, there is no threshold term.

**Hypothesis (H) (Hirota form).** For each channel $b$ there are solutions $\psi_b^\mp(\cdot,\lambda)$ of $(\mathcal A-\lambda)u=0$, meromorphic in $\lambda$ on the physical sheet, of the form

$$
\psi_b^\mp(x,\lambda)=e^{\pm\kappa_b x}R_b^\mp(x,\lambda).
$$

Here each $R_b^\mp$ is a rational function, with vector coefficients, of finitely many exponentials $E_i=e^{\mu_ix+\xi_i}$ with $\mu_i>0$. Its limits are as follows: $R_b^-\to u_b$ as $x\to-\infty$; $R_b^-\to X_b(\lambda)\,u_b$ as $x\to+\infty$; and $R_b^+\to u_b$ as $x\to+\infty$. The channel vectors are normalized so that $u_b^{\mathsf T}u_c=\delta_{c\bar b}$.

This is the form of the fluctuation solutions obtained by linearizing a Hirota two-soliton solution in the amplitude of the second soliton, continued to a linear wave. $\psi_b^+$ is the solution at $-k$, divided by $X_b(-k)$. Sections 8.6 and 8.7 show that (H) holds for $a_n^{(1)}$ and $c_n^{(1)}$.

**Theorem G [Theorem].** Assume (H). Then:

1. The Wronskian matrix is $W(\psi_b^-,\psi_d^+)=-2\kappa_bX_b(\lambda)\,\delta_{d\bar b}$; that is, $\mathcal W=X\mathcal W_0$.
2. The factorization (F) holds: $\det(I+K(\lambda))=\prod_bX_b(\lambda)$.

Consequently, part 1 of Proposition 8.5 holds unconditionally for such solitons. The algebraic multiplicity of every isolated fluctuation eigenvalue is the net order (zeros minus poles) of $\prod_bX_b$.

*Proof.*

1. The Wronskian $W(\psi_b^-,\psi_d^+)$ is independent of $x$. For large $x$ it equals $e^{(\kappa_b-\kappa_d)x}F(x)$. Here $F(x)=\sum_rf_re^{-rx}$ is the convergent expansion of a rational function of the $E_i$, running over a discrete set of exponents $r\ge0$, and $f_0=-(\kappa_b+\kappa_d)X_b\,u_b^{\mathsf T}u_d$.
2. Suppose $m_b\ne m_d$. For $\lambda$ off a discrete set, $\kappa_b-\kappa_d$ is non-real or avoids every $r$. Since exponentials with distinct exponents are linearly independent, constancy in $x$ forces every term to vanish, so $W=0$. Suppose instead $m_b=m_d$. Then only the $r=0$ term can survive, and $W=f_0=-2\kappa_bX_b\,u_b^{\mathsf T}u_d$. Both identities extend to all $\lambda$ by meromorphy. This proves item 1.
3. By the matrix Jost–Pais formula [Jost–Pais 1951; Gesztesy–Makarov 2003; Gesztesy–Latushkin–Makarov 2007], $\det(I+K)=\det\mathcal W^{\rm J}/\det\mathcal W_0$. Here $\mathcal W^{\rm J}$ is the Wronskian matrix of the Jost solutions $\psi_b^{{\rm J},\mp}=e^{\pm\kappa_bx}(u_b+o(1))$ at $\mp\infty$.
4. Take $\lambda$ off a null set, where the real parts $\operatorname{Re}\kappa_c$ are distinct. The difference $\psi_b^--\psi_b^{{\rm J},-}$ decays at $-\infty$ and is $o(e^{\operatorname{Re}\kappa_bx})$ there. So it is a combination of the $\psi_c^{{\rm J},-}$ with $\operatorname{Re}\kappa_c>\operatorname{Re}\kappa_b$.
5. Order the channels by $\operatorname{Re}\kappa$. Then the change of basis from Jost to Hirota solutions is unit triangular, and the same holds at $+\infty$. Hence $\det\mathcal W=\det\mathcal W^{\rm J}$, and item 2 follows from item 1. $\square$

*Numerical check.* For the $a_4^{(1)}$ soliton $a=1$ and the $c_2^{(1)}$ and $c_3^{(1)}$ solitons, the Hirota Wronskian matrix is diagonal to $10^{-8}$ at complex $\lambda$, with diagonal entries $X_b$. The ratio $\det(\mathcal A_N-\lambda)/\det(\mathcal A_{0,N}-\lambda)$ on a Fourier grid, corrected for the momentum cutoff by the factor $\exp\!\big(-\int\operatorname{tr}(M-M_0)\,dx/\pi k_{\max}\big)$, equals $\prod_bX_b(\lambda)$ to about $10^{-4}$. It converges to it as the grid is refined: for $a_4^{(1)}$ the ratio is $0.99999$ at $N=640$.


**Check.** For the sine-Gordon kink, $\mathcal A=-\partial^2+1-2\operatorname{sech}^2x$ and $X(k)=(k+i)/(k-i)$. Part 2 gives $1-\operatorname{erfc}\sqrt t=\operatorname{erf}\sqrt t$. This has the correct small-$t$ behaviour, $2\sqrt{t/\pi}$, and tends to $1$ as $t\to\infty$, as the zero mode requires.

**Consequence for (8.1).** Pass to (8.1) through the Mellin representation of $\sqrt\lambda$; the counterterm removes the small-$t$ logarithmic divergence. Each zero $\zeta$ of $X_b$ then contributes $\operatorname{ord}_\zeta X_b\cdot\tfrac12\omega_\zeta$, with $\omega_\zeta=\sqrt{m_b^2+\zeta^2}$, and each pole subtracts in the same way. Where (R) holds, this weight is exactly the channel's contribution to $\dim\mathcal R$ (Theorem F), whether the eigenvalue is isolated or embedded. At threshold resonances it can be larger (Proposition 8.7).

### 8.6 Example: $a_n^{(1)}$, and why zeros alone over-count

Let $h=n+1$, $m_b=2m\sin(\pi b/h)$, $A=\pi a/h$ and $B=\pi b/h$. Take the static soliton of species $a$, $\tau_j=1+\omega^{ja}E$ with $\omega=e^{2\pi i/h}$ and $E=e^{m_ax+\xi}$. Linearizing the two-soliton Hirota solution in the amplitude of the second soliton, with that soliton continued to a linear wave $e^{ikx-i\omega t}$ of species $b$, gives exact fluctuation solutions. Their transmission factor is the two-soliton interaction coefficient [Hollowood 1992; Fring–Johnson–Kneipp–Olive 1994] evaluated at $\cosh\theta=ik/m_b$:

$$
X_{ab}(k)=\frac{ik/m_b-\cos(A-B)}{ik/m_b-\cos(A+B)}.
$$

So $X_{ab}$ has a zero in the bound-state half-plane when $\cos(A-B)>0$, at $\kappa=m_b\cos(A-B)$ and $\lambda=m_b^2\sin^2(A-B)$. It has a pole there when $\cos(A+B)>0$.

Explicitly, the channel-$b$ solution has components $f_j=e^{ikx}\big(\omega^{jb}+\omega^{j(a+b)}X_{ab}E\big)/(1+\omega^{ja}E)$ along $\vec\alpha_j$. This is a rational function of $E$ times $e^{ikx}$, with limits $u_b\propto\sum_j\vec\alpha_j\omega^{jb}$ at $-\infty$ and $X_{ab}u_b$ at $+\infty$. So these solutions have the Hirota form (H), and Theorem G applies.

**Proposition 8.6.**

1. **[Theorem]** Take a zero of $X_{ab}$ with $b\ne a$ and $\cos(A-B)>0$. The exact solution at the zero is in $L^2$ if and only if $\cos B\,\sin(A-B)>0$. Otherwise it grows at $+\infty$ in channel $c\equiv b-a\pmod h$, and $X_{ac}$ has a pole at the same energy. (For $b=a$ the zero is the translation mode at $\lambda=0$.)
2. **[Theorem, by Theorem G and Proposition 8.5; checked numerically for $n\le6$]** For every $n$ and every soliton, the isolated eigenvalues of $\mathcal A$, counted with algebraic multiplicity, are exactly the net count (zeros minus poles) of $\prod_bX_{ab}$. Counting zeros alone over-counts. For example, in $a_4^{(1)}$ with $a=1$, $X_{12}$ has a simple zero at $\lambda=1.25\,m^2$, below the lowest threshold $m_1^2=1.382\,m^2$. But $X_{11}$ has a pole there, and $\mathcal A$ has no eigenvalue there: below threshold it has only the zero mode.

*Proof of part 1.*

1. At the zero, the channel-$b$ solution is $f_j=e^{\kappa x}\omega^{jb}/(1+\omega^{ja}E)$.
2. As $x\to+\infty$ its leading term is $e^{(\kappa-m_a)x}\omega^{j(b-a)}$, and $\kappa-m_a=-2m\cos B\sin(A-B)$. So the solution decays if and only if $\cos B\sin(A-B)>0$.
3. At the zero, $\kappa_c(\lambda)=m_c\lvert\cos B\rvert$ and $\cos(A+C)=\pm\cos B$, where $C=\pi c/h$. In the growing case these agree, which is exactly the pole condition for $X_{ac}$. $\square$

*Numerics.* The operator $\mathcal A$ was discretized on a periodic Fourier grid, with $\operatorname{Im}\xi$ chosen to keep every $\tau_j$ away from zero on the real line. Two grids were used, up to 768 points on periods 32–56, and all below-threshold eigenvalues were compared with the net count. The agreement covers both which eigenvalues occur and their multiplicities, for example $\{0,\,0.4603\,(\times2),\,0.7157\,(\times2)\}$ for $a_6^{(1)}$, $a=3$. All the multiplicities seen are consistent with semisimple eigenvalues, as Theorem F predicts for simple zeros in distinct channels.

### 8.7 Example: the $c_n^{(1)}$ Jordan blocks are threshold blocks

Realize $c_n^{(1)}$ as the $\sigma$-invariant sector of $a_{2n-1}^{(1)}$, where $\sigma:\vec\alpha_j\mapsto\vec\alpha_{-j}$. Then $h=2n$, $m_b=2m\sin(\pi b/2n)$, $A=\pi a/2n$ and $B=\pi b/2n$. The soliton $a<n$ is the coincident static pair $(a,\,2n-a)$, with

$$
\tau_j=1+2\cos(\pi ja/n)\,E+\cos^2\!A\;E^2,\qquad E=e^{m_ax+\xi}.
$$

Channel $b$ of $c_n^{(1)}$ is the $\sigma$-invariant combination of the $a_{2n-1}^{(1)}$ channels $b$ and $2n-b$. Its fluctuation passes through both solitons, so its transmission factor is the product of the two $a_{2n-1}^{(1)}$ factors of Section 8.6. The fluctuation solutions come from the linearized three-soliton solution. They are rational in $E$ times $e^{ikx}$, so they have the Hirota form (H), and Theorem G applies. With $z=\kappa_b/m_b$,

$$
X_b(z)=\frac{\big(z-\cos(A-B)\big)\big(z+\cos(A+B)\big)}{\big(z-\cos(A+B)\big)\big(z+\cos(A-B)\big)},\qquad X_n(z)=\Big(\frac{z-\sin A}{z+\sin A}\Big)^2 .
$$

**Proposition 8.7.**

1. **[Theorem]** $X_n$ has a double zero at $z=\sin A$, that is at $\lambda_0=m_n^2\cos^2\!A=m_{n-a}^2$, the threshold of channel $n-a$. So $\lambda_0$ is never an isolated eigenvalue. It is embedded in the continua of channels $1,\dots,n-a-1$ when $n-a\ge2$, and sits at the bottom of the continuum when $a=n-1$.
2. **[Theorem]** The eigenvector $\psi_n(\lambda_0)$ is in $L^2$ and self-orthogonal. The generalized eigenvector $\partial_\lambda\psi_n(\lambda_0)$ satisfies $(\mathcal A-\lambda_0)\partial_\lambda\psi_n=\psi_n$ and decays at $-\infty$. At $+\infty$, however, it tends to a nonzero constant vector in channel $n-a$. No homogeneous solution can correct this, so the $L^2$ root space contributed by channel $n$ at $\lambda_0$ is one-dimensional. The second member of the chain is a bounded threshold resonance of channel $n-a$. Hypothesis (R) fails at $\lambda_0$, and Theorem F does not apply.
3. **[Numerical; (F) itself holds by Theorem G]** For every soliton with $n\le4$, the heat trace $\operatorname{Tr}(e^{-t\mathcal A}-e^{-t\mathcal A_0})$ in the $\sigma$-invariant sector agrees with Proposition 8.5, part 2, to four decimals at $t=0.3,0.6,1,1.5$, with weight 2 at $\lambda_0$. With weight 1 instead, the disagreement is $e^{-t\lambda_0}$. The check includes the embedded blocks and, for $c_3^{(1)}$ with $a=1$, a zero of $X_2$ at $\lambda=\tfrac34$ that is cancelled by a pole of $X_1$.

*Proof of part 1.* At $B=\pi/2$, $\cos(A\mp B)=\pm\sin A$, so the two zeros of $X_n$ coincide. Also $m_n^2\cos^2\!A=4m^2\cos^2(\pi a/2n)=m_{n-a}^2$. $\square$

*Proof of part 2.*

1. In both $a_{2n-1}^{(1)}$ factors the zero at $z_0=\sin A$ is simple. At $z_0$ only the first term of the linearized three-soliton numerator survives, so $\psi_n=e^{\kappa_0x}\sum_j\vec\alpha_j(-1)^j/\tau_j$.
2. This decays like $e^{\kappa_0x}$ at $-\infty$. At $+\infty$ it decays like $e^{(\kappa_0-2m_a)x}$, and $\kappa_0=m_n\sin A=m_a$. Self-orthogonality follows from the identity $\langle\psi,\psi\rangle=\lim_{L\to\infty}W(\partial_\lambda\psi_n,\psi)(L)$, because the tail of $\partial_\lambda\psi_n$ is bounded while $\psi$ decays.
3. Differentiating at $z_0$ gives $\partial_\kappa\psi_n=x\,\psi_n+e^{\kappa_0x}\sum_j\vec\alpha_j(-1)^j\big(\omega^{ja}A_1'+\omega^{-ja}A_2'\big)E/\tau_j$, with $A_1'=A_2'=1/(2\sin A)\ne0$.
4. As $x\to+\infty$, $E/\tau_j\sim e^{-m_ax}/\cos^2\!A$. So the second term tends to a constant multiple of $\sum_j\vec\alpha_j(-1)^j\cos(\pi ja/n)$, which is the $\sigma$-invariant channel-$(n-a)$ vector.
5. Any $L^2(\mathbb R_-)$ solution at $\lambda_0$ is a combination of closed-channel solutions $\psi_c(\lambda_0)$. Their $+\infty$ exponents are $\kappa_c-jm_a$ with $j\in\{0,1,2\}$. An exponent of zero requires $\kappa_c(\lambda_0)=m_a$, which means $m_c^2=m_a^2+m_{n-a}^2=m_n^2$, so $c=n$. But the rate-zero coefficient of $\psi_n(\lambda_0)$ is proportional to $A_1(z_0)$ and $A_2(z_0)$, which vanish. So the constant tail cannot be removed. $\square$

*Numerics.* The soliton and every fluctuation solution were checked against the field and fluctuation equations, with residuals at the level of the finite-difference error. The block structure was checked for all $a<n\le4$. Numerically $\lvert\langle\psi,\psi\rangle\rvert\sim10^{-15}$ against $\lVert\psi\rVert^2\sim10$–$40$. The tail of $\partial_\lambda\psi_n$ is constant to six digits from $x=10$ to $x=40$, and lies along channel $n-a$ to six digits. Heat traces were computed on periodic Fourier grids of periods 36–72. Because the background is reflectionless, they showed no dependence on the box.

**What this means for the counting rule.** The weight with which the one-loop trace (8.1) counts an energy is the net order of the product of all transmission factors (Proposition 8.5). This weight can differ from the dimension of the $L^2$ root space in both directions. It is smaller when a zero is cancelled by a pole in another channel (Proposition 8.6). It is larger at the $c_n^{(1)}$ blocks, where the second unit of weight is carried by a bounded threshold resonance (Proposition 8.7). On a lattice or in a box the threshold resonance becomes normalizable, the block becomes a genuine finite-dimensional Jordan block, and Proposition 8.1 counts it twice. This is why the lattice evaluation of (8.1) and the exact S-matrices both agree with Delius–Grisaru. It also fits the exact $c_2^{(1)}$ result of Section 8.3: quantum corrections pull both members of the threshold block below threshold, as two real bound states.

---

## 9. Scattering without unitarity

### 9.1 What survives, and why

Every known imaginary-coupling soliton S-matrix satisfies:

- factorization and the Yang–Baxter equation;
- crossing;
- braiding (R-matrix) unitarity, $\check S_{ab}(\theta)\check S_{ba}(-\theta)=1$;
- the bootstrap.

These are not independent assumptions. Delius [1995] showed that all four follow from the defining properties of the universal R-matrix, for any theory whose conserved charges generate $U_q(\hat g)$ and have definite Lorentz spins:

- the intertwining property $R\Delta=\Delta^{\mathsf T}R$ gives braiding unitarity and factorization;
- $(\mathcal S\otimes1)R=R^{-1}$ gives crossing;
- $(\Delta\otimes1)R=R_{13}R_{12}$ gives the bootstrap.

Two further consequences of his analysis matter here.

- **Crossing constrains the gradation.** Crossing is possible only for gradations satisfying a constraint that introduces a coupling-dependent "quantum" dual Coxeter number. The particle poles, and hence the quantum mass ratios, are therefore not rigid.
- **Real masses need $\lvert q\rvert=1$.** Particle poles lie on the imaginary rapidity axis, as stable particles require, only when $q$ is a pure phase. In imaginary-coupling ATFT, $q$ is a pure phase for real $\beta$ [Bernard–LeClair 1991], so this condition holds.

Given the quantum affine symmetry, these statements are [Established]. What they do not give is unitarity in the quantum-field-theory sense, $S(\theta)S(\theta)^\dagger=1$ for real $\theta$. Delius states explicitly that braiding unitarity is unrelated to whether the field theory is unitary. The S-matrices do not in general satisfy parity [Takács–Watts 1999], and braiding unitarity plus real analyticity plus parity would be needed to give QFT unitarity. The question is what replaces it.

### 9.2 Krein-unitarity

**Proposition 9.1 [Proposition, sketch].** Suppose the theory is realized on a Krein space with fundamental symmetry $P$ (Theorem C, in the continuum limit), its Hamiltonian is $P$-self-adjoint, and asymptotic completeness holds in the Krein sense. Then the S-matrix is Krein-unitary:

$$
S^{-1}=PS^\dagger P .
$$

Consequently the eigenvalues of $S_{ab}(\theta)$ for real $\theta$, and of every multi-particle Bethe–Yang transfer matrix $T_A$, are either pure phases or occur in pairs $(s,1/\bar s)$.

*Sketch.* If $H$ is $P$-self-adjoint, then $e^{-iHt}$ is Krein-unitary, and so are the Møller operators and $S$. For a Krein-unitary operator, $\sigma(S^{-1})=\sigma(PS^\dagger P)=\overline{\sigma(S)}$. The transfer matrices are products of two-body S-matrices acting on tensor products of Krein spaces, so they are Krein-unitary as well. $\square$

**Check against known results.** Takács and Watts found that the $U_q(a_2^{(1)})$ R-matrix eigenvalues

$$
e_1=\frac{1+q\sqrt x}{q+\sqrt x},\qquad e_2=\frac{1-q\sqrt x}{q-\sqrt x}
$$

are not phases for $\lvert q\rvert=1$ and $x<0$. Writing $\sqrt x=iy$, a short computation gives $1/\bar e_2=e_1$ exactly. The non-phase eigenvalues are precisely a pair $(s,1/\bar s)$, as Proposition 9.1 requires.

**Proposition 9.2 (an exactly Krein-unitary integrable regularization) [Theorem].** Let $q=e^{i\mu}$ with $0<\mu<\pi$, and let $x>0$. Let $R(x)$ be Jimbo's $U_q(\widehat{sl}_n)$ R-matrix on $\mathbb C^n\otimes\mathbb C^n$ in the homogeneous gradation,

$$
R(x)=(xq-x^{-1}q^{-1})\sum_iE_{ii}\otimes E_{ii}+(x-x^{-1})\sum_{i\ne j}E_{ii}\otimes E_{jj}+(q-q^{-1})\sum_{i\ne j}x^{\operatorname{sgn}(j-i)}E_{ij}\otimes E_{ji},
$$

let $S$ be the swap of the two factors, and let $G(x)=S\,R(x)/(xq-x^{-1}q^{-1})$. On $2M$ sites with periodic conditions, the light-cone (brickwork) evolution is $U=U_eU_o$. Here $U_o$ is the product of $G$ over the bonds $(1,2),(3,4),\dots$ and $U_e$ the product over $(2,3),\dots,(2M,1)$. The light-cone cutoff rapidity $\Lambda$ enters through $x=e^{2\Lambda}$ [Destri–de Vega 1987]. Let $\Pi$ be the reflection $j\mapsto2M+1-j$ about a bond centre, $\Pi_s$ the reflection $j\mapsto2-j$ (mod $2M$) about a site, $C$ the colour reversal $i\mapsto n+1-i$ on every site, and $K$ complex conjugation in the standard basis. Then:

1. $G(x)^{\mathsf T}=G(x)$, $G(x)G(x^{-1})=1$, and $\overline{G(x)}=S\,G(x)^{-1}S$.
2. $U^\dagger\Pi U=\Pi$, so $U$ is unitary for the Krein form $[\psi,\chi]=\langle\psi,\Pi\chi\rangle$.
3. $\Pi_s\overline U\,\Pi_s=U^{-1}$, so the antilinear operator $\Theta_{\rm lat}=\Pi_sK$ satisfies $\Theta_{\rm lat}U\Theta_{\rm lat}^{-1}=U^{-1}$ and commutes with $H=\frac ia\log U$.
4. $C\Pi$ commutes with $U$.

Consequently the spectrum of $U$ is invariant under $\lambda\mapsto1/\bar\lambda$, and quasi-energies that are not real come in complex-conjugate pairs.

*Proof.*

1. The diagonal parts of $R(x)$ are symmetric and swap-invariant, and transposition and conjugation by $S$ both send $E_{ij}\otimes E_{ji}$ to $E_{ji}\otimes E_{ij}$. So $R^{\mathsf T}=SRS$, and $G^{\mathsf T}=R^{\mathsf T}S/(xq-x^{-1}q^{-1})=G$. The second identity is the unitarity relation $R_{12}(x)R_{21}(x^{-1})=(xq-x^{-1}q^{-1})(x^{-1}q-xq^{-1})$ [Established]. For the third, an entrywise check using $\bar q=q^{-1}$ and $x$ real gives $R(x)^\dagger=-R(x^{-1})$ and $\overline{xq-x^{-1}q^{-1}}=-(x^{-1}q-xq^{-1})$. Hence $\overline{G(x)}=G(x)^\dagger=R(x^{-1})S/(x^{-1}q-xq^{-1})=S\,G(x^{-1})\,S$, and $G(x^{-1})=G(x)^{-1}$.
2. $\Pi$ maps each layer to itself and conjugates the gate on a bond $b$ into $SGS$ on the mirror bond $b'$. By part 1, $G_b^\dagger=\Pi\,G_{b'}^{-1}\Pi$. Gates in one layer commute, so $U_o^\dagger=\Pi U_o^{-1}\Pi$ and $U_e^\dagger=\Pi U_e^{-1}\Pi$, and $U^\dagger=U_o^\dagger U_e^\dagger=\Pi U^{-1}\Pi$.
3. $\Pi_s$ exchanges the two layers. The same computation with $\overline G=SG^{-1}S$ gives $\Pi_s\overline{U_e}\Pi_s=U_o^{-1}$ and $\Pi_s\overline{U_o}\Pi_s=U_e^{-1}$, hence $\Pi_s\overline U\Pi_s=U_o^{-1}U_e^{-1}=U^{-1}$.
4. $C$ reverses $\operatorname{sgn}(j-i)$, so $(C\otimes C)R(C\otimes C)=SRS$ and $(C\otimes C)G(C\otimes C)=SGS$. Combined with step 2, $C\Pi$ maps each gate to the same gate on the mirror bond. $\square$

*Checks.* The identities of part 1 hold to machine precision for $n=2,3,4$, and parts 2–4 for $n=3$ on 6 sites. For $n=3$ neither $\Pi$ nor $C$ alone commutes with $U$. For $n=2$ every eigenvalue of $U$ lies on the unit circle, as for the unitary Destri–de Vega lattice. For $n=3$, on 6 and 8 sites, pairs $(\lambda,1/\bar\lambda)$ off the unit circle occur only in sectors that contain all three colours; sectors with two colours reduce to $n=2$. The pairs appear as soon as $\Lambda>0$, so at these sizes they are lattice-scale effects. Whether the low-lying ones survive the scaling limit is a question for the non-linear integral equations of open problem 8. Unlike Theorem B, the Euclidean weights ($\lvert x\rvert=1$) carry the phases $x^{\pm1}$. In this gauge, on 6 sites with $\mu=1$ and $x=e^{i(\pi-\mu)/2}$, the ratio of $\operatorname{Tr}V^m$ to the same trace with the phases removed is $0.78$, $0.66$, $0.49$ and $0.18$ for $m=1,2,4,8$.

Heuristically, $C$ acts on the weights of the vector representation as the longest Weyl element $w_0$, and $-w_0$ is the diagram automorphism. Then $C\Pi$ corresponds to spatial parity, and $\Pi$ to the field parity $P$ of Theorem C composed with the diagram automorphism and spatial parity, both symmetries of the theory. In this sense Proposition 9.2 is a lattice-exact, real-time counterpart of Theorem C, and its Krein-unitarity is the structure that Proposition 9.1 assumes in the continuum.

**Hermitian analyticity.** Miramontes [1999] showed that braiding unitarity becomes QFT unitarity exactly when, in addition, the S-matrix is Hermitian analytic. He also showed that Hermitian analyticity is consistent with the bootstrap but depends on the basis of states. Together with the results above, this gives a ladder of conditions, from weakest to strongest:

1. braiding unitarity, automatic for quantum-group S-matrices [Delius 1995];
2. Krein-unitarity (Proposition 9.1), which allows eigenvalue pairs $(s,1/\bar s)$;
3. unimodular eigenvalues of every multi-particle transfer matrix, which is what a real finite-volume spectrum requires (to leading order in $e^{-mL}$) [Takács–Watts 1999];
4. Hermitian analyticity in some basis of states, which together with 1 is QFT unitarity [Miramontes 1999].

**Unbroken and broken sectors.** We call a sector **$\Theta$-unbroken** when condition 3 holds in it, and **$\Theta$-broken** otherwise.

- In a $\Theta$-broken sector there are eigenvalues with $\lvert s\rvert\ne1$, occurring in pairs $(s,1/\bar s)$. Eigenvalues do not depend on the basis, so condition 4 fails in every basis. Non-unitarity there appears as complex-conjugate pairs of finite-volume levels, as in non-unitary perturbed minimal models.
- In a $\Theta$-unbroken sector, condition 4 may or may not hold. Takács and Watts stress that the required change of basis must act on one-particle states only, and that this can fail even when every two-particle S-matrix has phase eigenvalues. Even when condition 4 does hold, the positive inner product need not be one under which the local fields are Hermitian (Section 10).

This reclassifies H9. QFT unitarity is not a property of the unrestricted theories and cannot be imposed. Krein-unitarity is the property that does hold. Section 9.3 classifies the sectors in which more can be recovered, and Section 10 lists what cannot be recovered anywhere.

### 9.3 Sector-by-sector status

| Sector | Status | Evidence |
|---|---|---|
| real coupling, all $\hat g$ | Hermitian and unitary | [Established] |
| breather–breather, all $\hat g$ | $\Theta$-unbroken | the amplitudes are real-coupling S-matrices continued to $B<0$; their blocks are phases for real $\theta$ [Established] |
| all sectors of $a_1^{(1)}$ (sine-Gordon) | unbroken; real finite-volume spectra | [Takács–Watts 2002] |
| $a_2^{(2)}$, soliton–soliton | unbroken (numerical, three particles) | [Takács–Watts 1999] |
| $a_2^{(2)}$, soliton–excited soliton | broken; spectrum real only in the repulsive regime | [Takács–Watts 1999] |
| $a_n^{(1)}$ ($n\ge2$), $d_{2n+1}^{(1)}$, $e_6^{(1)}$: soliton–antisoliton and soliton–breather | broken for generic $\beta$; mirrored by complex classical time delays | [Takács–Watts 1999] for $a_2^{(1)}$; [Conjecture] in general |
| self-conjugate theories, soliton–breather | unbroken | [Conjecture]. Evidence: one-loop transmission factors are pure phases for real $k$ (shown for $c_n^{(1)}$); the transmission factors fail parity only when both species are non-self-conjugate [Takács–Watts 1999] |
| self-conjugate theories, soliton–excited soliton | to be checked case by case | $a_2^{(2)}$ shows these sectors can be broken in attractive regimes |

### 9.4 Unitary subtheories

At $q$ a root of unity, RSOS restriction replaces the Krein space by the space of states of a perturbed minimal model, which carries a positive inner product when the minimal model is unitary. The transfer-matrix eigenvalues after restriction are a *different* set from those before [Takács–Watts 1999]. So restriction can cure a broken sector, but it cannot touch quantum-group singlets such as breathers.

**Conjecture 9.3 [Conjecture].** An RSOS-restricted ATFT is unitary if and only if its restricted (IRF) S-matrix is Hermitian analytic [Miramontes 1999]. We conjecture that this holds exactly when (i) the underlying minimal CFT is unitary and (ii) every sector surviving the restriction is $\Theta$-unbroken.

Evidence:

- the de Vega–Fateev unitary $W_3$ series, at couplings where breathers decouple;
- $M_{6,7}+\Phi_{1,2}$, unitary;
- $M_{3,14}+\Phi_{1,5}$, where TCSA confirms a complex spectrum [Takács–Kausch–Watts 1997; Takács–Watts 1999];
- the $U_q(su(2)^{(1)})$-based S-matrices of Hoare, Hollowood and Miramontes [2013], whose IRF form at $q$ a root of unity is manifestly Hermitian analytic and hence unitary. This is the mechanism by which restriction restores unitarity. To our knowledge it has not been shown for the RSOS restrictions of higher-rank ATFTs.

This reduces H12 to a criterion that can be checked from the S-matrix.

---

## 10. What unitarity would give, and what cannot be rescued

In a unitary quantum field theory, unitarity is the source of a list of physical properties. This section sorts those properties, for imaginary-coupling ATFTs, into three groups:

- properties recovered in every theory;
- properties recoverable only in some sectors or after restriction;
- properties that **cannot be rescued in any way**.

**Scope.** The negative statements concern the *unrestricted* imaginary-coupling theories of rank $r\ge2$. Real-coupling ATFTs are unitary. Sine-Gordon theory, $a_1^{(1)}$, is unitary and has real finite-volume spectra in all sectors [Takács–Watts 2002]. Particular restricted models can be unitary (Section 9.4).

**What "cannot be rescued" means.** A property is unrescuable when its failure is detected by a quantity that no admissible change of description alters: an eigenvalue of the Hamiltonian or of a transfer matrix, the modulus of such an eigenvalue, the existence of a Jordan block, or the residue of a scalar amplitude. Such quantities do not depend on the inner product, on the basis of states (including rapidity-dependent bases, so regradation does not help either), or on how the Krein form is chosen. For amplitudes of quantum-group singlets such as breathers, they are also unchanged by RSOS restriction [Takács–Watts 1999].

### 10.1 A no-go lemma

**Lemma 10.1 [Theorem].** Let $H$ be a linear operator on a vector space equipped with some positive-definite inner product $\langle\cdot,\cdot\rangle$ with respect to which $H$ is symmetric. Then:

- every eigenvalue of $H$ is real;
- if the space is finite-dimensional, $H$ is diagonalizable.

Likewise, if $U$ is unitary with respect to some positive-definite inner product, every eigenvalue of $U$ has modulus 1. Consequently:

- a non-real energy rules out every positive inner product making $H$ self-adjoint;
- a Jordan block rules it out in any finite-dimensional regularization;
- a transfer-matrix eigenvalue with $\lvert s\rvert\ne1$ rules out every positive inner product making the S-matrix unitary.

*Proof.* If $H\psi=E\psi$, then $E\langle\psi,\psi\rangle=\langle\psi,H\psi\rangle=\langle H\psi,\psi\rangle=\bar E\langle\psi,\psi\rangle$, and $\langle\psi,\psi\rangle>0$, so $E$ is real. A symmetric operator on a finite-dimensional inner-product space is diagonalizable by the spectral theorem. If $U\psi=s\psi$, then $\lvert s\rvert^2\langle\psi,\psi\rangle=\langle U\psi,U\psi\rangle=\langle\psi,\psi\rangle$, so $\lvert s\rvert=1$. $\square$

The lemma covers non-local metric operators as well. No choice of metric $\eta>0$ in the sense of quasi-Hermiticity [Scholtz–Geyer–Hahne 1992; Mostafazadeh 2002] can make an operator with a non-real eigenvalue self-adjoint.

### 10.2 The properties that cannot be rescued

Each item below gives the property as it holds in a unitary theory, then how it fails in imaginary-coupling ATFT, then the reason it cannot be rescued and the scope.

**U1. Real energies and stationary states in finite volume.**

- *In a unitary theory:* all energy levels are real.
- *In imaginary-coupling ATFT:* complex-conjugate pairs of levels occur in $\Theta$-broken sectors. This is predicted at leading Bethe–Yang order for $a_2^{(1)}$ whenever solitons and breathers are both present, and for $a_2^{(2)}$ in the attractive regime [Takács–Watts 1999]. It is confirmed by TCSA in $M_{3,14}+\Phi_{1,5}$ [Takács–Kausch–Watts 1997]. The zero-mode reductions of $a_2$ and $c_2$ show complex pairs among excited levels at all couplings (Section 5.6).
- *Why it cannot be rescued:* energies are eigenvalues and do not depend on the basis.
- *Scope:* every theory with a $\Theta$-broken sector. The zero-mode results suggest that complex finite-volume levels may be generic for $r\ge2$, but the minisuperspace approximation does not prove this.

**U2. A global probability interpretation (Born rule, conserved positive total probability).**

- *In a unitary theory:* the state space carries a positive inner product preserved by time evolution.
- *In imaginary-coupling ATFT:* no positive inner product makes $H$ self-adjoint.
- *Why it cannot be rescued:* this follows from Lemma 10.1 and U1. The only conserved form is the indefinite Krein form (Theorem C, Proposition 9.1).
- *Scope:* every theory with a non-real level or a Jordan block anywhere in its finite-volume spectrum.

**U3. Probability conservation in scattering ($\lvert s\rvert=1$; transmission and reflection probabilities in $[0,1]$).**

- *In a unitary theory:* every eigenvalue of the S-matrix is a pure phase.
- *In imaginary-coupling ATFT:* soliton–antisoliton and soliton–breather amplitudes in $a_2^{(1)}$, and soliton–excited-soliton amplitudes in $a_2^{(2)}$, have eigenvalue pairs $(s,1/\bar s)$ with $\lvert s\rvert\ne1$ [Takács–Watts 1999; Proposition 9.1]. One member of each pair has $\lvert s\rvert>1$.
- *Why it cannot be rescued:* eigenvalue moduli do not depend on the basis, so Hermitian analyticity [Miramontes 1999] cannot hold in any basis. Breather amplitudes are quantum-group singlets and are untouched by RSOS restriction.
- *Scope:* $a_n^{(1)}$ ($n\ge2$), $d_{2n+1}^{(1)}$ and $e_6^{(1)}$ for generic $\beta$ (established for $a_2^{(1)}$, conjectured for the others); $a_2^{(2)}$ in the attractive regime.

**U4. Bounded real-time evolution.**

- *In a unitary theory:* the norm of a state is constant in time.
- *In imaginary-coupling ATFT:* in a $\Theta$-broken sector one member of each complex pair grows like $e^{\lvert\operatorname{Im}E\rvert t}$. At exceptional points, evolution grows polynomially.
- *Why it cannot be rescued:* the growth rates are $\operatorname{Im}E$, which is basis-independent.
- *Scope:* as U1. The Euclidean, statistical interpretation is unaffected. In Euclidean time a complex pair gives oscillating exponential decay, $e^{-\operatorname{Re}E\,\tau}\cos(\operatorname{Im}E\,\tau)$, which is a real and measurable feature of correlation functions.

**U5. Positive residues of bound-state poles (real on-shell three-point couplings, positive norm of bound states).**

- *In a unitary theory:* each forward-channel bound-state pole has a residue of definite sign, set by the square of a real coupling.
- *In imaginary-coupling ATFT:* continuing $\beta\to i\beta$ makes the cubic couplings of the real-coupling theory imaginary. Some physical-strip poles of the breather amplitudes then change sign, which Hollowood called the hallmark of non-unitarity [Hollowood 1991]. The scaling Lee–Yang model, itself an RSOS restriction of sine-Gordon, is the classic example [Cardy–Mussardo 1989].
- *Why it cannot be rescued:* the residue of a scalar (singlet) amplitude is basis-independent and unchanged by restriction.
- *Scope:* [Established] for $a_n^{(1)}$; [Conjecture] for every $\hat g$ whose real-coupling theory has cubic couplings, which is all except $a_1^{(1)}$. That exception is consistent with the unitarity of sine-Gordon theory.

**U6. Reflection positivity and positive spectral densities of local fields.**

- *In a unitary theory:* Osterwalder–Schrader positivity gives positive Källén–Lehmann spectral densities for local fields.
- *In imaginary-coupling ATFT:* only twisted reflection positivity holds (Section 5.4). The Hamiltonian is self-adjoint only for the indefinite Krein form, under which the field $\vec\phi$ is $P$-odd. Through the form-factor expansion, the wrong-sign residues of U5 give negative contributions to spectral densities.
- *Why it cannot be rescued:* a positive quasi-Hermitian metric may exist on a $\Theta$-unbroken sector, but it is in general non-local, and the local fields are not Hermitian with respect to it. Positivity of *local* spectral densities is not restored.
- *Scope:* the structural statement is a [Theorem] (Section 5.4). Negative spectral weight in specific channels follows from U5 at the level of [Conjecture], beyond the cases where U5 is established.

**U7. Monotonicity results that rely on reflection positivity.**

- *In a unitary theory:* Zamolodchikov's $c$-theorem holds, and the ground-state scaling function measures the Virasoro central charge, $c_{\rm eff}=c$.
- *In imaginary-coupling ATFT:* RSOS restrictions include non-unitary minimal models, for example $M_{3,14}$ and $M_{2,5}$ with $c<0$, where $c_{\rm eff}\ne c$ and the $c$-theorem's positivity input is absent.
- *Why it cannot be rescued:* these models are established non-unitary CFTs.
- *Scope:* the quantity that remains meaningful is $c_{\rm eff}$, from finite-size scaling or the thermodynamic Bethe ansatz.

**U8. Diagonalizability (a complete set of stationary states or normal modes).**

- *In a unitary theory:* the Hamiltonian is diagonalizable.
- *In imaginary-coupling ATFT:* it fails at exceptional points. Semiclassically these are the Jordan blocks of the soliton fluctuation operators (Section 8.3). In finite volume they are expected wherever a real pair of levels turns into a complex pair as $L$ or $\beta$ varies.
- *Why it cannot be rescued:* the existence of a Jordan block is basis-independent (Lemma 10.1).
- *Scope:* this failure is harmless for partition functions and one-loop masses, which count by algebraic multiplicity (Section 8.2).

### 10.3 What is recovered

| Property | Status | Scope |
|---|---|---|
| real vacuum energy, of minimal real part | Theorem D | every theory, finite lattice volume |
| positive partition function; real, non-negative vertex-operator correlation sums | Theorem B | every theory |
| conjugation-symmetric spectrum; conserved Krein form | Theorem C, Proposition 9.1 | every theory |
| real masses of all asymptotic particles | Proposition 7.1 (classical); one-loop (Section 8); $\lvert q\rvert=1$ in every known exact S-matrix [Delius 1995] | every theory |
| factorization, Yang–Baxter, crossing, braiding unitarity, bootstrap | universal R-matrix [Delius 1995] | wherever the quantum affine symmetry holds |
| real finite-volume spectra; unimodular S-matrix eigenvalues | Bethe–Yang [Takács–Watts 1999] | $\Theta$-unbroken sectors only: all of sine-Gordon; breather–breather sectors; $a_2^{(2)}$ repulsive regime; conjecturally soliton–breather sectors of self-conjugate theories |
| QFT unitarity | Hermitian analyticity [Miramontes 1999] | specific RSOS restrictions (de Vega–Fateev $W_3$ series, $M_{6,7}+\Phi_{1,2}$); mechanism shown for $U_q(su(2)^{(1)})$-based S-matrices [Hoare–Hollowood–Miramontes 2013] |

In short, imaginary-coupling ATFTs with $r\ge2$ have a protected vacuum, real particle masses and a consistent S-matrix. They do not have, and cannot be given, a global probability interpretation, probability-conserving scattering, bounded real-time evolution, positive local spectral densities, or positivity-based monotonicity theorems. All of these are lost exactly where $\Theta$ is broken, or, for U5 to U7, wherever imaginary cubic couplings or non-unitary restrictions arise.

---

## 11. Verdict

| Criterion | Status |
|---|---|
| (C1) definition | lattice: Theorems A and B. Continuum, finite volume: Theorem E for $\beta^2<2\pi$; for $2\pi\le\beta^2<\beta^2_{\rm UV}$ every term converges (Proposition 6.13, sketch), and Proposition 6.1 is conditional on Conjecture 6.14; further counterterms needed above (Proposition 6.7). Infinite volume: Conjecture 6.2 |
| (C2) Euclidean axioms | twisted reflection positivity, giving a Krein-space OS reconstruction (Section 5.4); other axioms as for Coulomb gases |
| (C3) Hamiltonian | $P$-self-adjoint transfer matrix, $\Theta$-symmetry, real vacuum energy: Theorems C and D |
| (C4) real masses | classical: Proposition 7.1. One-loop: $\Theta$ plus complex translation, with lattice evaluation. Exact: real in every known S-matrix |
| (C5) scattering | factorization, Yang–Baxter, crossing, braiding unitarity and bootstrap follow from the universal R-matrix [Delius 1995] where S-matrices exist; Krein-unitarity: Proposition 9.1 (sketch), exact on the integrable light-cone lattice (Proposition 9.2). QFT unitarity needs Hermitian analyticity [Miramontes 1999] and fails in $\Theta$-broken sectors. **Incomplete for four families (H10)** |
| (C6) predictivity | TBA and TCSA agree with restricted S-matrices, including complex spectra; semiclassics agrees with exact S-matrices once eigenvalues are counted by algebraic multiplicity, which is the net order of the transmission factors (Theorem F, Propositions 8.5 and 8.6) |

**Conclusion.** Imaginary-coupling ATFTs are valid **non-unitary** quantum field theories.

- Their Euclidean theory is a positive multi-component Coulomb gas.
- Their Hamiltonian theory lives on a Krein space with a $PT$-type symmetry, whose vacuum energy is provably real.
- Their asymptotic particles have real masses.
- Their S-matrices are Krein-unitary and fully consistent.

They are not unitary. In $\Theta$-unbroken sectors their finite-volume spectra are real and their S-matrix eigenvalues are phases. QFT unitarity additionally requires Hermitian analyticity, which RSOS restriction can supply in favourable cases. The properties listed in Section 10 cannot be rescued by any choice of inner product, basis or restriction: a global probability interpretation, probability-conserving scattering, bounded real-time evolution in broken sectors, positive bound-state residues, positive local spectral densities, and positivity-based monotonicity theorems. Where $\Theta$ is broken, the complex finite-volume levels are physical predictions, confirmed numerically, and not inconsistencies. In the sense that matters for statistical field theory, and for the description of non-unitary perturbed CFTs, these are as legitimate as the scaling Lee–Yang model.

---

## 12. Open problems

1. Conjecture 6.14, the bound of Theorem 6.9 with the intra-unit factors of (6.5) included. By Lemmas 6.10–6.12 and Proposition 6.13 it completes Proposition 6.1 in the window $2\pi\le\beta^2<\beta^2_{\rm UV}$, where every term of the expansion is already shown, at sketch level, to converge. By Lemma 6.15, the bound (S1) that it would give is equivalent to a stability bound $\sup_a\log\Xi_a(z)\le A(1+z^{1/(1-\kappa)})$ at real fugacity. Also open here: a full write-up of the power counting sketched in Proposition 6.13, and the infinite-volume limit and mass gap (Conjecture 6.2), using multi-species Coulomb-gas methods, including the claim that only vacuum-energy counterterms are needed in the bulk for simply-laced algebras. Theorem E settles Proposition 6.1 for $\beta^2<2\pi$. Theorem 6.9 proves the bound without the unit factors. Proposition 6.8 shows that per-charge (Onsager or moment) bounds cannot go further.
2. Exact soliton S-matrices for $c_n^{(1)}$ ($n\ge3$), $f_4^{(1)}$, $e_6^{(2)}$ and $a_{2n-1}^{(2)}$. These require R-matrices in representations that are not multiplicity-free. Our one-loop results constrain them: masses $M_a\propto\epsilon_a\sin(a\pi/H)$ with $H=2n-\beta^2/4\pi+O(\beta^4)$ for $c_n^{(1)}$, and the exceptional points predict pairs of nearby poles. Two routes avoid the missing R-matrices. The nested Bethe ansatz of the lattice models of open problem 8 gives finite-volume levels and multi-particle transfer-matrix eigenvalues directly, which is all the $\Theta$-breaking criterion of Section 9.2 uses. For $a_{2n-1}^{(2)}$, the lattice mass ratios of [Vernier–Jacobsen–Saleur 2016] (Section 3.3), converted to the normalization of Section 2.1 (their longest roots have squared length 4), give a direct check of the one-loop counting rule of Section 8.
3. Conjecture 8.2 in general. The counting rule is now Theorems F and G and Proposition 8.5, including for embedded eigenvalues and eigenvalues at thresholds. What remains is fourfold. First, check the Hirota form (H) for the $d_n^{(1)}$ and $e_n^{(1)}$ solitons and their foldings ($b_n^{(1)}$, $g_2^{(1)}$, $f_4^{(1)}$ and the twisted algebras). By Theorem G this gives the factorization (F) for them. (H) is proved for $a_n^{(1)}$ and $c_n^{(1)}$ (Propositions 8.6 and 8.7). Second, give a general theory of the trace weight at threshold resonances, which can exceed the dimension of the $L^2$ root space (Proposition 8.7). Third, treat spectral singularities (zeros of $X_b$ at real $k$), which can occur in $\Theta$-broken sectors. Fourth, prove the mass-ratio statement.
4. Computation of the thimble intersection numbers, and Conjecture 8.3 (the filling of quantum multiplets).
5. A proof of Krein asymptotic completeness, completing Proposition 9.1.
6. The sector classification for self-conjugate theories, in particular soliton–excited-soliton sectors in attractive regimes.
7. Whether $\Theta$-unbroken sectors of higher-rank RSOS restrictions admit Hermitian-analytic (IRF) bases, extending [Hoare–Hollowood–Miramontes 2013] beyond $U_q(su(2)^{(1)})$.
8. Integrable lattice regularizations: light-cone (staggered) $U_q$ vertex models [Destri–de Vega 1987; Reshetikhin–Saleur 1994], with Bethe-ansatz control of the scaling limit through non-linear integral equations [Destri–de Vega 1992, 1995; Zinn-Justin 1998 for untwisted simply-laced $\hat g$]. Doikou and Nepomechie [1999] identify the hole S-matrix of the critical $A_{N-1}^{(1)}$ chain with the $a_{N-1}^{(1)}$ soliton S-matrix. Matching their S-matrix parameter $e^{i\pi/(\nu-1)}$, $\nu=\pi/\mu$, with Section 2.3 shows that for $a_{N-1}^{(1)}$ the lattice anisotropy $q=e^{i\mu}$ corresponds to $\beta^2=4(\pi-\mu)$; for $N=2$ this is the Destri–de Vega relation $\beta_{\rm SG}^2=8(\pi-\mu)$. Theorem E's range is then $\mu>\pi/2$, the attractive regime. The window is $\pi/N<\mu\le\pi/2$, and the neutral thresholds of Section 6.1 sit at $q^{kN}=-1$. The light-cone evolution is exactly Krein-unitary (Proposition 9.2).

   Three questions precede any use as a definition. (i) Which algebra and which regime. [Vernier–Jacobsen–Saleur 2016] obtain imaginary $a_{2n-1}^{(2)}$ Toda from the $U_q(a_{2n-1}^{(2)})$ chain. This suggests that the lattice carries $U_q(\hat g)$ and that $\hat g^\vee$ appears only in the continuum S-matrix, with a renormalized $q$ (Section 3.3). The same R-matrix also has non-Toda continuum limits in other regimes. (ii) The ground-state structure, which [Zinn-Justin 1998; Doikou–Nepomechie 1999] assume and which is subtle for twisted algebras. Rigorous Bethe-ansatz tools exist so far only for six-vertex models [Duminil-Copin–Kozlowski–Krachun–Manolescu–Tikhonovskaia 2022], that is for $a_1^{(1)}$, where the window is empty. The Euclidean weights carry phases, so Theorems B and D have no evident counterpart (Proposition 9.2, checks). (iii) The comparison with Section 6. The lattice has periodic boundary conditions and constant coupling, while Theorem E's theory has the coupling profile of Section 6.2. The natural target is a bridge theorem: for $0<\beta^2<2\pi$, the vacuum energy given by the non-linear integral equations on a circle of length $L_t$ equals $-\lim_{L_x\to\infty}L_x^{-1}\log\Xi$, where far from the walls the coupling is constant. This needs a crossed-channel version of Section 5.4, with the spatial zero mode compactified as in Appendix A.

   This route does not bypass Conjecture 6.14. What it can supply is the bulk energy and its singularities across the window (Section 6.1), finite-volume form factors through the quantum inverse problem [Hegedűs 2017], and the data missing in problem 2.

---

## Appendix A. The zero-mode model and the sign of the fugacity

**Model.** The zero-mode model is $H=\tfrac12\lvert\vec p\rvert^2-g\sum_jn_je^{i\vec\alpha_j\cdot\vec x}$ on $\mathbb R^r/(2\pi\Lambda^\vee)$, in a plane-wave basis $\vec p=\sum_ic_i\vec\alpha_i$ with $\lvert c_i\rvert\le P$. For $a_2$ we use roots of squared length 4; for $c_2$, long roots of squared length 4. The reported levels are converged for $P$ from 8 to 20.

**Sign equivalence.** A shift $\vec x\mapsto\vec x+\vec s$ multiplies each term by $e^{i\vec\alpha_j\cdot\vec s}$. One can choose $\vec s$ with $\vec\alpha_j\cdot\vec s\equiv\pi\pmod{2\pi}$ for all $j$ if and only if $\sum_jn_j$ is even. When it is, $g$ and $-g$ give the same spectrum: this holds for $c_2$ ($\sum n_j=4$) but not for $a_2$ ($\sum n_j=3$). For $a_2$, only the physical sign gives a positive $Z(T)$ and a real ground state.

## Appendix B. Lattice evaluation of the one-loop trace

Formula (8.1) is evaluated with Dirichlet conditions on $[-L/2,L/2]$ and second-order finite differences. The counterterm uses the lattice propagator $\langle\phi_i\phi_i^{\mathsf T}\rangle=(2h)^{-1}[\mathcal A_0^{-1/2}]_{ii}$. The full complex spectrum is used, the lattice-broken zero mode is removed, and results are extrapolated in $h^2$. Non-simply-laced theories are realized as folding-invariant sectors of their simply-laced parents.

Calibration:

- sine-Gordon: $-2m/\pi$ reproduced to $10^{-4}$;
- Hollowood's $a_{N-1}^{(1)}$ masses: reproduced for $N\le6$;
- results independent of $L$ and of the topological sector.

---

## References

- H. Aratyn, C.P. Constantinidis, L.A. Ferreira, Nucl. Phys. B406 (1993) 727.
- J. Alexander, R. Gardner, C. Jones, J. Reine Angew. Math. 410 (1990) 167.
- S. Albeverio and R. Høegh-Krohn, J. Funct. Anal. 16 (1974) 39.
- A. Alexandru, G. Başar, P.F. Bedaque, N.C. Warrington, *Complex paths around the sign problem*, Rev. Mod. Phys. 94 (2022) 015006.
- Y. Ashida, Z. Gong, M. Ueda, *Non-Hermitian physics*, Adv. Phys. 69 (2020) 249, arXiv:2006.01837.
- G. Benfatto, G. Gallavotti, F. Nicolò, Commun. Math. Phys. 83 (1982) 387.
- C.M. Bender and S. Boettcher, Phys. Rev. Lett. 80 (1998) 5243; C.M. Bender, Rep. Prog. Phys. 70 (2007) 947.
- L. Bercini, M. Fabri, A. Homrich, P. Vieira, *S-matrix bootstrap: supersymmetry, $\mathbb Z_2$, and $\mathbb Z_4$ symmetry*, arXiv:1909.06453.
- D. Bernard and A. LeClair, Commun. Math. Phys. 142 (1991) 99.
- H.W. Braden, E. Corrigan, P.E. Dorey, R. Sasaki, Nucl. Phys. B338 (1990) 689; B356 (1991) 469.
- J.L. Cardy and G. Mussardo, Phys. Lett. B225 (1989) 275.
- E.A. Coddington and N. Levinson, *Theory of Ordinary Differential Equations*, McGraw–Hill (1955).
- P. Christe and G. Mussardo, Nucl. Phys. B330 (1990) 465.
- E. Corrigan, P.E. Dorey, R. Sasaki, Nucl. Phys. B408 (1993) 579.
- G.W. Delius, *Exact S-matrices with affine quantum group symmetry*, Nucl. Phys. B451 (1995) 445, hep-th/9503079.
- G.W. Delius and M.T. Grisaru, Nucl. Phys. B441 (1995) 259.
- G.W. Delius, M.D. Gould, Y.-Z. Zhang, q-alg/9508012.
- G.W. Delius, M.T. Grisaru, D. Zanon, Nucl. Phys. B382 (1992) 365.
- C. Destri and H.J. de Vega, Nucl. Phys. B290 (1987) 363.
- C. Destri and H.J. de Vega, Phys. Rev. Lett. 69 (1992) 2313; Nucl. Phys. B438 (1995) 413 (hep-th/9407117).
- H.J. de Vega and V.A. Fateev, Int. J. Mod. Phys. A6 (1991) 3221.
- A. Doikou and R.I. Nepomechie, *Soliton S matrices for the critical $A^{(1)}_{N-1}$ chain*, hep-th/9906069.
- P.E. Dorey, Nucl. Phys. B358 (1991) 654; B374 (1992) 741; Phys. Lett. B312 (1993) 291.
- P. Dorey, S. Faldella, S. Negro, R. Tateo, Phil. Trans. R. Soc. A371 (2013) 20120052.
- P. Dorey and D. Polvara, JHEP 02 (2022) 199.
- G.V. Dunne and M. Ünsal, *New nonperturbative methods in quantum field theory: from large-N orbifold equivalence to bions and resurgence*, arXiv:1511.05977.
- G.V. Dunne, lectures on resurgence (CERN, 2024), arXiv:2511.15528.
- H. Duminil-Copin, K.K. Kozlowski, D. Krachun, I. Manolescu, T. Tikhonovskaia, *On the six-vertex model's free energy*, Commun. Math. Phys. 395 (2022) 1383.
- C.J. Efthimiou, Nucl. Phys. B398 (1993) 697.
- T. Eguchi and S.-K. Yang, Phys. Lett. B224 (1989) 373.
- V.A. Fateev, Phys. Lett. B324 (1994) 45.
- F. Fabri and D. Polvara, arXiv:2402.12087.
- G. Felder and A. LeClair, Int. J. Mod. Phys. A7 (1992) 239.
- B. Feigin and E. Frenkel, in *Integrable Systems, Quantum Groups and Quantum Field Theories* (1993), hep-th/9310022.
- M.E. Fisher, Phys. Rev. Lett. 40 (1978) 1610.
- M.D. Freeman, Phys. Lett. B261 (1991) 57; Nucl. Phys. B433 (1995) 657.
- A. Fring, P.R. Johnson, M.A.C. Kneipp, D.I. Olive, Nucl. Phys. B430 (1994) 597.
- A. Fring, G. Mussardo, P. Simonetti, Nucl. Phys. B393 (1993) 413.
- J. Fröhlich, Commun. Math. Phys. 47 (1976) 233.
- G.M. Gandenberger, Nucl. Phys. B449 (1995) 375; Int. J. Mod. Phys. A13 (1998) 4553 (hep-th/9703158).
- G.M. Gandenberger and N.J. MacKay, Nucl. Phys. B457 (1995) 240; Phys. Lett. B390 (1997) 185.
- G.M. Gandenberger, N.J. MacKay, G.M.T. Watts (GMW), Nucl. Phys. B465 (1996) 329.
- L. Ge, Y.D. Chong, A.D. Stone, Phys. Rev. A85 (2012) 023802.
- F. Gesztesy and K.A. Makarov, Integral Equations Operator Theory 47 (2003) 457.
- F. Gesztesy, Y. Latushkin, K.A. Makarov, Arch. Ration. Mech. Anal. 186 (2007) 361.
- I.C. Gohberg and E.I. Sigal, Math. USSR Sb. 13 (1971) 603.
- V. Gorbenko, S. Rychkov, B. Zan, *Walking, weak first-order transitions, and complex CFTs*, JHEP 10 (2018) 108, arXiv:1807.11512; part II, SciPost Phys. 5 (2018) 050, arXiv:1808.04380.
- C. Guillarmou, A. Kupiainen, R. Rhodes, *Compactified imaginary Liouville theory*, Commun. Am. Math. Soc. 5 (2025).
- U. Harder, A.A. Iskandar, W.A. McGhee, Int. J. Mod. Phys. A10 (1995) 1879.
- Á. Hegedűs, *Lattice approach to finite volume form-factors of the massive Thirring (sine-Gordon) model*, JHEP 08 (2017) 059, arXiv:1705.00319.
- T.J. Hollowood, Nucl. Phys. B384 (1992) 523; Phys. Lett. B300 (1993) 73; Int. J. Mod. Phys. A8 (1993) 947.
- T.J. Hollowood and P. Mansfield, Phys. Lett. B226 (1989) 73.
- K. Ito and C. Locke, Nucl. Phys. B885 (2014) 600.
- K. Ito and H. Shu, arXiv:1805.08062 (massive ODE/IM for $A_r^{(1)}$-type modified affine Toda).
- K. Ito and H. Shu, *ODE/IM Correspondence and Quantum Periods*, Springer (2025).
- P.R. Johnson, *On soliton quantum S-matrices in simply laced affine Toda field theories*, Nucl. Phys. B (1997).
- J. Junnila, E. Saksman, C. Webb, Ann. Appl. Probab. 30 (2020) 2099.
- S.P. Khastgir and R. Sasaki, Prog. Theor. Phys. 95 (1996) 485.
- T.R. Klassen and E. Melzer, Nucl. Phys. B338 (1990) 485.
- M.A.C. Kneipp and D.I. Olive, Commun. Math. Phys. 177 (1996) 561.
- R. Jost and A. Pais, Phys. Rev. 82 (1951) 840.
- R. Konno and S.T. Kuroda, J. Fac. Sci. Univ. Tokyo Sect. I 13 (1966) 55.
- H. Lacoin, R. Rhodes, V. Vargas, Commun. Math. Phys. 337 (2015) 569.
- G.F. Lawler and V. Limic, *Random Walk: A Modern Introduction*, Cambridge University Press (2010).
- S.L. Lukyanov and A.B. Zamolodchikov, JHEP 07 (2010) 008.
- N.J. MacKay and W.A. McGhee, Int. J. Mod. Phys. A8 (1993) 2791.
- N.J. MacKay and G.M.T. Watts, Nucl. Phys. B441 (1995) 277.
- W.A. McGhee, Int. J. Mod. Phys. A9 (1994) 2645.
- A. Mostafazadeh, J. Math. Phys. 43 (2002) 205; Int. J. Geom. Meth. Mod. Phys. 7 (2010) 1191.
- M.R. Niedermaier, Nucl. Phys. B424 (1994) 184.
- D.I. Olive and N. Turok, Nucl. Phys. B215 (1983) 470; B257 (1985) 277.
- D.I. Olive, N. Turok, J.W.R. Underwood, Nucl. Phys. B401 (1993) 663; B409 (1993) 509.
- K. Osterwalder and R. Schrader, Commun. Math. Phys. 31 (1973) 83; 42 (1975) 281.
- M.F. Paulos, J. Penedones, J. Toledo, B.C. van Rees, P. Vieira, *The S-matrix bootstrap II: two dimensional amplitudes*, JHEP 11 (2017) 143, arXiv:1607.06110.
- F. Pham, Proc. Symp. Pure Math. 40 (1983) 319.
- D. Polvara, JHEP 04 (2023) 020.
- N.Yu. Reshetikhin and F.A. Smirnov, Commun. Math. Phys. 131 (1990) 157.
- N.Yu. Reshetikhin and H. Saleur, Nucl. Phys. B419 (1994) 507.
- F.A. Smirnov, Int. J. Mod. Phys. A6 (1991) 1407.
- G. Takács, Nucl. Phys. B489 (1997) 532; B501 (1997) 711; B502 (1997) 629.
- G. Takács, H.G. Kausch, G.M.T. Watts, Nucl. Phys. B489 (1997) 557.
- G. Takács and G.M.T. Watts, hep-th/9810006 (*Non-unitarity in quantum affine Toda theory and perturbed conformal field theory*); hep-th/0203073 (*RSOS revisited*).
- B. Hoare, T.J. Hollowood, J.L. Miramontes, *Restoring unitarity in the q-deformed world-sheet S-matrix*, JHEP 10 (2013) 050, arXiv:1303.1447.
- T.J. Hollowood, *Quantum solitons in affine Toda field theories*, hep-th/9110010.
- J.L. Miramontes, *Hermitian analyticity versus real analyticity in two-dimensional factorised S-matrix theories*, Phys. Lett. B455 (1999) 231, hep-th/9901145.
- F.G. Scholtz, H.B. Geyer, F.J.W. Hahne, Ann. Phys. 213 (1992) 74.
- E. Vernier, J.L. Jacobsen, H. Saleur, *The continuum limit of $a^{(2)}_{N-1}$ spin chains*, arXiv:1601.01559.
- E. Witten, arXiv:1001.2933.
- Yao, *Compactified imaginary Toda theory*, arXiv:2605.25494 (2026).
- A.B. Zamolodchikov, Adv. Stud. Pure Math. 19 (1989) 641.
- A.B. Zamolodchikov and Al.B. Zamolodchikov, Ann. Phys. 120 (1979) 253.
- Al.B. Zamolodchikov, Nucl. Phys. B342 (1990) 695.
- Al.B. Zamolodchikov, *Mass scale in the sine-Gordon model and its reductions*, Int. J. Mod. Phys. A10 (1995) 1125.
- Z. Zhu and D.G. Caldi, *Reality of complex affine Toda solitons*, J. Math. Phys. 36 (1995).
- P. Zinn-Justin, *Non-linear integral equations for complex affine Toda models associated to simply laced Lie algebras*, J. Phys. A31 (1998) 6747, hep-th/9712222.
