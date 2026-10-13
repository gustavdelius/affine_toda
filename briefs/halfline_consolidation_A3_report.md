<!-- Consolidation agent report A3, 2026-10-13. Model: Claude Opus 5.5, effort max. Prompt: briefs/halfline_consolidation_common.md + briefs/halfline_consolidation_A3_prompt.md. Interrupted several times by usage limits and resumed. Saved verbatim by the coordinator (extracted from the agent transcript). Scripts, logs copied to briefs/halfline_consolidation_A3_scripts/ (explore/, logs/); the downloaded papers are not copied. -->

# A3 report: the two classical claims of problem 11

Both claims are **VERIFIED**, and claim 5 holds in a stronger form than the worker stated.
- **Claim 4 (coweights outside $2\Lambda_W^\vee$).** These are integrable with $K_+=e^{i\pi\vec\nu\cdot\vec H}$. The book's conclusion (ch. 30 image-form paragraph, chapter notes, Appendix D) is wrong, and the reason it gives ("not a vacuum") does not hold.
- **Claim 5 (Dirichlet on all components).** It is not integrable for any untwisted $\hat g\ne a_1^{(1)}$. This holds for every linear σ, in every representation, at real coupling and for complex fields at imaginary coupling. A separate test, independent of Sklyanin's formalism (the spin-±2 charges), agrees.
- **Gemini's argument for claim 5** is not valid as written, although its conclusion is right.

All checks use the book's own code (`boundary_classical`, `boundary_images`, `boundary_common`) or code written from scratch; none of the worker's code was used. Papers were read at the source.

---

## 1. Claim 4: half coweights in the soliton-preserving condition. Verdict: VERIFIED

**By hand.** Take σ = Hermitian conjugation at imaginary coupling. The book's Lax pair (eq-boundary-lax-at) gives:
- **Cartan part.** For general $K_+$, the Cartan part of $K_+a_t(\lambda)+a_t(1/\lambda)^\dagger K_+$ is
  $\frac{i\beta}2[K_+,\partial_x\operatorname{Re}\vec\phi\cdot\vec H]-\frac\beta2\{K_+,\partial_x\operatorname{Im}\vec\phi\cdot\vec H\}$.
- **Consequences of Dirichlet on Re φ.** $\partial_x\operatorname{Re}\vec\phi$ is free there, so $K_+$ commutes with $\vec H$. In the defining representation of $a_n^{(1)}$ it is therefore diagonal, and Neumann on $\operatorname{Im}\vec\phi$ follows.
- **Step terms.** With $K_+E_{\pm j}K_+^{-1}=\kappa_j^{\pm1}E_{\pm j}$ and $w_j=e^{i\beta\vec\alpha_j\cdot\vec\phi/2}$, they require $\kappa_j=\kappa_j^{-1}=\bar w_j/w_j=e^{-i\beta\vec\alpha_j\cdot\operatorname{Re}\vec\phi}$.
- **Result.** $\kappa_j=\pm1$, $\operatorname{Re}\vec\phi(0)=\frac\pi\beta\vec\nu$ with $\vec\nu$ any coweight, and $K_+\propto e^{i\pi\vec\nu\cdot\vec H}$. The worker's "$K_{i+1}/K_i=e^{\pm i\beta\alpha_j\cdot\operatorname{Re}\phi}$" is right if read as "both signs".
- **The "if" direction holds for every algebra.** $e^{i\pi\nu\cdot H}$ gives $\kappa_j=(-1)^{\alpha_j\cdot\nu}=\bar w_j/w_j$.
- **Reflection equation.** Every constant diagonal $K_+$ solves eq-classical-reflection at real spectral parameters. $D\otimes D$ commutes with every term of $r$ and $r^{t_1t_2}$, and the $e_{ik}\otimes e_{ik}$ terms of $r^{t_2}$ and $r^{t_1}$ pick up the same factor $d_id_k$. So that equation does not select the signs; the compatibility condition does.

**What the book claims, taken apart:**
- **prp-soliton-preserving** ($K_+=1$ iff $\operatorname{Re}\vec\phi(0)\in\frac{2\pi}\beta\Lambda_W^\vee$ and Neumann on Im) is **correct as stated**.
- **"The condition has no free parameter, apart from the discrete choice of the vacuum value"** (ch. 30, line 240) is **incomplete**. With $K_+\ne1$ there are $2^n$ classes in $\Lambda_W^\vee/2\Lambda_W^\vee$ for $a_n^{(1)}$, not one.
- **Image-form paragraph** (line 248): "a coweight not in $2\Lambda_W^\vee$ gives $\operatorname{Re}\vec\phi(0)=-\frac\pi\beta\vec\lambda$, which is not a vacuum, and the condition fails". This is true only for $K_+=1$. "Not a vacuum" is not what makes it fail: a Dirichlet value need not be a vacuum.
- **Notes** (line 328, "In the image form (1.4) the coweight must lie in $2\Lambda_W^\vee$") and the **Appendix D item** are wrong as general statements; see below.

Three questions the brief asks to keep separate:
- **Classically integrable?** Yes, for every coweight class.
- **Is the boundary value a vacuum?** Only for $\vec\nu\in2\Lambda_W^\vee$. For $h$ even and all $\vec\alpha_j\cdot\vec\nu$ odd, it is a stationary point of $V$ (a "maximum"); otherwise it is not even stationary.
- **Is the field next to the wall a constant vacuum?** Only in the trivial class. In the other classes no constant field of finite energy satisfies the condition. The static solutions are halves of static solitons:
  - A species-$b$ soliton at rest with $\operatorname{Re}\xi=0$ satisfies $\vec\phi(x)+\overline{\vec\phi(-x)}=2\operatorname{Re}\vec\phi(0)$, where $\frac\beta\pi\operatorname{Re}\vec\phi(0)$ is its topological charge.
  - $\operatorname{Re}E_{x<0}=M_b/2$.
  - Each class contains a weight of representation $b$ and one of $h-b$, so there are two degenerate half-solitons. Their vacua at $-\infty$ differ by one soliton charge.

**Sine-Gordon analogue: yes.** For $n=1$ and $\nu=\omega^\vee$, $\operatorname{Re}\phi(0)=\pi/\beta_{\rm SG}$, the maximum of the potential, with $K_+\propto\mathrm{diag}(1,-1)$. The ground states are the half-kink and half-antikink, which is Ghoshal–Zamolodchikov Dirichlet at the maximum.

**Numerical evidence** (book code: `boundary_classical.a_t/algebra/RE_plus/r_toda`, `boundary_images.Sol/interaction/roots_an`; written from scratch: $a_x$ from eq-boundary-lax-at, a 4th-order Magnus transport, and τ):
- **Compatibility.** Residual $\le1.7\times10^{-16}$ with $D=e^{i\pi\nu\cdot H}$ for all $2^n$ classes, $n=1,2,3$. $K_+=1$ gives $\ge0.34$ off $2\Lambda$.
- **Uniqueness.** The general $h\times h$ $K$ is unique up to scale and equals $D$ (to $<10^{-12}$). There is no $K$ if $\operatorname{Re}\phi(0)$ is moved off $(\pi/\beta)\Lambda$, or if Neumann on Im is dropped.
- **Dynamics on exact solutions.** Three-soliton image solutions: a moving soliton $(a,\theta,\xi_1)$, its mirror with $\xi_2=-\bar\xi_1-\ln A_{aa}(2\theta)-\ln A_{ab}(\theta)$, and a half-soliton $(b,0,\xi_3)$ with $\operatorname{Re}\xi_3=-\ln A_{ab}(\theta)$.
  - They satisfy the new condition at all times (residual $3\times10^{-14}$).
  - **Main evidence:** $\tau(\lambda)=\mathrm{tr}(K_+T(\lambda)T(1/\lambda)^\dagger)$ is conserved to $2.8\times10^{-8}$ with $K_+=D$, and varies by $\ge3.8\times10^{-2}$ with $K_+=1$.
  - Control: the book's ΘP pair (eq-thetap-mirror) with $K_+=1$ conserves τ to $3.8\times10^{-10}$.
  - Cases: $a_1$ (including the real kink reflected at the Dirichlet-at-maximum wall), $a_2$, $a_3$.
  - Some phase choices are singular: one of my first configurations came within $4\times10^{-4}$ of a zero of $\tau_j$ at $x=0$, which made branch jumps. Regular phases were found by a scan, and only those are used.

**Appendix D (Delius 1998b, item (1.4)): right correction, wrong reason, overstated conclusion.** I read hep-th/9809140 at the source:
- **What the paper says.**
  - (1.4) is the image condition $\phi^*(-x)=\phi(x)+2\pi i\lambda$ with λ "in the coweight lattice". Its justification is "adding 2πi times a coweight to a soliton solution gives an equivalent soliton solution".
  - The paper then states "The boundary condition corresponding to (1.4) is (1.5)".
  - §5, eqs. (5.1)–(5.2): "If one imposes our new boundary condition (1.5) then these equations are solved by $K_+=1$".
- **What is wrong.**
  - (1.4) gives $\operatorname{Im}\phi_D(0)=-\pi\lambda$, so (1.5) follows only for $\lambda\in2\Lambda_W^\vee$.
  - Delius's own justification produces only $2\Lambda_W^\vee$ ($\lambda\to\lambda-2\mu$).
  - For the other coweights, his own (5.1)–(5.2) are solved by $K_+=e^{i\pi\lambda\cdot H}$ (residual $4\times10^{-15}$; check 11), not by 1.
  - So the real error is the claimed correspondence of (1.4) and (1.5); (1.4) itself defines integrable boundaries. Corrected item in §4(G).
- The original brief's "known trap" (`halfline_open_11_soliton_preserving_scalar.md`, line 42: "The Dirichlet value must be a vacuum") is likewise wrong.

**Still open (claim 4):**
- The quantum reflection matrices for the new classes when $n\ge2$. For $n=1$ it is GZ Dirichlet at $\phi_0=\pi/\beta_{\rm SG}$.
- The involution of the charges (σ antilinear; open as for $K_+=1$).
- Whether the half-solitons are ground states for $n\ge2$.

---

## 2. Claim 5: Dirichlet on all components. Verdict: VERIFIED and strengthened

**Theorem (proof by hand).** Let $\hat g=g^{(1)}\ne a_1^{(1)}$ and impose $\vec\phi(0,t)=\vec\phi_0$ on all components, at real coupling or for complex fields at imaginary coupling. Then no $K_+(\lambda)\ne0$ that is independent of the fields satisfies eq-boundary-compatibility at generic λ. This holds for every σ that is the antipode or the transposition composed with a conjugation, and in every finite-dimensional representation without trivial summands. At imaginary coupling, antilinear σ fails too.

*Proof.*
1. Absorb $\mathrm{Ad}_A$ into $K_+$. Using $\rho(X)^T=-\rho^*(X)$, both cases become $K\rho(a_t(\lambda))=\rho'(a_t(1/\lambda))K$, with $\rho'=\rho$ or $\rho^*$.
2. The Cartan part makes $K$ map each weight space to the weight space of the same weight.
3. The step terms of node $j$ shift weights by $\pm\alpha_j$. These $2(r+1)$ shifts are distinct unless $\hat g=a_1^{(1)}$ (where $\alpha_0=-\alpha_1$), so each node cancels by itself.
4. The factor $c_je^{\beta\alpha_j\cdot\phi_0/2}$ is common to both terms of a node and **drops out**. This leaves $K\rho(E_{\pm j})=\lambda^{\mp2}\rho'(E_{\pm j})K$ for all $j$.
5. $\ker K$ is a subrepresentation, so on its image $\rho'\simeq\rho\circ\gamma$, where $\gamma(E_{\pm j})=\lambda^{\pm2}E_{\pm j}$.
6. $E_{\alpha_0}\propto E_{-\theta}$ is a commutator of $h-1$ of the $E_{-\alpha_j}$. Consistency then forces $\lambda^2=\lambda^{-2(h-1)}$, i.e. $\lambda^{2h}=1$.
7. Antilinear σ with complex $\partial_x\phi$: $K$ must both commute and anticommute with $\vec H$, and the step terms then give $K=0$.
8. At real coupling with real fields and real λ, antilinear σ acts as a linear one.

**Numerics** (book `a_t`/`algebra`/`d4_rep`; representations built from scratch and validated):
- **No $K$:**
  - $a_2,a_3,a_4$ in the defining representation, for every combination: transposition or antipode, alone or composed with every $\mathbb Z_h$ rotation or the reflection. Also †/−conj on real data, generic $\phi_0$ and $\phi_0=0$.
  - Adjoint of $a_2$ and $a_3$, $\Lambda^2$ of $a_3$, $V\oplus\bar V$ of $a_2$, and the vector representation of $d_4^{(1)}$.
  - Complex Dirichlet at imaginary coupling (all four σ, $a_2$ and $a_3$).
- **Controls:**
  - Sinh-Gordon has a $K$ in every case: for linear σ only in the complex case; the adjoint of $a_1$ too.
  - The mechanism is confirmed: at $\lambda=e^{i\pi/h}$ a $K$ appears (antipode in the defining representation; both types in the self-conjugate adjoint and $\Lambda^2$). At $\lambda=e^{i\pi/2h}$ there is none.

**Independent of Sklyanin's formalism: the spin-±2 charges** (the method of Corrigan–Dorey–Rietdijk–Sasaki, hep-th/9404108 §2):
- **The current.** For $a_{2,3,4}$ the current $T=A(p,p,p)+B(q,p)$, $\Theta=-B(p,\nabla V)$ is unique. It is verified on an exact two-soliton solution, with 4th-order convergence down to relative $5\times10^{-8}$.
- **At a Dirichlet wall** ($u=\partial_x\phi(0)$):
  - The energy-like flux is $2A(u,u,u)-4B(u,\nabla V)$. It has no $\dot u$ term and is nonzero, so it cannot be a time derivative.
  - The momentum-like flux is $4B(\dot u,u)$. It is not a time derivative because $B$ is antisymmetric and nonzero.
- So no spin-±2 charge survives.

**Loopholes:**
- **$K$ depending on boundary data:** needs $\dot K=Ka_t+\sigma(a_t)K$. The right side contains no $\dot u$, so $K$ is constant. No loophole.
- **Other representations:** covered by the theorem.
- **A nonlocal $K$** (built from the time history): trivially makes τ constant but gives no local charges.
- **A dynamical boundary** with extra degrees of freedom: changes the boundary condition.
- **"Dirichlet at infinity"** (field at infinity at the wall, so some exponentials vanish): breaks the cycle, but is degenerate.
- **What survives: partial Dirichlet.** Take an involutive diagram automorphism $M$ and the images $\phi(x)=M\phi(-x)$: Dirichlet (=0) on the $M$-odd part, Neumann on the $M$-even part. These are integrable by the image argument. Full Dirichlet would need $M=-1$, which is a diagram symmetry only for $a_1^{(1)}$.
- **Sklyanin $K$ for those partial conditions:**
  - $a_3$: exists (half-turn: transposition; reflections: antipode).
  - $a_2$: none for the reflections in $V$, the adjoint or $V\oplus\bar V$. The reflection fixing node 0 is not an automorphism of the evaluated algebra for odd $h$, so the pairing $\lambda\leftrightarrow1/\lambda$ would have to be twisted.

**Gemini's derivation:**
- It assumes a λ-independent $K$ when matching powers of λ, and drops the constants.
- For λ-dependent $K$ the bracket $[E,E_-]$ produces no sign conflict. Its argument would also exclude sinh-Gordon Dirichlet, whose $K$ is λ-dependent.
- The worker's argument (ψ_j = λ² at every node, product $\lambda^{2N}=1$) is right.
- The coordinator's worry about $c_jc_j'=-1$ does not apply: after absorbing $A$ the constants cancel node by node.

**Prior art** (read at the source):
- **BCDR hep-th/9501098 and CDRS hep-th/9404108** assume $\partial_x\phi=-\partial\mathcal B/\partial\phi$: BCDR (1.2)–(1.4), CDRS (1.5)–(1.6) by count. Dirichlet is outside their ansatz, and the rigidity $|A_i|=2\sqrt{n_i}$ excludes it as a limit.
- **CDR hep-th/9407148:** no mention of Dirichlet.
- **Corrigan, hep-th/9612138** (Urumqi 1996 lectures; World Scientific 1998, pp. 9–32), §4: "setting the field to a specific value at the boundary will not be compatible with integrability… Dirichlet boundary conditions are not generally permitted", and "the Neumann condition is possible in all cases (but not Dirichlet)". This is the claim, stated without proof.
- **Doikou, JHEP 05 (2008) 091 (0803.0943), §4:** the antipode with the diagonal de Vega–González-Ruiz matrices for real $a_2$; the boundary conditions are deferred.
- **Avan–Doikou, NPB 821 (2009) 481 (0809.2734), §5:** for the diagonal $K$, Neumann (5.60) plus the singular constraint (5.61) $\Omega_2^{-2}\Omega_3^2(0)\to0$, i.e. the field at infinity. This is consistent with the obstruction.
- Nothing contradicts claim 5. The statement is known (Corrigan), but I found no published proof.

**Open problem 7(c), "Dirichlet walls":**
- A wall fixes the boundary value. Saddle points of the Euclidean integral keep it fixed, so the relevant complex classical solutions obey $\phi(0)=0$ in both real and imaginary parts. That is full Dirichlet, not the soliton-preserving condition, which leaves $\operatorname{Im}\phi(0)$ free.
- Real configurations satisfy both conditions. But for $n\ge2$ the equation of motion has a non-real force on real fields: $\operatorname{Im}=(m^2/\beta)\sum_j\alpha_j\cos(\beta\alpha_j\cdot\phi)\ne0$.
- **Implication:** classically the walls are not integrable for $a_n$, $n\ge2$. They are integrable for sine-Gordon (GZ Dirichlet at 0). Whether the quantum theory changes this is open.

---

## 3. Check script

`/tmp/claude-1000/-home-gustav-Git-affine-toda/b441568d-def3-432f-aa17-ea3eda7016d3/scratchpad/a3/boundary_dirichlet.py`, proposed as `textbook_code/boundary_dirichlet.py`.
- It has a docstring with book ids, conventions and 11 numbered checks, and reports through `boundary_common.report`.
- It execs the definition blocks of `boundary_classical.py` and `boundary_images.py`, which run their checks on import. Those modules are found beside the script, or via `PYTHONPATH`.
- Runtime about 2 min (2m16s).
- It also passed 11/11 in a simulated `textbook_code/` layout; that run used the revision before check 3 was extended, and the path logic is unchanged.

Output:
```
PASS sigma = dagger, Re phi(0) = (pi/beta) nu (mod (2pi/beta) Lambda), Neumann on Im phi: K_+ = e^{i pi nu.H} is compatible for all 2^n classes, K_+ = 1 only for nu in 2 Lambda (n = 1, 2, 3)  [residuals: D 1.7e-16; K = 1 off 2 Lambda >= 0.34]
PASS classification (sigma = dagger, these boundary data, n = 1, 2, 3): the constant K is unique up to scale and equals e^{i pi nu.H} for Re phi(0) in (pi/beta) Lambda; no K for Re phi(0) off (pi/beta) Lambda or without Neumann on Im phi
PASS e^{i pi nu.H} (all classes), and indeed any constant diagonal K_+, solves eq-classical-reflection with the book r-matrix at real spectral parameters, N = 2, 3, 4; a random non-diagonal constant K does not  [6.5e-16; random 1.17]
PASS static self-mirror solitons: phi(x) + conj phi(-x) = 2 Re phi(0) with beta Re phi(0)/pi a coweight outside 2 Lambda, d_x Im phi(0) = 0; e^{i pi nu.H} compatible, K = 1 not; Re E(x<0) = E/2, E real (half-soliton at the wall; a1 = sine-Gordon half-kink)  [a1 b=1 y=1.57: nu=[-0.5  0.5], E=9.877, Im E(x<0)=-0.000; a1 b=1 y=0.70: nu=[-0.5  0.5], E=9.877, Im E(x<0)=-5.863; a2 b=1 y=0.90: nu=[-0.333 -0.333  0.667], E=12.830, Im E(x<0)=28.576; a2 b=2 y=-1.30: nu=[ 0.333 -0.667  0.333], E=12.830, Im E(x<0)=16.097; a3 b=1 y=0.40: nu=[-0.25 -0.25  0.75 -0.25], E=13.968, Im E(x<0)=-6.783; a3 b=2 y=0.80: nu=[-0.5  0.5 -0.5  0.5], E=19.753, Im E(x<0)=-9.592]
PASS exact three-soliton image solutions (soliton + mirror + half-soliton; a1, a2, a3) obey the new condition at all times, and tau(lam) is conserved with K_+ = e^{i pi nu.H}, not with K_+ = 1; control: Theta P pair, K_+ = 1  [image 3.2e-14; rel. variation D 2.8e-08, K=1 >= 3.8e-02; control 3.8e-10]
PASS real coupling, Dirichlet on all components (generic phi_0 and phi_0 = 0): an invertible K exists for a1 and for no linear sigma (transposition, antipode, composed with all Z_h rotations and the reflection; dagger, -conj on real data) for a2, a3, a4  [nullities a1: [1]; a2: [0]; a3: [0]; a4: [0]]
PASS other representations, Dirichlet on all components: no K (transposition, antipode) in adjoint(a2, a3), Lambda^2(a3), V+Vbar(a2), vector of d4^(1); K exists for a1 adjoint  [a2 adjoint (0, 0); a3 adjoint (0, 0); a3 Lambda^2 (0, 0); a2 V+Vbar (0, 0); d4 vector (0, 0); a1 adjoint (1, 1)]
PASS imaginary coupling, complex Dirichlet phi(0) = phi_0 on all components: no K for any sigma (a2, a3); a1: linear sigma only  [a1 {'T': 1, 'S': 1, 'dag': 0, 'nconj': 0}; a2 {'T': 0, 'S': 0, 'dag': 0, 'nconj': 0}; a3 {'T': 0, 'S': 0, 'dag': 0, 'nconj': 0}]
PASS the obstruction is the cycle condition: K exists at lam = e^(i pi/h) (lam^(2h) = 1), not at lam = e^(i pi/2h)  [a2 defining ('S',): [1] vs [0]; a3 defining ('S',): [1] vs [0]; a2 adjoint ('T', 'S'): [1, 1] vs [0, 0]; a3 Lambda^2 ('T', 'S'): [1, 1] vs [0, 0]]
PASS spin +-2 charges (no Sklyanin): unique current for a2, a3, a4, conserved on an exact solution [rel. 2.1e-07]; at a Dirichlet wall the energy-like flux has no u-dot term and is nonzero, the momentum-like flux 4B(u-dot, u) is not a time derivative (B antisymmetric)  [a1: none; a2: dim 1, max |energy flux| 58.4, |A| 0.68, |B| 1.24; a3: dim 1, max |energy flux| 80.1, |A| 0.60, |B| 1.25; a4: dim 1, max |energy flux| 79.7, |A| 0.54, |B| 1.27]
PASS Delius hep-th/9809140 (5.1)-(5.2), with his image condition (1.4) for every coweight lambda: solved by K_+ = e^{i pi lambda.H}; K_+ = 1 (his proof) only for lambda in 2 Lambda, the case in which (1.4) implies his (1.5)  [3.6e-15; K = 1 off 2 Lambda >= 2.57]
```
The full-line energies equal $2hm_b/\beta^2$ to $10^{-6}$ (this is asserted in check 4).

---

## 4. Proposed book text

**(A) ch. 30, `sec-soliton-preserving`.** Replace the paragraph "The condition has no free parameter, apart from … imaginary part." with:

The condition has no free parameter. For sine-Gordon, where the physical field is real, it is the Dirichlet condition $\phi(0)\in\frac{2\pi}\beta\Lambda_W^\vee$ at a vacuum. For higher theories, where solitons are complex, it is the natural generalization of Dirichlet: Dirichlet on the real part, Neumann on the imaginary part. With $K_+\ne1$ the real part may also sit halfway between vacua.

::: {#prp-soliton-preserving-half}
**[Theorem]** Let $\vec\nu$ be a coweight. With $\sigma$ = Hermitian conjugation and $K_+=e^{i\pi\vec\nu\cdot\vec H}$, the compatibility condition (-@eq-boundary-compatibility) holds at imaginary coupling if

$$
\operatorname{Re}\vec\phi(0,t)\in\frac\pi\beta\vec\nu+\frac{2\pi}\beta\Lambda_W^\vee,\qquad\partial_x\operatorname{Im}\vec\phi(0,t)=0 .
$$ {#eq-soliton-preserving-half}

For $a_n^{(1)}$ in the defining representation, a field-independent $K_+$ compatible with Dirichlet on $\operatorname{Re}\vec\phi$ and Neumann on $\operatorname{Im}\vec\phi$ exists only if $\operatorname{Re}\vec\phi(0)\in\frac\pi\beta\Lambda_W^\vee$, and it is then a multiple of $e^{i\pi\vec\nu\cdot\vec H}$. Every constant diagonal $K_+$ solves (-@eq-classical-reflection) at real spectral parameters. As for $K_+=1$, the involution of the charges is not proved.
:::

::: {.proof}
The Cartan part of $K_+a_t(\lambda)+a_t(1/\lambda)^\dagger K_+$ is $\frac{i\beta}2[K_+,\partial_x\operatorname{Re}\vec\phi\cdot\vec H]-\frac\beta2\{K_+,\partial_x\operatorname{Im}\vec\phi\cdot\vec H\}$. At a Dirichlet boundary $\partial_x\operatorname{Re}\vec\phi$ is free, so $K_+$ commutes with $\vec H$; for $a_n^{(1)}$ in the defining representation it is diagonal, and the anticommutator vanishes iff $\partial_x\operatorname{Im}\vec\phi=0$. Write $K_+E_{\pm\vec\alpha_j}K_+^{-1}=\kappa_j^{\pm1}E_{\pm\vec\alpha_j}$ and $w_j=e^{i\beta\vec\alpha_j\cdot\vec\phi/2}$. The terms in $\lambda E_{\vec\alpha_j}$ and $\lambda^{-1}E_{-\vec\alpha_j}$ vanish iff $\kappa_j=\kappa_j^{-1}=\bar w_j/w_j=e^{-i\beta\vec\alpha_j\cdot\operatorname{Re}\vec\phi}$. So $\kappa_j=\pm1$ and $\operatorname{Re}\vec\phi\in\frac\pi\beta\Lambda_W^\vee$, and $\kappa_j=(-1)^{\vec\alpha_j\cdot\vec\nu}$ fixes $K_+$ up to a factor. In any algebra, $e^{i\pi\vec\nu\cdot\vec H}$ commutes with $\vec H$ and has $\kappa_j=(-1)^{\vec\alpha_j\cdot\vec\nu}=\bar w_j/w_j$. A constant diagonal $K_+$ commutes with the terms $e_{ik}\otimes e_{ki}$ and $e_{ii}\otimes e_{kk}$ of $r$ and $r^{t_1t_2}$, and it multiplies the matching terms of $r^{t_2}$ and $r^{t_1}$ by the same factor $k_ik_k$. So (-@eq-classical-reflection) reduces to the case $K_+=1$. (All three statements are checked for $n\le3$ in `textbook_code/boundary_dirichlet.py`.)
:::

**Half-solitons at the wall.** For $a_n^{(1)}$ the coweights modulo $2\Lambda_W^\vee$ label $2^n$ distinct boundaries, and $\vec\nu\in2\Lambda_W^\vee$ is (-@eq-soliton-preserving). In the other classes the boundary value is not a vacuum, and no constant field of finite energy satisfies the condition. The static solutions are halves of static solitons. A soliton of species $b$ at rest with $\operatorname{Re}\xi=0$ satisfies $\vec\phi(x)+\overline{\vec\phi(-x)}=2\operatorname{Re}\vec\phi(0)$. Here $\frac\beta\pi\operatorname{Re}\vec\phi(0)$ is its topological charge, a weight of the $b$-th fundamental representation, and the real part of its energy on $x\le0$ is half the soliton mass. Each class contains a weight of the representation $b$ and one of $h-b$. So there are two such half-solitons of equal mass, whose vacua at $x\to-\infty$ differ by the charge of a soliton. For sine-Gordon the class of the fundamental coweight is Ghoshal and Zamolodchikov's Dirichlet boundary at the maximum of the potential, $\phi_0=\pi/\beta_{\rm SG}$, with $K_+\propto\mathrm{diag}(1,-1)$; its two ground states are the half-kink and the half-antikink (@sec-boundary-sine-gordon). A soliton scattered off such a wall is a three-soliton image solution: the soliton, its $\Theta P$ mirror and the half-soliton (@sec-images-preserving). Along it $\tau(\lambda)$ is conserved with $K_+=e^{i\pi\vec\nu\cdot\vec H}$, and not with $K_+=1$ (checked for $a_1^{(1)}$, $a_2^{(1)}$ and $a_3^{(1)}$ in `textbook_code/boundary_dirichlet.py`).

**(B) ch. 30, after eq-theta-p-reflection.** Replace "In this form $\vec\lambda=-2\vec v\beta/2\pi$ must lie in … (checked in `textbook_code/boundary_classical.py`)." with:

At $x=0$ it gives $\operatorname{Re}\vec\phi(0)=-\frac\pi\beta\vec\lambda$. This is (-@eq-soliton-preserving) only for $\vec\lambda\in2\Lambda_W^\vee$. For the other coweights it is (-@eq-soliton-preserving-half) with $\vec\nu=-\vec\lambda$, which is integrable with $K_+=e^{i\pi\vec\lambda\cdot\vec H}$ and not with $K_+=1$ (`textbook_code/boundary_classical.py`, `textbook_code/boundary_dirichlet.py`). The image form (-@eq-theta-p-reflection) with $\vec v\in\frac\pi\beta\Lambda_W^\vee$ covers both.

**(C) ch. 30, `sec-durham-classification`, after the paragraph "Why are the $C_j$ so rigid?".**

**Dirichlet conditions.** For sine-Gordon the Dirichlet condition $\phi(0,t)=\phi_0$ is the limit $M_0\to\infty$ of the Ghoshal–Zamolodchikov boundary, and it is integrable for every $\phi_0$ (@sec-boundary-sine-gordon). For the other theories the $C_j$ cannot be sent to infinity, and Dirichlet conditions are not integrable at all, as Corrigan remarked without proof [@corrigan1996].

::: {#prp-no-dirichlet}
**[Theorem]** Let $\hat g=g^{(1)}\ne a_1^{(1)}$, and impose $\vec\phi(0,t)=\vec\phi_0$ on all components, at real coupling or, for complex fields, at imaginary coupling. For generic $\lambda$ no field-independent $K_+(\lambda)\ne0$ satisfies (-@eq-boundary-compatibility), when $\sigma$ is the antipode or the transposition composed with any conjugation, in any finite-dimensional representation without trivial summands. At imaginary coupling no antilinear $\sigma$ works either.
:::

::: {.proof}
At $x=0$ only $\partial_x\vec\phi$ is free. The conjugation in $\sigma$ can be absorbed into $K_+$, and since $\rho(X)^{\rm T}=-\rho^*(X)$, both cases read $K_+\rho(a_t(\lambda))=\rho'(a_t(1/\lambda))K_+$ with $\rho'=\rho$ or $\rho^*$. The Cartan part makes $K_+$ map weight spaces to weight spaces of the same weight. The step terms of node $j$ shift weights by $\pm\vec\alpha_j$, and for $\hat g\ne a_1^{(1)}$ these $2(r+1)$ shifts are distinct, so each cancels by itself. The factor $c_je^{\beta\vec\alpha_j\cdot\vec\phi_0/2}$ is common to both terms of a node and drops out, leaving $K_+\rho(E_{\pm\vec\alpha_j})=\lambda^{\mp2}\rho'(E_{\pm\vec\alpha_j})K_+$ for $j=0,\dots,r$. So the kernel of $K_+$ is a subrepresentation. On its image $\rho'$ is equivalent to $\rho\circ\gamma$, where $\gamma(E_{\pm\vec\alpha_j})=\lambda^{\pm2}E_{\pm\vec\alpha_j}$. Since $E_{\vec\alpha_0}\propto E_{-\vec\theta}$ is a commutator of $h-1$ generators $E_{-\vec\alpha_j}$, $j\ge1$, this requires $\lambda^{2h}=1$. For antilinear $\sigma$ and complex $\partial_x\vec\phi$, $K_+$ must both commute and anticommute with $\vec H$, and the step terms then give $K_+=0$. (Checked numerically for $a_2^{(1)}$, $a_3^{(1)}$, $a_4^{(1)}$ with all such $\sigma$, including compositions with the diagram automorphisms; in the adjoint, $\Lambda^2$ and $V\oplus\bar V$ representations; and in the vector representation of $d_4^{(1)}$. For sinh-Gordon a $K_+$ exists in every case; `textbook_code/boundary_dirichlet.py`.)
:::

The method of [@cdrs1994] gives an independent check that uses no $K_+$. For $a_n^{(1)}$, $n\ge2$, the lowest charges have spin $\pm2$, and no combination of them survives $\vec\phi(0,t)=\vec\phi_0$. The energy-like flux at $x=0$ is a nonzero cubic in $\partial_x\vec\phi(0)$, and the momentum-like flux is not a time derivative. [Optional, by the image argument: Dirichlet on part of the components can survive. If $M$ is an involutive symmetry of the affine Dynkin diagram, the image condition $\vec\phi(x)=M\vec\phi(-x)$, Dirichlet on the $M$-odd and Neumann on the $M$-even components, is integrable. Full Dirichlet would need $M=-1$, which is a diagram symmetry only for $a_1^{(1)}$.]

**(D) ch. 30, notes, last bullet.** Replace "In the image form (1.4) the coweight must lie in $2\Lambda_W^\vee$." with:

Its image form (1.4) allows every coweight, but (1.4) implies (1.5), and the proof with $K_+=1$ applies, only for $\vec\lambda\in2\Lambda_W^\vee$. For the other coweights (1.4) is integrable with $K_+=e^{i\pi\vec\lambda\cdot\vec H}$ (@prp-soliton-preserving-half), which also solves the paper's (5.1)–(5.2).

**(E) ch. 31, `sec-images-preserving`.** Add after "…with shifts as in (-@eq-multi-pair-shift).":

A soliton and its mirror carry the same charge, so these pairs keep $\vec v$ in $\frac{2\pi}\beta\Lambda_W^\vee$. For the boundaries of @prp-soliton-preserving-half with $\vec\nu\notin2\Lambda_W^\vee$, the image solution contains a third, self-mirror soliton at rest, $(b,0,\xi_3)$ with $\operatorname{Re}\xi_3=-\ln A_{ab}(\theta)$, and the mirror shift becomes $\xi_2=-\bar\xi_1-\ln A_{aa}(2\theta)-\ln A_{ab}(\theta)$. At phases where it is regular, it satisfies (-@eq-soliton-preserving-half) at all times (`textbook_code/boundary_dirichlet.py`).

In `tbl-boundary-types`, soliton-preserving row, vacuum column: "$\operatorname{Re}\vec\phi(0)$ at a vacuum; off the vacua (@prp-soliton-preserving-half) a static half-soliton".

**(F) ch. 35, open problem 7.** Replace "And are the Dirichlet walls … is open." with:

And are the Dirichlet walls of the constructive continuum limit (@sec-collapse-thresholds, @thm-window) integrable boundaries? A wall fixes the boundary value, and saddle points of the Euclidean integral keep it fixed. So its classical limit is $\vec\phi(0)=0$ on both the real and the imaginary part, not the soliton-preserving condition (-@eq-soliton-preserving), which leaves $\operatorname{Im}\vec\phi(0)$ free. Real configurations satisfy both, but for $n\ge2$ the imaginary-coupling field equation does not keep the field real. Dirichlet on all components is not classically integrable for $n\ge2$ (@prp-no-dirichlet). So classically the walls are integrable only for sine-Gordon; the quantum question remains open.

Optional addition to open problem 11: "The boundaries of @prp-soliton-preserving-half need reflection matrices as well; for $n=1$ the half class is Ghoshal–Zamolodchikov Dirichlet at $\phi_0=\pi/\beta_{\rm SG}$."

**(G) Appendix D, Delius 1998b.** Replace the (1.4) item with:

- **(1.4)–(1.5), the image form and the boundary condition.** The image condition (1.4), $\vec\phi(x)=-\overline{\vec\phi(-x)}-\frac{2\pi}\beta\vec\lambda$ with $\vec\lambda$ in the coweight lattice, gives $\operatorname{Re}\vec\phi(0)=-\frac\pi\beta\vec\lambda$. The paper calls (1.5) the corresponding boundary condition and proves it integrable with $K_+=1$. But (1.4) implies (1.5) only for $\vec\lambda\in2\Lambda_W^\vee$, which are also the only coweights produced by the reason the paper gives for the freedom, a shift by a vacuum. For the other coweights (1.4) is still integrable: the paper's own conditions (5.1)–(5.2) are solved by $K_+=e^{i\pi\vec\lambda\cdot\vec H}$, not by $K_+=1$. The boundary value is then not a vacuum, and the static solutions at the wall are half-solitons, as for sine-Gordon Dirichlet at the maximum of the potential. *Error* in the stated correspondence of (1.4) and (1.5); it could also be classed *Incomplete* (@prp-soliton-preserving-half; `textbook_code/boundary_dirichlet.py`).

**(H) `references.bib`:**
```
@incollection{corrigan1996, author = {Corrigan, E.}, title = {Integrable field theory with boundary conditions},
  booktitle = {Frontiers in Quantum Field Theory (Urumqi, 1996)}, editor = {Zha, C.-Z. and Wu, K.},
  publisher = {World Scientific}, address = {Singapore}, year = {1998}, pages = {9--32},
  eprint = {hep-th/9612138}, archivePrefix = {arXiv}, note = {arXiv:hep-th/9612138}}
```
If cited, also add Doikou, JHEP 05 (2008) 091 (doi 10.1088/1126-6708/2008/05/091) and Avan–Doikou, Nucl. Phys. B 821 (2009) 481 (doi 10.1016/j.nuclphysb.2009.05.010).

**(I) `textbook_code`:**
- README row: `boundary_dirichlet.py` | `prp-soliton-preserving-half`, `eq-soliton-preserving-half`, `prp-no-dirichlet`, `sec-images-preserving` | soliton-preserving boundaries at half coweights with $K_+=e^{i\pi\nu\cdot H}$; no full Dirichlet for $n\ge2$; spin-±2 test; Delius (5.1)–(5.2) | about two minutes.
- Suggested relabel in `boundary_classical.py`: "(pi/beta) w (half a coweight) fails *with K = 1*".

---

## 5. Book statements my results contradict

1. ch. 30, line 248: "must lie in $2\Lambda_W^\vee$: … which is not a vacuum, and the condition fails". This holds only for $K_+=1$; the stated reason is invalid.
2. ch. 30, line 240: "no free parameter, apart from the discrete choice of the vacuum value". Incomplete.
3. ch. 30, notes line 328: "the coweight must lie in $2\Lambda_W^\vee$".
4. Appendix D, Delius 1998b item (1.4): wrong reason and overstated conclusion (corrected item in §4(G)).
5. ch. 35, open problem 7: "a Dirichlet wall $\vec\phi=0$ satisfies the soliton-preserving condition entirely…". Literally true for real configurations, but misleading for the dynamics; the classical limit of a wall is full Dirichlet.
6. ch. 31, `tbl-boundary-types`, soliton-preserving row: "at a vacuum" is incomplete.
7. Outside the book: the brief's "known trap" (line 42), and Gemini's argument for claim 5.

Not contradicted: prp-soliton-preserving, ch. 33/35 "neither Dirichlet nor Neumann", `tbl-durham`.

**Lead for open problem 11 (not verified).** The reflection-type partial-Dirichlet images:
- preserve solitons and reflect particles $a\to h-a$, exactly as the soliton-preserving matrices do;
- correspond to the worker's coideal involution τ: j → l − j;
- reduce for $n=1$, $l=1$ to Dirichlet at 0.

They are candidates for the classical boundary of the de Vega–González-Ruiz family.

---

## 6. Files
All paths are under `/tmp/claude-1000/-home-gustav-Git-affine-toda/b441568d-def3-432f-aa17-ea3eda7016d3/scratchpad/a3/`.

**Deliverable:**
- `boundary_dirichlet.py`

**Exploration scripts** (`explore/`):
- `booklib.py`
- `e1_claim4.py`
- `e2_solitons.py`
- `e2b_singular.py`
- `e2c_regular_tau.py`
- `e3_claim5.py`
- `e4_spin2.py`
- `e4b_spin2_convergence.py`
- `e5_partial.py`
- `e5b_debug.py`
- `e5c_partial_reps.py`
- `e6_diag_reflection.py`

The last four were run inline with identical code and saved afterwards.

**Layout test:**
- `sim_textbook_code/` (copies of the script and the three book modules)

**Logs** (`logs/`):
- `boundary_dirichlet_final.log` (final run, 11 PASS)
- `boundary_dirichlet_run1.log`, `boundary_dirichlet_run2.log`
- `boundary_dirichlet_sim_textbook_code.log`
- `boundary_classical_run.log` (reference run of the book script, all PASS, 1m44s)
- `e1_claim4.log`
- `e2_solitons.log`
- `e2b_singular.log`
- `e2c_regular_tau.log`
- `e3_claim5.log`
- `e4_spin2.log`
- `e4b_spin2_convergence.log`
- `e5_partial.log`
- `e5b_debug.log`
- `e5c_partial_reps.log`
- `e6_diag_reflection.log`

**Papers** (`papers/`):
- `hep-th_9809140.tex` (Delius 1998b)
- `hep-th_9501098.tex` (BCDR)
- `hep-th_9404108.tex` (CDRS)
- `hep-th_9407148.tex` (CDR)
- `9612138/urumqir.tex` (Corrigan review)
- `arXiv_0803.0943.tex` (Doikou 2008)
- `ad0809/main.tex` (Avan–Doikou 2009)
- INSPIRE json files
