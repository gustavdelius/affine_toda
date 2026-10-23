# A6 report: problem 2 leftovers (the classical boundary of the d_n ε̂=0 amplitude, and the d_n s=±1 scalar factors)

Scratch directory: `/tmp/claude-1000/-home-gustav-Git-affine-toda/b441568d-def3-432f-aa17-ea3eda7016d3/scratchpad/a6/` (written `$A6/` below). I edited nothing in the repository and ran no git commands.

**Summary.**
- **Task 1 is done.** It is verified and extended.
  - **ε̂=0 is the uniform boundary C≡+1.** I showed this for d_4, d_5 and d_6. The classical limit of the book's K_1, and of CAF's \hat R_a for every particle, is the linearized reflection about the C≡+1 vacuum. Agreement is 3·10⁻¹¹.
  - **The vacuum is unique, and this is now proved.** For C_j≥0 the energy is strictly convex.
- **Task 2 is done for n=4, 5, 6 (scalar factors and breathers) and n=4..7 (continued amplitudes).**
  - The s=±1 matrices have exact unitarity and crossing functions, Gamma-product scalar factors and closed breather amplitudes.
  - The relevant sign is ς=(−1)ⁿs, not s.
  - **ς=+1 is the Neumann amplitude, but only with a CDD factor.** After the CDD factor that Kim–Cho's d_4 one-loop result requires, it is exactly the B→2−B dual of the uniform amplitude.
  - ς=−1 equals K_1/K_1^cl.
  - This is the a_n pattern.
- **Task 3: minimal.** The folding argument holds, and the classical c_2 check is consistent.

---

## 1. Verdicts

### Task 1: classical identification

Own code throughout. The d_n roots are those of `boundary_common.dn_vector`, and the potentials are the book's eq-boundary-potential with C_j=+1. The comparison S-matrix is the book's `smatrix_simply_laced.setup/S`.

**(1a) The invariant static kinks: VERIFIED.**
- **d_4.** ψ=f(e₂−e₃). This gives Tzitzéica, f''=2m²(e^{2f}−e^{−f}), with u=e^f, u'=m(u−1)√(2u(u+2)) and u(0)=9/8.
  - Closed form: u=2/(s²−1) with s=√3(1−E)/(1+E), E=E₀e^{√6 mx} and E₀=(3√3−5)/(3√3+5).
- **d_5.** ψ=f(e₂−e₄). This gives sinh-Gordon, f''=4m² sinh f, with tanh(f/4)=e^{2m(x−x₀)} and e^{f(0)}=4−2√2.
- Both satisfy the field equation and the boundary condition to 10⁻³⁰ (mpmath). Both are regular on x≤0, with f increasing from 0 to f(0).
- **Energies, now in closed form:**
  - d_4: E=(16√2−8−6√6) m/β² = −0.0695214587.
  - d_5: E=(32(1+√2)sin(π/8)−24−4√2) m/β² = −0.0927092091.
  - Both agree with the worker's numbers.
- Logs: `$A6/t1_vac.log`; the check is in the extended `dn_uniform.py`, check 7.

**(1b) "Is it the minimum among the candidates": PROVED, plus numerics.**
- **Proof.** For real C_j≥0 every term of V and of the boundary potential is a positive multiple of exp(linear function of φ), and the α_j span the root space. So the static energy ∫(½φ'²+V)+𝓑(φ(0)) is strictly convex on the convex space of finite-energy fields.
- Hence a static solution satisfying the boundary condition is the unique global minimum. A diagram symmetry preserving C≡+1 must fix it, which justifies the invariant ansatz.
- This also holds for a_n with C≡+1 or Neumann, where it says φ=0 is the unique vacuum.
- **Numerics.** I minimized the discretized energy (graded grid, analytic gradient, L-BFGS) in the full field space, without assuming the symmetry, from 6 random starts each.
  - d_4: −0.06952138. d_5: −0.09270913 (discretization error about 10⁻⁷).
  - The boundary values match the kinks.
  - Files: `$A6/t22_min.py`, `$A6/t22_min_d45.log`.
- The worker's "poorly converged lattice minimization" is superseded.

**(1c) Linearized reflection matrix about the vacuum: VERIFIED and EXTENDED.**
- **Method.** I integrate the ODE ε''+Ω²ε=M²(x)ε in matrix form from x=−L, with free waves in the mass eigenbasis; closed channels start with their decaying exponential, rescaled to O(1) at x=0. The boundary condition is ε'(0)=−Rε(0), with R=½Σ√n_j e^{α_j·ψ(0)/2}α_jα_jᵀ.
- **Results** (d_4, d_5; 5 values of Ω covering all channel openings):
  - K is diagonal in the particle basis: off-diagonal ≤5·10⁻¹⁵, including inside the degenerate d_4 triplet and the d_5 spinor pair.
  - For **every** particle, spinors included, K_a=∏_b S_ab(θ,B=1) to 3·10⁻¹².
  - For a=1 this equals the B→0 limit of eq-dn-uniform-amplitude, −(h/2−1)(h/2+1)/((1)(h−1)(h/2)²), to 3·10⁻¹².
  - Since CAF's \hat R_a(θ,0)=∏_b S_ab(θ,1) (their text after (Rhat); \tilde R→1), this identifies the B→0 limit of CAF's whole \hat R family.
- **Independent second route.** The vacuum was also found by shooting, with no closed form: same K to 1.4·10⁻¹¹.
- **d_6 (new).** The invariant sector is two-dimensional. The vacuum has tails in both invariant modes, m_2=1.6625 and m_4=2.6900 (m=1). It was found by shooting, seeded by the minimization.
  - Energy −0.1032515 m/β². Minimization gives −0.10325144, and the full-space minimization agrees.
  - Linearized reflection: diagonal, K_a=∏_b S_ab(θ,1) to 2·10⁻¹¹.
- Files: `$A6/cl_lin.py`, `$A6/t2_lin.py`, `$A6/t2_lin.log`, `$A6/t21_dn_vac.py`, `$A6/t21_d4.log`, `$A6/t21_d5.log`, `$A6/t21_d6.log`.

**(1d) New: no classical boundary bound states, and redundant poles.**
- K_1^cl has a physical-strip pole at θ=iπ(½−1/h). It is **not** a bound state.
  - The bound-state determinant has no zero in (0, m_min) for d_4 or d_5, and none at Ω=m_1 sin(π/h).
  - There are no negative modes, as convexity requires.
- Explanation: at that pole p=i m_2/2, where m_2 is exactly the decay rate of the vacuum's tail. It is Ma's redundant pole of an exponential tail (the Jost solution e^{ipx}(1+a e^{m_2x}) has a=V₁/(m_2²+2im_2p)).
- All physical-strip poles of ∏_bS_ab(θ,1) for d_4 and d_5 sit at p=i m_2/2.
- **Quantum consequence.** The coupling-independent pole of the book's K_1 (block (h/2−1)) is presumably not a boundary bound state. Its explanation (particle 1 absorbing a particle 2 from the condensate?) is open.
- Files: `$A6/t3_bound.py`, `$A6/t3_bound.log`; extended `dn_uniform.py`, check 8.

**(1e) "ε̂=0 is C≡+1 for d_4, d_5": VERIFIED (and for d_6), at 3·10⁻¹¹ in the classical limit.**
- What remains:
  - CDD factors that tend to 1 as B→0 are not excluded, as for a_n.
  - There are no one-loop data for this boundary.
  - n≥7 is not done (for n≥6 the vacuum has several tails and no closed form).
  - Other sign patterns are not done: for C_j<0 the convexity argument fails.

**(1f) Round-1 claims now computed rather than argued or assumed:**
- (i) At ε̂=0 the matrix part of crossing-unitarity is exactly 1: Ř_{V̄V̄}(2θ) fixes the boundary-state vector, to 3·10⁻¹⁴. The construction is that of `an_nondiagonal_scalar.py` check 6.
- (ii) (1⊗J)Ř(1⊗J) maps the B_1 singlet (rapidity difference i(π−π/T)) to itself with coefficient exactly 1, for n=4, 5 (n=6 in `$A6/t9_proj_d6.log`).
- Round-1 integration claim 2 ("verified by argument") and the assumption in prp-dn-scalar ("with the matrix part acting trivially") are now checks.
- The worker's claims that the d_5 spinors reflect diagonally with \hat R_s(B=0), and that K_1^cl=∏S_1j(θ,1): VERIFIED.

### Task 2: the d_n s=±1 scalar factors (new)

**Setup.**
- ε̂₀=sε_*, ε̂_{j≥1}=ε_*, with ε_*=1/(2cos(πω/2)) on the analytic branch (the book's K_raw-style principal square root flips all signs for ω>1, which is a gauge).
- A torus element can flip ε̂_j on any node with n_j=2 and on an even number of end nodes. So s is the only invariant.

**(2a) Exact K and the sign ς: VERIFIED.**
- I solved the book's intertwining equations symbolically in w=(−q)^{1/2} and y (sympy nullspace) for n=4, 5, 6. They agree with `boundary_common.k_intertwiner` to 10⁻¹⁴.
- In polynomial normalization: K₁₁=y^h(w²−1)(w^h−1) for both s, and the corner entry ∝(y^h−ς)(w^{h−2}y^h−ς) with **ς=(−1)ⁿs**.
- **Every scalar-factor datum depends on ς, not s.** For odd n the "all ε̂_j equal" matrix behaves like a_n's s=−1.
- Normalization used: K₁₁=−sin(πT/2)/sin πμ for ς=+1, and −i sin(πT/2)/cos πμ for ς=−1.
- Files: `$A6/t5_symK.py`; `$A6/t5_symK_d4.log`, `$A6/t5_symK_d5.log`, `$A6/t5_symK_d6.log`; `$A6/symK_d*.pkl`; `$A6/t6_unit.log`.

**(2b) Matrix-level data: VERIFIED for n=4, 5, 6 at 10⁻¹³–10⁻¹⁵** (book code: `k_intertwiner`, the `rcheck_w` solver, the Part VI coproduct):
- **Unitarity:** K(−θ)K(θ)=U_ς𝟙, with
  U_ς = ς(cos2πμ−ς cosπT)(cos2πμ−ς cosπω) / (2 sin²(πω/2)(cos2πμ−ς)).
  This is the a_n U_s times (cos2πμ−ς cosπω)/(cos2πμ−ς).
- **Crossing:** k(θ)=Ψ_ς Ř_{V̄V̄}(2θ)k(−θ), with the a_n functions Ψ_± unchanged. For ε̂=0, Ψ=1.
- **Singlet:** c_ς(μ)=c^{a_n}_ς(μ)·e_ς(μ), where
  - e_+ = cos π(μ−T/2+ω/2) cos π(μ−T/2−ω/2) / cos² π(μ−T/2);
  - e_− is the same with sin.
  - I found the zero and pole positions numerically (orders from winding numbers) and then checked them at random μ.
- Files: `$A6/t7_cross.py`, `$A6/t8_psi.py`, `$A6/t9_proj.py`, `$A6/t12_orders.py`, `$A6/t13_scan2.py`, `$A6/t14_cform.py`, with logs `t7_cross_d{4,5,6}.log`, `t8_psi_d{4,5,6}.log`, `t9_proj_d{4,5,6}.log`, `t14_d{4,5,6}.log`, `t12_d4_sp.log`, `t13_d4_sm.log`.

**(2c) Minimal Gamma-product solutions: DERIVED and VERIFIED to 10⁻²⁶–10⁻²⁹ (n=4, 5, 6).**
- A_ς = A₀(μ)·ρ_ς(μ)·g_ς(μ), where A₀ is the book's eq-dn-scalar.
- ρ_ς is the a_n factor of prp-an-nondiagonal-scalar evaluated at T=(n−1)ω (it depends only on T and ω):
  - ρ_+ = (−T/2)_I sin(πω/2)/sin π(μ−T/2) ∏_j G_j(μ,ω/2−3T/2)G_j(μ,1−ω/2−T/2) / (G_j(μ,ω/2−T/2)G_j(μ,1−ω/2−3T/2));
  - ρ_− = −i sin(πω/2)/cos π(μ−T/2) ∏_j[the same with ½ shifts];
  - G_j(μ,a)=Γ(μ+2jT+a)/Γ(−μ+2jT+a).
- g_ς(μ) = v_ς(μ)∏_{k≥1}[v_ς(μ+kT)/v_ς(−μ+kT)]^{(−1)^k}, with
  - v_+(x)=Γ(x+ω/2)Γ(1+x−ω/2)/(Γ(x)Γ(1+x));
  - v_−(x)=Γ(½+x+ω/2)Γ(½+x−ω/2)/Γ(½+x)².
- I derived g from g(μ)g(−μ)=U^{a_n}/U and g(T−μ)=g(μ). g has no poles in the strip for ω<2.

**(2d) Breather amplitudes, closed forms: VERIFIED to 10⁻²⁶ (n=4, 5, 6).** Here μ_{a,b}=μ∓(T/2−½).
- ς=−1:
  K_B = −(ω/2)_I(½−T)_I(−ω/2−½)_I(T/2+ω/2+½)_I(T/2)_I / ((T/2+½)_I(T/2+ω/2)_I).
- ς=+1 (minimal):
  K_B = (ω/2−T)_I(½−T)_I²(−ω/2−½)_I(T/2−ω/2)_I(T/2+ω/2+½)_I(T/2−½)_I / [(T/2)_I(½)_I(½−ω/2)_I(T+ω/2−½)_I(T/2−½+ω/2)_I(T/2+½−ω/2)_I].
- ς=+1 with the CDD factor Ξ=(T/2)_I²(T/2−ω/2)_I(T/2+ω/2)_I(ω/2)_I(T−ω/2)_I:
  K_B^N = (ω/2−T)_I(T/2−ω/2)_I(T/2−½)_I(½)_I(T−ω/2−½)_I / ((T/2)_I(T/2−ω/2−½)_I).
  This is the a_n Neumann form times (T/2−ω/2)_I(T/2−½)_I/((T/2)_I(T/2−ω/2−½)_I).
- Files: `$A6/dn_scal.py`, `$A6/t15_fe.log`, `$A6/t16_breather.log`, `$A6/t25_kbn.log`.

**(2e) Continuation ω→−2/B: VERIFIED for n=4..7.**
- K^{(−)}_1 = −(1)(h/2)(h/2+1−B/2) / [(1−B/2)(h+B/2)(h/2−B/2)(h/2+1)] = **K_1/K_1^cl**.
- K^N_1 = (h/2−1)(h/2+B/2)(h−1+B/2) / ((h−1)(h/2)(h/2−1+B/2)(B/2)) = **K_1|_{B→2−B}**, to 10⁻³⁰.
- The minimal ς=+1 amplitude without Ξ is in `$A6/t17_cont.py`.
- **Axioms:**
  - unitarity, and crossing-unitarity with the book's S_11, to 10⁻¹⁵;
  - bootstrap closure over all vector-like particles, for n=4..7, with and without Ξ, using my own fusing finder (forward poles with −i Res>0), validated on the book's uniform K_1;
  - all three tend to 1 as B→0.
- Files: `$A6/t17_cont.log`, `$A6/t19_boot.log`, `$A6/t24_dual.log`.

**(2f) Kim's d_4 one-loop amplitude.**
- **Source correction.** The d_4 one-loop Neumann calculation is **not** in Kim hep-th/9506031 (the book's kim1995). It is in **Kim and Cho, hep-th/9505138** (Phys. Rev. D 53 (1996) 4441). I read it at the arXiv TeX: eq. (29) for K_1, eq. (40) for K_2, eq. (42) for the exact self-dual conjecture; numbering checked against Kim's (33), (39), which the book already cites. hep-th/9506031 contains only the self-dual d_n conjectures (eq. (20) with fig. 7 for the vector particle).
- **Results at O(β²)** (B=β²/2π):
  - ς=+1 **with Ξ** reproduces Kim–Cho (29), and with K_2=K_1(θ+iπ/6)K_1(θ−iπ/6)S_11(2θ) also (40), to 2·10⁻¹⁷. Without Ξ it misses by 0.6.
  - ς=−1 (minimal) reproduces Kim–Cho at **θ+iπ**, to 6·10⁻¹⁸.
  - The a_n analogue holds too: the book's a_3 K^{(−)}_1(θ) equals Kim's (33) at θ+iπ (10⁻²⁷).
  - Against Kim's d_n vector conjecture: the same holds for n=5, 6, 7.
  - That conjecture equals CAF's \tilde R_1(θ+iπ) (CAF (dn1), (KimCF)) to 10⁻⁴⁰ for n=4..7, which confirms CAF's statement.
- **Beyond one loop**, K^N_1 ≠ Kim's self-dual amplitude (difference 8–15% at B=0.6), as for a_n.
- **The CDD multiplicities are fitted.** Within the family P_x=(T/2−x)_I(T/2+x)_I, the one-loop data fix the multiplicities (one each of x=0, ω/2, (T−ω)/2).
  - They do not distinguish P_{kω/2} from P_{½+kω/2}, which give the same O(β²).
  - So Ξ is fixed by Kim–Cho plus uniformity, not derived.
  - What favours it: with Ξ the result is exactly the dual of the uniform amplitude for every n, the d_n version of DG's (3.20).
- **Semiclassics does not discriminate.** Both ς give (β²/4h)d(−i log A)/dθ → ½ of the same quantity for F(2θ) (−0.3381 at θ=0.6, Richardson-extrapolated), i.e. the image picture of a Neumann wall.
  - The book's a_n A^{(−)} also reproduces the Neumann time advance at leading order.
- Files: `$A6/t17_cont.log`, `$A6/t18_semi.log`, `$A6/t20_kim.log`.

**Which s gives Neumann?** On the evidence, ς=+1 (s=+1 for d_4; s=−1, i.e. ε̂₀=−ε_*, for d_5), with Ξ, given cnj-reflection-continuation. ς=−1 is the analogue of a_n's s=−1 (K/K^cl), and no classical boundary is known for it.

**n=5 status.**
- Matrix data and scalar factors are verified.
- The one-loop test for n≥5 uses only Kim's conjecture, which he says was checked perturbatively but did not publish.

### Task 3: c_2

- **The folding argument is CONFIRMED.** c_2's S_22 (book eq-cn-snn, H=4+B) differs from a_3's {1}{3}: the ratio is 0.84−0.62i at B=0.4. The mass ratio floats (0.6616 at B=0.4, against 0.7071). So the a_3 bootstrap and crossing-unitarity cannot hold for folded amplitudes.
- **Classical check: done.** In c_2 itself (long roots of length² 2, short of length² 1, Kac labels (1,2,1)), the boundary matrix with √(2n_j/α_j²) is exactly M²/2.
  - The linearized factors at φ=0 are −[(1)(3)]^{−C} and −[(2)²]^{−C} with h=4, to 4·10⁻¹⁶. These are the folded a_3 targets.
  - This is a consistency check of the normalization in sec-images-cn. No c_2 quantum amplitude exists to compare with.
- Files: `$A6/t23_c2.py`, `$A6/t23_c2.log`.

---

## 2. Check scripts

### (A) Extended `dn_uniform.py`

File: `$A6/code_new/dn_uniform.py`.
- Book checks 1–5 are unchanged, apart from location-independent loading of `smatrix_simply_laced`.
- New checks 6–9: ε̂=0 matrix crossing factor and singlet coefficient; C≡+1 vacua of d_4, d_5; linearized reflection and absence of bound states; d_6 by shooting.
- Runtime about 3 min. Run with `PYTHONPATH=code` from any directory. Log: `$A6/dn_uniform_ext.log`.

```
PASS d_n, eps-hat = 0: K = J antidiagonal with entries +-1, independent of y, J^2 = 1 (n = 4, 5) [7e-15]
PASS d_n vector scalar factor F: unitarity and crossing F(T - mu) = F(mu) c(-mu) (n = 4, 5, 6) [8.4e-28]
PASS d_n boundary scalar factor A: unitarity, crossing A(T/2-mu)/A(T/2+mu) = F(2mu), A(0) = 1; breather amplitude = closed form (n = 4, 5, 6) [7.7e-27]
PASS d_n: continuation gives K_1, which satisfies unitarity and crossing-unitarity with the book d_n S-matrix and tends to the classical factor (n = 4..7) [5.1e-13]
PASS K_1 = Castro-Alvaredo-Fring's R^_1, and their spinor amplitude satisfies their crossing-unitarity with S_ss for n = 4..7 (no error for odd n) [2.3e-15]
PASS d_n, eps-hat = 0: R-check_(Vbar Vbar)(2 theta) fixes the boundary-state vector (crossing factor 1), and (1 x J) R-check (1 x J) maps the B_1 singlet (u = pi - pi/T) to itself with coefficient 1 (n = 4, 5) [3e-14]
PASS C = +1 vacuum: the invariant kinks of d_4 (Tzitzeica, e^f(0) = 9/8) and d_5 (sinh-Gordon, e^f(0) = 4 - 2 sqrt2) solve the static field equation and the boundary condition, are regular on x <= 0, and have energies 16 sqrt2 - 8 - 6 sqrt6 and -24 - 4 sqrt2 + 32(1 + sqrt2) sin(pi/8) (m/beta^2) [4.1e-30]
PASS linearized reflection about the C = +1 vacuum (d_4, d_5): diagonal in the particle basis, K_a = prod_j S_aj(theta, B=1) (book S-matrix) for every particle, spinors included, and K_1 = the B -> 0 limit of eq-dn-uniform-amplitude [3e-11; off-diagonal 2e-14]; no boundary bound state below threshold
PASS d_6, C = +1: the invariant static solution (two tail modes, m_2 and m_4) satisfies the boundary condition, and the linearized reflection about it is diagonal with K_a = prod_j S_aj(theta, B=1) [2e-11]; energy -0.103252 m/beta^2
```

### (B) New `dn_nondiagonal_scalar.py`

File: `$A6/code_new/dn_nondiagonal_scalar.py`.
- 9 checks, as listed in its docstring. Imports only from `code`: `boundary_common`, plus definitions loaded from `reflection_matrices.py` and `smatrix_simply_laced.py`.
- Runtime about 2.5 min. Log: `$A6/dn_nondiagonal_scalar.log`. The trailing "FIN" line in the log is my end marker; the script does not print it.

```
PASS d_n, eps-hat_j^2 = eps_*^2: K(-theta) K(theta) = U_varsigma 1, varsigma = (-1)^n s (n = 4, 5, 6) [5e-14]
PASS crossing-unitarity of the matrix part: k(theta) = Psi_varsigma R-check_(Vbar Vbar)(2 theta) k(-theta), Psi_+- as for a_n (n = 4, 5) [2e-14]
PASS B_1 singlet mapped to c_varsigma(mu) = c^an_varsigma(mu) e_varsigma(mu) times itself (n = 4, 5) [2e-13]
PASS scalar factors A_+-: unitarity A(mu)A(-mu)U = 1 and crossing A(T/2-mu)/A(T/2+mu) = F(2mu)/Psi (n = 4, 5, 6) [1.5e-26]
PASS breather amplitudes: A(mu_a) F(2mu) A(mu_b) c(mu) = closed block forms (n = 4, 5, 6) [1.4e-26]
PASS omega -> -2/B: K^(+-)_1 = closed forms, unitary, crossing-unitary with the book S_11 (n = 4..7), bootstrap closure over the vector-like particles (n = 5, 6), K -> 1 as B -> 0 [2.4e-15]
PASS one loop, d_4 (Kim-Cho (K1), (K2)): varsigma = +1 with the CDD factor agrees [1.9e-17] (without it: 0.61); varsigma = -1 agrees at theta + i pi [5.6e-18] (at theta: 0.31); the same against Kim's d_n conjecture = CAF R~_1(theta + i pi), n = 4..7 [1.6e-26]; a_3: K^(-)_1(theta) = Kim (33)(theta + i pi) [1.7e-27]
PASS semiclassics (d_4): both varsigma = +-1 give half the time delay of F(2 theta), as a Neumann wall by images [-0.338016, -0.338097, -0.338162]
PASS varsigma = +1 with the CDD factor is the dual K_1(theta, 2 - B) of the uniform amplitude, = K^N_1 in closed form; varsigma = -1 is K_1/K_1^cl (n = 4..7) [1.3e-30]
```

(The script's labels "(K1), (K2)" are Kim–Cho's TeX labels, printed as eqs. (29), (40).)

### Proposed README rows (`code/README.md`)

- `dn_uniform.py`: add ids `prp-dn-uniform-vacuum` and `prp-dn-uniform-classical`. Add to "what it checks": "ε̂=0 crossing factor and singlet coefficient 1; C≡+1 vacua of d_4, d_5 (closed form) and d_6 (shooting); linearized reflection = ∏_b S_ab(θ,1); no classical bound states". Runtime about 3 min.
- New row: `dn_nondiagonal_scalar.py` | `prp-dn-nondiagonal-scalar`, `eq-dn-nondiagonal-scalar`, `eq-dn-neumann-amplitude`, `sec-dn-uniform` | "d_n s=±1 matrices: U_ς, Ψ_ς, singlet c_ς (ς=(−1)ⁿs); scalar factors and breathers to 30 digits; continued amplitudes, axioms, bootstrap; Kim–Cho one loop; K^N_1 = dual of uniform; semiclassics" | about 2.5 min.

---

## 3. Proposed book text

**(a) ch. 30, `sec-boundary-potential`, line 53.** Replace "For $a_n^{(1)}$ such solutions are constructed in @sec-real-coupling-vacua; for $d_n^{(1)}$ none is known." with:

> For $a_n^{(1)}$ such solutions are constructed in @sec-real-coupling-vacua. For $d_n^{(1)}$ with $C\equiv+1$ the vacuum is a static kink, unique because the energy is convex (@prp-dn-uniform-vacuum).

**(b) ch. 33, `sec-dn-vector-k`.**

*Structural fix.* Lines 210–232: the paragraph "**The scalar factor at $\hat\epsilon=0$.**" through "…continues to is in @sec-dn-uniform." sits between `::: {#prp-dn-r}` and that proposition's title ("## The $d_n$ vector R-matrix"), so prp-dn-scalar is nested inside prp-dn-r. Move the paragraph after the paragraph that follows @prp-dn-k, because it uses $K=AJ$ before the matrices are introduced. Restore prp-dn-r as `::: {#prp-dn-r}` / `## The $d_n$ vector R-matrix` / `**[Numerical]** …` / `:::`.

*In prp-dn-scalar.* Replace "With the matrix part acting trivially on the singlet, the lowest-breather amplitude is" by "The matrix part maps the singlet of $V\otimes V$ at the pole $\langle n-1\rangle$ to itself with coefficient 1, and the lowest-breather amplitude is". After "contributes no factor." add "(computed in `code/dn_uniform.py`)".

*New, after it:*

> **The matrices with $\hat\epsilon_j^2=\epsilon_*^2$.** Fix the gauge as for $a_n$: $\hat\epsilon_0=s\epsilon_*$, $\hat\epsilon_{j\ge1}=\epsilon_*$, $\epsilon_*=1/(2\cos\frac{\pi\omega}2)$. A torus element flips $\hat\epsilon_j$ on any node with $n_j=2$, and on an even number of end nodes, so $s$ is the only invariant. The scalar factor depends on $\varsigma=(-1)^ns$: for odd $n$ the matrix with all $\hat\epsilon_j$ equal behaves like the $a_n$ matrix with $s=-1$. Normalize $K$ by $K_{11}=-\sin\frac{\pi T}2/\sin\pi\mu$ for $\varsigma=+1$ and $-i\sin\frac{\pi T}2/\cos\pi\mu$ for $\varsigma=-1$.
>
> ::: {#prp-dn-nondiagonal-scalar}
> ## Scalar factors of the $d_n$ matrices with $\hat\epsilon_j^2=\epsilon_*^2$
> **[Numerical]** For $n=4,5,6$: (i) $K(-\theta)K(\theta)=U_\varsigma\mathbf 1$ with
> $$U_\varsigma=\varsigma\,\frac{(\cos2\pi\mu-\varsigma\cos\pi T)(\cos2\pi\mu-\varsigma\cos\pi\omega)}{2\sin^2\frac{\pi\omega}2\,(\cos2\pi\mu-\varsigma)},$$
> the $a_n$ factor of @sec-an-nondiagonal-k times $(\cos2\pi\mu-\varsigma\cos\pi\omega)/(\cos2\pi\mu-\varsigma)$; (ii) the matrix part of crossing-unitarity is the $a_n$ function $\Psi_\varsigma$; (iii)
> $$A_\varsigma(\mu)=A(\mu)\,\rho_\varsigma(\mu)\,g_\varsigma(\mu)$$ {#eq-dn-nondiagonal-scalar}
> with $A$ of (-@eq-dn-scalar) solves $A_\varsigma(\mu)A_\varsigma(-\mu)U_\varsigma=1$ and $A_\varsigma(\frac T2-\mu)/A_\varsigma(\frac T2+\mu)=F(2\mu)/\Psi_\varsigma(\mu)$, where
> $$\rho_+=\Big(-\frac T2\Big)_I\frac{\sin\frac{\pi\omega}2}{\sin\pi(\mu-\frac T2)}\prod_{j\ge1}\frac{G_j(\mu,\frac\omega2-\frac{3T}2)\,G_j(\mu,1-\frac\omega2-\frac T2)}{G_j(\mu,\frac\omega2-\frac T2)\,G_j(\mu,1-\frac\omega2-\frac{3T}2)},\qquad \rho_-=\frac{-i\sin\frac{\pi\omega}2}{\cos\pi(\mu-\frac T2)}\prod_{j\ge1}\frac{G_j(\mu,\frac12+\frac\omega2-\frac{3T}2)\,G_j(\mu,\frac12-\frac\omega2-\frac T2)}{G_j(\mu,\frac12+\frac\omega2-\frac T2)\,G_j(\mu,\frac12-\frac\omega2-\frac{3T}2)},$$
> $G_j(\mu,a)=\Gamma(\mu+2jT+a)/\Gamma(-\mu+2jT+a)$, and
> $$g_\varsigma(\mu)=v_\varsigma(\mu)\prod_{k\ge1}\Big(\frac{v_\varsigma(\mu+kT)}{v_\varsigma(-\mu+kT)}\Big)^{(-1)^k},\quad v_+(x)=\frac{\Gamma(x+\frac\omega2)\Gamma(1+x-\frac\omega2)}{\Gamma(x)\Gamma(1+x)},\quad v_-(x)=\frac{\Gamma(\frac12+x+\frac\omega2)\Gamma(\frac12+x-\frac\omega2)}{\Gamma(\frac12+x)^2};$$
> (iv) the singlet is mapped to $c_\varsigma$ times itself, with $c_+$ the $c_+$ of @prp-an-nondiagonal-scalar times $\cos\pi(\mu-\frac T2+\frac\omega2)\cos\pi(\mu-\frac T2-\frac\omega2)/\cos^2\pi(\mu-\frac T2)$, and $c_-$ the $c_-$ there times the same with sines; (v) with $\mu_{a,b}=\mu\mp(\frac T2-\frac12)$,
> $$A_-(\mu_a)F(2\mu)A_-(\mu_b)c_-=-\Big(\frac\omega2\Big)_I\Big(\frac12-T\Big)_I\Big(-\frac\omega2-\frac12\Big)_I\frac{(\frac T2+\frac\omega2+\frac12)_I(\frac T2)_I}{(\frac T2+\frac12)_I(\frac T2+\frac\omega2)_I},$$
> and with $A_+$ multiplied by $\Xi=(\frac T2)_I^2(\frac T2-\frac\omega2)_I(\frac T2+\frac\omega2)_I(\frac\omega2)_I(T-\frac\omega2)_I$,
> $$\Xi A_+(\mu_a)F(2\mu)\,\Xi A_+(\mu_b)\,c_+=\Big(\frac\omega2-T\Big)_I\Big(\frac12\Big)_I\Big(T-\frac\omega2-\frac12\Big)_I\frac{(\frac T2-\frac\omega2)_I(\frac T2-\frac12)_I}{(\frac T2)_I(\frac T2-\frac\omega2-\frac12)_I}.$$
> :::
>
> Checked in `code/dn_nondiagonal_scalar.py`: (i), (ii), (iv) at the matrix level to $10^{-13}$; (iii), (v) to $10^{-26}$. The CDD factor $\Xi$ is unitary and crossing-symmetric by itself. It is not fixed by the functional equations; it is the factor that Kim and Cho's one-loop $d_4$ amplitudes require (@sec-dn-uniform). It plays the role of the factor $(\frac T2\mp\frac\omega2)_I$ that Kim's $a_3$ amplitudes require (@sec-an-nondiagonal-k).

*Replace the last paragraph (line 255):*

> So $d_n$ has the same three solutions as $a_n$, and they describe the same boundaries. The $\hat\epsilon=0$ matrix is the uniform boundary $C\equiv+1$, whose vacuum is a static kink rather than $\vec\phi=0$ (@sec-dn-uniform). The matrix with $\varsigma=+1$, with the CDD factor fixed at one loop, gives the dual of the uniform amplitude, the $d_n$ form of the Neumann amplitude. The matrix with $\varsigma=-1$ has no known classical boundary.

**(c) ch. 35, `sec-dn-uniform`.** Replace the last two sentences of line 302 ("Which classical boundary … should be $C\equiv+1$."):

> The vacuum of the uniform boundary can be found exactly.
>
> ::: {#prp-dn-uniform-vacuum}
> ## The uniform vacuum of $d_n^{(1)}$
> **[Theorem]** At real coupling with real $C_j\ge0$, a static field of finite energy that satisfies (-@eq-boundary-condition) is unique and is the vacuum. For $C\equiv+1$ it is invariant under the symmetries of the extended Dynkin diagram. For $d_4^{(1)}$ it is $\beta\vec\phi=f(\vec e_2-\vec e_3)$ with $u=e^f$, $u'=m(u-1)\sqrt{2u(u+2)}$ and $u(0)=\frac98$, a Tzitzéica kink. For $d_5^{(1)}$ it is $\beta\vec\phi=f(\vec e_2-\vec e_4)$ with $\tanh\frac f4=e^{2m(x-x_0)}$ and $e^{f(0)}=4-2\sqrt2$, a sinh-Gordon kink. Both are regular on $x\le0$, with energies $(16\sqrt2-8-6\sqrt6)\frac m{\beta^2}$ and $\big(32(1+\sqrt2)\sin\frac\pi8-24-4\sqrt2\big)\frac m{\beta^2}$.
> :::
>
> ::: {.proof}
> Each term of $V$ is a positive multiple of the exponential of a linear function of $\vec\phi$, and the $\vec\alpha_j$ span the root space. For $C_j\ge0$ the same holds for $\mathcal B$. So the static energy $\int(\frac12\vec\phi'^2+V)dx+\mathcal B(\vec\phi(0))$ is strictly convex on the convex set of finite-energy fields. Its critical points, the static solutions, are therefore global minima, and there is at most one. A diagram symmetry maps the problem with $C\equiv+1$ to itself, so it fixes the solution. On the invariant fields of $d_4^{(1)}$ the exponents are $-f$ on the four end nodes and $2f$ on the centre, and the field equation becomes $f''=2m^2(e^{2f}-e^{-f})$. For $d_5^{(1)}$ the exponents are $\mp f$ and the equation is $f''=4m^2\sinh f$. The first integral and the boundary condition fix $f(0)$, and the energies follow by quadrature.
> :::
>
> Checked in `code/dn_uniform.py`. A direct minimization of the discretized energy, from random fields and without the symmetry, reaches the same values. For $d_6^{(1)}$ the invariant sector is two-dimensional, and the vacuum, found by shooting, has tails in both invariant directions, $m_2$ and $m_4$.
>
> ::: {#prp-dn-uniform-classical}
> **[Numerical]** For $d_4^{(1)}$, $d_5^{(1)}$ and $d_6^{(1)}$, the linearized reflection about the $C\equiv+1$ vacuum is diagonal in the particle basis, and for every particle $a$, spinors included, $K_a^{\rm cl}(\theta)=\prod_bS_{ab}(\theta)\big|_{B=1}$. For $a=1$ this is the limit above. The linear problem has no bound state.
> :::
>
> Checked in `code/dn_uniform.py` to $10^{-11}$. So the $\hat\epsilon=0$ matrix is the uniform boundary, as for $a_n^{(1)}$; the kink adds the factor $(\frac h2-1)(\frac h2+1)/(\frac h2)^2$. The right-hand side is Castro-Alvaredo and Fring's $\hat R_a$ at $B=0$, so the classical limit of their whole family is the uniform boundary. The pole of $K_1$ at $\theta=i\pi(\frac12-\frac1h)$ is present at every $B$. Classically it is not a bound state: there $p=\frac i2m_2$, and the pole is the redundant pole of the $e^{m_2x}$ tail of the vacuum.
>
> **The other two matrices.** Continued by $\omega\to-2/B$, the breather amplitudes of @prp-dn-nondiagonal-scalar become
> $$K^N_1=\frac{(\frac h2-1)(\frac h2+\frac B2)(h-1+\frac B2)}{(h-1)(\frac h2)(\frac h2-1+\frac B2)(\frac B2)}=K_1\big|_{B\to2-B},\qquad K^{(-)}_1=\frac{K_1}{K_1^{\rm cl}} .$$ {#eq-dn-neumann-amplitude}
> Both are unitary and crossing-unitary with $S_{11}$, close the bootstrap over the vector-like particles, and tend to 1, the Neumann value, as $B\to0$ ($n\le7$). At $O(\beta^2)$, $K^N_1$ and $K^N_2(\theta)=K^N_1(\theta+\frac{i\pi}6)K^N_1(\theta-\frac{i\pi}6)S_{11}(2\theta)$ of $d_4^{(1)}$ are Kim and Cho's one-loop Neumann amplitudes [@kim_cho1996, eqs. (29), (40)]. $K^{(-)}_1$ agrees with them only after $\theta\to\theta+i\pi$. The same holds against Kim's self-dual conjecture for $d_5$, $d_6$, $d_7$ [@kim1995, eq. (20)], which is Castro-Alvaredo and Fring's $\tilde R_1(\theta+i\pi)$. This is the $a_n$ pattern: (-@eq-neumann-amplitude) is the dual of the uniform amplitude, $K^{(-)}_a=K_a/K_a^{\rm cl}$ in @sec-neumann-duality, and $K^{(-)}_1$ of $a_3$ is Kim's amplitude at $\theta+i\pi$. So $\varsigma=+1$ is the Neumann boundary, given @cnj-reflection-continuation and the CDD factor $\Xi$ that one loop fixes. Neither the classical limit nor the semiclassical soliton time delay distinguishes the two matrices: both give half the phase of $F(2\theta)$, as a Neumann wall does by images (`code/dn_nondiagonal_scalar.py`).

**(d) Open problem 2 (line 342).** Replace the d_n sentences:

> For $d_n^{(1)}$ all three vector reflection matrices now have scalar factors and breather amplitudes (@prp-dn-scalar, @prp-dn-nondiagonal-scalar). $\hat\epsilon=0$ is the uniform boundary $C\equiv+1$, whose vacuum is a static kink (@prp-dn-uniform-vacuum, @prp-dn-uniform-classical, checked for $n\le6$), and its continuation is Castro-Alvaredo and Fring's $\hat R_1$. The matrix with $(-1)^ns=+1$ gives the dual of the uniform amplitude, which reproduces Kim and Cho's one-loop $d_4^{(1)}$ Neumann amplitudes [@kim_cho1996] with a CDD factor chosen at one loop (@sec-dn-uniform). Still open for $d_n^{(1)}$: a derivation of that CDD factor; the vacua for $n\ge7$ and for the other sign patterns; the spinor amplitudes; a classical boundary for $(-1)^ns=-1$; and the quantum origin of the coupling-independent pole of $K_1$ at $\theta=i\pi(\frac12-\frac1h)$, which classically is a redundant pole of the vacuum's tail.

Keep the c_n and CAF sentences. The CAF clause becomes "…apart from the $\hat R_a$ of $d_n^{(1)}$ (uniform), which boundary conditions they describe is not known."

**(e) `references.bib`.** Add:

```
@article{kim_cho1996,
  author = {Kim, J.D. and Cho, H.S.},
  title = {Boundary reflection matrix for {$D_4^{(1)}$} affine {Toda} field theory},
  journal = {Phys. Rev. D}, volume = {53}, pages = {4441}, year = {1996},
  eprint = {hep-th/9505138}, archivePrefix = {arXiv}
}
```

The DOI is 10.1103/PhysRevD.53.4441, from arXiv metadata. The equation numbers are from the arXiv version.

**(f) Optional, ch. 35 line 166.** After "The scalar factor reproduces the classical Neumann time advance (@sec-boundary-semiclassics)" add "(as does $A^{(-)}$, so this does not single out $s=+1$)".

---

## 4. Statements in the book or Appendix D that my results contradict or supersede

1. **ch. 33, lines 210–232: markup error.** prp-dn-scalar and its paragraph are inside the `prp-dn-r` div, before its title, and they use K=AJ before prp-dn-k is introduced.
2. **ch. 33, line 255; ch. 35, line 302; ch. 30, line 53; open problem 2:** "which classical d_n boundary … is open / should be C≡+1"; "for d_n none is known"; "scalar factors of the s=±1 matrices" open. All are superseded.
3. **Open problem 2 cites [@kim1995] for "Kim's one-loop d_4 Neumann amplitude".** That calculation is Kim and Cho, hep-th/9505138. kim1995 (hep-th/9506031) has only the self-dual d_n conjectures.
4. **prp-dn-scalar's hedge** "with the matrix part acting trivially on the singlet" is now computed: the coefficient is exactly 1.
5. **Implicit in "same three solutions … s=±1" (ch. 33):** for odd n, the a_n label s corresponds to ς=−s. The all-equal matrix is the a_n-s=−1 analogue.
6. **ch. 35, line 166, and the book's semiclassical evidence for s=+1:** a_n's A^{(−)} has the same leading semiclassical limit (`$A6/t18_semi.log`: a_2 and a_3 A^{(−)} → the Neumann value). This is not wrong as stated, but it does not discriminate between the two matrices.
7. **Appendix D:** nothing contradicted.
   - CAF's (KimCF), Kim = \tilde R(θ+iπ), is confirmed for d_n vector particles, n=4..7.
   - The appendix's a_n remark "Kim's self-dual amplitudes agree at one loop, not beyond" extends to d_n (n=4..7, against K^N_1).
   - Round-1's CAF-odd-n refutation still stands (book check 5 passes).

---

## 5. Files (all under `$A6/`)

**Check scripts:**
- `code_new/dn_uniform.py`, log `dn_uniform_ext.log`
- `code_new/dn_nondiagonal_scalar.py`, log `dn_nondiagonal_scalar.log`

**Task 1:**
- Classical machinery: `cl_lin.py`
- Vacua: `t1_vac.py` + `t1_vac.log`
- Linearized reflection: `t2_lin.py` + `t2_lin.log`, `t2_lin_d5.log`
- Bound states: `t3_bound.py` + `t3_bound.log`
- Shooting vacua: `t21_dn_vac.py` + `t21_d4.log`, `t21_d5.log`, `t21_d6.log`
- Direct minimization: `t22_min.py` + `t22_min_d45.log`, `t22_min_d6.log`, `min_d*.npy`, `minprof_d6_*.npy`

**Task 2:**
- Book-code helpers: `qk.py`
- Symbolic K: `t5_symK.py` + `t5_symK_d4.log`, `t5_symK_d5.log`, `t5_symK_d6.log`, `symK_d4.pkl`, `symK_d5.pkl`, `symK_d6.pkl`
- Unitarity: `t6_unit.py` + `t6_unit.log`
- Crossing: `t7_cross.py` + `t7_cross_d4.log`, `t7_cross_d5.log`, `t7_cross_d6.log`; `t8_psi.py` + `t8_psi_d4.log`, `t8_psi_d5.log`, `t8_psi_d6.log`
- Singlet: `t9_proj.py` + `t9_proj_d4.log`, `t9_proj_d5.log`, `t9_proj_d6.log`; `t10_fitc.py`; `t11_scan.py`; `t12_orders.py` + `t12_d4_sp.log`; `t13_scan2.py` + `t13_d4_sm.log`; `t14_cform.py` + `t14_d4.log`, `t14_d5.log`, `t14_d6.log`
- Scalar factors: `dn_scal.py`; `t15_fe.py` + `t15_fe.log`
- Breathers: `t16_breather.py` + `t16_breather.log`
- Continuation and one loop: `t17_cont.py` + `t17_cont.log`
- Semiclassics: `t18_semi.py` + `t18_semi.log`
- Bootstrap: `t19_boot.py` + `t19_boot.log`
- Kim and CAF comparison: `t20_kim.py` + `t20_kim.log`
- Duality relations: `t24_dual.py` + `t24_dual.log`
- Dressed breather: `t25_kbn.py` + `t25_kbn.log`

**Task 3:** `t23_c2.py` + `t23_c2.log`

**Sources read (arXiv TeX):**
- `dl_kim/kim.tex` (hep-th/9506031)
- `dl_kimcho/kc.tex` (hep-th/9505138)
- `dl_caf/x/Matrixreloaded.tex` (hep-th/0308004)

**Scratch debugging, not needed:** `dbg*.py`, `dbg8.out`, `dbg9.out`.

**Process note.** One early background job was killed by the system for low memory. I did not restart it in the background; I reran the work in the foreground after fixing an ODE-tolerance inefficiency, and all logs listed are from completed runs.
