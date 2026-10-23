# Ch. 42 drafting checkpoint (defects meet boundaries)

Scratch dir: /tmp/claude-1000/-home-gustav-Git-affine-toda--claude-worktrees-integrable-defects-framework-a54a84/cfddf2dd-7e99-4c1c-b8f5-c799ea33aced/scratchpad/ch42/

## Status (2026-10-12)
- [x] Read conventions, status, ch. 36, 40, 34 (model), 29 (GZ, shG), 37 fusion section, sources digest, cz2012 text, Bajnok-Simon text, Bao-Wang Rem. 4.10(2), Wang ICM survey.
- [x] Fetched SY 2504.08375 and AV 2402.16676 PDFs (WebFetch saved them; pdftotext in scratch ch42/sy2025.txt, av2024.txt).
- [x] Scratch numerics (a)-(d) done; (e) BDG by hand + to script; (f) tensor K to do in script.
- [x] Script code/defect_doubling.py: 17 PASS, 0 FAIL, ~30 s (output copy: scratch ch42/script_output.txt).
- [x] Chapter text part8-defects/42-defects-and-boundaries.qmd (241 lines; all refs and bib keys resolve).
- [x] Report briefs/draft_ch42_report.md (final).  If resumed: only proofreading left.

## Decisions / plan
- Reflected copy: V'(u) := V(-u) as a rep with U acting through Delta^op (U^cop). PRcP intertwines Delta^op.
  Via omega_q (Hopf iso U^cop -> U), V(1/y) o omega_q = Vbar(y) exactly (matrices), so the reflected copy is the
  conjugate multiplet with the standard coproduct; PRcP ~ Rc_{Vbar Vbar}. For a_2: reflected soliton = antisoliton.
- Doubled algebra A = U (x) U^cop, Delta_A(g(x)h) = (g1 (x) h2) (x) (g2 (x) h1).  Delta(U) is a LEFT coideal of A
  (coassociativity), = Bao-Wang's (omega (x) 1)Delta(U) in U(x)U after 1(x)omega and factor swap.
  Delta(U_q(b^-+)) is also a left coideal (Hopf subalgebra) = "diagonal Borel".
- Action of g(x)h on M = V + V': diag(pi(g) eps(h), eps(g) pi'(h)); implement sigma1 = pi + eps, sigma2 = eps + pi'.
- Transparent defect test (b): search space X = span{g(x)h}, g,h in {k^a, e k^a, f k^a}, a in {-1,0,1}, total e/f degree <= 1
  (45 elements per node). Impose single-particle K-commutation and two-particle doubled reflection intertwining;
  compare operator image dimension with that of Delta(U) cap X (expect 7: 1, k(x)k, k^-1(x)k^-1, Delta e, Delta(k^-1 e),
  Delta f, Delta(f k)).  Single-particle condition alone admits e(x)1 + 1(x)e (fails on mixed V(x)V' states).
- Jump defect test (c): diagonal Borel holds by construction (= eq-borel-intertwiner).  Full twisted diagonal
  (1(x)psi)Delta(U), psi = torus and/or diagram swap tau: least squares over free theta-independent Y = chi(x) and
  kappa = chi(k(x)psi k) on W window:  T(pi(e)(x)1 + pi(k)(x)Y) = (Y(x)1 + kappa(x)pi(psi e)) T.  Expect O(1) residual;
  controls: transparent defect and T = Rc_{V,W}(y/y_d) (residual 0); limit theta - theta_d -> -inf (residual -> 0).
- Fusion (d): K_f = Tbar(theta) (1_W (x) K_D(theta)) T(theta) : V(theta)(x)W -> V(-theta)(x)W, Tbar(theta) = T(-theta)^-1.
  cz2012 (2.18): type I + Dirichlet = GZ with alpha-independent entries after a similarity transformation.
  Expect GZ with eta = eta_D, vartheta = +-lambda theta_d -> lambda theta_d = (lambda+1) theta_d^cl (cnj-defect-parameter).
  Try the scalar identity too (Phi, R_0 sigma of boundary_sine_gordon.py).
- BDG (e): defect potential = H_osc with p = beta psi/2 - i theta_d, q = beta psi/2 + i theta_d (book normalizations).
  Dirichlet: frozen oscillator at complex point -> static GZ.  Neumann (Bajnok-Simon (5)): one boundary variable psi,
  conjugate to phi(0); boundary condition has extra psi_t term; BDG flow does not preserve p - q = const.
  So the fused system is NOT the BDG boundary.
- Sources: Bao-Wang Rem 4.10(2) is a finite-type remark ("one checks"); Wang ICM Ex. 0.1 and Sec. 2 "quasi K for
  diagonal type reduces to quasi R"; Balagovic-Kolb do not state diagonal type (they identify X = Upsilon);
  Kolb 2014 assumes indecomposable A (excludes g^ + g^).  Kolb-Yakimov "diagonal type" = diagonal braiding (different!).

## FINDINGS (with scratch scripts in ch42/)
(a) s1_doubled.py: doubled Rc (Rc on VV, P Rc P on V'V', flip mixed) intertwines Delta_A of U(x)U^cop to 8e-16 (sG, a2);
    with U(x)U (Delta on 2nd factor) fails 1.7.  V(1/y) o omega = Vbar(y) exactly (a2: 2e-16; sG: Vbar = C V C);
    P Rc P = Rc_{Vbar Vbar} up to scalar (4e-16, 6e-16).
(b) dd_engine.py + s7_coideal.py + s8_twisted.py: search space X = span{g(x)h}, g,h in {k^a, e, f} (deg<=1, 21 elts per node);
    one-particle commutant of the doubled K: 17-dim; adding two-particle doubled RE: 7 (spin-1 probe) or 12 (spin-1/2, a2 vector)
    -- spurious symmetric Cartan combos (artifacts of <=2 particles).  The MAXIMAL LEFT COIDEAL inside the commutant
    (iterate S -> {x in S: Delta_A(x) in A (x) S}) is 5-dim = span{1, Delta(k^+-1), Delta(e), Delta(f)} = Delta(U) cap X
    exactly (sG spin 1/2 and spin 1 probes, a2 vector; all nodes).  Analogue of Appel-Vlaar Thm 4.1.3 (B_xi maximal coideal).
    Invertible lines: twist e^{i alpha mu} -> (1(x)t)Delta(U), t_1 = e^{2i alpha}, t_0 = e^{-2i alpha}; charge conjugation C ->
    (1(x)tau)Delta(U), tau = node swap.  All exact (1e-15).
(c) s9_jump.py, s11_gen_test.py, s12_limits.py: type I_- (Weyl, b^-), window N=5,7,9, interior margin 2:
    - Delta(b^-) (f's, k's) with W = Weyl module: holds (6e-16).  [I_+: Delta(b^+) holds 1e-15]
    - least squares over theta-independent weight-graded W-operators for generators of (1(x)psi)Delta(U):
      I_-: e/id 0.95-0.96 FAIL; e/tau 1e-15 PASS; f/id PASS; f/tau 0.96 FAIL.   I_+: e/id PASS; e/tau 0.93 FAIL; f/id 0.97 FAIL; f/tau PASS.
      => no twisted diagonal (any torus x {id,tau}) is a symmetry: cnj-twisted-diagonal REFUTED as stated.
      => NEW: I_- also has the tau-twisted OPPOSITE Borel (1(x)tau)Delta(b^+), with a second b^+-module structure rho~ on W
         (by hand, sg matrix normalization A=q^-n, s): rho~(e_1)|n> = a q^{-2n}|n+1>, a = i y_d q^{-s}/(q^2-1);
         rho~(k_1)|n> = q^{-2n-2-s}|n+2>; rho~(e_0)|n> = b q^{2n}|n-1>, b = i y_d q^s/(q^2-1); rho~(k_0)|n> = q^{2n-2+s}|n-2>.
         k_0 k_1 = 1, k e k^-1 ok, e_1 e_0 = e_0 e_1 = ab q^2 central (q-Serre).  Free constant c in Y = c + a q^-2n is
         a homogeneous solution; c = 0 forced by k e k^-1 = q^2 e.  Interpretation: = b^+ symmetry of the I_+ defect obtained by
         psi -> -psi (prp-two-borels (3)): psi -> -psi is C on the reflected copy = tau.
      - Combination impossible: rho(k_1) (diag q^{2n+s}) and rho~(k_1)-type kappa (shift by 2) q^4-commute, while k_1(x)k_1 and
        k_1(x)k_0 commute in A.  Abstractly: any coideal containing Delta(k_1) and e_1(x)1 + k_1(x)e_0 contains e_1(x)1
        (q^4 != 1), and e_1(x)1 is not a symmetry (y^2 C_n term) -> no module of any kind (incl. BCZ two-label space).
    - Limits (s12_limits.py, N=6): theta_d -> +inf: e/id residual 0.44, 0.063, 0.008 at theta_d = 1.5, 3, 4.5 (-> twist line,
      torus-twisted diagonal restored, rate ~ e^{-lambda theta_d}); theta_d -> -inf: f/tau 0.27, 0.036, 0.0046 (-> charge
      conjugation line, tau-twisted diagonal restored).  e/tau and f/id stay exact throughout.
(d) s13_fusion.py, s14_fusion_re.py, s15_scalar.py: K_f = Tbar (1 (x) K_D) T, T = I_- matrix:
    - gauge invariants of the 2x2 charge-sector blocks n-independent (8e-16) [cz2012 (2.18) alpha-independence];
    - equal to GZ eq-bsg-k-matrix with (eta, vartheta) = (eta_D, lambda theta_d) exactly (1e-15); alternatives (-lam th_d, eta->-eta,
      vartheta = theta_d) fail 0.02-0.3;
    - RE on windows N=5,7,9: 6e-16; perturbed control 0.03-0.07;
    - FULL scalar identity: Tbar K_D T (with Phi of eq-sgd-scalar, K_D = R_0 sigma(eta,u)/(2 cos eta) diag(ye^{i eta}+..)) equals
      BPTT K_full(theta; eta_D, lambda theta_d) entry by entry: ratio 1 to 5e-16 at three theta.  (sigma(i vt,u) -> 1 as vt->inf.)
    - With the classical limit of cnj-bsg-uv-ir: vartheta = (lambda+1) theta_d^cl => lambda theta_d = (lambda+1) theta_d^cl
      = cnj-defect-parameter (ch. 37), now from the soliton sector.
(e) BDG by hand: B_d(phi, psi) = -(2m0/b^2)[e^{i b phi/2} cos(b psi/2 + i theta_d) + e^{-i b phi/2} cos(b psi/2 - i theta_d)]
    = H_osc with q^ = b psi/2 + i theta_d, p^ = b psi/2 - i theta_d (b = beta_SG).  Dirichlet psi = phi_D: frozen oscillator at a
    complex point = static GZ with C_{0,1} = cosh(theta_d +- i b phi_D/2) (ch. 37).  Neumann: one boundary variable psi(t),
    conjugate to phi(0); BC phi_x = psi_t - d_phi B_d has the extra psi_t; BDG flow does not preserve p^ - q^ = -2 i theta_d.
    Classical dressing: K_phi(lam) = sigma(L_d^{-1}(1/lam)) K_psi(lam) L_d^{-1}(lam) satisfies the dynamical compatibility
    d_t K = K a_t(lam) + sigma(a_t(1/lam)) K if K_psi does (two-line proof). To check numerically in script.

## Dead ends
- Commutant of K inside a monomial search space with only spin-1/2 probes and <=2 particles: many spurious solutions
  (ek^a degenerate on spin 1/2; symmetric Cartan combos).  Fixed by the maximal-coideal reduction.
- boundary_classical, boundary_sine_gordon, qsp_checks, boundary_coideal, defect_axioms run their checks at import: do not
  import them; only boundary_common and defect_common are import-safe.  Reimplement small helpers (R_0, sigma) in the script.
