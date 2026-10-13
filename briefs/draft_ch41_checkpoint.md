# Chapter 41 drafting checkpoint (2026-10-12)

State: exploration done in the scratchpad; NOTHING written yet to `part8-defects/41-soliton-transmission.qmd`,
`textbook_code/defect_an_transmission.py` or `textbook_code/defect_common.py`. Scratchpad (session-specific, may vanish):
`SP=/tmp/claude-1000/-home-gustav-Git-affine-toda--claude-worktrees-integrable-defects-framework-a54a84/cfddf2dd-7e99-4c1c-b8f5-c799ea33aced/scratchpad`.
Key helper draft: `$SP/arhelp.py` (an_root_window, an_weyl_module, solve_intertwiner, T_vector, block_inverse,
crossing_residual, crossing_residual_bar, tensor_rep, fusion_map, solve_intertwiner2). a_2^(2): `$SP/a22.py`, `a22d.py`, `a22e.py`, `a22j.py`.
All tests used LAM = omega = gamma = 1.37, q = q_phys(LAM) = -exp(-i pi LAM), y = e^{LAM theta}.

## 1. Findings (all numerical, residuals ~1e-14 unless stated)

**a_r Borel (Weyl) module, r >= 2** (book conventions, Part VI coproduct). Labels N in Z^h, sum N = 0 (root lattice) + offset s
(class p: s = lambda_p = e_1+..+e_p - p/h). alpha_j = e_j - e_{j+1} (j=1..r), alpha_0 = e_h - e_1 (matches boundary_common.an_vector).
  k_j|N> = q^{alpha_j.(N+s)}|N>;  b^-: f_j|N> = (mu/y_d) q^{(N+s)_{j+1}} |N - alpha_j>  (component j+1 mod h, 1-based; j=0 -> comp 1);
  b^+: e_j|N> = mu y_d q^{(N+s)_{j+1}} |N + alpha_j>.  Then f_j f_{j+1} = q^{-1} f_{j+1} f_j, others commute => q-Serre holds.
  Control: plain shifts (no q-power) violate q-Serre (0.55). Script: `$SP/explore1.py`.
  Normalization chosen: mu = -e^{-i pi omega/2}/(q - q^{-1}), so X := -(q-q^-1) mu y/y_d = e^{omega(theta - theta_d - i pi/2)} = CZ's xhat
  with Delta = theta_d + i pi/2 (cz2009a (6.12)-(6.13)). Equivalent: xhat = (-q)^{1/2} y/y_d with (-q)^{1/2} := e^{-i pi omega/2}.
**Orientation**: c=+1 (the q-power above) <-> type I_- with the book's pi (b^-), and type I_+ (b^+); c=-1 gives pi^{-1} (cz2007's convention).
**Intertwiners V (x) W -> W (x) V, b^- (type I_-)**: unique up to scalar (a_2: deep window rank 1; a_3: rank 1, ~35 s), closed forms
(universal-R normalization, = CZ 2009a (5.8) and (5.11)/(6.5) in the SAME gauge, Q = 1/q):
  V:    T(v_i (x)|N>) = q^{-(N+s)_i}|N>(x)v_i + X |N + e_i - e_{i-1}>(x)v_{i-1}      [only two kinds of entries: quantum selection rule, moves by pi]
  Vbar: T(vb_i(x)|N>) = sum_{j=0}^{h-1} X^j q^{(N+s).(e_i+...+e_{i+j})} |N - e_i + e_{i+j}>(x)vb_{i+j}   [all h entries; tunnelling]
  Scripts `$SP/explore3.py`, `explore6.py`, `explore8.py` (closed-form residual 2e-15..7e-10 on deep labels).
**b^+ (type I_+)**: unique; Vbar is the linear one, moving vb_i -> vb_{i-1} (= pi on weights -e_i): A_i = q^{-(N+s)_i},
  C_i = kappa mu (y_d/y) q^{-1 + (N+s)_i - (N+s)_{i-1}} (fitted for s=0, a_2; verify for s != 0); V is full. `$SP/explore4.py`, `explore14/15.py`.
  Rank needs tol 1e-6 (Gram-matrix conditioning; tol 1e-8 gave a false rank 0 for generic s).
**Hypothesis (ii) of thm-t-triangle**: Hom_{b^-}(V(y1)(x)V(y2)(x)W, W(x)V(y2)(x)V(y1)) one-dimensional for a_2, R2=14, generic s (25 s). `$SP/explore22.py`.
**Triangle relation** with boundary_common.rcheck(an_vector): NOT yet run for a_r (expected to pass: thm-t-triangle).
**Conjugation C**: dual(V(y)) = Vbar(-y/q) with C = identity (v_i^* -> vb_i), both signs of q; for h even also Vbar(y/q) with C = diag((-1)^i).
  `$SP/explore5.py`. Physical shift: y(theta - i pi) = -q_phys y.
**Crossing (I)** (eq-defect-crossing for V): M_V(theta) = c_I Cross[M_Vbar(theta - i pi)^{-1}], c_I = 1 + (-q)^r X^h, s-independent (s=0, lambda_1).
**Crossing (II)** (for Vbar): M_Vbar(theta) = c_II Cross[M_V(theta - i pi)^{-1}], c_II = 1 - q X^h (= CZ (5.11) factor 1 - xhat^h Q^{-1}).
  r = 2, 3, 4. `$SP/explore9.py`, `explore10.py`, `explore16.py`.
**q-sign**: at -q_phys crossing FAILS for h odd (a_2, a_4) with every C (crossing point x_c = (-1/q)^h flips sign; dual(V) only = Vbar(y/q_phys));
  for h even (a_3) it fails with C = 1 and HOLDS with C' = diag((-1)^i) (graded conjugation, analogue of ch.37's [[0,1],[-1,0]]).
**Bootstrap**: fusing r V's at theta + i pi(2m - r + 1)/h, m = 0..r-1 (left particle lowest) into Vbar = Lambda^r V gives exactly
  M_Vbar (ratio 1) for a_2 and a_3 => gbar(theta) = prod_m g(theta + i pi(2m-r+1)/h)  (= cz2009a (6.6)-(6.7)). `$SP/explore17.py`, `explore18.py`.
**Scalar factor** (book form, convergent): g = c_r x^{-1/2} G(z), z = i h gamma (theta - Delta)/(2 pi), x = xhat,
  G = prod_{k>=0} Gam(1/2+r g/2+z+khg) Gam(1/2+r g/2-z+khg) / [Gam(1/2+(r/2+1)g+z+khg) Gam(1/2+(r/2+1)g-z+khg)], each factor divided by its z=0 value
  (g = gamma). Satisfies cz2009a (6.8) g(theta+2 pi i)/g = (1 + (-Q)^r X^h)/(1 + (-Q)^{r+2} X^h), (II) and (I) to 1e-9 with c_r real > 0 fixed by (II)
  (r = 2,3,4; gamma = 1.37, 0.6). `$SP/explore19.py`. CZ's (6.10) is the SAME product re-paired (one numerator shifted by one k):
  as printed its partial products diverge like K^{-z(r+1)/h} (theta-dependent phase drift) -> correction, analogous to BCZ (3.47).
**Poles / unstable states**: g pole at theta_1 = Delta - i(pi r/h + pi/(h gamma)); there det of each h x h charge block of M_V vanishes
  (det = q^r [1 + (-Q)^r X^h]) => residue rank h-1 per sector (not 1!). gbar pole at Delta - i(pi/h + pi/(h gamma)).
  a_2 in BOOK labels: species 1: theta_d - i pi/6 - i pi/(3 gamma), unstable for gamma > 1 (includes classical limit, -> ch.39 absorption pole);
  species 2: theta_d + i pi/6 - i pi/(3 gamma), unstable window 1/2 < gamma <= 2 (= cz2007's window; cz2007 calls it species 1 because their D = D^T).
  General: window max(0, a - h/2) < 1/gamma < a (cz2009a text "1/gamma < h/2 - 1 or (h-1)/2" is a misprint for 1/gamma < a).
**Charged classes (kernel bootstrap)**: kernel of M_V(theta_1) in V(theta_1)(x)W is 1-dim per sector; T_23 Rc_12 restricted to it = class-1 module
  (offset lambda_1, A-ratios exactly q^{-Q_i}) with X^h -> q^2 X^h, i.e. theta_d -> theta_d + 2 pi i/h (mod 2 pi i/(h gamma)): quantum eq-ds-class-shift (I_-).
  a_2, `$SP/explore20.py`. Consequence (by hand): classical real absorptions a + 2p = h/2 mod 2h become resonances at theta_d - i pi/(h gamma).
**Semiclassics**: (2 pi/(h(gamma+1))) d arg g/dtheta -> ln|z_1(theta - theta_d)| (eq-dc-ar-delay), dev 0.092/0.025/0.0063 (a_2) and
  0.067/0.018/0.0046 (a_3) at gamma = 10/40/160; same for gbar (|z_r| = |z_1|). `$SP/explore21.py`. arg z: arg z_1(theta_d) = pi/h = half a
  charge sector, so a soliton at the centre of its sector changes charge iff theta > theta_d, the gamma -> inf limit of |X|^2/(1+|X|^2) (by hand).
**a_2^(2)** (CZ 2011 (4.1)-(4.10) algebra converted to Part VI coproduct by e_i = K_i X_i^+, f_i = X_i^- K_i^{-1}, k_i = K_i^2; q_i = q^{d_i},
  d = (4,1), Cartan [[2,-1],[-4,2]]; homogeneous x on node 0; generic q = e^{-0.77i}): 3-dim rep relations 2e-15; Borel module X1+ = a^dag,
  X0+ = aa, F(N) = f1 q^{-2N-1} + f2 q^{2N+1} passes (incl. 5th-order Serre); F = 1 and F = (-1)^N FAIL (no Weyl/type-I module, cf. prp-ar-only).
  Intertwiners unique for (f1,f2) = (1,.7),(1,0),(0,1),(.4+.3i,-1.1). Triangle with rcheck(vec3) holds 5e-15 (control 0.1).
  Crossing point: dual(V(x)) = V(x x_c), x_c = -q^{-6}; C antidiagonal [[0,0,1],[0,c,0],[c',0,0]]. Crossing holds EXACTLY with
  x(theta - i pi) = x/x_c and the TRANSPOSED placement Cm[d,b] = C[b,d] (a_r and SG insensitive since C symmetric). At -q: holds with C(-q)
  (= C with middle entry sign-flipped), fails (0.99) with C(q): q-sign seen only relative to the fixed conjugation (like SG, h even).
  Scripts `$SP/a22.py` (rep, module, solve_T22), `a22b.py` (triangle), `a22c.py` (x_c), `a22j.py` (crossing, both q).
  NOT done for a_2^(2): gauge invariants vs CZ (3.11)/(4.14); scalar factor; physical q/x dictionary; unstable windows (quote CZ (3.12)-(3.13)).

## 3. Conventions / decisions
- Book type I_- = cz2009a's B (2.38) = clockwise; cz2007 uses D^T: their species a = book h-a (already in tbl-defect-dictionary).
- Universal-R normalization = CZ's normalization (diag q^{-(N+s)_i}); its crossing scalars are s-independent -> one g for all classes.
- cz2007 (6.14) prints xhat = (-q)^{-1/2} e^{gamma(theta-eta)}; with Delta = eta + i pi/2 it should be (-q)^{+1/2} (minor; verify before stating).
- Gram-matrix null spaces: use tol 1e-6 and report the deep-label rank.

## 4. Dead ends
- Cube windows / 'deep' with ball R2 small: a_3 needs R2 = 14 (135 labels), a_4 uniqueness too expensive (R2=8 gives deep=1, 4275 unknowns) -> a_4 by closed form only.
- a_2^(2) crossing with C (untransposed) or x*x_c fails; first a22g run had an argument-order bug (ignore its output).
- Image-of-residue bootstrap for a_r: image is (h-1)-dim per sector, not a single Weyl module (only the kernel is).
- Robertson 1408.1960 not in ~/referee-copy/papers: fused a_3 matrix not checked.

## 5. Remaining, in order
1. Append helpers to `textbook_code/defect_common.py` (an_root_window, an_weyl_module, an_solve_intertwiner, an_T_vector, an_block_inverse,
   an_crossing_residual, an_scalar_g) -- append only.
2. Write `textbook_code/defect_an_transmission.py`: checks (a) module relations + control; (b) uniqueness r=2,3 + closed forms r=2,3,4 + b^+;
   (c) triangle with rcheck r=2,3 (both V, Vbar), hypothesis (ii) a_2; (d) crossing I/II, s-independence, q-sign; scalar factor eqs;
   bootstrap a_2 (a_3 if time); (e) poles, residue rank, windows, kernel bootstrap; (f) Wigner; (g) a_2^(2); (h) selection rule & tunnelling
   suppression |X|^j/max ratio -> 0. Keep runtime < 3 min (a_3 uniqueness ~35 s, hyp. (ii) ~25 s).
3. Chapter 41 (~300 lines, mirror ch. 33): learning goals; sec-t-procedure; sec-an-soliton-transmission (prp-an-transmission);
   sec-transmission-unitarity-crossing (prp-transmission-crossing, eq for g); sec-a22-transmission (prp-a22-transmission);
   sec-an-defect-bound-states (prp-a2-unstable); sec-transmission-semiclassics (prp-transmission-semiclassical); 7-8 exr-dt-*; Notes with corrections:
   cz2009a (6.10) pairing; cz2009a window misprint; cz2007 species labels (book a <-> h-a) for the window; possibly cz2007 (6.14).
4. Report `briefs/draft_ch41_report.md` and SubagentHandback.

## Progress (second agent)
- 2026-10-12: read checkpoint; helpers appended by agent 1 are complete (an_root_window .. an_scalar_constant, defect_crossing_residual); testing next.
- helpers fixed: an_solve_intertwiner vectorized and drops unknowns in no equation (then the FULL null space is 1-dim, gap ~0.02-0.05; a_3 R2=18 in 10 s);
  defect_crossing_residual vectorized (einsum); added an_tensor_rep, an_triangle_residual (sparse). Script sections (a)-(e) written and PASS (55 s).
  Findings: orient=-1 module gives V FULL and Vbar linear vb_i -> vb_{i+1} (pi^-1) = charge-conjugate defect (cz2007's D^T), not "V moves by pi^-1";
  crossing fit a = 1 exactly (c_I = 1 + (-q)^r X^h, c_II = 1 - q X^h), s = 0, lambda_1, generic; c_2,3,4 = 0.3983, 0.3933, 0.3910;
  bootstrap ratio 1 for a_2 AND a_3; kernel bootstrap PASS (sectors must be "in window", not "interior").
- script complete: 21 checks, all PASS, ~100 s (load ~2). Kernel bootstrap also from class 1 -> 2 (q^2 again). Added non-unitarity check (r >= 2).
  Next: write the chapter (ids eq-dt-weyl, eq-dt-tv, eq-dt-tvbar, eq-dt-tplus, eq-dt-crossing, eq-dt-g, eq-dt-bootstrap, eq-dt-pole, eq-dt-kernel), then fix the script docstring ids.
- chapter 41 drafted in full (all six sections, 8 exercises exr-dt-*, Notes); now reviewing/fixing text, then docstring ids and report.
- DONE 2026-10-12: chapter reviewed (cross-refs and bib keys verified by grep), script docstring updated (21 PASS, 97-132 s), report written to briefs/draft_ch41_report.md.
