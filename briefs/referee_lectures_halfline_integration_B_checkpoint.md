# Integration agent B checkpoint (Lectures 3 and 4)

Items: L3-1..L3-10, L4-1..L4-16, G1-G4 in `lecture-notes/half-line/03-*.qmd`, `04-*.qmd`, `code/lectures/halfline_l3_f1_mirror_images.py`, `code/lectures/halfline_l4_f1_coideal_sg.py` and their SVGs.
Scratch: /tmp/claude-1000/-home-gustav-Git-affine-toda/af9e27ec-122e-43f3-a547-2dbab370a65e/scratchpad/integ_B/

## Log

- (start 2026-10-25) Reading reports and lectures.
- Read P1, P2, combined checklist, L3, L4, both scripts, book ch. 29-33 sections. Plan: scratch checks integ_B_l3_checks.py and integ_B_l4_checks.py, then edit L3, L4, scripts, rerun. Next: write scratch L3 checks.
- Scratch integ_B_l3_checks.py: ALL PASS. Verified: static-pair pattern C_j=sgn(c_{j+1}c_{j-1}) with integer indices (h=3..8, all a), unique, equals book (-1)^{alpha_j.lambda} with lambda from continuous log, never uniform, product +1; a_2 table (-pi,-pi/3):(-,+,-), (-pi/3,pi/3):(+,-,-), (pi/3,pi):(-,-,+); chi=0.25 -> (+,-,-); E_vac a_2 and a_3; inertia beta^2 m_1 I = 6; lower-bound positivity (old root claim false); corrected derivative 2(1-Cs), k s = 2m sin^2, ratio k-mC, branch = integer-index product; K_cl=[-(a)(h-a)]^{-C}; a=h/2,C=-1 pole at i pi/2 is DOUBLE in theta (strip edge); breather C=+1 wall zero, C=-1 interior zeros. No lecture edits yet. Next: L4 scratch checks.
- (resumed after usage stop) Scratch integ_B_l4_checks.py: ALL PASS. Verified: eps=0 nonzero K iff c0c1=q^-2 for generic y (y^4=1 exceptional), c0=c1=-1/q also allowed; braid R solved from coproduct: R(-q)=(Pi x 1)R(q)(Pi x 1); K(eps;-q)=K(-eps;q); RE at q and -q; matrix crossing k(theta) ~ R(2theta)k(-theta): q plain C holds, -q plain fails (0.773), -q twisted C holds; y_n=i(-1)^n e^{-i eta}, v_n=(1,(-1)^n e^{-vartheta}), beta_{1,0} correct; common eigenvectors exactly at y^2 in {-e^{+-2i eta},-e^{+-2 vartheta}}; det(b_1-beta k_1) y-independent; P^pm swap; K/(-4i)=[[P+,Q],[Q,P-]]; Neumann eps_* and (1-q)(1-1/q)=4cos^2; UV-IR combination; Neumann poles 0<n<lam/2, E=m_{2n}/2; contact identity sign. Also by hand: SG current weight (8pi/beta)^2/(8pi)=lambda+1; Kolb (1207.6036 eq. 3.3-3.4, Prop. 5.2) has the SAME coproduct and right coideals B_i=F_i+c_i theta_q(F_iK_i)K_i^-1+s_iK_i^-1; omega: e<->f, k<->k^-1 is an algebra automorphism reversing the coproduct and maps b_j to Kolb's form. dVGR abstract (six/eight-vertex K-matrices) and MN abstract (free-fermion IM) read. Next: edit L3.
- L3 edits done in 03-mirror-images-and-vacua.qmd (not yet logged per item; items below). Next: L3 script rewrite (interval table check), then L4 edits.
- L3-1 verified (wording own): intro now says L2 treated a_1,a_2, links prp-durham; mirror enters from x>0, collision behind/at/in front; antisoliton return scoped (parity + book's numerical check species 1 a_2,a_3); real-coupling vacua "believed, not proved"; chi "not fixed by the shift", regular except isolated chi; solitonic = product +1 and not all +1, construction misses some patterns (link sec-boundary-missing-charges). File: 03. Evidence: integ_B_l3_checks.py + book ch. 32 read.
- L3-2 verified: k s = 2m sin^2(pi a/h); derivative 2(1-Cs); constant-ratio step k-mC written out in text with branch (vacuum v_a, integer indices) and in exr-lec3-mirror-square (new part 4 on index reduction). File: 03. Evidence: integ_B_l3_checks.py (branch from unwrapped log = 4c_{j+1}c_{j-1}/(1+Cs), h=3,4,5).
- L3-3 verified: zeta^{-ja} in dominant-term factor. File: 03. By hand.
- L3-5 verified: positivity via lower bound ((1-u)^2+s(1-u^2))/(1+s); old root claim removed (roots at cos=1 are negative). File: 03. integ_B_l3_checks.py.
- L3-6 modified: rotor sentence now: inertia I=6/(beta^2 m_1) (checked by quadrature), diagonal tower b_{n,-n} rotor-like at weak coupling, global rotor through three intervals is an interpretation, fixed BC = one interval of length 2pi/3. Matches coordinator's wording for agent C. File: 03. integ_B_l3_checks.py.
- L3-7 verified: [-(a)(h-a)]^{-C}, block (x) defined, K=1 at C=0 stated. File: 03. integ_B_l3_checks.py. Book note: eq-linear-reflection prints -[(a)(h-a)]^{-C} (wrong at C=0, though next sentence says K=1) -> pending.
- L3-8 verified: exr-lec3-mirror-square heading now links book:prp-images (C=+-1 mirror shift), prefix removed. File: 03.
- L3-9 verified: Exercise 3.2 part 2 adds wall zeros (singular on wall, regular only on open half-line); part 3 precise cases (s<1 factorization + rho!=0 negative at wall; s>1 tau(0)<0; s=1 separate). Title no longer claims binding. File: 03. integ_B_l3_checks.py + book prp-boundary-breathers proof.
- L3-10 modified: Exercise 3.3 pole at theta=i pi/2 on the strip edge (stated); energy argument via eq-lec3-de-dchi hint. Found that the pole is double in theta (simple in p); not stated in lecture, only "edge". File: 03. integ_B_l3_checks.py.
- G3 (L3): all typed numbers in link text and exercise prefixes removed. G1 (L3): "exact"/"decisive benchmark" -> "proposed amplitudes ... compared". Also removed "!" inside three displays (read as factorials) and "Wait". File: 03.
- L3-4 verified (recomputed from the BC on the real field, h=3..8; agrees with book prp-solitonic-vacua C_j=(-1)^{alpha_j.lambda}, lambda from continuous log): text now derives C_j=eps_{j+1}eps_{j-1} with integer indices, eps_{k+h}=(-1)^a eps_k, product argument with wraparound, corrected a_2 table (-pi,-pi/3):(-,+,-); (-pi/3,pi/3):(+,-,-); (pi/3,pi):(-,-,+); caption (+1,-1,-1) at chi=0.25. Script halfline_l3_f1_mirror_images.py rewritten: tests the printed table interval by interval (3 chi per interval), all 7 other patterns fail, formula at integer indices reproduces table, positivity, EOM, energy; new check of the C=+-1 blackboard steps; K_cl block form incl. C=0. Legend fixed (C_0=+1,C_1=-1,C_2=-1), SVG regenerated (real change, kept). 7 PASS, ~1.6 s. Files: 03, script, l3 SVG. Evidence: integ_B_l3_checks.py.
- G2 (L3): L3:127 claim now "at sampled parameters ... to 1e-13 with v_a computed from the solution"; caption says what the script checks. G4 (L3): legend/caption fixed.
- L4 text edits applied in 04-coideals-and-sine-gordon.qmd (items logged below after script update). Next: L4 script (add c0c1 rank check, twisted-C crossing at -q, explicit v_n and branch set, P/Q parametrization), rerun, caption already fixed; then check-refs, copy scratch scripts, final table.
- L4 script updated: 7 checks, all PASS, ~1 s (coideal coproduct; intertwiner + RE; c0c1=q^-2 rank test incl. c0=c1=-1/q; R(-q)=(Pi x1)R(q)(Pi x1) and K(eps;-q)=K(-eps;q); matrix crossing plain/wrong-sign/twisted; v_n eigenvector + resultant branch set; (eta,vartheta) form + Neumann point/poles). Figure code unchanged; SVG diff was ids/dates only -> git restore.
- L4-1 verified: callout (matrix part up to scalar; first-level states), L4:253 scalar has CDD freedom, UV-IR "Assuming the conjectured relation", zero-coideal point "approaches" uniform, links cnj-neumann-coideal; L4:310 first-level states, higher levels need bootstrap. Neumann point also attributed to GZ free boundary (read GZ 9306002.tex lines 1645-1662: eta=pi(lam+1)/2, vartheta=0 from C-symmetry and expected pole). File: 04. Evidence: book ch.30/33, GZ source, integ_B_l4_checks.py.
- L4-2 verified: right-moving excitation hits wall, returns left-moving. File: 04. By hand.
- L4-3 verified: bar = parity transform, not complex conjugation; scalar case phi(x)+phi(-x), Neumann derivative zero. File: 04. By hand + book eq-parity-image.
- L4-4 verified: "spin-1" removed; weight coupling-dependent (lambda+1 in SG: (8pi/beta)^2/(8pi)); pole order coupling-independent (exponent product). File: 04. By hand.
- L4-5 verified: first order, +O(lambda_b^2), eps_j ~ lambda_b eps_j, "suggests"/"motivates, not proves exact conservation"; simply-laced restriction stated. DM source read (0112023 source.tex ~800-850: exponent e^{i(1/b-b)phi} = q^T in their normalization, "conserved to first order"). File: 04.
- L4-6 verified: MN constructed boundary IM at free-fermion point (abstract read, 9709078); DM developed coideal method. File: 04.
- L4-7 verified: normalization hypothesis spelled out (scalar endomorphisms -> proportional at each rapidity pair; factor may depend on rapidities; needs compatible normalization agreeing on a fixed nonzero matrix element for every generic pair; highest-weight normalization not automatic for mu != nu). File: 04. By hand. Book G5c remains pending (coordinator).
- L4-8 verified: dVGR general six/eight-vertex K-matrices (abstract read 9306089), GZ physical amplitude; "eight linear equations" now derive the matrix part. File: 04.
- L4-9 verified: "physical parameters (eta,vartheta) of BPTT, written in the book's orientation and node conventions". File: 04. K/(-4i)=[[P+,Q],[Q,P-]] checked.
- L4-10 verified: "if and only if" -> generic, matrix part; "shown in @fig-lec4-coideal" deleted; twisted C restores at wrong sign (new script check); script scope (proportionality, no scalar) stated. File: 04, L4 script.
- L4-11 modified: replaced wrong component display by the book's vector form k(theta)=S(2theta)k(-theta), k=sum K_{abar}^b(i pi/2-theta) v_a x v_b, linked to Lecture 1's eq-lec1-crossing-unitarity (agent A's current wording uses the same K_{abar}^b convention). Not the referee's reflection-rapidity component form, to avoid a second index convention. File: 04. Script check tests exactly this vector form.
- L4-12 verified: determinant detour removed; quadratics for alpha, v_n=(1,(-1)^n e^{-vartheta}), all four branches y^2 in {-e^{+-2i eta},-e^{+-2 vartheta}}, scalar poles select y^2=-e^{-2i eta}. File: 04 + exr-lec4-1d-submodule part 4. integ_B_l4_checks.py (sympy resultant).
- L4-13 verified: caption no longer claims crossing residuals or energies; UV-IR labels flagged. File: 04.
- L4-14 verified: generic y (y^4!=1 exceptional), product fixed, parity allows c0=c1=+-q^-1, choose q^-1 (also from charge normalization, book eq-coideal-generators). Title no longer says "forced". File: 04. Script rank check.
- L4-15 verified: GZ free point / UV-IR at mu=0; analytic branch 2cos(pi lam/2); 0<n<lam/2 interior = B_{2n} existence range; n=0 strip edge zero energy not a bound state; mirror argument heuristic only. File: 04. integ_B_l4_checks.py + script check 5.
- L4-16 verified (optional item done): parenthesis in callout: Kolb uses same coproduct (1207.6036 (3.3)-(3.4)), generators F_i + c_i theta_q(F_iK_i)K_i^-1 + s_iK_i^-1, right coideals (Prop. 5.2); omega: e<->f, k<->k^-1 reverses coproduct (checked by hand). File: 04. Note: P2 said "opposite coproduct"; Kolb's coproduct is the same, the difference is the generators -> wording modified accordingly.
- Node-label paragraph (ed0585a): agrees with current book ch.30 "Node labels"; only link texts changed (G3).
- G1/G2/G3 (L4): all typed numbers in link text removed; "Hold on!", "wait:" removed; script claims at L4:263, 312, 353 now state what the script tests.
- Final: both lecture scripts print only PASS (7 + 7); check-refs.py lecture-notes/half-line PASS; check-refs.py PASS except pre-existing WARN citation.qmd:11 @sec-about (not my file). Scratch scripts copied to briefs/referee_lectures_halfline_scripts/integration_B_l3_checks.py and integration_B_l4_checks.py (both ALL PASS). L4 SVG restored (ids/dates only); L3 SVG kept (legend fix).

## Table

| Item | Verdict | How checked | Files changed |
|---|---|---|---|
| L3-1 | verified (own wording) | book ch. 32/31 read; integration_B_l3_checks.py | 03 |
| L3-2 | verified | by hand + integration_B_l3_checks.py (derivative, ks, ratio, branch vs v_a) | 03, L3 script |
| L3-3 | verified | by hand | 03 |
| L3-4 | verified | BC on real field recomputed, h=3..8; equals book (-1)^{alpha.lambda}; a_2 table; script now tests printed table | 03, L3 script, l3 SVG |
| L3-5 | verified | lower bound by hand; root counterexample numerically | 03 |
| L3-6 | modified (coordinator's wording; inertia 6/(beta^2 m_1) recomputed) | quadrature in integration_B_l3_checks.py | 03 |
| L3-7 | verified | integration_B_l3_checks.py (blocks, C=-1,0,1); new script check | 03, L3 script |
| L3-8 | verified | link target read (prp-images) | 03 |
| L3-9 | verified | factorization/sign cases numerically; book proof read | 03 |
| L3-10 | modified (strip edge stated; energy via variational hint; double pole in theta noted only here) | numerically | 03 |
| L4-1 | verified | book ch. 30/33 labels; GZ source for Neumann point | 04 |
| L4-2 | verified | by hand | 04 |
| L4-3 | verified | by hand, book eq-parity-image | 04 |
| L4-4 | verified | by hand (weight lambda+1, pole order) | 04 |
| L4-5 | verified | DM source (first order, q^T exponent) | 04 |
| L4-6 | verified | MN abstract | 04 |
| L4-7 | verified (own wording) | by hand (rescaling argument) | 04 |
| L4-8 | verified | dVGR abstract | 04 |
| L4-9 | verified | K/(-4i) form checked | 04, L4 script |
| L4-10 | verified | integration_B_l4_checks.py (R solved from coproduct); new script check (twisted C) | 04, L4 script |
| L4-11 | modified (book's vector form instead of component form) | script tests that form | 04 |
| L4-12 | verified | sympy resultant, numerics | 04, L4 script |
| L4-13 | verified | SVG content | 04 |
| L4-14 | verified | rank test (also y^4=1) | 04, L4 script |
| L4-15 | verified | by hand + script check 5 | 04, L4 script |
| L4-16 | modified (Kolb has the same coproduct; difference is the generators) | Kolb 1207.6036 (3.3)-(3.4), Prop. 5.2; omega by hand | 04 |
| G1 | verified in 03, 04 | book labels | 03, 04 |
| G2 | verified in 03, 04 | script contents | 03, 04 |
| G3 | verified in 03, 04 | grep for typed numbers; check-refs | 03, 04 |
| G4 | verified | l3 legend fixed; l4 caption fixed, figure unchanged | L3 script, l3 SVG, 04 |

Counts (26 numbered items): 22 verified, 4 modified (L3-6, L3-10, L4-11, L4-16), 0 rejected. L3: 8 verified, 2 modified; L4: 14 verified, 2 modified. G1-G4 applied in both files.

## Pending

Book findings (for the coordinator's book PR; not edited):
- `part7-half-line/32-half-line-solitons.qmd` eq-linear-reflection (~228): block form printed as -[(a)(h-a)]^{-C}; correct for C=+-1 but gives -1 at C=0, while the next sentence says K=1 for Neumann. Suggest [-(a)(h-a)]^{-C}, which holds for C in {0,+-1}.
- G5c (coordinator's): thm-k-reflection hypothesis (iii) "agree at one point" should mean a fixed nonzero matrix element for every generic rapidity pair; the lecture now states this.
- Minor: eq-sg-neumann-eps writes 1/sqrt((1-q)(1-q^{-1})) without the branch; the value is the analytic branch 2cos(pi lambda/2) (may be negative).
- Minor, informational: for a=h/2, C=-1 the classical factor (sinh theta + i)/(sinh theta - i) has a double pole in theta at i pi/2 (simple in p); the book's "pole sits at theta = i pi/2" is correct as stated.
- For agent C / coordinator (L5-6): P2 says Kolb's coideal convention needs the "opposite coproduct". Kolb 1207.6036 (3.3)-(3.4) uses the same coproduct as the book; his B_i = F_i + c_i theta_q(F_iK_i)K_i^-1 + s_iK_i^-1 are right coideals because of the generator choice; e<->f, k<->k^-1 maps one to the other.

Not finished: nothing.

Scripts (for code/lectures/README.md):
- `code/lectures/halfline_l3_f1_mirror_images.py`: 7 PASS checks, ~1.6 s: (1) mirror pair BC at imaginary coupling for C=0,+-1 relative to v_a; (2) C=+-1 blackboard steps (derivative with 2(1-Cs), ks=2m sin^2(pi a/h), ratio k-mC with integer-index square root); (3) Delta t_N average and virtual collision d; (4) stationary vacuum positivity and static EOM; (5) BC with the lecture's interval table (-pi,-pi/3):(-,+,-), (-pi/3,pi/3):(+,-,-), (pi/3,pi):(-,-,+), all other patterns fail, eq-lec3-cj-from-chi at integer indices reproduces it; (6) E_vac = 2(4-3 sqrt3) in all intervals; (7) linearized K_cl and [-(a)(h-a)]^{-C}, C=+1,0,-1. Figure legend now (C_0,C_1,C_2)=(+1,-1,-1) at chi=0.25.
- `code/lectures/halfline_l4_f1_coideal_sg.py`: 7 PASS checks, ~1 s: (1) left coideal coproduct; (2) SG intertwiner and braid reflection equation; (3) eps=0 rank test c0c1=q^-2 (and c0=c1=-1/q); (4) R(-q)=(Pi x 1)R(q)(Pi x 1), K(eps;-q)=K(-eps;q); (5) matrix part of boundary crossing-unitarity: plain C at physical q, fails at -q, twisted C restores it; (6) common eigenvector v_n at y_n with beta_{1,0}, and all branches y^2 in {-e^{+-2i eta}, -e^{+-2 vartheta}}; (7) (eta,vartheta) form of K, Neumann eps_*, Neumann poles and energies, interior range 0<n<lambda/2.
