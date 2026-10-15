# Codex phase-2 continuation checkpoint

Date: 2026-10-15. Worktree: `/home/gustav/Git/affine_toda/.claude/worktrees/referee-briefs-part-iii`; reviewed HEAD bcd82b9ee9667975eed02109f08a257d50ba38b1. Main checkout is a03edd4 and differs; do not switch reviewed version.

## Authorization and restrictions

User requested phase 2 per `briefs/referee_classical_09_11_items.md` (read in full). Report must be `briefs/referee_classical_09_11_codex_report.md`, scratch `referee_scratch/classical_09_11/codex/`. Every created filename includes codex. Do not edit any book file or git state, do not render Quarto, do not run author's scripts or spawn agents. Preserve frozen `briefs/referee_classical_09_11_codex_phase1_report.md`. DO NOT READ other referee reports `briefs/referee_classical_09_11_phase1_report.md` or `briefs/referee_classical_09_11_report.md` (user explicitly forbids main copies; we avoid all copies). Earlier Claude reports are allowed by brief but have not been read and are unnecessary. Before changing own existing files make a Codex copy, or write a new version. No need to reread all phase-1 context: own report contains full derivations and ranked 20-item list.

## Phase-2 report status

`briefs/referee_classical_09_11_codex_report.md` currently has header, frozen phase-1 hash and coordinator-list hash only; no check sections yet. Required final format: sections 1–9 following coordinator items, explicit overlap with phase-1 ranks, evidence/replacements, overlap table with one row per coordinator item AND per phase-1 entry, errors found only by either/both lists, unchecked claims/reasons, unverified observations. Phase 1 should be cited instead of repeated when fully covering a subitem. Need all 17 exercises in outline and a citation-by-citation audit (unavailable sources explicitly unclear).

## Independent new computations already created

All under `referee_scratch/classical_09_11/codex/`:

1. `codex_phase2_currents_v3.py`, output `codex_phase2_currents_v3_output.txt`, generated `codex_phase2_currents_v3.json`.
   Constructs A2 spin-2 and spin-4 local currents from scratch by solving polynomial conservation identities. Fields U=i beta phi in R3 sum-zero coordinates (u,v,-u-v), beta=m=1. D=partial_t-partial_x (gives spin +s in book rapidity convention). Exponentials z0=exp(u-v), z1=exp(u+2v), z2=exp(-2u-v); mixed equations F=(-z0+z2,z0-z1). Solves partial_+ T + D X=0; charge integral is T+X.
   Spin 2 representative:
   T=-u1^2 v1-u1 v1^2-u2 v1-2 v1 v2+v3; X=3 z0 v1.
   Spin 4 representative printed in output/JSON; T has 18 monomials, X=3 z0(u1 v2+v1^3-v1 v2).
   Current nullspaces dimensions 6 and 21, total-derivative ranks 5 and 20, one primitive per spin. Coefficients rational/exact. v1/u1 denote D derivatives, etc.
   Preliminary scripts v1 and v2 kept. v2 failed because old SymPy DomainMatrix.nullspace requires .to_field(); v3 fixed it. v1 is a slow plain SymPy solve, exec session 6827 may still be running; no PID available across per-command process namespaces. v3 finished successfully. Can poll v1 with write_stdin, but do not depend on it.

2. `codex_phase2_integrals.py`, output `codex_phase2_integrals_output.txt`.
   Integrates actual canonical energy, momentum and the above current densities on one solitons and breathers, analytic exponential jets, no asymptotic/additivity substitution. mpmath 32 digits, x range +-40.
   Tau roots alpha0=e3-e1, alpha1=e1-e2, alpha2=e2-e3, so U=(log tau0-log tau1,log tau1-log tau2,log tau2-log tau0). D log derivatives from cumulants of exact monomial derivatives. Checks all three canonical stress components against T00=2 sum logtau_xx, T01=2 sum logtau_xt, T11=2 sum logtau_tt.
   Numerical results:
   - species1 theta0 xi=.3i: E=10.392304845413263761165, P=0, Q2raw=+5.19615242270663188058 i=3sqrt3 i, Q4raw=-7.79422863405994782087 i=-(9sqrt3/2)i.
   - species2 theta0: q2=-q1 for both currents.
   - complex single theta=.3+.2i, xi=.4i, t=.1: E=10.646930242183003067385+.628722418343127100433i; P=3.101585063440180077815+2.158239607451275474235i. Canonical E,P match formulas to 8e-34, stress to 5e-32. Current integrals match q_a exp(s theta).
   - regular breather theta0=.23, u=pi/3, both xi=i*pi/2, t=.17: E=10.668395194383174381312; P=2.411359787244120294605; Q2raw=-14.25666586495033597 (~1e-14 imaginary quadrature error), Q4raw=-33.87542026414001856. Normalized Q2/q1=2.74369662495854 i, predicted 2i exp(.46)sin(2pi/3). Q4/q1=-4.346218446313987 i, predicted 2i exp(.92)sin(4pi/3).
   - shifted breather xi1=.2+i*pi/2,xi2=-.2+i*pi/2,t=.43 gives same charges (Q2 ~8e-15 quadrature error).
   IMPORTANT: current script reports FAIL for Q2 breather under tolerance 1e-16, error ~1e-14 whereas Q4 and E/P agree ~1e-29. Need diagnose/refine Q2 quadrature, not misreport passes. Likely lambdify constants / quadrature? Raw T has coefficients integers (v3 JSON), no obvious float except Python 2 division? Check script. Can increase mp.dps and use subintervals/quadgl or remove total derivative terms to improve cancellation. Do not simply declare additivity proved generally. Final test at reasonable verified numerical tolerance plus convergence comparison acceptable. User asked actual integration; done.
   Theta in U coordinates is U -> conjugate(U). Raw real-coefficient currents are Theta-even. Normalizing q1 to real (=1) divides by pure imaginary q1, hence normalized Q2,Q4 are Theta-odd. Not EVERY current is definite parity: (1+i)*raw is a concrete counterexample, or arbitrary mixtures. An antilinear involution allows a real basis; a chosen complex multiple need not have +/- parity. Phase independence also needs regularity/analytic continuation; additivity alone is insufficient. Phase1 already gives topological-charge counterexample to universal rule.

3. `codex_phase2_centralizers.py`, `codex_phase2_centralizers_output.txt`: finished all PASS.
   Independently builds defining matrix loop algebras and solves [E1,X_s]=0 plus twisting constraints by SVD. Numerical residuals <=4e-15, nonzero singular gaps >=.77, so robust finite-dimensional checks, not claimed exact general proof.
   A_(h-1) untwisted, h3,4,5: dimensions1 at all residues except0 modh.
   A4^(2): residues1,3,7,9 modulo10 (odd not divisible5).
   C2^(1): residues1,3 mod4.
   D3^(2)=A3^(2): residues1,3,5 mod6 (same set odd).
   G2^(1): residues1,5 mod6.
   D4^(3): residues1,5,7,11 mod12 (same set as G2).
   Construction: type A superdiagonal plus z E_(h1); automorphism sigma(X)=-J X^T J^-1, alternating anti-diagonal J. Grading j-i+Lk. D3 twisted E1=E12+E23+E34+z(E31+E42), L3. D4 triality built by closure of Chevalley brackets in so8 defining matrices, cyclic e1->e3->e4, e2 fixed. G2 fixed sector with affine negative highest root, L6. D4^3 uses z times eigenvalue-zeta vector of negative height3, L4, sigma(Xk)=zeta^k Xk.
   General A derivation: centralizer polynomials of regular cyclic E, E^h=zI; traceless excludes h multiples. A_2n twisted anti-transpose forces odd degrees. Can explain analytic support plus numerical requested dual pairs.

4. `codex_phase2_couplings.py`, `codex_phase2_couplings_output.txt`: launched session11748; need read output/poll.
   A4 Fourier charge basis tests all35 triples and area rule only nonzero couplings, then real eigenbasis counterexamples. E6 explicit Cartan Cholesky roots, highest marks (1,2,2,3,2,1), chain1-3-4-5-6 with2 attached4. Affine diagram centre rotation 0->1->6->0, 2->3->5->2. Node masses from Coxeter projection, characters from lambda1 dot lambda_a, simultaneous eigenvectors of mass and centre; tests all56 triples against root orbits and area. Real eigh basis counterexamples. Potential mapping issues must inspect, don't assume passes.

## Source files downloaded, read status

New PDFs and pdftotext files all `codex_source_<id>.pdf/.txt` in scratch:
- 1602.05465 (Cen/Fring, already substantial phase1 reading, full paper now downloaded).
- 1608.01691 (Cen/Correa/Fring, phase1 read reality section; full text now available).
- 2007.15425 (Fring/Taira, read intro and section2 pp1–4 so far).
- 2110.06825 (Correa/Fring/Taira, read intro/general reality argument, model and type I solitons pp1–7 so far).
- hep-th_9203076 (Hollowood): read section5 semiclassical pp26–29 and eq4.5/4.6, see sign issue below.
Need finish relevant source reading (user says read four papers) and search further work; targeted web search done (phase1 also searched Delius). Can read all in chunks without other reports.
Existing sources under `/home/gustav/referee-copy/papers/` (read only) include all phase1-listed arXiv papers. No need download those.

Key source findings:
- OTU hep-th/9305160 section8.3 pp26–27, local text lines1328–1355: additive light-cone momenta hold even for complex e^eta PROVIDED real part positive; explicitly independent of Q_i and sign epsilon_i. This directly supports phase1's missing-domain correction. Quote <=25 words/source in output reports if sourced via web (safer keep all direct excerpts short). Full exact stress-C formula not located; need locate original OTU1993a or Zhu/Caldi.
- Khastgir/Sasaki hep-th/9507001 eq3.10–3.12 (text255–300): tau=1+omega^(aj) exp(sigma x-lambda t+x0), sigma²-lambda²=m_a². Their c-solitons correspond theta complex, sigma=m_a cosh theta,lambda=m_a sinh theta. Our numerical complex single above constructs one. At tau zero solve two real linear equations; need calculate explicit first singular (x,t) and show regular initial slice. Read section5 for negative-energy examples too.
- H6 at `part9-foundations/44-twelve-gaps.qmd:33` asks physical-solution criterion, not only reality. Table line81 says closed. Phase1 already refutes sufficient-reality -> physical selection. Need concrete replacement in both chapter11 and H6 table.
- Fring/Taira 2007.15425 eq2.7 proves E=conjE by x->-x and density conjugation; eq2.8 also allows antilinear map between equal-energy BPS/anti-BPS partners. This already contains equal-energy rather than fixed-field logic. It does NOT justify general inequality E>=|Q| for arbitrary complex fields; avoid endorsing that source overclaim, describe precise sufficient identity.
- 2110.06825 section2 p2 says when symmetry maps field to itself OR second degenerate solution, energy real and reasoning extends to other conserved quantities. Model L=.5(dphi)^2-exp(phi)-.5exp(-2phi)+1.5 (real-coupling BD analytically complexified); maps to imaginary ATFT via U=i beta phi and normalization, energy sign differs under field scaling. TypeI solutions log[(cosh(beta+kx+sqrt(k²-3)t) +/-2)/(cosh(...) -/+1)], real energies E=-6|k| independent of complex beta. PT is spacetime reversal+conjugation, not book's field-only Theta. TypeII and stability sections still to read. Broad novelty clearly prior art; exact reduced species data phrasing may be presentation but no mathematical new theorem established.

## Important unfinished check: quantum scattering phase sign

Hollowood hep-th/9203076 uses n=h (SL(n)). Section5 eq5.8 shows log S leading derivative proportional to log[sin(theta/2i+pi/h) sin(theta/2i-pi/h)/sin²(theta/2i)] = -log A11. With lambda0=4pi, log S = i*(2h/beta²) integral log(1/A11), i.e. derivative of phase NEGATIVE of book's (2h/beta²) log A11. It also has time-delay conventions reversed in his eq5.3 (positive for gamma=logA<0). Need carefully map S orientation/rapidity/inverse to book, not announce error prematurely. Phase1 verified book SG S0 has derivative +16/bSG² log tanh and is correct in book convention. For h=2 Hollowood's formula may give inverse SG S0. Need compare explicit gamma product eq4.5/4.6 and book physical coproduct. May be physical S inverse/rapidity orientation. Need derive independently and document convention. Book chapter21 `sec-an-example` does not actually print scalar Gamma product; it refers to Hollowood. Full source eq3.12 etc available. General fused S_ab from product S11 shifted constituent rapidities could telescope to log A_ab; finite channel factors may only O(1) derivatives at generic real theta. Need calculate A2/A3 same and distinct species, with careful analytic branches and poles. If unable resolve must mark unclear with concrete obstruction, not invent match.

## Other remaining computations

- Explicit a2^(2) Lax 3x3 matrices. Use arbitrary root alpha1=a,alpha0=-2a, marks n=(1,2), H=a diag(1,0,-1), e1=sqrt2(E12+E23),f1=transpose, e0=E31,f0=E13. Then c1=a,c0=sqrt2 a. For reviewed finite-root norm2 set a=sqrt2; if affine-long norm2 set a=1/sqrt2. Abstract curvature should hold; state scaling explicitly to avoid book normalization conflict. Need symbolic check plus homogeneous gradation boost gauge (principal scale requires spectral power h and Cartan conjugation).
- Current spin3 one-particle value: phase1 proves displayed partial_+ current scales e^-3theta (not +). Linearized free scalar with p_-=E-p=m exp(-theta): quadratic T=(partial_+²phi)², X=m²(partial_+phi)² gives one-particle normal-ordered eigenvalue proportional to m³ e^-3theta (normalization 2 depending conventions). Derive properly and note X sign convention: displayed identity partial_-T+partial_+X=0 vs generic partial_-T=partial_+X requires Xgeneric=-Xprinted.
- g2 folded soliton degrees from scratch. Fold D4 triality, tauF at affine0, tauG common three outer nodes, tauH central. Bilinear equations (m1): F Fxx-Fx²=F²-H; G same=G²-H; H Hxx-Hx²=2(H²-F G³). Invariant central species mass sqrt6 solution F=G=1+E,H=1-4E+E² verifies degrees (1,1,2)=n_j^vee. Triple outer orbit with k=sqrt2, initial coefficients F1=3,G1=-1,H1=0 should produce degrees(3,3,6). Solve coefficient recursion up to6; nonlinear polynomial exact. This checks |O| n_j^vee. Need derive from field ansatz factor and degree scope. c2 explicit polynomial already verified phase1.
- Embedded SG n2,n3: h=2n, a=n tau even1+iE,odd1-iE. Sum even roots v, odd=-v, phi= (i/beta)v log((1+iE)/(1-iE))= -(2/beta)v arctanE (sign depends roots); v real with norm²=2n, real field. k=2m; substitution roots alpha_j·v=+/-2 leads SG along v with b_eff=2beta/|v|. Check numeric or symbolic n2/3 and canonical E.
- Source audit McGhee count; Kneipp/Olive actually state untwisted Lie duality in intro? Need rg/read. Fring time delays and Hollowood interaction coeff source check. PartVIII prp-stored-charge and prp-bianchi local references read and algebraically verify. Citations quantum conservation need read feigin9310022, dgz9201067,niedermaier9401078; distinguish perturbative explicit from general theorem.
- Sources pre-arXiv braden1990,freeman1991,dorey1991,fring1991,mikhailov1981,olive1985,Kac need attempt primary sources online; if unavailable explicit unclear, not claim read. Dorey sequel9110058 was read phase1 and verifies area rule with nonzero qualifier and fusing orbits.
- Full exercise outlines can mostly refer phase1 derivations; chapter11 phi->i phi sign directly from eom is +i (not -i generally). Theta interaction A conjugacy elementary cosine. Analytic tau polynomial continues through zero but logfield and ordinary charges do not, and physical-Theta reading of TW is interpretation, not established.

## Process/tool notes

Last active sessions: 16385 integrals (may finished; output has two numerical FAILs as above);11748 couplings;6827 old slow current solver;5668 trivial DomainMatrix test (probably fails NotAField, irrelevant). Read/poll as needed. Tools process namespaces isolate ps, so do not use ps in another exec to conclude sessions stopped. Final usable currents v3 already finished; don't repeat tests that already pass unless changed.

User asked whether to compact. Assistant answered yes, saving this checkpoint. Continue phase2 autonomously after compaction; do not mark report complete until remaining required checks attempted and report delivered. Keep frequent concise commentary.

## Completion checkpoint (supersedes the pending-task sections above)

2026-10-15: Phase-2 report completed at `briefs/referee_classical_09_11_codex_report.md`; final evidence archived in `briefs/referee_classical_09_11_codex_scripts/`.
- Item 1 finished: Theta proof/counterexamples, prior art and H6 replacements.
- Item 2 finished: Canonical energy/momentum and all stress components.
- Item 3 finished: Exact A2 local currents and converged breather integrations.
- Item 4 finished: Both time shifts; SG phase; unresolved higher-rank convention explicitly documented.
- Item 5 finished: Twisted Lax, boosts, one-particle current, centralizers and quantum-source scopes.
- Item 6 finished: A4/E6 cubic tensors and original Fring-Liao-Olive source.
- Item 7 finished: Embedded SG, C2/G2 degrees, charge count and Lie-duality source.
- Item 8 finished: Other reality claims, singularities and missing-saddle scope.
- Item 9 finished: Every citing sentence audited and all 17 exercise outlines.

All 29 overlap rows and the remaining-uncertainties section are present. Final integrity audit: `referee_scratch/classical_09_11/codex/codex_phase2_final_audit.txt`. Reviewed chapter/Appendix E snapshots and frozen phase-1/list hashes match. No book edits, git mutations, other referee reports, author scripts or subagents. The obsolete slow preliminary current solver was stopped; all final numerical jobs completed. No work remains on this referee report unless the editor requests follow-up or integration.
