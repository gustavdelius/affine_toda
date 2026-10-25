# Half-line lecture notes: combined checklist of the two referee reports

The two-phase Codex pass on `lecture-notes/half-line/` produced `referee_lectures_halfline_phase1_report.md` (**P1**, 25 ranked risks) and `referee_lectures_halfline_report.md` (**P2**, 12 coordinator items and 18 exercises). Both were committed as received (35ffaf7, 1460bd5). The referee's scripts are in `referee_lectures_halfline_scripts/`; rerunning them reproduces the logged numbers (phase-2 logs byte-identical).

This file merges the unique findings into one checklist, so that no passage is verified or edited twice. Line numbers refer to the lectures at 035351d (unchanged at the time of writing); match passages by content. Every item must be verified independently before it is edited (Rule 1 of `AGENTS.md`); the references to P1/P2 say where the referee's evidence and proposed text are, not that the item is right.

Agents: **A** = `index.qmd`, Lectures 1 and 2; **B** = Lectures 3 and 4; **C** = Lectures 5 and 6. Each agent owns the lecture files, the lecture's check/figure script `code/lectures/halfline_lN_f1_*.py` and its SVG. The coordinator owns `code/lectures/README.md`, the book findings (G5) and the record.

## General items (each agent, in its own files)

- **G1 Status language.** Apply the rows of P2 §1.1 that fall in your files: sentences that state as derived, exact, resolved or discovered what the book labels [Conjecture] or [Numerical], or that state a result more broadly than the book. Keep the mathematics; adapt the wording to the lecture's voice. P1 2.17.
- **G2 Script claims.** Apply the rows of P2 §12 in your files: say what the script actually tests. P1 2.25.
- **G3 Typed book numbers.** Link text such as `[Exercise 29.1](book:…)`, `[Chapter 31](book:…)`, `[Section 29.1](book:…)`, `[eq. 33.5](book:…)`, `[Proposition 33.1](book:…)`, `[Table 31.1](book:…)` breaks the rule in `AGENTS.md` ("never type a book chapter or statement number"). Replace with descriptive link text. Remove the typed prefixes "Exercise 1.1:" etc. from exercise titles (Quarto numbers the `exr-` divs). Where a heading's text names two book objects but links one (P2 §11: Exercises 1.3, 2.2), link each object or name only the one linked.
- **G4 Figures.** Where a caption or the drawn annotation/legend is wrong, fix the figure script, rerun it (it overwrites its SVG in place in this worktree) and `git restore` any other SVG whose diff is only ids and timestamps.
- **G5 Book findings (coordinator; separate book PR, not edited here).** (a) `part7-half-line/35-quantum-symmetric-pairs.qmd` ≈158: the Serre-obstruction proof uses the inverse of $k_i^2k_j$ inside the coideal without justification (P1 2.8). (b) Chapter 31 ≈145 and `exr-node-local`: omitting $\sqrt{n_j}$ keeps the equations node-local; "couple the nodes" is misleading (P2 §8). (c) `thm-k-reflection`: "agree at one point" should mean a fixed nonzero matrix element for every generic rapidity pair (P2 §5). (d) Chapter 36 ≈282: "So the rotor must run through all three … patterns" needs the qualification of P1 2.2.

## Index and Lecture 1 (agent A)

- **I1** Index:26 rotation $(x,t_E)\mapsto(-t_E,x)$ sends $x\le0$ to negative Euclidean time; use $(\tilde x,\tilde t_E)=(t_E,-x)$, as in the book's figure. Also the garbled clause "spatial translations conservation of momentum". P2 §3.
- **I2** Index:26–27 "can only emit pairs": qualify (zero-momentum one-particle terms). P1 2.21.
- **I3** Index:12, 30, 34, 46: rigidity "governs all higher-rank", "cannot be tuned continuously at all", Dirichlet "never integrable", antisoliton return and backgrounds "needed", "resolve both puzzles", WKB checks (the lectures do none), "literally exchange places", "only half … admit … vacua". P2 §1.1 rows; P2 §4 (mirror position is not universal).
- **I4** Index:46 attribution of the duality to Corrigan and Delius–Gandenberger (they proposed/conjectured it); image citations. P2 §9.1, §9.6.
- **I5** Index:69 "Every non-trivial calculation … verifies … from first principles". P1 2.25 replacement.
- **I6** Index:11 "solitons are their own antiparticles ($V\simeq\bar V$)": check and reword if it conflates the representation with the states (P1 §4, unverified lead).
- **L1-1** Continuity equations L1:20–51: the printed $(t,x)$ form has the densities and fluxes interchanged; with the book's normalization $Q_1+Q_{-1}=2E$ and $\Sigma_1=2\mathcal B$. Match the book's convention in `sec-half-line-kinematics`. P1 2.20, P2 §3.
- **L1-2** L1:62 "the boundary term … does not act": on asymptotic states it is a constant. P2 §3.
- **L1-3** L1:73–77 elasticity for arbitrary $\gamma_b^{(s)}$, and Exercise 1.1(2) (multiplicity of the maximal rapidity; "same mass" is not enough for several species). P2 §3, §11.
- **L1-4** L1:157 braid conversion: the display has $K$ in the **second** slot; the conversion needs $\check S=S_{12}P=PS_{21}$. P1 2.11.
- **L1-5** `eq-lec1-reflection-twisted` RHS: $\check R_{\bar\mu,\nu}\to\check R_{\nu,\bar\mu}$, also in `exr-lec1-twisted`. P1 2.11.
- **L1-6** L1:188–202 and the figure caption: the map is a reflection, and $\tilde E=-ip,\tilde p=-iE$ does not give $\theta=i\pi/2-\tilde\theta$. P1 2.20, P2 §3.
- **L1-7** L1:211–215 boundary state: one-particle term, and the measure's normalization. P1 2.21, P2 §3.
- **L1-8** L1:232, 238: $S^{ba}_{dc}\to S^{ba}_{cd}$ in the crossing derivation. Check against `eq-boundary-crossing-unitarity`. P2 §3 (new).
- **L1-9** L1:255 attributes the centred product form to CDRS, which prints the shifted form; Corrigan–Delius (`corrigan_delius1999`, (2.7)) has it. The exercise calls it Corrigan–Delius. P2 §3.
- **L1-10** `eq-lec1-excited-bootstrap` and `eq-lec1-bulk-bootstrap`: $S_{b\bar a}\to S_{ba}$ (diagonal reflection keeps the species). Check against `eq-boundary-bootstrap`, `eq-boundary-bootstrap-bulk`. P2 §3 (new).
- **L1-11** L1:279, 293 pole interpretation ("signals", fixed vs moving poles). P1 2.21.
- **L1-12** `eq-lec1-one-point-residue`: the sign is $-\frac i2\Gamma g$ at the reflection pole (Jacobian of $t=i\pi/2-\theta$). Check against `eq-one-point-residue` and GZ (3.45). P2 §2 (new).
- **L1-13** L1:369–403 Lee–Yang doubling: `eq-lec1-s2theta` is false ($S_{\rm LY}(2\theta)=-(1)(2)/((4)(5))$ in $h=6$ blocks); abandoned steps; replace by a complete derivation. P1 2.4.
- **L1-14** L1:452 $\Gamma^2=-2\sqrt3$ (convention $\operatorname{Res}S=i\Gamma^2$); "cardy" → "Cardy". CFT data and $g$ correct. P1 2.22, P2 §2.
- **L1-15** L1:353, 404, 464–485 status of the Lee–Yang factors (book `prp-ly-reflection` is [Numerical]); the shift $K_{b+2}$ is correct (P1 2.3). P2 §1.1.
- **L1-16** Exercise 1.3(4): $f_x=-(6x)(6-6x)$; the ratio is right. Heading link (G3). P1 2.4, P2 §11.
- **L1-17** Figure 1: "b=0: identity wall" → equal amplitudes, not equal boundary conditions; tick labels vs axis title; rotation map if drawn. P1 2.25.
- **L1-18** Script claims L1:162, 264, 266, 404. P2 §12.
- **L1-19** Typos: "inits" (L1:83), "therotated" (L1:289).

## Lecture 2 (agent A)

- **L2-1** Status/scope at L2:4, 73, 176, 313, 325. P2 §1.1.
- **L2-2** L2:126–130: $K_-=1$ uses the transpose and the vacuum asymptotics; involution also needs the bulk bracket and both reflection equations. P2 §8.
- **L2-3** L2:238 "invertible for arbitrary $(C_0,C_1)$" → generic spectral parameter (the matrix vanishes at $C=0,\lambda=1$). P1 2.19.
- **L2-4** L2:286 Branch I: supply the mixed-branch propagation. P1 2.19.
- **L2-5** Figure 2 caption: the orange slice vanishes only at $C_0=\pm1$; it is not Neumann at $C_0=0$. P1 2.25.
- **L2-6** L2:327 Dirichlet: the Cartan condition is the anticommutator $\{K,H\}=0$; state the scope of `prp-no-dirichlet`. P1 2.5.
- **L2-7** L2:354 script claim for the coweight extension. P2 §12.
- **L2-8** Exercise 2.1: name the $d_4$ vector representation and its spectator subspaces; the branch is $C_2^2=2$ in the altered normalization. P2 §8, §11.
- **L2-9** Exercise 2.2 heading links one of the two objects it names. P2 §11.
- **L2-10** Exercise 2.3: $D$ flips $C_j$ (not $\lambda$); coweight flips preserve $s$ (two classes of four); inversion flips $s$; together one orbit of eight; the $s=+1$ class is uniform plus three solitonic patterns. P2 §8.

## Lecture 3 (agent B)

- **L3-1** Status/scope at L3:4, 55, 118–127, 191–243, 238 ("solitonic" must exclude the uniform pattern); "mirror behind the wall" is not universal. P2 §1.1, §4.
- **L3-2** L3:127 and `exr-lec3-mirror-square`: $k=m_a\cosh\theta=2m\sin^2(\pi a/h)/s$; the derivative at L3:298 has $2(1-Cs)$; supply the constant-ratio step with its branch qualification. P1 2.10.
- **L3-3** L3:161: $1+\zeta^{-ja}AE_2$, not $\zeta^{-2ja}$. P1 2.10.
- **L3-4** L3:229–243 and Figure 3: extended-index signs ($\epsilon_{j+h}=(-1)^a\epsilon_j$); the interval table is wrong for two intervals; caption and legend at $\chi=0.25$ give $(+1,-1,-1)$. Make the L3 script test the stated assignment directly, not choose the best-fitting pattern. P1 2.1, P2 §4.
- **L3-5** L3:211–222 positivity: the roots need not lie at $u\ge1$; use the lower bound $\big((1-u)^2+s(1-u^2)\big)/(1+s)$ and the wall condition. P1 2.23.
- **L3-6** L3:256 rotor sentence: qualify as in L6-12 (consistent wording with agent C: "rotor-like tower", an interpretation). P1 2.2.
- **L3-7** `eq-lec3-linear-reflection`: $[-(a)(h-a)]^{-C}$ (the printed form gives $-1$ at $C=0$). P1 2.15.
- **L3-8** Exercise 3.1 heading links `exr-mirror-neumann` (the Neumann exercise); link the $C=\pm1$ mirror construction instead. P2 §11.
- **L3-9** Exercise 3.2: open-half-line positivity holds, but the $\rho_-=0$ breather has wall zeros; correct the $C=-1$ cases ($s<1$, $s>1$, $s=1$). P1 2.23.
- **L3-10** Exercise 3.3: the pole at $i\pi/2$ is at the strip edge. P1 2.23, P2 §4.

## Lecture 4 (agent B)

- **L4-1** Status/scope at L4:8, 10, 253 (scalar has CDD freedom), 267–280 ("Assuming the conjectured UV–IR relation"; the zero-coideal point approaches the uniform boundary), 310 (first-level states only). P2 §1.1.
- **L4-2** L4:24 direction ("a right-moving excitation hits the wall and returns left-moving"). P2 §5.
- **L4-3** After L4:39 define the bar on $\Phi$ (parity image, not complex conjugation); reconcile $\tilde\Phi$ and $\tilde\phi$. P2 §5.
- **L4-4** L4:56–58 "spin-1 non-local currents": weights depend on the coupling ($\lambda+1$ in sine-Gordon). P2 §5.
- **L4-5** L4:64–75: first order in $\lambda_b$, $+O(\lambda_b^2)$; the first-order computation motivates, not proves, exact conservation; restrict to simply-laced unless dual factors are given. P2 §5.
- **L4-6** L4:28 attribution: Mezincescu–Nepomechie (free-fermion integrals of motion), Delius–MacKay (coideal method). P2 §9.3.
- **L4-7** L4:149–150 "matching highest weights fixes the constant to 1": state the normalization hypothesis (book finding G5c). P2 §5.
- **L4-8** L4:158 de Vega–González-Ruiz solved the six-vertex reflection-matrix problem; GZ the physical amplitude. P2 §10.
- **L4-9** L4:238 BPTT parametrization wording (book's orientation and node conventions). P2 §10.
- **L4-10** L4:257–263 sign of $q$: "if and only if" → generic, matrix part; delete "and shown in @fig-lec4-coideal". P1 2.16, 2.25.
- **L4-11** L4:261 crossing display has the wrong argument and undefined indices; use $K_{\bar a}^{b}(\theta)=\sum S^{ba}_{cd}(i\pi-2\theta)K^d_{\bar c}(i\pi-\theta)$ or the vector form. P2 §5 (new).
- **L4-12** L4:299–304 abandoned determinant paragraph; eigenvector $(1,(-1)^ne^{-\vartheta})$; other reducibility branches exist. P1 2.12.
- **L4-13** Figure 4 caption: it shows neither crossing residuals nor energies. P1 2.25.
- **L4-14** Exercise 4.1 last sentence (parity allows $c_0=c_1=\pm q^{-1}$; genericity in $y$). P2 §11.
- **L4-15** Exercise 4.2: $0<n<\lambda/2$; $n=0$ is a strip-edge pole; analytic branch of $\epsilon_*$; input is the UV–IR conjecture. P2 §5, §11.
- **L4-16** (optional) Coproduct-side convention of "left coideal" vs Kolb's $\Delta(B)\subset B\otimes U$. P2 §10.

## Lecture 5 (agent C)

- **L5-1** Status/scope at L5:6, 7, 10, 58, 97–108, 114–131, 232–251. P2 §1.1, P1 2.17.
- **L5-2** After L5:108: mixed-branch propagation for $a_2$. P2 §6.
- **L5-3** L5:96, 209: $\epsilon_*$ and the generalized $q$-Onsager relations were published by Baseilhac–Belliard. L5:114: cite Gandenberger's explicit $a_n^{(1)}$ solutions (`gandenberger1999b`, hep-th/9911178); the referee did not establish formula-level equivalence, so check before saying more than "found explicit non-diagonal solutions". P2 §9.5, §9.
- **L5-4** L5:194–196 obstruction argument: the scalar $d_{ij}$ is missing, and the inverse of $k_i^2k_j$ is not known to lie in the coideal; argue in the represented algebra (or by the $a_2$ spectra). P1 2.8.
- **L5-5** L5:207 rank one: word reversal alone is not the reason; expand and cite the $q$-Dolan–Grady relation. P2 §6.
- **L5-6** L5:208 Kolb's parameter sets ("for the split diagram used here"); universal $K$: Balagović–Kolb finite type, Appel–Vlaar Kac–Moody with a generalized reflection equation. P2 §6.
- **L5-7** L5:209 "defining … obtained": the embedded coideal realizes the relations and is the image of $B_0$ under the character twist. P2 §6.
- **L5-8** L5:218–228: the classical $s=\prod C_j$ and quantum $s=\prod\hat\epsilon_j/\epsilon_*$ are different; the book pairs $\hat\epsilon=0$ with the product-$+1$ class and $s=+1$ with Neumann, and leaves $s=-1$ unpaired. P2 §6.
- **L5-9** Figure 5: the zero is at $\mu_0=\cos[\pi/(2(\lambda+1))]$, not $\mu\to1$; label the curve conditional on UV–IR. P1 2.17.
- **L5-10** Exercise 5.1(1): node 0 gives $K_{32},K_{21}$; $K_{31}$ needs node 1. P2 §11.
- **L5-11** Exercise 5.3: $Z^T=Z^{-1}$. P1 2.9.
- **L5-12** L5:131 script claim. P2 §12.

## Lecture 6 (agent C)

- **L6-1** Status/scope at L6:6–9, 22, 34, 53, 81, 101, 111, 119–134, 151, 165–181, 206–210, 229, 254. P2 §1.1.
- **L6-2** L6:20 history (Sasaki, Fring–Köberle, Kim; Kim assumed self-duality and says the earlier factors disagree with perturbation theory). P2 §7, §10.
- **L6-3** L6:22 and callout: Corrigan (rank one), Gandenberger and Delius–Gandenberger proposed/conjectured the identifications and duality. P2 §9.1.
- **L6-4** L6:48: the crossed boundary vector is $\sum_i|\bar v_i\rangle\otimes|\bar v_i\rangle$. P2 §7 (new).
- **L6-5** L6:53, 80–81 "telescopes": cite the book for the Gamma-product reduction; CDD freedom remains. P1 2.18.
- **L6-6** L6:94: $(-h-B/2)=(h-B/2)$ by the period $2h$. P1 2.7.
- **L6-7** L6:123–124: the CDD factor $P_{\omega/2}$ is chosen for exact duality; other choices share the classical and one-loop limits. P2 §1.2.
- **L6-8** L6:119 $c_+$ is tested for $a_2$ only. P2 §12.
- **L6-9** L6:151 Kim parenthesis ($a_2$, $a_3$ checked; does not fix higher orders). P2 §9.7.
- **L6-10** L6:165–167 solitonic background (imaginary-coupling stationary soliton with the vacuum shift, continued). P2 §7.
- **L6-11** L6:181: physical ranges $|n|,|m|<1/(2B)+1/2$ for $n+m>0$, diagonal $0\le n<3/(2B)+1/2$; negative labels are bootstrap labels. P1 2.24, P2 §7.
- **L6-12** L6:206–210 rotor: inertia correct, global $2\pi$ rotor not derived (interpretation). P1 2.2.
- **L6-13** L6:214–227 Coleman–Thun: $5-2\cdot2-1=0$; DTW do not print $P-2L-Z$; "machine precision" → tolerance $10^{-8}$; the residue identity and the $b_{-1,2}$ sentence are correct. P1 2.6, P2 §7, §12.
- **L6-14** Figure 6: "Quantized Vacuum Rotor" → rotor-like diagonal tower; duality of the proposed amplitudes. P1 2.25.
- **L6-15** Exercise 6.1: $u_{11}^2=2\pi/3$, $\bar u_{11}^2=\pi/3$. P1 2.14.
- **L6-16** Exercise 6.2: delete the false second expression for $K_1^N$; replace part 3 (e.g. compare Kim's $[1/2][3/2]$); "true" → "proposed", no two-loop claim. P1 2.13, P2 §11.
- **L6-17** Exercise 6.3: state labels, ranges, "same rapidities" meaning; degeneracy at $B=\frac12$ is correct. P2 §11.

## Coordinator

- **C1** `code/lectures/README.md` half-line rows: stale ids, $K_{-b}$ (actually $K_{b+2}$), vacuum energy $(3-\sqrt3)$ (actually $(4-3\sqrt3)$), attraction labels, and the new checks the agents add. P1 2.25.
- **C2** G5 book findings, recorded as pending for a separate book PR.
- **C3** Render, `check-refs.py`, the `## Integration` record below, `briefs/README.md`, and unfreezing the course.

## Integration (2026-10-25)

Integrated on branch `claude/referee-halfline-lectures` by three Opus integration agents, A (index, Lectures 1–2), B (Lectures 3–4) and C (Lectures 5–6), from `referee_lectures_halfline_integration_prompt.md`. Each agent verified every item before editing, with its own scripts, written from the lectures' and the book's definitions (`referee_lectures_halfline_scripts/integration_*`; none imports the referee's scripts). Their checkpoints, `referee_lectures_halfline_integration_{A,B,C}_checkpoint.md`, give one line per item with its evidence, a `## Table` and a `## Pending` section. The coordinator reviewed every diff, re-derived the items listed below, and wrote the README rows and this record.

### Verdicts

| Group | Items | Verified | Modified | Rejected |
|---|---|---|---|---|
| A: I1–I6, L1-1–L1-19, L2-1–L2-10, G1–G4 | 41 | 37 | 4 | 0 |
| B: L3-1–L3-10, L4-1–L4-16, G1–G4 | 30 | 26 | 4 | 0 |
| C: L5-1–L5-12, L6-1–L6-17, G1–G4 | 33 | 27 | 6 | 0 |
| **Total** | **104** | **90** | **14** | **0** |

No finding was rejected. The modified items, where the edit departs from the report's proposal:

- **L1-1** (continuity equations): the repair follows the book's spin labels (`sec-factorization`): a current $(T_{s+1},\mathcal X_{s-1})$ with $\partial_-T=\partial_+\mathcal X$ has spin $-s$, density $T-\mathcal X$ and flux $T+\mathcal X$; $Q_{\pm1}=E\pm P$, so $\widehat Q_1=2(E_{\rm bulk}+\mathcal B)$. The referee's counterexample ($T=t+x$, $\Theta=t-x$) satisfies both the printed and the corrected form and refutes nothing; $T=t-x$, $\Theta=t+x$ refutes the printed form. The finding itself is right.
- **L1-6** (rotation): derived as the continuation of a boost with rapidity $-i\pi/2$, which gives both rapidity assignments, instead of the report's prose.
- **L2-4** (mixed branches): a shorter propagation argument; node 2's spectator equations kill $k_{12},k_{21}$ for any $C_2$.
- **I3**, **L3-6**, **L6-12**: the coordinator's common wording for the rotor (inertia $I=6/(\beta^2m_1)$, a rotor-like tower, a global rotor is an interpretation) and for the status sentences of the index.
- **L3-10**: the pole of the self-conjugate exercise is on the strip edge; it is double in $\theta$ and simple in $p$.
- **L4-11**: the book's vector form of crossing-unitarity instead of the report's component display, to avoid a second index convention.
- **L4-16**: the report's claim that Kolb uses the opposite coproduct is **wrong**. Kolb (arXiv:1207.6036, eqs. (3.3)–(3.4) and Prop. 5.2) uses the book's coproduct; his coideals are right coideals because their generators are built around $f_j$. The lecture says so in one parenthesis.
- **L5-3**: Gandenberger hep-th/9911178 is cited for explicit non-diagonal $a_n^{(1)}$ solutions, without a claim of formula-level equivalence with the Toeplitz matrix.
- **L5-4**: for $a_2$ the obstruction is proved from the disjoint spectra of $k_i^2k_j$ on $V$ and $\bar V$ (no inverse needed); for general $n$ the conjugation argument is done in the represented algebra.
- **L6-7**: the CDD statement was verified; the alternative factor is called admissible only for unitarity and crossing (see book finding B8).
- **L6-13**: the Coleman–Thun count is made explicit (four bulk lines and two boundary segments, two loops, order 2, lowered to 1 by the zero), from the topology of Dorey, Tateo and Watts's Lee–Yang diagram.
- **L6-16**: Exercise 6.2(3) uses Kim's actual $a_2$ proposal $[\frac12][\frac32]$ instead of the report's $B(2-B)/2$ construction.

### Coordinator's checks

- Re-derived by hand: the $a_2$ interval-to-pattern table (literal integer indices), $ks=2m\sin^2\frac{\pi a}h$ and the constant ratio $k-mC$, the rotation as a continued boost, the Lee–Yang doubling identity and the bootstrap sign, the $C=-1$ cases of the breather exercise, the node-locality of `eq-boundary-k-equation` without $n_j$, and the soliton-bootstrap test of the CDD factors (`integration_coord_cdd_bootstrap.py`).
- Read the source behind the self-duality history: Delius–Gandenberger write that the earlier papers "had assumed that the reflection amplitudes should be invariant under $B\rightarrow 2-B$".
- Fixed two small points the agents left: a non-sequitur about the asymptotics of $E$ in Lecture 3, and the description of the $q$-Dolan–Grady check in Lecture 5.
- Reran all six lecture scripts: 10 + 9 + 7 + 7 + 8 + 11 = 52 PASS, each in at most about 8 s. `check-refs.py lecture-notes/half-line` passes. Rewrote the six half-line rows of `code/lectures/README.md` (item C1): current ids, the shift $K_b^{(1)}=K_{b+2}$, $E_{\rm vac}=\frac{2m}{\beta^2}(4-3\sqrt3)$, and the new checks. Every id in those rows exists.
- Rendered the course, HTML and PDF, in a scratch copy with `render-scratch.sh --pdf --project lecture-notes/half-line`: no warnings, no unresolved cross-references in the HTML, and none in the 97-page PDF. Checked the corrected legend of `fig-lec3-mirror-images` by eye.

### Phase 1 against phase 2

Every major mathematical error was already in the blind phase-1 list: the vacuum sign table, the Lee–Yang doubling step and $\Gamma^2$, the continuity equations and the Euclidean map, the twisted-equation domains, the Dirichlet argument, the mirror $k$ and positivity proof, the Robin factor at $C=0$, $Z^T$, the $a_2$ fusion angle, the Coleman–Thun arithmetic, and the conjectures presented as settled. Phase 1 alone also found the stale script README.

The directed phase 2 added:
- index and sign errors inside derivations: the order $S_{cd}^{ba}$ in the crossing derivation, the sign of the one-point residue, the bootstrap channels $S_{ba}$, the factor 2 in $Q_1+Q_{-1}$, the space of the crossed boundary vector, and the crossing display of Lecture 4;
- exercise defects (1.1, 2.3, 4.1, 4.2, 5.1) and the coupling-dependent weight of the "spin-1" currents;
- provenance and scope: Kim's self-duality history, the Dorey–Tateo–Watts formula, the Mezincescu–Nepomechie and Baseilhac–Belliard credits, Kolb's parameter sets and universal $K$-matrices, and the new citation of Gandenberger hep-th/9911178.

As in the Part III experiment, the blind phase found the substance and the directed phase found the provenance and the smaller derivation slips. The referee made two errors that the integration caught: the Kolb coproduct (L4-16) and the void counterexample (L1-1). P1 2.15's remark that the book "separates the Neumann case correctly" is too generous (B5).

### Book findings (pending; for a separate book PR, not edited here)

- **B1** `part7-half-line/35-quantum-symmetric-pairs.qmd` ≈158: the proof after `prp-serre-obstruction` conjugates by $\kappa^{-r}$ without showing $\kappa^{-1}\in B_{\hat\epsilon}$. Repair in the represented algebra (Cayley–Hamilton), or by the disjoint spectra for $a_2$. Confirmed by the coordinator and agent C.
- **B2** Chapter 31 ≈145 and `exr-node-local`: without $n_j$ the equations stay node-local, with $C_j/\sqrt{n_j}$; for $d_4^{(1)}$ the central node then needs $C_2^2=2$. "Couple the nodes" and "no longer decouple" are wrong. Confirmed by the coordinator and agent A.
- **B3** `thm-k-reflection`, hypothesis (iii): "agree at one point" should read "agree on a fixed nonzero matrix element for every generic pair of rapidities".
- **B4** Chapter 36 ≈282: "So the rotor must run through all three … patterns" is part of an interpretation; "would have to".
- **B5** `eq-linear-reflection` (ch. 32): $-[(a)(h-a)]^{-C}$ gives $-1$ at $C=0$; write $[-(a)(h-a)]^{-C}$.
- **B6** `eq-sg-neumann-eps`: name the branch $2\cos(\pi\lambda/2)$ of the square root.
- **B7** Chapter 36 ≈288: Dorey, Tateo and Watts do not print a general $P-2L$ rule lowered by one per zero; "as applied by".
- **B8** `sec-an-nondiagonal-k` (ch. 34): "The replacement also satisfies the soliton-bootstrap condition" holds only for odd $n$. For $n=2,4$ both $P_{\omega/2}$ and $P_{(\omega+1)/2}$ fail $\sigma(\mu)=\prod_{k=1}^n\sigma(\mu+T-k\omega)$ (`integration_C_cdd_bootstrap.py`, `integration_coord_cdd_bootstrap.py`; the earlier check covered $a_3$ only). Whether that diagonal condition is the right test for $A^{(+)}$ at all is open; for even $n$, including Gandenberger's $a_2$, the book does not say whether $A^{(+)}$ satisfies a soliton bootstrap.
- **B9** Chapter 36 ≈158: "The same applies to the other self-dual proposals" can be read as saying that they agree at one loop; Kim says the Fring–Köberle and Sasaki $a_3$ factors do not.
- **B10** `sec-boundary-lee-yang`: "Unitarity and crossing-unitarity are linear in $K$"; unitarity is quadratic (the conclusion, invariance under $K\to-K$, holds).
- **B11** (optional) Chapter 33: for $\lambda>2$ the UV–IR curve has several zeros of $\hat\epsilon$; $\cos\frac{\pi}{2(\lambda+1)}$ is the one nearest the uniform point.

### Not settled by this pass

The global quantization of the rotor, a derivation (rather than a count) of the Coleman–Thun diagram, an all-rank proof of the scalar-factor telescoping, the four pre-arXiv sources (Cherednik, Sklyanin, Bernard–LeClair, Letzter), and the precise scope of Bernard–LeClair's $U_q(\hat{\mathfrak g}^\vee)$ remain as the report leaves them. The lectures now state these points with the corresponding qualification. No status label changed: the lectures carry none, and the book is not edited here.
