# The open problems of Part VII: briefs and feasibility ranking

Written 2026-10-11. Section `sec-boundary-open-problems` (§35.7, `part7-half-line/35-particle-reflection.qmd`) lists twelve open problems for the half-line. Each problem except 10 (defects, which goes into the planned defects Part) has an agent-ready brief here. Every brief contains the problem as stated in the book, the relevant ids and scripts, a literature section, a plan of attack, known traps, required checks, a deliverable, and a feasibility assessment.

**Caveat on the literature.** The searches were done on 2026-10-11. Most entries were checked only at the level of the arXiv abstract; some rest on a machine summary of the full text (marked in the briefs). One equation was read directly: Ahn–Kim–Rim eq. (10). New references are marked "not in `references.bib`". Following the repository rule, a worker must re-derive or recompute every claim before the book is changed.

## Ranking by feasibility

Feasibility here means the likelihood that one focused agent session (or a few) produces a result that can go into the book, either a solution or a sharp partial answer. It is not the importance of the problem.

| Rank | Problem | Brief | Feasibility | Model | Effort | Why |
|---|---|---|---|---|---|---|
| 1 | 3. Neumann from the soliton side | `halfline_open_03_neumann_soliton_route.md` | High | Opus | xhigh | Every ingredient is in the repository. The matrix part of crossing-unitarity already factors (`prp-an-nondiagonal` (iii)). The scalar factor is a Gamma-product computation with exact checks at $n=1$ (GZ Neumann) and $n=2$ (Gandenberger). The outcome confirms or breaks `cnj-neumann-amplitude` and `cnj-neumann-coideal`. |
| 2 | 12. Rank one | `halfline_open_12_rank_one.md` | High (parts b and the lead), medium (a), low (c as a proof) | Opus for the lead, (a) and (c); Sonnet for (b) | xhigh (Opus), high (Sonnet) | Extending the excited-state checks of boundary sine-Gordon is routine numerics. The lead found while writing: AKR's normalization $\mu_B^2=\frac\mu2\cot\frac{\pi b^2}2$ of the $(+)$ boundary appears to continue to the $\hat\epsilon=0$ point $\cos\frac\pi{2(\lambda+1)}$. If that holds, the $E=\frac B2$ versus $E=0$ puzzle has independent non-perturbative support without a two-loop calculation. |
| 3 | 6. Missing vacua | `halfline_open_06_missing_vacua.md` | High | Sonnet; with problem 1 (classical) in one agent | high | Classical numerics: energy minimization on the half-line, plus fits to Hirota forms. Bowcock–Perkins (2003) classified vacua for $A_r^{(1)}$, $r\le5$, including multi-soliton and "exceptional" backgrounds that the book's table excludes, so part of the answer may already exist. The link to `cnj-missing-charge` is the research part. |
| 4 | 1. $\prod_jC_j=-1$ | `halfline_open_01_negative_product_vacua.md` | Medium-high (classical), low (quantum) | Sonnet for the classical part (with problem 6); Opus for the amplitudes, after problem 3 | high (Sonnet), xhigh (Opus) | Same numerics as problem 6, shared code, and Bowcock–Perkins's tables again. Amplitudes would need a new soliton reflection matrix or the $s=-1$ scalar factor (a by-product of problem 3). |
| 5 | 11. Soliton-preserving scalar factors | `halfline_open_11_soliton_preserving_scalar.md` | Medium-high (scalar factor, breathers), medium (interpretation) | Opus | xhigh | Standard machinery with an exact $n=1$ (Dirichlet) check. Regelskis–Vlaar 2018 classify all solutions of both reflection equations for $A^{(1)}_{N-1}$, which removes the guesswork. A tension surfaced while writing: classically the soliton-preserving condition at real coupling is Neumann on real fields, but for $n=1$ the continuation gives Dirichlet. This needs resolving and may yield a second route to the Neumann amplitudes. |
| 6 | 7. Non-unitary half-line | `halfline_open_07_nonunitary_half_line.md` | Medium-high (a), medium (b), low (c) | Opus for (a) and (b); (c) deferred until problem 11 is done | xhigh | Boundary Krein-unitarity is a small extension of `prp-krein-fundamental`. Boundary Bethe–Yang numerics for $a_2$ are standard. Recent exactly solved $\mathcal{PT}$-symmetric boundary spin chains are a guide. Whether Dirichlet walls are integrable touches the constructive theory. |
| 7 | 8. Masked bound states | `halfline_open_08_masked_bound_states.md` | Medium | Opus | max | A concrete first target exists: the unexplained $a_2$ pole $\frac12+(2n-1)\frac B2$. The bulk residue technology is now mature (Dorey–Polvara 2022, Homrich et al. 2022), and Bajnok–Böhm–Takács give the boundary Landau analysis. The boundary residue formula with reflections at zeros must first be calibrated on boundary Lee–Yang and sine-Gordon. |
| 8 | 2. Beyond $a_n^{(1)}$ | `halfline_open_02_beyond_an.md` | Medium (first steps), low (all algebras) | Opus | xhigh | The $d_n$ $\hat\epsilon=0$ scalar factor and breather bootstrap are routine. Identifying the classical boundary is not, since $\vec\phi=0$ is not a vacuum for $|C_j|=1$. For $c_n$, folding plus Fateev–Onofri (2003), who treat boundary $B_r$, $C_r$, $BC_r$ theories with continuous parameters, gives an independent comparison. |
| 9 | 9. Dynamical boundaries | `halfline_open_09_dynamical_boundaries.md` | Medium | Opus | xhigh | Belliard–Fomin (2012) built dynamical K-matrices for $\mathcal O_q(a_2^{(1)})$, and Hills (York thesis, 2016) built $a_2^{(1)}$ boundaries by dressing with defects. The work is translation, verification and classical–quantum matching. Best run with the defects Part. |
| 10 | 5. Completeness of coideals | `halfline_open_05_coideal_completeness.md` | Medium ($a_n$ vector), low (general) | Opus | max | Regelskis–Vlaar 2018 solved the reflection equations completely for the vector representation, so the quantum census is translation. The classical census with automorphism-twisted $\sigma$ is finite. "Exhaust all integrable boundaries" in general is a programme. |
| 11 | 4. Duality on the half-line | `halfline_open_04_half_line_duality.md` | Medium-low | Opus | max | The CFT mechanism (FZZ, AKR, Fateev–Onofri) and a boundary TBA comparison are concrete. Identifying the duals of the solitonic boundaries may or may not succeed, and there is no algebraic framework for a coideal-level statement. |

## Dependencies and suggested order

- **Problem 3 first.** Its scalar-factor technique transfers to problems 11 and 2 ($d_n$). Its $s=-1$ by-product feeds problems 1, 4 and 5.
- **Problems 6 and 1 together.** One energy-minimization code, one reading of Bowcock–Perkins.
- **Problem 12 early.** The AKR lead is cheap to check, and its result changes the wording of `sec-uniform-boundary` and of open problem 12 itself.
- **Problems 11 and 7(c) are linked.** The soliton-preserving amplitudes are the input for testing whether Dirichlet walls are integrable.
- **Problem 9 with the defects Part.** It shares the dressing construction and the anomaly caveat (Penati–Polvara).
- Per the subagent lessons, run at most two or three agents at a time, and follow every result with an adversarial referee pass before any status label changes.

## Status (2026-10-13): consolidation (section A of `halfline_open_continuation.md`)

Branch `half-line-consolidation`, from the unmerged PR #12 branch. The user chose consolidation only (no new attacks). Prompts: `briefs/halfline_consolidation_common.md` and `halfline_consolidation_A<k>_prompt.md`; agent reports `halfline_consolidation_A<k>_report.md` with scripts in `halfline_consolidation_A<k>_scripts/`. Verdicts are recorded in the Integration tables of the round-1 reports (03, 07, 06_01, 11, 02).

| Item | Who | Outcome | Book |
|---|---|---|---|
| A1, $s=-1$ scalar factor | coordinator | verified; **new**: the duals of $K^{(-)}_a=-(a)(h-a)K_a$ are CDRS (5.16) for all $n\le7$ | `prp-an-nondiagonal-scalar` extended; ch. 35 "The third matrix"; `an_nondiagonal_scalar.py` checks 8–10, `particle_reflection.py` check 5 |
| A2, problem 11 scalar factors | Opus xhigh + coordinator | verified with corrections ($k$ closed form, $n=1$ sign, breather phase is a convention); **new**: SP soliton bootstrap; $(-\frac T2)_I^2$ excluded for even $n$ | ch. 33 `prp-sp-scalar`; ch. 35 `eq-sp-breather`; open problem 11; `sp_scalar.py` |
| A3, problem 11 classical claims | Opus max + coordinator | both verified; claim 5 strengthened to [Theorem] `prp-no-dirichlet` (all untwisted $\hat g\ne a_1$); Appendix D Delius 1998b (1.4) item rewritten (the book's earlier reason was wrong) | ch. 30 `prp-soliton-preserving-half`, `prp-no-dirichlet`; ch. 31; open problem 7(c); `boundary_dirichlet.py` |
| A4, problem 7(b) thresholds | coordinator | verified (UN VV 0.698, UN $V\bar V$ 0.666, NN VV 0.536, NN $V\bar V$ 0.501; three solitons broken); "NN VV unbroken at 0.67" refuted | `sec-boundary-krein` "Two and three solitons"; `boundary_krein.py` check 7 |
| A5, problem 6 numerical $K$ | Opus high + coordinator (own ODE) | DG (5.7) verified about every regular vacuum ($a_2$–$a_5$); singular backgrounds: path-dependent, so round 1's "BP Table 2 corrected" is **withdrawn** | ch. 31 `prp-solitonic-linear`, `sec-boundary-missing-charges`; `solitonic_linear.py` |
| A6, problem 2 leftovers | Opus xhigh + coordinator (own BVP vacuum and linearization, own one-loop and functional-equation checks) | verified: $\hat\epsilon=0$ is $C\equiv+1$ for $d_4,d_5,d_6$ (convex energy ⇒ unique vacuum, closed-form kinks, $K_a^{\rm cl}=\prod_bS_{ab}(\theta,1)$ for every particle); **new**: $d_n$ $s=\pm1$ scalar factors, the relevant sign is $\varsigma=(-1)^ns$; $\varsigma=+1$ with a one-loop CDD factor is the dual of the uniform amplitude and reproduces Kim–Cho's $d_4$ one loop; $\varsigma=-1$ is $K_1/K_1^{\rm cl}$; the book's "Kim's one-loop $d_4$" citation was wrong (it is Kim–Cho hep-th/9505138); ch. 33 `prp-dn-r` nesting fixed; c_2 step: folding argument confirmed, nothing new for the book | ch. 30 sec-boundary-potential; ch. 33 `sec-dn-vector-k` restructured, new `prp-dn-nondiagonal-scalar`; ch. 35 `sec-dn-uniform` (`prp-dn-uniform-vacuum` [Theorem], `prp-dn-uniform-classical`, `eq-dn-neumann-amplitude`), open problem 2; Appendix D self-dual item; `dn_uniform.py` checks 6–9, new `dn_nondiagonal_scalar.py` |
| A7, problem 8 coincident diagrams | — | not attempted (research) | — |

The three new [Theorem]s (`prp-soliton-preserving-half`, `prp-no-dirichlet`, `prp-dn-uniform-vacuum`) and the Appendix D changes were checked by an Opus agent and by the coordinator, but have had no Gemini referee pass.

## Status (2026-10-12): solving paused, integration under way

On 2026-10-12 the user asked to stop attacking further open problems and to consolidate and integrate the solutions into the book.

**Reports in hand:**

| Tag | Problem | Report |
|---|---|---|
| 03 | Neumann from the soliton side | `briefs/halfline_open_03_report.md` |
| 12a | rank one, lead and part (a) | `briefs/halfline_open_12a_report.md` |
| 12b | rank one, part (b) | `briefs/halfline_open_12b_report.md` |
| 06_01 | missing vacua; classical part of problem 1 | `briefs/halfline_open_06_01_report.md` |
| 07 | non-unitary half-line, parts (a) and (b) | `briefs/halfline_open_07_report.md` |
| 11 | soliton-preserving scalar factors | `briefs/halfline_open_11_report.md` |
| 02 | beyond $a_n$: $d_n$, $\hat\epsilon=0$, and the map | `briefs/halfline_open_02_report.md` |
| 08 | masked bound states: the $a_2$ pole | `briefs/halfline_open_08_report.md` |

Each report has a Gemini verdict (`*_referee.md`), except 08, whose prompt is pending.

**Opus follow-ups:** run for 03 + 12a and for 12b + 06_01. The follow-up prompts for 07, 11 and 02 (`*_opus_prompt.md`) were not launched. Their items are verified by the integrating agent instead.

**Not started:** problems 5, 4, 9, 7(c) and the problem 1 amplitudes.

**Integration rule** (see the repo's practice): re-derive or recompute every finding before it enters the book, and record verdicts in an "Integration" section appended to each report.

## Instructions for a coordinating agent

You coordinate work on these briefs. You delegate; you do not solve the problems yourself. The repository is `/home/gustav/Git/affine_toda`.

**Delegation.**
- Start each worker with the Agent tool, `subagent_type: "general-purpose"`, and give it the brief file's full contents verbatim as its prompt.
- Set `model` and `effort` as in the table above. The user explicitly asks for these effort levels. Use `"opus"` or `"sonnet"` for `model`; Fable is not available to the user.
- Never use `subagent_type: "fork"`: a fork ignores the `model` setting.
- Every brief already tells its worker not to spawn subagents and not to edit the repository.

**Schedule.** Check the remaining usage before each launch, and launch only agents that can finish before the quota runs out. Run at most three agents at once; a stopped agent can be resumed with SendMessage after a usage reset. Suggested batches:
1. Problem 3 (Opus). Problem 12, lead and (a) (Opus). Problems 6 and 1, classical part, as one Sonnet agent.
2. Problem 12 (b) (Sonnet). Problem 11 (Opus). Problem 7 (a), (b) (Opus).
3. Problem 8 (Opus, `max`). Problem 2 (Opus; pass it problem 3's scalar-factor method). Problem 1, amplitudes (Opus, only if problem 3 produced the $s=-1$ scalar factor).
4. Problems 5 and 4 (Opus, `max`). Problem 9 (Opus), unless the user says the defects Part should take it.

If an Opus worker reports that it is blocked on an idea rather than a computation, start a second, fresh Opus agent at `max`. Give it the brief, the first worker's report, and the instruction to attack the blocking point from a different angle. Do not retry for runtime or convention errors.

**Reports.** Workers report in their final message. Save each report verbatim as `briefs/halfline_open_<tag>_report.md`, with the date, model and effort in a header line. `<tag>` is the zero-padded problem number, with a suffix when a problem is split or two are combined (`03`, `12a`, `12b`, `06_01`). Large reports may arrive as a persisted tool-results file; read them from there. Copy the worker's scripts from its scratch directory to `briefs/halfline_open_<tag>_scripts/`, because scratch directories do not last and the referees must rerun the scripts.

**Referee passes (changed 2026-10-11: Gemini referees, to save Claude quota for the workers).** The referees run on Gemini (High) in Antigravity, under a separate Gemini coordinator that follows `briefs/halfline_referee_coordinator.md`. The two coordinators communicate only through files; that brief defines the protocol. For every report that claims a solved problem, a new [Theorem] or [Numerical] statement, or a correction to the literature:
- Write `briefs/halfline_open_<tag>_referee_prompt.md` in the format given in the Gemini brief: tier, brief, report and scripts paths, the claims as a numbered list with their claimed labels, and coordinator notes on suspected weak points. Do not start a Claude referee or Sonnet rerun agent.
- **Tier A** (no [Theorem] claim, no correction to the literature): the Gemini verdict is the referee pass.
- **Tier B** (a [Theorem] or a literature correction): after the Gemini verdict `briefs/halfline_open_<tag>_referee.md` appears, read it critically and choose the Opus items yourself. Gemini's "Items for an Opus pass" is input, not the list (on 12a it missed the real weak points). Look especially for claims it confirmed although its description differs from the worker's, for qualitative claims it confirmed loosely, and for label upgrades. Write the items to `briefs/halfline_open_<tag>_opus_prompt.md` and start a fresh Opus agent at `max` on them. Give it the brief, the report, the Gemini verdict, and the instruction to attack those items, to distrust steps that merely look plausible, and not to spawn subagents. Save its verdicts as `briefs/halfline_open_<tag>_referee_opus.md`.
- A Gemini verdict file is complete when it exists; it ends with "Referee pass complete." When no more prompts will come, create the empty file `briefs/halfline_referees_done` so that the Gemini coordinator stops.
- The Gemini model is newer than Claude's training data and has not yet been calibrated on this project. Read its verdicts critically too. Where it and the worker disagree, settle the point before anything moves toward the book.

**What you do not do.**
- Do not edit book chapters, `references.bib` or the scripts. Integration into the book happens later, with each finding verified first.
- Do not commit.

**Final message.** One table: problem, model and effort used, outcome (solved / partial / blocked / not started), referee verdict, and the paths of the report files. Then a short list of the findings that would change the book, each with its evidence and referee status.

## Literature found (not yet in `references.bib`)

| Reference | arXiv | Relevant to |
|---|---|---|
| Fateev–Onofri, Boundary one-point functions, scattering and background vacuum solutions in Toda theories, IJMPA 18 (2003) 879 | hep-th/0207152 | 1, 2, 4, 6 |
| Fateev, Normalization factors, reflection amplitudes and integrable systems | hep-th/0103014 | 1, 2, 4 |
| Fateev–Zamolodchikov–Zamolodchikov, Boundary Liouville field theory I | hep-th/0001012 | 4, 12 |
| Riva, Boundary bootstrap principle in two-dimensional integrable QFTs, NPB 604 (2001) 511; thesis | hep-th/0102163, hep-th/0106268 | 2, 8 |
| Belliard–Fomin, Generalized q-Onsager algebras and dynamical K-matrices, J. Phys. A 45 (2012) 025201 | 1106.1317 | 9 |
| Regelskis–Vlaar, Solutions of the $U_q(\widehat{sl}_N)$ reflection equations, J. Phys. A (2018) | 1803.06491 | 5, 11 |
| Appel–Vlaar, Tensor K-matrices for quantum symmetric pairs | 2402.16676 | 5 |
| Avan–Doikou, Boundary Lax pairs for the $A_n^{(1)}$ Toda field theories, NPB 821 (2009) 481 | 0809.2734 | 5, 9 |
| Arnaudon–Avan–Crampé–Doikou–Frappat–Ragoucy, General boundary conditions for the $sl(N)$ and $sl(M\vert N)$ open spin chains | math-ph/0406021 | 5, 11 |
| Lima-Santos–Malara, regular reflection K-matrices for non-exceptional vertex models | nlin/0412058 | 5, 11 |
| Bajnok–Böhm–Takács, On perturbative QFT with boundary, NPB 682 (2004) 585; Boundary reduction formula | hep-th/0309119, hep-th/0207079 | 7, 8, 12 |
| Dorey–Polvara, From tree- to loop-simplicity in affine Toda theories I, JHEP 09 (2022) 220 | 2206.09368 | 8 |
| Homrich et al., Nonperturbative anomalous thresholds | 2212.06157 | 8 |
| Chenaghlou, Second order quantum corrections to the classical reflection factor of the sinh-Gordon model, IJMPA 16 (2001) 4613 | hep-th/0108077 | 12 |
| On the quantum reflection factor for the sinh-Gordon model with general boundary conditions | hep-th/0004121 | 12 |
| Kormos–Takács, boundary form factors of exponential operators in sinh-Gordon | 0801.0962 | 12 |
| Bajnok–Palla–Takács, Finite size effects in boundary sine-Gordon theory; Boundary sine-Gordon model (review) | hep-th/0108157, hep-th/0211132 | 12 |
| Corrigan–Zambon, Infinite dimension reflection matrices in the sine-Gordon model with a boundary, JHEP 06 (2012) 050 | 1202.6016 | 9 |
| Bajnok–Simon, Solving topological defects via fusion | 0712.4292 | 9 |
| Hills, Generating boundary conditions for integrable field theories using defects, PhD thesis, York (2016) | — | 9 |
| Perkins, Affine Toda field theories on a half-line, PhD thesis, Durham | — | 1, 6 |
| Beggs–Johnson, Inverse scattering and solitons in $A_{n-1}$ affine Toda field theories I, II | hep-th/9610104, hep-th/9803248 | 6 |
| Exact solution of a non-Hermitian $\mathcal{PT}$-symmetric spin chain; a spin chain with $\mathcal{PT}$-symmetric boundary couplings | 2301.06004, 2406.10334 | 7 |
| W. Wang, Quantum symmetric pairs (survey) | 2112.10911 | 5 |

Author lists marked "et al." or omitted were not confirmed in the search; check them before adding bib entries.
