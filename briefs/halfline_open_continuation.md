# Brief: continue the open problems of Part VII (coordinator, second round)

Written 2026-10-12 by the coordinator of the first round, so that a fresh coordinator can take over without the earlier conversation. Read this brief, then the documents in "Where things are", before doing anything.

## Where things stand

The book is a Quarto textbook on affine Toda field theory, in the repository `/home/gustav/Git/affine_toda`. Part VII (`part7-half-line/`, ch. 28–35) treats the theory on the half-line. Its open problems are in §35.7, `sec-boundary-open-problems` in `part7-half-line/35-particle-reflection.qmd`.

**Round 1 (2026-10-11/12).**
- Agents attacked problems 3, 12, 6+1 (classical part), 7(a, b), 11, 2 and 8.
- Gemini refereed every report, and Opus did follow-up passes on the claims that needed them.
- The coordinator re-verified each result and integrated it into the book.
- All of this is on branch `half-line-open-problems-integration`, **PR #12** (https://github.com/gustavdelius/affine_toda/pull/12). Check whether it has been merged before you branch.

**What the book now contains**, by problem:

| Problem | Result in the book |
|---|---|
| 3 | `prp-an-nondiagonal-scalar`: the $s=+1$ scalar factor. Neumann amplitudes are derived for $n=2,3,4$, with antisoliton $K^{\rm T}$, crossing function $\Psi_s$, CDD factor (Kim at $n=3$) and semiclassics. |
| 12 | The uniform sinh-Gordon boundary is $M/M_0=\cos\frac{\pi b^2}2$ ($E=\frac B2$); `prp-bsg-excited` holds for all $n$; a paragraph on the unexplained $\nu_{-j}$ poles. |
| 7 | `prp-boundary-krein` [Theorem]; $\mathbb Z_{n+1}$ sector pairing; the ch. 40 bulk remark. |
| 11 | Soliton-preserving matrices in book conventions; they reflect particles $a\to h-a$, so they are neither Dirichlet nor Neumann. |
| 2 | $d_n$ $\hat\epsilon=0$ scalar factor, breather amplitude and $K_1$ (= Castro-Alvaredo–Fring's $\hat R_1$). |
| 8 | The $a_2$ "unexplained" pole is a boundary Coleman–Thun pole. |
| 6+1 | Bowcock–Perkins boundary-singular vacua as energy infima; the $a_3$ bound. |
| Appendix D | New entries: AKR, BRZ, Mattsson–Dorey, BPT. |

**On 2026-10-12 the user paused new attacks to consolidate.** Ask the user whether to resume attacks (section B) or only finish consolidation (section A).

## Where things are

| What | Where |
|---|---|
| Ranking, original briefs, coordinator instructions | `briefs/halfline_open_problems.md`; the coordinator section and its "Status (2026-10-12)" block apply, except where this brief overrides them |
| One brief per problem | `briefs/halfline_open_NN_<name>.md` |
| Worker reports | `briefs/halfline_open_<tag>_report.md`, tags `03`, `12a`, `12b`, `06_01`, `07`, `11`, `02`, `08`. **Each ends with an Integration table**: verified / modified / rejected / not integrated, with reasons. Section A comes from those tables. |
| Worker scripts | `briefs/halfline_open_<tag>_scripts/` (the scratch originals are gone or ephemeral) |
| Gemini verdicts | `briefs/halfline_open_<tag>_referee.md` |
| Opus follow-up prompts / verdicts | `_opus_prompt.md` / `_referee_opus.md`. Verdicts exist for 03, 12a, 12b, 06_01. Prompts for 07, 11, 02 were written but not run; the coordinator handled their items while integrating, except the ones listed in A. |
| Gemini referee-coordinator brief | `briefs/halfline_referee_coordinator.md`. The stop file `briefs/halfline_referees_done` was created on 2026-10-12. **Delete it before writing new referee prompts**, and tell the user to restart the Gemini coordinator in Antigravity. |
| Book check scripts | `textbook_code/` (index in `textbook_code/README.md`). Round 1 added `an_nondiagonal_scalar.py`, `boundary_krein.py`, `shg_boundary_tba.py`, `soliton_preserving.py`, `dn_uniform.py`, `boundary_sine_gordon.py` check 12b and `particle_reflection.py` section 9. |
| Memory | `/home/gustav/.claude/projects/-home-gustav-Git-affine-toda/memory/`. Read `MEMORY.md`, especially `half-line-open-problem-briefs`, `gemini-referee-workflow`, `subagent-workflow-lessons`, `verify-brief-findings`, `literature-corrections-appendix`, `book-structure`, `no-fable-access`. |

## The work

### A. Consolidation: integrate what round 1 left out (no new attacks needed)

Each item already has worker output. What is missing is an independent check (with book code, not the worker's) and, for anything that goes into the book, a check script in `textbook_code/`. Do these yourself or with one agent each; most are a session or less.

1. **The $s=-1$ scalar factor** (`03` report §3, §5; claims 6, 10).
   - Verify $A^{(-)}$ and the breather amplitude DG(3.11)$\cdot\sigma_1$, which continues to $K^{(-)}_a=-(a)(h-a)K_a^{\rm unif}$, against Gandenberger's $n=2$ (6.4).
   - Use the Stirling/Hurwitz product evaluation of `an_nondiagonal_scalar.py` and extend that script.
   - Then update "The third matrix" in ch. 35, open problem 3, and line 45 of ch. 35 ("For $s=-1$ … not yet established").
2. **Problem 11, scalar factors and breathers** (`11` report §3–4; claims 7, 8). The crossing function, the Gamma-product $a,b$, the CDD candidates $\sigma_A,\sigma_B$ and the breather amplitude $-K_0E_\eta$ are not independently checked. $K_0$ is a *fit*. Recompute, and settle the CDD factor with the soliton-preserving soliton bootstrap; the worker suggests adapting `briefs/halfline_open_03_scripts/sol_boot.py`.
3. **Problem 11, classical claims** (claims 4, 5; prompt items 1–2 in `halfline_open_11_opus_prompt.md`).
   - **Claim 4:** `sec-soliton-preserving` should allow coweights outside $2\Lambda_W^\vee$ with $K_+=e^{i\pi\mu\cdot H}$. This contradicts Appendix D's Delius 1998b (1.4) entry, so be careful.
   - **Claim 5:** real-coupling Dirichlet is not classically integrable for $n\ge2$. Gemini's proof drops the constants $e^{\beta\alpha_j\cdot\phi_0/2}$. Check prior art (Bowcock–Corrigan–Dorey–Rietdijk hep-th/9501098).
   - Claim 5 also bears on 7(c).
4. **Problem 7(b) thresholds** (`07` report Table 2). Two-soliton reality thresholds (UN VV $\omega\approx0.69$, NN VV $\approx0.535$) and broken three-soliton states. Gemini bracketed UN VV independently. Port a compact transfer-matrix check into `boundary_krein.py`, then add the numbers to `sec-boundary-krein`. The open question whether the $\zeta\ne1$ sectors are physical (topological charge on the interval) is the item from `halfline_open_07_opus_prompt.md` that was never done.
5. **Problem 6, numerical $K$-matrices** (`06_01` report §3; `06_01_referee_opus.md` item 1).
   - The DG (5.7) linearization check over all regular pair vacua: two independent codes agree, but there is no book script. `tbl-boundary-types` still says "not checked here".
   - The $K$ limits about singular backgrounds, e.g. $a_5$ `++----` non-diagonal against BP's Table 2 (computed by the Opus τ-method, rerun by the coordinator).
   - The τ-method script is in the Opus work directory, which may be gone; the method is described in the verdict.
6. **Problem 2 leftovers.**
   - The classical boundary of the $d_n$ amplitude: the $d_4$ Tzitzéica / $d_5$ sinh-Gordon kink vacua and the linearized reflection match `02` report §3. Only the worker's script exists.
   - The $s=\pm1$ $d_n$ scalar factors, with Kim's $d_4$ Neumann one loop as the test.
   - The $c_2$ step.
7. **Problem 8 leftovers**: the coincident-diagram cases ($a_4$ with $n_2=-n_3$; $a_3$ with $b=2$) and a field-theoretic residue rule at zeros. These are research, really part of B.

### B. New attacks (only if the user asks)

The original ranking and briefs still apply (`halfline_open_problems.md`, "Ranking"). Not started:

| Problem | Effort | Note |
|---|---|---|
| 5, coideal completeness | Opus `max` | Regelskis–Vlaar 2018 classify the reflection-equation solutions |
| 4, duality on the half-line | Opus `max` | |
| 9, dynamical boundaries | Opus | **Ask the user** whether the defects Part should take it |
| 7(c), Dirichlet walls | | Now informed by problem 11: the soliton-preserving matrices are not real Dirichlet; see A3 |
| 1, particle amplitudes | | Low value now: no $\prod C=-1$ boundary has classical $K\to1$, and most have non-diagonal $K$ |

When you launch a brief, add a coordinator note on top giving what round 1 learned for that problem (pointers to the relevant reports).

## How to work (lessons from round 1; these override older instructions where they differ)

**Delegation**
- Agent tool, `subagent_type: "general-purpose"`, the brief verbatim plus a short coordinator note, and model/effort as in the ranking table. **Never Fable** (the user has no access).
- Every prompt says: do not edit the repo, do not spawn subagents, deliver the full report in the final message, and list full script paths.
- At most **three agents at once**. The session limit was still hit four times. Stopped agents resume with SendMessage and keep their progress, so tell them to check the output files of their own background jobs.
- **Give each agent its own scratch subdirectory**, e.g. `<scratchpad>/pNN/`. The scratchpad is shared.
- **Copy an agent's scripts to `briefs/halfline_open_<tag>_scripts/` as soon as it reports.** Leave out unchanged copies of `textbook_code/*.py` (check with `cmp`) and `__pycache__`.

**Refereeing**
- Write `briefs/halfline_open_<tag>_referee_prompt.md` for the Gemini coordinator. Tier B means a [Theorem] or a literature correction.
- **Choose the Opus follow-up items yourself** (`_opus_prompt.md`); do not take Gemini's list.
- Gemini's numerics are good. It is unreliable in four ways:
  - it confirms by restating;
  - it twice presented a rerun of the worker's own script as independent (once "confirming" a false error in a published paper);
  - it misquoted a paper and built a "literature correction" on the misquote;
  - it upgrades labels to [Theorem] too easily.
- Opus follow-ups found real problems every time.

**Integration (the user's rule)**
- Never put a finding into the book without re-deriving or recomputing it yourself first.
- Prefer the **book's own code** (`boundary_common`, `particle_reflection`, `smatrix_simply_laced`) over the worker's for that check.
- Read quoted papers at the TeX line yourself before any Appendix D entry.
- Record verdicts in the report's Integration table.
- Numerical claims in the book need a passing check in `textbook_code/`, listed in its README.
- Corrections to published papers go into Appendix D in the same change (`appendices/D-corrections.qmd`, ordered by year).

**Technical pitfalls met in round 1**
- `smatrix_simply_laced.py`, `boundary_coideal.py` and `reflection_matrices.py` run all their checks on import, which is slow. To get only the definitions, `exec` the source up to the first check (e.g. split at `'\nB = 0.37'`).
- `boundary_common.rcheck` is the generic solver and is far too slow for $d_n$ tensor products. Use structural arguments, or the closed form `eq-dn-vector-r`.
- `boundary_coideal.k_closed` takes the wrong branch for odd $N$ when $\cos\frac{\pi\omega}2<0$. Use the analytic-branch `kmat` in `an_nondiagonal_scalar.py`.
- Avoid the special couplings $B=\frac1{2k}$ in residue checks (level crossings), and normalize eigenvectors and intertwiners by a fixed entry, not by `max|·|`.
- Never `pkill -f`/`pgrep -f` from Bash: it kills its own shell. Kill by PID.
- Book notation:
  - the book's $\Phi$ in ch. 33 is taken (the crossing function is $\Psi_s$);
  - write the fundamental range as $\frac\pi2(\lambda+1)$, not $\eta_N$;
  - no `[Status]` label inside prose with a cross-reference;
  - the status tags are bold labels at the start of a statement.
- Render in a scratch copy: rsync the repo without `.git`, `.claude`, `_site`, `.quarto`, then `quarto render --to html`, and grep `_site` for `?@`. The user usually runs `quarto preview`.

**Git:** branch from `main` (or from the PR branch if it is unmerged). Commit or open a PR only when the user asks. Commit messages end with the attribution line given in the session's system reminder.

## Deliverable of each session

Update the "Status" block of `briefs/halfline_open_problems.md` and the Integration tables. Finish with a short report to the user: what was verified and integrated, what was rejected (with the reason), what is still pending, and whether anything changes the book's status labels.
