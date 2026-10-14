# Referee brief: the half-line consolidation round (Part VII, 2026-10-13)

**Referee:** GPT6-Astra in Codex, started by hand.
**Branch:** `half-line-consolidation`.
**Report to write:** `briefs/halfline_consolidation_referee_report.md`.

## Rules

You are a referee. You do not draft book text.

- **Do not edit** anything under `part*/`, `appendices/`, `references.bib` or `textbook_code/`.
- **Where to write.** Write only your report, plus your own scripts in a new directory `briefs/halfline_consolidation_referee_scripts/`.
- **Replacement text.** Where you find an error, propose replacement text in the report. The coordinator re-derives every finding before it changes the book.
- **Work independently.**
  - Re-derive by hand, or recompute with your own code. Running our scripts is useful, but it is not a check: our scripts are what is being checked.
  - You may import the book's library modules from `textbook_code/`, for example `boundary_common.py`, `reflection_matrices.py` and `smatrix_simply_laced.py` (the S-matrix blocks).
  - Read the coordinator's own verification scripts (`briefs/halfline_consolidation_A6_scripts/coordinator/`) only after you have done your own check of an item.
- **Git and files.** Use git only read-only (`git diff`, `git log`, `git show`). Never checkout, stash, reset, add or commit: the worktree is shared with running Claude sessions. Use the directory `briefs/halfline_consolidation_referee_scripts/` for your own scripts and scratch files.
- **Literature.** Read every paper you rely on at the source, from the arXiv TeX or PDF, and quote equation numbers as printed. ArXiv ids are given below and in `references.bib`.
- **Running scripts.** Run the book's check scripts from `textbook_code/` (`cd textbook_code && python3 <script>.py`). Each prints one `PASS`/`FAIL` line per check, and each takes up to about 3 minutes.
- **Verdict for each item:** one of
  - **CONFIRMED** (re-derived or independently recomputed; say how);
  - **CONFIRMED WITH CHANGES** (give the change);
  - **WRONG** (give the counterexample or the error, and replacement text);
  - **UNVERIFIABLE** (say what is missing).

  Also check whether the status label ([Theorem], [Numerical], [Conjecture]) is justified, and whether the result is already in the literature (prior art).
- **Precision.** Numerical claims need agreement at the stated precision, at points you choose yourself (complex rapidities, several couplings), not at the points our scripts use.

## Background

- **The round.** It verified and integrated leftovers of the Part VII open problems (§35.7, `sec-boundary-open-problems` in `part7-half-line/35-particle-reflection.qmd`).
- **Diffs.** `git diff origin/main...half-line-consolidation -- part7-half-line appendices references.bib` shows every change to the book.
- **Records.** For each item, the agent report and the coordinator's verdicts are in:
  - `briefs/halfline_consolidation_A<k>_report.md`;
  - the round-1 report the item came from, in its "Integration" table (`briefs/halfline_open_NN_report.md`, sections dated 2026-10-13);
  - the summary in the status block of `briefs/halfline_open_problems.md`.

**Book conventions** (getting these wrong is the usual source of false alarms):
- Blocks:
  - $(x)_I=\sin\frac{\pi(\mu+x)}{2T}/\sin\frac{\pi(\mu-x)}{2T}$, with $\mu=-\frac{iT\theta}\pi$ and $T=\frac{h\omega}2$;
  - particle blocks $(x)=\sinh(\frac\theta2+\frac{i\pi x}{2h})/\sinh(\frac\theta2-\frac{i\pi x}{2h})$.
- Continuation $\omega\to-\frac2B$; $q=-e^{-i\pi\omega}$.
- Real-coupling boundary potential: `eq-boundary-potential` in ch. 30, with $C\equiv+1$ the "uniform" boundary.
- Bulk potential: $V=\frac{m^2}{\beta^2}\sum_jn_j(e^{\beta\vec\alpha_j\cdot\vec\phi}-1)$.

## Tier A: theorems and corrections to the literature (highest priority)

**A1. `prp-no-dirichlet`** (ch. 30, `30-boundary-conditions.qmd` around line 185), [Theorem].
- **Statement:** the real-coupling Dirichlet condition $\vec\phi(0)=\vec\phi_0$ is not classically integrable, in Sklyanin's sense with a field-independent $K$, for any untwisted $\hat g\ne a_1^{(1)}$.
- **Check the proof step by step,** in particular:
  - that a Dirichlet condition forces $\sigma=-\mathrm{Ad}_K$ on $\mathfrak h$;
  - the "cycle condition" $\lambda^{2h}=1$, and why it fails;
  - why sine-Gordon ($a_1$) escapes;
  - the spin $\pm2$ check;
  - the partial-Dirichlet remark.
- **Loopholes:** does the hypothesis (field-independent $K$, linear $\sigma$, which representation) leave any? Is the result stronger or weaker than what the proof shows?
- **Prior art:** Bowcock–Corrigan–Dorey–Rietdijk hep-th/9501098; Corrigan's review hep-th/9612138 (`corrigan1996`, which the book quotes); Delius hep-th/9809140. Is the statement known? Is it contradicted there?

**A2. `prp-soliton-preserving-half`** (ch. 30, around line 258), [Theorem], with `eq-soliton-preserving-half` and the "Half-solitons at the wall" paragraph.
- **Statement:** the soliton-preserving (image-form) boundary of Delius hep-th/9809140 is integrable for every coweight, not only for even ones, with $K_+=e^{i\pi\nu\cdot H}$.
- **Re-derive the proof.** Check that the boundary value is a legitimate boundary condition. Check the half-soliton three-soliton images (ch. 31, `sec-images-preserving`).

**A3. Appendix D, Delius 1998b entry** (`appendices/D-corrections.qmd`, item "(1.4)–(1.5), the image form and the boundary condition").
- **The claim:**
  - Delius's image condition (1.4) gives his boundary condition (1.5) only for even coweights;
  - the book's earlier reason (the boundary value "is not a vacuum") was wrong;
  - the other coweights are integrable with $K_+=e^{i\pi\lambda\cdot H}$.
- **Read Delius's (1.4), (1.5) and (5.1)–(5.2) at the source** and decide whether this is a fair statement of an error in the paper.

**A4. `prp-dn-uniform-vacuum`** (ch. 35, `35-particle-reflection.qmd` around line 304), [Theorem], with its proof.
- **Uniqueness by convexity.** For real $C_j\ge0$ the static energy is claimed strictly convex on finite-energy fields on $x\le0$, so a static solution satisfying the boundary condition is the unique vacuum. Check:
  - strict convexity, including the role of the boundary term;
  - that "critical point ⇒ global minimum" is legitimate on this function space;
  - the invariance argument.
- **The closed forms.** Check the $d_4$ Tzitzéica and $d_5$ sinh-Gordon kinks, $e^{f(0)}=\frac98$ and $4-2\sqrt2$, and the energies $(16\sqrt2-8-6\sqrt6)\frac m{\beta^2}$ and $\big(32(1+\sqrt2)\sin\frac\pi8-24-4\sqrt2\big)\frac m{\beta^2}$, in the book's root convention (stated in the proposition).
- **Literature.** Is the convexity argument standard (e.g. for Toda boundaries in Bowcock–Corrigan–Dorey–Rietdijk, or Fujii–Sasaki)? Then a citation is needed.

**A5. Literature attributions introduced in this round.** Check each against the source:
- **(a)** Kim–Cho hep-th/9505138 (`kim_cho1996`): their eqs. (29) and (40) are the one-loop $d_4^{(1)}$ Neumann amplitudes $K_1$, $K_2$. The book previously attributed this calculation to Kim hep-th/9506031 (`kim1995`).
- **(b)** Kim hep-th/9506031: the self-dual $d_n$ vector conjecture (his fig. 7) equals Castro-Alvaredo–Fring's $\tilde R_1(\theta+i\pi)$ (hep-th/0308004, their (KimCF)).
- **(c)** Appendix D "Self-dual reflection factors", the new last sentence:
  - Kim–Cho's self-dual $d_4$ amplitude, which they call uniquely determined by duality, agrees with `eq-dn-neumann-amplitude` at one loop but not beyond;
  - the same holds for Kim's $d_5$–$d_7$ vector amplitude.
- **(d)** ch. 33 `prp-an-nondiagonal-scalar` and ch. 35 "The third matrix": for $a_2^{(1)}$ the $s=-1$ breather amplitude equals Gandenberger's (5.17) (hep-th/9806003). The duals of $K^{(-)}_a=-(a)(h-a)K_a$ are the conjecture (5.16) of Corrigan–Dorey–Rietdijk–Sasaki (hep-th/9404108, `cdrs1994`) for $C\equiv-1$, for all $n\le7$. Is (5.16) stated there for general $a_n$, and is "their conjecture" the right description?

## Tier B: new numerical results

For each item, recompute at least one representative case with your own code.

**B1. `prp-dn-nondiagonal-scalar`** (ch. 33, around line 257) and `eq-dn-neumann-amplitude` (ch. 35, `sec-dn-uniform`, "The other two matrices").
- **(a) The sign rule.** The scalar data depend on $\varsigma=(-1)^ns$, not on $s$. Build the $d_n$ vector $K$ from the book's coideal (`boundary_common.coideal_gen`, `intertwiners`, `dn_vector`) for $n=4,5$.
- **(b) The data.** Check $U_\varsigma$, the crossing function $\Psi_\varsigma$ and the singlet coefficient $c^d_\varsigma$.
- **(c) The scalar factors.** Check that `eq-dn-nondiagonal-scalar` solves both functional equations. Use high precision and handle the infinite-product tails properly; a naive truncation converges like $1/J$.
- **(d) The breather amplitudes** `eq-dn-nondiagonal-breather`, including the CDD factor $\Xi$.
- **(e) The continued amplitudes.** Check the identities $K^N_1=K_1|_{B\to2-B}$ and $K^{(-)}_1=K_1/K_1^{\rm cl}$, and that both satisfy unitarity, crossing-unitarity with $S_{11}$ and the bootstrap.
- **(f) One loop.** Do your own $O(\beta^2)$ expansion against Kim–Cho (29), (40); note that their $B=\frac{\beta^2}{2\pi}/(1+\beta^2/4\pi)$.
- **(g) Assessment.** The CDD factor $\Xi$ is chosen to match one loop, not derived. Is the text honest about this? Is there a derivation, for example from a boundary bootstrap, a self-duality requirement or a semiclassical argument?

**B2. `prp-dn-uniform-classical`** (ch. 35, around line 320).
- **The claim:** for $d_4,d_5,d_6$ the linearized reflection about the $C\equiv+1$ vacuum is diagonal, with $K^{\rm cl}_a=\prod_bS_{ab}(\theta)|_{B=1}$ for every particle.
- **Recompute** with your own multichannel integration for at least $d_4$ and one other.
- **"No bound state."** This has been checked only by the agent's script (`dn_uniform.py` check 8). Check it independently: no normalizable solution with $0<\Omega<m_{\min}$.
- **The poles.** Check the interpretation of the physical-strip poles as redundant poles of the vacuum's tail ($p=\frac i2m_2$; for $d_6$ also $m_4$ and $m_2+m_4$).

**B3. `prp-sp-scalar`** (ch. 33, around line 342), with `eq-sp-scalar-eqs`, `eq-sp-scalar`, and ch. 35 `eq-sp-breather`: the soliton-preserving scalar factors and breather amplitude. Check:
- the crossing functions $\Phi_{1,2}$;
- the scalar factor at one $(n,l)$;
- the breather amplitude's block content (`textbook_code/sp_scalar.py` is the book's check).

**B4. `prp-solitonic-linear`** (ch. 31, around line 237).
- **The claim:** the linearized reflection about every regular continued vacuum ($a_2$–$a_5$) is DG's (5.7) (hep-th/9904002) at $B=0$, independent of $\chi$.
- **Also check** the three bullet cases at the end of `sec-boundary-missing-charges`, in particular "pair endpoints are path-dependent". This is why the round withdrew the claim that Bowcock–Perkins (hep-th/9909174) Table 2 is "corrected".

**B5. ch. 35 `sec-boundary-krein`, "Two and three solitons"** [Numerical], with `boundary_krein.py` check 7.
- **Thresholds:** two-soliton Krein-unitarity breaks at $\omega\approx0.698$ (UN VV), $0.666$ (UN $V\bar V$), $0.536$ (NN VV) and $0.501$ (NN $V\bar V$), and three solitons are broken.
- **The mirror convention:** check the wall convention that makes the double-row transfer matrices commute.

**B6. ch. 33 `prp-an-nondiagonal-scalar`** (the $s=-1$ part, `eq-an-nondiagonal-scalar-minus`, `eq-an-minus-breather`) and ch. 35 "The third matrix". Check:
- that $K^{(-)}_a=-(a)(h-a)K_a$ satisfies unitarity, crossing-unitarity and the bootstrap;
- the breather identity at one $n$.

## Report format

- **One section per item** (A1–A5, B1–B6), each with:
  - the verdict;
  - what you did, with script names;
  - numbers at your points;
  - prior art found;
  - proposed replacement text where needed.
- **Summary table at the top:** item | verdict | label OK? | action needed.
- **A final section** "Other problems noticed", for anything else wrong in the changed text: notation, cross-references, claims stronger than the evidence, missing citations.
