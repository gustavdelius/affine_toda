# Brief: find the missing real-coupling vacua of the solitonic boundaries (half-line open problem 6)

You are a mathematical physicist with expertise in classical integrable field theory, Hirota solutions of affine Toda theory and numerical analysis. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files. Do NOT spawn subagents. This brief was written from the book, from abstracts and from machine summaries of papers. Treat every claim in it as unverified until you have rechecked it.

## The problem

Open problem 6 of `sec-boundary-open-problems` (`part7-half-line/35-particle-reflection.qmd`):

> **Missing vacua and missing charges.** For the patterns of `tbl-missing-vacua` (for $a_3^{(1)}$ the alternating ones), find the classical vacua that Delius and Gandenberger's amplitudes presuppose, perhaps from the complex two-soliton configurations of `cnj-missing-charge`.

`tbl-missing-vacua` counts, for $a_1,\dots,a_7$, the solitonic patterns ($\prod_jC_j=1$, not all $+1$) that have a continued **single-soliton** vacuum: $1,3,5,10,10,21,21$ out of $1,3,7,15,31,63,127$. For $a_3^{(1)}$ the missing ones are $(+,-,+,-)$ and $(-,+,-,+)$. Their weights, for example $\epsilon_1+\epsilon_2\in\Lambda_2$, are carried by no single species-2 soliton. Delius and Gandenberger's amplitude formula nevertheless assigns them the amplitude of their class $\{b,h-b\}$.

## What the book already has

Read `sec-real-coupling-vacua`, `sec-boundary-missing-charges`, `sec-images-boundary` and `sec-images-cn` (ch. 31); `sec-solitonic-boundaries` (ch. 35); `cnj-missing-charge` and `prp-static-saddles` (`part5-semiclassics/17-thimbles.qmd`). Scripts: `textbook_code/boundary_vacua.py` (stationary pairs, `prp-dg-weights`, the table, the self-conjugate family) and `boundary_images.py`. Facts to use:

- A real-coupling static solution is $\vec\phi=-\frac1\beta\sum_j\vec\alpha_j\ln\tau_j$ with all $\tau_j>0$ on $x\le0$. Singularities may hide at $x>0$, and only static configurations qualify.
- For $a_3^{(1)}$ the four missing weights of $\Lambda_2$ are reached statically only with species $\{1,2,3\}$, at energy $(1+\sqrt2)M_2$. Bound states of two species-1 solitons at rapidities $\pm iu$, $u\to\pi/4$, have energy $\to M_2$ and cover all of $\Lambda_2$ at $t=0$ (`prp-static-saddles`, citing Harder et al.).
- The boundary energy (`eq-boundary-energy`) and the $\chi$-independence argument in the proof of `prp-self-conjugate-vacuum` ($dE/d\xi=0$ for a family obeying one boundary condition).

## Literature

- **Bowcock–Perkins, hep-th/9909174 (`bowcock_perkins2003`). Read this first.** According to a machine summary, they derive a Bogomolny bound on the half-line ($E\ge W|_{x=0}-W|_{-\infty}+\mathcal B|_{x=0}$, stability iff $A_i\le1$ in Bowcock's sign, i.e. $C_i\ge-1$). They also classify vacuum configurations and classical reflection factors for $A_r^{(1)}$ up to $r=5$, including multi-soliton and "exceptional" backgrounds (a 4-soliton $1\bar1 2\bar2$ background for one pattern, 6- and 8-soliton backgrounds in larger rank). Their tables 1–2 mark some patterns "not diagonal" (no consistent diagonal classical scattering). The book's table explicitly excludes these multi-soliton backgrounds. **Their tables may already contain vacua for the missing patterns.** Extract them and translate to book conventions ($A_i=-C_i$, `tbl-boundary-dictionary`) before doing anything else.
- Perkins's Durham PhD thesis, "Affine Toda field theories on a half-line" (etheses.durham.ac.uk/4560): the longer version of the same classification.
- Bowcock, hep-th/9609233 (`bowcock1998`), sec. 3.2: static configurations with several stationary solitons.
- Fateev–Onofri, hep-th/0207152 (IJMPA 18 (2003) 879; not in `references.bib`). They extract exact classical vacuum solutions from the $b\to0$ limit of boundary one-point functions, as $\tau$-functions summing over binary configurations $\sigma_j\in\{0,1\}$. Only all-positive boundary signs, and non-simply-laced algebras, are treated; the mixed-sign cases are "deferred to a separate publication" (which we did not find). The method transfers to $a_n$.
- Harder–Iskandar–McGhee and McGhee (`mcghee1994`); Beggs–Johnson, "Inverse scattering and solitons in $A_{n-1}$ affine Toda field theories" I and II (hep-th/9610104, hep-th/9803248): new single-soliton solutions with extra topological charges, which fill some holes but leave others.

## Suggested plan of attack

1. **Read and translate Bowcock–Perkins.** For $a_2,\dots,a_5$, list every sign pattern with its vacuum type (flat, single pair, multi-soliton, exceptional, none) and energy. Reconcile with `tbl-missing-vacua`. If they found vacua for the alternating $a_3$ patterns, recompute them (regularity on $x\le0$, the boundary condition, the energy) and stop the search for $a_3$ there.
2. **Direct numerical search.** Independently of the literature, minimize the static energy functional on $[-L,0]$ with the boundary potential at $x=0$ and $\vec\phi\to0$ (or Neumann) at $x=-L$. Use a discretized gradient flow or Newton continuation from many random starts. Do this for the missing patterns of $a_3,a_4,a_5$. Record whether the infimum is attained by a regular configuration, its energy, and whether it is a family (a modulus). Compare with $E_{\rm vac}$ of (`eq-vacuum-energy`) for the same class $b$.
3. **Fit to Hirota data.** If a numerical vacuum is found, fit $\tau_j$ to Hirota forms: multi-soliton static configurations with complex $\xi$, and the continued two-soliton bound states at the fusing angle (`cnj-missing-charge`) in their static limit. This tests Delius's suggestion that the mechanism that supplies the missing soliton states also supplies the missing vacua.
4. **Consistency with the amplitudes.** For each vacuum found, compute the linearized reflection factor (`sec-classical-reflection-factors`). Compare it with the $B\to0$ limit DG (5.7) of (`eq-solitonic-factor`) for its class (`prp-solitonic-amplitude` (iii)). Agreement would confirm DG's assignment for the missing patterns.
5. **If there is no vacuum.** Show it: for example, a Bogomolny-type bound that is not attained, or a minimizing sequence that runs into a singularity at $x=0$. Then state what that means for DG's prediction.

## Known traps

- At real coupling there is no shift symmetry. Patterns related by (`eq-coweight-shift-c`) are different theories, so each pattern needs its own search.
- Normalization of the boundary energy: the book uses $\mathcal B(0)=0$; Bowcock omits the constants $-1$. Energies differ by a pattern-dependent constant.
- Bowcock's $A_i=-C_i$; Delius 1998a's $\epsilon$ with $C\equiv-\epsilon$; `prp-dg-weights` maps weights to patterns two-to-one.
- A regular solution can have $\tau_j(0)=0$ at isolated parameter values; at those points the pattern changes (`prp-solitonic-vacua`).

## Deliverable

- The translated Bowcock–Perkins classification for $a_2$–$a_5$, with verdicts (recomputed or not).
- For each missing pattern of $a_3,a_4$ (and $a_5$ if time permits): a vacuum (formula or numerical profile, energy, moduli) or evidence that none exists.
- A revised `tbl-missing-vacua` with a column for multi-soliton vacua, and proposed text for `sec-boundary-missing-charges`.
- Scripts.

## Feasibility (assessment when the brief was written, 2026-10-11)

**High.** It is classical and numerical, uses existing scripts, and part of the answer may already be in Bowcock–Perkins. The link to `cnj-missing-charge` (step 3) is the research part and may stay open.
