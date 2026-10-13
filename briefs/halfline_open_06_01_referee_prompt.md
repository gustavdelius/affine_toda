# Referee prompt: half-line open problems 06 and 01 (classical part)
Tier: B
Brief: briefs/halfline_open_06_missing_vacua.md and briefs/halfline_open_01_negative_product_vacua.md (steps 1–4 and 6 only)
Report: briefs/halfline_open_06_01_report.md
Scripts: briefs/halfline_open_06_01_scripts/ (29 files)
- **Final code:** `redmin2.py` (`python3 redmin2.py <n> <nstarts>`), `tauvac.py`, `halfline_min.py`, `linrefl2.py`, `kcheck.py`, `bp_table.py`, and `exc4*.py`. The `*.out` files are the worker's outputs.
- **Superseded:** `redmin.py`, `run_red.py` and `slsq.py`.

Scope: all of problem 6 and the classical part of problem 1. The worker was Claude Sonnet. Particle amplitudes (problem 1 step 5) are not in scope. Tier B because of claim 4 (a proof) and because claim 7 says DG's published amplitude assignment does not fit the $a_4$ vacua.

## Claims to referee

1. [Literature translation] **Bowcock–Perkins (hep-th/9909174) tables 1–2 in book conventions.** They contain vacua for every sign pattern of $a_2$–$a_5$, including the missing patterns of `tbl-missing-vacua` and all $\prod C_j=-1$ patterns. The dictionary is $A_i=-C_i$, BP list one pattern per cyclic orbit, and $E_{\rm book}=-2(n+1)C_{\rm BP}-2\sum_jC_j$ in units of $m/\beta^2$.
2. [Numerical] **The energies.** The worker's reduced minimization over the static Hirota family reproduces BP's energy for every pattern class it ran (tables in report §1) and never goes below it. $a_5$ (+++−−−) is not converged.
3. [Claim] **The vacua are singular.** Apart from flat $\vec\phi=0$ and the regular intervals of the continued pair, the vacua have $\tau_j(0)=0$ on nodes with $C_j=-1$. $\phi$ diverges logarithmically at $x=0$, and the energy infimum is finite but not attained.
4. [Theorem] **The rigorous bound (report §2.3).**
   - $E_P=E_Q+4\sum_{j:C_j=+1}(e^{\psi_j(0)/2}-1)$, with $Q\equiv-1$, and BP's Bogomolny bound gives $E_Q\ge0$, so $E_P\ge-4\#\{C_j=+1\}$.
   - The bound is saturated by the $d\to1$ endpoint of the self-conjugate family for $a_3$ (+−+−) ($E=-8$) and (+−−−) ($E=-4$).
   - For $h$ odd, $\vec\phi=0$ is the unique minimum for $C\equiv-1$.
5. [Numerical] **The a_3 alternating patterns** (the missing ones of problem 6). The vacuum is the $d\to1$ endpoint of $\tau_j=1+(-1)^jde^{2x}$, i.e. species 2 only. Its classical $K$ equals the flat $C\equiv-1$ value and DG's class-$b=2$ limit. So DG's assignment is consistent, and `cnj-missing-charge` is not needed for vacua.
6. [Numerical] **DG (5.7) checked.** The linearized reflection about every regular continued-pair vacuum of $a_2$–$a_5$ (30 patterns) equals the $B\to0$ limit of DG (5.7) to $10^{-11}$, with off-diagonals ≤$10^{-11}$. The book says this was unchecked.
7. [Numerical, indicative] **The a_4 missing pattern (++−+−) and its rotations.**
   - The vacuum is BP's exceptional $1\bar12\bar2$ background ($E=-10.4162$; Hirota parameters $d_1=0.2579$, $d_2=0.5637$).
   - Its linearized $K$ is non-diagonal (off-diagonals ≈0.14), so it does not reproduce DG's diagonal $b=2$ amplitude.
   - The solver is unreliable below $\tau(0)=10^{-5}$.
8. [Numerical] **$\prod C=-1$ (problem 1).**
   - Every pattern has a BP vacuum that is singular at $x=0$, except $C\equiv-1$ with $h$ odd, which has flat $\vec\phi=0$ with no negative modes (boundary quadratic-form argument, report §5).
   - Classical $K$ is diagonal only for the middle-soliton type, where it equals the $C\equiv-1$ value; otherwise it is non-diagonal (the $a_2$ (++−) numbers are in §3).
   - $\prod C=-1$ vacua are not limits of regular solutions of the same pattern.
9. [Claim] **$s=-1$.** No $\prod C=-1$ boundary has classical $K\to1$, so classical data do not support assigning $s=-1$ to any of them.
10. [Counting] **The missing patterns.** $a_3$: 2 of 2, $a_4$: 5 of 5, and $a_5$: 21 of 21 have BP vacua. Of the 21 in $a_5$, 12 have diagonal $K$; this count is by row counting, "unverified by direct enumeration".
11. [Text] **The proposed book text** for `tbl-missing-vacua`, `sec-boundary-missing-charges`, the last row of `tbl-boundary-types`, and `sec-real-coupling-vacua`.

## Coordinator notes (suspected weak points)

- **Whether a singular infimum counts as a vacuum is the central interpretive question.** The worker says the field diverges at $x=0$ and the bulk and boundary energies diverge separately, cancelling only in the limit. Decide whether the proposed book text states this honestly. Calling such a configuration "the vacuum that DG's amplitudes presuppose" (claim 5) is a strong reading: is perturbation theory about it, and the linearized $K$ of claim 5, well defined? Check what BP themselves say about these singular backgrounds.
- **Claim 4.** Rederive $E_P=E_Q+4\sum(e^{\psi_j(0)/2}-1)$ from the book's boundary potential and its normalization $\mathcal B(0)=0$ (`eq-boundary-energy`, `tbl-boundary-dictionary`). Check that BP's bound $E_Q\ge0$ applies in the book's normalization, and that "saturated but not attained" is correct.
- **Claim 1.** Check the energy conversion formula against BP's text and against `eq-vacuum-energy` for one solitonic and one singular case. Check that cyclic rotation of the extended Dynkin diagram is a symmetry at real coupling (unlike the coweight shift), and how the book's `tbl-missing-vacua` counts patterns (individually or by orbit); the counts in claim 10 must use the book's convention.
- **Claim 6.** Write your own linearized-reflection solver, or a WKB-free shooting code, for at least one pattern each of $a_3$ and $a_4$, and compare with DG (5.7) at $B\to0$. Do not rely only on `kcheck.py`.
- **Claim 5, the singular limit.** The worker argues that the Hessian terms that differ between the alternating pattern and $C\equiv-1$ vanish as $d\to1$. Check which nodes have $\tau_j(0)\to0$, and that $e^{\psi_j(0)/2}\to0$ on exactly the nodes where the two patterns differ.
- **Claim 7** is marked indicative by the worker. Recommend a label no stronger than that, unless you can confirm the non-diagonal $K$ with a more stable computation.
- **Claim 10.** Enumerate the $a_5$ patterns directly; do not count rows.
- **Prior art.** Perkins's Durham thesis (etheses.durham.ac.uk/4560) may already contain claims 5–7 (the linearized $K$ about BP backgrounds). Check before the book credits anything as new.
