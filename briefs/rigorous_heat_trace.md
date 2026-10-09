# Brief: make Proposition `prp-net-count` part 2 and Proposition `prp-one-loop-mass` rigorous

You are a mathematical physicist with expertise in spectral theory of non-self-adjoint Schrödinger operators. The repository is `/home/gustav/Git/affine_toda`. Do not edit any file in it. Work in your own scratch directory and deliver your full report in your final message: subagents may be unable to write report files.

## The task

In Part I (`foundations/*.qmd`), two results that the one-loop soliton masses rest on are labelled [Proposition, sketch]:

- **Proposition `prp-net-count` part 2** (Section `sec-counting-rule-transmission-factors`): the heat-trace formula
  $$\operatorname{Tr}\big(e^{-t\mathcal A}-e^{-t\mathcal A_0}\big)=\sum_\mu\Big[\sum_\zeta\operatorname{ord}_\zeta X_\mu\,e^{-t(\mu^2+\zeta^2)}+\frac1{\pi i}\int_0^\infty e^{-t(\mu^2+k^2)}\,\partial_k\log X_\mu(k)\,dk\Big],$$
  with $\mu$ running over the distinct masses and $X_\mu=\prod_{c:\,m_c=\mu}X_c$.
- **Proposition `prp-one-loop-mass`** (Section `sec-one-loop-masses-folded`): the closed form of the one-loop mass (`eq-one-loop-mass`),
  $$\Delta M=\tfrac12\sum_bm_b\sum_w\operatorname{ord}_w(X_b)\,\psi(\arcsin w),\qquad\psi(t)=\frac{t\cos t-\sin t}{\pi}.$$

The goal is to prove both as theorems, under hypotheses that hold for every static Hirota soliton treated in the paper. If a complete proof is out of reach, isolate exactly what is missing.

## What is available

Read Sections `sec-lagrangians` and `sec-counting-rule` to `sec-one-loop-masses-folded` of Part I (`foundations/*.qmd`) in full. In particular:

- $\mathcal A=-\partial_x^2+M(x)$ on $L^2(\mathbb R,\mathbb C^r)$, with $M$ complex symmetric and converging exponentially to $M_0$ at both ends.
- `thm-hirota-factorization`: under the Hirota form (H), $\det(I+K(\lambda))=\prod_bX_b(\lambda)$ on the physical sheet.
- Proposition `prp-net-count` part 1 (proved): $d'/d=\operatorname{tr}(R_0-R)$.
- `lem-wronskian` (corrected version): $X_b(k)X_{\bar b}(-k)=1$, so threshold zeros of single channels cancel within a mass level.
- Proposition `prp-hirota-form`: (H) holds for the static solitons of $a_n^{(1)}$, $d_n^{(1)}$, $e_n^{(1)}$ and their foldings given a degree bound. The $X_{ab}$ are rational in $z=\kappa_b/m_b$, with all zeros and poles at real $z=\cos(\pi x/h)\in[-1,1]$, so there are no zeros or poles at real $k\ne0$ (part 4).

A natural route for the heat trace:

1. $\mathcal A$ is sectorial because $M$ is bounded, so $e^{-t\mathcal A}$ has a Dunford integral over a contour around the spectrum.
2. $R-R_0=-R_0V_1(I+K)^{-1}V_2R_0$ is trace class, so the trace of the difference is $\frac1{2\pi i}\oint e^{-t\lambda}\,\partial_\lambda\log d\,d\lambda$.
3. Split $\partial\log d$ by (`eq-net-factorization`) into mass levels. Deform each level's contour onto its own cut $[\mu^2,\infty)$ and around its zeros and poles on the physical sheet.
4. Control the large-$\lambda$ behaviour ($X_\mu\to1$ for rational $X$) and the threshold, where grouping by level removes the $1/k$ singularity.

For the mass formula: insert the Mellin representations of $\sqrt\lambda$ and $\lambda^{-1/2}$, justify exchanging the $t$- and $k$-integrals, fix the counterterm channel by channel with the Born approximation $\log X_b=-ic_b/2k+O(k^{-2})$, and pair $b$ with $\bar b$.

## Known traps (found by a referee)

- The closed form is valid only after pairing $b$ with $\bar b$. Unpaired, a root $w$ contributes $\tfrac12m_b\big[\psi(\arcsin w)+\tfrac12\sqrt{1-w^2}\big]$, and the second term cancels in pairs.
- Roots at $w=\pm1$ occur (the translation mode is at $w=1$), so the hypothesis must allow $[-1,1]$, not $(-1,1)$.
- In Proposition `prp-one-loop-mass`, $w$ runs over zeros and poles in both half-planes. Reconcile this carefully with the heat-trace formula, which has only roots with $\operatorname{Im}\zeta<0$ on the physical sheet.
- The "Consequence for (`eq-one-loop-mass`)" paragraph of Section `sec-counting-rule-transmission-factors` gives only the discrete part; the continuum also depends on the roots.
- Embedded eigenvalues and threshold resonances need no special treatment in the formula. Each level's contour wraps its own cut, which is why the trace weight is the net order (Section `sec-hirota-form-dn1-en1`, threshold resonances, and Proposition `prp-cn-threshold-blocks`).
- The naive energy-weighted phase-shift sum misses a term $c_b/4\pi$ per channel; for the sine-Gordon kink it gives 0 instead of $-m_1/\pi$.

## Required checks

Verify the final formulas independently against:

- the sine-Gordon kink: $X=(z-1)/(z+1)$; the heat trace is $\operatorname{erf}\sqrt t$ for $m_1=1$, and $\Delta M=-m_1/\pi$;
- Hollowood's $a_n^{(1)}$ masses, $\Delta M_a=-m_a\big[\tfrac h{2\pi}-\tfrac14\cot\tfrac\pi h\big]$;
- the $a_3^{(1)}$ soliton $a=1$, whose heat trace is $(1+e^{-2t})\operatorname{erf}\sqrt{2t}$;
- the shipped scripts `foundations_code/oneloop/massformula.py`, `heatcheck.py`, `calib_an_cn.py` and `foundations_code/hirota/heat_a3.py`. Copy them to your scratch directory before running.

## Deliverable

In your final message:

- full proofs in the paper's style, with [Theorem] only where the proof is complete;
- exact hypotheses, and which solitons of the paper satisfy them;
- proposed replacement text for Proposition `prp-net-count` part 2 and Proposition `prp-one-loop-mass`, with the sketches replaced by proofs (short declarative sentences, LaTeX in `$...$`, numbered proof steps);
- the numerical checks with their precision;
- anything that remains open, stated precisely.

Usage is limited: work efficiently, and prefer a complete proof of the heat-trace formula over a partial treatment of both.
