# Referee report: exchange statistics and Bloch sectors (PR #2)

Brief: `briefs/referee_exchange_statistics.md`. Date: 2026-10-10. One independent referee agent, after `main` (PR #4, the physical sign of $q$) had been merged into the branch. The referee edited no repository file. Its scripts and outputs are in `referee_exchange_statistics_scripts/`.

## Integration (2026-10-10)

The findings below were integrated into PR #2. Changes beyond the referee's replacement text:

- Findings 1 and 2 were re-verified independently before integration:
  - the rank-3 projectors are degenerate at $\omega=0.3$, $0.7$ (weight violation $0.15$) and clean at $0.3013$;
  - $\check R(x;-q)=(\Pi\otimes1)\check R(x;q)(\Pi\otimes1)$ holds to $10^{-14}$;
  - the per-sector phase $c_Q=(-1)^{nN}e^{i\pi\gamma\cdot Q}$ holds in every sector checked.
- Code:
  - `sectors/lib.py` now has `check_projectors`, which refuses degenerate projectors;
  - rank-3 runs moved from $\omega=0.3$ to $0.3013$;
  - new scripts `bloch/t15_qsign_identity.py` (the operator identities) and the reruns of `t5`, `t6`, `t14`, `window.py`, `graded.py 3` and `twist.py 3 c`.
- The $\theta_c$ range in chapter 32 (main, PR #4) is corrected to $0.003$–$0.22$.
- Twist sign: `prp-bloch-sectors` now labels the sector by $U_\lambda=e^{-i\alpha\cdot\lambda}$, so its Bethe–Yang factor is $e^{i\alpha\cdot\mu_j}=\Omega_\alpha$, as in the rest of the section.
- Bibliography:
  - The Klassen–Melzer title is "Sine-Gordon $\neq$ massive Thirring, and related heresies"; the report had "=".
  - Göhmann–Murakami is cited only for introducing Fermi operators into lattice Yang–Baxter solutions, which is what its abstract states.
  - All three entries were checked against arXiv.

Not done (open for a later pass):
- brief items 6 (beyond the $\Theta$ remark), 8, 9 (beyond the labels touched) and 10;
- the $\pi/N$ thresholds by sector for $N\le5$ with targeted sampling (the referee suspects the $N=5$ case may need it, as for the graded rule);
- the $\psi=\pi$ exception for $N\ge4$;
- the parity and coproduct-ordering variants in item 3;
- whether $a_2^{(2)}$ admits a fermionic global form.

## Verdicts

| Result | Verdict |
|---|---|
| `prp-exchange-twist` | **Keep [Theorem].** The first sentence and its proof are correct. **Fix the second sentence**: its hypothesis fails for odd-rank spinors (finding 3). Credit prior art. |
| `prp-bloch-sectors` | **Keep [Proposition, sketch]**, with small fixes: the sign of the twist, twisted boundary conditions, the spinor torus, and $\Theta$ flipping $Q$. |
| Convention argument | **The general claim is correct**, derived and checked numerically. **Fix "a change of exactly this type"** for the sign of $q$: the mechanism is a different one. |
| Sign of $q$ as a Bloch twist | **Upgrade.** It follows exactly from a one-line lemma, checked numerically, with a sketch of a proof. **Fix "the two-soliton spinor spectra are unchanged"**: in some sectors every eigenvalue changes sign. |
| Scope: "every global form" | **Narrow it** to "every bosonic global form". |
| Row: $a_2^{(2)}$ kinks at $\alpha=0$ | Keep, with "every bosonic global form". |
| Row: graded rule | **Fix.** Confirmed at $N=4$; at $N=5$ the onset found is about $0.41\pi$, not $2\pi/5$. |
| Row: two kinks, general $\alpha$ | **Keep [Numerical]**; confirmed at 14 more couplings. |
| Row: spinor Bloch torus | **Fix.** "Nowhere on the grid for $c_3$ at $\omega=0.3$" is a numerical artefact. |
| Krein-pairing argument | Not refereed (item 6); $\Theta$ maps $(Q,\alpha)$ to $(-Q,\alpha)$. |

## Findings, most severe first

**1. The rank-3 spinor R-matrices are wrong at $\omega=0.3$ and $0.7$.** There $q^{40}=1$ (and $q^{20}=1$ for the sign used), and the projector construction in `spin_common.make` breaks down for $n=3$ ($c_3^{(1)}$, $a_5^{(2)}$).
- At $\omega=0.3$: weight-conservation violation $0.15$; Yang–Baxter residual $7.4$ ($c_3$) and $9.5$ ($a_5^{(2)}$); $|\check R(x)\check R(1/x)-1|=4.3$ and $5.2$.
- At $\omega=0.7$: the same checks give $0.05$, $3.4$ and $0.9$.
- Clean (Yang–Baxter $\le6\times10^{-14}$) at $\omega=0.29,0.3013,0.31,0.45,0.73,1.2,1.61,2.37$. $c_2$ is clean everywhere checked. Script: `r7e.py`.
- With a correct R near $\omega=0.3$, the two-soliton sector $Q=0$ of $c_3$ is unbroken at $330/1728$ ($\omega=0.29$), $342/1728$ ($0.3013$, $0.31$) and $330/1728$ ($0.73$) grid points (12 steps, threshold $10^{-6}$; `r5_c3_near03.py`). So "nowhere on the grid for $c_3^{(1)}$ at $\omega=0.3$" was an artefact.
- Also affected: the `t14` entries at $\omega=0.3$ for rank 3, the $n=3$ lines of `t5` at $0.3$, and on `main` the $\omega=0.3$ windows of $c_3$ and $a_5^{(2)}$ ($l_c=0.001$) and the $\omega=0.3$ line of `graded_3_qm`.
- "0 of 512 bicharacter twists" survives: $0/512$ with $\omega\in\{2.37,1.61,0.3013\}$ at the physical $q$.
- Suggested: a Yang–Baxter and weight-conservation guard in the setup, and $\omega$ offset away from roots of unity.

**2. The sign of $q$: the mechanism is misstated, and "spectra unchanged" is false.**
- $\check R(x;-q)=(\Pi\otimes1)\check R(x;q)(\Pi\otimes1)$, $\Pi=(-1)^g$, $g$ the parity of the number of minus signs. The $\pm1$ matrix $L$ of `qsign_checks.py` equals $\Pi\otimes1$ times a function of the total weight. Residual 0 for $n=2$ at $\omega=0.3,0.7,1.61,2.37$; $\le2\times10^{-15}$ for $n=3$ at clean couplings.
- $L$ is not gauge-equivalent to any $\mathbb Z_2$-valued bilinear form in the sign labels (exhaustive search). For $n=2$ it is equivalent to the complex bicharacter $e^{2\pi i\mu_{a,1}\mu_{b,2}}$.
- Proof sketch: $\Pi$ commutes with the long-root generators and anticommutes with the short-root and affine-node generators, and $k_i(-q)/k_i(q)=-1$ exactly on those nodes. So $(\Pi\otimes1)\Delta_q(\Pi\otimes1)=\Delta_{-q}$, and uniqueness of the intertwiner with $\check R(1)=1$ gives the lemma.
- Pairing consecutive collisions: $T_1(-q)=S\,(-1)^G\,\Pi_1^N\,T_1(q)\,S^{-1}$, $S$ the product of the parities of the even-numbered solitons. Residual $\le10^{-32}$ ($c_2$, $a_3^{(2)}$, $N=2$–$5$, four couplings) and $\le7\times10^{-15}$ ($c_3$, $a_5^{(2)}$, $N=2,3$). Script: `r7_qsign_exact.py`.
- $S$ commutes with every $\Omega_\alpha$, so $(Q,\alpha)\mapsto(Q,\alpha+\pi N\gamma)$ with phase $c_Q=(-1)^{nN}e^{i\pi\gamma\cdot Q}$, exactly constant within each charge sector.
- For $N=2$ the eigenvalues in sectors with odd $\gamma\cdot Q$ change sign (`r7f_c2N3.py`). Moduli, windows, the broken/unbroken pattern and the Bloch tori are unchanged; finite-volume levels are not.
- Three $c_2$ solitons at $\alpha=0$: unbroken at the physical $q$ for $\omega=0.3$ and broken at $-q$ ($0.76$); at $\omega=0.7$ the reverse ($0.75$). Confirmed.
- The general bicharacter statement is right: for $\check R\to F\check RF^{-1}$, $F=e^{i\pi B}$, $T_1'=e^{i\pi\Phi}e^{i\pi A(Q,\mu_1^{\rm out})}T_1e^{-i\pi\Phi}$ with $\Phi=\sum_{i<k}B(w_i,w_k)$, $A=B-B^{\mathsf T}$; residuals $2\times10^{-15}$–$10^{-13}$ for $N=2,3,4$ (`r3d_bichar.py`).
- $q\to1/q$ gives the complex-conjugate spectrum in every sector (`r3c_qinv.py`), so $(Q,\alpha)\mapsto(Q,-\alpha)$ and $(0,0)$ is fixed. Parity and the coproduct ordering were not checked.

**3. The second sentence of `prp-exchange-twist` fails for odd-rank spinors.**
- For weights $(\pm\frac12)^n$ with $\gamma=(1,\dots,1)$, $(-1)^g=c\,e^{i\pi\gamma\cdot\mu}$ with $c=e^{-i\pi n/2}$: $c=-1$ for $n=2$, $\pm i$ for $n=3$.
- For odd $n$, $\gamma\cdot\mu$ is a half-integer, so $\alpha=\pi(G-1)\gamma$ depends on $G$ as an integer, which is not fixed in a charge sector (example $n=3$, $N=3$, $Q=(\frac12,\frac12,\frac12)$: $G=1$ and $G=3$ both occur).
- The phase is $c^{[G-1]}$, common within a sector, not across sectors. It does not change moduli but shifts the momentum quantization. For sine-Gordon (all solitons odd, $\gamma=0$, $c=-1$) this phase is the whole sine-Gordon/massive Thirring difference.
- The proof is correct: an independent derivation matches line by line. Random-matrix check with an embedding built from graded adjacent swaps: $|T^{\rm gr}-DT|=0$ exactly, $|T^{\rm gr}-TD|=O(1)$ (`r1_random_even.py`).
- Minor: $g_j$ was used but not defined; only $T_1$ is treated.

**4. "Exists in every global form" is overstated.** A Ramond-type fermionic form (graded exchange, periodic spin structure) puts even $Q$ at $\alpha=\pi$ and has no bosonic $(Q=0,\alpha=0)$ sector; there the two-kink sector $Q=0$ is unbroken at every coupling. Whether $a_2^{(2)}$ admits a fermionisation is not settled. "The breaking appears there" is true but not exclusive (three kinks at $\psi=-0.891\pi$ have $Q=1$ broken at $\alpha=0$; two $c_3$ spinors have $Q=e_i$ broken). The Dirichlet-strip lattice fixes no global form.

**5. The graded rule at $N=5$.**
- $N=4$: unbroken at $0.4895\pi$, $0.4970\pi$; broken at $0.5020\pi$ ($3\times10^{-4}$), $0.5095\pi$, $0.5195\pi$. Threshold $2\pi/4$ confirmed to $\pm0.003\pi$.
- $N=5$: unbroken at $0.3897\pi$–$0.4094\pi$ ($\le8\times10^{-12}$; 500 random sets including clustered ones, a line scan with $|t|$ from $10^{-5}$ to $2.5$, Nelder–Mead continuation); broken at $0.4122\pi$ ($7.1\times10^{-3}$), $0.4147\pi$, $0.4172\pi$, $0.4197\pi$. The `t4` sampling (40 random sets) misses even the $0.42\pi$ breaking. Scripts: `r5_*`.

**6. `prp-bloch-sectors`: small points.** Covariant one-soliton states are fine. $\Theta$ preserves $\alpha$ but maps $Q\to-Q$. The twist sign differed between the propositions. Twisted boundary conditions $\vec\phi(L)=\vec\phi(0)+\frac{2\pi}\beta\vec\lambda$ give $Q=\vec\lambda$ with every $\alpha$. The spinor charges generate $\mathbb Z^n+\mathbb Z(\frac12,\dots,\frac12)$, so the torus is twice $[0,2\pi)^n$, and a shift by $2\pi e_i$ multiplies the twist by $-1$: the `t6` grid suffices. The Ramond/NS remark is correct.

**7. Minor.** `t3` uses a $10^{-8}$ threshold while the README said $10^{-6}$.

## Re-run numbers

Two-kink threshold (`r5_threshold.py`): $\psi/\pi=+0.417$, $-0.455$, $\pm0.55$, $\pm0.62$, $\pm0.70$, $\pm0.78$, $\pm0.88$, $\pm0.97$; $|\log x|$ from $10^{-5}$ to $80$, both signs; 22-step bisection. $\alpha_c$ agrees with $\max(0,\min(2|\psi|-\pi,\pi-|\psi|))$ to $\le6.2\times10^{-7}$ in $\alpha/\pi$.

## Prior art

- Klassen and Melzer, "Sine-Gordon $\neq$ massive Thirring, and related heresies", Int. J. Mod. Phys. A 8 (1993) 4131, hep-th/9206114.
- Feverati, Ravanini and Takács, "Scaling functions in the odd charge sector of sine-Gordon/massive Thirring theory", Phys. Lett. B 444 (1998) 442, hep-th/9807160.
- Göhmann and Murakami, "Fermionic representations of integrable lattice systems", J. Phys. A 31 (1998) 7729, cond-mat/9805129.
- Lieb, Schultz and Mattis, Ann. Phys. 16 (1961) 407 (parity-dependent Jordan–Wigner twist).
- Bajnok, Palla, Takács and Wágner, Nucl. Phys. B 587 (2000) 585, hep-th/0004181: Bloch waves $\vartheta_n=2\pi n/k$, NLIE twist, topological charge in $k\mathbb Z$, Bethe–Yang per sector against TCSA. The book's attribution is accurate for real-coupling sine-Gordon (read through a summary).

The general $N$-body telescoping identity was not found stated elsewhere.

## Doubts stated by the referee

- The $\check R(-q)$ lemma is numerical; the proof sketch assumes a normalisation of $k_i$ not checked against `lib.spinor_rep`.
- The $N=5$ graded onset could be below $0.41\pi$ if the broken region is below the sampling resolution.
- Whether a fermionic form exists for $a_2^{(2)}$ is not established.
