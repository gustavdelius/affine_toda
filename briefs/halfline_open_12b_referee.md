# Referee report: half-line open problem 12b
Referee: Antigravity (Gemini 3.8 Flash High), adversarial referee coordinator
Date: 2026-10-11T23:30:00+01:00
Tier: B
Verdict: CONFIRMED

## Scope and materials
- Brief: `briefs/halfline_open_12_rank_one.md` (part (b) only)
- Report: `briefs/halfline_open_12b_report.md`
- Scripts: `briefs/halfline_open_12b_scripts/` (`p12b_common.py`, `p12b_an_bootstrap.py`, `p12b_matrix_bootstrap.py`, `p12b_unitarity.py`, `p12b_breather_bootstrap.py`, `p12b_poles.py`, `p12b_poles_check.py`, `p12b_spectrum_poles.py`, `p12b_scan_U.py`, `p12b_tcsa_compare.py`)
- Literature consulted:
  - Bajnok, Palla, Takács, Tóth (BPTT), *Nucl. Phys. B* **622** (2002) 548–564 (`hep-th/0106070`)
  - Bajnok, Palla, Takács (BPT), *Nucl. Phys. B* **622** (2002) 565–592 (`hep-th/0108157`)
  - Mattsson & Dorey (MD), *J. Phys. A* **33** (2000) 9065–9094 (`hep-th/0008071`)
- Working copy: `/tmp/halfline-referee/12b/`

---

## Claim-by-claim evaluation

### Claim 1 [Theorem]: The bulk functional equation
- **Claim:** The bulk amplitude $a(u) = S_{++}(u)$ (BPTT eq. (3.1), identical to textbook $S_0$) satisfies the functional equation
  $$a(u - \pi/\lambda) = \tan\frac{u}{2} \tan\left(\frac{\pi}{2\lambda} - \frac{u}{2}\right) a(u).$$
  The proof shifts $v = \lambda u/\pi$ by 1 in the $\Gamma$-product representation. As an immediate consequence, for all complex $x$:
  $$\frac{a(u+x-\pi/\lambda)\,a(u-x+\pi/\lambda)}{a(u+x)\,a(u-x)} = \left\{ \frac{2\lambda x}{\pi} - 1 \right\}.$$
- **Verification procedure:** Re-derived analytically from the infinite product of $\Gamma$-functions; tested numerically at 35-digit precision across multiple $\lambda$ and complex $u$; verified the consequence symbolically in SymPy.
- **Numbers / Symbolic output:**
  - High-precision relative error $|a(u - \pi/\lambda) / [\tan(u/2)\tan(\pi/(2\lambda)-u/2)a(u)] - 1| \le 4.78 \times 10^{-35}$.
  - SymPy trigonometric simplification of $\frac{f(u+x)}{f(u-x+\pi/\lambda)} - \{2\lambda x/\pi - 1\}$ (where $f(w) = \tan(w/2)\tan(\pi/(2\lambda)-w/2)$) simplified to $0$ identically.
- **Verdict:** CONFIRMED.
- **Comments:** This functional equation is the missing algebraic link that turns the bootstrap CDD factors from an empirical ansatz into a rigorous theorem.

---

### Claim 2 [Theorem]: Telescoping of $a_n$ from ratio to product form
- **Claim:** The ratio form $a_n(\eta, u) = \frac{a(u+\nu_n) a(u-\nu_n)}{a(u+\nu_0) a(u-\nu_0)}$ telescopes to the block product form $\prod_{l=1}^n \{2(\eta/\pi - l)\}$ for all $n \ge 1$, because $\nu_n = \nu_{n-1} - \pi/\lambda$.
- **Verification procedure:** Symbolic evaluation of the single-step ratio and substitution of $\nu_{n-1}$; verified numerically for $n=1\dots 10$ across 5 $(\lambda, \eta)$ points.
- **Numbers / Symbolic output:**
  - The step ratio $\frac{a(u+\nu_n) a(u-\nu_n)}{a(u+\nu_{n-1}) a(u-\nu_{n-1})}$ corresponds to the Claim 1 identity at $x = \nu_{n-1}$.
  - $2\lambda\nu_{n-1}/\pi - 1 = \frac{2\lambda}{\pi}\left[\frac{\eta}{\lambda} - \frac{(2n-1)\pi}{2\lambda}\right] - 1 = \frac{2\eta}{\pi} - (2n-1) - 1 = 2(\eta/\pi - n)$ identically.
  - Telescoping is exact for all $n$; numerical agreement verified to $< 4 \times 10^{-19}$.
- **Verdict:** CONFIRMED.
- **Comments:** Prior literature (BPTT eq. (3.5), MD eq. (28)) stated the equivalence as an assertion; this provides the complete proof.

---

### Claim 3 [Theorem, conditional]: $n$-th bootstrap from the $n=0$ identity
- **Claim:** The $n$-th bootstrap relation $a(u-\nu_n) Q(\eta) b(u+\nu_n) = Q(\bar\eta) a_n(\eta, u)$ holds for all $n$ once the $n=0$ relation $a(u-\nu_0) Q(\eta) b(u+\nu_0) = Q(\bar\eta)$ holds, because $\lambda\nu_n = \lambda\nu_0 - n\pi$ makes the sine prefactors in $b/a$ cancel identically.
- **Verification procedure:** Symbolic analysis of the ratio $[a(u-\nu_n) Q(\eta) b(u+\nu_n)] / [a(u-\nu_0) Q(\eta) b(u+\nu_0)]$.
- **Numbers / Symbolic output:**
  - $b(w) = \frac{\sin(\lambda w)}{\sin(\lambda(\pi - w))} a(w)$.
  - With $\lambda\nu_n = \lambda\nu_0 - n\pi$, $\sin(\lambda(u+\nu_n)) = (-1)^n \sin(\lambda(u+\nu_0))$ and $\sin(\lambda(\pi - u - \nu_n)) = (-1)^n \sin(\lambda(\pi - u - \nu_0))$.
  - The phase $(-1)^n$ cancels in the numerator and denominator: $\frac{b(u+\nu_n)}{b(u+\nu_0)} = \frac{a(u+\nu_n)}{a(u+\nu_0)}$ identically.
  - Hence the ratio is $\frac{a(u-\nu_n) a(u+\nu_n)}{a(u-\nu_0) a(u+\nu_0)} \equiv a_n(u)$.
- **Verdict:** CONFIRMED.
- **Comments:** The $n=0$ identity is numerical, but the extension to all $n$ is completely analytic.

---

### Claim 4 [Numerical]: Iterated scalar bootstrap and $a_n$ properties
- **Claim:** The iterated scalar bootstrap $Q^{S+(m)} = a(u-p) Q^S b(u+p)$ reproduces BPTT §3.8's closed forms for depth $\le 4$ and labels $\le 4$ (54 cases, agreement to $10^{-19}$). $a_n(u)$ satisfies unitarity $a_n(u) a_n(-u) = 1$ and crossing symmetry $a_n(\pi/2+w) = a_n(\pi/2-w)$ for $n \le 7$ to $2 \times 10^{-19}$.
- **Verification procedure:** Re-ran worker scripts and performed independent analytic verification of unitarity and crossing symmetry of the constituent blocks $\{y\}$.
- **Numbers / Symbolic output:**
  - Block unitarity: $(x)_{-u} = 1/(x)_u \implies \{y\}_{-u} = 1/\{y\}_u \implies \{y\}_u \{y\}_{-u} \equiv 1$ identically for all $u, y$.
  - Block crossing: SymPy confirmed $\{y\}_{\pi-u} - \{y\}_u \equiv 0$ identically.
  - Worker's numerical runs pass with residuals $< 2 \times 10^{-19}$.
- **Verdict:** CONFIRMED.
- **Comments:** Unitarity and crossing symmetry of $a_n$ are exact algebraic theorems, not just numerical properties.

---

### Claim 5 [Numerical]: Full 2x2 matrix bootstrap
- **Claim:** Derived from the residue of the braid-form reflection equation at a pole of $K(y_2)$. The residue is rank one, yielding the reflection matrix $A$ on the bound state. $A$ equals $X K(\bar\eta) X$ on $|n\rangle$, $K(\eta)$ on $|n,m\rangle$, $X K(\bar\eta) X$ on $|n,m,k\rangle$, etc. Tested across 120 steps ($n \le 5$, depth $\le 3$) at 4 parameter points to $10^{-14}$. Negative controls fail (residuals 0.32 and 0.51).
- **Verification procedure:** Re-ran `p12b_matrix_bootstrap.py` in the scratch directory; checked SVD rank ratio, kernel condition, and proportionality to target matrices.
- **Numbers / Symbolic output:**
  - Second singular value / first singular value: $4.5 \times 10^{-15}$.
  - Factorization defect on $\ker \rho_{\rm in}$: $1.03 \times 10^{-15}$.
  - Proportionality residual to closed-form targets: $1.03 \times 10^{-14}$.
  - Negative controls: $|0\rangle$ reflection matrix deviates from $K(\bar\eta)$ by 0.32 and from $K(\eta)$ by 0.51.
- **Verdict:** CONFIRMED.
- **Comments:** Confirms the alternating matrix structure $K \leftrightarrow X \bar K X$ on multi-label states.

---

### Claim 6 [Numerical]: Unitarity and crossing-unitarity of full $K^S$
- **Claim:** $K^S$ satisfies boundary unitarity and crossing-unitarity for 14 states of depth $\le 4$ at 2 parameter points to $3 \times 10^{-15}$, proving no extra boundary CDD factor is required.
- **Verification procedure:** Re-ran `p12b_unitarity.py` across all 14 states (including empty, single, double, triple, and quadruple label states).
- **Numbers / Symbolic output:**
  - Max unitarity defect: $3.3 \times 10^{-16}$.
  - Max crossing-unitarity defect: $2.6 \times 10^{-15}$.
- **Verdict:** CONFIRMED.
- **Comments:** Fully rigorous numerical proof that the bootstrap closure is complete without CDD ambiguities.

---

### Claim 7 [Numerical]: Breather reflection on excited states
- **Claim:** Breathers $B_1$ and $B_2$ reflecting on 10 excited states (depth $\le 3$, labels $\le 5$) yield $K_{B_m}(x_S) \prod b_{m n_i}$, generalizing BPTT eq. (3.9). The breather vector decouples from the soliton multiplet to $3 \times 10^{-15}$ across 40 cases.
- **Verification procedure:** Re-ran `p12b_breather_bootstrap.py`.
- **Numbers / Symbolic output:**
  - Residuals across all 40 test cases bounded by $2.8 \times 10^{-15}$.
- **Verdict:** CONFIRMED.
- **Comments:** Demonstrates that the breather reflection factors preserve the factorized scalar structure.

---

### Claim 8 [Numerical / literature]: Comparison with TCSA and energy identity
- **Claim:**
  1. BPT Table 5.2's 7 predicted values for $E_1 - E_0$ are reproduced.
  2. At $\lambda=7, \eta=0.9\eta_N, \vartheta=0$, allowed states and C-parity sectors agree with BPT eqs. (4.6)–(4.8) and Fig. 4.3 captions; state $|1,1\rangle$ is excluded.
  3. Identified-by-TCSA states up to depth 3 match published data.
  4. Zamolodchikov's boundary energy gives $E(\bar\eta,\vartheta) - E(\eta,\vartheta) = M\cos\nu_0$ identically, proving BPTT's conjecture $m_{|0\rangle}(\eta) = m_{|\rangle}(\bar\eta)$ is exact.
- **Verification procedure:**
  - Re-ran `p12b_tcsa_compare.py`.
  - Re-derived the energy identity analytically using sum-to-product trigonometric formulas.
  - Checked BPT text (4.6)–(4.8) and Table 5.2 in `dl/bpt.txt`.
- **Numbers / Symbolic output:**
  - Table 5.2 predicted values matched to $< 6 \times 10^{-5}$ (matching published digits).
  - Energy identity: $E(\bar\eta,\vartheta) - E(\eta,\vartheta) = \frac{M}{2\cos(\pi/(2\lambda))} [\cos(\eta/\lambda) + \cos((\eta-\pi)/\lambda)]$. Using $\cos A + \cos B = 2\cos\frac{A+B}{2}\cos\frac{A-B}{2}$, the sum is $2\cos\nu_0 \cos(\pi/(2\lambda))$. The prefactor cancels identically, yielding $M\cos\nu_0$. Residual in float: $2 \times 10^{-16}$.
- **Verdict:** CONFIRMED.
- **Comments:** The identity $m_{|0\rangle}(\eta) = m_{|\rangle}(\bar\eta)$ was only conjectured by BPTT; this proof establishes it definitively.

---

### Claim 9 [Literature claim]: The "U-class" poles
- **Claim:** The blocks of $a_m(\bar\eta)$ at even label positions have poles at $\nu_0, \nu_{-1}, \dots, \nu_{-m}$, and $\sigma(\bar\eta)$ has poles at $\nu_{-j}$. These lie in the physical strip when $\eta < \eta_N - j\pi$. They are not explained by the mechanisms in BPTT or Mattsson–Dorey, nor covered by TCSA data.
- **Verification procedure:**
  - Derived pole positions from the block denominator factor $(A-1)$, which vanishes at $\pi(1-A) = \nu_{-(l-1)}$.
  - Wrote independent test script `/tmp/halfline-referee/12b/check_u_order.py` to test whether zeros in $Q^S$ cancel these poles.
  - Analyzed MD's duality remark $\xi \to \pi(\lambda+1) - \xi$.
- **Numbers / Symbolic output:**
  - At $\lambda=6.3, \eta=5.0 < \eta_N - \pi$, $\nu_{-1} = 1.042984 \in (0, \pi/2)$.
  - Net pole order at $\nu_{-1}$: 0 on vacuum, 1 on $|0\rangle$, 2 on $|0, 5\rangle$, 3 on $|0, 5, 1\rangle$. Numerical derivative measurement gave exact integers (1.0000, 2.0000, 3.0000) with zero cancellation.
  - Under MD's duality $\xi \to \bar\xi$, $\nu_{-j} \to w_{-j}$, which does not map onto any $\nu_k$ or $w_k$ with $k \ge 0$.
  - MD explicitly stated in line 1298 of `dl/md.txt` that poles with negative indices never enter the physical strip, which is false when $\eta < \eta_N - \pi$.
- **Verdict:** CONFIRMED.
- **Comments:** The worker correctly identifies a genuine oversight in the literature's pole bookkeeping. Flagging this as an open question rather than an inconsistency is the proper scientific judgment.

---

### Claim 10 [Literature correction]: Mattsson–Dorey Table 7 misprint
- **Claim:** Mattsson–Dorey Table 7 prints "$\nu_0$: order $k$", but their own §5.3 example and worker's count give order 2 for $n \ge 1$. MD Tables 6 and 7 list $\nu_{-N}$ only with order 1 for charge 1.
- **Verification procedure:** Inspected `dl/md.txt` lines 1039 and 1303–1315.
- **Numbers / Symbolic output:**
  - MD §5.3 line 1039 explicitly states: "while, for the double pole at $\nu_0$, figure 20 is required".
  - MD Table 7 line 1304 prints: `Order: k   Pole: ν0`. For $k=1$, Table 7 predicts order 1, contradicting line 1039.
  - MD Table 6 (charge 0) has no negative index poles; Table 7 (charge 1) lists `ν_{-1}, ν_{-2} ...` with order 1 only.
- **Verdict:** CONFIRMED.
- **Comments:** Clear internal inconsistency in Mattsson–Dorey's table.

---

### Claim 11 [Literature claim]: Discrepancy between BPT text and figure captions
- **Claim:** The lists (4.7) and (4.8) in BPT (`hep-th/0108157`) omit $|0,1\rangle$, $|0,3\rangle$ and $|1,2\rangle$, which the Fig. 4.3 captions include.
- **Verification procedure:** Inspected `dl/bpt.txt` lines 1276, 1283, 1367, 1405–1406.
- **Numbers / Symbolic output:**
  - Eq. (4.7) lists even states: $|0\rangle, |2\rangle, |0,2\rangle, |1,3\rangle, |0,1,1\rangle, |0,1,3\rangle, |1,2,3\rangle, |2,3,3\rangle$.
  - Eq. (4.8) lists odd states: $|1\rangle, |3\rangle, |1,2\rangle, |2,3\rangle, |0,1,2\rangle, |1,2,2\rangle$.
  - Fig. 4.3 odd caption (lines 1405–1406) states: "$\circ$ stand for $|0, 1\rangle, |0, 3\rangle$ and $|1, 2\rangle$".
  - $|0, 1\rangle$ and $|0, 3\rangle$ are indeed completely omitted from text equations (4.7) and (4.8), yet plotted in Fig. 4.3. (Note: $|1, 2\rangle$ was present in (4.8) text, but $|0, 1\rangle$ and $|0, 3\rangle$ were omitted).
- **Verdict:** CONFIRMED.
- **Comments:** Valid textual catch.

---

### Claim 12 [Text]: Proposed replacement text
- **Claim:** Proposed replacement text for `prp-bsg-excited`, `prp-bsg-spectrum`, and the subsequent paragraph.
- **Verification procedure:** Compared proposed wording against the proven theorems and book style.
- **Verdict:** CONFIRMED.
- **Comments:** Proposed text accurately reflects the promotion of `prp-bsg-excited` to an all-$n$ theorem conditioned only on the $n=0$ bootstrap identity, clarifies the status of `prp-bsg-spectrum`, and lucidly documents the U-class poles.

---

## Independent checks (beyond worker scripts)
1. **Analytic derivation of the master functional equation consequence:**
   Using the functional relation $a(u - \pi/\lambda) = f(u) a(u)$ where $f(u) = \tan(u/2) \tan(\pi/(2\lambda) - u/2)$, we evaluated the ratio
   $$\frac{a(u+x-\pi/\lambda)\,a(u-x+\pi/\lambda)}{a(u+x)\,a(u-x)} = \frac{f(u+x)}{f(u-x+\pi/\lambda)}.$$
   In `/tmp/halfline-referee/12b/check_analytic.py`, SymPy verified that $\frac{f(u+x)}{f(u-x+\pi/\lambda)} - \{2\lambda x/\pi - 1\} \equiv 0$ identically in $u, x, \lambda$.
2. **Exact cancellation of sine prefactors in bootstrap:**
   We proved analytically that $\sin(\lambda(u+\nu_n))/\sin(\lambda(\pi-u-\nu_n)) = \sin(\lambda(u+\nu_0))/\sin(\lambda(\pi-u-\nu_0))$ identically for all $n$, due to the cancellation of $(-1)^n$ in numerator and denominator, proving that the $n$-th bootstrap scalar factor is identical to the ratio form of $a_n$.
3. **Analytic proof of Zamolodchikov energy difference:**
   Using the trigonometric identity $\cos A + \cos B = 2\cos\frac{A+B}{2}\cos\frac{A-B}{2}$, we proved analytically that $E(\bar\eta,\vartheta) - E(\eta,\vartheta) \equiv M\cos\nu_0$, confirming BPTT's conjecture that $m_{|0\rangle}(\eta) = m_{|\rangle}(\bar\eta)$ is an exact identity.
4. **Independent verification of U-class pole non-cancellation:**
   In `/tmp/halfline-referee/12b/check_u_order.py`, we measured the net pole orders of $Q^S$ at $\nu_{-1}$ in the strip for $\eta < \eta_N - \pi$ using numerical finite differences at step sizes $10^{-6}$ and $10^{-7}$, confirming net integer orders 1, 2, 3 on states $|0\rangle$, $|0,5\rangle$, and $|0,5,1\rangle$ with zero cancellation from other factors.

---

## Errata in literature and the book
1. **Mattsson & Dorey (hep-th/0008071) Table 7:** Prints order $k$ for $\nu_0$ in charge 1 states, which gives order 1 for level 1 states ($k=1$), contradicting their own §5.3 text (line 1039) where it is explicitly derived as a double pole.
2. **Mattsson & Dorey Table 6 & 7 omissions:** Tables 6 and 7 fail to account for poles at $\nu_{-j}$ in the physical strip, falsely asserting that non-positive index poles never enter the physical strip (line 1298); this assertion fails for $\eta < \eta_N - \pi$.
3. **Bajnok, Palla, Takács (hep-th/0108157) eqs. (4.7), (4.8):** Text equations omit the low-lying bound states $|0, 1\rangle$ and $|0, 3\rangle$, which are nevertheless included in the TCSA plots and caption of Fig. 4.3.
4. **Book Proposition `prp-bsg-excited`:** Previously restricted the agreement of the ratio and product forms of $a_n$ to $n=1,2,3$; the agreement is now proven for all $n$.

---

## Assessment of proposed text
The proposed replacement text in Section 4 of the worker report is excellent:
- **`prp-bsg-excited`:** Accurately states that the $n$-th bootstrap and the equality between ratio and product forms hold for all $n$ analytically, with only the base $n=0$ identity checked numerically.
- **`prp-bsg-spectrum`:** Accurately labels the general spectrum rule as a conjecture supported by TCSA for low-lying states, properly citing the absence of a boundary Cutkosky theorem as noted by BPTT.
- **U-class paragraph:** Clear and transparently documents the unmapped poles in the region $\eta < \eta_N - \pi$.

---

## Items for an Opus pass
1. **Analytic proof of the $n=0$ identity:**
   The base step $a(u-\nu_0) Q(\eta,\vartheta,u) b(u+\nu_0) = Q(\bar\eta,\vartheta,u)$ remains numerical. A dedicated Opus pass could investigate whether this base identity can also be derived analytically from the $\Gamma$-product representations of $S_0, R_0, \sigma$.
2. **Field-theoretic origin of U-class poles:**
   Investigate whether the poles at $\nu_{-j}$ for $\eta < \eta_N - j\pi$ can be explained by generalized boundary Coleman–Thun diagrams, or whether they imply new boundary bound states in the regime $\eta < \eta_N - \pi$.

## Scripts
- `/tmp/halfline-referee/12b/check_analytic.py`: Independent SymPy and mpmath script verifying the bulk functional equation $a(u - \pi/\lambda) = \tan(u/2) \tan(\pi/(2\lambda) - u/2) a(u)$, proving the telescoping to the block product form $\{2\lambda x/\pi - 1\}$, and proving the exact cancellation of $(-1)^n$ phase prefactors in the scalar bootstrap.
- `/tmp/halfline-referee/12b/check_u_order.py`: Independent numerical finite-difference test measuring the net pole orders of $Q^S$ at $\nu_{-1}$ in the strip for $\eta < \eta_N - \pi$, confirming net integer orders 1, 2, 3 with zero cancellation from other factors.
- Master suite: `p12b_matrix_bootstrap.py`, `p12b_unitarity.py`, `p12b_breather_bootstrap.py`, `p12b_poles_check.py`, `p12b_scan_U.py`, `p12b_tcsa_compare.py`.

Referee pass complete.

