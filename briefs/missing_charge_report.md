# Research and verification report: `prp-static-saddles` and `cnj-missing-charge`

Date: 2026-10-27. Model: Gemini. Commit reviewed: `4ae0994c57c729d2f214e9c4fa2198ddfd45f180`.

Reviewed files:
- `part5-semiclassics/17-thimbles.qmd` (lines 116–138, `sec-which-saddles-contribute`, `prp-static-saddles`, `cnj-missing-charge`)
- `part3-classical/10-solitons.qmd` and `part3-classical/11-reality.qmd` (`sec-topological-charges`, `sec-missing-charges`)
- `part9-foundations/44-twelve-gaps.qmd` (gap H11, line 90) and `part9-foundations/45-lattice.qmd` (`thm-coulomb-gas`)
- `outlook/50-outlook.qmd` (open problem 4, line 32)
- `appendices/E-original-results.qmd` (lines 12–13)
- `textbook_code/thimbles/` (`lie.py`, `hirota.py`, `model.py`, `newton_batch.py`, `flows.py`, `flows2.py`, `charges.py`, `charges_multi.py`)

Accompanying scripts in `briefs/missing_charge_scripts/`:
- `verify_missing_charge.py` (consolidated 6-check `PASS`/`FAIL` script ready to be adapted into `textbook_code/thimbles/`)
- `test_a3_thimbles.py` (exact Hirota/PDE residuals for the fusing-angle solution, exact Coulomb-gas $Z_1(\vec\lambda)$ for all 6 weights of $\Lambda^2$ in $a_3^{(1)}$, and $L=1$ critical points)
- `test_continuum_12.py` (exact $L=1$ intersection numbers $n_\sigma$ and 1-loop/2-loop saddle expansion in the missing sector $(1,2)$, plus continuum scaling)
- `test_v12_bvp.py` ($L=1, 2, 3$ critical points in sector $(1,2)$ and continuum BVP test on $\mathcal{V}_{(1,2)}$)
- `test_shooting_12.py` (1D shooting proof that no smooth symmetric continuum static solution exists on $\mathcal{V}_{(1,2)}$)
- `test_fusing_action.py` (spatial integrals of $(\partial_x\vec u)^2/2$, $(\partial_t\vec u)^2/2$, $W(\vec u)$, $E_{\text{Mink}}$, and $L_{\text{Eucl}}$ across chambers and as $\epsilon \to 0^+$)

---

## Executive summary

| Item | Target | Verdict / Status | Summary of findings |
|---|---|---|---|
| **1** | `prp-static-saddles` (`17-thimbles.qmd` L120–124) | **Correct** (with one precision fix on "non-resonant"); upgrade **`[Proposition, sketch]` $\to$ `[Theorem]`** | Complete analytical proof via Poincaré–Dulac non-resonance ($\sum_{a=1}^3 n_a m_a \notin \{m_1,m_2,m_3\}$ for $\sum n_a \ge 2$), exact static $\{1,2,3\}$ Hirota family ($A_{12}(0)=A_{23}(0)=(\sqrt{2}-1)^2, A_{13}(0)=1/2$), $\mathbb{Z}_4$ N-ality selection rule $\sum_{a\in S} a \equiv 2 \pmod 4$, and additivity of topological charge for spatially separated kinks ($x_1 \ll x_2 \ll x_3$). Removes the first `[Numerical]` tag in `prp-static-saddles`. |
| **2** | Fusing-angle bound states in `prp-static-saddles` (`17-thimbles.qmd` L122) | **Correct**; upgrade second **`[Numerical]` $\to$ `[Theorem]`** | The asymmetric fusing limit $e^{\xi_1} \to C_1 \in \mathbb{C}$, $e^{\xi_2} \sim 2(u-\pi/4)C_4/C_1 \to 0$ of two species-1 solitons at $\theta_{1,2} = \pm i u$ ($u \to \pi/4$) gives the exact $(1+1)$-dimensional Hirota solutions $\tau_j^\pm(x,t) = 1 + i^j C_1 e^{m(x\mp it)} + (-1)^j C_4 e^{2mx}$. Proved analytically that $\tau_j^\pm$ has exact rest energy $E = M_2$ and momentum $P = 0$, covers all 6 weights of $\Lambda^2\mathbb{C}^4$ in closed form, and for $C_4=-1, C_1=i\epsilon$ ($0<\epsilon<2$) has strictly real Euclidean Lagrangian $L_{\text{Eucl}} = M_2 - O(\epsilon^2) \to M_2$ as $\epsilon \to 0^+$ in the missing sector $(1,2)$. Generalizes to all $a_n^{(1)}$ and $\Lambda^a\mathbb{C}^h$. |
| **3** | Finite-lattice part of `cnj-missing-charge` (`17-thimbles.qmd` L128–131) | **New `[Theorem]`** for any finite $(1+1)$D lattice $\Lambda$ | On any $(1+1)$D lattice $\Lambda = \mathbb{Z}_{N_t} \times \{1,\dots,L\}$, the combined antilinear spatial-reflection symmetry $\Theta \circ R_x$ makes $Z_\Lambda(\vec\lambda)$ strictly real for every weight $\vec\lambda \in P(g)$, diagram automorphisms $D_h$ make $Z_\Lambda(\vec\lambda)$ equal across each orbit ($\{(1,2),(2,3),(3,4),(1,4)\}$ and $\{(1,3),(2,4)\}$), and the Coulomb-gas expansion proves $Z_\Lambda(\vec\lambda) \ne 0$ for generic $(\kappa,\hbar)$. By Picard–Lefschetz completeness $Z_\Lambda(\vec\lambda) = \sum_\sigma n_\sigma Z_\sigma$, every missing-charge sector has at least one complex saddle $\sigma$ with $n_\sigma \ne 0$. |
| **4** | Reduced 1D model in missing sector $(1,2)$ (`17-thimbles.qmd` L133–137) | **New `[Numerical]` check + analytical continuum-limit theorem** | At $L=1, \kappa=1$, $Z_1((1,2)) > Z_1((1,3)) > 0$, and sector $(1,2)$ has 54 critical points, of which 4 have $\|n_\sigma\| = 1$; their 2-loop sum reproduces $Z_1((1,2))$ across 35 orders of magnitude ($\hbar \in [0.15, 1.0]$, rel. err $1.7\times 10^{-3}$ at $\hbar=0.25$). However, as the lattice spacing $a \to 0$, the dominant 1D saddle in $\mathcal{V}_{(1,2)}$ is a **link-pinned lattice artifact** straddling a real zero of $\tau_1(x)=\tau_3(x)$, whose continuum action diverges as $S_{\text{cont}} \sim 4.93/a \to +\infty$. Turning on the $(1+1)$D fusing-angle mode $C_1 = i\epsilon$ splits this zero off the real axis to $x = \pm i\epsilon/(2m)$ with finite continuum action $M_2$. |
| **5** | Full continuum limit of `cnj-missing-charge` | **Remains `[Conjecture]`** | Identifying the finite-lattice saddles with $n_\sigma \ne 0$ on the $(1+1)$D lattice $\Lambda = \mathbb{Z}_{N_t} \times \{1,\dots,L\}$ as $a \to 0, N_t a \to \infty$ with the continuum fusing-angle bound states $\tau_j^\pm$, and matching their relative quantum weights across $D_h$-orbits to $U_q(\hat g^\vee)$, remains open. |

---

## Item 1: Complete analytical proof of `prp-static-saddles`

**Verdict: Correct statement; upgrade `**[Proposition, sketch]**` (and the first internal `[Numerical]` tag) to `**[Theorem]**`, with a precision fix on "non-resonant".**

### 1.1 Precision point on "In $a_3^{(1)}$ the masses are non-resonant"

In `part5-semiclassics/17-thimbles.qmd` line 122, `prp-static-saddles` currently says:
> "In $a_3^{(1)}$ the masses are non-resonant. So every complex static solution leaving the vacuum lies in the unstable manifold of the vacuum, which is parametrized by the coefficients of the $e^{m_ax}$ and is realized by static Hirota multi-solitons of distinct species."

In $a_3^{(1)}$ ($h=4$), the soliton inverse widths (in units where $m = \sqrt{\kappa} = 1$) are $m_a = 2\sin(\pi a/4)$:
$$
m_1 = \sqrt{2}, \qquad m_2 = 2, \qquad m_3 = \sqrt{2}.
$$
Thus $m_1 = m_3$ ($1:1$ equality between conjugate species $1$ and $3$), and $m_1 + m_3 = 2\sqrt{2}$. What is needed for the Poincaré–Dulac theorem on the unstable manifold $\mathcal{W}^u(0)$ of the vacuum $\vec u = 0$ as $x \to -\infty$ is that **no non-negative integer linear combination $\sum_{a=1}^3 n_a m_a$ of degree $n_1 + n_2 + n_3 \ge 2$ equals any $m_b \in \{\sqrt{2}, 2\}$**. Because $\sqrt{2}$ is irrational:
- $(n_1 + n_3)\sqrt{2} + 2 n_2 = \sqrt{2} \implies n_2 = 0$ and $n_1 + n_3 = 1$, contradicting $n_1 + n_2 + n_3 \ge 2$;
- $(n_1 + n_3)\sqrt{2} + 2 n_2 = 2 \implies n_1 + n_3 = 0$ and $n_2 = 1$, contradicting $n_1 + n_2 + n_3 \ge 2$.

Furthermore, in the Hirota bilinear equation $\tau_j \tau_j'' - (\tau_j')^2 = m^2(\tau_j^2 - \tau_{j-1}\tau_{j+1})$, the cross term between species $1$ ($i^j c_1 e^{\sqrt{2} m x}$) and species $3$ ($(-i)^j c_3 e^{\sqrt{2} m x}$) has $N$-ality $1 + 3 \equiv 0 \pmod 4$ ($j$-independent) and exponent $2\sqrt{2} m x$. At $\theta_1 = \theta_3 = 0$, the Hirota interaction coefficient (`@eq-hirota-Aab`) is completely regular and positive for all three pairs:
$$
A_{ab}(0) = \frac{1 - \cos(\pi(a-b)/4)}{1 - \cos(\pi(a+b)/4)}
\implies A_{12}(0) = A_{23}(0) = \frac{1 - 1/\sqrt{2}}{1 + 1/\sqrt{2}} = (\sqrt{2}-1)^2 > 0, \qquad A_{13}(0) = \frac{1 - 0}{1 - (-1)} = \frac{1}{2} > 0.
$$
Because $A_{12}(0), A_{23}(0), A_{13}(0)$ are all finite and non-zero, the static 3-soliton Hirota formula
$$
\tau_j(x; c_1, c_2, c_3) = 1 + \sum_{a=1}^3 i^{a j} c_a e^{m_a x} + \sum_{1 \le a < b \le 3} A_{ab}(0) i^{(a+b)j} c_a c_b e^{(m_a + m_b)x} + A_{123} (-1)^j c_1 c_2 c_3 e^{(m_1 + m_2 + m_3)x}
$$
(with $A_{123} = A_{12}(0)A_{23}(0)A_{13}(0) = \frac{1}{2}(\sqrt{2}-1)^4$) is an exact static solution for every $(c_1, c_2, c_3) \in \mathbb{C}^3$ whose asymptotic expansion as $x \to -\infty$ has the arbitrary leading coefficients $(c_1, c_2, c_3)$ of $e^{m_a x}$. By Poincaré–Dulac uniqueness, this 3-parameter Hirota family exhausts the unstable manifold $\mathcal{W}^u(0)$ of the vacuum as $x \to -\infty$.

### 1.2 Analytical proof that only $S = \{2\}$ and $S = \{1,2,3\}$ have N-ality $2 \pmod 4$, and that $S = \{1,2,3\}$ reaches all 6 weights of $\Lambda^2\mathbb{C}^4$

Let $S = \{a \in \{1, 2, 3\} : c_a \ne 0\}$ be the non-empty set of active species of a regular static trajectory in $\mathcal{W}^u(0)$ (where no $\tau_j(x)$ vanishes on $x \in \mathbb{R}$). Since $A_{ab}(0) > 0$ for all $1 \le a < b \le 3$ and $m_a > 0$ for all $a \in S$, the single dominant term in $\tau_j(x; c_1, c_2, c_3)$ as $x \to +\infty$ is the full-subset product over $S$:
$$
\tau_j(x) \sim \left(\prod_{a < b \in S} A_{ab}(0)\right)\left(\prod_{a \in S} c_a\right) i^{\left(\sum_{a\in S} a\right) j} \exp\left(\sum_{a\in S} m_a x\right) \qquad (x \to +\infty).
$$
From this single asymptotic relation, two exact consequences follow:
1. **Energy:** By `@prp-hirota-energy`, the classical energy is determined by the leading exponent $\sum_{a\in S} m_a$ as $x \to +\infty$ and equals $\sum_{a\in S} M_a$.
2. **N-ality selection rule:** As $x \to +\infty$, $\tau_j(+\infty)/\tau_0(+\infty) = i^{(\sum_{a\in S} a) j}$, so $\Delta\arg\tau_j - \Delta\arg\tau_0 \equiv \frac{2\pi j}{4}\sum_{a\in S} a \pmod{2\pi}$. Since $\vec Q = -\frac{1}{2\pi}\sum_{j=0}^3 \vec\alpha_j \Delta\arg\tau_j$ and $-\sum_{j=1}^3 \frac{j}{4}\vec\alpha_j = \varpi_1$ (the fundamental weight of $\Lambda^1\mathbb{C}^4$), the topological charge $\vec Q$ lies in the weight-lattice coset $\varpi_k + Q(A_3)$ with N-ality
   $$
   k \equiv \sum_{a \in S} a \pmod 4.
   $$
   Evaluating $\sum_{a\in S} a \bmod 4$ on all $2^3 - 1 = 7$ non-empty subsets $S \subseteq \{1, 2, 3\}$:
   - $S = \{1\}, \{2,3\}$: $\sum a \equiv 1 \pmod 4$ (coset of $\Lambda^1\mathbb{C}^4$),
   - $S = \{3\}, \{1,2\}$: $\sum a \equiv 3 \pmod 4$ (coset of $\Lambda^3\mathbb{C}^4$),
   - $S = \{1,3\}$: $\sum a \equiv 0 \pmod 4$ (root lattice $Q(A_3)$),
   - $S = \{2\}$ and $S = \{1,2,3\}$: $\sum a \equiv 2 \pmod 4$ (coset of $\Lambda^2\mathbb{C}^4$).

   Thus **no subset other than $S = \{2\}$ (energy $M_2$) and $S = \{1,2,3\}$ (energy $M_1 + M_2 + M_3 = (1+\sqrt{2})M_2$) can ever produce a charge in $\Lambda^2\mathbb{C}^4$**.

3. **Charges reached by $S = \{2\}$ and $S = \{1,2,3\}$:**
   - For $S = \{2\}$, $\tau_j(x) = 1 + (-1)^j c_2 e^{2mx}$ has $\tau_0(x) = \tau_2(x)$ and $\tau_1(x) = \tau_3(x)$ for all $x \in \mathbb{R}$. For $c_2 = |c_2| e^{i\psi}$ ($\psi \in (-\pi, \pi)\setminus\{0\}$), $\Delta\arg\tau_0 = \Delta\arg\tau_2 = \psi$ and $\Delta\arg\tau_1 = \Delta\arg\tau_3 = \psi - \pi\operatorname{sgn}(\psi)$, so
     $$
     \vec Q = \frac{\operatorname{sgn}(\psi)}{2}(\vec\alpha_1 + \vec\alpha_3) = \operatorname{sgn}(\psi)\frac{1}{2}(1, -1, 1, -1) \in \big\{\vec\lambda_{(1,3)},\ \vec\lambda_{(2,4)}\big\}.
     $$
   - For $S = \{1,2,3\}$, write $c_a = e^{m_a x_a + i\phi_a}$. Choosing the real centers $x_1 \ll x_2 \ll x_3$ (for example $x_1 = -L, x_2 = 0, x_3 = +L$ with $L \gg 1$), the solution decomposes on $[-2L, 2L]$ into three well-separated single-species kinks: a species-1 kink near $x_1$ (with effective phase $e^{i\phi_1}$), a species-2 kink near $x_2$ (with effective phase $A_{12}(0)e^{i\phi_2}$), and a species-3 kink near $x_3$ (with effective phase $A_{13}(0)A_{23}(0)e^{i\phi_3}$). Since $A_{ab}(0) > 0$ are positive real numbers, they do not shift the phases $\phi_a$. A single species-1 kink reaches all 4 weights $e_a - \frac{1}{4}\mathbf{1}$ of $\Lambda^1$, a single species-2 kink reaches $e_1 + e_3 - \frac{1}{2}\mathbf{1}$ and $e_2 + e_4 - \frac{1}{2}\mathbf{1}$ of $\Lambda^2$, and a single species-3 kink reaches all 4 weights $-e_b + \frac{1}{4}\mathbf{1}$ of $\Lambda^3$. By additivity of $\Delta\arg\tau_j$ across the three separated regions, $\vec Q = \vec Q_1 + \vec Q_2 + \vec Q_3$ realizes every weight of the form $(e_1 + e_3 - \frac{1}{2}\mathbf{1}) + e_a - e_b$ ($b \in \{1,3\}, a \in \{2,4\}$), which covers all four missing weights $\{(1,2), (2,3), (3,4), (1,4)\}$ as well as $\{(1,3), (2,4)\}$. (Verified numerically in Check 1 of `verify_missing_charge.py`.)

---

## Item 2: Exact $(1+1)$-dimensional fusing-angle bound states at energy $M_2$

**Verdict: Correct; the second `[Numerical]` claim in `prp-static-saddles` can also be upgraded to `**[Theorem]**` via the exact asymmetric fusing-angle solution $\tau_j^\pm(x, t)$.**

### 2.1 Derivation of the asymmetric fusing limit

For two species-1 solitons of $a_3^{(1)}$ at rapidities $\theta_1 = +i u$ and $\theta_2 = -i u$,
$$
\tau_j(x, t) = 1 + i^j (E_1 + E_2) + (-1)^j A_{11}(2iu) E_1 E_2, \qquad E_{1,2} = e^{\xi_{1,2}} \exp\big(\sqrt{2} m (x\cos u \mp i t\sin u)\big),
$$
the interaction coefficient $A_{11}(2iu) = 1 - \frac{1}{\cos(2u)}$ has a simple pole $A_{11}(2iu) \sim \frac{1}{2\delta}$ as $u = \frac{\pi}{4} + \delta \to \frac{\pi}{4}$ (fusing angle $U = 2u = \frac{\pi}{2}$), where $\sqrt{2}m\cos(\pi/4) = \sqrt{2}m\sin(\pi/4) = m = m_2/2$.
- Taking the symmetric limit $e^{\xi_1}, e^{\xi_2} = O(\sqrt{\delta}) \to 0$ with $A_{11}(2iu) e^{\xi_1+\xi_2} \to C_4$ kills the $i^j(E_1 + E_2)$ term and reduces to the static species-2 soliton $1 + (-1)^j C_4 e^{2mx}$.
- Taking instead the **asymmetric fusing limit** where one constituent amplitude remains finite ($e^{\xi_1} \to C_1 \in \mathbb{C}$) while the other scales as $e^{\xi_2} \sim \frac{2 C_4}{C_1}\delta \to 0$ (or vice versa, exchanging $E_1 \leftrightarrow E_2$) yields the two-parameter family of functions
  $$
  \tau_j^\pm(x, t) = 1 + i^j C_1 e^{m(x \mp i t)} + (-1)^j C_4 e^{2 m x} \qquad (j \in \mathbb{Z}_4).
  $$

### 2.2 Proof that $\tau_j^\pm(x,t)$ is an exact solution with $E = M_2, P = 0$ covering all of $\Lambda^2\mathbb{C}^4$

1. **Exact Hirota solution:** Set $A = C_1 e^{m(x \mp i t)}$ and $B = C_4 e^{2mx}$. Since $(\partial_x^2 - \partial_t^2)A = 2m^2 A$ and $(\partial_x A)^2 - (\partial_t A)^2 = m^2 A^2 - (-m^2 A^2) = 2m^2 A^2$, the $A^2$ terms cancel on the LHS of the Hirota bilinear equation (`@eq-hirota-bilinear`):
   $$
   \tau_j^\pm(\partial_x^2 - \partial_t^2)\tau_j^\pm - \big((\partial_x\tau_j^\pm)^2 - (\partial_t\tau_j^\pm)^2\big)
   = m^2\big[2 i^j A + 4(-1)^j B + 2 i^j(-1)^j A B\big].
   $$
   On the RHS, using $i^{j-1} + i^{j+1} = 0$, $i^{j-1}i^{j+1} = (-1)^j$, and $(-1)^{j\pm 1} = -(-1)^j$:
   $$
   (\tau_j^\pm)^2 - \tau_{j-1}^\pm \tau_{j+1}^\pm = 2 i^j A + 4(-1)^j B + 2 i^j(-1)^j A B.
   $$
   Thus `@eq-hirota-bilinear` holds identically for all $(C_1, C_4) \in \mathbb{C}^2$. (In `test_a3_thimbles.py`, the maximum numerical residual of `@eq-hirota-bilinear` and of $\vec u_{xx} - \vec u_{tt} = \nabla W(\vec u)$ is $2.4\times 10^{-13}$.)
2. **Exact rest energy $E = M_2$ and momentum $P = 0$:** On any regular trajectory (where no $\tau_j^\pm$ vanishes on $x \in \mathbb{R}$), $A = C_1 e^{m(x\mp it)}$ is subleading to $B = C_4 e^{2mx}$ as $x \to +\infty$ ($\partial_x \log\tau_j^\pm \to 2m$, $\partial_t\log\tau_j^\pm = O(e^{-mx}) \to 0$) and vanishes as $x \to -\infty$. By `@prp-hirota-energy`, the total energy and momentum depend only on the $x \to \pm\infty$ asymptotics of $\partial_\mu \log\tau_j^\pm$, so $E = M_2$ and $P = 0$ identically for all $(C_1, C_4)$.
3. **Closed-form coverage of all 6 weights of $\Lambda^2\mathbb{C}^4$ at $t = 0$:** Write $r = e^{mx} \in (0, \infty)$ at $t = 0$, so $\tau_j(r) = 1 + i^j C_1 r + (-1)^j C_4 r^2$:
   - **Reachable weights $(1,3)$ and $(2,4)$:** At $C_1 = 0, C_4 = e^{i\psi}$ ($\psi \in (-\pi,\pi)\setminus\{0\}$), $\vec Q = \operatorname{sgn}(\psi)\frac{1}{2}(1,-1,1,-1)$, giving $\vec\lambda_{(1,3)}$ for $\psi > 0$ and $\vec\lambda_{(2,4)}$ for $\psi < 0$.
   - **Missing weights $(1,2)$ and $(3,4)$:** Set $C_4 = -1$ and $C_1 = i\epsilon$ with $\epsilon \in (-2, 2)\setminus\{0\}$. Then for all $r \in (0, \infty)$:
     $$
     \tau_1(r) = 1 - \epsilon r + r^2 = \left(r - \frac{\epsilon}{2}\right)^2 + 1 - \frac{\epsilon^2}{4} > 0, \qquad
     \tau_3(r) = 1 + \epsilon r + r^2 = \left(r + \frac{\epsilon}{2}\right)^2 + 1 - \frac{\epsilon^2}{4} > 0,
     $$
     so $\tau_1, \tau_3$ are strictly positive real and $\Delta\arg\tau_1 = \Delta\arg\tau_3 = 0$. Meanwhile, $\tau_0(r) = (1 - r^2) + i\epsilon r$ and $\tau_2(r) = (1 - r^2) - i\epsilon r = \overline{\tau_0(r)}$ have non-vanishing imaginary parts $\pm\epsilon r \ne 0$ for all $r \in (0, \infty)$ while their real part $1 - r^2$ goes from $+1$ to $-\infty$, so $\Delta\arg\tau_0 = -\Delta\arg\tau_2 = \pi\operatorname{sgn}(\epsilon)$. Therefore:
     $$
     \vec Q = -\frac{1}{2\pi}\sum_{j=0}^3 \vec\alpha_j \Delta\arg\tau_j = \frac{\operatorname{sgn}(\epsilon)}{2}(\vec\alpha_2 - \vec\alpha_0) = \operatorname{sgn}(\epsilon)\frac{1}{2}(1, 1, -1, -1) = \begin{cases} \vec\lambda_{(1,2)} & (0 < \epsilon < 2), \\ \vec\lambda_{(3,4)} & (-2 < \epsilon < 0). \end{cases}
     $$
   - **Missing weights $(2,3)$ and $(1,4)$:** Applying the diagram rotation $j \mapsto j+1 \pmod 4$ ($C_1 \mapsto i C_1, C_4 \mapsto -C_4$), i.e. $C_4 = +1$ and $C_1 = \epsilon \in (-2, 2)\setminus\{0\}$, makes $\tau_0(r) = 1 + \epsilon r + r^2 > 0$ and $\tau_2(r) = 1 - \epsilon r + r^2 > 0$ strictly positive real ($\Delta\arg\tau_0 = \Delta\arg\tau_2 = 0$), while $\tau_1(r) = (1 - r^2) + i\epsilon r$ and $\tau_3(r) = \overline{\tau_1(r)}$ wind by $\Delta\arg\tau_1 = -\Delta\arg\tau_3 = \pi\operatorname{sgn}(\epsilon)$. Therefore:
     $$
     \vec Q = \frac{\operatorname{sgn}(\epsilon)}{2}(\vec\alpha_3 - \vec\alpha_1) = \operatorname{sgn}(\epsilon)\frac{1}{2}(-1, 1, 1, -1) = \begin{cases} \vec\lambda_{(2,3)} & (0 < \epsilon < 2), \\ \vec\lambda_{(1,4)} & (-2 < \epsilon < 0). \end{cases}
     $$
4. **Euclidean action and the $\epsilon \to 0^+$ limit:** In Euclidean time $\tau = i t$, $e^{m(x\mp it)} = e^{m(x\mp\tau)}$, so $\epsilon(\tau) = \epsilon e^{\mp m\tau}$ preserves its sign and purely imaginary $C_1 = i\epsilon(\tau)$ character for all $\tau$ in the slab where $0 < |\epsilon| e^{\mp m\tau} < 2$. In this regime, $\tau_2 = \overline{\tau_0}$ and $\tau_1, \tau_3 \in \mathbb{R}_+$, so $\vec u(x, \tau)$ is invariant under the antilinear involution $\sigma\circ\Theta$ ($\sigma = (1\leftrightarrow 2)(3\leftrightarrow 4)$) and its Euclidean Lagrangian $L_{\text{Eucl}}(\tau) = \int_{-\infty}^\infty dx \left[\frac{1}{2}(\partial_\tau\vec u)^2 + \frac{1}{2}(\partial_x\vec u)^2 + W(\vec u)\right]$ is strictly real. As $\epsilon \to 0^+$, the root of $\tau_0(x) = 1 - e^{2mx} + i\epsilon e^{mx}$ approaches $x = 0$ from the upper half-plane at $x_* = \frac{i\epsilon}{2m} + O(\epsilon^3)$. Near $x = 0$:
   - $\partial_x \log\tau_0(x) \approx \frac{1}{x - i\epsilon/(2m)}$ has a **simple pole** with residue $+1$, whose imaginary part integrates to $+\pi$, producing the $\vec\lambda_{(1,2)}$ topological charge for every $\epsilon > 0$;
   - $(\partial_\tau \log\tau_0(x))^2 \approx -\frac{\epsilon^2/4}{(x - i\epsilon/(2m))^2}$ has a **double pole** in the upper half-plane with **residue zero**, so its singular peak integrates to zero by Cauchy's theorem, leaving $\int_{-\infty}^\infty (\partial_\tau\vec u)^2 dx = O(\epsilon^2) \to 0$ and $\lim_{\epsilon\to 0^+} L_{\text{Eucl}} = M_2$.

   Numerical output from `test_fusing_action.py` (with $m = 1, M_2 = 16$):
   ```text
   === Testing C4 = -1.0, C1 = i * eps (missing charge 12, strictly real L_Eucl!) ===
   eps= 1.90 [charge 12]: int(u_x^2/2)= -6.14525 | int(u_t^2/2)=204.98553 | int(W)=-182.84028 | E_Mink= 16.00000 | L_Eucl=-393.97105 (Im=3.4e-15)
   eps= 1.50 [charge 12]: int(u_x^2/2)=  3.45253 | int(u_t^2/2)= 16.13246 | int(W)= -3.58499 | E_Mink= 16.00000 | L_Eucl=-16.26492 (Im=1.4e-14)
   eps= 1.00 [charge 12]: int(u_x^2/2)=  5.58160 | int(u_t^2/2)=  3.47280 | int(W)=  6.94560 | E_Mink= 16.00000 | L_Eucl=  9.05440 (Im=1.5e-14)
   eps= 0.50 [charge 12]: int(u_x^2/2)=  6.74243 | int(u_t^2/2)=  0.61717 | int(W)=  8.64040 | E_Mink= 16.00000 | L_Eucl= 14.76566 (Im=1.5e-14)
   eps= 0.20 [charge 12]: int(u_x^2/2)=  7.43781 | int(u_t^2/2)=  0.08649 | int(W)=  8.47570 | E_Mink= 16.00000 | L_Eucl= 15.82703 (Im=1.5e-14)
   eps= 0.10 [charge 12]: int(u_x^2/2)=  7.70400 | int(u_t^2/2)=  0.02079 | int(W)=  8.27521 | E_Mink= 16.00000 | L_Eucl= 15.95842 (Im=1.5e-14)
   eps= 0.05 [charge 12]: int(u_x^2/2)=  7.84768 | int(u_t^2/2)=  0.00510 | int(W)=  8.14722 | E_Mink= 16.00000 | L_Eucl= 15.98980 (Im=1.2e-14)
   ```

5. **Generalization to all $a_n^{(1)}$ ($h = n+1$) and $\Lambda^a\mathbb{C}^h$:** Fusing $a$ species-1 solitons at the static Bethe-string rapidities $\theta_k = \frac{i\pi}{h}(a + 1 - 2k)$ ($k = 1, \dots, a$, where adjacent rapidities differ by the fundamental fusing angle $\frac{2\pi i}{h}$) and taking the asymmetric fusing limit produces the exact $a$-parameter Hirota solution
   $$
   \tau_j(x, t) = 1 + \sum_{b=1}^{a-1} \omega^{b j} C_b e^{k_b x - i\omega_b t} + \omega^{a j} C_a e^{m_a x} \qquad (\omega = e^{2\pi i/h},\ k_b = m_b\cos\tfrac{\pi(a-b)}{h},\ \omega_b = \pm m_b\sin\tfrac{\pi(a-b)}{h}),
   $$
   which has exact rest energy $E = M_a$, momentum $P = 0$, and reaches all $\binom{h}{a}$ weights of $\Lambda^a\mathbb{C}^h$ (verified for $a_4^{(1)}, h=5, \Lambda^2\mathbb{C}^5$ in Check 3 of `verify_missing_charge.py`, where all $\binom{5}{2} = 10$ weights are reached with bilinear residual $< 10^{-12}$).

---

## Item 3: Finite-lattice theorem for topological sectors on the $(1+1)$-dimensional lattice

**Verdict: New `[Theorem]` proving the existence of complex saddles with $n_\sigma \ne 0$ in every missing-charge sector on any finite lattice.**

Consider the $(1+1)$-dimensional lattice $\Lambda = \mathbb{Z}_{N_t} \times \{1, \dots, L\}$ of `part9-foundations/45-lattice.qmd` (`@eq-lattice-action`) in topological sector $\vec\lambda \in P(g)$, defined by the spatial boundary conditions $\vec u(t, 0) = 0$ and $\vec u(t, L+1) = 2\pi\vec\lambda$ for all $t \in \mathbb{Z}_{N_t}$ (with $N_t = 1$ recovering the reduced 1D static chain of `sec-which-saddles-contribute`).

### 3.1 Proof of reality, diagram-orbit equality, and non-vanishing of $n_\sigma$

1. **Reality in every topological sector $\vec\lambda \in P(g)$:**
   Even when $\vec\lambda$ is not $\Theta$-invariant (so $\Theta: \vec u \mapsto -\bar{\vec u}$ does not preserve the boundary condition at $x = L+1$), combine $\Theta$ with spatial reflection $R_x: x \mapsto L + 1 - x$:
   $$
   (\Theta \circ R_x)(\vec u)(t, x) := 2\pi\vec\lambda - \overline{\vec u(t, L+1-x)}.
   $$
   Then $(\Theta \circ R_x)(\vec u)(t, 0) = 2\pi\vec\lambda - 2\pi\vec\lambda = 0$ and $(\Theta \circ R_x)(\vec u)(t, L+1) = 2\pi\vec\lambda - 0 = 2\pi\vec\lambda$. Since $\vec\alpha_j \cdot \vec\lambda \in \mathbb{Z}$ for every weight $\vec\lambda \in P(g)$ and every affine root $\vec\alpha_j$, we have $e^{i\vec\alpha_j \cdot (2\pi\vec\lambda)} = 1$, so
   $$
   e^{i\vec\alpha_j \cdot (\Theta \circ R_x)(\vec u)(t, x)} = e^{-i\vec\alpha_j \cdot \overline{\vec u(t, L+1-x)}} = \overline{e^{i\vec\alpha_j \cdot \vec u(t, L+1-x)}}.
   $$
   Summing over $x \in \{1, \dots, L\}$ and $t \in \mathbb{Z}_{N_t}$ gives $S_\Lambda[(\Theta \circ R_x)(\vec u)] = \overline{S_\Lambda[\vec u]}$. Since $\Theta \circ R_x$ is a real-volume-preserving involution on the real integration cycle $\Gamma_{\mathbb{R}} = \mathbb{R}^{r N_t L}$, the partition function $Z_\Lambda(\vec\lambda)$ is **strictly real** for all $\vec\lambda \in P(g), N_t, L \ge 1$, and $\kappa, \hbar > 0$.
   Equivalently, splitting $\vec u(t, x) = \frac{x}{L+1}(2\pi\vec\lambda) + \vec v(t, x)$ (where the harmonic background $\vec u_0(x) = \frac{x}{L+1}(2\pi\vec\lambda)$ is orthogonal to the Dirichlet fluctuation $\vec v$) and expanding in the fugacity $\kappa/\hbar$ as in `@thm-coulomb-gas` gives the exact sector Coulomb-gas series:
   $$
   Z_\Lambda(\vec\lambda) = Z_{\Lambda,0}\, e^{-\frac{2\pi^2 N_t |\vec\lambda|^2}{(L+1)\hbar} - \frac{\kappa|\Lambda|h}{\hbar}} \sum_{\{m_{u,j}\ge 0\}} \left[\prod_{u,j} \frac{(\kappa k_j/\hbar)^{m_{u,j}}}{m_{u,j}!}\right] \cos\!\left(\sum_{t,x,j} m_{(t,x),j} \frac{2\pi x}{L+1}(\vec\alpha_j \cdot \vec\lambda)\right) e^{-\frac{\hbar}{2}\sum_{u,u'} \vec Q_u \cdot \vec Q_{u'} C_\Lambda(u,u')}.
   $$
2. **Equality across diagram-automorphism orbits and non-vanishing of $n_\sigma$:**
   Any automorphism $\pi \in \operatorname{Aut}(\hat{\Gamma})$ of the affine Dynkin diagram permutes the simple roots $\{\vec\alpha_j\}$ and Kac labels $k_j$ isometrically, so $Z_\Lambda(\pi(\vec\lambda)) = Z_\Lambda(\vec\lambda)$. For $a_3^{(1)}$, the dihedral group $D_4 = \operatorname{Aut}(\hat{\Gamma})$ partitions the 6 weights of $\Lambda^2\mathbb{C}^4$ into the reachable orbit $\mathcal{O}_{\text{reach}} = \{(1,3), (2,4)\}$ and the missing-charge orbit $\mathcal{O}_{\text{miss}} = \{(1,2), (2,3), (3,4), (1,4)\}$, so
   $$
   Z_\Lambda((1,2)) = Z_\Lambda((2,3)) = Z_\Lambda((3,4)) = Z_\Lambda((1,4)), \qquad Z_\Lambda((1,3)) = Z_\Lambda((2,4)).
   $$
   Since the Coulomb-gas power series in $\kappa$ converges absolutely on $\mathbb{C}$ (hence is real-analytic on $\kappa \in (0, \infty)$) and its $m_{u,j} \equiv 0$ leading term is $1 \ne 0$, $Z_\Lambda(\vec\lambda) \ne 0$ for all but at most a discrete set of $\kappa > 0$. By the Picard–Lefschetz decomposition $Z_\Lambda(\vec\lambda) = \sum_\sigma n_\sigma Z_\sigma$, whenever $Z_\Lambda(\vec\lambda) \ne 0$, **at least one complex saddle $\sigma$ in each missing-charge sector has a nonzero intersection number $n_\sigma \ne 0$**.

---

## Item 4: Explicit Picard–Lefschetz decomposition in the missing sector $(1,2)$ of $a_3^{(1)}$ and continuum decoupling of the 1D static saddle

**Verdict: New `[Numerical]` check and analytical continuum-limit explanation for `sec-which-saddles-contribute`.**

### 4.1 Exact $L=1$ Coulomb-gas $Z_1(\vec\lambda)$ and saddle decomposition in sector $(1,2)$ vs. $(1,3)$

From `test_a3_thimbles.py` and `test_continuum_12.py`, for $a_3^{(1)}$ at $L=1, \kappa=1$:

```text
2. Exact Z_1(lambda) for all 6 weights of Lambda^2 in a_3^(1) (kappa=1.0, hbar=0.5):
   sector (1, 2): Z_1 = 5.599088489751e-15
   sector (1, 3): Z_1 = 7.306189799424e-16
   sector (1, 4): Z_1 = 5.599088489750e-15
   sector (2, 3): Z_1 = 5.599088489750e-15
   sector (2, 4): Z_1 = 7.306189799424e-16
   sector (3, 4): Z_1 = 5.599088489751e-15

=== 1. Sector missing (1,2): L=1, kappa=1.0, min_R Re(S) = 13.869604, Imax = 4.0000 ===
  cp  0: S= 14.86816 +0.00000j  n=-1  #hits=1  mu=[4.941 3.724 3.217]  u1/2pi=[ 0.3804+0.0395j  0.3804-0.0395j -0.3804+0.0395j -0.3804-0.0395j]
  cp  1: S= 14.86816 -0.00000j  n=-1  #hits=1  mu=[4.941 3.724 3.217]  u1/2pi=[ 0.1196+0.0395j  0.1196-0.0395j -0.1196+0.0395j -0.1196-0.0395j]
  cp  2: S= 15.88236 +3.00681j  n=-1  #hits=1  mu=[4.835 3.774 2.498]  u1/2pi=[ 0.1461+0.0504j  0.3539-0.0504j -0.3539+0.0504j -0.1461-0.0504j]
  cp  3: S= 15.88236 -3.00681j  n=-1  #hits=1  mu=[4.835 3.774 2.498]  u1/2pi=[ 0.3539+0.0504j  0.1461-0.0504j -0.1461+0.0504j -0.3539-0.0504j]
  cp  4: S= 22.30199 +3.24420j  n=+0  #hits=0  mu=[7.744 3.925 2.918]  u1/2pi=[-0.1052-0.1139j  0.6052+0.1139j -0.2687+0.0248j -0.2313-0.0248j]
  --- Saddle expansion vs exact Coulomb-gas Z_1 for missing (1,2) ---
  hbar=1.00: Z_exact=6.18242000e-09 | sum n Z_1=+5.93895126e-09 (rel.err 3.94e-02) | sum n Z_2=+6.39896232e-09 (rel.err 3.50e-02)
  hbar=0.50: Z_exact=5.59908849e-15 | sum n Z_1=+5.50821633e-15 (rel.err 1.62e-02) | sum n Z_2=+5.63593779e-15 (rel.err 6.58e-03)
  hbar=0.25: Z_exact=1.96380778e-27 | sum n Z_1=+1.94114341e-27 (rel.err 1.15e-02) | sum n Z_2=+1.96716233e-27 (rel.err 1.71e-03)
  hbar=0.15: Z_exact=2.57791503e-44 | sum n Z_1=+2.54830766e-44 (rel.err 1.15e-02) | sum n Z_2=+2.56905157e-44 (rel.err 3.44e-03)
```

- In the **reachable sector $(1,3)$**, the floor on $\mathbb{R}^3$ is $\min_{\mathbb{R}^3}\operatorname{Re} S = 16.186413$, which is attained at a $\Theta\circ R_x$-pair of **purely real** critical points $u_1/(2\pi) = \gamma_1(1, -1, 1, -1)$ with $\gamma_1 \in \{0.4008, 0.0992\}$ in the single-soliton subspace $\operatorname{span}(\vec\alpha_1 + \vec\alpha_3)$ (the 4 complex critical points with $S = 15.62961 < 16.18641$ lie below the real floor and have $n_\sigma = 0$). In the continuum limit $L \to \infty, a = \ell/(L+1) \to 0$ with $\kappa = a^2$, this real saddle converges smoothly to the static species-2 soliton $\gamma(x) = 2\arctan(e^{2x})$ with continuum action $S_{\text{cont}} = S_{\text{chain}}/a \to 16.0 = M_2$.
- In the **missing sector $(1,2)$**, there are **no real critical points** on $\mathbb{R}^3$. The 4 contributing critical points (`cp 0..3`, all with $|n_\sigma| = 1$) consist of:
  1. A $\Theta\circ R_x$-pair (`cp 0, 1`) with **strictly real action** $S = 14.86816$ lying in the 2-real-dimensional subspace
     $$
     \mathcal{V}_{(1,2)} = \big\{\vec u = \alpha(1, 1, -1, -1) + i\beta(1, -1, 1, -1) : \alpha, \beta \in \mathbb{R}\big\},
     $$
     which is fixed pointwise by the antilinear symmetry $\sigma\circ\Theta$ ($\sigma = (1\leftrightarrow 2)(3\leftrightarrow 4)$). Because $(1,1,-1,-1) \perp (1,-1,1,-1)$, the action on $\mathcal{V}_{(1,2)}^L$ is identically real for any $L$:
     $$
     S[\alpha, \beta] = \sum_{x=0}^L \big[2(\alpha_{x+1}-\alpha_x)^2 - 2(\beta_{x+1}-\beta_x)^2\big] + \kappa \sum_{x=1}^L \big[4 - 2e^{2\beta_x}\cos(2\alpha_x) - 2e^{-2\beta_x}\big].
     $$
  2. A complex-action conjugate pair (`cp 2, 3`, $S = 15.88236 \pm 3.00681i$) lying in the reflection-invariant subspace $\mathcal{W}_{(1,2)} = \{u_4 = -u_1, u_3 = -u_2\}$.

### 4.2 Why the 1D static lattice saddle on $\mathcal{V}_{(1,2)}$ is a link-pinned lattice artifact that decouples as $a \to 0$

Why does the 1D chain have a low-action saddle in $\mathcal{V}_{(1,2)}$ at $L=1,2,3$ when `prp-static-saddles` proves that no smooth continuum static saddle of energy $M_2$ exists in sector $(1,2)$?

In Hirota variables on $\mathcal{V}_{(1,2)}$, $\tau_2(x) = \overline{\tau_0(x)} = |\tau_0(x)|e^{-i\alpha(x)}$ and $\tau_1(x) = \tau_3(x) \in \mathbb{R}$ with $e^{2\beta(x)} = |\tau_1(x)|/|\tau_0(x)|$. Any static continuum solution in sector $(1,2)$ must go from $\tau_j(-\infty) = +1$ to $\tau_j(+\infty) \sim (-1)^j c_4 e^{2x}$ ($c_4 > 0$), so $\tau_1(x) = \tau_3(x) \in \mathbb{R}$ starts at $+1 > 0$ and ends at $-\infty < 0$. By the intermediate value theorem, $\tau_1(x)$ **must vanish at some real $x_* \in \mathbb{R}$**, where $\beta(x) \sim \frac{1}{2}\log|x - x_*| \to -\infty$. (This is confirmed directly by 1D shooting of the continuum ODEs $\alpha'' = e^{2\beta}\sin(2\alpha), \beta'' = e^{2\beta}\cos(2\alpha) - e^{-2\beta}$ in `test_shooting_12.py`: every trajectory on the stable manifold of $(\alpha, \beta) = (\pi, 0)$ arrives at $\alpha = \pi/2$ with strictly positive derivative $\beta' > 0$, and every symmetric trajectory with $\alpha(0) = \pi/2, \beta'(0) = 0$ hits a logarithmic singularity $\beta(x_*) \to -\infty$ at a finite $x_* \in (1.0, 2.3)$.)

On a discrete chain with lattice spacing $a > 0$, two adjacent lattice sites $x_k, x_{k+1}$ straddle the singularity $x_*$ without landing on it (in `test_v12_bvp.py`, $\alpha_x/(2\pi)$ jumps from $0.10$ to $0.40$ across a single bond at $L=1, 2, 3$). Across that bond, $\Delta\alpha = O(1)$, so the discrete chain action $S_{\text{chain}} = \sum \frac{1}{2}(\Delta u)^2 + a^2 \sum W(u)$ stays $O(1)$ ($\approx 4.93$) as $a \to 0$, and the physical continuum action
$$
S_{\text{cont}} = \frac{S_{\text{chain}}}{a} \sim \frac{4.93}{a} \to +\infty \qquad (a \to 0)
$$
diverges like $1/a$. Thus the 1D link-pinned saddle decouples as $e^{-c/(a\hbar)} \to 0$ in the continuum limit, whereas turning on the $(1+1)$-dimensional fusing-angle mode $C_1 = i\epsilon$ in Item 2 shifts the zero off the real axis to $x = \pm i\epsilon/(2m)$ and yields a smooth continuum saddle of action $M_2$.

---

## Item 5: Recommended edits to the book

### 5.1 `part5-semiclassics/17-thimbles.qmd` (lines 120–138)

1. **Upgrade `prp-static-saddles` (lines 120–124) from `**[Proposition, sketch]**` to `**[Theorem]**` and include its proof:**

   Replace lines 120–126 of `part5-semiclassics/17-thimbles.qmd`:
   ```markdown
   ::: {#prp-static-saddles}
   ## Static saddles cannot fill the multiplets
   **[Theorem]** A single static soliton of species $a$ of $a_n^{(1)}$ carries $h/\gcd(a,h)$ topological charges [@mcghee1994]; scanning $\operatorname{Im}\xi$ confirms this for $n\le4$. So the species-2 soliton of $a_3^{(1)}$ carries only the weights $13$ and $24$ of the six weights of $\Lambda^2$. In $a_3^{(1)}$ the spatial exponents $(m_1,m_2,m_3)=(\sqrt2,2,\sqrt2)m$ have no Poincaré–Dulac resonances of degree at least $2$. Every complex static solution leaving the vacuum therefore lies in the unstable manifold of the vacuum, which is parametrized by the coefficients of the $e^{m_ax}$ and is realized by static Hirota multi-solitons of distinct species $S\subseteq\{1,2,3\}$. Its energy is $\sum_{a\in S}M_a$ and its topological charge has $\mathbb Z_4$ N-ality $\sum_{a\in S}a\bmod 4$. The four missing weights of $\Lambda^2$ are reached statically only with $S=\{1,2,3\}$, at energy $(1+\sqrt2)M_2$. So no static saddle of energy $M_2$ carries them. Bound states of two species-1 solitons at rapidities $\pm iu$, by contrast, have energy $2M_1\cos u$, which tends to $M_2$ at $u=\pi/4$, so that the fusing angle is $U=2u=\pi/2$. In the asymmetric fusing limit as $u\to\pi/4$, they yield the exact solutions
   $$
   \tau_j^\pm(x,t)=1+i^jC_1e^{m(x\mp it)}+(-1)^jC_4e^{2mx}
   $$ {#eq-fusing-angle-soliton}
   of energy $M_2$ and momentum $0$, whose charges at $t=0$ cover all six weights of $\Lambda^2$, consistently with [@harder1995].
   :::

   ::: {.proof}
   For $a_3^{(1)}$, $(n_1+n_3)\sqrt2+2n_2\notin\{\sqrt2,2\}$ for all non-negative integers with $n_1+n_2+n_3\ge2$. By the Poincaré–Dulac theorem, every static trajectory leaving $\vec u=0$ as $x\to-\infty$ is uniquely determined by the coefficients $(c_1,c_2,c_3)\in\mathbb C^3$ of $e^{m_ax}$. Since $A_{12}(0)=A_{23}(0)=(\sqrt2-1)^2>0$ and $A_{13}(0)=\frac12>0$ are finite and non-zero, the static $\{1,2,3\}$ Hirota family is an exact solution for every $(c_1,c_2,c_3)\in\mathbb C^3$ and exhausts the unstable manifold. For any non-empty active subset $S=\{a:c_a\ne0\}$, the dominant term as $x\to+\infty$ has exponent $\sum_{a\in S}m_ax$ and phase $i^{(\sum_{a\in S}a)j}$, giving energy $\sum_{a\in S}M_a$ by @prp-hirota-energy and N-ality $\sum_{a\in S}a\bmod 4$. Among the seven non-empty subsets of $\{1,2,3\}$, only $S=\{2\}$ and $S=\{1,2,3\}$ have $\sum_{a\in S}a\equiv2\pmod4$. For $S=\{2\}$, $\tau_0=\tau_2$ and $\tau_1=\tau_3$, which gives only $\vec\lambda_{(1,3)}$ and $\vec\lambda_{(2,4)}$. For $S=\{1,2,3\}$, choosing the real centres $x_1\ll x_2\ll x_3$ separates the three kinks, and additivity of topological charge $\vec Q=\vec Q_1+\vec Q_2+\vec Q_3$ reaches all six weights of $\Lambda^2$.

   For the two-soliton bound state at $\theta_{1,2}=\pm iu$, $A_{11}(2iu)=1-1/\cos(2u)$ has a simple pole at $u=\pi/4$. Taking $e^{\xi_1}\to C_1$ and $e^{\xi_2}\sim 2(u-\pi/4)C_4/C_1\to0$ gives (-@eq-fusing-angle-soliton), which satisfies @eq-hirota-bilinear identically and has energy $M_2$ and momentum $0$ by @prp-hirota-energy because $e^{m(x\mp it)}$ is subleading to $e^{2mx}$ as $x\to+\infty$. At $t=0$ with $r=e^{mx}>0$: setting $C_1=0, C_4=e^{\pm i\psi}$ ($\psi\in(0,\pi)$) gives $(1,3)$ and $(2,4)$; setting $C_4=-1, C_1=\pm i\epsilon$ ($0<\epsilon<2$) keeps $\tau_1,\tau_3>0$ while $\Delta\arg\tau_0=-\Delta\arg\tau_2=\pm\pi$, giving $(1,2)$ and $(3,4)$; and rotating by $j\mapsto j+1$ ($C_4=+1, C_1=\pm\epsilon$) gives $(2,3)$ and $(1,4)$.
   :::
   ```

2. **Expand the reduced-model bullets (lines 135–137 of `part5-semiclassics/17-thimbles.qmd`) to record the finite-lattice theorem and the $a_3^{(1)}$ missing-sector $(1,2)$ check:**
   - Add a bullet noting that by $\Theta\circ R_x$ symmetry and `@thm-coulomb-gas`, $Z_\Lambda(\vec\lambda)$ on any $(1+1)$-dimensional lattice is real, invariant across diagram-automorphism orbits, and non-zero for generic $(\kappa,\hbar)$, so every missing-charge sector has at least one saddle with $\nu_\sigma\ne 0$.
   - In the final bullet (line 137), note that for $a_3^{(1)}$ at $L=1, \kappa=1$ the missing sector $(1,2)$ is reproduced by four complex saddles with $\lvert\nu_\sigma\rvert=1$, but as $a\to 0$ the real-action 1D saddle in $(1,2)$ is a link-pinned lattice artifact straddling a real zero of $\tau_j$ whose continuum action diverges as $O(1/a)$, whereas turning on $C_1=i\epsilon$ in (-@eq-fusing-angle-soliton) shifts the zero off the real axis with Euclidean Lagrangian $M_2-O(\epsilon^2)\to M_2$.

### 5.2 `appendices/E-original-results.qmd` (line 12)

Whenever `prp-static-saddles` is upgraded from `**[Proposition, sketch]**` to `**[Theorem]**`, update line 12 of `appendices/E-original-results.qmd` in the same commit so `check-refs.py` passes:
```markdown
- **[Theorem]** Picard–Lefschetz theory applied to the saddles of affine Toda theory. Static soliton saddles of energy $M_2$ cannot fill the quantum-group multiplets of topological charges, whereas the asymmetric fusing-angle limit of two-soliton bound states gives exact solutions of energy $M_2$ covering the full multiplet (@prp-static-saddles).
```
*(Reminder from `AGENTS.md`: any upgrade to `[Theorem]` made without an independent referee pass must be listed in the PR description.)*
